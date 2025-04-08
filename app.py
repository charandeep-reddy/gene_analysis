from flask import Flask, render_template, request, jsonify, flash, make_response
from Bio import SeqIO, Align
from Bio.Seq import Seq
import numpy as np
from collections import Counter
import io
import re
import json
import logging
from pprint import pformat
import traceback
import atexit
import multiprocessing
import os
import signal
import sys
from werkzeug.utils import secure_filename

# Setup basic logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Import our custom modules - with proper error handling for missing modules
try:
    from alignment import needleman_wunsch, smith_waterman
    from pattern_search import kmp_pattern_search
    from gene_network import create_gene_network, find_shortest_path
    from codon_optimization import calculate_cai
except ImportError as e:
    logger.error(f"Failed to import required module: {e}")
    logger.error("Make sure all required modules are in the same directory as app.py")
    sys.exit(1)

# Constants
MAX_CONTENT_LENGTH = 5 * 1024 * 1024  # 5MB file size limit
ALLOWED_EXTENSIONS = {'txt', 'fasta', 'fa', 'fna', 'ffn', 'faa', 'frn', 'sif', 'csv', 'tsv'}

app = Flask(__name__)
app.secret_key = os.environ.get('SECRET_KEY', 'dev-secret-key-change-in-production')
app.config['MAX_CONTENT_LENGTH'] = MAX_CONTENT_LENGTH
app.config['UPLOAD_FOLDER'] = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'tmp_uploads')

# Ensure upload folder exists
os.makedirs(app.config['UPLOAD_FOLDER'], exist_ok=True)

def allowed_file(filename):
    """Check if file extension is allowed"""
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

def sanitize_input(text):
    """Basic input sanitization"""
    if not text:
        return ""
    # Remove any non-alphanumeric characters except basic punctuation and biological sequence chars
    return re.sub(r'[^a-zA-Z0-9\s,._\-ACGTURYKMSWBDHVN*]', '', text)

@app.template_filter('pprint')
def pprint_filter(value):
    return pformat(value)

def parse_fasta_content(content):
    """Parse FASTA content and return the sequence."""
    try:
        if not content or not content.strip():
            return ""
        
        # Try multiple encodings if needed
        errors_to_try = ['utf-8', 'latin-1', 'iso-8859-1']
        for encoding in errors_to_try:
            try:
                if isinstance(content, bytes):
                    content = content.decode(encoding)
                break
            except UnicodeDecodeError:
                continue
        
        # Try to use BioPython's parser first
        with io.StringIO(content) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                return str(record.seq)
        
        # If no FASTA records found, try manual parsing
        lines = content.strip().split('\n')
        sequence = ''
        if lines[0].startswith('>'):
            for line in lines[1:]:
                if not line.startswith('>'):  # Skip header lines
                    sequence += line.strip()
        else:
            # Assume raw sequence
            for line in lines:
                # Remove whitespace and common non-sequence characters
                clean_line = re.sub(r'[\s0-9,\.;:\-_]', '', line)
                sequence += clean_line
                
        return sequence.upper()  # Return standardized uppercase sequence

    except Exception as e:
        logger.error(f"Error parsing FASTA content: {str(e)}")
        return ""

def read_file_content(file_path):
    """Read file content with encoding error handling"""
    for encoding in ['utf-8', 'latin-1', 'iso-8859-1']:
        try:
            with open(file_path, 'r', encoding=encoding) as f:
                return f.read()
        except UnicodeDecodeError:
            continue
    
    # If all encodings fail, try binary mode
    try:
        with open(file_path, 'rb') as f:
            content = f.read()
            # Try to decode with replacement character
            return content.decode('utf-8', errors='replace')
    except Exception as e:
        logger.error(f"Failed to read file with any encoding: {str(e)}")
        raise

def sequence_alignment(seq1, seq2):
    """Perform sequence alignment with basic output."""
    try:
        # Input validation
        if not seq1 or not seq2 or len(seq1) < 3 or len(seq2) < 3:
            return "Sequences too short or invalid for alignment"
            
        aligner = Align.PairwiseAligner()
        aligner.mode = 'global'
        aligner.match_score = 1.0
        aligner.mismatch_score = -1.0
        aligner.gap_score = -2.0
        
        # Get the alignments without using max_alignments
        alignments = list(aligner.align(seq1, seq2))
        
        if not alignments:
            return "No alignment could be generated"
            
        alignment = alignments[0]
        
        # Format aligned sequences
        aligned_str = str(alignment)
        aligned_seqs = aligned_str.split('\n')
        
        if len(aligned_seqs) < 3:
            return "Invalid alignment format returned"
            
        # Format the output exactly as requested
        output = []
        output.append("\nAlignment Score: {:.2f}".format(alignment.score))
        output.append("\nAligned Sequences:")
        output.append(aligned_seqs[0])  # First sequence
        output.append(aligned_seqs[2])  # Second sequence (skip the match line)
        
        return '\n'.join(output)
    except Exception as e:
        logger.error(f"Alignment error: {str(e)}")
        logger.error(traceback.format_exc())
        return f"Alignment error: {str(e)}"

def validate_sequence(sequence):
    """Validate if string contains a proper DNA/RNA/protein sequence."""
    if not sequence:
        return False
    # Expanded validation for more sequence types
    valid_chars = set('ACGTURYKMSWBDHVNacgturykmswbdhvn-*')
    sequence = sequence.strip()
    return bool(sequence) and all(c in valid_chars for c in sequence)

def kmer_frequency(sequence, k=3):
    """Calculate k-mer frequency in a sequence"""
    if not sequence or len(sequence) < k:
        return {}
    kmers = [sequence[i:i+k] for i in range(len(sequence)-k+1)]
    freq = Counter(kmers)
    return dict(freq)

def motif_detection(sequence, motif_length=4):
    """Detect motifs in a sequence"""
    if not sequence or len(sequence) < motif_length:
        return []
    motifs = []
    for i in range(len(sequence)-motif_length+1):
        motif = sequence[i:i+motif_length]
        if sequence.count(motif) > 1:
            motifs.append((motif, sequence.count(motif)))
    return list(set(motifs))

def get_file_type(filename):
    """Determine the type of file based on extension."""
    if not filename:
        return "unknown"
    
    ext = filename.rsplit('.', 1)[-1].lower() if '.' in filename else ""
    
    if ext in ['fasta', 'fa', 'fna', 'ffn', 'faa', 'frn']:
        return "fasta"
    elif ext in ['sif', 'csv', 'tsv', 'txt']:
        if ext == 'sif':
            return "network"
        elif ext in ['csv', 'tsv']:
            return "potential_network"
        else:
            # For txt files, try to determine content type
            return "text"
    
    return "unknown"

def is_network_data(content):
    """Try to determine if content represents network data."""
    if not content:
        return False
        
    lines = content.strip().split('\n')
    if not lines:
        return False

    # Check if lines follow network format (at least two columns)
    valid_lines = 0
    for line in lines[:10]:  # Check first 10 lines
        parts = line.strip().split()
        if len(parts) >= 2:
            try:
                if len(parts) >= 3:
                    float(parts[2])  # Check if third column is numeric
                valid_lines += 1
            except ValueError:
                pass
    return valid_lines > 0

def validate_file_content(content, expected_type="any"):
    """Validate file content matches expected type."""
    if not content or not content.strip():
        return False, "Empty file content"
    
    content = content.strip()
    
    # Check for FASTA format
    if content.startswith('>'):
        sequence = ''.join(line.strip() for line in content.split('\n') 
                         if not line.startswith('>'))
        if validate_sequence(sequence):
            return True, "fasta"
    
    # Check for network data
    if expected_type in ["any", "network"]:
        if is_network_data(content):
            return True, "network"
    
    # If we get here and expected_type is "any", assume it's raw sequence data
    if expected_type == "any":
        cleaned = ''.join(c for c in content if c.isalpha())
        if validate_sequence(cleaned):
            return True, "sequence"
    
    return False, "unknown"

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/analyze', methods=['POST'])
def analyze():
    """Main analysis endpoint for gene sequence processing"""
    file1_path = None
    file2_path = None
    
    try:
        # Input sanitization
        seq1 = request.form.get('sequence1', '').strip()
        seq2 = request.form.get('sequence2', '').strip()
        
        # Handle file uploads with size validation
        if 'file1' in request.files:
            file1 = request.files['file1']
            if file1 and file1.filename and allowed_file(file1.filename):
                try:
                    filename = secure_filename(file1.filename)
                    file1_path = os.path.join(app.config['UPLOAD_FOLDER'], filename)
                    file1.save(file1_path)
                    
                    # Use our improved file reading function
                    file1_content = read_file_content(file1_path)
                    
                    is_valid, content_type = validate_file_content(file1_content)
                    if not is_valid:
                        return render_template('index.html', 
                            error=f"Invalid content in file: {file1.filename}")
                            
                    if content_type in ["fasta", "sequence"]:
                        parsed_seq = parse_fasta_content(file1_content)
                        if parsed_seq:
                            seq1 = parsed_seq
                        else:
                            return render_template('index.html', 
                                error=f"Could not parse sequence from file: {file1.filename}")
                except Exception as e:
                    logger.error(f"Error processing file1: {file1.filename}")
                    logger.error(traceback.format_exc())
                    return render_template('index.html', 
                        error=f"Error reading file {file1.filename}: {str(e)}")

        # Similar handling for file2
        if 'file2' in request.files:
            file2 = request.files['file2']
            if file2 and file2.filename and allowed_file(file2.filename):
                try:
                    filename = secure_filename(file2.filename)
                    file2_path = os.path.join(app.config['UPLOAD_FOLDER'], filename)
                    file2.save(file2_path)
                    
                    # Use our improved file reading function
                    file2_content = read_file_content(file2_path)
                        
                    file2_type = get_file_type(file2.filename)
                    if file2_type in ["network", "potential_network"] or (file2_type == "text" and is_network_data(file2_content)):
                        analysis_method = request.form.get('analysis_method')
                        if analysis_method != 'gene_network':
                            return render_template('index.html', 
                                error=f"The uploaded file '{file2.filename}' appears to be a network file. Please upload it as file 1 and select Gene Network Analysis.")
                    else:
                        parsed_seq = parse_fasta_content(file2_content)
                        if parsed_seq:
                            seq2 = parsed_seq
                except Exception as e:
                    logger.error(f"Error processing file2: {file2.filename}")
                    logger.error(traceback.format_exc())
                    return render_template('index.html', 
                        error=f"Error reading file {file2.filename}: {str(e)}")

        # Get and validate analysis method
        analysis_method = sanitize_input(request.form.get('analysis_method', ''))
        if not analysis_method:
            return render_template('index.html', error="Invalid or missing analysis method")

        # Enhanced validation based on method and available data
        if analysis_method == 'sequence_alignment':
            if not seq2:
                return render_template('index.html', error="Sequence alignment requires two sequences.")
        
        elif analysis_method == 'gene_network':
            edge_data = request.form.get('edge_data', '').strip()
            if not edge_data:
                return render_template('index.html', 
                    error="Gene Network Analysis requires network data. Please upload a network file or enter edge data manually.")
        
        elif analysis_method in ['kmer_frequency', 'motif_detection', 'pattern_search', 'codon_optimization']:
            if not seq1:
                return render_template('index.html', 
                    error=f"{analysis_method.replace('_', ' ').title()} requires at least one sequence.")
        
        # Method-specific validation
        if analysis_method != 'gene_network':  # Skip sequence validation for network analysis
            if not seq1:
                return render_template('index.html', error="Please provide at least one sequence.")
                 
            if seq1 and not validate_sequence(seq1):
                return render_template('index.html', error="Sequence 1 contains invalid characters.")
            
            if seq2 and not validate_sequence(seq2):
                return render_template('index.html', error="Sequence 2 contains invalid characters.")
        
        # Continue with the regular analysis
        result = {}
        
        # Process based on the method
        if analysis_method == 'sequence_alignment':
            # Validate sequences have minimum viable length
            if len(seq1) < 3 or len(seq2) < 3:
                return render_template('index.html', error="Sequences must be at least 3 characters long for alignment.")
                
            # Choose alignment algorithm
            alignment_type = request.form.get('alignment_type', 'biopython')
            
            if alignment_type == 'needleman_wunsch':
                try:
                    score, aligned_seq1, aligned_seq2 = needleman_wunsch(seq1, seq2)
                    result['alignment'] = f"\nAlignment Score: {score}\n\nAligned Sequences:\n{aligned_seq1}\n{aligned_seq2}"
                except Exception as e:
                    logger.error(f"Needleman-Wunsch alignment error: {str(e)}")
                    return render_template('index.html', 
                        error=f"Alignment error: {str(e)}")
                        
            elif alignment_type == 'smith_waterman':
                try:
                    score, aligned_seq1, aligned_seq2 = smith_waterman(seq1, seq2)
                    result['alignment'] = f"\nLocal Alignment Score: {score}\n\nAligned Sequences:\n{aligned_seq1}\n{aligned_seq2}"
                except Exception as e:
                    logger.error(f"Smith-Waterman alignment error: {str(e)}")
                    return render_template('index.html', 
                        error=f"Alignment error: {str(e)}")
            else:
                # Use Bio.Align (default)
                alignment_result = sequence_alignment(seq1, seq2)
                if alignment_result.startswith("Alignment error"):
                    return render_template('index.html', error=alignment_result)
                result['alignment'] = alignment_result
            
        elif analysis_method == 'kmer_frequency':
            result['kmer_freq_seq1'] = kmer_frequency(seq1)
            if seq2:
                result['kmer_freq_seq2'] = kmer_frequency(seq2)
        
        elif analysis_method == 'motif_detection':
            result['motifs_seq1'] = motif_detection(seq1)
            if seq2:
                result['motifs_seq2'] = motif_detection(seq2)
                
        elif analysis_method == 'pattern_search':
            pattern = sanitize_input(request.form.get('pattern', '').strip())
            if not pattern:
                return render_template('index.html', error="Pattern search requires a pattern.")
            
            try:
                matches, frequency = kmp_pattern_search(seq1, pattern)
                result['pattern_search'] = {   
                    'pattern': pattern,
                    'matches': matches,
                    'frequency': frequency
                }
            except Exception as e:
                logger.error(f"Pattern search error: {str(e)}")
                return render_template('index.html', error=f"Pattern search error: {str(e)}")
        
        elif analysis_method == 'codon_optimization':
            if len(seq1) % 3 != 0:
                return render_template('index.html', error="Sequence length must be a multiple of 3 for codon optimization.")
                
            try:
                codon_counts, cai_score, optimized_seq = calculate_cai(seq1)
                result['codon_optimization'] = {
                    'codon_counts': dict(codon_counts),
                    'cai_score': cai_score,
                    'optimized_sequence': optimized_seq
                }
            except Exception as e:
                logger.error(f"Codon optimization error: {str(e)}")
                return render_template('index.html', error=f"Codon optimization error: {str(e)}")
            
        elif analysis_method == 'gene_network':
            # Parse edge data from form
            edge_data = request.form.get('edge_data', '').strip()
            disease_genes = sanitize_input(request.form.get('disease_genes', '').strip())
            
            edges = []
            if edge_data:
                for line in edge_data.split('\n'):
                    if line.strip():
                        parts = line.strip().split()
                        if len(parts) >= 3:
                            try:
                                edges.append((parts[0], parts[1], float(parts[2])))
                            except ValueError:
                                continue
            
            if not edges:
                return render_template('index.html', error="Please provide valid gene interaction data.")
                
            disease_gene_list = [g.strip() for g in disease_genes.split(',') if g.strip()]
            
            try:
                G, network_info = create_gene_network(edges, disease_gene_list)
                result['gene_network'] = network_info
                
                # Check for path finding
                start_gene = sanitize_input(request.form.get('start_gene', '').strip())
                end_gene = sanitize_input(request.form.get('end_gene', '').strip())
                
                if start_gene and end_gene:
                    path_info = find_shortest_path(G, start_gene, end_gene)
                    result['shortest_path'] = path_info
            except Exception as e:
                logger.error(f"Gene network analysis error: {str(e)}")
                return render_template('index.html', error=f"Gene network analysis error: {str(e)}")
        else:
            return render_template('index.html', error=f"Unknown analysis method: {analysis_method}")
        
        # Check if result.html exists, provide fallback rendering if not
        try:
            # Create response with appropriate headers
            response = make_response(render_template('result.html', result=result, method=analysis_method))
            response.headers['Cache-Control'] = 'no-cache, no-store, must-revalidate'
            return response
        except Exception as template_error:
            logger.error(f"Template error: {str(template_error)}")
            # Fallback to simple JSON response if template is missing
            return jsonify({
                'success': True, 
                'result': result, 
                'method': analysis_method
            })

    except Exception as e:
        logger.error(f"Error during analysis: {str(e)}")
        logger.error(traceback.format_exc())
        return render_template('index.html', error=f"Analysis error: {str(e)}")
    finally:
        # Final cleanup of any temporary files
        for path in [file1_path, file2_path]:
            if path and os.path.exists(path):
                try:
                    os.unlink(path)
                except Exception as e:
                    logger.error(f"Error cleaning up file {path}: {str(e)}")

def cleanup_resources():
    """Cleanup any remaining resources"""
    try:
        # Only cleanup if we're the main process
        if multiprocessing.current_process().name == 'MainProcess':
            # Cleanup temp directory
            if os.path.exists(app.config['UPLOAD_FOLDER']):
                for filename in os.listdir(app.config['UPLOAD_FOLDER']):
                    try:
                        filepath = os.path.join(app.config['UPLOAD_FOLDER'], filename)
                        if os.path.isfile(filepath):
                            os.unlink(filepath)
                    except Exception as e:
                        logger.error(f"Error cleaning up file {filepath}: {str(e)}")
    except Exception as e:
        logger.error(f"Error in cleanup: {str(e)}")

def signal_handler(signum, frame):
    """Handle shutdown signals"""
    logger.info("Received shutdown signal, cleaning up...")
    cleanup_resources()
    sys.exit(0)

# Register cleanup handlers
atexit.register(cleanup_resources)
signal.signal(signal.SIGINT, signal_handler)
signal.signal(signal.SIGTERM, signal_handler)

if __name__ == '__main__':
    try:
        # Create upload folder if it doesn't exist
        os.makedirs(app.config['UPLOAD_FOLDER'], exist_ok=True)
        
        # Start the Flask app
        app.run(
            debug=True, 
            host='0.0.0.0', 
            port=5000,
            use_reloader=True
        )
    except Exception as e:
        logger.error(f"Error starting application: {str(e)}")
        sys.exit(1)
    finally:
        cleanup_resources()