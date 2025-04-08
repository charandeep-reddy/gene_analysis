from flask import Flask, render_template, request, jsonify, flash
from Bio import SeqIO, Align
from Bio.Seq import Seq
import numpy as np
from collections import Counter
import io
import re
import json
import logging

# Setup basic logging
logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)

# Import our custom modules
try:
    from alignment import needleman_wunsch, smith_waterman
    from codon_optimization import calculate_cai, DEFAULT_CODON_USAGE
    from gene_network import create_gene_network, find_shortest_path
    from pattern_search import kmp_pattern_search
except ImportError as e:
    logger.error(f"Failed to import required modules: {str(e)}")
    raise

app = Flask(__name__)
app.secret_key = 'your-secret-key-here'  # Required for flash messages

def parse_fasta_content(content):
    """Parse FASTA content and return the sequence."""
    try:
        # Try to use BioPython's parser first
        with io.StringIO(content) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                return str(record.seq)
        
        # If no FASTA records found, try manual parsing
        lines = content.split('\n')
        sequence = ''
        for line in lines:
            if not line.startswith('>'):  # Skip header lines
                sequence += line.strip()
        return sequence
    except Exception as e:
        logger.error(f"Error parsing FASTA content: {str(e)}")
        return ""

def sequence_alignment(seq1, seq2):
    """Perform sequence alignment with basic output."""
    try:
        aligner = Align.PairwiseAligner()
        aligner.mode = 'global'
        aligner.match_score = 1.0
        aligner.mismatch_score = -1.0
        aligner.gap_score = -2.0
        
        # Get the alignment object
        alignments = aligner.align(seq1, seq2)
        alignment = alignments[0]
        
        # Format aligned sequences
        aligned_seqs = str(alignment).split('\n')
        
        # Format the output exactly as requested
        output = []
        output.append("\nAlignment Score: {}".format(alignment.score))
        output.append("\nAligned Sequences:")
        output.append(aligned_seqs[0])  # First sequence
        output.append(aligned_seqs[2])  # Second sequence (skip the match line)
        
        return '\n'.join(output)
    except Exception as e:
        logger.error(f"Alignment error: {str(e)}")
        return f"Alignment error: {str(e)}"

def validate_sequence(sequence):
    """Validate if string contains a proper DNA/RNA/protein sequence."""
    if not sequence:
        return False
    # Simple validation - check if sequence contains only valid characters
    valid_chars = set('ACGTURYKMSWBDHVN-acgturykmswbdhvn')
    return all(c in valid_chars for c in sequence)

def kmer_frequency(sequence, k=3):
    if len(sequence) < k:
        return {}
    kmers = [sequence[i:i+k] for i in range(len(sequence)-k+1)]
    freq = Counter(kmers)
    return dict(freq)

def motif_detection(sequence, motif_length=4):
    if len(sequence) < motif_length:
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
    
    ext = filename.rsplit('.', 1)[-1].lower()
    
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

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/analyze', methods=['POST'])
def analyze():
    try:
        # Get input sequences
        seq1 = request.form.get('sequence1', '').strip()
        seq2 = request.form.get('sequence2', '').strip()
        
        # Initialize file types
        file1_type = "none"
        file2_type = "none"
        file1_content = None
        file2_content = None
        
        # Handle file uploads for file1
        if 'file1' in request.files and request.files['file1'].filename:
            file1 = request.files['file1']
            try:
                file1_type = get_file_type(file1.filename)
                file1_content = file1.read().decode('utf-8')
                
                # If it's a network file, process accordingly
                if file1_type in ["network", "potential_network"] or (file1_type == "text" and is_network_data(file1_content)):
                    # Check if we're using the correct analysis method
                    analysis_method = request.form.get('analysis_method')
                    if analysis_method != 'gene_network':
                        return render_template('index.html', 
                                              error=f"The uploaded file '{file1.filename}' appears to be a network file. Please select Gene Network Analysis method.")
                    
                    # Store the content in the form for network processing
                    edge_data = []
                    for line in file1_content.strip().split('\n'):
                        parts = line.strip().split()
                        if len(parts) >= 2:
                            # Default weight to 1.0 if not provided
                            weight = 1.0
                            if len(parts) >= 3:
                                try:
                                    weight = float(parts[2])
                                except ValueError:
                                    pass
                            edge_data.append(f"{parts[0]} {parts[1]} {weight}")
                    
                    # Use werkzeug's MultiDict to update the form data
                    from werkzeug.datastructures import MultiDict
                    form_data = MultiDict(request.form)
                    form_data['edge_data'] = '\n'.join(edge_data)
                    request.form = form_data
                else:
                    # Regular sequence processing
                    parsed_seq = parse_fasta_content(file1_content)
                    if parsed_seq:
                        seq1 = parsed_seq
            except Exception as e:
                logger.exception(f"Error processing file1: {file1.filename}")
                return render_template('index.html', error=f"Error reading file 1: {str(e)}")
        
        # Handle file uploads for file2
        if 'file2' in request.files and request.files['file2'].filename:
            file2 = request.files['file2']
            try:
                file2_type = get_file_type(file2.filename)
                file2_content = file2.read().decode('utf-8')
                
                if file2_type in ["network", "potential_network"] or (file2_type == "text" and is_network_data(file2_content)):
                    # Network files in file2 slot are not commonly used - show warning
                    analysis_method = request.form.get('analysis_method')
                    if analysis_method != 'gene_network':
                        return render_template('index.html', 
                                              error=f"The uploaded file '{file2.filename}' appears to be a network file. Please upload it as file 1 and select Gene Network Analysis.")
                else:
                    # Regular sequence processing
                    parsed_seq = parse_fasta_content(file2_content)
                    if parsed_seq:
                        seq2 = parsed_seq
            except Exception as e:
                logger.exception(f"Error processing file2: {file2.filename}")
                return render_template('index.html', error=f"Error reading file 2: {str(e)}")

        # Get analysis method
        analysis_method = request.form.get('analysis_method')
        
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
            if not seq1 and analysis_method != 'gene_network':
                return render_template('index.html', error="Please provide at least one sequence.")
                
            if seq1 and not validate_sequence(seq1):
                return render_template('index.html', error="Sequence 1 contains invalid characters.")
                
            if seq2 and not validate_sequence(seq2):
                return render_template('index.html', error="Sequence 2 contains invalid characters.")

        # Continue with the regular analysis
        result = {}
        
        # Process based on the method (existing code)
        if analysis_method == 'sequence_alignment':
            if not seq2:
                return render_template('index.html', error="Sequence alignment requires two sequences.")
            
            # Choose alignment algorithm
            alignment_type = request.form.get('alignment_type', 'biopython')
            
            if alignment_type == 'needleman_wunsch':
                score, aligned_seq1, aligned_seq2 = needleman_wunsch(seq1, seq2)
                result['alignment'] = f"\nAlignment Score: {score}\n\nAligned Sequences:\n{aligned_seq1}\n{aligned_seq2}"
            elif alignment_type == 'smith_waterman':
                score, aligned_seq1, aligned_seq2 = smith_waterman(seq1, seq2)
                result['alignment'] = f"\nLocal Alignment Score: {score}\n\nAligned Sequences:\n{aligned_seq1}\n{aligned_seq2}"
            else:
                # Use Bio.Align (default)
                result['alignment'] = sequence_alignment(seq1, seq2)
            
        elif analysis_method == 'kmer_frequency':
            result['kmer_freq_seq1'] = kmer_frequency(seq1)
            if seq2:
                result['kmer_freq_seq2'] = kmer_frequency(seq2)
                
        elif analysis_method == 'motif_detection':
            result['motifs_seq1'] = motif_detection(seq1)
            if seq2:
                result['motifs_seq2'] = motif_detection(seq2)
                
        elif analysis_method == 'pattern_search':
            pattern = request.form.get('pattern', '').strip()
            if not pattern:
                return render_template('index.html', error="Pattern search requires a pattern.")
            
            matches, frequency = kmp_pattern_search(seq1, pattern)
            result['pattern_search'] = {
                'pattern': pattern,
                'matches': matches,
                'frequency': frequency
            }
        
        elif analysis_method == 'codon_optimization':
            if len(seq1) % 3 != 0:
                return render_template('index.html', error="Sequence length must be a multiple of 3 for codon optimization.")
                
            codon_counts, cai_score, optimized_seq = calculate_cai(seq1)
            result['codon_optimization'] = {
                'codon_counts': dict(codon_counts),
                'cai_score': cai_score,
                'optimized_sequence': optimized_seq
            }
            
        elif analysis_method == 'gene_network':
            # Parse edge data from form
            edge_data = request.form.get('edge_data', '').strip()
            disease_genes = request.form.get('disease_genes', '').strip()
            
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
            
            G, network_info = create_gene_network(edges, disease_gene_list)
            result['gene_network'] = network_info
            
            # Check for path finding
            start_gene = request.form.get('start_gene', '').strip()
            end_gene = request.form.get('end_gene', '').strip()
            
            if start_gene and end_gene:
                path_info = find_shortest_path(G, start_gene, end_gene)
                result['shortest_path'] = path_info
        else:
            return render_template('index.html', error=f"Unknown analysis method: {analysis_method}")

        return render_template('result.html', result=result, method=analysis_method)
    except Exception as e:
        logger.exception("Unexpected error in analyze route")
        return render_template('index.html', error=f"Analysis error: {str(e)}")

if __name__ == '__main__':
    app.run(debug=True)