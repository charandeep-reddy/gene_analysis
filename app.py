from flask import Flask, render_template, request, jsonify, flash
from Bio import SeqIO, Align
from Bio.Seq import Seq
import numpy as np
from collections import Counter
import io
import re
import json

# Import our custom modules
from alignment import needleman_wunsch, smith_waterman
from codon_optimization import calculate_cai, DEFAULT_CODON_USAGE
from gene_network import create_gene_network, find_shortest_path
from pattern_search import kmp_pattern_search

app = Flask(__name__)
app.secret_key = 'your-secret-key-here'  # Required for flash messages

def parse_fasta_content(content):
    """Parse FASTA content and return the sequence."""
    lines = content.split('\n')
    sequence = ''
    for line in lines:
        if not line.startswith('>'):  # Skip header lines
            sequence += line.strip()
    return sequence

def sequence_alignment(seq1, seq2):
    """Perform sequence alignment with basic output."""
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

def kmer_frequency(sequence, k=3):
    kmers = [sequence[i:i+k] for i in range(len(sequence)-k+1)]
    freq = Counter(kmers)
    return dict(freq)

def motif_detection(sequence, motif_length=4):
    motifs = []
    for i in range(len(sequence)-motif_length+1):
        motif = sequence[i:i+motif_length]
        if sequence.count(motif) > 1:
            motifs.append((motif, sequence.count(motif)))
    return list(set(motifs))

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/analyze', methods=['POST'])
def analyze():
    # Get input sequences
    seq1 = request.form.get('sequence1', '').strip()
    seq2 = request.form.get('sequence2', '').strip()
    
    # Handle file uploads if present
    if 'file1' in request.files:
        file1 = request.files['file1']
        if file1.filename:
            content = file1.read().decode('utf-8')
            seq1 = parse_fasta_content(content)
    
    if 'file2' in request.files:
        file2 = request.files['file2']
        if file2.filename:
            content = file2.read().decode('utf-8')
            seq2 = parse_fasta_content(content)

    # Validate input sequences
    if not seq1:
        return render_template('index.html', error="Please provide at least one sequence.")

    analysis_method = request.form.get('analysis_method')
    result = {}

    try:
        if analysis_method == 'sequence_alignment':
            if not seq2:
                return render_template('index.html', error="Sequence alignment requires two sequences.")
            
            # Choose alignment algorithm
            alignment_type = request.form.get('alignment_type', 'global')
            
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
            
            disease_gene_list = [g.strip() for g in disease_genes.split(',') if g.strip()]
            
            G, network_info = create_gene_network(edges, disease_gene_list)
            result['gene_network'] = network_info
            
            # Check for path finding
            start_gene = request.form.get('start_gene', '').strip()
            end_gene = request.form.get('end_gene', '').strip()
            
            if start_gene and end_gene:
                path_info = find_shortest_path(G, start_gene, end_gene)
                result['shortest_path'] = path_info

        return render_template('result.html', result=result, method=analysis_method)
    except Exception as e:
        return render_template('index.html', error=f"Analysis error: {str(e)}")

if __name__ == '__main__':
    app.run(debug=True)