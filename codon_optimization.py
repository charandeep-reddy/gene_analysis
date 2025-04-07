import numpy as np
from collections import Counter

# Default reference codon usage frequencies for common organisms (human)
DEFAULT_CODON_USAGE = {
    'ATC': 0.48, 'ATA': 0.18, 'ATT': 0.34,  # Isoleucine
    'GAA': 0.68, 'GAG': 0.32,  # Glutamic Acid
    'TTC': 0.55, 'TTT': 0.45,  # Phenylalanine
    'GGT': 0.36, 'GGC': 0.37, 'GGA': 0.15, 'GGG': 0.12,  # Glycine
    'CTT': 0.13, 'CTC': 0.20, 'CTA': 0.07, 'CTG': 0.60,  # Leucine
    'TAT': 0.57, 'TAC': 0.43,  # Tyrosine
    # Add more codons as needed
}

def calculate_cai(dna_seq, reference_codon_usage=DEFAULT_CODON_USAGE):
    """
    Calculate the Codon Adaptation Index (CAI) for a DNA sequence.
    
    Args:
        dna_seq: DNA sequence as a string
        reference_codon_usage: Dictionary mapping codons to their usage frequency
        
    Returns:
        Tuple of (codon_counts, CAI score, optimized DNA sequence)
    """
    # Ensure sequence length is a multiple of 3
    if len(dna_seq) % 3 != 0:
        dna_seq = dna_seq[:-(len(dna_seq) % 3)]
    
    # Split into codons
    codons = [dna_seq[i:i+3] for i in range(0, len(dna_seq), 3)]
    codon_counts = Counter(codons)
    
    # Calculate weights and create optimized sequence
    codon_weights = {}
    optimized_sequence = []
    
    for codon in codons:
        # Find synonymous codons (codons with same first two letters)
        synonymous_codons = {c: reference_codon_usage.get(c, 0) for c in reference_codon_usage 
                            if c[:2] == codon[:2]}
        
        if synonymous_codons:
            max_freq_codon = max(synonymous_codons, key=synonymous_codons.get)
            optimized_sequence.append(max_freq_codon)
            
            if codon in reference_codon_usage:
                max_freq = max(synonymous_codons.values(), default=1)
                if max_freq > 0:
                    codon_weights[codon] = reference_codon_usage[codon] / max_freq
                else:
                    codon_weights[codon] = 0
            else:
                codon_weights[codon] = 0
        else:
            optimized_sequence.append(codon)
            codon_weights[codon] = 0
    
    # Calculate CAI score
    non_zero_weights = [w for w in codon_weights.values() if w > 0]
    if not non_zero_weights:
        cai_score = 0.0
    else:
        cai_score = np.exp(np.mean(np.log(non_zero_weights)))
    
    optimized_dna_seq = ''.join(optimized_sequence)
    
    return codon_counts, cai_score, optimized_dna_seq