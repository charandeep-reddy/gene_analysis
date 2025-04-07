import numpy as np

def needleman_wunsch(seq1, seq2, match=1, mismatch=-1, gap=-2):
    """
    Performs global sequence alignment using the Needleman-Wunsch algorithm.
    
    Args:
        seq1: First sequence
        seq2: Second sequence
        match: Score for matching characters (default: 1)
        mismatch: Score for mismatched characters (default: -1)
        gap: Gap penalty (default: -2)
    
    Returns:
        Tuple of (alignment score, aligned sequence 1, aligned sequence 2)
    """
    m, n = len(seq1), len(seq2)
    dp = np.zeros((m+1, n+1))

    # Initialize first row and column
    for i in range(m+1):
        dp[i][0] = i * gap
    for j in range(n+1):
        dp[0][j] = j * gap

    # Fill the matrix
    for i in range(1, m+1):
        for j in range(1, n+1):
            match_score = match if seq1[i-1] == seq2[j-1] else mismatch
            dp[i][j] = max(
                dp[i-1][j-1] + match_score,
                dp[i-1][j] + gap,
                dp[i][j-1] + gap
            )
    
    # Traceback to find the alignment
    aligned_seq1, aligned_seq2 = "", ""
    i, j = m, n
    while i > 0 or j > 0:
        if i > 0 and j > 0 and (dp[i][j] == dp[i-1][j-1] + (match if seq1[i-1] == seq2[j-1] else mismatch)):
            aligned_seq1 = seq1[i-1] + aligned_seq1
            aligned_seq2 = seq2[j-1] + aligned_seq2
            i -= 1
            j -= 1
        elif i > 0 and (dp[i][j] == dp[i-1][j] + gap):
            aligned_seq1 = seq1[i-1] + aligned_seq1
            aligned_seq2 = "-" + aligned_seq2
            i -= 1
        else:
            aligned_seq1 = "-" + aligned_seq1
            aligned_seq2 = seq2[j-1] + aligned_seq2
            j -= 1

    return dp[m][n], aligned_seq1, aligned_seq2

def smith_waterman(seq1, seq2, match=1, mismatch=-1, gap=-2):
    """
    Performs local sequence alignment using the Smith-Waterman algorithm.
    
    Args:
        seq1: First sequence
        seq2: Second sequence
        match: Score for matching characters (default: 1)
        mismatch: Score for mismatched characters (default: -1)
        gap: Gap penalty (default: -2)
    
    Returns:
        Tuple of (alignment score, aligned sequence 1, aligned sequence 2)
    """
    m, n = len(seq1), len(seq2)
    dp = np.zeros((m+1, n+1))

    max_score = 0
    max_i, max_j = 0, 0

    # Fill the matrix
    for i in range(1, m+1):
        for j in range(1, n+1):
            match_score = match if seq1[i-1] == seq2[j-1] else mismatch
            dp[i][j] = max(
                0,
                dp[i-1][j-1] + match_score,
                dp[i-1][j] + gap,
                dp[i][j-1] + gap
            )
            if dp[i][j] > max_score:
                max_score = dp[i][j]
                max_i, max_j = i, j
    
    # Traceback to find the alignment
    aligned_seq1 = ""
    aligned_seq2 = ""
    i, j = max_i, max_j

    while i > 0 and j > 0 and dp[i][j] > 0:
        if dp[i][j] == dp[i-1][j-1] + (match if seq1[i-1] == seq2[j-1] else mismatch):
            aligned_seq1 = seq1[i-1] + aligned_seq1
            aligned_seq2 = seq2[j-1] + aligned_seq2
            i -= 1
            j -= 1
        elif dp[i][j] == dp[i-1][j] + gap:
            aligned_seq1 = seq1[i-1] + aligned_seq1
            aligned_seq2 = "-" + aligned_seq2
            i -= 1
        else:
            aligned_seq1 = "-" + aligned_seq1
            aligned_seq2 = seq2[j-1] + aligned_seq2
            j -= 1

    return max_score, aligned_seq1, aligned_seq2