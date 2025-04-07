def kmp_pattern_search(text, pattern):
    """
    Knuth-Morris-Pratt algorithm for efficient pattern matching in sequences.
    
    Args:
        text: The text to search in
        pattern: The pattern to search for
        
    Returns:
        Tuple of (list of match positions, frequency)
    """
    def compute_lps(pattern):
        """Compute Longest Prefix Suffix array for KMP algorithm"""
        lps = [0] * len(pattern)
        length = 0
        i = 1
        while i < len(pattern):
            if pattern[i] == pattern[length]:
                length += 1
                lps[i] = length
                i += 1
            else:
                if length:
                    length = lps[length - 1]
                else:
                    lps[i] = 0
                    i += 1
        return lps

    lps = compute_lps(pattern)
    i = j = 0
    matches = []
    
    while i < len(text):
        if pattern[j] == text[i]:
            i += 1
            j += 1
        if j == len(pattern):
            matches.append(i - j)
            j = lps[j - 1]
        elif i < len(text) and pattern[j] != text[i]:
            if j != 0:
                j = lps[j - 1]
            else:
                i += 1

    frequency = len(matches)
    return matches, frequency