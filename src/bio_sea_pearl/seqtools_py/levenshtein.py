"""Levenshtein distance implementation in pure Python.

This module defines a function `distance` which computes the Levenshtein
(edit) distance between two strings.  It is based on a dynamic
programming algorithm that runs in O(n*m) time and O(min(n,m)) space.
"""

from __future__ import annotations

def distance(s: str, t: str) -> int:
    """Compute the Levenshtein distance between two strings.

    Args:
        s: First sequence.
        t: Second sequence.

    Returns:
        The edit distance between `s` and `t`.
    """
    # Ensure s is the shorter string to reduce memory usage
    if len(s) > len(t):
        s, t = t, s
    n, m = len(s), len(t)
    # Previous and current row of the DP table
    prev = list(range(m + 1))
    curr = [0] * (m + 1)
    for i in range(1, n + 1):
        curr[0] = i
        sc = s[i - 1]
        for j in range(1, m + 1):
            cost = 0 if sc == t[j - 1] else 1
            ins = curr[j - 1] + 1
            delete = prev[j] + 1
            substitute = prev[j - 1] + cost
            # Choose minimum of insertion, deletion and substitution
            curr[j] = ins if ins <= delete and ins <= substitute else (delete if delete <= substitute else substitute)
        prev, curr = curr, prev
    return prev[m]