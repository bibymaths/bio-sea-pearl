"""K‑mer counting in pure Python.

Provides a function `counts` that computes the frequency of all k‑mers in
a given sequence.  The input sequence can be any iterable supporting
slicing; however, typical usage is with strings.
"""

from __future__ import annotations

from typing import Dict


def counts(seq: str, k: int) -> Dict[str, int]:
    """Return a dictionary mapping each k‑mer to its frequency.

    Args:
        seq: The input sequence.
        k: The length of k‑mers to count.  Must be ≤ len(seq).

    Raises:
        ValueError: If k is non‑positive or greater than the length of seq.

    Returns:
        A dictionary with k‑mers as keys and counts as values.
    """
    n = len(seq)
    if k <= 0:
        raise ValueError("k must be positive")
    if k > n:
        raise ValueError(f"k ({k}) must be ≤ sequence length ({n})")
    result: Dict[str, int] = {}
    end = n - k
    for i in range(end + 1):
        kmer = seq[i: i + k]
        result[kmer] = result.get(kmer, 0) + 1
    return result
