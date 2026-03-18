"""Hamming distance implementation in pure Python.

This module defines a single function, `distance`, which computes the
Hamming distance between two equal‑length strings.  It raises a
`ValueError` if the sequences have different lengths.
"""

from __future__ import annotations

def distance(a: str, b: str) -> int:
    """Return the number of positions at which the two strings differ.

    Args:
        a: First sequence (must be same length as `b`).
        b: Second sequence.

    Raises:
        ValueError: If the sequences have different lengths.

    Returns:
        The Hamming distance as an integer.
    """
    if len(a) != len(b):
        raise ValueError(f"Sequences must be of equal length: got {len(a)} and {len(b)}")
    return sum(ch1 != ch2 for ch1, ch2 in zip(a, b))