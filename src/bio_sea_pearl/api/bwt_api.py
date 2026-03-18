"""API for the Burrows–Wheeler Transform and FM‑index.

This module provides convenience functions to build an FM‑index from a
sequence and to search for patterns.  It reuses the classes defined in
`src/bio_sea_pearl/bwt/transform.py`.
"""

from __future__ import annotations

from typing import List, Tuple

from ..bwt.transform import FMIndex


def build_fm_index(sequence: str, sentinel: str = "$") -> FMIndex:
    """Construct an FM‑index for a given sequence.

    Args:
        sequence: The input string to index.
        sentinel: A sentinel character not present in the sequence.

    Returns:
        An `FMIndex` instance.
    """
    return FMIndex(sequence, sentinel=sentinel)


def search_fm_index(index: FMIndex, pattern: str) -> List[int]:
    """Search for a pattern using a prebuilt FM‑index.

    Args:
        index: An `FMIndex` object.
        pattern: The substring to search for.

    Returns:
        A list of starting positions where the pattern occurs in the indexed sequence.
    """
    l, r = index.backward_search(pattern)
    return list(index.sa[l:r])
