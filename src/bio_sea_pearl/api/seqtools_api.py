"""SeqTools API for the Bio Sea Pearl toolkit.

This module defines user-friendly functions for computing common sequence
metrics and k‑mer counts.  Whenever possible, it uses the native Python
implementations from the `seqtools_py` package.  When those are not
available or if an error occurs, it falls back to the Perl wrappers.
"""

from __future__ import annotations

from typing import Dict

try:
    # Try to import the Python implementations
    from ..seqtools_py.hamming import distance as _hamming_py
    from ..seqtools_py.levenshtein import distance as _levenshtein_py
    from ..seqtools_py.kmer import counts as _kmer_py
    _PY_PORT_AVAILABLE = True
except Exception:
    _PY_PORT_AVAILABLE = False

from ..perl_wrappers import (
    distance_hamming_perl,
    distance_levenshtein_perl,
    kmer_counts_perl,
)


def hamming_distance(a: str, b: str) -> int:
    """Compute the Hamming distance between two sequences.

    Falls back to the Perl implementation if the Python port is unavailable.
    """
    if _PY_PORT_AVAILABLE:
        try:
            return _hamming_py(a, b)
        except Exception:
            pass
    return distance_hamming_perl(a, b)


def levenshtein_distance(a: str, b: str) -> int:
    """Compute the Levenshtein distance between two sequences.

    Falls back to the Perl implementation if the Python port is unavailable.
    """
    if _PY_PORT_AVAILABLE:
        try:
            return _levenshtein_py(a, b)
        except Exception:
            pass
    return distance_levenshtein_perl(a, b)


def kmer_counts(seq: str, k: int) -> Dict[str, int]:
    """Compute k‑mer counts for a sequence.

    Falls back to the Perl implementation if the Python port is unavailable.
    """
    if _PY_PORT_AVAILABLE:
        try:
            return _kmer_py(seq, k)
        except Exception:
            pass
    return kmer_counts_perl(seq, k)