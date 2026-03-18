"""Alignment API for the Bio Sea Pearl toolkit.

Provides a function `align_sequences` that computes an alignment between two
FASTA files.  Under the hood it calls the alignment wrapper, which
delegates to either the Python or Perl implementation.
"""

from __future__ import annotations

from typing import Optional

from ..perl_wrappers import run_alignment


def align_sequences(
    fasta1: str,
    fasta2: str,
    matrix: Optional[str] = None,
    mode: str = "global",
) -> str:
    """Align two sequences stored in FASTA files.

    Args:
        fasta1: Path to the first FASTA file.
        fasta2: Path to the second FASTA file.
        matrix: Optional path to a substitution matrix.
        mode: Alignment mode (`global`, `local`, or `lcs`).

    Returns:
        Alignment result as a string.
    """
    return run_alignment(fasta1, fasta2, matrix=matrix, mode=mode)