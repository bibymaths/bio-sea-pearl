"""Markov API for the Bio Sea Pearl toolkit.

This module exposes a simple function for generating random walks from a
first- or higher-order Markov chain based on training sequences in a FASTA
file.  It calls the underlying Perl implementation through the wrapper.
"""

from __future__ import annotations

from typing import Optional

from ..perl_wrappers import run_markov_walk


def generate_walk(
    fasta: str,
    length: int,
    start: str = "A",
    order: int = 1,
    method: str = "alias",
    pseudocount: int = 0,
) -> str:
    """Generate a random walk from a Markov model learned from a FASTA file.

    Args:
        fasta: Path to the FASTA file with training sequences.
        length: Length of the random walk to generate.
        start: Starting state (usually a single character).
        order: Order of the Markov chain (1 for first‑order, 2 for second‑order, etc.).
        method: Sampling method ("alias" or "binsrch").
        pseudocount: Non‑negative integer pseudocount to add to transitions.

    Returns:
        A string representing the generated random sequence.
    """
    return run_markov_walk(fasta, length, start=start, order=order, method=method, pseudocount=pseudocount)