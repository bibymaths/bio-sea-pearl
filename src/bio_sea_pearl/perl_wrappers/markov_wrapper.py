"""Wrapper for Markov chain simulators.

This module exposes a Python function to generate random walks based on the
existing Perl Markov chain implementation.  It invokes the `randomwalk.pl`
script via `subprocess` and returns the generated sequence.
"""

from __future__ import annotations

import subprocess
from pathlib import Path
from typing import Optional


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[3]


def run_markov_walk(
    fasta: str,
    length: int,
    start: str = "A",
    order: int = 1,
    method: str = "alias",
    pseudocount: int = 0,
) -> str:
    """Generate a random walk from a Markov model built from a FASTA file.

    Args:
        fasta: Path to the input FASTA file containing training sequences.
        length: Number of characters in the generated walk.
        start: Starting state (single character).  Defaults to 'A'.
        order: Order of the Markov chain (1 for first‑order, 2 for second‑order, etc.).
        method: Sampling method ('alias' or 'binsrch').  Defaults to 'alias'.
        pseudocount: Pseudocount added to all transitions to avoid zeros.

    Returns:
        A string representing the random walk.
    """
    root = _repo_root()
    script = root / "markov" / "bin" / "randomwalk.pl"
    cmd = [
        "perl",
        str(script),
        "--fasta",
        fasta,
        "--length",
        str(length),
        "--start",
        start,
        "--order",
        str(order),
        "--method",
        method,
    ]
    if pseudocount:
        cmd += ["--pseudocount", str(pseudocount)]
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=root, check=True)
    return result.stdout