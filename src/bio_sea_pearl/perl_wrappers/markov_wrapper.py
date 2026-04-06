"""Wrapper for Markov chain simulators.

This module exposes a Python function to generate random walks based on the
existing Perl Markov chain implementation.  It invokes the ``randomwalk.pl``
script via ``subprocess`` and returns the generated sequence.
"""

from __future__ import annotations

import os
import subprocess
from pathlib import Path


def _repo_root() -> Path:
    """Return the repository root directory.

    Checks the ``BIOSEA_REPO_ROOT`` environment variable first so that
    containerised or installed deployments can override the default
    heuristic (three levels above this file in the source tree).
    """
    env = os.environ.get("BIOSEA_REPO_ROOT")
    if env:
        return Path(env)
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
        start: Starting state (single character).  Defaults to ``'A'``.
        order: Order of the Markov chain (1 for first-order, 2+ for higher).
        method: Sampling method (``'alias'`` or ``'binsrch'``).
        pseudocount: Pseudocount added to all transitions to avoid zeros.

    Returns:
        A string representing the random walk.

    Raises:
        subprocess.CalledProcessError: If the Perl command fails.
        FileNotFoundError: If the randomwalk.pl script is not found.
    """
    root = _repo_root()
    script = root / "markov" / "bin" / "randomwalk.pl"
    if not script.exists():
        raise FileNotFoundError(f"randomwalk.pl not found at {script}")
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
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=str(root), check=True)
    return result.stdout
