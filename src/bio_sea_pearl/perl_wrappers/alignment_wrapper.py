"""Wrapper for sequence alignment using existing scripts.

This module provides a Python function that invokes the existing alignment
implementations.  It prefers the Python implementation (`align.py`) because
it integrates cleanly with the Python runtime, but can fall back to the Perl
implementation if necessary.  The wrapper constructs an external command
using `subprocess.run` and returns the standard output.

Example:

    from bio_sea_pearl.perl_wrappers import run_alignment
    result = run_alignment("seq1.fasta", "seq2.fasta",
                          matrix="alignment/scoring/BLOSUM62.mat",
                          mode="global")
    print(result)
"""

from __future__ import annotations

import subprocess
from pathlib import Path
from typing import Optional


def _repo_root() -> Path:
    """Return the root directory of the repository (two levels up)."""
    return Path(__file__).resolve().parents[3]


def run_alignment(
    fasta1: str,
    fasta2: str,
    matrix: Optional[str] = None,
    mode: str = "global",
) -> str:
    """Run the sequence alignment on two FASTA files and return the result.

    Args:
        fasta1: Path to the first FASTA file.
        fasta2: Path to the second FASTA file.
        matrix: Optional path to a scoring matrix file.  If not provided,
            a default matrix will be used by the underlying script.
        mode: Alignment mode (`global`, `local` or `lcs`).

    Returns:
        The standard output of the alignment command.

    Raises:
        subprocess.CalledProcessError: If the underlying command exits with
            a non-zero status.
    """
    root = _repo_root()
    # Prefer the Python implementation if present
    align_py = root / "alignment" / "bin" / "align.py"
    if align_py.exists():
        cmd = [
            "python",
            str(align_py),
            "--mode",
            mode,
        ]
        if matrix:
            cmd += ["--matrix", matrix]
        cmd += [fasta1, fasta2]
    else:
        # Fallback to Perl script
        align_pl = root / "alignment" / "bin" / "align.pl"
        cmd = [
            "perl",
            str(align_pl),
            "--mode",
            mode,
        ]
        if matrix:
            cmd += ["--matrix", matrix]
        cmd += [fasta1, fasta2]
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=root, check=True)
    return result.stdout