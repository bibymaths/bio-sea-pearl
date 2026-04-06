"""Wrapper for sequence alignment using existing scripts.

This module provides a Python function that invokes the existing alignment
implementations.  It prefers the Python implementation (``align.py``) because
it integrates cleanly with the Python runtime, but can fall back to the Perl
implementation if necessary.  The wrapper constructs an external command
using ``subprocess.run`` and returns the standard output.

Example::

    from bio_sea_pearl.perl_wrappers import run_alignment
    result = run_alignment("seq1.fasta", "seq2.fasta",
                          matrix="alignment/scoring/blosum62.mat",
                          mode="global")
    print(result)
"""

from __future__ import annotations

import os
import subprocess
from pathlib import Path
from typing import Optional


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
        mode: Alignment mode (``global``, ``local`` or ``lcs``).

    Returns:
        The standard output of the alignment command.

    Raises:
        subprocess.CalledProcessError: If the underlying command exits with
            a non-zero status.
        FileNotFoundError: If neither the Python nor Perl align script exists,
            or if the FASTA input files do not exist.
    """
    root = _repo_root()

    # Validate that the input files exist to prevent arbitrary path injection.
    for label, path_str in (("fasta1", fasta1), ("fasta2", fasta2)):
        p = Path(path_str)
        if not p.is_file():
            raise FileNotFoundError(f"{label} file does not exist: {path_str}")

    # Resolve the matrix path; fall back to a bundled default when the
    # caller does not specify one.
    if matrix:
        mat_path = Path(matrix)
        if not mat_path.is_absolute():
            mat_path = root / mat_path
    else:
        mat_path = root / "alignment" / "scoring" / "blosum62.mat"

    align_py = root / "alignment" / "bin" / "align.py"
    if align_py.exists():
        cmd: list[str] = [
            "python",
            str(align_py),
            "--mode",
            mode,
            "--matrix",
            str(mat_path),
            fasta1,
            fasta2,
        ]
    else:
        align_pl = root / "alignment" / "bin" / "align.pl"
        if not align_pl.exists():
            raise FileNotFoundError(
                f"Neither align.py nor align.pl found under {root / 'alignment' / 'bin'}"
            )
        cmd = [
            "perl",
            str(align_pl),
            "--mode",
            mode,
            "--matrix",
            str(mat_path),
            fasta1,
            fasta2,
        ]

    result = subprocess.run(cmd, capture_output=True, text=True, cwd=str(root), check=True)
    return result.stdout
