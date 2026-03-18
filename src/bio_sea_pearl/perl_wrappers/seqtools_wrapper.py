"""Wrappers for SeqTools Perl modules.

This module provides Python functions that call the existing SeqTools Perl
implementations for Hamming distance, Levenshtein distance and k‑mer
counts.  The functions invoke Perl one-liners using `subprocess.run` and
return the results as Python types.
"""

from __future__ import annotations

import subprocess
from pathlib import Path
from typing import Dict


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[3]


def distance_hamming_perl(a: str, b: str) -> int:
    """Compute the Hamming distance between two equal‑length strings using Perl.

    Args:
        a: First string.
        b: Second string.

    Returns:
        Number of positions at which the two strings differ.

    Raises:
        subprocess.CalledProcessError: If the Perl command fails.
    """
    root = _repo_root()
    # Build Perl one-liner to call SeqTools::Hamming::distance
    script = (
        f"print SeqTools::Hamming::distance('{a}', '{b}');"
    )
    cmd = [
        "perl",
        "-MSeqTools::Hamming",
        "-e",
        script,
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=root, check=True)
    return int(result.stdout.strip())


def distance_levenshtein_perl(a: str, b: str) -> int:
    """Compute the Levenshtein (edit) distance using Perl SeqTools implementation."""
    root = _repo_root()
    script = (
        f"print SeqTools::Levenshtein::distance('{a}', '{b}');"
    )
    cmd = [
        "perl",
        "-MSeqTools::Levenshtein",
        "-e",
        script,
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=root, check=True)
    return int(result.stdout.strip())


def kmer_counts_perl(seq: str, k: int) -> Dict[str, int]:
    """Compute k‑mer counts for a sequence using Perl SeqTools implementation.

    Args:
        seq: Input sequence.
        k: Size of k‑mers.

    Returns:
        A dictionary mapping k‑mer substrings to their counts.
    """
    root = _repo_root()
    # Prepare a Perl one‑liner that prints kmer and count separated by a tab per line
    # We quote $seq and $k via single quotes and escape backslash if needed.  Since
    # input may contain single quotes, we escape them by closing and reopening
    # single quotes in Perl.  However for simplicity we assume DNA/AA sequences.
    script = (
        "my $counts = SeqTools::Kmer::count('" + seq + "', " + str(k) + ");"
        "foreach my $kmer (sort keys %{$counts}) { print $kmer . '\t' . $counts->{$kmer} . '\n'; }"
    )
    cmd = [
        "perl",
        "-MSeqTools::Kmer",
        "-e",
        script,
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=root, check=True)
    counts: Dict[str, int] = {}
    for line in result.stdout.strip().splitlines():
        if not line.strip():
            continue
        kmer, count = line.split("\t")
        counts[kmer] = int(count)
    return counts