"""Shared fixtures for the Bio Sea Pearl test suite."""

from __future__ import annotations

from pathlib import Path

import pytest


@pytest.fixture()
def tmp_fasta(tmp_path: Path) -> Path:
    """Create a small FASTA file and return its path."""
    fasta = tmp_path / "test.fa"
    fasta.write_text(">seq1\nACGTACGT\n")
    return fasta


@pytest.fixture()
def tmp_fasta_pair(tmp_path: Path) -> tuple[Path, Path]:
    """Create a pair of identical FASTA files for alignment tests."""
    content = ">seq\nACGT\n"
    f1 = tmp_path / "seq1.fa"
    f2 = tmp_path / "seq2.fa"
    f1.write_text(content)
    f2.write_text(content)
    return f1, f2
