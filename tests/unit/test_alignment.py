"""Tests for sequence alignment via the public API."""

from __future__ import annotations

from pathlib import Path

import pytest

from bio_sea_pearl.api import align_sequences

REPO_ROOT = Path(__file__).resolve().parents[2]
ALIGN_SCRIPT = REPO_ROOT / "alignment" / "bin" / "align.py"
BLOSUM62 = REPO_ROOT / "alignment" / "scoring" / "blosum62.mat"

_skip_no_align = pytest.mark.skipif(
    not ALIGN_SCRIPT.exists(),
    reason="alignment/bin/align.py not found",
)


@_skip_no_align
def test_align_identical_sequences(tmp_fasta_pair):
    """Identical sequences should produce non-empty alignment output."""
    f1, f2 = tmp_fasta_pair
    result = align_sequences(str(f1), str(f2))
    assert result, "Expected non-empty alignment output for identical sequences"


@_skip_no_align
def test_align_with_default_matrix(tmp_fasta_pair):
    """Alignment without explicit matrix should use the default blosum62."""
    f1, f2 = tmp_fasta_pair
    result = align_sequences(str(f1), str(f2))
    assert isinstance(result, str)
    assert len(result) > 0


@_skip_no_align
def test_align_with_explicit_matrix(tmp_fasta_pair):
    """Alignment with an explicit matrix path should succeed."""
    f1, f2 = tmp_fasta_pair
    result = align_sequences(str(f1), str(f2), matrix="alignment/scoring/blosum62.mat")
    assert isinstance(result, str)
    assert len(result) > 0


@_skip_no_align
def test_align_different_sequences(tmp_path):
    """Alignment of two distinct sequences should return output."""
    f1 = tmp_path / "a.fa"
    f2 = tmp_path / "b.fa"
    f1.write_text(">seqA\nACGTACGT\n")
    f2.write_text(">seqB\nTGCATGCA\n")
    result = align_sequences(str(f1), str(f2))
    assert result


@_skip_no_align
def test_align_local_mode(tmp_fasta_pair):
    """Local alignment mode should also return output."""
    f1, f2 = tmp_fasta_pair
    result = align_sequences(str(f1), str(f2), mode="local")
    assert isinstance(result, str)
