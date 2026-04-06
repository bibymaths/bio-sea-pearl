"""Tests for Markov chain random walk generation."""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest

from bio_sea_pearl.api import generate_walk

REPO_ROOT = Path(__file__).resolve().parents[2]
MARKOV_SCRIPT = REPO_ROOT / "markov" / "bin" / "randomwalk.pl"

_skip_no_perl = pytest.mark.skipif(
    shutil.which("perl") is None,
    reason="Perl interpreter not available",
)
_skip_no_script = pytest.mark.skipif(
    not MARKOV_SCRIPT.exists(),
    reason="markov/bin/randomwalk.pl not found",
)

pytestmark = [_skip_no_perl, _skip_no_script]


def test_generate_walk_basic(tmp_fasta):
    """A short walk should return a non-empty string."""
    walk = generate_walk(str(tmp_fasta), 10, start="A", order=1, method="alias")
    walk = walk.strip()
    assert len(walk) > 0


def test_generate_walk_length(tmp_fasta):
    """The generated walk should have approximately the requested length."""
    length = 20
    walk = generate_walk(str(tmp_fasta), length, start="A", order=1, method="alias")
    walk = walk.strip()
    assert length <= len(walk) <= length + 1


def test_generate_walk_characters(tmp_fasta):
    """The walk should only contain characters from the training alphabet."""
    walk = generate_walk(str(tmp_fasta), 50, start="A", order=1, method="alias")
    walk = walk.strip()
    allowed = set("ACGT")
    assert set(walk).issubset(allowed), f"Unexpected characters in walk: {set(walk) - allowed}"


def test_generate_walk_binsrch_method(tmp_fasta):
    """The binsrch sampling method should also produce a valid walk."""
    walk = generate_walk(str(tmp_fasta), 10, start="A", order=1, method="binsrch")
    walk = walk.strip()
    assert len(walk) > 0
