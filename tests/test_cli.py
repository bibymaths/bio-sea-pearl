"""CLI smoke tests – verify commands exit cleanly and produce expected output."""

from __future__ import annotations

import subprocess
import sys


def _run(*args: str) -> subprocess.CompletedProcess[str]:
    """Run a CLI command and return the completed process."""
    return subprocess.run(
        [sys.executable, "-m", "bio_sea_pearl.cli", *args],
        capture_output=True,
        text=True,
        timeout=60,
    )


def test_help_flag():
    result = _run("--help")
    assert result.returncode == 0
    assert "Bio Sea Pearl" in result.stdout or "Usage" in result.stdout


def test_seqtools_hamming():
    result = _run("seqtools", "hamming", "ACGT", "AGGT")
    assert result.returncode == 0
    assert result.stdout.strip() == "1"


def test_seqtools_hamming_identical():
    result = _run("seqtools", "hamming", "AAAA", "AAAA")
    assert result.returncode == 0
    assert result.stdout.strip() == "0"


def test_seqtools_levenshtein():
    result = _run("seqtools", "levenshtein", "kitten", "sitting")
    assert result.returncode == 0
    assert result.stdout.strip() == "3"


def test_seqtools_kmer():
    result = _run("seqtools", "kmer", "ACGTAC", "--k", "2")
    assert result.returncode == 0
    import json
    counts = json.loads(result.stdout)
    assert counts["AC"] == 2


def test_bwt_search():
    result = _run("bwt", "search", "--sequence", "ACGTACGT", "--pattern", "CGT")
    assert result.returncode == 0
    positions = list(map(int, result.stdout.strip().split()))
    assert sorted(positions) == [1, 5]


def test_bwt_search_not_found():
    result = _run("bwt", "search", "--sequence", "ACGTACGT", "--pattern", "ZZZ")
    assert result.returncode == 0
    assert result.stdout.strip() == ""
