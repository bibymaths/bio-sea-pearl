"""Integration tests for the FastAPI REST endpoints."""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

from api.server import app

REPO_ROOT = Path(__file__).resolve().parents[1]
ALIGN_SCRIPT = REPO_ROOT / "alignment" / "bin" / "align.py"
MARKOV_SCRIPT = REPO_ROOT / "markov" / "bin" / "randomwalk.pl"

client = TestClient(app)

_skip_no_align = pytest.mark.skipif(
    not ALIGN_SCRIPT.exists(),
    reason="alignment/bin/align.py not found",
)
_skip_no_markov = pytest.mark.skipif(
    shutil.which("perl") is None or not MARKOV_SCRIPT.exists(),
    reason="Perl or markov/bin/randomwalk.pl not available",
)


# ── Root endpoint ─────────────────────────────────────────────────

def test_root_returns_welcome():
    resp = client.get("/")
    assert resp.status_code == 200
    data = resp.json()
    assert "message" in data
    assert "Welcome" in data["message"]


# ── Distance endpoint ────────────────────────────────────────────

def test_distance_hamming():
    resp = client.post("/distance", json={"seq1": "ACGT", "seq2": "AGGT", "metric": "hamming"})
    assert resp.status_code == 200
    assert resp.json()["distance"] == 1


def test_distance_levenshtein():
    resp = client.post("/distance", json={"seq1": "kitten", "seq2": "sitting", "metric": "levenshtein"})
    assert resp.status_code == 200
    assert resp.json()["distance"] == 3


def test_distance_unsupported_metric():
    resp = client.post("/distance", json={"seq1": "A", "seq2": "A", "metric": "jaccard"})
    assert resp.status_code == 400
    assert "Unsupported metric" in resp.json()["detail"]


def test_distance_default_metric_is_hamming():
    resp = client.post("/distance", json={"seq1": "AAAA", "seq2": "AAAA"})
    assert resp.status_code == 200
    assert resp.json()["distance"] == 0


# ── K-mer endpoint ───────────────────────────────────────────────

def test_kmer_basic():
    resp = client.post("/kmer", json={"sequence": "ACGTAC", "k": 2})
    assert resp.status_code == 200
    counts = resp.json()["counts"]
    assert counts["AC"] == 2
    assert counts["CG"] == 1


def test_kmer_k_equals_length():
    resp = client.post("/kmer", json={"sequence": "ACG", "k": 3})
    assert resp.status_code == 200
    assert resp.json()["counts"] == {"ACG": 1}


def test_kmer_invalid_k():
    resp = client.post("/kmer", json={"sequence": "AC", "k": 5})
    assert resp.status_code == 400
    assert "detail" in resp.json()


# ── BWT / FM-index endpoint ──────────────────────────────────────

def test_bwt_search_banana_ana():
    resp = client.post("/bwt/search", json={"sequence": "banana", "pattern": "ana"})
    assert resp.status_code == 200
    assert sorted(resp.json()["positions"]) == [1, 3]


def test_bwt_search_not_found():
    resp = client.post("/bwt/search", json={"sequence": "banana", "pattern": "xyz"})
    assert resp.status_code == 200
    assert resp.json()["positions"] == []


def test_bwt_search_dna():
    resp = client.post("/bwt/search", json={"sequence": "ACGTACGT", "pattern": "CGT"})
    assert resp.status_code == 200
    assert sorted(resp.json()["positions"]) == [1, 5]


# ── Alignment endpoint ──────────────────────────────────────────

@_skip_no_align
def test_align_endpoint_identical():
    resp = client.post("/align", json={
        "fasta1": ">seq\nACGT\n",
        "fasta2": ">seq\nACGT\n",
        "mode": "global",
    })
    assert resp.status_code == 200
    assert "result" in resp.json()
    assert len(resp.json()["result"]) > 0


@_skip_no_align
def test_align_endpoint_different_sequences():
    resp = client.post("/align", json={
        "fasta1": ">seqA\nACGTACGT\n",
        "fasta2": ">seqB\nTGCATGCA\n",
        "mode": "global",
    })
    assert resp.status_code == 200
    assert resp.json()["result"]


@_skip_no_align
def test_align_endpoint_with_matrix():
    resp = client.post("/align", json={
        "fasta1": ">seq\nACGT\n",
        "fasta2": ">seq\nACGT\n",
        "matrix": "alignment/scoring/blosum62.mat",
        "mode": "global",
    })
    assert resp.status_code == 200


# ── Markov endpoint ──────────────────────────────────────────────

@_skip_no_markov
def test_markov_endpoint_basic():
    resp = client.post("/markov", json={
        "fasta": ">seq\nACGTACGT\n",
        "length": 10,
        "start": "A",
        "order": 1,
    })
    assert resp.status_code == 200
    walk = resp.json()["walk"]
    assert len(walk) > 0


@_skip_no_markov
def test_markov_endpoint_characters():
    resp = client.post("/markov", json={
        "fasta": ">seq\nACGTACGT\n",
        "length": 30,
        "start": "A",
        "order": 1,
    })
    assert resp.status_code == 200
    walk = resp.json()["walk"]
    assert set(walk).issubset(set("ACGT"))


# ── Error handling ───────────────────────────────────────────────

def test_missing_required_field_distance():
    resp = client.post("/distance", json={"seq1": "A"})
    assert resp.status_code == 422


def test_missing_required_field_kmer():
    resp = client.post("/kmer", json={"sequence": "ACGT"})
    assert resp.status_code == 422


def test_missing_required_field_bwt():
    resp = client.post("/bwt/search", json={"sequence": "ACGT"})
    assert resp.status_code == 422
