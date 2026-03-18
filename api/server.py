"""FastAPI service exposing Bio Sea Pearl functionality.

This module defines a FastAPI application that wraps the toolkit's API
functions in REST endpoints.  The server can be launched with Uvicorn,
for example:

    uvicorn api.server:app --reload

Endpoints:

* POST `/align` – align two sequences given FASTA file contents.
* POST `/markov` – generate a random walk from a Markov model.
* POST `/distance` – compute Hamming or Levenshtein distance between two strings.
* POST `/kmer` – compute k‑mer counts of a sequence.
* POST `/bwt/search` – search for a pattern within a sequence via the FM‑index.
"""

from __future__ import annotations

from fastapi import FastAPI, HTTPException
from pydantic import BaseModel

from bio_sea_pearl.api import (
    align_sequences,
    generate_walk,
    hamming_distance,
    levenshtein_distance,
    kmer_counts,
    build_fm_index,
    search_fm_index,
)

app = FastAPI(title="Bio Sea Pearl API", version="0.1.0")

@app.get("/")
def read_root():
    return {"message": "Welcome to the Bio Sea Pearl API! Go to /docs to test the endpoints."}

class AlignRequest(BaseModel):
    fasta1: str
    fasta2: str
    matrix: str | None = None
    mode: str = "global"


class MarkovRequest(BaseModel):
    fasta: str
    length: int
    start: str = "A"
    order: int = 1
    method: str = "alias"
    pseudocount: int = 0


class DistanceRequest(BaseModel):
    seq1: str
    seq2: str
    metric: str = "hamming"  # or 'levenshtein'


class KmerRequest(BaseModel):
    sequence: str
    k: int


class BWTRequest(BaseModel):
    sequence: str
    pattern: str


@app.post("/align")
def align_endpoint(req: AlignRequest) -> dict:
    try:
        result = align_sequences(req.fasta1, req.fasta2, matrix=req.matrix, mode=req.mode)
        return {"result": result.strip()}
    except Exception as exc:
        raise HTTPException(status_code=500, detail=str(exc))


@app.post("/markov")
def markov_endpoint(req: MarkovRequest) -> dict:
    try:
        walk = generate_walk(req.fasta, req.length, start=req.start, order=req.order, method=req.method,
                             pseudocount=req.pseudocount)
        return {"walk": walk.strip()}
    except Exception as exc:
        raise HTTPException(status_code=500, detail=str(exc))


@app.post("/distance")
def distance_endpoint(req: DistanceRequest) -> dict:
    metric = req.metric.lower()
    if metric == "hamming":
        func = hamming_distance
    elif metric == "levenshtein":
        func = levenshtein_distance
    else:
        raise HTTPException(status_code=400, detail=f"Unsupported metric '{req.metric}'")
    try:
        d = func(req.seq1, req.seq2)
        return {"distance": d}
    except Exception as exc:
        raise HTTPException(status_code=500, detail=str(exc))


@app.post("/kmer")
def kmer_endpoint(req: KmerRequest) -> dict:
    try:
        counts = kmer_counts(req.sequence, req.k)
        return {"counts": counts}
    except Exception as exc:
        raise HTTPException(status_code=500, detail=str(exc))


@app.post("/bwt/search")
def bwt_search_endpoint(req: BWTRequest) -> dict:
    try:
        index = build_fm_index(req.sequence)
        positions = search_fm_index(index, req.pattern)
        return {"positions": positions}
    except Exception as exc:
        raise HTTPException(status_code=500, detail=str(exc))
