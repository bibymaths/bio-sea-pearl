"""FastAPI service exposing Bio Sea Pearl functionality.

This module defines a FastAPI application that wraps the toolkit's API
functions in REST endpoints.  The server can be launched with Uvicorn,
for example:

    uvicorn api.server:app --reload

Endpoints:

* POST ``/align``  – align two sequences given FASTA file contents.
* POST ``/markov`` – generate a random walk from a Markov model.
* POST ``/distance`` – compute Hamming or Levenshtein distance between two strings.
* POST ``/kmer``   – compute k-mer counts of a sequence.
* POST ``/bwt/search`` – search for a pattern within a sequence via the FM-index.
"""

from __future__ import annotations

import logging
import tempfile
from pathlib import Path

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

logger = logging.getLogger(__name__)

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
    metric: str = "hamming"


class KmerRequest(BaseModel):
    sequence: str
    k: int


class BWTRequest(BaseModel):
    sequence: str
    pattern: str


def _write_temp_fasta(content: str, suffix: str = ".fa") -> Path:
    """Write FASTA content to a named temporary file and return its path.

    The caller is responsible for cleanup.
    """
    tmp = tempfile.NamedTemporaryFile(mode="w", suffix=suffix, delete=False)
    tmp.write(content)
    tmp.flush()
    tmp.close()
    return Path(tmp.name)


@app.post("/align")
def align_endpoint(req: AlignRequest) -> dict:
    tmp1: Path | None = None
    tmp2: Path | None = None
    try:
        if not req.fasta1.strip().startswith(">") or not req.fasta2.strip().startswith(">"):
            raise HTTPException(400, "Expected FASTA content (must start with '>')")
        tmp1 = _write_temp_fasta(req.fasta1)
        tmp2 = _write_temp_fasta(req.fasta2)
        result = align_sequences(str(tmp1), str(tmp2), matrix=req.matrix, mode=req.mode)
        return {"result": result.strip()}
    except HTTPException:
        raise
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc
    except FileNotFoundError as exc:
        raise HTTPException(status_code=404, detail=str(exc)) from exc
    except Exception as exc:
        logger.exception("Unexpected error in /align")
        raise HTTPException(status_code=500, detail="Internal server error") from exc
    finally:
        for p in (tmp1, tmp2):
            if p is not None:
                p.unlink(missing_ok=True)


@app.post("/markov")
def markov_endpoint(req: MarkovRequest) -> dict:
    tmp: Path | None = None
    try:
        tmp = _write_temp_fasta(req.fasta)
        if not req.fasta.strip().startswith(">"):
            raise HTTPException(400, "Expected FASTA content (must start with '>')")
        walk = generate_walk(str(tmp), req.length, start=req.start, order=req.order, method=req.method,
                             pseudocount=req.pseudocount)
        return {"walk": walk.strip()}
    except HTTPException:
        raise
    except (ValueError, TypeError) as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc
    except FileNotFoundError as exc:
        raise HTTPException(status_code=404, detail=str(exc)) from exc
    except Exception as exc:
        logger.exception("Unexpected error in /markov")
        raise HTTPException(status_code=500, detail="Internal server error") from exc
    finally:
        if tmp is not None:
            tmp.unlink(missing_ok=True)


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
    except HTTPException:
        raise
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc
    except Exception as exc:
        logger.exception("Unexpected error in /distance")
        raise HTTPException(status_code=500, detail="Internal server error") from exc


@app.post("/kmer")
def kmer_endpoint(req: KmerRequest) -> dict:
    try:
        counts = kmer_counts(req.sequence, req.k)
        return {"counts": counts}
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc
    except Exception as exc:
        logger.exception("Unexpected error in /kmer")
        raise HTTPException(status_code=500, detail="Internal server error") from exc


@app.post("/bwt/search")
def bwt_search_endpoint(req: BWTRequest) -> dict:
    try:
        index = build_fm_index(req.sequence)
        positions = search_fm_index(index, req.pattern)
        return {"positions": positions}
    except Exception as exc:
        logger.exception("Unexpected error in /bwt/search")
        raise HTTPException(status_code=500, detail="Internal server error") from exc
