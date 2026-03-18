---
title: REST API
---

# REST API

The project includes an optional FastAPI service which exposes the toolkit as
a web API.  The server is defined in `api/server.py` and can be launched
using Uvicorn:

```
uvicorn api.server:app --reload
```

## Endpoints

All endpoints accept JSON payloads and return JSON responses.  Errors are
reported with appropriate HTTP status codes.

### `POST /align`

Align two sequences from FASTA files.  The request body:

```json
{
  "fasta1": ">seq1\nACGT\n",
  "fasta2": ">seq2\nAGGT\n",
  "matrix": "alignment/scoring/BLOSUM62.mat",
  "mode": "global"
}
```

Returns:

```json
{
  "result": "<alignment output>"
}
```

### `POST /markov`

Generate a random walk from sequences in a FASTA file.  Request:

```json
{
  "fasta": ">train\nACGTACGT\n",
  "length": 50,
  "start": "A",
  "order": 1,
  "method": "alias",
  "pseudocount": 0
}
```

Returns:

```json
{
  "walk": "ACGTACGT..."
}
```

### `POST /distance`

Compute the Hamming or Levenshtein distance between two strings.  Set
`metric` to `hamming` or `levenshtein`.

```json
{
  "seq1": "kitten",
  "seq2": "sitting",
  "metric": "levenshtein"
}
```

Returns:

```json
{
  "distance": 3
}
```

### `POST /kmer`

Compute k‑mer counts for a sequence.  Request:

```json
{
  "sequence": "ACGTACGT",
  "k": 2
}
```

Returns a mapping of k‑mers to counts:

```json
{
  "counts": {
    "AC": 2,
    "CG": 1,
    "GT": 1,
    "TA": 1
  }
}
```

### `POST /bwt/search`

Search for a pattern using an FM‑index.  The sequence is indexed on the fly.

```json
{
  "sequence": "banana",
  "pattern": "ana"
}
```

Returns:

```json
{
  "positions": [1, 3]
}
```