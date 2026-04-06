---
title: REST API
---

# REST API

Bio Sea Pearl exposes its functionality as a REST API via [FastAPI](https://fastapi.tiangolo.com/). The server is defined in `api/server.py` and serves five endpoints that mirror the capabilities of the `biosea` CLI.

---

## Starting the server

=== "Direct"

    ```bash
    uvicorn api.server:app --host 0.0.0.0 --port 8000
    ```

=== "Docker"

    ```bash
    docker compose up -d
    ```

=== "Docker (pull)"

    ```bash
    docker run -p 8000:8000 ghcr.io/bibymaths/bio-sea-pearl:latest
    ```

Once running, the interactive Swagger UI is available at `http://localhost:8000/docs` and the ReDoc UI at `http://localhost:8000/redoc`.

---

## Interactive API documentation

<swagger-ui src="openapi.json" />

!!! info
    The OpenAPI schema above is generated from the FastAPI application during CI and embedded here. If you are viewing this locally and the schema does not render, start the server and visit `http://localhost:8000/docs` instead.

---

## Endpoints

### `GET /`

Health check / welcome endpoint.

**Response:**

```json
{"message": "Welcome to Bio Sea Pearl API"}
```

---

### `POST /align`

Align two sequences provided as FASTA-format strings.

**Request body:**

| Field   | Type   | Required | Default    | Description |
|---------|--------|----------|------------|-------------|
| `fasta1` | string | yes | — | FASTA content for the first sequence |
| `fasta2` | string | yes | — | FASTA content for the second sequence |
| `matrix` | string | no | *(built-in)* | Scoring matrix name (e.g. `BLOSUM62.mat`) |
| `mode`   | string | no | `"global"` | Alignment mode: `global`, `local`, or `lcs` |

**Example request:**

```bash
curl -s -X POST http://localhost:8000/align \
    -H "Content-Type: application/json" \
    -d '{
        "fasta1": ">seq1\nACGTACGT",
        "fasta2": ">seq2\nACGTACGA",
        "matrix": "BLOSUM62.mat",
        "mode": "global"
    }'
```

!!! note
    The server writes the FASTA content to temporary files, passes them to the alignment wrapper, and cleans up the files after the response is returned. This is a **content-based** endpoint — you send sequence data, not file paths.

---

### `POST /markov`

Generate a random walk from a Markov chain trained on the provided FASTA sequences.

**Request body:**

| Field        | Type    | Required | Default    | Description |
|--------------|---------|----------|------------|-------------|
| `fasta`      | string  | yes      | —          | FASTA content for training |
| `length`     | integer | yes      | —          | Walk length (characters to generate) |
| `start`      | string  | no       | `"A"`      | Starting state |
| `order`      | integer | no       | `1`        | Markov chain order |
| `method`     | string  | no       | `"alias"`  | Sampling method: `alias` or `binsrch` |
| `pseudocount`| number  | no       | `0`        | Smoothing pseudocount |

**Example request:**

```bash
curl -s -X POST http://localhost:8000/markov \
    -H "Content-Type: application/json" \
    -d '{
        "fasta": ">train\nACGTACGTACGTACGT",
        "length": 20,
        "start": "A",
        "order": 1,
        "method": "alias"
    }'
```

---

### `POST /distance`

Compute the Hamming or Levenshtein distance between two strings.

**Request body:**

| Field   | Type   | Required | Default      | Description |
|---------|--------|----------|--------------|-------------|
| `seq1`  | string | yes      | —            | First sequence |
| `seq2`  | string | yes      | —            | Second sequence |
| `metric`| string | no       | `"hamming"`  | Distance metric: `hamming` or `levenshtein` |

**Example request:**

```bash
curl -s -X POST http://localhost:8000/distance \
    -H "Content-Type: application/json" \
    -d '{"seq1": "kitten", "seq2": "sitting", "metric": "levenshtein"}'
```

**Example response:**

```json
{
    "metric": "levenshtein",
    "seq1": "kitten",
    "seq2": "sitting",
    "distance": 3
}
```

---

### `POST /kmer`

Count k-mer frequencies in a sequence.

**Request body:**

| Field     | Type    | Required | Description |
|-----------|---------|----------|-------------|
| `sequence`| string  | yes      | Input sequence |
| `k`       | integer | yes      | K-mer length |

**Example request:**

```bash
curl -s -X POST http://localhost:8000/kmer \
    -H "Content-Type: application/json" \
    -d '{"sequence": "ACGTACGT", "k": 3}'
```

**Example response:**

```json
{
    "sequence": "ACGTACGT",
    "k": 3,
    "kmers": {"ACG": 2, "CGT": 2, "GTA": 1, "TAC": 1}
}
```

---

### `POST /bwt/search`

Search for a pattern in a sequence using the FM-index.

**Request body:**

| Field     | Type   | Required | Description |
|-----------|--------|----------|-------------|
| `sequence`| string | yes      | Text to index and search |
| `pattern` | string | yes      | Pattern to find |

**Example request:**

```bash
curl -s -X POST http://localhost:8000/bwt/search \
    -H "Content-Type: application/json" \
    -d '{"sequence": "ACGTACGT", "pattern": "CGT"}'
```

**Example response:**

```json
{
    "sequence": "ACGTACGT",
    "pattern": "CGT",
    "positions": [1, 5]
}
```

Positions are zero-indexed. An empty list is returned if the pattern is not found.

---

## Error handling

The API returns standard HTTP status codes:

| Status | Meaning |
|--------|---------|
| `200`  | Success |
| `400`  | Validation error (e.g. sequences of different length for Hamming) |
| `404`  | Referenced file not found (e.g. missing scoring matrix) |
| `500`  | Internal server error (Perl not available, subprocess failure, etc.) |

Error responses include a `detail` field with a human-readable message:

```json
{"detail": "Hamming distance requires equal-length sequences"}
```

---

## Content-based vs path-based

All endpoints accept **content** (sequence strings or FASTA text), not file paths. The `/align` and `/markov` endpoints internally write FASTA content to temporary files because the underlying scripts expect file arguments. These temporary files are cleaned up after each request.

This design means:

- Clients do not need filesystem access to the server.
- The server is stateless between requests.
- Docker deployment works without volume mounts for sequence data.

---

## Related pages

- [CLI Reference](cli.md) — equivalent functionality from the command line
- [Architecture](architecture.md) — how the server dispatches to the API layer
- [Quickstart](getting-started/quickstart.md) — basic API usage examples