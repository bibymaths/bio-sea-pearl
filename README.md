<img src="assets/logo.svg" width="700">

A hybrid **Python + Perl bioinformatics toolkit** for sequence alignment, Markov modeling, sequence analysis, and
FM-index–based search — unified under a single Python CLI and API.

---

## Documentation

* 📘 **Live Docs (GitHub Pages)**
  <https://bibymaths.github.io/bio-sea-pearl/>

* 📂 Local Docs
  Located in `docs/`

---

## Quick Start

### 1. Clone

```bash
git clone https://github.com/bibymaths/bio-sea-pearl.git
cd bio-sea-pearl
```

---

### 2. Install Dependencies

This project uses **Python ≥ 3.10** and is built with [Hatch](https://hatch.pypa.io/).

#### With uv (recommended)

```bash
uv pip install -e ".[dev]"
```

#### With pip

```bash
pip install -e ".[dev]"
```

> **Note:** Some features (alignment, Markov generation) delegate to Perl
> scripts at runtime. Install Perl ≥ 5.26 if you need those features.

---

### 3. Run the CLI

The unified CLI entrypoint is:

```bash
biosea --help
```

---

## CLI Usage

### Alignment

```bash
biosea align seq1.fa seq2.fa \
  --matrix alignment/scoring/blosum62.mat \
  --mode global
```

⚠️ Notes:

* Matrix filenames are **case-sensitive**
* Use `blosum62.mat`, not `BLOSUM62.mat`

Optionally, you can generate dotplot in svg format:

```bash
perl alignment/bin/dotplot.pl align.matrix.tsv dotplot.svg
```

---

### Markov Chain

```bash
biosea markov \
  --fasta seq1.fa \
  --length 100 \
  --start A \
  --order 1 \
  --method alias
```

For higher-order models:

```bash
--order 2 --start AA
```

⚠️ Constraint:

* `start` length must equal `order`

---

### Sequence Utilities

```bash
# Hamming distance
biosea seqtools hamming ACGT AGGT

# Levenshtein distance
biosea seqtools levenshtein kitten sitting

# k-mer counts
biosea seqtools kmer ACGTACGT --k 3
```

---

### BWT / FM-Index Search

```bash
biosea bwt search \
  --sequence ACGTACGT \
  --pattern CGT
```

---

## REST API

Start the FastAPI server:

```bash
uvicorn api.server:app --reload
```

Endpoints:

* `POST /align`
* `POST /markov`
* `POST /distance`
* `POST /kmer`
* `POST /bwt/search`

Example:

```bash
curl -X POST http://localhost:8000/distance \
  -H "Content-Type: application/json" \
  -d '{"seq1": "kitten", "seq2": "sitting", "metric": "levenshtein"}'
```

Interactive API documentation is available at <http://localhost:8000/docs>.

---

## Docker

### Build and start

```bash
./docker_up.sh
```

### Stop and remove

```bash
./docker_down.sh
```

### Interactive shell

```bash
./docker_interactive.sh
```

### Manual Docker commands

```bash
docker compose up --build -d
docker compose exec biosea biosea --help
docker compose down
```

---

## Project Structure

```
src/bio_sea_pearl/
├── cli.py                 # Unified CLI
├── api/                   # Clean Python API layer
├── perl_wrappers/         # Bridge to legacy Perl scripts
├── seqtools_py/           # Python ports of core algorithms
└── bwt/                   # Native Python FM-index

alignment/                 # Legacy alignment tools
markov/                    # Perl Markov models
seqtools/                  # Perl sequence utilities
api/server.py              # FastAPI layer
docs/                      # MkDocs documentation
tests/                     # Unit + integration tests
```

---

## Architecture Overview

The system is layered:

```
CLI / API
   ↓
Python API Layer
   ↓
Wrappers (subprocess)
   ↓
Perl + Python legacy tools
```

This design:

* preserves legacy code
* enables gradual Python migration
* provides production-ready interfaces

---

## Running Tests

```bash
pytest
```

---

## Building

```bash
pip install build
python -m build
```

This produces a source distribution and wheel in `dist/`.

---

## Releasing

Releases are automated via GitHub Actions. Push a tag to trigger the workflow:

```bash
git tag v0.1.0
git push origin v0.1.0
```

This will:

1. Run the test suite
2. Create a GitHub release
3. Build and publish the package to PyPI
4. Build and push multi-arch Docker images to `ghcr.io/bibymaths/bio-sea-pearl`

---

## Troubleshooting

### Alignment fails

* Check matrix path:

  ```
  alignment/scoring/blosum62.mat
  ```
* Avoid uppercase filenames

---

### Markov fails

Error:

```
Start state must be length N
```

Fix:

```
--order N → start string length must be N
```

---

### CLI not found

```bash
pip install -e .
biosea --help
```

---

## License

This project is licensed under the MIT License. See `LICENSE` for details.
