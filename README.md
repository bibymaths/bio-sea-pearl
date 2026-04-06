<div align="center">

<img src="assets/logo.svg" width="300">

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19435099.svg)](https://doi.org/10.5281/zenodo.19435099)
[![PyPI version](https://img.shields.io/pypi/v/bio-sea-pearl.svg)](https://pypi.org/project/bio-sea-pearl/)
[![Python versions](https://img.shields.io/pypi/pyversions/bio-sea-pearl.svg)](https://pypi.org/project/bio-sea-pearl/)
[![PyPI downloads](https://img.shields.io/pypi/dm/bio-sea-pearl.svg)](https://pypi.org/project/bio-sea-pearl/)
[![GitHub release](https://img.shields.io/github/v/release/bibymaths/bio-sea-pearl)](https://github.com/bibymaths/bio-sea-pearl/releases)
[![Docs](https://img.shields.io/badge/docs-GitHub%20Pages-blue)](https://bibymaths.github.io/bio-sea-pearl/)
[![License](https://img.shields.io/github/license/bibymaths/bio-sea-pearl)](LICENSE)
[![CI](https://img.shields.io/github/actions/workflow/status/bibymaths/bio-sea-pearl/release.yml?branch=main&label=release)](https://github.com/bibymaths/bio-sea-pearl/actions/workflows/release.yml)
[![Docs workflow](https://img.shields.io/github/actions/workflow/status/bibymaths/bio-sea-pearl/docs.yml?branch=main&label=docs)](https://github.com/bibymaths/bio-sea-pearl/actions/workflows/docs.yml)
[![GHCR](https://img.shields.io/badge/container-ghcr.io-blue)](https://github.com/bibymaths/bio-sea-pearl/pkgs/container/bio-sea-pearl)
[![Docker pulls](https://img.shields.io/badge/docker-GHCR-blue)](https://github.com/bibymaths/bio-sea-pearl/pkgs/container/bio-sea-pearl)


</div>

***

A hybrid **Python + Perl bioinformatics toolkit** for sequence alignment, Markov modeling, sequence analysis, and FM-index-based search — unified under a single Python CLI and API.

***

[Read Documentation](https://bibymaths.github.io/bio-sea-pearl/)

***

## Quick Start

Requires **Python >= 3.10** and optionally **Perl >= 5.26** (for alignment and Markov features).

### 1. Install

```bash
git clone https://github.com/bibymaths/bio-sea-pearl.git
cd bio-sea-pearl
uv pip install -e ".[dev]"   # or: pip install -e ".[dev]"
```

### 2. Run

```bash
biosea --help
```

***

## CLI Usage

### Alignment

```bash
biosea align seq1.fa seq2.fa --matrix alignment/scoring/blosum62.mat --mode global
```

> Matrix filenames are case-sensitive (`blosum62.mat`, not `BLOSUM62.mat`). For a dotplot: `perl alignment/bin/dotplot.pl align.matrix.tsv dotplot.svg`

### Markov Chain

```bash
biosea markov --fasta seq1.fa --length 100 --start A --order 1 --method alias
```

> `--start` length must equal `--order` (e.g. `--order 2 --start AA`)

### Sequence Utilities

```bash
biosea seqtools hamming ACGT AGGT
biosea seqtools levenshtein kitten sitting
biosea seqtools kmer ACGTACGT --k 3
```

### BWT / FM-Index

```bash
biosea bwt search --sequence ACGTACGT --pattern CGT
```

***

## REST API

```bash
uvicorn api.server:app --reload
```

Endpoints: `POST /align`, `/markov`, `/distance`, `/kmer`, `/bwt/search`. Interactive docs at [http://localhost:8000/docs](http://localhost:8000/docs).

```bash
curl -X POST http://localhost:8000/distance \
  -H "Content-Type: application/json" \
  -d '{"seq1": "kitten", "seq2": "sitting", "metric": "levenshtein"}'
```

***

## Docker

```bash
./docker_up.sh          # Build and start
./docker_down.sh        # Stop and remove
./docker_interactive.sh # Interactive shell
```

***

## Project Structure

```
src/bio_sea_pearl/
├── cli.py               # Unified CLI
├── api/                 # Python API layer
├── perl_wrappers/       # Perl subprocess bridge
├── seqtools_py/         # Python algorithm ports
└── bwt/                 # FM-index (pure Python)

alignment/  markov/  seqtools/  api/server.py  docs/  tests/
```

The layered architecture (`CLI/API → Python layer → subprocess wrappers → Perl/Python tools`) preserves legacy code while enabling gradual migration.

***

## Build, Test & Release

```bash
pytest                        # Run tests
python -m build               # Build wheel + sdist into dist/
git tag v0.1.0 && git push origin v0.1.0  # Trigger release CI
```

Pushing a tag runs tests, creates a GitHub release, publishes to PyPI, and pushes multi-arch Docker images to `ghcr.io/bibymaths/bio-sea-pearl`.

***

## Troubleshooting

- **Alignment fails:** verify matrix path (`alignment/scoring/blosum62.mat`) and use lowercase filenames.
- **Markov fails** (`Start state must be length N`): ensure `--start` string length equals `--order`.
- **CLI not found:** run `pip install -e .` then retry.

***

## License

[MIT License](LICENSE)

***