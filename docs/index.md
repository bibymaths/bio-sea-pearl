---
title: Bio Sea Pearl
---

<div style="text-align: center;">

<p>
<a href="https://doi.org/10.5281/zenodo.19435099"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.19435099.svg" alt="DOI"></a>
<a href="https://pypi.org/project/bio-sea-pearl/"><img src="https://img.shields.io/pypi/v/bio-sea-pearl.svg" alt="PyPI version"></a>
<a href="https://pypi.org/project/bio-sea-pearl/"><img src="https://img.shields.io/pypi/pyversions/bio-sea-pearl.svg" alt="Python versions"></a>
<a href="https://github.com/bibymaths/bio-sea-pearl/releases"><img src="https://img.shields.io/github/v/release/bibymaths/bio-sea-pearl" alt="Release"></a>
<a href="https://github.com/bibymaths/bio-sea-pearl/blob/main/LICENSE"><img src="https://img.shields.io/github/license/bibymaths/bio-sea-pearl" alt="License"></a>
<a href="https://github.com/bibymaths/bio-sea-pearl/pkgs/container/bio-sea-pearl"><img src="https://img.shields.io/badge/container-ghcr.io-blue" alt="Container"></a>
</p>

</div>

# Bio Sea Pearl

**Sequence analysis and bioinformatics utilities in Python and Perl.**

Bio Sea Pearl is a dual-language bioinformatics toolkit that integrates mature Perl implementations with a modern Python interface. It provides pairwise sequence alignment, Markov chain simulation, sequence distance metrics, k-mer analysis, and full-text indexing via the Burrows–Wheeler Transform — accessible through a unified CLI, a REST API, or direct Python imports.

---

## Core capabilities

| Subsystem | What it does | Implementation |
|-----------|-------------|----------------|
| [**Alignment**](modules/alignment.md) | Global, local, and LCS pairwise alignment with affine gap penalties | Python + Perl (Gotoh algorithm) |
| [**Markov Chains**](modules/markov.md) | Train transition models on sequences and sample random walks | Perl (first-order and higher-order) |
| [**Sequence Tools**](modules/seqtools.md) | Hamming distance, Levenshtein distance, k-mer counting, pattern search | Python (native) + Perl (Inline C) |
| [**BWT & FM-index**](modules/bwt.md) | Suffix arrays, Burrows–Wheeler Transform, FM-index substring search | Pure Python |

---

## Three ways to use it

=== "CLI"

    ```bash
    biosea seqtools hamming ACGTACGT ACGTACGA
    biosea bwt search --sequence ACGTACGT --pattern CGT
    biosea align seq1.fa seq2.fa --mode global
    ```

=== "REST API"

    ```bash
    curl -X POST http://localhost:8000/distance \
        -H "Content-Type: application/json" \
        -d '{"seq1": "kitten", "seq2": "sitting", "metric": "levenshtein"}'
    ```

=== "Python"

    ```python
    from bio_sea_pearl.api import hamming_distance, build_fm_index, search_fm_index

    print(hamming_distance("ACGT", "AGGT"))  # 1

    idx = build_fm_index("ACGTACGT")
    print(search_fm_index(idx, "CGT"))       # [1, 5]
    ```

---

## Architecture at a glance

```
              ┌──────────┐    ┌──────────────┐
              │  biosea  │    │ FastAPI       │
              │   CLI    │    │ REST server   │
              └────┬─────┘    └──────┬────────┘
                   │                 │
                   └────────┬────────┘
                            │
                   ┌────────▼────────┐
                   │   API layer     │
                   │ (bio_sea_pearl) │
                   └───┬─────────┬───┘
                       │         │
              ┌────────▼──┐  ┌───▼──────────┐
              │  Wrappers │  │ Pure Python   │
              │ (Perl via │  │ (BWT, seqtools│
              │ subprocess)│  │  distances)   │
              └────┬──────┘  └──────────────┘
                   │
          ┌────────▼────────┐
          │  Perl scripts   │
          │  & modules      │
          │ (alignment,     │
          │  markov,        │
          │  seqtools)      │
          └─────────────────┘
```

The Python API layer is the single dispatch point. Whether a request arrives via the CLI, the REST API, or a direct import, it flows through the same API functions in `src/bio_sea_pearl/api/`. See [Architecture](architecture.md) for the full breakdown.

---

## Quick start

```bash
pip install bio-sea-pearl
biosea --help
```

:material-arrow-right: [Installation guide](getting-started/installation.md) · [Quickstart tutorial](getting-started/quickstart.md)

---

## Documentation map

| Section | Description |
|---------|-------------|
| [Installation](getting-started/installation.md) | Install from PyPI, source, or Docker |
| [Quickstart](getting-started/quickstart.md) | First commands across all subsystems |
| [CLI Reference](cli.md) | Complete `biosea` command documentation |
| [REST API](api.md) | Endpoint schemas and interactive docs |
| [Architecture](architecture.md) | Layers, dispatch, data flow, repo layout |
| [Alignment](modules/alignment.md) | Gotoh DP, scoring matrices, parallelisation |
| [Markov Chains](modules/markov.md) | Transition models, sampling, random walks |
| [Sequence Tools](modules/seqtools.md) | Distances, k-mers, Boyer–Moore search |
| [BWT & FM-index](modules/bwt.md) | Suffix arrays, BWT, backward search |
| [Extending](developer/extending.md) | Adding modules, porting Perl to Python |
| [Troubleshooting](developer/troubleshooting.md) | Common errors and failure modes |

