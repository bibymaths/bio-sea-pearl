---
title: Architecture
---

# Architecture

Bio Sea Pearl is a layered system that wraps legacy Perl bioinformatics code in a modern Python interface. This page describes the layers, how they interact, and where each component lives in the repository.

---

## High-level layers

```
┌─────────────────────────────────────────────────────────────────┐
│                      User-facing interfaces                     │
│                                                                 │
│  ┌──────────────┐   ┌────────────────┐   ┌──────────────────┐  │
│  │ biosea CLI   │   │ FastAPI server │   │ Python imports   │  │
│  │ (Typer)      │   │ (api/server.py)│   │ (library usage)  │  │
│  └──────┬───────┘   └───────┬────────┘   └────────┬─────────┘  │
│         │                   │                     │             │
│         └───────────────────┼─────────────────────┘             │
│                             │                                   │
│                    ┌────────▼─────────┐                         │
│                    │    API layer     │                         │
│                    │ src/bio_sea_pearl│                         │
│                    │    /api/         │                         │
│                    └───┬─────────┬───┘                         │
│                        │         │                              │
│               ┌────────▼──┐  ┌───▼──────────────┐              │
│               │ Perl      │  │ Python-native    │              │
│               │ wrappers  │  │ implementations  │              │
│               └─────┬─────┘  └──────────────────┘              │
│                     │                                           │
│            ┌────────▼──────────┐                                │
│            │  Legacy scripts   │                                │
│            │  alignment/       │                                │
│            │  markov/          │                                │
│            │  seqtools/        │                                │
│            └───────────────────┘                                │
└─────────────────────────────────────────────────────────────────┘
```

### Interface layer

Three entry points exist, and all converge on the same API functions:

| Entry point | Technology | Location |
|-------------|-----------|----------|
| `biosea` CLI | Typer | `src/bio_sea_pearl/cli.py` |
| REST API | FastAPI + Uvicorn | `api/server.py` |
| Library import | Direct Python | `from bio_sea_pearl.api import ...` |

### API layer

The `src/bio_sea_pearl/api/` package exports seven public functions:

```python
from bio_sea_pearl.api import (
    align_sequences,      # alignment_api.py
    generate_walk,        # markov_api.py
    hamming_distance,     # seqtools_api.py
    levenshtein_distance, # seqtools_api.py
    kmer_counts,          # seqtools_api.py
    build_fm_index,       # bwt_api.py
    search_fm_index,      # bwt_api.py
)
```

Each function encapsulates the decision of whether to call a pure-Python implementation or delegate to a Perl subprocess.

### Wrapper layer

`src/bio_sea_pearl/perl_wrappers/` contains bridge functions that invoke Perl scripts via `subprocess.run()`. The wrappers:

1. Locate the repository root using the `BIOSEA_REPO_ROOT` environment variable, or by traversing three parent directories above the wrapper file.
2. Construct the appropriate `perl` or `python` command line.
3. Execute the subprocess, capture stdout, and strip ANSI escape codes from the output.
4. Parse the result into Python data types.

### Native Python layer

Two subsystems have pure Python implementations that do not require Perl:

- **BWT / FM-index** (`src/bio_sea_pearl/bwt/`) — suffix arrays, BWT construction, and FM-index backward search.
- **SeqTools ports** (`src/bio_sea_pearl/seqtools_py/`) — Hamming distance, Levenshtein distance, and k-mer counting.

### Legacy Perl layer

The `alignment/`, `markov/`, and `seqtools/` directories contain the original Perl (and in alignment's case, Python) scripts. These are standalone programs with their own CLI arguments and module systems. The Python wrapper layer calls them as subprocesses.

---

## Repository layout

```
bio-sea-pearl/
├── src/bio_sea_pearl/           # Installable Python package
│   ├── cli.py                   # Typer CLI entrypoint
│   ├── api/                     # Public API functions
│   │   ├── alignment_api.py     #   → align_sequences()
│   │   ├── markov_api.py        #   → generate_walk()
│   │   ├── seqtools_api.py      #   → hamming, levenshtein, kmer
│   │   └── bwt_api.py           #   → build_fm_index, search_fm_index
│   ├── bwt/                     # FM-index (pure Python)
│   │   ├── transform.py         #   FMIndex class, suffix_array(), bwt()
│   │   └── find.py              #   Example search script
│   ├── seqtools_py/             # Python ports of Perl SeqTools
│   │   ├── hamming.py
│   │   ├── levenshtein.py
│   │   └── kmer.py
│   └── perl_wrappers/           # Subprocess bridges to Perl
│       ├── alignment_wrapper.py
│       ├── markov_wrapper.py
│       ├── seqtools_wrapper.py
│       └── term_utils.py        # ANSI escape stripping
├── api/
│   └── server.py                # FastAPI application
├── alignment/
│   ├── bin/                     # align.py, align.pl, dotplot.pl, etc.
│   └── scoring/                 # 80+ substitution matrices
├── markov/
│   ├── bin/                     # randomwalk.pl, validate.pl
│   └── lib/                     # MarkovChainFirstOrder.pm, HigherOrder.pm
├── seqtools/
│   ├── run.pl                   # Perl CLI dispatcher
│   └── SeqTools/                # Perl modules (hamming, levenshtein, etc.)
├── tests/                       # pytest test suite
├── config/
│   └── logging.yaml             # Logging configuration
├── Dockerfile                   # Python 3.12-slim + Perl
├── docker-compose.yml           # Single-service Compose file
└── pyproject.toml               # Hatchling build config
```

---

## Execution flow

### CLI path

```
User runs: biosea seqtools hamming ACGT AGGT
  │
  ▼
cli.py → seqtools_hamming_cmd(seq1, seq2)
  │
  ▼
api.seqtools_api.hamming_distance("ACGT", "AGGT")
  │
  ├─ Try: seqtools_py.hamming.distance()  ← pure Python
  │    └─ Success → return 1
  │
  └─ Fallback: perl_wrappers.distance_hamming_perl()
       └─ subprocess: perl -I seqtools -MSeqTools::hamming ...
```

### REST API path

```
Client sends: POST /distance {"seq1":"kitten","seq2":"sitting","metric":"levenshtein"}
  │
  ▼
api/server.py → validates request via Pydantic model
  │
  ▼
api.seqtools_api.levenshtein_distance("kitten", "sitting")
  │
  ├─ Try: seqtools_py.levenshtein.distance()  ← pure Python
  │    └─ Success → return 3
  │
  └─ Fallback: perl_wrappers.distance_levenshtein_perl()
       └─ subprocess: perl one-liner with SeqTools::levenshtein
```

### Alignment path (Perl-dependent)

```
User runs: biosea align seq1.fa seq2.fa --matrix BLOSUM62.mat --mode global
  │
  ▼
cli.py → align_cmd(fasta1, fasta2, matrix, mode)
  │
  ▼
api.alignment_api.align_sequences(fasta1, fasta2, matrix, mode)
  │
  ▼
perl_wrappers.alignment_wrapper.run_alignment(...)
  │
  ├─ Resolves repo root (BIOSEA_REPO_ROOT or parents[3])
  ├─ Prefers alignment/bin/align.py (Python Gotoh implementation)
  ├─ Falls back to alignment/bin/align.pl (Perl) if needed
  └─ Returns ANSI-stripped alignment output
```

---

## Python-native vs wrapper-based execution

The API layer in `seqtools_api.py` implements a **try-Python-first** strategy:

1. If the pure Python port is importable (`_PY_PORT_AVAILABLE` is `True`), call it.
2. If it raises a `ValueError`, re-raise immediately — this is a genuine validation error.
3. If it raises any other exception, fall back to the Perl wrapper.
4. If the Python port is not importable at all, go directly to Perl.

This means:

- **BWT / FM-index**: Always runs in pure Python. No Perl involved.
- **Hamming / Levenshtein / k-mer**: Runs in Python by default; Perl is the fallback.
- **Alignment**: Always delegates to the Perl wrapper layer (which internally runs `align.py`, a standalone Python DP script, or `align.pl`).
- **Markov**: Always delegates to the Perl wrapper layer (calls `randomwalk.pl`).

---

## Configuration

### Environment variables

| Variable | Purpose | Default |
|----------|---------|---------|
| `BIOSEA_REPO_ROOT` | Root directory for locating Perl scripts and scoring matrices | Three levels above `perl_wrappers/` |

### Logging

The `config/logging.yaml` file defines a `biosea` logger that writes to stdout at `INFO` level:

```yaml
loggers:
  biosea:
    level: INFO
    handlers: [console]
    propagate: no
```

### Build system

The project uses **Hatchling** as its build backend. The wheel bundles:

- `src/bio_sea_pearl/` — the Python package
- `alignment/`, `markov/`, `seqtools/` — legacy scripts, included via `force-include`
- `api/` — the FastAPI server module

See `pyproject.toml` for the full build configuration.

---

## Related pages

- [CLI Reference](cli.md) — command-level documentation
- [REST API](api.md) — endpoint schemas
- [Extending the Toolkit](developer/extending.md) — how to add new modules
