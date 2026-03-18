"""Public API layer for the Bio Sea Pearl toolkit.

This package exposes pure Python functions that wrap the underlying Perl
implementations or use native Python code.  The intent is to provide a
uniform interface for other modules (including the CLI and REST API) to
call without worrying about language interop.  All functions here should
only depend on Python types, should include type hints, and should not
perform any I/O beyond calling subprocesses via the wrappers.

Modules exported:

* `alignment_api` – sequence alignment operations
* `markov_api` – Markov chain generation
* `seqtools_api` – distance and k‑mer utilities
* `bwt_api` – Burrows–Wheeler Transform and FM‑index utilities
"""

from .alignment_api import align_sequences
from .markov_api import generate_walk
from .seqtools_api import (
    hamming_distance,
    levenshtein_distance,
    kmer_counts,
)
from .bwt_api import (
    build_fm_index,
    search_fm_index,
)

__all__ = [
    "align_sequences",
    "generate_walk",
    "hamming_distance",
    "levenshtein_distance",
    "kmer_counts",
    "build_fm_index",
    "search_fm_index",
]
