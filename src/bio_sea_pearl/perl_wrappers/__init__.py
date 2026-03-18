"""Perl wrapper package.

This package exposes Python functions that delegate to existing Perl scripts or
modules in the repository.  The wrappers use `subprocess` to execute the
underlying Perl programs and return results as native Python objects.

They are designed to insulate the rest of the codebase from the details of
invoking Perl.  Any errors in the external scripts are propagated as
exceptions.
"""

from .alignment_wrapper import run_alignment
from .markov_wrapper import run_markov_walk
from .seqtools_wrapper import (distance_hamming_perl,
                               distance_levenshtein_perl,
                               kmer_counts_perl)

__all__ = [
    "run_alignment",
    "run_markov_walk",
    "distance_hamming_perl",
    "distance_levenshtein_perl",
    "kmer_counts_perl",
]
