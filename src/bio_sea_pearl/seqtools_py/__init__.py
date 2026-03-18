"""Pure Python re‑implementations of select SeqTools algorithms.

These modules provide native Python implementations of the Hamming and
Levenshtein distance algorithms and k‑mer counting.  They are intended
to mirror the behaviour of the existing Perl modules while eliminating
cross-language overhead.  When performance is critical, consider
leveraging NumPy for further acceleration.
"""

from .hamming import distance as hamming_distance
from .levenshtein import distance as levenshtein_distance
from .kmer import counts as kmer_counts

__all__ = [
    "hamming_distance",
    "levenshtein_distance",
    "kmer_counts",
]