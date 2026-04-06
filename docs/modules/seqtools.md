---
title: Sequence Tools
---

# Sequence Tools

Sequence tools span two implementations: the original **Perl modules** in `seqtools/` and the **pure Python ports** in `src/bio_sea_pearl/seqtools_py/`. The API layer transparently selects the Python implementation when available, falling back to Perl if needed.

---

## Capabilities

| Function | Python port | Perl module | Description |
|----------|:-----------:|:-----------:|-------------|
| Hamming distance | :material-check: | :material-check: | Positional mismatch count for equal-length strings |
| Levenshtein distance | :material-check: | :material-check: | Edit distance (insertions, deletions, substitutions) |
| K-mer counting | :material-check: | :material-check: | Frequency table of all k-length substrings |
| Boyer–Moore search | — | :material-check: | Fast pattern matching with skip heuristics |
| FASTA parsing | — | :material-check: | Read FASTA files into sequence records |

---

## Python implementations

Located in `src/bio_sea_pearl/seqtools_py/`. These are the preferred runtime path.

### `hamming.py`

```python
def distance(a: str, b: str) -> int
```

Single-pass comparison. O(n) time, O(1) space.
Raises `ValueError` if the two strings have different lengths.

### `levenshtein.py`

```python
def distance(s: str, t: str) -> int
```

Classic dynamic programming with **space optimisation** — uses two rows instead of the full matrix. O(n×m) time, O(min(n,m)) space. Handles empty strings correctly.

### `kmer.py`

```python
def counts(seq: str, k: int) -> Dict[str, int]
```

Sliding-window hash-based counting. O(n) time.
Raises `ValueError` if `k ≤ 0` or `k > len(seq)`.

### Import paths

```python
# Via the public API (recommended)
from bio_sea_pearl.api import hamming_distance, levenshtein_distance, kmer_counts

# Direct access to the Python implementations
from bio_sea_pearl.seqtools_py import hamming_distance, levenshtein_distance, kmer_counts
```

---

## Perl implementations

Located in `seqtools/SeqTools/`. These are the original implementations and serve as fallbacks.

### `hamming.pm`

Computes Hamming distance using **Inline C** for performance. The C kernel compares characters in a tight loop, avoiding Perl interpreter overhead. Falls back to pure Perl if Inline C is not available.

### `levenshtein.pm`

Edit distance via **Inline C** with a space-optimised DP approach — O(min(n,m)) working space. Similar fallback behaviour as `hamming.pm`.

### `kmer.pm`

K-mer counting using Perl hash aggregation. Outputs tab-separated `kmer\tcount` pairs.

### `boyermoore.pm`

Boyer–Moore string matching algorithm. Uses the bad-character heuristic to skip ahead on mismatches, achieving sublinear average-case performance. This module is available only through Perl; there is no Python port.

### `fasta.pm`

Minimal FASTA parser. Reads header lines (starting with `>`) and concatenates sequence lines into single records.

### Perl CLI dispatcher

`seqtools/run.pl` is a unified entry point for Perl-based operations:

```bash
perl seqtools/run.pl kmer input.fa 3
perl seqtools/run.pl hamming seq1.fa seq2.fa
perl seqtools/run.pl levenshtein seq1.fa seq2.fa
perl seqtools/run.pl motif input.fa PATTERN
```

---

## Fallback strategy

The API layer in `seqtools_api.py` implements a **try-Python-first** pattern:

```
hamming_distance("ACGT", "AGGT")
  │
  ├─ _PY_PORT_AVAILABLE is True?
  │    ├─ Yes → call seqtools_py.hamming.distance()
  │    │         ├─ ValueError → re-raise (validation error)
  │    │         └─ Other exception → fall back to Perl
  │    │
  │    └─ No → call perl_wrappers.distance_hamming_perl()
  │
  └─ Return integer result
```

This means:

- **No Perl installed?** Hamming, Levenshtein, and k-mer functions still work via the Python ports.
- **Python port broken?** The system degrades gracefully to Perl.
- **Genuine input errors** (e.g. unequal-length strings for Hamming) are raised immediately regardless of which backend is used.

!!! tip
    If you are deploying without Perl (e.g. a minimal Python container), the seqtools commands will work correctly as long as the `seqtools_py` module is importable. Only alignment and Markov commands strictly require Perl.

---

## Usage

=== "CLI"

    ```bash
    biosea seqtools hamming ACGTACGT ACGTACGA
    # Hamming distance: 1

    biosea seqtools levenshtein kitten sitting
    # Levenshtein distance: 3

    biosea seqtools kmer ACGTACGT --k 3
    # {"ACG": 2, "CGT": 2, "GTA": 1, "TAC": 1}
    ```

=== "REST API"

    ```bash
    # Distance
    curl -X POST http://localhost:8000/distance \
        -H "Content-Type: application/json" \
        -d '{"seq1": "ACGT", "seq2": "AGGT", "metric": "hamming"}'

    # K-mer
    curl -X POST http://localhost:8000/kmer \
        -H "Content-Type: application/json" \
        -d '{"sequence": "ACGTACGT", "k": 3}'
    ```

=== "Python"

    ```python
    from bio_sea_pearl.api import hamming_distance, levenshtein_distance, kmer_counts

    hamming_distance("ACGT", "AGGT")       # 1
    levenshtein_distance("kitten", "sitting")  # 3
    kmer_counts("ACGTACGT", 3)             # {"ACG": 2, "CGT": 2, "GTA": 1, "TAC": 1}
    ```

---

## Related pages

- [CLI Reference — `biosea seqtools`](../cli.md#biosea-seqtools)
- [REST API — `POST /distance` and `POST /kmer`](../api.md#post-distance)
- [Architecture — Python-native vs wrapper-based](../architecture.md#python-native-vs-wrapper-based-execution)
