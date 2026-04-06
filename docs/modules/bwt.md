---
title: BWT & FM-index
---

# Burrows–Wheeler Transform & FM-index

The `src/bio_sea_pearl/bwt/` package provides pure Python implementations of the **suffix array**, **Burrows–Wheeler Transform (BWT)**, and **FM-index**. These data structures enable efficient exact substring search and form the foundation of modern read mappers and text compressors.

This module is entirely Python-native — no Perl dependency.

---

## Core data structures

### Suffix array

```python
from bio_sea_pearl.bwt.transform import suffix_array

sa = suffix_array("banana$")
# sa = [6, 5, 3, 1, 0, 4, 2]
# Suffixes in sorted order:
#   $, a$, ana$, anana$, banana$, na$, nana$
```

The `suffix_array(s)` function builds the suffix array using a **doubling algorithm**:

- **Time complexity:** O(n log² n)
- **Space complexity:** O(n log n)

The suffix array is a permutation of indices `[0, n)` such that `s[sa[i]:]` is the i-th suffix in lexicographic order.

### Burrows–Wheeler Transform

```python
from bio_sea_pearl.bwt.transform import bwt

transformed = bwt("banana", sentinel="$")
# transformed = "annb$aa"
```

The `bwt(s, sentinel='$')` function:

1. Appends the sentinel character to the input string.
2. Builds the suffix array.
3. Collects the character preceding each suffix (i.e. `s[sa[i] - 1]` for each position).

The BWT is reversible and tends to group identical characters together, making it useful for compression.

### FM-index

The `FMIndex` class combines the suffix array, BWT, a **C-table**, and an **Occ-matrix** to support O(m) exact pattern matching (where m is the pattern length).

```python
from bio_sea_pearl.bwt.transform import FMIndex

index = FMIndex("ACGTACGT")
l, r = index.backward_search("CGT")
# [l, r) is the half-open interval in the suffix array
# containing all suffixes that start with "CGT"
```

**Internal tables:**

| Table | Definition | Purpose |
|-------|------------|---------|
| `C[c]` | Count of characters strictly less than `c` in the BWT | Base offset for rank queries |
| `Occ[c][i]` | Number of occurrences of `c` in `BWT[0:i]` | Rank function for backward search |

**Backward search algorithm:**

The `backward_search(pattern)` method processes the pattern **in reverse**:

```
For each character c in pattern (right to left):
    l = C[c] + Occ[c][l]
    r = C[c] + Occ[c][r]
```

If `l < r` after processing all characters, the pattern exists in the text and the positions are `suffix_array[l:r]`. If `l ≥ r`, the pattern is not found.

- **Time complexity:** O(m) for a pattern of length m
- **Pattern length limit:** m ≤ 256 (enforced in the implementation)

---

## API functions

The `bwt_api` module in `src/bio_sea_pearl/api/` provides two high-level functions:

```python
from bio_sea_pearl.api import build_fm_index, search_fm_index

# Build an FM-index for a sequence
index = build_fm_index("ACGTACGT")

# Search for a pattern — returns sorted list of zero-indexed positions
positions = search_fm_index(index, "CGT")
# positions = [1, 5]
```

`build_fm_index(sequence, sentinel='$')` constructs and returns an `FMIndex` instance.

`search_fm_index(index, pattern)` performs backward search and converts the suffix array interval to a sorted list of positions.

---

## Utility scripts

### `transform.py`

Located at `src/bio_sea_pearl/bwt/transform.py`. Contains the `suffix_array()`, `bwt()`, and `FMIndex` class implementations. This is the core module.

### `find.py`

Located at `src/bio_sea_pearl/bwt/find.py`. An example script that loads a serialised `FMIndex` from a pickle file (`bwt.fmidx`) and performs a hard-coded search. Useful as a reference for programmatic usage.

---

## Usage

=== "CLI"

    ```bash
    biosea bwt search --sequence ACGTACGT --pattern CGT
    # FM-index search positions: [1, 5]
    ```

=== "REST API"

    ```bash
    curl -X POST http://localhost:8000/bwt/search \
        -H "Content-Type: application/json" \
        -d '{"sequence": "ACGTACGT", "pattern": "CGT"}'
    ```

    Response:

    ```json
    {"sequence": "ACGTACGT", "pattern": "CGT", "positions": [1, 5]}
    ```

=== "Python"

    ```python
    from bio_sea_pearl.api import build_fm_index, search_fm_index

    idx = build_fm_index("banana")
    print(search_fm_index(idx, "ana"))  # [1, 3]
    print(search_fm_index(idx, "xyz"))  # []
    ```

!!! info
    The CLI and REST API build a fresh FM-index for every request. For repeated queries against the same sequence, use the Python API to build the index once and search multiple patterns.

---

## Complexity summary

| Operation | Time | Space |
|-----------|------|-------|
| Build suffix array | O(n log² n) | O(n log n) |
| Build BWT | O(n) (given SA) | O(n) |
| Build FM-index (C + Occ) | O(n × |Σ|) | O(n × |Σ|) |
| Backward search | O(m) | O(1) |

Where n = text length, m = pattern length, |Σ| = alphabet size.

---

## Related pages

- [CLI Reference — `biosea bwt`](../cli.md#biosea-bwt)
- [REST API — `POST /bwt/search`](../api.md#post-bwtsearch)
- [Architecture](../architecture.md) — how BWT fits in the pure-Python layer
