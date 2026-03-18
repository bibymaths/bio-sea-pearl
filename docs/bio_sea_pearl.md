---
title: BWT and FM‑index
---

# Burrows–Wheeler Transform and FM‑index

The `src/bio_sea_pearl/bwt` package contains Python implementations of the Burrows–Wheeler Transform (BWT) and the
FM‑index. These data structures underlie modern read mappers and compressors because they enable efficient substring
queries while keeping memory use low.

## Components

* **`suffix_array(s)`** – builds the suffix array for a string `s` using a doubling algorithm. The suffix array is an
  array of indices into the string such that the suffixes appear in lexicographic order.
* **`bwt(s, sentinel='$')`** – computes the BWT of a string by sorting rotations via the suffix array and collecting the
  last characters. A sentinel character not present in the alphabet is appended automatically.
* **`FMIndex` class** – encapsulates the FM‑index. Upon construction it builds the suffix array, BWT, the `C` table (
  cumulative counts of characters) and the `Occ` matrix (rank queries for characters). The method
  `backward_search(pattern)` returns the half‑open interval `[l,r)` in the suffix array where suffixes start with
  `pattern`.
* **`transform.py`** – command‑line script that reads a FASTA file, constructs an `FMIndex` for each sequence and
  serialises it to a `.fmidx` file using `pickle`.
* **`find.py`** – example script that loads a serialised FM‑index (`bwt.fmidx`) and searches for a hard‑coded pattern
  using backward search.

## Example

To build FM‑indices for sequences in a FASTA file and perform a search:

```sh
# Build the index (outputs one .fmidx per sequence)
python src/bio_sea_pearl/bwt/transform.py sequences.fa

# Load the index and search for a pattern
python src/bio_sea_pearl/bwt/find.py
```

The FM‑index can be wrapped into a high‑level Python API or served via a web service to enable fast substring queries on
large genomes. Integrating this component with the alignment and Markov modules is part of the long‑term roadmap.
