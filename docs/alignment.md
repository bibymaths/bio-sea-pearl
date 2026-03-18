---
title: Alignment
---

# Alignment

The alignment module provides flexible dynamic‑programming algorithms for pairwise sequence alignment with separate gap
opening and extension penalties. It supports global (Needleman–Wunsch), local (Smith–Waterman) and longest common
subsequence (LCS) modes, and can utilise a variety of scoring matrices.

## Features

* **Affine gap penalties** – both the Python and Perl implementations use the Gotoh algorithm. Gap opening (`-g`) and
  gap extension (`-e`) penalties can be adjusted on the command line. The default values are 10 for opening and 1 for
  extension.
* **Multiple modes** – global alignment scores the entire sequences, local alignment finds the best matching
  subsequences and LCS computes the longest common subsequence length.
* **Scoring matrices** – the `alignment/scoring` directory contains numerous substitution matrices (BLOSUM, PAM, VTML
  etc.), enabling meaningful protein alignments. The Python script `read_matrix_file` parses these files into nested
  dictionaries for fast lookup.
* **Parallelisation** – `align.py` uses Python’s `multiprocessing` module to split large alignments into smaller tasks.
  A heuristic based on sequence lengths determines when to fork child processes for improved performance.
* **Perl utilities** – the Perl script `align.pl` implements similar alignment logic and is used by `dotplot.pl` and
  `traceplot.pl` to visualise alignment scores.  `dotplot.pl` produces dot plots of similarity scores between two
  sequences, while `traceplot.pl` draws the dynamic‑programming matrix and traceback arrows.

## Usage examples

Align two sequences using Python with the BLOSUM62 matrix:

```sh
poetry run python alignment/bin/align.py \
    --matrix alignment/scoring/BLOSUM62.mat \
    --mode global \
    seq1.fa seq2.fa
```

Generate a dot plot using the Perl scripts:

```sh
# Install Perl if necessary (on Ubuntu: sudo apt-get install -y perl)
perl alignment/bin/dotplot.pl seq1.fa seq2.fa > dotplot.svg
```

## Transformations and plotting

The auxiliary script `alignfm.py` writes DP tables to files for downstream visualisation.  `transform.py` performs
common conversions and can be adapted to pre‑process input sequences. For example, to view the scoring matrix as a
heatmap you can parse the `.mat` files with Python and plot them via `matplotlib`.
