---
title: Overview
---

# Bio Sea Pearl Toolkit Overview

This toolkit brings together a collection of Python and Perl utilities for sequence analysis and stochastic modelling.  It grew from practical needs in bioinformatics: aligning DNA or protein sequences, exploring Markov chain behaviour and building efficient text indexes.  As a result it includes both modern Python code and a mature body of Perl modules.

The project is organised into a few high‑level components:

* **Alignment** – dynamic‑programming algorithms for global, local and longest common subsequence (LCS) alignments.  Python scripts (`align.py`, `alignfm.py`, `transform.py`) implement affine gap penalties, parallelisation and various output modes.  Perl scripts (`align.pl`, `dotplot.pl`, `traceplot.pl`) provide plotting and compatibility with earlier versions.  Scoring matrices in `alignment/scoring/` include BLOSUM, PAM and VTML families for protein alignment.
* **Markov** – high‑precision Markov chain simulators implemented in Perl.  `MarkovChainFirstOrder.pm` and `MarkovChainHigherOrder.pm` build transition tables from FASTA sequences, support pseudocount smoothing and provide fast sampling via the alias method or binary search.  Utilities like `randomwalk.pl` and `validate.pl` generate random walks or evaluate hidden Markov models.
* **SeqTools** – a collection of Perl modules for common sequence operations.  Modules implement pattern matching (Boyer–Moore search), FASTA parsing, Hamming and Levenshtein distances, k‑mer generation and more.
* **Burrows–Wheeler Transform (BWT) and FM‑index** – pure Python implementations of suffix arrays, the BWT and the FM‑index are provided in `src/bio_sea_pearl/bwt`.  These enable fast substring queries and serve as a foundation for compressed text indices.

The long‑term goal is to unify these capabilities under a single Python interface while preserving the performance and features of the existing Perl code.  See the updated roadmap in `TODO.md` for a detailed plan.
