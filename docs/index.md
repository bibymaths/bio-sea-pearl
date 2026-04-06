---
title: Overview
---

<div style="text-align: center;">

<img src="../assets/logo.svg" width="300">

<p>
<a href="https://doi.org/10.5281/zenodo.19435099"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.19435099.svg"></a>
<a href="https://pypi.org/project/bio-sea-pearl/"><img src="https://img.shields.io/pypi/v/bio-sea-pearl.svg"></a>
<a href="https://pypi.org/project/bio-sea-pearl/"><img src="https://img.shields.io/pypi/pyversions/bio-sea-pearl.svg"></a>
<a href="https://pypi.org/project/bio-sea-pearl/"><img src="https://img.shields.io/pypi/dm/bio-sea-pearl.svg"></a>
<a href="https://github.com/bibymaths/bio-sea-pearl/releases"><img src="https://img.shields.io/github/v/release/bibymaths/bio-sea-pearl"></a>
<a href="https://bibymaths.github.io/bio-sea-pearl/"><img src="https://img.shields.io/badge/docs-GitHub%20Pages-blue"></a>
<a href="https://github.com/bibymaths/bio-sea-pearl/blob/main/LICENSE"><img src="https://img.shields.io/github/license/bibymaths/bio-sea-pearl"></a>
<a href="https://github.com/bibymaths/bio-sea-pearl/actions/workflows/release.yml"><img src="https://img.shields.io/github/actions/workflow/status/bibymaths/bio-sea-pearl/release.yml?branch=main&label=release"></a>
<a href="https://github.com/bibymaths/bio-sea-pearl/actions/workflows/docs.yml"><img src="https://img.shields.io/github/actions/workflow/status/bibymaths/bio-sea-pearl/docs.yml?branch=main&label=docs"></a>
<a href="https://github.com/bibymaths/bio-sea-pearl/pkgs/container/bio-sea-pearl"><img src="https://img.shields.io/badge/container-ghcr.io-blue"></a>
</p>

</div>

This toolkit brings together a collection of Python and Perl utilities for sequence analysis and stochastic modelling.
It grew from practical needs in bioinformatics: aligning DNA or protein sequences, exploring Markov chain behaviour and
building efficient text indexes. As a result it includes both modern Python code and a mature body of Perl modules.

The project is organised into a few high‑level components:

* **Alignment** – dynamic‑programming algorithms for global, local and longest common subsequence (LCS) alignments.
  Python scripts (`align.py`, `alignfm.py`, `transform.py`) implement affine gap penalties, parallelisation and various
  output modes. Perl scripts (`align.pl`, `dotplot.pl`, `traceplot.pl`) provide plotting and compatibility with earlier
  versions. Scoring matrices in `alignment/scoring/` include BLOSUM, PAM and VTML families for protein alignment.
* **Markov** – high‑precision Markov chain simulators implemented in Perl.  `MarkovChainFirstOrder.pm` and
  `MarkovChainHigherOrder.pm` build transition tables from FASTA sequences, support pseudocount smoothing and provide
  fast sampling via the alias method or binary search. Utilities like `randomwalk.pl` and `validate.pl` generate random
  walks or evaluate hidden Markov models.
* **SeqTools** – a collection of Perl modules for common sequence operations. Modules implement pattern matching (
  Boyer–Moore search), FASTA parsing, Hamming and Levenshtein distances, k‑mer generation and more.
* **Burrows–Wheeler Transform (BWT) and FM‑index** – pure Python implementations of suffix arrays, the BWT and the
  FM‑index are provided in `src/bio_sea_pearl/bwt`. These enable fast substring queries and serve as a foundation for
  compressed text indices.
