---
title: Markov Chains
---

# Markov Chains

The `markov` directory contains a set of Perl modules and scripts implementing first‑ and higher‑order Markov chain
simulators. These tools make it easy to learn transition probabilities from sequence data and generate random walks or
evaluate the fit of a model to observed data.

## Components

### `MarkovChainFirstOrder.pm`

This module exposes an object‑oriented interface for constructing and sampling from first‑order Markov chains. Key
features include:

* **Transition table construction** – the class method `build_transitions` reads one or more sequences (via
  `read_fasta`) and counts all adjacent pairs to obtain transitions.
* **Pseudocount support** – when creating a new chain, a pseudocount can be supplied to avoid zero‑probability events.
* **Sampling methods** – two sampling strategies are implemented:
    * *Alias method* – constant‑time sampling using a precomputed alias table.
    * *Binary search* – logarithmic‑time sampling by binary searching cumulative probabilities.
* **Random walks** – `generate(N, start)` returns a string of length `N` sampled from the chain starting at state
  `start`.

### `MarkovChainHigherOrder.pm`

Higher‑order chains treat k‑mers as states. This module extends the first‑order implementation to support an arbitrary
order by grouping adjacent symbols. It shares the same sampling and pseudocount features.

### Command‑line utilities

* `randomwalk.pl` – builds a Markov chain from an input FASTA file and emits random sequences of a specified length.
  Options include chain order, pseudocount and sampling method.
* `validate.pl` – compares a hidden Markov model against observed sequences, computing log‑likelihoods and summarising
  model fit. This can be used to tune the order and pseudocount parameters.

## Example usage

Generate a random DNA walk of length 100 starting from state `A` using the alias method:

```sh
perl markov/bin/randomwalk.pl \
    --fasta tests/data/sample.fasta \
    --length 100 \
    --start A \
    --method alias
```

The Markov modules can be accessed from Python via the `subprocess` module or ported to native Python for integration
into a unified CLI (see the roadmap for details).
