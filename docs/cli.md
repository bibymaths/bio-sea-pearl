---
title: Command‑Line Interface
---

# Command‑Line Interface

The `biosea` CLI provides a unified entry point to the toolkit.  It is
implemented using Typer and exposes several subcommands corresponding to
the major domains of the project.

## Usage

```
biosea [OPTIONS] COMMAND [ARGS]...
```

Use `biosea --help` to view global options and a list of subcommands.  Each
subcommand provides its own help message.

### `align`

Align two sequences stored in FASTA files.  Example:

```
biosea align seq1.fa seq2.fa --matrix alignment/scoring/BLOSUM62.mat --mode global
```

Supported modes are `global`, `local` and `lcs`.  If no matrix is
specified, the underlying script uses its default scoring system.

### `markov`

Generate a random walk from a Markov model built from an input FASTA file:

```
biosea markov --fasta train.fa --length 100 --start A --order 1 --method alias
```

You can adjust the starting state, chain order, sampling method
(`alias` or `binsrch`) and pseudocount.

### `seqtools`

Compute sequence statistics:

- `biosea seqtools hamming SEQ1 SEQ2` – Hamming distance.
- `biosea seqtools levenshtein SEQ1 SEQ2` – Levenshtein distance.
- `biosea seqtools kmer SEQUENCE --k 3` – k‑mer counts returned as JSON.

### `bwt`

Search for a pattern in a sequence using an FM‑index:

```
biosea bwt search --sequence ACGTACGT --pattern CGT
```

This builds a temporary FM‑index and returns the starting positions of the
pattern.