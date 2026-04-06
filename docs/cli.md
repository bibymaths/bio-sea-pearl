---
title: CLI Reference
---

# CLI Reference

The `biosea` command is the unified entry point to the Bio Sea Pearl toolkit. It is built with [Typer](https://typer.tiangolo.com/) and delegates all work to the API layer in `src/bio_sea_pearl/api/`.

```
biosea [OPTIONS] COMMAND [ARGS]...
```

---

## Global options

| Option   | Description                |
|----------|----------------------------|
| `--help` | Show help and exit         |

---

## `biosea align`

Align two sequences from FASTA files using the Gotoh dynamic programming algorithm.

```bash
biosea align FASTA1 FASTA2 [OPTIONS]
```

**Arguments:**

| Argument | Type | Description |
|----------|------|-------------|
| `FASTA1` | Path | Path to the first FASTA file |
| `FASTA2` | Path | Path to the second FASTA file |

**Options:**

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--matrix` | Path | *(built-in default)* | Path to a substitution matrix file (e.g. `alignment/scoring/BLOSUM62.mat`) |
| `--mode` | String | `global` | Alignment mode: `global`, `local`, or `lcs` |

**Examples:**

```bash
# Global alignment with BLOSUM62
biosea align seq1.fa seq2.fa --matrix alignment/scoring/BLOSUM62.mat --mode global

# Local alignment (Smith–Waterman)
biosea align protein_a.fa protein_b.fa --matrix alignment/scoring/PAM250.mat --mode local

# Longest common subsequence
biosea align seq1.fa seq2.fa --mode lcs
```

!!! warning
    This command requires Perl to be installed. The underlying `align.py` script (or `align.pl` fallback) is invoked via subprocess. If neither is found, the command will fail with an error.

**Execution path:**
`cli.py` → `align_sequences()` → `alignment_wrapper.run_alignment()` → `alignment/bin/align.py`

---

## `biosea markov`

Generate a random walk from a Markov chain trained on sequences in a FASTA file.

```bash
biosea markov [OPTIONS]
```

**Options:**

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--fasta` | Path | *(required)* | Path to the training FASTA file |
| `--length` | Integer | *(required)* | Number of characters to generate |
| `--start` | String | `A` | Starting state for the walk |
| `--order` | Integer | `1` | Markov chain order (1 = first-order, higher = k-mer states) |
| `--method` | String | `alias` | Sampling method: `alias` or `binsrch` |
| `--pseudocount` | Float | `0` | Smoothing pseudocount added to all transitions |

**Example:**

```bash
biosea markov --fasta training.fa --length 200 --start A --order 2 --method alias --pseudocount 0.01
```

!!! note
    The output string length is `len(start) + length`. The starting state is included as a prefix of the generated sequence.

!!! warning
    This command requires Perl. It invokes `markov/bin/randomwalk.pl` via subprocess.

**Execution path:**
`cli.py` → `generate_walk()` → `markov_wrapper.run_markov_walk()` → `markov/bin/randomwalk.pl`

---

## `biosea seqtools`

Sequence distance and k-mer utilities. This command group contains three subcommands.

### `biosea seqtools hamming`

Compute the Hamming distance between two equal-length sequences.

```bash
biosea seqtools hamming SEQ1 SEQ2
```

**Arguments:**

| Argument | Type | Description |
|----------|------|-------------|
| `SEQ1` | String | First sequence |
| `SEQ2` | String | Second sequence (must be same length as SEQ1) |

**Example:**

```bash
biosea seqtools hamming ACGTACGT ACGTACGA
# Hamming distance: 1
```

### `biosea seqtools levenshtein`

Compute the Levenshtein (edit) distance between two sequences. Handles sequences of different lengths.

```bash
biosea seqtools levenshtein SEQ1 SEQ2
```

**Example:**

```bash
biosea seqtools levenshtein kitten sitting
# Levenshtein distance: 3
```

### `biosea seqtools kmer`

Count k-mer frequencies in a sequence. Output is JSON.

```bash
biosea seqtools kmer SEQUENCE --k K
```

**Options:**

| Option | Type | Description |
|--------|------|-------------|
| `--k` | Integer | k-mer length (required) |

**Example:**

```bash
biosea seqtools kmer ACGTACGT --k 3
# {"ACG": 2, "CGT": 2, "GTA": 1, "TAC": 1}
```

!!! tip
    The seqtools commands use pure Python implementations by default and do **not** require Perl. They fall back to the Perl SeqTools modules only if the Python implementation fails.

**Execution path:**
`cli.py` → `hamming_distance()` / `levenshtein_distance()` / `kmer_counts()` → `seqtools_py` (Python) or `seqtools_wrapper` (Perl fallback)

---

## `biosea bwt`

Burrows–Wheeler Transform and FM-index tools.

### `biosea bwt search`

Build a transient FM-index for a sequence and search for a pattern.

```bash
biosea bwt search --sequence SEQUENCE --pattern PATTERN
```

**Options:**

| Option | Type | Description |
|--------|------|-------------|
| `--sequence` | String | The text to index |
| `--pattern` | String | The pattern to search for |

**Example:**

```bash
biosea bwt search --sequence ACGTACGT --pattern CGT
# FM-index search positions: [1, 5]
```

Positions are **zero-indexed** into the input sequence. If the pattern is not found, an empty list is returned.

!!! info
    The FM-index is built in-memory for each invocation. It is not persisted to disk. For repeated searches against the same sequence, use the Python API directly — see [BWT & FM-index](modules/bwt.md).

**Execution path:**
`cli.py` → `build_fm_index()` + `search_fm_index()` → `bwt.transform.FMIndex`

---

## Exit codes

| Code | Meaning |
|------|---------|
| `0`  | Success |
| `1`  | Runtime error (missing file, invalid input, Perl not found, etc.) |

---

## Related pages

- [Quickstart](getting-started/quickstart.md) — worked examples for each command
- [REST API](api.md) — equivalent functionality over HTTP
- [Architecture](architecture.md) — how the CLI dispatches to the API layer