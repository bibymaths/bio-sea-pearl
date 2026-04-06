---
title: Quickstart
---

# Quickstart

This page walks through representative tasks using the `biosea` CLI. Each example maps directly to a subsystem in the toolkit.

!!! note "Prerequisites"
    The examples below assume you have installed the package (see [Installation](installation.md)). Alignment and Markov commands require Perl. Sequence tools and BWT commands work with Python alone.

---

## Compute sequence distances

Hamming distance compares two equal-length strings position by position:

```bash
biosea seqtools hamming ACGTACGT ACGTACGA
```

```
Hamming distance: 1
```

Levenshtein distance handles strings of different lengths by counting insertions, deletions, and substitutions:

```bash
biosea seqtools levenshtein kitten sitting
```

```
Levenshtein distance: 3
```

---

## Count k-mers

Extract k-mer frequencies from a sequence:

```bash
biosea seqtools kmer ACGTACGT --k 3
```

```json
{"ACG": 2, "CGT": 2, "GTA": 1, "TAC": 1}
```

---

## Search with BWT / FM-index

Build a transient FM-index and search for a pattern:

```bash
biosea bwt search --sequence ACGTACGT --pattern CGT
```

```
FM-index search positions: [1, 5]
```

Positions are zero-indexed into the input sequence.

---

## Align two sequences

Prepare two FASTA files, then run a global alignment with a scoring matrix:

```bash
biosea align seq1.fa seq2.fa \
    --matrix alignment/scoring/BLOSUM62.mat \
    --mode global
```

Supported alignment modes:

| Mode     | Algorithm            | Description                              |
|----------|----------------------|------------------------------------------|
| `global` | Needleman–Wunsch     | End-to-end alignment of both sequences   |
| `local`  | Smith–Waterman       | Best-scoring local subsequence alignment |
| `lcs`    | Longest common subseq| Longest common subsequence length        |

!!! tip
    If you omit `--matrix`, the alignment script uses its built-in default scoring. For protein sequences, specifying a substitution matrix (BLOSUM, PAM, VTML) is strongly recommended.

---

## Generate a Markov random walk

Train a Markov chain on a FASTA file and sample a sequence:

```bash
biosea markov \
    --fasta training.fa \
    --length 100 \
    --start A \
    --order 1 \
    --method alias
```

This builds a first-order Markov chain from the sequences in `training.fa`, then generates a 100-character random walk starting from state `A` using the alias sampling method.

---

## Use the REST API

Start the server:

```bash
uvicorn api.server:app --host 0.0.0.0 --port 8000
```

Then call any endpoint. For example, compute a Levenshtein distance:

```bash
curl -s -X POST http://127.0.0.1:8000/distance \
    -H "Content-Type: application/json" \
    -d '{"seq1": "kitten", "seq2": "sitting", "metric": "levenshtein"}' | python -m json.tool
```

```json
{
    "metric": "levenshtein",
    "seq1": "kitten",
    "seq2": "sitting",
    "distance": 3
}
```

See the [REST API reference](../api.md) for all endpoints and request schemas.

---

## Where to go next

| Goal                          | Page                                            |
|-------------------------------|-------------------------------------------------|
| Full CLI command reference    | [CLI Reference](../cli.md)                      |
| REST API endpoints            | [REST API](../api.md)                           |
| Understand the architecture   | [Architecture](../architecture.md)              |
| Deep-dive into alignment      | [Alignment Module](../modules/alignment.md)     |
| Explore Markov chains         | [Markov Module](../modules/markov.md)           |
| BWT and FM-index internals    | [BWT & FM-index](../modules/bwt.md)             |
| Extend the toolkit            | [Developer Guide](../developer/extending.md)    |
