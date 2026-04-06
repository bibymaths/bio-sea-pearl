---
title: Alignment
---

# Alignment Module

The alignment subsystem implements pairwise sequence alignment using the Gotoh algorithm — a dynamic programming approach that supports affine gap penalties. It handles global (Needleman–Wunsch), local (Smith–Waterman), and longest common subsequence (LCS) alignment modes.

---

## Components

### `alignment/bin/align.py`

The primary alignment implementation (813 lines). This standalone Python script:

- Reads two FASTA files as input.
- Parses a substitution matrix from a `.mat` file.
- Constructs three DP matrices (match, insert, delete) for the Gotoh algorithm.
- Supports **multiprocessing**: a heuristic based on sequence lengths determines when to fork child processes to split the computation.
- Outputs the aligned sequences, score, and traceback to stdout.

**Command-line interface:**

```bash
python alignment/bin/align.py \
    --matrix alignment/scoring/BLOSUM62.mat \
    --mode global \
    seq1.fa seq2.fa
```

| Option     | Description |
|------------|-------------|
| `--matrix` | Path to a scoring matrix file |
| `--mode`   | `global`, `local`, or `lcs` |
| `-g`       | Gap opening penalty (default: 10) |
| `-e`       | Gap extension penalty (default: 1) |

### `alignment/bin/align.pl`

Perl implementation of the same algorithm. Used as a fallback if `align.py` is not available. The wrapper layer (`alignment_wrapper.py`) prefers the Python version.

### Auxiliary scripts

| Script | Language | Purpose |
|--------|----------|---------|
| `alignfm.py` | Python | Writes DP tables to files for downstream visualisation |
| `transform.py` | Python | Pre-processes and transforms input sequences |
| `dotplot.pl` | Perl | Generates dot-plot visualisations of pairwise similarity |
| `traceplot.pl` | Perl | Draws the DP matrix with traceback arrows |

---

## Scoring matrices

The `alignment/scoring/` directory contains **80+ substitution matrices** for protein and nucleotide alignment:

| Family | Matrices | Typical use |
|--------|----------|-------------|
| **BLOSUM** | BLOSUM30–BLOSUM100, BLOSUMN | Protein alignment (block substitution) |
| **PAM** | PAM10–PAM500 | Protein alignment (point accepted mutation) |
| **VTML** | VTML10, VTML20, VTML40, VTML80, VTML120, VTML160 | Protein alignment (variable time) |
| **BENNER** | BENNER6, BENNER22, BENNER74 | Protein alignment |
| **GONNET** | GONNET | General protein alignment |
| **DAYHOFF** | DAYHOFF | Classic protein matrix |
| **NUC.4.4** | NUC.4.4 | DNA/RNA nucleotide alignment |

Matrix files use a tabular format. The `read_matrix_file` function in `align.py` parses them into nested dictionaries for O(1) lookup during DP.

!!! tip
    For protein sequences, BLOSUM62 is the most commonly used matrix. For DNA, use `NUC.4.4` or omit the matrix to use the script's built-in default scoring. All matrices are located under `alignment/scoring/`.

---

## Affine gap model

The Gotoh algorithm extends basic Needleman–Wunsch by separating gap opening and gap extension penalties:

```
Gap cost = opening_penalty + (gap_length - 1) × extension_penalty
```

With default penalties of `g=10` (opening) and `e=1` (extension):

- A gap of length 1 costs 10
- A gap of length 5 costs 14
- A gap of length 10 costs 19

This model penalises opening new gaps more heavily than extending existing ones, which better reflects biological insertion/deletion events.

---

## Alignment modes

| Mode | Algorithm | Behaviour |
|------|-----------|-----------|
| `global` | Needleman–Wunsch with Gotoh extension | Aligns the full length of both sequences end-to-end |
| `local` | Smith–Waterman with Gotoh extension | Finds the highest-scoring local subsequence alignment |
| `lcs` | Longest common subsequence | Computes the longest common subsequence (match-only, no mismatches) |

---

## Wrapper integration

The Python API layer accesses alignment through:

```
bio_sea_pearl.api.align_sequences(fasta1, fasta2, matrix=None, mode="global")
  └─ perl_wrappers.alignment_wrapper.run_alignment(...)
       ├─ Locates repo root (BIOSEA_REPO_ROOT or path traversal)
       ├─ Checks for alignment/bin/align.py → preferred
       ├─ Falls back to alignment/bin/align.pl → Perl
       ├─ Resolves matrix path (defaults to blosum62.mat if unspecified)
       └─ Returns ANSI-stripped stdout
```

The wrapper validates that both FASTA files exist and that the matrix file (if specified) is present before spawning the subprocess.

---

## Usage via CLI and API

=== "CLI"

    ```bash
    biosea align seq1.fa seq2.fa --matrix alignment/scoring/BLOSUM62.mat --mode global
    ```

=== "REST API"

    ```bash
    curl -X POST http://localhost:8000/align \
        -H "Content-Type: application/json" \
        -d '{
            "fasta1": ">s1\nMKWVTFISL",
            "fasta2": ">s2\nMKWVTFISM",
            "matrix": "BLOSUM62.mat",
            "mode": "global"
        }'
    ```

=== "Python"

    ```python
    from bio_sea_pearl.api import align_sequences
    result = align_sequences("seq1.fa", "seq2.fa", matrix="BLOSUM62.mat", mode="global")
    print(result)
    ```

!!! warning
    Alignment requires Perl to be installed, even though the preferred backend is `align.py` (Python). The wrapper layer relies on the Perl ecosystem being available as a fallback.

---

## Related pages

- [CLI Reference — `biosea align`](../cli.md#biosea-align)
- [REST API — `POST /align`](../api.md#post-align)
- [Architecture — Alignment execution path](../architecture.md#alignment-path-perl-dependent)
