---
title: Markov Chains
---

# Markov Chain Module

The `markov/` directory contains Perl modules and scripts for training Markov chain models on sequence data and generating random walks. The implementation supports both first-order and arbitrary higher-order chains, with two sampling strategies and optional pseudocount smoothing.

---

## Components

### Perl modules

#### `MarkovChainFirstOrder.pm`

Object-oriented first-order Markov chain implementation (`markov/lib/`).

**Key methods:**

| Method | Description |
|--------|-------------|
| `new(pseudocount)` | Constructor; accepts an optional pseudocount for smoothing |
| `build_transitions(sequences)` | Reads sequences (via `read_fasta`) and counts all adjacent symbol pairs |
| `generate(N, start)` | Returns a string of length `len(start) + N` starting from `start` |

**Sampling methods:**

- **Alias method** — constant-time O(1) sampling using a precomputed alias table. Preferred for large generation tasks.
- **Binary search** — O(log k) sampling by binary searching cumulative probability distributions. Simpler to implement but slower per sample.

**Pseudocount smoothing:**

When a pseudocount is provided, it is added to every possible transition before normalising. This prevents zero-probability events and ensures the chain can always make a transition, even for states unseen in the training data.

#### `MarkovChainHigherOrder.pm`

Extends the first-order model to support **arbitrary-order** chains. Instead of single symbols, higher-order chains use **k-mers as states**. A chain of order `k` groups `k` consecutive symbols to form each state and predicts the next symbol based on that context.

The interface mirrors `MarkovChainFirstOrder.pm`:

- Same constructor, `build_transitions`, and `generate` signatures.
- Same sampling methods (alias, binary search).
- Same pseudocount support.

!!! note
    The `generate(N, start)` method returns a string of length `len(start) + N`. The start string is included as a prefix of the output. For a higher-order chain, the start string should be at least as long as the chain order.

---

### Command-line scripts

#### `randomwalk.pl`

The primary entry point for Markov chain generation (`markov/bin/`).

```bash
perl markov/bin/randomwalk.pl \
    --fasta INPUT.fa \
    --length 100 \
    --start A \
    --order 1 \
    --method alias \
    --pseudocount 0.01
```

| Option          | Type    | Default  | Description |
|-----------------|---------|----------|-------------|
| `--fasta`       | Path    | required | Training FASTA file |
| `--length`      | Integer | required | Characters to generate |
| `--start`       | String  | `A`      | Starting state |
| `--order`       | Integer | `1`      | Chain order |
| `--method`      | String  | `alias`  | Sampling method: `alias` or `binsrch` |
| `--pseudocount` | Float   | `0`      | Smoothing pseudocount |

Output: a single line containing the generated DNA sequence.

#### `validate.pl`

Compares a Markov model against observed sequences by computing log-likelihoods. Useful for evaluating model fit and tuning the order and pseudocount parameters.

---

## Wrapper integration

The Markov module is accessed through the Python wrapper layer:

```
bio_sea_pearl.api.generate_walk(fasta, length, start, order, method, pseudocount)
  └─ perl_wrappers.markov_wrapper.run_markov_walk(...)
       ├─ Locates repo root
       ├─ Constructs: perl markov/bin/randomwalk.pl --fasta ... --length ...
       ├─ Sets PERL5LIB to include markov/lib/
       ├─ Filters stdout for valid DNA sequences
       └─ Returns the generated sequence string
```

The wrapper filters the output to extract only valid DNA character lines, discarding any diagnostic messages or warnings from the Perl script.

---

## Usage

=== "CLI"

    ```bash
    biosea markov --fasta training.fa --length 200 --start A --order 2 --method alias
    ```

=== "REST API"

    ```bash
    curl -X POST http://localhost:8000/markov \
        -H "Content-Type: application/json" \
        -d '{
            "fasta": ">train\nACGTACGTACGT",
            "length": 50,
            "start": "AC",
            "order": 2,
            "method": "alias",
            "pseudocount": 0.01
        }'
    ```

=== "Python"

    ```python
    from bio_sea_pearl.api import generate_walk
    walk = generate_walk("training.fa", length=100, start="A", order=1)
    print(walk)
    ```

---

## Choosing parameters

| Parameter | Guidance |
|-----------|----------|
| **Order** | Higher orders capture longer-range dependencies but require exponentially more training data. Start with order 1 or 2. |
| **Pseudocount** | Use a small value (0.001–0.1) to smooth rare transitions. Set to 0 if your training data is large enough to cover all transitions. |
| **Method** | `alias` is faster for large generation tasks. `binsrch` is simpler and sufficient for small outputs. |
| **Start** | For higher-order chains, the start string must be at least `order` characters long. |

---

## Related pages

- [CLI Reference — `biosea markov`](../cli.md#biosea-markov)
- [REST API — `POST /markov`](../api.md#post-markov)
- [Architecture — Execution flow](../architecture.md#execution-flow)
