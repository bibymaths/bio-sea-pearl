# Project Roadmap

This repository combines legacy **Perl** scripts with newer **Python** modules to provide sequence alignment, Markov
chain modelling, sequence utilities, and FM‑index search. To evolve this toolkit into a unified, maintainable, and
Python‑centric platform, the following roadmap outlines actionable tasks and a clear strategy for bridging the language
gap.

## 1. Unify the Command‑Line Interface

* **Create a single CLI entry point** (`biosea`) with subcommands for `align`, `markov`, `seqtools` and `bwt`. Use a
  modern Python CLI framework such as Typer or argparse.
* **Wrap existing Perl scripts** using Python’s `subprocess` module so that users can invoke tools like `align.pl`,
  `randomwalk.pl`, `run.pl`, `dotplot.pl`, and `traceplot.pl` transparently. Standardise input and output formats (FASTA
  in, plain text or JSON out) so that users do not need to know which language is performing the work.
* **Document usage** with helpful `--help` messages and examples. Include CLI examples in the MkDocs documentation.

## 2. Bridge Perl and Python via Wrappers

* **Perl wrappers**: Implement a Python `perl_wrappers` package that provides functions such as `run_alignment`,
  `generate_walk`, `hamming_distance`, `levenshtein_distance`, `kmer_count`, and `boyer_moore_search`. These functions
  should call the appropriate Perl scripts/modules through `subprocess.run` and return parsed results as Python data
  structures. Capture both stdout and stderr to provide informative error messages.
* **Python API layer**: Build a Python `api` package that exposes high‑level functions (`align_sequences`,
  `generate_markov_walk`, `hamming_distance_py`, `levenshtein_distance_py`, `count_kmers`, `fm_search`) which delegate
  to either the wrappers (for Perl implementations) or pure‑Python modules when available. Use type hints for clarity
  and maintainability.
* **Swappable back‑ends**: Design the API so that the underlying implementation (Perl vs. Python) can be selected via
  configuration or dependency injection. This allows gradual migration away from Perl without breaking the public
  interface.

## 3. Gradually Port Perl Modules to Python

* **SeqTools port**: Re‑implement the `SeqTools` library in Python. Start with the most commonly used functions—Hamming
  distance, Levenshtein distance and k‑mer generation—using pure Python or NumPy for acceleration. Ensure that the
  Python implementations produce identical outputs to the Perl version by writing unit tests that compare results across
  languages.
* **Markov chains**: Port the Markov chain simulators (`MarkovChainFirstOrder.pm` and `MarkovChainHigherOrder.pm`) to
  Python. Preserve features such as pseudocounts and alias sampling. Validate parity with the Perl implementations via
  integration tests.
* **Deprecate Perl gradually**: As Python ports reach feature parity, update the wrappers to favour Python
  implementations while retaining the ability to fall back to Perl for backward compatibility.

## 4. Expose a RESTful API

* **Service layer**: Build a FastAPI server (`api/server.py`) that exposes alignment, Markov modelling, distance
  calculations, k‑mer counting, and BWT search as JSON‑driven HTTP endpoints. Keep route handlers light; delegate
  business logic to the API layer.
* **OpenAPI documentation**: Use FastAPI’s automatic OpenAPI generation to provide interactive API documentation.
  Incorporate examples into the MkDocs site so that users can try endpoints directly from the docs.
* **Asynchronous design**: Write endpoints as `async` functions and ensure that subprocess calls are non‑blocking (e.g.
  using `run_in_threadpool`) to allow concurrent processing.

## 5. Testing and Continuous Integration

* **Unit tests**: Expand the `tests/unit` suite to cover all Python modules and wrapper functions. Use `pytest` for
  Python and `prove`/`Test::More` for any remaining Perl tests. Include both positive and negative cases (empty inputs,
  invalid characters, mismatched lengths).
* **Integration tests**: Add `tests/integration` to verify that the CLI and REST API perform end‑to‑end operations
  correctly. Spin up the FastAPI server in the test context and make HTTP requests to assert expected responses.
* **GitHub Actions**: Ensure the CI pipeline runs both Python (`pytest`) and Perl syntax checks on each commit. Add
  coverage reporting and static analysis (e.g. `ruff` for Python, `perlcritic` for Perl) to enforce code quality.

## 6. Documentation Enhancements

* **MkDocs**: Build out the documentation in the `docs/` directory. Write tutorials for using each CLI subcommand, API
  endpoint, and library function. Include diagrams illustrating dynamic programming for alignment and the FM‑index
  structure. The Material for MkDocs configuration should specify `site_name`, `site_url`, and the Material theme as
  recommended in the official docs【190512879997651†L195-L204】.
* **Architecture diagram**: Provide a high‑level diagram showing how the CLI, wrappers, API layer, FastAPI server, and
  underlying Perl/Python implementations interact. Explain the rationale for the bridging layer and how the system can
  switch between back‑ends.
* **Setup instructions**: Document installation on Linux, macOS, and Windows. Describe how to install Perl dependencies
  and how to install the Python package via Poetry. Explain how to run the CLI, start the API server, and contribute to
  the codebase.

## 7. Configuration and Logging

* **Configuration files**: Add a `config/` directory with a `.env` file for environment variables and a `logging.yaml`
  file for structured logging. Document how to override configurations (e.g. select Perl vs. Python back‑end, set
  logging levels).
* **Logging**: Use human‑readable logs for CLI output and JSON‑formatted logs for the API server. Log key events (
  invocations, errors, runtime durations) to facilitate debugging.

## 8. Long‑Term and Advanced Goals

* **Performance optimisation**: Investigate using Cython or Rust to accelerate critical routines once the pure‑Python
  implementations are stable. Provide optional compiled extensions that users can install for large datasets.
* **Additional alphabets**: Extend sequence support to RNA (including ambiguous IUPAC codes) and other alphabets. Add
  weighted alignments based on amino‑acid properties or custom scoring matrices.
* **Real‑world applications**: Integrate Markov modelling and alignment with pharmacokinetic or network‑based models for
  cancer systems biology and precision medicine use‑cases. Build example workflows demonstrating the toolkit’s
  application to multi‑omics data.

This roadmap aims to guide development toward a cohesive, modern, and flexible toolkit. Contributions are welcome via
pull requests and issue discussions.
