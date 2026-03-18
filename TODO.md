# Project Roadmap

This repository currently combines Python and Perl code to deliver sequence alignment, Markov chain modelling, sequence utilities and text indexing.  To make the toolkit easier to use and extend, the following roadmap outlines concrete steps for future development.

## Unified Command‑Line Interface

* Develop a Python entry point (e.g. `cli.py`) that exposes subcommands for each major feature: `align`, `markov`, `seqtools` and `bwt`.  Each subcommand should accept consistent arguments and delegate work to the underlying implementation.
* Wrap existing Perl scripts using Python’s `subprocess` module so that users can invoke `dotplot.pl`, `randomwalk.pl` and other utilities without leaving Python.  Standardise input and output formats (e.g. FASTA, plain text).
* Provide helpful `--help` messages and examples for every subcommand.

## Port Perl Modules to Python

* Re‑implement the `SeqTools` library in Python.  The Hamming and Levenshtein distance functions can be implemented using pure Python or NumPy for acceleration.  K‑mer generation and Boyer–Moore search can be ported and optimised using built‑in libraries.
* Port the Markov chain simulators (`MarkovChainFirstOrder.pm` and `MarkovChainHigherOrder.pm`) to Python.  Maintain feature parity (alias and binary search sampling, pseudocounts) and validate behaviour against the Perl version using unit tests.
* Gradually phase out Perl dependencies in favour of pure‑Python implementations while maintaining backwards compatibility via the CLI wrapper.

## API Layer

* Build a RESTful API using FastAPI or Flask to expose alignment, Markov modelling and BWT search as HTTP endpoints.  This will allow the toolkit to serve as a backend service for web applications or command‑line clients.
* Design the API so that back‑ends can be swapped: initially call the existing Perl scripts via `subprocess`, then switch to Python implementations when available.  Use dependency injection to decouple the interface from the implementation language.
* Document the API endpoints with OpenAPI/Swagger and include usage examples in the MkDocs documentation.

## Continuous Integration and Testing

* Expand the minimal test data in `tests/data/` into comprehensive unit tests.  Use `pytest` for Python code and `Test::More`/`prove` for Perl scripts.  Tests should cover typical usage scenarios as well as edge cases (e.g. empty sequences, invalid characters, high orders for Markov models).
* Update the GitHub Actions workflow to run both the Python tests and the Perl test suite on each push and pull request.  Include code coverage reports to track progress.
* Add static analysis tools such as `ruff` for Python and `perlcritic` for Perl to enforce coding standards.

## Documentation

* Flesh out the MkDocs documentation (see the `docs/` directory) with detailed tutorials, API reference, and architectural diagrams.  Include step‑by‑step examples for each module.
* Add diagrams illustrating the dynamic‑programming matrices used in alignment and the data structures underlying the FM‑index.  Use `mkdocs-mermaid2-plugin` if diagrams need to be rendered from Mermaid syntax.
* Provide guidance on installing dependencies on various platforms (Linux, macOS, Windows) and using the CLI wrapper.

## Long‑Term Goals

* Investigate replacing performance‑critical routines with compiled extensions (e.g. Cython, Rust) for speed‑up, especially for large-scale alignments and Markov sampling.
* Extend support to additional alphabets such as RNA (including ambiguous IUPAC codes) and to amino‑acid properties for weighted alignments.
* Apply Markov modelling to real datasets (e.g. cancer genomic sequences) and integrate with pharmacokinetic or network‑based models to support precision medicine applications.

This roadmap is intended to guide development toward a cohesive, modern toolkit.  Community contributions and suggestions are welcome through pull requests and issue discussions.
