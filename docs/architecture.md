---
title: Architecture
---

# Architecture

The Bio Sea Pearl project is intentionally modular to accommodate both
Python and Perl code.  The layers of the architecture are outlined below.

## Repository layout

* **src/bio_sea_pearl** – Core Python package.
  * **bwt/** – Burrows–Wheeler Transform and FM‑index implementations.
  * **perl_wrappers/** – Functions that invoke existing Perl scripts via
    `subprocess`.  These isolate language interop details from the rest
    of the code.
  * **seqtools_py/** – Pure Python ports of select SeqTools modules
    (Hamming, Levenshtein, k‑mer counting).  These provide native
    implementations for performance and ease of deployment.
  * **api/** – High‑level Python interface.  Functions here call into
    the wrappers or native code and return Python data structures.
  * **cli.py** – Typer‑based command‑line interface exposing the API.

* **alignment/**, **markov/** and **seqtools/** – Legacy Perl and Python
  implementations.  These remain unchanged and are invoked via the wrappers.

* **api/server.py** – FastAPI application exposing toolkit functionality
  over HTTP.

## Data flow

1. **User input** – The user interacts with the CLI (`biosea`) or sends
   HTTP requests to the FastAPI server.
2. **API layer** – The CLI and server call functions in
   `bio_sea_pearl.api`.  This layer defines clear interfaces for
   alignment, Markov chain generation, sequence metrics and text search.
3. **Wrappers and native code** – The API layer either invokes
   pure‑Python code (e.g. the BWT or SeqTools ports) or calls the legacy
   Perl scripts using the wrappers.
4. **External scripts** – Perl scripts and modules perform heavy
   computation such as dynamic programming or Markov sampling.  Results
   are captured and returned through Python.

This separation allows the project to evolve towards a fully Python
implementation while preserving existing functionality.