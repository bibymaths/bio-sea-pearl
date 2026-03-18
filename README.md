To get the project running:

1. **Unpack** the ZIP you downloaded. You should see a `bio-sea-pearl-main` directory containing the source code, tests, documentation, and configuration files.

2. **Install dependencies**. The project uses Python≥3.13 and Poetry for dependency management. From inside the `bio-sea-pearl-main` directory:

   ```bash
   pip install poetry   # if you don’t already have Poetry
   poetry install       # installs the package and its dependencies
   ```

   If you prefer plain `pip`, you can run `pip install .`, but Poetry simplifies handling the lockfile.

3. **Run the command‑line tool**. Poetry installs an entry point named `biosea`. Get help:

   ```bash
   poetry run biosea --help
   ```

   Example commands:

   * Align two FASTA files:

     ```bash
     poetry run biosea align seq1.fa seq2.fa --matrix alignment/scoring/blosum62.mat --mode global
     ```
   * Generate a Markov random walk:

     ```bash
     poetry run biosea markov --fasta sample.fa --length 100 --start A --order 1 --method alias
     ```
   * Sequence utilities:

     ```bash
     poetry run biosea seqtools hamming ACGT AGGT       # Hamming distance
     poetry run biosea seqtools levenshtein kitten sitting  # Levenshtein distance
     poetry run biosea seqtools kmer ACGTACGT --k 3     # k‑mer counts in JSON
     ```
   * BWT/FM‑index search:

     ```bash
     poetry run biosea bwt search --sequence ACGTACGT --pattern CGT
     ```

4. **Run the REST API**. A FastAPI server lives in `api/server.py`. Launch it with Uvicorn:

   ```bash
   poetry run uvicorn api.server:app --reload
   ```

   The API exposes endpoints such as `/align`, `/markov`, `/distance`, `/kmer` and `/bwt/search`. Use a tool like `curl` or Postman to POST JSON to these endpoints; see `docs/api.md` for request/response formats.

5. **View the documentation**. The repository includes an MkDocs site. Serve it locally:

   ```bash
   pip install mkdocs-material
   mkdocs serve
   ```

   Then browse to the local URL (typically `http://localhost:8000`) to read the tutorials, architecture overview, and API documentation.
