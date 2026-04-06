---
title: Troubleshooting
---

# Troubleshooting

This page covers common issues encountered when installing, running, or developing with Bio Sea Pearl. Each section describes the symptom, cause, and resolution.

---

## Installation issues

### `biosea: command not found`

**Cause:** The package is not installed, or the install location is not in your `PATH`.

**Resolution:**

```bash
pip install bio-sea-pearl
# or, from a local clone:
pip install -e ".[dev]"
```

If you installed with `--user`, ensure `~/.local/bin` is in your `PATH`:

```bash
export PATH="$HOME/.local/bin:$PATH"
```

### `ModuleNotFoundError: No module named 'bio_sea_pearl'`

**Cause:** The package is not installed in the active Python environment, or you are using a different Python than expected.

**Resolution:**

```bash
# Verify which Python is active
which python
python -m pip show bio-sea-pearl
```

If the package is not listed, install it. If you are in a virtual environment, make sure it is activated.

---

## Perl-related failures

### `FileNotFoundError: [Errno 2] No such file or directory: 'perl'`

**Cause:** Perl is not installed on the system.

**Resolution:** Install Perl:

=== "Ubuntu / Debian"

    ```bash
    sudo apt-get install -y perl
    ```

=== "macOS"

    ```bash
    brew install perl
    ```

=== "Docker"

    Use the provided `Dockerfile`, which installs Perl automatically.

!!! info
    Perl is required for alignment and Markov commands. Sequence tools (hamming, levenshtein, k-mer) and BWT commands work without Perl.

### `Can't locate MarkovChainHigherOrder.pm in @INC`

**Cause:** The Perl module search path does not include `markov/lib/`.

**Resolution:** The wrapper sets `PERL5LIB` automatically based on the repository root. Verify that `BIOSEA_REPO_ROOT` is set correctly:

```bash
echo $BIOSEA_REPO_ROOT
```

If running outside Docker or a standard install, set it explicitly:

```bash
export BIOSEA_REPO_ROOT=/path/to/bio-sea-pearl
```

### `Can't locate Inline.pm` or `Inline C compilation failed`

**Cause:** The Perl `Inline::C` module is not installed, or no C compiler is available. This affects `hamming.pm` and `levenshtein.pm` in the Perl SeqTools.

**Resolution:** This only matters if you are using the Perl fallback path. The Python implementations (`seqtools_py`) do not require a C compiler. If the Python port is available, the system bypasses Perl automatically.

If you need the Perl path:

```bash
cpan Inline::C
# or
sudo apt-get install -y build-essential
```

---

## Alignment issues

### `align.py not found` or alignment returns empty output

**Cause:** The wrapper cannot locate `alignment/bin/align.py` relative to the repository root.

**Resolution:**

1. Check that `BIOSEA_REPO_ROOT` points to the correct directory:

    ```bash
    ls $BIOSEA_REPO_ROOT/alignment/bin/align.py
    ```

2. If running from a pip-installed wheel (not a local clone), ensure the wheel was built with `force-include` entries intact. The `alignment/` directory should be bundled.

### `Matrix file not found`

**Cause:** The specified scoring matrix path does not exist.

**Resolution:** Scoring matrices are in `alignment/scoring/`. Use the full relative path:

```bash
biosea align seq1.fa seq2.fa --matrix alignment/scoring/BLOSUM62.mat
```

If no matrix is specified, the alignment script uses its built-in default. The wrapper falls back to `blosum62.mat` in the scoring directory.

---

## REST API issues

### `Connection refused` on `localhost:8000`

**Cause:** The API server is not running.

**Resolution:**

```bash
uvicorn api.server:app --host 0.0.0.0 --port 8000
```

Or via Docker:

```bash
docker compose up -d
```

### `500 Internal Server Error` from `/align` or `/markov`

**Cause:** Usually a Perl subprocess failure. The error detail will include the stderr from the Perl script.

**Resolution:**

1. Check the response body for the `detail` field — it contains the error message.
2. Verify Perl is installed and `BIOSEA_REPO_ROOT` is set (especially in Docker).
3. Test the Perl script directly:

    ```bash
    perl markov/bin/randomwalk.pl --fasta test.fa --length 10 --start A
    ```

### `400 Bad Request` from `/distance`

**Cause:** Validation error — for example, requesting Hamming distance on strings of different lengths.

**Resolution:** Check the `detail` field in the response. For Hamming distance, both sequences must be the same length. Use `levenshtein` for different-length strings.

---

## Docker issues

### Container starts but API returns errors

**Cause:** The `BIOSEA_REPO_ROOT` environment variable may not match the actual file layout in the container.

**Resolution:** The default `docker-compose.yml` sets `BIOSEA_REPO_ROOT=/app`. Verify:

```bash
docker exec biosea ls /app/alignment/bin/align.py
docker exec biosea ls /app/markov/bin/randomwalk.pl
```

### Port conflict on 8000

**Cause:** Another service is using port 8000.

**Resolution:** Change the port mapping in `docker-compose.yml`:

```yaml
ports:
  - "9000:8000"  # Map host port 9000 to container port 8000
```

---

## Test failures

### Tests skip with `Perl not found`

**Cause:** Expected behaviour. Tests that require Perl are decorated with `@pytest.mark.skipif(shutil.which("perl") is None, ...)` and are skipped when Perl is not installed.

**Resolution:** Install Perl to run the full test suite. The BWT and seqtools Python tests run without Perl.

### `ImportError` in tests

**Cause:** The package is not installed in the test environment.

**Resolution:**

```bash
pip install -e ".[dev]"
pytest
```

The `pyproject.toml` sets `pythonpath = ["src", "."]` for pytest, so tests should find the package when installed in editable mode.

---

## Common environment variables

| Variable | Purpose | Default |
|----------|---------|---------|
| `BIOSEA_REPO_ROOT` | Root directory for locating scripts and scoring matrices | Auto-detected from wrapper file location |
| `PERL5LIB` | Perl module search path (set by wrappers internally) | Not required for user configuration |

---

## Getting help

- **GitHub Issues:** [github.com/bibymaths/bio-sea-pearl/issues](https://github.com/bibymaths/bio-sea-pearl/issues)
- **Source code:** [github.com/bibymaths/bio-sea-pearl](https://github.com/bibymaths/bio-sea-pearl)

---

## Related pages

- [Installation](../getting-started/installation.md) — setup instructions
- [Architecture](../architecture.md) — understanding the wrapper layer and dispatch
- [Extending the Toolkit](extending.md) — developer workflow
