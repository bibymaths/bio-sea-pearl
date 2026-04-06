---
title: Installation
---

# Installation

Bio Sea Pearl requires **Python 3.10+** and optionally **Perl 5** for legacy module support. The toolkit can be installed from PyPI, from source, or run as a Docker container.

---

## From PyPI

```bash
pip install bio-sea-pearl
```

This installs the `biosea` CLI, the Python library, and all Python dependencies. Perl scripts and scoring matrices are bundled in the wheel.

!!! note
    The PyPI package bundles the Perl scripts and scoring matrices inside the wheel. You do not need to clone the repository for basic usage.

---

## From source (development)

Clone the repository and install in editable mode with development dependencies:

```bash
git clone https://github.com/bibymaths/bio-sea-pearl.git
cd bio-sea-pearl
pip install -e ".[dev]"
```

This installs:

| Dependency  | Purpose                       |
|-------------|-------------------------------|
| `typer`     | CLI framework                 |
| `fastapi`   | REST API framework            |
| `uvicorn`   | ASGI server                   |
| `pydantic`  | Request/response validation   |
| `requests`  | HTTP client utilities         |
| `pytest`    | Test runner (dev only)        |
| `httpx`     | Async HTTP test client (dev)  |
| `ruff`      | Linter and formatter (dev)    |

---

## Perl runtime (optional)

Several modules — alignment, Markov chain generation, and the original SeqTools library — are implemented in Perl. If Perl is not available, the toolkit falls back to pure Python implementations where they exist (sequence distances, k-mer counting) or raises a clear error for Perl-only features (alignment, Markov).

=== "Ubuntu / Debian"

    ```bash
    sudo apt-get install -y perl
    ```

=== "macOS (Homebrew)"

    ```bash
    brew install perl
    ```

=== "Verify installation"

    ```bash
    perl --version
    ```

!!! tip
    The SeqTools Perl modules (`hamming.pm`, `levenshtein.pm`) use Inline C for performance. If you need the Perl path for these modules, ensure a C compiler is available. The Python fallbacks do not require a compiler.

---

## Docker

The simplest way to run the full stack (including Perl) is via Docker:

```bash
docker compose up --build -d
```

This builds an image from the included `Dockerfile`, installs Python and Perl inside a `python:3.12-slim` base, and starts the REST API on port **8000**.

To stop the container:

```bash
docker compose down
```

The container sets `BIOSEA_REPO_ROOT=/app`, which tells the wrapper layer where to find the Perl scripts and scoring matrices.

!!! info
    The Docker image is also published to `ghcr.io/bibymaths/bio-sea-pearl:latest`. You can pull it directly:

    ```bash
    docker pull ghcr.io/bibymaths/bio-sea-pearl:latest
    docker run -p 8000:8000 ghcr.io/bibymaths/bio-sea-pearl:latest
    ```

---

## Verify the installation

After installing, confirm the CLI is available:

```bash
biosea --help
```

Expected output:

```
Usage: biosea [OPTIONS] COMMAND [ARGS]...

Options:
  --help  Show this message and exit.

Commands:
  align     Align two FASTA sequences.
  bwt       Burrows–Wheeler / FM-index tools.
  markov    Generate a random walk from a Markov model.
  seqtools  Sequence distance and k-mer utilities.
```

To verify the REST API (requires uvicorn):

```bash
uvicorn api.server:app --host 127.0.0.1 --port 8000 &
curl http://127.0.0.1:8000/
```

---

## Next steps

Proceed to the [Quickstart](quickstart.md) to run your first commands, or jump to the [CLI Reference](../cli.md) for full command documentation.
