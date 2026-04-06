---
title: Extending the Toolkit
---

# Extending the Toolkit

This guide explains how to add new functionality to Bio Sea Pearl, whether you are porting an existing Perl module to Python, wrapping a new external tool, or adding an entirely new subsystem.

---

## Project conventions

Before writing code, understand the layered structure:

```
Interface layer    →  cli.py  /  api/server.py  /  direct imports
API layer          →  src/bio_sea_pearl/api/*.py
Implementation     →  src/bio_sea_pearl/bwt/  or  seqtools_py/  (Python-native)
                      src/bio_sea_pearl/perl_wrappers/*.py       (Perl bridge)
Legacy scripts     →  alignment/  markov/  seqtools/
```

Every new feature should be accessible through the **API layer** first. The CLI and REST API are thin wrappers that call API functions. Do not put logic directly in `cli.py` or `server.py`.

---

## Adding a new Python-native function

**Example:** Adding a GC-content calculator.

### 1. Create the implementation module

```
src/bio_sea_pearl/seqtools_py/gc_content.py
```

```python
def gc_content(seq: str) -> float:
    """Return the GC content of a DNA sequence as a fraction."""
    if not seq:
        raise ValueError("Empty sequence")
    seq = seq.upper()
    gc = sum(1 for c in seq if c in "GC")
    return gc / len(seq)
```

### 2. Export from the API layer

Add to `src/bio_sea_pearl/api/seqtools_api.py`:

```python
from bio_sea_pearl.seqtools_py.gc_content import gc_content as _gc_content_py

def gc_content(seq: str) -> float:
    return _gc_content_py(seq)
```

Export from `src/bio_sea_pearl/api/__init__.py`:

```python
from .seqtools_api import gc_content
```

### 3. Add a CLI subcommand

In `src/bio_sea_pearl/cli.py`, register a new command under the `seqtools` group:

```python
@seqtools_app.command()
def gc(sequence: str = typer.Argument(...)):
    """Compute GC content of a sequence."""
    from bio_sea_pearl.api import gc_content
    result = gc_content(sequence)
    typer.echo(f"GC content: {result:.4f}")
```

### 4. Add a REST API endpoint

In `api/server.py`, add a Pydantic model and endpoint:

```python
class GCRequest(BaseModel):
    sequence: str

@app.post("/gc")
def gc_endpoint(req: GCRequest):
    result = gc_content(req.sequence)
    return {"sequence": req.sequence, "gc_content": result}
```

### 5. Add tests

Create `tests/unit/test_gc_content.py`:

```python
from bio_sea_pearl.api import gc_content

def test_gc_content_basic():
    assert gc_content("ACGT") == 0.5

def test_gc_content_all_gc():
    assert gc_content("GCGC") == 1.0
```

---

## Wrapping a new Perl script

If you have a Perl script that you want to expose through the Python interface:

### 1. Place the script

Put it in the appropriate directory (e.g. `seqtools/` for sequence utilities, or create a new top-level directory for a new subsystem).

### 2. Create a wrapper

Add a new file in `src/bio_sea_pearl/perl_wrappers/`, following the pattern of existing wrappers:

```python
import os
import subprocess
from pathlib import Path

def _repo_root() -> Path:
    env = os.environ.get("BIOSEA_REPO_ROOT")
    if env:
        return Path(env)
    return Path(__file__).resolve().parents[3]

def run_my_tool(arg1: str, arg2: str) -> str:
    root = _repo_root()
    script = root / "my_tool" / "bin" / "my_script.pl"
    result = subprocess.run(
        ["perl", str(script), arg1, arg2],
        capture_output=True, text=True, cwd=str(root)
    )
    if result.returncode != 0:
        raise RuntimeError(result.stderr.strip())
    return result.stdout.strip()
```

### 3. Export and expose

Follow the same pattern as above: add an API function, export it, add CLI and REST endpoints.

!!! tip
    Always strip ANSI escape codes from Perl output using `term_utils.strip_ansi()`. Many Perl scripts emit coloured terminal output that corrupts JSON responses.

---

## Porting a Perl module to Python

The existing `seqtools_py/` directory demonstrates the porting pattern:

1. **Write a pure Python implementation** with the same semantics as the Perl module.
2. **Add it to `seqtools_py/`** (or a new package under `src/bio_sea_pearl/`).
3. **Update the API layer** to try the Python version first and fall back to Perl.

The fallback logic in `seqtools_api.py` handles this automatically:

```python
try:
    return _python_implementation(args)
except ValueError:
    raise  # Validation errors propagate immediately
except Exception:
    return _perl_fallback(args)  # Other failures fall back
```

This pattern means you can deploy without Perl once all required functions have Python ports.

---

## Build and packaging

The project uses **Hatchling** as its build backend. Key configuration in `pyproject.toml`:

```toml
[tool.hatch.build.targets.wheel]
packages = ["src/bio_sea_pearl"]

[tool.hatch.build.targets.wheel.force-include]
"api" = "api"
"alignment" = "alignment"
"markov" = "markov"
"seqtools" = "seqtools"
```

If you add a new top-level directory that needs to be included in the wheel, add a `force-include` entry.

---

## Running tests

```bash
# Install dev dependencies
pip install -e ".[dev]"

# Run the full test suite
pytest

# Run a specific test file
pytest tests/unit/test_bwt.py

# Run with verbose output
pytest -v
```

The test suite uses `@pytest.mark.skipif` to skip tests that require Perl when Perl is not installed. Follow this pattern for any new tests that depend on external tools.

---

## Repository root resolution

The wrapper layer locates the repository root in two ways:

1. **`BIOSEA_REPO_ROOT` environment variable** — set this in Docker containers or CI environments where the repo is at a non-standard path.
2. **Path traversal** — `Path(__file__).resolve().parents[3]` walks three levels up from the wrapper file, which reaches the repository root from `src/bio_sea_pearl/perl_wrappers/`.

If you move wrapper files to a different directory depth, update the `parents[3]` index accordingly.

---

## Related pages

- [Architecture](../architecture.md) — understand the full layer model
- [Troubleshooting](troubleshooting.md) — common issues during development
