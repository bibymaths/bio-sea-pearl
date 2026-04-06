"""Packaging and import smoke tests."""

from __future__ import annotations


def test_import_top_level():
    import bio_sea_pearl  # noqa: F401


def test_import_api_functions():
    from bio_sea_pearl.api import (  # noqa: F401
        align_sequences,
        generate_walk,
        hamming_distance,
        levenshtein_distance,
        kmer_counts,
        build_fm_index,
        search_fm_index,
    )


def test_import_server_app():
    from api.server import app  # noqa: F401


def test_app_is_fastapi_instance():
    from api.server import app
    from fastapi import FastAPI

    assert isinstance(app, FastAPI)


def test_app_has_title():
    from api.server import app

    assert app.title == "Bio Sea Pearl API"


def test_api_all_exports():
    """The __all__ list should contain exactly the expected public functions."""
    from bio_sea_pearl.api import __all__ as exported

    expected = {
        "align_sequences",
        "generate_walk",
        "hamming_distance",
        "levenshtein_distance",
        "kmer_counts",
        "build_fm_index",
        "search_fm_index",
    }
    assert set(exported) == expected
