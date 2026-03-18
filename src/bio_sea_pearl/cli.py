"""Command-line interface for the Bio Sea Pearl toolkit.

This CLI exposes the major functionalities of the project via distinct
subcommands.  It uses Typer to define commands, parse arguments and
provide help messages.  Each subcommand delegates to the API layer,
which in turn calls the appropriate wrappers or native Python code.

Usage examples:

    biosea align --matrix alignment/scoring/BLOSUM62.mat seq1.fa seq2.fa
    biosea markov --fasta data.fa --length 100 --start A
    biosea seqtools hamming ACCT AGGT
    biosea seqtools kmer --k 3 ACGTACGT
    biosea bwt search --sequence ACGTACGT --pattern GT
"""

from __future__ import annotations

import json
import typer

from .api import (
    align_sequences,
    generate_walk,
    hamming_distance,
    levenshtein_distance,
    kmer_counts,
    build_fm_index,
    search_fm_index,
)

app = typer.Typer(help="Bio Sea Pearl unified command-line interface")


@app.command()
def align(
        fasta1: str = typer.Argument(..., help="Path to the first FASTA file."),
        fasta2: str = typer.Argument(..., help="Path to the second FASTA file."),
        matrix: str = typer.Option(
            None,
            "--matrix",
            "-m",
            help="Optional substitution matrix file (e.g. alignment/scoring/BLOSUM62.mat).",
        ),
        mode: str = typer.Option(
            "global",
            "--mode",
            "-M",
            help="Alignment mode: global, local or lcs.",
        ),
) -> None:
    """Align two sequences stored in FASTA files."""
    result = align_sequences(fasta1, fasta2, matrix=matrix, mode=mode)
    typer.echo(result.strip())


@app.command()
def markov(
        fasta: str = typer.Option(..., help="Path to the FASTA file containing training sequences."),
        length: int = typer.Option(..., help="Length of the random walk to generate."),
        start: str = typer.Option("A", help="Starting state for the Markov chain."),
        order: int = typer.Option(1, help="Order of the Markov chain."),
        method: str = typer.Option("alias", help="Sampling method: alias or binsrch."),
        pseudocount: int = typer.Option(0, help="Pseudocount added to each transition."),
) -> None:
    """Generate a random sequence from a Markov model built from a FASTA file."""
    walk = generate_walk(fasta, length, start=start, order=order, method=method, pseudocount=pseudocount)
    typer.echo(walk.strip())


seqtools_app = typer.Typer(help="SeqTools commands")
app.add_typer(seqtools_app, name="seqtools")


@seqtools_app.command()
def hamming(seq1: str, seq2: str) -> None:
    """Compute the Hamming distance between two equal‑length sequences."""
    try:
        d = hamming_distance(seq1, seq2)
        typer.echo(str(d))
    except Exception as exc:
        typer.echo(f"Error computing Hamming distance: {exc}", err=True)


@seqtools_app.command()
def levenshtein(seq1: str, seq2: str) -> None:
    """Compute the Levenshtein distance between two sequences."""
    try:
        d = levenshtein_distance(seq1, seq2)
        typer.echo(str(d))
    except Exception as exc:
        typer.echo(f"Error computing Levenshtein distance: {exc}", err=True)


@seqtools_app.command()
def kmer(seq: str, k: int = typer.Option(..., help="Length of k‑mers to count.")) -> None:
    """Count k‑mers in the given sequence and output as JSON."""
    try:
        counts = kmer_counts(seq, k)
        typer.echo(json.dumps(counts))
    except Exception as exc:
        typer.echo(f"Error computing k‑mer counts: {exc}", err=True)


bwt_app = typer.Typer(help="BWT and FM‑index commands")
app.add_typer(bwt_app, name="bwt")


@bwt_app.command()
def search(
        sequence: str = typer.Option(..., help="Sequence to build the FM‑index from."),
        pattern: str = typer.Option(..., help="Pattern to search for."),
) -> None:
    """Search for a pattern in a sequence using the FM‑index and return positions."""
    index = build_fm_index(sequence)
    positions = search_fm_index(index, pattern)
    # Print positions as a space‑separated list
    typer.echo(" ".join(map(str, positions)))


if __name__ == "__main__":
    app()
