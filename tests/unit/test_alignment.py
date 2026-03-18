import os
from bio_sea_pearl.api import align_sequences


def test_align_sequences(tmp_path):
    """Basic check that alignment returns non-empty output for identical sequences."""
    # Create two identical small FASTA files
    fasta1 = tmp_path / "seq1.fa"
    fasta2 = tmp_path / "seq2.fa"
    content = ">seq\nACGT\n"
    fasta1.write_text(content)
    fasta2.write_text(content)
    result = align_sequences(str(fasta1), str(fasta2))
    # The output should contain at least one line of alignment
    assert result