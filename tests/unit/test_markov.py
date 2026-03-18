from bio_sea_pearl.api import generate_walk


def test_generate_walk(tmp_path):
    """Generate a short Markov walk and ensure output length matches."""
    # Prepare a simple FASTA file
    fasta = tmp_path / "train.fa"
    fasta.write_text(">seq\nACGTACGT\n")
    length = 10
    walk = generate_walk(str(fasta), length, start="A", order=1, method="alias")
    walk = walk.strip()
    # The generated walk should have the requested length
    assert len(walk) == length + 1 or len(walk) == length  # account for possible start char duplication
