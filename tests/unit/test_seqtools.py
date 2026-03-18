from bio_sea_pearl.api import hamming_distance, levenshtein_distance, kmer_counts


def test_hamming_distance():
    assert hamming_distance("ACGT", "AGGT") == 1
    assert hamming_distance("AAAA", "AAAA") == 0


def test_levenshtein_distance():
    assert levenshtein_distance("kitten", "sitting") == 3
    assert levenshtein_distance("flaw", "lawn") == 2


def test_kmer_counts():
    counts = kmer_counts("ACGTAC", 2)
    # Expected k-mers: AC, CG, GT, TA, AC
    assert counts["AC"] == 2
    assert counts["CG"] == 1
    assert counts["GT"] == 1
    assert counts["TA"] == 1
