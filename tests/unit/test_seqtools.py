"""Tests for sequence-tools API functions: Hamming, Levenshtein, k-mer counts."""

from __future__ import annotations

import pytest

from bio_sea_pearl.api import hamming_distance, levenshtein_distance, kmer_counts


# ── Hamming distance ───────────────────────────────────────────────

class TestHammingDistance:
    def test_identical_strings(self):
        assert hamming_distance("AAAA", "AAAA") == 0

    def test_single_mismatch(self):
        assert hamming_distance("ACGT", "AGGT") == 1

    def test_all_mismatches(self):
        assert hamming_distance("AAAA", "TTTT") == 4

    def test_single_char(self):
        assert hamming_distance("A", "A") == 0
        assert hamming_distance("A", "T") == 1

    def test_empty_strings(self):
        assert hamming_distance("", "") == 0

    def test_unequal_lengths_raises(self):
        with pytest.raises((ValueError, Exception)):
            hamming_distance("ABC", "AB")


# ── Levenshtein distance ──────────────────────────────────────────

class TestLevenshteinDistance:
    def test_kitten_sitting(self):
        assert levenshtein_distance("kitten", "sitting") == 3

    def test_flaw_lawn(self):
        assert levenshtein_distance("flaw", "lawn") == 2

    def test_identical_strings(self):
        assert levenshtein_distance("hello", "hello") == 0

    def test_empty_to_nonempty(self):
        assert levenshtein_distance("", "abc") == 3

    def test_nonempty_to_empty(self):
        assert levenshtein_distance("abc", "") == 3

    def test_both_empty(self):
        assert levenshtein_distance("", "") == 0

    def test_single_insert(self):
        assert levenshtein_distance("abc", "abcd") == 1

    def test_single_delete(self):
        assert levenshtein_distance("abcd", "abc") == 1

    def test_single_substitution(self):
        assert levenshtein_distance("abc", "aXc") == 1


# ── K-mer counts ─────────────────────────────────────────────────

class TestKmerCounts:
    def test_basic_2mers(self):
        counts = kmer_counts("ACGTAC", 2)
        assert counts["AC"] == 2
        assert counts["CG"] == 1
        assert counts["GT"] == 1
        assert counts["TA"] == 1

    def test_k_equals_seq_length(self):
        counts = kmer_counts("ACG", 3)
        assert counts == {"ACG": 1}

    def test_single_character_kmers(self):
        counts = kmer_counts("AABBA", 1)
        assert counts["A"] == 3
        assert counts["B"] == 2

    def test_homopolymer(self):
        counts = kmer_counts("AAAA", 2)
        assert counts == {"AA": 3}

    def test_k_greater_than_length_raises(self):
        with pytest.raises((ValueError, Exception)):
            kmer_counts("AC", 5)

    def test_k_zero_raises(self):
        with pytest.raises((ValueError, Exception)):
            kmer_counts("ACGT", 0)
