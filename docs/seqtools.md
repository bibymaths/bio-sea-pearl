---
title: SeqTools Library
---

# SeqTools Library

`seqtools` is a Perl collection of sequence and string utilities designed to support bioinformatics pipelines.  Rather than reinventing core algorithms in every script, it exposes reusable modules for pattern searching, distance calculations and k‑mer analysis.

## Modules

* **SeqTools/fasta.pm** – minimal FASTA parser.  Reads headers and sequence lines, returning an array of sequence records.
* **SeqTools/hamming.pm** – computes the Hamming distance between two equal‑length sequences.  Useful for quick comparison of SNP profiles.
* **SeqTools/levenshtein.pm** – calculates the Levenshtein (edit) distance allowing insertions, deletions and substitutions.  This is slower than Hamming but handles sequences of different lengths.
* **SeqTools/kmer.pm** – generates k‑mer frequency tables from sequences.  K‑mer counts are essential for genome assembly and complexity analysis.
* **SeqTools/boyermoore.pm** – implements the Boyer–Moore pattern search algorithm, which skips ahead based on mismatches for sublinear time complexity.

These modules are intended to be `use`d in Perl scripts:

```perl
use SeqTools::hamming qw(distance);
my $d = distance('ACGT', 'AGGT');
print "Hamming distance = $d\n";

use SeqTools::levenshtein qw(edit_distance);
print edit_distance('kitten', 'sitting');
```

Future versions of this toolkit may port many of these functions to Python to enable a unified command‑line interface and to leverage modern Python libraries such as NumPy for accelerated computation.
