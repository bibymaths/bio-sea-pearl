#!/usr/bin/env perl
use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin";

use SeqTools::Fasta qw(read_fasta);
use SeqTools::Kmer qw(count);
use SeqTools::Hamming qw(distance);
use SeqTools::Levenshtein qw(distance);
use SeqTools::BoyerMoore qw(search);

# Alias levenshtein distance as lev
*lev = \&SeqTools::Levenshtein::distance;

# Die with usage message if insufficient arguments
unless (@ARGV) {
    die <<"USAGE";
Usage:
  \$0 kmer <fasta_file> <k>
  \$0 motif <fasta_file> <pattern>
  \$0 hamming <fasta1> <fasta2>
  \$0 levenshtein <fasta1> <fasta2>
USAGE
}

my ($cmd, @args) = @ARGV;

if ($cmd eq 'kmer') {
    die "Need FASTA file & k
" unless @args == 2;
    my ($fasta, $k) = @args;
    my $seq = read_fasta($fasta);
    my $h   = count($seq, $k);
    print "$_	$h->{$_}
" for sort keys %$h;

} elsif ($cmd eq 'motif') {
    die "Need FASTA file & pattern
" unless @args == 2;
    my ($fasta, $pat) = @args;
    my $seq  = read_fasta($fasta);
    my $hits = search($seq, $pat);
    if (@$hits) {
        print "Matches at positions:\n", join("\n", @$hits), "
";
    } else {
        print "No matches found
";
    }

} elsif ($cmd eq 'hamming') {
    die "Need two FASTA files
" unless @args == 2;
    my ($fa1, $fa2) = @args;
    my $s1 = read_fasta($fa1);
    my $s2 = read_fasta($fa2);
    print distance($s1, $s2), "
";

} elsif ($cmd eq 'levenshtein') {
    die "Need two FASTA files
" unless @args == 2;
    my ($fa1, $fa2) = @args;
    my $s1 = read_fasta($fa1);
    my $s2 = read_fasta($fa2);
    print lev($s1, $s2), "
";

} else {
    die "Unknown command '$cmd'
";
}