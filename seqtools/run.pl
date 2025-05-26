#!/usr/bin/env perl
use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin";

use SeqTools::Kmer qw(count);
use SeqTools::Hamming qw(distance);
use SeqTools::Levenshtein qw(distance as lev);
use SeqTools::BoyerMoore qw(search);

die <<"USAGE" unless @ARGV;
Usage:
  $0 kmer <sequence> <k>
  $0 hamming <seq1> <seq2>
  $0 levenshtein <seq1> <seq2>
  $0 motif <sequence> <pattern>
USAGE

my ($cmd, @args) = @ARGV;

if ($cmd eq 'kmer') {
    die "Need seq & k\n" unless @args == 2;
    my ($s, $k) = @args;
    my $h = count($s, $k);
    print "$_\t$h->{$_}\n" for sort keys %$h;

} elsif ($cmd eq 'hamming') {
    die "Need two equal-length seqs\n" unless @args == 2;
    print distance(@args), "\n";

} elsif ($cmd eq 'levenshtein') {
    die "Need two seqs\n" unless @args == 2;
    print lev(@args), "\n";

} elsif ($cmd eq 'motif') {
    die "Need sequence & pattern\n" unless @args == 2;
    my $hits = search(@args);
    if (@$hits) {
        print "Matches at positions: ", join(", ", @$hits), "\n";
    } else {
        print "No matches found\n";
    }

} else {
    die "Unknown command '$cmd'\n";
}
