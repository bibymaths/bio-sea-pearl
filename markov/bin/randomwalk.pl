#!/usr/bin/env perl
use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin/../lib";
use Getopt::Long;
use MarkovChainHigherOrder;

# Default parameters
my $fasta;
my $length      = 1000;
my $start;
my $method      = 'alias';
my $pseudocount = 0;
my $order       = 1;
my $help;

# Command-line options
GetOptions(
    'fasta=s'       => \$fasta,
    'length=i'      => \$length,
    'start=s'       => \$start,
    'method=s'      => \$method,
    'pseudocount=f' => \$pseudocount,
    'order=i'       => \$order,
    'help|h'        => \$help,
) or usage();

sub usage {
    die <<"USAGE";
Usage: $0 --fasta FILE --length N --start X \
          [--method alias|binsrch] [--pseudocount FLOAT] [--order K]

Options:
  --fasta        Path to FASTA file containing sequences (DNA or protein)
  --length       Number of transitions to generate (default: 1000)
  --start        Starting k-mer of length K (default: length 1)
  --method       Sampling method: 'alias' or 'binsrch' (default: alias)
  --pseudocount  Pseudocount to add to counts (default: 0)
  --order        Markov model order K (integer ≥1, default: 1)
  --help         Show this help message
USAGE
}

# Validate inputs
usage() unless defined $fasta && defined $start;
die "Order must be ≥1\n" unless $order >= 1;
die "Start state must be length $order\n"
    unless length($start) == $order;

# Read sequences from FASTA
my $seqs = MarkovChainHigherOrder->read_fasta($fasta);

# Build transition counts for order = K
my ($T, $ord) = MarkovChainHigherOrder->build_transitions($seqs, $order);

# Create MarkovChain object with specified order
my $mc = MarkovChainHigherOrder->new(
    $T,
    method      => $method,
    pseudocount => $pseudocount,
    order       => $ord,
);

# Generate the random walk
my $walk = $mc->generate($length, $start);

# Output the result
print "$walk\n";