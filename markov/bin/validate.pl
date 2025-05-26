#!/usr/bin/env perl
use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin/../lib";
use Getopt::Long;
use MarkovChainHigherOrder;

# CLI options
my ($fasta, $order, $pseudocount, $help);
$order      = 1;
$pseudocount = 0;
GetOptions(
    'fasta=s'       => \$fasta,
    'order=i'       => \$order,
    'pseudocount=f' => \$pseudocount,
    'help|h'        => \$help,
) or usage();
usage() if $help || !$fasta;

sub usage {
    die <<"USAGE";
Usage: $0 --fasta FILE [--order K] [--pseudocount FLOAT]

  --fasta        Input FASTA file
  --order        Markov order K (integer ≥1; default: 1)
  --pseudocount  Smoothing pseudocount for each observed transition (default: 0)
  --help         Show this message
USAGE
}

# Read sequences & build counts
my $seqs = MarkovChainHigherOrder->read_fasta($fasta);
my ($T, $ord) = MarkovChainHigherOrder->build_transitions($seqs, $order);

# Compute smoothed transition probabilities P[state][next]
my %P;
for my $state (keys %$T) {
    # apply pseudocount only to observed transitions
    my %W = map { $_ => $T->{$state}{$_} + $pseudocount } keys %{ $T->{$state} };
    my $sum = 0; $sum += $_ for values %W;
    die "State '$state' has no outgoing weight\n" if $sum == 0;
    $P{$state}{$_} = $W{$_} / $sum for keys %W;
}

# Score held-in data
my $total_ll      = 0;
my $n_transitions = 0;
for my $seq (@$seqs) {
    my @a = split //, uc $seq;
    next if @a <= $ord;
    for my $i (0 .. $#a - $ord - 1) {
        my $state = join('', @a[$i .. $i+$ord-1]);
        my $next  = $a[$i+$ord];
        my $p = $P{$state}{$next} // 1e-12;  # floor for unseen
        $total_ll      += log($p);
        $n_transitions++;
    }
}

# Metrics
my $avg_ll     = $total_ll / $n_transitions;
my $perplexity = exp( - $avg_ll );

printf "Total transitions: %d\n", $n_transitions;
printf "Average log-likelihood per transition: %.4f nats\n", $avg_ll;
printf "Perplexity: %.4f\n", $perplexity;
