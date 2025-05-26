#!/usr/bin/perl
use strict;
use warnings;

package MarkovChainFirstOrder;
our $VERSION = '1.01';

# Constructor:
#   - $trans: hashref of { state => { next1 => count, … } }
#   - %opts:
#       method       => 'alias' (default) or 'binsrch'
#       pseudocount  => numeric ≥0 (default 0)
sub new {
    my ($class, $trans, %opts) = @_;
    my $method      = $opts{method}      || 'alias';
    my $pseudo      = $opts{pseudocount} || 0;
    die "Unknown method '$method'\n"
        unless $method =~ /^(?:alias|binsrch)$/;

    my $self = { method => $method, chain => {} };

    while (my ($state, $succ) = each %$trans) {
        # apply pseudocount
        my %W = map { $_ => ($succ->{$_}||0) + $pseudo } keys %$succ;
        my $total = 0; $total += $_ for values %W;
        die "State '$state' has no outgoing mass\n" if $total == 0;

        if ($method eq 'alias') {
            # build alias table
            $self->{chain}{$state} = _make_alias(\%W, $total);
        }
        else {
            # binsrch: build sorted pairs + cumulative
            my @pairs = sort keys %W;
            my @cum; my $acc = 0;
            foreach my $n (@pairs) {
                $acc += $W{$n} / $total;
                push @cum, [$n, $acc];
            }
            $self->{chain}{$state} = \@cum;
        }
    }

    return bless $self, $class;
}

# Build transitions automatically from arrayref of sequences
sub build_transitions {
    my ($class, $seqs) = @_;
    my %T;
    for my $seq (@$seqs) {
        my @a = split //, uc $seq;
        for my $i (0 .. $#a-1) {
            $T{ $a[$i] }{ $a[$i+1] }++;
        }
    }
    return \%T;
}

# Read all sequences from a FASTA file; return arrayref
sub read_fasta {
    my ($class, $file) = @_;
    open my $fh, '<', $file or die "Can't open '$file': $!\n";
    my @seqs; my $acc = '';
    while (<$fh>) {
        chomp;
        if (/^>/) {
            push @seqs, $acc if length $acc;
            $acc = '';
        } else {
            $acc .= $_;
        }
    }
    push @seqs, $acc if length $acc;
    return \@seqs;
}

# Draw next state
sub next_state {
    my ($self, $state) = @_;
    my $spec = $self->{chain}{$state}
        or die "Unknown state '$state'\n";
    return $self->{method} eq 'alias'
         ? _draw_alias($spec)
         : _draw_binsrch($spec);
}

# Generate a walk of N transitions from $start
sub generate {
    my ($self, $N, $start) = @_;
    my @walk = ($start);
    my $cur   = $start;
    for (1..$N) {
        $cur = $self->next_state($cur);
        push @walk, $cur;
    }
    # if single‐char states, return string
    return join('', @walk) if length($start)==1;
    return \@walk;
}

########################################
#— internal: alias‐method builder & draw
sub _make_alias {
    my ($W, $total) = @_;
    my @keys = keys %$W;
    my $n    = @keys;
    # scaled probs
    my @P = map { $W->{$_} / $total * $n } @keys;
    my (@small, @large);
    for my $i (0..$#P) {
        push @{ $P[$i] < 1 ? \@small : \@large }, $i;
    }
    my @alias = (0) x $n;
    while (@small && @large) {
        my $l = pop @small;
        my $g = pop @large;
        $alias[$l] = $g;
        $P[$g]    = $P[$g] - (1 - $P[$l]);
        push @{ $P[$g] < 1 ? \@small : \@large }, $g;
    }
    # all remaining have P ≈ 1
    $_ = 1 for @P;
    return { keys => \@keys, prob => \@P, alias => \@alias };
}

sub _draw_alias {
    my ($spec) = @_;
    my $n = @{ $spec->{keys} };
    my $i = int rand($n);
    # accept or alias
    return $spec->{keys}[ $i ]
        if rand() < $spec->{prob}[$i];
    return $spec->{keys}[ $spec->{alias}[$i] ];
}

########################################
#— internal: binary‐search over cumulative
sub _draw_binsrch {
    my ($cum) = @_;
    my $u = rand;
    # manual binary search
    my ($lo, $hi) = (0, $#$cum);
    while ($lo < $hi) {
        my $mid = int(($lo + $hi) / 2);
        if ($u <= $cum->[$mid][1]) {
            $hi = $mid;
        } else {
            $lo = $mid + 1;
        }
    }
    return $cum->[$lo][0];
}

1;
__END__

=pod

=head1 NAME

MarkovChain – zero‐dependency, high‐precision Markov‐chain simulator

=head1 SYNOPSIS

    use MarkovChain;
    # build from FASTA:
    my $dna_seqs = MarkovChain->read_fasta('genome.fa');
    my $T        = MarkovChain->build_transitions($dna_seqs);
    my $mc       = MarkovChain->new($T, method=>'alias', pseudocount=>1);
    my $walk     = $mc->generate(1000, 'A');

=head1 DESCRIPTION

This module supports both DNA (4 states) and protein (20 states) chains,
with optional pseudocount smoothing, and two very efficient sampling methods.

=over

=item * alias (O(1) per draw)
=item * binary‐search on cumulative (O(log σ) per draw)

=back

=head1 AUTHOR

Abhinav Mishra <mishraabhinav36@gmail.com>

=cut
