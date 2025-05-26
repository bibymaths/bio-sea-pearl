#!/usr/bin/perl
use strict;
use warnings;

package MarkovChainHigherOrder;
our $VERSION = '1.02';

# Constructor:
#   - $trans: hashref of { state => { next1 => count, … } }
#   - %opts:
#       method      => 'alias' or 'binsrch' (default: alias)
#       pseudocount => numeric ≥0 (default: 0)
#       order       => integer ≥1 (default: 1)
sub new {
    my ($class, $trans, %opts) = @_;
    my $method = $opts{method}      || 'alias';
    my $pseudo = $opts{pseudocount} || 0;
    my $order  = $opts{order}       || 1;

    die "Unknown method '$method'\n"
        unless $method =~ /^(?:alias|binsrch)$/;
    die "Order must be ≥1\n" unless $order >= 1;

    my $self = { method => $method, pseudocount => $pseudo, order => $order, chain => {} };

    while (my ($state, $succ) = each %$trans) {
        # apply pseudocount
        my %W = map { $_ => ($succ->{$_}||0) + $pseudo } keys %$succ;
        my $total = 0; $total += $_ for values %W;
        die "State '$state' has no outgoing mass\n" if $total == 0;

        if ($method eq 'alias') {
            $self->{chain}{$state} = _make_alias(\%W, $total);
        }
        else {
            my @keys = sort keys %W;
            my @cum; my $acc = 0;
            for my $k (@keys) {
                $acc += $W{$k} / $total;
                push @cum, [$k, $acc];
            }
            $self->{chain}{$state} = \@cum;
        }
    }

    return bless $self, $class;
}

# Build transitions from sequences for arbitrary order
# Returns \\%T and the order
# Usage: my ($T, $ord) = MarkovChain->build_transitions([@seqs], $order);
sub build_transitions {
    my ($class, $seqs, $order) = @_;
    $order ||= 1;
    die "Order must be ≥1\n" unless $order >=1;

    my %T;
    for my $seq (@$seqs) {
        my @a = split //, uc $seq;
        next if @a <= $order;
        for my $i (0 .. $#a - $order) {
            my $state = join('', @a[$i .. $i + $order - 1]);
            my $next  = $a[$i + $order];
            $T{$state}{$next}++;
        }
    }
    return (\%T, $order);
}

# Read FASTA file into arrayref of sequences
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
    close $fh;
    return \@seqs;
}

# Draw next symbol for a given state
sub next_state {
    my ($self, $state) = @_;
    my $spec = $self->{chain}{$state}
        or die "Unknown state '$state'\n";
    return $self->{method} eq 'alias'
         ? _draw_alias($spec)
         : _draw_binsrch($spec);
}

# Generate a walk of N transitions from a start state (length = start + N)
sub generate {
    my ($self, $N, $start) = @_;
    my $k = $self->{order};
    die "Start state must be length $k\n" unless length($start) == $k;

    my @out = (split //, $start);
    my $cur = $start;
    for (1..$N) {
        my $next = $self->next_state($cur);
        push @out, $next;
        $cur = substr($cur, 1) . $next;
    }
    return join('', @out);
}

### Internal: Alias method
sub _make_alias {
    my ($W, $total) = @_;
    my @keys = keys %$W;
    my $n = @keys;
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
    $_ = 1 for @P;
    return { keys => \@keys, prob => \@P, alias => \@alias };
}

sub _draw_alias {
    my ($spec) = @_;
    my $n = @{ $spec->{keys} };
    my $i = int rand($n);
    return $spec->{keys}[$i] if rand() < $spec->{prob}[$i];
    return $spec->{keys}[ $spec->{alias}[$i] ];
}

### Internal: Binary search on cumulative
sub _draw_binsrch {
    my ($cum) = @_;
    my $u = rand;
    my ($lo, $hi) = (0, $#$cum);
    while ($lo < $hi) {
        my $mid = int(($lo + $hi) / 2);
        $u <= $cum->[$mid][1] ? ($hi = $mid) : ($lo = $mid + 1);
    }
    return $cum->[$lo][0];
}

1;
__END__

=pod

=head1 NAME

MarkovChain – zero-dependency, order-k Markov-chain simulator

=head1 SYNOPSIS

    use MarkovChain;
    my $seqs = MarkovChain->read_fasta('data.fa');
    my ($T, $ord) = MarkovChain->build_transitions($seqs, 2);
    my $mc = MarkovChain->new(
        $T,
        method      => 'alias',
        pseudocount => 1,
        order       => $ord,
    );
    my $walk = $mc->generate(1000, 'AT');
    print "$walk\n";

=head1 DESCRIPTION

Supports arbitrary-order (k) Markov chains, alias or binary-search sampling, and pseudocount smoothing.

=cut
