package SeqTools::BoyerMoore;
use strict;
use warnings;
use Exporter 'import';
our @EXPORT_OK = qw(search);

sub search {
    my ($text, $pat) = @_;
    my @T = split //, $text;
    my @P = split //, $pat;
    my $n = @T;
    my $m = @P;
    return [] if $m == 0 || $n < $m;

    # bad-char table
    my @bad = (-1) x 256;
    $bad[ ord($P[$_]) ] = $_ for 0..$m-1;

    my @matches;
    my $i = 0;
    while ($i <= $n - $m) {
        my $j = $m - 1;
        while ($j >= 0 && $P[$j] eq $T[$i + $j]) {
            $j--;
        }
        if ($j < 0) {
            push @matches, $i;
            $i++;
        } else {
            my $lo = $bad[ ord $T[$i + $j] ];
            my $shift = $j - $lo;
            $i += $shift > 0 ? $shift : 1;
        }
    }
    return \@matches;
}

1;
