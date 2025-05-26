package SeqTools::Levenshtein;
use strict;
use warnings;
use Exporter 'import';
our @EXPORT_OK = qw(distance);

sub distance {
    my ($s, $t) = @_;
    # ensure $s is shorter
    ($s, $t) = ($t, $s) if length($s) > length($t);
    my @S = split //, $s;
    my @T = split //, $t;
    my $m = @S;
    my $n = @T;
    my @prev = (0 .. $n);
    my @curr = (0) x ($n + 1);

    for my $i (1 .. $m) {
        $curr[0] = $i;
        my $si = $S[$i - 1];
        for my $j (1 .. $n) {
            my $cost = $si eq $T[$j - 1] ? 0 : 1;
            my $ins = $curr[$j - 1] + 1;
            my $del = $prev[$j]       + 1;
            my $sub = $prev[$j - 1]   + $cost;
            my $min = $ins < $del ? $ins : $del;
            $curr[$j] = $min < $sub ? $min : $sub;
        }
        @prev = @curr;
    }
    return $prev[$n];
}

1;
