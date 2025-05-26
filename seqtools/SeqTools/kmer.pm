package SeqTools::Kmer;
use strict;
use warnings;
use Exporter 'import';
our @EXPORT_OK = qw(count);

sub count {
    my ($s, $k) = @_;
    my $len = length($s);
    die "k ($k) must be ≤ sequence length ($len)\n" if $k > $len;
    my %counts;
    my $end = $len - $k;
    for my $i (0 .. $end) {
        $counts{ substr($s, $i, $k) }++;
    }
    return \%counts;
}

1;
