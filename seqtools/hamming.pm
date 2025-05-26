package SeqTools::Hamming;
use strict;
use warnings;
use Exporter 'import';
our @EXPORT_OK = qw(distance);

sub distance {
    my ($a, $b) = @_;
    die "Lengths differ\n" unless length($a) == length($b);
    # count bytes where XOR =! NUL
    return ($a ^ $b) =~ tr/\0//c;
}

1;
