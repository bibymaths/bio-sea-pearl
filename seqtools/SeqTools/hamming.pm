package SeqTools::Hamming;
use strict;
use warnings;
use Carp qw(croak);

# Inline C for high-performance byte-wise comparison
# to implement Hamming distance in C-level loop.
use Inline C => <<'END_C';
#include <string.h>
#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"

int
hamming_c(const char *a, const char *b) {
    size_t len = strlen(a);
    if (len != strlen(b)) {
        croak("Lengths differ: %zu vs %zu", len, strlen(b));
    }
    size_t count = 0;
    for (size_t i = 0; i < len; i++) {
        if (a[i] != b[i]) {
            count++;
        }
    }
    return (int)count;
}
END_C

use Exporter 'import';
our @EXPORT_OK = qw(distance);

sub distance {
    my ($a, $b) = @_;
    # C implementation
    return hamming_c($a, $b);
}

1;
