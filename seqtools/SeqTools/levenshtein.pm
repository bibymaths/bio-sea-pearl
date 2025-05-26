package SeqTools::Levenshtein;
use strict;
use warnings;
use Carp;
use Exporter 'import';
our @EXPORT_OK = qw(distance);

# Inline C implementation for high-performance Levenshtein distance
use Inline C => <<'END_C';
#include <string.h>
#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"

int
levenshtein_c(const char *s, const char *t) {
    size_t n = strlen(s);
    size_t m = strlen(t);
    /* Ensure n <= m to use less memory */
    if (n > m) {
        const char *tmp = s; s = t; t = tmp;
        size_t tmp_len = n; n = m; m = tmp_len;
    }
    /* Allocate two rows */
    int *prev = (int *)malloc((m+1) * sizeof(int));
    int *curr = (int *)malloc((m+1) * sizeof(int));
    if (!prev || !curr) {
        croak("Memory allocation failed");
    }
    /* Initialize base case: distance from empty prefix */
    for (size_t j = 0; j <= m; j++) {
        prev[j] = j;
    }
    /* Compute DP rows */
    for (size_t i = 1; i <= n; i++) {
        curr[0] = i;
        for (size_t j = 1; j <= m; j++) {
            int cost = (s[i-1] == t[j-1]) ? 0 : 1;
            int ins = curr[j-1] + 1;
            int del = prev[j] + 1;
            int sub = prev[j-1] + cost;
            int min = ins < del ? ins : del;
            curr[j] = (min < sub ? min : sub);
        }
        /* Swap rows for next iteration */
        int *tmp_row = prev; prev = curr; curr = tmp_row;
    }
    int result = prev[m];
    free(prev);
    free(curr);
    return result;
}
END_C

sub distance {
    my ($s, $t) = @_;
    return levenshtein_c($s, $t);
}

1;
