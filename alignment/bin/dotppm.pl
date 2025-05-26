#!/usr/bin/env perl
use strict;
use warnings;

my ($tsv, $svg) = @ARGV;
die "Usage: $0 MATRIX.tsv OUTPUT.svg\n" unless $svg;

# 1) Read TSV, grab col labels and row labels+values
open my $fh, '<', $tsv or die $!;
chomp(my $hdr = <$fh>);
my @cols = split /\t/, $hdr;
# drop the leading empty cell
shift @cols;

my @rows;
my @mat;
while (<$fh>) {
  chomp;
  my @f = split /\t/;
  my $rlabel = shift @f;
  push @rows, $rlabel;
  push @mat, [@f];
}
close $fh;

my $R = @rows;
my $C = @cols;
die "Empty matrix\n" unless $R && $C;

# 2) find magnitude maxima for positive & negative scaling
my ($max_pos, $max_neg) = (0, 0);
for my $i (0..$R-1) {
  for my $j (0..$C-1) {
    my $v = defined $mat[$i][$j] ? $mat[$i][$j] : 0;
    if    ($v >  $max_pos) { $max_pos =  $v }
    elsif ($v < -$max_neg) { $max_neg = -$v }
  }
}
die "No non-zero values to plot\n" if $max_pos==0 && $max_neg==0;


# 3) layout parameters
my $cell = 20;                 # pixels per cell
my $rmax = int($cell * 0.4);   # max circle radius
my $ml   = 100;                # left margin for row labels
my $mt   = 100;                # top margin for col labels
my $w    = $ml + $C * $cell + 20;
my $h    = $mt + $R * $cell + 20;

open my $out, '>', $svg or die $!;
print $out <<"SVG_HDR";
<?xml version="1.0" encoding="UTF-8"?>
<svg width="$w" height="$h" xmlns="http://www.w3.org/2000/svg" version="1.1">
 <rect width="100%" height="100%" fill="white"/>
 <style> text { font-family: monospace; font-size:12px; } </style>
SVG_HDR

# 4) column labels (rotated)
for my $j (0..$C-1) {
  my $x = $ml + $j*$cell + $cell/2;
  my $y = $mt - 5;
  my $lbl = $cols[$j];
  print $out qq{ <text x="$x" y="$y" text-anchor="middle" transform="rotate(-45 $x,$y)">},
              $lbl,qq{</text>\n};
}

# 5) row labels
for my $i (0..$R-1) {
  my $x = $ml - 5;
  my $y = $mt + $i*$cell + $cell/2 + 4;
  my $lbl = $rows[$i];
  print $out qq{ <text x="$x" y="$y" text-anchor="end">},$lbl,qq{</text>\n};
}

# 6) dots for each non-zero matrix cell
for my $i (0..$R-1) {
  for my $j (0..$C-1) {
    my $v = defined $mat[$i][$j] ? $mat[$i][$j] : 0;
    next if $v == 0;

    # choose scale based on sign
    my $mag   = $v > 0 ? $max_pos : $max_neg;
    my $r     = int( $rmax * (abs($v) / $mag) + 0.5 ) || 1;
    my $cx    = $ml + $j*$cell + $cell/2;
    my $cy    = $mt + $i*$cell + $cell/2;
    my $color = $v > 0 ? 'black' : 'red';

    print $out qq{ <circle cx="$cx" cy="$cy" r="$r" fill="$color" />\n};
  }
}


print $out "</svg>\n";
close $out;
print "Wrote SVG to $svg ($w×$h px)\n";
