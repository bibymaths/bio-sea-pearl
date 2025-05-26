#!/usr/bin/env perl
use strict;
use warnings;

# usage:
#   perl dotplot.pl matrix.tsv dotplot.ppm [cell_scale] [max_dot_radius]
# example:
#   perl dotplot.pl myalign.matrix.tsv dotplot.ppm 4 10
#
# cell_scale: pixels per matrix‐cell (default 4)
# max_dot_radius: max circle radius in pixels for the largest score (default 10)

my ($tsv, $out, $scale, $rmax) = @ARGV;
die "Usage: $0 MATRIX.tsv OUT.ppm [cell_scale] [max_dot]\n"
  unless defined $tsv && defined $out;

$scale ||= 4;
$rmax  ||= 10;

# 1) Read the matrix and strip any header row/column
open my $fh, '<', $tsv or die "Cannot open $tsv: $!";
chomp(my $hdr = <$fh>);
# assume header row of labels; ignore it
my @mat;
while (<$fh>) {
  chomp;
  my @f = split /\t/;
  shift @f;               # drop first column (row label)
  push @mat, \@f;
}
close $fh;
my $rows = @mat;
my $cols = @{ $mat[0] // [] };
die "Empty or malformed matrix\n" unless $rows && $cols;

# 2) find the maximum positive score
my $max = 0;
for my $r (0..$rows-1) {
  for my $c (0..$cols-1) {
    my $v = $mat[$r][$c];
       $max = $v if $v > $max;
  }
}
die "No positive scores to plot\n" if $max <= 0;

# 3) prepare a white canvas of size (cols*scale × rows*scale)
my $W = $cols * $scale;
my $H = $rows * $scale;
my @img = map { [ (255,255,255) x $W ] } 0..$H-1;

# 4) draw a black circle for each positive‐score cell
for my $r (0..$rows-1) {
  for my $c (0..$cols-1) {
    my $v = $mat[$r][$c];
    next if $v <= 0;
    # radius scaled by value/max
    my $rad = int( $rmax * ($v/$max) + 0.5 );
    next unless $rad > 0;
    # center of this cell
    my $cx = int($c*$scale + $scale/2);
    my $cy = int($r*$scale + $scale/2);
    # draw filled circle via simple bounding‐box scan
    my $r2 = $rad*$rad;
    for my $dy (-$rad .. $rad) {
      my $y = $cy + $dy;
      next if $y < 0 || $y >= $H;
      my $dxmax = int( sqrt($r2 - $dy*$dy) );
      for my $dx (-$dxmax .. $dxmax) {
        my $x = $cx + $dx;
        next if $x < 0 || $x >= $W;
        # set pixel to black
        $img[$y][$x*3  ] = 0;
        $img[$y][$x*3+1] = 0;
        $img[$y][$x*3+2] = 0;
      }
    }
  }
}

# 5) write out a binary PPM (P6)
open my $outfh, '>', $out or die "Cannot create $out: $!";
binmode $outfh;
print $outfh "P6\n$W $H\n255\n";
for my $y (0..$H-1) {
  print $outfh pack("C*", @{ $img[$y] });
}
close $outfh;

print "Wrote dot-plot to $out ($W×$H pixels)\n";
