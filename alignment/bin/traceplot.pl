#!/usr/bin/env perl
use strict;
use warnings;

my ($tsv, $alnA_fa, $alnB_fa, $svg) = @ARGV;
die <<"USAGE" unless $svg;
Usage: $0 MATRIX.tsv ALIGN_A.fa ALIGN_B.fa OUTPUT.svg
USAGE

# 1) Read matrix
open my $fh, '<', $tsv or die "Cannot open $tsv: $!";
chomp( my $hdr = <$fh> );
my @cols = split /\t/, $hdr; shift @cols;
my (@rows, @mat);
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

# 2) Find pos/neg maxima
my ($max_pos, $max_neg) = (0,0);
for my $ri (0..$R-1) {
  for my $ci (0..$C-1) {
    my $v = defined $mat[$ri][$ci] ? $mat[$ri][$ci] : 0;
    $max_pos = $v       if $v >  $max_pos;
    $max_neg = -$v      if $v < -$max_neg;
  }
}
die "No non-zero values\n" unless $max_pos || $max_neg;

# 3) Read alignments and compute traceback path
open my $faA, '<', $alnA_fa or die $!;
<$faA>; chomp(my $alnA = <$faA>); close $faA;
open my $faB, '<', $alnB_fa or die $!;
<$faB>; chomp(my $alnB = <$faB>); close $faB;
die "Alignment lengths differ\n"
  unless length($alnA) == length($alnB);

my @trace;
my ($ti, $tj) = (0,0);
for my $k (0 .. length($alnA)-1) {
  my $ca = substr($alnA, $k, 1);
  my $cb = substr($alnB, $k, 1);
  if    ($ca ne '-' && $cb ne '-') { $ti++; $tj++ }
  elsif ($ca eq '-' && $cb ne '-') {      $tj++ }
  elsif ($ca ne '-' && $cb eq '-') { $ti++      }
  else                              { next }  # both gaps
  push @trace, [ $ti-1, $tj-1 ];
}

# 4) Layout
my $cell = 20;
my $rmax = int($cell * 0.4);
my $ml   = 100;
my $mt   = 100;
my $w    = $ml + $C*$cell + 20;
my $h    = $mt + $R*$cell + 20;

open my $out, '>', $svg or die "Cannot open $svg: $!";
print $out <<"SVG_HDR";
<?xml version="1.0" encoding="UTF-8"?>
<svg width="$w" height="$h" xmlns="http://www.w3.org/2000/svg" version="1.1">
 <rect width="100%" height="100%" fill="white"/>
 <style> text { font-family: monospace; font-size:12px } </style>
 <defs>
   <marker id="arrow" markerWidth="6" markerHeight="6"
           refX="5" refY="3" orient="auto">
     <path d="M0,0 L6,3 L0,6 Z" fill="green"/>
   </marker>
 </defs>
SVG_HDR

# 5) Column labels (rotated)
for my $ci (0..$C-1) {
  my $x   = $ml + $ci*$cell + $cell/2;
  my $y   = $mt - 5;
  my $clb = $cols[$ci];
  print $out qq{ <text x="$x" y="$y" text-anchor="middle" }
           .  qq{transform="rotate(-45 $x,$y)">}
           .  qq{$clb</text>\n};
}

# 6) Row labels
for my $ri (0..$R-1) {
  my $x   = $ml - 5;
  my $y   = $mt + $ri*$cell + $cell/2 + 4;
  my $rlb = $rows[$ri];
  print $out qq{ <text x="$x" y="$y" text-anchor="end">}
           .  qq{$rlb</text>\n};
}

# 7) Dots
for my $ri (0..$R-1) {
  for my $ci (0..$C-1) {
    my $v = defined $mat[$ri][$ci] ? $mat[$ri][$ci] : 0;
    next if $v == 0;
    my $mag   = $v > 0 ? $max_pos : $max_neg;
    my $r     = int($rmax * abs($v)/$mag + 0.5) || 1;
    my $cx    = $ml + $ci*$cell + $cell/2;
    my $cy    = $mt + $ri*$cell + $cell/2;
    my $col   = $v > 0 ? 'black' : 'red';
    print $out qq{ <circle cx="$cx" cy="$cy" r="$r" fill="$col" />\n};
  }
}

# 8) Traceback polyline
if (@trace > 1) {
  my @pts = map {
    my ($xr,$xc) = @$_;
    my $x = $ml + $xc*$cell + $cell/2;
    my $y = $mt + $xr*$cell + $cell/2;
    "$x,$y"
  } @trace;
  my $pts_str = join(' ', @pts);
  print $out qq{ <polyline points="$pts_str"\n}
           .  qq{            fill="none" stroke="green" stroke-width="2"\n}
           .  qq{            marker-end="url(#arrow)" />\n};
}

print $out "</svg>\n";
close $out;
print "Wrote SVG to $svg ($w×$h px)\n";
