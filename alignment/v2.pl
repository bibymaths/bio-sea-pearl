#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;

##
##  Sequence alignment (global / local / LCS)
##  – rolling two‐row DP (O(n) memory)
##  – full traceback via pointer matrix
##  – streams full score matrix to TSV as it goes
##  – no substr in the inner loop, no locks, no threads
##

#─ command‐line options ─────────────────────────────────────────────────────────
my $mode      = 'global';      # global | local | lcs
my $threads   = 1;             # unused (see note below)
my $out_pref  = 'align';
my $match     =  2;
my $mismatch  = -1;
my $gap       = -2;

GetOptions(
  'mode=s'     => \$mode,
  'threads=i'  => \$threads,
  'out=s'      => \$out_pref,
  'match=i'    => \$match,
  'mismatch=i' => \$mismatch,
  'gap=i'      => \$gap,
) or die "Usage: $0 --mode [global|local|lcs] --out PREFIX seqA.fasta seqB.fasta\n";

@ARGV == 2 or die "Please give two sequence files (FASTA or plain-text).\n";
my ($seqA, $seqB) = map { slurp_seq($_) } @ARGV;

#─ preprocess ──────────────────────────────────────────────────────────────────
my @A = split //, $seqA;
my @B = split //, $seqB;
my $m = @A;
my $n = @B;

#─ open the TSV matrix file and print header + row0 ────────────────────────────
open my $MF, '>', "$out_pref.matrix.tsv"
  or die "Cannot write matrix file: $!";
print $MF "\t", join("\t", @B), "\n";

# initial row 0 of scores
my @H_prev =
    $mode eq 'global'
  ? map { $gap * $_ } 0..$n
  : (0) x ($n+1);
print $MF join("\t", ('') , @H_prev), "\n";

#─ allocate full pointer matrix for traceback ──────────────────────────────────
#   0 = diag, 1 = up, 2 = left, 3 = stop (for local)
my @Trace;
# row 0 pointers
$Trace[0] = [ map { $mode eq 'global' ? 2 : 3 } 0..$n ];
# column 0 pointers for rows 1..m
for my $i (1..$m) {
  $Trace[$i] = [ (0) x ($n+1) ];
  $Trace[$i][0] = $mode eq 'global' ? 1 : 3;
}

#─ the main DP loop (rolling rows) ─────────────────────────────────────────────
my @H_cur = (0) x ($n+1);

for my $i (1..$m) {
  # initialize first column of this row
  $H_cur[0] = $mode eq 'global' ? $i * $gap : 0;

  for my $j (1..$n) {
    # compute the three scores
    my $diag = $H_prev[$j-1] + ($A[$i-1] eq $B[$j-1] ? $match : $mismatch);
    my $up   = $H_prev[$j]   + $gap;
    my $le   = $H_cur[$j-1]  + $gap;

    if ($mode eq 'local') {
      # Smith‐Waterman (zero clamp)
      my $best = 0;
      my $ptr  = 3;
      if ($diag > $best) { $best = $diag; $ptr = 0 }
      if ($up   > $best) { $best = $up;   $ptr = 1 }
      if ($le   > $best) { $best = $le;   $ptr = 2 }
      $H_cur[$j]       = $best;
      $Trace[$i][$j]   = $ptr;
    }
    elsif ($mode eq 'global') {
      # Needleman–Wunsch
      my $best = $diag;
      my $ptr  = 0;
      if ($up > $best) { $best = $up;  $ptr = 1 }
      if ($le > $best) { $best = $le;  $ptr = 2 }
      $H_cur[$j]       = $best;
      $Trace[$i][$j]   = $ptr;
    }
    else {
      # LCS
      if ($A[$i-1] eq $B[$j-1]) {
        $H_cur[$j]     = $H_prev[$j-1] + 1;
        $Trace[$i][$j] = 0;
      } else {
        if ($H_prev[$j] >= $H_cur[$j-1]) {
          $H_cur[$j]     = $H_prev[$j];
          $Trace[$i][$j] = 1;
        } else {
          $H_cur[$j]     = $H_cur[$j-1];
          $Trace[$i][$j] = 2;
        }
      }
    }
  }

  # stream this row to disk
  print $MF join("\t", $A[$i-1], @H_cur), "\n";

  # roll
  @H_prev = @H_cur;
}

close $MF;

#─ traceback ───────────────────────────────────────────────────────────────────
# find start point
my ($si, $sj) = $mode eq 'local'
  ? do {
      my ($bi, $bj, $bv) = (0,0,0);
      for my $i (0..$m) {
        for my $j (0..$n) {
          if ($Trace[$i][$j] != 3 && $H_prev[$j] > $bv) {
            ($bi,$bj,$bv) = ($i,$j,$H_prev[$j]);
          }
        }
      }
      ($bi,$bj);
    }
  : ($m,$n);

my ($Aa, $Ba) = ('','');
my ($i,$j) = ($si,$sj);

while ($i>0 || $j>0) {
  last if $mode eq 'local' && $Trace[$i][$j] == 3;
  my $p = $Trace[$i][$j];
  if    ($p == 0) { $Aa .= $A[$i-1];     $Ba .= $B[$j-1];     $i--; $j-- }
  elsif ($p == 1) { $Aa .= $A[$i-1];     $Ba .= '-';          $i-- }
  elsif ($p == 2) { $Aa .= '-';          $Ba .= $B[$j-1];     $j-- }
  else            { last }
}
$Aa = reverse $Aa;
$Ba = reverse $Ba;

#─ write alignments ────────────────────────────────────────────────────────────
open my $FA, '>', "$out_pref.A.fa" or die $!;
open my $FB, '>', "$out_pref.B.fa" or die $!;
print $FA ">A_$mode\n$Aa\n";
print $FB ">B_$mode\n$Ba\n";
close $FA;
close $FB;

#─ final report ────────────────────────────────────────────────────────────────
my $score = $mode eq 'local'
  ? do { my $mx=0; $mx = $_ > $mx ? $_ : $mx for @H_prev; $mx }
  : $H_prev[$n];

print <<"EOF";
Completed [$mode] alignment.
Score: $score
Alignment → $out_pref.A.fa, $out_pref.B.fa
Full DP matrix → $out_pref.matrix.tsv
EOF

#───────────────────────────────────────────────────────────────────────────────
sub slurp_seq {
  my $fn = shift;
  open my $fh, '<', $fn or die "Cannot open $fn: $!";
  local $/;
  my $t = <$fh>;
  close $fh;
  $t =~ s/>.*?\R//g;    # strip FASTA header lines
  $t =~ s/\s+//g;
  return uc $t;
}