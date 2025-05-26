#!/usr/bin/env perl
use strict;
use warnings;
no warnings 'recursion';      # silence deep‐recursion warnings

use threads;
use threads::shared;
use Getopt::Long;

### command‐line options ###
my ($mode, $out_pref, $matfile, $gapopen, $gapext, $threads) =
    ('global','align','',10,1,4);

GetOptions(
  'mode=s'     => \$mode,      # global | local | lcs
  'out=s'      => \$out_pref,
  'matrix=s'   => \$matfile,   # substitution matrix file
  'gapopen=i'  => \$gapopen,
  'gapext=i'   => \$gapext,
  'threads=i'  => \$threads,
) or die "Usage: $0 --mode [global|local|lcs] --matrix MATRIX_FILE \\\n"
       . "          --gapopen N --gapext M --threads T --out PREFIX seqA seqB\n";

@ARGV == 2 or die "Please provide two sequence files (FASTA or plain‐text) and --matrix\n";

#── helper: slurp and clean a sequence file ────────────────────────────────────
sub slurp_seq {
  my $fn = shift;
  open my $fh, '<', $fn or die "Cannot open sequence $fn: $!";
  local $/; my $t = <$fh>;
  close $fh;
  $t =~ s/>.*?\R//g;     # strip FASTA headers
  $t =~ s/\s+//g;
  return uc $t;
}

#── helper: read a whitespace‐delimited matrix with header row/col ─────────────
sub read_matrix {
  my ($file) = @_;
  open my $fh, '<', $file or die "Cannot open matrix file '$file': $!";

  # 1) find the header row (skip any C‐ or shell‐style comments)
  my @cols;
  while (<$fh>) {
    next if /^\s*(?:\/\*|\*|\/\/|#|;)/;   # skip comment lines
    chomp;
    my @tok = grep { length } split /\s+/, $_;
    next if @tok < 2;                     # not our header
    @cols = @tok;
    last;
  }
  die "No header row found in '$file'\n" unless @cols;

  # 2) read only rows with exactly (1 + #cols) tokens
  my %M;
  while (<$fh>) {
    next if /^\s*(?:\/\*|\*|\/\/|#|;)/;
    chomp;
    my @f = grep { length } split /\s+/, $_;
    next unless @f == @cols + 1;
    my $row = shift @f;
    for my $i (0..$#cols) {
      $M{$row}{$cols[$i]} = $f[$i];
    }
  }
  close $fh;
  return \%M;
}


#── dynamic scoring lookup ─────────────────────────────────────────────────────
die "--matrix is required\n" unless $matfile;
my $MAT = read_matrix($matfile);
sub score {
  my ($a,$b) = @_;
  return $MAT->{$a}{$b} // 0;
}

#── read sequences ─────────────────────────────────────────────────────────────
my ($seqA, $seqB) = map { slurp_seq($_) } @ARGV;
my @A = split //, $seqA;
my @B = split //, $seqB;
my $m = @A;
my $n = @B;

#── allocate shared DP & trace matrices: M, E, F + TM, TE, TF ────────────────
my (@SM, @SE, @SF, @TM, @TE, @TF) : shared;

sub shared_row {
  my @r : shared = map { undef } 0..$n;
  return \@r;
}

for my $i (0..$m) {
  $SM[$i] = shared_row();
  $SE[$i] = shared_row();
  $SF[$i] = shared_row();
  $TM[$i] = shared_row();
  $TE[$i] = shared_row();
  $TF[$i] = shared_row();
}

#── initialize base cases ─────────────────────────────────────────────────────
my $NEG_INF = 0;
if ($mode eq 'global') {
  # (0,0)
  $SM[0][0] = 0;
  $SE[0][0] = $SF[0][0] = $NEG_INF;
  $TM[0][0] = $TE[0][0] = $TF[0][0] = 3;
  # row 0
  for my $j (1..$n) {
    $SM[0][$j] = $NEG_INF;
    $SE[0][$j] = -$gapopen - ($j-1)*$gapext;
    $SF[0][$j] = $NEG_INF;
    $TM[0][$j] = 3;
    $TE[0][$j] = $j==1 ? 0 : 1;   # open vs extend
    $TF[0][$j] = 3;
  }
  # col 0
  for my $i (1..$m) {
    $SM[$i][0] = $NEG_INF;
    $SE[$i][0] = $NEG_INF;
    $SF[$i][0] = -$gapopen - ($i-1)*$gapext;
    $TM[$i][0] = 3;
    $TE[$i][0] = 3;
    $TF[$i][0] = $i==1 ? 0 : 1;
  }
}
else {
  # local or LCS: zero‐clamp all
  for my $i (0..$m) {
    for my $j (0..$n) {
      $SM[$i][$j] = 0;
      $SE[$i][$j] = 0 if $mode eq 'local';  # clamp E/F too for local
      $SF[$i][$j] = 0 if $mode eq 'local';
      $SE[$i][$j] = $SF[$i][$j] = $NEG_INF if $mode eq 'lcs';
      $TM[$i][$j] = $TE[$i][$j] = $TF[$i][$j] = 3;
    }
  }
}

#── open matrix TSV, print header + row0 ──────────────────────────────────────
open my $MF, '>', "$out_pref.matrix.tsv" or die $!;
print $MF "\t", join("\t", @B), "\n";
{
  my @r0 = map {
    my $v = $SM[0][$_];
    $v = 0 if $mode eq 'local' && $v < 0;
    $v
  } 0..$n;
  print $MF join("\t", '', @r0), "\n";
}

#── compute one diagonal’s cells in parallel ──────────────────────────────────
sub compute_cells {
  my ($cells) = @_;
  foreach my $ij (@$cells) {
    my ($i,$j) = @$ij;
    # 1) M‐matrix
    my $m1 = $SM[$i-1][$j-1];
    my $e1 = $SE[$i-1][$j-1];
    my $f1 = $SF[$i-1][$j-1];
    my $prev_best = $m1 >= $e1
                  ? ($m1 >= $f1 ? $m1 : $f1)
                  : ($e1 >= $f1 ? $e1 : $f1);
    my $sc = $mode eq 'lcs'
           ? ($A[$i-1] eq $B[$j-1] ? 1 : 0)
           : score($A[$i-1], $B[$j-1]);
    my $mval = $prev_best + $sc;
    my $pm   = $prev_best == $m1 ? 0 : ($prev_best == $e1 ? 1 : 2);
    if ($mode eq 'local' && $mval < 0) { $mval = 0; $pm = 3 }
    $SM[$i][$j] = $mval;
    $TM[$i][$j] = $pm;

    # 2) E‐matrix (gap in A; horizontal)
    my $openE = $SM[$i][$j-1] - $gapopen;
    my $extE  = $SE[$i][$j-1] - $gapext;
    my ($eval, $pe) = $openE >= $extE ? ($openE,0) : ($extE,1);
    if ($mode eq 'local' && $eval < 0) { $eval = 0; $pe = 3 }
    $SE[$i][$j] = $eval;
    $TE[$i][$j] = $pe;

    # 3) F‐matrix (gap in B; vertical)
    my $openF = $SM[$i-1][$j] - $gapopen;
    my $extF  = $SF[$i-1][$j] - $gapext;
    my ($fval, $pf) = $openF >= $extF ? ($openF,0) : ($extF,1);
    if ($mode eq 'local' && $fval < 0) { $fval = 0; $pf = 3 }
    $SF[$i][$j] = $fval;
    $TF[$i][$j] = $pf;
  }
}

#── recursive diagonal fill ───────────────────────────────────────────────────
sub fill_diag {
  my ($d) = @_;
  return if $d > $m + $n;
  my @cells;
  for my $i (1..$m) {
    my $j = $d - $i;
    next unless $j >= 1 && $j <= $n;
    push @cells, [$i,$j];
  }
  my $chunk = int(@cells / $threads) || 1;
  my @workers;
  while (@cells) {
    my @sl = splice(@cells, 0, $chunk);
    push @workers, threads->create(\&compute_cells, \@sl);
  }
  $_->join for @workers;
  fill_diag($d+1);
}

fill_diag(2);

#── stream rows 1..m to the TSV (best of M/E/F) ──────────────────────────────
for my $i (1..$m) {
  my @row = map {
    my $v = $SM[$i][$_];
    $v = 0 if $mode eq 'local' && $v < 0;
    $v
  } 0..$n;
  print $MF join("\t", $A[$i-1], @row), "\n";
}
close $MF;

#── locate best endpoint & starting matrix ───────────────────────────────────
my ($bi,$bj,$cur) = (0,0,'M');
if ($mode eq 'local') {
  my $best = 0;
  for my $i (0..$m) {
    for my $j (0..$n) {
      for my $mat ('M','E','F') {
        my $v = $mat eq 'M' ? $SM[$i][$j]
              : $mat eq 'E' ? $SE[$i][$j]
                             : $SF[$i][$j];
        if ($v > $best) {
          ($best,$bi,$bj,$cur) = ($v,$i,$j,$mat);
        }
      }
    }
  }
} else {
  # global or LCS: end at (m,n), pick best matrix
  my $vM = $SM[$m][$n];
  my $vE = $SE[$m][$n];
  my $vF = $SF[$m][$n];
  if    ($vM >= $vE && $vM >= $vF) { $cur='M' }
  elsif ($vE >= $vM && $vE >= $vF) { $cur='E' }
  else                              { $cur='F' }
  ($bi,$bj) = ($m,$n);
}

#── build the pointer path (single‐thread) ────────────────────────────────────
my @path;
{
  my ($i,$j,$mat) = ($bi,$bj,$cur);
  while ($i>0 || $j>0) {
    my $p;
    if ($mat eq 'M') {
      $p = $TM[$i][$j];
      last if $mode eq 'local' && $p == 3;
      push @path, [$i,$j,'M'];
      $i--; $j--;
      $mat = $p==0 ? 'M' : $p==1 ? 'E' : 'F';
    }
    elsif ($mat eq 'E') {
      $p = $TE[$i][$j];
      last if $mode eq 'local' && $p == 3;
      push @path, [$i,$j,'E'];
      $j--;
      $mat = $p==0 ? 'M' : 'E';
    }
    else {  # 'F'
      $p = $TF[$i][$j];
      last if $mode eq 'local' && $p == 3;
      push @path, [$i,$j,'F'];
      $i--;
      $mat = $p==0 ? 'M' : 'F';
    }
  }
}
@path = reverse @path;

#── parallel traceback to build aligned strings ───────────────────────────────
my $plen  = @path;
my $chunk = int($plen / $threads) || 1;
my @thr;
for my $t (0..$threads-1) {
  my $start = $t * $chunk;
  my $end   = $t == $threads-1 ? $plen-1 : $start + $chunk-1;
  last if $start > $end;
  push @thr, threads->create(sub {
    my ($as, $bs) = ('','');
    for my $k ($start..$end) {
      my ($ii,$jj,$mm) = @{ $path[$k] };
      if    ($mm eq 'M') { $as .= $A[$ii-1]; $bs .= $B[$jj-1] }
      elsif ($mm eq 'E') { $as .= '-';        $bs .= $B[$jj-1] }
      else               { $as .= $A[$ii-1]; $bs .= '-'        }
    }
    return ($as, $bs);
  });
}

my ($Aa, $Ba) = ('','');
for my $thr (@thr) {
  my ($pa,$pb) = $thr->join;
  $Aa .= $pa;
  $Ba .= $pb;
}

#── write FASTA alignments ───────────────────────────────────────────────────
open my $FA, '>', "$out_pref.A.fa" or die $!;
open my $FB, '>', "$out_pref.B.fa" or die $!;
print $FA ">A_$mode\n$Aa\n";
print $FB ">B_$mode\n$Ba\n";
close $FA; close $FB;

#── final report ─────────────────────────────────────────────────────────────
my $final_score = ($cur eq 'M' ? $SM[$bi][$bj]
                 : $cur eq 'E' ? $SE[$bi][$bj]
                                : $SF[$bi][$bj]);

print <<"EOT";
Completed [$mode] affine‐gap alignment.
Matrix loaded from: $matfile
Gap-open:   $gapopen
Gap-extend: $gapext
Score:      $final_score

Alignment → $out_pref.A.fa, $out_pref.B.fa
DP matrix  → $out_pref.matrix.tsv
EOT
