#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use IO::Handle;
STDOUT->autoflush(1);


# ANSI color constants
use constant {
  COL_GREEN => "\e[32m",
  COL_CYAN  => "\e[36m",
  COL_RED   => "\e[31m",
  COL_RESET => "\e[0m",
};

our @DP;
### command‐line options ###
my ($mode, $matfile, $gapopen, $gapext, $out_pref) =
    ('global','',10,1,'align');
GetOptions(
  'mode=s'    => \$mode,      # global | local | lcs
  'matrix=s'  => \$matfile,   # .mat file for scoring
  'gapopen=i' => \$gapopen,
  'gapext=i'  => \$gapext,
  'out=s'     => \$out_pref,
) or die "Usage: $0 --mode [global|local|lcs] --matrix MATRIX_FILE \\\n"
        . "           --gapopen N --gapext M --out PREFIX seqA seqB\n";
@ARGV==2 or die "Need exactly two sequence files\n";
die "--matrix is required\n" unless $matfile;

# Slurp a FASTA or plain‐text sequence
sub slurp_seq {
  my $fn = shift;
  open my $fh, '<', $fn or die "Cannot open $fn: $!";
  local $/; my $t = <$fh>;
  close $fh;
  $t =~ s/>.*?\R//g;
  $t =~ s/\s+//g;
  return uc $t;
}
my ($seqA, $seqB) = map { slurp_seq($_) } @ARGV;
my @A = split //, $seqA;
my @B = split //, $seqB;

# Read a whitespace‐delimited substitution matrix
sub read_matrix {
  my ($file) = @_;
  open my $fh, '<', $file or die "Cannot open matrix '$file': $!";
  my @cols;

  # find header
  while (<$fh>) {
    next if /^\s*(?:\/\*|\*|\/\/|#|;)/;
    chomp;
    my @tok = grep { length } split /\s+/, $_;
    next unless @tok >= 2;
    @cols = @tok;
    last;
  }
  die "No header found in '$file'\n" unless @cols;

  # read rows
  my %M;
  while (<$fh>) {
    next if /^\s*(?:\/\*|\*|\/\/|#|;)/;
    chomp;
    my @f = grep { length } split /\s+/, $_;
    next unless @f == @cols + 1;
    my $row = shift @f;
    $M{$row}{$cols[$_]} = $f[$_] for 0..$#cols;
  }
  close $fh;
  return \%M;
}

my $MAT = read_matrix($matfile);
sub score { $MAT->{ $_[0] }{ $_[1] } // 0 }

# Gotoh DP in linear space: returns \@best_scores for A --> B (or reversed)
sub NWScore {
  my ($Aref,$Bref,$rev) = @_;
  my @Av = @$Aref;  my @Bv = @$Bref;
  @Av = reverse @Av if $rev;
  @Bv = reverse @Bv if $rev;
  my ($m,$n) = (scalar @Av, scalar @Bv);
  my $NEG = -1e9;
  my (@Sprev,@Eprev,@Fprev,@Scur,@Ecur,@Fcur);

  # init row 0
  $Sprev[0] = 0;
  $Eprev[0] = $Fprev[0] = $NEG;
  for my $j (1..$n) {
    if ($mode eq 'global') {
      $Sprev[$j] = -$gapopen - ($j-1)*$gapext;
      $Eprev[$j] = $Sprev[$j] - $gapopen;
      $Fprev[$j] = $NEG;
    } else {
      $Sprev[$j] = 0;
      $Eprev[$j] = $Fprev[$j] = ($mode eq 'lcs' ? $NEG : 0);
    }
  }

  if (!$rev) {
    $DP[0] = [ @Sprev ];   # store scores for i=0
  }

  # fill rows
  for my $i (1..$m) {
    if ($mode eq 'global') {
      $Scur[0] = -$gapopen - ($i-1)*$gapext;
      $Ecur[0] = $NEG;
      $Fcur[0] = $Scur[0] - $gapopen;
    } else {
      $Scur[0] = 0;
      $Ecur[0] = $Fcur[0] = ($mode eq 'lcs' ? $NEG : 0);
    }

    for my $j (1..$n) {
      # match/mismatch or LCS
      my $sub = $mode eq 'lcs'
              ? ($Av[$i-1] eq $Bv[$j-1] ? 1 : 0)
              : score($Av[$i-1], $Bv[$j-1]);

      # best of M/E/F from (i-1,j-1)
      my $bp = $Sprev[$j-1];
      $bp = $Eprev[$j-1] if $Eprev[$j-1] > $bp;
      $bp = $Fprev[$j-1] if $Fprev[$j-1] > $bp;
      my $mval = $bp + $sub;
      $mval = 0 if $mode eq 'local' && $mval < 0;

      # E (gap in A; horiz)
      my $eopen = $Scur[$j-1] - $gapopen;
      my $eext  = $Ecur[$j-1] - $gapext;
      my $eval  = $eopen >= $eext ? $eopen : $eext;
      $eval = 0 if $mode eq 'local' && $eval < 0;

      # F (gap in B; vert)
      my $fopen = $Sprev[$j]   - $gapopen;
      my $fext  = $Fprev[$j]   - $gapext;
      my $fval  = $fopen >= $fext ? $fopen : $fext;
      $fval = 0 if $mode eq 'local' && $fval < 0;

      # best-of-three
      my $sv = $mval;
      $sv = $eval if $eval > $sv;
      $sv = $fval if $fval > $sv;

      $Scur[$j] = $sv;
      $Ecur[$j] = $eval;
      $Fcur[$j] = $fval;
    }
    # pointer matrix
    if (!$rev) {
      $DP[$i] = [ @Scur ];   # store the full row i
    }
    progress_bar($i, $m);
    # rotate in
    @Sprev = @Scur;
    @Eprev = @Ecur;
    @Fprev = @Fcur;
  }
  return [ @Sprev ];
}

### full-DP + traceback when one dimension == 1 ###
sub small_align {
  my ($Aref,$Bref) = @_;
  my @Av = @$Aref;  my @Bv = @$Bref;
  my ($m,$n) = (scalar @Av, scalar @Bv);
  my $NEG = -1e9;
  my (@M,@E,@F,@P);

  # initialize
  for my $i (0..$m) {
    for my $j (0..$n) {
      $M[$i][$j] = $E[$i][$j] = $F[$i][$j] = 0;
      $P[$i][$j] = '';
    }
  }
  $E[0][0] = $F[0][0] = $NEG;
  $P[0][0] = 'X';

  # edges
  for my $j (1..$n) {
    if ($mode eq 'global') {
      $M[0][$j] = $NEG;
      $E[0][$j] = -$gapopen - ($j-1)*$gapext;
      $F[0][$j] = $NEG;
      $P[0][$j] = 'E';
    } else {
      $M[0][$j] = 0;
      $E[0][$j] = $F[0][$j] = ($mode eq 'lcs' ? $NEG : 0);
      $P[0][$j] = 'X';
    }
  }
  for my $i (1..$m) {
    if ($mode eq 'global') {
      $M[$i][0] = $NEG;
      $E[$i][0] = $NEG;
      $F[$i][0] = -$gapopen - ($i-1)*$gapext;
      $P[$i][0] = 'F';
    } else {
      $M[$i][0] = 0;
      $E[$i][0] = $F[$i][0] = ($mode eq 'lcs' ? $NEG : 0);
      $P[$i][0] = 'X';
    }
  }

  # fill
  for my $i (1..$m) {
    for my $j (1..$n) {
      # M
      my $sub = $mode eq 'lcs'
              ? ($Av[$i-1] eq $Bv[$j-1] ? 1:0)
              : score($Av[$i-1], $Bv[$j-1]);
      my $sp = $M[$i-1][$j-1];
      $sp = $E[$i-1][$j-1] if $E[$i-1][$j-1] > $sp;
      $sp = $F[$i-1][$j-1] if $F[$i-1][$j-1] > $sp;
      my $mval = $sp + $sub;
      if ($mode eq 'local' && $mval < 0) { $mval = 0 }
      $M[$i][$j] = $mval;

      # E
      my $eopen = $M[$i][$j-1] - $gapopen;
      my $eext  = $E[$i][$j-1] - $gapext;
      my $eval  = $eopen >= $eext ? $eopen : $eext;
      if ($mode eq 'local' && $eval < 0) { $eval = 0 }
      $E[$i][$j] = $eval;

      # F
      my $fopen = $M[$i-1][$j] - $gapopen;
      my $fext  = $F[$i-1][$j] - $gapext;
      my $fval  = $fopen >= $fext ? $fopen : $fext;
      if ($mode eq 'local' && $fval < 0) { $fval = 0 }
      $F[$i][$j] = $fval;

      # pointer
      my $best = $M[$i][$j];
      my $mat  = 'M';
      if ($E[$i][$j] > $best) { $best = $E[$i][$j]; $mat = 'E' }
      if ($F[$i][$j] > $best) { $best = $F[$i][$j]; $mat = 'F' }
      $P[$i][$j] = $mat;
    }
   progress_bar($i, $m);
  }

  # endpoint
  my ($ei,$ej) = ($m,$n);
  if ($mode eq 'local') {
    my $mx = 0;
    ($ei,$ej) = (0,0);
    for my $i (0..$m) {
      for my $j (0..$n) {
        if ($M[$i][$j] > $mx) {
          $mx = $M[$i][$j];
          ($ei,$ej) = ($i,$j);
        }
      }
    }
  }

  # traceback
  my ($i,$j)=($ei,$ej);
  my ($Aout,$Bout)=('','');
  while ($i>0 || $j>0) {
    last if $mode eq 'local' && $P[$i][$j] eq 'X';
    if    ($P[$i][$j] eq 'M') { $Aout.= $Av[$i-1]; $Bout.= $Bv[$j-1]; $i--, $j-- }
    elsif ($P[$i][$j] eq 'E') { $Aout.= '-';       $Bout.= $Bv[$j-1];        $j--    }
    elsif ($P[$i][$j] eq 'F') { $Aout.= $Av[$i-1]; $Bout.= '-';               $i--    }
    else                       { last }
  }
  return (reverse $Aout, reverse $Bout);
}

### Hirschberg divide‐and‐conquer ###
sub hirschberg {
  my ($Aref,$Bref) = @_;
  my $m = @$Aref;
  my $n = @$Bref;
  return ('-' x $n, join('',@$Bref))        if $m==0;
  return (join('',@$Aref), '-' x $m)        if $n==0;
  return small_align($Aref,$Bref)           if $m==1 || $n==1;

  my $i_mid = int($m/2);
  my $L = NWScore([ @$Aref[0..$i_mid-1] ], $Bref, 0);
  my $R = NWScore([ @$Aref[$i_mid..$m-1] ], $Bref, 1);

  my ($best,$j_mid) = (-1e9,0);
  for my $j (0..$n) {
    my $s = $L->[$j] + $R->[$n-$j];
    if ($s > $best) { $best=$s; $j_mid=$j }
  }

  my ($A1,$B1) = hirschberg(
    [ @$Aref[0..$i_mid-1] ], [ @$Bref[0..$j_mid-1] ]
  );
  my ($A2,$B2) = hirschberg(
    [ @$Aref[$i_mid..$m-1] ], [ @$Bref[$j_mid..$n-1] ]
  );

  return ($A1.$A2, $B1.$B2);
}

# Build the full DP matrix for global/local/lcs alignment
sub build_full_DP {
    my ($Aref, $Bref) = @_;
    my @Av = @$Aref;
    my @Bv = @$Bref;
    my ($m, $n) = (scalar(@Av), scalar(@Bv));
    my $NEG = -1e9;

    # rows of best/E/F
    my (@Sprev, @Eprev, @Fprev, @Scur, @Ecur, @Fcur);
    # full DP storage
    my @DP;

    # init row 0
    $Sprev[0] = 0;
    $Eprev[0] = $Fprev[0] = $NEG;
    for my $j (1..$n) {
        if ($mode eq 'global') {
            $Sprev[$j] = -$gapopen - ($j-1)*$gapext;
            $Eprev[$j] = $Sprev[$j] - $gapopen;
            $Fprev[$j] = $NEG;
        } else {
            $Sprev[$j] = 0;
            $Eprev[$j] = $Fprev[$j] = ($mode eq 'lcs' ? $NEG : 0);
        }
    }
    # store row 0
    $DP[0] = [ @Sprev ];

    # fill rows 1..m
    for my $i (1..$m) {
        if ($mode eq 'global') {
            $Scur[0] = -$gapopen - ($i-1)*$gapext;
            $Ecur[0] = $NEG;
            $Fcur[0] = $Scur[0] - $gapopen;
        } else {
            $Scur[0] = 0;
            $Ecur[0] = $Fcur[0] = ($mode eq 'lcs' ? $NEG : 0);
        }

        for my $j (1..$n) {
            # substitution or LCS
            my $sub = $mode eq 'lcs'
                    ? ($Av[$i-1] eq $Bv[$j-1] ? 1 : 0)
                    : score($Av[$i-1], $Bv[$j-1]);

            # M from (i-1,j-1)
            my $bp = $Sprev[$j-1];
            $bp = $Eprev[$j-1] if $Eprev[$j-1] > $bp;
            $bp = $Fprev[$j-1] if $Fprev[$j-1] > $bp;
            my $mval = $bp + $sub;
            $mval = 0 if $mode eq 'local' && $mval < 0;

            # E (gap in A; horiz)
            my $eopen = $Scur[$j-1] - $gapopen;
            my $eext  = $Ecur[$j-1] - $gapext;
            my $eval  = $eopen >= $eext ? $eopen : $eext;
            $eval = 0 if $mode eq 'local' && $eval < 0;

            # F (gap in B; vert)
            my $fopen = $Sprev[$j]   - $gapopen;
            my $fext  = $Fprev[$j]   - $gapext;
            my $fval  = $fopen >= $fext ? $fopen : $fext;
            $fval = 0 if $mode eq 'local' && $fval < 0;

            # best-of-three
            my $sv = $mval;
            $sv = $eval if $eval > $sv;
            $sv = $fval if $fval > $sv;

            $Scur[$j] = $sv;
            $Ecur[$j] = $eval;
            $Fcur[$j] = $fval;
        }

        # store row i
        $DP[$i] = [ @Scur ];

        # rotate
        @Sprev = @Scur;
        @Eprev = @Ecur;
        @Fprev = @Fcur;
    }

    return \@DP;
}

# Longest Common Subsequence (LCS) extraction
sub extract_lcs {
  my ($a, $b) = @_;
  my @runs;
  my $cur = "";

  for my $i (0 .. length($a)-1) {
    my $ca = substr($a, $i, 1);
    my $cb = substr($b, $i, 1);

    if ($ca eq $cb and $ca ne '-') {
      # still matching, accumulate
      $cur .= $ca;
    }
    else {
      # break in the run
      push @runs, $cur if length $cur;
      $cur = "";
    }
  }
  push @runs, $cur if length $cur;   # final run
  return @runs;
}

# Write the DP matrix to a binary file
sub write_DP_binary {
    my ($file, $DP_ref) = @_;
    open my $out, '>:raw', $file
      or die "Cannot write $file: $!";
    my $nrows = scalar @$DP_ref;
    my $ncols = scalar @{ $DP_ref->[0] };
    # header
    print $out pack('N2', $nrows, $ncols);
    # each row
    for my $row (@$DP_ref) {
      die "DP row length mismatch\n" unless @$row == $ncols;
      print $out pack('l>*', @$row);
    }
    close $out;
    warn "Wrote DP binary: $file ($nrows×$ncols cells)\n";
}

# Print the alignment in a colored format
sub print_alignment {
  my ($seqA, $seqB) = @_;
  my $len = length $seqA;
  my ($posA, $posB) = (0,0);

  for (my $off = 0; $off < $len; $off += 80) {
    # extract up to 80‐column block
    my $blockA = substr($seqA, $off, 80);
    my $blockB = substr($seqB, $off, 80);
    my $blen   = length $blockA;

    # line for Seq A
    printf "%6d ", $posA+1;
    for (my $i=0; $i < $blen; $i++) {
      my $a = substr($blockA,$i,1);
      my $b = substr($blockB,$i,1);
      # choose color
      my $col = ($a eq '-' || $b eq '-') ? COL_RED
              : ($a eq $b)              ? COL_GREEN
                                        : COL_CYAN;
      print $col, $a, COL_RESET;
      $posA++ if $a ne '-';
    }
    print " ", $posA, "\n";

    # line for Seq B
    printf "%6d ", $posB+1;
    for (my $i=0; $i < $blen; $i++) {
      my $a = substr($blockA,$i,1);
      my $b = substr($blockB,$i,1);
      my $col = ($a eq '-' || $b eq '-') ? COL_RED
              : ($a eq $b)              ? COL_GREEN
                                        : COL_CYAN;
      print $col, $b, COL_RESET;
      $posB++ if $b ne '-';
    }
    print " ", $posB, "\n\n";
  }
}

# Progress bar for long operations
sub progress_bar {
  my ($done, $total, $width) = @_;
  $width ||= 40;
  my $pct    = $total ? $done/$total : 0;
  my $filled = int($pct * $width);
  my $bar    = "[" . ("#" x $filled) . (" " x ($width - $filled)) . "]";
  printf "\r%s %3.0f%%", $bar, $pct*100;
  if ($done >= $total) {
    print "\n";
  }
}

# Main execution starts here
my ($AlnA, $AlnB);

if ( $mode eq 'lcs' ) {
  # LCS: we want the full gotoh on the whole matrix
  ($AlnA, $AlnB) = small_align(\@A, \@B);
}
else {
  # both global and local go through Hirschberg --> which will
  # route down to small_align only when one seq has length 1
  ($AlnA, $AlnB) = hirschberg(\@A, \@B);
}

print "\n\n";
print_alignment($AlnA, $AlnB);

open my $FA, '>', "$out_pref.A.fa" or die $!;
print $FA ">A_$mode\n$AlnA\n";
close $FA;

open my $FB, '>', "$out_pref.B.fa" or die $!;
print $FB ">B_$mode\n$AlnB\n";
close $FB;

# build the full matrix (O(mn) memory!)
my $DP_full = build_full_DP(\@A, \@B);

# open TSV
open my $MF, '>', "$out_pref.matrix.tsv"
  or die "Cannot write $out_pref.matrix.tsv: $!";
# header row: a leading empty cell, then B[0],B[1],…B[n-1]
print $MF "\t", join("\t", @B), "\n";
# each DP_full->[$i] is an arrayref of length n+1
for my $i (0 .. $#{ $DP_full }) {
  # row‐label: empty for i=0, otherwise A[i-1]
  my $label = $i ? $A[$i-1] : '';
  # join the scores
  my @row = @{ $DP_full->[$i] };
  print $MF $label, "\t", join("\t", @row), "\n";
}
close $MF;
warn "Wrote full DP matrix to $out_pref.matrix.tsv\n";

if ($mode eq 'lcs') {
  my @runs = extract_lcs($AlnA, $AlnB);
  my ($longest) = sort { length $b <=> length $a } @runs;
  open my $LC, '>', "$out_pref.lcs.fa" or die $!;
  print $LC ">${out_pref}_LCS\n";
  print $LC "$longest\n";
  close $LC;
  warn "Wrote LCS runs to $out_pref.lcs.fa\n";
}

write_DP_binary("$out_pref.matrix.bin", $DP_full);

# compute final score correctly
my $final_score;
if ( $mode eq 'local' ) {
  # local: maximum over all cells
  $final_score = -1e9;
  for my $i (0..$#$DP_full) {
    for my $j (0..$#{ $DP_full->[0] }) {
      $final_score = $DP_full->[$i][$j]
        if $DP_full->[$i][$j] > $final_score;
    }
  }
} else {
  # global (or lcs): bottom‐right corner
  my $m = scalar(@A);
  my $n = scalar(@B);
  $final_score = $DP_full->[$m][$n];
}

print <<"EOF";
Score   :  $final_score
Output  :  $out_pref.A.fa, $out_pref.B.fa
Completed $mode alignment.
EOF
