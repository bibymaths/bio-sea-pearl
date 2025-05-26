#!/usr/bin/env perl
use strict;
use warnings;
no warnings 'recursion';    # disable “Deep recursion on subroutine” :contentReference[oaicite:0]{index=0}

# Optional: if you’d rather keep warnings but only every N levels, you can
# install the following helper (from PerlMonks) to raise the threshold:
#
#   use B       'svref_2object';
#   use Symbol 'qualify_to_ref';
#   sub change_recursion_limit {
#     my ($subname, $limit) = @_;
#     my $subref = \&$subname;
#     my $gv     = svref_2object($subref)->GV;
#     no warnings 'redefine';
#     *{ qualify_to_ref $subname } = sub {
#       if ( $gv->CV->DEPTH % $limit == 0 ) {
#         # emit warning only every $limit calls
#         warn sprintf
#           "Deep recursion on '%s' at %s line %d\n",
#           join('::',$gv->STASH->NAME,$gv->NAME), $0, (caller)[2];
#       }
#       goto &$subref;
#     };
#   }
#   # e.g. to warn only every 1000 levels:
#   change_recursion_limit('comp', 1000); :contentReference[oaicite:1]{index=1}

use Getopt::Long;

### options ###
my $mode     = 'global';    # global | local | lcs
my $out_pref = 'align';
my $match    =  2;
my $mismatch = -1;
my $gap      = -2;

GetOptions(
  'mode=s'     => \$mode,
  'out=s'      => \$out_pref,
  'match=i'    => \$match,
  'mismatch=i' => \$mismatch,
  'gap=i'      => \$gap,
) or die "Usage: $0 --mode [global|local|lcs] --out PREFIX seqA.fa seqB.fa\n";

@ARGV==2 or die "Provide two sequence files.\n";
my ($seqA, $seqB) = map { slurp_seq($_) } @ARGV;
my @A = split //, $seqA;
my @B = split //, $seqB;
my $m = @A;
my $n = @B;

# memo tables
my (@Score, @Trace);
my ($best_score, $best_i, $best_j) = (undef,0,0);

sub comp {
  my ($i,$j) = @_;
  return $Score[$i][$j] if defined $Score[$i][$j];

  my ($res,$ptr);
  if ($i==0 || $j==0) {
    if ($mode eq 'global') {
      ($res,$ptr) = $i==0 && $j==0 ? (0,3)
                   : $i==0         ? ($j*$gap,2)
                                   : ($i*$gap,1);
    } else {
      ($res,$ptr) = (0,3);
    }
  }
  else {
    my $d = comp($i-1,$j-1) + ($A[$i-1] eq $B[$j-1] ? $match : $mismatch);
    my $u = comp($i-1,$j)   + $gap;
    my $l = comp($i  ,$j-1) + $gap;

    if ($mode eq 'local') {
      ($res,$ptr) = (0,3);
      ($res,$ptr) = ($d,0) if $d > $res;
      ($res,$ptr) = ($u,1) if $u > $res;
      ($res,$ptr) = ($l,2) if $l > $res;
      if (!defined $best_score || $res > $best_score) {
        ($best_score,$best_i,$best_j) = ($res,$i,$j);
      }
    }
    elsif ($mode eq 'global') {
      ($res,$ptr) = ($d,0);
      ($res,$ptr) = ($u,1) if $u > $res;
      ($res,$ptr) = ($l,2) if $l > $res;
    }
    else {  # LCS
      if ($A[$i-1] eq $B[$j-1]) {
        $res = comp($i-1,$j-1) + 1; $ptr = 0;
      } else {
        if ($u >= $l) { $res=$u; $ptr=1 } else { $res=$l; $ptr=2 }
      }
    }
  }

  $Score[$i][$j] = $res;
  $Trace[$i][$j] = $ptr;
  return $res;
}

# run
if ($mode eq 'local') {
  comp($m,$n);                # best_* set in comp()
} else {
  $best_score = comp($m,$n);
  ($best_i,$best_j) = ($m,$n);
}

# dump matrix
open my $MF, '>', "$out_pref.matrix.tsv" or die $!;
print $MF "\t", join("\t",@B), "\n";
for my $i (0..$m) {
  my $r = defined $A[$i-1] ? $A[$i-1]: '';
  print $MF $r, "\t", join("\t", map {$Score[$i][$_]//0} 0..$n), "\n";
}
close $MF;

# traceback
my ($i,$j) = ($best_i,$best_j);
my ($Aa,$Ba) = ('','');
while ($i>0 || $j>0) {
  last if $mode eq 'local' && $Trace[$i][$j]==3;
  my $p = $Trace[$i][$j];
  if    ($p==0){ $Aa.= $A[$i-1];  $Ba.=$B[$j-1];  $i--; $j-- }
  elsif ($p==1){ $Aa.= $A[$i-1];  $Ba.='-';     $i-- }
  elsif ($p==2){ $Aa.='-';       $Ba.=$B[$j-1]; $j-- }
  else         { last }
}
$Aa=reverse $Aa; $Ba=reverse $Ba;

# write FASTA
open my $FA, '>', "$out_pref.A.fa" or die $!;
open my $FB, '>', "$out_pref.B.fa" or die $!;
print $FA ">A_$mode\n$Aa\n";
print $FB ">B_$mode\n$Ba\n";
close $FA; close $FB;

print <<"EOF";
Completed [$mode] alignment.
Score: $best_score
Alignment → $out_pref.A.fa, $out_pref.B.fa
Full DP matrix → $out_pref.matrix.tsv
EOF

sub slurp_seq {
  my $fn=shift; open my $fh,'<',$fn or die $!;
  local $/; my $t=<$fh>; close $fh;
  $t =~ s/>.*?\R//g; $t =~ s/\s+//g;
  return uc $t;
}
