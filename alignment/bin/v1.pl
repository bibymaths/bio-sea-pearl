#!/usr/bin/env perl
use strict;
use warnings;
use threads;
use threads::shared;
use Getopt::Long;

### parameters ###
my $match     =  2;
my $mismatch  = -1;
my $gap       = -2;
my $threads   =  4;
my $mode      = 'global';  # global | local | lcs
my $out_prefix= 'align';

GetOptions(
    'mode=s'   => \$mode,
    'threads=i'=> \$threads,
    'out=s'    => \$out_prefix,
    'match=i'  => \$match,
    'mismatch=i'=>\$mismatch,
    'gap=i'    => \$gap,
) or die "Usage: $0 --mode [global|local|lcs] --threads N --out PREFIX seqA.fasta seqB.fasta\n";

@ARGV==2 or die "Provide two sequence files or plain-text sequence files\n";
my ($seqA, $seqB) = map { slurp_seq($_) } @ARGV;
my $lenA = length $seqA;
my $lenB = length $seqB;

## allocate DP and trace matrices (shared for threads)
my @H :shared;       # score matrix
my @Trace :shared;   # pointer matrix: 0=diag,1=up,2=left,3=zero(for local)

# initialize row 0 and col 0
sub init_matrix {
    $H[0] = shared_clone([ map { $_ } 0 .. $lenB ]);
    $Trace[0] = shared_clone([ map { 3 } 0 .. $lenB ]);
    for my $i (1..$lenA) {
        $H[$i] = shared_clone([ 0 x ($lenB+1) ]);
        $Trace[$i] = shared_clone([ 0 x ($lenB+1) ]);
        if ($mode eq 'global' or $mode eq 'lcs') {
            $H[$i][0] = $mode eq 'lcs' ? 0 : $i * $gap;
            $Trace[$i][0] = 1;
        }
    }
}

# worker to fill a block of rows
sub worker_fill {
    my ($i_start, $i_end) = @_;
    for my $i ($i_start .. $i_end) {
        my @row = (0) x ($lenB+1);
        my @tptr= (0) x ($lenB+1);
        my $Ai = substr($seqA, $i-1, 1);
        for my $j (1..$lenB) {
            # SIMD trick: pack 3 candidates into binary (one cell at a time here)
            my $diag = $H[$i-1][$j-1] + ( $Ai eq substr($seqB,$j-1,1) ? $match : $mismatch );
            my $up   = $H[$i-1][$j] + $gap;
            my $left = $H[$i][$j-1] + $gap;
            if ($mode eq 'local') {
                my $best = $diag;
                my $ptr  = 0;
                if ($up > $best)    { $best = $up;    $ptr=1 }
                if ($left> $best)   { $best = $left;  $ptr=2 }
                if ($best < 0)      { $best = 0;      $ptr=3 }
                $row[$j] = $best;
                $tptr[$j]= $ptr;
            }
            elsif ($mode eq 'global') {
                my $best = $diag;
                my $ptr  = 0;
                if ($up > $best)    { $best = $up;    $ptr=1 }
                if ($left> $best)   { $best = $left;  $ptr=2 }
                $row[$j]   = $best;
                $tptr[$j]  = $ptr;
            }
            else {  # LCS
                if ($Ai eq substr($seqB,$j-1,1)) {
                    $row[$j] = $H[$i-1][$j-1] + 1;
                    $tptr[$j]= 0;
                } else {
                    if ($H[$i-1][$j] >= $H[$i][$j-1]) {
                        $row[$j] = $H[$i-1][$j];
                        $tptr[$j]= 1;
                    } else {
                        $row[$j] = $H[$i][$j-1];
                        $tptr[$j]= 2;
                    }
                }
            }
        }
        { lock(@H);     @{$H[$i]}     = @row; }
        { lock(@Trace); @{$Trace[$i]} = @tptr; }
    }
}

# run fills in parallel
sub fill_matrix {
    init_matrix();
    my $block = int($lenA / $threads) || 1;
    my @thr;
    for my $t (0..$threads-1) {
        my $i0 = 1 + $t*$block;
        last if $i0 > $lenA;
        my $i1 = ($t == $threads-1) ? $lenA : $i0 + $block -1;
        push @thr, threads->create(\&worker_fill, $i0, $i1);
    }
    $_->join for @thr;
}

# traceback from (i,j)
sub traceback {
    my ($i, $j) = @_;
    my ($Aout, $Bout) = ('','');
    while ($i>0 or $j>0) {
        my $ptr = $Trace[$i][$j];
        last if $mode eq 'local' && $ptr==3;
        if ($ptr==0) {
            $Aout .= substr($seqA,$i-1,1);
            $Bout .= substr($seqB,$j-1,1);
            $i--; $j--;
        }
        elsif ($ptr==1) {
            $Aout .= substr($seqA,$i-1,1);
            $Bout .= '-';
            $i--;
        }
        elsif ($ptr==2) {
            $Aout .= '-';
            $Bout .= substr($seqB,$j-1,1);
            $j--;
        }
    }
    return (reverse $Aout, reverse $Bout);
}

# save DP matrix
sub save_matrix {
    open my $fh, '>', "$out_prefix.matrix.tsv" or die $!;
    print $fh join("\t", "", split //, $seqB), "\n";
    for my $i (0..$lenA) {
        my $row = join("\t", ($i ? substr($seqA,$i-1,1) : ''), @{$H[$i]});
        print $fh "$row\n";
    }
    close $fh;
}

# save alignment
sub save_alignment {
    my ($Aal, $Bal) = @_;
    open my $fa, '>', "$out_prefix.A.fa" or die $!;
    open my $fb, '>', "$out_prefix.B.fa" or die $!;
    print $fa ">A\_$mode\n$Aal\n";
    print $fb ">B\_$mode\n$Bal\n";
    close $fa; close $fb;
}

### main ###
fill_matrix();

# decide starting point
my ($si,$sj) = ( $mode eq 'local'
                  ? do { my $best=0; my($bi,$bj)=(0,0);
                         for my $i(0..$lenA){ for my $j(0..$lenB){
                             if($H[$i][$j]>$best){ $best=$H[$i][$j]; $bi=$i;$bj=$j }
                         }}; ($bi,$bj) }
                  : ($lenA,$lenB)
               );

my ($Aal, $Bal) = traceback($si,$sj);
save_matrix();
save_alignment($Aal, $Bal);

print "Mode: $mode\n";
print "Score: $H[$si][$sj]\n";
print "Alignment written to $out_prefix.A.fa and $out_prefix.B.fa\n";
print "DP matrix to $out_prefix.matrix.tsv\n";

###############
### helpers ###
###############
sub slurp_seq {
    my $f=shift;
    open my $fh, '<', $f or die $!;
    local $/;
    my $t = <$fh>;
    close $fh;
    # strip FASTA header
    $t =~ s/>.*?\n//g;
    $t =~ s/\s+//g;
    return uc $t;
}

__END__
