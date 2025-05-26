use strict;
use warnings;
use GD;

sub plot_dp_heatmap {
    my ($tsv_file, $png_file, $scale) = @_;
    $scale ||= 2;  # pixels per cell

    # 1) read matrix (skip first header row/col)
    open my $fh, '<', $tsv_file or die $!;
    my $hdr = <$fh>;
    my @mat;
    while (<$fh>) {
        chomp;
        my @F = split /\t/;
        shift @F;           # row‐label
        push @mat, \@F;
    }
    close $fh;
    my $rows = @mat;
    my $cols = @{$mat[0]};

    # 2) find min/max for color scale
    my ($min,$max) = ( 1e9, -1e9 );
    for my $r (0..$rows-1) {
      for my $c (0..$cols-1) {
        my $v = $mat[$r][$c];
        $min = $v if $v < $min;
        $max = $v if $v > $max;
      }
    }
    my $range = $max - $min || 1;

    # 3) create image
    my $img = GD::Image->new($cols*$scale, $rows*$scale);
    my $white = $img->colorAllocate(255,255,255);

    # 4) draw each cell as a scaled‐grayscale rectangle
    for my $r (0..$rows-1) {
      for my $c (0..$cols-1) {
        my $v = $mat[$r][$c];
        my $gray = int(255 * ( ($v - $min)/$range ));
        # invert so high scores are dark
        $gray = 255 - $gray;
        my $color = $img->colorAllocate($gray,$gray,$gray);
        my $x1 = $c*$scale;
        my $y1 = $r*$scale;
        $img->filledRectangle($x1,$y1, $x1+$scale-1, $y1+$scale-1, $color);
      }
    }

    # 5) save PNG
    open my $out, '>', $png_file or die $!;
    binmode $out;
    print $out $img->png;
    close $out;
}

sub plot_dp_last_row {
    my ($tsv_file, $png_file, $height) = @_;
    $height ||= 200;

    # 1) read only last row of matrix
    open my $fh, '<', $tsv_file or die $!;
    my $hdr = <$fh>;
    my $last;
    while (<$fh>) {
        chomp;
        $last = [ split /\t/ ] if eof $fh;
    }
    close $fh;
    shift @$last;  # remove row‐label
    my @scores = @$last;
    my $n = @scores;

    # 2) find min/max
    my ($min,$max) = ( 1e9, -1e9 );
    $min = $_ if $_ < $min for @scores;
    $max = $_ if $_ > $max for @scores;
    my $range = $max - $min || 1;

    # 3) create image
    my $img = GD::Image->new($n, $height);
    my $white = $img->colorAllocate(255,255,255);
    my $black = $img->colorAllocate(  0,  0,  0);

    # 4) plot polyline
    my @points;
    for my $i (0..$n-1) {
      my $v = $scores[$i];
      my $y = $height - 1 - int( ($v - $min)/$range * ($height-1) );
      push @points, ($i, $y);
    }
    $img->line(@points, $black);

    # 5) save PNG
    open my $out, '>', $png_file or die $!;
    binmode $out;
    print $out $img->png;
    close $out;
}

plot_dp_heatmap("myalign.matrix.tsv", "myalign.heatmap.png", 3);
plot_dp_last_row("myalign.matrix.tsv", "myalign.lastrow.png", 250);
