package SeqTools::Fasta;
use strict;
use warnings;
use Exporter 'import';

our @EXPORT_OK = qw(read_fasta);

# Read the first sequence record from a FASTA file and return its sequence string.
sub read_fasta {
    my ($file) = @_;
    open my $fh, '<', $file
        or die "Cannot open FASTA file '$file': $!";
    local $/ = "\n>";
    my $rec = <$fh>;
    chomp $rec;
    $rec =~ s/^>//;
    my (undef, @lines) = split /\n/, $rec;
    return join('', @lines);
}

1;
