#!/usr/bin/perl -w

## formatMixedPopSequencing.pl by Rohan Maddamsetti.
## This script reformats MixedPopSequencingSNPs.csv so that
## 1) the names of genes are consistent with the names of SNPs in the
## Veracode project, and 2) so that the data is nicely readable in R.

use strict;
use Text::CSV;

##This stuff is boilerplate code from the Text::CSV documentation.
my @rows;
my $csv = Text::CSV->new ( { binary => 1 } ) or
  die "Cannot use CSV: ".Text::CSV->error_diag ();

open my $fh, "<:encoding(cp1252)", "MixedPopSequencingSNPs.csv" or die "WTF: $!";
open my $output_fh, ">:encoding(utf-8)", "EditedMixedPopSequencingSNPs.csv" or die "WTF: $!";

while (my $row = $csv->getline($fh)) {
  push @rows, $row;
}
$csv->eof or $csv->error_diag();
close $fh;

# Now, go through each row of data, and modify the gene names.
my $first_line = 1;
for my $rizzow (@rows) {
    if ($first_line) {
	$csv->print ($output_fh, $rizzow);
	$first_line = 0;
    }
    else {
    my $raw_gene = $rizzow->[7];
    $raw_gene =~ s/\//-/g;
    my $gene = "Ara-1-". $raw_gene . "-snp";
    $rizzow->[7] = $gene;
    $csv->print ($output_fh, $rizzow);
#    print "@$rizzow", "\n";
    }
}


