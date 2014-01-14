#!/usr/bin/perl -w

##cleanMuller.pl by Rohan Maddamsetti
##This script does some simple post-processing of the file created
##in CalibrateAndGraphSNPs_final.R in the makeMullerInput function
##that exists in GenotypingLibrary.R.
##Precisely, this script removes the column names in the header; it makes all
##values greater than 1 (i.e. 1.01 values) into 1s; and it removes prefixes and
##suffixes from the SNP Names, such as 'Ara-1-' and '-snp.'
##This also adds an 'anc' genotype line to the head of the file.

##usage: perl cleanMuller.pl > mullerinput.dat

use strict;

open ROUGHINPUT, "rough_mullerinput.dat" or die "cannot open file: $!";

# Make the ancestor line.
my @ancpoints = split(//, "1" x 62);
unshift @ancpoints, 'anc';
my $anc_genotype = join("\t", @ancpoints);
print $anc_genotype, "\n";

while (<ROUGHINPUT>){
# ignore the header.    
    next if m/^\d/;
# globally replace prefixes and suffixes with nothing.
    s/Ara-1-//g;
    s/-snp//g;
    s/1.01/1/g;
    s/1.02/1/g;
# If an allele frequency of 1 is reached, change the rest of the allele frequencies to 1.
    my @line = split(/\t/, $_);
    for my $counter (0..$#line){
	my $mark = 0; # This to make sure that the rest of the line is formatted only once.
	if ($line[$counter] eq "1" and $mark == 0) {
	    my @end = @line[$counter .. $#line];
	    my @rest = split(//, "1" x ($#end+1));
	    @line[$counter .. $#line] = @rest;
	    $mark = 1;
	    #print "length of array slices", scalar(@end), " ", scalar(@rest), "\n";
	}
	last if $mark;

##do a similar procedure, where all frequencies before the last zero value are set to zero.
	for my $counter2 (reverse (0..$#line)) {
	    my $mark2 = 0; #To make sure that the rest of the line is formatted only once.
#	    print $counter2, "\n";
	    if ($line[$counter2] eq "0" and $mark == 0){
		my @beginning = @line[1..$counter2];
		my @rest2 = split(//, "0" x ($counter2));
		@line[1..$counter2] = @rest2;
		$mark = 1;
		#print "length of beginning array: ", scalar(@beginning), "\n";
		#print "length of rest2 array: ", scalar(@rest2), "\n";
	    }
	}
    }
    my $formattedline = join("\t", @line);
# add the ancestral background to all genotype names.
    print "anc|". $formattedline."\n";
}
