#!/usr/bin/perl -w
use strict;

##fixMixedClones.pl by Rohan Maddamsetti
##This fixes the mixed clone SNP frequencies in 091130_VeraCode_A-1_assay.tab,
##which are currently reversed.

open ASSAYFILE, "input/091130_VeraCode_A-1_assay.tab" or
  die "cannot open file for input: $!";

open OUTPUT, ">091130_VeraCode_A-1_assay.tab" or
    die "cannot open file for output: $!";

while(<ASSAYFILE>){
    my @data = split(/\t/);
    print OUTPUT join("\t",@data) unless $data[6] eq 'MC';
    if ($data[6] eq 'MC'){
	my @slice_data =  @data[9..$#data];
	foreach (@slice_data){
	    next if $_ == 0;
	    $_ = 1 - $_; #take complement.
	}      
	my @fixed_data = @data[0..8];
	push(@fixed_data,@slice_data);
	print OUTPUT join("\t",@fixed_data),"\n";
    }
}
