#!/usr/bin/perl -w
use strict;
##fixIDs.pl by Rohan Maddamsetti

open INPUT, "input/091130_VeraCode_A-1_assay.tab"
  or die "cannot open file for input: $!\n";
open OUTPUT, ">input/091130_VeraCode_A-1_assay_rm.tab"
  or die "cannot create file for output:$!\n";

my %letter_hash = ('A' => 1,
		'B' => 2,
		'C' => 3,
		'D' => 4,
		'E' => 5,
		'F' => 6,
		'G' => 7,
		'H' => 8, );

my $count = 0;
foreach my $line (<INPUT>) {
    if ($line =~ /R00\d_C0\d{2}/){
	my $row_letter = substr($line,0,1);
	my $row_replacer = $letter_hash{$row_letter};
	my $column_number = substr($line,1,2);
	$column_number =~ s/\t//;
	$column_number = '0'.$column_number if int($column_number) < 10;
	#print "row_letter, row_replacer, and column_number: ",$row_letter,' ',$row_replacer,' ', $column_number,"\n";
	
	$line =~ s/R00\d_C0\d{2}/R00($row_replacer)_C0($column_number)/;
	$line =~ s/\(|\)//g;
	$count++;
	print OUTPUT $line;
    }
    else { print OUTPUT $line; }
}
print "$count lines replaced.\n";
