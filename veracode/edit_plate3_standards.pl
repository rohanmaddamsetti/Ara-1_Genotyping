#!/usr/bin/perl -w

# edit_plate3_standards.pl by Rohan Maddamsetti.
# This script makes a copy of the 091130_VeraCode_A-1_assay.tab file
# entitled 100823_VeraCode_A-1_assay3.tab. The only difference is that
# the Sample ID columns is changed to reflect the Sample ID columns in the
# data input file, "Plate_3_FinalReport_no_header.txt".

use strict;

open PROPER_HEADERS, "<input/Plate_3_FinalReport_no_header.txt" or die "cannot open file for input $!";

my %new_sample_id_column;
while (<PROPER_HEADERS>){
  chomp;
  my @line_data = split('\t');
  my $sample_id = $line_data[2];
  # the key is the last 9 characters of the string. Format: "R001_C001".
  my $str_length = length($sample_id);
  my $key = substr($sample_id, $str_length - 9 , 9);
#  print $key, "\n";
  $new_sample_id_column{$key} = $sample_id;
#print $sample_id, "\n";
#  print $_,"\n";
}
#print "@sample_id_column";

open NEW_STANDARDS, ">input/100823_VeraCode_A-1_assay3.tab" or die "cannot open file for output: $!";
open OLD_STANDARDS, "<input/091130_VeraCode_A-1_assay.tab" or die "cannot open file for input: $!";

while (<OLD_STANDARDS>){
  chomp;
#  print $_,"\n";
  my @standards_line = split('\t');
  my $old_sample_id = $standards_line[1];
  #  print $old_sample_id,"\n";
  my $old_str_length = length($old_sample_id);
  my $lookup = substr($old_sample_id, $old_str_length - 9, 9);
  $standards_line[1] = $new_sample_id_column{$lookup};
  my $output_line = join("\t",@standards_line);
#  print $output_line, "\n";
  print NEW_STANDARDS $output_line, "\n"; 
}
