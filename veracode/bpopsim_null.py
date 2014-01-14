#!/usr/bin/python

## bpopsim_null.py by Rohan Maddamsetti.
## This script produces a *.csv file suitable for analysis with the R script
## bpopsim_null_analysis.R.

## Usage: python bpopsim_null.py ~/Desktop/bpopsim_out > bpopsim_mutation_divergence.csv")

## Notes from an email from Jeff Barrick:
## Jeff will send me files, for looking at how often a clade X or greater
## mutations off the line of descent to the final dominant reaches Y% of
## the population. The series of files is named 
## "mutation_divergence_depth_frequency_X.tab. Notice that it is "X or greater"
## and the frequency is for the biggest clade in existence that is that distance
## (in mutations) from where it diverged from the line of descent to the final
## dominant. Columns are time in transfers (75 = 500 generations, so this is
## realistically coarse-grained.) Each row is a different simulation run.
## I think you will want to take the maximum frequency in each row and use that
## as the test statistic for that replicate. I [Jeff] guess the X=3 file is the
## one of most interest?

## One issue that might need to be dealt with is that if there is not a sweep
## near the end of the simulation, the final dominant might not be the ultimate
## winner (it may be on its way out). There seem to generally be higher numbers
## at the tail end of the simulations for that reason, but I [Jeff] haven't
## investigated it fully. Maybe you could look at the data to see if this is the
## case. We could extend the simulations a lot at the end of see which would
## have been the eventual winner and get rid of this artifact.

from sys import argv ## for simple command-line argument processing
from sys import exit
from os import listdir

## input: a directory containing bpopsim summary statistic files
def main():

	output_rows = [] ## This is going to be a list of rows (which are lists themselves).
	output_header = [] ## for scope reasons.
	## must take only one argument: a directory containing input files.
	if (len(argv) != 2):
		exit("Usage: python bpopsim_null.py ~/Desktop/bpopsim_out > bpopsim_mutation_divergence.csv")
	my_dir = argv[1]
	files = listdir(my_dir)
	for fname in files:
		output_header = ["X"] ## This will be added to later.
		##open each file for reading if it matches the proper naming convention.
		prefix = "mutation_divergence_depth_frequency_"
		suffix = ".tab"
		if fname.startswith(prefix) and fname.endswith(suffix):
			## x is the number of mutations off the line of descent--named in the filename.
			x = fname[len(prefix):len(fname) - len(suffix)]
			fpath = my_dir + "/" + fname
			fhandle = open(fpath)
			out_row = [x] ## This list contains the relevant info from the file.
			for lineno, line  in enumerate(fhandle, 1):
				## remove trailing newlines.
				line = line.strip()
				## skip the header.
				if lineno == 1:
					continue
				row = line.split("\t")
				##print row
				replicate = row[0]
				frequencies = row[1:]
				##print frequencies
				maxfreq = max(frequencies)
				#print replicate, maxfreq
				## this next line is a bit stupid, all files must have the same number
				## of replicates as the last file examined.
				output_header.append("Replicate" + replicate)
				out_row.append(maxfreq) 
			output_rows.append(out_row)
	print ", ".join(output_header)
	## sort the rows by X numerically, not alphanumerically.
	output_rows.sort(key=lambda this_row: int(this_row[0]))
	for r in output_rows:
		print ", ".join(r)
    
main()
