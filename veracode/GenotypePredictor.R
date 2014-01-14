# GenotypePredictor.R by Rohan Maddamsetti.

#############################################################################################
# NOTE!!!
# This script is unfinished. See my purple notebook for solution to intermediate alleles
# on the line of descent. This has been left unfinished since Jeff thinks it would be more
# valuable to graph multiple SNPs at the same time, and so infer genotypes from the order of
# fixations/changes in allele frequency, rather than figuring it out 'analytically' in this
# script.

# I've kept this script, since its current output is useful for keeping track of which SNPs
# fixed at what time.
#############################################################################################
# Input: predicted SNP trajectories, from MixedPop_SNP_FittedValues.tab,
# and a list of bad SNPs to be removed from the data.
# MixedPop_SNP_FittedValues.tab is an output file from CalibrateAndGraphSNPs.R.

# Output: Predicted frequencies of genotypes at various times in the LTEE.

# This script predicts LTEE population structure from mixed population
# Illumina genotyping data. From predicted allele frequencies, it
# attempts to figure out genotypes present at different timepoints in
# the LTEE, and the frequencies of these genotypes.

# A simple heuristic is used to figure out genotypes, as follows.
# If SNP (1) increases in frequency, and SNP (2) doesn't change in
# frequency, then either A) SNP (2) is fixed, and SNP (1) occurred
# on that background, or B) SNP (2) has not occurred yet.
# Mike Wiser has also pointed out that (1) might increase at the
# expense of (2) but not at the expense of (3), if neutral with
# regard to (3). I might then mistakenly put (1) on the background
# of (3).
# If SNP (1) increases in frequency, and SNP (2) increases concomitantly,
# then SNP (1) occurred on a SNP (2) background, and (2) has not fixed.
# If SNP (1) increases in frequency, and SNP (2) decreases concomitantly,
# then SNP (1) did not occur on a SNP (2) background, and (2) had not fixed.

# The most important error case is if the allele frequencies are reversed.
# This script does not do this error checking on the data.

# At any time interval, at most 1 unlinked SNP can go to fixation.

snp.trajectory.data <- read.table("MixedPop_SNP_FittedValues.tab",sep="\t",header=T)
#bad.snps <- c("Ara-1-manB-cpsG-xst-1","Ara-1-aroA-ycaL-del-u","Ara-1-atoC-snp",
#              "Ara-1-glmU-atpC-ins-d","Ara-1-manB-cpsG-xst-1","Ara-1-nmpC-ECB_00513-xst-2",
#              "")

# A genotype is a character vector of alleles, plus a genotype frequency.
setClass("Genotype",
         representation(
                        alleles="character",
                        frequency="numeric"
                        )
         )

generations <- unique(snp.trajectory.data$generation)
# The data structure for the genotype reconstructions is a list of key-value pairs,
# where the keys are values in the 'generations' vector, and the values are lists
# of Genotype objects.
#reconstructions <- list()
#for (gen in generations){
#  reconstructions[as.character(gen)] <- list()
#}

line.of.descent <- list()
# Go backwards in time, then I know what is on the eventual line of descent, i.e. what has been fixed.
# By definition, these are SNPs that have been calibrated.
for (gen in rev(generations)) {
  gen.data <- subset(snp.trajectory.data,generation==gen)
  fixed.snps <- gen.data[gen.data$predicted.b.allele.freq.fit > 0.9,]$SNP.Name
  # Now remove NA values from fixed.snps, and coerce to a vector of characters.
  # In the case of replicates from the same time point, remove duplicates, too.
  fixed.snps <- unique(as.character(fixed.snps[!is.na(fixed.snps)]))
  fixed.genotype <- new("Genotype",alleles=fixed.snps,frequency=1)
  print(gen)
  print(fixed.genotype)
  line.of.descent[as.character(gen)] <- fixed.genotype
}
line.of.descent <- rev(line.of.descent)
# Now handle missing data points for already fixed SNPs.
for (i in 2:length(generations)) {
  previous.gen <- as.character(generations[i-1])
  current.gen <- as.character(generations[i])
  the.previous <- line.of.descent[[previous.gen]]@alleles
  the.current <- line.of.descent[[current.gen]]@alleles
  # SNPs in the last timepoint that are not in the current timepoint.
  in.previous <- the.previous[!(the.previous %in% the.current)]
  # And add these to the current timepoint.
  line.of.descent[[current.gen]]@alleles <- sort(c(the.current,in.previous))
}


## Now reconstruct intermediate genotype frequencies on the line of descent.
##for (gen in rev(generations)) {
##  gen.data <- subset(snp.trajectory.data,generation==gen)
##  middle.snps <- gen.data[gen.data$predicted.b.allele.freq.fit < 0.9
##                          & gen.data$predicted.b.allele.freq.fit > 0,
##                          c("SNP.Name","predicted.b.allele.freq.fit")]
##  # Remove NA values.
##  middle.snps <- subset(middle.snps,!is.na(SNP.Name) & !is.na(predicted.b.allele.freq.fit))
##  # In the case of replicates from the same timepoint, average over SNPs.
##  test <- XXX
##  # Rank intermediate alleles based on frequency.
##  middle.snps <- middle.snps[order(middle.snps$predicted.b.allele.freq.fit),]
##  # Generate genotype frequencies.
##  middle.genotypes <- list()
##  subtract.me <- 0
##  for (i in length(middle.genotypes)) {
##    this.snp <- middle.genotypes
##    
##    a.middle.genotype <- new("Genotype",alleles=)
##  }
##  
##  print(middle.snps)
##}




