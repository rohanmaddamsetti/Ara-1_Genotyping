##findBadDataPredictors.R by Rohan Maddamsetti.
##Find a method to filter out the bad data:

##  Write a script that takes all the bad, spiked points, then graph the
##distribution of all of the original, raw data for these versus the rest of 
##the "good" data.
##I want to find some predictor for when the data is bad, and then I can use
## this to filter the data for the analysis.

source("GenotypingLibrary.R")

##Import all five data sets.
tenK.plate.data <- plate.data <- read.table("input/Plate_2_FinalReport_no_header_10K_Clones.txt",sep="\t",header=T)
seven.point.fiveK.plate.data <- read.table("input/Plate_5_FinalReport_no_header_7.5K_Clones.txt",sep="\t",header=T)

##I need the mixed.plate.data to look at the data for the ancestor and derived clones.
##This is to see whether the assay data needs to be reversed or not.
old.mixed.plate.data <- read.table("input/Plate_1_FinalReport_no_header.txt",sep="\t",header=T)
old.calibration.standards <- read.table("input/091130_VeraCode_A-1_assay.tab",sep="\t",header=T)
new.mixed.plate.data <- read.table("input/Plate_3_FinalReport_no_header.txt",sep="\t",header=T)
new.calibration.standards <- read.table("input/100823_VeraCode_A-1_assay3.tab",sep="\t",header=T)
newest.mixed.plate.data <- read.table("input/Plate_4_FinalReport_no_header.txt" ,sep="\t",header=T)
newest.calibration.standards <- read.table("input/110208_VeraCode_A-1_assay4.tab", sep="\t",header=T)

ancestors1 <- getAncestorData(old.mixed.plate.data,old.calibration.standards,FALSE)
ancestors2 <- getAncestorData(new.mixed.plate.data,new.calibration.standards,FALSE)
ancestors3 <- getAncestorData(newest.mixed.plate.data, newest.calibration.standards, FALSE)
all.ancestors <- data.frame(rbind(ancestors1,ancestors2, ancestors3))
snps.to.reverse <- reverseSNP(all.ancestors)

##This function maps time points to sample ID. should be helpful
##for finding bad datapoints from my graphs.
TimetoSampleID <- function(calibration.data) {
  time.to.sample = subset(calibration.data,select=c(Sample.ID,Gen))
  return(time.to.sample)
}

##This function makes a graph that replicates the normalized
##SNP graphs in the GenomeStudio software.
##Example: snp.name <- factor("Ara−1−gltI−snp"). Another reason why I hate R.
GenomeStudioGrapher <- function(snp.name, plate.data) {
  quartz()
  snp.data <- plate.data[plate.data$SNP.Name==snp.name,]
  ##This next line accurately recapitulated the normalized SNP graph in GenomeStudio.
  plot(x=snp.data$Theta, y=snp.data$R)
#  plot(x=snp.data$X, y=snp.data$Y)
#  plot(x=snp.data$X.Raw, y=snp.data$Y.Raw)
  dev.off()
}

