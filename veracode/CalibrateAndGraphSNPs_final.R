##CalibrateAndGraphSNPs_final.R by Rohan Maddamsetti

##This is the final, polished script that analyzes all three mixed pop. datasets.
##CalibrateAndGraphClones.R should be run before this script is run (for the code
##that compares clone and mixed. pop genotyping datasets).

##import functions from my library.
source("GenotypingLibrary.R")

##This is the first set of calibration standards.
first.calibration.standards <- read.table("input/091130_VeraCode_A-1_assay.tab",sep="\t",header=T)
##This is the first run of the mixed population data.
first.plate.data <- read.table("input/Plate_1_FinalReport_no_header.txt",sep="\t",header=T)
##This is the second set of calibration standards.
second.calibration.standards <- read.table("input/100823_VeraCode_A-1_assay3.tab",sep="\t",header=T)
## This is the second run of the mixed population data.
second.plate.data <- read.table("input/Plate_3_FinalReport_no_header.txt",sep="\t",header=T)
## This is the third set of calibration standards.
third.calibration.standards <- read.table("input/110208_VeraCode_A-1_assay4.tab",sep="\t",header=T)
## This is the third run of the mixed population data.
third.plate.data <- read.table("input/Plate_4_FinalReport_no_header.txt",sep="\t",header=T)

## Remove bad samples from the data and calibration standards:
## IMPORTANT! These datapoints are anomalous, but they may not actually be bad!
## Caution must be exercised in assuming that these are bad data.
first.bad.ones <- c("Lenski Plate 1_R005_C011")
first.calibration.standards <- removeBadStandards(first.bad.ones, first.calibration.standards)
first.plate.data <- removeBadData(first.bad.ones, first.plate.data)

## Remove the 21K datapoint from the second mixed pop. run.
second.bad.ones <- c("LTEE Genotyping II_R003_C002")
second.calibration.standards <- removeBadStandards(second.bad.ones, second.calibration.standards)
second.plate.data <- removeBadData(second.bad.ones, second.plate.data)

## Remove the 24.5K datapoint from the third mixed pop. run.
third.bad.ones <- c("Lenski Plate 4_R008_C008")
third.calibration.standards <- removeBadStandards(third.bad.ones, third.calibration.standards)
third.plate.data <- removeBadData(third.bad.ones, third.plate.data)

## Remove rearrangements from the plates.
first.plate.data <- removeRearrangements(first.plate.data)
second.plate.data <- removeRearrangements(second.plate.data)
third.plate.data <- removeRearrangements(third.plate.data)

## Remove all bad SNP Assays, as listed in BadSNPAssays.txt.
first.plate.data <- removeBadSNPAssays(first.plate.data)
second.plate.data <- removeBadSNPAssays(second.plate.data)
third.plate.data <- removeBadSNPAssays(third.plate.data)

##Now, weight the allele frequencies in the mixed clone data in the calibration
##standards by the ratio of 8593A DNA in the mixed clone DNA.
##This factor was calculated from flow cytometry data.
scaled.first.calibration.standards <- scaleMixedClones(first.calibration.standards, weight=1.67)
scaled.second.calibration.standards <- scaleMixedClones(second.calibration.standards, weight=1.67)
scaled.third.calibration.standards <- scaleMixedClones(third.calibration.standards, weight=1.67)
## Now calibrate the data based on the standards.
first.calibrated <- calibrateSNPs(first.plate.data, scaled.first.calibration.standards)
second.calibrated <- calibrateSNPs(second.plate.data, scaled.second.calibration.standards)
third.calibrated <- calibrateSNPs(third.plate.data, scaled.third.calibration.standards)
all.calibrated <- rbind(first.calibrated,second.calibrated,third.calibrated)
averaged.plates <- averageRelevantMixedPlateData(all.calibrated)
almost.final.data1 <- averaged.plates
## remove mutT and mutations with little data from the dataset. 
almost.final.data2 <- subset(almost.final.data1, !almost.final.data1$SNP.Name %in% c("Ara-1-mutT-ins-u", "Ara-1-tdcR-yhaB-snp", "Ara-1-mgrB-yobH-snp", "Ara-1-x-thiC-snp", "Ara-1-uvrB-del-d", "Ara-1-yaaH-snp", "Ara-1-ycbX-ycbY-snp", "Ara-1-nirC-snp", "Ara-1-iclR..2-snp", "Ara-1-ycgH-snp"))
## cut the data to 20,000 generations.
final.data <- subset(almost.final.data2, !almost.final.data2$generation %in% c(20500,21000,21500,22000,22500,23000,23500,24000,24500,25000,25500,26000,26500,27000,27500,28000,28500,29000,29500,30000))


#################################################################################

## Generate figures from calibrated data.

######## Figure 1. 

## Make plots for Steve to post-process into a Muller plot with Illustrator.

## First, I need to divide up the data into fixations and non-fixations.
## Both calibrated and non-calibrated data will be graphed.
## Save the movie graphs as a set of *.svg files for post-processing with
## Adobe Illustrator.

## Adjustable parameters for plotsForMuller:
## 1) length of plot
## 2) width of plot
## 3) with or without error bars
## 4) with or without gene name.

all.output <- "/Users/Rohandinho/Desktop/MullerMutations.pdf"
plotsForMuller(final.data, all.output, use.theta=FALSE, width.parameter=6, height.parameter=3,bottomline=TRUE)

## Make series without a bottom line for Figure S1.

S1.output <- "/Users/Rohandinho/Desktop/S1MullerMutations.pdf"
plotsForMuller(final.data, S1.output, use.theta=FALSE, width.parameter=6, height.parameter=3,bottomline=FALSE)

######## Supplementary Figures.

## Figure S4: Example graph of effect of calibration to fix calibration.
S4.output <- "/Users/Rohandinho/Desktop/FigureS4.pdf"
graphS4(final.data, S4.output)

##Plot the pyrosequencing results that Jeff Barrick got.
#pyrosequencing.data <- read.csv("/Users/Rohandinho/Desktop/Projects/Ara-1_Genotyping/veracode/input/PyrosequencingResultsSummary.csv", header=T)
#plotPyrosequencingResults(pyrosequencing.data,
#                          outfile="/Users/Rohandinho/Desktop/pyroPlot.pdf")
##Plot the same results, but from the Veracode Genotyping data.
#makeGenotypingPyroPlot(final.data, outfile="/Users/Rohandinho/Desktop/GenotypingPyroPlot.pdf")

##Make a scatterplot comparing the relevant mixed pop. genotyping data to
##the pyrosequencing results.
#compareSequencingToGenotyping(final.data,pyrosequencing.data,
#                              outfile="/Users/Rohandinho/Desktop/pyroCompPlot.pdf")

##Test code for makeTimeHistogram.
#time.hist1 <- makeTimeHistogram(final.data)
#time.hist2 <- makeTimeHistogram(final.data, omitfirst2K=TRUE)
#write.csv(time.hist1, file="/Users/Rohandinho/Desktop/with2K.csv")
#write.csv(time.hist2, file="/Users/Rohandinho/Desktop/without2K.csv")


##Compare 10K clone allele frequencies with the mixed 10K allele frequencies.
#clone.ten.k.freqs <- read.csv("10K_allele_frequencies.csv",header=TRUE)
#regression.10K <- compareCloneAndMixedFreqs(clone.ten.k.freqs, final.data, gen=10000,
#                          outfile="/Users/Rohandinho/Desktop/10KClone-Mixed-Comparison.pdf")

##Compare 7.5K clone allele frequencies with the mixed 7.5K allele frequencies.
#clone.7.5k.freqs <- read.csv("7.5K_allele_frequencies.csv", header=FALSE)
#regression.7.5k <- compareCloneAndMixedFreqs(clone.7.5k.freqs, final.data, gen=7500,
#                          outfile="/Users/Rohandinho/Desktop/7.5KClone-Mixed-Comparison.pdf")

##Superimpose graphs from all mixed plate datasets.
#compareMixedPopDatasets(list(first.calibrated,second.calibrated,third.calibrated))

## Graph all SNPs.
#graphAllSNPs(final.data, "/Users/Rohandinho/Desktop/finalfudge.pdf")

## Make a scatterplot comparing allele frequencies as measured in previous
## mixed pop. sequencing (Jeff's Cold Spring Harbor paper results).
#mixed.pop.sequencing <- read.csv("input/MixedPopSequencingComparison/FinalMixedPopSequencingSNPs.csv", header=T)
#compareCSHresultsToGenotyping(final.data, mixed.pop.sequencing)
