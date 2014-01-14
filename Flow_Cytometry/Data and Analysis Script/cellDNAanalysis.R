##cellDNAanalysis.R by Rohan Maddamsetti.
##This script analyzes the raw flow cytometry data in .fcs format
##that represents DNA content per cell in REL606 and in REL8593A, which is a
##clone isolated from Ara-1 at 20,000 generations.
##I want to find the fold-increase in DNA content per cell after 20,000
##generations of evolution.

##This script uses packages from Bioconductor that are amenable for flow
##cytometry analysis.

## Load packages.

library(flowCore)
library(flowStats)
library(flowQ)
library(flowViz) # for flow data visualization.

## DM25.flowData is the data for DM25. This is relevant for the LTEE, but not for the Veracode genotyping.
DM25.flowData <- read.flowSet(path="./FACS080311/80311FCS3", alter.names=TRUE, transformation=FALSE)
## flowData is the data from the DM100 experiment. This is what is relevant for the Veracode genotyping.
flowData <- read.flowSet(path="./FACS081911/RM081911", alter.names=TRUE, transformation=FALSE)
ancestral.data <- flowData[[1]]
evolved.data <- flowData[[3]]
##plot a histogram of the fluorescence.
quartz()
plot(ancestral.data, "FL.1.A", breaks = 1000)
plot(evolved.data, "FL.1.A", breaks = 1000, ylim=c(0,400))

##gate the data based on the FSC-SSC scatterplot.
##find a region that most resembles a
##bivariate Normal distribution.
n2filter <- norm2Filter("FSC.A", "SSC.A", scale=2)
anc.n2filter.results <- filter(ancestral.data,n2filter)
xyplot(FSC.A ~ SSC.A, data=ancestral.data, filter=anc.n2filter.results)
gated.anc <- Subset(ancestral.data, subset=anc.n2filter.results)
plot(gated.anc, "FL.1.A", breaks=1000)
evo.n2filter.results <- filter(evolved.data, n2filter)
gated.evo <- Subset(evolved.data, subset=evo.n2filter.results)
xyplot(FSC.A ~ SSC.A, data=evolved.data, filter=evo.n2filter.results)
plot(gated.evo, "FL.1.A", breaks=1000)

##integrate FL-1 data for both experiments.
##The rows of the data matrix in a flowFrame object
##represent the events, and the columns represent
##the various parameters of interest.
mean.anc.fluor <- each_col(gated.anc$FL.1.A, mean)
sd.anc.fluor <- each_col(gated.anc$FL.1.A, sd)
std.err.anc.fluor <- sd.anc.fluor/sqrt(nrow(gated.anc))
mean.evo.fluor <- each_col(gated.evo$FL.1.A, mean)
sd.evo.fluor <- each_col(gated.evo$FL.1.A, sd)
std.err.evo.fluor <- sd.evo.fluor/sqrt(nrow(gated.evo))
##divide 20K by 606, and calculate the standard error.
fold.increase <- mean.evo.fluor/mean.anc.fluor
##The fold increase is 1.796752 for DM25.
##The fold increase is 1.879978 for DM100.
##The standard error calculatation comes from
##the following discussion on error analysis in
##physics:
##http://teacher.pas.rochester.edu/PHY_LABS/AppendixB/AppendixB.html
std.error.fold.increase <- fold.increase*sqrt((std.err.anc.fluor/mean.anc.fluor)^2 + (std.err.evo.fluor/mean.evo.fluor)^2)
conf.interval <- c(fold.increase-1.96*std.error.fold.increase, fold.increase+1.96*std.error.fold.increase)
