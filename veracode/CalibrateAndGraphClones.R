##CalibrateAndGraphClones.R by Rohan Maddamsetti
##This script analyzes the Veracode 7.5K and 10K Clone data.
##This script replaces the separate scripts for the 7.5K and 10K Clone data.

##Import libraries.
source("GenotypingLibrary.R")

## Functions.

printTables <- function(genotype.table, path="./", filename.prefix) {
 ## Write the table to an output file.
  write.table(genotype.table,file=paste(path, filename.prefix,"clone_genotypes.tab",sep=""),
              quote=F, sep="\t",row.names=F)
  
  allele.frequency.table <- countAlleleFrequencies(genotype.table)
  ## Write the table to an output file.
  write.csv(allele.frequency.table,file=paste(path, filename.prefix,"allele_frequencies.csv", sep=""),quote=F, row.names=F)
            
  genotype.freqs <- countGenotypeFrequencies(genotype.table)
  genotype.names <- genotypeStringRepr(genotype.table)
  ##Write out the genotype names.
  write.table(genotype.names, file=paste(path, filename.prefix, "genotype_names.txt",sep=""),
              quote=F, row.names=F, sep="\t")
}

makeHeatmap <- function(genotype.table, control.wells,  path="./", filename.prefix) {
 
  ##remove control rows of the genotype table.
  genotype.table <- genotype.table[!(genotype.table$Well %in% control.wells),]
  
  row.names <- genotype.table$Well
  no.well.genotype.table <- subset(genotype.table,select = names(genotype.table) != "Well")

  genotype.matrix <- as.matrix(no.well.genotype.table)
  
  ##Remove snps that wont't appear in the 10K heatmap plot.
  keep.these.snps <- c("Ara.1.atoS.snp", "Ara.1.nuoG.snp", "Ara.1.leuO.ilvI.snp", "Ara.1.nadR..2.snp", "Ara.1.elaD.snp", "Ara.1.hsdM.snp", "Ara.1.maeB.talA.snp", "Ara.1.pflC.del.d", "Ara.1.araJ.snp", "Ara.1.yhdG.fis.snp","Ara.1.yghJ.snp", "Ara.1.rpsM.snp", "Ara.1.yedW.yedX.snp", "Ara.1.gltB.del.d", "Ara.1.nuoM.snp", "Ara.1.acs.nrfA.snp", "Ara.1.rpsA.snp", "Ara.1.topA.snp", "Ara.1.mrdA.snp", "Ara.1.infB.snp", "Ara.1.ybaL.ins.d", "Ara.1.yegI.snp", "Ara.1.malT.snp", "Ara.1.spoT.snp", "Ara.1.nagC.snp")

 genotype.matrix <- subset(genotype.matrix, select=keep.these.snps)  
  
  ##Remove degenerate columns of the genotype matrix (where no snps appear at all).
  genotype.matrix <- genotype.matrix[,colSums(genotype.matrix, na.rm=TRUE) > 0]

  ##This sets the labels for the heatmap diagram.
  rownames(genotype.matrix) <- row.names
  ##Remove the "Ara.1." prefix from the names of the snps.
  col.names <- colnames(genotype.matrix)
  col.names <- sapply(col.names, function (x) substring(x,7))
  colnames(genotype.matrix) <- col.names
  ## Get initial heatmap return values for manual reordering for better readability.
  heatmap.returnvals <- heatmap(genotype.matrix,col=topo.colors(2),cexRow=0.3, cexCol=0.8, Colv=NA,)
  dev.off()
  col.inds <- heatmap.returnvals$colInd
  row.inds <- heatmap.returnvals[1]$rowInd
 
  ##initialize this matrix.
  new.order.genotype.matrix <- genotype.matrix
  if (filename.prefix == "7.5K_") {
    ##This manually reorders the samples of the 7.5K samples, so the rows are more like the 10K heatmap.
    new.order <- c(row.inds[25:57], row.inds[72:90], row.inds[58:71], row.inds[1:24])
    new.order.names <- sapply(new.order, function (x) row.names[x])
    new.order.genotype.matrix <- t(cbind(sapply(new.order.names, function (x) genotype.matrix[x,])))
  }
  else { ##For 10K Clones--manually reorder clones with atoS and nuoG mutations.
    new.order <- c(row.inds[1:8],row.inds[9:90])
    new.order.names <- sapply(new.order, function (x) row.names[x])
    new.order.genotype.matrix <- t(cbind(sapply(new.order.names, function (x) genotype.matrix[x,])))
    ##put atoS and nuoG snps further down in the heatmap for readability.
    new.col.order <- c(col.inds[3:14],col.inds[2],col.inds[1],col.inds[15:25])
    new.col.order.names <- sapply(new.col.order, function (x) col.names[x])
    new.order.genotype.matrix <- cbind(sapply(as.vector(new.col.order.names), function (x) new.order.genotype.matrix[,x]))

}
  pdf(paste(path, filename.prefix, "heatmap.pdf",sep=""))
    ##Make the heatmap diagram.
  heatmap(new.order.genotype.matrix,col=c("grey90", "lightgreen"), cexRow=0.3, cexCol=0.75, Colv=NA , Rowv=NA)
  dev.off()
}

##Import data.
tenK.plate.data <- read.table("input/Plate_2_FinalReport_no_header_10K_Clones.txt",sep="\t",header=T)
seven.point.fiveK.plate.data <- read.table("input/Plate_5_FinalReport_no_header_7.5K_Clones.txt",sep="\t",header=T)

##Remove all bad SNPs, as listed in BadSNPAssays.txt.
tenK.plate.data <- removeBadSNPAssays(tenK.plate.data)
seven.point.fiveK.plate.data <- removeBadSNPAssays(seven.point.fiveK.plate.data)

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

tenK.control.wells <- c("R008_C006", "R005_C003", "R006_C001", "R006_C006", "R004_C001", "R002_C005")
seven.point.fiveK.control.wells <- c("R001_C011", "R003_C002", "R003_C008", "R005_C009", "R001_C005", "R004_C009")

tenK.genotype.table <- makeGenotypeTable(tenK.plate.data, snps.to.reverse, tenK.control.wells, FALSE)
seven.point.fiveK.genotype.table <- makeGenotypeTable(seven.point.fiveK.plate.data, snps.to.reverse, seven.point.fiveK.control.wells, FALSE)

printTables(tenK.genotype.table, path="/Users/Rohandinho/Desktop/", filename.prefix="10K_")
printTables(seven.point.fiveK.genotype.table, path="/Users/Rohandinho/Desktop/", filename.prefix="7.5K_")

##Don't forget to end the pathname with a slash!

makeHeatmap(tenK.genotype.table, tenK.control.wells, path="/Users/Rohandinho/Desktop/", filename.prefix="10K_")

makeHeatmap(seven.point.fiveK.genotype.table, seven.point.fiveK.control.wells, path="/Users/Rohandinho/Desktop/", filename.prefix="7.5K_")
