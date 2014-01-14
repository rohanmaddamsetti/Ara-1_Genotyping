## How to use this script:
##
## Two parameters:
##   snp=<SNP>, name of snp, file names "<SNP>.tab" expected
##   rev=[T/F], optional, ancestral/evolved order reversed?
##
## Examples:
## > R --vanilla snp=nagC < calibrate_snp.r
## > R --vanilla snp=mrdA rev=T < calibrate_snp.r

rev=F;

commandArgs()

for (e in commandArgs()) {
  ta = strsplit(e,"=",fixed=TRUE)
  if(! is.na(ta[[1]][2])) {
    temp = ta[[1]][2]
    if(substr(ta[[1]][1],nchar(ta[[1]][1]),nchar(ta[[1]][1])) == "I") {
      temp = as.integer(temp)
    }
    if(substr(ta[[1]][1],nchar(ta[[1]][1]),nchar(ta[[1]][1])) == "N") {
      temp = as.numeric(temp)
    }
    assign(ta[[1]][1],temp)
    cat("assigned ",ta[[1]][1]," the value of |",temp,"|\n")
  } else {
    assign(ta[[1]][1],TRUE)
    cat("assigned ",ta[[1]][1]," the value of TRUE\n")
  }
}

data<-read.table(paste(snp, ".tab", sep=""), header=T)
if(rev)
{
	data$measured <- 1 - data$measured;
}
##I (Rohan) should add code here to deal with down curves as well as up curves: see where reverse flag is used, or just look at endpoints of the curve.

standard<-subset(data, actual != -1) ##May not need this line. 
standard$residual <- standard$actual-standard$measured
fit<-lm(residual~0+I(measured-measured^2)+I(measured-measured^4)+I(measured-measured^6)+I(measured-measured^8)+I(measured-measured^10), standard)

standard$predicted <- predict(fit, standard, level=0.95, interval='confidence') + standard$measured

data$predicted <- predict(fit, data, level=0.95, interval='confidence') + data$measured
data.sorted<-data[order(data$order),] 

#graphs
pdf(paste(snp, "-output.pdf", sep=""))

##combine the following two graphs and use a different symbol to differentiate them.
plot(standard$actual,standard$measured, main="Mutation Frequency Response Curve - Raw Data", xlab="known frequency", ylab="predicted frequency")
lines(cbind(0,1), cbind(0,1))

plot(standard$actual,standard$predicted[,1], main="Mutation Frequency Response Curve - After Calibration", xlab="known frequency", ylab="predicted frequency")
lines(cbind(0,1), cbind(0,1))

#graph of all samples
mycolors = cbind("red", "blue", "orange", "green", "purple");
mycolorlist = c();
for (i in 1:length(data.sorted$measured))
{
	mycolorlist = rbind(mycolorlist, mycolors[as.numeric(data.sorted$color[i])+1]);
	data.sorted$generation[i]<-as.numeric(data.sorted$generation[i]);
}
data.sorted$graph_colors = mycolorlist;


#subset and sort by generation
data.mixed_pops = data.sorted[c(1,26:96),];
plot(data.mixed_pops$generation, data.mixed_pops$predicted[,1], main="Predicted Mutation Frequencies - Mixed Populations", xlab="generation", ylab="mutation frequency")
lines(data.mixed_pops$generation, data.mixed_pops$predicted[,1])


barplot(data.sorted$predicted[,1], col=data.sorted$graph_colors, main="Predicted Mutation Frequencies - All Samples", xlab="sample ID (1-96 = A1, A2 ... H11, H12)", ylab="mutation frequency", xaxt="n")

dev.off()