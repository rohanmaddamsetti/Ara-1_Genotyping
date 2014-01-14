##  bpopsim_time_null_analysis.R by Rohan Maddamsetti
##  Input: sorted_cat_longest_sweep_time.tab.


## calculate.p.value(input.data.frame, y)
## Input: the data from the csv file as a dataframe, and a threshold generation y.
## Output: the p-value of a sweep time as extreme as y.
## In the end, Y should come from the Veracode data (Y=9000)
calculate.p.value <- function(the.data, y=9000) {

  ## For each Generation, make a logical vector of Generations > y. Add the trues together,
  ## and divide by the number of replicates (the length) to get the p-value.
  ## Produce a vector as long as a column in the data frame.

  p.val <- sum(sapply(the.data$Generations, function (x) ifelse(x>=y, 1, 0)))/length(the.data$Generations)
  return(p.val)
}

## Function 2: plot the distribution of Generations, and color in the
## p-value tail of the distribution, given a Y.
plot.null.distribution <- function(the.data, y=9000) {

  library(ggplot2)
  
  the.data$Tail <- sapply(the.data$Generations, function (gen) ifelse(gen>=y, TRUE, FALSE))

  this.plot <- ggplot(the.data, aes(x=Generations, fill=Tail)) + geom_histogram(binwidth=500)
  plot.filename = "/Users/Rohandinho/Desktop/bpopsim_sweeptime_null.pdf"
  ggsave(plot.filename)
}

## This is where the analysis runs.

data <- read.csv("input/sorted_cat_longest_sweep_time.tab", header=FALSE)
names(data) <- c("Transfers")
data$Generations <- sapply(data$Transfers, function (x) return(x*6.6666667))

calculate.p.value(data)
plot.null.distribution(data)
