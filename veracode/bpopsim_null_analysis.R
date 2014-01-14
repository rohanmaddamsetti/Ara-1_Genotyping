## bpopsim_null_analysis.R by Rohan Maddamsetti.
## This script analyzes a *.csv file produced by the python script bpopsim_null.py.

## This file is called bpopsim_mutation_divergence.csv.

## calculate.p.values(input.data.frame, y)
## Input: the data from the csv file as a dataframe, and a threshold frequency y.
## Output: the p-value of a frequency as extreme as y, for all X in the table,
## where X is the number of mutations off the line of descent.
## In the end, Y should come from the Veracode data.
calculate.p.values <- function(the.data, y=0.90) {

  ## For each row, make a logical vector of frequencies > y. Add the trues together,
  ## and divide by the number of replicates (the length) to get the p-value.
  ## Produce a vector as long as a column in the data frame.
  
  just.frequencies <- the.data
  just.frequencies$X <- NULL
  
  p.vals <- apply(just.frequencies, 1, function (row) sum(sapply(row, function(x) ifelse(x >= y, 1, 0)))/length(row))
  
  return(p.vals)
}

## Function 2: plot the distribution of maximum frequencies, and color in the
## p-value tail of the distribution, for all X, given a Y.
plot.null.distributions <- function(the.data, y=0.90) {

  library(ggplot2)

  ## Iterate over all rows in the data.
  for (i in 1:nrow(the.data)) {
    current.row <- the.data[i,]
    ## From this row, make a dataframe suitable for ggplot2.
    x = as.vector(unlist(current.row[1]))
    row.without.x <- current.row[2:length(current.row)]
    data.to.plot <- data.frame(Replicate=1:length(row.without.x), Frequency=as.vector(unlist(row.without.x)),Tail=sapply(row.without.x, function(freq) ifelse(freq >= y, TRUE, FALSE)))
    this.plot <- ggplot(data.to.plot, aes(x=Frequency, fill=Tail)) + geom_histogram(binwidth=0.05)
    plot.filename <- sprintf("%sNullMutationDivergence%d.pdf", "/Users/Rohandinho/Desktop/bpopsim_plots/", x)
    ggsave(plot.filename)
  }
}

data <- read.csv("bpopsim_mutation_divergence.csv")
