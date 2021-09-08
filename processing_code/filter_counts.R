################################################################################
#### Project: TP Network
#### Title:   Function | Small function | Filter count matrix
#### Author:  Tom Walker (thomas.walker@usys.ethz.ch)
#### Date:    9 December 2020
#### ---------------------------------------------------------------------------

filter_counts <- function(x, threshold){
  # generate presence/absence data
  pa <- apply(x > 0, 2, as.numeric)
  # calculate coverage across samples
  cover <- rowSums(pa) / ncol(pa) * 100
  # change into index
  out <- cover > threshold
  # return
  return(out)
}
