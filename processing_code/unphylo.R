################################################################################
#### Project: TP Network
#### Title:   Small function | Unclass phyloseq object
#### Author:  Tom Walker (thomas.walker@usys.ethz.ch)
#### Date:    21 February 2021
#### ---------------------------------------------------------------------------

unphylo <- function(x){
  as.data.frame(unclass(x))
}
