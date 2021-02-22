################################################################################
#### Project: TP Network
#### Title:   Function | Analyse | Do NMDS
#### Author:  Tom Walker (thomas.walker@usys.ethz.ch)
#### Date:    16 February 2021
#### ---------------------------------------------------------------------------

do_nmds <- function(x){
  # calculate distance matrix
  d <- vegan::vegdist(t(x@otu_table))
  # do NMDS
  nmds <- vegan::metaMDS(d, k = 2, trymax = 200)
  # return
  return(nmds)
}
