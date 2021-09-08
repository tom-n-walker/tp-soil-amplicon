################################################################################
#### Project: TP Network
#### Title:   Function | Wrangle | Apply site-wise NMDS
#### Author:  Tom Walker (thomas.walker@usys.ethz.ch)
#### Date:    16 February 2021
#### ---------------------------------------------------------------------------

sitewise_nmds <- function(seqData){
  # Split OTU data by site
  splitData <- phyloseq_sep_variable(seqData, "site", drop_zeroes = F)
  # Apply NMDS to each split
  allNMDS <- lapply(splitData, do_nmds)
  # Extract scores
  allScores <- lapply(allNMDS, function(x) x$points)
  # Bind scores
  output <- do.call(rbind, allScores) %>%
    as.data.frame %>%
    rownames_to_column("ID")
  # return
  return(output)
}
