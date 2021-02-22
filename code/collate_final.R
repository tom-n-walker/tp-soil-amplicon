################################################################################
#### Project: TP Network
#### Title:   Function | Wrangle | Collapse final data
#### Author:  Tom Walker (thomas.walker@usys.ethz.ch)
#### Date:    16 February 2021
#### ---------------------------------------------------------------------------

collate_final <- function(funSeq, bacFullNMDS, funFullNMDS, bacSubNMDS, funSubNMDS){
  # organise sample data
  samples <- sample_data(funSeq) %>%
    unclass %>%
    as.data.frame %>%
    arrange(site)
  # format NMDS data
  allReady <- funFullNMDS %>%
    left_join(., funSubNMDS, by = "ID") %>%
    left_join(., bacFullNMDS, by = "ID") %>%
    left_join(., bacSubNMDS, by = "ID") %>%
    select(
      -"ID",
      funFull1 = MDS1.x,
      funFull2 = MDS2.x,
      funSub1 = MDS1.y,
      funSub2 = MDS2.y,
      bacFull1 = MDS1.x.x,
      bacFull2 = MDS2.x.x,
      bacSub1 = MDS1.y.y,
      bacSub2 = MDS2.y.y
    )
  # bind together
  output <- cbind(samples, allReady)
  # return
  return(output)
}
