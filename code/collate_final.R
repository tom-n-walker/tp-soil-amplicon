################################################################################
#### Project: TP Network
#### Title:   Function | Wrangle | Collapse final data
#### Author:  Tom Walker (thomas.walker@usys.ethz.ch)
#### Date:    16 February 2021
#### ---------------------------------------------------------------------------

collate_final <- function(soilData, fungi, bacFullNMDS, funFullNMDS, bacSubNMDS, funSubNMDS){
  # organise sample data
  samples <- sample_data(fungi) %>%
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
  samplesNMDS <- cbind(samples, allReady) %>%
    mutate(rep = as.numeric(rep_block)) %>%
    select(-rep_block)
  # add soil data
  output <- soilData %>%
    rename(treat_a = treatment) %>%
    left_join(., samplesNMDS) %>%
    select(
      site, gradient:site_id, 
           treat_a, treat_b:treat_c, rep, 
           elev_m:bacSub2, h2o_mgg:TRate_d
    )
  # return
  return(output)
}
