################################################################################
#### Project: TP Network
#### Title:   Function | Wrangle | Collapse final data
#### Author:  Tom Walker (thomas.walker@usys.ethz.ch)
#### Date:    16 February 2021
#### ---------------------------------------------------------------------------

collate_final <- function(siteData, sampleData, soilData, fungi, bacFullNMDS, funFullNMDS, bacSubNMDS, funSubNMDS){
  ## Format data ----
  # Subset and format OTU sample data
  otuSamples <- sample_data(fungi) %>%
    # unclass phylogseq object
    unclass %>%
    as.data.frame %>%
    # arrange by site to make compatible with NMDS objects
    arrange(site) %>%
    # rename columns and subset
    mutate(rep = as.numeric(rep_block)) %>%
    select(site, treat_a, rep)
  # Format NMDS data
  nmdsScores <- funFullNMDS %>%
    # subsequent joins through all NMDS objects
    left_join(., funSubNMDS, by = "ID") %>%
    left_join(., bacFullNMDS, by = "ID") %>%
    left_join(., bacSubNMDS, by = "ID") %>%
    # select columns
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
    ) %>%
    # bind to sample data
    cbind(otuSamples, .)
  # Format experiment sample data
  allSamples <- sampleData %>% 
    # rename columns and subset
    mutate(rep = as.numeric(rep_block)) %>%
    select(site, gradient:site_id, lon, lat, elev_cat, elev_m, treat_a:treat_c, rep, start_year:P_ann)
  # Subset and format soil data
  soil <- soilData %>%
    rename(treat_a = treatment) %>%
    select(site, treat_a, rep, h2o_mgg:rpa_mom)
  # Bind together
  output <- soil %>%
    left_join(allSamples, .) %>%
    left_join(., nmdsScores) %>%
    # remove YAN - mistake with sampling so not comparable
    filter(site != "YAN")
  
  ## Summarise soil responses with PCA ----
  # do
  soilPCA <- output %>%
    # get data
    select(pH:CmicNmic, R_ugCgh:rpa_mom, funFull1:funFull2, bacFull1:bacFull2) %>%
    # impute missing values
    apply_mice(., 5) %>%
    # do PCA
    prcomp(., scale = T, center = T)
  # get scores
  output$soilPC1 <- soilPCA$x[, 1]
  output$soilPC2 <- soilPCA$x[, 2]
  
  ## Spread treatments to make meta-analyses metrics ----
  listOut <- output %>%
    # select variables and split to bins in list
    select(h2o_mgg:soilPC2) %>%
    as.list %>%
    # bind responses to categorical variables of interest
    lapply(., function(x){
      cbind(
        select(output, site, treat_a, rep),
        data.frame(response = x)
      )
    }) %>%
    # pivot wider
    lapply(., function(x){
      pivot_wider(
        x, 
        values_from = response, 
        names_from = treat_a
      )
    }) %>%
    # calculate divergence, convergence, etc.
    lapply(., function(x){
      x %>%
        group_by(site) %>%
        summarise(
          .groups = "keep",
          # get original values
          HH = HH,
          HL = HL,
          LL = LL,
          # get differences based on means of controls
          HLvHHmean = (HL - mean(HH, na.rm = T))/mean(HH, na.rm = T) * 100,
          HLvLLmean = (HL - mean(LL, na.rm = T))/mean(LL, na.rm = T) * 100,
          HHvLLmean = (HH - mean(LL, na.rm = T))/mean(LL, na.rm = T) * 100,
        ) %>%
        ungroup %>%
        as.data.frame %>%
        right_join(siteData, .)
      })
  
  ## Site-wise treatment means ± SEs ----
  # separate data frames
  means <- output %>%
    group_by(site, treat_a) %>%
    summarise(across(h2o_mgg:bacSub2, mean, na.rm = T)) %>%
    ungroup %>%
    rename_with(~paste0("mean_", .x), h2o_mgg:bacSub2)
  SEs <- output %>%
    group_by(site, treat_a) %>%
    summarise(across(h2o_mgg:bacSub2, ~sd(.x, na.rm = T)/sqrt(n()))) %>%
    ungroup %>%
    rename_with(~paste0("se_", .x), h2o_mgg:bacSub2)
  # bind
  summaries <- left_join(means, SEs)
  
  ## Return ----
  finished <- list(
    fullDF = output,
    meansSEs = summaries,
    metaMetrics = listOut
  )
  return(finished)
}
