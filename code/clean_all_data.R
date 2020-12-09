################################################################################
#### Project: TP Network
#### Title:   Function | Wrangle | Format sequence data
#### Author:  Tom Walker (thomas.walker@usys.ethz.ch)
#### Date:    24 November 2020
#### ---------------------------------------------------------------------------

clean_all_data <- function(subSeqData, climData, soilData){
  ## Format count data ----
  # Index counts for non-missing taxa
  bacNotMiss <- rowSums(subSeqData$counts$bacteria) > 0
  funNotMiss <- rowSums(subSeqData$counts$fungi) > 0
  # Clean bacteria
  bacClean <- subSeqData$counts$bacteria %>%
    # remove missing
    filter(bacNotMiss) %>%
    # calculate relative abundance
    apply(
      2,
      function(x){
        x / sum(x)
      }
    ) %>%
    # transpose and make data frame
    t %>%
    as.data.frame
  # Clean fungi
  funClean <- subSeqData$counts$fungi %>%
    # remove missing
    filter(funNotMiss) %>%
    # calculate relative abundance
    apply(
      2,
      function(x){
        x / sum(x)
      }
    ) %>%
    # transpose and make data frame
    t %>%
    as.data.frame
  
  ## Format taxonomy data ----
  bacTaxPres <- filter(subSeqData$taxonomy$bacteria, bacNotMiss)
  funTaxPres <- filter(subSeqData$taxonomy$fungi, funNotMiss)
  
  ## Format metadata ----
  # Join datasets
  allMeta <- left_join(
    subSeqData$metadata,
    climData,
    by = c("site", "treat_c")
  )
  # Calculate site level metadata
  siteMeta <- allMeta %>%
    # group
    group_by(site) %>%
    # get site level information
    summarise(
      # difference between high and low elevation sites
      elev_range = diff(range(elev_m)),
      # all other variables take first (same across all elements)
      continent = first(continent),
      country = first(country),
      lat = first(lat),
      lon = first(long),
      mat = first(mat),
      map = first(map),
      age = first(year_range)
    ) %>%
    ungroup
  
  ## Format soil data ----
  soil <- soilData %>%
    # make rep column a character
    mutate(rep_block = as.character(rep)) %>%
    # join to meta data for order
    left_join(
      allMeta,
      ., 
      by = c("site", "treat_a" = "treatment", "rep_block")
    ) %>%
    # select only process rates
    select(h2o_mgg:TRate_d)
  
  ## Build output ----
  # microbes
  microbes <- list(
    counts = list(
      bacteria = bacClean,
      fungi = funClean
    ),
    taxonomy = list(
      bacteria = bacTaxPres,
      fungi = funTaxPres
    )
  )
  # combine
  out <- list(
    metadata = list(
      plots = allMeta,
      sites = siteMeta
    ),
    microbes = microbes,
    processes = soil
  )
  return(out)
}