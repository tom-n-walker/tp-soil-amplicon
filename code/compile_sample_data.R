################################################################################
#### Project: TP Network
#### Title:   Function | Wrangle | Compile sample data
#### Author:  Tom Walker (thomas.walker@usys.ethz.ch)
#### Date:    9 December 2020
#### ---------------------------------------------------------------------------

compile_sample_data <- function(metadata, climData, soilData){
  ## Format metadata ----
  # Join datasets
  allMeta <- left_join(
    metadata,
    climData,
    by = c("site", "treat_c")
  )
  ## Add soil data ----
  finished <- soilData %>%
    # make rep column a character
    mutate(rep_block = as.character(rep)) %>%
    # join to meta data for order
    left_join(
      allMeta,
      ., 
      by = c("site", "treat_a" = "treatment", "rep_block")
    ) %>%
    select(sample_code, site:map, h2o_mgg:TRate_d)
  ## Return ----
  return(finished)
}