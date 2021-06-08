################################################################################
#### Project: TP Network
#### Title:   Function | Wrangle | Compile sample data
#### Author:  Tom Walker (thomas.walker@usys.ethz.ch)
#### Date:    9 December 2020
#### ---------------------------------------------------------------------------

compile_sample_data <- function(metadata, climData){
  ## Format metadata ----
  # Join datasets
  allMeta <- left_join(
    metadata,
    climData,
    by = c("site", "treat_c")
  )
  ## Return ----
  return(allMeta)
}
