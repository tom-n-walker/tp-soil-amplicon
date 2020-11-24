################################################################################
#### Project: TP Network
#### Title:   Function | Wrangle | Format sequence data
#### Author:  Tom Walker (thomas.walker@usys.ethz.ch)
#### Date:    24 November 2020
#### ---------------------------------------------------------------------------

clean_sequences <- function(subData){
  ## Index for non-missing taxa ----
  bacNotMiss <- rowSums(subData$counts$bacteria) > 0
  funNotMiss <- rowSums(subData$counts$fungi) > 0
  ## Format count data ----
  # Bacteria
  bacClean <- subData$counts$bacteria %>%
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
  # Fungi
  funClean <- subData$counts$fungi %>%
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
  bacTaxPres <- filter(subData$taxonomy$bacteria, bacNotMiss)
  funTaxPres <- filter(subData$taxonomy$fungi, funNotMiss)
  ## Build output ----
  out <- list(
    metadata = subData$metadata,
    counts = list(
      bacteria = bacClean,
      fungi = funClean
    ),
    taxonomy = list(
      bacteria = bacTaxPres,
      fungi = funTaxPres
    )
  )
  return(out)
}