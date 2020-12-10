################################################################################
#### Project: TP Network
#### Title:   Function | Wrangle | Compile site metadata
#### Author:  Tom Walker (thomas.walker@usys.ethz.ch)
#### Date:    9 December 2020
#### ---------------------------------------------------------------------------

compile_site_data <- function(metadata, climData){
  ## Format metadata ----
  # Join datasets
  allMeta <- left_join(
    metadata,
    climData,
    by = c("site", "treat_c")
  )
  ## Format site-level data ----
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
  ## Return ----
  return(siteMeta)
}