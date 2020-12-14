################################################################################
#### Project: TP Network
#### Title:   Function | Wrangle | Load experiment information
#### Author:  Tom Walker (thomas.walker@usys.ethz.ch)
#### Date:    24 November 2020
#### ---------------------------------------------------------------------------

load_clim_data <- function(){
  # Load data
  all <- fread(
    file = paste0(
      "/Users/tomwalker/Dropbox/projects/2018_transplant/SNF_network/data",
      "/soil_amplicon/metadata_transplant.csv"
    )
  )
  # Format
  out <- all %>%
    filter(elev_cat != "low_2") %>%
    mutate(treat_c = ifelse(elev_cat == "high", "high_site", "low_site")) %>%
    select(
      site = site_code,
      continent:treat_c
    )
  # Return
  return(out)
}
