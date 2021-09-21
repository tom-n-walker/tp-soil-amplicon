################################################################################
#### Project: TP Network
#### Title:   Function | Wrangle | Load soil data
#### Author:  Tom Walker (thomas.walker@usys.ethz.ch)
#### Date:    24 November 2020
#### ---------------------------------------------------------------------------

load_soil_data <- function(){
  # load basic soil data
  basic <- fread(
    file = "./data/collated_soil_data.csv",
    data.table = F
  )
  # load drift data
  drift <- fread(
    file = "./data/drift_data.csv",
    data.table = F
  )
  # compile
  out <- drift %>%
    select(site, treatment = treat, rep = tw_rep, rpa_spom:rpa_mom) %>%
    left_join(basic, ., by = c("site", "treatment", "rep"))
  # return
  return(out)
}
