################################################################################
#### Project: TP Network
#### Title:   Function | Wrangle | Load soil data
#### Author:  Tom Walker (thomas.walker@usys.ethz.ch)
#### Date:    24 November 2020
#### ---------------------------------------------------------------------------

load_soil_data <- function(){
  out <- fread(
    file = paste0(
      "/Users/tomwalker/Dropbox/projects/2018_transplant/SNF_network/data/",
      "soil_vienna/tp_collated_data.csv"
    ),
    data.table = F
  )
  return(out)
}
