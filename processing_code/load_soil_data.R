################################################################################
#### Project: TP Network
#### Title:   Function | Wrangle | Load soil data
#### Author:  Tom Walker (thomas.walker@usys.ethz.ch)
#### Date:    24 November 2020
#### ---------------------------------------------------------------------------

load_soil_data <- function(){
  out <- fread(
    file = "./data/collated_soil_data.csv",
    data.table = F
  )
  return(out)
}
