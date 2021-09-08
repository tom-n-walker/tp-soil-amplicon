################################################################################
#### Project: TP Network
#### Title:   Function | Wrangle | Load experiment information
#### Author:  Tom Walker (thomas.walker@usys.ethz.ch)
#### Date:    24 November 2020
#### ---------------------------------------------------------------------------

load_clim_data <- function(){
  # directory path
  dir <- "/Users/tomwalker/Dropbox/projects/2018_transplant/SNF_network/data/soil_vienna/"
  # load climate and metadata
  meta <- fread(paste0(dir, "site_metadata.csv"), data.table = F)
  clim <- fread(paste0(dir, "TPclimatedat.csv"), data.table = F)
  # format
  out <- meta %>%
    left_join(., clim, by = c("gradient" = "Gradient", "site_id" = "destSiteID")) %>%
    mutate(year_range = 2018 - start_year) %>%
    mutate(treat_c = ifelse(elev_cat == "high", "high_site", "low_site")) %>%
    rename(site = site_code) %>%
    select(gradient:site_id, treat_c, elev_cat:start_year, year_range, plot_m2, T_ann_cor:P_ann)
  # return
  return(out)
}
