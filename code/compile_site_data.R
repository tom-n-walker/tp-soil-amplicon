################################################################################
#### Project: TP Network
#### Title:   Function | Wrangle | Compile site metadata
#### Author:  Tom Walker (thomas.walker@usys.ethz.ch)
#### Date:    9 December 2020
#### ---------------------------------------------------------------------------

compile_site_data <- function(climData){
  ## Format site-level data ----
  allData <- climData %>%
    # remove yan site (sampling problem led to them being incompatible)
    filter(site != "YAN") %>%
    # arrange data to make diff calculations work
    arrange(site, elev_cat) %>%
    # arrange data to make diff calculations work
    group_by(site) %>%
    # for differences, take range
    summarise(
      # categorical site level take first
      across(c(gradient:country, start_year:plot_m2), first),
      # mean climate information
      mat_mean = mean(T_ann_cor),
      mst_mean = mean(T_sum_cor),
      tap_mean = mean(P_ann),
      vpd_mean = mean(V_ann),
      # ranges (elevation reversed)
      elev_range = diff(range(elev_m)),
      mat_diff = diff(T_ann_cor),
      mst_diff = diff(T_sum_cor),
      tap_diff = diff(P_ann),
      vpd_diff = diff(V_ann),
      # re-express in terms of years of warming
      mat_cumdiff = mat_diff * first(year_range),
      mst_cumdiff = mst_diff * first(year_range),
      tap_cumdiff = tap_diff * first(year_range),
      vpd_cumdiff = vpd_diff * first(year_range)
    ) %>%
    # ungroup, order sites by cumulative summer warming and make data frame
    ungroup %>%
    mutate(site = fct_reorder(site, mst_cumdiff, min)) %>%
    as.data.frame
  ## Return ----
  return(allData)
}
