################################################################################
#### Project: TP Network
#### Title:   Analysis - mcc
#### Author:  Tom Walker (thomas.walker@usys.ethz.ch)
#### Date:    23 February 2021
#### ---------------------------------------------------------------------------


#### PROLOGUE ------------------------------------------------------------------

## Options ----
# remove objects from global environment
rm(list = ls())
# R session options (no factors, bias against scientific #s)
options(
stringsAsFactors = F,
scipen = 6
)

## Libraries ----
# standard library set
library(tidyverse)
library(data.table)

## Code ----
scripts <- list.files(
  "./code",
  full.names = T
)

sapply(
  scripts,
  source
)

## Data ----
# drake processed data
drake::loadd()
# climate data
clim <- fread(
  paste0(
    "../../../Dropbox/projects/2018_transplant/SNF_network/", 
    "data/climate_data/worlclim2_processedtemp.csv"
  ),
  data.table = F
)


#### FORMAT --------------------------------------------------------------------

## Site-level climate data ----
climSum <- clim %>%
  # arrange data to make diff calculations work
  arrange(site, elev_cat) %>%
  # bin useless columns
  select(site, elev_cat:elev, T_ann_cor:T_win_cor) %>%
  # bin china
  filter(site != "YAN") %>%
  # group by site
  group_by(site) %>%
  summarise(
    .groups = "keep",
    # calculate warming between high and low sites
    annWarm = diff(T_ann_cor),
    sumWarm = diff(T_sum_cor),
    winWarm = diff(T_win_cor),
    # re-express in terms of years of warming
    annCumWarm = annWarm * first(year_range),
    sumCumWarm = sumWarm * first(year_range),
    winCumWarm = winWarm * first(year_range),
    # calculate elevation and year ranges
    elevRange = diff(rev(elev)),
    yearRange = first(year_range)
  ) %>%
  ungroup %>%
  mutate(site = fct_reorder(site, sumCumWarm, min))

## Variable names ----
varNames <- c(
  # bulk soil pH and C/N pools
  "pH", "Soil C", "Soil N", "Soil C:N", 
  # dissolved C/N pools
  "DOC", "TDN", "NH4", "NO3", "DON",
  # microbial C/N pools
  "Microbial C", "Microbial N", "Microbial C:N", 
  # microbial processes
  "R", "G", "Mass-specific R", "Mass-specific G", 
  "C-specific R", "C-specific G", "CUE", "Turnover rate",
  # mcc
  "Fungi 1", "Fungi 2", "Bacteria 1", "Bacteria 2"
)

## Metadata ----
finalCats <- finalDF %>%
  # recode factors for full analysis
  mutate(treatment = factor(treat_a, c("HH", "HL", "LL"))) %>%
  mutate(site = factor(site, levels(climSum$site))) %>%
  # select some variables before joining to avoid duplicates
  select(site, continent, country, elev_cat, elev_m, treatment) %>%
  # join and subset variables again
  left_join(., clim, by = c("site", "elev_cat")) %>%
  select(site:elev_m, lat, lon, year_range, T_ann_cor:T_win_cum, treatment)

## Numeric data ----
finalNums <- finalDF %>%
  # create soil C corrected microbial data
  mutate(RC_ugSCh = R_ugCgh / soilC_mgCg * 1000,
         GC_ugSCh = G_ugCgh / soilC_mgCg * 1000) %>%
  # select numbers
  select(pH:CmicNmic, R_ugCgh:Gm_ugCmic, RC_ugSCh, GC_ugSCh, 
         CUE:funFull2, bacFull1:bacFull2) %>%
  # impute missing values
  apply_mice(., 5)

## Collate data ----
# classic format
allData <- bind_cols(finalCats, finalNums)
# nested DF
nestedCats <- finalCats %>% 
  group_by(site) %>% 
  nest
nestedNums <- finalNums %>% 
  mutate(site = finalCats$site) %>% 
  group_by(site) %>% 
  nest
nestData <- nestedCats %>%
  left_join(nestedNums, by = "site") %>%
  left_join(group_by(climSum, site), by = "site") %>%
  rename(cats = data.x, nums = data.y)


#### PLOT ----------------------------------------------------------------------

## Baseline plot ----
baseplot <- ggplot() +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_colour_manual(values = c("#666699", "#D9A464", "#996666")) +
  guides(col = "none") +
  aes(x = PC1, y = PC2, col = groupVar, xend = meanPC1, yend = meanPC2) +
  geom_segment() +
  geom_point() +
  geom_point(aes(x = meanPC1, y = meanPC2), shape = 21, fill = "white", size = 3) +
  labs(x = "Axis 1", y = "Axis 2")

## Make biplot data ----
biplotData <- nestData %>%
  mutate(all_data = map2(cats, nums, bind_cols)) %>%
  mutate(fun_plotdata = map(all_data, ~make_biplot_data(.x$funFull1, .x$funFull2, .x$treatment))) %>%
  mutate(bac_plotdata = map(all_data, ~make_biplot_data(.x$bacFull1, .x$bacFull2, .x$treatment))) %>%
  mutate(fun_plots = map(fun_plotdata, ~baseplot %+% .x)) %>%
  mutate(bac_plots = map(bac_plotdata, ~baseplot %+% .x))

## Collate ----
# arrange by warming
outData <- arrange(biplotData, sumCumWarm)
# fungi
postscript(
  file = "./plots/fungi.eps", 
  width = mm2in(180),
  height = mm2in(240)
)
cowplot::plot_grid(
  plotlist = outData$fun_plots,
  nrow = 4,
  align = "hv"
)
dev.off()
# bacteria
postscript(
  file = "./plots/bacteria.eps", 
  width = mm2in(160),
  height = mm2in(240)
)
cowplot::plot_grid(
  plotlist = outData$bac_plots,
  nrow = 4,
  align = "hv"
)
dev.off()



