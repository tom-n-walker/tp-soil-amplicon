################################################################################
#### Project: 
#### Title:   
#### Author:  Tom Walker (thomas.walker@usys.ethz.ch)
#### Date:    
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
library(mice)
library(metagMisc)
library(nlme)

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

## Climate data ----
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

## Process data ----
# generate factor order for plots
varNames <- c(
  # bulk soil pH and C/N pools
  "pH", "Soil C", "Soil N", "Soil C:N", 
  # dissolved C/N pools
  "DOC", "TDN", "NH4", "NO3", "DON",
  # microbial C/N pools
  "Microbial C", "Microbial N", "Microbial C:N", 
  # microbial processes
  "R", "G", "Mass-specific R", "Mass-specific G", "C-specific R", "C-specific G",
  "CUE", "Turnover rate",
  # mcc
  "Fungi 1", "Fungi 2", "Bacteria 1", "Bacteria 2"
)
# metadata
finalCats <- finalDF %>%
  # recode factors for full analysis
  mutate(treat_diff = factor(treat_a, c("HL", "HH", "LL"))) %>%
  mutate(treatment = factor(treat_a, c("HH", "HL", "LL"))) %>%
  mutate(site = factor(site, levels(climSum$site))) %>%
  select(site, treat_diff, treatment, rep_block, continent, country, elev_cat:year_range)
# all data
final <- finalDF %>%
  # create soil C corrected microbial data
  mutate(RC_ugSCh = R_ugCgh / soilC_mgCg * 1000,
         GC_ugSCh = G_ugCgh / soilC_mgCg * 1000) %>%
  # select numbers
  select(pH:CmicNmic, R_ugCgh:Gm_ugCmic, RC_ugSCh, GC_ugSCh, 
         CUE:funFull2, bacFull1:bacFull2) %>%
  # impute missing values
  apply_mice(., 5) %>%
  # bind to categorical data
  bind_cols(finalCats, .)


#### PCAs ----------------------------------------------------------------------

## Do PCAs ----
# group climate data
climGroup <- climSum %>%
  group_by(site)
# Do PCAs
pcas <- final %>%
  # group and setup nesting
  group_by(site) %>%
  nest() %>%
  left_join(climGroup, .) %>%
  # separate numbers and do PCA
  mutate(nums = map(data, ~select(.x, pH:bacFull2))) %>%
  mutate(pcas = map(nums, ~prcomp(.x, scale = T, center = T))) %>%
  # extract loadings and filter for PC1, create factor
  mutate(loadings = map(pcas, ~tidy(.x, matrix = "rotation"))) %>%
  mutate(loadings = map(loadings, ~filter(.x, PC == 1))) %>%
  # make axis directions the same
  mutate(loadings = map(loadings, ~flip_axes(.x))) %>%
  mutate(loadings = map(loadings, ~mutate(.x, variable = factor(varNames, varNames))))
# Unnest loadings
loadings <- pcas %>%
  unnest(loadings) %>%
  ungroup %>%
  mutate(site = fct_reorder(site, yearRange, min))

## Plot PCAs ----
ggplot(loadings) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  coord_flip() +
  facet_wrap(~ site) +
  aes(x = variable, y = value, yend = 0, xend = variable) +
  geom_hline(yintercept = 0) +
  geom_segment() +
  geom_point(shape = 21, fill = "white", size = 2) +
  labs(x = "", y = "PC1 loadings")


#### MODELS --------------------------------------------------------------------

## Format data ----
# Rescale within site
rescaled <- final %>%
  # group and setup nesting
  group_by(site) %>%
  nest %>%
  mutate(rescaled = map(data, ~mutate(.x, across(pH:bacFull2, arm::rescale)))) %>%
  unnest(rescaled)
# calculate mean ± se
summaries <- rescaled %>%
  group_by(site, treatment) %>%
  summarise(
    across(
      .cols = pH:bacFull2,
      .fns = list(
        mean = mean, 
        se = ~sd(.x)/sqrt(length(.x))
      )
    )
  )
# calculate response ratios (raw values OK here)
grouped <- final %>%
  group_by(site) %>%
  nest() %>%
  left_join(climGroup, .) %>%
  mutate(soilCrrs = map(data, ~rrs(.x, "soilC_mgCg"))) %>%
  mutate(respRRs = map(data, ~rrs(.x, "R_ugCgh"))) %>%
  mutate(rCrrs = map(data, ~rrs(.x, "RC_ugSCh"))) %>%
  mutate(rCmicRRs = map(data, ~rrs(.x, "Rm_ugCmic"))) %>%
  mutate(cmicRRs = map(data, ~rrs(.x, "Cmic_ugCg")))

## Plot ----
# soil C
soilCrrHH <- grouped %>%
  unnest(soilCrrs) %>%
  ggplot() +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  aes(x = sumCumWarm, y = hlhh_mean, ymin = hlhh_mean - hlhh_se, ymax = hlhh_mean + hlhh_se) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(width = 0.5) +
  geom_point(shape = 21, fill = "#D9A464") +
  geom_smooth(method = "lm", se = F, col = "#D9A464", size = 0.5) +
  labs(x = "Cum. warming (ºC)", y = "Difference from HH")
soilCrrLL <- grouped %>%
  unnest(soilCrrs) %>%
  ggplot() +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  aes(x = sumCumWarm, y = hlll_mean, ymin = hlll_mean - hlll_se, ymax = hlll_mean + hlll_se) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(width = 0.5) +
  geom_point(shape = 21, fill = "#D9A464") +
  geom_smooth(method = "lm", se = F, col = "#D9A464", size = 0.5) +
  labs(x = "Cum. warming (ºC)", y = "Difference from LL")
cowplot::plot_grid(soilCrrHH, soilCrrLL)
# respiration
resprrHH <- grouped %>%
  unnest(respRRs) %>%
  ggplot() +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  aes(x = sumCumWarm, y = hlhh_mean, ymin = hlhh_mean - hlhh_se, ymax = hlhh_mean + hlhh_se) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(width = 0.5) +
  geom_point(shape = 21, fill = "#D9A464") +
  geom_smooth(method = "lm", se = F, col = "#D9A464", size = 0.5) +
  labs(x = "Cum. warming (ºC)", y = "Difference from HH")
resprrLL <- grouped %>%
  unnest(respRRs) %>%
  ggplot() +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  aes(x = sumCumWarm, y = hlll_mean, ymin = hlll_mean - hlll_se, ymax = hlll_mean + hlll_se) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(width = 0.5) +
  geom_point(shape = 21, fill = "#D9A464") +
  geom_smooth(method = "lm", se = F, col = "#D9A464", size = 0.5) +
  labs(x = "Cum. warming (ºC)", y = "Difference from LL")
cowplot::plot_grid(resprrHH, resprrLL)
# Cmic
micCrrHH <- grouped %>%
  unnest(cmicRRs) %>%
  ggplot() +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  aes(x = sumCumWarm, y = hlhh_mean, ymin = hlhh_mean - hlhh_se, ymax = hlhh_mean + hlhh_se) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(width = 0.5) +
  geom_point(shape = 21, fill = "#D9A464") +
  geom_smooth(method = "lm", se = F, col = "#D9A464", size = 0.5) +
  labs(x = "Cum. warming (ºC)", y = "Difference from HH")
micCrrLL <- grouped %>%
  unnest(cmicRRs) %>%
  ggplot() +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  aes(x = sumCumWarm, y = hlll_mean, ymin = hlll_mean - hlll_se, ymax = hlll_mean + hlll_se) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(width = 0.5) +
  geom_point(shape = 21, fill = "#D9A464") +
  geom_smooth(method = "lm", se = F, col = "#D9A464", size = 0.5) +
  labs(x = "Cum. warming (ºC)", y = "Difference from LL")
cowplot::plot_grid(micCrrHH, micCrrLL)

# plot to file
postscript(file = "./plots/hh_divergence.eps", width = mm2in(180), height = mm2in(60))
cowplot::plot_grid(resprrHH, soilCrrHH, micCrrHH, nrow = 1)
dev.off()


## Plot data ----
# Soil C
soilCplot <- ggplot(summaries) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  guides(fill = "none") +
  aes(
    x = site, 
    y = soilC_mgCg_mean, 
    ymin = soilC_mgCg_mean - soilC_mgCg_se,
    ymax = soilC_mgCg_mean + soilC_mgCg_se,
    fill = treatment
  ) +
  scale_fill_manual(values = c("#666699", "#D9A464", "#996666")) +
  geom_errorbar(width = 0.25, position = position_dodge(0.5)) +
  geom_point(position = position_dodge(0.5), shape = 21) +
  labs(x = "", y = expression(paste("Soil C (mg ", g^{-1}, ")")))
# Respiration
rplot <- ggplot(summaries) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  guides(fill = "none") +
  aes(
    x = site, 
    y = R_ugCgh_mean, 
    ymin = R_ugCgh_mean - R_ugCgh_se,
    ymax = R_ugCgh_mean + R_ugCgh_se,
    fill = treatment
  ) +
  scale_fill_manual(values = c("#666699", "#D9A464", "#996666")) +
  geom_errorbar(width = 0.25, position = position_dodge(0.5)) +
  geom_point(position = position_dodge(0.5), shape = 21) +
  labs(x = "", y = expression(paste("R (µg ", g^{-1}, " ", h^{-1}, ")")))
# Respiration (per soil C)
cmicPlot <- ggplot(summaries) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  guides(fill = "none") +
  aes(
    x = site, 
    y = Cmic_ugCg_mean, 
    ymin = Cmic_ugCg_mean - Cmic_ugCg_se,
    ymax = Cmic_ugCg_mean + Cmic_ugCg_se,
    fill = treatment
  ) +
  scale_fill_manual(values = c("#666699", "#D9A464", "#996666")) +
  geom_errorbar(width = 0.25, position = position_dodge(0.5)) +
  geom_point(position = position_dodge(0.5), shape = 21) +
  labs(x = "", y = expression(paste(C[mic], " (µg C ", g^{-1}, ")")))
# plot to file
postscript(file = "./plots/cmic.eps", width = mm2in(150), height = mm2in(60))
print(cmicPlot)
dev.off()
postscript(file = "./plots/resp.eps", width = mm2in(150), height = mm2in(60))
print(rplot)
dev.off()
postscript(file = "./plots/soilC.eps", width = mm2in(150), height = mm2in(60))
print(soilCplot)
dev.off()



## Prepare data frame ----
nested <- final %>%
  # group and setup nesting
  group_by(site) %>%
  nest() %>%
  left_join(climSum, .)

## Do all models ----
# Using own function
soilCmodels <- nested_lms(nested, "soilC_mgCg")
micCmodels <- nested_lms(nested, "Cmic_ugCg")
rModels <- nested_lms(nested, "R_ugCgh")
rMassModels <- nested_lms(nested, "Rm_ugCmic")
gMassModels <- nested_lms(nested, "Gm_ugCmic")
rSoilCModels <- nested_lms(nested, "RC_ugSCh")
gSoilCModels <- nested_lms(nested, "GC_ugSCh")
gModels <- nested_lms(nested, "G_ugCgh")
cueModels <- nested_lms(nested, "CUE")
fun1Models <- nested_lms(nested, "funFull1")
fun2Models <- nested_lms(nested, "funFull2")
bac1Models <- nested_lms(nested, "bacFull1")
bac2Models <- nested_lms(nested, "bacFull2")






# Plot coefficients
test <- micCmodels %>%
  unnest(llhl_coefs) %>%
  ungroup
anova(lm(estimate ~ sumWarm * yearRange, test))

ggplot(test) +
  aes(x = sumCumWarm, y = estimate) +
  geom_point() +
  geom_smooth(method = "lm")










