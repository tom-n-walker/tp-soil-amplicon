################################################################################
#### Project: TP Network
#### Title:   Analysis | Main patterns and trends
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
library(mice)
library(nlme)
library(emmeans)

## Scripts ----
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

## Contributing variables ----
# plot colour scale
myCols <- c("#666699", "#D9A464", "#996666")


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
  group_by(treatment) %>% 
  nest
nestedNums <- finalNums %>% 
  mutate(treatment = finalCats$treatment) %>% 
  group_by(treatment) %>% 
  nest
nestData <- left_join(nestedCats, nestedNums, by = "treatment") %>%
  rename(cats = data.x, nums = data.y)


#### PCA -----------------------------------------------------------------------

## Do PCAs ----
# full data
fullPCA <- prcomp(finalNums, scale = T, center = T)
fullScores <- bind_cols(finalCats, as.data.frame(fullPCA$x[, 1:2]))
fullSummary <- fullScores %>% 
  group_by(treatment) %>%
  summarise(PC1mean = mean(PC1), PC1se = sd(PC1)/sqrt(n()),
            PC2mean = mean(PC2), PC2se = sd(PC2)/sqrt(n()))
fullLoadings <- fullPCA$rotation[, 1:2] %>%
  as.data.frame %>%
  rownames_to_column("variable") %>%
  mutate(plotName = varNames) %>%
  mutate(plotName = fct_reorder(plotName, PC1, min))
fullPlotData <- make_biplot_data(
  fullScores$PC1, 
  fullScores$PC2, 
  fullScores$treatment
)
# nested data
nestData <- nestData %>%
  mutate(pcas = map(nums, ~prcomp(.x, scale = T, center = T))) %>%
  mutate(scores = map(pcas, ~.x$x[, 1:2])) %>%
  mutate(loadings = map(pcas, ~.x$rotation[, 1:2]))

## Plot full PCA scores ----
# biplot
fullPlot <- ggplot(fullPlotData) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_colour_manual(values = myCols) +
  guides(col = "none") +
  aes(x = PC1, y = PC2, col = groupVar, xend = meanPC1, yend = meanPC2) +
  geom_segment() +
  geom_point() +
  geom_point(aes(x = meanPC1, y = meanPC2), shape = 21, fill = "white", size = 3)
# separate axes
fullPC1 <- ggplot(fullSummary) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_fill_manual(values = myCols) +
  guides(fill = "none") +
  aes(x = treatment, y = PC1mean, fill = treatment, 
      ymax = PC1mean + PC1se, ymin = PC1mean - PC1se) +
  geom_hline(yintercept = 0) +
  geom_errorbar(width = 0.2) +
  geom_bar(stat = "identity", col = "black") +
  labs(x = "", y = "PC1 score")
fullPC2 <- ggplot(fullSummary) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_fill_manual(values = myCols) +
  scale_y_continuous(limits = c(-1, 1)) +
  guides(fill = "none") +
  aes(x = treatment, y = PC2mean, fill = treatment, 
      ymax = PC2mean + PC2se, ymin = PC2mean - PC2se) +
  geom_hline(yintercept = 0) +
  geom_errorbar(width = 0.2) +
  geom_bar(stat = "identity", col = "black") +
  labs(x = "", y = "PC2 score")
# loadings
fullPC1loads <- ggplot(fullLoadings) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  coord_flip() +
  scale_y_continuous(limits = c(-0.5, 0.5), breaks = c(-0.5, 0, 0.5)) +
  aes(x = plotName, y = PC1, yend = 0, xend = plotName) +
  geom_hline(yintercept = 0) +
  geom_segment() +
  geom_point(shape = 21, fill = "white", size = 2) +
  labs(x = "", y = "PC1 loading")
fullPC2loads <- ggplot(fullLoadings) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  coord_flip() +
  scale_y_continuous(limits = c(-0.5, 0.5), breaks = c(-0.5, 0, 0.5)) +
  aes(x = plotName, y = PC2, yend = 0, xend = plotName) +
  geom_hline(yintercept = 0) +
  geom_segment() +
  geom_point(shape = 21, fill = "white", size = 2) +
  labs(x = "", y = "PC2 loading")
# print to file
postscript(
  file = "./plots/full_pca_scores.eps", 
  width = mm2in(150), 
  height = mm2in(75)
)
cowplot::plot_grid(
  fullPlot, fullPC1, fullPC2, 
  nrow = 1, 
  rel_widths = c(2, 1, 1),
  align = "hv"
)
dev.off()
postscript(
  file = "./plots/full_pca_loadings.eps", 
  width = mm2in(180), 
  height = mm2in(75)
)
cowplot::plot_grid(
  fullPC1loads, fullPC2loads, 
  nrow = 1, 
  align = "hv"
)
dev.off()

## Analyse ----
# Full PCA PERMANOVA
fullPCAmodel <- vegan::adonis(finalNums ~ treatment, finalCats)
fullPCAmodel
# PC scores
fullPC1model <- lme(PC1 ~ treatment, random = ~1 | site, fullScores, method = "ML")
fullPC2model <- lme(PC2 ~ treatment, random = ~1 | site, fullScores, method = "ML")
anova(fullPC1model, update(fullPC1model, ~.- treatment))
anova(fullPC2model, update(fullPC2model, ~.- treatment))
pairs(emmeans(fullPC1model, "treatment"))
pairs(emmeans(fullPC2model, "treatment"))



#### COUPLING ------------------------------------------------------------------

## Baseline plots ----
soilCr <- ggplot() +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  aes(x = soilC_mgCg, y = R_ugCgh, fill = T_sum_cum) +
  guides(fill = "none") +
  scale_fill_viridis_b() +
  scale_x_continuous(limits = c(10, 240)) +
  scale_y_continuous(limits = c(0, 13)) +
  geom_point(shape = 21, size = 1) +
  geom_smooth(col = "black", size = 0.5, se = F) +
  labs(
    x = expression(paste("Soil C (mg ", g^{-1}, ")")),
    y = expression(paste("R (µg C ", g^{-1}, " ", h^{-1}, ")"))
  )
soilCmicC <- ggplot() +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  aes(x = soilC_mgCg, y = Cmic_ugCg, fill = T_sum_cum) +
  guides(fill = "none") +
  scale_fill_viridis_b() +
  scale_x_continuous(limits = c(10, 240)) +
  scale_y_continuous(limits = c(20, 5500)) +
  geom_point(shape = 21, size = 1) +
  geom_smooth(col = "black", size = 0.5, se = F) +
  labs(
    x = expression(paste("Soil C (mg ", g^{-1}, ")")),
    y = expression(paste(C[mic], " (mg ", g^{-1}, ")"))
  )
micCr <- ggplot() +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  aes(x = Cmic_ugCg, y = R_ugCgh, fill = T_sum_cum) +
  guides(fill = "none") +
  scale_fill_viridis_b() +
  scale_x_continuous(limits = c(20, 5500)) +
  scale_y_continuous(limits = c(0, 13)) +
  geom_point(shape = 21, size = 1) +
  geom_smooth(col = "black", size = 0.5, se = F) +
  labs(
    x = expression(paste(C[mic], " (mg ", g^{-1}, ")")),
    y = expression(paste("R (µg C ", g^{-1}, " ", h^{-1}, ")"))
  )

## Add full data ----
scatterFull1 <- soilCr %+% allData
scatterFull2 <- soilCmicC %+% allData
scatterFull3 <- micCr %+% allData
cowplot::plot_grid(scatterFull1, scatterFull2, scatterFull3, nrow = 1)

## Add subsetted data ----
nestData <- nestData %>%
  mutate(all_data = map2(cats, nums, bind_cols)) %>%
  mutate(soilCr_plots = map(all_data, ~soilCr %+% .x)) %>%
  mutate(soilCmicC_plots = map(all_data, ~soilCmicC %+% .x)) %>%
  mutate(micCr_plots = map(all_data, ~micCr %+% .x))

## Combine ----
postscript(
  file = "./plots/coupling.eps", 
  width = mm2in(160), 
  height = mm2in(240)
)
cowplot::plot_grid(
  scatterFull1, scatterFull2, scatterFull3,
  nestData$soilCr_plots[[1]], nestData$soilCmicC_plots[[1]], nestData$micCr_plots[[1]],
  nestData$soilCr_plots[[2]], nestData$soilCmicC_plots[[2]], nestData$micCr_plots[[2]],
  nestData$soilCr_plots[[3]], nestData$soilCmicC_plots[[3]], nestData$micCr_plots[[3]],
  nrow = 4,
  align = "hv"
)
dev.off()
