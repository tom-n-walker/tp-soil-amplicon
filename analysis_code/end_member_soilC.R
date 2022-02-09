################################################################################
#### Project: TP Network Soil
#### Title:   Main effect - soil carbon content
#### Author:  Tom Walker (thomas.walker@usys.ethz.ch)
#### Date:    7 July 2021
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
library(nlme)
library(emmeans)
library(cowplot)
library(tidyverse)
source("./processing_code/mm2in.R")

## Data ----
# load
allData <- drake::readd(finalDF)
siteData <- drake::readd(siteData)
# colour palettes
myCols <- wesanderson::wes_palette("GrandBudapest1", n = 3)


#### FORMAT --------------------------------------------------------------------

## Assemble data on background context ----
# get meaningful context
drivers <- allData$metaMetrics$soilC_mgCg %>%
  select(site, year_range, elev_range,
         mat_mean, mst_mean, tap_mean, 
         mat_diff, mst_diff, tap_diff, 
         mat_cumdiff, mst_cumdiff, tap_cumdiff) %>%
  mutate(soilChh = allData$metaMetrics$soilC_mgCg$HH)

## Assemble responses from metrics list ----
# hl to hh
HLvHH <- lapply(allData$metaMetrics, function(x){x$HLvHHmean}) %>%
  do.call(cbind, .) %>%
  as.data.frame %>%
  bind_cols(drivers, .) %>%
  # remove where no soil C data
  filter(!is.na(soilChh))
# hl from ll
HLvLL <- lapply(allData$metaMetrics, function(x){x$HLvLLmean}) %>%
  do.call(cbind, .) %>%
  as.data.frame %>%
  bind_cols(drivers, .) %>%
  # remove where no soil C data
  filter(!is.na(soilChh))
# ll from hh
LLvHH <- lapply(allData$metaMetrics, function(x){x$HHvLLmean}) %>%
  do.call(cbind, .) %>%
  as.data.frame %>%
  bind_cols(drivers, .) %>%
  # remove where no soil C data
  filter(!is.na(soilChh))

# add mineral/organic soil split
HLvHH[which(HLvHH$soilChh <= 100), "soilCcat"] <- "mineral"
HLvHH[which(HLvHH$soilChh > 100), "soilCcat"] <- "organic"
HLvHH$soilCcat <- factor(HLvHH$soilCcat, levels = c("mineral", "organic"))
HLvLL$soilCcat <- LLvHH$soilCcat <- HLvHH$soilCcat


#### MODEL SOIL C --------------------------------------------------------------

# create weights for modelling
varSite <- varIdent(form = ~ 1 | site)

## HL versus HH ----
# model fit
m1 <- gls(
  soilC_mgCg ~ mst_cumdiff, 
  weights = varSite,
  HLvHH, 
  na.action = "na.exclude", 
  method = "ML"
)
# diagnose
r1 <- residuals(m1, type = "pearson")
par(mfrow = c(1, 3))
plot(r1 ~ fitted(m1))
boxplot(r1 ~ HLvHH$site)
hist(r1)
# effects
drop1(m1, test = "Chisq")
m1a <- update(m1, ~.- mst_cumdiff:tap_mean)
drop1(m1a, test = "Chisq")








## From HH control ----
# model fit
m1 <- gls(
  soilC_mgCg ~ soilChh, 
  weights = varSite,
  LLvHH, 
  na.action = "na.exclude", 
  method = "ML"
)
# diagnose
r1 <- residuals(m1, type = "pearson")
par(mfrow = c(1, 3))
plot(r1 ~ fitted(m1))
boxplot(r1 ~ HLvHH$site)
hist(r1)
# effects
drop1(m1, test = "Chisq")
m1a <- update(m1, ~.- mst_cumdiff:tap_mean)
drop1(m1a, test = "Chisq")

## To LL control ----
# model fit
m2 <- gls(
  soilC_mgCg ~ mst_cumdiff * soilChh + mst_cumdiff * tap_mean, 
  weights = varSite,
  HLvLL, 
  na.action = "na.exclude", 
  method = "ML"
)
# diagnose
r2 <- residuals(m2, type = "pearson")
par(mfrow = c(1, 3))
plot(r2 ~ fitted(m2))
boxplot(r2 ~ HLvLL$site)
hist(r2)
# effects
drop1(m2, test = "Chisq")
m2a <- update(m2, ~.- mst_cumdiff:tap_mean - mst_cumdiff:soilChh)
drop1(m2a, test = "Chisq")


#### PLOT ----------------------------------------------------------------------

## From HH ----
fromHH <- ggplot(HLvHH) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  guides(fill = "none", col = "none") +
  aes(x = mst_cumdiff, y = soilC_mgCg, col = soilCcat) +
  scale_x_continuous(limits = c(0, NA)) +
  scale_colour_manual(values = myCols) +
  geom_hline(yintercept = 0) +
  geom_point(size = 1.25) +
  geom_smooth(method = "lm", se = F) +
  xlab("Cumulative warming (ºC)") +
  ylab("Diff. from origin (%)")

## To LL ----
toLL <- ggplot(HLvLL) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  guides(fill = "none", col = "none") +
  aes(x = mst_cumdiff, y = soilC_mgCg, col = soilCcat) +
  scale_x_continuous(limits = c(0, NA)) +
  scale_colour_manual(values = myCols) +
  geom_hline(yintercept = 0) +
  geom_point(col = "black", size = 1.25) +
  geom_smooth(method = "lm", se = F) +
  geom_vline(xintercept = 70.8/2.35) +
  xlab("Cumulative warming (ºC)") +
  ylab("Diff. from destination (%)")

## Write ----
postscript(
  file = "./figure_builds/soilC_endmembers.eps", 
  width = mm2in(120),
  height = mm2in(60)
)
plot_grid(
  fromHH, toLL, 
  nrow = 1,
  align = "hv", axis = "tlbr"
  )
dev.off()
