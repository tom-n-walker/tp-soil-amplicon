################################################################################
#### Project: TP Network Soil
#### Title:   Is warming confounded with design
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
allData <- drake::readd(siteData)
plotData <- drake::readd(finalDF)$fullDF

# subset
contVars <- allData %>%
  select(year_range:tap_cumdiff) %>%
  as.matrix %>%
  cor %>%
  abs
# make labels
labels <- c(
  "Duration", "Plot size", 
  "MAT", "MST", "TAP",
  "δ ASL", "δ MAT", "δ MST", "δ TAP",
  "Σ MAT", "Σ MST", "Σ TAP"
)
colnames(contVars) <- rownames(contVars) <- labels


#### ELEVATION -----------------------------------------------------------------

## Format data ----
stData <- plotData %>%
  select(site, elev_cat, T_ann_cor:P_ann) %>%
  mutate(across(T_ann_cor:P_ann, arm::rescale))

## Models ----
# build
m1 <- lme(
  T_ann_cor ~ elev_cat,
  random = ~ 1 | site,
  method = "ML",
  data = stData,
  na.action = "na.exclude"
)
m2 <- lme(
  T_sum_cor ~ elev_cat,
  random = ~ 1 | site,
  method = "ML",
  data = stData,
  na.action = "na.exclude"
)
m3 <- lme(
  P_ann ~ elev_cat,
  random = ~ 1 | site,
  method = "ML",
  data = stData,
  na.action = "na.exclude"
)
# coefficients
m1Coef <- intervals(m1)$fixed[2, ]
m2Coef <- intervals(m2)$fixed[2, ]
m3Coef <- intervals(m3)$fixed[2, ]
allCoefs <- data.frame(
  response = c("MAT", "MST", "TAP"),
  lower = c(m1Coef[1], m2Coef[1], m3Coef[1]),
  mean = c(m1Coef[2], m2Coef[2], m3Coef[2]),
  upper = c(m1Coef[3], m2Coef[3], m3Coef[3])
)

## MAT ----
m1 <- lme(
  plot_m2 ~ elev_cat,
  random = ~ 1 | site,
  method = "ML",
  data = plotData,
  na.action = "na.exclude"
)
anova(m1, update(m1, ~.- elev_cat))

## MST ----
m2 <- lme(
  T_sum_cor ~ elev_cat,
  random = ~ 1 | site,
  method = "ML",
  data = plotData,
  na.action = "na.exclude"
)
anova(m2, update(m2, ~.- elev_cat))

## TAP ----
m3 <- lme(
  P_ann ~ elev_cat,
  random = ~ 1 | site,
  method = "ML",
  data = plotData,
  na.action = "na.exclude"
)
anova(m3, update(m3, ~.- elev_cat))


#### ANALYSE -------------------------------------------------------------------

## Heatmap of all context variables ----
# filter for correlations > 0.6
contVars[contVars < 0.65] <- 0
# make heatmap
postscript(
  file = "./figure_builds/design_heatmap.eps",
  width = mm2in(180),
  height = mm2in(180)
)
gplots::heatmap.2(
  x = contVars,
  trace = "none",
  density.info = "none",
  col = colorRampPalette(colors = c("white", "#313639"))(100)
)
dev.off()
