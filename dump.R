################################################################################
#### Project: TP Network soil
#### Title:   Analyse treatment effects on headline responses
#### Author:  Tom Walker (thomas.walker@usys.ethz.ch)
#### Date:    10 June 2021
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

## Scripts ----
source("./code/apply_mice.R")
source("./code/make_biplot_data.R")
source("./code/mm2in.R")

## Libraries ----
# standard library set
library(tidyverse)
library(data.table)
library(nlme)
library(mice)
library(emmeans)

## Data ----
# load
allData <- drake::readd(finalDF)
siteData <- drake::readd(siteData)

# plot colour scale
myCols <- c("#666699", "#D9A464", "#996666")


#### FORMAT --------------------------------------------------------------------

## Assemble data on background context ----
# get meaningful context
drivers <- allData$metaMetrics$soilC_mgCg %>%
  select(site, year_range, elev_range,
         mat_mean, mst_mean, tap_mean, 
         mat_diff, mst_diff, tap_diff, 
         mat_cumdiff, mst_cumdiff, tap_cumdiff) %>%
  mutate(soilChh = allData$metaMetrics$soilC_mgCg$HH)
# assemble responses
HLvHH <- lapply(allData$metaMetrics, function(x){x$HLvHHmean}) %>%
  do.call(cbind, .) %>%
  as.data.frame %>%
  bind_cols(drivers, .) %>%
  # make categorical soil C based on median of value (for interaction)
  mutate(soilCcat = ifelse(soilChh <= 94.56, "low", "high"))
LLvHH <- lapply(allData$metaMetrics, function(x){x$HHvLLmean}) %>%
  do.call(cbind, .) %>%
  as.data.frame %>%
  bind_cols(drivers, .) %>%
  # make categorical soil C based on median of value (for interaction)
  mutate(soilCcat = ifelse(soilChh <= 94.56, "low", "high"))


#### MODEL ---------------------------------------------------------------------

## Build dataset for predicting ----
predictHLvHH <- data.frame(
  mst_cumdiff = rep(0:30, 2), # create sequence of cumulative warming
  tap_mean = mean(drivers$tap_mean, na.rm = T), # use mean TAP for this
  soilCcat = rep(c("L", "H"), each = 31),
  soilChh = rep(c(68, 122), each = 31) # use 1st and 3rd quantiles for low/high
)

## Models ----
# build
m1 <- gls(
  soilC_mgCg ~ mst_cumdiff * soilChh + mst_cumdiff * tap_mean, 
  weights = varIdent(form = ~ 1 | site),
  HLvHH, 
  na.action = "na.exclude", 
  method = "ML"
  )


m2 <- lm(Cmic_ugCg ~ mst_cumdiff * soilChh, HLvHH)
m3 <- lm(R_ugCgh ~ mst_cumdiff * soilChh, HLvHH)
m4 <- lm(G_ugCgh ~ mst_cumdiff * soilChh, HLvHH)
m5 <- lm(Rm_ugCmic ~ mst_cumdiff * soilChh, HLvHH)
m6 <- lm(Gm_ugCmic ~ mst_cumdiff * soilChh, HLvHH)
m7 <- lm(CUE ~ mst_cumdiff * soilChh, HLvHH)
# test
drop1(m1, test = "Chisq")


anova(m1)
anova(m2)
anova(m2)
anova(m3)
anova(m5)
anova(m6)
anova(m7)
# add predicted
predictHLvHH <- predictHLvHH %>%
  mutate(soilC = predict(m1, newdata = predictHLvHH)) %>%
  mutate(Cmic = predict(m2, newdata = predictHLvHH)) %>%
  mutate(R = predict(m3, newdata = predictHLvHH)) %>%
  mutate(G = predict(m4, newdata = predictHLvHH)) %>%
  mutate(Rmic = predict(m5, newdata = predictHLvHH)) %>%
  mutate(Gmic = predict(m6, newdata = predictHLvHH)) %>%
  mutate(CUE = predict(m7, newdata = predictHLvHH))
# make longer for ease of plotting (only to interpret)
longPredicted <- pivot_longer(predictHLvHH, soilC:CUE)


#### PLOT ----------------------------------------------------------------------

# combined plot
ggplot(longPredicted) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  aes(x = mst_cumdiff, y = value, col = soilCcat) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line() +
  facet_wrap(~ name, scales = "free") +
  labs(x = "Cumulative MST (ºC)", y = "HL from HH (%)")




soilCmodel <- ggplot(predictHLvHH) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  guides(col = "none") +
  aes(x = mst_cumdiff, y = soilC, col = soilCcat) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line() +
  labs(x = "Cumulative MST (ºC)", y = "HL from HH (%)")
cmicModel <- ggplot(predictHLvHH) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  guides(col = "none") +
  aes(x = mst_cumdiff, y = Cmic, col = soilCcat) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line() +
  labs(x = "Cumulative MST (ºC)", y = "HL from HH (%)")
rModel <- ggplot(predictHLvHH) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  guides(col = "none") +
  aes(x = mst_cumdiff, y = R, col = soilCcat) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line() +
  labs(x = "Cumulative MST (ºC)", y = "HL from HH (%)")
gModel <- ggplot(predictHLvHH) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  guides(col = "none") +
  aes(x = mst_cumdiff, y = G, col = soilCcat) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line() +
  labs(x = "Cumulative MST (ºC)", y = "HL from HH (%)")
rmicModel <- ggplot(predictHLvHH) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  guides(col = "none") +
  aes(x = mst_cumdiff, y = Rmic, col = soilCcat) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line() +
  labs(x = "Cumulative MST (ºC)", y = "HL from HH (%)")
gmicModel <- ggplot(predictHLvHH) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  guides(col = "none") +
  aes(x = mst_cumdiff, y = Gmic, col = soilCcat) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line() +
  labs(x = "Cumulative MST (ºC)", y = "HL from HH (%)")
cueModel <- ggplot(predictHLvHH) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  guides(col = "none") +
  aes(x = mst_cumdiff, y = CUE, col = soilCcat) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line() +
  labs(x = "Cumulative MST (ºC)", y = "HL from HH (%)")






ggplot(drop_na(HLvHH)) +
  aes(x = mst_cumdiff, y = R_ugCgh, col = soilCcat) +
  geom_point() +
  geom_smooth(method = "lm")



plot_grid(soilCmodel, cmicModel, rmicModel, gmicModel, cueModel, nrow = 1)














