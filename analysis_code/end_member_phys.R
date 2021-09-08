################################################################################
#### Project: TP Network Soil
#### Title:   Analysis | End member physiology
#### Author:  Tom Walker (thomas.walker@usys.ethz.ch)
#### Date:    12 July 2021
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
library(nlme)
library(emmeans)
library(cowplot)
source("./code/mm2in.R")

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
  mutate(soilChh = allData$metaMetrics$soilC_mgCg$HH,
         soilCll = allData$metaMetrics$soilC_mgCg$LL)

## Assemble responses from metrics list ----
# hl to hh
HLvHH <- lapply(allData$metaMetrics, function(x){x$HLvHHmean}) %>%
  do.call(cbind, .) %>%
  as.data.frame %>%
  bind_cols(drivers, .) %>%
  # remove where no soil C data
  filter(!is.na(soilChh)) %>%
  filter(Rm_ugCmic < 1500)
# hl from ll
HLvLL <- lapply(allData$metaMetrics, function(x){x$HLvLLmean}) %>%
  do.call(cbind, .) %>%
  as.data.frame %>%
  bind_cols(drivers, .) %>%
  # remove where no soil C data
  filter(!is.na(soilChh)) %>%
  filter(Rm_ugCmic < 500)
# add quartile-based split in starting soil C for plotting
HLvHH[which(HLvHH$soilChh <= 58.7), "soilCcat"] <- "low"
HLvHH[which(HLvHH$soilChh > 58.7 & HLvHH$soilChh <= 87.8), "soilCcat"] <- "medium"
HLvHH[which(HLvHH$soilChh > 87.8), "soilCcat"] <- "high"
HLvHH$soilCcat <- factor(HLvHH$soilCcat, levels = c("low", "medium", "high"))
HLvLL$soilCcat <- HLvHH$soilCcat


#### MODEL FROM ORIGIN ---------------------------------------------------------

# create weights for modelling
varSite <- varIdent(form = ~ 1 | site)

## Rmic ----
# model fit
m1 <- gls(
  Rm_ugCmic ~ mst_cumdiff * soilChh, 
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
m1a <- update(m1, ~.- mst_cumdiff:soilChh)
drop1(m1a, test = "Chisq")

## Gmic ----
# model fit
m2 <- gls(
  Gm_ugCmic ~ mst_cumdiff * soilChh, 
  weights = varSite,
  HLvHH, 
  na.action = "na.exclude", 
  method = "ML"
)
# diagnose
r2 <- residuals(m2, type = "pearson")
par(mfrow = c(1, 3))
plot(r2 ~ fitted(m2))
boxplot(r2 ~ HLvHH$site)
hist(r2)
# effects
drop1(m2, test = "Chisq")
m2a <- update(m2, ~.- mst_cumdiff:soilChh)
drop1(m2a, test = "Chisq")

## Rtot ----
# model fit
m3 <- gls(
  R_ugCgh ~ mst_cumdiff * soilChh, 
  weights = varSite,
  HLvHH, 
  na.action = "na.exclude", 
  method = "ML"
)
# diagnose
r3 <- residuals(m3, type = "pearson")
par(mfrow = c(1, 3))
plot(r3 ~ fitted(m3))
boxplot(r3 ~ HLvHH$site)
hist(r3)
# effects
drop1(m3, test = "Chisq")
m3a <- update(m3, ~.- mst_cumdiff:soilChh)
drop1(m3a, test = "Chisq")

## Gtot ----
# model fit
m4 <- gls(
  G_ugCgh ~ mst_cumdiff * soilChh, 
  weights = varSite,
  HLvHH, 
  na.action = "na.exclude", 
  method = "ML"
)
# diagnose
r4 <- residuals(m4, type = "pearson")
par(mfrow = c(1, 3))
plot(r4 ~ fitted(m4))
boxplot(r4 ~ HLvHH$site)
hist(r4)
# effects
drop1(m4, test = "Chisq")
m4a <- update(m4, ~.- mst_cumdiff:soilChh)
drop1(m4a, test = "Chisq")

## CUE ----
# model fit
m5 <- gls(
  CUE ~ mst_cumdiff * soilChh, 
  weights = varSite,
  HLvHH, 
  na.action = "na.exclude", 
  method = "ML"
)
# diagnose
r5 <- residuals(m5, type = "pearson")
par(mfrow = c(1, 3))
plot(r5 ~ fitted(m5))
boxplot(r5 ~ HLvHH$site)
hist(r5)
# effects
drop1(m5, test = "Chisq")
m5a <- update(m5, ~.- mst_cumdiff:soilChh)
drop1(m5a, test = "Chisq")

## Cmic ----
# model fit
m6 <- gls(
  Cmic_ugCg ~ mst_cumdiff * soilChh, 
  weights = varSite,
  HLvHH, 
  na.action = "na.exclude", 
  method = "ML"
)
# diagnose
r6 <- residuals(m6, type = "pearson")
par(mfrow = c(1, 3))
plot(r6 ~ fitted(m6))
boxplot(r6 ~ HLvHH$site)
hist(r6)
# effects
drop1(m6, test = "Chisq")
m6a <- update(m6, ~.- mst_cumdiff:soilChh)
drop1(m6a, test = "Chisq")


#### COEFFICIENTS --------------------------------------------------------------

## Coefficients for separate soil C amounts (quartile split, see above) ----
coef(lm(R_ugCgh ~ mst_cumdiff, filter(HLvHH, soilCcat == "medium")))
coef(lm(R_ugCgh ~ mst_cumdiff, HLvHH))

#### PLOT ----------------------------------------------------------------------

## Build ----
Gtot <- ggplot(HLvHH) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  guides(col = "none") +
  aes(x = mst_cumdiff, y = G_ugCgh, col = soilCcat) +
  scale_x_continuous(limits = c(0, NA)) +
  scale_colour_manual(values = myCols) +
  geom_hline(yintercept = 0) +
  geom_point(size = 1.25) +
  geom_smooth(method = "lm", se = F) +
  xlab("Cumulative warming (ºC)") +
  ylab("Diff. from origin (%)")
Rtot <- ggplot(HLvHH) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  guides(col = "none") +
  aes(x = mst_cumdiff, y = R_ugCgh, col = soilCcat) +
  scale_x_continuous(limits = c(0, NA)) +
  scale_colour_manual(values = myCols) +
  geom_hline(yintercept = 0) +
  geom_point(size = 1.25) +
  geom_smooth(aes(col = NULL), method = "lm", se = F, col = "black") +
  xlab("Cumulative warming (ºC)") +
  ylab("Diff. from origin (%)")
Cmic <- ggplot(HLvHH) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  guides(col = "none") +
  aes(x = mst_cumdiff, y = Cmic_ugCg, col = soilCcat) +
  scale_x_continuous(limits = c(0, NA)) +
  scale_colour_manual(values = myCols) +
  geom_hline(yintercept = 0) +
  geom_point(size = 1.25) +
  geom_smooth(method = "lm", se = F) +
  xlab("Cumulative warming (ºC)") +
  ylab("Diff. from origin (%)")
Gmic <- ggplot(HLvHH) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  guides(col = "none") +
  aes(x = mst_cumdiff, y = Gm_ugCmic, col = soilCcat) +
  scale_x_continuous(limits = c(0, NA)) +
  scale_colour_manual(values = myCols) +
  geom_hline(yintercept = 0) +
  geom_point(size = 1.25) +
  geom_smooth(method = "lm", se = F) +
  xlab("Cumulative warming (ºC)") +
  ylab("Diff. from origin (%)")
Rmic <- ggplot(HLvHH) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  guides(col = "none") +
  aes(x = mst_cumdiff, y = Rm_ugCmic, col = soilCcat) +
  scale_x_continuous(limits = c(0, NA)) +
  scale_colour_manual(values = myCols) +
  geom_hline(yintercept = 0) +
  geom_point(size = 1.25) +
  geom_smooth(aes(col = NULL), method = "lm", se = F, col = "black") +
  xlab("Cumulative warming (ºC)") +
  ylab("Diff. from origin (%)")
cue <- ggplot(HLvHH) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  guides(col = "none") +
  aes(x = mst_cumdiff, y = CUE, col = soilCcat) +
  scale_x_continuous(limits = c(0, NA)) +
  scale_colour_manual(values = myCols) +
  geom_hline(yintercept = 0) +
  geom_point(size = 1.25) +
  xlab("Cumulative warming (ºC)") +
  ylab("Diff. from origin (%)")

## Combine ----
postscript(
  file = "./plots/phys_endmembers.eps",
  width = mm2in(180), height = mm2in(120),
)
plot_grid(
  Rmic, Gmic, cue, 
  Rtot, Gtot, Cmic,
  align = "tlbr",
  nrow = 2)
dev.off()
