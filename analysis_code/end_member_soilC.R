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
myCols <- wes_palette("GrandBudapest1", n = 3)


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
# add quartile-based split in starting soil C for plotting
HLvHH[which(HLvHH$soilChh <= 58.7), "soilCcat"] <- "low"
HLvHH[which(HLvHH$soilChh > 58.7 & HLvHH$soilChh <= 87.8), "soilCcat"] <- "medium"
HLvHH[which(HLvHH$soilChh > 87.8), "soilCcat"] <- "high"
HLvHH$soilCcat <- factor(HLvHH$soilCcat, levels = c("low", "medium", "high"))
HLvLL$soilCcat <- HLvHH$soilCcat

## Create data split by soil level ----
lowSoil <- filter(HLvHH, mst_cumdiff <= 5)
highSoil <- filter(HLvHH, mst_cumdiff >= 20)


#### MODEL SOIL C --------------------------------------------------------------

## From HH control ----
# create weights for modelling
varSite <- varIdent(form = ~ 1 | site)
# model fit
m1 <- gls(
  soilC_mgCg ~ mst_cumdiff * soilChh + mst_cumdiff * tap_mean, 
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

## Low warming categorical ----
# model fit
m3 <- gls(
  soilC_mgCg ~ soilCcat,
  weights = varSite,
  lowSoil,
  method = "ML"
)
# diagnose
r3 <- residuals(m3, type = "pearson")
par(mfrow = c(1, 3))
plot(r3 ~ fitted(m3))
boxplot(r3 ~ lowSoil$site)
hist(r3)
# effects
drop1(m3, test = "Chisq")
emmeans(m3, pairwise ~ soilCcat)

## High warming categorical ----
# model fit
m4 <- lm(
  soilC_mgCg ~ soilCcat,
  highSoil
)
# diagnose
r4 <- residuals(m4, type = "pearson")
par(mfrow = c(1, 3))
plot(r4 ~ fitted(m4))
boxplot(r4 ~ highSoil$site)
hist(r4)
# effects
drop1(m4, test = "Chisq")
emmeans(m4, pairwise ~ soilCcat)


#### PLOT ----------------------------------------------------------------------

## Coefficients for separate soil C amounts (quartile split, see above) ----
coef(lm(soilC_mgCg ~ mst_cumdiff, filter(HLvHH, soilCcat == "low")))
coef(lm(soilC_mgCg ~ mst_cumdiff, filter(HLvHH, soilCcat == "medium")))
coef(lm(soilC_mgCg ~ mst_cumdiff, filter(HLvHH, soilCcat == "high")))
coef(lm(soilC_mgCg ~ mst_cumdiff, HLvLL))

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
  geom_vline(xintercept = 16.97/1.2) +
  geom_vline(xintercept = 2.47/0.46) +
  geom_vline(xintercept = 7.41/0.88) +
  xlab("Cumulative warming (ºC)") +
  ylab("Diff. from origin (%)")

## To LL ----
toLL <- ggplot(HLvLL) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  guides(fill = "none", col = "none") +
  aes(x = mst_cumdiff, y = soilC_mgCg) +
  scale_x_continuous(limits = c(0, NA)) +
  geom_hline(yintercept = 0) +
  geom_point(col = "black", size = 1.25) +
  geom_smooth(method = "lm", se = F, col = "black") +
  geom_vline(xintercept = 70.8/2.35) +
  xlab("Cumulative warming (ºC)") +
  ylab("Diff. from destination (%)")

## Split by section to illustrate differences ----
# create data
lowSoilSum <- lowSoil %>%
  group_by(soilCcat) %>%
  summarise(
    mean = mean(soilC_mgCg, na.rm = T),
    se = sd(soilC_mgCg, na.rm = T)/sqrt(n())
  ) %>%
  ungroup
highSoilSum <- highSoil %>%
  group_by(soilCcat) %>%
  summarise(
    mean = mean(soilC_mgCg, na.rm = T),
    se = sd(soilC_mgCg, na.rm = T)/sqrt(n())
  ) %>%
  ungroup
# plot
lowSoilPlot <- ggplot(lowSoilSum) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  guides(fill = "none", col = "none") +
  scale_fill_manual(values = myCols)  +
  scale_y_continuous(limits = c(-22, 50)) +
  aes(x = soilCcat, fill = soilCcat, y = mean) +
  geom_hline(yintercept = 0) +
  geom_errorbar(aes(ymax = mean + se, ymin = mean - se), width = 0.4) +
  geom_bar(col = "black", stat = "identity") +
  xlab("Soil C content") +
  ylab("Diff. from origin (%)")
highSoilPlot <- ggplot(highSoilSum) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  guides(fill = "none", col = "none") +
  scale_y_continuous(limits = c(-22, 50)) +
  scale_fill_manual(values = myCols)  +
  aes(x = soilCcat, fill = soilCcat, y = mean) +
  geom_hline(yintercept = 0) +
  geom_errorbar(aes(ymax = mean + se, ymin = mean - se), width = 0.4) +
  geom_bar(col = "black", stat = "identity") +
  xlab("Soil C content") +
  ylab("Diff. from origin (%)")




## Write ----
postscript(
  file = "./plots/soilC_endmembers.eps", 
  width = mm2in(180),
  height = mm2in(60)
)
plot_grid(
  fromHH, lowSoilPlot, highSoilPlot, toLL, 
  rel_widths = c(2, 1, 1, 2), nrow = 1,
  align = "hv", axis = "tlbr"
  )
dev.off()





