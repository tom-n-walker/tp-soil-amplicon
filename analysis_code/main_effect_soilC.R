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

## Split data ----
fullData <- allData$fullDF

## Assemble data on background context ----
# get meaningful context
soilContext <- data.frame(
  soilChh = allData$metaMetrics$soilC_mgCg$HH,
  soilPC1hh = allData$metaMetrics$soilPC1$HH,
  soilPC2hh = allData$metaMetrics$soilPC2$HH,
  soilPC1hl = allData$metaMetrics$soilPC1$HL,
  soilPC2hl = allData$metaMetrics$soilPC2$HL
)
# assemble drivers
drivers <- allData$metaMetrics$soilC_mgCg %>%
  select(site, year_range, elev_range,
         mat_mean, mst_mean, tap_mean, 
         mat_diff, mst_diff, tap_diff, 
         mat_cumdiff, mst_cumdiff, tap_cumdiff) %>%
  bind_cols(., soilContext)
# assemble responses
responses <- allData$metaMetrics$soilC_mgCg %>%
  select(HLvHHmean:HHvLLmean)
# bind all
causeEffect <- bind_cols(drivers, responses) %>%
  mutate(tap_mean = tap_mean / 1000)

## Build factors ----
fullData$site <- factor(fullData$site, levels = siteData$site)
causeEffect$site <- factor(causeEffect$site, levels = siteData$site)


#### MODEL ---------------------------------------------------------------------

## Model treatments ----
# build
m1 <- lm(
  log10(soilC_mgCg) ~ treat_a * site,
  data = fullData,
  na.action = na.exclude
)
# diagnose
r1 <- residuals(m1)
par(mfrow = c(1, 4))
plot(r1 ~ fitted(m1))
boxplot(r1 ~ fullData$treat_a)
boxplot(r1 ~ fullData$site)
hist(r1)
# main effects
anova(m1)
# post-hoc
m1ref <- regrid(ref_grid(m1))
emmeans(m1ref, pairwise ~ treat_a | site)$contrasts

## Model change from HH end member ----
# build model
m2 <- lm(HLvHHmean ~ mst_cumdiff * soilChh + tap_mean + soilChh, causeEffect)
m2 <- lm(HLvHHmean ~ mst_cumdiff * soilChh + tap_mean + soilChh, filter(causeEffect, tap_mean < 6))
anova(m2)
# build dataset for plotting
predictHHvHL <- data.frame(
  mst_cumdiff = rep(0:30, 2), # create sequence of cumulative warming
  tap_mean = mean(causeEffect$tap_mean, na.rm = T), # use mean TAP for this
  soilCcat = rep(c("L", "H"), each = 31),
  soilChh = rep(c(68, 122), each = 31) # use 1st and 3rd quantiles for low/high
)
predictHHvHL$predicted <- predict(m2, newdata = predictHHvHL)
# plot
HLvHHtap <- ggplot(causeEffect) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  aes(x = tap_mean, y = HLvHHmean) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(shape = 21, fill = "#666699") +
  geom_smooth(method = "lm", col = "#666699") +
  labs(x = "TAP (m)", y = "HL from HH (%)")
HLvHHwarm <- ggplot(causeEffect) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  aes(x = mst_cumdiff, y = HLvHHmean) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(shape = 21, fill = "#666699") +
  geom_smooth(method = "lm", col = "#666699") +
  labs(x = "Cumulative MST (ºC)", y = "HL from HH (%)")
HLvHHsoilC <- ggplot(causeEffect) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  aes(x = soilChh, y = HLvHHmean) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(shape = 21, fill = "#666699") +
  geom_smooth(method = "lm", col = "#666699") +
  labs(x = "Initial soil C (mg/g)", y = "HL from HH (%)")
HLvHHmodel <- ggplot(predictHHvHL) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  guides(col = "none") +
  aes(x = mst_cumdiff, y = predicted, col = soilCcat) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line() +
  labs(x = "Cumulative MST (ºC)", y = "HL from HH (%)")

## Model change from LL end member ----
# build model
m3 <- lm(HLvLLmean ~ mst_cumdiff + soilChh + tap_mean + soilChh, causeEffect)
m3 <- lm(HLvLLmean ~ mst_cumdiff * soilChh + tap_mean * soilChh, filter(causeEffect, tap_mean < 6))
anova(m3)
# plot
HLvLLtap <- ggplot(causeEffect) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  aes(x = tap_mean, y = HLvLLmean) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(shape = 21, fill = "#996666") +
  geom_smooth(data = filter(causeEffect, tap_mean < 6), method = "lm", col = "#996666") +
  labs(x = "TAP (m)", y = "HL from LL (%)")
HLvLLwarm <- ggplot(causeEffect) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  aes(x = mst_cumdiff, y = HLvLLmean) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(shape = 21, fill = "#996666") +
  geom_smooth(method = "lm", col = "#996666") +
  labs(x = "Cumulative MST (ºC)", y = "HL from LL (%)")

## Plot together ----
postscript(file = "./plots/soilC_HHvHLvLL.eps", width = mm2in(110), height = mm2in(240))
plot_grid(HLvHHwarm, HLvLLwarm, 
          HLvHHtap, HLvLLtap,
          HLvHHsoilC, NULL,
          HLvHHmodel, NULL, 
          nrow = 4)
dev.off()





#### PLOT MAIN -----------------------------------------------------------------

## Organise data ----
meansSEs <- allData$meansSEs %>%
  mutate(site = factor(site, levels = levels(siteData$site))) %>%
  select(site, treat_a, mean_soilC_mgCg, se_soilC_mgCg)

## Plot ----
treatPlot <- ggplot(meansSEs) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  aes(x = treat_a, 
      y = mean_soilC_mgCg, 
      ymax = mean_soilC_mgCg + se_soilC_mgCg,
      ymin = mean_soilC_mgCg - se_soilC_mgCg,
      fill = treat_a
      ) +
  scale_fill_manual(values = myCols) +
  guides(fill = "none") +
  labs(x = "", y = expression(paste("Soil C (mg ", g^{-1}, ")"))) +
  geom_errorbar(width = 0.2, position = position_dodge(0.5)) +
  geom_point(size = 3, shape = 21, position = position_dodge(0.5)) +
  facet_wrap(~ site, nrow = 1, scales = "free")
