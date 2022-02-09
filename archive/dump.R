################################################################################
#### Project: TP Network
#### Title:   Main effects
#### Author:  Tom Walker (thomas.walker@usys.ethz.ch)
#### Date:    09 June 2021
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

## Custom functions ----
source("./code/rrs.R")

## Data ----
# load
allData <- drake::readd(finalDF) %>%
  filter(site != "YAN")
siteData <- drake::readd(siteData) %>%
  filter(site != "YAN")


#### FORMAT --------------------------------------------------------------------

## Categorical TAP data ----
siteData[siteData$tap_mean < 337, "tap_cat"] <- "low"
siteData[siteData$tap_mean >= 337 & siteData$tap_mean < 666, "tap_cat"] <- "mid"
siteData[siteData$tap_mean >= 666 & siteData$tap_mean < 1068, "tap_cat"] <- "high"
siteData[siteData$tap_mean >= 1068, "tap_cat"] <- "very high"
  

## Create factor from site information ----
allData$site <- factor(
  allData$site, 
  levels = levels(siteData$site)
)

## Site-wise treatment means ± SEs ----
# separate data frames
means <- allData %>%
  group_by(site, treat_a) %>%
  summarise(across(funFull1:TRate_d, mean, na.rm = T)) %>%
  ungroup %>%
  rename_with(~paste0("mean_", .x), funFull1:TRate_d)
SEs <- allData %>%
  group_by(site, treat_a) %>%
  summarise(across(funFull1:TRate_d, ~sd(.x, na.rm = T)/sqrt(n()))) %>%
  ungroup %>%
  rename_with(~paste0("se_", .x), funFull1:TRate_d)
# bind
summaries <- left_join(means, SEs)

## Nested data for meta-analysis metrics ----
# build pipeline
nestData <- allData %>%
  # group and nest by site
  group_by(site) %>%
  nest() %>%
  # separate drivers and responses and remove full dataset
  mutate(drivers = map(data, ~select(.x, gradient:P_ann))) %>%
  mutate(responses = map(data, ~select(.x, funFull1:TRate_d))) %>%
  select(-data) %>%
  # filter for different treatments
  mutate(HH = map2(drivers, responses, ~filter(.y, .x$treat_a == "HH"))) %>%
  mutate(HL = map2(drivers, responses, ~filter(.y, .x$treat_a == "HL"))) %>%
  mutate(LL = map2(drivers, responses, ~filter(.y, .x$treat_a == "LL"))) %>%
  # calculate site-wise treatment means
  mutate(HHmean = map(HH, ~summarise(.x, across(everything(), mean, na.rm = T)))) %>%
  mutate(LLmean = map(LL, ~summarise(.x, across(everything(), mean, na.rm = T)))) %>%
  # calculate plot-wise response of HL from both end members
  mutate(HLvHH = map2(HHmean, HL, ~do.call(rbind, apply(.y, 1, function(z) (z - .x)/.x * 100)))) %>%
  mutate(HLvLL = map2(LLmean, HL, ~do.call(rbind, apply(.y, 1, function(z) (z - .x)/.x * 100)))) %>%
  mutate(HHvLL = map2(LLmean, HH, ~do.call(rbind, apply(.y, 1, function(z) (z - .x)/.x * 100))))
# separate
HLHH <- nestData %>%
  select(HLvHH) %>%
  unnest(HLvHH) %>%
  as.data.frame %>%
  left_join(., siteData)
HLLL <- nestData %>%
  select(HLvLL) %>%
  unnest(HLvLL) %>%
  as.data.frame %>%
  left_join(., siteData)
HHLL <- nestData %>%
  select(HHvLL) %>%
  unnest(HHvLL) %>%
  as.data.frame %>%
  left_join(., siteData)


hlTurfs <- allData %>%
  filter(treat_a == "HL") %>%
  filter(Rm_ugCmic < 1500) %>%
  left_join(., siteData)




#### PLOT CHANGE ---------------------------------------------------------------

mA <- lm(soilC_mgCg ~ 1, HLLL)
mB <- update(mA, ~.+ mst_mean * tap_mean)

anova(mA, mB)

m1 <- lm(Cmic_ugCg ~ mst_mean * tap_mean + mst_cumdiff * tap_mean, HHLL)
m2 <- lm(Cmic_ugCg ~ mst_mean * tap_mean + mst_cumdiff * tap_mean, HLHH)
m3 <- lm(Cmic_ugCg ~ mst_mean * tap_mean + mst_cumdiff * tap_mean, HLLL)

anova(m1)
anova(m2)
anova(m3)


modData <- data.frame(
  tap_mean = rep(c(100, 500, 1000, 5000), each = 21),
  mst_cumdiff = rep(0:20, 4)
)
modData$HHLL <- predict(m1, newdata = modData)
modData$HLHH <- predict(m2, newdata = modData)
modData$HLLL <- predict(m3, newdata = modData)

ggplot(modData) +
  aes(x = mst_cumdiff, y = HLLL) +
  geom_line() +
  facet_wrap(~tap_mean, nrow = 1)
  



anova(m1)
anova(m2)
anova(m3)

ggplot(HHLL) +
  aes(x = mst_cumdiff, col = tap_cat, y = soilC_mgCg) +
  geom_smooth() +
  geom_point()

anova(lm(soilC_mgCg ~ mst_cumdiff * tap_mean + mst_cumdiff * tap_cumdiff, HHLL))
anova(lm(soilC_mgCg ~ mst_cumdiff * tap_mean, HLHH))
summary(lm(soilC_mgCg ~ mst_cumdiff * tap_mean, HLLL))


soilCtest <- ggplot(HHLL) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  aes(x = tap_mean, y = soilC_mgCg) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.5) +
  geom_point(shape = 21, fill = "#D9A464") +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2)) +
  labs(y = "HH from LL (%)")
soilCdiv <- ggplot(HLHH) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  aes(x = tap_mean, y = soilC_mgCg) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.5) +
  geom_point(shape = 21, fill = "#D9A464") +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2)) +
  labs(y = "HL from HH (%)")
soilCcon <- ggplot(HLLL) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  aes(x = tap_mean, y = soilC_mgCg) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.5) +
  geom_point(shape = 21, fill = "#D9A464") +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2)) +
  labs(y = "HL from LL (%)")

cowplot::plot_grid(soilCtest, soilCdiv, soilCcon, nrow = 1)


#### SOIL C --------------------------------------------------------------------

## Model fit ----
# fit
m1 <- gls(
  soilC_mgCg ~ treat_a * site,
  method = "ML",
  na.action = "na.exclude",
  data = allData
)
# diagnose
r1 <- residuals(m1, type = "pearson")
par(mfrow = c(1, 3))
plot(r1 ~ fitted(m1))
boxplot(r1 ~ allData$treat_a)
hist(r1)
# main effects
anova(m1, update(m1, ~.- treat_a:site))
# post-hoc
emmeans(m1, pairwise ~ treat_a | site)

## Divergence from/to end members ----
rrs(allData)

## Plot ----
# plot
ggplot(summaries) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  guides(fill = "none") +
  scale_fill_manual(values = c("#666699", "#D9A464", "#996666")) +
  aes(
    x = treat_a,
    fill = treat_a,
    y = mean_soilC_mgCg,
    ymin = mean_soilC_mgCg - se_soilC_mgCg,
    ymax = mean_soilC_mgCg + se_soilC_mgCg,
  ) +
  geom_errorbar(position = position_dodge(0.9), width = 0.2) +
  geom_bar(position = position_dodge(0.9), col = "black", stat = "identity") +
  facet_wrap(~site, nrow = 1, scales = "free")









