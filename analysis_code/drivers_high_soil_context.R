################################################################################
#### Project: TP Network Soil
#### Title:   Why does high site soil C vary?
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


#### FORMAT --------------------------------------------------------------------

## Subset and filter ----
highSoil <- allData$fullDF %>%
  # filter for high elevation control
  filter(treat_a == "HH") %>%
  # generate categorical soil C 
  # add mineral/organic soil split
  mutate(soilCcat = ifelse(soilC_mgCg <= 100, "mineral", "organic")) %>%
  # select variables of interest
  select(
    # geographical/experimental
    site, continent, elev_m, plot_m2,
    # climate
    T_ann_cor:P_ann,
    # soil system
    pH, soilCcat, soilC_mgCg:rpa_mom
  ) %>%
  drop_na

## Check which  DRIFT variables antagonistic (relative) ----
# plot
par(mfrow = c(1, 3))
plot(rpa_spom ~ rpa_cpom, highSoil)
abline(lm(rpa_spom ~ rpa_cpom, highSoil), col = "red")
plot(rpa_spom ~ rpa_mom, highSoil)
plot(rpa_cpom ~ rpa_mom, highSoil)


#### ANALYSE -------------------------------------------------------------------

## Continuous model build ----
# build models
m1 <- lme(
  fixed = soilC_mgCg ~ 
    continent + elev_m + plot_m2 + T_sum_cor + P_ann + pH + rpa_spom + rpa_mom,
  random = ~ 1 | site,
  data = highSoil,
  method = "ML",
  na.action = "na.exclude"
)
# diagnose
r1 <- residuals(m1, type = "pearson")
par(mfrow = c(1, 3))
plot(r1 ~ fitted(m1))
hist(r1)
qqnorm(r1)
# significance
drop1(m1, test = "Chisq")


## Categorical soil C on SPOM ----
m3 <- lme(
  fixed = rpa_spom ~ soilCcat,
  random = ~ 1 | site,
  data = highSoil,
  method = "ML",
  na.action = "na.exclude"
)
# diagnose
r3 <- residuals(m3, type = "pearson")
par(mfrow = c(1, 3))
plot(r3 ~ fitted(m3))
hist(r3)
qqnorm(r3)
# significance
drop1(m3, test = "Chisq")

## Categorical soil C on MOM ----
m4 <- lme(
  fixed = rpa_mom ~ soilCcat,
  random = ~ 1 | site,
  data = highSoil,
  method = "ML",
  na.action = "na.exclude"
)
# diagnose
r4 <- residuals(m4, type = "pearson")
par(mfrow = c(1, 3))
plot(r4 ~ fitted(m4))
hist(r4)
qqnorm(r4)
# significance
drop1(m4, test = "Chisq")


#### PLOT ----------------------------------------------------------------------

## Build plots ----
spom_soc <- ggplot(highSoil) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  guides(fill = "none", col = "none") +
  aes(x = rpa_spom, y = soilC_mgCg) +
  geom_point(col = "black", size = 1.25) +
  geom_smooth(method = "lm", se = F, col = "black") +
  ylab(expression(paste("HH soil C (mg C ", g^{-1}, " dm)"))) +
  xlab("sPOM (ratio peak area)")
cpom_soc <- ggplot(highSoil) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  guides(fill = "none", col = "none") +
  aes(x = rpa_cpom, y = soilC_mgCg) +
  geom_point(col = "black", size = 1.25) +
  geom_smooth(method = "lm", se = F, col = "black") +
  ylab(expression(paste("HH soil C (mg C ", g^{-1}, " dm)"))) +
  xlab("cPOM (ratio peak area)")
mom_soc <- ggplot(highSoil) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  guides(fill = "none", col = "none") +
  aes(x = rpa_mom, y = soilC_mgCg) +
  geom_point(col = "black", size = 1.25) +
  ylab(expression(paste("HH soil C (mg C ", g^{-1}, " dm)"))) +
  xlab("MOM (ratio peak area)")

## Combine ----
postscript(
  file = "./figure_builds/soilC_som_context.eps", 
  width = mm2in(180),
  height = mm2in(60)
)
plot_grid(
  spom_soc, cpom_soc, mom_soc,
  nrow = 1, 
  align = "hv", axis = "tlbr"
)
dev.off()





