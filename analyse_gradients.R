################################################################################
#### Project: TP Network soil
#### Title:   What are gradients doing
#### Author:  Tom Walker (thomas.walker@usys.ethz.ch)
#### Date:    11 June 2021
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
library(cowplot)
# my scripts
source("./code/make_biplot_data.R")
source("./code/mm2in.R")

## Data ----
siteData <- drake::readd(siteData) %>%
  mutate(continent = toupper(substr(continent, 1, 2))) %>%
  mutate(continent = factor(continent, levels = c("US", "SC", "AL", "AS")))


#### DATA ----------------------------------------------------------------------

## PCA ----
# construct
sitePCA <- prcomp(
  select(siteData, mat_mean:tap_mean, elev_range:tap_diff, mat_cumdiff:tap_cumdiff),
  center = T, scale = T
)
# inspect
summary(sitePCA)
sitePCA$rotation[, 1:2]
# plot
par(mfrow = c(1, 1))
biplot(sitePCA)
# extract scores
siteData$PC1 <- sitePCA$x[, 1]
siteData$PC2 <- sitePCA$x[, 2]
# extract loadings
siteLoadings <- sitePCA$rotation[, 1:2] %>%
  as.data.frame %>%
  rownames_to_column("variable") %>%
  mutate(variable = fct_reorder(variable, PC1, min))


#### CONTEXT -------------------------------------------------------------------

## PCA plots ----
# make biplot data
biPlotData <- make_biplot_data(
  siteData$PC1, 
  siteData$PC2, 
  siteData$continent
)
# biplot
pcaBiPlot <- ggplot(biPlotData) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  guides(col = "none") +
  aes(x = PC1, 
      y = PC2, 
      col = groupVar, 
      xend = meanPC1, 
      yend = meanPC2) +
  geom_segment() +
  geom_point() +
  geom_point(aes(x = meanPC1, y = meanPC2), 
             shape = 21, 
             fill = "white", 
             size = 3)
# loadings plots
sitePC1loads <- ggplot(siteLoadings) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  coord_flip() +
  scale_y_continuous(limits = c(-0.6, 0.6), breaks = c(-0.5, 0, 0.5)) +
  aes(x = variable, y = PC1, yend = 0, xend = variable) +
  geom_hline(yintercept = 0) +
  geom_segment() +
  geom_point(shape = 21, fill = "white", size = 2) +
  labs(x = "", y = "PC1 loading")
sitePC2loads <- ggplot(siteLoadings) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  coord_flip() +
  scale_y_continuous(limits = c(-0.6, 0.6), breaks = c(-0.5, 0, 0.5)) +
  aes(x = variable, y = PC2, yend = 0, xend = variable) +
  geom_hline(yintercept = 0) +
  geom_segment() +
  geom_point(shape = 21, fill = "white", size = 2) +
  labs(x = "", y = "PC2 loading")

## Continent plots ----
mstMean <- ggplot(siteData) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  aes(x = continent, fill = continent, y = mst_mean) +
  guides(fill = "none") +
  geom_boxplot() +
  labs(x = "", y = "MST (ºC)")
matMean <- ggplot(siteData) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  aes(x = continent, fill = continent, y = mat_mean) +
  guides(fill = "none") +
  geom_boxplot() +
  labs(x = "", y = "MAT (ºC)")
tapMean <- ggplot(siteData) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  aes(x = continent, fill = continent, y = tap_mean) +
  guides(fill = "none") +
  geom_boxplot() +
  labs(x = "", y = "TAP (mm)")
yearRange <- ggplot(siteData) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  aes(x = continent, fill = continent, y = year_range) +
  guides(fill = "none") +
  geom_boxplot() +
  labs(x = "", y = "Year range")
mstWarm <- ggplot(siteData) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  aes(x = continent, fill = continent, y = mst_cumdiff) +
  guides(fill = "none") +
  geom_boxplot() +
  labs(x = "", y = "Cumulative warming (ºC)")
pc1 <- ggplot(siteData) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  aes(x = continent, fill = continent, y = PC1) +
  guides(fill = "none") +
  geom_boxplot() +
  labs(x = "", y = "PC1 (52%)")
pc2 <- ggplot(siteData) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  aes(x = continent, fill = continent, y = PC2) +
  guides(fill = "none") +
  geom_boxplot() +
  labs(x = "", y = "PC2 (19%)")

## Combine ----
contextPlot1 <- plot_grid(pcaBiPlot, pc1, pc2, nrow = 1, rel_widths = c(2, 1, 1), axis = "tblr", align = "hv")
contextPlot2 <- plot_grid(sitePC1loads, sitePC2loads, nrow = 1, axis = "tblr")
contextPlot3 <- plot_grid(mstMean, matMean, tapMean, yearRange, mstWarm, nrow = 1, axis = "tblr", align = "hv")
postscript(file = "./plots/gradients_context.eps", width = mm2in(120), height = mm2in(180))
plot_grid(contextPlot1, contextPlot2, contextPlot3, nrow = 3, align = "hv", axis = "tblr")
dev.off()

 #### YEAR RANGE ---------------------------------------------------------------

## Is year range confounding ----
# model
anova(lm(mat_mean ~ year_range, siteData))
anova(lm(mst_mean ~ year_range, siteData))
anova(lm(tap_mean ~ year_range, siteData))
anova(lm(elev_range ~ year_range, siteData))
anova(lm(mat_diff ~ year_range, siteData))
anova(lm(mst_diff ~ year_range, siteData))
anova(lm(tap_diff ~ year_range, siteData))
anova(lm(mat_cumdiff ~ year_range, siteData))
anova(lm(mst_cumdiff ~ year_range, siteData)) # only significant relationship!
anova(lm(tap_cumdiff ~ year_range, siteData))
anova(lm(PC1 ~ year_range, siteData))
anova(lm(PC2 ~ year_range, siteData))
# plot
p0 <- ggplot(siteData) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  geom_point(shape = 21, fill = "grey")
p1 <- p0 + aes(x = year_range, y = elev_range)
p2 <- p0 + aes(x = year_range, y = mat_mean)
p3 <- p0 + aes(x = year_range, y = mst_mean)
p4 <- p0 + aes(x = year_range, y = tap_mean)
p5 <- p0 + aes(x = year_range, y = mst_diff)
p6 <- p0 + aes(x = year_range, y = tap_diff)
p7 <- p0 + aes(x = year_range, y = mst_cumdiff) + geom_smooth(method = "lm")
p8 <- p0 + aes(x = year_range, y = tap_cumdiff)
postscript(file = "./plots/gradients_yearrange.eps", width = mm2in(240), height = mm2in(110))
plot_grid(p1, p2, p3, p4, p5, p6, p7, p8, nrow = 2, align = "hv", axis = "tblr")
dev.off()
