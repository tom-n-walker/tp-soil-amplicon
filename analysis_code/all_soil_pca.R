################################################################################
#### Project: TP Network soil
#### Title:   Main analysis - PCA on all variables
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

## Scrips ----
source("./processing_code/apply_mice.R")
source("./processing_code/make_biplot_data.R")
source("./processing_code/mm2in.R")

## Libraries ----
# standard library set
library(tidyverse)
library(data.table)
library(nlme)
library(mice)
library(emmeans)

## Data ----
allData <- drake::readd(finalDF)
siteData <- drake::readd(siteData)

## Contributing variables ----
# plot colour scale
myCols <- c("#666699", "#D9A464", "#996666")

#### DATA ----------------------------------------------------------------------

## Format numeric data ----
responses <- allData$fullDF %>%
  # get data
  select(
    pH:CmicNmic, 
    R_ugCgh:rpa_mom, 
    funFull1:funFull2, 
    bacFull1:bacFull2
    ) %>%
  # impute missing values
  apply_mice(., 5)

## Format drivers ----
drivers <- allData$fullDF %>%
  select(site:P_ann)

## Variable names ----
varNames <- c(
  # bulk soil pH and C/N pools
  "pH", "Soil C", "Soil N", "Soil C:N", 
  # dissolved C/N pools
  "DOC", "TDN", "NH4", "NO3", "DON",
  # microbial C/N pools
  "Microbial C", "Microbial N", "Microbial C:N", 
  # microbial processes
  "R", "G", "Mass-specific R", "Mass-specific G", "CUE", "Turnover rate",
  # spectroscopy
  "sPOM", "cPOM", "MOM",
  # mcc
  "Fungi 1", "Fungi 2", "Bacteria 1", "Bacteria 2"
)

#### DO PCA --------------------------------------------------------------------

## Do PCA ----
# full data
fullPCA <- prcomp(responses, scale = T, center = T)
fullScores <- bind_cols(drivers, as.data.frame(fullPCA$x[, 1:2]))
fullSummary <- fullScores %>% 
  group_by(treat_a) %>%
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
  fullScores$treat_a
)


#### PLOT PCA ------------------------------------------------------------------

## Plot full PCA scores ----
# biplot
fullPlot <- ggplot(fullPlotData) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_colour_manual(values = myCols) +
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
# separate axes
fullPC1 <- ggplot(fullSummary) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_fill_manual(values = myCols) +
  guides(fill = "none") +
  aes(x = treat_a, y = PC1mean, fill = treat_a, 
      ymax = PC1mean + PC1se, ymin = PC1mean - PC1se) +
  geom_hline(yintercept = 0) +
  geom_errorbar(width = 0.2) +
  geom_bar(stat = "identity", col = "black") +
  labs(x = "", y = "PC1 score")
fullPC2 <- ggplot(fullSummary) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_fill_manual(values = myCols) +
  guides(fill = "none") +
  aes(x = treat_a, y = PC2mean, fill = treat_a, 
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
  file = "./figure_builds/full_pca_scores.eps", 
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
  file = "./figure_builds/full_pca_loadings.eps", 
  width = mm2in(180), 
  height = mm2in(75)
)
cowplot::plot_grid(
  fullPC1loads, fullPC2loads, 
  nrow = 1, 
  align = "hv"
)
dev.off()


#### ANALYSE PCA ---------------------------------------------------------------

## Permanova ----
fullPCAmodel <- vegan::adonis(responses ~ treat_a * site, drivers)
fullPCAmodel

## PC scores ----
# build models
fullPC1model <- lme(
  PC1 ~ treat_a, 
  random = ~1 | site,
  fullScores, 
  method = "ML"
)
fullPC2model <- lme(
  PC2 ~ treat_a, 
  random = ~1 | site,
  fullScores, 
  method = "ML"
)
# main effects
anova(fullPC1model, update(fullPC1model, ~.- treat_a))
anova(fullPC2model, update(fullPC2model, ~.- treat_a))
# post-hoc
emmeans(fullPC1model, pairwise ~ treat_a)
emmeans(fullPC2model, pairwise ~ treat_a)
