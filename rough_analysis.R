################################################################################
#### Project: TP Network
#### Title:   Rough analysis
#### Author:  Tom Walker (thomas.walker@usys.ethz.ch)
#### Date:    24 November 2020
#### ---------------------------------------------------------------------------


#### PROLOGUE ------------------------------------------------------------------

## Options ----
# remove objects from global environment
rm(list = ls())

# configure default R session options (no factors, bias against scientific #s)
options(stringsAsFactors = F,
        scipen = 6)


## Libraries ----
# standard library set
source("packages.R")

## Code ----
scripts <- list.files(
  "./code",
  full.names = T
)

sapply(
  scripts,
  source
)

## Load data ----
drake::loadd()


#### NMDS ----------------------------------------------------------------------

## Build ----
funNMDS <- cmdscale(
  d = vegdist(cleanData$microbes$counts$fungi),
  k = 2
) %>%
  as.data.frame
bacNMDS <- cmdscale(
  d = vegdist(cleanData$microbes$counts$bacteria),
  k = 2
) %>% 
  as.data.frame

## Create biplot data ----
funContinent <- make_biplot_data(
  PC1 = funNMDS[, 1],
  PC2 = funNMDS[, 2],
  groupVar = cleanData$metadata$plots$continent
)
bacContinent <- make_biplot_data(
  PC1 = bacNMDS[, 1],
  PC2 = bacNMDS[, 2],
  groupVar = cleanData$metadata$plots$continent
)
funSite <- make_biplot_data(
  PC1 = funNMDS[, 1],
  PC2 = funNMDS[, 2],
  groupVar = cleanData$metadata$plots$site
)
bacSite <- make_biplot_data(
  PC1 = bacNMDS[, 1],
  PC2 = bacNMDS[, 2],
  groupVar = cleanData$metadata$plots$site
)

## Plot ----
baseBiPlot <- ggplot(bacContinent) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  aes(x = PC1, y = PC2, col = groupVar) +
  geom_point() +
  geom_segment(aes(x = meanPC1, y = meanPC2, xend = PC1, yend = PC2)) +
  geom_point(aes(x = meanPC1, y = meanPC2), size = 5, shape = 21, fill = "white") +
  xlab("NMDS 1") + ylab("NMDS 2")

cowplot::plot_grid(
  baseBiPlot,
  baseBiPlot %+% bacSite,
  baseBiPlot %+% funContinent,
  baseBiPlot %+% funSite
)


#### WARMING -------------------------------------------------------------------

## Bind ----
funAll <- bind_cols(cleanData$metadata$plots, funNMDS)
bacAll <- bind_cols(cleanData$metadata$plots, bacNMDS)

## Plot ----
bacResponse <- ggplot(bacAll) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  aes(x = treat_a, y = V2) +
  stat_boxplot(geom = "errorbar", width = 0.2) +
  geom_boxplot(outlier.shape = NA, fill = "grey") +
  facet_wrap(~ site, scales = "free", nrow = 1) +
  ylab("NMDS 2") + xlab("")

cowplot::plot_grid(
  bacResponse,
  bacResponse %+% funAll,
  nrow = 2,
  align = "hv"
)


#### PROCESSES -----------------------------------------------------------------

## Impute ----
nums <- cleanData$processes %>%
  select(h2o_mgg:TRate_d) %>%
  mutate(funPC1 = funNMDS$V1,
         funPC2 = funNMDS$V2,
         bacPC1 = bacNMDS$V1,
         bacPC2 = bacNMDS$V2) %>%
  filter(complete.cases(soilC_mgCg)) %>%
  filter(complete.cases(G_ugCgh)) %>%
  filter(complete.cases(CmicNmic)) %>%
  filter(R_ugCgh < 10)

numsRF <- nums %>%
  select(h2o_mgg:R_ugCgh, -DNA_ugg, funPC1:bacPC2)


set.seed(4)
rf <- randomForest(
  R_ugCgh ~ .,
  data = numsRF
)
rf
impRF <- importance(rf)[, 1]
plotRF <- data.frame(variable = names(impRF),
                     importance = impRF,
                     relImportance = impRF/sum(impRF) * 100)
plotRF$plotVariable <- as.character(plotRF$variable)
plotRF$plotVariable <- factor(plotRF$plotVariable,
                              levels = plotRF$plotVariable[order(plotRF$importance)])

p4 <- ggplot(plotRF) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  coord_flip() +
  aes(x = plotVariable, y = relImportance) +
  geom_bar(stat = "identity", fill = "goldenrod", col = "black") +
  labs(x = "", y = "Relative importance (%)")

## Random forest ----


## Analyse ----
numsAll <- bind_cols(cleanData$processes, cleanData$metadata$plots)

ggplot(numsAll) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  aes(x = treat_a, y = R_ugCgh) +
  stat_boxplot(geom = "errorbar", width = 0.2) +
  geom_boxplot(outlier.shape = NA, fill = "grey") +
  facet_wrap(~ site, scales = "free", nrow = 1) +
  xlab("")











