################################################################################
#### Project: TP Network Soil
#### Title:   Main effect - predict soil carbon loss
#### Author:  Tom Walker (thomas.walker@usys.ethz.ch)
#### Date:    10 February 2022
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
library(randomForest)
library(caret)
library(nlme)
library(emmeans)
library(cowplot)
library(tidyverse)
source("./processing_code/mm2in.R")

## Data ----
# load
allData <- drake::readd(finalDF)
siteData <- drake::readd(siteData)


#### FORMAT --------------------------------------------------------------------

## Get data about end member controls ----
# high site
hhVars <- allData$metaMetrics %>%
  # gather from list object
  lapply(function(x){x$HH}) %>%
  # bind into data frame
  do.call(cbind, .) %>%
  as.data.frame %>%
  # select and rename columns with qualifier
  select(
    pH, soilC_mgCg:CmicNmic, R_ugCgh:TRate_d, 
    funSub1, funSub2, bacSub1, bacSub2
  ) %>%
  rename_with(~ paste0("HH_", .x))
# low site
llVars <- allData$metaMetrics %>%
  # gather from list object
  lapply(function(x){x$LL}) %>%
  # bind into data frame
  do.call(cbind, .) %>%
  as.data.frame %>%
  # select and rename columns with qualifier
  select(
    pH, soilC_mgCg:CmicNmic, R_ugCgh:TRate_d, 
    funSub1, funSub2, bacSub1, bacSub2
  ) %>%
  rename_with(~ paste0("LL_", .x))

## Assemble predictive dataset ----
allVars <- allData$metaMetrics$soilC_mgCg %>%
  # select abiotic context and HL soil C as response
  select(site, year_range, mat_mean:tap_cumdiff, HL) %>%
  # bind HH and LL context to this dataset, drop NAs
  bind_cols(., hhVars, llVars) %>%
  drop_na

## Assemble grouped dataset ----
meanVars <- allVars %>%
  select(site:mat_mean, mst_diff, HL, HH_soilC_mgCg, LL_soilC_mgCg) %>%
  group_by(site, year_range, mst_diff, mat_mean) %>%
  summarise(across(everything(), mean, na.rm = T))
seVars <- allVars %>%
  select(site:mat_mean, mst_diff, HL, HH_soilC_mgCg, LL_soilC_mgCg) %>%
  group_by(site, year_range, mst_diff, mat_mean) %>%
  summarise(across(everything(), function(x){sd(x, na.rm = T)}))
sumVars <- left_join(
  meanVars, seVars,
  by = c("site", "year_range", "mst_diff", "mat_mean"),
  suffix = c("_mean", "_se")
) %>%
  ungroup %>%
  mutate(yearFac = as.factor(as.numeric(cut(year_range, 2)))) %>%
  mutate(warmFac = as.factor(as.numeric(cut(mst_diff, 2))))

#### MODEL ---------------------------------------------------------------------

## Check covariance between experimental design predictors ----
psych::pairs.panels(
  select(allVars, year_range:tap_cumdiff), 
  method = "pearson",
  hist.col = "grey",
  density = T,
  ellipses = T
)

## Tiers 1-3 - starting soil C, treatments, background climate ----
# tier 1: we only know soil C content
m1 <- lm(HL ~ year_range * mst_diff, allVars)
summary(m1)
# tier 2: we also have a duration and warming magnitude in mind
m2a <- update(m1, ~.+ HH_soilC_mgCg)
m2b <- update(m1, ~.+ HH_soilC_mgCg + LL_soilC_mgCg)
summary(m2a)
summary(m2b)
anova(m1, m2a)
anova(m1, m2b)
m2 <- m2b
# tier 3: we also know something about the background climate
m3 <- update(m2, ~.+ mat_mean)
m3a <- update(m2, ~.+ mat_mean + tap_mean)
summary(m3)
summary(m3a)
anova(m2, m3, m3a)

## Tier 4 - soil C/N, pH ----
# covariances between easy remaining variables
psych::pairs.panels(
  select(
    allVars, HH_soilC_mgCg, LL_soilC_mgCg, 
    HH_pH, HH_soilN_mgCg, HH_soil_CN, 
    LL_pH, LL_soilN_mgCg, LL_soil_CN
  ),
  method = "pearson",
  hist.col = "grey",
  density = T,
  ellipses = T
)
# model
m4 <- update(m3, ~.+ HH_pH + HH_soil_CN + LL_pH + LL_soil_CN)
summary(m4)
anova(m3, m4)
m4 <- m3 # doesn't help, so retain previous best model

## Tier 5 - soil dissolved C/N ----
# covariances between next remaining variables
psych::pairs.panels(
  select(
    allVars, HH_soilC_mgCg, LL_soilC_mgCg, 
    HH_DOC_ugCg, HH_DN_ugNg, HH_NH4_ugNg, HH_NO3_ugNg,
    LL_DOC_ugCg, LL_DN_ugNg, LL_NH4_ugNg, LL_NO3_ugNg
  ),
  method = "pearson",
  hist.col = "grey",
  density = T,
  ellipses = T
)
# model
m5 <- update(m4, ~.+ HH_DOC_ugCg + LL_DOC_ugCg + HH_DN_ugNg + LL_DN_ugNg)
summary(m5)
anova(m4, m5)
m5 <- m4 # doesn't help, so retain previous
m6 <- update(m5, ~.+ HH_NH4_ugNg + HH_NO3_ugNg + LL_NH4_ugNg + LL_NO3_ugNg)
summary(m6)
anova(m5, m6)
m6 <- m5 # doesn't help, so retain previous

## Tier 6 - microbial C/N ----
# covariances between next remaining variables
psych::pairs.panels(
  select(
    allVars, HH_soilC_mgCg, LL_soilC_mgCg, 
    HH_Cmic_ugCg, HH_Nmic_ugNg, HH_CmicNmic,
    LL_Cmic_ugCg, LL_Nmic_ugNg, LL_CmicNmic
  ),
  method = "pearson",
  hist.col = "grey",
  density = T,
  ellipses = T
)
# model
m7 <- update(m6, ~.+ HH_CmicNmic + LL_CmicNmic)
summary(m7)
anova(m6, m7)
m7 <- m6 # doesn't help

## Tier 7 - microbial respiration ----
# covariances between next remaining variables
psych::pairs.panels(
  select(
    allVars, 
    HH_soilC_mgCg, LL_soilC_mgCg, 
    HH_R_ugCgh, HH_Rm_ugCmic, HH_Rc_ugSOC,
    LL_R_ugCgh, LL_Rm_ugCmic, LL_Rc_ugSOC
  ),
  method = "pearson",
  hist.col = "grey",
  density = T,
  ellipses = T
)
# model
m8a <- update(m7, ~.+ HH_Rm_ugCmic + LL_Rm_ugCmic)
m8b <- update(m7, ~.+ HH_Rc_ugSOC + LL_Rc_ugSOC)
summary(m8a)
summary(m8b)
anova(m7, m8a)
anova(m7, m8b)
m8 <- m7 # doesn't help

## Tier 8 - microbial physiology ----
# covariances between next remaining variables
psych::pairs.panels(
  select(
    allVars, 
    HH_soilC_mgCg, LL_soilC_mgCg, 
    HH_CUE, HH_TRate_d,
    LL_CUE, LL_TRate_d
  ),
  method = "pearson",
  hist.col = "grey",
  density = T,
  ellipses = T
)
# model
m9a <- update(m8, ~.+ HH_CUE + LL_CUE + HH_TRate_d + LL_TRate_d)
m9b <- update(m8, ~.+ HH_CUE + LL_CUE)
m9c <- update(m8, ~.+ HH_CUE)
summary(m9a)
summary(m9b)
summary(m9c)
anova(m8, m9a)
anova(m8, m9b)
anova(m8, m9c)
m9 <- m8

## Tier 9 - microbial community ----
psych::pairs.panels(
  select(
    allVars, 
    HH_soilC_mgCg, LL_soilC_mgCg, 
    HH_funSub1, HH_funSub2, HH_bacSub1, HH_bacSub2,
    LL_funSub1, LL_funSub2, LL_bacSub1, LL_bacSub2
  ),
  method = "pearson",
  hist.col = "grey",
  density = T,
  ellipses = T
)
# model
m10a <- update(m9, ~.+ HH_funSub1 + HH_funSub2 + LL_funSub1 + LL_funSub2)
m10b <- update(m9, ~.+ HH_bacSub1 + HH_bacSub2 + LL_bacSub1 + LL_bacSub2)
summary(m10a)
summary(m10b)
anova(m9, m10a)
anova(m9, m10b)
m10 <- m9


## Test loss of explanatory power without end member ----
m11 <- update(m10, ~.- LL_soilC_mgCg)
summary(m11)


#### RANDOM FOREST -------------------------------------------------------------

## Create validation dataset ----
validation <- data.frame(
  rf_train = NA, rf_test = NA, 
  lm_train = NA, lm_test = NA
)

## Set up dataset ----
rfHL <- allVars$HL
rfVars <- allVars %>% select(-site, -HL)

# ## Create tuning parameters ----
# # basic settings
# seed <- 7
# metric <- "RMSE"
# control <- trainControl(
#   method = "repeatedcv", 
#   number = 10, 
#   repeats = 3, 
#   search = 'grid'
# )
# tunegrid <- expand.grid(
#   .maxnodes=c(20, 30,40,50), 
#   .ntree=c(900, 1000, 1100)
# )
# # create RF list
# customRF <- list(
#   type = "Regression", 
#   library = "randomForest", 
#   loop = NULL,
#   parameters = data.frame(
#     parameter = c("maxnodes", "ntree"), 
#     class = rep("numeric", 2), 
#     label = c("maxnodes", "ntree")
#   ),
#   grid = function(x, y, len = NULL, search = "grid"){},
#   fit = function(x, y, wts, param, lev, last, weights, classProbs, ...){
#     randomForest(
#       x, y, 
#       maxnodes = param$maxnodes, 
#       ntree = param$ntree, 
#       ...
#     )
#   },
#   predict = function(modelFit, newdata, preProc = NULL, submodels = NULL){
#     predict(modelFit, newdata)
#   },
#   prob = function(modelFit, newdata, preProc = NULL, submodels = NULL){
#     predict(modelFit, newdata, type = "prob")
#   },
#   sort = function(x){x[order(x[,1]),]},
#   levels = function(x){x$classes}
# )
# 
# ## Run Grid Search Model ----
# set.seed(seed)
# rf_gridsearch <- train(
#   x = rfVars, 
#   y = rfHL, 
#   method = customRF, 
#   metric = metric, 
#   tuneGrid = tunegrid, 
#   trControl = control
# )

## Build final RF model ----
rfFinal <- randomForest(
  x = rfVars,
  y = rfHL,
  maxnodes = 40,
  ntree = 1100
)


#### VALIDATE ------------------------------------------------------------------

## Build validation dataset ----
# create base dataset
validation <- allVars %>%
  select(year_range, mst_diff, mat_mean, HH_soilC_mgCg, LL_soilC_mgCg, HL) %>%
  rename(HH = HH_soilC_mgCg, LL = LL_soilC_mgCg) %>%
  mutate(rf_test = NA, lm_test = NA)
# row by row predictions
for(i in 1:nrow(allVars)){
  # setup training and test data
  train <- allVars[-i, ]
  test <- allVars[i, ]
  # do models
  modLM <- update(m10, data = train)
  modRF <- update(rfFinal, x = select(train, -HL), y = train$HL)
  # predict
  validation[i, "lm_test"] <- predict(modLM, test)
  validation[i, "rf_test"] <- predict(modRF, test)
}
# model for r-squared
lmLM <- lm(lm_test ~ HL, validation)
lmRF <- lm(rf_test ~ HL, validation)
# add residuals and other error estimate to model
validation <- validation %>%
  mutate(lm_residuals = residuals(lmLM)) %>%
  mutate(lm_error = lm_test - HL)

## Explore prediction accuracy ----
# solve points where accuracy models intercept
cmLM <- rbind(c(0, 1),coef(lmLM))
solveLM <- unique(c(-solve(cbind(cmLM[,2],-1)) %*% cmLM[,1]))
cmRF <- rbind(c(0, 1),coef(lmRF))
solveRF <- unique(c(-solve(cbind(cmRF[,2],-1)) %*% cmRF[,1]))
# model what factors control magnitude of error
mX1 <- lm(lm_error ~ year_range + mst_diff + HH + LL + HL, validation)
anova(mX1)
mX2 <- update(mX1, ~.- year_range - mst_diff - HH - LL)
coef(mX2)
# plot accuracy against 1:1 line
lmTestPlot <- ggplot(validation) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_y_continuous(limits = c(24, 165)) +
  scale_x_continuous(limits = c(24, 165)) +
  aes(x = HL, y = lm_test) +
  geom_point() +
  geom_segment(
    y = 0, 
    yend = solveLM, 
    x = solveLM, 
    xend = solveLM, 
    size = 0.5, 
    col = "red"
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_abline(slope = coef(lmLM)[2], intercept = coef(lmLM)[1]) +
  xlab(expression(paste("Measured soil C (mg C ", g^{-1}, ")"))) +
  ylab(expression(paste("Predicted soil C (mg C ", g^{-1}, ")")))
# plot prediction error against HL concentration
errorHL <- ggplot(validation) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  aes(x = HL, y = lm_error, col = outlier) +
  scale_color_manual(values = c("black", "grey")) +
  guides(color = "none") +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_abline(slope = -0.20, intercept = 19.3) +
  xlab(expression(paste("Measured soil C (mg C ", g^{-1}, ")"))) +
  ylab(expression(paste("Prediction error (mg C ", g^{-1}, ")")))
# combine
cowplot::plot_grid(lmTestPlot, errorHL)