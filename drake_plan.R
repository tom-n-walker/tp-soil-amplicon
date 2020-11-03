################################################################################
#### Project: TP Network
#### Title:   Drake plan
#### Author:  Tom Walker (thomas.walker@usys.ethz.ch)
#### Date:    03 November 2020
#### ---------------------------------------------------------------------------


#### PROLOGUE ------------------------------------------------------------------

## Options ----
# remove objects from global environment
rm(list = ls())

# configure default R session options (no factors, bias against scientific #s)
options(stringsAsFactors = F,
scipen = 6)

## Libraries ----
source("packages.r")

## Code ----
scripts <- list.files(
  "./code",
  full.names = T
)

sapply(
  scripts,
  source
)


#### PLANS ---------------------------------------------------------------------

## Format sequence data ----
sequencePlan <- drake_plan(
  ## Load and split raw data ----
  rawData = load_raw_data(),
  project = "MT",
  subData = split_raw_data(
    input = rawData,
    project = project
  )
)


#### MAKE ----------------------------------------------------------------------

## Bind plans ----
thePlan <- bind_rows(
  sequencePlan
)

## Make ----
make(thePlan)

## Visualise targets ----
vis_drake_graph(thePlan)