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

## Load and format data ----
formatPlan <- drake_plan(
  ## Load and split raw data ----
  seqData = load_raw_data(),
  climData = load_clim_data(),
  project = "MT",
  subSeqData = split_raw_data(
    input = seqData,
    project = project
  ),
  soilData = load_soil_data(),
  ## Format data ----
  cleanData = clean_all_data(
    subSeqData = subSeqData,
    climData = climData,
    soilData = soilData
  )
)

## Format other data ----


#### MAKE ----------------------------------------------------------------------

## Bind plans ----
thePlan <- bind_rows(
  formatPlan
)

## Make ----
make(thePlan)

## Visualise targets ----
vis_drake_graph(thePlan)