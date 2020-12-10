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
loadPlan <- drake_plan(
  # Load all data
  seqData = load_raw_data(),
  climData = load_clim_data(),
  soilData = load_soil_data(),
  # Subset ASV data
  subSeqData = split_raw_data(
    input = seqData,
    project = "MT"
  )
)
formatPlan <- drake_plan(
  # Sample and metadata
  sampleData = compile_sample_data(
    metadata = subSeqData$metadata,
    climData = climData,
    soilData = soilData
  ),
  siteData = compile_site_data(
    metadata = subSeqData$metadata,
    climData = climData
  ),
  # Filter and normalise ASV data
  bacClean = clean_asv_data(
    metadata = sampleData,
    counts = subSeqData$counts$bacteria,
    taxonomy = subSeqData$taxonomy$bacteria,
    prev_cutoff1 = 2,
    prev_cutoff2 = 0.1,
    perc_cutoff = 0.01,
    count_cutoff = 2500
  ),
  funClean = clean_asv_data(
    metadata = sampleData,
    counts = subSeqData$counts$fungi,
    taxonomy = subSeqData$taxonomy$fungi,
    prev_cutoff1 = 2,
    prev_cutoff2 = 0.1,
    perc_cutoff = 0.01,
    count_cutoff = 0
  )
)

## Format other data ----


#### MAKE ----------------------------------------------------------------------

## Bind plans ----
thePlan <- bind_rows(
  loadPlan,
  formatPlan
)

## Make ----
make(thePlan)

## Visualise targets ----
vis_drake_graph(thePlan)