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
  "./processing_code/",
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
  # Subset ASV data for this project code
  subSeqData = split_raw_data(
    input = seqData,
    project = "MT"
  )
)
formatPlan <- drake_plan(
  # Sample and metadata
  sampleData = compile_sample_data(
    metadata = subSeqData$metadata,
    climData = climData
  ),
  siteData = compile_site_data(
    climData = climData
  ),
  # Filter and normalise ASV data
  bacFull = clean_asv_data(
    metadata = sampleData,
    counts = subSeqData$counts$bacteria,
    taxonomy = subSeqData$taxonomy$bacteria,
    perc.cutoff = 0.0001, 
    prev.cutoff1 = 100,
    prev.cutoff2 = 0.1,
    count.cutoff = 2500
  ),
  funFull = clean_asv_data(
    metadata = sampleData,
    counts = subSeqData$counts$fungi,
    taxonomy = subSeqData$taxonomy$fungi,
    perc.cutoff = 0.0001, 
    prev.cutoff1 = 100,
    prev.cutoff2 = 0.1,
    count.cutoff = 0
  )
)

## QC on ASV data ----
qcPlan <- drake_plan(
  bacFilterQC = plot_filter_bac(
    seqData = bacFull
  ),
  funFilterQC = plot_filter_fun(
    seqData = funFull
  )
)

## Analysis ----
analysePlan <- drake_plan(
  # Select filtering option
  fungi = funFull$percGMPR,
  bacteria = bacFull$percGMPR,
  # Subset for simplicity
  bacSub = tax_glom(bacteria, taxrank = "Phylum"),
  funSub = tax_glom(fungi, taxrank = "Class"),
  # NMDS on full data
  funFullNMDS = sitewise_nmds(seqData = fungi),
  bacFullNMDS = sitewise_nmds(seqData = bacteria),
  # NMDS on subset data
  funSubNMDS = sitewise_nmds(seqData = funSub),
  bacSubNMDS = sitewise_nmds(seqData = bacSub),
  # Collate data
  finalDF = collate_final(
    siteData = siteData,
    sampleData = sampleData,
    soilData = soilData,
    fungi = fungi,
    funFullNMDS = funFullNMDS,
    funSubNMDS = funSubNMDS,
    bacFullNMDS = bacFullNMDS,
    bacSubNMDS = bacSubNMDS
  )
)


#### MAKE ----------------------------------------------------------------------

## Bind plans ----
thePlan <- bind_rows(
  loadPlan,
  formatPlan,
  qcPlan,
  analysePlan
)

## Make ----
make(thePlan)

## Visualise targets ----
vis_drake_graph(thePlan, targets_only = T)
