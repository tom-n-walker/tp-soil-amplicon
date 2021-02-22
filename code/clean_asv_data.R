################################################################################
#### Project: TP Network
#### Title:   Function | Wrangle | Format ASV data
#### Author:  Tom Walker (thomas.walker@usys.ethz.ch)
#### Date:    9 December 2020
#### ---------------------------------------------------------------------------

clean_asv_data <- function(
  metadata, counts, taxonomy, 
  perc.cutoff = 0.0001, 
  prev.cutoff1 = 100,
  prev.cutoff2 = 0.1,
  count.cutoff = 0
){
  ## Generate phyloseq object ----
  # Add ASV information to count data rownames
  countTable <- counts %>%
    mutate(sample_name = taxonomy$V1) %>%
    column_to_rownames("sample_name") %>%
    otu_table(., taxa_are_rows = T)
  # Move ASV names to rownames taxonomy
  taxoTable <- taxonomy %>%
    column_to_rownames("V1") %>%
    as.matrix %>%
    tax_table(.)
  # Make sample_code rownames metadata
  metaTable <- metadata %>%
    column_to_rownames("sample_code") %>%
    sample_data(.)
  # Combine into phyloseq object
  rawSeq <- phyloseq(countTable, taxoTable, metaTable)
  
  ## Basic formatting ----
  # filter against poor read depth
  allSeq <- prune_samples(
    samples = sample_sums(rawSeq) > count.cutoff,
    x = rawSeq
  )
  
  ## Filtering ----
  # Percent cut-off
  percASV <- phyloseq_filter_sample_wise_abund_trim(
    physeq = allSeq,
    relabund = T,
    minabund = perc.cutoff
  )
  # Prevalence cut-off
  prevASV <- phyloseq_filter_prevalence(
    allSeq, 
    prev.trh = prev.cutoff2, 
    abund.trh = prev.cutoff1, 
    threshold_condition = "OR"
  )
  
  ## Standardisation ----
  # Relative abundance
  rawRA <- phyloseq_standardize_otu_abundance(
    physeq = allSeq,
    method = "total"
  )
  percRA <- phyloseq_standardize_otu_abundance(
    physeq = percASV,
    method = "total"
  )
  prevRA <- phyloseq_standardize_otu_abundance(
    physeq = prevASV,
    method = "total"
  )
  # GMPR (my function)
  rawGMPR <- apply_gmpr(allSeq)
  percGMPR <- apply_gmpr(percASV)
  prevGMPR <- apply_gmpr(prevASV)

  ## Generate output ----
  out <- list(
    rawCounts = allSeq,
    rawRA = rawRA,
    rawGMPR = rawGMPR,
    percRA = percRA,
    percGMPR = percGMPR,
    prevRA = prevRA,
    prevGMPR = prevGMPR
  )
  return(out)
}
