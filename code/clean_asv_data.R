################################################################################
#### Project: TP Network
#### Title:   Function | Wrangle | Format ASV data
#### Author:  Tom Walker (thomas.walker@usys.ethz.ch)
#### Date:    9 December 2020
#### ---------------------------------------------------------------------------

clean_asv_data <- function(
  metadata, counts, taxonomy, 
  perc_cutoff = 0.01, 
  prev_cutoff1 = 3, 
  prev_cutoff2 = 0.2,
  count_cutoff = 0
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
  
  ## Filter data against poor sequence depth ----
  allSeq <- prune_samples(
    samples = sample_sums(rawSeq) > count_cutoff,
    x = rawSeq
  )
  
  ## Filter out rare taxa with relative abundance cutoff ----
  # Calculate relative abundances (for now, only for filtering)
  relAbunSeq <- transform_sample_counts(
    physeq = allSeq,
    fun = function(x) x / sum(x) * 100
  )
  # Use this to generate percent filter
  percFilter <- filter_taxa(
    physeq = relAbunSeq,
    # filter where rel abun of ASV (x) is more than [perc_cutoff] % in samples
    flist = function(x) sum(x) > perc_cutoff
  )
  # Apply percent filter to count data
  percASV <- prune_taxa(
    taxa = percFilter,
    x = allSeq
  )
  
  ## Filter out rare taxa with prevalance cut-off ----
  prevASV <- filter_taxa(
    physeq = allSeq, 
    # filter for ASVs with counts > [prev_cutoff1] in min. [prev_cutoff2] data
    flist = function(x) sum(x > prev_cutoff1) > prev_cutoff2 * length(x), 
    prune = TRUE
  )
  
  ## Normalize filtered datasets ----
  # Relative abundance
  percRelAbun <- transform_sample_counts(
    physeq = percASV,
    fun = function(x) x / sum(x) * 100
  )
  prevRelAbun <- transform_sample_counts(
    physeq = prevASV,
    fun = function(x) x / sum(x) * 100
  )
  # GMPR (my function)
  percGMPR <- apply_gmpr(percASV)
  prevGMPR <- apply_gmpr(prevASV)
  
  ## Generate output ----
  out <- list(
    raw = allSeq,
    percFilter = list(
      relAbun = percRelAbun,
      gmpr = percGMPR
    ),
    prevFilter = list(
      relAbun = prevRelAbun,
      gmpr = prevGMPR
    )
  )
  return(out)
}