################################################################################
#### Project: TP Network
#### Title:   Function | Wrangle | Apply GMPR
#### Author:  Tom Walker (thomas.walker@usys.ethz.ch)
#### Date:    10 December 2020
#### ---------------------------------------------------------------------------

apply_gmpr <- function(physeq){
  # Calculate GMPR correction (otu-wise)
  gmprFactor <- GMPR(physeq@otu_table@.Data) # my function
  # Transform otu table
  okCounts <- apply(
    physeq@otu_table, 
    1, 
    function(x) x * gmprFactor
  ) %>% t
  # Rebuild phyloseq object
  out <- phyloseq(
    otu_table(okCounts, taxa_are_rows = T),
    physeq@tax_table,
    physeq@sam_data
  )
  # Return
  return(out)
}