################################################################################
#### Project: TP Network
#### Title:   Function | Wrangle | Load raw sequence data
#### Author:  Tom Walker (thomas.walker@usys.ethz.ch)
#### Date:    03 November 2020
#### ---------------------------------------------------------------------------

load_raw_data <- function(){
  ## Load and format metadata ----
  # Load sample metadata
  meta <- fread(
    "./data/metadata.csv", 
    data.table = F
  )
  ## Load count data
  bac_nums <- fread(
    file = "./data/MT_16S_matrix.csv",
    drop = 1,
    data.table = F
  )
  fun_nums <- fread(
    file = "./data/MT_ITS1_countmatrix.csv",
    drop = 1,
    data.table = F
  )
  ## Load taxonomy data
  bac_tax <- fread(
    file = "./data/MT_16S_taxonomy.csv", 
    header = T,
    data.table = F
  )
  fun_tax <- fread(
    file = "./data/MT_ITS1_taxonomy.csv", 
    header = T,
    data.table = F
  )
  # construct list and return
  out <- list(
    metadata = meta,
    counts = list(bacteria = bac_nums,
                  fungi = fun_nums),
    taxonomy = list(bacteria = bac_tax,
                    fungi = fun_tax)
  )
  return(out)
}