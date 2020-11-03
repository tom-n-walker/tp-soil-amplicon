################################################################################
#### Project: TP Network
#### Title:   Function | Wrangle | Load raw sequence data
#### Author:  Tom Walker (thomas.walker@usys.ethz.ch)
#### Date:    03 November 2020
#### ---------------------------------------------------------------------------

load_raw_data <- function(){
  # specify directories
  meta_dir <- paste0(getwd(), "/from_hannes/")
  data_dir <- paste0(meta_dir, "/processed files/")
  # load metadata
  meta <- fread(paste0(meta_dir, "metadata.csv"))
  # load count data - drop row names (sample ID)
  bac_nums <- fread(
    file = paste0(data_dir, "MT_16S_matrix.csv"), 
    drop = 1
  )
  fun_nums <- fread(
    file = paste0(data_dir, "MT_ITS1_countmatrix.csv"), 
    drop = 1
  )
  # load taxonomy data - specify header (sample ID column no header)
  bac_tax <- fread(
    file = paste0(data_dir, "MT_16S_taxonomy.csv"), 
    header = T
  )
  fun_tax <- fread(
    file = paste0(data_dir, "MT_ITS1_taxonomy.csv"), 
    header = T
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