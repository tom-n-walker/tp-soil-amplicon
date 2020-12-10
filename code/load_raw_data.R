################################################################################
#### Project: TP Network
#### Title:   Function | Wrangle | Load raw sequence data
#### Author:  Tom Walker (thomas.walker@usys.ethz.ch)
#### Date:    03 November 2020
#### ---------------------------------------------------------------------------

load_raw_data <- function(){
  # project directory
  projectDir <- paste0(
    "/Users/tomwalker/Dropbox/projects/2018_transplant/",
    "SNF_network/data/soil_amplicon/processing/"
  )
  # Load and format metadata ----
  meta <- fread(
    file = paste0(
      projectDir, 
      "from_hannes/metadata.csv"
    ),
    data.table = F
  )
  ## Load count data ----
  bac_nums <- fread(
    file = paste0(
      projectDir, 
      "from_hannes/processed_files/MT_16S_matrix.csv"
    ),
    drop = 1,
    data.table = F
  )
  fun_nums <- fread(
    file = paste0(
      projectDir, 
      "from_hannes/processed_files/MT_ITS1_countmatrix.csv"
    ),
    drop = 1,
    data.table = F
  )
  ## Load taxonomy data ----
  bac_tax <- fread(
    file = paste0(
      projectDir,
      "from_hannes/processed_files/MT_16S_taxonomy.csv"
    ), 
    header = T,
    data.table = F
  )
  fun_tax <- fread(
    file = paste0(
      projectDir,
      "from_hannes/processed_files/MT_ITS1_taxonomy.csv"
    ), 
    header = T,
    data.table = F
  )
  ## Construct list and return ----
  out <- list(
    metadata = meta,
    counts = list(bacteria = bac_nums,
                  fungi = fun_nums),
    taxonomy = list(bacteria = bac_tax,
                    fungi = fun_tax)
  )
  return(out)
}