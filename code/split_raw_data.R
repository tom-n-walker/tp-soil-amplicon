## -------------------------------------------------
##
## Title: Functions to split data between projects
##
## Author: Tom WN Walker
## Email: thomas.walker@usys.ethz.ch
##
## Created: 2018-04-11
##
## -------------------------------------------------


### Function to load data ------------------------------------------------------

load_data <- function () {
  # specify directories
  meta_dir <- paste0(getwd(), "/from_hannes/")
  data_dir <- paste0(meta_dir, "/processed files/")
  # load metadata
  meta <- fread(paste0(meta_dir, "metadata.csv"))
  # load count data - drop row names (sample ID)
  bac_nums <- fread(file = paste0(data_dir, "MT_16S_matrix.csv"), 
                    drop = 1)
  fun_nums <- fread(file = paste0(data_dir, "MT_ITS1_countmatrix.csv"), 
                    drop = 1)
  # load taxonomy data - specify header (sample ID column no header)
  bac_tax <- fread(file = paste0(data_dir, "MT_16S_taxonomy.csv"), 
                   header = T)
  fun_tax <- fread(file = paste0(data_dir, "MT_ITS1_taxonomy.csv"), 
                   header = T)
  # construct list and return
  out <- list(metadata = meta,
              counts = list(bacteria = bac_nums,
                            fungi = fun_nums),
              taxonomy = list(bacteria = bac_tax,
                              fungi = fun_tax))
  return(out)
}


### Function to split data between projects ------------------------------------

split_data <- function (input, labels = c("MT", "F", "S")) {
  # split count data by label names, subset and return
  my_split <- lapply(labels, function (x) {
    # index (positions shared by all numeric and metadata)
    position <- grep(pattern = x, 
                     x = colnames(input$counts$bacteria))
    # subset count data using this index (cols)
    bac_sub <- input$counts$bacteria %>% 
      select(all_of(position))
    fun_sub <- input$counts$fungi %>% 
      select(all_of(position))
    # subset metadata using this index (rows)
    meta_sub <- input$metadata %>%
      slice(position)
    # construct list and return
    out <- list(metadata = meta_sub,
                counts = list(bacteria = bac_sub,
                              ffungi = fun_sub),
                taxonomy = input$taxonomy)
    return(out)
  })
  # add labels as bin names
  names(my_split) <- labels
  # return split data
  return(my_split)
}


### Pull functions into convenience function -----------------------------------

load_split_raw <- function(){
  raw <- load_data()
  split <- split_data(input = raw)
  return(split)
}
