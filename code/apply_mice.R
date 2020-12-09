################################################################################
#### Project: Arctic priming
#### Title:   Function | Small function | Apply MICE imputation
#### Author:  Tom Walker (thomas.walker@usys.ethz.ch)
#### Date:    15 July 2020
#### ---------------------------------------------------------------------------

apply_mice <- function(matrix, iter){
  # get column names
  cols <- colnames(matrix)
  colnames(matrix) <- paste0("C", 1:ncol(matrix))
  # apply mice imputation
  imp <- mice(
    matrix, 
    iter
  )
  # generate list and add each imputation round to a bin
  imp_list <- list()
  for(i in 1:iter){
    imp_list[[i]] <- as.matrix(
      complete(imp, i)
    )
  }
  # take means of these iterations
  out <- apply(
    simplify2array(imp_list), 
    1:2, 
    mean, 
    na.rm = T
  ) %>% as.data.frame
  colnames(out) <- cols
  # return
  return(out)
}