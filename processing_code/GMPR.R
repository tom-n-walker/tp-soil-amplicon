################################################################################
#### Project: TP Network
#### Title:   Function | Small function | GMPR function
#### Author:  https://github.com/jchen1981/GMPR
#### Date:    10 December 2020
#### ---------------------------------------------------------------------------


GMPR <- function(comm, intersect.no = 5, ct.min = 1, trace = TRUE){
  # Args:
  #   comm: a matrix of counts, row - OTUs, column - sample
  #   intersect.no: min. no. shared OTUs between pairs (ratio is calculated)
  #   ct.min: min. no. counts required to calculate ratios
  #
  # Returns:
  #   vector of size factors with attribute 'NSS'. Samples with distinct sets 
  #   of features will be given NA.
  #   NSS: no. samples with significant sharing (> intersect.no) including self
  
  # mask counts < ct.min
  comm[comm < ct.min] <- 0
  # make colum names if necessary
  if (is.null(colnames(comm))) {
    colnames(comm) <- paste0('S', 1:ncol(comm))
  }
  # verbose output
  if (trace) cat('Begin GMPR size factor calculation ...\n')
  # number of samples
  comm.no <- numeric(ncol(comm))
  # apply gmpr function
  gmpr <- sapply(
    1:ncol(comm), 
    function(i){		
      if (i %% 50 == 0) {
        cat(i, '\n')
      }
      x <- comm[, i]
      # Compute the pairwise ratio
      pr <- x / comm
      # Handling of the NA, NaN, Inf
      pr[is.nan(pr) | !is.finite(pr) | pr == 0] <- NA
      # Counting the number of non-NA, NaN, Inf
      incl.no <- colSums(!is.na(pr))		
      # Calculate the median of PR
      pr.median <- colMedians(pr, na.rm=TRUE)
      # Record the number of samples used for calculating the GMPR
      comm.no[i] <<- sum(incl.no >= intersect.no)
      # Geometric mean of PR median
      if (comm.no[i] > 1) {
        return(exp(mean(log(pr.median[incl.no >= intersect.no]))))
      } else {
        return(NA)
      }
    }
  )
  #  messages
  if (sum(is.na(gmpr))) {
    warning(paste0(
      'The following samples\n ', 
      paste(colnames(comm)[is.na(gmpr)], collapse='\n'), 
      '\ndo not share at least ', intersect.no, 
      ' common taxa with the rest samples! ',
      'For these samples, their size factors are set to be NA! \n', 
      'You may consider removing these samples',
      'since they are potentially outliers or negative controls!\n',
      'You may also consider decreasing the minimum',
      'number of intersecting taxa and rerun the procedure!\n'
    )
    )
  }
  if (trace) cat('Completed!\n')
  if (trace) cat(
    'Please watch for the samples with limited sharing with',
    'other samples based on NSS! They may be outliers! \n'
  )
  # give names and attributes
  names(gmpr) <- names(comm.no) <- colnames(comm)
  attr(gmpr, 'NSS') <- comm.no
  # return
  return(gmpr)
}
