################################################################################
#### Project: TP Network
#### Title:   Small Function | Generate plot data for biplots
#### Author:  Tom Walker (thomas.walker@usys.ethz.ch)
#### Date:    24 November 2020
#### ---------------------------------------------------------------------------

make_biplot_data <- function(PC1, PC2, groupVar){
  ## Bind data ----
  allData <- data.frame(PC1, PC2, groupVar)
  xyData <- cbind(meanPC1 = PC1, meanPC2 = PC2)
  ## Generate centroids ----
  centroids <- aggregate(
    xyData ~ groupVar,
    allData,
    mean
  )
  ## Merge with raw data ----
  out <- merge(
    allData,
    centroids,
    by = "groupVar"
  )
  ## Return ----
  return(out)
}
