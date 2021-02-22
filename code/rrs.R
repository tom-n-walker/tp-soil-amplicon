################################################################################
#### Project: TP Network
#### Title:   Small function | Apply RR calculation
#### Author:  Tom Walker (thomas.walker@usys.ethz.ch)
#### Date:    22 February 2021
#### ---------------------------------------------------------------------------

rrs <- function(x, variable){
  # set up data frames
  hh <- filter(x, treatment == "HH")
  hl <- filter(x, treatment == "HL")
  ll <- filter(x, treatment == "LL")
  # sample
  sampleHH <- sample(hh[, variable, drop = T], size = 100, replace = T)
  sampleHL <- sample(hl[, variable, drop = T], size = 100, replace = T)
  sampleLL <- sample(ll[, variable, drop = T], size = 100, replace = T)
  # calculate
  hlhh <- (sampleHL - sampleHH) / sampleHH
  hlll <- (sampleHL - sampleLL) / sampleLL
  # return
  out <- data.frame(
    hlhh_mean = mean(hlhh),
    hlhh_se = se(hlhh),
    hlll_mean = mean(hlll),
    hlll_se = se(hlll)
  )
  # return
  return(out)
}
