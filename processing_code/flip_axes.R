################################################################################
#### Project: TP Network
#### Title:   Small function | flip PC axes if negative
#### Author:  Tom Walker (thomas.walker@usys.ethz.ch)
#### Date:    22 February 2021
#### ---------------------------------------------------------------------------

flip_axes <- function(x){
  if(x[2, 3] < 0){
    x$value <- x$value * -1
  }
  return(x)
}
