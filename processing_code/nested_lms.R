################################################################################
#### Project: TP Network
#### Title:   Small function | apply LMs to nested data frame
#### Author:  Tom Walker (thomas.walker@usys.ethz.ch)
#### Date:    22 February 2021
#### ---------------------------------------------------------------------------

nested_lms <- function(nd, response){
  # build
  out <- nd %>%
    # model
    mutate(lms = map(
      data, 
      ~lm(arm::rescale(.x[, response, drop = T]) ~ treat_diff, .x))
    ) %>%
    # get coefficients
    mutate(hhhl_coefs = map(lms, ~tidy(.x)[2, ])) %>%
    mutate(llhl_coefs = map(lms, ~tidy(.x)[3, ]))
  # return
  return(out)
}
