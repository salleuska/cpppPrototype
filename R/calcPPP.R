#' Calculate
#'
#' simulate samples from the posterior predictive and calculate summaries
#'
#' @param modelInfo
#'
#'
#' @return ??
#'

## matrix of posterior samples
samples <- matrix(c(1,2,1,3), 2,2)

## function that can sample from the posterior predictive
## number of samples
## other arguments
## output is a list containing each posterior predictive sample
posterior_predictive <- function(ndraws, ... ) {
  output <- list()
  for(i in 1:ndraws) output[[i]] <- rnorm(10)
  output
}
