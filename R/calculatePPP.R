#' function that simulate samples from the posterior predictive and calculate summaries
#'
#' @param data: original data
#' @param samples: matrix of posterior samples
#' @param ndraws: number of samples to simulate
#' @param sampleFromPP: function that can sample new data from the posterior predictive using the provided samples
#' @param discrepancy: a function that calculates the discrepancy (statistics) for the data (original or sampled)
#' @param ...: additional arguments to be passed
#'
#' @return list containing the p-value, sampled discrepancies and original one
#'

calculatePPP <- function(data, sample, ndraws, sampleFromPP, discrepancy, ...) {
  ## simulate samples from the posterior predictive
  ppSamples <- sampleFromPP(ndraws, samples)

  ## calculate the discrepancy for each sample
  discPP <- sapply(ppSamples, function(x) discrepancy(x))

  ## calculate the discrepancy for the original data
  discData <- discrepancy(data, ...)

  ## calculate the p-value
  pvalue <- sum(discPP - discData >= 0) / ndraws

  ## return the p-value, sampled discrepancies and original one
  list(pvalue = pvalue, sampledDisc = discPP, discData = discData)

}

####################
## Silly testing ##
####################
## matrix of posterior samples
## 2 samples for 2 variables
samples <- matrix(c(1,2,1,3), 2,2)

## function that can sample from the posterior predictive
## number of samples
## other arguments
## output is a list containing each posterior predictive sample
sampleFromPP <- function(ndraws, samples, ... ) {
  output <- list()
  for(i in 1:ndraws) output[[i]] <- rnorm(10)
  output
}

## function that calculates the discrepancy
## input is a vector of data
## output is a scalar
discrepancy <- function(data) {
  mean(data)
}

## test

data <- rnorm(10)
res <- calculatePPP(data, samples, 100, sampleFromPP, discrepancy)
hist(res$sampledDisc, breaks= 20)
abline(v = res$discData, col = "red")


#' function that run the calibration and compute the cppp samples from the posterior predictive
#' and calculate the cppp
#'
#' @param data: original data
#' @param nCalibrationReplicates: number of calibration replicates
#' @param MCMCcontrol: list with number of posterior samples per replicates (and optional thinning and burnin)
#' @param sampleFromPosterior: function that can sample from the posterior distribution
#' @param sampleFromPP: function that can sample new data from the posterior predictive using some provided samples
#' @param discrepancy: a function that calculates the discrepancy (statistics) for the data (original or sampled)
#' @param origMCMCSamples: matrix of posterior samples from the original data
#' @param ...: additional arguments to be passed
#'
#'

# calculateCPPP <- function(data,
#                           nCalibrationReplicates = 100,
#                           MCMCcontrol = list(niter = 200,
#                                              thin = 1,
#                                              nburnin = 0),
#                           sampleFromPosterior,
#                           sampleFromPP,
#                           discrepancy,
#                           origMCMCSamples,
#                           ...) {
#
#   #### Resumee here
#   ## calculatePPP on the original data if not provided
#   res <- calculatePPP(data, origMCMCSamples,ndraws = dim(origMCMCSamples)[1], sampleFromPP, discrepancy)
#
#
#   ## run the calibration
#   calibration <- runCalibration(nCalibrationReplicates, MCMCcontrol, sampleFromPosterior, ...)
#
#   ## calculate the cppp
#
#   ## return the cppp
#   cppp
# }

# runCalibration <- function(nCalibrationReplicates, MCMCcontrol, sampleFromPosterior, ...) {
#   ## run the calibration:
#   ## 1. sample new data from the posterior predictive
#   ## 2. treat the new data as the observed data and get the posterior distribution
#   ## 3. repeat the above steps nCalibrationReplicates times
#
#   ## initialize the list to store the results
#   calibration <- list()
#
#   ## loop over the replicates
#   for(i in 1:nCalibrationReplicates) {
#     ## sample new data from the posterior predictive
#     newData <- sampleFromPP(1, ...)
#
#     ## treat the new data as the observed data and get the posterior distribution
#     posterior <- sampleFromPosterior(newData, ...)
#
#     ## store the results
#     calibration[[i]] <- list(data = newData, samples = posterior)
#   }
#
#   ## return the results
#   calibration
# }

