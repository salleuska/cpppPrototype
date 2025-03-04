#' function that run the calibration and compute the cppp samples from the posterior predictive
#' and calculate the cppp
#'
#' @param data: original data
#' @param nCalibrationReplicates: number of calibration replicates
#' @param MCMCcontrol: list with number of posterior samples per replicates (and optional thinning and burnin)
#' @param sampleFromPosterior: function that can sample from the posterior distribution
#' @param sampleFromPP: function that can sample new data from the posterior predictive using some provided samples
#' @param discrepancy: a function that calculates the discrepancy (statistics) for the data (original or sampled)
#' @param ...: additional arguments to be passed
#'
#'

calculateCPPP <- function(data,
                          nCalibrationReplicates = 100,
                          MCMCcontrol = list(niter = 200,
                                             thin = 1,
                                             nburnin = 0),
                          sampleFromPosterior,
                          sampleFromPP, discrepancy, ...) {

  ## calculatePPP on the original data if not provided

  ## run the calibration
  calibration <- runCalibration(nCalibrationReplicates, MCMCcontrol, sampleFromPosterior, ...)

  ## calculate the cppp

  ## return the cppp
  cppp
}

runCalibration <- function(nCalibrationReplicates, MCMCcontrol, sampleFromPosterior, ...) {
  ## run the calibration:
  ## 1. sample new data from the posterior predictive
  ## 2. treat the new data as the observed data and get the posterior distribution
  ## 3. repeat the above steps nCalibrationReplicates times

  ## initialize the list to store the results
  calibration <- list()

  ## loop over the replicates
  for(i in 1:nCalibrationReplicates) {
    ## sample new data from the posterior predictive
    newData <- sampleFromPP(1, ...)

    ## treat the new data as the observed data and get the posterior distribution
    posterior <- sampleFromPosterior(newData, ...)

    ## store the results
    calibration[[i]] <- list(data = newData, samples = posterior)
  }

  ## return the results
  calibration
}
