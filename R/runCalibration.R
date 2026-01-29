############################################################
## Generic calibration engine (backend-agnostic)
############################################################

#' Run posterior predictive calibration
#'
#' @param MCMCSamples Matrix of posterior draws from the observed-data fit.
#' @param observedData Observed dataset (any R object).
#' @param MCMCFun Function `function(newData, control)` that runs a short MCMC  on `newData` and returns posterior samples.
#' @param simulateNewDataFun Function `function(thetaRow, observedData, control)` that  simulates one replicated dataset from the posterior predictive. SP: We assume that new data is sampled from the posterior predictive of the model. In principle we may want to consider sampling from the prior predictive.
#' @param discFun Function `function(MCMCSamples, data, control)` that returns  a list with at least a numeric vector `D` of discrepancy values, one per row
#'   of `MCMCSamples`.
#' @param nReps Number of calibration replications.
#' @param drawIndexSelector Optional function `function(MCMCSamples, nReps, control)`
#'   returning the indices of rows to use as seeds for calibration.
#' @param control List of additional backend-specific arguments.
#' @param ... Not used currently.
#'
#' @return a list (future `S3` class `cpppResult` objects) containing the CPPP, observed and replicated repPPP, observed discrepancies and replicated discrepancies.
#' @export

runCalibration <- function(
    MCMCSamples,
    observedData,
    MCMCFun,
    simulateNewDataFun,
    discFun,
    nReps,
    drawIndexSelector = NULL,
    control = list(),
    ...
) {
  MCMCSamples <- as.matrix(MCMCSamples)
  nDraws      <- nrow(MCMCSamples)

  if (nDraws < 1L) stop("MCMCSamples must contain at least one row.")

  ##  Choose rows in MCMCSamples to simulate data for calibration
  if (is.null(drawIndexSelector)) {
    drawnIndices <- floor(seq(1, nDraws, length.out = nReps))
  } else {
    drawnIndices <- drawIndexSelector(MCMCSamples, nReps, control)
  }
  if (length(drawnIndices) != nReps) {
    stop("drawIndexSelector must return exactly nReps indices.")
  }

  ## 2. Discrepancies + PPP for the observed data
  obsDisc <- discFun(MCMCSamples  = MCMCSamples,
                       newData     = observedData,
                       control     = control)
  if (!all(c("obs", "sim") %in% names(obsDisc))) {
    stop("discFun must return a list with components 'obs' and 'sim'.")
  }
  ## discrepancy computed on real data
  origObsDisc <- as.numeric(obsDisc$obs)
  origSimDisc <- as.numeric(obsDisc$sim)
  if (length(origObsDisc) != length(origSimDisc)) {
    stop("'obs' and 'sim' must have the same length.")
  }

  # scalar PPP for the observed data
  obsPPP <- mean(origSimDisc >= origObsDisc)

  ## 3. Calibration worlds: PPP for each replicates
  repPPP <- numeric(nReps)
  repDiscList <- vector("list", nReps)

  for (r in seq_len(nReps)) {
    thetaRow <- MCMCSamples[drawnIndices[r], , drop = FALSE]

    # 3a. simulate a new dataset y^(r) from posterior predictive of the original model
    newData <- simulateNewDataFun(thetaRow = thetaRow,
                             observedData = observedData,
                             control = control)

    # 3b. fit model on y^(r) (short chain)
    repMCMC <- MCMCFun(newData = newData,
                         control  = control)

    # 3c. compute discrepancies + PPP in this world
    repDisc <- discFun(MCMCSamples = repMCMC,
                         newData     = newData,
                         control      = control)

    if (!all(c("obs", "sim") %in% names(repDisc))) {
      stop("discFun must return a list with components 'obs' and 'sim'.")
    }
    repObsDisc <- as.numeric(repDisc$obs)
    repSimDisc <- as.numeric(repDisc$sim)
    if (length(repObsDisc) != length(repSimDisc)) {
      stop("'obs' and 'sim' must have the same length in each calibration world.")
    }

    repPPP[r]         <- mean(repSimDisc >= repObsDisc)
    repDiscList[[r]] <- repDisc

  }

  ## 4. CPPP: how extreme obsPPP is under the calibration distribution
  CPPP <- mean(repPPP <= obsPPP)

  ## 5. Collect all discrepancies
  discrepancies <- list(
    obs = list(
      obs = origObsDisc,
      sim = origSimDisc
    ),
    rep = lapply(repDiscList, function(d) {
      list(
        obs = as.numeric(d$obs),
        sim = as.numeric(d$sim)
      )
    })
  )

  ## 6. Return cpppResult object
  newCpppResult(
    CPPP         = CPPP,
    repPPP          = repPPP,
    obsPPP      = obsPPP,
    discrepancies = discrepancies,
    drawnIndices   = drawnIndices
  )
}




