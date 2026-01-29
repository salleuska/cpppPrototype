#' Discrepancy extractor (online)
#'
#' Returns a `discFun` that reads already-computed discrepancy columns
#' (observed vs replicated) from the MCMC output and returns them for use
#' in `runCalibration()`, which computes the PPP.
#'
#' @param colObs Character. Column name for observed discrepancies.
#' @param colSim Character. Column name for replicated discrepancies.
#'
#' @return A function `function(MCMCSamples, newData, ...)` that returns
#'   `list(obs = <numeric>, sim = <numeric>)` with equal lengths.
#'
#' @export

makeColDiscFun <- function(colObs = "discrepancy_model",
                             colSim = "discrepancy_simulated") {
  function(MCMCSamples, ...) {
    ## check that the discrepancy columns are in the samples
    if (!all(c(colObs, colSim) %in% colnames(MCMCSamples))) {
      stop(sprintf("Columns '%s' and/or '%s' not found in new_samples.", colObs, colSim))
    }
    list(
      obs = as.numeric(MCMCSamples[, colObs]),
      sim = as.numeric(MCMCSamples[, colSim])
    )
  }
}

#' Make offline discrepancy function
#'
#' Creates a `discFun` that computes observed and replicated discrepancies
#' for each posterior draw. This is the "offline" version used when the
#' MCMC output does not already include discrepancy columns.
#'
#' The returned function evaluates discrepancies conditional on `targetData`,
#' which may be the original observed dataset or a replicated dataset from a
#' calibration world. It expects posterior samples (`MCMCSamples`) and the
#' corresponding `targetData`, and uses the model-specific functions in
#' `control` to generate posterior predictive replicates and compute the
#' discrepancy values.
#'
#' @param control A list containing at least two functions:
#'   \describe{
#'     \item{simulateNewDataFun}{A function `function(thetaRow, control)` that simulates one
#'       replicated dataset `y*` given one draws of the parameters}
#'     \item{discrepancy}{A function `function(data, thetaRow, ...)` that returns
#'       the scalar discrepancy `D(data, θ)` for one draw.}
#'   }
#'
#' @return A function of the form `function(MCMCSamples, targetData, ...)`
#'   that returns a list with two numeric vectors:
#'   \itemize{
#'     \item `obs`: discrepancies for the observed (replicate) data.
#'     \item `sim`: discrepancies for the replicated data.
#'   }
#'
#' @examples
#' \dontrun{
#' control <- list(
#'   simulateNewDataFun = function(theta) rnorm(10, theta),
#'   discrepancy  = function(data, theta) mean((data - theta)^2)
#' )
#' discFun <- makeOfflineDiscFun(control)
#' }
#'
#' @export

makeOfflineDiscFun <- function(control) {
  function(MCMCSamples, targetData, ...) {
    ## Offline discrepancy calculator: compute D(data, theta) and D(y*, theta)
    if (!is.list(control) || !all(c("simulateNewDataFun", "discrepancy") %in% names(control))) {
      stop("control must be a list with elements 'simulateNewDataFun' and 'discrepancy'.")
    }
    simulateNewDataFun <- control$simulateNewDataFun
    discrepancy  <- control$discrepancy

    dObs <- apply(MCMCSamples, 1, function(thetaRow) discrepancy(targetData, thetaRow, ...))
    dSim <- apply(MCMCSamples, 1, function(thetaRow) {
      y <- simulateNewDataFun(thetaRow, control)
      discrepancy(y, thetaRow, ...)
    })
    list(obs = dObs, sim = dSim)
  }
}
