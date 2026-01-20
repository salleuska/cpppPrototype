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

make_col_discFun <- function(colObs = "discrepancy_model",
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
#' The returned function expects to receive MCMC samples (`MCMCSamples`)
#' and the corresponding dataset (`newData`). It uses the model-specific
#' functions provided in `control` to simulate replicated data and compute
#' the discrepancy values.
#'
#' @param control A list containing at least two functions:
#'   \describe{
#'     \item{simulateNewDataFun}{A function `function(theta_row, ...)` that simulates one
#'       replicated dataset `y*` given a single parameter draw.}
#'     \item{discrepancy}{A function `function(data, theta_row, ...)` that returns
#'       the scalar discrepancy `D(data, θ)` for one draw.}
#'   }
#'
#' @return A function of the form `function(MCMCSamples, newData, ...)`
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
#' discFun <- make_offline_discFun(control)
#' }
#'
#' @export

make_offline_discFun <- function(control) {
  function(MCMCSamples, newData, ...) {
    ## Offline discrepancy calculator: compute D(data, theta) and D(y*, theta)
    if (!is.list(control) || !all(c("simulateNewDataFun", "discrepancy") %in% names(control))) {
      stop("control must be a list with elements 'simulateNewDataFun' and 'discrepancy'.")
    }
    simulateNewDataFun <- control$simulateNewDataFun
    discrepancy  <- control$discrepancy

    dObs <- apply(MCMCSamples, 1, function(th) discrepancy(newData, th, ...))
    dSim <- apply(MCMCSamples, 1, function(th) {
      y <- simulateNewDataFun(th, ...)
      discrepancy(y, th, ...)
    })
    list(obs = dObs, sim = dSim)
  }
}
