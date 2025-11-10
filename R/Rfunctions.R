# R/Rfunctions.R

#' Compute calibrated posterior predictive p-value (placeholder)
#'
#' @param pppObs Numeric scalar in [0,1]; observed PPP.
#' @param pppRep Numeric vector in [0,1]; calibration PPPs.
#' @return Numeric scalar; calibrated PPP.
#' @export
computeCppp <- function(pppObs, pppRep) {
  stop("Not implemented", call. = FALSE)
}


#' Estimate transfer autocorrelation (placeholder)
#'
#' Placeholder for the function that will estimate the transfer
#' autocorrelation from the observed Delta-chain
#'
#' @param deltaChain Data or structure containing \eqn{D(y^{rep}, \theta) - D(y^{obs}, \theta)} values.
#' @param pppObs Numeric scalar in [0,1].
#' @param pppRep Numeric vector in [0,1].
#' @param mTilde Numeric vector or scalar; effective Monte Carlo sample size per replication.
#' @param ... Additional arguments (future use).
#'
#' @return List with components to be defined (e.g., `iact`, `essTransfer`).
#' @export
transferAutocorrelation <- function(deltaChain, pppObs, pppRep, mTilde, ...) {
  stop("Not implemented", call. = FALSE)
}
