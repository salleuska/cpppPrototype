# R/Rfunctions.R

#' Compute the calibrated posterior predictive p-value
#'
#' @param pppObs Numeric scalar in \eqn{[0,1]}. The observed PPP.
#' @param pppCal Numeric vector in \eqn{[0,1]}. PPPs from calibration datasets.
#'
#' @return Numeric scalar in \eqn{[0,1]}, the calibrated PPP.
#' @export
computeCppp <- function(pppObs, pppCal) {
  # basic validation
  if (!is.numeric(pppObs) || length(pppObs) != 1L)
    stop("pppObs must be a numeric scalar.", call. = FALSE)
  if (!is.numeric(pppCal))
    stop("pppCal must be numeric.", call. = FALSE)

  # handle empty vector as NA ??
  if (length(pppCal) == 0L)
    return(NA_real_)

  mean(pppCal <= pppObs)
}


#' Estimate transfer autocorrelation (placeholder)
#'
#' Placeholder for the function that will estimate the transfer
#' autocorrelation from the observed Delta-chain
#'
#' @param deltaChain Data or structure containing \eqn{D(y^{rep}, \theta) - D(y^{obs}, \theta)} values.
#' @param pppObs Numeric scalar in  \eqn{[0,1]}.
#' @param pppRep Numeric vector in  \eqn{[0,1]}.
#' @param mTilde Numeric vector or scalar; effective Monte Carlo sample size per replication.
#' @param ... Additional arguments (future use).
#'
#' @return List with components to be defined (e.g., `iact`, `essTransfer`).
#' @export
transferAutocorrelation <- function(deltaChain, pppObs, pppRep, mTilde, ...) {
  stop("Not implemented", call. = FALSE)
}
