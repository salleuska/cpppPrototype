#' Create a cpppResult object
#'
#' Constructs an S3 object of class `cpppResult`, which stores the results
#' related to the calibration of the predictive p-value.
#'
#' @param cppp Numeric scalar. The approximated calibrated posterior predictive p-value.
#' @param ppp Numeric vector. The posterior predictive p-values across calibration datasets.
#' @param obs_ppp Numeric vector. The observed posterior predictive p-values for the observed data.
#
#' @param discrepancies Optional list or array containing discrepancy values for
#'   each replication and dataset (observed and simulated).

#' @param ... Additional elements to store in the object.
#' SP E.g., do we want to save in this object functions used in the pipelie?
#'
#' @return An object of class `cppresults`.
#' @export
new_cppresults <- function(cppp = NA_real_,
                           ppp = numeric(),
                           obs_ppp = numeric(),
                           discrepancies = NULL,
                           ...) {
  x <- list(
    cppp = cppp,
    ppp = ppp,
    obs_ppp = obs_ppp,
    discrepancies = discrepancies,
    ...
  )
  class(x) <- "cpppResult"
  x
}
