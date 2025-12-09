#' Create a cpppResult object
#'
#' Constructs an S3 object of class `cpppResult`, which stores the results
#' related to the calibration of the predictive p-value.
#'
#' @param cppp Numeric scalar. The approximated calibrated posterior predictive p-value.
#' @param ppp Numeric vector. The posterior predictive p-values across calibration datasets.
#' @param obs_ppp Numeric vector. The observed posterior predictive p-values for the observed data.
#'
#' @param discrepancies Optional list or array containing discrepancy values for
#'   each replication and dataset (observed and simulated).
#'
#' @param ... Additional elements to store in the object.
#' SP E.g., do we want to save in this object functions used in the pipeline?
#'
#' @return An object of class `cpppResults`.
#' @export
new_cpppResults <- function(cppp = NA_real_,
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
  class(x) <- c("cpppResult", "list")
  validate_cpppResult(x)
}

# Minimal internal validator: keeps inputs sane but stays permissive

#' @keywords internal
validate_cpppResult <- function(x) {
  stopifnot(inherits(x, "cpppResult"))

  # cppp: allow NA for now; if finite, must be a scalar in [0,1]
  if (!(length(x$cppp) == 1L && is.numeric(x$cppp))) {
    stop("`cppp` must be a numeric scalar (NA allowed).", call. = FALSE)
  }
  if (is.finite(x$cppp) && (x$cppp < 0 || x$cppp > 1)) {
    stop("`cppp` must be in [0,1] when finite.", call. = FALSE)
  }

  # ppp: numeric vector, all finite in [0,1] if provided
  if (!is.numeric(x$ppp)) {
    stop("`ppp` must be a numeric vector.", call. = FALSE)
  }
  if (length(x$ppp) && (!all(is.finite(x$ppp)) || any(x$ppp < 0 | x$ppp > 1))) {
    stop("All `ppp` values must be finite and in [0,1].", call. = FALSE)
  }

  # obs_ppp: numeric vector, all finite in [0,1] if provided
  if (!is.numeric(x$obs_ppp)) {
    stop("`obs_ppp` must be a numeric vector.", call. = FALSE)
  }
  if (length(x$obs_ppp) && (!all(is.finite(x$obs_ppp)) || any(x$obs_ppp < 0 | x$obs_ppp > 1))) {
    stop("All `obs_ppp` values must be finite and in [0,1].", call. = FALSE)
  }

  # discrepancies: optional; no strict contract yet
  # (kept permissive to allow list/array/tibble as you iterate)

  invisible(x)
}
