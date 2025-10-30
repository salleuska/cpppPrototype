#' Discrepancy extractor (online)
#'
#' Returns a `disc_fun` that reads already-computed discrepancy columns
#' (observed vs replicated) from the MCMC output and returns them for use
#' in `runCalibration()`, which computes the PPP.
#'
#' @param col_obs Character. Column name for observed discrepancies.
#' @param col_rep Character. Column name for replicated discrepancies.
#'
#' @return A function `function(new_samples, new_data, ...)` that returns
#'   `list(obs = <numeric>, rep = <numeric>)` with equal lengths.
#'
#' @export
make_col_disc_fun <- function(col_obs = "disc_obs",
                              col_rep = "disc_rep") {
  function(MCMC_samples, ...) {
    ## check that the discrepancy columns are in the samples
    if (!all(c(col_obs, col_rep) %in% colnames(new_samples))) {
      stop(sprintf("Columns '%s' and/or '%s' not found in new_samples.", col_obs, col_rep))
    }
    list(
      obs = as.numeric(MCMC_samples[, col_obs]),
      rep = as.numeric(MCMC_samples[, col_rep])
    )
  }
}

#' Make offline discrepancy function
#'
#' Creates a `disc_fun` that computes observed and replicated discrepancies
#' for each posterior draw. This is the "offline" version used when the
#' MCMC output does not already include discrepancy columns.
#'
#' The returned function expects to receive MCMC samples (`MCMC_samples`)
#' and the corresponding dataset (`new_data`). It uses the model-specific
#' functions provided in `control` to simulate replicated data and compute
#' the discrepancy values.
#'
#' @param control A list containing at least two functions:
#'   \describe{
#'     \item{new_data_fun}{A function `function(theta_row, ...)` that simulates one
#'       replicated dataset `y*` given a single parameter draw.}
#'     \item{discrepancy}{A function `function(data, theta_row, ...)` that returns
#'       the scalar discrepancy `D(data, θ)` for one draw.}
#'   }
#'
#' @return A function of the form `function(MCMC_samples, new_data, ...)`
#'   that returns a list with two numeric vectors:
#'   \itemize{
#'     \item `obs`: discrepancies for the observed (replicate) data.
#'     \item `rep`: discrepancies for the replicated data.
#'   }
#'
#' @examples
#' \dontrun{
#' control <- list(
#'   new_data_fun = function(theta) rnorm(10, theta),
#'   discrepancy  = function(data, theta) mean((data - theta)^2)
#' )
#' disc_fun <- make_offline_disc_fun(control)
#' }
#'
#' @export
make_offline_disc_fun <- function(control) {
  function(MCMC_samples, new_data, ...) {
    ## Offline discrepancy calculator: compute D(data, θ) and D(y*, θ)
    if (!is.list(control) || !all(c("new_data_fun", "discrepancy") %in% names(control))) {
      stop("control must be a list with elements 'new_data_fun' and 'discrepancy'.")
    }
    new_data_fun <- control$new_data_fun
    discrepancy  <- control$discrepancy

    d_obs <- apply(MCMC_samples, 1, function(th) discrepancy(new_data, th, ...))
    d_rep <- apply(MCMC_samples, 1, function(th) {
      y_rep <- new_data_fun(th, ...)
      discrepancy(y_rep, th, ...)
    })
    list(obs = as.numeric(d_obs), rep = as.numeric(d_rep))
  }
}
