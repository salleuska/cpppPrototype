###################################################
#' Run calibration
#'
#' This function runs the calibration to obtain the calibrated posterior predictive p-value.
#'
#' @param MCMC_samples A matrix or list containing posterior samples from the abserved-data (long) run. Passed to `new_data_fun()` as context.
#' @param MCMC_fun A function of the form `function(new_data, control, ...)` that runs an MCMC algorithm and returns posterior samples.
#' @param new_data_fun A function of the form `function(MCMC_samples, ...)` that generates a new synthetic dataset for each calibration replicate. SP: We assume that new data is sampled from the posterior predictive of the model. In principle we may want to consider sampling from the prior predictive.
#' @param disc_fun A function of the form `function(mcmc_samples, ...)`
#'   that computes and returns a list with components `rep` and `obs`,
#'   the replicated and observed discrepancies for that replicate.
#' @param num_reps Integer. Number of calibration replicates to run.
#' @param control Optional list of controls passed to `MCMC_fun`
#'   (e.g., `niter`, seeds, etc.).
#' @param ... Additional arguments passed through to the user-supplied functions.
#' @details
#' This function is fully engine-agnostic and contains no NIMBLE-specific code. It provides the scaffolding for the cppp calibration procedure.
#'
#' @return Numeric vector of length `num_reps` containing the posterior predictive p-values for each calibration replicate.
#' @export
runCalibration <- function(MCMC_samples,
                           MCMC_fun,
                           new_data_fun,
                           disc_fun,
                           num_reps,
                           control = NULL,
                           ...) {
  ## No nimble-specific concepts or code in this function.
  ## Ignore where the "main" run happens

  if (!is.numeric(num_reps) || length(num_reps) != 1L || num_reps < 1L) {
    stop("`num_reps` must be a positive integer.", call. = FALSE)
  }
  num_reps <- as.integer(num_reps)

  ## SP: here may be where we want the parallel option to work
  ppp <- numeric(num_reps)

  for (i in seq_len(num_reps)) {
    new_data    <- new_data_fun(MCMC_samples, ...)
    new_samples <- MCMC_fun(new_data, control, ...)
    ## disc_fun() will return replicated discrepancies and observed one
    discrepancies <- disc_fun(new_samples, ...)

    if (!is.list(discrepancies) ||
        is.null(discrepancies$rep) ||
        is.null(discrepancies$obs)) {
      stop("`disc_fun` must return a list with components `rep` and `obs`.",
           call. = FALSE)
    }

    ppp[i] <- mean(discrepancies$rep >= discrepancies$obs, na.rm = TRUE)
  }

  ## return vector of ppp
  ppp
}
