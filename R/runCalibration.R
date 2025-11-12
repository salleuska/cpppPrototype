############################################################
## Generic calibration engine (backend-agnostic)
############################################################

#' Run posterior predictive calibration
#'
#' @param MCMC_samples Matrix/data.frame of posterior draws from the observed-data fit.
#' @param observed_data Observed dataset (any R object).
#' @param MCMC_fun Function `function(new_data, control)` that runs a short MCMC  on `new_data` and returns posterior samples.
#' @param new_data_fun Function `function(theta_row, observed_data, control)` that  simulates one replicated dataset from the posterior predictive. SP: We assume that new data is sampled from the posterior predictive of the model. In principle we may want to consider sampling from the prior predictive.
#' @param disc_fun Function `function(MCMC_samples, data, control)` that returns  a list with at least a numeric vector `D` of discrepancy values, one per row
#'   of `MCMC_samples`.
#' @param n_reps Number of calibration replications.
#' @param row_selector Optional function `function(MCMC_samples, n_reps, control)`
#'   returning the indices of rows to use as seeds for calibration.
#' @param control List of additional backend-specific arguments.
#' @param ... Not used currently.
#'
#' @return A list (later your cpppResults) with observed discrepancies and replicated discrepancies / PPP skeleton.
runCalibration <- function(
    MCMC_samples,
    observed_data,
    MCMC_fun,
    new_data_fun,
    disc_fun,
    n_reps,
    row_selector = NULL,
    control = list(),
    ...
) {

  MCMC_samples <- as.matrix(MCMC_samples)
  n_draws      <- nrow(MCMC_samples)

  ## Check that MCMC samples has at leaset one row
  if (n_draws < 1L) {
    stop("MCMC_samples must contain at least one row.")
  }

  ## 1. Choose which posterior draws seed the calibration replications
  if (is.null(row_selector)) {
    # Default: evenly spaced draws across the chain
    row_indices <- floor(seq(1, n_draws, length.out = n_reps))
  } else {
    row_indices <- row_selector(MCMC_samples, n_reps, control)
  }

  if (length(row_indices) != n_reps) {
    stop("row_selector must return exactly n_reps indices.")
  }

  ## 2. Discrepancy for observed data
  ## Note - this function either calculate the discrepancy or
  ## extract the right MCMC output
  obs_disc <- disc_fun(MCMC_samples = MCMC_samples,
                       data         = observed_data,
                       control      = control)

  if (is.null(obs_disc$D)) {
    stop("disc_fun must return a list with numeric vector 'D'.")
  }
  D_obs <- as.numeric(obs_disc$D)
  n_disc <- length(D_obs)

  ## Storage for replicated discrepancies
  D_rep <- matrix(NA_real_, nrow = n_reps, ncol = n_disc)

  ## 3. Calibration replications
  for (r in seq_len(n_reps)) {
    theta_row <- MCMC_samples[row_indices[r], , drop = FALSE]

    # 3a. Simulate replicated data from posterior predictive
    new_data <- new_data_fun(theta_row = theta_row,
                             observed_data = observed_data,
                             control = control)

    # 3b. Fit model on replicated data (short chain)
    MCMC_rep <- MCMC_fun(new_data = new_data,
                         control  = control)

    # 3c. Compute discrepancies for replicated dataset
    rep_disc <- disc_fun(MCMC_samples = MCMC_rep,
                         data         = new_data,
                         control      = control)

    D_rep[r, ] <- as.numeric(rep_disc$D)
  }

  colnames(D_rep) <- names(D_obs)

  ## 4. Placeholder: PPP / CPPP calculation
  # Example scalar PPP per discrepancy: P(D_rep >= D_obs)
  PPP_obs <- vapply(
    X   = seq_len(n_disc),
    FUN = function(j) mean(D_rep[, j] >= D_obs[j]),
    FUN.VALUE = numeric(1)
  )

  ## Example "calibrated" p-value placeholder: use PPP distribution vs PPP_obs
  # This is intentionally minimal; replace with your exact CPPP definition.
  CPPP <- vapply(
    X   = seq_len(n_disc),
    FUN = function(j) mean(PPP_obs[j] <= PPP_obs[j]),  # dummy; replace later
    FUN.VALUE = numeric(1)
  )

  ## 5. Return structure – later you can wrap this into cpppResults S3
  list(
    observed_discrepancy   = D_obs,
    replicated_discrepancy = D_rep,
    PPP_obs                = PPP_obs,
    CPPP                   = CPPP,
    row_indices            = row_indices
  )
}





