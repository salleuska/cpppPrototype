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
#' @export

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

  if (n_draws < 1L) stop("MCMC_samples must contain at least one row.")

  ## 1. Choose rows to seed calibration worlds
  if (is.null(row_selector)) {
    row_indices <- floor(seq(1, n_draws, length.out = n_reps))
  } else {
    row_indices <- row_selector(MCMC_samples, n_reps, control)
  }
  if (length(row_indices) != n_reps) {
    stop("row_selector must return exactly n_reps indices.")
  }

  ## 2. Discrepancies + PPP for the observed data
  obs_disc <- disc_fun(MCMC_samples = MCMC_samples,
                       new_data     = observed_data,
                       control      = control)
  if (!all(c("obs", "sim") %in% names(obs_disc))) {
    stop("disc_fun must return a list with components 'obs' and 'sim'.")
  }
  obs_obs <- as.numeric(obs_disc$obs)
  obs_sim <- as.numeric(obs_disc$sim)
  if (length(obs_obs) != length(obs_sim)) {
    stop("'obs' and 'sim' must have the same length.")
  }

  # scalar PPP for the observed data
  PPP_obs <- mean(obs_sim >= obs_obs)

  ## 3. Calibration worlds: PPP in each replicated world
  PPP_rep <- numeric(n_reps)

  for (r in seq_len(n_reps)) {
    theta_row <- MCMC_samples[row_indices[r], , drop = FALSE]

    # 3a. simulate a new dataset y^(r) from posterior predictive of the original model
    new_data <- new_data_fun(theta_row = theta_row,
                             observed_data = observed_data,
                             control = control)

    # 3b. fit model on y^(r) (short chain)
    MCMC_rep <- MCMC_fun(new_data = new_data,
                         control  = control)

    # 3c. compute discrepancies + PPP in this world
    rep_disc <- disc_fun(MCMC_samples = MCMC_rep,
                         new_data     = new_data,
                         control      = control)

    if (!all(c("obs", "sim") %in% names(rep_disc))) {
      stop("disc_fun must return a list with components 'obs' and 'sim'.")
    }
    rep_obs <- as.numeric(rep_disc$obs)
    rep_sim <- as.numeric(rep_disc$sim)
    if (length(rep_obs) != length(rep_sim)) {
      stop("'obs' and 'sim' must have the same length in each calibration world.")
    }

    PPP_rep[r] <- mean(rep_sim >= rep_obs)
  }

  ## 4. CPPP: how extreme PPP_obs is under the calibration distribution
  CPPP <- mean(PPP_rep <= PPP_obs)

  list(
    PPP_obs    = PPP_obs,
    PPP_rep    = PPP_rep,
    CPPP       = CPPP,
    row_indices = row_indices
  )
}




