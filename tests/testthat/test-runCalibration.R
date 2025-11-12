# tests/testthat/test-runCalibration.R

test_that("runCalibration returns numeric PPP vector of correct length", {
  skip_if_not("runCalibration" %in% ls("package:cpppPrototype"),
              "runCalibration() not exported yet")

  set.seed(1)

  # Fake "long-run" MCMC samples
  MCMC_samples <- matrix(0, nrow = 10, ncol = 1)

  # Toy new_data_fun: ignore MCMC_samples, just simulate data
  new_data_fun <- function(MCMC_samples, ...) {
    rnorm(20)
  }

  # Toy MCMC_fun: given data, return a "chain" whose mean depends on the data
  MCMC_fun <- function(new_data, control, ...) {
    theta <- rnorm(100, mean(new_data), sd = 1 / sqrt(length(new_data)))
    cbind(theta = theta)
  }

  # Toy disc_fun: discrepancy is |mean(y)|; PPP via comparing |theta| to |mean(y)|
  disc_fun <- function(mcmc_samples, ...) {
    d_obs <- abs(mean(mcmc_samples[, "theta"]))
    d_rep <- abs(mcmc_samples[, "theta"])  # pretend these are replicated discrepancies
    list(
      rep = d_rep,
      obs = d_obs
    )
  }

  num_reps <- 5

  ppp <- runCalibration(
    MCMC_samples = MCMC_samples,
    MCMC_fun     = MCMC_fun,
    new_data_fun = new_data_fun,
    disc_fun     = disc_fun,
    num_reps     = num_reps
  )

  expect_type(ppp, "double")
  expect_length(ppp, num_reps)
  expect_true(all(is.finite(ppp)))
  expect_true(all(ppp >= 0 & ppp <= 1))
})

test_that("runCalibration checks num_reps and disc_fun output", {
  skip_if_not("runCalibration" %in% ls("package:cpppPrototype"),
              "runCalibration() not exported yet")

  MCMC_samples <- matrix(0, nrow = 5, ncol = 1)

  new_data_fun <- function(MCMC_samples, ...) rnorm(10)
  MCMC_fun     <- function(new_data, control, ...) matrix(new_data, ncol = 1)

  bad_disc_fun <- function(mcmc_samples, ...) {
    0.5  # numeric scalar, no $rep / $obs
  }

  # num_reps must be positive integer
  expect_error(
    runCalibration(MCMC_samples, MCMC_fun, new_data_fun, bad_disc_fun, num_reps = 0),
    "positive integer"
  )

  # disc_fun must return list with $rep and $obs
  expect_error(
    runCalibration(MCMC_samples, MCMC_fun, new_data_fun, bad_disc_fun, num_reps = 1),
    "list with components `rep` and `obs`"
  )
})
