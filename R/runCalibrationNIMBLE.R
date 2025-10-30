runCalibrationNIMBLE <- function(nimble_specific_args) {
  # Use all the nimble-specific steps to create the generic inputs
  #  that are a set of functions for runCalibration
  # build and compile the model and mcmc pieces
  # Set up "MCMC_fun"
  MCMCfun <- make_MCMCfun(calib_niter, simulated_data_nodes, cmodel, cmcmc)
  # Set up disc_fun
  # Set up the data_sim_fun
  disc_fun <- make_col_disc_fun("some_disc_output")
  results <- runCalibration(MCMCfun, disc_fun)
  results

}
