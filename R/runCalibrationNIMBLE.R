
runCalibrationNIMBLE <- function(
    model,
    dataNames    = NULL,
    paramNames,
    mcmcConfFun  = NULL,
    MCMCcontrolMain = list(niter = 5000, thin = 1, nburnin = 1000),
    MCMCcontrol     = list(niter = 500,  thin = 1, nburnin = 0),
    nCalibrationReplicates = 100,
    new_data_fun = NULL,
    discrepancy  = NULL,
    disc_control = list(),
    ...
) {
  ## 0. Default/validate dataNames, mcmcConfFun, nCalibrationReplicates, etc.
  # (reuse your existing checks almost verbatim)


  ## 1. Build Nimble MCMC config and compiled objects
  mcmcConf <- if (is.null(mcmcConfFun)) {
    configureMCMC(model, monitors = paramNames, print = FALSE)
  } else {
    mcmcConfFun(model)
  }

  # Optional: if you keep derived discrepancy machinery
  # addDiscrepancies(mcmcConf, discrepancyFunctions, discrepancyFunctionsArgs)

  mcmcUncompiled <- buildMCMC(mcmcConf)
  cmcmc          <- compileNimble(mcmcUncompiled, project = model, resetFunctions = TRUE)
  cmodel         <- model  # if already compiled; otherwise, result of compileNimble(model)

  ## 2. Run long chain on observed data
  originalOutput <- runMCMC(cmcmc,
                            niter   = MCMCcontrolMain$niter,
                            nburnin = MCMCcontrolMain$nburnin,
                            thin    = MCMCcontrolMain$thin)
  MCMC_samples   <- originalOutput$samples

  ## 3. Build new_data_fun (if not supplied) using setAndSimNodes
  if (is.null(new_data_fun)) {
    simNodes <- unique(c(
      model$expandNodeNames(dataNames),
      model$getDependencies(paramNames, includeData = FALSE, self = FALSE)
    ))
    setAndSimPP   <- setAndSimNodes(model = model, nodes = paramNames, simNodes = simNodes)
    cSetAndSimPP  <- compileNimble(setAndSimPP, project = model)

    new_data_fun <- make_nimble_new_data_fun(cSetAndSimPP, dataNames)
  }

  ## 4. Build MCMC_fun
  MCMC_fun <- make_MCMCfun(
    calib_niter          = MCMCcontrol$niter,
    simulated_data_nodes = dataNames,
    cmodel               = cmodel,
    cmcmc                = cmcmc,
    MCMCcontrol          = MCMCcontrol
  )

  ## 5. Build disc_fun: either Nimble-derived via calculatePPP or offline
  if (is.null(discrepancy)) {
    # option A: Nimble-derived approach
    disc_fun <- make_nimble_disc_fun(...)
  } else {
    # option B: offline discrepancy; build disc_fun around user discrepancy
    disc_fun <- make_offline_disc_fun(discrepancy, disc_control)
  }

  ## 7. Call the generic engine
  results <- runCalibration(
    MCMC_samples = MCMC_samples,
    MCMC_fun     = MCMC_fun,
    new_data_fun = new_data_fun,
    disc_fun     = disc_fun,
    control      = control,
    num_reps    = nCalibrationReplicates,
    ...
  )

}
