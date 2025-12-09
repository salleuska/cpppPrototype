#' Run calibration using a NIMBLE model
#'
#' @param model either an uncompiled or compiled nimbleModel with observed data set.
#' @param dataNames Optional character vector of data node names. If NULL,
#'   nodes flagged as data in the model are used.
#' @param paramNames Character vector of parameter node names to monitor.
#' @param disc_fun Function `function(data, theta_row, control)` returning a
#'   scalar or vector discrepancy for one posterior draw.
#' @param new_data_fun Function `function(theta_row, observed_data, control)` that  simulates one replicated dataset from the posterior predictive. SP: We assume that new data is sampled from the posterior predictive of the model. In principle we may want to consider sampling from the prior predictive.
#' @param n_reps Number of calibration replications.
#' @param MCMCcontrolMain List with `niter`, `nburnin`, `thin` for main chain.
#' @param MCMCcontrolRep List with `niter`, `nburnin`, `thin` for calibration chains.
#' @param mcmcConfFun Optional function `function(model)` returning an MCMC configuration.
#' @param row_selector Optional row selector (see runCalibration()).
#' @param control List of additional options passed to discrepancy / helpers.
#' @param ... Not used currently.
#' @export

runCalibrationNIMBLE <- function(
    model,
    dataNames    = NULL,
    paramNames,
    disc_fun,
    new_data_fun,
    n_reps       = 100,
    MCMCcontrolMain = list(niter = 5000, nburnin = 1000, thin = 1),
    MCMCcontrolRep  = list(niter = 500,  nburnin = 0,    thin = 1),
    mcmcConfFun = NULL,
    row_selector = NULL,
    control = list(),
    ...
) {
  ## 0. Data names and checks
  if (is.null(dataNames)) {
    dataNames <- model$getNodeNames(dataOnly = TRUE)
  }
  ## expand to nodes
  dataNodes <- model$expandNodeNames(dataNames)
  # ensure dataNames correspond to stochastic nodes
  testDataNames <- all(dataNodes %in%
                         model$getNodeNames(stochOnly = TRUE))
  if (!testDataNames) {
    stop("All dataNames must be stochastic nodes in the model.")
  }

  ## check if the model is compiled model
  if (inherits(model, "CmodelBaseClass")) {
    ## Model is already compiled
    cmodel <- model
  } else if (inherits(model, "RmodelBaseClass")) {
    ## Model is not compiled yet
    cmodel <- compileNimble(model)
  } else {
    stop("Argument 'model' must be a nimbleModel or a compiled nimble model.")
  }

  ## 1. Configure and compile MCMC for main chain
  if (is.null(mcmcConfFun)) {
    mcmcConfFun <- function(model) {
      configureMCMC(model, monitors = paramNames, print = FALSE)
    }
  }
  mcmcConf       <- mcmcConfFun(model)
  mcmcUncompiled <- buildMCMC(mcmcConf)
  cmcmc          <- compileNimble(mcmcUncompiled, project = model, resetFunctions = TRUE)

  ## 2. Run main chain on observed data
  main_out <- runMCMC(
    cmcmc,
    niter   = MCMCcontrolMain$niter,
    nburnin = MCMCcontrolMain$nburnin,
    thin    = MCMCcontrolMain$thin
  )
  MCMC_samples <- as.matrix(main_out)

  ## Extract observed data from the model
  ## SP: need to think better if data is a vector/matrix/array - something else?
  observed_data <- cmodel[[dataNames]]


  ## 4. Build MCMC_fun for replicated datasets
  ## SP: do we want a make_MCMCfun?
  MCMC_fun <- function(new_data,
                       control) {
    # assign new data into cmodel
    cmodel[[dataNames]] <- new_data
    # run short chain
    rep_out <- runMCMC(
      cmcmc,
      niter   = MCMCcontrolRep$niter,
      nburnin = MCMCcontrolRep$nburnin,
      thin    = MCMCcontrolRep$thin
    )
    as.matrix(rep_out)
  }

  ## 5. Call generic engine
  runCalibration(
    MCMC_samples  = MCMC_samples,
    observed_data = observed_data,
    MCMC_fun      = MCMC_fun,
    new_data_fun  = new_data_fun,
    disc_fun      = disc_fun,
    n_reps        = n_reps,
    row_selector  = row_selector,
    control       = control
  )

}

