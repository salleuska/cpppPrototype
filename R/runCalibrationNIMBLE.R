#' Run calibration using a NIMBLE model
#'
#' @param model either an uncompiled or compiled nimbleModel with observed data set.
#' @param dataNames Optional character vector of data node names. If NULL,
#'   nodes flagged as data in the model are used.
#' @param paramNames Character vector of parameter node names to monitor. (SP: if NULL some default?)
#' @param discFun Function `function(data, theta_row, control)` returning a
#'   scalar or vector discrepancy for one posterior draw.
#' @param simulateNewDataFun Function `function(theta_row, observedData, control)` that  simulates one replicated dataset from the posterior predictive. SP: We assume that new data is sampled from the posterior predictive of the model. In principle we may want to consider sampling from the prior predictive.
#' @param nReps Number of calibration replications.
#' @param MCMCcontrolMain List with `niter`, `nburnin`, `thin` for main chain.
#' @param MCMCcontrolRep List with `niter`, `nburnin`, `thin` for calibration chains.
#' @param mcmcConfFun Optional function `function(model)` returning an MCMC configuration.
#' @param drawIndices Optional row selector (see runCalibration()).
#' @param control List of additional options passed to discrepancy / helpers.
#' @param ... Not used currently.
#' @export

runCalibrationNIMBLE <- function(
    model,
    dataNames    = NULL,
    paramNames,
    discFun,
    simulateNewDataFun,
    nReps       = 100,
    MCMCcontrolMain = list(niter = 5000, nburnin = 1000, thin = 1),
    MCMCcontrolRep  = list(niter = 500,  nburnin = 0,    thin = 1),
    mcmcConfFun = NULL,
    drawIndices = NULL,
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
  obsMCMC <- runMCMC(
    cmcmc,
    niter   = MCMCcontrolMain$niter,
    nburnin = MCMCcontrolMain$nburnin,
    thin    = MCMCcontrolMain$thin
  )
  MCMCSamples <- as.matrix(obsMCMC)

  ## Extract observed data from the model
  ## SP: need to think better if data is a vector/matrix/array - something else?
  observedData <- cmodel[[dataNames]]


  ## 4. Build MCMCFun for replicated datasets
  ## SP: do we want a function makeMCMCfun?
  MCMCFun <- function(newData,
                       control) {
    # assign new data into cmodel
    cmodel[[dataNames]] <- newData
    # run short chain
    repMCMC <- runMCMC(
      cmcmc,
      niter   = MCMCcontrolRep$niter,
      nburnin = MCMCcontrolRep$nburnin,
      thin    = MCMCcontrolRep$thin
    )
    as.matrix(repMCMC)
  }

  ## 5. Call generic engine
  runCalibration(
    MCMCSamples  = MCMCSamples,
    observedData = observedData,
    MCMCFun      = MCMCFun,
    simulateNewDataFun  = simulateNewDataFun,
    discFun      = discFun,
    nReps        = nReps,
    drawIndices  = drawIndices,
    control       = control
  )

}

