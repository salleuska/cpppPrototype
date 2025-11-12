############################################################
## Newcomb data example: NIMBLE + offline discrepancy + CPPP
############################################################

## 1) Packages ------------------------------------------------
library(nimble)
library(cpppPrototype)   # or whatever your package is called


## 2) Data: Newcomb light-speed measurements ------------------
## Assuming light.txt is in the working directory
newcomb_data <- list(
  y = read.table("inst/examples/light.txt")$V1
)
newcomb_constants <- list(
  n = length(newcomb_data$y)
)


## 3) NIMBLE model -------------------------------------------
newcomb_code <- nimbleCode({
  for (i in 1:n) {
    y[i] ~ dnorm(mu, sd = sigma)
  }
  mu ~ dflat()
  log(sigma) ~ dflat()
})

newcomb_inits <- list(mu = 0, log_sigma = 2)

newcomb_model <- nimbleModel(
  code      = newcomb_code,
  data      = newcomb_data,
  inits     = newcomb_inits,
  constants = newcomb_constants
)

dataNames  <- "y"
paramNames <- c("mu", "log_sigma")


## 4) Offline discrepancy pieces ------------------------------
## Discrepancy: data-only, min(y)
min_disc <- function(data, theta_row, control) {
  min(data$y)
}

## Inner posterior predictive simulator used by make_offline_disc_fun
## This uses the *uncompiled* model to simulate y* given theta_row.
newcomb_pp_sim_inner <- local({
  model      <- newcomb_model
  paramNames <- paramNames  # captured from outer scope

  function(theta_row, ...) {
    theta_vec <- as.numeric(theta_row)
    names(theta_vec) <- paramNames

    ## write parameters into the model
    for (nm in paramNames) {
      model[[nm]] <<- theta_vec[[nm]]
    }

    ## simulate downstream including data
    model$simulate(nodes = "y", includeData = TRUE)

    list(y = as.numeric(model[["y"]]))
  }
})

## Build disc_fun via the package helper
disc_control <- list(
  new_data_fun = newcomb_pp_sim_inner,
  discrepancy  = min_disc
)
disc_fun <- make_offline_disc_fun(disc_control)


## 5) Run calibration via the package wrapper -----------------
## Note: this uses runCalibrationNIMBLE from your package, which internally
## handles the main chain and short chains. The outer new_data_fun used for
## calibration worlds is whatever implementation you currently have in the
## package (you can later upgrade it to use a posterior predictive simulator).

set.seed(1)
res_newcomb <- runCalibrationNIMBLE(
  model      = newcomb_model,
  dataNames  = dataNames,
  paramNames = paramNames,
  disc_fun   = disc_fun,
  n_reps     = 50,  # smallish number for quick runs
  MCMCcontrolMain = list(niter = 5000, nburnin = 1000, thin = 1),
  MCMCcontrolRep  = list(niter = 1000, nburnin = 200,  thin = 1)
)

print(res_newcomb$PPP_obs)
print(res_newcomb$CPPP)
