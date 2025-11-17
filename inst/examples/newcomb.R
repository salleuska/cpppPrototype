############################################################
## Newcomb data example: NIMBLE + offline discrepancy + CPPP
############################################################

## 1) Packages
library(nimble)
library(cpppPrototype)


## 2) Data: Newcomb light-speed measurements
## Assuming light.txt is in the working directory
newcombData <- list(
  y = read.table("inst/examples/light.txt")$V1
)
constants <- list(
  n = length(newcombData$y)
)


## 3) NIMBLE model
newcombCode <- nimbleCode({
  for (i in 1:n) {
    y[i] ~ dnorm(mu, sd = sigma)
  }
  mu ~ dflat()
  log(sigma) ~ dflat()
})


inits <- list(mu = 0, log_sigma = 2)

newcomb_model <- nimbleModel(
  code      = newcombCode,
  data      = newcombData,
  inits     = inits,
  constants = constants
)

dataNames  <- "y"
paramNames <- c("mu", "log_sigma")

## 4) Offline discrepancy
## Discrepancy: data-only, min(y)
min_disc <- function(data, theta_row, control) {
  min(data$y)
}

## Uses newcomb_model and paramNames directly from the enclosing environment
newcomb_newData <- function(theta_row, ...) {
  theta_vec <- as.numeric(theta_row)
  names(theta_vec) <- paramNames  # paramNames defined earlier

  ## write parameters into the model
  for (nm in paramNames) {
    newcomb_model[[nm]] <<- theta_vec[[nm]]
  }

  ## simulate downstream including data
  newcomb_model$simulate(nodes = "y", includeData = TRUE)

  list(y = as.numeric(newcomb_model[["y"]]))
}

## Build disc_fun via the package helper
disc_control <- list(
  new_data_fun = newcomb_newData,
  discrepancy  = min_disc
)
disc_fun <- make_offline_disc_fun(disc_control)

# #####
# ## Test disc_fun
# ## fake θ draws, just to test disc_fun mechanics
# MCMC_samples_test <- matrix(
#   c(0, 2,   # mu = 0, log_sigma = 2
#     5, 1),  # mu = 5, log_sigma = 1
#   ncol = 2,
#   byrow = TRUE
# )
# colnames(MCMC_samples_test) <- paramNames  # c("mu", "log_sigma")
# new_data_test <- list(y = newcombData$y)
# disc_fun(MCMC_samples_test, new_data_test)
##############

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
