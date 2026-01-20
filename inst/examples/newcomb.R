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
  min(data)
}

## asymmetry discrepancy using R
asymm_disc <- function(data, theta_row, control) {

  mu <- theta_row["mu"]
  dataSorted <- sort(data)

  abs(dataSorted[61] - mu) - abs(dataSorted[6] - mu)
}


## function that generates new data

newcomb_newData <- function(theta_row, paramNames, dataNames) {
  theta_vec <- as.numeric(theta_row)

  ## write parameters into the model
  for (nm in paramNames) {
    newcomb_model[[nm]] <<- theta_vec[[nm]]
  }

  ## simulate downstream including data
  newcomb_model$simulate(nodes = dataNames, includeData = TRUE)

  newdata <- newcomb_model[[dataNames]]
  newdata
}

## Build disc_fun via the package helper
# disc_control <- list(
#   new_data_fun = newcomb_newData,
#   discrepancy  = min_disc
# )

disc_control <- list(
  new_data_fun = newcomb_newData,
  discrepancy  = asymm_disc
)

disc_fun <- make_offline_disc_fun(disc_control)



set.seed(1)
res_newcomb <- runCalibrationNIMBLE(
  model      = newcomb_model,
  dataNames  = dataNames,
  paramNames = paramNames,
  disc_fun   = disc_fun,
  new_data_fun = newcomb_newData,
  n_reps       = 100,  # smallish number for quick runs
  MCMCcontrolMain = list(niter = 10000, nburnin = 5000, thin = 1),
  MCMCcontrolRep  = list(niter = 500, nburnin = 0,  thin = 1)
)


print(res_newcomb$cppp)
print(res_newcomb$ppp)

# obsDisc <- res_newcomb$discrepancies$obs
# plot(obsDisc$obs, obsDisc$sim)
# abline()
############################




######
samples <- nimbleMCMC(newcomb_model, niter = 5000, nburn = 1000)
newcomb_model$y <- newcombData$y

# # Test disc_fun
# new_data_test <- newcombData$y
# res <- disc_fun(samples, new_data_test)
# str(res)

## 2. Discrepancies + PPP for the observed data
obs_disc <- disc_fun(MCMC_samples = samples,
                     new_data     = newcomb_model$y)
plot(obsDisc$obs, obsDisc$sim)
mean(obsDisc$sim > obsDisc$obs)
abline()
####
## line by line check
dObs <- apply(samples, 1, function(th) asymm_disc(newcomb_model$y, th))
hist(dObs)
sim_data <- list()
disc_rep <- list()
for(i in seq_along(NROW(samples))) {
  sim_data[[i]]  <- newcomb_newData(samples[i, ], paramNames, dataNames)
  disc_rep[[i]] <-   asymm_disc(samples[i, ], i)
}
# dSim <- apply(samples, 1, function(th) {
#   y <- newcomb_newData(th)
#   asymm_disc(y, th)
# })
hist(dSim)
plot(dObs, dSim)
abline()
##############
