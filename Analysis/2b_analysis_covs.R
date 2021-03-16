# Analyzing the birds and plants data with a model that uses the covariates.

# Where the processed data are saved:
data_path <- '~/Github/Birds_and_plants/Application/Data/Aves_analysis/'
# Where you want to save MCMC results:
save_path <- '~/Github/Birds_and_plants/Application/Results/'
# Where the functions are available:
source_path <- '~/Github/BiasedNetwork/R/'


# ------ STEP 0: Some functions. --------- #

source(paste0(source_path, 'useful_functions.R'))
source(paste0(source_path, 'UpdProbObs_function.R'))
source('~/Github/Birds_and_plants/Simulations/functions/ModelCovariates_function.R')
source('~/Github/Birds_and_plants/Simulations/functions/GenDataCovariates_function.R')
source('~/Github/Birds_and_plants/Simulations/functions/PredPower_function.R')
source('~/Github/Birds_and_plants/Simulations/functions/AllPredPower_function.R')

# --------------------------------------------------------------- #

# Loading the data:
load(paste0(data_path, 'Cu.dat'))
load(paste0(data_path, 'Cv.dat'))
load(paste0(data_path, 'obs_A.dat'))
load(paste0(data_path, 'obs_n.dat'))
load(paste0(data_path, 'obs_W.dat'))
load(paste0(data_path, 'obs_X.dat'))

# Excluding the two variables with the highest missingness.
# Convergence was poor when these covariates were incldued:
obs_W <- obs_W[, - c(4, 12)]

nB <- nrow(Cu)
nP <- nrow(Cv)

# -------------- STEP 1: Specifications. ------------ #

bias_cor <- TRUE

Nsims <- 1000
burn <- 40000
thin <- 40
mh_n_pis <- 70  # Parameter for proposal in Metropolis-Hastings for pi update.
mh_n_pjs <- 70

# Prior distributions:
prior_mu0 <- 0
prior_sigmasq0 <- 10
prior_sigmasq <- c(1, 1)

start_values <- NULL
sampling <- NULL


# ------------ STEP 2: MCMC -------------- #

# We recommend running the three chains in parallel:

for (cc in 1 : 3) {
  
  mcmc <- ModelCovariates(obs_A = obs_A, obs_n = obs_n, obs_X = obs_X, obs_W = obs_W,
                          Nsims = Nsims, burn = burn, thin = thin, bias_cor = bias_cor,
                          mh_n_pis = mh_n_pis, mh_n_pjs = mh_n_pjs, prior_mu0 = prior_mu0,
                          prior_sigmasq0 = prior_sigmasq0, prior_sigmasq = prior_sigmasq,
                          start_values = start_values, sampling = sampling)
  
  all_pred <- abind::abind(pred_L = mcmc$Ls, probL = mcmc$mod_pL1s, pL1s = mcmc$pL1s, along = 4)
  
  res <- list(all_pred = all_pred)
  save(res, file = paste0(save_path, 'alt_res_', cc, '.dat'))
  rm(res)
  
}

