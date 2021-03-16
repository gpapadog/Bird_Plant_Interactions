# Performing cross-validation of bird-plant interaction models based on the
# model that depends directly on covariates.
#
# Created: October 12, 2020.

rm(list = ls())
dev.off()

# ------ STEP 0: Some functions. --------- #

source('~/Github/Birds_and_plants/functions/useful_functions.R')
source('~/Github/Birds_and_plants/functions/UpdProbObs_function.R')
source('~/Github/Birds_and_plants/Simulations/functions/ModelCovariates_function.R')
source('~/Github/Birds_and_plants/Simulations/functions/GenDataCovariates_function.R')
source('~/Github/Birds_and_plants/Simulations/functions/PredPower_function.R')
source('~/Github/Birds_and_plants/Simulations/functions/AllPredPower_function.R')

# --------------------------------------------------------------- #

load('~/Github/Birds_and_plants/Application/Data/Aves_analysis/Cu.dat')
load('~/Github/Birds_and_plants/Application/Data/Aves_analysis/Cv.dat')
load('~/Github/Birds_and_plants/Application/Data/Aves_analysis/obs_A.dat')
load('~/Github/Birds_and_plants/Application/Data/Aves_analysis/obs_n.dat')
load('~/Github/Birds_and_plants/Application/Data/Aves_analysis/obs_W.dat')
load('~/Github/Birds_and_plants/Application/Data/Aves_analysis/obs_X.dat')

nB <- nrow(Cu)
nP <- nrow(Cv)

# -------------- STEP 1: Specifications. ------------ #

bias_cor <- TRUE


# ------------- STEP 2: Setting some interactions to out of sample ------------ #

set.seed(1234)
set_out <- matrix(0, nrow = nB, ncol = nP)
set_out[sample(which(obs_A == 1), 100)] <- 1

# Zero-ing out corresponding entries in A.
use_A <- obs_A
use_A[which(set_out == 1)] <- 0


# ------------- STEP 3: MCMC specifications. ------------ #

Nsims <- 1000
burn <- 300
thin <- 2
mh_n_pis <- 100  # Parameter for proposal in Metropolis-Hastings for pi update.
mh_n_pjs <- 100

# Prior distributions:
prior_mu0 <- 0
prior_sigmasq0 <- 10
prior_sigmasq <- c(1, 1)

start_values <- NULL
sampling <- NULL


# ------------- STEP 4: MCMC. ------------ #

mcmc <- ModelCovariates(obs_A = use_A, obs_n = obs_n, obs_X = obs_X, obs_W = obs_W,
                        Nsims = Nsims, burn = burn, thin = thin, bias_cor = bias_cor,
                        mh_n_pis = mh_n_pis, mh_n_pjs = mh_n_pjs, prior_mu0 = prior_mu0,
                        prior_sigmasq0 = prior_sigmasq0, prior_sigmasq = prior_sigmasq,
                        start_values = start_values, sampling = sampling)



pred <- matrix(NA, nrow = Nsims, ncol = 100)
colnames(pred) <- 1 : 100
wh <- which(set_out == 1)

for (ii in 1 : 100) {
  row_ii <- wh[ii] %% nB
  row_ii <- ifelse(row_ii == 0, nB, row_ii)
  col_ii <- ceiling(wh[ii] / nB)
  if (set_out[row_ii, col_ii] == 0) {
    warning('something is wrong here')
  }
  pred[, ii] <- mcmc$Ls[, row_ii, col_ii]
  colnames(pred)[ii] <- paste0(row_ii, ', ', col_ii)
}


