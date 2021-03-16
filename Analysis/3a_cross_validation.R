# Performing cross-validation of bird-plant interaction models.
# Created: October 12, 2020.

rm(list = ls())
dev.off()

# ------ STEP 0: Some functions. --------- #

source('~/Github/Birds_and_plants/functions/UpdExtraVar_function.R')
source('~/Github/Birds_and_plants/functions/UpdTraitCoef_function.R')
source('~/Github/Birds_and_plants/functions/UpdLatFac_function.R')
source('~/Github/Birds_and_plants/functions/UpdProbObs_function.R')
source('~/Github/Birds_and_plants/functions/UpdRho_function.R')
source('~/Github/Birds_and_plants/functions/OmegaFromV_function.R')
source('~/Github/Birds_and_plants/functions/useful_functions.R')
source('~/Github/Birds_and_plants/functions/CorrMat_function.R')
source('~/Github/Birds_and_plants/functions/MCMC_function.R')
source('~/Github/Birds_and_plants/functions/PredictInteractions_function.R')
source('~/Github/Birds_and_plants/functions/GetPredLatFac_function.R')
source('~/Github/Birds_and_plants/functions/GetPredWeights_function.R')
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

in_sample_mcmc <- FALSE
bias_cor <- TRUE


# ------------- STEP 2: Setting some interactions to out of sample ------------ #

set.seed(1234)
set_out <- matrix(0, nrow = nB, ncol = nP)
set_out[sample(which(obs_A == 1), 100)] <- 1

# Zero-ing out corresponding entries in A.
use_A <- obs_A
use_A[which(set_out == 1)] <- 0


# STEP 3: MCMC specifications.

Nsims <- 1000
burn <- 3
thin <- 2
use_H <- 5
theta_inf <- 0.01
mh_n_pis <- 100  # Parameter for proposal in Metropolis-Hastings for pi update.
mh_n_pjs <- 100
mh_n_rho <- 100

# Prior distributions:
stick_alpha <- 5
prior_theta <- c(1, 1)
prior_tau <- c(5, 5)
prior_rho <- c(5, 5)  # I do not update this for now.
prior_mu0 <- 0
prior_sigmasq0 <- 10
prior_sigmasq <- c(1, 1)


start_values <- NULL
sampling <- NULL


# STEP 4: MCMC.

mcmc <- MCMC(obs_A = use_A, obs_n = obs_n, obs_X = obs_X, obs_W = obs_W,
             Cu = Cu, Cv = Cv, Nsims = Nsims, burn = burn, thin = thin,
             use_H = use_H, bias_cor = bias_cor,
             theta_inf = theta_inf, mh_n_pis = mh_n_pis,
             mh_n_pjs = mh_n_pjs, mh_n_rho = mh_n_rho,
             stick_alpha = stick_alpha, prior_theta = prior_theta,
             prior_tau = prior_tau, prior_rho = prior_rho,
             prior_mu0 = prior_mu0, prior_sigmasq0 = prior_sigmasq0,
             prior_sigmasq = prior_sigmasq, start_values = start_values,
             sampling = sampling)


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


