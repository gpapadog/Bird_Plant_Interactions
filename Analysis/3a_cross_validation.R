# Performing cross-validation using the proposed bird-plant interaction model.

# The directory where the analysis is performed:
wd_path <- 'Bird_Plant_Interactions/'
# Where the processed data are saved:
data_path <- 'Data/'
# Where you want to save MCMC results:
result_path <- 'Results/'
# Where the functions are available:
source_path <- 'HelperScripts/'


# ------ STEP 0: Some functions. --------- #

setwd(wd_path)
source(paste0(source_path, 'UpdExtraVar_function.R'))
source(paste0(source_path, 'UpdTraitCoef_function.R'))
source(paste0(source_path, 'UpdLatFac_function.R'))
source(paste0(source_path, 'UpdProbObs_function.R'))
source(paste0(source_path, 'UpdRho_function.R'))
source(paste0(source_path, 'OmegaFromV_function.R'))
source(paste0(source_path, 'useful_functions.R'))
source(paste0(source_path, 'CorrMat_function.R'))
source(paste0(source_path, 'MCMC_function.R'))
source(paste0(source_path, 'PredictInteractions_function.R'))
source(paste0(source_path, 'GetPredLatFac_function.R'))
source(paste0(source_path, 'GetPredWeights_function.R'))


# --------------------------------------------------------------- #

# Loading the data:
load(paste0(data_path, 'Cu_phylo.dat'))
load(paste0(data_path, 'Cv_phylo.dat'))
load(paste0(data_path, 'obs_A.dat'))
load(paste0(data_path, 'obs_W.dat'))
load(paste0(data_path, 'obs_X.dat'))
load(paste0(data_path, 'obs_F.dat'))
load(paste0(data_path, 'obs_OB.dat'))
load(paste0(data_path, 'obs_OP.dat'))
load(paste0(data_path, 'birds_232.dat'))

Cu <- Cu_phylo
Cv <- Cv_phylo

# Restricting to the birds with species information:
wh_keep <- which(rownames(obs_A) %in% birds_232)
obs_A <- obs_A[wh_keep, , ]
obs_F <- obs_F[wh_keep, , ]
obs_OB <- obs_OB[wh_keep, ]
obs_X <- obs_X[wh_keep, ]

# Sample sizes of the two sets of species:
nB <- nrow(Cu)
nP <- nrow(Cv)
nS <- dim(obs_A)[3]
print(c(nB, nP, nS))


# Getting the combined network for the interactions recorded in any study
comb_A <- apply(obs_A, c(1, 2), sum)
comb_A <- (comb_A > 0) * 1



# -------------- STEP 1: Specifications. ------------ #

bias_cor <- TRUE

Nsims <- 600
burn <- 20000
thin <- 20
use_H <- 10
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


# ------------- STEP 2: Setting some interactions to out of sample ------------ #

# We highly recommend running the following code in parallel on 30 machines.

repetitions <- 30


for (rr in 1  : repetitions) {
  
  set.seed(rr)
  
  # Matrix that chooses 100 recorded interactions to be held-out.
  set_out <- matrix(0, nrow = nB, ncol = nP)
  set_out[sample(which(comb_A == 1), 100)] <- 1
  
  # Zero-ing out corresponding entries in A.
  use_A <- obs_A
  use_A[which(set_out == 1)] <- 0  
  
  # Getting the indices that were zero-ed out.
  cv_indices <- matrix(NA, nrow = 100, ncol = 2)
  wh <- which(set_out == 1)
  for (ii in 1 : 100) {
    row_ii <- wh[ii] %% nB
    row_ii <- ifelse(row_ii == 0, nB, row_ii)
    col_ii <- ceiling(wh[ii] / nB)
    cv_indices[ii, ] <- c(row_ii, col_ii)
  }
  
  # Zero-ing out corresponding entries in A.
  use_A <- obs_A
  for (ii in 1 : 100) {
    use_A[cv_indices[ii, 1], cv_indices[ii, 2], ] <- 0
  }
  
  # Running the MCMC with the new recorded interaction matrix:
  set.seed(rr)
  mcmc <- MCMC(obs_A = use_A, focus = obs_F, occur_B = obs_OB, occur_P = obs_OP,
               obs_X = obs_X, obs_W = obs_W, Cu = Cu, Cv = Cv,
               Nsims = Nsims, burn = burn, thin = thin,
               use_H = use_H, bias_cor = bias_cor,
               theta_inf = theta_inf, mh_n_pis = mh_n_pis,
               mh_n_pjs = mh_n_pjs, mh_n_rho = mh_n_rho,
               stick_alpha = stick_alpha, prior_theta = prior_theta,
               prior_tau = prior_tau, prior_rho = prior_rho,
               prior_mu0 = prior_mu0, prior_sigmasq0 = prior_sigmasq0,
               prior_sigmasq = prior_sigmasq, start_values = start_values,
               sampling = sampling)
  
  # Saving the predictions:
  pred <- apply(mcmc$Ls, c(2, 3), mean)
  save(pred, file = paste0(result_path, 'pred_', rr, '.dat'))
  rm(pred)
  
  
  save(cv_indices, file = paste0(result_path, 'cv_indices_', rr, '.dat'))
  rm(cv_indices)
}
  

