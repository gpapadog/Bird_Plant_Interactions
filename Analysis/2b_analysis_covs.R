# Analyzing the birds and plants data with a model that uses the covariates.

# The directory where the analysis is performed:
wd_path <- 'Bird_Plant_Interactions/'
# Where the processed data are saved:
data_path <- 'Data/'
# Where you want to save MCMC results:
save_path <- 'Results/'
# Where the functions are available:
source_path <- 'HelperScripts/'


# ------ STEP 0: Some functions. --------- #

setwd(wd_path)
source(paste0(source_path, 'useful_functions.R'))
source(paste0(source_path, 'UpdProbObs_function.R'))
source(paste0(source_path, 'ModelCovariates_function.R'))


# --------------------------------------------------------------- #

# Loading the data:
load(paste0(data_path, 'obs_A.dat'))
load(paste0(data_path, 'obs_W.dat'))
load(paste0(data_path, 'obs_X.dat'))
load(paste0(data_path, 'obs_F.dat'))
load(paste0(data_path, 'obs_OB.dat'))
load(paste0(data_path, 'obs_OP.dat'))
load(paste0(data_path, 'birds_232.dat'))

# Excluding the two variables with the highest missingness.
# Convergence was poor when these covariates were incldued:
obs_W <- obs_W[, - c(4, 12)]

# Restricting to the analysis of the 232 bird species:
wh_keep <- which(rownames(obs_A) %in% birds_232)
obs_A <- obs_A[wh_keep, , ]
obs_F <- obs_F[wh_keep, , ]
obs_OB <- obs_OB[wh_keep, ]
obs_X <- obs_X[wh_keep, ]


nB <- nrow(obs_A)
nP <- nrow(obs_A)

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
  
  mcmc <- ModelCovariates(obs_A = obs_A, focus = obs_F, occur_B = obs_OB, occur_P = obs_OP,
                          obs_X = obs_X, obs_W = obs_W,
                          Nsims = Nsims, burn = burn, thin = thin, bias_cor = bias_cor,
                          mh_n_pis = mh_n_pis, mh_n_pjs = mh_n_pjs, prior_mu0 = prior_mu0,
                          prior_sigmasq0 = prior_sigmasq0, prior_sigmasq = prior_sigmasq,
                          start_values = start_values, sampling = sampling)
  
  all_pred <- abind::abind(pred_L = mcmc$Ls, probL = mcmc$mod_pL1s, pL1s = mcmc$pL1s, along = 4)
  
  res <- list(all_pred = all_pred)
  save(res, file = paste0(save_path, 'alt_res_', cc, '.dat'))
  rm(res)
  
}

