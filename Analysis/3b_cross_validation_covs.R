# Performing cross-validation of bird-plant interaction models based on the
# model that depends directly on covariates.

# The directory where the analysis is performed:
wd_path <- 'Bird_Plant_Interactions/'
# Where the processed data are saved:
data_path <- 'Data/'
# Where the functions are available:
source_path <- 'HelperScripts/'
# Where the results should be saved:
result_path <- 'Results/'

# ------ STEP 0: Some functions. --------- #

setwd(wd_path)
source(paste0(source_path, 'useful_functions.R'))
source(paste0(source_path, 'UpdProbObs_function.R'))
source(paste0(source_path, 'ModelCovariates_function.R'))

# --------------------------------------------------------------- #

# Loading the data:
load(paste0(data_path, 'Cu.dat'))
load(paste0(data_path, 'Cv.dat'))
load(paste0(data_path, 'obs_A.dat'))
load(paste0(data_path, 'obs_n.dat'))
load(paste0(data_path, 'obs_W.dat'))
load(paste0(data_path, 'obs_X.dat'))

# Excluding the two variables with the highest missingness
obs_W <- obs_W[, - c(4, 12)]

# Sample sizes of the two sets of species:
nB <- nrow(Cu)
nP <- nrow(Cv)

# -------------- STEP 1: Specifications. ------------ #

bias_cor <- TRUE

Nsims <- 600
burn <- 20000
thin <- 20
mh_n_pis <- 100  # Parameter for proposal in Metropolis-Hastings for pi update.
mh_n_pjs <- 100

# Prior distributions:
prior_mu0 <- 0
prior_sigmasq0 <- 10
prior_sigmasq <- c(1, 1)

start_values <- NULL
sampling <- NULL


# ------------- STEP 2: Setting some interactions to out of sample ------------ #

# We highly recommend running the following code in parallel on 20 machines.

repetitions <- 20


for (rr in 1  : repetitions) {
  
  set.seed(rr)
  
  # Matrix that chooses 100 recorded interactions:
  set_out <- matrix(0, nrow = nB, ncol = nP)
  set_out[sample(which(obs_A == 1), 100)] <- 1  
  
  # Zero-ing out corresponding entries in A.
  use_A <- obs_A
  use_A[which(set_out == 1)] <- 0  
  
  
  # Running the MCMC with the new recorded interaction matrix:
  set.seed(rr)
  mcmc <- ModelCovariates(obs_A = use_A, obs_n = obs_n, obs_X = obs_X, obs_W = obs_W,
                          Nsims = Nsims, burn = burn, thin = thin, bias_cor = bias_cor,
                          mh_n_pis = mh_n_pis, mh_n_pjs = mh_n_pjs, prior_mu0 = prior_mu0,
                          prior_sigmasq0 = prior_sigmasq0, prior_sigmasq = prior_sigmasq,
                          start_values = start_values, sampling = sampling)
  
  alt_pred <- apply(mcmc$Ls, c(2, 3), mean)
  
  save(alt_pred, file = paste0(result_path, 'alt_pred_', rr, '.dat'))
  rm(alt_pred)
  
}

