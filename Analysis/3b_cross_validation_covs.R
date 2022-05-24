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
load(paste0(data_path, 'obs_A.dat'))
load(paste0(data_path, 'obs_W.dat'))
load(paste0(data_path, 'obs_X.dat'))
load(paste0(data_path, 'obs_F.dat'))
load(paste0(data_path, 'obs_OB.dat'))
load(paste0(data_path, 'obs_OP.dat'))
load(paste0(data_path, 'birds_232.dat'))

Cu <- Cu_phylo
Cv <- Cv_phylo

# Excluding the two variables with the highest missingness
obs_W <- obs_W[, - c(4, 12)]

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

# Getting the combined network:
comb_A <- apply(obs_A, c(1, 2), sum)
comb_A <- (comb_A > 0) * 1


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
  mcmc <- ModelCovariates(obs_A = use_A, focus = obs_F, occur_B = obs_OB, occur_P = obs_OP,
                          obs_X = obs_X, obs_W = obs_W,
                          Nsims = Nsims, burn = burn, thin = thin, bias_cor = bias_cor,
                          mh_n_pis = mh_n_pis, mh_n_pjs = mh_n_pjs, prior_mu0 = prior_mu0,
                          prior_sigmasq0 = prior_sigmasq0, prior_sigmasq = prior_sigmasq,
                          start_values = start_values, sampling = sampling)
  
  alt_pred <- apply(mcmc$Ls, c(2, 3), mean)
  
  save(alt_pred, file = paste0(result_path, 'alt_pred_', rr, '.dat'))
  rm(alt_pred)
  
}

