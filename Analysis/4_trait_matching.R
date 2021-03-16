# Performing cross-validation of bird-plant interaction models.

# The directory where the analysis is performed:
wd_path <- 'Bird_Plant_Interactions/'
# Where the processed data are saved:
data_path <- '~/Github/Birds_and_plants/Application/Data/Aves_analysis/'
# Where the MCMC results are saved and the trait matching will be saved:
save_path <- '~/Github/Birds_and_plants/Application/Results/'
# Where the functions are available:
source_path <- '~/Github/BiasedNetwork/R/'


# ------ STEP 0: Some functions. --------- #

setwd(wd_path)
source(paste0(source_path, 'TraitMatching2_function.R'))
source(paste0(source_path, 'useful_functions.R'))

# --------------------------------------------------------------- #

# Loading the data:
load(paste0(data_path, 'obs_W.dat'))
load(paste0(data_path, 'obs_X.dat'))

# Getting the sample sizes:
nB <- nrow(obs_X)
nP <- nrow(obs_W)



# --------------- STEP 1: Getting the results together ----------------- #

# MCMC chains saved:
nchains <- 3

# Putting together the predictions from the chains:
all_pred <- NULL
for (ii in 1 : nchains) {
  load(paste0(save_path, 'res_', ii, '.dat'))
  all_pred[[ii]] <- res$all_pred
}

# Number of posterior samples by chain:
Nsims <- dim(all_pred[[1]])[1]

# Using the linear predictor of the interaction model:
mod_pL1s <- array(NA, dim = c(nchains * Nsims, nB, nP))
for (cc in 1 : 3) {
  wh_entries <- Nsims * (cc - 1) + 1 : Nsims
  mod_pL1s[wh_entries, , ] <- all_pred[[cc]][, , , 2]
}



# --------------- STEP 2: Performing trait matching ----------------- #

# Dealing with extreme values for which logit(x) is infinite.
mod_pL1s[mod_pL1s > 1 - 10^{-10}] <- 1 - 10^{-10}

trait_match <- TraitMatching2(B = 500, mod_pL1s = use_mod_pL1s,
                              Xs = NULL, Ws = NULL,  # Imputed values not used.
                              obs_X = obs_X, obs_W = obs_W, obs_only = TRUE)

rsq_resampling_X <- trait_match$rsq_resampling_X
rsq_resampling_W <- trait_match$rsq_resampling_W
rsq_obs_X <- trait_match$rsq_obs_X
rsq_obs_W <- trait_match$rsq_obs_W

save(rsq_obs_X, file = paste0(save_path,  'rsq_obs_X.dat'))
save(rsq_obs_W, file = paste0(save_path, 'rsq_obs_W.dat'))
save(rsq_resampling_X, file = paste0(save_path, 'rsq_resampling_X.dat'))
save(rsq_resampling_W, file = paste0(save_path, 'rsq_resampling_W.dat'))


