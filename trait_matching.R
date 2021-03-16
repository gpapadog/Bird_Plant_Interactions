# Performing cross-validation of bird-plant interaction models.
# Created: October 12, 2020.

rm(list = ls())
dev.off()

# ------ STEP 0: Some functions. --------- #

setwd('~/Documents/Research/Birds_and_plants/Analysis/')
source('~/Github/Birds_and_plants/functions/TraitMatching3_function.R')
source('~/Github/Birds_and_plants/functions/useful_functions.R')

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

wh_results <- 3
nchains <- 3


# --------------- STEP 1: Getting the results together ----------------- #

all_pred <- NULL
for (ii in 1 : nchains) {
  load(paste0('Output/1results/results', wh_results, '/res_', ii, '.dat'))
  all_pred[[ii]] <- res$all_pred
}

all_imputed <- NULL
for (ii in 1 : nchains) {
  load(paste0('Output/1results/results', wh_results, '/imputed_', ii, '.dat'))
  all_imputed[[ii]] <- imputed
}



Nsims <- dim(all_pred[[1]])[1]


mod_pL1s <- array(NA, dim = c(nchains * Nsims, nB, nP))
for (cc in 1 : 3) {
  wh_entries <- Nsims * (cc - 1) + 1 : Nsims
  mod_pL1s[wh_entries, , ] <- all_pred[[cc]][, , , 2]
}

Xs <- all_imputed[[1]]$Xs
for (cc in 2 : 3) {
  for (cov in 1 : 5) {
    Xs[[cov]] <- rbind(Xs[[cov]], all_imputed[[cc]]$Xs[[cov]])
  }
}

Ws <- all_imputed[[1]]$Ws
for (cc in 2 : 3) {
  for (cov in 1 : 12) {
    Ws[[cov]] <- rbind(Ws[[cov]], all_imputed[[cc]]$Ws[[cov]])
  }
}



# --------------- STEP 2: Performing trait matching ----------------- #


use_mod_pL1s <- mod_pL1s
use_mod_pL1s[use_mod_pL1s > 1 - 10^{-10}] <- 1 - 10^{-10}


trait_match3 <- TraitMatching3(B = 4, mod_pL1s = use_mod_pL1s,
                               Xs = Xs, Ws = Ws,
                               obs_X = obs_X, obs_W = obs_W)



