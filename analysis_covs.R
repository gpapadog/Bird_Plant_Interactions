# Analyzing the birds and plants data
#
# Created: August 5, 2020.

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


# --------- In sample MCMC or including out of sample species ----------- #

use_A <- obs_A
use_n <- obs_n
use_X <- obs_X
use_W <- obs_W
use_Cu <- Cu
use_Cv <- Cv
use_nB <- nB
use_nP <- nP



# STEP 3: MCMC specifications.

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


# STEP 4: MCMC.


mcmc <- ModelCovariates(obs_A = use_A, obs_n = use_n, obs_X = use_X, obs_W = use_W,
                        Nsims = Nsims, burn = burn, thin = thin, bias_cor = bias_cor,
                        mh_n_pis = mh_n_pis, mh_n_pjs = mh_n_pjs, prior_mu0 = prior_mu0,
                        prior_sigmasq0 = prior_sigmasq0, prior_sigmasq = prior_sigmasq,
                        start_values = start_values, sampling = sampling)

attach(mcmc)

mean_pred <- apply(mcmc$Ls, c(2, 3), mean)
est_probL <- apply(mcmc$mod_pL1s, c(2, 3), mean)
mean_pL1s <- apply(pL1s, c(2, 3), mean)


gplots::heatmap.2(apply(Ls, c(2, 3), mean), trace = 'none',
                  dendrogram = 'none', Rowv = FALSE, Colv = FALSE)



plot(mean_pred[obs_A != 1], est_probL[obs_A != 1], pch = 16, cex = 0.3)
abline(a = 0, b = 1, col = 'red')
par(mfrow = c(1, 2))
hist(mean_pred[obs_A != 1] - est_probL[obs_A != 1], breaks = 100)
hist(est_probL[obs_A == 1], breaks = 100)


plot(mean_pred[obs_n == 0], est_probL[obs_n == 0], pch = 16, cex = 0.3)
abline(a = 0, b = 1, col = 'red')






# ---------- For the probability of interactions ----------- #

par(mfrow = c(1, 1))
hist(apply(logit(est_probL), c(2, 3), mean), main = 'Logit scale')









# ---------------- Predictions ----------------- #

pairs <- list(out_of_sample = list((nB - pred_nB + 1) : nB, (nP - pred_nP + 1) : nP),
              in_sample = list(1 : (nB - pred_nB), 1 : (nP - pred_nP)),
              half_in_sample = list(1 : nB, (nP - pred_nP + 1) : nP))

pred_values <- list(mean_pred = mean_pred,
                    mean_est_probL = apply(est_probL, c(2, 3), mean),
                    mean_pL1s = apply(pL1s, c(2, 3), mean),
                    true_probL = true_probL)

pred_power <- AllPredPower(obs_n = obs_n, pairs = pairs,
                           pred_values = pred_values, obs_A = obs_A,
                           true_L = true_L)$pred_power

pred_power[, , 2, 1, 1]


# Prediction of missing covariates.
miss_X_ind <- apply(use_X, 2, function(x) which(is.na(x)))
par(mfrow = c(2, floor((sum(pB) + 1) / 2)), mar = rep(2, 4))
for (pp in 1 : sum(pB)) {
  if (length(miss_X_ind[[pp]]) > 0) {
    plot(dta$X_lat[miss_X_ind[[pp]], pp],
         apply(Xs[[pp]], 2, mean, na.rm = TRUE), main = pp)
    abline(a = 0, b = 1)
  }
}

miss_W_ind <- apply(use_W, 2, function(x) which(is.na(x)))
par(mfrow = c(2, floor((sum(pP) + 1) / 2)), mar = rep(2, 4))
for (pp in 1 : sum(pB)) {
  if (length(miss_W_ind[[pp]]) > 0) {
    plot(dta$W_lat[miss_W_ind[[pp]], pp], apply(Ws[[pp]], 2, mean, na.rm = TRUE), main = pp)
    abline(a = 0, b = 1)
  }
}





# ---- For the probability of observing a bird -------- #

par(mfrow = c(1, 1))
est <- sweep(abind::abind(array(1, dim = c(dim(Us)[1 : 2], 1)), Us, along = 3),
             MARGIN = c(1, 3), FUN = '*', deltas)
plot(true_pis, apply(expit(apply(est, c(1, 2), sum)), 2, mean), xlim = c(0, 1), ylim = c(0, 1))
abline(a = 0, b = 1)


par(mfrow = c(1, 1))
plot(sigmasq_pB, type = 'l')
abline(h = true_sigmasq_pB, col = 'red')


# ---- For the probability of observing a plant -------- #

par(mfrow = c(1, 1))
est <- sweep(abind::abind(array(1, dim = c(dim(Vs)[1 : 2], 1)), Vs, along = 3),
             MARGIN = c(1, 3), FUN = '*', zetas)
plot(true_pjs, apply(expit(apply(est, c(1, 2), sum)), 2, mean), xlim = c(0, 1), ylim = c(0, 1))
abline(a = 0, b = 1)


par(mfrow = c(1, 1))
plot(sigmasq_pP, type = 'l')
abline(h = true_sigmasq_pP, col = 'red')


par(mfrow = c(2, 2), mar = rep(2, 4))
plot(true_pis[1 : use_nB], apply(pis, 2, mean), xlim = c(0, 1), ylim = c(0, 1))
abline(a = 0, b = 1)
plot(true_pjs[1 : use_nP], apply(pjs, 2, mean), xlim = c(0, 1), ylim = c(0, 1))
abline(a = 0, b = 1)

hist(pi_accepted, main = 'pi acceptance')
hist(pj_accepted, main = 'pj acceptance')


# ---------- For the correlation matrix of latent factors:
par(mfrow = c(2, 1), mar = rep(2, 4))
plot(rU, type = 'l')
abline(h = true_ru, col = 'red')
plot(rV, type = 'l')
abline(h = true_rv, col = 'red')

ru_accepted
rv_accepted

# --------------------------------------------------------------- #
# gplots::heatmap.2(apply(Ls, c(2, 3), mean), trace = 'none',
#                   dendrogram = 'none', Rowv = FALSE, Colv = FALSE)


pairs <- list(out_of_sample = list((nB - pred_nB + 1) : nB, (nP - pred_nP + 1) : nP),
              in_sample = list(1 : (nB - pred_nB), 1 : (nP - pred_nP)),
              half_in_sample = list(1 : (nB - pred_nB), (nP - pred_nP + 1) : nP))

pred_values <- list(predictive = mean_pred,
                    network_model = apply(est_probL, c(2, 3), mean),
                    mcmc_probs = apply(pL1s, c(2, 3), mean),
                    truth = true_probL)

pred_power <- AllPredPower(obs_n = obs_n, pairs = pairs,
                           pred_values = pred_values, obs_A = obs_A,
                           true_L = true_L)$pred_power

wh_pp <- c(1, 2, 3)
wh_aa <- 2
wh_pp_nB <- 1:5
wh_pp_nP <- 1

pred_power[wh_pp, , wh_aa, wh_pp_nB, wh_pp_nP]



par(mfrow = c(2, 2))
wh_out <- 4
for (hh in 1 : use_H) {
  plot(taus_beta[, wh_out, hh], type = 'l', main = hh)
  abline(h = 1, col = 'red')
}



# For the shrinkage prior:


par(mfrow = c(2, use_H / 2))
for (hh in 1 : use_H) {
  plot(thetas[, hh], main = hh, type = 'l')
}

round(apply(thetas, 2, mean), 4)

table(apply(zs, 1, function(x) sum(x <= 1 : use_H)))


