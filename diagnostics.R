# preparing diagnostic plots for convergence.

setwd('~/Documents/Research/Birds_and_plants/')
out_path <- 'Analysis/Output/1results/results3/'
nchains <- 3

load('~/Github/Birds_and_plants/Application/Data/Aves_analysis/obs_A.dat')
load('~/Github/Birds_and_plants/Application/Data/Aves_analysis/obs_W.dat')
source('~/Github/Birds_and_plants/functions/useful_functions.R')

library(ggplot2)
library(reshape2)

nB <- nrow(obs_A)
nP <- ncol(obs_A)

choose_pairs <- function(how_many) {

  chosen_pairs <- NULL
  chosen_pairs.df <- data.frame(bird = NULL, plant = NULL)
  while (length(chosen_pairs) < how_many) {
    chosen_bird <- sample(1 : nB, 1)
    chosen_plant <- sample(1 : nP, 1)
    if (obs_A[chosen_bird, chosen_plant] == 0) {
      name_pair <- paste0(chosen_bird, ', ', chosen_plant)
      if (!(name_pair %in% chosen_pairs)) {
        chosen_pairs <- c(chosen_pairs, name_pair)
        chosen_pairs.df <- rbind(chosen_pairs.df, data.frame(bird = chosen_bird, plant = chosen_plant))
      }
    }
  }
  return(list(pairs = chosen_pairs, df = chosen_pairs.df))
}


all_res <- NULL
for (ii in 1 : nchains) {
  print(ii)
  load(paste0(out_path, 'res_', ii, '.dat'))
  all_res[[ii]] <- res
  rm(res)
}

Nsims <- length(all_res[[1]]$sigmasq_pB)
pB <- dim(all_res[[1]]$linpred_X)[3]
pP <- dim(all_res[[1]]$linpred_W)[3]

# ---------------------------------------------------- #

# DIAGNOSTIC PLOT FOR PREDICTIONS
#
# I will plot the traceplot for the probability of interaction
# according to the model, and according to the model with bias
# correction

# ggplot was taking too long, so I will use base-R:


# For probL: Traceplots of randomly chosen pairs:

chosen_pairs <- choose_pairs(6)

png(file = paste0(out_path, 'Diagnostics/diag_probL.png'), width = 700, height = 400)
par(mfrow = c(2, 3), mar = c(1, 1, 2, 1), oma = c(0, 0, 2, 0))
for (ii in 1 : 6) {
  plot(1, xlim = c(1, Nsims), ylim = c(0, 1), type = 'n', main = chosen_pairs$pairs[ii], axes = F)
  for (cc in 1 : nchains) {
    lines(1 : Nsims, all_res[[cc]]$all_pred[, chosen_pairs$df[ii, 1], chosen_pairs$df[ii, 2], 2], col = cc)
  }
}
title(main = 'Probability of interaction', outer = TRUE, line = 0)
dev.off()


chosen_pairs <- choose_pairs(6)

png(file = paste0(out_path, 'Diagnostics/diag_pL1.png'), width = 700, height = 400)
par(mfrow = c(2, 3), mar = c(1, 1, 2, 1), oma = c(0, 0, 2, 0))
for (ii in 1 : 6) {
  plot(1, xlim = c(1, Nsims), ylim = c(0, 1), type = 'n', main = chosen_pairs$pairs[ii], axes = F)
  axis(2)
  for (cc in 1 : nchains) {
    lines(1 : Nsims, all_res[[cc]]$all_pred[, chosen_pairs$df[ii, 1], chosen_pairs$df[ii, 2], 3], col = cc)
  }
}
title(main = 'Probability of interaction with bias correction', outer = TRUE, line = 0)
dev.off()


# ---------------------------------------------------- #

mean_pred1 <- sapply(all_res, function(x) apply(x$all_pred[, , , 1], c(2, 3), mean))
mean_pred1 <- array(mean_pred1, dim = c(nB, nP, nchains))

mean_pred2 <- sapply(all_res, function(x) apply(x$all_pred[, , , 2], c(2, 3), mean))
mean_pred2 <- array(mean_pred2, dim = c(nB, nP, nchains))

mean_pred3 <- sapply(all_res, function(x) apply(x$all_pred[, , , 3], c(2, 3), mean))
mean_pred3 <- array(mean_pred3, dim = c(nB, nP, nchains))

plot(mean_pred1[obs_A == 0], mean_pred2[obs_A == 0], col = rep(1 : 5, each = nB * nP), pch = 16, cex = 0.1)
plot(mean_pred1[obs_A == 0], mean_pred3[obs_A == 0], col = rep(1 : 5, each = nB * nP), pch = 16, cex = 0.1)
plot(mean_pred2[obs_A == 0], mean_pred3[obs_A == 0], col = rep(1 : 5, each = nB * nP), pch = 16, cex = 0.1)


png(file = paste0(out_path, 'Diagnostics/diag_mean_pL1.png'), width = 700, height = 400)
plot(mean_pred1[, , 1], mean_pred1[, , 2], pch = 16, cex = 0.1)
for (cc in 3 : nchains) {
  points(mean_pred1[, , 1], mean_pred1[, , cc], col = cc, pch = 16, cex = 0.1)
  Sys.sleep(3)
}
dev.off()



pvals <- apply(mean_pred1, c(1, 2), function(x) chisq.test(x = rbind(x * Nsims, (1 - x) * Nsims))$p.value)
pvals <- as.numeric(pvals[obs_A == 0])

png(file = paste0(out_path, 'Diagnostics/diag_pvals_QQplot.png'), width = 800, height = 500)
EnvStats::qqPlot(pvals, distribution = 'unif', param.list = list(min = 0, max = 1))
abline(a = 0, b = 1, col = 'red')
dev.off()



# ---------------------------------------------------- #

# DIAGNOSTIC PLOT FOR Probability of observation
#

# For BIRDS:

# Merging predictions from the different chains:
plot_dta <- array(NA, dim = c(nchains, dim(all_res[[1]]$prob_obsB)))
for (ii in 1 : nchains) {
  plot_dta[ii, , , ] <- all_res[[ii]]$prob_obsB
}
dimnames(plot_dta) <- list(chain = 1 : nchains,
                           sim = 1 : dim(plot_dta)[2],
                           bird = 1 : nB,
                           quant = c('Model probability', 'Samples'))
names(dimnames(plot_dta)) <- c('chain', 'sim', 'bird', 'quant')

chosen_birds <- sample(1 : nB, 4)
plot_dta <- plot_dta[, , chosen_birds, ]
plot_dta <- reshape2::melt(plot_dta)
plot_dta$chain <- factor(plot_dta$chain)

png(file = paste0(out_path, 'Diagnostics/diag_prob_obsB.png'), width = 600, height = 300)
ggplot() +
  geom_line(aes(x = sim, y = value, group = chain, color = chain), data = plot_dta, alpha = 0.5) +
  ylab('Probability of observing a bird') +  xlab('') +
  theme_bw() +
  facet_grid(quant ~ bird) 
dev.off()


# FOR PLANTS:

plot_dta <- array(NA, dim = c(nchains, dim(all_res[[1]]$prob_obsP)))
for (ii in 1 : nchains) {
  plot_dta[ii, , , ] <- all_res[[ii]]$prob_obsP
}
dimnames(plot_dta) <- list(chain = 1 : nchains,
                           sim = 1 : dim(plot_dta)[2],
                           plant = 1 : nP,
                           quant = c('Model probability', 'Samples'))
names(dimnames(plot_dta)) <- c('chain', 'sim', 'plant', 'quant')

chosen_plants <- sample(1 : nP, 4)
plot_dta <- plot_dta[, , chosen_plants, ]
plot_dta <- reshape2::melt(plot_dta)
plot_dta$chain <- factor(plot_dta$chain)


png(file = paste0(out_path, 'Diagnostics/diag_prob_obsP.png'), width = 600, height = 300)
ggplot() +
  geom_line(aes(x = sim, y = value, group = chain, color = chain), data = plot_dta, alpha = 0.5) +
  ylab('Probability of observing a plant') + xlab('') +
  theme_bw() +
  facet_grid(quant ~ plant)
dev.off()


# RESIDUAL VARIANCES

plot_dta <- array(NA, dim = c(nchains, length(all_res[[1]]$sigmasq_pB), 2))
for (cc in 1 : nchains) {
  plot_dta[cc, , 1] <- all_res[[cc]]$sigmasq_pB
  plot_dta[cc, , 2] <- all_res[[cc]]$sigmasq_pP
}
dimnames(plot_dta) <- list(1 : nchains, 1 : Nsims, c('birds', 'plants'))
names(dimnames(plot_dta)) <- c('chain', 'sim', 'species')
plot_dta <- reshape2::melt(plot_dta)
plot_dta$chain <- factor(plot_dta$chain)

png(file = paste0(out_path, 'Diagnostics/diag_sigmasq_p.png'), width = 500, height = 200)
ggplot() +
  geom_line(aes(x = sim, y = value, group = chain, color = chain), alpha = 0.5, data = plot_dta) +
  ylab('Residual variance for\ndetection probability model') + xlab('') +
  theme_bw() +
  facet_wrap(~ species, nrow = 1, scales = 'free_y')
dev.off()

# ---------------------------------------------------- #



# ---------------------------------------------------- #

# DIAGNOSTIC PLOT FOR COVARIATE MODELS

# FOR BIRDS:

# Linear predictor:
chosen_birds <- sample(1 : nB, 4)

plot_dta <- array(NA, dim = c(nchains, Nsims, length(chosen_birds), pB))
for (ii in 1 : nchains) {
  plot_dta[ii, , , ] <- all_res[[ii]]$linpred_X[, chosen_birds, ]
}
dimnames(plot_dta) <- list(1 : nchains, 1 : Nsims, chosen_birds, 1 : dim(plot_dta)[4])
names(dimnames(plot_dta)) <- c('chain', 'sim', 'bird', 'cov')
plot_dta <- reshape2::melt(plot_dta)
plot_dta$chain <- factor(plot_dta$chain)
plot_dta$use_facet <- paste(plot_dta$bird, '-', plot_dta$cov)
plot_dta <- plot_dta[order(plot_dta$bird), ]
plot_dta$use_facet <- factor(plot_dta$use_facet, levels = unique(plot_dta$use_facet))

png(file = paste0(out_path, 'Diagnostics/diag_linpred_X.png'), width = 1000, height = 500)
ggplot() +
  geom_line(aes(x = sim, y = value, group = chain, color = chain), data = plot_dta, alpha = 0.5) +
  ylab('Linear predictor for bird covariate model') + xlab('') +
  theme_bw() +
  facet_wrap(~ use_facet, scale = 'free_y', nrow = length(chosen_plants)) +
  theme(legend.position = 'none')
dev.off()

# Residual variance:
plot_dta <- array(NA, dim = c(nchains, dim(all_res[[1]]$sigmasq_m)))
for (ii in 1 : nchains) {
  plot_dta[ii, , ] <- all_res[[ii]]$sigmasq_m
}
dimnames(plot_dta) <- list(1 : nchains, 1 : Nsims, cov = 1 : dim(plot_dta)[3])
names(dimnames(plot_dta)) <- c('chain', 'sim', 'cov')
plot_dta <- reshape2::melt(plot_dta)
plot_dta$chain <- factor(plot_dta$chain)


png(file = paste0(out_path, 'Diagnostics/diag_sigmasq_m.png'), width = 600, height = 200)
ggplot() +
  geom_line(aes(x = sim, y = value, group = chain, color = chain),
            data = plot_dta, alpha = 0.5) +
  ylab('Residual variance\ntrait models for birds') + xlab('') +
  theme_bw() +
  facet_wrap(~ cov, scales = 'free_y')
dev.off()



# FOR PLANTS:

# Linear predictor:
chosen_plants <- sample(1 : nP, 4)

plot_dta <- array(NA, dim = c(nchains, Nsims, length(chosen_plants), pP))
for (ii in 1 : nchains) {
  plot_dta[ii, , , ] <- all_res[[ii]]$linpred_W[, chosen_plants, ]
}
dimnames(plot_dta) <- list(1 : nchains, 1 : Nsims, chosen_plants, 1 : dim(plot_dta)[4])
names(dimnames(plot_dta)) <- c('chain', 'sim', 'plant', 'cov')
plot_dta <- reshape2::melt(plot_dta)
plot_dta$chain <- factor(plot_dta$chain)
plot_dta$use_facet <- paste(plot_dta$plant, '-', plot_dta$cov)
plot_dta <- plot_dta[order(plot_dta$plant), ]
plot_dta$use_facet <- factor(plot_dta$use_facet, levels = unique(plot_dta$use_facet))

# plot_dta <- subset(plot_dta, cov %% 2 == 1)
plot_dta <- subset(plot_dta, cov %% 2 == 0)
# 
# png(file = paste0(out_path, 'Diagnostics/diag_linpred_W.png'), width = 950, height = 500)
ggplot() +
  geom_line(aes(x = sim, y = value, group = chain, color = chain), data = plot_dta, alpha = 0.5) +
  ylab('Linear predictor for plant covariate model') + xlab('') +
  theme_bw() +
  facet_wrap(~ use_facet, scale = 'free_y', nrow = length(chosen_plants)) +
  theme(legend.position = 'none')
dev.off()



plot_dta <- array(NA, dim = c(nchains, dim(all_res[[1]]$sigmasq_l)))
for (ii in 1 : nchains) {
  plot_dta[ii, , ] <- all_res[[ii]]$sigmasq_l
}
dimnames(plot_dta) <- list(1 : nchains, 1 : Nsims, cov = 1 : dim(plot_dta)[3])
names(dimnames(plot_dta)) <- c('chain', 'sim', 'cov')
plot_dta <- reshape2::melt(plot_dta)
plot_dta$chain <- factor(plot_dta$chain)

png(file = paste0(out_path, 'Diagnostics/diag_sigmasq_l.png'), width = 600, height = 150)
ggplot() +
  geom_line(aes(x = sim, y = value, group = chain, color = chain),
            data = plot_dta, alpha = 0.5) +
  ylab('Residual variance for\ntrait models for plants') + xlab('') +
  theme_bw() +
  theme(legend.position = 'none') +
  facet_wrap(~ cov, nrow = 1)
dev.off()



# ---------------------------------------------------- #

# DIAGNOSTIC PLOT FOR CORRELATIONS

plot_dta <- array(NA, dim = c(nchains, Nsims, 2))
dimnames(plot_dta) <- list(1 : nchains, 1 : Nsims, c('birds', 'plants'))
names(dimnames(plot_dta)) <- c('chain', 'sim', 'species')
for (ii in 1 : nchains) {
  plot_dta[ii, , ] <- all_res[[ii]]$correlations
}
plot_dta <- reshape2::melt(plot_dta)
plot_dta$chain <- factor(plot_dta$chain)

png(file = paste0(out_path, 'Diagnostics/diag_correlations.png'), width = 600, height = 200)
ggplot() + 
  geom_line(aes(x = sim, y = value, group = chain, color = chain), data = plot_dta, alpha = 0.5) +
  ylab('Correlation in latent factors') + xlab('') +
  theme_bw() +
  theme(legend.position = 'none') +
  facet_wrap(~ species)
dev.off()



# ----------- RUNNING MEANS FOR THE PROBABILITY OF INTERACTION -------- #

chosen_pairs <- choose_pairs(40)
run_means <- array(NA, dim = c(Nsims, nB, nP, nchains))

for (ii in 1 : nrow(chosen_pairs$df)) {
  bird_index <- chosen_pairs$df[ii, 1]
  plant_index <- chosen_pairs$df[ii, 2]
  
  for (cc in 1 : nchains) {
    run_means[, bird_index, plant_index, cc] <- cumsum(all_res[[cc]]$all_pred[, bird_index, plant_index, 1]) / (1 : Nsims)
  }
}



plot_num <- 9
reordering <- sample(1 : nrow(chosen_pairs$df), nrow(chosen_pairs$df))
chosen_pairs$pairs <- chosen_pairs$pairs[reordering]
chosen_pairs$df <- chosen_pairs$df[reordering, ]


png(file = paste0(out_path, 'Diagnostics/diag_prob_inter.png'), width = 500, height = 350)
par(mfrow = rep(sqrt(plot_num), 2), mar = rep(2, 4), oma = c(0, 0, 2, 0))
tt <- 1
for (ii in 1 : 40) {
  bird_index <- chosen_pairs$df[ii, 1]
  plant_index <- chosen_pairs$df[ii, 2]

  if (obs_A[bird_index, plant_index] == 0 & tt <= plot_num) {
    use_title <- paste0('(', chosen_pairs$pairs[ii], ')')
    plot(1, type = 'n',xlim = c(1, Nsims), ylim = c(0, 1), main = use_title)
    for (cc in 1 : nchains) {
      lines(1 : Nsims, run_means[, bird_index, plant_index, cc], col = cc)
    }
    tt <- tt + 1
  }
}
title(main = 'Probability of interaction', line = 0.5, outer = TRUE)
dev.off()




