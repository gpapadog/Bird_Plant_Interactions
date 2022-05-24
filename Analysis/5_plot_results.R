# The directory where the analysis is performed:
wd_path <- 'Bird_Plant_Interactions/'
# Where the processed data are saved:
data_path <- 'Data/'
# Where the results are saved:
result_path <- 'Results/'


library(ggplot2)
library(RColorBrewer)
library(gplots)
library(superheat)
library(abind)
library(gridExtra)
library(grid)

setwd(wd_path)

# Loading the data:
load(paste0(data_path, 'bird_order_232.dat'))
load(paste0(data_path, 'bird_order_info_232.dat'))
load(paste0(data_path, 'plant_order.dat'))
load(paste0(data_path, 'plant_order_info.dat'))
load(paste0(data_path, 'obs_A.dat'))
load(paste0(data_path, 'obs_X.dat'))
load(paste0(data_path, 'obs_W.dat'))
load(paste0(data_path, 'birds_232.dat'))

# Restricting to the analysis of the 232 bird species:
wh_keep <- which(rownames(obs_A) %in% birds_232)
obs_A <- obs_A[wh_keep, , ]
obs_X <- obs_X[wh_keep, ]

# The combined network:
comb_A <- apply(obs_A, c(1, 2), sum)
comb_A <- (comb_A > 0) * 1

nB <- nrow(obs_A)
nP <- ncol(obs_A)

# Number of MCMC chains for our method and for the alternative method:
nchains <- 4
nchains_alt <- 3
# Number of cross validation repetitions:
repetitions <- 30

# Covariate names that are nicer for plotting:
good_namesX <- c('Body Mass', 'Gape Size', 'Large*', 'Fruit\nDependent*', 'Endangered*')
good_namesW <- c('Fruit\nDiameter', 'Fruit\nLength', 'Seed\nDiameter', 'Seed\nLength', 'Native*',
                 'Tree*', 'Black\nFruit*', 'Red\nFruit*', 'Yellow/Orange\nFruit*', 'Green\nFruit*',
                 'Lipid*', 'Endangered*')


# --------------- STEP 1: Getting the results together ----------------- #

all_res <- NULL
for (ii in 1 : nchains) {
  load(paste0(result_path, 'res_', ii, '.dat'))
  all_res[[ii]] <- res
}

alt_res <- NULL
for (ii in 1 : nchains_alt) {
  load(paste0(result_path, 'alt_res_', ii, '.dat'))
  alt_res[[ii]] <- res
}



# --------------- STEP 2: Plotting our analysis results ----------------- #

# Binding the posterior samples for the interactions from across chains.


# ----- Based on our model:

# Number of posterior samples used:
use_Nsims <- dim(all_res[[1]]$all_pred)[1]

# Creating an array to bind results across chains:
pred_ours <- array(NA, dim = c(nchains * use_Nsims, nB, nP))
for (ii in 1 : nchains) {
  # Using the posterior samples of the L matrix:
  pred_ours[1 : use_Nsims + use_Nsims * (ii - 1), , ] <- all_res[[ii]]$all_pred[, , , 1]
}
dimnames(pred_ours)[2 : 3] <- list(bird = rownames(obs_A), plant = colnames(obs_A))

# Calculating the posterior probability of interaction by averaging across
# posterior samples:
pred_ours <- apply(pred_ours, c(2, 3), mean)



# ----- Based on the alternative model:

# Number of posterior samples used:
use_Nsims <- dim(alt_res[[1]]$all_pred)[1]

# Creating an array to bind results across chains:
pred_alt <- array(NA, dim = c(nchains_alt * use_Nsims, nB, nP))
for (ii in 1 : nchains_alt) {
  # Using the posterior samples of the L matrix:
  pred_alt[1 : use_Nsims + use_Nsims * (ii - 1), , ] <- alt_res[[ii]]$all_pred[, , , 1]
}
dimnames(pred_alt)[2 : 3] <- list(bird = rownames(obs_A), plant = colnames(obs_A))
# Posterior probability of interaction:
pred_alt <- apply(pred_alt, c(2, 3), mean)


# Setting the recorded interactions to NA (so that they dont overpower the colors)
pred_ours[comb_A == 1] <- NA
pred_alt[comb_A == 1] <- NA

# Re-ordering according to the bird and plant order (in order to plot along
# taxonomic information)
pred_ours <- pred_ours[bird_order, plant_order]
pred_alt <- pred_alt[bird_order, plant_order]



# ----------- PART A: PLOTTING THE HEATMAP:

# Creating the clusters that will be used
# The following two lines specify that horizontal and vertical lines in our
# plot will separate species by taxonomic families:
bird_group <- bird_order_info$Frug_Family
plant_group <- as.character(plant_order_info$Plant_family)

# Calculating the size of each cluster, will be used when plotting results
# for families of certain size:
bird_size_cluster <- sapply(unique(bird_group), function(x) sum(bird_group == x))
plant_size_cluster <- sapply(unique(plant_group), function(x) sum(plant_group == x))

# Set plot_pred to pred_ours for results based on our method and to
# pred_alt for results based on the alternative method:
plot_pred <- pred_ours

# Set the minimum cluster size that should be plotted. For the results of the
# manuscript, we set min_bird_size to 10, and min_plant_size to 20. Setting both
# to 0 will produce the full results.
min_bird_size <- 10
min_plant_size <- 20

keep_bird_groups <- names(which(bird_size_cluster >= min_bird_size))
keep_plant_groups <- names(which(plant_size_cluster >= min_plant_size))

keep_bird_index <- which(bird_group %in% keep_bird_groups)
keep_plant_index <- which(plant_group %in% keep_plant_groups)


# Plotting those with minimum size as specified:
superheat(plot_pred[keep_bird_index, keep_plant_index],
          membership.rows = bird_group[keep_bird_index],
          membership.cols = plant_group[keep_plant_index],
          grid.hline.col = "#00257D", grid.vline.col = '#00257D',
          grid.hline.size = 0.3, grid.vline.size = 0.3,
          bottom.label.text.angle = 90,
          left.label.text.size = 3,
          bottom.label.text.size = 3,
          bottom.label.size = 0.2, left.label.size = 0.12,
          legend.breaks = seq(0, 1, by = 0.2),
          legend.vspace = 0.05,
          heat.col.scheme = "grey", heat.na.col = 'black',
          heat.pal.values = seq(0, 1, by = 0.05))



# ------------- PART B: Comparison of predictions across the two models:

# Combining the results in one array:
comb_pred <- abind::abind(pred_ours, pred_alt, along = 3)
dimnames(comb_pred)[[3]] <- c('Latent Factors', 'Covariates')

# Plotting histograms of predicted probability of interaction, grouped by the
# method used:
g1 <- ggplot(data = reshape2::melt(comb_pred),
             aes(x = value, group = Var3, fill = as.factor(Var3))) +
  geom_histogram(alpha = 0.6, position = 'identity', breaks = seq(0, 1, by = 0.05)) +
  theme_bw() +
  scale_fill_manual(name = 'Method', values = c('#F5352A', '#0C40A8')) +
  theme(legend.position = 'top') +
  xlab('Posterior probability of interaction\n')

# Plotting the predicted interaction probabilities against each other:
x <- as.data.frame(cbind(as.numeric(comb_pred[, , 1]), as.numeric(comb_pred[, , 2])))
names(x) <- c('Latent_Factors', 'Covariates')

g2 <- ggplot(data = x, aes(x = Latent_Factors, y = Covariates)) +
  geom_point(size = 0.001, alpha = 0.1) +
  theme_bw() +
  xlab('Posterior probability of interaction\nbased on latent factor model') +
  ylab('Posterior probability of interaction\nbased on covariate model')

# Any missing values exist because we have set the recorded interactions to NA,
# so that they don't skew our understanding of interaction prevalence.
gridExtra::grid.arrange(g1, g2, nrow = 1)



# --------------- STEP 3: Taxonomic correlation of latent factors ----------------- #

all_cor <- abind::abind(all_res[[1]]$correlations, all_res[[2]]$correlations, along = 3)
for (cc in 3 : nchains) {
  all_cor <- abind::abind(all_cor, all_res[[cc]]$correlations, along = 3)
}

# Posterior means and 95% credible intervals for the rho parameters in the
# latent factors for bird and plant species:
apply(all_cor, 2, mean)
apply(all_cor, 2, quantile, probs = c(0.025, 0.975))




# --------------- STEP 4: Cross validation results ----------------- #

# Getting the results together (held out indicies and predictions)
all_indices <- array(NA, dim = c(repetitions, 100, 2))
our_preds <- array(NA, dim = c(repetitions, nB, nP))
alt_preds <- array(NA, dim = c(repetitions, nB, nP))

for (rr in 1 : repetitions) {
  load(paste0(result_path, 'cv_indices_', rr, '.dat'))
  load(paste0(result_path, 'pred_', rr, '.dat'))
  load(paste0(result_path, 'alt_pred_', rr, '.dat'))
  all_indices[rr, , ] <- cv_indices
  our_preds[rr, , ] <- pred
  alt_preds[rr, , ] <- alt_pred
}

# Predictions of the held out data from both models:
pred <- array(NA, dim = c(repetitions, 100, 2))
for (rr in 1 : repetitions) {
  for (ii in 1 : 100) {
    pred[rr, ii, 1] <- our_preds[rr, all_indices[rr, ii, 1], all_indices[rr, ii, 2]]
    pred[rr, ii, 2] <- alt_preds[rr, all_indices[rr, ii, 1], all_indices[rr, ii, 2]]
  }
}

# Average and median probability of interaction based on the two models in the
# overall data:
overall_mean <- cbind(apply(our_preds, 1, mean), apply(alt_preds, 1, mean))
overall_median <- cbind(apply(our_preds, 1, median), apply(alt_preds, 1, median))

# Average and median in the held out data.
pred_mean <- apply(pred, c(1, 3), mean)
pred_median <- apply(pred, c(1, 3), median)

# Creating the data frame we will plot:
plot_dta <- abind::abind(pred_mean / overall_mean, pred_median / overall_median, along = 3)
dimnames(plot_dta)[2 : 3] <- list(method = c('Latent Factors', 'Covariates'),
                                  stat = c('Prediction mean / Overall mean',
                                           'Prediction median / Overall median'))
names(dimnames(plot_dta)) <- c('Iteration', 'Method', 'Statistic')
plot_dta <- reshape2::melt(plot_dta)

# Plotting cross validation results:
ggplot(data = plot_dta) +
  geom_boxplot(aes(x = Method, y = value)) +
  facet_wrap(~ Statistic, scales = 'free_y') +
  theme_bw() +
  ylab('') +
  xlab('') +
  ggtitle('Out of sample performance', subtitle = 'Predictions for held-out recorded interactions') +
  theme(legend.position = 'none') +
  scale_y_continuous(limits = function(x) c(1, x[2]), n.breaks = 6)



# --------------- STEP 5: Variable importance measure ----------------- #


# ------- PART A: Plotting the variables in order of importance:

load(paste0(result_path, 'rsq_obs_X.dat'))
load(paste0(result_path, 'rsq_obs_W.dat'))
load(paste0(result_path, 'rsq_resampling_X.dat'))
load(paste0(result_path, 'rsq_resampling_W.dat'))


# Calculating the number of permuted standard deviations away from the mean.

# Starting from the bird covariates:
wh_obs <- rsq_obs_X
wh_resampling <- rsq_resampling_X
sd_awayX <- rep(NA,  length(wh_obs))
names(sd_awayX) <- good_namesX
for  (cc in 1 : length(wh_obs)) {
  sd_awayX[cc] <- (wh_obs[cc] - mean(wh_resampling[, cc])) / sd(wh_resampling[, cc])
}

# And for the plant covariates:
wh_obs <- rsq_obs_W
wh_resampling <- rsq_resampling_W
sd_awayW <- rep(NA,  length(wh_obs))
names(sd_awayW) <- good_namesW
for  (cc in 1 : length(wh_obs)) {
  sd_awayW[cc] <- (wh_obs[cc] - mean(wh_resampling[, cc])) / sd(wh_resampling[, cc])
}


# Plotting the tiles of variable importance ordering the variables in
# decreasing importance:

# For the bird species:
xx <- data.frame(value = sd_awayX, covariate = names(sd_awayX), y = 1)
xx <- xx[order(- xx$value), ]
xx$covariate <- factor(xx$covariate, levels = xx$covariate)

ggplot() + geom_tile(aes(x = covariate, y = y, fill = value), color = 'white', data = xx) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid = element_blank()) +
  theme(axis.text.y = element_blank(), axis.ticks = element_blank()) +
  theme(legend.position = 'none', axis.title = element_blank()) +
  scale_fill_gradient(low = '#BFF0B6', high = '#3B6E32') +
  theme(axis.text = element_text(angle = 0, hjust = 0.5, vjust = 0, size = 8))

# For the plant species:
ww <- data.frame(value = sd_awayW, covariate = names(sd_awayW), y = 1)
ww <- ww[order(- ww$value), ]
ww$covariate <- factor(ww$covariate, levels = ww$covariate)

ggplot() + geom_tile(aes(x = covariate, y = y, fill = value), color = 'white', data = ww) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid = element_blank()) +
  theme(axis.text.y = element_blank(), axis.ticks = element_blank()) +
  theme(legend.position = 'none', axis.title = element_blank()) +
  scale_fill_gradient(low = '#D4DEF7', high = '#425075') +
  theme(axis.text = element_text(angle = 0, hjust = 0.5, vjust = 0, size = 8))




# ------- PART B: Posterior probabilities based on the important covariates.

# Taking the posterior probabilities of interaction across the three chains,
# and setting the recorded interactions to NA:
use_Ls <- do.call(abind, c(lapply(all_res, function(x) x$all_pred[, , , 1]), along = 1))
use_mean_Ls <- apply(use_Ls, c(2, 3), mean)
use_mean_Ls[comb_A == 1] <- NA

# Which covariate is to be plotted. The ones we want are listed first.
wh_X <- 1
wh_W <- 1

# Showing only the species that have the covariate measured.
keep_birds <- which(!is.na(obs_X[, wh_X]))
keep_plants <- which(!is.na(obs_W[, wh_W]))
use_out <- use_mean_Ls[keep_birds, keep_plants]

# Because some species have identical values for the covariate, in order for
# plot to show all of them, we need to slightly pertube their values. That way,
# the increasing or decreasing order is not altered, but there is no overlap in
# the covariate values:
bird_cov <- obs_X[keep_birds, wh_X]
bird_cov <- bird_cov + rnorm(length(keep_birds), sd = sd(bird_cov) * 0.0001)
plant_cov <- obs_W[keep_plants, wh_W]
plant_cov <- plant_cov + rnorm(length(plant_cov), sd = sd(plant_cov) * 0.0001)

# Creating a data frame in which the species are ordered by their covariate
# values. This will allow us to plot the probability of interaction across the
# covariates in an interpretable way. We also note that we need to turn the
# covariates to factors in order for them to be plotted in the correct order.
plot_dta <- data.frame(cov_bird = rep(bird_cov, length(keep_plants)),
                       cov_plant = rep(plant_cov, each = length(keep_birds)),
                       probability = as.numeric(use_out))
plot_dta$use_cov_bird <- factor(as.numeric(factor(plot_dta$cov_bird)))
plot_dta$use_cov_plant <- factor(as.numeric(factor(plot_dta$cov_plant)))

g <- ggplot() +
  geom_raster(aes(x = use_cov_plant, y = use_cov_bird, fill = probability), data = plot_dta) +
  scale_fill_gradient(low = "#F5D4C7", high = "#02A65F", na.value = '#016B3B',
                      name = 'Posterior\ninteraction\nprobability\n', limits = c(0, 1)) +
  theme(axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_text(size = 30),
        axis.title.y = element_text(vjust = - 1)) +
  ylab(expression(symbol('\256'))) + xlab(expression(symbol('\256')))

gridExtra::grid.arrange(g, left = textGrob("Bird information: Increasing Body Mass", rot = 90,
                                           x = 1.3, y = 0.57, gp = gpar(fontsize = 12)),
                        bottom = textGrob("Plant information: Increasing Fruit Diameter", 
                                          x = 0.435, y = 1.3, gp = gpar(fontsize = 12)),
                        vp=viewport(width=0.5, height=0.6))

