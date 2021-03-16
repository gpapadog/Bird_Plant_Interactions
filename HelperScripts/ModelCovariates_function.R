# Alternative estimation approaches for birds and plants.
#
# This code perfoms MCMC for a model that uses the covariates directly in the
# interaction and probability of observing models.
# This model DOES NOT use latent factors.
#
#' @param obs_A Observed interaction matrix. Contains values of 0 and 1.
#' @param obs_n Matrix with the same dimensions as obs_A. Includes the number
#' of studies in which both species were observed.
#' @param obs_X Matrix of observed covariates for the first set of species.
#' Rows correspond to species, columns to covariates. Continuous covariates
#' should be first, and binary covariates should follow.
#' @param obs_W Same structure as obs_X but for the second set of species.
#' @param Nsims Number of posterior samples to be kept.
#' @param burn Number of samples to be burnt in the beginning of the MCMC.
#' @param thin MCMC sampling thinning.
#' @param bias_cor Logical. Whether the model should aim to correct for
#' geographical and taxonomical bias. If set to FALSE, the MCMC is simplied
#' since the observed matrix of interactions is considered the same as the
#' true matrix of interactions. Defaults to TRUE.
#' @param mh_n_pis Parameter n in the Beta proposal for updating pi. Defaults
#' to 100.
#' @param mh_n_pjs Same as mh_n_pis but for the update of pj.
#' @param prior_mu0 Mean of the normal prior distribution for all coefficients.
#' Defaults to 0.
#' @param prior_sigmasq0 Variance of the normal prior distribution for all
#' coefficients. Defaults to 10.
#' @param prior_sigmasq Hyperparameters of the inverse gamma prior on the
#' variance terms. Defaults to (1, 1).
#' @param start_values List that can include starting values for any subset of
#' parameters. Starting values for parameters that are not specified in
#' start_values will be sampled. Defaults to NULL.
#' @param sampling List specifying which parameters should be sampled by
#' setting the value to TRUE, and which should not and should be kept at
#' their starting value by setting the corresponding list element to FALSE.
#' If sampling is set to FALSE for a parameter, it is recommended that the
#' corresponding start_value is specified to be set to the parameter's true
#' value. Defaults to NULL, and when set to NULL all parameters are sampled.
#' 
ModelCovariates <- function(obs_A, obs_n, obs_X, obs_W, Nsims, burn, thin,
                            bias_cor = TRUE, mh_n_pis = 100, mh_n_pjs = 100,
                            prior_mu0 = 0, prior_sigmasq0 = 10,
                            prior_sigmasq = c(1, 1),
                            start_values = NULL, sampling = NULL) {
  
  
  if (is.null(sampling)) {
    sampling <- list(L = TRUE, inter_coef = TRUE, sigmasq_p = TRUE,
                     delta = TRUE, zeta = TRUE, pis = TRUE, pjs = TRUE,
                     miss_X = TRUE, miss_W = TRUE)
  }
  if (!bias_cor) {
    cat('Without bias correction a number of parameters will not be sampled.', fill = TRUE)
    sampling$L <- FALSE
    sampling$sigmasq_p <- FALSE
    sampling$delta <- FALSE
    sampling$zeta <- FALSE
    sampling$pis <- FALSE
    sampling$pjs <- FALSE
  }
  
  # Sample size
  nB <- nrow(obs_A)
  nP <- ncol(obs_A)
  
  # Which A indices are equal to 0.
  index_A0 <- which(obs_A == 0)
  quant_A0 <- length(index_A0)
  
  # Continuous and binary covariates.
  entries_X <- apply(obs_X, 2, function(x) length(unique(x[!is.na(x)])))
  pB <- c(sum(entries_X > 2), sum(entries_X == 2))
  entries_W <- apply(obs_W, 2, function(x) length(unique(x[!is.na(x)])))
  pP <- c(sum(entries_W > 2), sum(entries_W == 2))
  # Making sure covariates are ordered as continuous first.
  if (any(c(entries_X[1 : pB[1]] == 2, entries_W[1 : pP[1]] == 2))) {
    stop('Reorder covariates')
  }
  
  # If the covariates have missing values, we check for it here.
  any_X_miss <- any(apply(obs_X, 2, function(x) sum(is.na(x))) > 0)
  any_W_miss <- any(apply(obs_W, 2, function(x) sum(is.na(x))) > 0)
  
  sampling$miss_X <- sampling$miss_X & any_X_miss
  sampling$miss_W <- sampling$miss_W & any_W_miss
  
  # Indices with missing values:
  if (any_X_miss) {
    miss_X_ind <- apply(obs_X, 2, function(x) which(is.na(x)))
  }
  if (any_W_miss) {
    miss_W_ind <- apply(obs_W, 2, function(x) which(is.na(x)))
  }
  
  
  
  # -------- PART 2 --------- #
  # Creating the elements where the MCMC samples will be saved.
  
  Ls <- array(NA, dim = c(Nsims, nB, nP))  # True interaction matrix.
  mod_pL1s <- array(NA, dim = c(Nsims, nB, nP))  # Probability of interaction without correction.
  pL1s <- array(NA, dim = c(Nsims, nB, nP))  # Probability of interaction.
  inter_coefs <- array(NA, dim = c(Nsims, sum(pB) + sum(pP) + 1))  # Coefficients, network model.
  deltas <- array(NA, dim = c(Nsims, sum(pB) + 1))  # Coefficients, observing bird.
  zetas <- array(NA, dim = c(Nsims, sum(pP) + 1))  # Coefficients, observing plant.
  mod_pis <- array(NA, dim = c(Nsims, nB))  # Probabilities of observing bird according to mean.
  pis <- array(NA, dim = c(Nsims, nB))  # Probabilities of observing bird.
  pi_accepted <- rep(0, nB)  # Number of MCMC iterations for which the bird proposal was accepted.
  mod_pjs <- array(NA, dim = c(Nsims, nP))  # Probabilities of observing bird according to mean.
  pjs <- array(NA, dim = c(Nsims, nP))  # Probabilities of observing plant.
  pj_accepted <- rep(0, nP)  # Number of MCMC iterations for which the plant proposal was accepted.
  sigmasq_pB <- rep(NA, Nsims)  # Residual variance, observing bird.
  sigmasq_pP <- rep(NA, Nsims)  # Residual variance, observing plant.
  
  if (any_X_miss) {
    beta0s <- array(NA, dim = c(Nsims, sum(pB)))
    sigmasq_m <- array(NA, dim = c(Nsims, pB[1]))
    Xs <- lapply(miss_X_ind, function(x) matrix(NA, nrow = Nsims, ncol = length(x)))
  }
  if (any_W_miss) {
    gamma0s <- array(NA, dim = c(Nsims, sum(pP)))
    sigmasq_l <- array(NA, dim = c(Nsims, pP[1]))
    Ws <- lapply(miss_W_ind, function(x) matrix(NA, nrow = Nsims, ncol = length(x)))
  }
  
  
  # --------- PART 3 --------- #
  # Setting starting values for the parameters.
  
  this_L <- matrix(rbinom(nB * nP, 1, 1 / 2), nB, nP)
  this_L[obs_A == 1] <- 1
  this_inter_coef <- rnorm(sum(pB) + sum(pP) + 1, 0, 1)
  this_delta <- rnorm(sum(pB) + 1, 0, 1)
  this_zeta <- rnorm(sum(pP) + 1, 0, 1)
  this_sigmasq_pB <- 1 / rgamma(1, 10, 10)
  this_sigmasq_pP <- 1 / rgamma(1, 10, 10)
  this_mod_pi <- rep(NA, nB)
  this_mod_pj <- rep(NA, nP)
  this_pi <- runif(nB, 0, 1)
  this_pj <- runif(nP, 0, 1)
  this_X <- obs_X
  this_W <- obs_W
  this_mod_pL1 <- matrix(NA, nrow = nB, ncol = nP)
  this_pL1 <- matrix(NA, nrow = nB, ncol = nP)
  
  # Draw starting values for missing covariates from observed distribution.
  if (any_X_miss) {
    
    this_beta0 <- rep(NA, sum(pB))
    this_sigmasq_m <- rep(NA, pB[1])
    
    for (jj in 1 : sum(pB)) {
      this_covar <- obs_X[, jj]
      this_covar <- this_covar[!is.na(this_covar)]
      use_stats <- c(mean(this_covar), sd(this_covar))
      
      if (jj > pB[1]) {  # Binary covariate
        this_X[miss_X_ind[[jj]], jj] <- rbinom(length(miss_X_ind[[jj]]), 1, use_stats[1])
        this_beta0[jj] <- logit(rbeta(1, use_stats[1] * 100, (1 - use_stats[1]) * 100))
      } else {  # Continuous covariate
        this_X[miss_X_ind[[jj]], jj] <- rnorm(length(miss_X_ind[[jj]]),
                                              use_stats[1], use_stats[2])
        this_beta0[jj] <- rnorm(1, mean = use_stats[1], sd = 0.5)
        this_sigmasq_m[jj] <- 1 / rgamma(1, 30, 30 * use_stats[2] ^ 2)
      }
    }
  }
  
  if (any_W_miss) {
    
    this_gamma0 <- rep(NA, sum(pP))
    this_sigmasq_l <- rep(NA, pP[1])
    
    for (jj in 1 : sum(pP)) {
      this_covar <- obs_W[, jj]
      this_covar <- this_covar[!is.na(this_covar)]
      use_stats <- c(mean(this_covar), sd(this_covar))
      
      if (jj > pP[1]) {
        this_W[miss_W_ind[[jj]], jj] <- rbinom(length(miss_W_ind[[jj]]), 1, use_stats[1])
        this_gamma0[jj] <- logit(rbeta(1, use_stats[1] * 100, (1 - use_stats[1]) * 100))
      } else {
        this_W[miss_W_ind[[jj]], jj] <- rnorm(length(miss_W_ind[[jj]]),
                                              use_stats[1], use_stats[2])
        this_gamma0[jj] <- rnorm(1, mean = use_stats[1], sd = 0.5)
        this_sigmasq_l[jj] <- 1 / rgamma(1, 30, 30 * use_stats[2] ^ 2)
      }
    }
  }
  
  # If the starting values have been specified in start_values, use those.
  if (!is.null(start_values)) {
    for (pp in 1 : length(start_values)) {
      assign(x = names(start_values)[pp], start_values[[pp]])
    }
  }
  
  if (length(this_inter_coef) != (sum(pB) + sum(pP) + 1)) {
    stop('Wrong length of interaction coefficient vector.')
  }
  
  # If we do not perform bias correction, change the starting values:
  if (!bias_cor) {
    this_L <- obs_A
    this_delta <- rep(NA, sum(pB) + 1)
    this_zeta <- rep(NA, sum(pP) + 1)
    this_sigmasq_pB <- NA
    this_sigmasq_pP <- NA
    this_pi <- rep(1, nB)
    this_pj <- rep(1, nP)
    obs_n <- matrix(NA, nB, nP)  # Does not inform anything.
  }
  
  
  # --------- PART 4 ----------- #
  # Performing the MCMC.
  
  # Which index we are saving at. Increased by 1 every thin iterations.
  keep_index <- 1
  
  cat('Total number of iterations:', Nsims * thin + burn, fill = TRUE)
  
  for (ss in 1 : (Nsims * thin + burn)) {
    
    if (ss %% 100 == 0) print(ss)
    
    
    # ------ L: Update true interactions. -------- #
    
    if (sampling$L | sampling$inter_coef | sampling$miss_X | sampling$miss_W) {
      
      bX <- this_X %*% matrix(this_inter_coef[2 : (sum(pB) + 1)], ncol = 1)
      gW <- this_W %*% matrix(this_inter_coef[(sum(pB) + 2) : (sum(pB) + sum(pP) + 1)])
      logit_pijL <- matrix(bX, nB, nP) + matrix(gW, nB, nP, byrow = TRUE)
      logit_pijL <- logit_pijL + this_inter_coef[1]
      this_mod_pL1 <- expit(logit_pijL)
      
    }
    
    if (sampling$L) {
      
      # Probability of observing interaction (i, j).
      pipj <- outer(this_pi, this_pj)
      
      # Probability of L = 1, or L = 0, when A = 0.
      pL1 <- this_mod_pL1 * (1 - pipj) ^ obs_n
      pL0 <- 1 - this_mod_pL1
      this_pL1 <- pL1 / (pL0 + pL1)  # Standardize.
      
      # Update L.
      this_L <- matrix(1, nB, nP)
      this_L[index_A0] <- rbinom(n = quant_A0, 1, prob = this_pL1[index_A0])
      
    }
    
    
    
    # ------- inter_coef: Update the coefficients in network model ------- #
    
    if (sampling$inter_coef | sampling$miss_X | sampling$miss_W) {
      # Sample the Polya-Gamma random variables:
      omega_L <- BayesLogit::rpg(nB * nP, h = rep(1, nB * nP), z = as.numeric(logit_pijL))
    }
    
    
    if (sampling$inter_coef) {
      des_mat <- matrix(1, nB * nP)
      for (jj in 1 : sum(pB)) {
        des_mat <- cbind(des_mat, rep(this_X[, jj], nP))
      }
      for (jj in 1 : sum(pP)) {
        des_mat <- cbind(des_mat, rep(this_W[, jj], each = nB))
      }
      
      prior_m <- rep(prior_mu0, sum(pB) + sum(pP) + 1)
      prior_S_inv <- diag(sum(pB) + sum(pP) + 1) / prior_sigmasq0
      
      # Following line is faster but identical to
      # t(des_mat) %*% diag(omega_L) %*% des_mat
      new_S <- sweep(t(des_mat), 2, FUN = '*', omega_L) %*% des_mat
      new_S <- chol2inv(chol(new_S + prior_S_inv))
      new_m <- t(des_mat) %*% matrix(as.numeric(this_L - 1 / 2), ncol = 1)
      new_m <- new_m + prior_S_inv %*% prior_m
      new_m <- new_S %*% new_m
      
      # Draw new values for the coefficients of the interaction model:
      this_inter_coef <- mvnfast::rmvn(1, new_m, new_S)
      
    }
    
    
    
    # ----- sigmasq_p: Residual variance for probability of observing ------- #
    
    if (sampling$sigmasq_p) {
      new_a <- prior_sigmasq[1] + nB / 2
      resid <- logit(this_pi) - cbind(1, this_X) %*% matrix(this_delta, ncol = 1)
      new_b <- prior_sigmasq[2] + sum(resid ^ 2) / 2
      this_sigmasq_pB <- 1 / rgamma(1, shape = new_a, rate = new_b)
      
      new_a <- prior_sigmasq[1] + nP / 2
      resid <- logit(this_pj) - cbind(1, this_W) %*% matrix(this_zeta, ncol = 1)
      new_b <- prior_sigmasq[2] + sum(resid ^ 2) / 2
      this_sigmasq_pP <- 1 / rgamma(1, shape = new_a, rate = new_b)  
    }
    
    
    # ------- deltas, zetas: Coefficients of probability of observing ------- #
    
    if (sampling$delta) {
      # deltas.
      des_mat <- cbind(1, this_X)
      prior_m <- matrix(rep(prior_mu0, sum(pB) + 1), ncol = 1)
      prior_S_inv <- diag(sum(pB) + 1) / prior_sigmasq0
      
      new_S <- chol2inv(chol(t(des_mat) %*% des_mat / this_sigmasq_pB + prior_S_inv))
      new_m <- t(des_mat) %*% matrix(logit(this_pi), ncol = 1) / this_sigmasq_pB
      new_m <- new_S %*% (new_m + prior_S_inv %*% prior_m)
      this_delta <- mvnfast::rmvn(1, new_m, new_S)
    }
    
    if (sampling$zeta) {
      # zetas:
      des_mat <- cbind(1, this_W)
      prior_m <- matrix(rep(prior_mu0, sum(pP) + 1), ncol = 1)
      prior_S_inv <- diag(sum(pP) + 1) / prior_sigmasq0
      
      new_S <- chol2inv(chol(t(des_mat) %*% des_mat / this_sigmasq_pP + prior_S_inv))
      new_m <- t(des_mat) %*% matrix(logit(this_pj), ncol = 1) / this_sigmasq_pP
      new_m <- new_S %*% (new_m + prior_S_inv %*% prior_m)
      this_zeta <- mvnfast::rmvn(1, new_m, new_S)
    }
    
    
    # ------- pis, pjs: Probability of observing for birds and plants.
    
    if (sampling$pis) {
      upd_prob_obs <- UpdProbObs(probobs_curr = this_pi,
                                 probobs_others = this_pj,
                                 curr_inter = this_L,
                                 obs_inter = obs_A,
                                 mh_n = mh_n_pis,
                                 n_studies = obs_n,
                                 latfac = this_X,  # Latent factors are substituted by covariates.
                                 coefs_probobs = this_delta,
                                 var_probobs = this_sigmasq_pB)
      this_mod_pi <- upd_prob_obs$model_values
      this_pi <- upd_prob_obs$new_values
      pi_accepted <- pi_accepted + upd_prob_obs$accepted
    }
    
    if (sampling$pjs) {
      upd_prob_obs <- UpdProbObs(probobs_curr = this_pj,
                                 probobs_others = this_pi,
                                 curr_inter = t(this_L),
                                 obs_inter = t(obs_A),
                                 mh_n = mh_n_pjs,
                                 n_studies = t(obs_n),
                                 latfac = this_W,
                                 coefs_probobs = this_zeta,
                                 var_probobs = this_sigmasq_pP)
      this_mod_pj <- upd_prob_obs$model_values
      this_pj <- upd_prob_obs$new_values
      pj_accepted <- pj_accepted + upd_prob_obs$accepted
    }
    
    
    
    # -------- Missing covariate values for the first set of species --------- #
    
    omega_L <- matrix(omega_L, nrow = nB, ncol = nP)
    
    if (sampling$miss_X) {
      
      # For the continuous covariates:
      for (jj in 1 : pB[1]) {
        
        # Updating the mean and variance of the covariate:
        
        s_new <- 1 / (nB / this_sigmasq_m[jj] + 1 / prior_sigmasq0)
        m_new <- sum(this_X[, jj]) / this_sigmasq_m[jj] + prior_mu0 / prior_sigmasq0
        m_new <- s_new * m_new
        this_beta0[jj] <- rnorm(1, mean = m_new, sd = sqrt(s_new))
        
        a_new <- prior_sigmasq[1] + nB / 2
        b_new <- prior_sigmasq[2] + sum((this_X[, jj] - this_beta0[jj]) ^ 2) / 2
        this_sigmasq_m[jj] <- 1 / rgamma(1, shape = a_new, rate = b_new)
        
        # Imputing the missing values
        
        this_miss <- miss_X_ind[[jj]]
        if (length(this_miss) > 0) {
          
          # Getting the residual variance.
          s_new <- (this_inter_coef[jj + 1] ^ 2 * apply(omega_L[this_miss, ], 1, sum) +
                      1 / this_sigmasq_m[jj])
          if (bias_cor) {
            s_new <- s_new + this_delta[jj + 1] ^ 2 / this_sigmasq_pB
          }
          s_new <- 1 / s_new
          m_new <- rep(NA, length(this_miss))
          
          # Getting the mean for each missing observation:
          for (ii in 1 : length(this_miss)) {
            this_ii <- this_miss[ii]
            
            # For the interaction model:
            lin_pred <- sum(this_inter_coef[setdiff(1 : (sum(pB) + 1), jj + 1)] *
                              c(1, this_X[this_ii, - jj]))
            lin_pred <- lin_pred + this_W %*% matrix(this_inter_coef[- c(1 : (sum(pB) + 1))], ncol = 1)
            A <- this_inter_coef[jj + 1] * sum(this_L[this_ii, ] - 1 / 2 -
                                                 omega_L[this_ii, ] * lin_pred)
            
            # For the covariate model:
            A <- A + this_beta0[jj] / this_sigmasq_m[jj]
            
            # For the probability of observing (if used):
            if (bias_cor) {
              lin_pred <- sum(this_delta[- (jj + 1)] * c(1, this_X[this_ii, - jj]))
              A <- A + this_delta[jj + 1] * (logit(this_pi[this_ii]) - lin_pred) / this_sigmasq_pB
            }
            
            m_new[ii] <- A
          }
          m_new <- s_new * m_new
          this_X[this_miss, jj] <- rnorm(length(this_miss), mean = m_new, sd = sqrt(s_new))
        }
      }
      
      
      # For the binary covariates:
      for (jj in (pB[1] + 1) : sum(pB)) {
        
        this_miss <- miss_X_ind[[jj]]
        if (length(this_miss) > 0) {
          
          # Updating the parameters (logit of probability):
          
          omega_m <- BayesLogit::rpg(nB, h = rep(1, nB), z = as.numeric(this_beta0[jj]))
          s_new <- 1 / (sum(omega_m) + 1 / prior_sigmasq0)
          m_new <- s_new * (sum(this_X[, jj] - 1 / 2) + prior_mu0 / prior_sigmasq0)
          this_beta0[jj] <- rnorm(1, mean = m_new, sd = sqrt(s_new))
          
          # Imputing the missing values:
          
          p_new <- rep(NA, length(this_miss))
          
          for (ii in 1 : length(this_miss)) {
            
            this_ii <- this_miss[ii]
            x0 <- this_X[this_ii, ]
            x0[jj] <- 0
            x1 <- this_X[this_ii, ]
            x1[jj] <- 1
            
            # Part about probability of interaction:
            lin_pred_add <- this_W %*% matrix(this_inter_coef[- c( 1 : (sum(pB) + 1))])
            lin_pred0 <- sum(this_inter_coef[1 : (sum(pB) + 1)] * c(1, x0)) + as.numeric(lin_pred_add)
            lin_pred1 <- sum(this_inter_coef[1 : (sum(pB) + 1)] * c(1, x1)) + as.numeric(lin_pred_add)
            
            p_pred0 <- expit(lin_pred0)
            p_pred1 <- expit(lin_pred1)
            
            p_pred0[this_L[this_ii, ] == 0] <- 1 - p_pred0[this_L[this_ii, ] == 0]
            p_pred1[this_L[this_ii, ] == 0] <- 1 - p_pred1[this_L[this_ii, ] == 0]
            
            # Before multiplying all the probabilities, we perform some sort of standardization:
            standard <- 1 / sort(c(p_pred0, p_pred1))[1 : floor(length(p_pred0) / 4)]
            
            p0 <- prod(c(standard, p_pred0))
            p1 <- prod(c(standard, p_pred1))
            
            # Part about model for the covariate:
            p0 <- p0 * (1 - expit(this_beta0[jj]))
            p1 <- p1 * expit(this_beta0[jj])
            
            # Part about probability of observing (if used):
            if (bias_cor) {
              lin_pred0 <- sum(this_delta * c(1, x0))
              lin_pred1 <- sum(this_delta * c(1, x1))
              
              p0 <- p0 * dnorm(logit(this_pi[this_ii]), mean = lin_pred0,
                               sd = sqrt(this_sigmasq_pB))
              p1 <- p1 * dnorm(logit(this_pi[this_ii]), mean = lin_pred1,
                               sd = sqrt(this_sigmasq_pB))
            }
            p_new[ii] <- p1 / (p0 + p1)
          }
          this_X[this_miss, jj] <- rbinom(length(this_miss), 1, prob = p_new)
        }
      }
    }
    
    
    
    # -------- Missing covariate values for the second set of species --------- #
    
    
    if (sampling$miss_W) {
      
      # For the continuous covariates:
      
      for (jj in 1 : pP[1]) {
        
        # Updating the mean and variance of the covariate:
        
        s_new <- 1 / (nP / this_sigmasq_l[jj] + 1 / prior_sigmasq0)
        m_new <- sum(this_W[, jj]) / this_sigmasq_l[jj] + prior_mu0 / prior_sigmasq0
        m_new <- s_new * m_new
        this_gamma0[jj] <- rnorm(1, mean = m_new, sd = sqrt(s_new))
        
        a_new <- prior_sigmasq[1] + nP / 2
        b_new <- prior_sigmasq[2] + sum((this_W[, jj] - this_gamma0[jj]) ^ 2) / 2
        this_sigmasq_l[jj] <- 1 / rgamma(1, shape = a_new, rate = b_new)
        
        # Imputing the missing values
        
        this_miss <- miss_W_ind[[jj]]
        if (length(this_miss) > 0) {
          
          # Getting the residual variance.
          s_new <- (this_inter_coef[1 + sum(pB) + jj] ^ 2 * apply(omega_L[, this_miss], 2, sum) +
                      1 / this_sigmasq_l[jj])
          if (bias_cor) {
            s_new <- s_new + this_zeta[jj + 1] ^ 2 / this_sigmasq_pP
          }
          s_new <- 1 / s_new
          m_new <- rep(NA, length(this_miss))
          
          # Getting the mean for each missing observation:
          for (ii in 1 : length(this_miss)) {
            this_ii <- this_miss[ii]
            
            # For the interaction model:
            lin_pred <- sum(this_inter_coef[1 + sum(pB) + setdiff(1 : sum(pP), jj)] * this_W[this_ii, - jj])
            lin_pred <- lin_pred + cbind(1, this_X) %*% matrix(this_inter_coef[1 : (sum(pB) + 1)], ncol = 1)
            A <- this_inter_coef[1 + sum(pB) + jj] * sum(this_L[, this_ii] - 1 / 2 -
                                                           omega_L[, this_ii] * lin_pred)
            
            # For the covariate model:
            A <- A + this_gamma0[jj] / this_sigmasq_l[jj]
            
            # For the probability of observing (if used):
            if (bias_cor) {
              lin_pred <- sum(this_zeta[- (jj + 1)] * c(1, this_W[this_ii, - jj]))
              A <- A + this_zeta[jj + 1] * (logit(this_pj[this_ii]) - lin_pred) / this_sigmasq_pP
            }
            
            m_new[ii] <- A
          }
          m_new <- s_new * m_new
          this_W[this_miss, jj] <- rnorm(length(this_miss), mean = m_new, sd = sqrt(s_new))
        }
      }
      
      
      # For the binary covariates:
      
      for (jj in (pP[1] + 1) : sum(pP)) {
        
        this_miss <- miss_W_ind[[jj]]
        if (length(this_miss) > 0) {
          
          # Updating the parameters:
          
          omega_m <- BayesLogit::rpg(nP, h = rep(1, nP), z = as.numeric(this_gamma0[jj]))
          s_new <- 1 / (sum(omega_m) + 1 / prior_sigmasq0)
          m_new <- s_new * (sum(this_W[, jj] - 1 / 2) + prior_mu0 / prior_sigmasq0)
          this_gamma0[jj] <- rnorm(1, mean = m_new, sd = sqrt(s_new))
          
          # Imputing the missing values:
          
          p_new <- rep(NA, length(this_miss))
          
          for (ii in 1 : length(this_miss)) {
            
            this_ii <- this_miss[ii]
            w0 <- this_W[this_ii, ]
            w0[jj] <- 0
            w1 <- this_W[this_ii, ]
            w1[jj] <- 1
            
            # Part about probability of interaction:
            lin_pred_add <- cbind(1, this_X) %*% matrix(this_inter_coef[1 : (sum(pB) + 1)], ncol = 1)
            lin_pred0 <- sum(this_inter_coef[1 + sum(pB) + 1 : sum(pP)] * w0) + as.numeric(lin_pred_add)
            lin_pred1 <- sum(this_inter_coef[1 + sum(pB) + 1 : sum(pP)] * w1) + as.numeric(lin_pred_add)
            
            p_pred0 <- expit(lin_pred0)
            p_pred1 <- expit(lin_pred1)
            
            p_pred0[this_L[, this_ii] == 0] <- 1 - p_pred0[this_L[, this_ii] == 0]
            p_pred1[this_L[, this_ii] == 0] <- 1 - p_pred1[this_L[, this_ii] == 0]
            
            # Before multiplying all the probabilities, we perform some sort of standardization:
            standard <- 1 / sort(c(p_pred0, p_pred1))[1 : floor(length(p_pred0) / 4)]
            
            p0 <- prod(c(standard, p_pred0))
            p1 <- prod(c(standard, p_pred1))
            
            # Part about model for the covariate:
            p0 <- p0 * (1 - expit(this_gamma0[jj]))
            p1 <- p1 * expit(this_gamma0[jj])
            
            # Part about probability of observing (if used):
            if (bias_cor) {
              lin_pred0 <- sum(this_zeta * c(1, w0))
              lin_pred1 <- sum(this_zeta * c(1, w1))
              
              p0 <- p0 * dnorm(logit(this_pj[this_ii]), mean = lin_pred0,
                               sd = sqrt(this_sigmasq_pP))
              p1 <- p1 * dnorm(logit(this_pj[this_ii]), mean = lin_pred1,
                               sd = sqrt(this_sigmasq_pP))
            }
            p_new[ii] <- p1 / (p0 + p1)
          }
          this_W[this_miss, jj] <- rbinom(length(this_miss), 1, prob = p_new)
        }
      }
    }
    
    
    # -------------- END OF MCMC UPDATES ------------- #
    
    
    # ------ Saving the results every thin iteration after burn in:
    if ((ss - burn) %% thin == 0 & ss > burn) {
      
      Ls[keep_index, , ] <- this_L
      pL1s[keep_index, , ] <- this_pL1
      mod_pL1s[keep_index, , ] <- this_mod_pL1
      inter_coefs[keep_index, ] <- this_inter_coef
      sigmasq_pB[keep_index] <- this_sigmasq_pB
      sigmasq_pP[keep_index] <- this_sigmasq_pP
      deltas[keep_index, ] <- this_delta
      zetas[keep_index, ] <- this_zeta
      mod_pis[keep_index, ] <- this_mod_pi
      pis[keep_index, ] <- this_pi
      mod_pjs[keep_index, ] <- this_mod_pj
      pjs[keep_index, ] <- this_pj
      
      if (any_X_miss) {
        beta0s[keep_index, ] <- this_beta0
        sigmasq_m[keep_index, ] <- this_sigmasq_m
        for (jj in 1 : sum(pB)) {
          Xs[[jj]][keep_index, ] <- this_X[miss_X_ind[[jj]], jj]
        }
      }
      if (any_W_miss) {
        gamma0s[keep_index, ] <- this_gamma0
        sigmasq_l[keep_index, ] <- this_sigmasq_l
        for (jj in 1 : sum(pP)) {
          Ws[[jj]][keep_index, ] <- this_W[miss_W_ind[[jj]], jj]
        }
      }
      
      # Increasing the index by 1.
      keep_index <- keep_index + 1
    }
    
  }
  
  
  # ---------- PART 5 ---------- #
  # Returning the results:
  
  r <- list(Ls = Ls, pL1s = pL1s, mod_pL1s = mod_pL1s,
            inter_coef = inter_coefs, sigmasq_pB = sigmasq_pB,
            sigmasq_pP = sigmasq_pP, deltas = deltas, zetas = zetas,
            pis = pis, pjs = pjs, mod_pis = mod_pis, mod_pjs = mod_pjs,
            pi_accepted = pi_accepted / (Nsims * thin + burn),
            pj_accepted = pj_accepted / (Nsims * thin + burn))
  
  if (any_X_miss) {
    r$Xs <- Xs
    r$beta0s <- beta0s
    r$sigmasq_m <- sigmasq_m
  }
  if (any_W_miss) {
    r$Ws <- Ws
    r$gamma0s <- gamma0s
    r$sigmasq_l <- sigmasq_l
  }
  
  return(r)
  
  
}