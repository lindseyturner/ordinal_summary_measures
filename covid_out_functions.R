# functions for covid out data analysis 

# same as the simulation function, but output data in a list instead of a vector 
# new function for weighted risk difference 
# assumes a partial proportional odds model 
# assumes the treatment is coded as -0.5 for control and 0.5 for treatment 
order_weighted_RD2 <- function(b = b, iprior = 1, dat, trt, outcome) { 
  # first calculate P(y >= i | ctl) and P(y >= i | trt) for all i 
  PPO_ctl_geq <- NA 
  PPO_trt_geq <- NA 
  PPO_ctl_leq <- NA 
  PPO_trt_leq <- NA 
  p_y_ctl <- NA
  p_y_trt <- NA
  p_y_ctl_est <- NA
  p_y_trt_est <- NA
  p_y <- NA
  c <- length(b$ylevels)
  
  if (iprior == 1) {
    new_coef <- c(-coef(b)[1:(c-1)], coef(b)[c:length(coef(b))])
  } else {
    new_coef <- coef(b)
  }
  
  for (i in 2:c) {
    # for the first level there is not a partial term, so only need first term and metformin term 
    if (i == 2) {
      PPO_ctl_geq[i] <- plogis(new_coef[1] - 0.5 * new_coef[c]) 
      PPO_trt_geq[i] <- plogis(new_coef[1] + 0.5 * new_coef[c])
      PPO_ctl_leq[i - 1] <- 1 - PPO_ctl_geq[i]
      PPO_trt_leq[i - 1] <- 1 - PPO_trt_geq[i]
    } else{
      phrase <- paste0("y>=", i)
      coefs <- new_coef[str_detect(names(new_coef), phrase)]
      PPO_ctl_geq[i] <- plogis(coefs[1] - 0.5 * coefs[2] - 0.5 * coef(b)[c]) 
      PPO_trt_geq[i] <- plogis(coefs[1] + 0.5 * coefs[2] + 0.5 * coef(b)[c])
      PPO_ctl_leq[i - 1] <- 1 - PPO_ctl_geq[i]
      PPO_trt_leq[i - 1] <- 1 - PPO_trt_geq[i]
    }
  }
  
  # calculate p(Y = i | control) for i = 1, ... , c
  for (i in 1:c) {
    if (i == 1) {
      p_y_ctl_est[1] <-  PPO_ctl_leq[1]
      p_y_trt_est[1] <-  PPO_trt_leq[1]
    } else if (i == c) {
      p_y_ctl_est[c] <- 1 - sum(p_y_ctl_est[1:(c-1)])
      p_y_trt_est[c] <- 1 - sum(p_y_trt_est[1:(c-1)])
    } else {
      p_y_ctl_est[i] <-  PPO_ctl_leq[i] - sum(p_y_ctl_est[1:(i-1)])
      p_y_trt_est[i] <-  PPO_trt_leq[i] - sum(p_y_trt_est[1:(i-1)])
    }
  }
  
  # correct for any negative probabilities 
  if (any(p_y_ctl_est < 0)) {
    p_y_ctl_est[which(p_y_ctl_est < 0)] <- 0 
    p_y_ctl_est <- p_y_ctl_est/sum(p_y_ctl_est)
  }
  
  if (any(p_y_trt_est < 0)) {
    p_y_trt_est[which(p_y_trt_est < 0)] <- 0 
    p_y_trt_est <- p_y_trt_est/sum(p_y_trt_est)
  }
  
  p_y_est_overall <- (p_y_ctl_est + p_y_trt_est)/2
  
  
  # get weights from the data 
  p_y_ctl <- unname(table(dat[[trt]], dat[[outcome]])[1,])/sum(dat[[trt]] == -0.5)
  p_y_trt <- unname(table(dat[[trt]], dat[[outcome]])[2,])/sum(dat[[trt]] == 0.5)
  p_y <- unname(table(dat[[outcome]]))/nrow(dat)
  
  
  RD_leq <- PPO_trt_leq - PPO_ctl_leq
  RD_geq <- PPO_ctl_geq - PPO_trt_geq
  
  # put control in the numerator so the OR is greater than one if metformin has a positive effect
  OR_le <- (PPO_trt_leq / (1- PPO_trt_leq)) / (PPO_ctl_leq / (1- PPO_ctl_leq))
  log_OR_le <- log((PPO_trt_leq / (1- PPO_trt_leq)) / (PPO_ctl_leq / (1- PPO_ctl_leq)))
  OR_ge <- (PPO_ctl_geq / (1- PPO_ctl_geq)) / (PPO_trt_geq / (1- PPO_trt_geq))
  # use this for the iprior = 1: (PPO_ctl_geq / (1- PPO_ctl_geq))/(PPO_trt_geq / (1 - PPO_trt_geq))
  
  # calculate RR 
  log_RR_le <- log(PPO_trt_leq/PPO_ctl_leq)
  log_RR_ge <- log(PPO_ctl_geq/PPO_trt_geq)
  
  j <- 1:(c-1)
  k <- 2:c
  # numerator <- sum((p_y_ctl[j] + p_y_ctl[k]) * (PPO_trt_leq[j] - PPO_ctl_leq[j]))
  # denominator <- sum(p_y_ctl[j]) + sum(p_y_ctl[k])
  # weighted_RD <- numerator/denominator
  # 
  # 
  # numerator <- sum((p_y[j] + p_y[k]) * (PPO_trt_leq[j] - PPO_ctl_leq[j]))
  # denominator <- sum(p_y[j]) + sum(p_y[k])
  # weighted_RD_overall <- numerator/denominator
  
  #fully bayesian approach 
  numerator <- sum((p_y_ctl_est[j] + p_y_ctl_est[k]) * (PPO_trt_leq[j] - PPO_ctl_leq[j]))
  denominator <- sum(p_y_ctl_est[j]) + sum(p_y_ctl_est[k])
  weighted_RD_B <- numerator/denominator
  
  # fully bayesian with only cumulative weights
  numerator <- sum((PPO_ctl_leq[j] * (1-PPO_ctl_leq[j])) * (PPO_trt_leq[j] - PPO_ctl_leq[j]))
  denominator <- sum((PPO_ctl_leq[j] * (1-PPO_ctl_leq[j])))
  weighted_RD_cum <- numerator/denominator
  
  # fully bayesian with both cumulative and split weights
  numerator <- sum((p_y_ctl_est[j] + p_y_ctl_est[k]) * (PPO_ctl_leq[j] * (1-PPO_ctl_leq[j])) * (PPO_trt_leq[j] - PPO_ctl_leq[j]))
  denominator <- sum((p_y_ctl_est[j] + p_y_ctl_est[k]) * (PPO_ctl_leq[j] * (1-PPO_ctl_leq[j])))
  weighted_RD_both <- numerator/denominator
  
  # numerator <- sum((p_y_est_overall[j] + p_y_est_overall[k]) * (PPO_trt_leq[j] - PPO_ctl_leq[j]))
  # denominator <- sum(p_y_est_overall[j]) + sum(p_y_est_overall[k])
  # weighted_RD_Boverall <- numerator/denominator
  
  # unweighted RD 
  unweighted_RD <- mean(PPO_trt_leq - PPO_ctl_leq)
  
  # # calculate weighted OR 
  # numerator <- sum((p_y_ctl[j] + p_y_ctl[k])  * log(OR_le[j]))
  # denominator <- sum(p_y_ctl[j]) + sum(p_y_ctl[k])
  # weighted_OR <- exp(numerator/denominator)
  # weighted_log_OR <- numerator/denominator
  # 
  # #calculate weighted OR using overall probability as weights
  # numerator <- sum((p_y[j] + p_y[k]) * log(OR_le[j]))
  # denominator <- sum(p_y[j]) + sum(p_y[k])
  # weighted_OR_overall <- exp(numerator/denominator)
  # weighted_log_OR_overall <- numerator/denominator
  
  #bayesian approach for OR
  #calculate weighted OR 
  numerator <- sum((p_y_ctl_est[j] + p_y_ctl_est[k])  * log(OR_le[j]))
  denominator <- sum(p_y_ctl_est[j]) + sum(p_y_ctl_est[k])
  weighted_OR_B <- exp(numerator/denominator)
  weighted_log_OR_B <- numerator/denominator
  
  #calculate weighted OR using overall probability as weights
  # numerator <- sum((p_y_est_overall[j] + p_y_est_overall[k]) * log(OR_le[j]))
  # denominator <- sum(p_y_est_overall[j]) + sum(p_y_est_overall[k])
  # weighted_OR_Boverall <- exp(numerator/denominator)
  # weighted_log_OR_Boverall <- numerator/denominator
  
  # fully bayesian with only cumulative weights
  numerator <- sum((PPO_ctl_leq[j] * (1-PPO_ctl_leq[j])) * log(OR_le[j]))
  denominator <- sum((PPO_ctl_leq[j] * (1-PPO_ctl_leq[j])))
  weighted_OR_cum <- exp(numerator/denominator)
  weighted_log_OR_cum <- numerator/denominator
  
  # fully bayesian with both cumulative and split weights
  numerator <- sum((p_y_ctl_est[j] + p_y_ctl_est[k]) * (PPO_ctl_leq[j] * (1-PPO_ctl_leq[j])) * log(OR_le[j]))
  denominator <- sum((p_y_ctl_est[j] + p_y_ctl_est[k]) * (PPO_ctl_leq[j] * (1-PPO_ctl_leq[j])))
  weighted_OR_both <- exp(numerator/denominator)
  weighted_log_OR_both <- numerator/denominator
  
  # unweighted OR 
  unweighted_log_OR <- mean(log(OR_le))
  unweighted_OR <- exp(mean(log(OR_le)))
  
  # # calculate weighted RR with control probabilities
  # numerator <- sum((p_y_ctl[j] + p_y_ctl[k]) * log_RR_le[j])
  # denominator <- sum(p_y_ctl[j]) + sum(p_y_ctl[k])
  # weighted_log_RR_control <- numerator/denominator
  # 
  # # calculate weighted RR with overall probabilities
  # numerator <- sum((p_y[j] + p_y[k]) * log_RR_le[j])
  # denominator <- sum(p_y[j]) + sum(p_y[k])
  # weighted_log_RR_overall <- numerator/denominator
  
  # bayesian approach RR 
  # calculate weighted RR with control probabilities 
  numerator <- sum((p_y_ctl_est[j] + p_y_ctl_est[k]) * log_RR_le[j])
  denominator <- sum(p_y_ctl_est[j]) + sum(p_y_ctl_est[k])
  weighted_log_RR_Bcontrol <- numerator/denominator
  
  # fully bayesian with only cumulative weights
  numerator <- sum((PPO_ctl_leq[j] * (1-PPO_ctl_leq[j])) * log_RR_le[j])
  denominator <- sum((PPO_ctl_leq[j] * (1-PPO_ctl_leq[j])))
  weighted_RR_cum <- exp(numerator/denominator)
  weighted_log_RR_cum <- numerator/denominator
  
  # fully bayesian with both cumulative and split weights
  numerator <- sum((p_y_ctl_est[j] + p_y_ctl_est[k]) * (PPO_ctl_leq[j] * (1-PPO_ctl_leq[j])) * log_RR_le[j])
  denominator <- sum((p_y_ctl_est[j] + p_y_ctl_est[k]) * (PPO_ctl_leq[j] * (1-PPO_ctl_leq[j])))
  weighted_RR_both <- exp(numerator/denominator)
  weighted_log_RR_both <- numerator/denominator
  
  # fully bayesian with both cumulative and split weights
  numerator <- sum((p_y_ctl_est[j] + p_y_ctl_est[k]) * (PPO_ctl_leq[j] * (1-PPO_ctl_leq[j])) * log_RR_ge[k])
  denominator <- sum((p_y_ctl_est[j] + p_y_ctl_est[k]) * (PPO_ctl_leq[j] * (1-PPO_ctl_leq[j])))
  weighted_RR_ge_both <- exp(numerator/denominator)
  weighted_log_RR_ge_both <- numerator/denominator
  
  # # calculate weighted RR with overall probabilities 
  # numerator <- sum((p_y_est_overall[j] + p_y_est_overall[k]) * log_RR_le[j]) 
  # denominator <- sum(p_y_est_overall[j]) + sum(p_y_est_overall[k])
  # weighted_log_RR_Boverall <- numerator/denominator
  
  # unweighted RR 
  unweighted_log_RR <- mean(log_RR_le)
  unweighted_RR <- exp(mean(log_RR_le))
  
  # unweighted RR ge
  unweighted_log_RR_ge <- mean(log_RR_ge, na.rm = T)
  unweighted_RR_ge <- exp(mean(log_RR_ge, na.rm = T))
  
  # calculate NB 
  NB <- sum(p_y_trt[j] * PPO_ctl_geq[j+1]) - (sum(p_y_ctl[j] * PPO_trt_geq[j+1]))
  
  weights_split <- (p_y_ctl_est[j] + p_y_ctl_est[k]) / sum((p_y_ctl_est[j] + p_y_ctl_est[k]))
  weights_cum <- (PPO_ctl_leq[j] * (1-PPO_ctl_leq[j])) / sum((PPO_ctl_leq[j] * (1-PPO_ctl_leq[j])))
  weights_both <- (p_y_ctl_est[j] + p_y_ctl_est[k]) * (PPO_ctl_leq[j] * (1-PPO_ctl_leq[j])) / 
    sum((p_y_ctl_est[j] + p_y_ctl_est[k]) * (PPO_ctl_leq[j] * (1-PPO_ctl_leq[j])))
  
  
  return(list(p_y_ctl = p_y_ctl,
              p_y_trt = p_y_trt,
              RD_leq = RD_leq,
              RD_geq = RD_geq,
              OR = OR_le,
              log_OR = log_OR_le,
              log_RR_leq = log_RR_le,
              log_RR_geq = log_RR_ge,
              # weighted_RD_control = weighted_RD, 
              # weighted_RD_overall = weighted_RD_overall,
              weighted_RD_Bcontrol = weighted_RD_B, 
              weighted_RD_cum = weighted_RD_cum,
              weighted_RD_both = weighted_RD_both,
              # weighted_RD_Boverall = weighted_RD_Boverall,
              unweighted_RD = unweighted_RD,
              # weighted_OR_control = weighted_OR,
              # weighted_OR_overall = weighted_OR_overall,
              weighted_OR_Bcontrol = weighted_OR_B,
              weighted_OR_cum = weighted_OR_cum,
              weighted_OR_both = weighted_OR_both,
              # weighted_OR_Boverall = weighted_OR_Boverall,
              unweighted_OR = unweighted_OR,
              # weighted_log_OR_overall = weighted_log_OR_overall,
              # weighted_log_OR_control = weighted_log_OR,
              # weighted_log_OR_Boverall = weighted_log_OR_Boverall,
              weighted_log_OR_cum = weighted_log_OR_cum,
              weighted_log_OR_both = weighted_log_OR_both,
              weighted_log_OR_Bcontrol = weighted_log_OR_B,
              unweighted_log_OR = unweighted_log_OR,
              # weighted_log_RR_control = weighted_log_RR_control, 
              # weighted_log_RR_overall = weighted_log_RR_overall,
              weighted_log_RR_Bcontrol = weighted_log_RR_Bcontrol, 
              weighted_log_RR_cum = weighted_log_RR_cum,
              weighted_log_RR_both = weighted_log_RR_both,
              weighted_log_RR_ge_both = weighted_log_RR_ge_both,
              # weighted_log_RR_Boverall = weighted_log_RR_Boverall,
              unweighted_log_RR = unweighted_log_RR,
              unweighted_RR = unweighted_RR,
              unweighted_log_RR_ge = unweighted_log_RR_ge,
              unweighted_RR_ge = unweighted_RR_ge,
              #WR = WR, 
              NB = NB,
              p_y_trt_est = p_y_trt_est,
              p_y_ctl_est = p_y_ctl_est,
              weights_split = weights_split,
              weights_cum = weights_cum,
              weights_both = weights_both))
  
}


# calculate the credible interval 
# assumes a partial proportional odds model 
# assumes the treatment is coded as -0.5 for control and 0.5 for treatment 
# b is the model 
# NOTE: this only works for the iprior = 0 because I didn't negate the alphas
order_weighted_RD2_CI <- function(b, n_post_draws = 4000, dat, trt, outcome) { 
  # first calculate P(y >= i | ctl) and P(y >= i | trt) for all i
  posterior_draws <- as.matrix(b$draws)
  c <- length(b$ylevels)
  PPO_ctl_geq <- matrix(data = NA, nrow = n_post_draws, ncol = c) 
  PPO_trt_geq <- matrix(data = NA, nrow = n_post_draws, ncol = c) 
  PPO_ctl_leq <- matrix(data = NA, nrow = n_post_draws, ncol = c) 
  PPO_trt_leq <- matrix(data = NA, nrow = n_post_draws, ncol = c) 
  p_y_ctl_est <- matrix(data = NA, nrow = n_post_draws, ncol = c)
  p_y_trt_est <- matrix(data = NA, nrow = n_post_draws, ncol = c)
  
  
  for (i in 2:c) {
    # for the first level there is not a partial term, so only need first term and metformin term 
    if (i == 2) {
      PPO_ctl_geq[,i] <- apply(posterior_draws, 1, function (b) plogis(b[1] - 0.5 * b[c]))
      PPO_trt_geq[,i] <- apply(posterior_draws, 1, function (b) plogis(b[1] + 0.5 * b[c]))
      PPO_ctl_leq[,(i - 1)] <- 1 - PPO_ctl_geq[,i]
      PPO_trt_leq[,(i - 1)] <- 1 - PPO_trt_geq[,i]
    } else{
      phrase <- paste0("y>=", i)
      coefs <- cbind(posterior_draws[,str_detect(colnames(posterior_draws), phrase)], posterior_draws[,c])
      PPO_ctl_geq[,i] <- apply(coefs, 1, 
                               function (b) plogis(b[1] - 0.5 * b[2] - 0.5 * b[3]))
      PPO_trt_geq[,i] <- apply(coefs, 1, 
                               function (b) plogis(b[1] + 0.5 * b[2] + 0.5 * b[3]))
      PPO_ctl_leq[,(i - 1)] <- 1 - PPO_ctl_geq[,i]
      PPO_trt_leq[,(i - 1)] <- 1 - PPO_trt_geq[,i]
    }
  }
  
  # calculate p(Y = i | control) using the model for i = 1, ... , c
  for (i in 1:c) {
    if (i == 1) {
      p_y_ctl_est[,1] <-  PPO_ctl_leq[,1]
      p_y_trt_est[,1] <-  PPO_trt_leq[,1]
    } else if (i == 2){
      p_y_ctl_est[,2] <- PPO_ctl_leq[,i] - PPO_ctl_leq[,1]
      p_y_trt_est[,2] <- PPO_trt_leq[,i] - PPO_trt_leq[,1]
    } else if (i == c) {
      p_y_ctl_est[,c] <- 1 - rowSums(p_y_ctl_est[,1:(c-1)])
      p_y_trt_est[,c] <- 1 - rowSums(p_y_trt_est[,1:(c-1)])
    } else {
      p_y_ctl_est[,i] <-  PPO_ctl_leq[,i] - apply(p_y_ctl_est[,1:(i-1)], 1, sum)
      p_y_trt_est[,i] <-  PPO_trt_leq[,i] - apply(p_y_trt_est[,1:(i-1)], 1, sum)
    }
  }
  
  # correct for any negative probabilities 
  if (any(p_y_ctl_est < 0)) {
    p_y_ctl_est[which(p_y_ctl_est < 0)] <- 0 
    p_y_ctl_est <- p_y_ctl_est/rowSums(p_y_ctl_est)
  }
  
  if (any(p_y_trt_est < 0)) {
    p_y_trt_est[which(p_y_trt_est < 0)] <- 0 
    p_y_trt_est <- p_y_trt_est/rowSums(p_y_trt_est)
  }
  
  
  p_y_est <- (p_y_ctl_est + p_y_trt_est)/2
  
  # get weights from the data 
  p_y_ctl <- unname(table(dat[[trt]], dat[[outcome]])[1,])/sum(dat[[trt]] == -0.5)
  p_y_trt <- unname(table(dat[[trt]], dat[[outcome]])[2,])/sum(dat[[trt]] == 0.5)
  p_y <- as.vector(unname(table(dat[[outcome]])))/nrow(dat)
  
  RD_leq <- PPO_trt_leq - PPO_ctl_leq
  RD_leq_95CI <- apply(RD_leq[,], 2, quantile, probs = c(0.025, 0.975), na.rm = T)
  RD_geq <- PPO_ctl_geq - PPO_trt_geq
  RD_geq_95CI <- apply(RD_geq[,], 2, quantile, probs = c(0.025, 0.975), na.rm = T)
  
  # put control in the numerator so the OR is greater than one if metformin has a positive effect
  OR_le <- (PPO_trt_leq / (1- PPO_trt_leq)) / (PPO_ctl_leq / (1- PPO_ctl_leq))
  OR_le_95CI <- apply(OR_le[, 1:(c-1)], 2, quantile, probs = c(0.025, 0.975), na.rm = T)
  
  log_OR_le <- log((PPO_trt_leq / (1- PPO_trt_leq)) / (PPO_ctl_leq / (1- PPO_ctl_leq)))
  log_OR_le_95CI <- apply(log_OR_le[, 1:(c-1)], 2, quantile, probs = c(0.025, 0.975), na.rm = T)
  
  
  OR_ge <- (PPO_ctl_geq / (1- PPO_ctl_geq)) / (PPO_trt_geq / (1- PPO_trt_geq))
  OR_ge_95CI <- apply(OR_ge[, 1:(c-1)], 2, quantile, probs = c(0.025, 0.975), na.rm = T)
  #for iprior = 1: (PPO_ctl_geq / (1- PPO_ctl_geq))/(PPO_trt_geq / (1 - PPO_trt_geq))
  
  log_RR_le <- log(PPO_trt_leq / PPO_ctl_leq)
  log_RR_le_95CI <- apply(log_RR_le[, 1:(c-1)], 2, quantile, probs = c(0.025, 0.975), na.rm = T)
  
  log_RR_ge <- log(PPO_ctl_geq / PPO_trt_geq)
  log_RR_ge_95CI <- apply(log_RR_ge[, 2:(c)], 2, quantile, probs = c(0.025, 0.975), na.rm = T)
  
  j <- 1:(c-1)
  k <- 2:c
  # numerator <- rowSums((p_y_ctl[j] + p_y_ctl[k]) * (PPO_trt_leq[,j] - PPO_ctl_leq[,j]))
  # denominator <- sum(p_y_ctl[j]) + sum(p_y_ctl[k])
  # weighted_RD <- numerator/denominator
  # weighted_RD_95CI <- quantile(weighted_RD, probs = c(0.025, 0.975))
  # weighted_RD_p <- 1 - sum(weighted_RD >= 0)/n_post_draws
  # weighted_RD_sd <- sd(weighted_RD)
  
  # numerator <- rowSums((p_y[j] + p_y[k]) * (PPO_trt_leq[,j] - PPO_ctl_leq[,j]))
  # denominator <- sum(p_y[j]) + sum(p_y[k])
  # weighted_RD_overall <- numerator/denominator
  # weighted_RD_95CI_overall <- quantile(weighted_RD_overall, probs = c(0.025, 0.975))
  # weighted_RD_p_overall <- 1 - sum(weighted_RD_overall >= 0)/n_post_draws
  # weighted_RD_sd_overall <- sd(weighted_RD_overall)
  
  # fully bayesian weighted RD 
  numerator <- rowSums((p_y_ctl_est[,j] + p_y_ctl_est[,k]) * (PPO_trt_leq[,j] - PPO_ctl_leq[,j]))
  denominator <- rowSums(p_y_ctl_est[,j] + p_y_ctl_est[,k])
  weighted_RD_B <- numerator/denominator
  weighted_RD_B95CI <- quantile(weighted_RD_B, probs = c(0.025, 0.975))
  weighted_RD_Bp <- 1 - sum(weighted_RD_B >= 0)/n_post_draws
  weighted_RD_Bsd <- sd(weighted_RD_B)
  
  numerator <- rowSums((PPO_ctl_leq[,j] * (1-PPO_ctl_leq[,j])) * (PPO_trt_leq[,j] - PPO_ctl_leq[,j]))
  denominator <- rowSums((PPO_ctl_leq[,j] * (1-PPO_ctl_leq[,j])))
  weighted_RD_cum <- numerator/denominator
  weighted_RD_cum95CI <- quantile(weighted_RD_cum, probs = c(0.025, 0.975))
  weighted_RD_cump <- 1 - sum(weighted_RD_cum >= 0)/n_post_draws
  weighted_RD_cumsd <- sd(weighted_RD_cum)
  
  numerator <- rowSums((p_y_ctl_est[,j] + p_y_ctl_est[,k]) * (PPO_ctl_leq[,j] * (1-PPO_ctl_leq[,j])) * (PPO_trt_leq[,j] - PPO_ctl_leq[,j]))
  denominator <- rowSums((p_y_ctl_est[,j] + p_y_ctl_est[,k]) * (PPO_ctl_leq[,j] * (1-PPO_ctl_leq[,j])))
  weighted_RD_both <- numerator/denominator
  weighted_RD_both95CI <- quantile(weighted_RD_both, probs = c(0.025, 0.975))
  weighted_RD_bothp <- 1 - sum(weighted_RD_both >= 0)/n_post_draws
  weighted_RD_bothsd <- sd(weighted_RD_both)
  
  # numerator <- rowSums((p_y_est[j] + p_y_est[k]) * (PPO_trt_leq[,j] - PPO_ctl_leq[,j]))
  # denominator <- sum(p_y_est[j]) + sum(p_y_est[k])
  # weighted_RD_Boverall <- numerator/denominator
  # weighted_RD_B95CI_overall <- quantile(weighted_RD_overall, probs = c(0.025, 0.975))
  # weighted_RD_Bp_overall <- 1 - sum(weighted_RD_overall >= 0)/n_post_draws
  # weighted_RD_Bsd_overall <- sd(weighted_RD_overall)
  
  #unweighted RD
  unweighted_RD <- rowMeans(RD_leq[,1:(c-1)])
  unweighted_RD_95CI <- quantile(unweighted_RD, probs = c(0.025, 0.975))
  unweighted_RD_p <- 1 - sum(unweighted_RD >= 0)/n_post_draws
  unweighted_RD_sd <- sd(unweighted_RD)
  
  
  #calculate weighted OR 
  # numerator <- rowSums((p_y_ctl[j] + p_y_ctl[k]) * log(OR_le[,j])) 
  # denominator <- sum(p_y_ctl[j]) + sum(p_y_ctl[k])
  # weighted_OR <- exp(numerator/denominator)
  # weighted_OR_95CI <- quantile(weighted_OR, probs = c(0.025, 0.975))
  # weighted_OR_p <- 1 - sum(weighted_OR >= 1)/n_post_draws
  # weighted_OR_sd <- sd(weighted_OR)
  
  # weighted_log_OR <- numerator/denominator
  # weighted_log_OR_95CI <- quantile(weighted_log_OR, probs = c(0.025, 0.975))
  # weighted_log_OR_p <- 1 - sum(weighted_log_OR >= 1)/n_post_draws
  # weighted_log_OR_sd <- sd(weighted_log_OR)
  
  # #calculate weighted OR with overall probability as weights 
  # numerator <- rowSums((p_y[j] + p_y[k]) * log(OR_le[,j])) 
  # denominator <- sum(p_y[j]) + sum(p_y[k])
  # weighted_OR_overall <- exp(numerator/denominator)
  # weighted_OR_95CI_overall <- quantile(weighted_OR_overall, probs = c(0.025, 0.975))
  # weighted_OR_p_overall <- 1 - sum(weighted_OR_overall >= 1)/n_post_draws
  # weighted_OR_sd_overall <- sd(weighted_OR_overall)
  # 
  # weighted_log_OR_overall <- numerator/denominator
  # weighted_log_OR_95CI_overall <- quantile(weighted_log_OR_overall, probs = c(0.025, 0.975))
  # weighted_log_OR_p_overall <- 1 - sum(weighted_log_OR_overall >= 1)/n_post_draws
  # weighted_log_OR_sd_overall <- sd(weighted_log_OR_overall)
  
  # fully bayesian estimates
  #calculate weighted OR 
  numerator <- rowSums((p_y_ctl_est[,j] + p_y_ctl_est[,k]) * log(OR_le[,j])) 
  denominator <- rowSums(p_y_ctl_est[,j] + p_y_ctl_est[,k])
  weighted_OR_B <- exp(numerator/denominator)
  weighted_OR_B95CI <- quantile(weighted_OR_B, probs = c(0.025, 0.975))
  weighted_OR_Bp <- 1 - sum(weighted_OR_B >= 1)/n_post_draws
  weighted_OR_Bsd <- sd(weighted_OR_B)
  
  weighted_log_OR_B <- numerator/denominator
  weighted_log_OR_B95CI <- quantile(weighted_log_OR_B, probs = c(0.025, 0.975))
  weighted_log_OR_Bp <- 1 - sum(weighted_log_OR_B >= 1)/n_post_draws
  weighted_log_OR_Bsd <- sd(weighted_log_OR_B)
  
  # cumulative weights only
  numerator <- rowSums((PPO_ctl_leq[,j] * (1-PPO_ctl_leq[,j])) * log(OR_le[,j]))
  denominator <- rowSums((PPO_ctl_leq[,j] * (1-PPO_ctl_leq[,j])))
  weighted_OR_cum <- exp(numerator/denominator)
  weighted_OR_cum95CI <- quantile(weighted_OR_cum, probs = c(0.025, 0.975))
  weighted_OR_cump <- 1 - sum(weighted_OR_cum >= 0)/n_post_draws
  weighted_OR_cumsd <- sd(weighted_OR_cum)
  
  weighted_log_OR_cum <- numerator/denominator
  weighted_log_OR_cum95CI <- quantile(weighted_log_OR_cum, probs = c(0.025, 0.975))
  weighted_log_OR_cump <- 1 - sum(weighted_log_OR_cum >= 0)/n_post_draws
  weighted_log_OR_cumsd <- sd(weighted_log_OR_cum)
  
  # both weights
  numerator <- rowSums((p_y_ctl_est[,j] + p_y_ctl_est[,k]) * (PPO_ctl_leq[,j] * (1-PPO_ctl_leq[,j])) * log(OR_le[,j]))
  denominator <- rowSums((p_y_ctl_est[,j] + p_y_ctl_est[,k]) * (PPO_ctl_leq[,j] * (1-PPO_ctl_leq[,j])))
  weighted_OR_both <- exp(numerator/denominator)
  weighted_OR_both95CI <- quantile(weighted_OR_both, probs = c(0.025, 0.975))
  weighted_OR_bothp <- 1 - sum(weighted_OR_both >= 0)/n_post_draws
  weighted_OR_bothsd <- sd(weighted_OR_both)
  
  weighted_log_OR_both <- numerator/denominator
  weighted_log_OR_both95CI <- quantile(weighted_log_OR_both, probs = c(0.025, 0.975))
  weighted_log_OR_bothp <- 1 - sum(weighted_log_OR_both >= 0)/n_post_draws
  weighted_log_OR_bothsd <- sd(weighted_log_OR_both)
  
  #calculate weighted OR with overall probability as weights 
  # numerator <- rowSums((p_y_est[j] + p_y_est[k]) * log(OR_le[,j])) 
  # denominator <- sum(p_y_est[j]) + sum(p_y_est[k])
  # weighted_OR_Boverall <- exp(numerator/denominator)
  # weighted_OR_B95CI_overall <- quantile(weighted_OR_overall, probs = c(0.025, 0.975))
  # weighted_OR_Bp_overall <- 1 - sum(weighted_OR_overall >= 1)/n_post_draws
  # weighted_OR_Bsd_overall <- sd(weighted_OR_overall)
  # 
  # weighted_log_OR_Boverall <- numerator/denominator
  # weighted_log_OR_B95CI_overall <- quantile(weighted_log_OR_overall, probs = c(0.025, 0.975))
  # weighted_log_OR_Bp_overall <- 1 - sum(weighted_log_OR_overall >= 1)/n_post_draws
  # weighted_log_OR_Bsd_overall <- sd(weighted_log_OR_overall)
  
  #unweighted OR 
  unweighted_log_OR <- rowMeans(log_OR_le[,1:(c-1)])
  unweighted_log_OR_95CI <- quantile(unweighted_log_OR, probs = c(0.025, 0.975))
  unweighted_log_OR_p <- 1 - sum(unweighted_log_OR >= 0)/n_post_draws
  unweighted_log_OR_sd <- sd(unweighted_log_OR)
  
  unweighted_OR_95CI <- quantile(exp(unweighted_log_OR), probs = c(0.025, 0.975))
  unweighted_OR_p <- 1 - sum(exp(unweighted_log_OR) >= 0)/n_post_draws
  unweighted_OR_sd <- sd(exp(unweighted_log_OR))
  
  # # calculate weighted RR with control probabilities 
  # numerator <- rowSums((p_y_ctl[j] + p_y_ctl[k]) * log_RR_le[,j])
  # denominator <- sum(p_y_ctl[j]) + sum(p_y_ctl[k])
  # weighted_log_RR_control <- numerator/denominator
  # weighted_log_RR_95CI_control <- quantile(weighted_log_RR_control, probs = c(0.025, 0.975))
  # weighted_log_RR_p_control <- 1 - sum(weighted_log_RR_control >= 0)/n_post_draws
  # weighted_log_RR_sd_control <- sd(weighted_log_RR_control)
  
  # # calculate weighted RR with overall probabilities 
  # numerator <- rowSums((p_y[j] + p_y[k]) * log_RR_le[,j])
  # denominator <- sum(p_y[j]) + sum(p_y[k])
  # weighted_log_RR_overall <- numerator/denominator
  # weighted_log_RR_95CI_overall <- quantile(weighted_log_RR_overall, probs = c(0.025, 0.975))
  # weighted_log_RR_p_overall <- 1 - sum(weighted_log_RR_overall >= 0)/n_post_draws
  # weighted_log_RR_sd_overall <- sd(weighted_log_RR_overall)
  
  # fully bayesian weighted estimate 
  # calculate weighted RR with control probabilities 
  numerator <- rowSums((p_y_ctl_est[,j] + p_y_ctl_est[,k]) * log_RR_le[,j])
  denominator <- rowSums(p_y_ctl_est[,j] + p_y_ctl_est[,k])
  weighted_log_RR_Bcontrol <- numerator/denominator
  weighted_log_RR_B95CI_control <- quantile(weighted_log_RR_Bcontrol, probs = c(0.025, 0.975))
  weighted_log_RR_Bp_control <- 1 - sum(weighted_log_RR_Bcontrol >= 0)/n_post_draws
  weighted_log_RR_Bsd_control <- sd(weighted_log_RR_Bcontrol)
  
  # cumulative weights only
  numerator <- rowSums((PPO_ctl_leq[,j] * (1-PPO_ctl_leq[,j])) * log_RR_le[,j])
  denominator <- rowSums((PPO_ctl_leq[,j] * (1-PPO_ctl_leq[,j])))
  weighted_log_RR_cum <- numerator/denominator
  weighted_log_RR_cum95CI <- quantile(weighted_log_RR_cum, probs = c(0.025, 0.975))
  weighted_log_RR_cump <- 1 - sum(weighted_log_RR_cum >= 0)/n_post_draws
  weighted_log_RR_cumsd <- sd(weighted_log_RR_cum)
  
  # both weights
  numerator <- rowSums((p_y_ctl_est[,j] + p_y_ctl_est[,k]) * (PPO_ctl_leq[,j] * (1-PPO_ctl_leq[,j])) * log_RR_le[,j])
  denominator <- rowSums((p_y_ctl_est[,j] + p_y_ctl_est[,k]) * (PPO_ctl_leq[,j] * (1-PPO_ctl_leq[,j])))
  weighted_log_RR_both <- numerator/denominator
  weighted_log_RR_both95CI <- quantile(weighted_log_RR_both, probs = c(0.025, 0.975))
  weighted_log_RR_bothp <- 1 - sum(weighted_log_RR_both >= 0)/n_post_draws
  weighted_log_RR_bothsd <- sd(weighted_log_RR_both)
  
  # both weights RR greater than or equal to 
  numerator <- rowSums((p_y_ctl_est[,j] + p_y_ctl_est[,k]) * (PPO_ctl_leq[,j] * (1-PPO_ctl_leq[,j])) * log_RR_ge[,k])
  denominator <- rowSums((p_y_ctl_est[,j] + p_y_ctl_est[,k]) * (PPO_ctl_leq[,j] * (1-PPO_ctl_leq[,j])))
  weighted_log_RR_ge_both <- numerator/denominator
  weighted_log_RR_ge_both95CI <- quantile(weighted_log_RR_ge_both, probs = c(0.025, 0.975))
  weighted_log_RR_ge_bothp <- 1 - sum(weighted_log_RR_ge_both >= 0)/n_post_draws
  weighted_log_RR_ge_bothsd <- sd(weighted_log_RR_ge_both)
  
  # # calculate weighted RR with overall probabilities 
  # numerator <- rowSums((p_y_est[j] + p_y_est[k]) * log_RR_le[,j])
  # denominator <- sum(p_y_est[j]) + sum(p_y_est[k])
  # weighted_log_RR_Boverall <- numerator/denominator
  # weighted_log_RR_B95CI_overall <- quantile(weighted_log_RR_overall, probs = c(0.025, 0.975))
  # weighted_log_RR_Bp_overall <- 1 - sum(weighted_log_RR_overall >= 0)/n_post_draws
  # weighted_log_RR_Bsd_overall <- sd(weighted_log_RR_overall)
  
  #unweighted RR
  unweighted_log_RR <- rowMeans(log_RR_le[, 1:(c-1)])
  unweighted_log_RR_95CI <- quantile(unweighted_log_RR, probs = c(0.025, 0.975))
  unweighted_log_RR_p <- 1 - sum(unweighted_log_RR >= 0)/n_post_draws
  unweighted_log_RR_sd <- sd(unweighted_log_RR)
  
  unweighted_RR_95CI <- quantile(exp(unweighted_log_RR), probs = c(0.025, 0.975))
  unweighted_RR_p <- 1 - sum(exp(unweighted_log_RR) >= 0)/n_post_draws
  unweighted_RR_sd <- sd(exp(unweighted_log_RR))
  
  #unweighted RR greater than or equal to
  unweighted_log_RR_ge <- rowMeans(log_RR_ge[, 2:(c)])
  unweighted_log_RR_ge_95CI <- quantile(unweighted_log_RR_ge, probs = c(0.025, 0.975))
  unweighted_log_RR_ge_p <- 1 - sum(unweighted_log_RR_ge >= 0)/n_post_draws
  unweighted_log_RR_ge_sd <- sd(unweighted_log_RR_ge)
  
  unweighted_RR_ge_95CI <- quantile(exp(unweighted_log_RR_ge), probs = c(0.025, 0.975))
  unweighted_RR_ge_p <- 1 - sum(exp(unweighted_log_RR_ge) >= 0)/n_post_draws
  unweighted_RR_ge_sd <- sd(exp(unweighted_log_RR_ge))
  
  # calculate NB 
  NB <- rowSums(p_y_trt[j] * PPO_ctl_geq[,(j+1)]) - (rowSums(p_y_ctl[j] * PPO_trt_geq[,(j+1)]))
  NB_95CI <- quantile(NB, probs = c(0.025, 0.975))
  NB_p <- 1 - sum(NB >= 0)/n_post_draws
  NB_sd <- sd(NB)
  
  #weights 
  weights_split <- (p_y_ctl_est[,j] + p_y_ctl_est[,k]) / rowSums(p_y_ctl_est[,j] + p_y_ctl_est[,k])
  weights_split <- apply(weights_split, 2, quantile, probs = c(0.025, 0.975))
  
  weights_cum <- (PPO_ctl_leq[,j] * (1-PPO_ctl_leq[,j])) / rowSums((PPO_ctl_leq[,j] * (1-PPO_ctl_leq[,j])))
  weights_cum <- apply(weights_cum, 2, quantile, probs = c(0.025, 0.975))
  
  weights_both <- ((p_y_ctl_est[,j] + p_y_ctl_est[,k]) * (PPO_ctl_leq[,j] * (1-PPO_ctl_leq[,j]))) /
    rowSums(((p_y_ctl_est[,j] + p_y_ctl_est[,k]) * (PPO_ctl_leq[,j] * (1-PPO_ctl_leq[,j]))))
  weights_both <- apply(weights_both, 2, quantile, probs = c(0.025, 0.975))
  
  return(list(RD_leq = RD_leq_95CI,
              RD_geq = RD_geq_95CI,
              OR = OR_le_95CI,
              log_OR = log_OR_le_95CI,
              log_RR_leq = log_RR_le_95CI,
              log_RR_geq = log_RR_ge_95CI,
              # weighted_RD_control_CI = weighted_RD_95CI,
              # weighted_RD_p_control = weighted_RD_p,
              # weighted_RD_sd = weighted_RD_sd, 
              # weighted_RD_sd_overall = weighted_RD_sd_overall,
              # weighted_RD_overall_CI = weighted_RD_95CI_overall,
              # weighted_RD_p_overall = weighted_RD_p_overall,
              weighted_RD_Bcontrol_CI = weighted_RD_B95CI,
              weighted_RD_Bp_control = weighted_RD_Bp,
              weighted_RD_Bsd = weighted_RD_Bsd, 
              weighted_RD_cum_CI = weighted_RD_cum95CI,
              weighted_RD_cum_p = weighted_RD_cump,
              weighted_RD_cum_sd = weighted_RD_cumsd,
              weighted_RD_both_CI = weighted_RD_both95CI,
              weighted_RD_both_p = weighted_RD_bothp,
              weighted_RD_both_sd = weighted_RD_bothsd, 
              # weighted_RD_Bsd_overall = weighted_RD_Bsd_overall,
              # weighted_RD_Boverall_CI = weighted_RD_B95CI_overall,
              # weighted_RD_Bp_overall = weighted_RD_Bp_overall,
              # weighted_OR_CI = weighted_OR_95CI, 
              # weighted_OR_p = weighted_OR_p,
              # weighted_OR_sd = weighted_OR_sd,
              # weighted_OR_overall_CI = weighted_OR_95CI_overall,
              # weighted_OR_p_overall = weighted_OR_p_overall,
              # weighted_OR_sd_overall = weighted_OR_sd_overall,
              # weighted_log_OR_control_CI = weighted_log_OR_95CI,
              # weighted_log_OR_p_control = weighted_log_OR_p,
              # weighted_log_OR_sd_control = weighted_log_OR_sd,
              # weighted_log_OR_overall_CI = weighted_log_OR_95CI_overall,
              # weighted_log_OR_p_overall = weighted_log_OR_p_overall,
              # weighted_log_OR_sd_overall = weighted_log_OR_sd_overall,
              weighted_log_OR_Bcontrol_CI = weighted_log_OR_B95CI,
              weighted_log_OR_Bp_control = weighted_log_OR_Bp,
              weighted_log_OR_Bsd_control = weighted_log_OR_Bsd,
              weighted_log_OR_cum_CI = weighted_log_OR_cum95CI,
              weighted_log_OR_cum_p = weighted_log_OR_cump,
              weighted_log_OR_cum_sd = weighted_log_OR_cumsd,
              weighted_log_OR_both_CI = weighted_log_OR_both95CI,
              weighted_log_OR_both_p = weighted_log_OR_bothp,
              weighted_log_OR_both_sd = weighted_log_OR_bothsd,
              # weighted_log_OR_Boverall_CI = weighted_log_OR_B95CI_overall,
              # weighted_log_OR_Bp_overall = weighted_log_OR_Bp_overall,
              # weighted_log_OR_Bsd_overall = weighted_log_OR_Bsd_overall,
              # weighted_log_RR_95CI_control = weighted_log_RR_95CI_control,
              # weighted_log_RR_p_control = weighted_log_RR_p_control, 
              # weighted_log_RR_sd_control = weighted_log_RR_sd_control,
              # weighted_log_RR_95CI_overall = weighted_log_RR_95CI_overall,
              # weighted_log_RR_p_overall = weighted_log_RR_p_overall, 
              # weighted_log_RR_sd_overall = weighted_log_RR_sd_overall,
              weighted_log_RR_B95CI_control = weighted_log_RR_B95CI_control,
              weighted_log_RR_Bp_control = weighted_log_RR_Bp_control, 
              weighted_log_RR_Bsd_control = weighted_log_RR_Bsd_control,
              weighted_log_RR_cum_CI = weighted_log_RR_cum95CI,
              weighted_log_RR_cum_p = weighted_log_RR_cump,
              weighted_log_RR_cum_sd = weighted_log_RR_cumsd,
              weighted_log_RR_both_CI = weighted_log_RR_both95CI,
              weighted_log_RR_both_p = weighted_log_RR_bothp,
              weighted_log_RR_both_sd = weighted_log_RR_bothsd,
              weighted_log_RR_ge_both_CI = weighted_log_RR_ge_both95CI,
              weighted_log_RR_ge_both_p = weighted_log_RR_ge_bothp,
              weighted_log_RR_ge_both_sd = weighted_log_RR_ge_bothsd,
              # weighted_log_RR_B95CI_overall = weighted_log_RR_B95CI_overall,
              # weighted_log_RR_Bp_overall = weighted_log_RR_Bp_overall, 
              # weighted_log_RR_Bsd_overall = weighted_log_RR_Bsd_overall,
              
              unweighted_RD_CI = unweighted_RD_95CI,
              unweighted_RD_p = unweighted_RD_p,
              unweighted_RD_sd = unweighted_RD_sd, 
              unweighted_OR_CI = unweighted_OR_95CI,
              unweighted_OR_p = unweighted_OR_p,
              unweighted_OR_sd = unweighted_OR_sd,  
              unweighted_log_OR_CI = unweighted_log_OR_95CI,
              unweighted_log_OR_p = unweighted_log_OR_p,
              unweighted_log_OR_sd = unweighted_log_OR_sd,
              unweighted_RR_CI = unweighted_RR_95CI,
              unweighted_RR_p = unweighted_RR_p,
              unweighted_RR_sd = unweighted_RR_sd,  
              unweighted_log_RR_CI = unweighted_log_RR_95CI,
              unweighted_log_RR_p = unweighted_log_RR_p,
              unweighted_log_RR_sd = unweighted_log_RR_sd,
              unweighted_RR_ge_CI = unweighted_RR_ge_95CI,
              unweighted_RR_ge_p = unweighted_RR_ge_p,
              unweighted_RR_ge_sd = unweighted_RR_ge_sd,  
              unweighted_log_RR_ge_CI = unweighted_log_RR_ge_95CI,
              unweighted_log_RR_ge_p = unweighted_log_RR_ge_p,
              unweighted_log_RR_ge_sd = unweighted_log_RR_ge_sd,
              
              #WR = WR_95CI,
              NB_CI = NB_95CI,
              NB_p = NB_p,
              NB_sd = NB_sd,
              
              weights_split = weights_split,
              weights_cum = weights_cum, 
              weights_both = weights_both))
  
}



# Adjusted Model

# modify adjusted PPO function to use Bayesian Bootstrap to get average 

order_weighted_RD2_PPO_adj_bb <- function(b = b, dat, iprior,
                                          outcome, trt, covs = 1) { 
  c <- length(b$ylevels)
  
  predicted_trt_geq <- matrix(data = NA, nrow = nrow(dat), ncol = c-1) 
  predicted_ctl_geq <- matrix(data = NA, nrow = nrow(dat), ncol = c-1) 
  predicted_trt_leq <- matrix(data = NA, nrow = nrow(dat), ncol = c-1) 
  predicted_ctl_leq <- matrix(data = NA, nrow = nrow(dat), ncol = c-1) 
  p_y_ctl <- matrix(data = NA, nrow = nrow(dat), ncol = c) 
  p_y_trt <- matrix(data = NA, nrow = nrow(dat), ncol = c) 
  p_y_ctl_est <- matrix(data = NA, nrow = nrow(dat), ncol = c) 
  p_y_trt_est <- matrix(data = NA, nrow = nrow(dat), ncol = c) 
  
  if (iprior == 0){
    coefs <- coef(b)
  } else {
    coefs <- c(-coef(b)[1:3], coef(b)[4:9])
  }
  
  # covariates <- unlist(str_split(deparse(b_ppo_model_adj$call$formula[[3]]), " \\+ "))
  # covariates <- covariates[-which(covariates == trt)]
  # num_covs <- length(covariates)
  
  covs <- b$Design$colnames
  num_covs <- length(covs)
  covariate_mat <- as.data.frame(matrix(data = NA, nrow = nrow(dat), ncol = num_covs - 1))
  
  # assume the first covariate is the treatment arm 
  for(i in 2:num_covs) {
    if (is.numeric(dat[[b$Design$name[i]]]) == T) {
      covariate_mat[,i-1] <- coefs[[covs[i]]] * (dat[[b$Design$name[i]]])
    } else {
      # only works for numeric and binary
      covariate_mat[,i-1] <- coefs[[covs[i]]] * (dat[[b$Design$name[i]]] == b$Design$parms[[b$Design$name[i]]][2])
    }
  }
  
  for (i in 2:c) {
    if (i == 2) {
      predicted_trt_geq[,i-1] <- plogis(coefs[i-1] + coefs[[trt]] * 0.5 + rowSums(covariate_mat))
      predicted_ctl_geq[,i-1] <- plogis(coefs[i-1] + coefs[[trt]] * -0.5 + rowSums(covariate_mat))
    } else {
      phrase <- paste0("y>=", i)
      coefs_new <- coefs[str_detect(names(coefs), phrase)]
      predicted_trt_geq[,i-1] <- plogis(coefs_new[1] + coefs_new[2] * 0.5 + coefs[[trt]] * 0.5 + rowSums(covariate_mat))
      predicted_ctl_geq[,i-1] <- plogis(coefs_new[1] + coefs_new[2] * -0.5 + coefs[[trt]] * -0.5 + rowSums(covariate_mat))
    }
  }
  
  
  predicted_trt_leq <-  1 - predicted_trt_geq
  predicted_ctl_leq <-  1 - predicted_ctl_geq
  
  # calculate p(Y = i | control) for i = 1, ... , c
  for (i in 1:c) {
    if (i == 1) {
      p_y_ctl_est[,1] <-  predicted_ctl_leq[,1]
      p_y_trt_est[,1] <-  predicted_trt_leq[,1]
    } else if (i == 2) {
      p_y_ctl_est[,i] <-  predicted_ctl_leq[,i] - p_y_ctl_est[,1]
      p_y_trt_est[,i] <-  predicted_trt_leq[,i] - p_y_trt_est[,1]
    } else if (i == c) {
      p_y_ctl_est[,c] <- 1 - rowSums(p_y_ctl_est[,1:(c-1)])
      p_y_trt_est[,c] <- 1 - rowSums(p_y_trt_est[,1:(c-1)])
    } else {
      p_y_ctl_est[,i] <-  predicted_ctl_leq[,i] - rowSums(p_y_ctl_est[,1:(i-1)])
      p_y_trt_est[,i] <-  predicted_trt_leq[,i] - rowSums(p_y_trt_est[,1:(i-1)])
    }
  }
  
  # correct for any negative probabilities 
  if (any(p_y_ctl_est < 0)) {
    p_y_ctl_est[which(p_y_ctl_est < 0)] <- 0 
    p_y_ctl_est <- p_y_ctl_est/rowSums(p_y_ctl_est)
  }
  
  if (any(p_y_trt_est < 0)) {
    p_y_trt_est[which(p_y_trt_est < 0)] <- 0 
    p_y_trt_est <- p_y_trt_est/rowSums(p_y_trt_est)
  }
  
  
  p_y_est <- (p_y_ctl_est + p_y_trt_est)/2
  
  # get weights from the data 
  p_y_ctl <- unname(table(dat[[trt]], dat[[outcome]])[1,])/sum(dat[[trt]] == -0.5)
  p_y_trt <- unname(table(dat[[trt]], dat[[outcome]])[2,])/sum(dat[[trt]] == 0.5)
  p_y <- unname(table(dat[[outcome]]))/nrow(dat)
  
  # use Bayesian bootstrap to marginalize the probabilities
  weights <- t(rdirichlet(1, rep(1, nrow(dat))))
  p_y_ctl_est_bb <- apply(p_y_ctl_est, 2, weighted.mean, w = weights)
  p_y_trt_est_bb <- apply(p_y_trt_est, 2, weighted.mean, w = weights)
  p_y_est_bb <- (p_y_ctl_est_bb + p_y_trt_est_bb) /2
  
  p_y_ctl_leq_average <- apply(predicted_ctl_leq, 2, weighted.mean, w = weights)
  p_y_trt_leq_average <- apply(predicted_trt_leq, 2, weighted.mean, w = weights)
  p_y_ctl_geq_average <- apply(predicted_ctl_geq, 2, weighted.mean, w = weights)
  p_y_trt_geq_average <- apply(predicted_trt_geq, 2, weighted.mean, w = weights)
  
  #p_y_est_average <- apply(p_y_est, 2, weighted.mean, w = weights)
  #p_y_ctl_est_average <- apply(p_y_ctl_est, 2, weighted.mean, w = weights)
  
  
  RD_leq <- p_y_trt_leq_average - p_y_ctl_leq_average
  RD_geq <- p_y_ctl_geq_average - p_y_trt_geq_average
  
  
  OR <- (p_y_trt_leq_average / (1- p_y_trt_leq_average)) / (p_y_ctl_leq_average / (1- p_y_ctl_leq_average))
  
  RR_leq <- p_y_trt_leq_average / p_y_ctl_leq_average
  RR_geq <- p_y_ctl_geq_average / p_y_trt_geq_average
  
  j <- 1:(c-1)
  k <- 2:c
  # numerator <- sum((p_y_ctl[j] + p_y_ctl[k]) * (p_y_trt_leq_average[j] - p_y_ctl_leq_average[j]))
  # denominator <- sum(p_y_ctl[j]) + sum(p_y_ctl[k])
  # weighted_RD <- numerator/denominator
  # 
  # #calculated weighted RD with overall probability as weights
  # numerator <- sum((p_y[j] + p_y[k]) * (p_y_trt_leq_average[j] - p_y_ctl_leq_average[j]))
  # denominator <- sum(p_y[j]) + sum(p_y[k])
  # weighted_RD_overall <- numerator/denominator
  
  # fully bayesian approach 
  # control weighted
  numerator <- sum((p_y_ctl_est_bb[j] + p_y_ctl_est_bb[k]) * (p_y_trt_leq_average[j] - p_y_ctl_leq_average[j]))
  denominator <- sum(p_y_ctl_est_bb[j]) + sum(p_y_ctl_est_bb[k])
  weighted_RD_B <- numerator/denominator
  
  numerator <- sum((p_y_ctl_leq_average[j] * (1 - p_y_ctl_leq_average[j])) * (p_y_trt_leq_average[j] - p_y_ctl_leq_average[j]))
  denominator <- sum((p_y_ctl_leq_average[j] * (1 - p_y_ctl_leq_average[j])))
  weighted_RD_cum <- numerator/denominator
  
  numerator <- sum((p_y_ctl_est_bb[j] + p_y_ctl_est_bb[k]) * (p_y_ctl_leq_average[j] * (1 - p_y_ctl_leq_average[j])) * (p_y_trt_leq_average[j] - p_y_ctl_leq_average[j]))
  denominator <- sum((p_y_ctl_est_bb[j] + p_y_ctl_est_bb[k]) * (p_y_ctl_leq_average[j] * (1 - p_y_ctl_leq_average[j])))
  weighted_RD_both <- numerator/denominator
  
  # # calculated weighted RD with overall probability as weights
  # numerator <- sum((p_y_est_bb[j] + p_y_est_bb[k]) * (p_y_trt_leq_average[j] - p_y_ctl_leq_average[j]))
  # denominator <- sum(p_y_est_bb[j]) + sum(p_y_est_bb[k])
  # weighted_RD_Boverall <- numerator/denominator
  
  # unweighted RD 
  unweighted_RD <- mean(p_y_trt_leq_average - p_y_ctl_leq_average)
  
  # #calculate weighted OR with control probability as weights
  # numerator <- sum((p_y_ctl[j] + p_y_ctl[k]) * log(OR[j]))
  # denominator <- sum(p_y_ctl[j]) + sum(p_y_ctl[k])
  # weighted_OR <- exp(numerator/denominator)
  # log_weighted_OR <- numerator/denominator
  # 
  # #calculate weighted OR with overall probability as weights
  # numerator <- sum((p_y[j] + p_y[k]) * log(OR[j]))
  # denominator <- sum(p_y[j]) + sum(p_y[k])
  # log_weighted_OR_overall <- numerator/denominator
  # weighted_OR_overall <- exp(numerator/denominator)
  
  # bayesian approach
  #calculate weighted OR with control probability as weights
  numerator <- sum((p_y_ctl_est_bb[j] + p_y_ctl_est_bb[k]) * log(OR[j]))
  denominator <- sum(p_y_ctl_est_bb[j]) + sum(p_y_ctl_est_bb[k])
  weighted_OR_B <- exp(numerator/denominator)
  log_weighted_OR_B <- numerator/denominator
  
  numerator <- sum((p_y_ctl_leq_average[j] * (1 - p_y_ctl_leq_average[j])) * log(OR[j]))
  denominator <- sum((p_y_ctl_leq_average[j] * (1 - p_y_ctl_leq_average[j])))
  weighted_OR_cum <- exp(numerator/denominator)
  log_weighted_OR_cum <- numerator/denominator
  
  numerator <- sum((p_y_ctl_est_bb[j] + p_y_ctl_est_bb[k]) * (p_y_ctl_leq_average[j] * (1 - p_y_ctl_leq_average[j])) * log(OR[j]))
  denominator <- sum((p_y_ctl_est_bb[j] + p_y_ctl_est_bb[k]) * (p_y_ctl_leq_average[j] * (1 - p_y_ctl_leq_average[j])))
  weighted_OR_both <- exp(numerator/denominator)
  log_weighted_OR_both <- numerator/denominator
  
  # #calculate weighted OR with overall probability as weights
  # numerator <- sum((p_y_est_bb[j] + p_y_est_bb[k]) * log(OR[j]))
  # denominator <- sum(p_y_est_bb[j]) + sum(p_y_est_bb[k])
  # log_weighted_OR_Boverall <- numerator/denominator
  # weighted_OR_Boverall <- exp(numerator/denominator)
  
  # unweighted OR 
  unweighted_log_OR <- mean(log(OR))
  unweighted_OR <- exp(mean(log(OR)))
  
  # #calculate weighted RR with control probability as weights
  # numerator <- sum((p_y_ctl[j] + p_y_ctl[k]) * log(RR_leq[j]))
  # denominator <- sum(p_y_ctl[j]) + sum(p_y_ctl[k])
  # weighted_RR <- exp(numerator/denominator)
  # log_weighted_RR <- numerator/denominator
  # 
  # #calculate weighted RR with overall probability as weights
  # numerator <- sum((p_y[j] + p_y[k]) * log(RR_leq[j]))
  # denominator <- sum(p_y[j]) + sum(p_y[k])
  # log_weighted_RR_overall <- numerator/denominator
  # weighted_RR_overall <- exp(numerator/denominator)
  
  # bayesian approach
  #calculate weighted RR with control probability as weights
  numerator <- sum((p_y_ctl_est_bb[j] + p_y_ctl_est_bb[k]) * log(RR_leq[j]))
  denominator <- sum(p_y_ctl_est_bb[j]) + sum(p_y_ctl_est_bb[k])
  weighted_RR_B <- exp(numerator/denominator)
  log_weighted_RR_B <- numerator/denominator
  
  # cumulative weights 
  numerator <- sum((p_y_ctl_leq_average[j] * (1 - p_y_ctl_leq_average[j])) * log(RR_leq[j]))
  denominator <- sum((p_y_ctl_leq_average[j] * (1 - p_y_ctl_leq_average[j])))
  weighted_RR_cum <- exp(numerator/denominator)
  log_weighted_RR_cum <- numerator/denominator
  
  # both weights
  numerator <- sum((p_y_ctl_est_bb[j] + p_y_ctl_est_bb[k]) * (p_y_ctl_leq_average[j] * (1 - p_y_ctl_leq_average[j])) * log(RR_leq[j]))
  denominator <- sum((p_y_ctl_est_bb[j] + p_y_ctl_est_bb[k]) * (p_y_ctl_leq_average[j] * (1 - p_y_ctl_leq_average[j])))
  weighted_RR_both <- exp(numerator/denominator)
  log_weighted_RR_both <- numerator/denominator
  
  # both weights RR greater than or equal to
  numerator <- sum((p_y_ctl_est_bb[j] + p_y_ctl_est_bb[k]) * (p_y_ctl_leq_average[j] * (1 - p_y_ctl_leq_average[j])) * log(RR_geq[j]))
  denominator <- sum((p_y_ctl_est_bb[j] + p_y_ctl_est_bb[k]) * (p_y_ctl_leq_average[j] * (1 - p_y_ctl_leq_average[j])))
  weighted_RR_ge_both <- exp(numerator/denominator)
  log_weighted_RR_ge_both <- numerator/denominator
  
  # #calculate weighted RR with overall probability as weights
  # numerator <- sum((p_y_est_bb[j] + p_y_est_bb[k]) * log(RR_leq[j]))
  # denominator <- sum(p_y_est_bb[j]) + sum(p_y_est_bb[k])
  # log_weighted_RR_Boverall <- numerator/denominator
  # weighted_RR_Boverall <- exp(numerator/denominator)
  
  # unweighted RR 
  unweighted_log_RR <- mean(log(RR_leq))
  unweighted_RR <- exp(mean(log(RR_leq)))
  
  # unweighted RR greater than or equal to  
  unweighted_log_RR_ge <- mean(log(RR_geq))
  unweighted_RR_ge <- exp(mean(log(RR_geq)))
  
  #calculate Net Benefit    
  l <- 1:c
  NB <- sum(p_y_trt[j] * p_y_ctl_geq_average[j]) - (sum(p_y_ctl[j] * p_y_trt_geq_average[j]))
  
  #save the weights 
  weights_split <- (p_y_ctl_est_bb[j] + p_y_ctl_est_bb[k]) / sum((p_y_ctl_est_bb[j] + p_y_ctl_est_bb[k]))
  weights_cum <- (p_y_ctl_leq_average[j] * (1 - p_y_ctl_leq_average[j])) / sum((p_y_ctl_leq_average[j] * (1 - p_y_ctl_leq_average[j])))
  weights_both <- ((p_y_ctl_est_bb[j] + p_y_ctl_est_bb[k]) * (p_y_ctl_leq_average[j] * (1 - p_y_ctl_leq_average[j]))) /
    sum((p_y_ctl_est_bb[j] + p_y_ctl_est_bb[k]) * (p_y_ctl_leq_average[j] * (1 - p_y_ctl_leq_average[j])))
  
  return(list(p_y_ctl = p_y_ctl, 
              p_y_trt = p_y_trt,
              RD_leq = RD_leq, 
              RD_geq = RD_geq, 
              RR_leq = RR_leq,
              RR_geq = RR_geq,
              OR = OR,
              log_RR_leq = log(RR_leq),
              log_RR_geq = log(RR_geq),
              log_OR = log(OR), 
              # weighted_RD = weighted_RD, 
              # weighted_RD_overall = weighted_RD_overall,
              weighted_RD_B = weighted_RD_B, 
              weighted_RD_cum = weighted_RD_cum, 
              weighted_RD_both = weighted_RD_both, 
              # weighted_RD_Boverall = weighted_RD_Boverall,
              unweighted_RD = unweighted_RD,
              # weighted_OR = weighted_OR,
              # log_weighted_OR = log_weighted_OR,
              # weighted_OR_overall = weighted_OR_overall,
              # log_weighted_OR_overall = log_weighted_OR_overall,
              weighted_OR_B = weighted_OR_B,
              log_weighted_OR_B = log_weighted_OR_B,
              weighted_OR_cum = weighted_OR_cum,
              log_weighted_OR_cum = log_weighted_OR_cum,
              weighted_OR_both = weighted_OR_both,
              log_weighted_OR_both = log_weighted_OR_both,
              # weighted_OR_Boverall = weighted_OR_Boverall,
              # log_weighted_OR_Boverall = log_weighted_OR_Boverall,
              unweighted_log_OR = unweighted_log_OR,
              unweighted_OR = unweighted_OR,
              # weighted_RR = weighted_RR,
              # log_weighted_RR = log_weighted_RR,
              # weighted_RR_overall = weighted_RR_overall,
              # log_weighted_RR_overall = log_weighted_RR_overall,
              weighted_RR_B = weighted_RR_B,
              log_weighted_RR_B = log_weighted_RR_B,
              weighted_RR_cum = weighted_RR_cum,
              log_weighted_RR_cum = log_weighted_RR_cum,
              weighted_RR_both = weighted_RR_both,
              log_weighted_RR_both = log_weighted_RR_both,
              weighted_RR_ge_both = weighted_RR_ge_both,
              log_weighted_RR_ge_both = log_weighted_RR_ge_both,
              # weighted_RR_Boverall = weighted_RR_Boverall,
              # log_weighted_RR_Boverall = log_weighted_RR_Boverall,
              unweighted_log_RR = unweighted_log_RR,
              unweighted_RR = unweighted_RR,
              unweighted_log_RR_ge = unweighted_log_RR_ge,
              unweighted_RR_ge = unweighted_RR_ge,
              NB = NB,
              weights_split = weights_split,
              weights_cum = weights_cum,
              weights_both = weights_both))
  
}


# 95% CI calculations 
# new function that uses bayesian bootstrap to average over population
# new function for adjusted CI for PPO model 

order_weighted_RD2_CI_adj_bb <- function(b, n_post_draws = 4000, dat, iprior, 
                                         outcome, trt, covs = 1) { 
  
  # first calculate P(y >= i | ctl) and P(y >= i | trt) for all i
  posterior_draws <- as.matrix(b$draws)
  c <- length(b$ylevels)
  PPO_ctl_geq_averages <- matrix(data = NA, nrow = n_post_draws, ncol = c) 
  PPO_trt_geq_averages <- matrix(data = NA, nrow = n_post_draws, ncol = c) 
  PPO_ctl_leq_averages <- matrix(data = NA, nrow = n_post_draws, ncol = c) 
  PPO_trt_leq_averages <- matrix(data = NA, nrow = n_post_draws, ncol = c) 
  
  wRD_control <- matrix(data = NA, nrow = n_post_draws, ncol = 1)
  wRD_overall <- matrix(data = NA, nrow = n_post_draws, ncol = 1)
  wOR_control <- matrix(data = NA, nrow = n_post_draws, ncol = 1)
  wOR_overall <- matrix(data = NA, nrow = n_post_draws, ncol = 1)
  log_wOR_control <- matrix(data = NA, nrow = n_post_draws, ncol = 1)
  log_wOR_overall <- matrix(data = NA, nrow = n_post_draws, ncol = 1)
  wRR_control <- matrix(data = NA, nrow = n_post_draws, ncol = 1)
  wRR_overall <- matrix(data = NA, nrow = n_post_draws, ncol = 1)
  log_wRR_control <- matrix(data = NA, nrow = n_post_draws, ncol = 1)
  log_wRR_overall <- matrix(data = NA, nrow = n_post_draws, ncol = 1)
  
  wRD_bcontrol <- matrix(data = NA, nrow = n_post_draws, ncol = 1)
  wRD_boverall <- matrix(data = NA, nrow = n_post_draws, ncol = 1)
  wRD_cum <- matrix(data = NA, nrow = n_post_draws, ncol = 1)
  wRD_both <- matrix(data = NA, nrow = n_post_draws, ncol = 1)
  wOR_bcontrol <- matrix(data = NA, nrow = n_post_draws, ncol = 1)
  wOR_boverall <- matrix(data = NA, nrow = n_post_draws, ncol = 1)
  wOR_cum <- matrix(data = NA, nrow = n_post_draws, ncol = 1)
  wOR_both <- matrix(data = NA, nrow = n_post_draws, ncol = 1)
  log_wOR_bcontrol <- matrix(data = NA, nrow = n_post_draws, ncol = 1)
  log_wOR_boverall <- matrix(data = NA, nrow = n_post_draws, ncol = 1)
  log_wOR_cum <- matrix(data = NA, nrow = n_post_draws, ncol = 1)
  log_wOR_both <- matrix(data = NA, nrow = n_post_draws, ncol = 1)
  wRR_bcontrol <- matrix(data = NA, nrow = n_post_draws, ncol = 1)
  wRR_boverall <- matrix(data = NA, nrow = n_post_draws, ncol = 1)
  wRR_cum <- matrix(data = NA, nrow = n_post_draws, ncol = 1)
  wRR_both <- matrix(data = NA, nrow = n_post_draws, ncol = 1)
  wRR_ge_both <- matrix(data = NA, nrow = n_post_draws, ncol = 1)
  log_wRR_bcontrol <- matrix(data = NA, nrow = n_post_draws, ncol = 1)
  log_wRR_boverall <- matrix(data = NA, nrow = n_post_draws, ncol = 1)
  log_wRR_cum <- matrix(data = NA, nrow = n_post_draws, ncol = 1)
  log_wRR_both <- matrix(data = NA, nrow = n_post_draws, ncol = 1)
  log_wRR_ge_both <- matrix(data = NA, nrow = n_post_draws, ncol = 1)
  
  weights_split <- matrix(data = NA, nrow = n_post_draws, ncol = c-1)
  weights_cum <- matrix(data = NA, nrow = n_post_draws, ncol = c-1)
  weights_both <- matrix(data = NA, nrow = n_post_draws, ncol = c-1)
  
  NB <- matrix(data = NA, nrow = n_post_draws, ncol = 1)
  RD_leq_save <- matrix(data = NA, nrow = n_post_draws, ncol = c-1)
  OR_save <- matrix(data = NA, nrow = n_post_draws, ncol = c-1)
  RR_save <- matrix(data = NA, nrow = n_post_draws, ncol = c-1)
  RR_save_geq <- matrix(data = NA, nrow = n_post_draws, ncol = c-1)
  log_OR_save <- matrix(data = NA, nrow = n_post_draws, ncol = c-1)
  log_RR_save <- matrix(data = NA, nrow = n_post_draws, ncol = c-1)
  log_RR_ge_save <- matrix(data = NA, nrow = n_post_draws, ncol = c-1)
  
  foreach (s = 1:n_post_draws) %do% {
    
    predicted_trt_geq <- matrix(data = NA, nrow = nrow(dat), ncol = c-1) 
    predicted_ctl_geq <- matrix(data = NA, nrow = nrow(dat), ncol = c-1) 
    predicted_trt_leq <- matrix(data = NA, nrow = nrow(dat), ncol = c-1) 
    predicted_ctl_leq <- matrix(data = NA, nrow = nrow(dat), ncol = c-1) 
    
    # if (iprior == 0){
    #   coefs <- posterior_draws[i, ]
    # } else {
    #   coefs <- c(-coef(b)[1:3], coef(b)[4:9])
    # }
    
    coefs <- posterior_draws[s, ]
    covs <- b$Design$colnames
    num_covs <- length(covs)
    covariate_mat <- as.data.frame(matrix(data = NA, nrow = nrow(dat), ncol = num_covs - 1))
    
    # assume the first covariate is the treatment arm 
    for(i in 2:num_covs) {
      if (is.numeric(dat[[b$Design$name[i]]]) == T) {
        covariate_mat[,i-1] <- coefs[[covs[i]]] * (dat[[b$Design$name[i]]])
      } else {
        # only works for numeric and binary
        covariate_mat[,i-1] <- coefs[[covs[i]]] * (dat[[b$Design$name[i]]] == b$Design$parms[[b$Design$name[i]]][2])
      }
    }
    
    for (i in 2:c) {
      if (i == 2) {
        predicted_trt_geq[,i-1] <- plogis(coefs[i-1] + coefs[[trt]] * 0.5 + rowSums(covariate_mat))
        predicted_ctl_geq[,i-1] <- plogis(coefs[i-1] + coefs[[trt]] * -0.5 + rowSums(covariate_mat))
      } else {
        phrase <- paste0("y>=", i)
        coefs_new <- coefs[str_detect(names(coefs), phrase)]
        predicted_trt_geq[,i-1] <- plogis(coefs_new[1] + coefs_new[2] * 0.5 + coefs[[trt]] * 0.5 + rowSums(covariate_mat))
        predicted_ctl_geq[,i-1] <- plogis(coefs_new[1] + coefs_new[2] * -0.5 + coefs[[trt]] * -0.5 + rowSums(covariate_mat))
      }
    }
    
    
    predicted_trt_leq <-  1 - predicted_trt_geq
    predicted_ctl_leq <-  1 - predicted_ctl_geq
    
    p_y_trt_est <- matrix(data = NA, nrow = nrow(dat), ncol = c) 
    p_y_ctl_est <- matrix(data = NA, nrow = nrow(dat), ncol = c) 
    for (i in 1:c) {
      if (i == 1) {
        p_y_ctl_est[,1] <-  predicted_ctl_leq[,1]
        p_y_trt_est[,1] <-  predicted_trt_leq[,1]
      } else if (i == 2){
        p_y_ctl_est[,2] <- predicted_ctl_leq[,i] - predicted_ctl_leq[,1]
        p_y_trt_est[,2] <- predicted_trt_leq[,i] - predicted_trt_leq[,1]
      } else if (i == c) {
        p_y_ctl_est[,c] <- 1 - rowSums(p_y_ctl_est[,1:(c-1)])
        p_y_trt_est[,c] <- 1 - rowSums(p_y_trt_est[,1:(c-1)])
      } else {
        p_y_ctl_est[,i] <-  predicted_ctl_leq[,i] - apply(p_y_ctl_est[,1:(i-1)], 1, sum)
        p_y_trt_est[,i] <-  predicted_trt_leq[,i] - apply(p_y_trt_est[,1:(i-1)], 1, sum)
      }
    }
    
    # correct for any negative probabilities 
    if (any(p_y_ctl_est < 0)) {
      p_y_ctl_est[which(p_y_ctl_est < 0)] <- 0 
      p_y_ctl_est <- p_y_ctl_est/rowSums(p_y_ctl_est)
    }
    
    if (any(p_y_trt_est < 0)) {
      p_y_trt_est[which(p_y_trt_est < 0)] <- 0 
      p_y_trt_est <- p_y_trt_est/rowSums(p_y_trt_est)
    }
    
    p_y_est <- (p_y_ctl_est + p_y_trt_est)/2
    
    # get weights from the data 
    p_y_ctl <- unname(table(dat[[trt]], dat[[outcome]])[1,])/sum(dat[[trt]] == -0.5)
    p_y_trt <- unname(table(dat[[trt]], dat[[outcome]])[2,])/sum(dat[[trt]] == 0.5)
    p_y <- as.vector(unname(table(dat[[outcome]]))/nrow(dat))
    
    weights <- t(rdirichlet(1, rep(1, nrow(dat))))
    
    p_y_ctl_leq_average <- apply(predicted_ctl_leq, 2, weighted.mean, w = weights)
    p_y_trt_leq_average <- apply(predicted_trt_leq, 2, weighted.mean, w = weights)
    p_y_ctl_geq_average <- apply(predicted_ctl_geq, 2, weighted.mean, w = weights)
    p_y_trt_geq_average <- apply(predicted_trt_geq, 2, weighted.mean, w = weights)
    
    p_y_est_average <- apply(p_y_est, 2, weighted.mean, w = weights)
    p_y_ctl_est_average <- apply(p_y_ctl_est, 2, weighted.mean, w = weights)
    
    RD_leq <- p_y_trt_leq_average - p_y_ctl_leq_average
    RD_leq_save[s,] <- RD_leq
    RD_geq <- p_y_ctl_geq_average - p_y_trt_geq_average
    
    
    OR <- (p_y_trt_leq_average/ (1- p_y_trt_leq_average)) / (p_y_ctl_leq_average / (1- p_y_ctl_leq_average))
    OR_save[s,] <- OR
    log_OR_save[s,] <- log(OR)
    
    RR_leq <-  p_y_trt_leq_average / p_y_ctl_leq_average
    RR_save[s,] <- RR_leq
    log_RR_save[s,] <- log(RR_leq)
    RR_geq <- p_y_ctl_geq_average / p_y_trt_geq_average
    RR_save_geq[s,] <- RR_geq
    log_RR_ge_save[s,] <- log(RR_geq)
    
    j <- 1:(c-1)
    k <- 2:c
    # numerator <- sum((p_y_ctl[j] + p_y_ctl[k]) * (p_y_trt_leq_average[j] - p_y_ctl_leq_average[j])) 
    # denominator <- sum(p_y_ctl[j]) + sum(p_y_ctl[k])
    # wRD_control[s,1] <- numerator/denominator
    # 
    # #calculated weighted RD with overall probability as weights
    # numerator <- sum((p_y[j] + p_y[k]) * (p_y_trt_leq_average[j] - p_y_ctl_leq_average[j]))
    # denominator <- sum(p_y[j]) + sum(p_y[k])
    # wRD_overall[s,1] <- numerator/denominator
    
    # bayesian approach 
    numerator <- sum((p_y_ctl_est_average[j] + p_y_ctl_est_average[k]) * (p_y_trt_leq_average[j] - p_y_ctl_leq_average[j])) 
    denominator <- sum(p_y_ctl_est_average[j]) + sum(p_y_ctl_est_average[k])
    wRD_bcontrol[s,1] <- numerator/denominator
    
    #cumulative weights 
    numerator <- sum((p_y_ctl_leq_average[j] * (1 - p_y_ctl_leq_average[j])) * (p_y_trt_leq_average[j] - p_y_ctl_leq_average[j])) 
    denominator <- sum((p_y_ctl_leq_average[j] * (1 - p_y_ctl_leq_average[j])))
    wRD_cum[s,1] <- numerator/denominator
    
    #both weights 
    numerator <- sum((p_y_ctl_est_average[j] + p_y_ctl_est_average[k]) * (p_y_ctl_leq_average[j] * (1 - p_y_ctl_leq_average[j])) * (p_y_trt_leq_average[j] - p_y_ctl_leq_average[j])) 
    denominator <- sum((p_y_ctl_est_average[j] + p_y_ctl_est_average[k]) * (p_y_ctl_leq_average[j] * (1 - p_y_ctl_leq_average[j])))
    wRD_both[s,1] <- numerator/denominator
    
    # #calculated weighted RD with overall probability as weights
    # numerator <- sum((p_y_est_average[j] + p_y_est_average[k]) * (p_y_trt_leq_average[j] - p_y_ctl_leq_average[j]))
    # denominator <- sum(p_y_est_average[j]) + sum(p_y_est_average[k])
    # wRD_boverall[s,1] <- numerator/denominator
    
    # #calculate log weighted OR with control probability as weights
    # numerator <- sum((p_y_ctl[j] + p_y_ctl[k]) * log(OR[j]))
    # denominator <- sum(p_y_ctl[j]) + sum(p_y_ctl[k])
    # log_wOR_control[s,1] <- numerator/denominator
    # wOR_control[s,1] <- exp(numerator/denominator)
    # 
    # #calculate log weighted OR with overall probability as weights
    # numerator <- sum((p_y[j] + p_y[k])  * log(OR[j])) 
    # denominator <- sum(p_y[j]) + sum(p_y[k])
    # log_wOR_overall[s,1] <- numerator/denominator
    # wOR_overall[s,1] <- exp(numerator/denominator)
    
    # bayesian approach 
    #calculate log weighted OR with control probability as weights
    numerator <- sum((p_y_ctl_est_average[j] + p_y_ctl_est_average[k]) * log(OR[j]))
    denominator <- sum(p_y_ctl_est_average[j]) + sum(p_y_ctl_est_average[k])
    log_wOR_bcontrol[s,1] <- numerator/denominator
    wOR_bcontrol[s,1] <- exp(numerator/denominator)
    
    # cumulative weights 
    numerator <- sum((p_y_ctl_leq_average[j] * (1 - p_y_ctl_leq_average[j])) * log(OR[j]))
    denominator <- sum((p_y_ctl_leq_average[j] * (1 - p_y_ctl_leq_average[j])))
    log_wOR_cum[s,1] <- numerator/denominator
    wOR_cum[s,1] <- exp(numerator/denominator)
    
    # both weights 
    numerator <- sum((p_y_ctl_est_average[j] + p_y_ctl_est_average[k]) * (p_y_ctl_leq_average[j] * (1 - p_y_ctl_leq_average[j])) * log(OR[j]))
    denominator <- sum((p_y_ctl_est_average[j] + p_y_ctl_est_average[k]) * (p_y_ctl_leq_average[j] * (1 - p_y_ctl_leq_average[j])))
    log_wOR_both[s,1] <- numerator/denominator
    wOR_both[s,1] <- exp(numerator/denominator)
    
    # #calculate log weighted OR with overall probability as weights
    # numerator <- sum((p_y_est_average[j] + p_y_est_average[k]) * log(OR[j])) 
    # denominator <- sum(p_y_est_average[j]) + sum(p_y_est_average[k])
    # log_wOR_boverall[s,1] <- numerator/denominator
    # wOR_boverall[s,1] <- exp(numerator/denominator)
    
    # #calculate log weighted RR with control probability as weights
    # numerator <- sum((p_y_ctl[j] + p_y_ctl[k]) * log(RR_leq[j]))
    # denominator <- sum(p_y_ctl[j]) + sum(p_y_ctl[k])
    # log_wRR_control[s,1] <- numerator/denominator
    # wRR_control[s,1] <- exp(numerator/denominator)
    # 
    # #calculate log weighted RR with overall probability as weights
    # numerator <- sum((p_y[j] + p_y[k])  * log(RR_leq[j]))
    # denominator <- sum(p_y[j]) + sum(p_y[k])
    # log_wRR_overall[s,1] <- numerator/denominator
    # wRR_overall[s,1] <- exp(numerator/denominator)
    
    # bayesian approach
    #calculate log weighted RR with control probability as weights
    numerator <- sum((p_y_ctl_est_average[j] + p_y_ctl_est_average[k]) * log(RR_leq[j]))
    denominator <- sum(p_y_ctl_est_average[j]) + sum(p_y_ctl_est_average[k])
    log_wRR_bcontrol[s,1] <- numerator/denominator
    wRR_bcontrol[s,1] <- exp(numerator/denominator)
    
    # cumulative weights 
    numerator <- sum((p_y_ctl_leq_average[j] * (1 - p_y_ctl_leq_average[j])) * log(RR_leq[j]))
    denominator <- sum((p_y_ctl_leq_average[j] * (1 - p_y_ctl_leq_average[j])))
    log_wRR_cum[s,1] <- numerator/denominator
    wRR_cum[s,1] <- exp(numerator/denominator)
    
    # both weights 
    numerator <- sum((p_y_ctl_est_average[j] + p_y_ctl_est_average[k]) * (p_y_ctl_leq_average[j] * (1 - p_y_ctl_leq_average[j])) * log(RR_leq[j]))
    denominator <- sum((p_y_ctl_est_average[j] + p_y_ctl_est_average[k]) * (p_y_ctl_leq_average[j] * (1 - p_y_ctl_leq_average[j])))
    log_wRR_both[s,1] <- numerator/denominator
    wRR_both[s,1] <- exp(numerator/denominator)
    
    # both weights 
    numerator <- sum((p_y_ctl_est_average[j] + p_y_ctl_est_average[k]) * (p_y_ctl_leq_average[j] * (1 - p_y_ctl_leq_average[j])) * log(RR_geq[j]))
    denominator <- sum((p_y_ctl_est_average[j] + p_y_ctl_est_average[k]) * (p_y_ctl_leq_average[j] * (1 - p_y_ctl_leq_average[j])))
    log_wRR_ge_both[s,1] <- numerator/denominator
    wRR_ge_both[s,1] <- exp(numerator/denominator)
    
    # #calculate log weighted RR with overall probability as weights
    # numerator <- sum((p_y_est_average[j] + p_y_est_average[k]) * log(RR_leq[j]))
    # denominator <- sum(p_y_est_average[j]) + sum(p_y_est_average[k])
    # log_wRR_boverall[s,1] <- numerator/denominator
    # wRR_boverall[s,1] <- exp(numerator/denominator)
    
    #calculate WR    
    l <- 1:c
    NB[s,1] <- sum(p_y_trt[j] * p_y_ctl_geq_average[j]) - (sum(p_y_ctl[j] * p_y_trt_geq_average[j]))
    
    weights_split[s,] <- (p_y_ctl_est_average[j] + p_y_ctl_est_average[k]) / sum(p_y_ctl_est_average[j] + p_y_ctl_est_average[k])
    weights_cum[s,] <- (p_y_ctl_leq_average[j] * (1 - p_y_ctl_leq_average[j])) / sum((p_y_ctl_leq_average[j] * (1 - p_y_ctl_leq_average[j])))
    weights_both[s,] <- (p_y_ctl_est_average[j] + p_y_ctl_est_average[k]) * (p_y_ctl_leq_average[j] * (1 - p_y_ctl_leq_average[j])) / 
      sum((p_y_ctl_est_average[j] + p_y_ctl_est_average[k]) * (p_y_ctl_leq_average[j] * (1 - p_y_ctl_leq_average[j])))
    
  }
  
  RD_leq_95CI <- apply(RD_leq_save, 2, function (x) quantile(x, probs = c(0.025, 0.975)))
  RD_leq_postprob <- colSums(RD_leq_save > 0)/n_post_draws
  
  #unweighted_RD
  unweighted_RD <- rowMeans(RD_leq_save)
  unweighted_RD_95CI <- quantile(unweighted_RD, probs = c(0.025, 0.975))
  
  OR_95CI <- apply(OR_save, 2, function (x) quantile(x, probs = c(0.025, 0.975)))
  OR_postprob <- colSums(OR_save > 1)/n_post_draws
  
  log_OR_95CI <- apply(log_OR_save, 2, function (x) quantile(x, probs = c(0.025, 0.975)))
  
  #unweighted_OR
  unweighted_log_OR <- rowMeans(log_OR_save)
  unweighted_log_OR_95CI <- quantile(unweighted_log_OR, probs = c(0.025, 0.975))
  unweighted_OR_95CI <- quantile(exp(unweighted_log_OR), probs = c(0.025, 0.975))
  
  RR_leq_95CI <- apply(RR_save, 2, function (x) quantile(x, probs = c(0.025, 0.975)))
  RR_leq_postprob <- colSums(RR_save > 1)/n_post_draws
  
  log_RR_leq_95CI <- apply(log_RR_save, 2, function (x) quantile(x, probs = c(0.025, 0.975)))
  
  #unweighted_RR
  unweighted_log_RR <- rowMeans(log_RR_save)
  unweighted_log_RR_95CI <- quantile(unweighted_log_RR, probs = c(0.025, 0.975))
  unweighted_RR_95CI <- quantile(exp(unweighted_log_RR), probs = c(0.025, 0.975))
  
  #unweighted_RR
  unweighted_log_RR_ge <- rowMeans(log_RR_ge_save)
  unweighted_log_RR_ge_95CI <- quantile(unweighted_log_RR_ge, probs = c(0.025, 0.975))
  unweighted_RR_ge_95CI <- quantile(exp(unweighted_log_RR_ge), probs = c(0.025, 0.975))
  
  RR_geq_95CI <- apply(RR_save_geq, 2, function (x) quantile(x, probs = c(0.025, 0.975)))
  RR_geq_postprob <- colSums(RR_save_geq > 1)/n_post_draws
  
  # wRD_control_95CI <- quantile(wRD_control, probs = c(0.025, 0.975))
  # wRD_control_postprob <- sum(wRD_control > 0)/n_post_draws
  # 
  # wRD_overall_95CI <- quantile(wRD_overall, probs = c(0.025, 0.975))
  # wRD_overall_postprob <- sum(wRD_overall > 0)/n_post_draws
  
  wRD_bcontrol_95CI <- quantile(wRD_bcontrol, probs = c(0.025, 0.975))
  wRD_bcontrol_postprob <- sum(wRD_bcontrol > 0)/n_post_draws
  
  wRD_cum_95CI <- quantile(wRD_cum, probs = c(0.025, 0.975))
  wRD_cum_postprob <- sum(wRD_cum > 0)/n_post_draws
  
  wRD_both_95CI <- quantile(wRD_both, probs = c(0.025, 0.975))
  wRD_both_postprob <- sum(wRD_both > 0)/n_post_draws
  
  # wRD_boverall_95CI <- quantile(wRD_boverall, probs = c(0.025, 0.975))
  # wRD_boverall_postprob <- sum(wRD_boverall > 0)/n_post_draws
  # 
  # wOR_control_95CI <- quantile(wOR_control, probs = c(0.025, 0.975))
  # wOR_control_postprob <- sum(wOR_control > 1)/n_post_draws
  # 
  # wOR_overall_95CI <- quantile(wOR_overall, probs = c(0.025, 0.975))
  # wOR_overall_postprob <- sum(wOR_overall > 1)/n_post_draws
  
  wOR_bcontrol_95CI <- quantile(wOR_bcontrol, probs = c(0.025, 0.975))
  wOR_bcontrol_postprob <- sum(wOR_bcontrol > 1)/n_post_draws
  
  wOR_cum_95CI <- quantile(wOR_cum, probs = c(0.025, 0.975))
  wOR_cum_postprob <- sum(wOR_cum > 1)/n_post_draws
  
  wOR_both_95CI <- quantile(wOR_both, probs = c(0.025, 0.975))
  wOR_both_postprob <- sum(wOR_both > 1)/n_post_draws
  
  # wOR_boverall_95CI <- quantile(wOR_boverall, probs = c(0.025, 0.975))
  # wOR_boverall_postprob <- sum(wOR_boverall > 1)/n_post_draws
  # 
  # log_wOR_control_95CI <- quantile(log_wOR_control, probs = c(0.025, 0.975))
  # log_wOR_control_postprob <- sum(log_wOR_control > 0)/n_post_draws
  # 
  # log_wOR_overall_95CI <- quantile(log_wOR_overall, probs = c(0.025, 0.975))
  # log_wOR_overall_postprob <- sum(log_wOR_overall > 0)/n_post_draws
  
  log_wOR_bcontrol_95CI <- quantile(log_wOR_bcontrol, probs = c(0.025, 0.975))
  log_wOR_bcontrol_postprob <- sum(log_wOR_bcontrol > 0)/n_post_draws
  
  log_wOR_cum_95CI <- quantile(log_wOR_cum, probs = c(0.025, 0.975))
  log_wOR_cum_postprob <- sum(log_wOR_cum > 0)/n_post_draws
  
  log_wOR_both_95CI <- quantile(log_wOR_both, probs = c(0.025, 0.975))
  log_wOR_both_postprob <- sum(log_wOR_both > 0)/n_post_draws
  
  # log_wOR_boverall_95CI <- quantile(log_wOR_boverall, probs = c(0.025, 0.975))
  # log_wOR_boverall_postprob <- sum(log_wOR_boverall > 0)/n_post_draws
  # 
  # wRR_control_95CI <- quantile(wRR_control, probs = c(0.025, 0.975))
  # wRR_control_postprob <- sum(wRR_control > 1)/n_post_draws
  # 
  # wRR_overall_95CI <- quantile(wRR_overall, probs = c(0.025, 0.975))
  # wRR_overall_postprob <- sum(wRR_overall > 1)/n_post_draws
  # 
  # log_wRR_control_95CI <- quantile(log_wRR_control, probs = c(0.025, 0.975))
  # log_wRR_control_postprob <- sum(log_wRR_control > 0)/n_post_draws
  # 
  # log_wRR_overall_95CI <- quantile(log_wRR_overall, probs = c(0.025, 0.975))
  # log_wRR_overall_postprob <- sum(log_wRR_overall > 0)/n_post_draws
  
  wRR_bcontrol_95CI <- quantile(wRR_bcontrol, probs = c(0.025, 0.975))
  wRR_bcontrol_postprob <- sum(wRR_bcontrol > 1)/n_post_draws
  
  wRR_cum_95CI <- quantile(wRR_cum, probs = c(0.025, 0.975))
  wRR_cum_postprob <- sum(wRR_cum > 1)/n_post_draws
  
  wRR_both_95CI <- quantile(wRR_both, probs = c(0.025, 0.975))
  wRR_both_postprob <- sum(wRR_both > 1)/n_post_draws
  
  wRR_ge_both_95CI <- quantile(wRR_ge_both, probs = c(0.025, 0.975))
  wRR_ge_both_postprob <- sum(wRR_ge_both > 1)/n_post_draws
  
  # wRR_boverall_95CI <- quantile(wRR_boverall, probs = c(0.025, 0.975))
  # wRR_boverall_postprob <- sum(wRR_boverall > 1)/n_post_draws
  
  log_wRR_bcontrol_95CI <- quantile(log_wRR_bcontrol, probs = c(0.025, 0.975))
  log_wRR_bcontrol_postprob <- sum(log_wRR_bcontrol > 0)/n_post_draws
  
  log_wRR_cum_95CI <- quantile(log_wRR_cum, probs = c(0.025, 0.975))
  log_wRR_cum_postprob <- sum(log_wRR_cum > 0)/n_post_draws
  
  log_wRR_both_95CI <- quantile(log_wRR_both, probs = c(0.025, 0.975))
  log_wRR_both_postprob <- sum(log_wRR_both > 0)/n_post_draws
  
  log_wRR_ge_both_95CI <- quantile(log_wRR_ge_both, probs = c(0.025, 0.975))
  log_wRR_ge_both_postprob <- sum(log_wRR_ge_both > 0)/n_post_draws
  
  # log_wRR_boverall_95CI <- quantile(log_wRR_boverall, probs = c(0.025, 0.975))
  # log_wRR_boverall_postprob <- sum(log_wRR_boverall > 0)/n_post_draws
  
  NB_95CI <- quantile(NB, probs = c(0.025, 0.975))
  NB_postprob <- sum(NB > 0)/n_post_draws
  
  weights_split_95CI <- apply(weights_split, 2, quantile, probs = c(0.025, 0.975))
  weights_cum_95CI <- apply(weights_cum, 2, quantile, probs = c(0.025, 0.975))
  weights_both_95CI <- apply(weights_both, 2, quantile, probs = c(0.025, 0.975))
  
  
  
  return(list(RD_leq_95CI = RD_leq_95CI, 
              RD_leq_postprob = RD_leq_postprob, 
              OR_95CI = OR_95CI, 
              log_OR_95CI = log_OR_95CI, 
              OR_postprob = OR_postprob, 
              RR_leq_95CI = RR_leq_95CI, 
              log_RR_leq_95CI = log_RR_leq_95CI, 
              RR_leq_postprob = RR_leq_postprob, 
              RR_geq_95CI = RR_geq_95CI,
              # wRD_control_95CI = wRD_control_95CI, 
              # wRD_control_postprob = wRD_control_postprob, 
              # wRD_overall_95CI = wRD_overall_95CI, 
              # wRD_overall_postprob = wRD_overall_postprob, 
              wRD_bcontrol_95CI = wRD_bcontrol_95CI, 
              wRD_bcontrol_postprob = wRD_bcontrol_postprob, 
              wRD_cum_95CI = wRD_cum_95CI, 
              wRD_cum_postprob = wRD_cum_postprob,
              wRD_both_95CI = wRD_both_95CI, 
              wRD_both_postprob = wRD_both_postprob,
              # wRD_boverall_95CI = wRD_boverall_95CI, 
              # wRD_boverall_postprob = wRD_boverall_postprob, 
              # wOR_control_95CI = wOR_control_95CI, 
              # wOR_control_postprob = wOR_control_postprob, 
              # wOR_overall_95CI = wOR_overall_95CI, 
              # wOR_overall_postprob = wOR_overall_postprob, 
              # log_wOR_control_95CI = log_wOR_control_95CI, 
              # log_wOR_control_postprob = log_wOR_control_postprob, 
              # log_wOR_overall_95CI = log_wOR_overall_95CI, 
              # log_wOR_overall_postprob = log_wOR_overall_postprob, 
              wOR_bcontrol_95CI = wOR_bcontrol_95CI,
              wOR_bcontrol_postprob = wOR_bcontrol_postprob,
              wOR_cum_95CI = wOR_cum_95CI,
              wOR_cum_postprob = wOR_cum_postprob,
              wOR_both_95CI = wOR_both_95CI,
              wOR_both_postprob = wOR_both_postprob,
              # wOR_boverall_95CI = wOR_boverall_95CI, 
              # wOR_boverall_postprob = wOR_boverall_postprob, 
              log_wOR_bcontrol_95CI = log_wOR_bcontrol_95CI, 
              log_wOR_bcontrol_postprob = log_wOR_bcontrol_postprob, 
              log_wOR_cum_95CI = log_wOR_cum_95CI, 
              log_wOR_cum_postprob = log_wOR_cum_postprob, 
              log_wOR_both_95CI = log_wOR_both_95CI, 
              log_wOR_both_postprob = log_wOR_both_postprob, 
              # log_wOR_boverall_95CI = log_wOR_boverall_95CI, 
              # log_wOR_boverall_postprob = log_wOR_boverall_postprob, 
              # wRR_control_95CI = wRR_control_95CI, 
              # wRR_control_postprob = wRR_control_postprob, 
              # wRR_overall_95CI = wRR_overall_95CI, 
              # wRR_overall_postprob = wRR_overall_postprob, 
              # log_wRR_control_95CI = log_wRR_control_95CI, 
              # log_wRR_control_postprob = log_wRR_control_postprob, 
              # log_wRR_overall_95CI = log_wRR_overall_95CI, 
              # log_wRR_overall_postprob = log_wRR_overall_postprob, 
              wRR_bcontrol_95CI = wRR_bcontrol_95CI, 
              wRR_bcontrol_postprob = wRR_bcontrol_postprob, 
              wRR_cum_95CI = wRR_cum_95CI, 
              wRR_cum_postprob = wRR_cum_postprob,
              wRR_both_95CI = wRR_both_95CI, 
              wRR_both_postprob = wRR_both_postprob,
              wRR_ge_both_95CI = wRR_ge_both_95CI, 
              wRR_ge_both_postprob = wRR_ge_both_postprob,
              # wRR_boverall_95CI = wRR_boverall_95CI, 
              # wRR_boverall_postprob = wRR_boverall_postprob, 
              log_wRR_bcontrol_95CI = log_wRR_bcontrol_95CI, 
              log_wRR_bcontrol_postprob = log_wRR_bcontrol_postprob, 
              log_wRR_cum_95CI = log_wRR_cum_95CI, 
              log_wRR_cum_postprob = log_wRR_cum_postprob, 
              log_wRR_both_95CI = log_wRR_both_95CI, 
              log_wRR_both_postprob = log_wRR_both_postprob, 
              log_wRR_ge_both_95CI = log_wRR_ge_both_95CI, 
              log_wRR_ge_both_postprob = log_wRR_ge_both_postprob, 
              # log_wRR_boverall_95CI = log_wRR_boverall_95CI, 
              # log_wRR_boverall_postprob = log_wRR_boverall_postprob, 
              unweighted_RD_95CI = unweighted_RD_95CI,
              unweighted_log_OR_95CI = unweighted_log_OR_95CI,
              unweighted_OR_95CI = unweighted_OR_95CI,
              unweighted_log_RR_95CI = unweighted_log_RR_95CI,
              unweighted_RR_95CI = unweighted_RR_95CI,
              unweighted_log_RR_ge_95CI = unweighted_log_RR_ge_95CI,
              unweighted_RR_ge_95CI = unweighted_RR_ge_95CI,
              NB_95CI = NB_95CI, 
              NB_postprob = NB_postprob,
              weights_split_95CI = weights_split_95CI,
              weights_cum_95CI = weights_cum_95CI,
              weights_both_95CI = weights_both_95CI))
  
}




# proportional odds models 
# new function for weighted risk difference 
# assumes a proportional odds model 
# assumes the treatment is coded as -0.5 for control and 0.5 for treatment 
order_weighted_RD2_PO <- function(b = b, dat, trt, outcome) { 
  # first calculate P(y >= i | ctl) and P(y >= i | trt) for all i 
  PO_ctl_geq <- NA 
  PO_trt_geq <- NA 
  PO_ctl_leq <- NA 
  PO_trt_leq <- NA 
  p_y_ctl_est <- NA
  p_y_trt_est <- NA
  c <- length(b$ylevels)
  
  for (i in 2:c) {
    # for the first level there is not a partial term, so only need first term and metformin term 
    PO_ctl_geq[i] <- plogis(coef(b)[i - 1] - 0.5 * coef(b)[c]) 
    PO_trt_geq[i] <- plogis(coef(b)[i - 1] + 0.5 * coef(b)[c])
    PO_ctl_leq[i - 1] <- 1 - PO_ctl_geq[i]
    PO_trt_leq[i - 1] <- 1 - PO_trt_geq[i]
  }
  
  # get weights from the data 
  p_y_ctl <- unname(table(dat[[trt]], dat[[outcome]])[1,])/sum(dat[[trt]] == -0.5)
  p_y_trt <- unname(table(dat[[trt]], dat[[outcome]])[2,])/sum(dat[[trt]] == 0.5)
  p_y <- as.vector(unname(table(dat[[outcome]])))/nrow(dat)
  
  # calculate p(Y = i | control) for i = 1, ... , c
  for (i in 1:c) {
    if (i == 1) {
      p_y_ctl_est[1] <-  PO_ctl_leq[1]
      p_y_trt_est[1] <-  PO_trt_leq[1]
    } else if (i == c) {
      p_y_ctl_est[c] <- 1 - sum(p_y_ctl_est[1:(c-1)])
      p_y_trt_est[c] <- 1 - sum(p_y_trt_est[1:(c-1)])
    } else {
      p_y_ctl_est[i] <-  PO_ctl_leq[i] - sum(p_y_ctl_est[1:(i-1)])
      p_y_trt_est[i] <-  PO_trt_leq[i] - sum(p_y_trt_est[1:(i-1)])
    }
  }
  
  p_y_est <- (p_y_ctl_est + p_y_trt_est)/2
  
  RD_leq <- PO_trt_leq - PO_ctl_leq
  unweighted_RD <- mean(RD_leq)
  RD_geq <- PO_ctl_geq - PO_trt_geq
  
  OR <- (PO_trt_leq / (1- PO_trt_leq)) / (PO_ctl_leq / (1- PO_ctl_leq))
  log_OR <- log(OR)
  unweighted_log_OR <- mean(log_OR)
  unweighted_OR <- exp(mean(log_OR))
  
  log_RR_le <- log(PO_trt_leq/PO_ctl_leq)
  unweighted_log_RR <- mean(log_RR_le)
  unweighted_RR <- exp(mean(log_RR_le))
  log_RR_ge <- log(PO_ctl_geq/PO_trt_geq)
  unweighted_log_RR_ge <- mean(log_RR_ge, na.rm = T)
  unweighted_RR_ge <- exp(mean(log_RR_ge, na.rm = T))
  
  j <- 1:(c-1)
  k <- 2:c
  # numerator <- sum((p_y_ctl[j] + p_y_ctl[k]) * (PO_trt_leq[j] - PO_ctl_leq[j]))
  # denominator <- sum(p_y_ctl[j]) + sum(p_y_ctl[k])
  # weighted_RD <- numerator/denominator
  
  # #calculated weighted RD with overall probability as weights
  # p_y <- (p_y_trt + p_y_ctl)/2
  # numerator <- sum((p_y[j] + p_y[k]) * (PO_trt_leq[j] - PO_ctl_leq[j]))
  # denominator <- sum(p_y[j]) + sum(p_y[k])
  # weighted_RD_overall <- numerator/denominator
  
  # fully bayesian approach RD
  numerator <- sum((p_y_ctl_est[j] + p_y_ctl_est[k]) * (PO_trt_leq[j] - PO_ctl_leq[j]))
  denominator <- sum(p_y_ctl_est[j]) + sum(p_y_ctl_est[k])
  weighted_RD_b <- numerator/denominator
  
  # cumulative weights 
  numerator <- sum((PO_ctl_leq[j] * (1 - PO_ctl_leq[j])) * (PO_trt_leq[j] - PO_ctl_leq[j]))
  denominator <- sum((PO_ctl_leq[j] * (1 - PO_ctl_leq[j])))
  weighted_RD_cum <- numerator/denominator
  
  # both weights 
  numerator <- sum((p_y_ctl_est[j] + p_y_ctl_est[k]) * (PO_ctl_leq[j] * (1 - PO_ctl_leq[j])) * (PO_trt_leq[j] - PO_ctl_leq[j]))
  denominator <- sum((p_y_ctl_est[j] + p_y_ctl_est[k]) * (PO_ctl_leq[j] * (1 - PO_ctl_leq[j])))
  weighted_RD_both <- numerator/denominator 
  
  #calculated weighted RD with overall probability as weights
  # p_y <- (p_y_trt + p_y_ctl)/2
  # numerator <- sum((p_y_est[j] + p_y_est[k]) * (PO_trt_leq[j] - PO_ctl_leq[j]))
  # denominator <- sum(p_y_est[j]) + sum(p_y_est[k])
  # weighted_RD_boverall <- numerator/denominator
  
  # # calculate weighted RR with control probabilities 
  # numerator <- sum((p_y_ctl[j] + p_y_ctl[k]) * log_RR_le[j])
  # denominator <- sum(p_y_ctl[j]) + sum(p_y_ctl[k])
  # weighted_log_RR_control <- numerator/denominator
  # 
  # # calculate weighted RR with overall probabilities 
  # numerator <- sum((p_y[j] + p_y[k]) * log_RR_le[j]) 
  # denominator <- sum(p_y[j]) + sum(p_y[k])
  # weighted_log_RR_overall <- numerator/denominator
  
  # fully bayesian approach
  # calculate weighted RR with control probabilities 
  numerator <- sum((p_y_ctl_est[j] + p_y_ctl_est[k]) * log_RR_le[j])
  denominator <- sum(p_y_ctl_est[j]) + sum(p_y_ctl_est[k])
  weighted_log_RR_bcontrol <- numerator/denominator
  
  numerator <- sum((PO_ctl_leq[j] * (1 - PO_ctl_leq[j])) * log_RR_le[j])
  denominator <- sum((PO_ctl_leq[j] * (1 - PO_ctl_leq[j])))
  weighted_log_RR_cum <- numerator/denominator
  
  numerator <- sum((p_y_ctl_est[j] + p_y_ctl_est[k]) * (PO_ctl_leq[j] * (1 - PO_ctl_leq[j])) * log_RR_le[j])
  denominator <- sum((p_y_ctl_est[j] + p_y_ctl_est[k]) * (PO_ctl_leq[j] * (1 - PO_ctl_leq[j])))
  weighted_log_RR_both <- numerator/denominator
  
  numerator <- sum((p_y_ctl_est[j] + p_y_ctl_est[k]) * (PO_ctl_leq[j] * (1 - PO_ctl_leq[j])) * log_RR_ge[k])
  denominator <- sum((p_y_ctl_est[j] + p_y_ctl_est[k]) * (PO_ctl_leq[j] * (1 - PO_ctl_leq[j])))
  weighted_log_RR_ge_both <- numerator/denominator
  
  # # calculate weighted RR with overall probabilities 
  # numerator <- sum((p_y_est[j] + p_y_est[k]) * log_RR_le[j]) 
  # denominator <- sum(p_y_est[j]) + sum(p_y_est[k])
  # weighted_log_RR_boverall <- numerator/denominator
  
  
  # calculate NB 
  NB <- sum(p_y_trt[j] * PO_ctl_geq[j+1]) - (sum(p_y_ctl[j] * PO_trt_geq[j+1]))
  
  # weights 
  weights_split <- (p_y_ctl_est[j] + p_y_ctl_est[k]) / sum(p_y_ctl_est[j] + p_y_ctl_est[k])
  weights_cum <- (PO_ctl_leq[j] * (1 - PO_ctl_leq[j])) / sum(PO_ctl_leq[j] * (1 - PO_ctl_leq[j]))
  weights_both <- ((p_y_ctl_est[j] + p_y_ctl_est[k]) * (PO_ctl_leq[j] * (1 - PO_ctl_leq[j]))) /
    sum((p_y_ctl_est[j] + p_y_ctl_est[k]) * (PO_ctl_leq[j] * (1 - PO_ctl_leq[j])))
  
  return(list(#p_y_ctl = p_y_ctl, 
    # p_y_trt = p_y_trt,
    RD_leq = RD_leq,
    RD_geq = RD_geq,
    OR = OR[1],
    log_OR = log_OR[1],
    # weighted_RD_control = weighted_RD, 
    # weighted_RD_overall = weighted_RD_overall,
    # weighted_log_RR_control = weighted_log_RR_control,
    # weighted_log_RR_overall = weighted_log_RR_overall,
    weighted_RD_bcontrol = weighted_RD_b,
    weighted_RD_cum = weighted_RD_cum,
    weighted_RD_both = weighted_RD_both,
    # weighted_RD_boverall = weighted_RD_boverall,
    weighted_log_RR_bcontrol = weighted_log_RR_bcontrol,
    weighted_log_RR_cum = weighted_log_RR_cum,
    weighted_log_RR_both = weighted_log_RR_both,
    weighted_log_RR_ge_both = weighted_log_RR_ge_both,
    # weighted_log_RR_boverall = weighted_log_RR_boverall,
    unweighted_RD = unweighted_RD,
    unweighted_log_OR = unweighted_log_OR,
    unweighted_OR = unweighted_OR,
    unweighted_log_RR = unweighted_log_RR,
    unweighted_RR = unweighted_RR,
    unweighted_log_RR_ge = unweighted_log_RR_ge,
    unweighted_RR_ge = unweighted_RR_ge,
    #WR = WR,
    NB = NB, 
    weights_split = weights_split, 
    weights_cum = weights_cum,
    weights_both = weights_both))
  
}

# calculate the credible interval 
# assumes a partial proportional odds model 
# assumes the treatment is coded as -0.5 for control and 0.5 for treatment 
# b is the model 
order_weighted_RD2_PO_CI <- function(b, n_post_draws = 4000, dat, trt, outcome) { 
  # first calculate P(y >= i | ctl) and P(y >= i | trt) for all i
  posterior_draws <- as.matrix(b$draws)
  c <- length(b$ylevels)
  PO_ctl_geq <- matrix(data = NA, nrow = n_post_draws, ncol = c) 
  PO_trt_geq <- matrix(data = NA, nrow = n_post_draws, ncol = c) 
  PO_ctl_leq <- matrix(data = NA, nrow = n_post_draws, ncol = c) 
  PO_trt_leq <- matrix(data = NA, nrow = n_post_draws, ncol = c) 
  
  p_y_ctl_est <- matrix(data = NA, nrow = n_post_draws, ncol = c) 
  p_y_trt_est <- matrix(data = NA, nrow = n_post_draws, ncol = c) 
  
  for (i in 2:c) {
    # for the first level there is not a partial term, so only need first term and metformin term 
    PO_ctl_geq[,i] <- apply(posterior_draws, 1, function (b) plogis(b[i-1] - 0.5 * b[c]))
    PO_trt_geq[,i] <- apply(posterior_draws, 1, function (b) plogis(b[i-1] + 0.5 * b[c]))
    PO_ctl_leq[,(i - 1)] <- 1 - PO_ctl_geq[,i]
    PO_trt_leq[,(i - 1)] <- 1 - PO_trt_geq[,i]
  }
  
  # calculate p(Y = i | control) for i = 1, ... , c
  for (i in 1:c) {
    if (i == 1) {
      p_y_ctl_est[,1] <-  PO_ctl_leq[,1]
      p_y_trt_est[,1] <-  PO_trt_leq[,1]
    } else if (i == 2){
      p_y_ctl_est[,2] <- PO_ctl_leq[,i] - PO_ctl_leq[,1]
      p_y_trt_est[,2] <- PO_trt_leq[,i] - PO_trt_leq[,1]
    } else if (i == c) {
      p_y_ctl_est[,c] <- 1 - rowSums(p_y_ctl_est[,1:(c-1)])
      p_y_trt_est[,c] <- 1 - rowSums(p_y_trt_est[,1:(c-1)])
    } else {
      p_y_ctl_est[,i] <-  PO_ctl_leq[,i] - apply(p_y_ctl_est[,1:(i-1)], 1, sum)
      p_y_trt_est[,i] <-  PO_trt_leq[,i] - apply(p_y_trt_est[,1:(i-1)], 1, sum)
    }
  }
  
  p_y_est <- (p_y_ctl_est + p_y_trt_est)/2
  
  # use data to properly weight each levels 
  p_y_ctl <- unname(table(dat[[trt]], dat[[outcome]])[1,])/sum(dat[[trt]] == -0.5)
  p_y_trt <- unname(table(dat[[trt]], dat[[outcome]])[2,])/sum(dat[[trt]] == 0.5)
  p_y <- as.vector(unname(table(dat[[outcome]])))/nrow(dat)
  
  RD_leq <- PO_trt_leq - PO_ctl_leq
  RD_leq_95CI <- apply(RD_leq[,], 2, quantile, probs = c(0.025, 0.975), na.rm = T)
  unweighted_RD <- rowMeans(RD_leq[,1:(c-1)])
  unweighted_RD_95CI <- quantile(unweighted_RD, probs = c(0.025, 0.975))
  
  RD_geq <- PO_ctl_geq - PO_trt_geq
  RD_geq_95CI <- apply(RD_geq[,], 2, quantile, probs = c(0.025, 0.975), na.rm = T)
  
  # put control in the numerator so the OR is greater than one if metformin has a positive effect
  OR <- (PO_ctl_geq / (1- PO_ctl_geq))/(PO_trt_geq / (1 - PO_trt_geq))
  OR_95CI <- apply(OR[, 2:4], 2, quantile, probs = c(0.025, 0.975), na.rm = T)[1:2]
  OR_p <- 1 - sum(OR[,2] >= 1)/n_post_draws
  OR_sd <- sd(OR[,2])
  #(PPO_trt_leq / (1- PPO_trt_leq)) / (PPO_ctl_leq / (1- PPO_ctl_leq))
  log_OR <- log(OR)
  log_OR_95CI <- apply(log_OR[, 2:4], 2, quantile, probs = c(0.025, 0.975), na.rm = T)[1:2]
  log_OR_p <- 1 - sum(log_OR[,2] >= 0)/n_post_draws
  log_OR_sd <- sd(log_OR[,2])
  
  unweighted_log_OR <- rowMeans(log(OR[,2:c]))
  unweighted_log_OR_95CI <- quantile(unweighted_log_OR, probs = c(0.025, 0.975))
  unweighted_OR_95CI <- quantile(exp(unweighted_log_OR), probs = c(0.025, 0.975))
  
  log_RR_le <- log(PO_trt_leq / PO_ctl_leq)
  log_RR_le_95CI <- apply(log_RR_le[, 1:3], 2, quantile, probs = c(0.025, 0.975), na.rm = T)
  
  unweighted_log_RR <- rowMeans(log_RR_le[,1:(c-1)])
  unweighted_log_RR_95CI <- quantile(unweighted_log_RR, probs = c(0.025, 0.975))
  unweighted_RR_95CI <- quantile(exp(unweighted_log_RR), probs = c(0.025, 0.975))
  
  log_RR_ge <- log(PO_ctl_geq / PO_trt_geq)
  log_RR_ge_95CI <- apply(log_RR_ge[, 2:4], 2, quantile, probs = c(0.025, 0.975), na.rm = T)
  
  unweighted_log_RR_ge <- rowMeans(log_RR_ge[,2:c])
  unweighted_log_RR_ge_95CI <- quantile(unweighted_log_RR_ge, probs = c(0.025, 0.975))
  unweighted_RR_ge_95CI <- quantile(exp(unweighted_log_RR_ge), probs = c(0.025, 0.975))
  
  # calculate weighted RD with control probability as weights
  j <- 1:(c-1)
  k <- 2:c
  # numerator <- rowSums((p_y_ctl[j] + p_y_ctl[k]) * (PO_trt_leq[,j] - PO_ctl_leq[,j])) 
  # denominator <- sum(p_y_ctl[j]) + sum(p_y_ctl[k])
  # weighted_RD <- numerator/denominator
  # weighted_RD_95CI <- quantile(weighted_RD, probs = c(0.025, 0.975))
  # weighted_RD_p <- 1 - sum(weighted_RD >= 0)/n_post_draws
  # weighted_RD_sd <- sd(weighted_RD)
  
  #calculate weighted risk difference with overall probability as weights
  # numerator <- rowSums((p_y[j] + p_y[k]) * (PO_trt_leq[,j] - PO_ctl_leq[,j]))
  # denominator <- sum(p_y[j]) + sum(p_y[k])
  # weighted_RD_overall <- numerator/denominator
  # weighted_RD_95CI_overall <- quantile(weighted_RD_overall, probs = c(0.025, 0.975))
  # weighted_RD_p_overall <- 1 - sum(weighted_RD_overall >= 0)/n_post_draws
  # weighted_RD_sd_overall <- sd(weighted_RD_overall)
  
  # fully bayesian approach 
  numerator <- rowSums((p_y_ctl_est[,j] + p_y_ctl_est[,k]) * (PO_trt_leq[,j] - PO_ctl_leq[,j])) 
  denominator <- rowSums(p_y_ctl_est[,j]) + rowSums(p_y_ctl_est[,k])
  weighted_RD_b <- numerator/denominator
  weighted_RD_95CI_b <- quantile(weighted_RD_b, probs = c(0.025, 0.975))
  weighted_RD_p_b <- 1 - sum(weighted_RD_b >= 0)/n_post_draws
  weighted_RD_sd_b <- sd(weighted_RD_b)
  
  numerator <- rowSums((PO_ctl_leq[,j] * (1 - PO_ctl_leq[,j])) * (PO_trt_leq[,j] - PO_ctl_leq[,j])) 
  denominator <- rowSums((PO_ctl_leq[,j] * (1 - PO_ctl_leq[,j])))
  weighted_RD_cum <- numerator/denominator
  weighted_RD_95CI_cum <- quantile(weighted_RD_cum, probs = c(0.025, 0.975))
  weighted_RD_p_cum <- 1 - sum(weighted_RD_cum >= 0)/n_post_draws
  weighted_RD_sd_cum <- sd(weighted_RD_cum)
  
  numerator <- rowSums((p_y_ctl_est[,j] + p_y_ctl_est[,k]) * (PO_ctl_leq[,j] * (1 - PO_ctl_leq[,j])) * (PO_trt_leq[,j] - PO_ctl_leq[,j])) 
  denominator <- rowSums((p_y_ctl_est[,j] + p_y_ctl_est[,k]) * (PO_ctl_leq[,j] * (1 - PO_ctl_leq[,j])))
  weighted_RD_both <- numerator/denominator
  weighted_RD_95CI_both <- quantile(weighted_RD_both, probs = c(0.025, 0.975))
  weighted_RD_p_both <- 1 - sum(weighted_RD_both >= 0)/n_post_draws
  weighted_RD_sd_both <- sd(weighted_RD_both)
  
  #calculate weighted risk difference with overall probability as weights
  # numerator <- rowSums((p_y_est[,j] + p_y_est[,k]) * (PO_trt_leq[,j] - PO_ctl_leq[,j]))
  # denominator <- rowSums(p_y_est[,j]) + rowSums(p_y_est[,k])
  # weighted_RD_boverall <- numerator/denominator
  # weighted_RD_95CI_boverall <- quantile(weighted_RD_boverall, probs = c(0.025, 0.975))
  # weighted_RD_p_boverall <- 1 - sum(weighted_RD_boverall >= 0)/n_post_draws
  # weighted_RD_sd_boverall <- sd(weighted_RD_boverall)
  
  # calculate weighted RR with control probabilities 
  # numerator <- rowSums((p_y_ctl[j] + p_y_ctl[k]) * log_RR_le[,j]) 
  # denominator <- sum(p_y_ctl[j]) + sum(p_y_ctl[k])
  # weighted_log_RR_control <- numerator/denominator
  # weighted_log_RR_95CI_control <- quantile(weighted_log_RR_control, probs = c(0.025, 0.975))
  # weighted_log_RR_p_control <- 1 - sum(weighted_log_RR_control >= 0)/n_post_draws
  # weighted_log_RR_sd_control <- sd(weighted_log_RR_control)
  
  # calculate weighted RR with overall probabilities 
  # numerator <- rowSums((p_y[j] + p_y[k]) * log_RR_le[,j])
  # denominator <- sum(p_y[j]) + sum(p_y[k])
  # weighted_log_RR_overall <- numerator/denominator
  # weighted_log_RR_95CI_overall <- quantile(weighted_log_RR_overall, probs = c(0.025, 0.975))
  # weighted_log_RR_p_overall <- 1 - sum(weighted_log_RR_overall >= 0)/n_post_draws
  # weighted_log_RR_sd_overall <- sd(weighted_log_RR_overall)
  
  # fully bayesian approach
  # calculate weighted RR with control probabilities 
  numerator <- rowSums((p_y_ctl_est[,j] + p_y_ctl_est[,k]) * log_RR_le[,j]) 
  denominator <- rowSums(p_y_ctl_est[,j]) + rowSums(p_y_ctl_est[,k])
  weighted_log_RR_bcontrol <- numerator/denominator
  weighted_log_RR_95CI_bcontrol <- quantile(weighted_log_RR_bcontrol, probs = c(0.025, 0.975))
  weighted_log_RR_p_bcontrol <- 1 - sum(weighted_log_RR_bcontrol >= 0)/n_post_draws
  weighted_log_RR_sd_bcontrol <- sd(weighted_log_RR_bcontrol)
  
  numerator <- rowSums((PO_ctl_leq[,j] * (1 - PO_ctl_leq[,j])) * log_RR_le[,j]) 
  denominator <- rowSums((PO_ctl_leq[,j] * (1 - PO_ctl_leq[,j])))
  weighted_log_RR_cum <- numerator/denominator
  weighted_log_RR_95CI_cum <- quantile(weighted_log_RR_cum, probs = c(0.025, 0.975))
  weighted_log_RR_p_cum <- 1 - sum(weighted_log_RR_cum >= 0)/n_post_draws
  weighted_log_RR_sd_cum <- sd(weighted_log_RR_cum)
  
  numerator <- rowSums((p_y_ctl_est[,j] + p_y_ctl_est[,k]) * (PO_ctl_leq[,j] * (1 - PO_ctl_leq[,j])) * log_RR_le[,j]) 
  denominator <- rowSums((p_y_ctl_est[,j] + p_y_ctl_est[,k]) * (PO_ctl_leq[,j] * (1 - PO_ctl_leq[,j])))
  weighted_log_RR_both <- numerator/denominator
  weighted_log_RR_95CI_both <- quantile(weighted_log_RR_both, probs = c(0.025, 0.975))
  weighted_log_RR_p_both <- 1 - sum(weighted_log_RR_both >= 0)/n_post_draws
  weighted_log_RR_sd_both <- sd(weighted_log_RR_both)
  
  numerator <- rowSums((p_y_ctl_est[,j] + p_y_ctl_est[,k]) * (PO_ctl_leq[,j] * (1 - PO_ctl_leq[,j])) * log_RR_ge[,k]) 
  denominator <- rowSums((p_y_ctl_est[,j] + p_y_ctl_est[,k]) * (PO_ctl_leq[,j] * (1 - PO_ctl_leq[,j])))
  weighted_log_RR_ge_both <- numerator/denominator
  weighted_log_RR_ge_95CI_both <- quantile(weighted_log_RR_ge_both, probs = c(0.025, 0.975))
  weighted_log_RR_ge_p_both <- 1 - sum(weighted_log_RR_ge_both >= 0)/n_post_draws
  weighted_log_RR_ge_sd_both <- sd(weighted_log_RR_ge_both)
  
  # calculate weighted RR with overall probabilities 
  # numerator <- rowSums((p_y_est[,j] + p_y_est[,k]) * log_RR_le[,j])
  # denominator <- rowSums(p_y_est[,j]) + rowSums(p_y_est[,k])
  # weighted_log_RR_boverall <- numerator/denominator
  # weighted_log_RR_95CI_boverall <- quantile(weighted_log_RR_boverall, probs = c(0.025, 0.975))
  # weighted_log_RR_p_boverall <- 1 - sum(weighted_log_RR_boverall >= 0)/n_post_draws
  # weighted_log_RR_sd_boverall <- sd(weighted_log_RR_boverall)
  
  
  # calculate NB 
  NB <- rowSums(p_y_trt[j] * PO_ctl_geq[,(j+1)]) - (rowSums(p_y_ctl[j] * PO_trt_geq[,(j+1)]))
  NB_95CI <- quantile(NB, probs = c(0.025, 0.975))
  NB_p <- 1 - sum(NB >= 0)/n_post_draws 
  NB_sd <- sd(NB)
  
  #weights 
  weights_split <- (p_y_ctl_est[,j] + p_y_ctl_est[,k]) / rowSums((p_y_ctl_est[,j] + p_y_ctl_est[,k]))
  weights_split_95CI <- apply(weights_split, 2, quantile, probs = c(0.025, 0.975))
  weights_cum <- (PO_ctl_leq[,j] * (1 - PO_ctl_leq[,j]))  / rowSums((PO_ctl_leq[,j] * (1 - PO_ctl_leq[,j])) )
  weights_cum_95CI <- apply(weights_cum, 2, quantile, probs = c(0.025, 0.975))
  weights_both <- (p_y_ctl_est[,j] + p_y_ctl_est[,k]) * (PO_ctl_leq[,j] * (1 - PO_ctl_leq[,j])) / 
    rowSums((p_y_ctl_est[,j] + p_y_ctl_est[,k]) * (PO_ctl_leq[,j] * (1 - PO_ctl_leq[,j])))
  weights_both_95CI <- apply(weights_both, 2, quantile, probs = c(0.025, 0.975))
  
  return(list(RD_leq = RD_leq_95CI, 
              RD_geq = RD_geq_95CI, 
              OR_CI = OR_95CI, 
              OR_p = OR_p,
              OR_sd = OR_sd,
              log_OR_CI = log_OR_95CI, 
              log_OR_p = log_OR_p,
              log_OR_sd = log_OR_sd,
              log_RR_le_95CI = log_RR_le_95CI,
              log_RR_ge_95CI = log_RR_ge_95CI,
              # weighted_RD_control_CI = weighted_RD_95CI, 
              # weighted_RD_p_control = weighted_RD_p,
              # weighted_RD_sd = weighted_RD_sd,
              # weighted_RD_overall_CI = weighted_RD_95CI_overall, 
              # weighted_RD_p_overall = weighted_RD_p_overall,
              # weighted_RD_sd_overall = weighted_RD_sd_overall,
              # weighted_log_RR_95CI_control = weighted_log_RR_95CI_control,
              # weighted_log_RR_p_control = weighted_log_RR_p_control, 
              # weighted_log_RR_sd_control = weighted_log_RR_sd_control,
              # weighted_log_RR_95CI_overall = weighted_log_RR_95CI_overall,
              # weighted_log_RR_p_overall = weighted_log_RR_p_overall, 
              # weighted_log_RR_sd_overall = weighted_log_RR_sd_overall,
              weighted_RD_bcontrol_CI = weighted_RD_95CI_b, 
              weighted_RD_p_bcontrol = weighted_RD_p_b,
              weighted_RD_sd_b = weighted_RD_sd_b,
              weighted_RD_cum_CI = weighted_RD_95CI_cum, 
              weighted_RD_p_cum = weighted_RD_p_cum,
              weighted_RD_sd_cum = weighted_RD_sd_cum,
              weighted_RD_both_CI = weighted_RD_95CI_both, 
              weighted_RD_p_both = weighted_RD_p_both,
              weighted_RD_sd_both = weighted_RD_sd_both,
              # weighted_RD_boverall_CI = weighted_RD_95CI_boverall, 
              # weighted_RD_p_boverall = weighted_RD_p_boverall,
              # weighted_RD_sd_boverall = weighted_RD_sd_boverall,
              weighted_log_RR_95CI_bcontrol = weighted_log_RR_95CI_bcontrol,
              weighted_log_RR_p_bcontrol = weighted_log_RR_p_bcontrol, 
              weighted_log_RR_sd_bcontrol = weighted_log_RR_sd_bcontrol,
              weighted_log_RR_95CI_cum = weighted_log_RR_95CI_cum,
              weighted_log_RR_p_cum = weighted_log_RR_p_cum, 
              weighted_log_RR_sd_cum = weighted_log_RR_sd_cum,
              weighted_log_RR_95CI_both = weighted_log_RR_95CI_both,
              weighted_log_RR_p_both = weighted_log_RR_p_both, 
              weighted_log_RR_sd_both = weighted_log_RR_sd_both,
              weighted_log_RR_ge_95CI_both = weighted_log_RR_ge_95CI_both,
              weighted_log_RR_ge_p_both = weighted_log_RR_ge_p_both, 
              weighted_log_RR_ge_sd_both = weighted_log_RR_ge_sd_both,
              # weighted_log_RR_95CI_boverall = weighted_log_RR_95CI_boverall,
              # weighted_log_RR_p_boverall = weighted_log_RR_p_boverall, 
              # weighted_log_RR_sd_boverall = weighted_log_RR_sd_boverall,
              unweighted_RD_95CI = unweighted_RD_95CI,
              unweighted_log_OR_95CI = unweighted_log_OR_95CI,
              unweighted_OR_95CI = unweighted_OR_95CI, 
              unweighted_log_RR_95CI = unweighted_log_RR_95CI,
              unweighted_RR_95CI = unweighted_RR_95CI,
              unweighted_log_RR_ge_95CI = unweighted_log_RR_ge_95CI,
              unweighted_RR_ge_95CI = unweighted_RR_ge_95CI,
              #WR = WR_95CI, 
              NB_CI = NB_95CI, 
              NB_p = NB_p, 
              NB_sd = NB_sd,
              weights_split_95CI = weights_split_95CI,
              weights_cum_95CI = weights_cum_95CI,
              weights_both_95CI = weights_both_95CI))
  
}


# proportional odds model adjusting for other covariates
# use bayesian bootstrap to combine over individuals
order_weighted_RD2_PO_adj_bb <- function(b = b, c, dat, trt_variable,
                                         outcome) { 
  
  dat_trt <- dat 
  dat_trt[[trt_variable]] <- 0.5
  
  dat_ctl <- dat 
  dat_ctl[[trt_variable]] <- -0.5
  
  predicted_trt_geq <- matrix(data = NA, nrow = nrow(dat), ncol = c-1) 
  predicted_ctl_geq <- matrix(data = NA, nrow = nrow(dat), ncol = c-1) 
  predicted_trt_leq <- matrix(data = NA, nrow = nrow(dat), ncol = c-1) 
  predicted_ctl_leq <- matrix(data = NA, nrow = nrow(dat), ncol = c-1) 
  
  p_y_ctl_est <- matrix(data = NA, nrow = nrow(dat), ncol = c)
  p_y_trt_est <- matrix(data = NA, nrow = nrow(dat), ncol = c)
  
  #predicted_trt_geq and predicted_ctl_geq are in the format where the first column is the P(y >= 1 |trt, covariates)
  #predicted_trt_leq and predicted_ctl_leq are in the format where the first column is the P(y <= 1 |trt, covariates)
  for (i in 1:(c-1)){
    predicted_trt_geq[,i] <-  plogis(predict(b, dat_trt, kint = i)$linear.predictors)
    predicted_ctl_geq[,i] <-  plogis(predict(b, dat_ctl, kint = i)$linear.predictors)
  }
  
  predicted_trt_leq <-  1 - predicted_trt_geq
  predicted_ctl_leq <-  1 - predicted_ctl_geq
  
  # calculate p(Y = i | control) for i = 1, ... , c
  for (i in 1:c) {
    if (i == 1) {
      p_y_ctl_est[,1] <-  predicted_ctl_leq[,1]
      p_y_trt_est[,1] <-  predicted_trt_leq[,1]
    } else if (i == 2) {
      p_y_ctl_est[,i] <-  predicted_ctl_leq[,i] - p_y_ctl_est[,1]
      p_y_trt_est[,i] <-  predicted_trt_leq[,i] - p_y_trt_est[,1]
    } else if (i == c) {
      p_y_ctl_est[,c] <- 1 - rowSums(p_y_ctl_est[,1:(c-1)])
      p_y_trt_est[,c] <- 1 - rowSums(p_y_trt_est[,1:(c-1)])
    } else {
      p_y_ctl_est[,i] <-  predicted_ctl_leq[,i] - rowSums(p_y_ctl_est[,1:(i-1)])
      p_y_trt_est[,i] <-  predicted_trt_leq[,i] - rowSums(p_y_trt_est[,1:(i-1)])
    }
  }
  
  p_y_est <- (p_y_ctl_est + p_y_trt_est)/2
  
  
  # get weights from the data 
  p_y_ctl <- unname(table(dat[[trt_variable]], dat[[outcome]])[1,])/sum(dat[[trt_variable]] == -0.5)
  p_y_trt <- unname(table(dat[[trt_variable]], dat[[outcome]])[2,])/sum(dat[[trt_variable]] == 0.5)
  p_y <- unname(table(dat[[outcome]]))/nrow(dat)
  
  weights <- t(rdirichlet(1, rep(1, nrow(dat))))
  
  p_y_ctl_leq_average <- apply(predicted_ctl_leq, 2, weighted.mean, w = weights)
  p_y_trt_leq_average <- apply(predicted_trt_leq, 2, weighted.mean, w = weights)
  p_y_ctl_geq_average <- apply(predicted_ctl_geq, 2, weighted.mean, w = weights)
  p_y_trt_geq_average <- apply(predicted_trt_geq, 2, weighted.mean, w = weights)
  
  p_y_ctl_est_average <- apply(p_y_ctl_est, 2, weighted.mean, w = weights)
  p_y_est_average <- apply(p_y_est, 2, weighted.mean, w = weights)
  
  RD_leq <- p_y_trt_leq_average - p_y_ctl_leq_average
  unweighted_RD <- mean(RD_leq)
  #unweighted_RD <- weighted.mean(unweighted_RD_all, w = weights)
  RD_geq <- p_y_ctl_geq_average - p_y_trt_geq_average
  
  
  OR <- (p_y_trt_leq_average / (1- p_y_trt_leq_average)) / (p_y_ctl_leq_average / (1- p_y_ctl_leq_average))
  unweighted_log_OR <- mean(log(OR))
  unweighted_OR <- exp(mean(log(OR)))
  
  RR_leq <- p_y_trt_leq_average / p_y_ctl_leq_average
  unweighted_log_RR <- mean(log(RR_leq))
  unweighted_RR <- exp(mean(log(RR_leq)))
  RR_geq <- p_y_ctl_geq_average / p_y_trt_geq_average
  unweighted_log_RR_ge <- mean(log(RR_geq))
  unweighted_RR_ge <- exp(mean(log(RR_geq)))
  
  j <- 1:(c-1)
  k <- 2:c
  # numerator <- sum((p_y_ctl[j] + p_y_ctl[k])  * (p_y_trt_leq_average[j] - p_y_ctl_leq_average[j]))
  # denominator <- sum(p_y_ctl[j]) + sum(p_y_ctl[k])
  # weighted_RD <- numerator/denominator
  
  
  #calculated weighted RD with overall probability as weights
  # numerator <- sum((p_y[j] + p_y[k]) * (p_y_trt_leq_average[j] - p_y_ctl_leq_average[j])) 
  # denominator <- sum(p_y[j]) + sum(p_y[k])
  # weighted_RD_overall <- numerator/denominator
  
  # fully bayesian approach
  numerator <- sum((p_y_ctl_est_average[j] + p_y_ctl_est_average[k]) * (p_y_trt_leq_average[j] - p_y_ctl_leq_average[j]))
  denominator <- sum(p_y_ctl_est_average[j]) + sum(p_y_ctl_est_average[k])
  weighted_RD_b <- numerator/denominator
  
  numerator <- sum((p_y_ctl_leq_average[j] * (1 - p_y_ctl_leq_average[j])) * (p_y_trt_leq_average[j] - p_y_ctl_leq_average[j]))
  denominator <- sum((p_y_ctl_leq_average[j] * (1 - p_y_ctl_leq_average[j])))
  weighted_RD_cum <- numerator/denominator
  
  numerator <- sum((p_y_ctl_est_average[j] + p_y_ctl_est_average[k]) * (p_y_ctl_leq_average[j] * (1 - p_y_ctl_leq_average[j])) * (p_y_trt_leq_average[j] - p_y_ctl_leq_average[j]))
  denominator <- sum((p_y_ctl_est_average[j] + p_y_ctl_est_average[k]) * (p_y_ctl_leq_average[j] * (1 - p_y_ctl_leq_average[j])))
  weighted_RD_both <- numerator/denominator 
  
  #calculated weighted RD with overall probability as weights
  # numerator <- sum((p_y_est_average[j] + p_y_est_average[k]) * (p_y_trt_leq_average[j] - p_y_ctl_leq_average[j])) 
  # denominator <- sum(p_y_est_average[j]) + sum(p_y_est_average[k])
  # weighted_RD_boverall <- numerator/denominator
  
  #calculate WR    
  l <- 1:c
  WP <- sum(p_y_trt[j] * p_y_ctl_geq_average[j]) + 0.5* sum(p_y_trt[l] * p_y_ctl[l])
  WR <- WP/(1-WP)
  
  
  
  #calculate weighted OR with control probability as weights
  # numerator <- sum((p_y_ctl[j] + p_y_ctl[k])  * log(OR[j]))
  # denominator <- sum(p_y_ctl[j]) + sum(p_y_ctl[k])
  # log_weighted_OR <- numerator/denominator
  # weighted_OR <- exp(numerator/denominator)
  # 
  # #calculate weighted OR with overall probability as weights
  # numerator <- sum((p_y[j] + p_y[k])  * log(OR[j])) 
  # denominator <- sum(p_y[j]) + sum(p_y[k])
  # log_weighted_OR_overall <- numerator/denominator
  # weighted_OR_overall <- exp(numerator/denominator)
  
  #fully bayesian approach
  #calculate weighted OR with control probability as weights
  numerator <- sum((p_y_ctl_est_average[j] + p_y_ctl_est_average[k]) * log(OR[j]))
  denominator <- sum(p_y_ctl_est_average[j]) + sum(p_y_ctl_est_average[k])
  log_weighted_OR_b <- numerator/denominator
  weighted_OR_b <- exp(numerator/denominator)
  
  # cumulative weights
  numerator <- sum((p_y_ctl_leq_average[j] * (1 - p_y_ctl_leq_average[j])) * log(OR[j]))
  denominator <- sum((p_y_ctl_leq_average[j] * (1 - p_y_ctl_leq_average[j])))
  log_weighted_OR_cum <- numerator/denominator
  weighted_OR_cum <- exp(numerator/denominator)
  
  # both weights
  numerator <- sum((p_y_ctl_est_average[j] + p_y_ctl_est_average[k]) * (p_y_ctl_leq_average[j] * (1 - p_y_ctl_leq_average[j])) * log(OR[j]))
  denominator <- sum((p_y_ctl_est_average[j] + p_y_ctl_est_average[k]) * (p_y_ctl_leq_average[j] * (1 - p_y_ctl_leq_average[j])))
  log_weighted_OR_both <- numerator/denominator
  weighted_OR_both <- exp(numerator/denominator)
  
  
  # #calculate weighted OR with overall probability as weights
  # numerator <- sum((p_y_est_average[j] + p_y_est_average[k]) * log(OR[j])) 
  # denominator <- sum(p_y_est_average[j]) + sum(p_y_est_average[k])
  # log_weighted_OR_boverall <- numerator/denominator
  # weighted_OR_boverall <- exp(numerator/denominator)
  
  
  # #calculate weighted RR with control probability as weights
  # numerator <- sum((p_y_ctl[j] + p_y_ctl[k]) * log(RR_leq[j]))
  # denominator <- sum(p_y_ctl[j]) + sum(p_y_ctl[k])
  # log_weighted_RR <- numerator/denominator
  # weighted_RR <- exp(numerator/denominator)
  
  #calculate weighted RR with overall probability as weights
  # numerator <- sum((p_y[j] + p_y[k])  * log(RR_leq[j]))
  # denominator <- sum(p_y[j]) + sum(p_y[k])
  # log_weighted_RR_overall <- numerator/denominator
  # weighted_RR_overall <- exp(numerator/denominator)
  
  # fully bayesian approach
  #calculate weighted RR with control probability as weights
  numerator <- sum((p_y_ctl_est_average[j] + p_y_ctl_est_average[k]) * log(RR_leq[j]))
  denominator <- sum(p_y_ctl_est_average[j]) + sum(p_y_ctl_est_average[k])
  log_weighted_RR_b <- numerator/denominator
  weighted_RR_b <- exp(numerator/denominator)
  
  # cumulative weights 
  numerator <- sum((p_y_ctl_leq_average[j] * (1 - p_y_ctl_leq_average[j])) * log(RR_leq[j]))
  denominator <- sum((p_y_ctl_leq_average[j] * (1 - p_y_ctl_leq_average[j])))
  log_weighted_RR_cum <- numerator/denominator
  weighted_RR_cum <- exp(numerator/denominator)
  
  # both weights
  numerator <- sum((p_y_ctl_est_average[j] + p_y_ctl_est_average[k]) * (p_y_ctl_leq_average[j] * (1 - p_y_ctl_leq_average[j])) * log(RR_leq[j]))
  denominator <- sum((p_y_ctl_est_average[j] + p_y_ctl_est_average[k]) * (p_y_ctl_leq_average[j] * (1 - p_y_ctl_leq_average[j])))
  log_weighted_RR_both <- numerator/denominator
  weighted_RR_both <- exp(numerator/denominator)
  
  # both weights RR greater than or equal to 
  numerator <- sum((p_y_ctl_est_average[j] + p_y_ctl_est_average[k]) * (p_y_ctl_leq_average[j] * (1 - p_y_ctl_leq_average[j])) * log(RR_geq[j]))
  denominator <- sum((p_y_ctl_est_average[j] + p_y_ctl_est_average[k]) * (p_y_ctl_leq_average[j] * (1 - p_y_ctl_leq_average[j])))
  log_weighted_RR_ge_both <- numerator/denominator
  weighted_RR_ge_both <- exp(numerator/denominator)
  
  #calculate weighted RR with overall probability as weights
  # numerator <- sum((p_y_est[j] + p_y_est[k]) * log(RR_leq[j]))
  # denominator <- sum(p_y_est[j]) + sum(p_y_est[k])
  # log_weighted_RR_boverall <- numerator/denominator
  # weighted_RR_boverall <- exp(numerator/denominator)
  
  #calculate WR    
  l <- 1:c
  NB <- sum(p_y_trt[j] * p_y_ctl_geq_average[j]) - (sum(p_y_ctl[j] * p_y_trt_geq_average[j]))
  
  #weights 
  weights_split <- (p_y_ctl_est_average[j] + p_y_ctl_est_average[k]) / sum((p_y_ctl_est_average[j] + p_y_ctl_est_average[k]))
  weights_cum <- (p_y_ctl_leq_average[j] * (1 - p_y_ctl_leq_average[j])) / sum((p_y_ctl_leq_average[j] * (1 - p_y_ctl_leq_average[j])))
  weights_both <- ((p_y_ctl_est_average[j] + p_y_ctl_est_average[k]) * (p_y_ctl_leq_average[j] * (1 - p_y_ctl_leq_average[j]))) / 
    sum((p_y_ctl_est_average[j] + p_y_ctl_est_average[k]) * (p_y_ctl_leq_average[j] * (1 - p_y_ctl_leq_average[j])))
  
  return(list(p_y_ctl = p_y_ctl, 
              p_y_trt = p_y_trt,
              RD_leq = RD_leq,
              RD_geq = RD_geq,
              OR = OR,
              RR_leq = RR_leq,
              RR_geq = RR_geq,
              log_OR = log(OR),
              log_RR_leq = log(RR_leq),
              # weighted_RD = weighted_RD, 
              # weighted_RD_overall = weighted_RD_overall,
              # weighted_OR = weighted_OR,
              # log_weighted_OR = log_weighted_OR,
              # weighted_OR_overall =  weighted_OR_overall,
              # log_weighted_OR_overall = log_weighted_OR_overall,
              # weighted_RR = weighted_RR,
              # log_weighted_RR = log_weighted_RR,
              # weighted_RR_overall =  weighted_RR_overall,
              # log_weighted_RR_overall = log_weighted_RR_overall,
              weighted_RD_b = weighted_RD_b, 
              weighted_RD_cum = weighted_RD_cum, 
              weighted_RD_both = weighted_RD_both, 
              # weighted_RD_boverall = weighted_RD_boverall,
              weighted_OR_b = weighted_OR_b,
              log_weighted_OR_b = log_weighted_OR_b,
              weighted_OR_cum = weighted_OR_cum,
              log_weighted_OR_cum = log_weighted_OR_cum,
              weighted_OR_both = weighted_OR_both,
              log_weighted_OR_both = log_weighted_OR_both,
              # weighted_OR_boverall =  weighted_OR_boverall,
              # log_weighted_OR_boverall = log_weighted_OR_boverall,
              weighted_RR_b = weighted_RR_b,
              log_weighted_RR_b = log_weighted_RR_b,
              weighted_RR_cum = weighted_RR_cum,
              log_weighted_RR_cum = log_weighted_RR_cum,
              weighted_RR_both = weighted_RR_both,
              log_weighted_RR_both = log_weighted_RR_both,
              weighted_RR_ge_both = weighted_RR_ge_both,
              log_weighted_RR_ge_both = log_weighted_RR_ge_both,
              # weighted_RR_boverall =  weighted_RR_boverall,
              # log_weighted_RR_boverall = log_weighted_RR_boverall,
              unweighted_RD = unweighted_RD,
              unweighted_log_OR = unweighted_log_OR,
              unweighted_OR = unweighted_OR,
              unweighted_log_RR = unweighted_log_RR,
              unweighted_RR = unweighted_RR,
              unweighted_log_RR_ge = unweighted_log_RR_ge,
              unweighted_RR_ge = unweighted_RR_ge,
              NB = NB,
              WR = WR,
              weights_split = weights_split,
              weights_cum = weights_cum,
              weights_both = weights_both))
  
}

# new function for adjusted CI for proportional odds model 
# using bayesian bootstrap to combine over individuals
order_weighted_RD2_CI_PO_adj_bb <- function(b, c, n_post_draws = 4000, dat, iprior,
                                            outcome, trt, covs = 1) { 
  
  # first calculate P(y >= i | ctl) and P(y >= i | trt) for all i
  posterior_draws <- as.matrix(b$draws)
  PPO_ctl_geq_averages <- matrix(data = NA, nrow = n_post_draws, ncol = c) 
  PPO_trt_geq_averages <- matrix(data = NA, nrow = n_post_draws, ncol = c) 
  PPO_ctl_leq_averages <- matrix(data = NA, nrow = n_post_draws, ncol = c) 
  PPO_trt_leq_averages <- matrix(data = NA, nrow = n_post_draws, ncol = c) 
  
  p_y_ctl_est_averages <- matrix(data = NA, nrow = n_post_draws, ncol = c)
  p_y_trt_est_averages <- matrix(data = NA, nrow = n_post_draws, ncol = c)
  p_y_est_averages <- matrix(data = NA, nrow = n_post_draws, ncol = c)
  
  wRD_control <- matrix(data = NA, nrow = n_post_draws, ncol = 1)
  wRD_overall <- matrix(data = NA, nrow = n_post_draws, ncol = 1)
  wOR_control <- matrix(data = NA, nrow = n_post_draws, ncol = 1)
  wOR_overall <- matrix(data = NA, nrow = n_post_draws, ncol = 1)
  log_wOR_control <- matrix(data = NA, nrow = n_post_draws, ncol = 1)
  log_wOR_overall <- matrix(data = NA, nrow = n_post_draws, ncol = 1)
  wRR_control <- matrix(data = NA, nrow = n_post_draws, ncol = 1)
  wRR_overall <- matrix(data = NA, nrow = n_post_draws, ncol = 1)
  log_wRR_control <- matrix(data = NA, nrow = n_post_draws, ncol = 1)
  log_wRR_overall <- matrix(data = NA, nrow = n_post_draws, ncol = 1)
  
  wRD_bcontrol <- matrix(data = NA, nrow = n_post_draws, ncol = 1)
  wRD_boverall <- matrix(data = NA, nrow = n_post_draws, ncol = 1)
  wOR_bcontrol <- matrix(data = NA, nrow = n_post_draws, ncol = 1)
  wOR_boverall <- matrix(data = NA, nrow = n_post_draws, ncol = 1)
  log_wOR_bcontrol <- matrix(data = NA, nrow = n_post_draws, ncol = 1)
  log_wOR_boverall <- matrix(data = NA, nrow = n_post_draws, ncol = 1)
  wRR_bcontrol <- matrix(data = NA, nrow = n_post_draws, ncol = 1)
  wRR_boverall <- matrix(data = NA, nrow = n_post_draws, ncol = 1)
  log_wRR_bcontrol <- matrix(data = NA, nrow = n_post_draws, ncol = 1)
  log_wRR_boverall <- matrix(data = NA, nrow = n_post_draws, ncol = 1)
  
  wRD_cum <- matrix(data = NA, nrow = n_post_draws, ncol = 1)
  wRD_both <- matrix(data = NA, nrow = n_post_draws, ncol = 1)
  wOR_cum <- matrix(data = NA, nrow = n_post_draws, ncol = 1)
  log_wOR_cum <- matrix(data = NA, nrow = n_post_draws, ncol = 1)
  wOR_both <- matrix(data = NA, nrow = n_post_draws, ncol = 1)
  log_wOR_both <- matrix(data = NA, nrow = n_post_draws, ncol = 1)
  wRR_cum <- matrix(data = NA, nrow = n_post_draws, ncol = 1)
  log_wRR_cum <- matrix(data = NA, nrow = n_post_draws, ncol = 1)
  wRR_both <- matrix(data = NA, nrow = n_post_draws, ncol = 1)
  log_wRR_both <- matrix(data = NA, nrow = n_post_draws, ncol = 1)
  wRR_ge_both <- matrix(data = NA, nrow = n_post_draws, ncol = 1)
  log_wRR_ge_both <- matrix(data = NA, nrow = n_post_draws, ncol = 1)
  
  weights_split <- matrix(data = NA, nrow = n_post_draws, ncol = c-1)
  weights_cum <- matrix(data = NA, nrow = n_post_draws, ncol = c-1)
  weights_both <- matrix(data = NA, nrow = n_post_draws, ncol = c-1)
  
  NB <- matrix(data = NA, nrow = n_post_draws, ncol = 1)
  RD_leq_save <-  matrix(data = NA, nrow = n_post_draws, ncol = c-1)
  OR_save <-  matrix(data = NA, nrow = n_post_draws, ncol = c-1)
  RR_leq_save <-  matrix(data = NA, nrow = n_post_draws, ncol = c-1)
  RR_geq_save <-  matrix(data = NA, nrow = n_post_draws, ncol = c-1)
  
  log_OR_save <-  matrix(data = NA, nrow = n_post_draws, ncol = c-1)
  log_RR_leq_save <-  matrix(data = NA, nrow = n_post_draws, ncol = c-1)
  log_RR_geq_save <-  matrix(data = NA, nrow = n_post_draws, ncol = c-1)
  
  foreach (s = 1:n_post_draws) %do% {
    
    predicted_trt_geq <- matrix(data = NA, nrow = nrow(dat), ncol = c-1) 
    predicted_ctl_geq <- matrix(data = NA, nrow = nrow(dat), ncol = c-1) 
    predicted_trt_leq <- matrix(data = NA, nrow = nrow(dat), ncol = c-1) 
    predicted_ctl_leq <- matrix(data = NA, nrow = nrow(dat), ncol = c-1) 
    
    # if (iprior == 0){
    #   coefs <- posterior_draws[i, ]
    # } else {
    #   coefs <- c(-coef(b)[1:3], coef(b)[4:9])
    # }
    
    coefs <- posterior_draws[s, ]
    covs <- b$Design$colnames
    num_covs <- length(covs)
    covariate_mat <- as.data.frame(matrix(data = NA, nrow = nrow(dat), ncol = num_covs - 1))
    
    # assume the first covariate is the treatment arm 
    for(i in 2:num_covs) {
      if (is.numeric(dat[[b$Design$name[i]]]) == T) {
        covariate_mat[,i-1] <- coefs[[covs[i]]] * (dat[[b$Design$name[i]]])
      } else {
        # only works for numeric and binary
        covariate_mat[,i-1] <- coefs[[covs[i]]] * (dat[[b$Design$name[i]]] == b$Design$parms[[b$Design$name[i]]][2])
      }
    }
    
    for (i in 2:c) {
      if (i == 2) {
        predicted_trt_geq[,i-1] <- plogis(coefs[i-1] + coefs[[trt]] * 0.5 + rowSums(covariate_mat))
        predicted_ctl_geq[,i-1] <- plogis(coefs[i-1] + coefs[[trt]] * -0.5 + rowSums(covariate_mat))
      } else {
        phrase <- paste0("y>=", i)
        coefs_new <- coefs[str_detect(names(coefs), phrase)]
        predicted_trt_geq[,i-1] <- plogis(coefs_new[1] + coefs[[trt]] * 0.5 + rowSums(covariate_mat))
        predicted_ctl_geq[,i-1] <- plogis(coefs_new[1] + coefs[[trt]] * -0.5 + rowSums(covariate_mat))
      }
    }
    
    predicted_trt_leq <-  1 - predicted_trt_geq
    predicted_ctl_leq <-  1 - predicted_ctl_geq
    
    p_y_trt_est <- matrix(data = NA, nrow = nrow(dat), ncol = c) 
    p_y_ctl_est <- matrix(data = NA, nrow = nrow(dat), ncol = c) 
    for (i in 1:c) {
      if (i == 1) {
        p_y_ctl_est[,1] <-  predicted_ctl_leq[,1]
        p_y_trt_est[,1] <-  predicted_trt_leq[,1]
      } else if (i == 2){
        p_y_ctl_est[,2] <- predicted_ctl_leq[,i] - predicted_ctl_leq[,1]
        p_y_trt_est[,2] <- predicted_trt_leq[,i] - predicted_trt_leq[,1]
      } else if (i == c) {
        p_y_ctl_est[,c] <- 1 - rowSums(p_y_ctl_est[,1:(c-1)])
        p_y_trt_est[,c] <- 1 - rowSums(p_y_trt_est[,1:(c-1)])
      } else {
        p_y_ctl_est[,i] <-  predicted_ctl_leq[,i] - apply(p_y_ctl_est[,1:(i-1)], 1, sum)
        p_y_trt_est[,i] <-  predicted_trt_leq[,i] - apply(p_y_trt_est[,1:(i-1)], 1, sum)
      }
    }
    
    
    p_y_est <- (p_y_ctl_est + p_y_trt_est)/2
    
    
    # get weights from the data 
    p_y_ctl <- unname(table(dat[[trt]], dat[[outcome]])[1,])/sum(dat[[trt]] == -0.5)
    p_y_trt <- unname(table(dat[[trt]], dat[[outcome]])[2,])/sum(dat[[trt]] == 0.5)
    p_y <- as.vector(unname(table(dat[[outcome]])))/nrow(dat)
    
    weights <- t(rdirichlet(1, rep(1, nrow(dat))))
    
    p_y_ctl_leq_average <- apply(predicted_ctl_leq, 2, weighted.mean, w = weights)
    p_y_trt_leq_average <- apply(predicted_trt_leq, 2, weighted.mean, w = weights)
    p_y_ctl_geq_average <- apply(predicted_ctl_geq, 2, weighted.mean, w = weights)
    p_y_trt_geq_average <- apply(predicted_trt_geq, 2, weighted.mean, w = weights)
    
    p_y_ctl_est_average <- apply(p_y_ctl_est, 2, weighted.mean, w = weights)
    p_y_est_average <- apply(p_y_est, 2, weighted.mean, w = weights)
    
    
    RD_leq <- p_y_trt_leq_average - p_y_ctl_leq_average
    RD_leq_save[s,] <- RD_leq
    RD_geq <- p_y_ctl_geq_average - p_y_trt_geq_average
    
    
    OR <- (p_y_trt_leq_average / (1- p_y_trt_leq_average)) / (p_y_ctl_leq_average / (1- p_y_ctl_leq_average))
    OR_save[s,] <- OR
    log_OR_save[s,] <- log(OR)
    
    RR_leq <- p_y_trt_leq_average / p_y_ctl_leq_average
    RR_leq_save[s,] <- RR_leq
    log_RR_leq_save[s,] <- log(RR_leq)
    RR_geq <- p_y_ctl_geq_average / p_y_trt_geq_average
    RR_geq_save[s,] <- RR_geq
    log_RR_geq_save[s,] <- log(RR_geq)
    
    j <- 1:(c-1)
    k <- 2:c
    # numerator <- sum((p_y_ctl[j] + p_y_ctl[k]) * (p_y_trt_leq_average[j] - p_y_ctl_leq_average[j]))
    # denominator <- sum(p_y_ctl[j]) + sum(p_y_ctl[k])
    # wRD_control[s,1] <- numerator/denominator
    
    #calculated weighted RD with overall probability as weights
    # numerator <- sum((p_y[j] + p_y[k]) * (p_y_trt_leq_average[j] - p_y_ctl_leq_average[j]))
    # denominator <- sum(p_y[j]) + sum(p_y[k])
    # wRD_overall[s,1] <- numerator/denominator
    
    # fully bayesian approach
    numerator <- sum((p_y_ctl_est_average[j] + p_y_ctl_est_average[k]) * (p_y_trt_leq_average[j] - p_y_ctl_leq_average[j]))
    denominator <- sum(p_y_ctl_est_average[j]) + sum(p_y_ctl_est_average[k])
    wRD_bcontrol[s,1] <- numerator/denominator
    
    # cumulative weights
    numerator <- sum((p_y_ctl_leq_average[j] * (1 - p_y_ctl_leq_average[j])) * (p_y_trt_leq_average[j] - p_y_ctl_leq_average[j]))
    denominator <- sum((p_y_ctl_leq_average[j] * (1 - p_y_ctl_leq_average[j])))
    wRD_cum[s,1] <- numerator/denominator
    
    # both weights
    numerator <- sum((p_y_ctl_est_average[j] + p_y_ctl_est_average[k]) * (p_y_ctl_leq_average[j] * (1 - p_y_ctl_leq_average[j])) * (p_y_trt_leq_average[j] - p_y_ctl_leq_average[j]))
    denominator <- sum((p_y_ctl_est_average[j] + p_y_ctl_est_average[k]) * (p_y_ctl_leq_average[j] * (1 - p_y_ctl_leq_average[j])))
    wRD_both[s,1] <- numerator/denominator
    
    #calculated weighted RD with overall probability as weights
    # numerator <- sum((p_y_est_average[j] + p_y_est_average[k]) * (p_y_trt_leq_average[j] - p_y_ctl_leq_average[j]))
    # denominator <- sum(p_y_est_average[j]) + sum(p_y_est_average[k])
    # wRD_boverall[s,1] <- numerator/denominator
    
    #calculate log weighted OR with control probability as weights
    # numerator <- sum((p_y_ctl[j] + p_y_ctl[k]) * log(OR[j]))
    # denominator <- sum(p_y_ctl[j]) + sum(p_y_ctl[k])
    # log_wOR_control[s,1] <- numerator/denominator
    # wOR_control[s,1] <- exp(numerator/denominator)
    
    #calculate log weighted OR with overall probability as weights
    # numerator <- sum((p_y[j] + p_y[k]) * log(OR[j]))
    # denominator <- sum(p_y[j]) + sum(p_y[k])
    # log_wOR_overall[s,1] <- numerator/denominator
    # wOR_overall[s,1] <- exp(numerator/denominator)
    
    #fully bayesian approach
    #calculate log weighted OR with control probability as weights
    numerator <- sum((p_y_ctl_est_average[j] + p_y_ctl_est_average[k]) * log(OR[j]))
    denominator <- sum(p_y_ctl_est_average[j]) + sum(p_y_ctl_est_average[k])
    log_wOR_bcontrol[s,1] <- numerator/denominator
    wOR_bcontrol[s,1] <- exp(numerator/denominator)
    
    #calculate log weighted OR with cumulative
    numerator <- sum((p_y_ctl_leq_average[j] * (1 - p_y_ctl_leq_average[j])) * log(OR[j]))
    denominator <- sum((p_y_ctl_leq_average[j] * (1 - p_y_ctl_leq_average[j])))
    log_wOR_cum[s,1] <- numerator/denominator
    wOR_cum[s,1] <- exp(numerator/denominator)
    
    #calculate log weighted OR with both weights
    numerator <- sum((p_y_ctl_est_average[j] + p_y_ctl_est_average[k]) * (p_y_ctl_leq_average[j] * (1 - p_y_ctl_leq_average[j])) * log(OR[j]))
    denominator <- sum((p_y_ctl_est_average[j] + p_y_ctl_est_average[k]) * (p_y_ctl_leq_average[j] * (1 - p_y_ctl_leq_average[j])))
    log_wOR_both[s,1] <- numerator/denominator
    wOR_both[s,1] <- exp(numerator/denominator)
    
    #calculate log weighted OR with overall probability as weights
    # numerator <- sum((p_y_est_average[j] + p_y_est_average[k]) * log(OR[j]))
    # denominator <- sum(p_y_est_average[j]) + sum(p_y_est_average[k])
    # log_wOR_boverall[s,1] <- numerator/denominator
    # wOR_boverall[s,1] <- exp(numerator/denominator)
    
    #calculate log weighted RR with control probability as weights
    # numerator <- sum((p_y_ctl[j] + p_y_ctl[k]) * log(RR_leq[j]))
    # denominator <- sum(p_y_ctl[j]) + sum(p_y_ctl[k])
    # log_wRR_control[s,1] <- numerator/denominator
    # wRR_control[s,1] <- exp(numerator/denominator)
    
    #calculate log weighted RR with overall probability as weights
    # numerator <- sum((p_y[j] + p_y[k]) * log(RR_leq[j]))
    # denominator <- sum(p_y[j]) + sum(p_y[k])
    # log_wRR_overall[s,1] <- numerator/denominator
    # wRR_overall[s,1] <- exp(numerator/denominator)
    
    #bayesian
    #calculate log weighted RR with control probability as weights
    numerator <- sum((p_y_ctl_est_average[j] + p_y_ctl_est_average[k]) * log(RR_leq[j]))
    denominator <- sum(p_y_ctl_est_average[j]) + sum(p_y_ctl_est_average[k])
    log_wRR_bcontrol[s,1] <- numerator/denominator
    wRR_bcontrol[s,1] <- exp(numerator/denominator)
    
    #cumulative weights 
    numerator <- sum((p_y_ctl_leq_average[j] * (1 - p_y_ctl_leq_average[j])) * log(RR_leq[j]))
    denominator <- sum((p_y_ctl_leq_average[j] * (1 - p_y_ctl_leq_average[j])))
    log_wRR_cum[s,1] <- numerator/denominator
    wRR_cum[s,1] <- exp(numerator/denominator)
    
    #both weights 
    numerator <- sum((p_y_ctl_est_average[j] + p_y_ctl_est_average[k]) * (p_y_ctl_leq_average[j] * (1 - p_y_ctl_leq_average[j])) * log(RR_leq[j]))
    denominator <- sum((p_y_ctl_est_average[j] + p_y_ctl_est_average[k]) * (p_y_ctl_leq_average[j] * (1 - p_y_ctl_leq_average[j])))
    log_wRR_both[s,1] <- numerator/denominator
    wRR_both[s,1] <- exp(numerator/denominator)
    
    #both weights 
    numerator <- sum((p_y_ctl_est_average[j] + p_y_ctl_est_average[k]) * (p_y_ctl_leq_average[j] * (1 - p_y_ctl_leq_average[j])) * log(RR_geq[j]))
    denominator <- sum((p_y_ctl_est_average[j] + p_y_ctl_est_average[k]) * (p_y_ctl_leq_average[j] * (1 - p_y_ctl_leq_average[j])))
    log_wRR_ge_both[s,1] <- numerator/denominator
    wRR_ge_both[s,1] <- exp(numerator/denominator)
    
    #calculate log weighted RR with overall probability as weights
    # numerator <- sum((p_y_est_average[j] + p_y_est_average[k]) * log(RR_leq[j]))
    # denominator <- sum(p_y_est_average[j]) + sum(p_y_est_average[k])
    # log_wRR_boverall[s,1] <- numerator/denominator
    # wRR_boverall[s,1] <- exp(numerator/denominator)
    
    
    #calculate WR    
    l <- 1:c
    NB[s,1] <- sum(p_y_trt[j] * p_y_ctl_geq_average[j]) - (sum(p_y_ctl[j] * p_y_trt_geq_average[j]))
    
    weights_split[s,] <- (p_y_ctl_est_average[j] + p_y_ctl_est_average[k]) / sum(p_y_ctl_est_average[j] + p_y_ctl_est_average[k])
    weights_cum[s,] <- (p_y_ctl_leq_average[j] * (1 - p_y_ctl_leq_average[j])) / sum((p_y_ctl_leq_average[j] * (1 - p_y_ctl_leq_average[j])))
    weights_both[s,] <- (p_y_ctl_est_average[j] + p_y_ctl_est_average[k]) * (p_y_ctl_leq_average[j] * (1 - p_y_ctl_leq_average[j])) / 
      sum((p_y_ctl_est_average[j] + p_y_ctl_est_average[k]) * (p_y_ctl_leq_average[j] * (1 - p_y_ctl_leq_average[j])))
    
    
  }
  
  RD_leq_95CI <- apply(RD_leq_save, 2, function (x) quantile(x, probs = c(0.025, 0.975)))
  RD_leq_postprob <- colSums(RD_leq_save > 0)/n_post_draws
  
  unweighted_RD <- rowMeans(RD_leq_save)
  unweighted_RD_95CI <- quantile(unweighted_RD, probs = c(0.025, 0.975))
  
  OR_95CI <- apply(OR_save, 2, function (x) quantile(x, probs = c(0.025, 0.975)))
  log_OR_95CI <- apply(log_OR_save, 2, function (x) quantile(x, probs = c(0.025, 0.975)))
  OR_postprob <- colSums(OR_save > 1)/n_post_draws
  
  unweighted_log_OR <- rowMeans(log_OR_save)
  unweighted_log_OR_95CI <- quantile(unweighted_log_OR, probs = c(0.025, 0.975))
  unweighted_OR_95CI <- quantile(exp(unweighted_log_OR), probs = c(0.025, 0.975))
  
  RR_leq_95CI <- apply(RR_leq_save, 2, function (x) quantile(x, probs = c(0.025, 0.975)))
  log_RR_leq_95CI <- apply(log_RR_leq_save, 2, function (x) quantile(x, probs = c(0.025, 0.975)))
  RR_leq_postprob <- colSums(RR_leq_save > 0)/n_post_draws
  
  unweighted_log_RR <- rowMeans(log_RR_leq_save)
  unweighted_log_RR_95CI <- quantile(unweighted_log_RR, probs = c(0.025, 0.975))
  unweighted_RR_95CI <- quantile(exp(unweighted_log_RR), probs = c(0.025, 0.975))
  
  RR_geq_95CI <- apply(RR_geq_save, 2, function (x) quantile(x, probs = c(0.025, 0.975)))
  RR_geq_postprob <- colSums(RR_geq_save > 0)/n_post_draws
  
  unweighted_log_RR_ge <- rowMeans(log_RR_geq_save)
  unweighted_log_RR_ge_95CI <- quantile(unweighted_log_RR_ge, probs = c(0.025, 0.975))
  unweighted_RR_ge_95CI <- quantile(exp(unweighted_log_RR_ge), probs = c(0.025, 0.975))
  
  
  # wRD_control_95CI <- quantile(wRD_control, probs = c(0.025, 0.975))
  # wRD_control_postprob <- sum(wRD_control > 0)/n_post_draws
  # 
  # wRD_overall_95CI <- quantile(wRD_overall, probs = c(0.025, 0.975))
  # wRD_overall_postprob <- sum(wRD_overall > 0)/n_post_draws
  # 
  # wOR_control_95CI <- quantile(wOR_control, probs = c(0.025, 0.975))
  # wOR_control_postprob <- sum(wOR_control > 1)/n_post_draws
  # 
  # wOR_overall_95CI <- quantile(wOR_overall, probs = c(0.025, 0.975))
  # wOR_overall_postprob <- sum(wOR_overall > 1)/n_post_draws
  # 
  # log_wOR_control_95CI <- quantile(log_wOR_control, probs = c(0.025, 0.975))
  # log_wOR_control_postprob <- sum(log_wOR_control > 0)/n_post_draws
  # 
  # log_wOR_overall_95CI <- quantile(log_wOR_overall, probs = c(0.025, 0.975))
  # log_wOR_overall_postprob <- sum(log_wOR_overall > 0)/n_post_draws
  # 
  # wRR_control_95CI <- quantile(wRR_control, probs = c(0.025, 0.975))
  # wRR_control_postprob <- sum(wRR_control > 1)/n_post_draws
  # 
  # wRR_overall_95CI <- quantile(wRR_overall, probs = c(0.025, 0.975))
  # wRR_overall_postprob <- sum(wRR_overall > 1)/n_post_draws
  # 
  # log_wRR_control_95CI <- quantile(log_wRR_control, probs = c(0.025, 0.975))
  # log_wRR_control_postprob <- sum(log_wRR_control > 0)/n_post_draws
  # 
  # log_wRR_overall_95CI <- quantile(log_wRR_overall, probs = c(0.025, 0.975))
  # log_wRR_overall_postprob <- sum(log_wRR_overall > 0)/n_post_draws
  
  # fully bayesian approach 
  wRD_bcontrol_95CI <- quantile(wRD_bcontrol, probs = c(0.025, 0.975))
  wRD_bcontrol_postprob <- sum(wRD_bcontrol > 0)/n_post_draws
  
  wRD_cum_95CI <- quantile(wRD_cum, probs = c(0.025, 0.975))
  wRD_cum_postprob <- sum(wRD_cum > 0)/n_post_draws
  
  wRD_both_95CI <- quantile(wRD_both, probs = c(0.025, 0.975))
  wRD_both_postprob <- sum(wRD_both > 0)/n_post_draws
  # wRD_boverall_95CI <- quantile(wRD_boverall, probs = c(0.025, 0.975))
  # wRD_boverall_postprob <- sum(wRD_boverall > 0)/n_post_draws
  
  wOR_bcontrol_95CI <- quantile(wOR_bcontrol, probs = c(0.025, 0.975))
  wOR_bcontrol_postprob <- sum(wOR_bcontrol > 1)/n_post_draws
  
  wOR_cum_95CI <- quantile(wOR_cum, probs = c(0.025, 0.975))
  wOR_cum_postprob <- sum(wOR_cum > 1)/n_post_draws
  
  wOR_both_95CI <- quantile(wOR_both, probs = c(0.025, 0.975))
  wOR_both_postprob <- sum(wOR_both > 1)/n_post_draws
  
  # wOR_boverall_95CI <- quantile(wOR_boverall, probs = c(0.025, 0.975))
  # wOR_boverall_postprob <- sum(wOR_boverall > 1)/n_post_draws
  # 
  log_wOR_bcontrol_95CI <- quantile(log_wOR_bcontrol, probs = c(0.025, 0.975))
  log_wOR_bcontrol_postprob <- sum(log_wOR_bcontrol > 0)/n_post_draws
  
  log_wOR_cum_95CI <- quantile(log_wOR_cum, probs = c(0.025, 0.975))
  log_wOR_cum_postprob <- sum(log_wOR_cum > 0)/n_post_draws
  
  log_wOR_both_95CI <- quantile(log_wOR_both, probs = c(0.025, 0.975))
  log_wOR_both_postprob <- sum(log_wOR_both > 0)/n_post_draws
  
  # log_wOR_boverall_95CI <- quantile(log_wOR_boverall, probs = c(0.025, 0.975))
  # log_wOR_boverall_postprob <- sum(log_wOR_boverall > 0)/n_post_draws
  
  wRR_bcontrol_95CI <- quantile(wRR_bcontrol, probs = c(0.025, 0.975))
  wRR_bcontrol_postprob <- sum(wRR_bcontrol > 1)/n_post_draws
  
  wRR_cum_95CI <- quantile(wRR_cum, probs = c(0.025, 0.975))
  wRR_cum_postprob <- sum(wRR_cum > 1)/n_post_draws
  
  wRR_both_95CI <- quantile(wRR_both, probs = c(0.025, 0.975))
  wRR_both_postprob <- sum(wRR_both > 1)/n_post_draws
  
  wRR_ge_both_95CI <- quantile(wRR_ge_both, probs = c(0.025, 0.975))
  wRR_ge_both_postprob <- sum(wRR_ge_both > 1)/n_post_draws
  
  # wRR_boverall_95CI <- quantile(wRR_boverall, probs = c(0.025, 0.975))
  # wRR_boverall_postprob <- sum(wRR_boverall > 1)/n_post_draws
  
  log_wRR_bcontrol_95CI <- quantile(log_wRR_bcontrol, probs = c(0.025, 0.975))
  log_wRR_bcontrol_postprob <- sum(log_wRR_bcontrol > 0)/n_post_draws
  
  log_wRR_cum_95CI <- quantile(log_wRR_cum, probs = c(0.025, 0.975))
  log_wRR_cum_postprob <- sum(log_wRR_cum > 0)/n_post_draws
  
  log_wRR_both_95CI <- quantile(log_wRR_both, probs = c(0.025, 0.975))
  log_wRR_both_postprob <- sum(log_wRR_both > 0)/n_post_draws
  
  log_wRR_ge_both_95CI <- quantile(log_wRR_ge_both, probs = c(0.025, 0.975))
  log_wRR_ge_both_postprob <- sum(log_wRR_ge_both > 0)/n_post_draws
  
  # log_wRR_boverall_95CI <- quantile(log_wRR_boverall, probs = c(0.025, 0.975))
  # log_wRR_boverall_postprob <- sum(log_wRR_boverall > 0)/n_post_draws
  
  NB_95CI <- quantile(NB, probs = c(0.025, 0.975))
  NB_postprob <- sum(NB > 0)/n_post_draws
  
  weights_split_95CI <- apply(weights_split, 2, quantile, probs = c(0.025, 0.975))
  weights_cum_95CI <- apply(weights_cum, 2, quantile, probs = c(0.025, 0.975))
  weights_both_95CI <- apply(weights_both, 2, quantile, probs = c(0.025, 0.975))
  
  
  
  return(list(RD_leq_95CI = RD_leq_95CI,
              RD_leq_postprob = RD_leq_postprob,
              OR_95CI = OR_95CI,
              log_OR_95CI = log_OR_95CI,
              OR_postprob = OR_postprob,
              RR_leq_95CI = RR_leq_95CI,
              log_RR_leq_95CI = log_RR_leq_95CI,
              RR_leq_postprob = RR_leq_postprob,
              RR_geq_95CI = RR_geq_95CI,
              # wRD_control_95CI = wRD_control_95CI, 
              # wRD_control_postprob = wRD_control_postprob, 
              # wRD_overall_95CI = wRD_overall_95CI, 
              # wRD_overall_postprob = wRD_overall_postprob, 
              # wOR_control_95CI = wOR_control_95CI, 
              # wOR_control_postprob = wOR_control_postprob, 
              # wOR_overall_95CI = wOR_overall_95CI, 
              # wOR_overall_postprob = wOR_overall_postprob, 
              # log_wOR_control_95CI = log_wOR_control_95CI, 
              # log_wOR_control_postprob = log_wOR_control_postprob, 
              # log_wOR_overall_95CI = log_wOR_overall_95CI, 
              # log_wOR_overall_postprob = log_wOR_overall_postprob, 
              # wRR_control_95CI = wRR_control_95CI, 
              # wRR_control_postprob = wRR_control_postprob, 
              # wRR_overall_95CI = wRR_overall_95CI, 
              # wRR_overall_postprob = wRR_overall_postprob, 
              # log_wRR_control_95CI = log_wRR_control_95CI, 
              # log_wRR_control_postprob = log_wRR_control_postprob, 
              # log_wRR_overall_95CI = log_wRR_overall_95CI, 
              # log_wRR_overall_postprob = log_wRR_overall_postprob, 
              
              wRD_bcontrol_95CI = wRD_bcontrol_95CI, 
              wRD_bcontrol_postprob = wRD_bcontrol_postprob, 
              wRD_cum_95CI = wRD_cum_95CI, 
              wRD_cum_postprob = wRD_cum_postprob, 
              wRD_both_95CI = wRD_both_95CI, 
              wRD_both_postprob = wRD_both_postprob, 
              # wRD_boverall_95CI = wRD_boverall_95CI, 
              # wRD_boverall_postprob = wRD_boverall_postprob, 
              wOR_bcontrol_95CI = wOR_bcontrol_95CI, 
              wOR_bcontrol_postprob = wOR_bcontrol_postprob, 
              wOR_cum_95CI = wOR_cum_95CI, 
              wOR_cum_postprob = wOR_cum_postprob, 
              wOR_both_95CI = wOR_both_95CI, 
              wOR_both_postprob = wOR_both_postprob, 
              # wOR_boverall_95CI = wOR_boverall_95CI, 
              # wOR_boverall_postprob = wOR_boverall_postprob, 
              log_wOR_bcontrol_95CI = log_wOR_bcontrol_95CI, 
              log_wOR_bcontrol_postprob = log_wOR_bcontrol_postprob, 
              log_wOR_cum_95CI = log_wOR_cum_95CI, 
              log_wOR_cum_postprob = log_wOR_cum_postprob, 
              log_wOR_both_95CI = log_wOR_both_95CI, 
              log_wOR_both_postprob = log_wOR_both_postprob, 
              # log_wOR_boverall_95CI = log_wOR_boverall_95CI, 
              # log_wOR_boverall_postprob = log_wOR_boverall_postprob, 
              wRR_bcontrol_95CI = wRR_bcontrol_95CI, 
              wRR_bcontrol_postprob = wRR_bcontrol_postprob, 
              wRR_cum_95CI = wRR_cum_95CI, 
              wRR_cum_postprob = wRR_cum_postprob, 
              wRR_both_95CI = wRR_both_95CI, 
              wRR_both_postprob = wRR_both_postprob, 
              wRR_ge_both_95CI = wRR_ge_both_95CI, 
              wRR_ge_both_postprob = wRR_ge_both_postprob, 
              # wRR_boverall_95CI = wRR_boverall_95CI, 
              # wRR_boverall_postprob = wRR_boverall_postprob, 
              log_wRR_bcontrol_95CI = log_wRR_bcontrol_95CI, 
              log_wRR_bcontrol_postprob = log_wRR_bcontrol_postprob, 
              log_wRR_cum_95CI = log_wRR_cum_95CI, 
              log_wRR_cum_postprob = log_wRR_cum_postprob, 
              log_wRR_both_95CI = log_wRR_both_95CI, 
              log_wRR_both_postprob = log_wRR_both_postprob, 
              log_wRR_ge_both_95CI = log_wRR_ge_both_95CI, 
              log_wRR_ge_both_postprob = log_wRR_ge_both_postprob, 
              # log_wRR_boverall_95CI = log_wRR_boverall_95CI, 
              # log_wRR_boverall_postprob = log_wRR_boverall_postprob, 
              
              unweighted_RD_95CI = unweighted_RD_95CI,
              unweighted_log_OR_95CI = unweighted_log_OR_95CI,
              unweighted_OR_95CI = unweighted_OR_95CI,
              unweighted_log_RR_95CI = unweighted_log_RR_95CI,
              unweighted_RR_95CI = unweighted_RR_95CI,
              unweighted_log_RR_ge_95CI = unweighted_log_RR_ge_95CI,
              unweighted_RR_ge_95CI = unweighted_RR_ge_95CI,
              
              NB_95CI = NB_95CI, 
              NB_postprob = NB_postprob,
              
              weights_split_95CI = weights_split_95CI, 
              weights_cum_95CI = weights_cum_95CI,
              weights_both_95CI = weights_both_95CI))
  
}

# write function to calculate order invariant weighted risk difference using polr 
polr_order_weighted_RD <- function(b = b, c) { 
  p_ctl_more <- NA
  p_trt_more <- NA
  p_ctl_less <- NA
  p_trt_less <- NA
  p_y_ctl <- NA
  
  for (i in 2:c) {
    # calculate p(y >= i | control) for i = 2, ... c
    p_ctl_more[i] <- 1 - plogis(b$zeta[i-1] + 0.5 * b$coefficients[1])
    # calculate p(y >= i | treatment) for i = 2, ... c
    # get all coefficients that have y>=i
    p_trt_more[i] <- 1 - plogis(b$zeta[i-1] - 0.5 * b$coefficients[1])
    
    # calculate p(y <= i | control) for i = 1, ..., c-1
    p_ctl_less[i-1] <- plogis(b$zeta[i-1] + 0.5 * b$coefficients[1]) 
    # calculate p(y <= i | treatment) for i = 2, ... c
    # get all coefficients that have y>=i
    p_trt_less[i-1] <- plogis(b$zeta[i-1] - 0.5 * b$coefficients[1])
  }
  
  # calculate p(Y = i) for i = 1, ... , c
  for (i in 1:c) {
    if (i == 1) {
      p_y_ctl[1] <-  p_ctl_less[1]
    } else if (i == c) {
      p_y_ctl[c] <- 1 - sum(p_y_ctl[1:c-1])
    } else {
      p_y_ctl[i] <-  p_ctl_less[i] - sum(p_y_ctl[1:i-1])
    }
  }
  
  RD_leq <- p_trt_less - p_ctl_less
  RD_geq <- p_ctl_more - p_trt_more
  
  OR <- (p_trt_less / (1- p_trt_less)) / (p_ctl_less / (1- p_ctl_less))
  
  j <- 1:c-1
  k <- 2:c
  numerator <- sum(p_y_ctl[j] * (p_trt_less[j] - p_ctl_less[j])) + sum(p_y_ctl[k] * (p_ctl_more[k] - p_trt_more[k]))
  denominator <- sum(p_y_ctl[j]) + sum(p_y_ctl[k])
  weighted_RD <- numerator/denominator
  return(list(RD_leq = RD_leq, RD_geq = RD_geq, OR = OR, weighted_RD = weighted_RD))
  
}
