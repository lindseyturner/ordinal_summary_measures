# simulation setup 
# data generation scenarios 
# two different data generations each with 5 categories
# the first scenario will mimic the covid out data
# 70%, 18%, 9%, 2%, 1%
# second scenario will have equal distribution to outcomes 
# 20%, 20%, 20%, 20%, 20%

# different treatment effects 
# null scenario - treatment and control are the same 
# proportional odds scenario - 2 different settings - one small effect and one large effect
# non proportional odds scenario - relatively constant risk difference

# what to output: 
# point estimates and credible intervals for each model (PPO, PO, polr, mlr) and 
# summary measure (RD, weighted RD, OR, weighted OR, WR)

# what to measure: 
# power, average z score (standardized effect size), bias, rMSE, coverage rates

# to do: 
# calculate the small and large effects for the proportional odds scenario 
# determine the non proportional odds scenario
# generate data 
# calculate point estimates and credible intervals 
# calculate measures



# generate data 
# function to generate data 
# start with 50% to treatment and 50% to control 
#

determine_start_i <- function(val) {
  if (val%%5 == 1){
    start <- 1
  } else if (val%%5 == 2) {
    start <- 401
  } else if (val%%5 == 3) {
    start <- 801
  } else if (val%%5 == 4) {
    start <- 1201
  } else if (val%%5 == 0) {
    start <- 1601
  }
  
  return(start)
}

determine_start_i2 <- function(val) {
  if (val%%2 == 1){
    start <- 1
  } else if (val%%2 == 0) {
    start <- 1001
  }
  
  return(start)
}

DataGen <- function(n) {
  pc1 <- c(0.2, 0.2, 0.2, 0.2, 0.2)
  pt_1_RD_constant <- c(0.25, 0.2, 0.2, 0.2, 0.15)
  beta1 <- 0.1
  pt_1_PO_low <- c(invlogit(logit(pc1[1]) + beta1),
                   invlogit(logit(pc1[1] + pc1[2]) + beta1) - invlogit(logit(pc1[1]) + beta1),
                   invlogit(logit(pc1[1] + pc1[2] + pc1[3]) + beta1) - invlogit(logit(pc1[1] + pc1[2]) + beta1),
                   invlogit(logit(pc1[1] + pc1[2] + pc1[3] + pc1[4]) + beta1) - invlogit(logit(pc1[1] + pc1[2] + pc1[3]) + beta1),
                   1 - invlogit(logit(pc1[1] + pc1[2] + pc1[3] + pc1[4]) + beta1))
  beta2 <- 0.5
  pt_1_PO_high <- c(invlogit(logit(pc1[1]) + beta2),
                    invlogit(logit(pc1[1] + pc1[2]) + beta2) - invlogit(logit(pc1[1]) + beta2),
                    invlogit(logit(pc1[1] + pc1[2] + pc1[3]) + beta2) - invlogit(logit(pc1[1] + pc1[2]) + beta2),
                    invlogit(logit(pc1[1] + pc1[2] + pc1[3] + pc1[4]) + beta2) - invlogit(logit(pc1[1] + pc1[2] + pc1[3]) + beta2),
                    1 - invlogit(logit(pc1[1] + pc1[2] + pc1[3] + pc1[4]) + beta2))
  
  # scenario 1, null treatment effect
  dat_1_null <- data.frame(matrix(NA, nrow = n, ncol = 3))
  colnames(dat_1_null) <- c("y", "trt", "trt_05")
  dat_1_null$trt <- rbinom(n, 1, 0.5)
  dat_1_null$trt_05 <- case_when(dat_1_null$trt == 0 ~ -0.5, 
                                 dat_1_null$trt == 1 ~ 0.5)
  dat_1_null$y <- sample(1:5, n, replace = TRUE, prob = pc1)
  
  #scenario 1, constant RD 
  dat_1_RD_constant <- data.frame(matrix(NA, nrow = n, ncol = 3))
  colnames(dat_1_RD_constant) <- c("y", "trt", "trt_05")
  dat_1_RD_constant$trt <- rbinom(n, 1, 0.5)
  dat_1_RD_constant$trt_05 <- case_when(dat_1_RD_constant$trt == 0 ~ -0.5, 
                                        dat_1_RD_constant$trt == 1 ~ 0.5)
  dat_1_RD_constant$y <- case_when(dat_1_RD_constant$trt == 0 ~ sample(1:5,n,replace=TRUE,prob=pc1), 
                                   dat_1_RD_constant$trt == 1 ~ sample(1:5,n,replace=TRUE,prob=pt_1_RD_constant))
  
  # scenario 1, small PO treatment effect 
  dat_1_RD_PO_low <- data.frame(matrix(NA, nrow = n, ncol = 3))
  colnames(dat_1_RD_PO_low) <- c("y", "trt", "trt_05")
  dat_1_RD_PO_low$trt <- rbinom(n, 1, 0.5)
  dat_1_RD_PO_low$trt_05 <- case_when(dat_1_RD_PO_low$trt == 0 ~ -0.5, 
                                      dat_1_RD_PO_low$trt == 1 ~ 0.5)
  dat_1_RD_PO_low$y <- case_when(dat_1_RD_PO_low$trt == 0 ~ sample(1:5,n,replace=TRUE,prob=pc1), 
                                 dat_1_RD_PO_low$trt == 1 ~ sample(1:5,n,replace=TRUE,prob=pt_1_PO_low))
  
  
  # scenario 1, large PO treatment effect 
  dat_1_RD_PO_high <- data.frame(matrix(NA, nrow = n, ncol = 3))
  colnames(dat_1_RD_PO_high) <- c("y", "trt", "trt_05")
  dat_1_RD_PO_high$trt <- rbinom(n, 1, 0.5)
  dat_1_RD_PO_high$trt_05 <- case_when(dat_1_RD_PO_high$trt == 0 ~ -0.5, 
                                       dat_1_RD_PO_high$trt == 1 ~ 0.5)
  dat_1_RD_PO_high$y <- case_when(dat_1_RD_PO_high$trt == 0 ~ sample(1:5,n,replace=TRUE,prob=pc1), 
                                  dat_1_RD_PO_high$trt == 1 ~ sample(1:5,n,replace=TRUE,prob=pt_1_PO_high))
  
  
  # scenario 2 probabilities
  pc2 <- c(0.70, 0.18, 0.09, 0.02, 0.01)
  pt_2_RD_constant <- c(0.71, 0.18, 0.09, 0.015, 0.005)
  pt_2_PO_low <- c(invlogit(logit(pc2[1]) + beta1),
                   invlogit(logit(pc2[1] + pc2[2]) + beta1) - invlogit(logit(pc2[1]) + beta1),
                   invlogit(logit(pc2[1] + pc2[2] + pc2[3]) + beta1) - invlogit(logit(pc2[1] + pc2[2]) + beta1),
                   invlogit(logit(pc2[1] + pc2[2] + pc2[3] + pc2[4]) + beta1) - invlogit(logit(pc2[1] + pc2[2] + pc2[3]) + beta1),
                   1 - invlogit(logit(pc2[1] + pc2[2] + pc2[3] + pc2[4]) + beta1))
  pt_2_PO_high <- c(invlogit(logit(pc2[1]) + beta2),
                    invlogit(logit(pc2[1] + pc2[2]) + beta2) - invlogit(logit(pc2[1]) + beta2),
                    invlogit(logit(pc2[1] + pc2[2] + pc2[3]) + beta2) - invlogit(logit(pc2[1] + pc2[2]) + beta2),
                    invlogit(logit(pc2[1] + pc2[2] + pc2[3] + pc2[4]) + beta2) - invlogit(logit(pc2[1] + pc2[2] + pc2[3]) + beta2),
                    1 - invlogit(logit(pc2[1] + pc2[2] + pc2[3] + pc2[4]) + beta2))
  
  # scenario 1, null treatment effect
  dat_2_null <- data.frame(matrix(NA, nrow = n, ncol = 3))
  colnames(dat_2_null) <- c("y", "trt", "trt_05")
  dat_2_null$trt <- rbinom(n, 1, 0.5)
  dat_2_null$trt_05 <- case_when(dat_2_null$trt == 0 ~ -0.5, 
                                 dat_2_null$trt == 1 ~ 0.5)
  dat_2_null$y <- sample(1:5, n, replace = TRUE, prob = pc2)
  
  #scenario 1, constant RD 
  dat_2_RD_constant <- data.frame(matrix(NA, nrow = n, ncol = 3))
  colnames(dat_2_RD_constant) <- c("y", "trt", "trt_05")
  dat_2_RD_constant$trt <- rbinom(n, 1, 0.5)
  dat_2_RD_constant$trt_05 <- case_when(dat_2_RD_constant$trt == 0 ~ -0.5, 
                                        dat_2_RD_constant$trt == 1 ~ 0.5)
  dat_2_RD_constant$y <- case_when(dat_2_RD_constant$trt == 0 ~ sample(1:5,n,replace=TRUE,prob=pc2), 
                                   dat_2_RD_constant$trt == 1 ~ sample(1:5,n,replace=TRUE,prob=pt_2_RD_constant))
  
  # scenario 2, small PO treatment effect 
  dat_2_RD_PO_low <- data.frame(matrix(NA, nrow = n, ncol = 3))
  colnames(dat_2_RD_PO_low) <- c("y", "trt", "trt_05")
  dat_2_RD_PO_low$trt <- rbinom(n, 1, 0.5)
  dat_2_RD_PO_low$trt_05 <- case_when(dat_2_RD_PO_low$trt == 0 ~ -0.5, 
                                      dat_2_RD_PO_low$trt == 1 ~ 0.5)
  dat_2_RD_PO_low$y <- case_when(dat_2_RD_PO_low$trt == 0 ~ sample(1:5,n,replace=TRUE,prob=pc2), 
                                 dat_2_RD_PO_low$trt == 1 ~ sample(1:5,n,replace=TRUE,prob=pt_2_PO_low))
  
  
  # scenario 2, large PO treatment effect 
  dat_2_RD_PO_high <- data.frame(matrix(NA, nrow = n, ncol = 3))
  colnames(dat_2_RD_PO_high) <- c("y", "trt", "trt_05")
  dat_2_RD_PO_high$trt <- rbinom(n, 1, 0.5)
  dat_2_RD_PO_high$trt_05 <- case_when(dat_2_RD_PO_high$trt == 0 ~ -0.5, 
                                       dat_2_RD_PO_high$trt == 1 ~ 0.5)
  dat_2_RD_PO_high$y <- case_when(dat_2_RD_PO_high$trt == 0 ~ sample(1:5,n,replace=TRUE,prob=pc2), 
                                  dat_2_RD_PO_high$trt == 1 ~ sample(1:5,n,replace=TRUE,prob=pt_2_PO_high))
  
  
  DFs <- list(
    dat_1_null = dat_1_null,
    dat_1_RD_constant = dat_1_RD_constant,
    dat_1_RD_PO_low = dat_1_RD_PO_low,
    dat_1_RD_PO_high = dat_1_RD_PO_high,
    dat_2_null = dat_2_null,
    dat_2_RD_constant = dat_2_RD_constant,
    dat_2_RD_PO_low = dat_2_RD_PO_low,
    dat_2_RD_PO_high = dat_2_RD_PO_high
  )
  
  return(DFs)
  
}


# data generation by array number 
DataGen_array <- function(n, array) {
  pc1 <- c(0.2, 0.2, 0.2, 0.2, 0.2)
  pc2 <- c(0.70, 0.18, 0.09, 0.02, 0.01)
  
  if (array %in% 1:2) {
    pc <- pc1
    pt <- pc1
  } else if (array %in% 3:4) {
    # RD low  is 0.06
    pc <- pc1
    pt <- c(0.26, 0.2, 0.2, 0.2, 0.14)
  } else if (array %in% 5:6) {
    # RD high  is 0.08
    pc <- pc1
    # was 0.28 and 0.12
    pt <- c(0.27, 0.2, 0.2, 0.2, 0.13)
  } else if (array %in% 7:8) {
    #RR low 
    pc <- pc1
    pt <- c(0.225, 0.225, 0.225, 0.225, 0.10)
  } else if (array %in% 9:10) {
    #RR high is 1.2
    pc <- pc1
    # old: c(0.34, 0.165, 0.165, 0.165, 0.165)
    pt <- c(0.23, 0.23, 0.23, 0.23, 0.08)
  } else if (array %in% 11:12) {
    #PO low is OR = 1.37
    pc <- pc1
    beta1 <- log(1.37)
    pt <- c(invlogit(logit(pc1[1]) + beta1),
                     invlogit(logit(pc1[1] + pc1[2]) + beta1) - invlogit(logit(pc1[1]) + beta1),
                     invlogit(logit(pc1[1] + pc1[2] + pc1[3]) + beta1) - invlogit(logit(pc1[1] + pc1[2]) + beta1),
                     invlogit(logit(pc1[1] + pc1[2] + pc1[3] + pc1[4]) + beta1) - invlogit(logit(pc1[1] + pc1[2] + pc1[3]) + beta1),
                     1 - invlogit(logit(pc1[1] + pc1[2] + pc1[3] + pc1[4]) + beta1))
  } else if (array %in% 13:14) {
    #PO high is OR = 1.47
    pc <- pc1
    beta2 <- log(1.47)
    pt <- c(invlogit(logit(pc1[1]) + beta2),
                      invlogit(logit(pc1[1] + pc1[2]) + beta2) - invlogit(logit(pc1[1]) + beta2),
                      invlogit(logit(pc1[1] + pc1[2] + pc1[3]) + beta2) - invlogit(logit(pc1[1] + pc1[2]) + beta2),
                      invlogit(logit(pc1[1] + pc1[2] + pc1[3] + pc1[4]) + beta2) - invlogit(logit(pc1[1] + pc1[2] + pc1[3]) + beta2),
                      1 - invlogit(logit(pc1[1] + pc1[2] + pc1[3] + pc1[4]) + beta2))
  } else if (array %in% 15:16) {
    # 20% reduction in the last 3 levels and then rest is put in level 1
    pc <- pc1
    # was c(0.32, 0.2, 0.16, 0.16, 0.16)
    pt <- c(0.29, 0.2, 0.17, 0.17, 0.17)
  } else if (array %in% 17:18) {
    # 25% reduction in the last 3 levels and then rest is put in level 2
    pc <- pc1
    pt <- c(0.2, 0.35, 0.15, 0.15, 0.15)
  } else if (array %in% 19:20) {
    pc <- pc2
    pt <- pc2
  } else if (array %in% 21:22) {
    #RD low is 0.03ish
    pc <- pc2
    pt <- c(0.76, 0.18, 0.05, 0.005, 0.005)
  } else if (array %in% 23:24) {
    #RD high is 0.05ish
    pc <- pc2
    pt <- c(0.78, 0.17, 0.04, 0.005, 0.005)
  } else if (array %in% 25:26) {
    #RR low is 0.9
    pc <- pc2
    #pt <- c(0.73, 0.162, 0.081, 0.018, 0.009)
    # 0.75
    pt <- c(0.78, 0.185, 0.025, 0.0075, 0.0025)
  } else if (array %in% 27:28) {
    #RR high is 0.8
    pc <- pc2
    #pt <- c(0.76, 0.144, 0.072, 0.016, 0.008)
    #try 0.6
    pt <- c(0.77, 0.18, 0.025, 0.0175, 0.0075)
  } else if (array %in% 29:30) {
    #OR low is 1.47
    pc <- pc2
    beta1_2 <- log(1.47)
    
    pt <- c(invlogit(logit(pc2[1]) + beta1_2),
                     invlogit(logit(pc2[1] + pc2[2]) + beta1_2) - invlogit(logit(pc2[1]) + beta1_2),
                     invlogit(logit(pc2[1] + pc2[2] + pc2[3]) + beta1_2) - invlogit(logit(pc2[1] + pc2[2]) + beta1_2),
                     invlogit(logit(pc2[1] + pc2[2] + pc2[3] + pc2[4]) + beta1_2) - invlogit(logit(pc2[1] + pc2[2] + pc2[3]) + beta1_2),
                     1 - invlogit(logit(pc2[1] + pc2[2] + pc2[3] + pc2[4]) + beta1_2))

  } else if (array %in% 31:32) {
    #OR high is 1.60
    pc <- pc2
    beta2_2 <- log(1.60)
    pt <- c(invlogit(logit(pc2[1]) + beta2_2),
                      invlogit(logit(pc2[1] + pc2[2]) + beta2_2) - invlogit(logit(pc2[1]) + beta2_2),
                      invlogit(logit(pc2[1] + pc2[2] + pc2[3]) + beta2_2) - invlogit(logit(pc2[1] + pc2[2]) + beta2_2),
                      invlogit(logit(pc2[1] + pc2[2] + pc2[3] + pc2[4]) + beta2_2) - invlogit(logit(pc2[1] + pc2[2] + pc2[3]) + beta2_2),
                      1 - invlogit(logit(pc2[1] + pc2[2] + pc2[3] + pc2[4]) + beta2_2))
  } else if (array %in% 33:34) {
    # 50% reduction in the last 3 levels and then rest is put in level 1
    pc <- pc2
    pt <- c(0.772, 0.18, 0.036, 0.008, 0.004)
  } else if (array %in% 35:36) {
    # 50% reduction in the last 3 levels and then rest is put in level 2
    pc <- pc2
    #pt <- c(0.70, 0.252, 0.036, 0.008, 0.004)
    #try 75% reduction 
    pt <- c(0.70, 0.27, 0.0225, 0.0050, 0.0025)
    
  }

  df <- data.frame(matrix(NA, nrow = n, ncol = 3))
  colnames(df) <- c("y", "trt", "trt_05")
  df$trt <- rbinom(n, 1, 0.5)
  df$trt_05 <- case_when(df$trt == 0 ~ -0.5, 
                         df$trt == 1 ~ 0.5)
  df$y <- case_when(df$trt == 0 ~ sample(1:5,n,replace=TRUE,prob=pc), 
                    df$trt == 1 ~ sample(1:5,n,replace=TRUE,prob=pt))
  
  return(df)
  
}

# functions for partial proportional odds model, proportional odds model, polr model
# functions for both point estimates and confidence intervals

# new function for weighted risk difference 
# assumes a partial proportional odds model 
# assumes the treatment is coded as -0.5 for control and 0.5 for treatment 
order_weighted_RD2 <- function(b = b, iprior = 1, dat) { 
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
  
  PPO_leq <- (PPO_ctl_leq + PPO_trt_leq)/2
  
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
  
  p_y_est <- (p_y_ctl_est + p_y_trt_est)/2
  
  # get weights from the data 
  p_y_ctl <- unname(table(dat$trt, dat$y)[1,])/500
  p_y_trt <- unname(table(dat$trt, dat$y)[2,])/500
  p_y <- unname(table(dat$y))/1000
  
  
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
  numerator <- sum((p_y_ctl[j] + p_y_ctl[k]) * (PPO_trt_leq[j] - PPO_ctl_leq[j]))
  denominator <- sum(p_y_ctl[j]) + sum(p_y_ctl[k])
  weighted_RD <- numerator/denominator
  
  
  numerator <- sum((p_y[j] + p_y[k]) * (PPO_trt_leq[j] - PPO_ctl_leq[j]))
  denominator <- sum(p_y[j]) + sum(p_y[k])
  weighted_RD_overall <- numerator/denominator
  
  # fully bayesian approach
  numerator <- sum((p_y_ctl_est[j] + p_y_ctl_est[k]) * (PPO_trt_leq[j] - PPO_ctl_leq[j]))
  denominator <- sum(p_y_ctl_est[j]) + sum(p_y_ctl_est[k])
  weighted_RD_b <- numerator/denominator
  
  
  numerator <- sum((p_y_est[j] + p_y_est[k]) * (PPO_trt_leq[j] - PPO_ctl_leq[j]))
  denominator <- sum(p_y_est[j]) + sum(p_y_est[k])
  weighted_RD_boverall <- numerator/denominator
  
  unweighted_RD <- mean(PPO_trt_leq - PPO_ctl_leq)
  
  # fully bayesian approach with cumulative weights 
  numerator <- sum((PPO_ctl_leq[j] * (1-PPO_ctl_leq[j])) * (PPO_trt_leq[j] - PPO_ctl_leq[j]))
  denominator <- sum((PPO_ctl_leq[j] * (1-PPO_ctl_leq[j])))
  weighted_RD_cum <- numerator/denominator
  
  # fully bayesian approach with both weights 
  numerator <- sum((p_y_ctl_est[j] + p_y_ctl_est[k]) * (PPO_ctl_leq[j] * (1-PPO_ctl_leq[j])) * (PPO_trt_leq[j] - PPO_ctl_leq[j]))
  denominator <- sum((p_y_ctl_est[j] + p_y_ctl_est[k]) * (PPO_ctl_leq[j] * (1-PPO_ctl_leq[j])))
  weighted_RD_both <- numerator/denominator
  
  # fully bayesian approach with both weights with overall prob
  numerator <- sum((p_y_est[j] + p_y_est[k]) * (PPO_leq[j] * (1-PPO_leq[j])) * (PPO_trt_leq[j] - PPO_ctl_leq[j]))
  denominator <- sum((p_y_est[j] + p_y_est[k]) * (PPO_leq[j] * (1-PPO_leq[j])))
  weighted_RD_both_ov <- numerator/denominator

  #calculate weighted OR 
  numerator <- sum((p_y_ctl[j] + p_y_ctl[k])  * log(OR_le[j]))
  denominator <- sum(p_y_ctl[j]) + sum(p_y_ctl[k])
  weighted_OR <- exp(numerator/denominator)
  weighted_log_OR <- numerator/denominator
  
  #calculate weighted OR using overall probability as weights
  numerator <- sum((p_y[j] + p_y[k]) * log(OR_le[j]))
  denominator <- sum(p_y[j]) + sum(p_y[k])
  weighted_OR_overall <- exp(numerator/denominator)
  weighted_log_OR_overall <- numerator/denominator
  
  # fully bayesian approach 
  #calculate weighted OR 
  numerator <- sum((p_y_ctl_est[j] + p_y_ctl_est[k])  * log(OR_le[j]))
  denominator <- sum(p_y_ctl_est[j]) + sum(p_y_ctl_est[k])
  weighted_OR_b <- exp(numerator/denominator)
  weighted_log_OR_b <- numerator/denominator
  
  #calculate weighted OR using overall probability as weights
  numerator <- sum((p_y_est[j] + p_y_est[k]) * log(OR_le[j]))
  denominator <- sum(p_y_est[j]) + sum(p_y_est[k])
  weighted_OR_boverall <- exp(numerator/denominator)
  weighted_log_OR_boverall <- numerator/denominator
  
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
  
  # fully bayesian with both cumulative and split weights
  numerator <- sum((p_y_est[j] + p_y_est[k]) * (PPO_leq[j] * (1-PPO_leq[j])) * log(OR_le[j]))
  denominator <- sum((p_y_est[j] + p_y_est[k]) * (PPO_leq[j] * (1-PPO_leq[j])))
  weighted_OR_both_ov <- exp(numerator/denominator)
  weighted_log_OR_both_ov <- numerator/denominator
  
  unweighted_log_OR <- mean(log(OR_le))
  unweighted_OR <- exp(mean(log(OR_le)))
  
  # calculate weighted RR with control probabilities 
  numerator <- sum((p_y_ctl[j] + p_y_ctl[k]) * log_RR_le[j])
  denominator <- sum(p_y_ctl[j]) + sum(p_y_ctl[k])
  weighted_log_RR_control <- numerator/denominator
  
  # calculate weighted RR with overall probabilities 
  numerator <- sum((p_y[j] + p_y[k]) * log_RR_le[j]) 
  denominator <- sum(p_y[j]) + sum(p_y[k])
  weighted_log_RR_overall <- numerator/denominator
  
  # fully bayesian approach 
  # calculate weighted RR with control probabilities 
  numerator <- sum((p_y_ctl_est[j] + p_y_ctl_est[k]) * log_RR_le[j])
  denominator <- sum(p_y_ctl_est[j]) + sum(p_y_ctl_est[k])
  weighted_log_RR_bcontrol <- numerator/denominator
  
  # calculate weighted RR with overall probabilities 
  numerator <- sum((p_y_est[j] + p_y_est[k]) * log_RR_le[j]) 
  denominator <- sum(p_y_est[j]) + sum(p_y_est[k])
  weighted_log_RR_boverall <- numerator/denominator
  
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
  numerator <- sum((p_y_est[j] + p_y_est[k]) * (PPO_leq[j] * (1-PPO_leq[j])) * log_RR_le[j])
  denominator <- sum((p_y_est[j] + p_y_est[k]) * (PPO_leq[j] * (1-PPO_leq[j])))
  weighted_RR_both_ov <- exp(numerator/denominator)
  weighted_log_RR_both_ov <- numerator/denominator
  
  ## for relative risk also look at using log_RR_ge 
  # fully bayesian approach 
  # calculate weighted RR with control probabilities 
  numerator <- sum((p_y_ctl_est[j] + p_y_ctl_est[k]) * log_RR_ge[k])
  denominator <- sum(p_y_ctl_est[j]) + sum(p_y_ctl_est[k])
  weighted_log_RR_bcontrol_ge <- numerator/denominator
  
  # calculate weighted RR with overall probabilities 
  numerator <- sum((p_y_est[j] + p_y_est[k]) * log_RR_ge[k]) 
  denominator <- sum(p_y_est[j]) + sum(p_y_est[k])
  weighted_log_RR_boverall_ge <- numerator/denominator
  
  # fully bayesian with only cumulative weights
  numerator <- sum((PPO_ctl_leq[j] * (1-PPO_ctl_leq[j])) * log_RR_ge[k])
  denominator <- sum((PPO_ctl_leq[j] * (1-PPO_ctl_leq[j])))
  weighted_RR_cum_ge <- exp(numerator/denominator)
  weighted_log_RR_cum_ge <- numerator/denominator
  
  # fully bayesian with both cumulative and split weights
  numerator <- sum((p_y_ctl_est[j] + p_y_ctl_est[k]) * (PPO_ctl_leq[j] * (1-PPO_ctl_leq[j])) * log_RR_ge[k])
  denominator <- sum((p_y_ctl_est[j] + p_y_ctl_est[k]) * (PPO_ctl_leq[j] * (1-PPO_ctl_leq[j])))
  weighted_RR_both_ge <- exp(numerator/denominator)
  weighted_log_RR_both_ge <- numerator/denominator
  
  # fully bayesian with both cumulative and split weights
  numerator <- sum((p_y_est[j] + p_y_est[k]) * (PPO_leq[j] * (1-PPO_leq[j])) * log_RR_ge[k])
  denominator <- sum((p_y_est[j] + p_y_est[k]) * (PPO_leq[j] * (1-PPO_leq[j])))
  weighted_RR_both_ov_ge <- exp(numerator/denominator)
  weighted_log_RR_both_ov_ge <- numerator/denominator
  
  unweighted_log_RR <- mean(log_RR_le)
  unweighted_RR <- exp(mean(log_RR_le))
  
  unweighted_log_RR_ge <- mean(log_RR_ge, na.rm = T)
  unweighted_RR_ge <- exp(mean(log_RR_ge, na.rm = T))
  
  
  # calculate NB 
  NB <- sum(p_y_trt_est[j] * PPO_ctl_geq[j+1]) - (sum(p_y_ctl_est[j] * PPO_trt_geq[j+1]))
  
  
   return(c(#p_y_ctl = p_y_ctl,
  #             p_y_trt = p_y_trt,
  #             RD_leq = RD_leq, 
  #             RD_geq = RD_geq, 
  #             OR = OR_le, 
  #             log_OR = log_OR_le,
              weighted_RD_control_PPO = weighted_RD, 
              weighted_RD_overall_PPO = weighted_RD_overall,
              weighted_RD_cum_PPO = weighted_RD_cum,
              weighted_RD_both_PPO = weighted_RD_both,
              weighted_RD_both_ov_PPO = weighted_RD_both_ov,
              unweighted_RD_PPO = unweighted_RD,
              weighted_OR_control_PPO = weighted_OR,
              weighted_OR_overall_PPO = weighted_OR_overall,
              weighted_log_OR_overall_PPO = weighted_log_OR_overall,
              weighted_log_OR_control_PPO = weighted_log_OR,
              unweighted_log_OR_PPO = unweighted_log_OR,
              unweighted_OR_PPO = unweighted_OR,
              weighted_log_RR_control_PPO = weighted_log_RR_control, 
              weighted_log_RR_overall_PPO = weighted_log_RR_overall,
              unweighted_log_RR_PPO = unweighted_log_RR,
              unweighted_RR_PPO = unweighted_RR,
              weighted_RD_bcontrol_PPO = weighted_RD_b, 
              weighted_RD_boverall_PPO = weighted_RD_boverall,
              weighted_OR_bcontrol_PPO = weighted_OR_b,
              weighted_OR_boverall_PPO = weighted_OR_boverall,
              weighted_OR_cum_PPO = weighted_OR_cum,
              weighted_OR_both_PPO = weighted_OR_both,
              weighted_OR_both_ov_PPO = weighted_OR_both_ov,
              weighted_log_OR_boverall_PPO = weighted_log_OR_boverall,
              weighted_log_OR_bcontrol_PPO = weighted_log_OR_b,
              weighted_log_OR_cum_PPO = weighted_log_OR_cum,
              weighted_log_OR_both_PPO = weighted_log_OR_both,
              weighted_log_OR_both_ov_PPO = weighted_log_OR_both_ov,
              weighted_log_RR_bcontrol_PPO = weighted_log_RR_bcontrol, 
              weighted_log_RR_boverall_PPO = weighted_log_RR_boverall,
              weighted_log_RR_cum_PPO = weighted_log_RR_cum,
              weighted_log_RR_both_PPO = weighted_log_RR_both,
              weighted_log_RR_both_ov_PPO = weighted_log_RR_both_ov,
              
              
              weighted_log_RR_bcontrol_ge_PPO = weighted_log_RR_bcontrol_ge, 
              weighted_log_RR_boverall_ge_PPO = weighted_log_RR_boverall_ge,
              weighted_log_RR_cum_ge_PPO = weighted_log_RR_cum_ge,
              weighted_log_RR_both_ge_PPO = weighted_log_RR_both_ge,
              weighted_log_RR_both_ov_ge_PPO = weighted_log_RR_both_ov_ge,
              unweighted_log_RR_ge_PPO = unweighted_log_RR_ge,
              unweighted_RR_ge_PPO = unweighted_RR_ge,
              #WR = WR, 
              NB_PPO = NB))
  
}


# calculate the credible interval 
# assumes a partial proportional odds model 
# assumes the treatment is coded as -0.5 for control and 0.5 for treatment 
# b is the model 
# NOTE: this only works for the iprior = 1 because I negated the alphas in rows 566, 567, 574, 576
order_weighted_RD2_CI <- function(b, n_post_draws = 4000, dat) { 
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
      PPO_ctl_leq[,(i - 1)] <-  1 - PPO_ctl_geq[,i]
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
  
  PPO_leq <- (PPO_ctl_leq + PPO_trt_leq)/2
  
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
  
  p_y_est <- (p_y_trt_est + p_y_ctl_est)/2
  
  # use data to properly weight each levels 
  p_y_ctl <- unname(table(dat$trt, dat$y)[1,])/sum(table(dat$trt, dat$y)[1,])
  p_y_trt <- unname(table(dat$trt, dat$y)[2,])/sum(table(dat$trt, dat$y)[2,])
  p_y <- as.vector(unname(table(dat$y))/sum(table(dat$y)))
  
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
  OR_ge_95CI <- apply(OR_ge[, 2:c], 2, quantile, probs = c(0.025, 0.975), na.rm = T)
  #for iprior = 1: (PPO_ctl_geq / (1- PPO_ctl_geq))/(PPO_trt_geq / (1 - PPO_trt_geq))
  
  log_RR_le <- log(PPO_trt_leq / PPO_ctl_leq)
  log_RR_le_95CI <- apply(log_RR_le[, 1:(c-1)], 2, quantile, probs = c(0.025, 0.975), na.rm = T)
  
  log_RR_ge <- log(PPO_ctl_geq / PPO_trt_geq)
  log_RR_ge_95CI <- apply(log_RR_ge[, 1:(c-1)], 2, quantile, probs = c(0.025, 0.975), na.rm = T)
  
  j <- 1:(c-1)
  k <- 2:c
  numerator <- rowSums((t(replicate(n_post_draws, p_y_ctl))[,j] + t(replicate(n_post_draws, p_y_ctl))[,k]) * 
                         (PPO_trt_leq[,j] - PPO_ctl_leq[,j]))
  denominator <- sum(p_y_ctl[j]) + sum(p_y_ctl[k])
  weighted_RD <- numerator/denominator
  weighted_RD_95CI <- quantile(weighted_RD, probs = c(0.025, 0.975))
  weighted_RD_p <- 1 - sum(weighted_RD >= 0)/n_post_draws
  weighted_RD_sd <- sd(weighted_RD)
  
  numerator <- rowSums((t(replicate(n_post_draws, p_y))[,j] + t(replicate(n_post_draws, p_y))[,k]) * 
                         (PPO_trt_leq[,j] - PPO_ctl_leq[,j]))
  denominator <- sum(p_y[j]) + sum(p_y[k])
  weighted_RD_overall <- numerator/denominator
  weighted_RD_95CI_overall <- quantile(weighted_RD_overall, probs = c(0.025, 0.975))
  weighted_RD_p_overall <- 1 - sum(weighted_RD_overall >= 0)/n_post_draws
  weighted_RD_sd_overall <- sd(weighted_RD_overall)
  
  # fully bayesian approach
  numerator <- rowSums((p_y_ctl_est[,j] + p_y_ctl_est[,k]) * 
                         (PPO_trt_leq[,j] - PPO_ctl_leq[,j]))
  denominator <- rowSums(p_y_ctl_est[,j]) + rowSums(p_y_ctl_est[,k])
  weighted_RD_b <- numerator/denominator
  weighted_RD_95CI_b <- quantile(weighted_RD_b, probs = c(0.025, 0.975))
  weighted_RD_p_b <- 1 - sum(weighted_RD_b >= 0)/n_post_draws
  weighted_RD_sd_b <- sd(weighted_RD_b)
  
  numerator <- rowSums((p_y_est[,j] + p_y_est[,k]) * 
                         (PPO_trt_leq[,j] - PPO_ctl_leq[,j]))
  denominator <- rowSums(p_y_est[,j]) + rowSums(p_y_est[,k])
  weighted_RD_boverall <- numerator/denominator
  weighted_RD_95CI_boverall <- quantile(weighted_RD_boverall, probs = c(0.025, 0.975))
  weighted_RD_p_boverall <- 1 - sum(weighted_RD_boverall >= 0)/n_post_draws
  weighted_RD_sd_boverall <- sd(weighted_RD_boverall)
  
  # cumulative weights
  numerator <- rowSums((PPO_ctl_leq[,j] * (1-PPO_ctl_leq[,j])) * (PPO_trt_leq[,j] - PPO_ctl_leq[,j]))
  denominator <- rowSums((PPO_ctl_leq[,j] * (1-PPO_ctl_leq[,j])))
  weighted_RD_cum <- numerator/denominator
  weighted_RD_cum95CI <- quantile(weighted_RD_cum, probs = c(0.025, 0.975))
  weighted_RD_cump <- 1 - sum(weighted_RD_cum >= 0)/n_post_draws
  weighted_RD_cumsd <- sd(weighted_RD_cum)
  
  # both weights
  numerator <- rowSums((p_y_ctl_est[,j] + p_y_ctl_est[,k]) * (PPO_ctl_leq[,j] * (1-PPO_ctl_leq[,j])) * (PPO_trt_leq[,j] - PPO_ctl_leq[,j]))
  denominator <- rowSums((p_y_ctl_est[,j] + p_y_ctl_est[,k]) * (PPO_ctl_leq[,j] * (1-PPO_ctl_leq[,j])))
  weighted_RD_both <- numerator/denominator
  weighted_RD_both95CI <- quantile(weighted_RD_both, probs = c(0.025, 0.975))
  weighted_RD_bothp <- 1 - sum(weighted_RD_both >= 0)/n_post_draws
  weighted_RD_bothsd <- sd(weighted_RD_both)
  
  # both weights
  numerator <- rowSums((p_y_est[,j] + p_y_est[,k]) * (PPO_leq[,j] * (1-PPO_leq[,j])) * (PPO_trt_leq[,j] - PPO_ctl_leq[,j]))
  denominator <- rowSums((p_y_est[,j] + p_y_est[,k]) * (PPO_leq[,j] * (1-PPO_leq[,j])))
  weighted_RD_both_ov <- numerator/denominator
  weighted_RD_both_ov95CI <- quantile(weighted_RD_both_ov, probs = c(0.025, 0.975))
  weighted_RD_both_ovp <- 1 - sum(weighted_RD_both_ov >= 0)/n_post_draws
  weighted_RD_both_ovsd <- sd(weighted_RD_both_ov)
  
  #unweighted RD
  unweighted_RD <- rowMeans(RD_leq, na.rm = T)
  unweighted_RD_95CI <- quantile(unweighted_RD, probs = c(0.025, 0.975))
  unweighted_RD_p <- 1 - sum(unweighted_RD >= 0)/n_post_draws
  unweighted_RD_sd <- sd(unweighted_RD)
  
  #calculate weighted OR 
  numerator <- rowSums((t(replicate(n_post_draws, p_y_ctl))[,j] + t(replicate(n_post_draws, p_y_ctl))[,k]) * 
                         log(OR_le[,j]))
  denominator <- sum(p_y_ctl[j]) + sum(p_y_ctl[k])
  weighted_OR <- exp(numerator/denominator)
  weighted_OR_95CI <- quantile(weighted_OR, probs = c(0.025, 0.975))
  weighted_OR_p <- 1 - sum(weighted_OR >= 1)/n_post_draws
  weighted_OR_sd <- sd(weighted_OR)
  
  weighted_log_OR <- numerator/denominator
  weighted_log_OR_95CI <- quantile(weighted_log_OR, probs = c(0.025, 0.975))
  weighted_log_OR_p <- 1 - sum(weighted_log_OR >= 0)/n_post_draws
  weighted_log_OR_sd <- sd(weighted_log_OR)
  
  #calculate weighted OR with overall probability as weights 
  numerator <- rowSums((t(replicate(n_post_draws, p_y))[,j] + t(replicate(n_post_draws, p_y))[,k]) * 
                         log(OR_le[,j]))
  denominator <- sum(p_y[j]) + sum(p_y[k])
  weighted_OR_overall <- exp(numerator/denominator)
  weighted_OR_95CI_overall <- quantile(weighted_OR_overall, probs = c(0.025, 0.975))
  weighted_OR_p_overall <- 1 - sum(weighted_OR_overall >= 1)/n_post_draws
  weighted_OR_sd_overall <- sd(weighted_OR_overall)

  weighted_log_OR_overall <- numerator/denominator
  weighted_log_OR_95CI_overall <- quantile(weighted_log_OR_overall, probs = c(0.025, 0.975))
  weighted_log_OR_p_overall <- 1 - sum(weighted_log_OR_overall >= 0)/n_post_draws
  weighted_log_OR_sd_overall <- sd(weighted_log_OR_overall)
  
  # fully bayesian approach 
  #calculate weighted OR 
  numerator <- rowSums((p_y_ctl_est[,j] + p_y_ctl_est[,k]) * 
                         log(OR_le[,j]))
  denominator <- rowSums(p_y_ctl_est[,j]) + rowSums(p_y_ctl_est[,k])
  weighted_OR_b <- exp(numerator/denominator)
  weighted_OR_95CI_b <- quantile(weighted_OR_b, probs = c(0.025, 0.975))
  weighted_OR_p_b <- 1 - sum(weighted_OR_b >= 1)/n_post_draws
  weighted_OR_sd_b <- sd(weighted_OR_b)
  
  weighted_log_OR_b <- numerator/denominator
  weighted_log_OR_95CI_b <- quantile(weighted_log_OR_b, probs = c(0.025, 0.975))
  weighted_log_OR_p_b <- 1 - sum(weighted_log_OR_b >= 0)/n_post_draws
  weighted_log_OR_sd_b <- sd(weighted_log_OR_b)
  
  
  #calculate weighted OR with overall probability as weights 
  numerator <- rowSums((p_y_est[,j] + p_y_est[,k]) * 
                         log(OR_le[,j]))
  denominator <- rowSums(p_y_est[,j]) + rowSums(p_y_est[,k])
  weighted_OR_boverall <- exp(numerator/denominator)
  weighted_OR_95CI_boverall <- quantile(weighted_OR_boverall, probs = c(0.025, 0.975))
  weighted_OR_p_boverall <- 1 - sum(weighted_OR_boverall >= 1)/n_post_draws
  weighted_OR_sd_boverall <- sd(weighted_OR_boverall)
  
  weighted_log_OR_boverall <- numerator/denominator
  weighted_log_OR_95CI_boverall <- quantile(weighted_log_OR_boverall, probs = c(0.025, 0.975))
  weighted_log_OR_p_boverall <- 1 - sum(weighted_log_OR_boverall >= 0)/n_post_draws
  weighted_log_OR_sd_boverall <- sd(weighted_log_OR_boverall)
  
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
  
  # both weights with overall probability
  numerator <- rowSums((p_y_est[,j] + p_y_est[,k]) * (PPO_leq[,j] * (1-PPO_leq[,j])) * log(OR_le[,j]))
  denominator <- rowSums((p_y_est[,j] + p_y_est[,k]) * (PPO_leq[,j] * (1-PPO_leq[,j])))
  weighted_OR_both_ov <- exp(numerator/denominator)
  weighted_OR_both_ov95CI <- quantile(weighted_OR_both_ov, probs = c(0.025, 0.975))
  weighted_OR_both_ovp <- 1 - sum(weighted_OR_both_ov >= 0)/n_post_draws
  weighted_OR_both_ovsd <- sd(weighted_OR_both_ov)
  
  weighted_log_OR_both_ov <- numerator/denominator
  weighted_log_OR_both_ov95CI <- quantile(weighted_log_OR_both_ov, probs = c(0.025, 0.975))
  weighted_log_OR_both_ovp <- 1 - sum(weighted_log_OR_both_ov >= 0)/n_post_draws
  weighted_log_OR_both_ovsd <- sd(weighted_log_OR_both_ov)
  
  #unweighted OR 
  unweighted_OR <- rowMeans(log_OR_le, na.rm = T)
  unweighted_log_OR_95CI <- quantile(unweighted_OR, probs = c(0.025, 0.975))
  unweighted_log_OR_p <- 1 - sum(unweighted_OR >= 0)/n_post_draws
  unweighted_log_OR_sd <- sd(unweighted_OR)
  
  
  unweighted_OR_95CI <- quantile(exp(unweighted_OR), probs = c(0.025, 0.975))
  unweighted_OR_p <- 1 - sum(exp(unweighted_OR) >= 0)/n_post_draws
  unweighted_OR_sd <- sd(exp(unweighted_OR))
  
  # calculate weighted RR with control probabilities 
  numerator <- rowSums((t(replicate(n_post_draws, p_y_ctl))[,j] + t(replicate(n_post_draws, p_y_ctl))[,k]) * 
                         log_RR_le[,j])
  denominator <- sum(p_y_ctl[j]) + sum(p_y_ctl[k])
  weighted_log_RR_control <- numerator/denominator
  weighted_log_RR_95CI_control <- quantile(weighted_log_RR_control, probs = c(0.025, 0.975))
  weighted_log_RR_p_control <- 1 - sum(weighted_log_RR_control >= 0)/n_post_draws
  weighted_log_RR_sd_control <- sd(weighted_log_RR_control)
  
  # calculate weighted RR with overall probabilities 
  numerator <- rowSums((t(replicate(n_post_draws, p_y))[,j] + t(replicate(n_post_draws, p_y))[,k]) * 
                         log_RR_le[,j])
  denominator <- sum(p_y[j]) + sum(p_y[k])
  weighted_log_RR_overall <- numerator/denominator
  weighted_log_RR_95CI_overall <- quantile(weighted_log_RR_overall, probs = c(0.025, 0.975))
  weighted_log_RR_p_overall <- 1 - sum(weighted_log_RR_overall >= 0)/n_post_draws
  weighted_log_RR_sd_overall <- sd(weighted_log_RR_overall)
  
  # fully bayesian 
  # calculate weighted RR with control probabilities 
  numerator <- rowSums((p_y_ctl_est[,j] + p_y_ctl_est[,k]) * 
                         log_RR_le[,j])
  denominator <- rowSums(p_y_ctl_est[,j]) + rowSums(p_y_ctl_est[,k])
  weighted_log_RR_bcontrol <- numerator/denominator
  weighted_log_RR_95CI_bcontrol <- quantile(weighted_log_RR_bcontrol, probs = c(0.025, 0.975))
  weighted_log_RR_p_bcontrol <- 1 - sum(weighted_log_RR_bcontrol >= 0)/n_post_draws
  weighted_log_RR_sd_bcontrol <- sd(weighted_log_RR_bcontrol)
  
  # calculate weighted RR with overall probabilities 
  numerator <- rowSums((p_y_est[,j] + p_y_est[,k]) * 
                         log_RR_le[,j])
  denominator <- rowSums(p_y_est[,j]) + rowSums(p_y_est[,k])
  weighted_log_RR_boverall <- numerator/denominator
  weighted_log_RR_95CI_boverall <- quantile(weighted_log_RR_boverall, probs = c(0.025, 0.975))
  weighted_log_RR_p_boverall <- 1 - sum(weighted_log_RR_boverall >= 0)/n_post_draws
  weighted_log_RR_sd_boverall <- sd(weighted_log_RR_boverall)
  
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
  
  # both weights
  numerator <- rowSums((p_y_est[,j] + p_y_est[,k]) * (PPO_leq[,j] * (1-PPO_leq[,j])) * log_RR_le[,j])
  denominator <- rowSums((p_y_est[,j] + p_y_est[,k]) * (PPO_leq[,j] * (1-PPO_leq[,j])))
  weighted_log_RR_both_ov <- numerator/denominator
  weighted_log_RR_both_ov95CI <- quantile(weighted_log_RR_both_ov, probs = c(0.025, 0.975))
  weighted_log_RR_both_ovp <- 1 - sum(weighted_log_RR_both_ov >= 0)/n_post_draws
  weighted_log_RR_both_ovsd <- sd(weighted_log_RR_both_ov)
  
  #unweighted RR
  unweighted_RR <- rowMeans(log_RR_le, na.rm = T)
  unweighted_log_RR_95CI <- quantile(unweighted_RR, probs = c(0.025, 0.975))
  unweighted_log_RR_p <- 1 - sum(unweighted_RR >= 0)/n_post_draws
  unweighted_log_RR_sd <- sd(unweighted_RR)
  
  
  unweighted_RR_95CI <- quantile(exp(unweighted_RR), probs = c(0.025, 0.975))
  unweighted_RR_p <- 1 - sum(exp(unweighted_RR) >= 0)/n_post_draws
  unweighted_RR_sd <- sd(exp(unweighted_RR))
  
  # fully bayesian RR greater than or equal to 
  # calculate weighted RR with control probabilities 
  numerator <- rowSums((p_y_ctl_est[,j] + p_y_ctl_est[,k]) * 
                         log_RR_ge[,k])
  denominator <- rowSums(p_y_ctl_est[,j]) + rowSums(p_y_ctl_est[,k])
  weighted_log_RR_ge_bcontrol <- numerator/denominator
  weighted_log_RR_ge_95CI_bcontrol <- quantile(weighted_log_RR_ge_bcontrol, probs = c(0.025, 0.975))
  weighted_log_RR_ge_p_bcontrol <- 1 - sum(weighted_log_RR_ge_bcontrol >= 0)/n_post_draws
  weighted_log_RR_ge_sd_bcontrol <- sd(weighted_log_RR_ge_bcontrol)
  
  # calculate weighted RR with overall probabilities 
  numerator <- rowSums((p_y_est[,j] + p_y_est[,k]) * 
                         log_RR_ge[,k])
  denominator <- rowSums(p_y_est[,j]) + rowSums(p_y_est[,k])
  weighted_log_RR_ge_boverall <- numerator/denominator
  weighted_log_RR_ge_95CI_boverall <- quantile(weighted_log_RR_ge_boverall, probs = c(0.025, 0.975))
  weighted_log_RR_ge_p_boverall <- 1 - sum(weighted_log_RR_ge_boverall >= 0)/n_post_draws
  weighted_log_RR_ge_sd_boverall <- sd(weighted_log_RR_ge_boverall)
  
  # cumulative weights only
  numerator <- rowSums((PPO_ctl_leq[,j] * (1-PPO_ctl_leq[,j])) * log_RR_ge[,k])
  denominator <- rowSums((PPO_ctl_leq[,j] * (1-PPO_ctl_leq[,j])))
  weighted_log_RR_ge_cum <- numerator/denominator
  weighted_log_RR_ge_cum95CI <- quantile(weighted_log_RR_ge_cum, probs = c(0.025, 0.975))
  weighted_log_RR_ge_cump <- 1 - sum(weighted_log_RR_ge_cum >= 0)/n_post_draws
  weighted_log_RR_ge_cumsd <- sd(weighted_log_RR_ge_cum)
  
  # both weights
  numerator <- rowSums((p_y_ctl_est[,j] + p_y_ctl_est[,k]) * (PPO_ctl_leq[,j] * (1-PPO_ctl_leq[,j])) * log_RR_ge[,k])
  denominator <- rowSums((p_y_ctl_est[,j] + p_y_ctl_est[,k]) * (PPO_ctl_leq[,j] * (1-PPO_ctl_leq[,j])))
  weighted_log_RR_ge_both <- numerator/denominator
  weighted_log_RR_ge_both95CI <- quantile(weighted_log_RR_ge_both, probs = c(0.025, 0.975))
  weighted_log_RR_ge_bothp <- 1 - sum(weighted_log_RR_ge_both >= 0)/n_post_draws
  weighted_log_RR_ge_bothsd <- sd(weighted_log_RR_ge_both)
  
  # both weights overall weights
  numerator <- rowSums((p_y_est[,j] + p_y_est[,k]) * (PPO_leq[,j] * (1-PPO_leq[,j])) * log_RR_ge[,k])
  denominator <- rowSums((p_y_est[,j] + p_y_est[,k]) * (PPO_leq[,j] * (1-PPO_leq[,j])))
  weighted_log_RR_ge_both_ov <- numerator/denominator
  weighted_log_RR_ge_both_ov95CI <- quantile(weighted_log_RR_ge_both_ov, probs = c(0.025, 0.975))
  weighted_log_RR_ge_both_ovp <- 1 - sum(weighted_log_RR_ge_both_ov >= 0)/n_post_draws
  weighted_log_RR_ge_both_ovsd <- sd(weighted_log_RR_ge_both_ov)
  
  #unweighted RR
  unweighted_RR_ge <- rowMeans(log_RR_ge, na.rm = T)
  unweighted_log_RR_ge_95CI <- quantile(unweighted_RR_ge, probs = c(0.025, 0.975))
  unweighted_log_RR_ge_p <- 1 - sum(unweighted_RR_ge >= 0)/n_post_draws
  unweighted_log_RR_ge_sd <- sd(unweighted_RR_ge)
  
  # calculate NB 
  NB <- rowSums(p_y_trt_est[,j] * PPO_ctl_geq[,(j+1)]) - (rowSums(p_y_ctl_est[,j] * PPO_trt_geq[,(j+1)]))
  NB_95CI <- quantile(NB, probs = c(0.025, 0.975))
  NB_p <- 1 - sum(NB >= 0)/n_post_draws
  NB_sd <- sd(NB)
  
  return(c(#p_y_ctl_CI = p_y_ctl_CI, 
              #p_y_trt_CI = p_y_trt_CI, 
              #RD_leq = RD_leq_95CI, 
              #RD_geq = RD_geq_95CI, 
              #OR = OR_le_95CI, 
              #log_OR = log_OR_le_95CI,
              weighted_RD_control_CI_PPO = weighted_RD_95CI,
              weighted_RD_p_control_PPO = weighted_RD_p,
              weighted_RD_sd_PPO = weighted_RD_sd, 
              weighted_RD_sd_overall_PPO = weighted_RD_sd_overall,
              weighted_RD_overall_CI_PPO = weighted_RD_95CI_overall,
              weighted_RD_p_overall_PPO = weighted_RD_p_overall,
              weighted_OR_CI_PPO = weighted_OR_95CI, 
              weighted_OR_p_PPO = weighted_OR_p,
              weighted_OR_sd_PPO = weighted_OR_sd,
              weighted_OR_overall_CI_PPO = weighted_OR_95CI_overall,
              weighted_OR_p_overall_PPO = weighted_OR_p_overall,
              weighted_OR_sd_overall_PPO = weighted_OR_sd_overall,
              weighted_log_OR_control_CI_PPO = weighted_log_OR_95CI,
              weighted_log_OR_p_control_PPO = weighted_log_OR_p,
              weighted_log_OR_sd_control_PPO = weighted_log_OR_sd,
              weighted_log_OR_overall_CI_PPO = weighted_log_OR_95CI_overall,
              weighted_log_OR_p_overall_PPO = weighted_log_OR_p_overall,
              weighted_log_OR_sd_overall_PPO = weighted_log_OR_sd_overall,
              weighted_log_RR_95CI_control_PPO = weighted_log_RR_95CI_control,
              weighted_log_RR_p_control_PPO = weighted_log_RR_p_control, 
              weighted_log_RR_sd_control_PPO = weighted_log_RR_sd_control,
              weighted_log_RR_95CI_overall_PPO = weighted_log_RR_95CI_overall,
              weighted_log_RR_p_overall_PPO = weighted_log_RR_p_overall, 
              weighted_log_RR_sd_overall_PPO = weighted_log_RR_sd_overall,
              
              weighted_RD_bcontrol_CI_PPO = weighted_RD_95CI_b,
              weighted_RD_p_bcontrol_PPO = weighted_RD_p_b,
              weighted_RD_sd_b_PPO = weighted_RD_sd_b, 
              weighted_RD_sd_boverall_PPO = weighted_RD_sd_boverall,
              weighted_RD_boverall_CI_PPO = weighted_RD_95CI_boverall,
              weighted_RD_p_boverall_PPO = weighted_RD_p_boverall,
              weighted_OR_CI_b_PPO = weighted_OR_95CI_b, 
              weighted_OR_p_b_PPO = weighted_OR_p_b,
              weighted_OR_sd_b_PPO = weighted_OR_sd_b,
              weighted_OR_boverall_CI_PPO = weighted_OR_95CI_boverall,
              weighted_OR_p_boverall_PPO = weighted_OR_p_boverall,
              weighted_OR_sd_boverall_PPO = weighted_OR_sd_boverall,
              weighted_log_OR_bcontrol_CI_PPO = weighted_log_OR_95CI_b,
              weighted_log_OR_p_bcontrol_PPO = weighted_log_OR_p_b,
              weighted_log_OR_sd_bcontrol_PPO = weighted_log_OR_sd_b,
              weighted_log_OR_boverall_CI_PPO = weighted_log_OR_95CI_boverall,
              weighted_log_OR_p_boverall_PPO = weighted_log_OR_p_boverall,
              weighted_log_OR_sd_boverall_PPO = weighted_log_OR_sd_boverall,
              weighted_log_RR_95CI_bcontrol_PPO = weighted_log_RR_95CI_bcontrol,
              weighted_log_RR_p_bcontrol_PPO = weighted_log_RR_p_bcontrol, 
              weighted_log_RR_sd_bcontrol_PPO = weighted_log_RR_sd_bcontrol,
              weighted_log_RR_95CI_boverall_PPO = weighted_log_RR_95CI_boverall,
              weighted_log_RR_p_boverall_PPO = weighted_log_RR_p_boverall, 
              weighted_log_RR_sd_boverall_PPO = weighted_log_RR_sd_boverall,
              weighted_log_RR_ge_95CI_bcontrol_PPO = weighted_log_RR_ge_95CI_bcontrol,
              weighted_log_RR_ge_p_bcontrol_PPO = weighted_log_RR_ge_p_bcontrol, 
              weighted_log_RR_ge_sd_bcontrol_PPO = weighted_log_RR_ge_sd_bcontrol,
              weighted_log_RR_ge_95CI_boverall_PPO = weighted_log_RR_ge_95CI_boverall,
              weighted_log_RR_ge_p_boverall_PPO = weighted_log_RR_ge_p_boverall, 
              weighted_log_RR_ge_sd_boverall_PPO = weighted_log_RR_ge_sd_boverall,
              
              unweighted_RD_CI_PPO = unweighted_RD_95CI,
              unweighted_RD_p_PPO = unweighted_RD_p,
              unweighted_RD_sd_PPO = unweighted_RD_sd, 
              unweighted_OR_CI_PPO = unweighted_OR_95CI,
              unweighted_OR_p_PPO = unweighted_OR_p,
              unweighted_OR_sd_PPO = unweighted_OR_sd,  
              unweighted_log_OR_CI_PPO = unweighted_log_OR_95CI,
              unweighted_log_OR_p_PPO = unweighted_log_OR_p,
              unweighted_log_OR_sd_PPO = unweighted_log_OR_sd,
              unweighted_RR_CI_PPO = unweighted_RR_95CI,
              unweighted_RR_p_PPO = unweighted_RR_p,
              unweighted_RR_sd_PPO = unweighted_RR_sd,  
              unweighted_log_RR_CI_PPO = unweighted_log_RR_95CI,
              unweighted_log_RR_p_PPO = unweighted_log_RR_p,
              unweighted_log_RR_sd_PPO = unweighted_log_RR_sd,
              unweighted_log_RR_ge_CI_PPO = unweighted_log_RR_ge_95CI,
              unweighted_log_RR_ge_p_PPO = unweighted_log_RR_ge_p,
              unweighted_log_RR_ge_sd_PPO = unweighted_log_RR_ge_sd,
              
              
              weighted_RD_cum_CI_PPO = weighted_RD_cum95CI,
              weighted_RD_cum_p_PPO = weighted_RD_cump,
              weighted_RD_cum_sd_PPO = weighted_RD_cumsd,
              weighted_RD_both_CI_PPO = weighted_RD_both95CI,
              weighted_RD_both_p_PPO = weighted_RD_bothp,
              weighted_RD_both_sd_PPO = weighted_RD_bothsd, 
              weighted_RD_both_ov_CI_PPO = weighted_RD_both_ov95CI,
              weighted_RD_both_ov_p_PPO = weighted_RD_both_ovp,
              weighted_RD_both_ov_sd_PPO = weighted_RD_both_ovsd, 
              weighted_log_OR_cum_CI_PPO = weighted_log_OR_cum95CI,
              weighted_log_OR_cum_p_PPO = weighted_log_OR_cump,
              weighted_log_OR_cum_sd_PPO = weighted_log_OR_cumsd,
              weighted_log_OR_both_CI_PPO = weighted_log_OR_both95CI,
              weighted_log_OR_both_p_PPO = weighted_log_OR_bothp,
              weighted_log_OR_both_sd_PPO = weighted_log_OR_bothsd,
              weighted_log_OR_both_ov_CI_PPO = weighted_log_OR_both_ov95CI,
              weighted_log_OR_both_ov_p_PPO = weighted_log_OR_both_ovp,
              weighted_log_OR_both_ov_sd_PPO = weighted_log_OR_both_ovsd,
              weighted_log_RR_cum_CI_PPO = weighted_log_RR_cum95CI,
              weighted_log_RR_cum_p_PPO = weighted_log_RR_cump,
              weighted_log_RR_cum_sd_PPO = weighted_log_RR_cumsd,
              weighted_log_RR_both_CI_PPO = weighted_log_RR_both95CI,
              weighted_log_RR_both_p_PPO = weighted_log_RR_bothp,
              weighted_log_RR_both_sd_PPO = weighted_log_RR_bothsd,
              weighted_log_RR_both_ov_CI_PPO = weighted_log_RR_both_ov95CI,
              weighted_log_RR_both_ov_p_PPO = weighted_log_RR_both_ovp,
              weighted_log_RR_both_ov_sd_PPO = weighted_log_RR_both_ovsd,
              weighted_log_RR_ge_cum_CI_PPO = weighted_log_RR_ge_cum95CI,
              weighted_log_RR_ge_cum_p_PPO = weighted_log_RR_ge_cump,
              weighted_log_RR_ge_cum_sd_PPO = weighted_log_RR_ge_cumsd,
              weighted_log_RR_ge_both_CI_PPO = weighted_log_RR_ge_both95CI,
              weighted_log_RR_ge_both_p_PPO = weighted_log_RR_ge_bothp,
              weighted_log_RR_ge_both_sd_PPO = weighted_log_RR_ge_bothsd,
              weighted_log_RR_ge_both_ov_CI_PPO = weighted_log_RR_ge_both_ov95CI,
              weighted_log_RR_ge_both_ov_p_PPO = weighted_log_RR_ge_both_ovp,
              weighted_log_RR_ge_both_ov_sd_PPO = weighted_log_RR_ge_both_ovsd,
              
              
              
              #WR = WR_95CI,
              NB_CI_PPO = NB_95CI,
              NB_p_PPO = NB_p,
              NB_sd_PPO = NB_sd))
  
}


# new function for weighted risk difference 
# assumes a proportional odds model 
# assumes the treatment is coded as -0.5 for control and 0.5 for treatment 
order_weighted_RD2_PO <- function(b = b, dat) { 
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
  
  # get weights from the data 
  p_y_ctl <- unname(table(dat$trt, dat$y)[1,])/500
  p_y_trt <- unname(table(dat$trt, dat$y)[2,])/500
  p_y <- unname(table(dat$y))/1000
  
  
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
  
  j <- 1:(c-1)
  k <- 2:c
  numerator <- sum((p_y_ctl[j] + p_y_ctl[k]) * (PO_trt_leq[j] - PO_ctl_leq[j]))
  denominator <- sum(p_y_ctl[j]) + sum(p_y_ctl[k])
  weighted_RD <- numerator/denominator
  
  #calculated weighted RD with overall probability as weights
  p_y <- (p_y_trt + p_y_ctl)/2
  numerator <- sum((p_y[j] + p_y[k]) * (PO_trt_leq[j] - PO_ctl_leq[j])) 
  denominator <- sum(p_y[j]) + sum(p_y[k])
  weighted_RD_overall <- numerator/denominator
  
  # calculate weighted RR with control probabilities 
  numerator <- sum((p_y_ctl[j] + p_y_ctl[k]) * log_RR_le[j]) 
  denominator <- sum(p_y_ctl[j]) + sum(p_y_ctl[k])
  weighted_log_RR_control <- numerator/denominator
  
  # calculate weighted RR with overall probabilities 
  numerator <- sum((p_y[j] + p_y[k]) * log_RR_le[j])
  denominator <- sum(p_y[j]) + sum(p_y[k])
  weighted_log_RR_overall <- numerator/denominator
  
  # calculate NB 
  NB <- sum(p_y_trt_est[j] * PO_ctl_geq[j+1]) - (sum(p_y_ctl_est[j] * PO_trt_geq[j+1]))
  
  return(c(#p_y_ctl = p_y_ctl, 
              #p_y_trt = p_y_trt,
              #RD_leq = RD_leq, 
              #RD_geq = RD_geq, 
              OR_PO = OR[1],
              log_OR_PO = log_OR[1],
              weighted_RD_control_PO = weighted_RD, 
              weighted_RD_overall_PO = weighted_RD_overall,
              weighted_log_RR_control_PO = weighted_log_RR_control,
              weighted_log_RR_overall_PO = weighted_log_RR_control,
              unweighted_log_RR = unweighted_log_RR,
              unweighted_RR = unweighted_RR, 
              unweighted_RD = unweighted_RD,
              unweighted_log_OR = unweighted_log_OR,
              unweighted_OR = unweighted_OR, 
              #WR = WR,
              NB_PO = NB))
  
}

# calculate the credible interval 
# assumes a partial proportional odds model 
# assumes the treatment is coded as -0.5 for control and 0.5 for treatment 
# b is the model 
order_weighted_RD2_PO_CI <- function(b, n_post_draws = 4000, dat) { 
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
  
  # use data to properly weight each levels 
  p_y_ctl <- unname(table(dat$trt, dat$y)[1,])/500
  p_y_trt <- unname(table(dat$trt, dat$y)[2,])/500
  p_y <- as.vector(unname(table(dat$y))/1000)
  
  RD_leq <- PO_trt_leq - PO_ctl_leq
  RD_leq_95CI <- apply(RD_leq[,], 2, quantile, probs = c(0.025, 0.975), na.rm = T)
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
  
  log_RR_le <- log(PO_trt_leq / PO_ctl_leq)
  log_RR_le_95CI <- apply(log_RR_le[, 1:3], 2, quantile, probs = c(0.025, 0.975), na.rm = T)
  
  log_RR_ge <- log(PO_ctl_geq / PO_trt_geq)
  log_RR_ge_95CI <- apply(log_RR_ge[, 1:3], 2, quantile, probs = c(0.025, 0.975), na.rm = T)
  
  
  # calculate weighted RD with control probability as weights
  j <- 1:(c-1)
  k <- 2:c
  numerator <- rowSums((t(replicate(n_post_draws, p_y_ctl))[,j] + t(replicate(n_post_draws, p_y_ctl))[,k]) *
                         (PO_trt_leq[,j] - PO_ctl_leq[,j])) 
  denominator <- sum(p_y_ctl[j]) + sum(p_y_ctl[k])
  weighted_RD <- numerator/denominator
  weighted_RD_95CI <- quantile(weighted_RD, probs = c(0.025, 0.975))
  weighted_RD_p <- 1 - sum(weighted_RD >= 0)/n_post_draws
  weighted_RD_sd <- sd(weighted_RD)
  
  #calculate weighted risk difference with overall probability as weights
  numerator <- rowSums((t(replicate(n_post_draws, p_y))[,j] + t(replicate(n_post_draws, p_y))[,k]) * 
                         (PO_trt_leq[,j] - PO_ctl_leq[,j]))
  denominator <- sum(p_y[j]) + sum(p_y[k])
  weighted_RD_overall <- numerator/denominator
  weighted_RD_95CI_overall <- quantile(weighted_RD_overall, probs = c(0.025, 0.975))
  weighted_RD_p_overall <- 1 - sum(weighted_RD_overall >= 0)/n_post_draws
  weighted_RD_sd_overall <- sd(weighted_RD_overall)
  
  # calculate weighted RR with control probabilities 
  numerator <- rowSums((t(replicate(n_post_draws, p_y_ctl))[,j] + t(replicate(n_post_draws, p_y_ctl))[,k]) * 
                         log_RR_le[,j]) 
  denominator <- sum(p_y_ctl[j]) + sum(p_y_ctl[k])
  weighted_log_RR_control <- numerator/denominator
  weighted_log_RR_95CI_control <- quantile(weighted_log_RR_control, probs = c(0.025, 0.975))
  weighted_log_RR_p_control <- 1 - sum(weighted_log_RR_control >= 0)/n_post_draws
  weighted_log_RR_sd_control <- sd(weighted_log_RR_control)
  
  # calculate weighted RR with overall probabilities 
  numerator <- rowSums((t(replicate(n_post_draws, p_y))[,j] + t(replicate(n_post_draws, p_y))[,k]) * 
                         log_RR_le[,j])
  denominator <- sum(p_y[j]) + sum(p_y[k])
  weighted_log_RR_overall <- numerator/denominator
  weighted_log_RR_95CI_overall <- quantile(weighted_log_RR_overall, probs = c(0.025, 0.975))
  weighted_log_RR_p_overall <- 1 - sum(weighted_log_RR_overall >= 0)/n_post_draws
  weighted_log_RR_sd_overall <- sd(weighted_log_RR_overall)

  
  # calculate NB 
  NB <- rowSums(p_y_trt_est[j] * PO_ctl_geq[,(j+1)]) - (rowSums(p_y_ctl_est[j] * PO_trt_geq[,(j+1)]))
  NB_95CI <- quantile(NB, probs = c(0.025, 0.975))
  NB_p <- 1 - sum(NB >= 0)/n_post_draws 
  NB_sd <- sd(NB)
  
  return(c(#RD_leq = RD_leq_95CI, 
              #RD_geq = RD_geq_95CI, 
              OR_CI_PO = OR_95CI, 
              OR_p_PO = OR_p,
              OR_sd_PO = OR_sd,
              log_OR_CI_PO = log_OR_95CI, 
              log_OR_p_PO = log_OR_p,
              log_OR_sd_PO = log_OR_sd,
              weighted_RD_control_CI_PO = weighted_RD_95CI, 
              weighted_RD_p_control_PO = weighted_RD_p,
              weighted_RD_sd_PO = weighted_RD_sd,
              weighted_RD_overall_CI_PO = weighted_RD_95CI_overall, 
              weighted_RD_p_overall_PO = weighted_RD_p_overall,
              weighted_RD_sd_overall_PO = weighted_RD_sd_overall,
              weighted_log_RR_95CI_control_PO = weighted_log_RR_95CI_control,
              weighted_log_RR_p_control_PO = weighted_log_RR_p_control, 
              weighted_log_RR_sd_control_PO = weighted_log_RR_sd_control,
              weighted_log_RR_95CI_overall_PO = weighted_log_RR_95CI_overall,
              weighted_log_RR_p_overall_PO = weighted_log_RR_p_overall, 
              weighted_log_RR_sd_overall_PO = weighted_log_RR_sd_overall,
              #WR = WR_95CI, 
              NB_CI_PO = NB_95CI, 
              NB_p_PO = NB_p, 
              NB_sd_PO = NB_sd))
  
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
      p_y_ctl[c] <- 1 - sum(p_y_ctl[1:(c-1)])
    } else {
      p_y_ctl[i] <-  p_ctl_less[i] - sum(p_y_ctl[1:(i-1)])
    }
  }
  
  RD_leq <- p_trt_less - p_ctl_less
  RD_geq <- p_ctl_more - p_trt_more
  
  OR <- (p_trt_less / (1- p_trt_less)) / (p_ctl_less / (1- p_ctl_less))
  
  j <- 1:(c-1)
  k <- 2:c
  numerator <- sum(p_y_ctl[j] * (p_trt_less[j] - p_ctl_less[j])) + sum(p_y_ctl[k] * (p_ctl_more[k] - p_trt_more[k]))
  denominator <- sum(p_y_ctl[j]) + sum(p_y_ctl[k])
  weighted_RD <- numerator/denominator
  return(list(RD_leq = RD_leq, RD_geq = RD_geq, OR = OR, weighted_RD = weighted_RD))
  
}