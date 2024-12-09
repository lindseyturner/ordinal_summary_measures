# array simulation code 
source('simulation_functions.R')
library(rmsb)
library(LaplacesDemon)
library(tidyverse)
library(parallel)
library(foreach)
library(BuyseTest)
library(ordinal)
library(doParallel)
library(doRNG)


# this is pulling the array number
# for the complete simulation, val should be all values from 1 to 36
val <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

start <- determine_start_i2(val = val)

myCluster <- makeCluster(7)
registerDoParallel(myCluster)



all_results <- foreach (i = start:(start + 1000), .combine = 'rbind',
                        .packages = c("rmsb", "LaplacesDemon", "tidyverse",
                                      "parallel", "foreach", "BuyseTest", "ordinal",
                                      "doParallel")) %dopar% {
 source('simulation_functions.R')
 set.seed(12345 + i)
 dat <- DataGen_array(n = 1000, array = val)
                                        
 b_ppo_model <- blrm(formula = y ~ trt_05, 
                                                            ppo = ~ trt_05, 
                                                            data = dat, 
                                                            method = "both", 
                                                            iter = 2000,
                                                            warmup = 1000, 
                                                            iprior = 0)
                                        
                                        results <- order_weighted_RD2(b_ppo_model, iprior = 0, dat = dat)
                                        results_CI <- order_weighted_RD2_CI(b_ppo_model, n_post_draws = 4000, dat = dat)
                                        
                                        b_po_model <- blrm(formula = y ~ trt_05,
                                                           data = dat, 
                                                           method = "both", 
                                                           iter = 2000,
                                                           warmup = 1000, 
                                                           iprior = 0)
                                        
                                        results_po <- order_weighted_RD2_PO(b_po_model, dat = dat)
                                        results_po_CI <- order_weighted_RD2_PO_CI(b_po_model, n_post_draws = 4000, dat = dat)
                                        
                                        freq_mod <- clm(as.factor(y) ~ trt_05, data = dat, link = "logit")
                                        freq_pval <- coef(summary(freq_mod))[nrow(coef(summary(freq_mod))), 4]
                                        
                                        log_reg_hypox <- glm(I(y >= 2) ~ trt_05, data = dat, family = "binomial")
                                        RD_hypox <- c(coef(log_reg_hypox)[2], confint(log_reg_hypox)[2,])
                                        
                                        log_reg_ED <- glm(I(y >= 3) ~ trt_05, data = dat, family = "binomial")
                                        RD_ED <- c(coef(log_reg_ED)[2], confint(log_reg_ED)[2,])
                                        
                                        # need to make a catch if no patients die 
                                        if (sum(dat$y >= 4) == 0) {
                                          # no patients die
                                          RD_hosp <- c(NA, NA, NA)
                                        } else if  (sum(dat$y >= 4 & dat$trt == 1) == 0 | sum(dat$y >=4 & dat$trt == 0) == 0)  {
                                          # no patients die in one group
                                          RD_hosp <- c(NA, NA, NA)
                                        } else {
                                          log_reg_hosp <- glm(I(y >= 4) ~ trt_05, data = dat, family = "binomial")
                                          RD_hosp <- c(coef(log_reg_hosp)[2], confint(log_reg_hosp)[2,])
                                        }
                                        
                                        
                                        # need to make a catch if no patients die 
                                        if (sum(dat$y == 5) == 0) {
                                          # no patients die
                                          RD_death <- c(NA, NA, NA)
                                        } else if  (sum(dat$y == 5 & dat$trt == 1) == 0 | sum(dat$y == 5 & dat$trt == 0) == 0)  {
                                          # no patients die in one group
                                          RD_death <- c(NA, NA, NA)
                                        } else {
                                          log_reg_death <- glm(I(y >= 5) ~ trt_05, data = dat, family = "binomial")
                                          RD_death <- c(coef(log_reg_death)[2], confint(log_reg_death)[2,])
                                        }
                                        
                                        BuyseTest_summary <- summary(BuyseTest(trt_05 ~ c(y, operator = "<0"), data = dat), statistic = "netBenefit")
                                        Frequentist_NB <- c(BuyseTest_summary$Delta, BuyseTest_summary$`CI [2.5% ; 97.5%]`)
                                        names(Frequentist_NB) <- c("Frequentist_NB", "Frequentist_NB_CI")
                                        names(freq_pval) <- "freq_pval"
                                        names(RD_hypox) <- c("RD_hypox", "RD_hypox_lower", "RD_hypox_upper")
                                        names(RD_ED) <- c("RD_ED", "RD_ED_lower", "RD_ED_upper")
                                        names(RD_hosp) <- c("RD_hosp", "RD_hosp_lower", "RD_hosp_upper")
                                        names(RD_death) <- c("RD_death", "RD_death_lower", "RD_death_upper")
                                        
                                        c(i, results, results_CI, results_po, results_po_CI, 
                                          Frequentist_NB, freq_pval, RD_hypox, RD_ED, RD_hosp, RD_death)
                                        
                                      }
                                     
stopImplicitCluster()
write_rds(all_results, paste0('results_array_',  val, '.rds'))


