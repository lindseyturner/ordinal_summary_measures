# example of covid out analysis 
# load all necessary packages 
library(tidyverse)
library(gtsummary)
library(rmsb)
library(MASS)
library(mice)
library(ggplot2)
library(kableExtra)
library(spatstat.utils)
library(brms)
library(WinRatio)
library(BuyseTest)
library(drord)
library(parallel)
library(foreach)
library(RColorBrewer)
library(scales)
library(cowplot)
set.seed(12345)

# set working directory 
setwd("/Users/lindseyturner/Library/CloudStorage/Box-Box/COVID_OUT_Data/")
# load in the file with the functions
source('covid_out_functions_github_050525.R')

# read in data 
covid_out <- read.csv("covid_out_data_2023-06-27.csv") 

# data cleaning 
mitt_covid <- covid_out %>% filter(in_MITT_analysis == "Yes")

mitt_covid$composite_primary_14 <- case_when(is.na(mitt_covid$screen_covid_vaccine) ~ NA, 
                                             mitt_covid$outcome_1_by_d14 == "Yes" | 
                                               mitt_covid$outcome_3_by_d14 == "Yes" |
                                               mitt_covid$outcome_4_by_d14 == "Yes" |
                                               mitt_covid$outcome_7_by_d14 == "Yes" 
                                             ~ "Yes", 
                                             mitt_covid$outcome_1_by_d14 == "No" & 
                                               mitt_covid$outcome_3_by_d14 == "No" &
                                               mitt_covid$outcome_4_by_d14 == "No" &
                                               mitt_covid$outcome_7_by_d14 == "No"  ~ "No",
                                             TRUE ~ NA)


mitt_covid$ED_hosp_dead_14 <- case_when(is.na(mitt_covid$screen_covid_vaccine) ~ NA, 
                                        mitt_covid$outcome_3_by_d14 == "Yes" |
                                          mitt_covid$outcome_4_by_d14 == "Yes" |
                                          mitt_covid$outcome_7_by_d14 == "Yes" 
                                        ~ "Yes", 
                                        mitt_covid$outcome_3_by_d14 == "No" &
                                          mitt_covid$outcome_4_by_d14 == "No" &
                                          mitt_covid$outcome_7_by_d14 == "No"  ~ "No",
                                        TRUE ~ NA)



mitt_covid$hosp_dead_14 <- case_when(is.na(mitt_covid$screen_covid_vaccine) ~ NA, 
                                     mitt_covid$outcome_4_by_d14 == "Yes" |
                                       mitt_covid$outcome_7_by_d14 == "Yes" 
                                     ~ "Yes", 
                                     mitt_covid$outcome_4_by_d14 == "No" &
                                       mitt_covid$outcome_7_by_d14 == "No"  ~ "No",
                                     TRUE ~ NA)

mitt_covid$hypoxemia_missingvax <- case_when(is.na(mitt_covid$screen_covid_vaccine) ~ NA, 
                                             TRUE ~ mitt_covid$outcome_1_by_d14)

mitt_covid$death_missingvax <- case_when(is.na(mitt_covid$screen_covid_vaccine) ~ NA, 
                                         TRUE ~ mitt_covid$outcome_7_by_d14)

mitt_covid$ordinal_outcome <- case_when(is.na(mitt_covid$outcome_1_by_d14) &
                                          is.na(mitt_covid$outcome_3_by_d14) &
                                          is.na(mitt_covid$outcome_4_by_d14) &
                                          is.na(mitt_covid$outcome_7_by_d14) ~ NA,
                                        mitt_covid$outcome_7_by_d14 == "Yes" ~ 4, 
                                        mitt_covid$outcome_4_by_d14 == "Yes" &
                                          is.na(mitt_covid$outcome_7_by_d14) == F ~ 4,
                                        mitt_covid$outcome_3_by_d14 == "Yes" &
                                          is.na(mitt_covid$outcome_4_by_d14) == F &
                                          is.na(mitt_covid$outcome_7_by_d14) == F ~ 3,
                                        mitt_covid$outcome_1_by_d14 == "Yes" &
                                          is.na(mitt_covid$outcome_3_by_d14) == F &
                                          is.na(mitt_covid$outcome_4_by_d14) == F &
                                          is.na(mitt_covid$outcome_7_by_d14) == F ~ 2,
                                        mitt_covid$outcome_1_by_d14 == "No" & 
                                          mitt_covid$outcome_3_by_d14 == "No" & 
                                          mitt_covid$outcome_4_by_d14 == "No" & 
                                          mitt_covid$outcome_7_by_d14 == "No" ~ 1,
                                        TRUE ~ NA)

covid_dropped <- mitt_covid %>% drop_na(ordinal_outcome)
covid_dropped$metformin_5 <- case_when(covid_dropped$metformin == "Yes" ~ 0.5,
                                       covid_dropped$metformin == "No" ~ -0.5)

mitt_covid$ordinal_outcome_28 <- case_when(is.na(mitt_covid$outcome_1_by_d28) &
                                             is.na(mitt_covid$outcome_3_by_d28) &
                                             is.na(mitt_covid$outcome_4_by_d28) &
                                             is.na(mitt_covid$outcome_7_by_d28) ~ NA,
                                           mitt_covid$outcome_7_by_d28 == "Yes" ~ 4, 
                                           mitt_covid$outcome_4_by_d28 == "Yes" &
                                             is.na(mitt_covid$outcome_7_by_d28) == F ~ 4,
                                           mitt_covid$outcome_3_by_d28 == "Yes" &
                                             is.na(mitt_covid$outcome_4_by_d28) == F &
                                             is.na(mitt_covid$outcome_7_by_d28) == F ~ 3,
                                           mitt_covid$outcome_1_by_d28 == "Yes" &
                                             is.na(mitt_covid$outcome_3_by_d28) == F &
                                             is.na(mitt_covid$outcome_4_by_d28) == F &
                                             is.na(mitt_covid$outcome_7_by_d28) == F ~ 2,
                                           mitt_covid$outcome_1_by_d28 == "No" & 
                                             mitt_covid$outcome_3_by_d28 == "No" & 
                                             mitt_covid$outcome_4_by_d28 == "No" & 
                                             mitt_covid$outcome_7_by_d28 == "No" ~ 1,
                                           TRUE ~ NA)


covid_dropped2 <- mitt_covid %>% drop_na(ordinal_outcome_28)
covid_dropped2$metformin_5 <- case_when(covid_dropped2$metformin == "Yes" ~ 0.5,
                                        covid_dropped2$metformin == "No" ~ -0.5)

##############################################################################################
# run unadjusted analysis 
b_ppo_model <- blrm(formula = ordinal_outcome ~ metformin_5, 
                    ppo = ~ metformin_5, 
                    data = covid_dropped, 
                    method = "both", 
                    iter = 2000,
                    warmup = 1000, 
                    iprior = 0)

#calculate order invariant weighted risk difference
order_weighted_RD_PPO <- order_weighted_RD2(b = b_ppo_model, iprior = 0, 
                                            dat = covid_dropped, trt = "metformin_5",
                                            outcome = "ordinal_outcome")
# get CI from  as.matrix(b_ppo_model$rstan) and do for each iteration and take quantiles
# check each iteration doesn't get negative value - look at probability of death in each arm 
order_weighted_RD_PPO_CI <- order_weighted_RD2_CI(b = b_ppo_model, dat = covid_dropped,
                                                  trt = "metformin_5", outcome = "ordinal_outcome")

# day 28
# fit bayesian partial proportional odds model 
b_ppo_model_28 <- blrm(formula = ordinal_outcome_28 ~ metformin_5, 
                       ppo = ~ metformin_5, 
                       data = covid_dropped2, 
                       method = "both", 
                       iter = 2000,
                       warmup = 1000, 
                       iprior = 0)

#calculate order invariant weighted risk difference
order_weighted_RD_PPO_28 <- order_weighted_RD2(b_ppo_model_28, iprior = 0, dat = covid_dropped2,
                                               trt = "metformin_5", outcome = "ordinal_outcome_28")
order_weighted_RD_PPO_28_CI <- order_weighted_RD2_CI(b = b_ppo_model_28, dat = covid_dropped2,
                                                     trt = "metformin_5", outcome = "ordinal_outcome_28")

##############################################################################################
# run adjusted analysis for fluvoxamine, ivermectin and vaccination status
### adjusted PPO model 
covid_dropped_dropvax <- covid_dropped %>% drop_na(screen_covid_vaccine)
b_ppo_model_adj <- blrm(formula = ordinal_outcome ~ metformin_5 + fluvoxamine + ivermectin + screen_covid_vaccine,  
                        ppo = ~ metformin_5, 
                        data = covid_dropped_dropvax, 
                        method = "both", 
                        iter = 2000,
                        warmup = 1000, 
                        iprior = 0)

b_ppo_model_adj_results_bb <- order_weighted_RD2_PPO_adj_bb(b_ppo_model_adj,
                                                            dat = covid_dropped_dropvax, iprior = 0,
                                                            outcome = "ordinal_outcome", trt = "metformin_5")


b_ppo_model_adj_results_CI_bb <- order_weighted_RD2_CI_adj_bb(b_ppo_model_adj, 
                                                              dat = covid_dropped_dropvax, iprior = 0,
                                                              outcome = "ordinal_outcome", trt = "metformin_5")

### adjusted PPO model 
covid_dropped_dropvax_28 <- covid_dropped2 %>% drop_na(screen_covid_vaccine)
b_ppo_model_adj_28 <- blrm(formula = ordinal_outcome_28 ~ metformin_5 + fluvoxamine + 
                             ivermectin + screen_covid_vaccine,  
                           ppo = ~ metformin_5, 
                           data = covid_dropped_dropvax_28, 
                           method = "both", 
                           iter = 2000,
                           warmup = 1000, 
                           iprior = 0)


b_ppo_model_adj_results_28_bb <- order_weighted_RD2_PPO_adj_bb(b_ppo_model_adj_28, 
                                                               dat = covid_dropped_dropvax_28, iprior = 0,
                                                               outcome = "ordinal_outcome_28", trt = "metformin_5")

b_ppo_model_adj_results_CI_28_bb <- order_weighted_RD2_CI_adj_bb(b_ppo_model_adj_28, 
                                                                 dat = covid_dropped_dropvax_28, iprior = 0,
                                                                 outcome = "ordinal_outcome_28", trt = "metformin_5")

##############################################################################################
# run adjusted analysis with ivermectin, fluvoxamine and vaccination status
# ppo model with hypoxemia and nothing bad combined 
covid_dropped_dropvax$ordinal_outcome_nohypo <- case_when(covid_dropped_dropvax$ordinal_outcome <= 2 ~ 1,
                                                          covid_dropped_dropvax$ordinal_outcome == 3 ~ 2,
                                                          covid_dropped_dropvax$ordinal_outcome == 4 ~ 3, 
                                                          TRUE ~ NA)

b_ppo_model_adj_nohyp <- blrm(formula = ordinal_outcome_nohypo ~ metformin_5 + fluvoxamine + 
                                ivermectin + screen_covid_vaccine,
                              ppo = ~ metformin_5, 
                              data = covid_dropped_dropvax, 
                              method = "both", 
                              iter = 2000,
                              warmup = 1000, 
                              iprior = 0)


b_ppo_model_adj_results_nohyp_bb <- order_weighted_RD2_PPO_adj_bb(b_ppo_model_adj_nohyp, 
                                                                  dat = covid_dropped_dropvax, iprior = 0,
                                                                  outcome = "ordinal_outcome_nohypo", trt = "metformin_5")
b_ppo_model_adj_results_CI_nohyp_bb <- order_weighted_RD2_CI_adj_bb(b_ppo_model_adj_nohyp, 
                                                                    dat = covid_dropped_dropvax, iprior = 0,
                                                                    outcome = "ordinal_outcome_nohypo", trt = "metformin_5")

covid_dropped_dropvax_28$ordinal_outcome_nohypo <- case_when(covid_dropped_dropvax_28$ordinal_outcome <= 2 ~ 1,
                                                             covid_dropped_dropvax_28$ordinal_outcome == 3 ~ 2,
                                                             covid_dropped_dropvax_28$ordinal_outcome == 4 ~ 3)

b_ppo_model_adj28_nohyp <- blrm(formula = ordinal_outcome_nohypo ~ metformin_5 + fluvoxamine + 
                                  ivermectin + screen_covid_vaccine,
                                ppo = ~ metformin_5, 
                                data = covid_dropped_dropvax_28, 
                                method = "both", 
                                iter = 2000,
                                warmup = 1000, 
                                iprior = 0)

b_ppo_model_adj_results_28_nohyp_bb <- order_weighted_RD2_PPO_adj_bb(b_ppo_model_adj28_nohyp, 
                                                                     dat = covid_dropped_dropvax_28, iprior = 0,
                                                                     outcome = "ordinal_outcome_nohypo", trt = "metformin_5")
b_ppo_model_adj_results_28_CI_nohyp_bb <- order_weighted_RD2_CI_adj_bb(b_ppo_model_adj28_nohyp, 
                                                                       dat = covid_dropped_dropvax_28, iprior = 0,
                                                                       outcome = "ordinal_outcome_nohypo", trt = "metformin_5")

##############################################################################################
# run adjusted analysis with ivermectin, fluvoxaminem, vaccination status, sex, diabetes and BMI
### adjusted PPO model adjusting for additional vars
b_ppo_model_adj_add <- blrm(formula = ordinal_outcome ~ metformin_5 + fluvoxamine + ivermectin + screen_covid_vaccine +
                              sex + pd_comorbidities___1 + bmi,
                            ppo = ~ metformin_5, 
                            data = covid_dropped_dropvax, 
                            method = "both", 
                            iter = 2000,
                            warmup = 1000, 
                            iprior = 0)


b_ppo_model_adj_results_additionalvars_bb <- order_weighted_RD2_PPO_adj_bb(b_ppo_model_adj_add, 
                                                                           dat = covid_dropped_dropvax, iprior = 0,
                                                                           outcome = "ordinal_outcome", trt = "metformin_5")
b_ppo_model_adj_results_CI_additionalvars_bb <- order_weighted_RD2_CI_adj_bb(b_ppo_model_adj_add, 
                                                                             dat = covid_dropped_dropvax, iprior = 0,
                                                                             outcome = "ordinal_outcome", trt = "metformin_5")

#adjusted PPO model at day 28 with additional vars
b_ppo_model_adj28_addition <- blrm(formula = ordinal_outcome_28 ~ metformin_5 + fluvoxamine + 
                                     ivermectin + screen_covid_vaccine + sex + pd_comorbidities___1 + bmi,
                                   ppo = ~ metformin_5, 
                                   data = covid_dropped_dropvax_28, 
                                   method = "both", 
                                   iter = 2000,
                                   warmup = 1000, 
                                   iprior = 0)


b_ppo_model_adj_results_28_additionalvars_bb <- order_weighted_RD2_PPO_adj_bb(b_ppo_model_adj28_addition,
                                                                              dat = covid_dropped_dropvax_28, iprior = 0,
                                                                              outcome = "ordinal_outcome_28", trt = "metformin_5")
b_ppo_model_adj_results_28_CI_additionalvars_bb <- order_weighted_RD2_CI_adj_bb(b_ppo_model_adj28_addition, 
                                                                                dat = covid_dropped_dropvax_28, iprior = 0,
                                                                                outcome = "ordinal_outcome_28", trt = "metformin_5")

##############################################################################################
# run unadjusted analysis 
# fit bayesian proportional odds model for day 14
b_po_model <- blrm(formula = ordinal_outcome ~ metformin_5,  
                   data = covid_dropped, 
                   method = "both", 
                   iter = 2000,
                   warmup = 1000, 
                   iprior = 0)

#calculate order invariant weighted risk difference
order_weighted_RD_PO <- order_weighted_RD2_PO(b_po_model, dat = covid_dropped,
                                              trt = "metformin_5", outcome = "ordinal_outcome")
order_weighted_RD_PO_CI <- order_weighted_RD2_PO_CI(b_po_model, dat = covid_dropped,
                                                    trt = "metformin_5", outcome = "ordinal_outcome")

# fit bayesian proportional odds model for day 28
b_po_model_28 <- blrm(formula = ordinal_outcome_28 ~ metformin_5,  
                      data = covid_dropped2, 
                      method = "both", 
                      iter = 2000,
                      warmup = 1000,
                      iprior = 0)

#calculate order invariant weighted risk difference
order_weighted_RD_PO_28 <- order_weighted_RD2_PO(b_po_model_28, dat = covid_dropped2,
                                                 trt = "metformin_5", outcome = "ordinal_outcome_28")
order_weighted_RD_PO_28_CI <- order_weighted_RD2_PO_CI(b_po_model_28, dat = covid_dropped2,
                                                       trt = "metformin_5", outcome = "ordinal_outcome_28")

##############################################################################################
# run PO adjusted analysis with fluvoxamine, ivermectin and vaccination status
# fit bayesian proportional odds model for day 14 adjusted
# drop the missing for covid vaccine 
covid_dropped_dropvax <- covid_dropped %>% drop_na(screen_covid_vaccine)
b_po_model_adj <- blrm(formula = ordinal_outcome ~ metformin_5 + fluvoxamine + ivermectin + screen_covid_vaccine,  
                       data = covid_dropped_dropvax, 
                       method = "both", 
                       iter = 2000,
                       warmup = 1000, 
                       iprior = 0)

order_weighted_RD_PO_adj_bb <- order_weighted_RD2_PO_adj_bb(b = b_po_model_adj, c = 4,
                                                            dat = covid_dropped_dropvax, 
                                                            trt_variable = "metformin_5",
                                                            outcome = "ordinal_outcome")
order_weighted_RD_PO_adj_CI_bb <- order_weighted_RD2_CI_PO_adj_bb(b = b_po_model_adj,
                                                                  c = 4, n_post_draws = 4000,
                                                                  iprior = 0,
                                                                  dat = covid_dropped_dropvax,
                                                                  outcome = "ordinal_outcome",
                                                                  trt = "metformin_5")

covid_dropped_dropvax_28 <- covid_dropped2 %>% drop_na(screen_covid_vaccine)
b_po_model_adj_28 <- blrm(formula = ordinal_outcome_28 ~ metformin_5 + fluvoxamine + 
                            ivermectin + screen_covid_vaccine,  
                          data = covid_dropped_dropvax_28, 
                          method = "both", 
                          iter = 2000,
                          warmup = 1000, 
                          iprior = 0)

order_weighted_RD_PO_28_adj_bb <- order_weighted_RD2_PO_adj_bb(b = b_po_model_adj_28, 
                                                               dat = covid_dropped_dropvax_28, c = 4,
                                                               trt_variable = "metformin_5",
                                                               outcome = "ordinal_outcome_28")

order_weighted_RD_PO_28_adj_CI_bb <- order_weighted_RD2_CI_PO_adj_bb(b = b_po_model_adj_28,
                                                                     c = 4, n_post_draws = 4000,
                                                                     iprior = 0,
                                                                     dat = covid_dropped_dropvax_28,
                                                                     outcome = "ordinal_outcome_28",
                                                                     trt = "metformin_5")

##############################################################################################
# run PO adjusted analysis with fluvoxamine, ivermectin and vaccination status
#ppo model with hypoxemia and nothing bad combined 
covid_dropped_dropvax$ordinal_outcome_nohypo <- case_when(covid_dropped_dropvax$ordinal_outcome <= 2 ~ 1,
                                                          covid_dropped_dropvax$ordinal_outcome == 3 ~ 2,
                                                          covid_dropped_dropvax$ordinal_outcome == 4 ~ 3, 
                                                          TRUE ~ NA)

b_po_model_adj_nohyp <- blrm(formula = ordinal_outcome_nohypo ~ metformin_5 + fluvoxamine + 
                               ivermectin + screen_covid_vaccine,
                             data = covid_dropped_dropvax, 
                             method = "both", 
                             iter = 2000,
                             warmup = 1000, 
                             iprior = 0)


b_po_model_adj_results_nohyp_bb <- order_weighted_RD2_PO_adj_bb(b_po_model_adj_nohyp, c= 3,
                                                                dat = covid_dropped_dropvax,
                                                                outcome = "ordinal_outcome_nohypo", trt = "metformin_5")
b_po_model_adj_results_CI_nohyp_bb <- order_weighted_RD2_CI_PO_adj_bb(b_po_model_adj_nohyp, c = 3,
                                                                      dat = covid_dropped_dropvax, 
                                                                      outcome = "ordinal_outcome_nohypo", trt = "metformin_5")

covid_dropped_dropvax_28$ordinal_outcome_nohypo <- case_when(covid_dropped_dropvax_28$ordinal_outcome <= 2 ~ 1,
                                                             covid_dropped_dropvax_28$ordinal_outcome == 3 ~ 2,
                                                             covid_dropped_dropvax_28$ordinal_outcome == 4 ~ 3)

b_po_model_adj28_nohypox <- blrm(formula = ordinal_outcome_nohypo ~ metformin_5 + fluvoxamine + 
                                   ivermectin + screen_covid_vaccine,
                                 data = covid_dropped_dropvax_28, 
                                 method = "both", 
                                 iter = 2000,
                                 warmup = 1000, 
                                 iprior = 0)

b_po_model_adj_results_28_nohyp_bb <- order_weighted_RD2_PO_adj_bb(b_po_model_adj28_nohypox,
                                                                   dat = covid_dropped_dropvax_28,c = 3,
                                                                   outcome = "ordinal_outcome_nohypo", trt = "metformin_5")
b_po_model_adj_results_CI_28_nohyp_bb <- order_weighted_RD2_CI_PO_adj_bb(b_po_model_adj28_nohypox, c = 3, 
                                                                         dat = covid_dropped_dropvax_28,
                                                                         outcome = "ordinal_outcome_nohypo", trt = "metformin_5")

##############################################################################################
# run PO adjusted analysis with fluvoxamine, ivermectin, vaccination status, sex, diabetes, and BMI
### adjusted PO model adjusting for additional vars
b_po_model_adj_add <- blrm(formula = ordinal_outcome ~ metformin_5 + fluvoxamine + ivermectin + screen_covid_vaccine +
                             sex + pd_comorbidities___1 + bmi,
                           data = covid_dropped_dropvax, 
                           method = "both", 
                           iter = 2000,
                           warmup = 1000, 
                           iprior = 0)

b_po_model_adj_results_additionalvars_bb <- order_weighted_RD2_PO_adj_bb(b_po_model_adj_add, c = 4,
                                                                         dat = covid_dropped_dropvax,
                                                                         outcome = "ordinal_outcome", trt = "metformin_5")
b_po_model_adj_results_CI_additionalvars_bb <- order_weighted_RD2_CI_PO_adj_bb(b_po_model_adj_add, c = 4, 
                                                                               dat = covid_dropped_dropvax,
                                                                               outcome = "ordinal_outcome", trt = "metformin_5")

#adjusted PO model at day 28 with additional vars
b_po_model_adj28_addition <- blrm(formula = ordinal_outcome_28 ~ metformin_5 + fluvoxamine + 
                                    ivermectin + screen_covid_vaccine + sex + pd_comorbidities___1 + bmi,
                                  data = covid_dropped_dropvax_28, 
                                  method = "both", 
                                  iter = 2000,
                                  warmup = 1000, 
                                  iprior = 0)

b_po_model_adj_results_28_additionalvars_bb <- order_weighted_RD2_PO_adj_bb(b_po_model_adj28_addition, 
                                                                            dat = covid_dropped_dropvax_28, c = 4,
                                                                            outcome = "ordinal_outcome_28", trt = "metformin_5")
b_po_model_adj_results_28_CI_additionalvars_bb <- order_weighted_RD2_CI_PO_adj_bb(b_po_model_adj28_addition, 
                                                                                  c = 4,
                                                                                  dat = covid_dropped_dropvax_28, 
                                                                                  outcome = "ordinal_outcome_28", trt = "metformin_5")


####################################################################################################
### Figures 

# set up data for plot
RD <- matrix(NA, nrow = 14, ncol = 6)
RD[,1] <- c(b_ppo_model_adj_results_bb$RD_leq, b_ppo_model_adj_results_bb$unweighted_RD,
            b_ppo_model_adj_results_bb$weighted_RD_B, b_ppo_model_adj_results_bb$weighted_RD_cum,
            b_ppo_model_adj_results_bb$weighted_RD_both,
            order_weighted_RD_PO_adj_bb$RD_leq, order_weighted_RD_PO_adj_bb$unweighted_RD,
            order_weighted_RD_PO_adj_bb$weighted_RD_b, order_weighted_RD_PO_adj_bb$weighted_RD_cum,
            order_weighted_RD_PO_adj_bb$weighted_RD_both)
RD[,2] <- c(rep("Bayesian PPO Model", 7),
            rep("Bayesian PO Model", 7))
RD[,3] <- c("Nothing Bad", "Hypoxemia or Better", "ED Visit or Better", "ARD", "wRD", "wRD_cumulative", "wRD_both",
            "Nothing Bad", "Hypoxemia or Better", "ED Visit or Better", "ARD", "wRD", "wRD_cumulative", "wRD_both")
RD[,4] <- c(b_ppo_model_adj_results_CI_bb$RD_leq_95CI[1,1:3],  b_ppo_model_adj_results_CI_bb$unweighted_RD_95CI[1],
            b_ppo_model_adj_results_CI_bb$wRD_bcontrol_95CI[1], b_ppo_model_adj_results_CI_bb$wRD_cum_95CI[1],
            b_ppo_model_adj_results_CI_bb$wRD_both_95CI[1],
            order_weighted_RD_PO_adj_CI_bb$RD_leq_95CI[1,1:3],
            order_weighted_RD_PO_adj_CI_bb$unweighted_RD_95CI[1],
            order_weighted_RD_PO_adj_CI_bb$wRD_bcontrol_95CI[1],
            order_weighted_RD_PO_adj_CI_bb$wRD_cum_95CI[1],
            order_weighted_RD_PO_adj_CI_bb$wRD_both_95CI[1])
RD[,5] <- c(b_ppo_model_adj_results_CI_bb$RD_leq_95CI[2,1:3], b_ppo_model_adj_results_CI_bb$unweighted_RD_95CI[2],
            b_ppo_model_adj_results_CI_bb$wRD_bcontrol_95CI[2],
            b_ppo_model_adj_results_CI_bb$wRD_cum_95CI[2], b_ppo_model_adj_results_CI_bb$wRD_both_95CI[2],
            order_weighted_RD_PO_adj_CI_bb$RD_leq_95CI[2,1:3], 
            order_weighted_RD_PO_adj_CI_bb$unweighted_RD_95CI[2], 
            order_weighted_RD_PO_adj_CI_bb$wRD_bcontrol_95CI[2],
            order_weighted_RD_PO_adj_CI_bb$wRD_cum_95CI[2],
            order_weighted_RD_PO_adj_CI_bb$wRD_both_95CI[2])
RD[,6] <- c(b_ppo_model_adj_results_bb$weights_split, mean(b_ppo_model_adj_results_bb$weights_split),
            mean(b_ppo_model_adj_results_bb$weights_split), mean(b_ppo_model_adj_results_bb$weights_split),
            mean(b_ppo_model_adj_results_bb$weights_split),
            order_weighted_RD_PO_adj_bb$weights_split, mean(order_weighted_RD_PO_adj_bb$weights_split),
            mean(order_weighted_RD_PO_adj_bb$weights_split), mean(order_weighted_RD_PO_adj_bb$weights_split),
            mean(order_weighted_RD_PO_adj_bb$weights_split)) * 2
RD <- as.data.frame(RD)
colnames(RD) <- c("RD", "Type", "Time", "Lower", "Upper", "Weights")
RD$Type <- factor(RD$Type, levels = c("Bayesian PPO Model", "Bayesian PO Model"))
RD$Time <- factor(RD$Time, levels = c("Nothing Bad", "Hypoxemia or Better", "ED Visit or Better", "ARD", "wRD", "wRD_cumulative", "wRD_both"))



RD_adj_14_pres3 <- RD %>% 
  filter(Time %in% c("Nothing Bad", "Hypoxemia or Better", "ED Visit or Better", "ARD", "wRD_both")) %>%
  mutate(Time2 = factor(Time, labels = c("Nothing Bad", "Hypoxemia or Better", "ED Visit or Better", "ARD", "wRD"))) %>%
  ggplot(aes(x=Time2, y=as.numeric(RD), colour=Type)) + 
  geom_errorbar(aes(ymin=as.numeric(Lower), ymax=as.numeric(Upper)), width=.2, 
                position=position_dodge(width=0.5)) +
  geom_point(position=position_dodge(width=0.5)) + 
  geom_hline(yintercept = 0, linetype = "dashed") + 
  geom_vline(xintercept = 3.5) +
  scale_color_manual(values = c(brewer.pal(8, "Paired")[2], brewer.pal(8, "Paired")[4])) + 
  #scale_x_discrete(labels = label_wrap(10)) +
  ylab("Risk Difference") + 
  xlab("Event") + theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_text(size = 9),
        axis.title = element_text(size = 14), legend.text = element_text(size = 9),
        legend.title = element_text(size = 14)) + 
  scale_x_discrete(labels=c("Nothing Bad" = "1|2", 
                            "Hypoxemia or Better" = "2|3", 
                            "ED Visit or Better" = "3|4", 
                            "ARD" = "ARD", 
                            "wRD" = "wRD"))



# set up data for plot
RD <- matrix(NA, nrow = 10, ncol = 5)
RD[,1] <- c(b_ppo_model_adj_results_28_bb$RD_leq, b_ppo_model_adj_results_28_bb$unweighted_RD,
            b_ppo_model_adj_results_28_bb$weighted_RD_both,
            order_weighted_RD_PO_28_adj_bb$RD_leq, order_weighted_RD_PO_28_adj_bb$unweighted_RD,
            order_weighted_RD_PO_28_adj_bb$weighted_RD_both)
RD[,2] <- c(rep("Bayesian PPO Model", 5),
            rep("Bayesian PO Model", 5))
RD[,3] <- c("Nothing Bad", "Hypoxemia or Better", "ED Visit or Better", "ARD", "wRD_both",
            "Nothing Bad", "Hypoxemia or Better", "ED Visit or Better", "ARD", "wRD_both")
RD[,4] <- c(b_ppo_model_adj_results_CI_28_bb$RD_leq_95CI[1,1:3], b_ppo_model_adj_results_CI_28_bb$unweighted_RD_95CI[1],
            b_ppo_model_adj_results_CI_28_bb$wRD_both_95CI[1],
            order_weighted_RD_PO_28_adj_CI_bb$RD_leq_95CI[1,1:3],
            order_weighted_RD_PO_28_adj_CI_bb$unweighted_RD_95CI[1],
            order_weighted_RD_PO_28_adj_CI_bb$wRD_both_95CI[1])
RD[,5] <- c(b_ppo_model_adj_results_CI_28_bb$RD_leq_95CI[2,1:3], b_ppo_model_adj_results_CI_28_bb$unweighted_RD_95CI[2],
            b_ppo_model_adj_results_CI_28_bb$wRD_both_95CI[2],
            order_weighted_RD_PO_28_adj_CI_bb$RD_leq_95CI[2,1:3], 
            order_weighted_RD_PO_28_adj_CI_bb$unweighted_RD_95CI[2], 
            order_weighted_RD_PO_28_adj_CI_bb$wRD_both_95CI[2])
RD <- as.data.frame(RD)
colnames(RD) <- c("RD", "Type", "Time", "Lower", "Upper")
RD$Type <- factor(RD$Type, levels = c("Bayesian PPO Model", "Bayesian PO Model"))
RD$Time <- factor(RD$Time, levels = c("Nothing Bad", "Hypoxemia or Better", "ED Visit or Better", "ARD", "wRD_both"))

RD_adj_28_pres3 <- RD %>% 
  filter(Time %in% c("Nothing Bad", "Hypoxemia or Better", "ED Visit or Better", "ARD", "wRD_both")) %>%
  mutate(Time2 = factor(Time, labels = c("Nothing Bad", "Hypoxemia or Better", "ED Visit or Better", "ARD", "wRD"))) %>%
  ggplot(aes(x=Time2, y=as.numeric(RD), colour=Type)) + 
  geom_errorbar(aes(ymin=as.numeric(Lower), ymax=as.numeric(Upper)), width=.2, 
                position=position_dodge(width=0.5)) +
  geom_point(position=position_dodge(width=0.5)) + 
  geom_hline(yintercept = 0, linetype = "dashed") + 
  geom_vline(xintercept = 3.5) +
  scale_color_manual(values = c(brewer.pal(8, "Paired")[2], brewer.pal(8, "Paired")[4])) + 
  #scale_x_discrete(labels = label_wrap(10)) +
  ylab("Risk Difference") + 
  xlab("Event") + theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_text(size = 9),
        axis.title = element_text(size = 14), legend.text = element_text(size = 9),
        legend.title = element_text(size = 14))  + 
  scale_x_discrete(labels=c("Nothing Bad" = "1|2", 
                            "Hypoxemia or Better" = "2|3", 
                            "ED Visit or Better" = "3|4", 
                            "ARD" = "ARD", 
                            "wRD" = "wRD"))

# set up OR for plot
OR <- matrix(NA, nrow = 14, ncol = 6)
OR[,1] <- c(b_ppo_model_adj_results_bb$log_OR, b_ppo_model_adj_results_bb$unweighted_log_OR,
            b_ppo_model_adj_results_bb$log_weighted_OR_B, b_ppo_model_adj_results_bb$log_weighted_OR_cum,
            b_ppo_model_adj_results_bb$log_weighted_OR_both,
            order_weighted_RD_PO_adj_bb$log_OR, order_weighted_RD_PO_adj_bb$unweighted_log_OR,
            order_weighted_RD_PO_adj_bb$log_weighted_OR_b, order_weighted_RD_PO_adj_bb$log_weighted_OR_cum,
            order_weighted_RD_PO_adj_bb$log_weighted_OR_both)
OR[,2] <- c(rep("Bayesian PPO Model", 7), rep("Bayesian PO Model", 7))
OR[,3] <- c("Nothing Bad", "Hypoxemia or Better", "ED Visit or Better", "AOR", "wOR", "wOR_cumulative", "wOR_both",
            "Nothing Bad", "Hypoxemia or Better", "ED Visit or Better", "AOR", "wOR", "wOR_cumulative", "wOR_both")
OR[,4] <- c(b_ppo_model_adj_results_CI_bb$log_OR_95CI[1,1:3], b_ppo_model_adj_results_CI_bb$unweighted_log_OR_95CI[1],
            b_ppo_model_adj_results_CI_bb$log_wOR_bcontrol_95CI[1],
            b_ppo_model_adj_results_CI_bb$log_wOR_cum_95CI[1],
            b_ppo_model_adj_results_CI_bb$log_wOR_both_95CI[1],
            order_weighted_RD_PO_adj_CI_bb$log_OR_95CI[1,1:3],
            order_weighted_RD_PO_adj_CI_bb$unweighted_log_OR_95CI[1],
            order_weighted_RD_PO_adj_CI_bb$log_wOR_bcontrol_95CI[1],
            order_weighted_RD_PO_adj_CI_bb$log_wOR_cum_95CI[1],
            order_weighted_RD_PO_adj_CI_bb$log_wOR_both_95CI[1])
OR[,5] <- c(b_ppo_model_adj_results_CI_bb$log_OR_95CI[2,1:3],  b_ppo_model_adj_results_CI_bb$unweighted_log_OR_95CI[2],
            b_ppo_model_adj_results_CI_bb$log_wOR_bcontrol_95CI[2],
            b_ppo_model_adj_results_CI_bb$log_wOR_cum_95CI[2],
            b_ppo_model_adj_results_CI_bb$log_wOR_both_95CI[2],
            order_weighted_RD_PO_adj_CI_bb$log_OR_95CI[2,1:3],
            order_weighted_RD_PO_adj_CI_bb$unweighted_log_OR_95CI[2],
            order_weighted_RD_PO_adj_CI_bb$log_wOR_bcontrol_95CI[2],
            order_weighted_RD_PO_adj_CI_bb$log_wOR_cum_95CI[2],
            order_weighted_RD_PO_adj_CI_bb$log_wOR_both_95CI[2])
OR[,6] <- c(b_ppo_model_adj_results_bb$weights_split, mean(b_ppo_model_adj_results_bb$weights_split),
            mean(b_ppo_model_adj_results_bb$weights_split), mean(b_ppo_model_adj_results_bb$weights_split),
            mean(b_ppo_model_adj_results_bb$weights_split),
            order_weighted_RD_PO_adj_bb$weights_split, mean(order_weighted_RD_PO_adj_bb$weights_split),
            mean(order_weighted_RD_PO_adj_bb$weights_split), mean(order_weighted_RD_PO_adj_bb$weights_split),
            mean(order_weighted_RD_PO_adj_bb$weights_split)) * 2
OR <- as.data.frame(OR)
colnames(OR) <- c("OR", "Type", "Time", "Lower", "Upper", "Weights")
OR$Type <- factor(OR$Type, levels = c( "Bayesian PPO Model", 
                                       "Bayesian PO Model"))
OR$Time <- factor(OR$Time, levels = c("Nothing Bad", "Hypoxemia or Better", "ED Visit or Better", "AOR", "wOR",
                                      "wOR_cumulative", "wOR_both"))


OR_adj_14_pres3 <- OR %>% 
  filter(Time %in% c("Nothing Bad", "Hypoxemia or Better", "ED Visit or Better", "AOR", "wOR_both")) %>%
  mutate(Time2 = factor(Time, labels = c("Nothing Bad", "Hypoxemia or Better", "ED Visit or Better", "AOR", "wOR"))) %>%
  ggplot(aes(x=Time2, y=as.numeric(OR), colour=Type)) + 
  geom_errorbar(aes(ymin=as.numeric(Lower), ymax=as.numeric(Upper)), width=.2, 
                position=position_dodge(width=0.5)) +
  geom_point(position=position_dodge(width=0.5)) + 
  geom_hline(yintercept = 0, linetype = "dashed") + 
  geom_vline(xintercept = 3.5) +
  scale_color_manual(values = c(brewer.pal(8, "Paired")[2], brewer.pal(8, "Paired")[4])) + 
  # scale_x_discrete(labels = label_wrap(10)) +
  ylab("log(Odds Ratio)") + 
  xlab("Event") + theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_text(size = 9),
        axis.title = element_text(size = 14), legend.text = element_text(size = 9),
        legend.title = element_text(size = 14))  + 
  scale_x_discrete(labels=c("Nothing Bad" = "1|2", 
                            "Hypoxemia or Better" = "2|3", 
                            "ED Visit or Better" = "3|4", 
                            "ARD" = "ARD", 
                            "wRD" = "wRD"))


# set up data for plot
OR <- matrix(NA, nrow = 10, ncol = 5)
OR[,1] <- c(b_ppo_model_adj_results_28_bb$log_OR, b_ppo_model_adj_results_28_bb$unweighted_log_OR,
            b_ppo_model_adj_results_28_bb$log_weighted_OR_both,
            order_weighted_RD_PO_28_adj_bb$log_OR, order_weighted_RD_PO_28_adj_bb$unweighted_log_OR,
            order_weighted_RD_PO_28_adj_bb$log_weighted_OR_both)
OR[,2] <- c(rep("Bayesian PPO Model", 5), rep("Bayesian PO Model", 5))
OR[,3] <- c("Nothing Bad", "Hypoxemia or Better", "ED Visit or Better", "AOR", "wOR_both",
            "Nothing Bad", "Hypoxemia or Better", "ED Visit or Better", "AOR", "wOR_both")
OR[,4] <- c(b_ppo_model_adj_results_CI_28_bb$log_OR_95CI[1,1:3], 
            b_ppo_model_adj_results_CI_28_bb$unweighted_log_OR_95CI[1],
            b_ppo_model_adj_results_CI_28_bb$log_wOR_both_95CI[1],
            order_weighted_RD_PO_28_adj_CI_bb$log_OR_95CI[1,1:3],
            order_weighted_RD_PO_28_adj_CI_bb$unweighted_log_OR_95CI[1],
            order_weighted_RD_PO_28_adj_CI_bb$log_wOR_both_95CI[1])
OR[,5] <- c(b_ppo_model_adj_results_CI_28_bb$log_OR_95CI[2,1:3], b_ppo_model_adj_results_CI_28_bb$unweighted_log_OR_95CI[2],
            b_ppo_model_adj_results_CI_28_bb$log_wOR_both_95CI[2],
            order_weighted_RD_PO_28_adj_CI_bb$log_OR_95CI[2,1:3],
            order_weighted_RD_PO_28_adj_CI_bb$unweighted_log_OR_95CI[2],
            order_weighted_RD_PO_28_adj_CI_bb$log_wOR_both_95CI[2])
OR <- as.data.frame(OR)
colnames(OR) <- c("OR", "Type", "Time", "Lower", "Upper")
OR$Type <- factor(OR$Type, levels = c( "Bayesian PPO Model", 
                                       "Bayesian PO Model"))
OR$Time <- factor(OR$Time, levels = c("Nothing Bad", "Hypoxemia or Better", "ED Visit or Better", "AOR", "wOR_both"))

OR_adj_28_pres3 <- OR %>% 
  filter(Time %in% c("Nothing Bad", "Hypoxemia or Better", "ED Visit or Better", "AOR", "wOR_both")) %>%
  mutate(Time2 = factor(Time, labels = c("Nothing Bad", "Hypoxemia or Better", "ED Visit or Better", "AOR", "wOR"))) %>%
  ggplot(aes(x=Time2, y=as.numeric(OR), colour=Type)) + 
  geom_errorbar(aes(ymin=as.numeric(Lower), ymax=as.numeric(Upper)), width=.2, 
                position=position_dodge(width=0.5)) +
  geom_point(position=position_dodge(width=0.5)) + 
  geom_hline(yintercept = 0, linetype = "dashed") + 
  geom_vline(xintercept = 3.5) +
  scale_color_manual(values = c(brewer.pal(8, "Paired")[2], brewer.pal(8, "Paired")[4])) + 
  # scale_x_discrete(labels = label_wrap(10)) +
  ylab("log(Odds Ratio)") + 
  xlab("Event") + theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_text(size = 9),
        axis.title = element_text(size = 14), legend.text = element_text(size = 9),
        legend.title = element_text(size = 14)) + 
  scale_x_discrete(labels=c("Nothing Bad" = "1|2", 
                            "Hypoxemia or Better" = "2|3", 
                            "ED Visit or Better" = "3|4", 
                            "ARD" = "ARD", 
                            "wRD" = "wRD"))


# set up data for plot
RR <- matrix(NA, nrow = 14, ncol = 6)
RR[,1] <- c(b_ppo_model_adj_results_bb$log_RR_leq,
            b_ppo_model_adj_results_bb$unweighted_log_RR,
            b_ppo_model_adj_results_bb$log_weighted_RR_B,
            b_ppo_model_adj_results_bb$log_weighted_RR_cum,
            b_ppo_model_adj_results_bb$log_weighted_RR_both,
            order_weighted_RD_PO_adj_bb$log_RR_leq, order_weighted_RD_PO_adj_bb$unweighted_log_RR,
            order_weighted_RD_PO_adj_bb$log_weighted_RR_b,
            order_weighted_RD_PO_adj_bb$log_weighted_RR_cum,
            order_weighted_RD_PO_adj_bb$log_weighted_RR_both)
RR[,2] <- c(rep("Bayesian PPO Model", 7), rep("Bayesian PO Model", 7))
RR[,3] <- c("Nothing Bad", "Hypoxemia or Better", "ED Visit or Better", "ARR", "wRR", "wRR_cumulative", "wRR_both",
            "Nothing Bad", "Hypoxemia or Better", "ED Visit or Better", "ARR", "wRR", "wRR_cumulative", "wRR_both")
RR[,4] <- c(b_ppo_model_adj_results_CI_bb$log_RR_leq_95CI[1,1:3], 
            b_ppo_model_adj_results_CI_bb$unweighted_log_RR_95CI[1],
            b_ppo_model_adj_results_CI_bb$log_wRR_bcontrol_95CI[1],
            b_ppo_model_adj_results_CI_bb$log_wRR_cum_95CI[1],
            b_ppo_model_adj_results_CI_bb$log_wRR_both_95CI[1],
            order_weighted_RD_PO_adj_CI_bb$log_RR_leq_95CI[1,1:3],
            order_weighted_RD_PO_adj_CI_bb$unweighted_log_RR_95CI[1],
            order_weighted_RD_PO_adj_CI_bb$log_wRR_bcontrol_95CI[1],
            order_weighted_RD_PO_adj_CI_bb$log_wRR_cum_95CI[1],
            order_weighted_RD_PO_adj_CI_bb$log_wRR_both_95CI[1])
RR[,5] <- c(b_ppo_model_adj_results_CI_bb$log_RR_leq_95CI[2,1:3],  
            b_ppo_model_adj_results_CI_bb$unweighted_log_RR_95CI[2],
            b_ppo_model_adj_results_CI_bb$log_wRR_bcontrol_95CI[2],
            b_ppo_model_adj_results_CI_bb$log_wRR_cum_95CI[2],
            b_ppo_model_adj_results_CI_bb$log_wRR_both_95CI[2],
            order_weighted_RD_PO_adj_CI_bb$log_RR_leq_95CI[2,1:3],
            order_weighted_RD_PO_adj_CI_bb$unweighted_log_RR_95CI[2],
            order_weighted_RD_PO_adj_CI_bb$log_wRR_bcontrol_95CI[2],
            order_weighted_RD_PO_adj_CI_bb$log_wRR_cum_95CI[2],
            order_weighted_RD_PO_adj_CI_bb$log_wRR_both_95CI[2])
RR[,6] <- c(b_ppo_model_adj_results_bb$weights_split, mean(b_ppo_model_adj_results_bb$weights_split),
            mean(b_ppo_model_adj_results_bb$weights_split), mean(b_ppo_model_adj_results_bb$weights_split),
            mean(b_ppo_model_adj_results_bb$weights_split),
            order_weighted_RD_PO_adj_bb$weights_split, mean(order_weighted_RD_PO_adj_bb$weights_split),
            mean(order_weighted_RD_PO_adj_bb$weights_split), mean(order_weighted_RD_PO_adj_bb$weights_split),
            mean(order_weighted_RD_PO_adj_bb$weights_split)) * 2
RR <- as.data.frame(RR)
colnames(RR) <- c("RR", "Type", "Time", "Lower", "Upper", "Weights")
RR$Type <- factor(RR$Type, levels = c("Bayesian PPO Model", "Bayesian PO Model"))
RR$Time <- factor(RR$Time, levels = c("Nothing Bad", "Hypoxemia or Better", "ED Visit or Better", "ARR", "wRR", "wRR_cumulative", "wRR_both"))


RR_adj_14_pres3 <- RR %>% 
  filter(Time %in% c("Nothing Bad", "Hypoxemia or Better", "ED Visit or Better", "ARR", "wRR_both")) %>%
  mutate(Time2 = factor(Time, labels = c("Nothing Bad", "Hypoxemia or Better", "ED Visit or Better", "ARR", "wRR"))) %>%
  ggplot(aes(x=Time2, y=as.numeric(RR), colour=Type)) + 
  geom_errorbar(aes(ymin=as.numeric(Lower), ymax=as.numeric(Upper)), width=.2, 
                position=position_dodge(width=0.5)) +
  geom_point(position=position_dodge(width=0.5)) +  
  geom_hline(yintercept = 0, linetype = "dashed") + 
  geom_vline(xintercept = 3.5) +
  scale_color_manual(values = c(brewer.pal(8, "Paired")[2], brewer.pal(8, "Paired")[4])) +
  #scale_x_discrete(labels = label_wrap(10)) +
  ylab("log(Relative Risk)") + 
  xlab("Event") + theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_text(size = 9),
        axis.title = element_text(size = 14), legend.text = element_text(size = 9),
        legend.title = element_text(size = 14)) + 
  scale_x_discrete(labels=c("Nothing Bad" = "\u2264 1", 
                            "Hypoxemia or Better" = "\u2264 2", 
                            "ED Visit or Better" = "\u2264 3", 
                            "ARR" = "ARR", 
                            "wRR" = "wRR"))

# set up data for plot
RR <- matrix(NA, nrow = 10, ncol = 6)
RR[,1] <- c(b_ppo_model_adj_results_bb$log_RR_geq,
            b_ppo_model_adj_results_bb$unweighted_log_RR_ge,
            b_ppo_model_adj_results_bb$log_weighted_RR_ge_both,
            log(order_weighted_RD_PO_adj_bb$RR_geq), order_weighted_RD_PO_adj_bb$unweighted_log_RR_ge,
            order_weighted_RD_PO_adj_bb$log_weighted_RR_ge_both)
RR[,2] <- c(rep("Bayesian PPO Model", 5), rep("Bayesian PO Model", 5))
RR[,3] <- c("Hypoxemia or Worse", "ED Visit or Worse", "Hospitalization/Death", "ARR", "wRR_both",
            "Hypoxemia or Worse", "ED Visit or Worse", "Hospitalization/Death", "ARR", "wRR_both")
RR[,4] <- c(log(b_ppo_model_adj_results_CI_bb$RR_geq_95CI[1,1:3]), 
            b_ppo_model_adj_results_CI_bb$unweighted_log_RR_ge_95CI[1],
            b_ppo_model_adj_results_CI_bb$log_wRR_ge_both_95CI[1],
            log(order_weighted_RD_PO_adj_CI_bb$RR_geq_95CI[1,1:3]),
            order_weighted_RD_PO_adj_CI_bb$unweighted_log_RR_ge_95CI[1],
            order_weighted_RD_PO_adj_CI_bb$log_wRR_ge_both_95CI[1])
RR[,5] <- c(log(b_ppo_model_adj_results_CI_bb$RR_geq_95CI[2,1:3]),  
            b_ppo_model_adj_results_CI_bb$unweighted_log_RR_ge_95CI[2],
            b_ppo_model_adj_results_CI_bb$log_wRR_ge_both_95CI[2],
            log(order_weighted_RD_PO_adj_CI_bb$RR_geq_95CI[2,1:3]),
            order_weighted_RD_PO_adj_CI_bb$unweighted_log_RR_ge_95CI[2],
            order_weighted_RD_PO_adj_CI_bb$log_wRR_ge_both_95CI[2])
RR[,6] <- c(b_ppo_model_adj_results_bb$weights_split, mean(b_ppo_model_adj_results_bb$weights_split),
            mean(b_ppo_model_adj_results_bb$weights_split),
            order_weighted_RD_PO_adj_bb$weights_split, mean(order_weighted_RD_PO_adj_bb$weights_split),
            mean(order_weighted_RD_PO_adj_bb$weights_split)) * 2
RR <- as.data.frame(RR)
colnames(RR) <- c("RR", "Type", "Time", "Lower", "Upper", "Weights")
RR$Type <- factor(RR$Type, levels = c("Bayesian PPO Model", "Bayesian PO Model"))
RR$Time <- factor(RR$Time, levels = c("Hypoxemia or Worse", "ED Visit or Worse", "Hospitalization/Death", "ARR", "wRR_both"))



RR_ge_adj_14_pres3 <- RR %>% 
  filter(Time %in% c("Hypoxemia or Worse", "ED Visit or Worse", "Hospitalization/Death", "ARR", "wRR_both")) %>%
  mutate(Time2 = factor(Time, labels = c("Hypoxemia or Worse", "ED Visit or Worse", "Hospitalization/Death", "ARR", "wRR"))) %>%
  ggplot(aes(x=Time2, y=as.numeric(RR), colour=Type)) + 
  geom_errorbar(aes(ymin=as.numeric(Lower), ymax=as.numeric(Upper)), width=.2, 
                position=position_dodge(width=0.5)) +
  geom_point(position=position_dodge(width=0.5)) +  
  geom_hline(yintercept = 0, linetype = "dashed") + 
  geom_vline(xintercept = 3.5) +
  scale_color_manual(values = c(brewer.pal(8, "Paired")[2], brewer.pal(8, "Paired")[4])) +
  #scale_x_discrete(labels = label_wrap(10)) +
  ylab("log(Relative Risk)") + 
  xlab("Event") + theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_text(size = 9),
        axis.title = element_text(size = 14), legend.text = element_text(size = 9),
        legend.title = element_text(size = 14)) + 
  scale_x_discrete(labels=c("Hypoxemia or Worse" = "\u2265 2", 
                            "ED Visit or Worse" = "\u2265 3", 
                            "Hospitalization/Death" = "\u2265 4", 
                            "ARR" = "ARR", 
                            "wRR" = "wRR"))

RR <- matrix(NA, nrow = 10, ncol = 5)
RR[,1] <- c(b_ppo_model_adj_results_28_bb$log_RR_leq,
            b_ppo_model_adj_results_28_bb$unweighted_log_RR,
            b_ppo_model_adj_results_28_bb$log_weighted_RR_both,
            order_weighted_RD_PO_28_adj_bb$log_RR_leq, order_weighted_RD_PO_28_adj_bb$unweighted_log_RR, 
            order_weighted_RD_PO_28_adj_bb$log_weighted_RR_both)
RR[,2] <- c(rep("Bayesian PPO Model", 5), rep("Bayesian PO Model", 5))
RR[,3] <- c("Nothing Bad", "Hypoxemia or Better", "ED Visit or Better", "ARR", "wRR_both",
            "Nothing Bad", "Hypoxemia or Better", "ED Visit or Better", "ARR", "wRR_both")
RR[,4] <- c(b_ppo_model_adj_results_CI_28_bb$log_RR_leq_95CI[1,1:3], 
            b_ppo_model_adj_results_CI_28_bb$unweighted_log_RR_95CI[1],
            b_ppo_model_adj_results_CI_28_bb$log_wRR_both_95CI[1],
            order_weighted_RD_PO_28_adj_CI_bb$log_RR_leq_95CI[1,1:3],
            order_weighted_RD_PO_28_adj_CI_bb$unweighted_log_RR_95CI[1],
            order_weighted_RD_PO_28_adj_CI_bb$log_wRR_both_95CI[1])
RR[,5] <- c(b_ppo_model_adj_results_CI_28_bb$log_RR_leq_95CI[2,1:3],  
            b_ppo_model_adj_results_CI_28_bb$unweighted_log_RR_95CI[2],
            b_ppo_model_adj_results_CI_28_bb$log_wRR_both_95CI[2],
            order_weighted_RD_PO_28_adj_CI_bb$log_RR_leq_95CI[2,1:3],
            order_weighted_RD_PO_28_adj_CI_bb$unweighted_log_RR_95CI[2],
            order_weighted_RD_PO_28_adj_CI_bb$log_wRR_both_95CI[2])
RR <- as.data.frame(RR)
colnames(RR) <- c("RR", "Type", "Time", "Lower", "Upper")
RR$Type <- factor(RR$Type, levels = c("Bayesian PPO Model", "Bayesian PO Model"))
RR$Time <- factor(RR$Time, levels = c("Nothing Bad", "Hypoxemia or Better", "ED Visit or Better", "ARR", "wRR_both"))

RR_adj_28_pres3 <- RR %>% 
  filter(Time %in% c("Nothing Bad", "Hypoxemia or Better", "ED Visit or Better", "ARR", "wRR_both")) %>%
  mutate(Time2 = factor(Time, labels = c("Nothing Bad", "Hypoxemia or Better", "ED Visit or Better", "ARR", "wRR"))) %>%
  ggplot(aes(x=Time2, y=as.numeric(RR), colour=Type)) + 
  geom_errorbar(aes(ymin=as.numeric(Lower), ymax=as.numeric(Upper)), width=.2, 
                position=position_dodge(width=0.5)) +
  geom_point(position=position_dodge(width=0.5)) +  
  geom_hline(yintercept = 0, linetype = "dashed") + 
  geom_vline(xintercept = 3.5) +
  scale_color_manual(values = c(brewer.pal(8, "Paired")[2], brewer.pal(8, "Paired")[4])) +
  #scale_x_discrete(labels = label_wrap(10)) +
  ylab("log(Relative Risk)") + 
  xlab("Event") + theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_text(size = 9),
        axis.title = element_text(size = 14), legend.text = element_text(size = 9),
        legend.title = element_text(size = 14))  + 
  scale_x_discrete(labels=c("Nothing Bad" = "\u2264 1", 
                            "Hypoxemia or Better" = "\u2264 2", 
                            "ED Visit or Better" = "\u2264 3", 
                            "ARR" = "ARR", 
                            "wRR" = "wRR"))


# set up data for plot
RR <- matrix(NA, nrow = 10, ncol = 6)
RR[,1] <- c(b_ppo_model_adj_results_28_bb$log_RR_geq,
            b_ppo_model_adj_results_28_bb$unweighted_log_RR_ge,
            b_ppo_model_adj_results_28_bb$log_weighted_RR_ge_both,
            log(order_weighted_RD_PO_28_adj_bb$RR_geq),
            order_weighted_RD_PO_28_adj_bb$unweighted_log_RR_ge,
            order_weighted_RD_PO_28_adj_bb$log_weighted_RR_ge_both)
RR[,2] <- c(rep("Bayesian PPO Model", 5), rep("Bayesian PO Model", 5))
RR[,3] <- c("Hypoxemia or Worse", "ED Visit or Worse", "Hospitalization/Death", "ARR", "wRR_both",
            "Hypoxemia or Worse", "ED Visit or Worse", "Hospitalization/Death", "ARR", "wRR_both")
RR[,4] <- c(log(b_ppo_model_adj_results_CI_28_bb$RR_geq_95CI[1,1:3]), 
            b_ppo_model_adj_results_CI_28_bb$unweighted_log_RR_ge_95CI[1],
            b_ppo_model_adj_results_CI_28_bb$log_wRR_ge_both_95CI[1],
            log(order_weighted_RD_PO_28_adj_CI_bb$RR_geq_95CI[1,1:3]),
            order_weighted_RD_PO_28_adj_CI_bb$unweighted_log_RR_ge_95CI[1],
            order_weighted_RD_PO_28_adj_CI_bb$log_wRR_ge_both_95CI[1])
RR[,5] <- c(log(b_ppo_model_adj_results_CI_28_bb$RR_geq_95CI[2,1:3]),  
            b_ppo_model_adj_results_CI_28_bb$unweighted_log_RR_ge_95CI[2],
            b_ppo_model_adj_results_CI_28_bb$log_wRR_ge_both_95CI[2],
            log(order_weighted_RD_PO_28_adj_CI_bb$RR_geq_95CI[2,1:3]),
            order_weighted_RD_PO_28_adj_CI_bb$unweighted_log_RR_ge_95CI[2],
            order_weighted_RD_PO_28_adj_CI_bb$log_wRR_ge_both_95CI[2])
RR[,6] <- c(b_ppo_model_adj_results_28_bb$weights_split,
            mean(b_ppo_model_adj_results_28_bb$weights_split),
            mean(b_ppo_model_adj_results_28_bb$weights_split),
            order_weighted_RD_PO_28_adj_bb$weights_split,
            mean(order_weighted_RD_PO_28_adj_bb$weights_split),
            mean(order_weighted_RD_PO_28_adj_bb$weights_split)) * 2
RR <- as.data.frame(RR)
colnames(RR) <- c("RR", "Type", "Time", "Lower", "Upper", "Weights")
RR$Type <- factor(RR$Type, levels = c("Bayesian PPO Model", "Bayesian PO Model"))
RR$Time <- factor(RR$Time, levels = c("Hypoxemia or Worse", "ED Visit or Worse", "Hospitalization/Death", "ARR", "wRR_both"))



RR_ge_adj_28_pres3 <- RR %>% 
  filter(Time %in% c("Hypoxemia or Worse", "ED Visit or Worse", "Hospitalization/Death", "ARR", "wRR_both")) %>%
  mutate(Time2 = factor(Time, labels = c("Hypoxemia or Worse", "ED Visit or Worse", "Hospitalization/Death", "ARR", "wRR"))) %>%
  ggplot(aes(x=Time2, y=as.numeric(RR), colour=Type)) + 
  geom_errorbar(aes(ymin=as.numeric(Lower), ymax=as.numeric(Upper)), width=.2, 
                position=position_dodge(width=0.5)) +
  geom_point(position=position_dodge(width=0.5)) +  
  geom_hline(yintercept = 0, linetype = "dashed") + 
  geom_vline(xintercept = 3.5) +
  scale_color_manual(values = c(brewer.pal(8, "Paired")[2], brewer.pal(8, "Paired")[4])) +
  #scale_x_discrete(labels = label_wrap(10)) +
  ylab("log(Relative Risk)") + 
  xlab("Event") + theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_text(size = 9),
        axis.title = element_text(size = 14), legend.text = element_text(size = 9),
        legend.title = element_text(size = 14))  + 
  scale_x_discrete(labels=c("Hypoxemia or Worse" = "\u2265 2", 
                            "ED Visit or Worse" = "\u2265 3", 
                            "Hospitalization/Death" = "\u2265 4", 
                            "ARR" = "ARR", 
                            "wRR" = "wRR"))


#stacked bar plot for paper

p1_day14 <- covid_dropped_dropvax %>% 
  ggplot(aes(x = metformin, fill = as.factor(ordinal_outcome))) + 
  geom_bar(position = "fill")  + scale_fill_brewer(palette = "Blues", name = "Outcome", labels = c("Nothing Bad", "Hypoxemia", "ED Visit", "Hospitalization/Death")) + theme_bw() + ylab("Proportion") + 
  xlab("Metformin") + theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14), legend.text = element_text(size = 9),
        legend.title = element_text(size = 14), 
        legend.position = "bottom") 

p1_day28 <- covid_dropped_dropvax_28 %>% 
  ggplot(aes(x = metformin, fill = as.factor(ordinal_outcome))) + 
  geom_bar(position = "fill")  + scale_fill_brewer(palette = "Blues", name = "Outcome", labels = c("Nothing Bad", "Hypoxemia", "ED Visit", "Hospitalization/Death")) + theme_bw() + ylab("Proportion") + 
  xlab("Metformin") + theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14), legend.text = element_text(size = 9),
        legend.title = element_text(size = 14)) 



# bar with whisker for the weights 
weights <- matrix(NA, nrow = 6, ncol = 5)
weights[,1] <- c(b_ppo_model_adj_results_bb$weights_both,
                 order_weighted_RD_PO_adj_bb$weights_both)
weights[,2] <- c(b_ppo_model_adj_results_CI_bb$weights_both_95CI[1,],
                 order_weighted_RD_PO_adj_CI_bb$weights_both_95CI[1,])
weights[,3] <- c(b_ppo_model_adj_results_CI_bb$weights_both_95CI[2,],
                 order_weighted_RD_PO_adj_CI_bb$weights_both_95CI[2,])
weights[,4] <- c(rep("Bayesian PPO Model", 3), rep("Bayesian PO Model", 3))
weights[,5] <- c("1|2", "2|3", "3|4", "1|2", "2|3", "3|4")
weights <- as.data.frame(weights)
colnames(weights) <- c("Value", "Low", "High", "Model", "Levels")



weights_d14 <- weights %>%
  ggplot(aes(y = as.numeric(Value), x = Levels, fill = Model)) + 
  geom_bar(position = position_dodge(width = 0.9), stat = "identity") + 
  geom_errorbar(aes(x = Levels, ymin = as.numeric(Low), ymax = as.numeric(High)),
                width = 0.2, 
                position = position_dodge(width = 0.9)) +
  theme_bw() + ylab("Weights") + 
  xlab("Event") + theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14), legend.text = element_text(size = 9),
        legend.title = element_text(size = 14)) +
  scale_fill_manual(values = c(brewer.pal(8, "Paired")[2], brewer.pal(8, "Paired")[4]))

# bar with whisker for the weights 
weights <- matrix(NA, nrow = 6, ncol = 5)
weights[,1] <- c(b_ppo_model_adj_results_28_bb$weights_both,
                 order_weighted_RD_PO_28_adj_bb$weights_both)
weights[,2] <- c(b_ppo_model_adj_results_CI_28_bb$weights_both_95CI[1,],
                 order_weighted_RD_PO_28_adj_CI_bb$weights_both_95CI[1,])
weights[,3] <- c(b_ppo_model_adj_results_CI_28_bb$weights_both_95CI[2,],
                 order_weighted_RD_PO_28_adj_CI_bb$weights_both_95CI[2,])
weights[,4] <- c(rep("Bayesian PPO Model", 3), rep("Bayesian PO Model", 3))
weights[,5] <- c("1|2", "2|3", "3|4", "1|2", "2|3", "3|4")
weights <- as.data.frame(weights)
colnames(weights) <- c("Value", "Low", "High", "Model", "Levels")

weights_d28 <- weights %>%
  ggplot(aes(y = as.numeric(Value), x = Levels, fill = Model)) + 
  geom_bar(position = position_dodge(width = 0.9), stat = "identity") + 
  geom_errorbar(aes(x = Levels, ymin = as.numeric(Low), ymax = as.numeric(High)),
                width = 0.2, 
                position = position_dodge(width = 0.9)) +
  theme_bw() + ylab("Weights") + 
  xlab("Event") + theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14), legend.text = element_text(size = 9),
        legend.title = element_text(size = 14)) +
  scale_fill_manual(values = c(brewer.pal(8, "Paired")[2], brewer.pal(8, "Paired")[4]))

# final figures saving 
legend1 <- cowplot::get_plot_component(p1_day14, 'guide-box-bottom')

legend <- get_legend(
  # create some space to the left of the legend
  OR_adj_14_paper + theme(legend.box.margin = margin(0, 0, 0, 5))
)

plots <- plot_grid(p1_day14 + theme(legend.position = "none"), weights_d14 + theme(legend.position = "none"), OR_adj_14_pres3+ theme(legend.position="none"), labels = c("A", "B", "C"), nrow = 1)
legends <- plot_grid(legend1, legend, nrow = 1, rel_widths = c(1.75,1))
plot_grid(plots, legends, nrow =2,  rel_heights = c(3, 0.5))
#ggsave("adj_day14_040725.png", width = 9, height = 6)

# day 28 figure 
plots <- plot_grid(p1_day28 + theme(legend.position = "none"), weights_d28 + theme(legend.position = "none"), OR_adj_28_pres3+ theme(legend.position="none"), labels = c("A", "B", "C"), nrow = 1)
legends <- plot_grid(legend1, legend, nrow = 1, rel_widths = c(1.75,1))
plot_grid(plots, legends, nrow =2,  rel_heights = c(3, 0.5))
#ggsave("adj_day28_040725.png", width = 9, height = 6)


# supplement day 14 figure 
plots <- plot_grid(RD_adj_14_pres3+ theme(legend.position="none"),  
                   RR_adj_14_pres3 + theme(legend.position="none"),
                   RR_ge_adj_14_pres3 + theme(legend.position = "none"), 
                   labels = c("A", "B", "C"), 
                   nrow = 1)

legend <- cowplot::get_plot_component(RD_adj_14_paper, 'guide-box-bottom')

plot_grid(plots, legend, nrow =2,  rel_heights = c(3, 0.5))
#ggsave("adj_day14_040725_supplement.png", width = 9, height = 6)


# supplement day 28 figure 
plots <- plot_grid(RD_adj_28_pres3+ theme(legend.position="none"),  
                   RR_adj_28_pres3 + theme(legend.position="none"),
                   RR_ge_adj_28_pres3 + theme(legend.position = "none"), 
                   labels = c("A", "B", "C"), 
                   nrow = 1)

legend <- cowplot::get_plot_component(RD_adj_14_paper, 'guide-box-bottom')

plot_grid(plots, legend, nrow =2,  rel_heights = c(3, 0.5))
#ggsave("adj_day28_040725_supplement.png", width = 9, height = 6)


