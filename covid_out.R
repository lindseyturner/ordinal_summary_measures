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
set.seed(12345)

# set working directory 
setwd("/Users/lindseyturner/Library/CloudStorage/Box-Box/COVID_OUT_Data/")
# load in the file with the functions
source('covid_out_functions.R')

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
