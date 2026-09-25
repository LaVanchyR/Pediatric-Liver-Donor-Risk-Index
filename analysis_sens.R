#title: Analysis for Pediatric Donor Risk Index
#author: Ryan LaVanchy
#date: 8/13/2026
#purpose: Perform the analysis for the pediatric donor risk index data.


#Libraries
library(tidyverse)
library(performance)
library(naniar)
library(mice)
library(glmnet)
library(haven)
library(survival)
library(riskRegression)
library(survAUC)
library(SurvMetrics)
library(timeROC)
library(ggsurvfit)
library(gtsummary)
library(tidyr)
library(ggsci)
library(gt)
library(patchwork)
library(randomForest)
library(pROC)
library(flowchart)
library(rpsftm)
library(rpsftmPlus)



#load the data
load("L:/Projects/Investigator/Jonathan Merola/data/analysis_sens.RData")

analysis <- analysis_sens

options(scipen = 999)

##################### partition the data for training and testing ########################

#create a random row index to select from analytic dataset
set.seed(3)
random_index <- sample(
  x = 1:nrow(analysis),
  size = (0.75 * nrow(analysis)),
  replace = F
)


#select rows from analytic file based on random index to establish each dataset
data_train <- analysis[random_index, ]
data_test <- analysis[-random_index, ]




#create a predictor matrix and outcome vector for train dataset
data_train_x <- model.matrix( ~ ., select(data_train, age_at_transplant_cat, donor_age_at_tx_cat, blood_type, donor_blood_type, albumin, bmi_cat, donor_bmi_cat, cold_ischemic_time, donor_creat, creat, donor_ddavp, donor_death_circum, donor_death_mech, donor_hist_diab, diagnosis, dialysis_within_last_week, ebv_status, last_enceph, functional_status, sex, donor_sex, donor_hist_hyperten, donor_inotrop, inr, life_support, prev_malig, med_condition, ventilator, pvt, prev_ab_surg, donor_sgot, bilirubin, tipss, donor_org_shared, exception, donor_cmv_status, prev_tx, donor_arginine, donor_blood_infect, donor_bun, donor_cdc_high_risk, donor_clinical_infect, donor_cod, donor_infect_other, donor_protein_urine, donor_diuretics, donor_steroids, donor_t4, donor_infect_lung, donor_bilirubin, donor_infect_urine, donor_vasodil, waitlist_time_cat, growth_failure, tx_year, donor_type))


data_train_y <- Surv(data_train$gs_time, data_train$gs_outcome)




#create a predictor matrix and outcome vector for test dataset
data_test_x <- model.matrix( ~ ., select(data_test, age_at_transplant_cat, donor_age_at_tx_cat, blood_type, donor_blood_type, albumin, bmi_cat, donor_bmi_cat, cold_ischemic_time, donor_creat, creat, donor_ddavp, donor_death_circum, donor_death_mech, donor_hist_diab, diagnosis, dialysis_within_last_week, ebv_status, last_enceph, functional_status, sex, donor_sex, donor_hist_hyperten, donor_inotrop, inr, life_support, prev_malig, med_condition, ventilator, pvt, prev_ab_surg, donor_sgot, bilirubin, tipss, donor_org_shared, exception, donor_cmv_status, prev_tx, donor_arginine, donor_blood_infect, donor_bun, donor_cdc_high_risk, donor_clinical_infect, donor_cod, donor_infect_other, donor_protein_urine, donor_diuretics, donor_steroids, donor_t4, donor_infect_lung, donor_bilirubin, donor_infect_urine, donor_vasodil, waitlist_time_cat, growth_failure, tx_year, donor_type))


data_test_y <- Surv(data_test$gs_time, data_test$gs_outcome)




#############
#LASSO Model#
#############

model_lasso <- cv.glmnet(x = data_train_x, 
                               y = data_train_y, 
                               nfolds = 10,
                               type.measure = "deviance", # Use deviance for cross-validation
                               alpha = 1,                 # Specify LASSO penalty
                               family = "cox", 
                               nlambda = 100)

#extract min lambda
lambda_min_lasso <- model_lasso$lambda.min

#extract coefficients
coef_model_lasso <- coef(model_lasso, s = "lambda.min")

#rerun the model with the lambda min value
model_lasso_min <- glmnet(
  x = data_train_x, 
  y = data_train_y, 
  alpha = 1,                 
  family = "cox", 
  lambda = lambda_min_lasso,
  nlambda = 1
)


# Use best lambda to predict TEST data
pred_lasso <- predict(model_lasso_min, 
                            s = "lambda.min",
                            type = 'response',
                            newx = data_test_x)

#get the linear predictor from the LASSO model
lp_lasso <- predict(model_lasso_min, newx = data_test_x, s = lambda_min_lasso, type = "link")
lp_lasso <- as.numeric(lp_lasso)


#calculate the area under the curve at five years
auc_results_lasso <- timeROC(
  T = data_test$gs_time,
  delta = data_test$gs_outcome,
  marker = lp_lasso,
  cause = 1,
  weighting = "marginal",
  times = c(1825),
  iid = T
)

plot(auc_results_lasso, time = 1825, col = "blue", title = F)

lasso_auc <- 0.622



#aggregate survival probabilities from the linear predictor
lasso_cox_lp <- coxph(data_test_y ~ offset(lp_lasso))
base_surv_lasso <- survfit(lasso_cox_lp)


lasso_surv_probs <- sapply(lp_lasso, function(lp_i) {
  base_surv_lasso$surv ^ exp(lp_i)
})

lasso_surv_probs <- as.data.frame(lasso_surv_probs)

data_test_times <- as.data.frame(unique(data_test$gs_time))

data_test_times <- data_test_times %>%
  rename(time = `unique(data_test$gs_time)`) %>%
  arrange(time)

lasso_surv_probs <- lasso_surv_probs %>%
  mutate(time = data_test_times$time) %>%
  filter(time == 1826) %>%
  pivot_longer(cols = c(V1:V1559))




#C-index
model_lasso_min_cindex <- glmnet::Cindex(pred = pred_lasso, y = data_test_y)


#deviance ratio
dr_lasso <- model_lasso_min$dev.ratio


#Brier score
model_lasso_min_brier <- Brier(
  object = data_test_y,
  pre_sp = lasso_surv_probs$value,
  t_star = 1825
)




#calibration curve - calculate the observed and expected survival probability, and plot them against each other


test <- survfit(model_lasso_min, s = 0.05, x = data_test_x, y = data_test_y) 
  
pred_surv_lasso <- data.frame(
  time_pred = test$time,
  estimate_pred = test$surv
)

  
lasso_obs_surv_prob <- survfit(data_test_y ~ 1)

lasso_obs_surv_prob <- data.frame(
  time_obs = lasso_obs_surv_prob$time,
  estimate_obs = lasso_obs_surv_prob$surv
)




lasso_cc_data <- cbind(lasso_obs_surv_prob, pred_surv_lasso)

ggplot(data = lasso_cc_data) +
  geom_step(aes(x = time_pred, y = estimate_pred), color = "red") + 
  geom_step(aes(x = time_obs, y = estimate_obs), color = "blue") +
  labs(x = "Time (Days)",
       y = "Survival Probability") +
  scale_ggsurvfit(y_scales = list(breaks = seq(0, 1, by = 0.1), limits = c(0,1)),
                  x_scales = list(breaks = seq(0, 1850, by = 300), limits = c(0,1850))) +
  theme_ggsurvfit_default() 


#calculate the Kolmogorov-Smirnov Test
ks.test(lasso_cc_data$estimate_obs, lasso_cc_data$estimate_pred)
lasso_ks <- 0.01

#survival at five years
lasso_5yr_surv <- 0.883




##############
#Ridge Model#
##############


model_ridge <- cv.glmnet(x = data_train_x, 
                         y = data_train_y, 
                         nfolds = 10,
                         type.measure = "deviance", # Use deviance for cross-validation
                         alpha = 0,                 # Specify LASSO penalty
                         family = "cox", 
                         nlambda = 100)

#extract min lambda
lambda_min_ridge <- model_ridge$lambda.min

#extract coefficients
coef_model_ridge <- coef(model_ridge, s = "lambda.min")

#rerun the model with the lambda min value
model_ridge_min <- glmnet(
  x = data_train_x, 
  y = data_train_y, 
  alpha = 0,                 
  family = "cox", 
  lambda = lambda_min_ridge,
  nlambda = 1
)



# Use best lambda to predict TEST data
pred_ridge <- predict(model_ridge_min, 
                      s = "lambda.min",
                      type = 'response',
                      newx = data_test_x)




#C-index
model_ridge_min_cindex <- glmnet::Cindex(pred = pred_ridge, y = data_test_y)


#deviance ratio
dr_ridge <- model_ridge_min$dev.ratio


#Brier score
lp_ridge <- predict(model_ridge_min, newx = data_test_x, s = lambda_min_ridge, type = "link")
lp_ridge <- as.numeric(lp_ridge)


ridge_cox_lp <- coxph(data_test_y ~ offset(lp_ridge))
base_surv_ridge <- survfit(ridge_cox_lp)


ridge_surv_probs <- sapply(lp_ridge, function(lp_i) {
  base_surv_ridge$surv ^ exp(lp_i)
})

ridge_surv_probs <- as.data.frame(ridge_surv_probs)


ridge_surv_probs <- ridge_surv_probs %>%
  mutate(time = data_test_times$time) %>%
  filter(time == 1826) %>%
  pivot_longer(cols = c(V1:V1559))


model_ridge_min_brier <- Brier(
  object = data_test_y,
  pre_sp = ridge_surv_probs$value,
  t_star = 1825
)





#calculate the area under the curve at five years
auc_results_ridge <- timeROC(
  T = data_test$gs_time,
  delta = data_test$gs_outcome,
  marker = lp_ridge,
  cause = 1,
  weighting = "marginal",
  times = c(1825),
  iid = T
)

plot(auc_results_ridge, time = 1825, col = "blue", title = F)

ridge_auc <- 0.614


#calibration curve - calculate the observed and expected survival probability, and plot them against each other


pred_surv_ridge <- survfit(model_ridge_min, s = 0.05, x = data_test_x, y = data_test_y) 

pred_surv_ridge <- data.frame(
  time_pred = pred_surv_ridge$time,
  estimate_pred = pred_surv_ridge$surv
)


ridge_obs_surv_prob <- survfit(data_test_y ~ 1)

ridge_obs_surv_prob <- data.frame(
  time_obs = ridge_obs_surv_prob$time,
  estimate_obs = ridge_obs_surv_prob$surv
)




ridge_cc_data <- cbind(ridge_obs_surv_prob, pred_surv_ridge)

ggplot(data = ridge_cc_data) +
  geom_step(aes(x = time_pred, y = estimate_pred), color = "red") + 
  geom_step(aes(x = time_obs, y = estimate_obs), color = "blue") +
  labs(x = "Time (Days)",
       y = "Survival Probability") +
  scale_ggsurvfit(y_scales = list(breaks = seq(0, 1, by = 0.1), limits = c(0,1)),
                  x_scales = list(breaks = seq(0, 1850, by = 300), limits = c(0,1850))) +
  theme_ggsurvfit_default()

#calculate the Kolmogorov-Smirnov Test
ks.test(ridge_cc_data$estimate_obs, ridge_cc_data$estimate_pred)
ridge_ks <- 0.177

#survival at five years
ridge_5yr_surv <- 0.876





###################
#Elastic Net Model#
###################

model_EL <- cv.glmnet(x = data_train_x, 
                         y = data_train_y, 
                         nfolds = 10,
                         type.measure = "deviance", # Use deviance for cross-validation
                         alpha = 0.25,                 # Specify LASSO penalty
                         family = "cox", 
                         nlambda = 100)

#extract the min lambda
lambda_min_EL <- model_EL$lambda.min

#extract coefficients
coef_model_EL <- coef(model_EL, s = "lambda.min")

#rerun the model with the lambda min value
model_EL_min <- glmnet(
  x = data_train_x, 
  y = data_train_y, 
  alpha = 0.25,                 
  family = "cox", 
  lambda = lambda_min_EL,
  nlambda = 1
)


# Use best lambda to predict TEST data
pred_EL <- predict(model_EL_min, 
                      s = "lambda.min",
                      type = 'response',
                      newx = data_test_x)




#C-index
model_EL_min_cindex <- glmnet::Cindex(pred = pred_EL, y = data_test_y)


#deviance ratio
dr_EL <- model_EL_min$dev.ratio

#Brier score
lp_EL <- predict(model_EL_min, newx = data_test_x, s = lambda_min_ridge, type = "link")
lp_EL <- as.numeric(lp_EL)


EL_cox_lp <- coxph(data_test_y ~ offset(lp_EL))
base_surv_EL <- survfit(EL_cox_lp)


EL_surv_probs <- sapply(lp_EL, function(lp_i) {
  base_surv_EL$surv ^ exp(lp_i)
})

EL_surv_probs <- as.data.frame(EL_surv_probs)


EL_surv_probs <- EL_surv_probs %>%
  mutate(time = data_test_times$time) %>%
  filter(time == 1826) %>%
  pivot_longer(cols = c(V1:V1559))


model_EL_min_brier <- Brier(
  object = data_test_y,
  pre_sp = EL_surv_probs$value,
  t_star = 1825
)




#calculate the area under the curve at five years
auc_results_EL <- timeROC(
  T = data_test$gs_time,
  delta = data_test$gs_outcome,
  marker = lp_EL,
  cause = 1,
  weighting = "marginal",
  times = c(1825),
  iid = T
)

plot(auc_results_EL, time = 1825, col = "blue", title = F)

EL_auc <- 0.622


#calibration curve - calculate the observed and expected survival probability, and plot them against each other


pred_surv_EL <- survfit(model_EL_min, s = 0.05, x = data_test_x, y = data_test_y) 

pred_surv_EL <- data.frame(
  time_pred = pred_surv_EL$time,
  estimate_pred = pred_surv_EL$surv
)


EL_obs_surv_prob <- survfit(data_test_y ~ 1)

EL_obs_surv_prob <- data.frame(
  time_obs = EL_obs_surv_prob$time,
  estimate_obs = EL_obs_surv_prob$surv
)




EL_cc_data <- cbind(EL_obs_surv_prob, pred_surv_EL)

ggplot(data = EL_cc_data) +
  geom_step(aes(x = time_pred, y = estimate_pred), color = "red") + 
  geom_step(aes(x = time_obs, y = estimate_obs), color = "blue") +
  labs(x = "Time (Days)",
       y = "Survival Probability") +
  scale_ggsurvfit(y_scales = list(breaks = seq(0, 1, by = 0.1), limits = c(0,1)),
                  x_scales = list(breaks = seq(0, 1850, by = 300), limits = c(0,1850))) +
  theme_ggsurvfit_default()

#calculate the Kolmogorov-Smirnov Test
ks.test(EL_cc_data$estimate_obs, EL_cc_data$estimate_pred)
EL_ks <- 0.01

#survival at five years
EL_5yr_surv <- 0.882



######################
#Adaptive LASSO Model#
######################

best_ridge_coef <- as.numeric(coef_model_ridge)

model_adapt <- cv.glmnet(x = data_train_x, 
                      y = data_train_y, 
                      nfolds = 10,
                      type.measure = "deviance", # Use deviance for cross-validation
                      alpha = 1,                 # Specify LASSO penalty
                      family = "cox",
                      penalty.factor = 1 / abs(best_ridge_coef))

#extract the min lambda
lambda_min_adapt <- model_adapt$lambda.min

#extract coefficients
coef_model_adapt <- coef(model_adapt, s = "lambda.min")

#rerun the model with the lambda min value
model_adapt_min <- glmnet(
  x = data_train_x, 
  y = data_train_y, 
  alpha = 1,                 
  family = "cox", 
  lambda = lambda_min_adapt,
  nlambda = 1,
  penalty.factor = 1 / abs(best_ridge_coef)
)



# Use best lambda to predict TEST data
pred_adapt <- predict(model_adapt_min, 
                   s = "lambda.min",
                   type = 'response',
                   newx = data_test_x)




#C-index
model_adapt_min_cindex <- glmnet::Cindex(pred = pred_adapt, y = data_test_y)


#deviance ratio
dr_adapt <- model_adapt_min$dev.ratio

#Brier score
lp_adapt <- predict(model_adapt_min, newx = data_test_x, s = lambda_min_ridge, type = "link")
lp_adapt <- as.numeric(lp_adapt)


adapt_cox_lp <- coxph(data_test_y ~ offset(lp_adapt))
base_surv_adapt <- survfit(adapt_cox_lp)


adapt_surv_probs <- sapply(lp_adapt, function(lp_i) {
  base_surv_adapt$surv ^ exp(lp_i)
})

adapt_surv_probs <- as.data.frame(adapt_surv_probs)


adapt_surv_probs <- adapt_surv_probs %>%
  mutate(time = data_test_times$time) %>%
  filter(time == 1826) %>%
  pivot_longer(cols = c(V1:V1559))


model_adapt_min_brier <- Brier(
  object = data_test_y,
  pre_sp = adapt_surv_probs$value,
  t_star = 1825
)


#calculate the area under the curve at five years
auc_results_adapt <- timeROC(
  T = data_test$gs_time,
  delta = data_test$gs_outcome,
  marker = lp_adapt,
  cause = 1,
  weighting = "marginal",
  times = c(1825),
  iid = T
)

plot(auc_results_adapt, time = 1825, col = "blue", title = F)

adapt_auc <- 0.625


#calibration curve - calculate the observed and expected survival probability, and plot them against each other


pred_surv_adapt <- survfit(model_adapt_min, s = 0.05, x = data_test_x, y = data_test_y) 

pred_surv_adapt <- data.frame(
  time_pred = pred_surv_adapt$time,
  estimate_pred = pred_surv_adapt$surv
)


adapt_obs_surv_prob <- survfit(data_test_y ~ 1)

adapt_obs_surv_prob <- data.frame(
  time_obs = adapt_obs_surv_prob$time,
  estimate_obs = adapt_obs_surv_prob$surv
)




adapt_cc_data <- cbind(adapt_obs_surv_prob, pred_surv_adapt)

ggplot(data = adapt_cc_data) +
  geom_step(aes(x = time_pred, y = estimate_pred), color = "red") + 
  geom_step(aes(x = time_obs, y = estimate_obs), color = "blue") +
  labs(x = "Time (Days)",
       y = "Survival Probability") +
  scale_ggsurvfit(y_scales = list(breaks = seq(0, 1, by = 0.1), limits = c(0,1)),
                  x_scales = list(breaks = seq(0, 1850, by = 300), limits = c(0,1850))) +
  theme_ggsurvfit_default() 

#calculate the Kolmogorov-Smirnov Test
ks.test(adapt_cc_data$estimate_obs, adapt_cc_data$estimate_pred)
adapt_ks <- 0.007

#survival at five years
adapt_5yr_surv <- 0.883



#initial model metrics to decide the best fitted model
initial_model_metrics <- data.frame(
  model = c("LASSO","Ridge","Elastic Net","Adaptive LASSO"),
  c_index = c(model_lasso_min_cindex, model_ridge_min_cindex, model_EL_min_cindex, model_adapt_min_cindex),
  brier_score = c(model_lasso_min_brier, model_ridge_min_brier, model_EL_min_brier, model_adapt_min_brier),
  dr = c(dr_lasso, dr_ridge, dr_EL, dr_adapt),
  auc = c(lasso_auc, ridge_auc, EL_auc, adapt_auc),
  ks = c(lasso_ks, ridge_ks, EL_ks, adapt_ks),
  surv_5yr = c(lasso_5yr_surv, ridge_5yr_surv, EL_5yr_surv, adapt_5yr_surv)
)

#Based on the Brier score, the elastic net model is the best model


###################################################################
#create a non-regularized cox model based on the best fitted model#
###################################################################

#Note - If one dummy variable from the regularized model was included, the entire factor was included in this model. 

coef_model_EL

cox_mle_model <- coxph(Surv(gs_time, gs_outcome) ~ age_at_transplant_cat + donor_age_at_tx_cat + blood_type + donor_blood_type + albumin + bmi_cat + donor_bmi_cat + cold_ischemic_time + donor_creat + creat + donor_ddavp + donor_death_circum + donor_death_mech + donor_hist_diab + diagnosis + dialysis_within_last_week + ebv_status + last_enceph + functional_status + donor_hist_hyperten + inr + prev_malig + med_condition + pvt + prev_ab_surg + donor_sgot + bilirubin + tipss + donor_org_shared + exception + prev_tx + donor_arginine + donor_blood_infect + donor_bun + donor_cdc_high_risk + donor_cod + donor_infect_other + donor_protein_urine + donor_diuretics + donor_steroids + donor_t4 + donor_infect_lung + donor_bilirubin + donor_infect_urine + waitlist_time_cat + growth_failure + tx_year + donor_type, data = data_train)


pretend.lm <- lm(gs_time ~ age_at_transplant_cat + donor_age_at_tx_cat + blood_type + donor_blood_type + albumin + bmi_cat + donor_bmi_cat + cold_ischemic_time + donor_creat + creat + donor_ddavp + donor_death_circum + donor_death_mech + donor_hist_diab + diagnosis + dialysis_within_last_week + ebv_status + last_enceph + functional_status + donor_hist_hyperten + inr + life_support + prev_malig + med_condition + ventilator + pvt + prev_ab_surg + donor_sgot + bilirubin + tipss + donor_org_shared + exception + prev_tx + donor_arginine + donor_blood_infect + donor_bun + donor_cdc_high_risk + donor_cod + donor_infect_other + donor_protein_urine + donor_diuretics + donor_steroids + donor_t4 + donor_infect_lung + donor_bilirubin + donor_infect_urine + waitlist_time_cat + growth_failure + tx_year + donor_type, data = data_train)

vif_lm <- car::vif(pretend.lm)



saveRDS(cox_mle_model, "L:/Projects/Investigator/Jonathan Merola/data/cox_mle_model_sens")

cox_mle_model_coefs <- as.data.frame(coef(cox_mle_model))

summary(cox_mle_model)


#predict using the test data set
pred_mle <- survfit(cox_mle_model, newdata = data_test)

#create a survival prob dataset
surv_prob_1825 <- as.data.frame(pred_mle$surv)

surv_prob_1825 <- surv_prob_1825 %>%
  slice(305) %>%
  pivot_longer(cols = everything(), names_to = c("time")) %>%
  rename(surv_prob_1825 = value) %>%
  select(surv_prob_1825)
  

#C-index
mle_cindex <- concordance(cox_mle_model, newdata = data_test)$concordance



#Brier score

mle_brier_score <- Brier(object = data_test_y,
              pre_sp = surv_prob_1825$surv_prob_1825,
              t_star = 1825)



#calculate the area under the curve at five years
mle_lp <- predict(cox_mle_model, newdata = data_test, type = "lp")

auc_results_mle <- timeROC(
  T = data_test$gs_time,
  delta = data_test$gs_outcome,
  marker = mle_lp,
  cause = 1,
  weighting = "marginal",
  times = c(1825),
  iid = T
)

plot(auc_results_mle, time = 1825, col = "blue", title = F)

mle_auc <- 0.602



########################################
#create a marginalized donor risk model#
########################################

#select just the donor factors
mle_cox_vars <- names(cox_mle_model$xlevels)
mle_cox_donor_factors <- c(grep("donor", mle_cox_vars, value = T), "cold_ischemic_time")

mle_cox_recip_factors <- names(cox_mle_model$xlevels)[!names(cox_mle_model$xlevels) %in% c(mle_cox_donor_factors)]


#obtain vector of means of expected values for each observation in the training dataset
t0 <- 365 * 5
base_sf <- survfit(cox_mle_model)
S0_t0 <- summary(base_sf, times = t0, extend = TRUE)$surv


#function to calculate average 5-year predicted survival probability from the Cox model holding the donor factors fixed to a patients values while allowing recipient factors to vary.
standardized_risk_5yr <- function(data,
                                  model,
                                  donor_vars,
                                  time){
  
  # Baseline survival at specified time
  base_sf <- survfit(model)
  
  S0 <- summary(
    base_sf,
    times = time,
    extend = TRUE
  )$surv
  
  std_risk <- numeric(nrow(data))
  
  for(i in seq_len(nrow(data))){
    
    cat("Row", i, "of", nrow(data), "\n")
    
    data_i <- data
    
    # Fix donor factors to subject i
    data_i[donor_vars] <-
      data_i[rep(i, nrow(data)), donor_vars]
    
    # Linear predictors
    lp <- predict(
      model,
      newdata = data_i,
      type = "lp"
    )
    
    # Survival probabilities at 5 years
    surv_prob <- S0 ^ exp(lp)
    
    # Event probabilities
    event_prob <- 1 - surv_prob
    
    # Standardized risk
    std_risk[i] <- mean(event_prob, na.rm = TRUE)
  }
  
  std_risk
}

std_risk_5yrs <- standardized_risk_5yr(
  data = data_train,
  model = cox_mle_model,
  donor_vars = mle_cox_donor_factors,
  time = 1825
)


data_train_std <- data_train %>%
  mutate(
    std_risk_5yr = std_risk_5yrs
  )

eps <- 1e-6


#logit transform the marginalized donor coefficients 
data_train_std <- data_train_std %>%
  mutate(logit_std_risk = qlogis(pmin(pmax(std_risk_5yr, eps), 1 - eps)))


#regress the donor factors on the logit transformed standardized risk to obtain marginalized donor coefficients.
model_std_lm <- lm(logit_std_risk ~ donor_age_at_tx_cat + donor_blood_type + donor_bmi_cat + cold_ischemic_time + donor_creat + donor_ddavp + donor_death_circum + donor_death_mech + donor_hist_diab + donor_hist_hyperten + donor_sgot + donor_org_shared + donor_arginine + donor_blood_infect + donor_bun + donor_cdc_high_risk + donor_cod + donor_infect_other + donor_protein_urine + donor_diuretics + donor_t4 + donor_infect_lung + donor_bilirubin + donor_infect_urine + donor_type, data = data_train_std)

test <- plogis(predict(model_std_lm, newdata = data_test))

summary(model_std_lm)


coef_model_std_lm <- as.data.frame(coef(model_std_lm))

coef_model_std_lm_CI <- confint(model_std_lm)






