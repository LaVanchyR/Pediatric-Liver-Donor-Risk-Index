#title: Analysis for Pediatric Donor Risk Index
#author: Ryan LaVanchy
#date: 8/13/2026
#purpose: Perform the analysis for the pediatric donor risk index project.

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
library(geomtextpath)
library(forcats)



#load the data
load("L:/Projects/Investigator/Jonathan Merola/data/analysis3.RData")

analysis <- analysis3

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


#write out the data_train data for the Pseudo R-squared analysis
save(
  data_train,
  file = "L:/Projects/Investigator/Jonathan Merola/data/data_train.RData"
)


#create a predictor matrix and outcome vector for train dataset
data_train_x <- model.matrix(
  ~.,
  select(
    data_train,
    age_at_transplant_cat,
    donor_age_at_tx_cat,
    blood_type,
    donor_blood_type,
    albumin,
    bmi_cat,
    donor_bmi_cat,
    cold_ischemic_time,
    donor_creat,
    creat,
    donor_ddavp,
    donor_death_circum,
    donor_death_mech,
    donor_hist_diab,
    diagnosis,
    dialysis_within_last_week,
    ebv_status,
    last_enceph,
    functional_status,
    sex,
    donor_sex,
    donor_hist_hyperten,
    donor_inotrop,
    inr,
    life_support,
    prev_malig,
    med_condition,
    ventilator,
    pvt,
    prev_ab_surg,
    donor_sgot,
    bilirubin,
    tipss,
    donor_org_shared,
    exception,
    donor_cmv_status,
    prev_tx,
    donor_arginine,
    donor_blood_infect,
    donor_bun,
    donor_cdc_high_risk,
    donor_clinical_infect,
    donor_cod,
    donor_infect_other,
    donor_protein_urine,
    donor_diuretics,
    donor_steroids,
    donor_t3,
    donor_t4,
    donor_infect_lung,
    donor_bilirubin,
    donor_infect_urine,
    donor_vasodil,
    waitlist_time_cat,
    growth_failure,
    tx_year,
    donor_type
  )
)


data_train_y <- Surv(data_train$gs_time, data_train$gs_outcome)


#create a predictor matrix and outcome vector for test dataset
data_test_x <- model.matrix(
  ~.,
  select(
    data_test,
    age_at_transplant_cat,
    donor_age_at_tx_cat,
    blood_type,
    donor_blood_type,
    albumin,
    bmi_cat,
    donor_bmi_cat,
    cold_ischemic_time,
    donor_creat,
    creat,
    donor_ddavp,
    donor_death_circum,
    donor_death_mech,
    donor_hist_diab,
    diagnosis,
    dialysis_within_last_week,
    ebv_status,
    last_enceph,
    functional_status,
    sex,
    donor_sex,
    donor_hist_hyperten,
    donor_inotrop,
    inr,
    life_support,
    prev_malig,
    med_condition,
    ventilator,
    pvt,
    prev_ab_surg,
    donor_sgot,
    bilirubin,
    tipss,
    donor_org_shared,
    exception,
    donor_cmv_status,
    prev_tx,
    donor_arginine,
    donor_blood_infect,
    donor_bun,
    donor_cdc_high_risk,
    donor_clinical_infect,
    donor_cod,
    donor_infect_other,
    donor_protein_urine,
    donor_diuretics,
    donor_steroids,
    donor_t3,
    donor_t4,
    donor_infect_lung,
    donor_bilirubin,
    donor_infect_urine,
    donor_vasodil,
    waitlist_time_cat,
    growth_failure,
    tx_year,
    donor_type
  )
)


data_test_y <- Surv(data_test$gs_time, data_test$gs_outcome)


#############
#LASSO Model#
#############

#run the LASSO model
model_lasso <- cv.glmnet(
  x = data_train_x,
  y = data_train_y,
  nfolds = 10,
  type.measure = "deviance", # Use deviance for cross-validation
  alpha = 1, # Specify LASSO penalty
  family = "cox",
  nlambda = 100
)

#extract the minimum Lamda
lambda_min_lasso <- model_lasso$lambda.min

#extract coefficients
coef_model_lasso <- coef(model_lasso, s = "lambda.min")

#rerun the model with the lambda min lasso value
model_lasso_min <- glmnet(
  x = data_train_x,
  y = data_train_y,
  alpha = 1,
  family = "cox",
  lambda = lambda_min_lasso,
  nlambda = 1
)



# Use best lambda to predict TEST data
pred_lasso <- predict(
  model_lasso_min,
  s = "lambda.min",
  type = 'response',
  newx = data_test_x
)


#get the linear predictor from the LASSO model
lp_lasso <- predict(
  model_lasso_min,
  newx = data_test_x,
  s = lambda_min_lasso,
  type = "link"
)

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

lasso_auc <- 0.630


#generate survival probabilities from the linear predictor
lasso_cox_lp <- coxph(data_test_y ~ offset(lp_lasso))
base_surv_lasso <- survfit(lasso_cox_lp)


lasso_surv_probs <- sapply(lp_lasso, function(lp_i) {
  base_surv_lasso$surv^exp(lp_i)
})

lasso_surv_probs <- as.data.frame(lasso_surv_probs)

data_test_times <- as.data.frame(unique(data_test$gs_time))

data_test_times <- data_test_times %>%
  rename(time = `unique(data_test$gs_time)`) %>%
  arrange(time)

lasso_surv_probs <- lasso_surv_probs %>%
  mutate(time = data_test_times$time) %>%
  filter(time == 1826) %>%
  pivot_longer(cols = c(V1:V2104))




#C-index
model_lasso_min_cindex <- glmnet::Cindex(pred = pred_lasso, y = data_test_y)


#Deviance ratio
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
  labs(x = "Time (Days)", y = "Survival Probability") +
  scale_ggsurvfit(
    y_scales = list(breaks = seq(0, 1, by = 0.1), limits = c(0, 1)),
    x_scales = list(breaks = seq(0, 1850, by = 300), limits = c(0, 1850))
  ) +
  theme_ggsurvfit_default() +
  geom_label(
    data = data.frame(
      x = 750,
      y = 0.40,
      label = "1-Year Observed: 90.3%\n 1-Year Predicted: 91.1%\n 3-Year Observed: 86.6%\n 3-Year Predicted: 87.6%\n 5-Year Observed: 84.7%\n 5-Year Expected: 86.6%"
    ),
    aes(x = x, y = y, label = label),
    hjust = 0,
    size = 3.5
  )


#calculate the Kolmogorov-Smirnov Test
ks.test(lasso_cc_data$estimate_obs, lasso_cc_data$estimate_pred)
lasso_ks <- 0.006

#survival at five years
lasso_5yr_surv <- 0.872




##############
#Ridge Model#
##############

#run the ridge model
model_ridge <- cv.glmnet(
  x = data_train_x,
  y = data_train_y,
  nfolds = 10,
  type.measure = "deviance", # Use deviance for cross-validation
  alpha = 0, # Specify LASSO penalty
  family = "cox",
  nlambda = 100
)

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
pred_ridge <- predict(
  model_ridge_min,
  s = "lambda.min",
  type = 'response',
  newx = data_test_x
)


#C-Index
model_ridge_min_cindex <- glmnet::Cindex(pred = pred_ridge, y = data_test_y)


#Deviance ratio
dr_ridge <- model_ridge_min$dev.ratio


#Brier score
lp_ridge <- predict(
  model_ridge_min,
  newx = data_test_x,
  s = lambda_min_ridge,
  type = "link"
)
lp_ridge <- as.numeric(lp_ridge)


ridge_cox_lp <- coxph(data_test_y ~ offset(lp_ridge))
base_surv_ridge <- survfit(ridge_cox_lp)


ridge_surv_probs <- sapply(lp_ridge, function(lp_i) {
  base_surv_ridge$surv^exp(lp_i)
})

ridge_surv_probs <- as.data.frame(ridge_surv_probs)


ridge_surv_probs <- ridge_surv_probs %>%
  mutate(time = data_test_times$time) %>%
  filter(time == 1826) %>%
  pivot_longer(cols = c(V1:V2104))


model_ridge_min_brier <- Brier(
  object = data_test_y,
  pre_sp = ridge_surv_probs$value,
  t_star = 1825
)


#AUROC
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

ridge_auc <- 0.640


#calibration curve - calculate the observed and expected survival probability, and plot them against each other

pred_surv_ridge <- survfit(
  model_ridge_min,
  s = 0.05,
  x = data_test_x,
  y = data_test_y
)

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
  labs(x = "Time (Days)", y = "Survival Probability") +
  scale_ggsurvfit(
    y_scales = list(breaks = seq(0, 1, by = 0.1), limits = c(0, 1)),
    x_scales = list(breaks = seq(0, 1850, by = 300), limits = c(0, 1850))
  ) +
  theme_ggsurvfit_default() +
  geom_label(
    data = data.frame(
      x = 750,
      y = 0.40,
      label = "1-Year Observed: 90.3%\n 1-Year Predicted: 90.6%\n 3-Year Observed: 86.6%\n 3-Year Predicted: 87%\n 5-Year Observed: 84.7%\n 5-Year Expected: 85.2%"
    ),
    aes(x = x, y = y, label = label),
    hjust = 0,
    size = 3.5
  )

#calculate the Kolmogorov-Smirnov Test
ks.test(ridge_cc_data$estimate_obs, ridge_cc_data$estimate_pred)
ridge_ks <- 0.039

#survival at five years
ridge_5yr_surv <- 0.868


###################
#Elastic Net Model#
###################

#run the elastic net model
model_EL <- cv.glmnet(
  x = data_train_x,
  y = data_train_y,
  nfolds = 10,
  type.measure = "deviance", # Use deviance for cross-validation
  alpha = 0.25, # Specify LASSO penalty
  family = "cox",
  nlambda = 100
)

#extract the minimum lambda
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
pred_EL <- predict(
  model_EL_min,
  s = "lambda.min",
  type = 'response',
  newx = data_test_x
)


#C-Index
model_EL_min_cindex <- glmnet::Cindex(pred = pred_EL, y = data_test_y)

#Deviance ratio
dr_EL <- model_EL_min$dev.ratio

#Brier score
lp_EL <- predict(
  model_EL_min,
  newx = data_test_x,
  s = lambda_min_ridge,
  type = "link"
)
lp_EL <- as.numeric(lp_EL)


EL_cox_lp <- coxph(data_test_y ~ offset(lp_EL))
base_surv_EL <- survfit(EL_cox_lp)


EL_surv_probs <- sapply(lp_EL, function(lp_i) {
  base_surv_EL$surv^exp(lp_i)
})

EL_surv_probs <- as.data.frame(EL_surv_probs)


EL_surv_probs <- EL_surv_probs %>%
  mutate(time = data_test_times$time) %>%
  filter(time == 1826) %>%
  pivot_longer(cols = c(V1:V2104))


model_EL_min_brier <- Brier(
  object = data_test_y,
  pre_sp = EL_surv_probs$value,
  t_star = 1825
)


#
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

EL_auc <- 0.632


#calibration curve - calculate the observed and expected survival probability, and plot them against each other

pred_surv_EL <- survfit(
  model_EL_min,
  s = 0.05,
  x = data_test_x,
  y = data_test_y
)

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
  labs(x = "Time (Days)", y = "Survival Probability") +
  scale_ggsurvfit(
    y_scales = list(breaks = seq(0, 1, by = 0.1), limits = c(0, 1)),
    x_scales = list(breaks = seq(0, 1850, by = 300), limits = c(0, 1850))
  ) +
  theme_ggsurvfit_default() +
  geom_label(
    data = data.frame(
      x = 750,
      y = 0.40,
      label = "1-Year Observed: 90.3%\n 1-Year Predicted: 91%\n 3-Year Observed: 86.6%\n 3-Year Predicted: 87.5%\n 5-Year Observed: 84.7%\n 5-Year Expected: 85.9%"
    ),
    aes(x = x, y = y, label = label),
    hjust = 0,
    size = 3.5
  )

#calculate the Kolmogorov-Smirnov Test
ks.test(EL_cc_data$estimate_obs, EL_cc_data$estimate_pred)
EL_ks <- 0.029

#survival at five years
EL_5yr_surv <- 0.869


######################
#Adaptive LASSO Model#
######################

best_ridge_coef <- as.numeric(coef_model_ridge)

#run the adaptive LASSO model
model_adapt <- cv.glmnet(
  x = data_train_x,
  y = data_train_y,
  nfolds = 10,
  type.measure = "deviance", # Use deviance for cross-validation
  alpha = 1, # Specify LASSO penalty
  family = "cox",
  penalty.factor = 1 / abs(best_ridge_coef)
)

#extract the minimum lambda
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
pred_adapt <- predict(
  model_adapt_min,
  s = "lambda.min",
  type = 'response',
  newx = data_test_x
)


#C-Index
model_adapt_min_cindex <- glmnet::Cindex(pred = pred_adapt, y = data_test_y)


#Deviance ratio
dr_adapt <- model_adapt_min$dev.ratio

#Brier score
lp_adapt <- predict(
  model_adapt_min,
  newx = data_test_x,
  s = lambda_min_ridge,
  type = "link"
)
lp_adapt <- as.numeric(lp_adapt)


adapt_cox_lp <- coxph(data_test_y ~ offset(lp_adapt))
base_surv_adapt <- survfit(adapt_cox_lp)


adapt_surv_probs <- sapply(lp_adapt, function(lp_i) {
  base_surv_adapt$surv^exp(lp_i)
})

adapt_surv_probs <- as.data.frame(adapt_surv_probs)


adapt_surv_probs <- adapt_surv_probs %>%
  mutate(time = data_test_times$time) %>%
  filter(time == 1826) %>%
  pivot_longer(cols = c(V1:V2104))


model_adapt_min_brier <- Brier(
  object = data_test_y,
  pre_sp = adapt_surv_probs$value,
  t_star = 1825
)

#AUROC
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

adapt_auc <- 0.626


#calibration curve - calculate the observed and expected survival probability, and plot them against each other

pred_surv_adapt <- survfit(
  model_adapt_min,
  s = 0.05,
  x = data_test_x,
  y = data_test_y
)

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
  labs(x = "Time (Days)", y = "Survival Probability") +
  scale_ggsurvfit(
    y_scales = list(breaks = seq(0, 1, by = 0.1), limits = c(0, 1)),
    x_scales = list(breaks = seq(0, 1850, by = 300), limits = c(0, 1850))
  ) +
  theme_ggsurvfit_default() +
  geom_label(
    data = data.frame(
      x = 750,
      y = 0.40,
      label = "1-Year Observed: 90.3%\n 1-Year Predicted: 91.6%\n 3-Year Observed: 86.6%\n 3-Year Predicted: 88.3%\n 5-Year Observed: 84.7%\n 5-Year Expected: 86.8%"
    ),
    aes(x = x, y = y, label = label),
    hjust = 0,
    size = 3.5
  )

#calculate the Kolmogorov-Smirnov Test
ks.test(adapt_cc_data$estimate_obs, adapt_cc_data$estimate_pred)
adapt_ks <- 0.000

#survival at five years
adapt_5yr_surv <- 0.877


#initial model metrics to decide the best fitted model
initial_model_metrics <- data.frame(
  model = c("LASSO", "Ridge", "Elastic Net", "Adaptive LASSO"),
  c_index = c(
    model_lasso_min_cindex,
    model_ridge_min_cindex,
    model_EL_min_cindex,
    model_adapt_min_cindex
  ),
  brier_score = c(
    model_lasso_min_brier,
    model_ridge_min_brier,
    model_EL_min_brier,
    model_adapt_min_brier
  ),
  dr = c(dr_lasso, dr_ridge, dr_EL, dr_adapt),
  auc = c(lasso_auc, ridge_auc, EL_auc, adapt_auc),
  ks = c(lasso_ks, ridge_ks, EL_ks, adapt_ks),
  surv_5yr = c(lasso_5yr_surv, ridge_5yr_surv, EL_5yr_surv, adapt_5yr_surv)
)

#Based on the Brier score, the elastic net model is the best model

###################################################################
#create a non-regularized cox model based on the best fitted model#
###################################################################

coef_model_EL

#test collinearity - you can use the vif command for cox regression since the concern is with the relationship among independent variables, the functional form of the model for the dependent variable is irrelevant to the estimate of collinearity. Here I use linear regression so that it works well with the vif function with the car package.
pretend.lm <- lm(
  gs_time ~ age_at_transplant_cat +
    donor_age_at_tx_cat +
    blood_type +
    albumin +
    bmi_cat +
    donor_bmi_cat +
    cold_ischemic_time +
    donor_creat +
    creat +
    donor_death_circum +
    donor_death_mech +
    donor_hist_diab +
    diagnosis +
    dialysis_within_last_week +
    ebv_status +
    functional_status +
    inr +
    life_support +
    prev_malig +
    med_condition +
    ventilator +
    pvt +
    bilirubin +
    tipss +
    donor_org_shared +
    prev_tx +
    donor_arginine +
    donor_cdc_high_risk +
    donor_cod +
    donor_protein_urine +
    donor_bilirubin +
    donor_infect_urine +
    growth_failure +
    tx_year +
    donor_type,
  data = data_train
)

vif_lm <- car::vif(pretend.lm)


#remove ventilator, life support, and donor cause of death
#Note - If one dummy variable from the regularized model was included, the entire factor was included in this model.
cox_mle_model <- coxph(
  Surv(gs_time, gs_outcome) ~ age_at_transplant_cat +
    donor_age_at_tx_cat +
    blood_type +
    albumin +
    bmi_cat +
    donor_bmi_cat +
    cold_ischemic_time +
    donor_creat +
    creat +
    donor_death_circum +
    donor_death_mech +
    donor_hist_diab +
    diagnosis +
    dialysis_within_last_week +
    ebv_status +
    functional_status +
    inr +
    prev_malig +
    med_condition +
    pvt +
    bilirubin +
    tipss +
    donor_org_shared +
    prev_tx +
    donor_arginine +
    donor_cdc_high_risk +
    donor_protein_urine +
    donor_bilirubin +
    donor_infect_urine +
    growth_failure +
    tx_year +
    donor_type,
  data = data_train
)


test <- car::vif(cox_mle_model)

#test the proportional hazards assumption
prop_hazards_assmpt <- cox.zph(cox_mle_model)
plot(prop_hazards_assmpt)
#although the overall global test says that there is a violation of the proportional hazards assumption the Schoenfeld residuals reveal that there it is only a minor violation at best

#calculate the median follow-up time in both the train and test datasets
summary(data_train$gs_time)
summary(data_test$gs_time)
# the median follow-up time for both datasets was five years.

saveRDS(
  cox_mle_model,
  "L:/Projects/Investigator/Jonathan Merola/data/cox_mle_model"
)

cox_mle_model_coefs <- as.data.frame(coef(cox_mle_model))

summary(cox_mle_model)


#predict using the test data set
pred_mle <- survfit(cox_mle_model, newdata = data_test)

#create a survival prob dataset
surv_prob_1825 <- as.data.frame(pred_mle$surv)

surv_prob_1825 <- surv_prob_1825 %>%
  slice(446) %>%
  pivot_longer(cols = everything(), names_to = c("time")) %>%
  rename(surv_prob_1825 = value) %>%
  select(surv_prob_1825)


#C-Index
mle_cindex <- concordance(cox_mle_model, newdata = data_test)$concordance




#Brier score

mle_brier_score <- Brier(
  object = data_test_y,
  pre_sp = surv_prob_1825$surv_prob_1825,
  t_star = 1825
)


#AUROC
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

mle_auc <- 0.624


########################################
#create a marginalized donor risk model#
########################################

#select just the donor factors
mle_cox_vars <- names(cox_mle_model$xlevels)
mle_cox_donor_factors <- c(
  grep("donor", mle_cox_vars, value = T),
  "cold_ischemic_time"
)

mle_cox_recip_factors <- names(cox_mle_model$xlevels)[
  !names(cox_mle_model$xlevels) %in% c(mle_cox_donor_factors)
]


#obtain vector of means of expected values for each observation in the training dataset
t0 <- 365 * 5
base_sf <- survfit(cox_mle_model)
S0_t0 <- summary(base_sf, times = t0, extend = TRUE)$surv


#function to calculate average 5-year predicted survival probability from the Cox model holding the donor factors fixed to a patients values while allowing recipient factors to vary.
standardized_risk_5yr <- function(data, model, donor_vars, time) {
  # Baseline survival at specified time
  base_sf <- survfit(model)

  S0 <- summary(
    base_sf,
    times = time,
    extend = TRUE
  )$surv

  std_risk <- numeric(nrow(data))

  for (i in seq_len(nrow(data))) {
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
    surv_prob <- S0^exp(lp)

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
model_std_lm <- lm(
  logit_std_risk ~ donor_age_at_tx_cat +
    donor_bmi_cat +
    cold_ischemic_time +
    donor_creat +
    donor_death_circum +
    donor_death_mech +
    donor_hist_diab +
    donor_org_shared +
    donor_arginine +
    donor_cdc_high_risk +
    donor_protein_urine +
    donor_bilirubin +
    donor_infect_urine +
    donor_type,
  data = data_train_std
)

test <- plogis(predict(model_std_lm, newdata = data_test))

summary(model_std_lm)


coef_model_std_lm <- as.data.frame(coef(model_std_lm))

coef_model_std_lm_CI <- confint(model_std_lm)


######################################################
#plot the time-dependent AUC curve for the best model#
######################################################

EL_auc_data <- data.frame(
  TP = auc_results_EL$TP,
  FP = auc_results_EL$FP
)

mle_auc_data <- data.frame(
  TP = auc_results_mle$TP,
  FP = auc_results_mle$FP
)

fig_auc_final_mle <- ggplot(
  data = mle_auc_data,
  aes(x = FP.t.1825, y = TP.t.1825)
) +
  geom_smooth(size = 0.5, col = 'gray', se = FALSE) +
  geom_point(size = .2, alpha = 0.5) +
  theme_bw() +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  geom_abline(col = 'red') +
  labs(x = "1 - Specificity", y = "Sensitivity") +
  geom_label(
    data = data.frame(x = 0.1, y = 0.95, label = "AUC at 5-Years: 0.624"),
    aes(x = x, y = y, label = label),
    hjust = 0,
    size = 4.5
  )


ggsave(
  filename = "L:/Projects/Investigator/Jonathan Merola/graphs and tables/fig_auc_final_mle.png",
  plot = fig_auc_final_mle,
  device = "png",
  dpi = "print",
  height = 6,
  width = 9
)

##################################################################################
#create a table of the baseline characteristics for the training vs. test dataset#
##################################################################################

baseline_train <- data_train %>%
  mutate(group = "Training Data")

baseline_test <- data_test %>%
  mutate(group = "Test Data")

baseline <- rbind(baseline_train, baseline_test)

baseline$waitlist_time_cat <- factor(
  baseline$waitlist_time_cat,
  levels = c("0-6", "6-12", "12-24", ">=24")
)
baseline$tx_year <- factor(
  baseline$tx_year,
  levels = c(
    "2004",
    "2005",
    "2006",
    "2007",
    "2008",
    "2009",
    "2010",
    "2011",
    "2012",
    "2013",
    "2014",
    "2015",
    "2016",
    "2017",
    "2018",
    "2019",
    "2020",
    "2021",
    "2022",
    "2023",
    "2024",
    "2025"
  )
)

theme_gtsummary_compact()
baseline_table <- tbl_summary(
  data = baseline,
  by = group,
  include = c(
    age_at_transplant_cat,
    donor_age_at_tx_cat,
    donor_bmi_cat,
    blood_type,
    donor_blood_type,
    albumin,
    bmi_cat,
    cold_ischemic_time,
    donor_creat,
    creat,
    donor_ddavp,
    donor_death_circum,
    donor_death_mech,
    donor_hist_diab,
    diagnosis,
    dialysis_within_last_week,
    ebv_status,
    last_enceph,
    functional_status,
    sex,
    donor_sex,
    donor_hist_hyperten,
    donor_inotrop,
    inr,
    life_support,
    prev_malig,
    med_condition,
    ventilator,
    pvt,
    prev_ab_surg,
    donor_sgot,
    bilirubin,
    tipss,
    donor_org_shared,
    exception,
    donor_cmv_status,
    prev_tx,
    donor_arginine,
    donor_blood_infect,
    donor_bun,
    donor_cdc_high_risk,
    donor_clinical_infect,
    donor_cod,
    donor_infect_other,
    donor_protein_urine,
    donor_diuretics,
    donor_steroids,
    donor_t3,
    donor_t4,
    donor_infect_lung,
    donor_bilirubin,
    donor_infect_urine,
    donor_vasodil,
    waitlist_time_cat,
    growth_failure,
    group,
    tx_year,
    donor_type,
    gs_time
  ),
  missing = "ifany",
  missing_text = "Missing",
  label = list(
    age_at_transplant_cat = "Age at Transplant",
    donor_age_at_tx_cat = "Donor Age",
    donor_bmi_cat = "Donor BMI (kg/m2)",
    blood_type = "Recipient Blood Type",
    donor_blood_type = "Donor Blood Type",
    albumin = "Recipient Last Albumin",
    bmi_cat = "Recipient BMI (kg/m2)",
    cold_ischemic_time = "Cold Ischemic Time (Hours)",
    donor_creat = "Donor Terminal Creatinine (mg/dl)",
    creat = "Recipient Last Creatinine (mg/dl)",
    donor_ddavp = "Meds Given to Donor: DDAVP",
    donor_death_circum = "Donor Death Circumstance",
    donor_death_mech = "Donor Mechanism of Death",
    donor_hist_diab = "Donor History of Diabetes",
    diagnosis = "Recipient Primary Diagnosis",
    dialysis_within_last_week = "Recipient Dialysis within Prior Week",
    ebv_status = "Recipient Epstein-Barr Virus (EBV) Status",
    last_enceph = "Recipient Last Encephalopathy",
    functional_status = "Recipient Functional Status",
    sex = "Recipient Sex",
    donor_sex = "Donor Sex",
    donor_hist_hyperten = "Donor History of Hypertension",
    donor_inotrop = "Donor on Inotropic Support",
    inr = "Recipient INR",
    life_support = "Recipient on Life Support",
    prev_malig = "Recipient Previous Malignancy",
    med_condition = "Recipient Medical Condition",
    ventilator = "Recipient on Ventilator",
    pvt = "Recipient History of Portal Vein Thrombosis",
    prev_ab_surg = "Recipient Previous Abdominal Surgery",
    donor_sgot = "Donor SGOT/AST (U/L)",
    bilirubin = "Recipient Final Bilirubin (mg/dl)",
    tipss = "Recipient History of TIPSS",
    donor_org_shared = "Donor Organ Shared",
    exception = "Recipient Ever Approved for a MELD/PELD Exception",
    donor_cmv_status = "Donor CMV Status",
    prev_tx = "Previous Transplant",
    donor_arginine = "Donor on Arginine Vasopressin",
    donor_blood_infect = "Donor Blood Infection",
    donor_bun = "Donor BUN (mg/dl)",
    donor_cdc_high_risk = "Donor CDC High Risk",
    donor_clinical_infect = "Donor Clinical Infection",
    donor_cod = "Donor Cause of Death",
    donor_infect_other = "Donor Other Infection",
    donor_protein_urine = "Donor Protein in Urine",
    donor_diuretics = "Donor on Pre-Recovery Meds: Diuretics",
    donor_steroids = "Donor on Pre-Recovery Meds: Steroids",
    donor_t3 = "Donor on Pre-Recovery Meds: T3",
    donor_t4 = "Donor on Pre-Recovery Meds: T4",
    donor_infect_lung = "Donor Lung Infection",
    donor_bilirubin = "Donor Total Bilirubin (mg/dl)",
    donor_infect_urine = "Donor Urine Infection",
    donor_vasodil = "Donor on Vasodilators",
    waitlist_time_cat = "Time on the Waitlist (Months)",
    growth_failure = "Recipient Growth Failure",
    tx_year = "Transplant Year",
    donor_type = "Donor Type",
    gs_time = "Follow-Up Time"
  )
) %>%
  add_p(pvalue_fun = label_style_pvalue(digits = 2)) %>%
  modify_fmt_fun(all_categorical() ~ function(x) style_percent(digits = 2)) %>%
  bold_labels() %>%
  add_overall()

baseline_table

saveRDS(
  baseline_table,
  "L:/Projects/Investigator/Jonathan Merola/data/baseline_table"
)


#########################################
#create a donor risk index point system#
#########################################

pred_surv_cox_tbl <- survfit(cox_mle_model, newdata = data_test)
pred_surv_cox_tbl <- as.data.frame(pred_surv_cox_tbl$surv)


pred_surv_cox_tbl <- pred_surv_cox_tbl %>%
  slice(446) %>%
  pivot_longer(cols = everything(), names_to = c("time")) %>%
  rename(pred_surv_prob_1825 = value) %>%
  select(pred_surv_prob_1825) %>%
  mutate(patient_id = c(1:2104))


data_test$age_at_transplant <- as_factor(data_test$age_at_transplant)


#find the optimal factor to be divided by the model coefficients
risk_lookup_list <- lapply(seq(0.005, 0.6, by = 0.005), function(i) {
  dri <- coef_model_std_lm %>%
    mutate(dri = round(`coef(model_std_lm)` / i))

  dri <- dri %>%
    rownames_to_column() %>%
    rename(coef = `coef(model_std_lm)`, variable = rowname)

  dri <- dri %>%
    rowwise() %>%
    mutate(
      var = names(data_test)[
        sapply(names(data_test), function(v) str_starts(variable, v))
      ][1]
    ) %>%
    ungroup() %>%
    mutate(level = str_remove(variable, paste0("^", var))) %>%
    select(variable = var, dri, level)

  dri$level <- as_factor(dri$level)

  dri <- dri %>%
    filter(!is.na(variable))

  data_test_long <- data_test %>%
    mutate(patient_id = c(1:2104)) %>%
    select(-c(gs_time, gs_outcome, tx_year)) %>%
    pivot_longer(
      cols = c(
        age_at_transplant_cat,
        age_at_transplant,
        donor_age_at_tx_cat,
        blood_type,
        donor_blood_type,
        albumin,
        bmi_cat,
        donor_bmi_cat,
        cold_ischemic_time,
        donor_creat,
        creat,
        donor_ddavp,
        donor_death_circum,
        donor_death_mech,
        donor_hist_diab,
        diagnosis,
        dialysis_within_last_week,
        ebv_status,
        last_enceph,
        functional_status,
        sex,
        donor_sex,
        donor_hist_hyperten,
        donor_inotrop,
        inr,
        life_support,
        prev_malig,
        med_condition,
        ventilator,
        pvt,
        prev_ab_surg,
        donor_sgot,
        bilirubin,
        tipss,
        donor_org_shared,
        exception,
        donor_cmv_status,
        prev_tx,
        donor_arginine,
        donor_blood_infect,
        donor_bun,
        donor_cdc_high_risk,
        donor_clinical_infect,
        donor_cod,
        donor_infect_other,
        donor_protein_urine,
        donor_diuretics,
        donor_steroids,
        donor_t3,
        donor_t4,
        donor_infect_lung,
        donor_bilirubin,
        donor_infect_urine,
        donor_vasodil,
        waitlist_time_cat,
        growth_failure,
        donor_type
      ),
      names_to = "variable",
      values_to = "level"
    )

  dri_test <- left_join(data_test_long, dri, by = c("variable", "level"))

  dri_scores <- dri_test %>%
    group_by(patient_id) %>%
    summarise(donor_risk_score = sum(dri, na.rm = T))

  pred_test <- left_join(pred_surv_cox_tbl, dri_scores, by = "patient_id")

  risk_lookup <- pred_test %>%
    group_by(donor_risk_score) %>%
    summarise(
      mean_5yr_risk = mean(pred_surv_prob_1825),
      sd_5yr_risk = sd(pred_surv_prob_1825),
      num_patients = n()
    )

  risk_lookup$mean_5yr_risk <- round(risk_lookup$mean_5yr_risk, 3)

  risk_lookup
})


names(risk_lookup_list) <- sprintf("%.3f", seq(0.005, 0.6, by = 0.005))
#0.335 is the best

#calculate the DRI
dri <- coef_model_std_lm %>%
  mutate(dri = round(`coef(model_std_lm)` / 0.335))


dri <- dri %>%
  rownames_to_column() %>%
  rename(coef = `coef(model_std_lm)`, variable = rowname)

dri <- dri %>%
  rowwise() %>%
  mutate(
    var = names(data_test)[
      sapply(names(data_test), function(v) str_starts(variable, v))
    ][1]
  ) %>%
  ungroup() %>%
  mutate(level = str_remove(variable, paste0("^", var))) %>%
  select(variable = var, dri, level)


dri$level <- as_factor(dri$level)

dri <- dri %>%
  filter(!is.na(variable))

data_test$age_at_transplant <- as_factor(data_test$age_at_transplant)

data_test_long <- data_test %>%
  mutate(patient_id = c(1:2104)) %>%
  select(-c(gs_time, gs_outcome, tx_year)) %>%
  pivot_longer(
    cols = c(
      age_at_transplant_cat,
      age_at_transplant,
      donor_age_at_tx_cat,
      blood_type,
      donor_blood_type,
      albumin,
      bmi_cat,
      donor_bmi_cat,
      cold_ischemic_time,
      donor_creat,
      creat,
      donor_ddavp,
      donor_death_circum,
      donor_death_mech,
      donor_hist_diab,
      diagnosis,
      dialysis_within_last_week,
      ebv_status,
      last_enceph,
      functional_status,
      sex,
      donor_sex,
      donor_hist_hyperten,
      donor_inotrop,
      inr,
      life_support,
      prev_malig,
      med_condition,
      ventilator,
      pvt,
      prev_ab_surg,
      donor_sgot,
      bilirubin,
      tipss,
      donor_org_shared,
      exception,
      donor_cmv_status,
      prev_tx,
      donor_arginine,
      donor_blood_infect,
      donor_bun,
      donor_cdc_high_risk,
      donor_clinical_infect,
      donor_cod,
      donor_infect_other,
      donor_protein_urine,
      donor_diuretics,
      donor_steroids,
      donor_t3,
      donor_t4,
      donor_infect_lung,
      donor_bilirubin,
      donor_infect_urine,
      donor_vasodil,
      waitlist_time_cat,
      growth_failure,
      donor_type
    ),
    names_to = "variable",
    values_to = "level"
  )


dri_test <- left_join(data_test_long, dri, by = c("variable", "level"))

#create dri quintiles
dri_scores <- dri_test %>%
  group_by(patient_id) %>%
  summarise(donor_risk_score = sum(dri, na.rm = T)) %>%
  mutate(
    donor_risk_score_quintile = cut(
      donor_risk_score,
      breaks = quantile(
        donor_risk_score,
        probs = c(0, 0.25, 0.50, 0.75, 1),
        na.rm = T
      ),
      include.lowest = T,
      labels = c("PLDRI Q1", "PLDRI Q2", "PLDRI Q3", "PLDRI Q4")
    )
  )



summary(dri_scores$donor_risk_score)


##########################################################################
#create a graph of the 5-year graft survival by donor risk index quartile#
##########################################################################

data_test_dri <- data_test %>%
  mutate(patient_id = c(1:2104))

data_test_dri <- left_join(data_test_dri, dri_scores, by = "patient_id")

# data_test_dri$donor_risk_score <- as_factor(data_test_dri$donor_risk_score)
data_test_dri$donor_risk_score_quintile <- as_factor(
  data_test_dri$donor_risk_score_quintile
)


test <- data_test_dri %>%
  tbl_cross(col = donor_risk_score_quintile, row = age_at_transplant)
test


data_test_dri_1 <- data_test_dri %>%
  filter(donor_risk_score_quintile == "PLDRI Q1")


pred_surv_cox_1 <- survfit(cox_mle_model, newdata = data_test_dri_1)

test <- as.data.frame(pred_surv_cox_1)

test <- pred_surv_cox_1$surv %>%
  as.data.frame() %>%
  slice(446) %>%
  pivot_longer(cols = everything(), names_to = c("time"))
hist(test$value)


pred_surv_cox_1_est <- rowMeans(pred_surv_cox_1$surv)
pred_surv_cox_1_lower <- rowMeans(pred_surv_cox_1$lower)
pred_surv_cox_1_upper <- rowMeans(pred_surv_cox_1$upper)


pred_surv_cox_1 <- data.frame(
  time_pred = pred_surv_cox_1$time,
  estimate_pred = pred_surv_cox_1_est,
  lower_pred = pred_surv_cox_1_lower,
  upper_pred = pred_surv_cox_1_upper,
  donor_risk_score_quintile = "PLDRI Q1"
)


data_test_dri_2 <- data_test_dri %>%
  filter(donor_risk_score_quintile == "PLDRI Q2")

pred_surv_cox_2 <- survfit(cox_mle_model, newdata = data_test_dri_2)


pred_surv_cox_2_est <- rowMeans(pred_surv_cox_2$surv)
pred_surv_cox_2_lower <- rowMeans(pred_surv_cox_2$lower)
pred_surv_cox_2_upper <- rowMeans(pred_surv_cox_2$upper)


pred_surv_cox_2 <- data.frame(
  time_pred = pred_surv_cox_2$time,
  estimate_pred = pred_surv_cox_2_est,
  lower_pred = pred_surv_cox_2_lower,
  upper_pred = pred_surv_cox_2_upper,
  donor_risk_score_quintile = "PLDRI Q2"
)


data_test_dri_3 <- data_test_dri %>%
  filter(donor_risk_score_quintile == "PLDRI Q3")

pred_surv_cox_3 <- survfit(cox_mle_model, newdata = data_test_dri_3)

pred_surv_cox_3_est <- rowMeans(pred_surv_cox_3$surv)
pred_surv_cox_3_lower <- rowMeans(pred_surv_cox_3$lower)
pred_surv_cox_3_upper <- rowMeans(pred_surv_cox_3$upper)


pred_surv_cox_3 <- data.frame(
  time_pred = pred_surv_cox_3$time,
  estimate_pred = pred_surv_cox_3_est,
  lower_pred = pred_surv_cox_3_lower,
  upper_pred = pred_surv_cox_3_upper,
  donor_risk_score_quintile = "PLDRI Q3"
)


data_test_dri_4 <- data_test_dri %>%
  filter(donor_risk_score_quintile == "PLDRI Q4")

pred_surv_cox_4 <- survfit(cox_mle_model, newdata = data_test_dri_4)

pred_surv_cox_4_est <- rowMeans(pred_surv_cox_4$surv)
pred_surv_cox_4_lower <- rowMeans(pred_surv_cox_4$lower)
pred_surv_cox_4_upper <- rowMeans(pred_surv_cox_4$upper)


pred_surv_cox_4 <- data.frame(
  time_pred = pred_surv_cox_4$time,
  estimate_pred = pred_surv_cox_4_est,
  lower_pred = pred_surv_cox_4_lower,
  upper_pred = pred_surv_cox_4_upper,
  donor_risk_score_quintile = "PLDRI Q4"
)



pred_surv_cox_cc <- rbind(
  pred_surv_cox_1,
  pred_surv_cox_2,
  pred_surv_cox_3,
  pred_surv_cox_4
)

pred_surv_cox_cc$donor_risk_score_quintile <- as_factor(
  pred_surv_cox_cc$donor_risk_score_quintile
)

five_year_gs <- pred_surv_cox_cc %>%
  filter(time_pred == 1826)

pred_surv_cox_cc <- pred_surv_cox_cc %>%
  mutate(strata = donor_risk_score_quintile) %>%
  select(time_pred, estimate_pred, strata)

load("L:/Projects/Investigator/Jonathan Merola/data/liver_surv_df.RData")

dri_km_data <- rbind(pred_surv_cox_cc, liver_surv_df)

legend_df <- data.frame(
  strata = c("PLDRI Q1", "PLDRI Q2", "PLDRI Q3", "PLDRI Q4", "Living Donor"),
  surv = c("90.2%", "86.0%", "82.9%", "76.8%", "89.1%"),
  x = 1775,
  y = c(0.918, 0.86, 0.829, 0.768, 0.893)
)


dri_km <- ggplot(data = dri_km_data) +
  geom_step(
    aes(x = time_pred, y = estimate_pred, color = strata),
    show.legend = FALSE
  ) +
  labs(
    x = "Time (Days)",
    y = "Survival Probability",
    color = "",
    title = "5A - All Recipients"
  ) +
  scale_ggsurvfit(
    y_scales = list(breaks = seq(0.6, 1, by = 0.1), limits = c(0.6, 1)),
    x_scales = list(breaks = seq(0, 2750, by = 250), limits = c(0, 2750))
  ) +
  theme_ggsurvfit_default() +
  scale_fill_jama() +
  scale_color_jama() +
  geom_text(
    data = legend_df,
    aes(x = x + 75, y = y, label = paste0(strata, ": ", surv), color = strata),
    hjust = 0,
    size = 3.5,
    show.legend = FALSE
  )


dri_km

data_test_dri_count <- data_test_dri %>%
  group_by(donor_risk_score_quintile) %>%
  summarise(n = n())


#########################################################################################################
#create a graph of the 5-year graft survival by donor risk index quartile stratified by recipient age 13-18 and < 6
#########################################################################################################

#######
#13-18#
#######

data_test_dri_13_18 <- data_test_dri %>%
  filter(age_at_transplant %in% c("13", "14", "15", "16", "17"))
summary(data_test_dri_13_18$age_at_transplant)


#1
data_test_dri_1_13_18 <- data_test_dri_13_18 %>%
  filter(donor_risk_score_quintile == "PLDRI Q1")

pred_surv_cox_1_13_18 <- survfit(cox_mle_model, newdata = data_test_dri_1_13_18)

pred_surv_cox_1_13_18_time <- pred_surv_cox_1_13_18$time
pred_surv_cox_1_13_18_est <- rowMeans(pred_surv_cox_1_13_18$surv)


pred_surv_cox_1_13_18 <- data.frame(
  time_pred = pred_surv_cox_1_13_18_time,
  estimate_pred = pred_surv_cox_1_13_18_est,
  donor_risk_score_quintile = "PLDRI Q1"
)


#2
data_test_dri_2_13_18 <- data_test_dri_13_18 %>%
  filter(donor_risk_score_quintile == "PLDRI Q2")

pred_surv_cox_2_13_18 <- survfit(cox_mle_model, newdata = data_test_dri_2_13_18)

pred_surv_cox_2_13_18_time <- pred_surv_cox_2_13_18$time
pred_surv_cox_2_13_18_est <- rowMeans(pred_surv_cox_2_13_18$surv)


pred_surv_cox_2_13_18 <- data.frame(
  time_pred = pred_surv_cox_2_13_18_time,
  estimate_pred = pred_surv_cox_2_13_18_est,
  donor_risk_score_quintile = "PLDRI Q2"
)


#3
data_test_dri_3_13_18 <- data_test_dri_13_18 %>%
  filter(donor_risk_score_quintile == "PLDRI Q3")

pred_surv_cox_3_13_18 <- survfit(cox_mle_model, newdata = data_test_dri_3_13_18)

pred_surv_cox_3_13_18_time <- pred_surv_cox_3_13_18$time
pred_surv_cox_3_13_18_est <- rowMeans(pred_surv_cox_3_13_18$surv)


pred_surv_cox_3_13_18 <- data.frame(
  time_pred = pred_surv_cox_3_13_18_time,
  estimate_pred = pred_surv_cox_3_13_18_est,
  donor_risk_score_quintile = "PLDRI Q3"
)


#4
data_test_dri_4_13_18 <- data_test_dri_13_18 %>%
  filter(donor_risk_score_quintile == "PLDRI Q4")

pred_surv_cox_4_13_18 <- survfit(cox_mle_model, newdata = data_test_dri_4_13_18)

pred_surv_cox_4_13_18_time <- pred_surv_cox_4_13_18$time
pred_surv_cox_4_13_18_est <- rowMeans(pred_surv_cox_4_13_18$surv)


pred_surv_cox_4_13_18 <- data.frame(
  time_pred = pred_surv_cox_4_13_18_time,
  estimate_pred = pred_surv_cox_4_13_18_est,
  donor_risk_score_quintile = "PLDRI Q4"
)



pred_surv_cox_cc_13_18 <- rbind(
  pred_surv_cox_1_13_18,
  pred_surv_cox_2_13_18,
  pred_surv_cox_3_13_18,
  pred_surv_cox_4_13_18
)

pred_surv_cox_cc_13_18$donor_risk_score_quintile <- as_factor(
  pred_surv_cox_cc_13_18$donor_risk_score_quintile
)

five_year_gs_13_18 <- pred_surv_cox_cc_13_18 %>%
  filter(time_pred == 1826)

pred_surv_cox_cc_13_18 <- pred_surv_cox_cc_13_18 %>%
  mutate(strata = donor_risk_score_quintile) %>%
  select(time_pred, estimate_pred, strata)

load("L:/Projects/Investigator/Jonathan Merola/data/liver_surv_13_18_df.RData")

dri_km_13_18_data <- rbind(pred_surv_cox_cc_13_18, liver_surv_df_13_18)


legend_df_13_18 <- data.frame(
  strata = c("PLDRI Q1", "PLDRI Q2", "PLDRI Q3", "PLDRI Q4", "Living Donor"),
  surv = c("88.3%", "80.7%", "74.3%", "68.2%", "85.9%"),
  x = 1775,
  y = c(0.890, 0.807, 0.743, 0.682, 0.859)
)

dri_km_13_18 <- ggplot(data = dri_km_13_18_data) +
  geom_step(
    aes(x = time_pred, y = estimate_pred, color = strata),
    show.legend = FALSE
  ) +
  labs(
    x = "Time (Days)",
    y = "Survival Probability",
    color = "",
    title = "5D - Recipient Ages 13-17"
  ) +
  scale_ggsurvfit(
    y_scales = list(breaks = seq(0.6, 1, by = 0.1), limits = c(0.6, 1)),
    x_scales = list(breaks = seq(0, 2750, by = 250), limits = c(0, 2750))
  ) +
  theme_ggsurvfit_default() +
  scale_fill_jama() +
  scale_color_jama() +
  geom_text(
    data = legend_df_13_18,
    aes(x = x + 75, y = y, label = paste0(strata, ": ", surv), color = strata),
    hjust = 0,
    size = 3.5,
    show.legend = FALSE
  )

dri_km_13_18


######
# < 6#
######

data_test_dri_6 <- data_test_dri %>%
  filter(age_at_transplant %in% c("5", "4", "3", "2", "1", "0"))


#1
data_test_dri_1_6 <- data_test_dri_6 %>%
  filter(donor_risk_score_quintile == "PLDRI Q1")

pred_surv_cox_1_6 <- survfit(cox_mle_model, newdata = data_test_dri_1_6)

pred_surv_cox_1_6_time <- pred_surv_cox_1_6$time
pred_surv_cox_1_6_est <- rowMeans(pred_surv_cox_1_6$surv)


pred_surv_cox_1_6 <- data.frame(
  time_pred = pred_surv_cox_1_6_time,
  estimate_pred = pred_surv_cox_1_6_est,
  donor_risk_score_quintile = "PLDRI Q1"
)


#2
data_test_dri_2_6 <- data_test_dri_6 %>%
  filter(donor_risk_score_quintile == "PLDRI Q2")

pred_surv_cox_2_6 <- survfit(cox_mle_model, newdata = data_test_dri_2_6)

pred_surv_cox_2_6_time <- pred_surv_cox_2_6$time
pred_surv_cox_2_6_est <- rowMeans(pred_surv_cox_2_6$surv)


pred_surv_cox_2_6 <- data.frame(
  time_pred = pred_surv_cox_2_6_time,
  estimate_pred = pred_surv_cox_2_6_est,
  donor_risk_score_quintile = "PLDRI Q2"
)


#3
data_test_dri_3_6 <- data_test_dri_6 %>%
  filter(donor_risk_score_quintile == "PLDRI Q3")

pred_surv_cox_3_6 <- survfit(cox_mle_model, newdata = data_test_dri_3_6)

pred_surv_cox_3_6_time <- pred_surv_cox_3_6$time
pred_surv_cox_3_6_est <- rowMeans(pred_surv_cox_3_6$surv)


pred_surv_cox_3_6 <- data.frame(
  time_pred = pred_surv_cox_3_6_time,
  estimate_pred = pred_surv_cox_3_6_est,
  donor_risk_score_quintile = "PLDRI Q3"
)


#4
data_test_dri_4_6 <- data_test_dri_6 %>%
  filter(donor_risk_score_quintile == "PLDRI Q4")

pred_surv_cox_4_6 <- survfit(cox_mle_model, newdata = data_test_dri_4_6)

pred_surv_cox_4_6_time <- pred_surv_cox_4_6$time
pred_surv_cox_4_6_est <- rowMeans(pred_surv_cox_4_6$surv)


pred_surv_cox_4_6 <- data.frame(
  time_pred = pred_surv_cox_4_6_time,
  estimate_pred = pred_surv_cox_4_6_est,
  donor_risk_score_quintile = "PLDRI Q4"
)




pred_surv_cox_cc_6 <- rbind(
  pred_surv_cox_1_6,
  pred_surv_cox_2_6,
  pred_surv_cox_3_6,
  pred_surv_cox_4_6
)

pred_surv_cox_cc_6$donor_risk_score_quintile <- as_factor(
  pred_surv_cox_cc_6$donor_risk_score_quintile
)

five_year_gs_6 <- pred_surv_cox_cc_6 %>%
  filter(time_pred == 1826)

pred_surv_cox_cc_6 <- pred_surv_cox_cc_6 %>%
  mutate(strata = donor_risk_score_quintile) %>%
  select(time_pred, estimate_pred, strata)

load("L:/Projects/Investigator/Jonathan Merola/data/liver_surv_6_df.RData")

dri_km_6_data <- rbind(pred_surv_cox_cc_6, liver_surv_df_6)


legend_df_6 <- data.frame(
  strata = c("PLDRI Q1", "PLDRI Q2", "PLDRI Q3", "PLDRI Q4", "Living Donor"),
  surv = c("90.6%", "87.8%", "83.2%", "77%", "90.6%"),
  x = 1775,
  y = c(0.928, 0.878, 0.832, 0.77, 0.905)
)

dri_km_6 <- ggplot(data = dri_km_6_data) +
  geom_step(
    aes(x = time_pred, y = estimate_pred, color = strata),
    show.legend = FALSE
  ) +
  labs(
    x = "Time (Days)",
    y = "Survival Probability",
    color = "",
    title = "5B - Recipient Ages < 6"
  ) +
  scale_ggsurvfit(
    y_scales = list(breaks = seq(0.6, 1, by = 0.1), limits = c(0.6, 1)),
    x_scales = list(breaks = seq(0, 2750, by = 250), limits = c(0, 2750))
  ) +
  theme_ggsurvfit_default() +
  scale_fill_jama() +
  scale_color_jama() +
  geom_text(
    data = legend_df_6,
    aes(x = x + 75, y = y, label = paste0(strata, ": ", surv), color = strata),
    hjust = 0,
    size = 3.5,
    show.legend = FALSE
  )

dri_km_6


######
#6-12#
######

data_test_dri_6_12 <- data_test_dri %>%
  filter(age_at_transplant %in% c("6", "7", "8", "9", "10", "11", "12"))


#1
data_test_dri_1_6_12 <- data_test_dri_6_12 %>%
  filter(donor_risk_score_quintile == "PLDRI Q1")

pred_surv_cox_1_6_12 <- survfit(cox_mle_model, newdata = data_test_dri_1_6_12)

pred_surv_cox_1_6_12_time <- pred_surv_cox_1_6_12$time
pred_surv_cox_1_6_12_est <- rowMeans(pred_surv_cox_1_6_12$surv)


pred_surv_cox_1_6_12 <- data.frame(
  time_pred = pred_surv_cox_1_6_12_time,
  estimate_pred = pred_surv_cox_1_6_12_est,
  donor_risk_score_quintile = "PLDRI Q1"
)


#2
data_test_dri_2_6_12 <- data_test_dri_6_12 %>%
  filter(donor_risk_score_quintile == "PLDRI Q2")

pred_surv_cox_2_6_12 <- survfit(cox_mle_model, newdata = data_test_dri_2_6_12)

pred_surv_cox_2_6_12_time <- pred_surv_cox_2_6_12$time
pred_surv_cox_2_6_12_est <- rowMeans(pred_surv_cox_2_6_12$surv)


pred_surv_cox_2_6_12 <- data.frame(
  time_pred = pred_surv_cox_2_6_12_time,
  estimate_pred = pred_surv_cox_2_6_12_est,
  donor_risk_score_quintile = "PLDRI Q2"
)


#3
data_test_dri_3_6_12 <- data_test_dri_6_12 %>%
  filter(donor_risk_score_quintile == "PLDRI Q3")

pred_surv_cox_3_6_12 <- survfit(cox_mle_model, newdata = data_test_dri_3_6_12)

pred_surv_cox_3_6_12_time <- pred_surv_cox_3_6_12$time
pred_surv_cox_3_6_12_est <- rowMeans(pred_surv_cox_3_6_12$surv)


pred_surv_cox_3_6_12 <- data.frame(
  time_pred = pred_surv_cox_3_6_12_time,
  estimate_pred = pred_surv_cox_3_6_12_est,
  donor_risk_score_quintile = "PLDRI Q3"
)


#4
data_test_dri_4_6_12 <- data_test_dri_6_12 %>%
  filter(donor_risk_score_quintile == "PLDRI Q4")

pred_surv_cox_4_6_12 <- survfit(cox_mle_model, newdata = data_test_dri_4_6_12)

pred_surv_cox_4_6_12_time <- pred_surv_cox_4_6_12$time
pred_surv_cox_4_6_12_est <- rowMeans(pred_surv_cox_4_6_12$surv)


pred_surv_cox_4_6_12 <- data.frame(
  time_pred = pred_surv_cox_4_6_12_time,
  estimate_pred = pred_surv_cox_4_6_12_est,
  donor_risk_score_quintile = "PLDRI Q4"
)




pred_surv_cox_cc_6_12 <- rbind(
  pred_surv_cox_1_6_12,
  pred_surv_cox_2_6_12,
  pred_surv_cox_3_6_12,
  pred_surv_cox_4_6_12
)

pred_surv_cox_cc_6_12$donor_risk_score_quintile <- as_factor(
  pred_surv_cox_cc_6_12$donor_risk_score_quintile
)

five_year_gs_6_12 <- pred_surv_cox_cc_6_12 %>%
  filter(time_pred == 1826)

pred_surv_cox_cc_6_12 <- pred_surv_cox_cc_6_12 %>%
  mutate(strata = donor_risk_score_quintile) %>%
  select(time_pred, estimate_pred, strata)

load("L:/Projects/Investigator/Jonathan Merola/data/liver_surv_6_12_df.RData")

dri_km_6_12_data <- rbind(pred_surv_cox_cc_6_12, liver_surv_df_6_12)


legend_df_6_12 <- data.frame(
  strata = c("PLDRI Q1", "PLDRI Q2", "PLDRI Q3", "PLDRI Q4", "Living Donor"),
  surv = c("91.9%", "86.7%", "84.4%", "81%", "81.6%"),
  x = 1775,
  y = c(0.919, 0.867, 0.844, 0.785, 0.816)
)

dri_km_6_12 <- ggplot(data = dri_km_6_12_data) +
  geom_step(
    aes(x = time_pred, y = estimate_pred, color = strata),
    show.legend = FALSE
  ) +
  labs(
    x = "Time (Days)",
    y = "Survival Probability",
    color = "",
    title = "5C - Recipient Ages 6-12"
  ) +
  scale_ggsurvfit(
    y_scales = list(breaks = seq(0.6, 1, by = 0.1), limits = c(0.6, 1)),
    x_scales = list(breaks = seq(0, 2750, by = 250), limits = c(0, 2750))
  ) +
  theme_ggsurvfit_default() +
  scale_fill_jama() +
  scale_color_jama() +
  geom_text(
    data = legend_df_6_12,
    aes(x = x + 75, y = y, label = paste0(strata, ": ", surv), color = strata),
    hjust = 0,
    size = 3.5,
    show.legend = FALSE
  )

dri_km_6_12


dri_km_comp_plot <- (dri_km + dri_km_6 + dri_km_6_12 + dri_km_13_18) +
  plot_layout(ncol = 2, nrow = 2, axis_titles = "collect")

ggsave(
  filename = "L:/Projects/Investigator/Jonathan Merola/graphs and tables/dri_km_comp_plot.png",
  plot = dri_km_comp_plot,
  device = "png",
  dpi = "print",
  height = 6,
  width = 10
)




###########################################################################################
#create a plot of the observed survival probability and the predicted survival probability#
###########################################################################################

km_test <- survfit(Surv(gs_time, gs_outcome) ~ 1, data = data_test)

pred_surv_cox <- survfit(cox_mle_model, newdata = data_test)
pred_surv_cox_surv <- as.data.frame(pred_surv_cox$surv)

pred_surv_cox_surv <- rowMeans(pred_surv_cox_surv)

pred_surv_df <- data.frame(
  surv = pred_surv_cox_surv,
  time = pred_surv_cox$time,
  strata = "Mean Predicted"
)

obs_surv_df <- data.frame(
  surv = km_test$surv,
  time = km_test$time,
  strata = "Observed"
)

#calculate the Kolmogorov-Smirnov Test
ks.test(obs_surv_df$surv, pred_surv_df$surv)
mle_ks <- 0.007


pred_obs_df <- rbind(pred_surv_df, obs_surv_df)

pred_obs_1825 <- pred_obs_df %>%
  filter(time == 1826)


pred_obs_km <- ggplot(data = pred_obs_df) +
  geom_step(aes(x = time, y = surv, color = strata)) +
  labs(x = "Time (Days)", y = "Survival Probability", color = "") +
  scale_ggsurvfit(
    y_scales = list(breaks = seq(0.8, 1, by = 0.02), limits = c(0.8, 1)),
    x_scales = list(breaks = seq(0, 1850, by = 300), limits = c(0, 1850))
  ) +
  theme_ggsurvfit_default() +
  scale_fill_jama() +
  scale_color_jama() +
  geom_label(
    data = data.frame(
      x = 10,
      y = 0.82,
      label = "Observed Survival Probability at 5-Years: 85.4%\n Mean Predicted Survival Probability at 5-Years: 85.1%"
    ),
    aes(x = x, y = y, label = label),
    hjust = 0,
    size = 3.5
  )

ggsave(
  filename = "L:/Projects/Investigator/Jonathan Merola/graphs and tables/pred_obs_KM.png",
  plot = pred_obs_km,
  device = "png",
  dpi = "print",
  height = 6,
  width = 10
)


#survival at five years
mle_5yr_surv <- 0.851


###############################################################
#create a figure with the distribution of the donor risk score#
###############################################################

theoretical_dri_max <- dri %>%
  group_by(variable) %>%
  slice_max(dri) %>%
  filter(dri > 0) %>%
  distinct(variable, .keep_all = T)
sum(theoretical_dri_max$dri)


theoretical_dri_min <- dri %>%
  group_by(variable) %>%
  slice_min(dri) %>%
  filter(dri < 0) %>%
  distinct(variable, .keep_all = T)
sum(theoretical_dri_min$dri)

#10,-3

dri_scores <- dri_scores %>%
  arrange(donor_risk_score) %>%
  mutate(
    donor_risk_score_updated = case_when(
      donor_risk_score == -3 ~ 1,
      donor_risk_score == -2 ~ 2,
      donor_risk_score == -1 ~ 3,
      donor_risk_score == 0 ~ 4,
      donor_risk_score == 1 ~ 5,
      donor_risk_score == 2 ~ 6,
      donor_risk_score == 3 ~ 7,
      donor_risk_score == 4 ~ 8,
      donor_risk_score == 5 ~ 9,
      donor_risk_score == 6 ~ 10,
      donor_risk_score == 7 ~ 11,
      donor_risk_score == 8 ~ 12,
      donor_risk_score == 9 ~ 13,
      donor_risk_score == 10 ~ 14,
      TRUE ~ NA
    )
  ) %>%
  select(-c(donor_risk_score)) %>%
  rename(donor_risk_score = donor_risk_score_updated)

table(dri_scores$donor_risk_score)
hist(dri_scores$donor_risk_score)
summary(dri_scores$donor_risk_score)

#figure out the range for each dri quartile
tbl_summary(
  data = dri_scores,
  include = c(donor_risk_score, donor_risk_score_quintile),
  by = donor_risk_score_quintile,
  type = list(donor_risk_score = "continuous"),
  statistic = all_continuous() ~ "{median} ({min}, {max})"
)

test <- dri_scores %>%
  group_by(donor_risk_score) %>%
  summarise(n = n())


dri_hist <- ggplot(dri_scores) +
  aes(x = donor_risk_score) +
  geom_histogram(binwidth = 0.5) +
  labs(x = "Pediatric Liver Donor Risk Score", y = "Frequency") +
  theme_bw() +
  # xlim(0,100) +
  geom_label(
    data = data.frame(
      x = 10,
      y = 300,
      label = "PLDRI Q1: 2-5\n PLDRI Q2: 6\n PLDRI Q3: 7\n PLDRI Q4: 8-10"
    ),
    aes(x = x, y = y, label = label),
    hjust = 0,
    size = 4.5
  ) +
  scale_x_continuous(
    breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14),
    limits = c(1, 14)
  ) +
  scale_y_continuous(breaks = c(100, 200, 300, 400, 500), limits = c(0, 500))

dri_hist


ggsave(
  filename = "L:/Projects/Investigator/Jonathan Merola/graphs and tables/dri_hist.png",
  plot = dri_hist,
  device = "png",
  dpi = "print",
  height = 6,
  width = 9
)


###########################################################################################
#create a table with the predicted survival probability for each possible donor risk score#
###########################################################################################

risk_lookup_data <- risk_lookup_list$`0.335` %>%
  mutate(
    donor_risk_score = case_when(
      donor_risk_score == -3 ~ 1,
      donor_risk_score == -2 ~ 2,
      donor_risk_score == -1 ~ 3,
      donor_risk_score == 0 ~ 4,
      donor_risk_score == 1 ~ 5,
      donor_risk_score == 2 ~ 6,
      donor_risk_score == 3 ~ 7,
      donor_risk_score == 4 ~ 8,
      donor_risk_score == 5 ~ 9,
      donor_risk_score == 6 ~ 10,
      donor_risk_score == 7 ~ 11,
      donor_risk_score == 8 ~ 12,
      donor_risk_score == 9 ~ 13,
      donor_risk_score == 10 ~ 14,
      TRUE ~ NA
    )
  )

risk_lookup_data$mean_5yr_risk <- as.character(
  risk_lookup_data$mean_5yr_risk * 100
)
risk_lookup_data$sd_5yr_risk <- as.character(
  round(risk_lookup_data$sd_5yr_risk, 3) * 100
)

risk_lookup_data <- risk_lookup_data %>%
  mutate(
    mean_5yr_risk = paste0(mean_5yr_risk, "%"),
    sd_5yr_risk = paste0(sd_5yr_risk, "%"),
    five_year_risk = paste0((mean_5yr_risk), " (", (sd_5yr_risk), ")"),
    percent_patients = paste0(
      round(num_patients / (sum(num_patients)) * 100, 1),
      "%"
    )
  )


obs_surv_probs <- survfit(
  Surv(gs_time, gs_outcome) ~ donor_risk_score,
  data = data_test_dri
) %>%
  tbl_survfit(
    times = 1825,
    label_header = "**5-Year Survival (95% CI)**"
  )

obs_surv_probs_data <- obs_surv_probs$table_body

obs_surv_probs_data <- obs_surv_probs_data %>%
  filter(label != "donor_risk_score") %>%
  select(obs_surv_probs = stat_1)

risk_lookup_data <- cbind(risk_lookup_data, obs_surv_probs_data)


risk_lookup_tbl <- risk_lookup_data %>%
  select(donor_risk_score, five_year_risk, num_patients, percent_patients) %>%
  gt() %>%
  cols_label(
    donor_risk_score = "PLDRI Score",
    five_year_risk = "Average Predicted Survival Probability at 5-Years (%)",
    num_patients = "Number of Patients in the Test Dataset",
    percent_patients = "% of Patient in the Test Dataset"
  ) %>%
  tab_footnote(
    footnote = c(
      "The 5-Year Graft Loss Survival Probability for Living Donors was 89.1%"
    ),
    locations = cells_column_labels(columns = c(five_year_risk))
  ) %>%
  tab_footnote(
    footnote = c("Mean (SD)"),
    locations = cells_column_labels(columns = c(five_year_risk))
  )


risk_lookup_tbl

saveRDS(
  risk_lookup_tbl,
  "L:/Projects/Investigator/Jonathan Merola/data/risk_lookup_tbl"
)


###########################################################################################
#create a table of the model variables included, coefficients, and donor risk index points#
###########################################################################################

var_coef_dri_data <- dri

include_vars <- dri %>%
  group_by(variable) %>%
  summarise(keep = any(dri != 0)) %>%
  filter(keep == TRUE)

var_coef_dri_data <- var_coef_dri_data %>%
  filter(variable %in% include_vars$variable)

var_coef_dri_data <- var_coef_dri_data %>%
  select(variable, level, dri)


var_coef_dri_tbl <- var_coef_dri_data %>%
  mutate(variable = "") %>%
  gt() %>%
  cols_label(
    variable = "Donor Factor",
    dri = "Points",
    level = "Level"
  ) %>%
  tab_row_group(
    label = "Donor Age",
    rows = 1:14
  ) %>%
  tab_row_group(
    label = "Donor Terminal Creatinine (mg/dl)",
    rows = 15:17
  ) %>%
  tab_row_group(
    label = "Donor Death Circumstance",
    rows = 18:21
  ) %>%
  tab_row_group(
    label = "Donor Death Mechanism",
    rows = 22:30
  ) %>%
  tab_row_group(
    label = "Donor History of Diabetes",
    rows = 31
  ) %>%
  tab_row_group(
    label = "Donor Total Bilirubin (mg/dl)",
    rows = 32:33
  ) %>%
  tab_row_group(
    label = "Donor Type",
    rows = 34:35
  ) %>%
  tab_footnote(
    footnote = c(
      "Donor factors removed due to lack of contribution to the DRI: Cold-Ischemic Time, BMI, arginine vasopressin, CDC high risk, organ share type, cause of death, protein in urine, and urine infection."
    ),
    locations = cells_column_labels(columns = c(variable))
  )

var_coef_dri_tbl

saveRDS(
  var_coef_dri_tbl,
  "L:/Projects/Investigator/Jonathan Merola/data/var_coef_dri_tbl"
)

##########################################################################################
#create a figure for the model variables included, coefficients, and donor risk index points
##########################################################################################

var_coef_dri_fig_data <- var_coef_dri_data %>%
  mutate(
    variable = case_when(
      variable == "donor_age_at_tx_cat" ~ "Donor Age at Transplant:",
      variable == "donor_creat" ~ "Donor Terminal Creatinine (mg/dl):",
      variable == "donor_death_circum" ~ "Donor Death Circumstance:",
      variable == "donor_death_mech" ~ "Donor Death Mechanism:",
      variable == "donor_hist_diab" ~ "Donor History of Diabetes:",
      variable == "donor_bilirubin" ~ "Donor Total Bilirubin (mg/dl):",
      variable == "donor_type" ~ "Donor Type:",
      TRUE ~ NA
    ),
    label = paste(variable, "", level),
    label = fct_reorder(label, dri)
  ) %>%
  filter(
    dri != 0
  )

var_coef_dri_fig_data <- var_coef_dri_fig_data %>%
  mutate(label= case_when(label == "Donor Age at Transplant:  >=50" ~ "Donor Age at Transplant: ≥50",
                           label == "Donor Total Bilirubin (mg/dl):  >=2" ~ "Donor Total Bilirubin (mg/dl): ≥2",
                           TRUE ~ label))


dri_factors <- ggplot(
  var_coef_dri_fig_data,
  aes(x = dri, y = fct_reorder(label, dri), fill = dri > 0)
) +
  geom_col(width = 0.7) +
  geom_vline(xintercept = 0, linewidth = 0.8, color = "black") +
  facet_grid(variable ~ ., scales = "free_y", space = "free_y") +
  scale_x_continuous(
    breaks = c(-1, 0, 1, 2, 3, 4, 5),
    limits = c(-1, 5)
  ) +
  scale_fill_manual(
    values = c("TRUE" = "#D73027", "FALSE" = "#4575B4"),
    guide = "none"
  ) +
labs(
  x = "PLDRI Points",
  y = NULL
) +
  theme_bw() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    axis.text.y = element_text(size = 13L, face = "bold"),
    strip.text = element_blank()
  )

dri_factors

  ggsave(
  filename = "L:/Projects/Investigator/Jonathan Merola/graphs and tables/dri_factors.png",
  plot = dri_factors,
  device = "png",
  dpi = "print",
  height = 6,
  width = 9
)

#######################################################
#create a table of all of the final model coefficients#
#######################################################

final_model_coefficients <- as.data.frame(as.matrix(coef_model_std_lm)) %>%
  rownames_to_column() %>%
  rename(coef = `coef(model_std_lm)`, variable = rowname) %>%
  filter(variable != "(Intercept)") %>%
  rowwise() %>%
  mutate(
    var = names(data_test)[
      sapply(names(data_test), function(v) str_starts(variable, v))
    ][1]
  ) %>%
  ungroup() %>%
  mutate(Level = str_remove(variable, paste0("^", var)))

final_model_coefficients$coef <- round(final_model_coefficients$coef, 3)




final_model_coef_tbl <- final_model_coefficients %>%
  select(-c(variable)) %>%
  select(var, Level, coef) %>%
  mutate(var = "") %>%
  gt() %>%
  cols_label(
    var = "Donor Factor",
    Level = "Level",
    coef = "Marginalized 5-Year Coefficient",
  ) %>%
  tab_row_group(
    label = "Donor Age at Transplant",
    rows = 1:14
  ) %>%
  tab_row_group(
    label = "Donor BMI (kg/m2)",
    rows = 15:16
  ) %>%
  tab_row_group(
    label = "Cold-Ischemic Time (Hours)",
    rows = 17:19
  ) %>%
  tab_row_group(
    label = "Donor Terminal Creatinine (mg/dl)",
    rows = 20:22
  ) %>%
  tab_row_group(
    label = "Donor Death Circumstance",
    rows = 23:26
  ) %>%
  tab_row_group(
    label = "Donor Death Mechanism",
    rows = 27:35
  ) %>%
  tab_row_group(
    label = "Donor History of Diabetes",
    rows = 36
  ) %>%
  tab_row_group(
    label = "Donor Organ Shared",
    rows = 37
  ) %>%
  tab_row_group(
    label = "Donor on Arginine Vasopressin",
    rows = 38
  ) %>%
  tab_row_group(
    label = "Donor CDC High Risk",
    rows = 39
  ) %>%
  tab_row_group(
    label = "Donor Protein in Urine",
    rows = 40
  ) %>%
  tab_row_group(
    label = "Donor Total Bilirubin (mg/dl)",
    rows = 41:42
  ) %>%
  tab_row_group(
    label = "Donor Urine Infection",
    rows = 43
  ) %>%
  tab_row_group(
    label = "Donor Type",
    rows = 44:45
  ) %>%
  tab_footnote(
    footnote = c(
      "In the 2010-2026 cohort, the most predictive factors of a decreased survival probability were still donor age at transplant, split donor type (0.924), partial donor type (0.623), and history of diabetes (0.378)."
    ),
    locations = cells_column_labels(columns = c(coef))
  ) %>%
  tab_footnote(
    footnote = c(
      "In the 2010-2026 cohort, the most predictive factors of an increased survival probability were other cause of death (-0.635), asphyxiation death mechanism (-0.532), seizure death mechanism (-0.531), and CDC high risk (-0.437). "
    ),
    locations = cells_column_labels(columns = c(coef))
  )

final_model_coef_tbl

saveRDS(
  final_model_coef_tbl,
  "L:/Projects/Investigator/Jonathan Merola/data/final_model_coef_tbl"
)


load("L:/Projects/Investigator/Jonathan Merola/data/final_coefs_1yr.RData")

#############################################
#construct a table for overall model metrics#
#############################################

model_metrics <- data.frame(
  model = c(
    "LASSO",
    "Ridge",
    "Elastic Net",
    "Adaptive LASSO",
    "Cox-Proportional Hazards"
  ),
  c_index = c(
    model_lasso_min_cindex,
    model_ridge_min_cindex,
    model_EL_min_cindex,
    model_adapt_min_cindex,
    mle_cindex
  ),
  brier_score = c(
    model_lasso_min_brier,
    model_ridge_min_brier,
    model_EL_min_brier,
    model_adapt_min_brier,
    mle_brier_score
  ),
  dr = c(dr_lasso, dr_ridge, dr_EL, dr_adapt, NA),
  auc = c(lasso_auc, ridge_auc, EL_auc, adapt_auc, mle_auc),
  ks = c(lasso_ks, ridge_ks, EL_ks, adapt_ks, mle_ks),
  surv_5yr = c(
    lasso_5yr_surv,
    ridge_5yr_surv,
    EL_5yr_surv,
    adapt_5yr_surv,
    mle_5yr_surv
  )
)

model_metrics$c_index <- round(model_metrics$c_index, 3)
model_metrics$brier_score <- round(model_metrics$brier_score, 3)
model_metrics$dr <- round(model_metrics$dr, 3)

model_metrics_tbl <- model_metrics %>%
  gt() %>%
  cols_label(
    model = "Model",
    c_index = "Harrel's C-Index",
    brier_score = "Brier Score",
    dr = "Deviance Ratio",
    auc = "AUROC",
    ks = "Kolmogorov-Smirnov Test",
    surv_5yr = "Predicted 5-Year Survival (%)"
  ) %>%
  tab_footnote(
    footnote = c("Brier Score and AUROC were calculated at the five years"),
    locations = cells_column_labels(columns = c(brier_score, auc))
  ) %>%
  tab_footnote(
    footnote = c("AUROC = Area under the receiver operating curve."),
    locations = cells_column_labels(columns = c(auc))
  ) %>%
  tab_footnote(
    footnote = c("The observed 5-Year Survival was 85.4%"),
    locations = cells_column_labels(columns = c(surv_5yr))
  ) %>%
  tab_footnote(
    footnote = c(
      "For the Cox-Proportional Hazards model, the predicted 5-year survival is the mean average."
    ),
    locations = cells_column_labels(columns = c(surv_5yr))
  ) %>%
  tab_footnote(
    footnote = c("Values shown are p-values."),
    locations = cells_column_labels(columns = c(ks))
  ) %>%
  fmt_percent(
    columns = surv_5yr,
    decimals = 1
  )

model_metrics_tbl

saveRDS(
  model_metrics_tbl,
  "L:/Projects/Investigator/Jonathan Merola/data/model_metrics_tbl"
)


######################################
#create a figure 1 exclusion criteria#
######################################

exclusion_text <- c(
  "Total excluded = 1,429 (14.5%)\n 1,140 Multi-Organ Transplants\n 289 No follow-up"
)

training_text <- c("5,447 (86.3%) Censored")
test_text <- c("1,810 (86%) Censored")


flow_chart <- as_fc(
  N = 9845,
  label = "Pediatric Liver Recipients from 2004-06-30 to 2026-01-01",
  text_fs = 11,
  text_pattern = "{label}\n N = {N}"
) %>%
  fc_filter(
    N = 8416,
    label = "N for Analysis",
    show_exc = T,
    perc_total = F,
    text_pattern = "{label}\n N = {n}",
    text_fs_exc = 10,
    text_fs = 11
  ) %>%
  fc_split(
    N = c(6312, 2104),
    label = c("Training Data", "Testing Data"),
    text_fs = 11,
    text_pattern = "{label}\n N = {n}"
  ) %>%
  fc_filter(
    N = c(865,294),
    label = c("Graft Losses within 5-Years"),
    text_fs = 11,
    text_pattern = "{label}\n N = {n}",
    show_exc = T,
    text_fs_exc = 10
  ) %>%
  fc_modify(
    ~ . |>
      mutate(
        text = case_when(id == 3 ~ exclusion_text,
                         id == 6 ~ "Graft Losses within 5-Years\n N = 865 (13.7%)",
                         id == 7 ~ training_text,
                         id == 8 ~ "Graft Losses within 5-Years\n N = 294 (14%)",
                         id == 9 ~ test_text,
                         TRUE ~ text),
        just = ifelse(id == 3, "left", just),
        x = case_when(
          id == 3 ~ 0.80,
          id %in% c(1, 2) ~ 0.42,
          id == 4 ~ 0.25,
          id == 5 ~ 0.60,
          id == 6 ~ 0.25,
          id == 7 ~ 0.45,
          id == 8 ~ 0.60,
          id == 9 ~ 0.78,
          TRUE ~ x
        ),
        y = case_when(
          id == 1 ~ 0.90,
          id == 2 ~ 0.67,
          id == 3 ~ 0.76,
          id == 2 ~ 0.40,
          id == 4 ~ 0.35,
          id == 5 ~ 0.35,
          id == 6 ~ 0.12,
          id == 7 ~ 0.24,
          id == 8 ~ 0.12,
          id == 9 ~ 0.24,
          TRUE ~ y
        )
      )
  )

flow_chart <- flow_chart %>%
  fc_draw(big.mark = ",") %>%
  fc_export(
    "L:/Projects/Investigator/Jonathan Merola/graphs and tables/fig_1.png",
    width = 8,
    height = 6,
    res = 300,
    units = "in"
  )


##################################
#calculate Somers D (w/censoring)#
##################################

data_test_transformed_dri <- data_test_dri %>%
  mutate(
    donor_risk_score = case_when(
      donor_risk_score == -3 ~ 1,
      donor_risk_score == -2 ~ 2,
      donor_risk_score == -1 ~ 3,
      donor_risk_score == 0 ~ 4,
      donor_risk_score == 1 ~ 5,
      donor_risk_score == 2 ~ 6,
      donor_risk_score == 3 ~ 7,
      donor_risk_score == 4 ~ 8,
      donor_risk_score == 5 ~ 9,
      donor_risk_score == 6 ~ 10,
      donor_risk_score == 7 ~ 11,
      donor_risk_score == 8 ~ 12,
      donor_risk_score == 9 ~ 13,
      donor_risk_score == 10 ~ 14,
      TRUE ~ NA
    )
  )


dri_model <- coxph(
  Surv(gs_time, gs_outcome) ~ donor_risk_score,
  data = data_test_transformed_dri
)
summary(dri_model)

dri_model_cindex <- concordance(dri_model)$concordance

somers_d <- (2 * dri_model_cindex) - 1


dri_model_tbl <- dri_model %>%
  tbl_regression(
    exponentiate = T,
    label = list(donor_risk_score = "PLDRI Score"),
    pvalue_fun = label_style_pvalue(digits = 2),
    estimate_fun = label_style_sigfig(digits = 3)
  ) %>%
  as_gt() %>%
  tab_footnote(
    footnote = "Harrell's C-Index: 0.570, Somers D: 0.139",
    locations = cells_column_labels(columns = c(estimate))
  )

dri_model_tbl


saveRDS(
  dri_model_tbl,
  "L:/Projects/Investigator/Jonathan Merola/data/dri_model_tbl"
)


