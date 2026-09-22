#title: Calculate and Plot Psuedo R-Squared values for Donor Variables
#date: 8/18/2026
#author: Ryan LaVanchy
#purpose: Calculate and Plot Psuedo R-Squared values for Donor Variables for the pediatric liver donor risk index project

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
library(performance)

#read in the training data
load("L:/Projects/Investigator/Jonathan Merola/data/data_train.RData")


#read in the main cox model
cox_mle_model <- readRDS("L:/Projects/Investigator/Jonathan Merola/data/cox_mle_model")

# cox_mle_model <- coxph(Surv(gs_time, gs_outcome) ~ age_at_transplant_cat + donor_age_at_tx_cat + blood_type + albumin + bmi_cat + donor_bmi_cat + cold_ischemic_time + donor_creat + creat + donor_death_circum + donor_death_mech + donor_hist_diab + diagnosis + dialysis_within_last_week + ebv_status + functional_status + inr + prev_malig + med_condition + pvt + bilirubin + tipss + donor_org_shared + prev_tx + donor_arginine + donor_cdc_high_risk + donor_protein_urine + donor_bilirubin + donor_infect_urine + growth_failure + tx_year + donor_type, data = data_train)

#calculate overall Kent & O-Quigley psuedo R-Squared
r2_full <- royston(cox_mle_model, adjust = T)[4]


#calculate Kent & O'Quigley psuedo R-Squared for each variable in the model as well as all donor and all recipient factors

#########################
#donor age at transplant#
#########################


no_donor_age_at_tx_model <- coxph(Surv(gs_time, gs_outcome) ~ age_at_transplant_cat + blood_type + albumin + bmi_cat + donor_bmi_cat + cold_ischemic_time + donor_creat + creat + donor_death_circum + donor_death_mech + donor_hist_diab + diagnosis + dialysis_within_last_week + ebv_status + functional_status + inr + prev_malig + med_condition + pvt + bilirubin + tipss + donor_org_shared + prev_tx + donor_arginine + donor_cdc_high_risk + donor_protein_urine + donor_bilirubin + donor_infect_urine + growth_failure + tx_year + donor_type, data = data_train)

r2_no_donor_age_at_tx <- royston(no_donor_age_at_tx_model, adjust = T)[4]

r2_no_donor_age_at_tx <- r2_full - r2_no_donor_age_at_tx




###################
#age at transplant#
###################


no_age_at_tx_model <- coxph(Surv(gs_time, gs_outcome) ~ donor_age_at_tx_cat + blood_type + albumin + bmi_cat + donor_bmi_cat + cold_ischemic_time + donor_creat + creat + donor_death_circum + donor_death_mech + donor_hist_diab + diagnosis + dialysis_within_last_week + ebv_status + functional_status + inr + prev_malig + med_condition + pvt + bilirubin + tipss + donor_org_shared + prev_tx + donor_arginine + donor_cdc_high_risk + donor_protein_urine + donor_bilirubin + donor_infect_urine + growth_failure + tx_year + donor_type, data = data_train)

r2_no_age_at_tx <- royston(no_age_at_tx_model, adjust = T)[4]

r2_no_age_at_tx <- r2_full - r2_no_age_at_tx


############
#blood type#
############

no_blood_type_model <- coxph(Surv(gs_time, gs_outcome) ~ age_at_transplant_cat + donor_age_at_tx_cat + albumin + bmi_cat + donor_bmi_cat + cold_ischemic_time + donor_creat + creat + donor_death_circum + donor_death_mech + donor_hist_diab + diagnosis + dialysis_within_last_week + ebv_status + functional_status + inr + prev_malig + med_condition + pvt + bilirubin + tipss + donor_org_shared + prev_tx + donor_arginine + donor_cdc_high_risk + donor_protein_urine + donor_bilirubin + donor_infect_urine + growth_failure + tx_year + donor_type, data = data_train)

r2_no_blood_type <- royston(no_blood_type_model, adjust = T)[4]

r2_no_blood_type <- r2_full - r2_no_blood_type


#########
#albumin#
#########

no_albumin_model <- coxph(Surv(gs_time, gs_outcome) ~ age_at_transplant_cat + donor_age_at_tx_cat + blood_type + bmi_cat + donor_bmi_cat + cold_ischemic_time + donor_creat + creat + donor_death_circum + donor_death_mech + donor_hist_diab + diagnosis + dialysis_within_last_week + ebv_status + functional_status + inr + prev_malig + med_condition + pvt + bilirubin + tipss + donor_org_shared + prev_tx + donor_arginine + donor_cdc_high_risk + donor_protein_urine + donor_bilirubin + donor_infect_urine + growth_failure + tx_year + donor_type, data = data_train)

r2_no_albumin <- royston(no_albumin_model, adjust = T)[4]

r2_no_albumin <- r2_full - r2_no_albumin


#####
#bmi#
#####

no_bmi_model <- coxph(Surv(gs_time, gs_outcome) ~ age_at_transplant_cat + donor_age_at_tx_cat + blood_type + albumin + donor_bmi_cat + cold_ischemic_time + donor_creat + creat + donor_death_circum + donor_death_mech + donor_hist_diab + diagnosis + dialysis_within_last_week + ebv_status + functional_status + inr + prev_malig + med_condition + pvt + bilirubin + tipss + donor_org_shared + prev_tx + donor_arginine + donor_cdc_high_risk + donor_protein_urine + donor_bilirubin + donor_infect_urine + growth_failure + tx_year + donor_type, data = data_train)

r2_no_bmi <- royston(no_bmi_model, adjust = T)[4]

r2_no_bmi <- r2_full - r2_no_bmi


###########
#donor bmi#
###########

no_donor_bmi_model <- coxph(Surv(gs_time, gs_outcome) ~ age_at_transplant_cat + donor_age_at_tx_cat + blood_type + albumin + bmi_cat + cold_ischemic_time + donor_creat + creat + donor_death_circum + donor_death_mech + donor_hist_diab + diagnosis + dialysis_within_last_week + ebv_status + functional_status + inr + prev_malig + med_condition + pvt + bilirubin + tipss + donor_org_shared + prev_tx + donor_arginine + donor_cdc_high_risk + donor_protein_urine + donor_bilirubin + donor_infect_urine + growth_failure + tx_year + donor_type, data = data_train)

r2_no_donor_bmi <- royston(no_donor_bmi_model, adjust = T)[4]

r2_no_donor_bmi <- r2_full - r2_no_donor_bmi


####################
#cold-ischemic time#
####################

no_cit_model <- coxph(Surv(gs_time, gs_outcome) ~ age_at_transplant_cat + donor_age_at_tx_cat + blood_type + albumin + bmi_cat + donor_bmi_cat + donor_creat + creat + donor_death_circum + donor_death_mech + donor_hist_diab + diagnosis + dialysis_within_last_week + ebv_status + functional_status + inr + prev_malig + med_condition + pvt + bilirubin + tipss + donor_org_shared + prev_tx + donor_arginine + donor_cdc_high_risk + donor_protein_urine + donor_bilirubin + donor_infect_urine + growth_failure + tx_year + donor_type, data = data_train)

r2_no_cit <- royston(no_cit_model, adjust = T)[4]

r2_no_cit <- r2_full - r2_no_cit


##################
#donor creatinine#
##################

no_donor_creat_model <- coxph(Surv(gs_time, gs_outcome) ~ age_at_transplant_cat + donor_age_at_tx_cat + blood_type + albumin + bmi_cat + donor_bmi_cat + cold_ischemic_time + creat + donor_death_circum + donor_death_mech + donor_hist_diab + diagnosis + dialysis_within_last_week + ebv_status + functional_status + inr + prev_malig + med_condition + pvt + bilirubin + tipss + donor_org_shared + prev_tx + donor_arginine + donor_cdc_high_risk + donor_protein_urine + donor_bilirubin + donor_infect_urine + growth_failure + tx_year + donor_type, data = data_train)

r2_no_donor_creat <- royston(no_donor_creat_model, adjust = T)[4]

r2_no_donor_creat <- r2_full - r2_no_donor_creat

######################
#recipient creatinine#
######################

no_creat_model <- coxph(Surv(gs_time, gs_outcome) ~ age_at_transplant_cat + donor_age_at_tx_cat + blood_type + albumin + bmi_cat + donor_bmi_cat + cold_ischemic_time + donor_creat + donor_death_circum + donor_death_mech + donor_hist_diab + diagnosis + dialysis_within_last_week + ebv_status + functional_status + inr + prev_malig + med_condition + pvt + bilirubin + tipss + donor_org_shared + prev_tx + donor_arginine + donor_cdc_high_risk + donor_protein_urine + donor_bilirubin + donor_infect_urine + growth_failure + tx_year + donor_type, data = data_train)

r2_no_creat <- royston(no_creat_model, adjust = T)[4]

r2_no_creat <- r2_full - r2_no_creat



##########################
#donor death circumstance#
##########################

no_donor_death_cir_model <- coxph(Surv(gs_time, gs_outcome) ~ age_at_transplant_cat + donor_age_at_tx_cat + blood_type + albumin + bmi_cat + donor_bmi_cat + cold_ischemic_time + donor_creat + creat + donor_death_mech + donor_hist_diab + diagnosis + dialysis_within_last_week + ebv_status + functional_status + inr + prev_malig + med_condition + pvt + bilirubin + tipss + donor_org_shared + prev_tx + donor_arginine + donor_cdc_high_risk + donor_protein_urine + donor_bilirubin + donor_infect_urine + growth_failure + tx_year + donor_type, data = data_train)

r2_no_donor_death_cir <- royston(no_donor_death_cir_model, adjust = T)[4]

r2_no_donor_death_cir <- r2_full - r2_no_donor_death_cir


#######################
#donor death mechanism#
#######################

no_donor_death_mech_model <- coxph(Surv(gs_time, gs_outcome) ~ age_at_transplant_cat + donor_age_at_tx_cat + blood_type + albumin + bmi_cat + donor_bmi_cat + cold_ischemic_time + donor_creat + creat + donor_death_circum + donor_hist_diab + diagnosis + dialysis_within_last_week + ebv_status + functional_status + inr + prev_malig + med_condition + pvt + bilirubin + tipss + donor_org_shared + prev_tx + donor_arginine + donor_cdc_high_risk + donor_protein_urine + donor_bilirubin + donor_infect_urine + growth_failure + tx_year + donor_type, data = data_train)

r2_no_donor_death_mech <- royston(no_donor_death_mech_model, adjust = T)[4]

r2_no_donor_death_mech <- r2_full - r2_no_donor_death_mech


###########################
#donor history of diabetes#
###########################

no_donor_hist_diab_model <- coxph(Surv(gs_time, gs_outcome) ~ age_at_transplant_cat + donor_age_at_tx_cat + blood_type + albumin + bmi_cat + donor_bmi_cat + cold_ischemic_time + donor_creat + creat + donor_death_circum + donor_death_mech + diagnosis + dialysis_within_last_week + ebv_status + functional_status + inr + prev_malig + med_condition + pvt + bilirubin + tipss + donor_org_shared + prev_tx + donor_arginine + donor_cdc_high_risk + donor_protein_urine + donor_bilirubin + donor_infect_urine + growth_failure + tx_year + donor_type, data = data_train)

r2_no_donor_hist_diab <- royston(no_donor_hist_diab_model, adjust = T)[4]

r2_no_donor_hist_diab <- r2_full - r2_no_donor_hist_diab


###########
#diagnosis#
###########

no_diagnosis_model <- coxph(Surv(gs_time, gs_outcome) ~ age_at_transplant_cat + donor_age_at_tx_cat + blood_type + albumin + bmi_cat + donor_bmi_cat + cold_ischemic_time + donor_creat + creat + donor_death_circum + donor_death_mech + donor_hist_diab + dialysis_within_last_week + ebv_status + functional_status + inr + prev_malig + med_condition + pvt + bilirubin + tipss + donor_org_shared + prev_tx + donor_arginine + donor_cdc_high_risk + donor_protein_urine + donor_bilirubin + donor_infect_urine + growth_failure + tx_year + donor_type, data = data_train)

r2_no_diagnosis <- royston(no_diagnosis_model, adjust = T)[4]

r2_no_diagnosis <- r2_full - r2_no_diagnosis


###########################
#dialysis within last week#
###########################

no_dial_last_week_model <- coxph(Surv(gs_time, gs_outcome) ~ age_at_transplant_cat + donor_age_at_tx_cat + blood_type + albumin + bmi_cat + donor_bmi_cat + cold_ischemic_time + donor_creat + creat + donor_death_circum + donor_death_mech + donor_hist_diab + diagnosis + ebv_status + functional_status + inr + prev_malig + med_condition + pvt + bilirubin + tipss + donor_org_shared + prev_tx + donor_arginine + donor_cdc_high_risk + donor_protein_urine + donor_bilirubin + donor_infect_urine + growth_failure + tx_year + donor_type, data = data_train)

r2_no_dial_last_week <- royston(no_dial_last_week_model, adjust = T)[4]

r2_no_dial_last_week <- r2_full - r2_no_dial_last_week


############
#EBV status#
############

no_ebv_status_model <- coxph(Surv(gs_time, gs_outcome) ~ age_at_transplant_cat + donor_age_at_tx_cat + blood_type + albumin + bmi_cat + donor_bmi_cat + cold_ischemic_time + donor_creat + creat + donor_death_circum + donor_death_mech + donor_hist_diab + diagnosis + dialysis_within_last_week + functional_status + inr + prev_malig + med_condition + pvt + bilirubin + tipss + donor_org_shared + prev_tx + donor_arginine + donor_cdc_high_risk + donor_protein_urine + donor_bilirubin + donor_infect_urine + growth_failure + tx_year + donor_type, data = data_train)

r2_no_ebv_status <- royston(no_ebv_status_model, adjust = T)[4]

r2_no_ebv_status <- r2_full - r2_no_ebv_status


###################
#functional status#
###################

no_func_status_model <- coxph(Surv(gs_time, gs_outcome) ~ age_at_transplant_cat + donor_age_at_tx_cat + blood_type + albumin + bmi_cat + donor_bmi_cat + cold_ischemic_time + donor_creat + creat + donor_death_circum + donor_death_mech + donor_hist_diab + diagnosis + dialysis_within_last_week + ebv_status + inr + prev_malig + med_condition + pvt + bilirubin + tipss + donor_org_shared + prev_tx + donor_arginine + donor_cdc_high_risk + donor_protein_urine + donor_bilirubin + donor_infect_urine + growth_failure + tx_year + donor_type, data = data_train)

r2_no_func_status <- royston(no_func_status_model, adjust = T)[4]

r2_no_func_status <- r2_full - r2_no_func_status



#####
#inr#
#####

no_inr_model <- coxph(Surv(gs_time, gs_outcome) ~ age_at_transplant_cat + donor_age_at_tx_cat + blood_type + albumin + bmi_cat + donor_bmi_cat + cold_ischemic_time + donor_creat + creat + donor_death_circum + donor_death_mech + donor_hist_diab + diagnosis + dialysis_within_last_week + ebv_status + functional_status + prev_malig + med_condition + pvt + bilirubin + tipss + donor_org_shared + prev_tx + donor_arginine + donor_cdc_high_risk + donor_protein_urine + donor_bilirubin + donor_infect_urine + growth_failure + tx_year + donor_type, data = data_train)

r2_no_inr <- royston(no_inr_model, adjust = T)[4]

r2_no_inr <- r2_full - r2_no_inr



#####################
#previous malignancy#
#####################

no_prev_malig_model <- coxph(Surv(gs_time, gs_outcome) ~ age_at_transplant_cat + donor_age_at_tx_cat + blood_type + albumin + bmi_cat + donor_bmi_cat + cold_ischemic_time + donor_creat + creat + donor_death_circum + donor_death_mech + donor_hist_diab + diagnosis + dialysis_within_last_week + ebv_status + functional_status + inr + med_condition + pvt + bilirubin + tipss + donor_org_shared + prev_tx + donor_arginine + donor_cdc_high_risk + donor_protein_urine + donor_bilirubin + donor_infect_urine + growth_failure + tx_year + donor_type, data = data_train)

r2_no_prev_malig <- royston(no_prev_malig_model, adjust = T)[4]

r2_no_prev_malig <- r2_full - r2_no_prev_malig


###################
#medical condition#
###################

no_med_cond_model <- coxph(Surv(gs_time, gs_outcome) ~ age_at_transplant_cat + donor_age_at_tx_cat + blood_type + albumin + bmi_cat + donor_bmi_cat + cold_ischemic_time + donor_creat + creat + donor_death_circum + donor_death_mech + donor_hist_diab + diagnosis + dialysis_within_last_week + ebv_status + functional_status + inr + prev_malig + pvt + bilirubin + tipss + donor_org_shared + prev_tx + donor_arginine + donor_cdc_high_risk + donor_protein_urine + donor_bilirubin + donor_infect_urine + growth_failure + tx_year + donor_type, data = data_train)

r2_no_med_cond <- royston(no_med_cond_model, adjust = T)[4]

r2_no_med_cond <- r2_full - r2_no_med_cond


#####
#pvt#
#####

no_pvt_model <- coxph(Surv(gs_time, gs_outcome) ~ age_at_transplant_cat + donor_age_at_tx_cat + blood_type + albumin + bmi_cat + donor_bmi_cat + cold_ischemic_time + donor_creat + creat + donor_death_circum + donor_death_mech + donor_hist_diab + diagnosis + dialysis_within_last_week + ebv_status + functional_status + inr + prev_malig + med_condition + bilirubin + tipss + donor_org_shared + prev_tx + donor_arginine + donor_cdc_high_risk + donor_protein_urine + donor_bilirubin + donor_infect_urine + growth_failure + tx_year + donor_type, data = data_train)

r2_no_pvt <- royston(no_pvt_model, adjust = T)[4]

r2_no_pvt <- r2_full - r2_no_pvt


###########
#bilirubin#
###########

no_bili_model <- coxph(Surv(gs_time, gs_outcome) ~ age_at_transplant_cat + donor_age_at_tx_cat + blood_type + albumin + bmi_cat + donor_bmi_cat + cold_ischemic_time + donor_creat + creat + donor_death_circum + donor_death_mech + donor_hist_diab + diagnosis + dialysis_within_last_week + ebv_status + functional_status + inr + prev_malig + med_condition + pvt + tipss + donor_org_shared + prev_tx + donor_arginine + donor_cdc_high_risk + donor_protein_urine + donor_bilirubin + donor_infect_urine + growth_failure + tx_year + donor_type, data = data_train)

r2_no_bili <- royston(no_bili_model, adjust = T)[4]

r2_no_bili <- r2_full - r2_no_bili


#######
#tipss#
#######

no_tipss_model <- coxph(Surv(gs_time, gs_outcome) ~ age_at_transplant_cat + donor_age_at_tx_cat + blood_type + albumin + bmi_cat + donor_bmi_cat + cold_ischemic_time + donor_creat + creat + donor_death_circum + donor_death_mech + donor_hist_diab + diagnosis + dialysis_within_last_week + ebv_status + functional_status + inr + prev_malig + med_condition + pvt + bilirubin + donor_org_shared + prev_tx + donor_arginine + donor_cdc_high_risk + donor_protein_urine + donor_bilirubin + donor_infect_urine + growth_failure + tx_year + donor_type, data = data_train)

r2_no_tipss <- royston(no_tipss_model, adjust = T)[4]

r2_no_tipss <- r2_full - r2_no_tipss


####################
#donor organ shared#
####################

no_donor_org_shared_model <- coxph(Surv(gs_time, gs_outcome) ~ age_at_transplant_cat + donor_age_at_tx_cat + blood_type + albumin + bmi_cat + donor_bmi_cat + cold_ischemic_time + donor_creat + creat + donor_death_circum + donor_death_mech + donor_hist_diab + diagnosis + dialysis_within_last_week + ebv_status + functional_status + inr + prev_malig + med_condition + pvt + bilirubin + tipss + prev_tx + donor_arginine + donor_cdc_high_risk + donor_protein_urine + donor_bilirubin + donor_infect_urine + growth_failure + tx_year + donor_type, data = data_train)

r2_no_donor_org_shared <- royston(no_donor_org_shared_model, adjust = T)[4]

r2_no_donor_org_shared <- r2_full - r2_no_donor_org_shared


#####################
#previous transplant#
#####################

no_prev_tx_model <- coxph(Surv(gs_time, gs_outcome) ~ age_at_transplant_cat + donor_age_at_tx_cat + blood_type + albumin + bmi_cat + donor_bmi_cat + cold_ischemic_time + donor_creat + creat + donor_death_circum + donor_death_mech + donor_hist_diab + diagnosis + dialysis_within_last_week + ebv_status + functional_status + inr + prev_malig + med_condition + pvt + bilirubin + tipss + donor_org_shared + donor_arginine + donor_cdc_high_risk + donor_protein_urine + donor_bilirubin + donor_infect_urine + growth_failure + tx_year + donor_type, data = data_train)

r2_no_prev_tx <- royston(no_prev_tx_model, adjust = T)[4]

r2_no_prev_tx <- r2_full - r2_no_prev_tx


################
#donor arginine#
################

no_donor_arginine_model <- coxph(Surv(gs_time, gs_outcome) ~ age_at_transplant_cat + donor_age_at_tx_cat + blood_type + albumin + bmi_cat + donor_bmi_cat + cold_ischemic_time + donor_creat + creat + donor_death_circum + donor_death_mech + donor_hist_diab + diagnosis + dialysis_within_last_week + ebv_status + functional_status + inr + prev_malig + med_condition + pvt + bilirubin + tipss + donor_org_shared + prev_tx + donor_cdc_high_risk + donor_protein_urine + donor_bilirubin + donor_infect_urine + growth_failure + tx_year + donor_type, data = data_train)

r2_no_donor_arginine <- royston(no_donor_arginine_model, adjust = T)[4]

r2_no_donor_arginine <- r2_full - r2_no_donor_arginine

#####################
#donor CDC high risk#
#####################

no_donor_cdc_hr_model <- coxph(Surv(gs_time, gs_outcome) ~ age_at_transplant_cat + donor_age_at_tx_cat + blood_type + albumin + bmi_cat + donor_bmi_cat + cold_ischemic_time + donor_creat + creat + donor_death_circum + donor_death_mech + donor_hist_diab + diagnosis + dialysis_within_last_week + ebv_status + functional_status + inr + prev_malig + med_condition + pvt + bilirubin + tipss + donor_org_shared + prev_tx + donor_arginine + donor_protein_urine + donor_bilirubin + donor_infect_urine + growth_failure + tx_year + donor_type, data = data_train)

r2_no_donor_cdc_hr <- royston(no_donor_cdc_hr_model, adjust = T)[4]

r2_no_donor_cdc_hr <- r2_full - r2_no_donor_cdc_hr


#####################
#donor protein urine#
#####################

no_donor_protein_urine_model <- coxph(Surv(gs_time, gs_outcome) ~ age_at_transplant_cat + donor_age_at_tx_cat + blood_type + albumin + bmi_cat + donor_bmi_cat + cold_ischemic_time + donor_creat + creat + donor_death_circum + donor_death_mech + donor_hist_diab + diagnosis + dialysis_within_last_week + ebv_status + functional_status + inr + prev_malig + med_condition + pvt + bilirubin + tipss + donor_org_shared + prev_tx + donor_arginine + donor_cdc_high_risk + donor_bilirubin + donor_infect_urine + growth_failure + tx_year + donor_type, data = data_train)

r2_no_donor_protein_urine <- royston(no_donor_protein_urine_model, adjust = T)[4]

r2_no_donor_protein_urine <- r2_full - r2_no_donor_protein_urine


#################
#donor bilirubin#
#################

no_donor_bili_model <- coxph(Surv(gs_time, gs_outcome) ~ age_at_transplant_cat + donor_age_at_tx_cat + blood_type + albumin + bmi_cat + donor_bmi_cat + cold_ischemic_time + donor_creat + creat + donor_death_circum + donor_death_mech + donor_hist_diab + diagnosis + dialysis_within_last_week + ebv_status + functional_status + inr + prev_malig + med_condition + pvt + bilirubin + tipss + donor_org_shared + prev_tx + donor_arginine + donor_cdc_high_risk + donor_protein_urine + donor_infect_urine + growth_failure + tx_year + donor_type, data = data_train)

r2_no_donor_bili <- royston(no_donor_bili_model, adjust = T)[4]

r2_no_donor_bili <- r2_full - r2_no_donor_bili


#######################
#donor urine infection#
#######################

no_donor_urine_infect_model <- coxph(Surv(gs_time, gs_outcome) ~ age_at_transplant_cat + donor_age_at_tx_cat + blood_type + albumin + bmi_cat + donor_bmi_cat + cold_ischemic_time + donor_creat + creat + donor_death_circum + donor_death_mech + donor_hist_diab + diagnosis + dialysis_within_last_week + ebv_status + functional_status + inr + prev_malig + med_condition + pvt + bilirubin + tipss + donor_org_shared + prev_tx + donor_arginine + donor_cdc_high_risk + donor_protein_urine + donor_bilirubin + growth_failure + tx_year + donor_type, data = data_train)

r2_no_donor_urine_infect <- royston(no_donor_urine_infect_model, adjust = T)[4]

r2_no_donor_urine_infect <- r2_full - r2_no_donor_urine_infect


################
#growth failure#
################

no_growth_failure_model <- coxph(Surv(gs_time, gs_outcome) ~ age_at_transplant_cat + donor_age_at_tx_cat + blood_type + albumin + bmi_cat + donor_bmi_cat + cold_ischemic_time + donor_creat + creat + donor_death_circum + donor_death_mech + donor_hist_diab + diagnosis + dialysis_within_last_week + ebv_status + functional_status + inr + prev_malig + med_condition + pvt + bilirubin + tipss + donor_org_shared + prev_tx + donor_arginine + donor_cdc_high_risk + donor_protein_urine + donor_bilirubin + donor_infect_urine + tx_year + donor_type, data = data_train)

r2_no_growth_failure <- royston(no_growth_failure_model, adjust = T)[4]

r2_no_growth_failure <- r2_full - r2_no_growth_failure


#########
#tx_year#
#########

no_tx_year_model <- coxph(Surv(gs_time, gs_outcome) ~ age_at_transplant_cat + donor_age_at_tx_cat + blood_type + albumin + bmi_cat + donor_bmi_cat + cold_ischemic_time + donor_creat + creat + donor_death_circum + donor_death_mech + donor_hist_diab + diagnosis + dialysis_within_last_week + ebv_status + functional_status + inr + prev_malig + med_condition + pvt + bilirubin + tipss + donor_org_shared + prev_tx + donor_arginine + donor_cdc_high_risk + donor_protein_urine + donor_bilirubin + donor_infect_urine + growth_failure + donor_type, data = data_train)

r2_no_tx_year <- royston(no_tx_year_model, adjust = T)[4]

r2_no_tx_year <- r2_full - r2_no_tx_year


############
#donor type#
############

no_donor_type_model <- coxph(Surv(gs_time, gs_outcome) ~ age_at_transplant_cat + donor_age_at_tx_cat + blood_type + albumin + bmi_cat + donor_bmi_cat + cold_ischemic_time + donor_creat + creat + donor_death_circum + donor_death_mech + donor_hist_diab + diagnosis + dialysis_within_last_week + ebv_status + functional_status + inr + prev_malig + med_condition + pvt + bilirubin + tipss + donor_org_shared + prev_tx + donor_arginine + donor_cdc_high_risk + donor_protein_urine + donor_bilirubin + donor_infect_urine + growth_failure + tx_year, data = data_train)

r2_no_donor_type <- royston(no_donor_type_model, adjust = T)[4]

r2_no_donor_type <- r2_full - r2_no_donor_type


###################
#recipient factors#
###################

no_recip_factor_model <- coxph(Surv(gs_time, gs_outcome) ~ donor_age_at_tx_cat + donor_bmi_cat + cold_ischemic_time + donor_creat + donor_death_circum + donor_death_mech + donor_hist_diab + donor_org_shared + donor_arginine + donor_cdc_high_risk + donor_protein_urine + donor_bilirubin + donor_infect_urine + donor_type, data = data_train)

r2_no_recip_factor <- royston(no_recip_factor_model, adjust = T)[4]

r2_no_recip_factor <- r2_full - r2_no_recip_factor


###############
#donor factors#
###############

no_donor_factors_model <- coxph(Surv(gs_time, gs_outcome) ~ age_at_transplant_cat + blood_type + albumin + bmi_cat + creat + diagnosis + dialysis_within_last_week + ebv_status + functional_status + inr + prev_malig + med_condition + pvt + bilirubin + tipss + prev_tx + growth_failure + tx_year, data = data_train)

r2_no_donor_factors <- royston(no_donor_factors_model, adjust = T)[4]

r2_no_donor_factors <- r2_full - r2_no_donor_factors






#create a data frame of the partial R2 and their factors
r2_data <- data.frame(
  partial_r2 = c(r2_no_recip_factor,r2_no_donor_factors, r2_no_donor_age_at_tx, r2_no_donor_bmi, r2_no_cit,
             r2_no_donor_creat, r2_no_donor_death_cir, r2_no_donor_death_mech,
             r2_no_donor_hist_diab, r2_no_donor_org_shared,
             r2_no_donor_arginine,
             r2_no_donor_cdc_hr,
             r2_no_donor_protein_urine, r2_no_donor_bili, r2_no_donor_urine_infect,
             r2_no_donor_type),
  factor = c("Recipient Factors", "Donor Factors","Donor Age at Transplant", "Donor BMI", "Cold-Ischemic Time", "Donor Terminal Creatinine", "Donor Death Circumstance", "Donor Mechanism of Death", "Donor History of Diabetes", "Donor Organ Shared", "Donor on Arginine Vasopressin", "Donor CDC High Risk", "Donor Protein in Urine", "Donor Total Bilirubin", "Donor Urine Infection", "Donor Type"),
  strata = c("Total Recipient Factors","Total Donor Factors","Individual Donor Factors","Individual Donor Factors","Individual Donor Factors","Individual Donor Factors","Individual Donor Factors","Individual Donor Factors","Individual Donor Factors","Individual Donor Factors","Individual Donor Factors","Individual Donor Factors","Individual Donor Factors","Individual Donor Factors","Individual Donor Factors","Individual Donor Factors")
)



#round the partial R2 variable
r2_data <- r2_data %>%
  mutate(partial_r2 = round(partial_r2, 5))


#create a pie chart of recipient and donor factors
pie_data <- r2_data %>%
  filter(factor %in% c("Recipient Factors","Donor Factors"))

pie_data$partial_r2 <- round(pie_data$partial_r2,2)

pos <- cumsum(pie_data$partial_r2) - 0.57 * pie_data$partial_r2


pie_chart <- ggplot(pie_data, aes(x = "", y = partial_r2, fill = factor)) +
  geom_col(width = 1, color = "white", show.legend = FALSE) +
  coord_polar(theta = "y") +
  geom_text(aes(y = pos, label = c("Recipient Factors:\n 0.14","Donor Factors:\n 0.06"), colour = "white"), fontface = "bold", size = 3.6) +
  labs(fill = "",
       title = "A. Donor and Recipient Factors") +
  scale_fill_jama() +
  scale_colour_identity() +
  theme_void() 

pie_chart

#create data and a graph for individual donor factors.

r2_individual_donor_factors_data <- r2_data %>%
  filter(strata == "Individual Donor Factors")



r2_individual_donor_factors_plot <- r2_individual_donor_factors_data %>%
  mutate(factor = reorder(factor, partial_r2, decreasing = T)) %>%
  ggplot() +
  aes(x = factor, y = partial_r2) +
  geom_col(fill = "#374e55") +
  labs(x = "Factor", y = "Partial Kent & O'Quigley Pseudo R-Squared", title = "B. Individual Donor Factors") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60L, hjust = 1L)) +
  scale_y_continuous(limits = c(0, 0.05), breaks = seq(0, 0.05, 0.01)) +
  geom_label(data = data.frame(x = "Donor Type", y = 0.045, label = "Overall Kent & O'Quigley Pseudo R-Squared: 0.251"), aes(x = x, y = y, label = label), hjust = 0, size = 4.5) 

r2_individual_donor_factors_plot


comp_r2_plot <- (pie_chart + r2_individual_donor_factors_plot) +
  plot_layout(ncol = 2, nrow = 1)

comp_r2_plot

#save the composite plot
ggsave(filename = "L:/Projects/Investigator/Jonathan Merola/graphs and tables/comp_r2_plot.png", plot = comp_r2_plot, device = "png", dpi = "print", height = 6, width = 10)
