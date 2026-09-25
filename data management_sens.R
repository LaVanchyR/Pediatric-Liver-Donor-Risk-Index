#title: Data Management for Pediatric Donor Risk Index
#author: Ryan LaVanchy
#date: 8/31/2026
#purpose: Manage data for the pediatric donor risk index data. Only include patients from 2010 onward.


#Libraries
library(tidyverse)
library(performance)
library(naniar)
library(mice)
library(glmnet)
library(haven)
library(gtsummary)


#import data
donor_deceased <- read_sas("S:/CORP/data/srtr/srtr_202603/pubsaf2603/donor_deceased.sas7bdat")
tx_li <- read_sas("S:/CORP/data/srtr/srtr_202603/pubsaf2603/tx_li.sas7bdat")
mp_except <- read_sas("S:/CORP/data/srtr/srtr_202603/pubsaf2603/mpexcept.sas7bdat")


#select relevant variables
liver_cohort <- tx_li %>%
  select(PERS_ID, PX_ID, DONOR_ID, DON_TY, REC_TX_DT, REC_HGT_CM, REC_WGT_KG, CAN_LISTING_DT, PERS_RETX, TFL_DEATH_DT, PERS_OPTN_DEATH_DT, PERS_SSA_DEATH_DT, TFL_ENDTXFU, REC_AGE_AT_TX, REC_AGE_IN_MONTHS_AT_TX, DON_AGE, REC_TX_ORG_TY, ORG_TY, CAN_ABO, DON_ABO, CAN_LAST_ALBUMIN, REC_ARTIFICIAL_LI, REC_ASCITES, REC_BMI, REC_COLD_ISCH_TM, DON_CREAT, CAN_LAST_SERUM_CREAT, DON_DDAVP, DON_DEATH_CIRCUM, DON_DEATH_MECH, CAN_DIAB, DON_HIST_DIAB, REC_DGN, CAN_LAST_DIAL_PRIOR_WEEK, REC_EBV_STAT, CAN_LAST_ENCEPH, CAN_ETHNICITY_SRTR, DON_ETHNICITY_SRTR, DON_RACE_SRTR, CAN_RACE_SRTR, REC_FUNCTN_STAT, CAN_GENDER, DON_HIST_CANCER, DON_HIST_CIGARETTE_GT20_PKYR, DON_HIST_COCAINE, DON_HIST_DIAB, DON_HIST_HYPERTEN, DON_HIST_OTHER_DRUG, CAN_LAST_INR, DON_INSULIN, REC_LIFE_SUPPORT, CAN_MALIG, REC_MED_COND, CAN_LAST_SRTR_LAB_MELD, DON_NON_HR_BEAT, DON_GENDER, REC_VENTILATOR, REC_PORTAL_VEIN_THROMB, REC_PORTAL_VEIN, REC_PREV_ABDOM_SURG, DON_SGOT, DON_SGPT_PREOP, DON_SGOT_PREOP, CAN_LAST_BILI, REC_TIPSS, REC_PROCEDURE_TY_LI, DON_ORG_SHARED, CANHX_MPXCPT_HCC_APPROVE_IND, DON_ANTI_CMV, CAN_PREV_TX, REC_ACADEMIC_LEVEL, REC_ACADEMIC_PROGRESS, REC_LI_TY, DON_AGE_IN_MONTHS, DON_WGT_KG, DON_HGT_CM)


liver_donor <- donor_deceased %>%
  select(DONOR_ID, DON_HEAVY_ALCOHOL, DON_ARGININE, DON_INFECT_BLOOD_CONFIRM, DON_INFECT_BLOOD, DON_BUN, DON_CARDIAC_ARREST_AFTER_DEATH, DON_MEET_CDC_HIGH_RISK, DON_CLINICAL_INFECT, DON_CAD_DON_COD, DON_EBV_IGG, DON_EBV_IGM, DON_HIST_INSULIN_DEPND, DON_HIST_PREV_MI, DON_INOTROP_SUPPORT, DON_INTRACRANIAL_CANCER, DON_EJECT_FRACT, DON_INFECT_OTHER, DON_INFECT_BLOOD_CONFIRM, DON_PROTEIN_URINE, DON_HYPERTEN_DIURETICS, DON_PRERECOV_DIURETICS, DON_PRERECOV_STEROIDS, DON_PRERECOV_T3, DON_PRERECOV_T4, DON_INFECT_LU, DON_INFECT_LU_CONFIRM, DON_RECOV_OUT_US, DON_RESUSCIT_DURATION, DON_TATTOOS, DON_TOT_BILI, DON_INFECT_URINE_CONFIRM, DON_INFECT_URINE, DON_VASODIL)


#examine missing data
liver_cohort[liver_cohort == ""] <- NA
liver_donor[liver_donor == ""] <- NA


liver_cohort_missing <- miss_var_summary(liver_cohort)
#these variables are missing more than 30% of data: DON_SGOT_PREOP, DON_SGPT_PREOP, REC_PORTAL_VEIN_THROMB, CAN_DIAB, REC_ASCITES, REC_PROCEDURE_TY_LI, DON_INSULIN, REC_ACADEMIC_LEVEL, REC_ACADEMIC_PROGESS


liver_donor_missing <- miss_var_summary(liver_donor)
#these variables are missing more than 30% of data: DON_INFECT_BLOOD_CONFIRM, DON_RESUSCIT_DURATION, DON_INFECT_URINE_CONFIRM, DON_HIST_INSULIN_DEPND, DON_INFECT_LU_CONFIRM, DON_HYPERTEN_DIURETICS, DON_EJECT_FRACT,DON_EBV_IGM, DON_EBV_IGG


#remove the variables with greater than 30% of missingness
liver_cohort <- liver_cohort %>%
  select(-c(DON_SGOT_PREOP, DON_SGPT_PREOP, REC_PORTAL_VEIN_THROMB, CAN_DIAB, REC_ASCITES, REC_PROCEDURE_TY_LI, DON_INSULIN, REC_ACADEMIC_LEVEL, REC_ACADEMIC_PROGRESS))

liver_donor <- liver_donor %>%
  select(-c(DON_INFECT_BLOOD_CONFIRM, DON_RESUSCIT_DURATION, DON_INFECT_URINE_CONFIRM, DON_HIST_INSULIN_DEPND, DON_INFECT_LU_CONFIRM, DON_HYPERTEN_DIURETICS, DON_EJECT_FRACT,DON_EBV_IGM, DON_EBV_IGG))


summary(liver_cohort$REC_TX_DT)
#filter dataset to be between 1/1/2010 and 2026-01-01, and be a deceased donor recipient, have an age less than 18 years, and signal their initial transplant. 
liver_cohort <- liver_cohort %>%
  filter(REC_TX_DT >= "2010-01-01" & REC_TX_DT <= "2026-01-01" & DON_TY == "C" & REC_AGE_AT_TX < 18) %>%
  group_by(PERS_ID) %>%
  slice_min(REC_TX_DT) %>%
  ungroup()


n_distinct(liver_cohort$PERS_ID)


#remove multiorgan recipients
multi_org_liver <- liver_cohort %>%
  filter(REC_TX_ORG_TY != ORG_TY) %>%
  select(PERS_ID, REC_TX_ORG_TY, ORG_TY)

liver_cohort <- liver_cohort %>%
  filter(!PERS_ID %in% multi_org_liver$PERS_ID)
#703





#join the liver cohort and the liver donor datasets
liver_cohort <- left_join(liver_cohort, liver_donor, by = "DONOR_ID")



#examine missing data over time
liver_cohort_missing_over_time <- liver_cohort %>%
  mutate(year = substr(REC_TX_DT,1,4)) %>%
  group_by(year) %>%
  miss_var_summary() %>%
  arrange(variable, year)
# remove donor history of cocaine, donor history of other drug, and donor tattoos. They have a lot or are completely missing for the last few years of data 

liver_cohort <- liver_cohort %>%
  select(-c(DON_HIST_COCAINE, DON_HIST_OTHER_DRUG, DON_TATTOOS))


#code categorical and continuous variables for analysis
liver_cohort_clean <- liver_cohort %>%
  mutate(death_date = case_when(!is.na(TFL_DEATH_DT) ~ TFL_DEATH_DT,
                                !is.na(PERS_OPTN_DEATH_DT) & is.na(TFL_DEATH_DT) ~ PERS_OPTN_DEATH_DT,
                                !is.na(PERS_SSA_DEATH_DT) & is.na(TFL_DEATH_DT) & is.na(PERS_OPTN_DEATH_DT) ~ PERS_SSA_DEATH_DT),
         graft_failure_date = pmin(death_date, PERS_RETX, na.rm = T),
         five_year_censor = as_date(REC_TX_DT %m+% months(60)),
         last_fu_date = pmin(graft_failure_date, TFL_ENDTXFU, five_year_censor, na.rm = T),
         gs_time = round(as.numeric(difftime(last_fu_date, REC_TX_DT, units = "days")),0),
         gs_outcome = case_when(last_fu_date == graft_failure_date ~ 1,
                                last_fu_date == TFL_ENDTXFU | last_fu_date == five_year_censor ~ 0,
                                TRUE ~ NA),
         age_at_transplant = REC_AGE_AT_TX,
         age_at_transplant_cat = factor(case_when(REC_AGE_AT_TX %in% c(0:2) ~ "0-2",
                                                  REC_AGE_AT_TX %in% c(3:5) ~ "3-5",
                                                  REC_AGE_AT_TX %in% c(6:8) ~ "6-8",
                                                  REC_AGE_AT_TX %in% c(9:11) ~ "9-11",
                                                  REC_AGE_AT_TX %in% c(12:14) ~ "12-14",
                                                  REC_AGE_AT_TX %in% c(15:17) ~ "15-17",
                                                  TRUE ~ NA),
                                        levels = c("0-2","3-5","6-8","9-11","12-14","15-17")),
         donor_age_at_tx = DON_AGE,
         donor_age_at_tx_cat = factor(case_when(DON_AGE_IN_MONTHS >= 0 & DON_AGE_IN_MONTHS < 6 ~ "0-6mo",
                                                DON_AGE_IN_MONTHS >= 6 & DON_AGE_IN_MONTHS < 12 ~ "6-12mo",
                                                DON_AGE_IN_MONTHS >= 12 & DON_AGE_IN_MONTHS < 24 ~ "1-2",
                                                DON_AGE_IN_MONTHS >= 24 & DON_AGE_IN_MONTHS < 48 ~ "2-4",
                                                DON_AGE_IN_MONTHS >= 48 & DON_AGE_IN_MONTHS < 72 ~ "4-6",
                                                DON_AGE_IN_MONTHS >= 72 & DON_AGE_IN_MONTHS < 96 ~ "6-8",
                                                DON_AGE_IN_MONTHS >= 96 & DON_AGE_IN_MONTHS < 120 ~ "8-10",
                                                DON_AGE_IN_MONTHS >= 120 & DON_AGE_IN_MONTHS < 144 ~ "10-12",
                                                DON_AGE_IN_MONTHS >= 144 & DON_AGE_IN_MONTHS < 168 ~ "12-14",
                                                DON_AGE_IN_MONTHS >= 168 & DON_AGE_IN_MONTHS < 192 ~ "14-16",
                                                DON_AGE_IN_MONTHS >= 192 & DON_AGE_IN_MONTHS < 216 ~ "16-18",
                                                DON_AGE_IN_MONTHS >= 216 & DON_AGE_IN_MONTHS < 360 ~ "18-30",
                                                DON_AGE_IN_MONTHS >= 360 & DON_AGE_IN_MONTHS < 480 ~ "30-40",
                                                DON_AGE_IN_MONTHS >= 480 & DON_AGE_IN_MONTHS < 600 ~ "40-50",
                                                DON_AGE_IN_MONTHS >= 600 ~ ">=50", TRUE ~ NA),
                                      levels = c("10-12","0-6mo","6-12mo","1-2","2-4","4-6","6-8","8-10",
                                                 "12-14","14-16","16-18","18-30","30-40","40-50",">=50")),
         blood_type = factor(case_when(CAN_ABO %in% c("A","A1","A2") ~ "A",
                                                       CAN_ABO %in% c("AB","A1B","A2B") ~ "AB",
                                                       CAN_ABO == "B" ~ "B",
                                                       CAN_ABO == "O" ~ "O",
                                                       TRUE ~ NA),
                             levels = c("A","AB","B","O")),
         donor_blood_type = factor(case_when(DON_ABO %in% c("A","A1","A2") ~ "A",
                                             DON_ABO %in% c("AB","A1B","A2B") ~ "AB",
                                             DON_ABO == "B" ~ "B",
                                             DON_ABO == "O" ~ "O",
                                             TRUE ~ NA),
                                   levels = c("A","AB","B","O")),
         albumin = factor(case_when(CAN_LAST_ALBUMIN < 2.5 ~ "<2.5",
                                    CAN_LAST_ALBUMIN >= 2.5 & CAN_LAST_ALBUMIN <= 3.5 ~ "2.5-3.5",
                                    CAN_LAST_ALBUMIN > 3.5 ~ ">3.5",
                                    TRUE ~ NA),
                          levels = c("<2.5", "2.5-3.5", ">3.5")),
         artificial_liver = factor(case_when(REC_ARTIFICIAL_LI == 1 ~ "Yes",
                                             REC_ARTIFICIAL_LI == 0 ~ "No",
                                             TRUE ~ NA),
                                   levels = c("No","Yes")),
         bmi_cat = factor(case_when(REC_BMI < 18.5 ~ "<18.5",
                                REC_BMI >= 18.5 & REC_BMI <= 25 ~ "18.5-25",
                                REC_BMI >= 25 ~ ">25",
                                REC_BMI %in% c(7621.9196,661.3042,235.2941,181.1111,115.6345,114.1080,
                                               108.2045,106.0281,102.1849,96.5706,78.9157,76.7849,
                                               66.9754) ~ NA,
                                TRUE ~ NA),
                      levels = c("<18.5","18.5-25",">25")),
         donor_bmi = ((DON_WGT_KG) / ((DON_HGT_CM / 100) * (DON_HGT_CM / 100))),
         donor_bmi_cat = factor(case_when(donor_bmi < 18.5 ~ "<18.5",
                                          donor_bmi >= 18.5 & donor_bmi <= 25 ~ "18.5-25",
                                          donor_bmi >= 25 & donor_bmi < 60 ~ ">25",
                                          TRUE ~ NA),
                                levels = c("<18.5","18.5-25",">25")),
         cold_ischemic_time = factor(case_when(REC_COLD_ISCH_TM >= 0 & REC_COLD_ISCH_TM < 4 ~ "0-4",
                                               REC_COLD_ISCH_TM >= 4 & REC_COLD_ISCH_TM < 6 ~ "4-6",
                                               REC_COLD_ISCH_TM >= 6 & REC_COLD_ISCH_TM < 8 ~ "6-8",
                                               REC_COLD_ISCH_TM >= 8 ~ ">=8",
                                               TRUE ~ NA),
                                     levels = c("0-4","4-6","6-8",">=8")),
         donor_creat = factor(case_when(DON_CREAT < 0.3 ~ "<0.3",
                                        DON_CREAT >= 0.3 & DON_CREAT < 0.7 ~ "0.3-0.7",
                                        DON_CREAT >= 0.7 & DON_CREAT < 1 ~ "0.7-1",
                                        DON_CREAT >= 1 ~ ">=1",
                                        TRUE ~ NA),
                              levels = c("<0.3","0.3-0.7","0.7-1",">=1")),
         creat = factor(case_when(CAN_LAST_SERUM_CREAT < 0.3 ~ "<0.3",
                                  CAN_LAST_SERUM_CREAT >= 0.3 & CAN_LAST_SERUM_CREAT < 0.7 ~ "0.3-0.7",
                                  CAN_LAST_SERUM_CREAT >= 0.7 & CAN_LAST_SERUM_CREAT < 1 ~ "0.7-1",
                                  CAN_LAST_SERUM_CREAT >= 1 ~ ">=1",
                                  TRUE ~ NA),
                        levels = c("<0.3","0.3-0.7","0.7-1",">=1")),
         donor_ddavp = factor(case_when(DON_DDAVP == "Y" ~ "Yes",
                                        DON_DDAVP == "N" ~ "No",
                                      TRUE ~ NA),
                              levels = c("No","Yes")),
         donor_death_circum = factor(case_when(DON_DEATH_CIRCUM == 1 ~ "Motor Vehicle Accident",
                                               DON_DEATH_CIRCUM == 3 ~ "Homicide",
                                               DON_DEATH_CIRCUM == 5 ~ "Accident (Non-MVA)",
                                               DON_DEATH_CIRCUM == 6 ~ "Natural Causes",
                                               DON_DEATH_CIRCUM %in% c(997,4,2) ~ "Other",
                                               TRUE ~ NA),
                                     levels = c("Motor Vehicle Accident","Homicide",
                                                "Accident (Non-MVA)","Natural Causes","Other")),
         donor_death_mech = factor(case_when(DON_DEATH_MECH == 1 ~ "Drowning",
                                             DON_DEATH_MECH == 2 ~ "Seizure",
                                             DON_DEATH_MECH == 3 ~ "Drug Intoxication",
                                             DON_DEATH_MECH == 4 ~ "Asphyxiation",
                                             DON_DEATH_MECH == 5 ~ "Cardiovascular",
                                             DON_DEATH_MECH %in% c(995, 7, 8) ~ "Gunshot or Stab Wound",
                                             DON_DEATH_MECH == 9 ~ "Blunt Injury",
                                             DON_DEATH_MECH == 11 ~ "Intracranial Hemorrhage",
                                             DON_DEATH_MECH == 12 ~ "Natural Causes",
                                             DON_DEATH_MECH %in% c(6, 997, 10) ~ "Other Causes",
                                             TRUE ~ NA),
                                   levels = c("Blunt Injury","Drowning","Seizure",
                                              "Drug Intoxication","Asphyxiation",
                                              "Cardiovascular","Gunshot or Stab Wound",
                                              "Intracranial Hemorrhage","Natural Causes",
                                              "Other Causes")),
         donor_hist_diab = factor(case_when(DON_HIST_DIAB %in% c("2", "3", "4", "5") ~ "Yes",
                                            DON_HIST_DIAB == "1" ~ "No",
                                                TRUE ~ NA),
                                  levels = c("No","Yes")),
         diagnosis = factor(case_when(REC_DGN %in% c("4100", "4101", "4102", 
                                                      "4103", "4104", "4105", 
                                                      "4106", "4107", "4108", 
                                                      "4110", "4200") ~ "Acute Liver Failure",
                                      REC_DGN %in% c("4300", "4301", "4302", 
                                                      "4303", "4304", "4305", 
                                                      "4306", "4307", "4308", 
                                                      "4315") ~ "Metabolic Liver Disease",
                                      REC_DGN %in% c("4230", "4231", "4235", 
                                                      "4240", "4241", "4242", 
                                                      "4245", "4250", "4255", 
                                                      "4260", "4264") ~ "Cholestatic Liver Disease",
                                      REC_DGN %in% c("4270", "4271", "4272", 
                                                      "4275") ~ "Biliary Atresia",
                                      REC_DGN %in% c("4400", "4401", "4402", 
                                                      "4403", "4404", "4405", 
                                                      "4410", "4420", "4430", 
                                                      "4450", "4451", "4455") ~ "Liver Cancers",
                                      REC_DGN %in% c("4201", "4202", "4203", 
                                                      "4204", "4205", "4206", 
                                                      "4207", "4208", "4209", 
                                                      "4210", "4212", "4213", 
                                                      "4214", "4215", "4216", 
                                                      "4217", "4218", "4219", 
                                                      "4220", "4230", "4231", 
                                                      "4235", "4240", "4241", 
                                                      "4242", "4245", "4250", 
                                                      "4255", "4260", "4264", 
                                                      "4265", "4270", "4271", 
                                                      "4272", "4275", "4280", 
                                                      "4285", "4290", "4300", 
                                                      "4301", "4302", "4303", 
                                                      "4304", "4305", "4306", 
                                                      "4307", "4308", "4315", 
                                                      "4400", "4401", "4402", 
                                                      "4403", "4404", "4405", 
                                                      "4410", "4420", "4430", 
                                                      "4450", "4451", "4455", 
                                                      "4500", "4510", "4520", 
                                                      "4592", "4593") ~ "Viral Hepatitis",
                                      REC_DGN %in% c("4210", "4212", "4213", 
                                                      "4214", "4215", "4216", 
                                                      "4217", "4218", "4219", 
                                                      "4220", "4230", "4231", 
                                                      "4235", "4240", "4241", 
                                                      "4242", "4245") ~ "Autoimmune Disease",
                                      REC_DGN %in% c("4214", "4215", "4216", 
                                                      "4217", "4218", "4219") ~ "NASH/Alcoholic Liver Disease",
                                      REC_DGN %in% c("4208", "4213", "4265", 
                                                      "4280", "4285", "4290", 
                                                      "4500", "4510", "4520", 
                                                      "4597", "4598", "999") ~ "Other")),
         dialysis_within_last_week = factor(case_when(CAN_LAST_DIAL_PRIOR_WEEK == "Y" ~ "Yes",
                                                      CAN_LAST_DIAL_PRIOR_WEEK == "N" ~ "No",
                                            TRUE ~ NA),
                                            levels = c("No","Yes")),
         ebv_status = factor(case_when(REC_EBV_STAT == "P" ~ "Positive",
                                       REC_EBV_STAT == "N" ~ "Negative",
                                       REC_EBV_STAT == "ND" ~ "Not Determined",
                                       TRUE ~ NA),
                             levels = c("Positive","Negative","Not Determined")),
         last_enceph = factor(case_when(CAN_LAST_ENCEPH == 1 ~ "None",
                                        CAN_LAST_ENCEPH == 2 ~ "1-2",
                                        CAN_LAST_ENCEPH == 3 ~ "3-4",
                                        TRUE ~ NA),
                              levels = c("None", "1-2", "3-4")),
         race = factor(case_when(CAN_RACE_SRTR == "ASIAN" ~ "Asian",
                                 CAN_RACE_SRTR == "BLACK" ~ "Black",
                                 CAN_RACE_SRTR == "MULTI" ~ "Multi-Racial",
                                 CAN_RACE_SRTR == "NATIVE" ~ "Native American",
                                 CAN_RACE_SRTR == "PACIFIC" ~ "Pacific Islander",
                                 CAN_RACE_SRTR == "WHITE" ~ "White",
                                 TRUE ~ NA),
                       levels = c("Asian","Black","Multi-Racial","Native American","Pacific Islander",
                                  "White")),
         ethnicity = factor(case_when(CAN_ETHNICITY_SRTR == "LATINO" ~ "Latino",
                                      CAN_ETHNICITY_SRTR == "NLATIN" ~ "Non-Latino or Unknown",
                                      TRUE ~ NA),
                            levels = c("Non-Latino or Unknown","Latino")),
         donor_race = factor(case_when(DON_RACE_SRTR == "ASIAN" ~ "Asian",
                                       DON_RACE_SRTR == "BLACK" ~ "Black",
                                       DON_RACE_SRTR == "MULTI" ~ "Multi-Racial",
                                       DON_RACE_SRTR == "NATIVE" ~ "Native American",
                                       DON_RACE_SRTR == "PACIFIC" ~ "Pacific Islander",
                                       DON_RACE_SRTR == "WHITE" ~ "White",
                                 TRUE ~ NA),
                       levels = c("Asian","Black","Multi-Racial","Native American","Pacific Islander",
                                  "White")),
         donor_ethnicity = factor(case_when(DON_ETHNICITY_SRTR == "LATINO" ~ "Latino",
                                      DON_ETHNICITY_SRTR == "NLATIN" ~ "Non-Latino or Unknown",
                                      TRUE ~ NA),
                            levels = c("Non-Latino or Unknown","Latino")),
         functional_status = factor(case_when(REC_FUNCTN_STAT == 996 | REC_AGE_AT_TX < 1 ~ "N/A (age < 1 year)",
                                              REC_FUNCTN_STAT %in% c(3, 4010, 4020, 4030) ~ "10%-30%",
                                              REC_FUNCTN_STAT %in% c(2, 4040, 4050) ~ "40%-60%",
                                              REC_FUNCTN_STAT %in% c(1, 4060, 4070, 4080, 4090,
                                                                       4100) ~ "70%-100%",
                                              TRUE ~ NA),
                                    levels = c("10%-30%",
                                               "40%-60%",
                                               "70%-100%",
                                               "N/A (age < 1 year)")),
         sex = factor(case_when(CAN_GENDER == "M" ~ "Male",
                                CAN_GENDER == "F" ~ "Female",
                                TRUE ~ NA),
                      levels = c("Male","Female")),
         donor_sex = factor(case_when(DON_GENDER == "M" ~ "Male",
                                      DON_GENDER == "F" ~ "Female",
                                      TRUE ~ NA),
                            levels = c("Male","Female")),
         donor_hist_cancer = factor(case_when(DON_HIST_CANCER == 1 ~ "No",
                                              DON_HIST_CANCER %in% c(2,4,6,9,12,13,14,18,21,29,998,999)~ "Yes",
                                              TRUE ~ NA),
                                    levels = c("No","Yes")),
         donor_hist_cig = factor(case_when(DON_HIST_CIGARETTE_GT20_PKYR == "Y" ~ "Yes",
                                           DON_HIST_CIGARETTE_GT20_PKYR == "N" ~ "No",
                                           TRUE ~ NA),
                                 levels = c("No","Yes")),
         donor_hist_hyperten = factor(case_when(DON_HIST_HYPERTEN == 1 ~ "No",
                                                DON_HIST_HYPERTEN %in% c(2,3,4,5) ~ "Yes",
                                                TRUE ~ NA),
                                      levels = c("No","Yes")),
         donor_hist_mi = factor(case_when(DON_HIST_PREV_MI == "N" ~ "No",
                                          DON_HIST_PREV_MI == "Y" ~ "Yes",
                                          TRUE ~ NA),
                                levels = c("No","Yes")),
         donor_inotrop = factor(case_when(DON_INOTROP_SUPPORT == "Y" ~ "Yes",
                                          DON_INOTROP_SUPPORT == "N" ~ "No",
                                          TRUE ~ NA),
                                levels = c("No","Yes")),
         inr = factor(case_when(CAN_LAST_INR < 2 ~ "<2",
                                CAN_LAST_INR >= 2 & CAN_LAST_INR <= 3 ~ "2-3",
                                CAN_LAST_INR > 3 ~ ">3",
                                TRUE ~ NA),
                      levels = c("<2", "2-3",">3")),
         donor_intracranial_cancer = factor(case_when(DON_INTRACRANIAL_CANCER == "Y" ~ "Yes",
                                                      DON_INTRACRANIAL_CANCER == "N" ~ "No",
                                                      TRUE ~ NA),
                                            levels = c("No","Yes")),
         life_support = factor(case_when(REC_LIFE_SUPPORT == "N" ~ "No",
                                         REC_LIFE_SUPPORT == "Y" ~ "Yes",
                                         TRUE ~ NA),
                               levels = c("No","Yes")),
         prev_malig = factor(case_when(CAN_MALIG == "N" ~ "No",
                                       CAN_MALIG == "Y" ~ "Yes",
                                       TRUE ~ NA),
                             levels = c("No","Yes")),
         med_condition = factor(case_when(REC_MED_COND == 1 ~ "Hospitalized - ICU",
                                          REC_MED_COND == 2 ~ "Hospitalized - Non-ICU",
                                          REC_MED_COND == 3 ~ "Not Hospitalized",
                                          TRUE ~ NA),
                                levels = c("Hospitalized - ICU","Hospitalized - Non-ICU","Not Hospitalized")),
         donor_dcd = factor(case_when(DON_NON_HR_BEAT == "Y" ~ "Yes",
                                      DON_NON_HR_BEAT == "N" ~ "No",
                                      TRUE ~ NA),
                            levels = c("No","Yes")),
         ventilator = factor(case_when(REC_VENTILATOR == 1 ~ "Yes",
                                        REC_VENTILATOR == 0 ~ "No",
                                        TRUE ~ NA),
                                       levels = c("No","Yes")),
         pvt = factor(case_when(REC_PORTAL_VEIN == "Y" ~ "Yes",
                                REC_PORTAL_VEIN == "N" ~ "No",
                                TRUE ~ NA),
                      levels = c("No","Yes")),
         prev_ab_surg = factor(case_when(REC_PREV_ABDOM_SURG == "Y" ~ "Yes",
                                         REC_PREV_ABDOM_SURG == "N" ~ "No",
                                         TRUE ~ NA),
                               levels = c("No","Yes")),
         donor_sgot = factor(case_when(DON_SGOT <= 40 ~ "≤40",
                                       DON_SGOT > 40 & DON_SGOT <= 80 ~ "41-80",
                                       DON_SGOT > 80 & DON_SGOT <= 120 ~ "81-120",
                                       DON_SGOT > 120 & DON_SGOT <= 160 ~ "121-160",
                                       DON_SGOT > 161 ~ ">160",
                                       TRUE ~ NA),
                             levels = c("≤40","41-80","81-120","121-160",">160")),
         bilirubin = factor(case_when(CAN_LAST_BILI >= 0 & CAN_LAST_BILI < 2 ~ "0-2",
                                      CAN_LAST_BILI >= 2 & CAN_LAST_BILI < 6 ~ "2-6",
                                      CAN_LAST_BILI >= 6 & CAN_LAST_BILI < 10 ~ "6-10",
                                      CAN_LAST_BILI >= 10 & CAN_LAST_BILI < 15 ~ "10-15",
                                      CAN_LAST_BILI >= 15 & CAN_LAST_BILI < 20 ~ "15-20",
                                      CAN_LAST_BILI > 20 ~ ">20",
                                      TRUE ~ NA),
                            levels = c("0-2","2-6","6-10","10-15","15-20",">20")),
        tipss = factor(case_when(REC_TIPSS == "N" ~ "No",
                                 REC_TIPSS == "Y" ~ "Yes",
                                 TRUE ~ NA),
                       levels = c("No","Yes")),
        donor_org_shared = factor(case_when(DON_ORG_SHARED == 1 ~ "Shared",
                                            DON_ORG_SHARED == 0 ~ "Local",
                                            TRUE ~ NA),
                                  levels = c("Shared","Local")),
        donor_cmv_status = factor(case_when(DON_ANTI_CMV == "P" ~ "Positive",
                                            DON_ANTI_CMV == "N" ~ "Negative",
                                            TRUE ~ NA),
                                  levels = c("Positive","Negative")),
        prev_tx = factor(case_when(CAN_PREV_TX == 1 ~ "Yes",
                                   CAN_PREV_TX == 0 ~ "No",
                                   TRUE ~ NA),
                         levels = c("No","Yes")),
        donor_heavy_alc = factor(case_when(DON_HEAVY_ALCOHOL == "Y" ~ "Yes",
                                           DON_HEAVY_ALCOHOL == "N" ~ "No",
                                           TRUE ~ NA),
                                 levels = c("No","Yes")),
        donor_arginine = factor(case_when(DON_ARGININE == "Y" ~ "Yes",
                                          DON_ARGININE == "N" ~ "No",
                                          TRUE ~ NA),
                                levels = c("No","Yes")),
        donor_blood_infect = factor(case_when(DON_INFECT_BLOOD == 1 ~ "Yes",
                                              DON_INFECT_BLOOD == 0 ~ "No",
                                              TRUE ~ NA),
                                    levels = c("No","Yes")),
        donor_bun = factor(case_when(DON_BUN >= 0 & DON_BUN < 6 ~ "0-6",
                                     DON_BUN >= 6 & DON_BUN < 12 ~ "6-12",
                                     DON_BUN >= 12 & DON_BUN < 20 ~ "12-20",
                                     DON_BUN >= 20 ~ ">=20",
                                     TRUE ~ NA),
                           levels = c("0-6","6-12","12-20",">=20")),
        donor_cdc_high_risk = factor(case_when(DON_MEET_CDC_HIGH_RISK == "N" ~ "No",
                                               DON_MEET_CDC_HIGH_RISK == "Y" ~ "Yes",
                                               TRUE ~ NA),
                                     levels = c("No","Yes")),
        donor_clinical_infect = factor(case_when(DON_CLINICAL_INFECT == "Y" ~ "Yes",
                                                 DON_CLINICAL_INFECT == "N" ~ "No",
                                                 TRUE ~ NA),
                                       levels = c("No","Yes")),
        donor_cod = factor(case_when(DON_CAD_DON_COD == "1" ~ "Anoxia",
                                     DON_CAD_DON_COD == "2" ~ "Stroke",
                                     DON_CAD_DON_COD == "3" ~ "Head Trauma",
                                     DON_CAD_DON_COD == "4" ~ "CNS Tumor",
                                     DON_CAD_DON_COD %in% c("999","4") ~ "Other",
                                     TRUE ~ NA),
                           levels = c("Head Trauma","Anoxia","Stroke","Other")),
        donor_infect_other = factor(case_when(DON_INFECT_OTHER == 1 ~ "Yes",
                                              DON_INFECT_OTHER == 0 ~ "No",
                                              TRUE ~ NA),
                                    levels = c("No","Yes")),
        donor_protein_urine = factor(case_when(DON_PROTEIN_URINE == "Y" ~ "Yes",
                                               DON_PROTEIN_URINE == "N" ~ "No",
                                               TRUE ~ NA),
                                     levels = c("No","Yes")),
        donor_diuretics = factor(case_when(DON_PRERECOV_DIURETICS == "Y" ~ "Yes",
                                           DON_PRERECOV_DIURETICS == "N" ~ "No",
                                           TRUE ~ NA),
                                 levels = c("No","Yes")),
        donor_steroids = factor(case_when(DON_PRERECOV_STEROIDS == "Y" ~ "Yes",
                                          DON_PRERECOV_STEROIDS == "N" ~ "No",
                                          TRUE ~ NA),
                                levels = c("No","Yes")),
        donor_t3 = factor(case_when(DON_PRERECOV_T3 == "Y" ~ "Yes",
                                    DON_PRERECOV_T3 == "N" ~ "No",
                                    TRUE ~ NA),
                          levels = c("No","Yes")),
        donor_t4 = factor(case_when(DON_PRERECOV_T4 == "Y" ~ "Yes",
                                    DON_PRERECOV_T4 == "N" ~ "No",
                                    TRUE ~ NA),
                          levels = c("No","Yes")),
        donor_infect_lung = factor(case_when(DON_INFECT_LU == 1 ~ "Yes",
                                             DON_INFECT_LU == 0 ~ "No",
                                             TRUE ~ NA),
                                   levels = c("No","Yes")),
        donor_recov_out_us = factor(case_when(DON_RECOV_OUT_US == "Y" ~ "Yes",
                                              DON_RECOV_OUT_US == "N" ~ "No",
                                              TRUE ~ NA),
                                    levels = c("No","Yes")),
        donor_bilirubin = factor(case_when(DON_TOT_BILI >= 0 & DON_TOT_BILI < 1.2 ~ "0-1.2",
                                           DON_TOT_BILI >= 1.2 & DON_TOT_BILI < 2 ~ "1.2-2",
                                           DON_TOT_BILI >= 2 ~ ">=2",
                                           TRUE ~ NA),
                                 levels = c("0-1.2","1.2-2",">=2")),
        donor_infect_urine = factor(case_when(DON_INFECT_URINE == 1 ~ "Yes",
                                              DON_INFECT_URINE == 0 ~ "No",
                                              TRUE ~ NA),
                                    levels = c("No","Yes")),
        donor_vasodil = factor(case_when(DON_VASODIL == "Y" ~ "Yes",
                                         DON_VASODIL == "N" ~ "No",
                                         TRUE ~ NA),
                               levels = c("No","Yes")),
        waitlist_time = as.numeric(round(difftime(REC_TX_DT, CAN_LISTING_DT, units = "days") / 30.45,0.1)),
        waitlist_time_cat = as_factor(case_when(waitlist_time >= 0 & waitlist_time < 6 ~ "0-6",
                                                waitlist_time >= 6 & waitlist_time < 12 ~ "6-12",
                                                waitlist_time >= 12 & waitlist_time < 24 ~ "12-24",
                                                waitlist_time >= 24  ~ ">=24",
                                                TRUE ~ NA)),
        tx_year = factor(substr(REC_TX_DT,1,4),
                         levels = c("2017","2010","2011","2012","2013",
                                    "2014","2015","2016","2018","2019","2020","2021","2022","2023",
                                    "2024","2025")),
        donor_type = factor(case_when(REC_LI_TY %in% c(1,2,3,4,5) ~ "Partial",
                                      REC_LI_TY %in% c(6,7,8,9,10,11,12,13,14,15) ~ "Split",
                                      is.na(REC_LI_TY) ~ "Whole",
                                      TRUE ~ NA),
                            levels = c("Whole","Partial","Split")))

#inspect the patients that had no follow-up time. remove patients with no follow-up time
no_fu_time <- liver_cohort_clean %>%
  select(gs_outcome, gs_time)

liver_cohort_clean <- liver_cohort_clean %>%
  filter(gs_time > 0)
#280


#add an exception variable from the mp_except data
exception_data <- mp_except %>%
  group_by(PX_ID) %>%
  summarise(exception = case_when(any(CANHX_MPXCPT_STAT %in% c(5,14,16,18,30)) ~ "Yes",
                                      TRUE ~ NA)) 
duplicates <- exception_data %>%
  group_by(PX_ID) %>%
  filter(n() > 1) %>%
  ungroup()


liver_cohort_clean <- left_join(liver_cohort_clean, exception_data, by = "PX_ID")
table(liver_cohort_clean$exception)

liver_cohort_clean <- liver_cohort_clean %>%
  mutate(exception = factor(case_when(is.na(exception) ~ "No",
                                      TRUE ~ exception),
                            levels = c("No","Yes")))
table(liver_cohort_clean$exception)



#create a growth failure variables

# age bands from the MELD/PELD table: [0,3), [3,6), ..., [216,219)
breaks <- seq(0, 219, by = 3)

# Male cutoffs from the documentation
male_height <- c(
  46.02, 55.80, 62.40, 67.02, 70.70, 73.70, 76.32, 78.64, 79.24, 81.38,
  83.36, 85.36, 87.28, 89.14, 91.00, 92.68, 94.38, 96.08, 97.72, 99.26,
  100.70, 102.16, 103.62, 105.00, 106.38, 107.66, 109.04, 110.22, 111.50,
  112.68, 113.96, 115.14, 116.30, 117.46, 118.60, 119.72, 120.84, 121.96,
  123.04, 124.12, 125.28, 126.40, 127.52, 128.72, 129.88, 131.02, 132.24,
  133.42, 134.68, 135.86, 137.16, 138.58, 139.94, 141.36, 142.88, 144.36,
  145.98, 147.68, 149.40, 151.12, 152.88, 154.60, 156.24, 157.78, 159.12,
  160.42, 161.58, 162.46, 163.06, 163.46, 163.66, 163.72, 163.60
)

male_weight <- c(
  2.48, 4.10, 5.84, 7.26, 8.16, 8.72, 9.14, 9.48, 9.64, 10.10,
  10.58, 11.08, 11.40, 11.82, 12.16, 12.52, 12.90, 13.28, 13.68, 14.06,
  14.46, 14.86, 15.26, 15.66, 16.04, 16.40, 16.76, 17.22, 17.64, 17.94,
  18.32, 18.78, 19.10, 19.50, 19.86, 20.20, 20.52, 20.90, 21.28, 21.74,
  22.10, 22.64, 23.08, 23.60, 24.14, 24.78, 25.40, 26.06, 26.82, 27.58,
  28.46, 29.38, 30.40, 31.44, 32.52, 33.72, 34.94, 36.18, 37.44, 38.62,
  39.94, 41.16, 42.40, 43.54, 44.72, 45.80, 46.80, 47.72, 48.54, 49.28,
  49.94, 50.48, 50.94
)

# Female cutoffs from the documentation
female_height <- c(
  45.56, 54.52, 60.62, 64.94, 68.62, 71.90, 74.76, 77.44, 78.00, 80.28,
  82.42, 84.52, 86.48, 88.32, 90.12, 91.90, 93.56, 95.10, 96.60, 98.10,
  99.58, 100.94, 102.28, 103.52, 104.76, 106.08, 107.20, 108.40, 109.62,
  110.84, 111.96, 113.08, 114.32, 115.46, 116.72, 118.00, 119.18, 120.48,
  121.82, 123.26, 124.64, 126.14, 127.66, 129.22, 130.92, 132.64, 134.40,
  136.18, 137.86, 139.54, 141.12, 142.50, 143.74, 144.86, 145.66, 146.46,
  147.02, 147.40, 147.76, 148.04, 148.32, 148.52, 148.64, 148.76, 149.08,
  149.40, 149.66, 150.02, 150.30, 150.72, 151.12, 151.50, 151.78
)

female_weight <- c(
  2.20, 3.92, 5.46, 6.66, 7.40, 7.96, 8.46, 8.98, 9.46, 9.88,
  10.32, 10.78, 11.14, 11.52, 11.92, 12.20, 12.60, 12.88, 13.18, 13.46,
  13.84, 14.10, 14.44, 14.68, 14.98, 15.28, 15.64, 15.98, 16.30, 16.68,
  17.12, 17.44, 17.84, 18.32, 18.78, 19.24, 19.76, 20.20, 20.82, 21.34,
  21.86, 22.48, 23.12, 23.76, 24.50, 25.18, 25.86, 26.66, 27.40, 28.26,
  29.02, 29.96, 30.80, 31.68, 32.58, 33.40, 34.26, 35.04, 35.88, 36.62,
  37.40, 38.12, 38.76, 39.34, 39.86, 40.30, 40.66, 41.06, 41.36, 41.56,
  41.78, 41.98, 42.08
)


# assign age interval index based on completed months
# intervals are [0,3), [3,6), ..., [216,219)
data_growthfailure <- liver_cohort_clean %>%
  mutate(
    age_band = cut(
      REC_AGE_IN_MONTHS_AT_TX,
      breaks = breaks,
      right = FALSE,
      include.lowest = TRUE,
      labels = FALSE
    ),
    low_height = case_when(
      sex == "Male"   ~ male_height[age_band],
      sex == "Female" ~ female_height[age_band],
      TRUE ~ NA_real_
    ),
    low_weight = case_when(
      sex == "Male"   ~ male_weight[age_band],
      sex == "Female" ~ female_weight[age_band],
      TRUE ~ NA_real_
    )
  )

# growth failure per MELD/PELD documentation:
# 1 if height <= low_height OR weight <= low_weight
# 0 if available values are above thresholds
# NA only if sex/age invalid or both height and weight missing
data_growthfailure <- data_growthfailure %>%
  mutate(
    growth_failure = case_when(
      is.na(sex) | is.na(age_band) ~ NA_real_,
      is.na(REC_HGT_CM) & is.na(REC_WGT_KG) ~ NA_real_,
      (!is.na(REC_HGT_CM) & REC_HGT_CM <= low_height) |
        (!is.na(REC_WGT_KG) & REC_WGT_KG <= low_weight) ~ 1,
      TRUE ~ 0
    )
  )


data_growthfailure <- data_growthfailure %>%
  select(PERS_ID, growth_failure)

liver_cohort_clean <- left_join(liver_cohort_clean, data_growthfailure, by = "PERS_ID")

liver_cohort_clean <- liver_cohort_clean %>%
  mutate(growth_failure = factor(case_when(growth_failure == 1 ~ "Yes",
                                           growth_failure == 0 ~ "No",
                                           TRUE ~ NA),
                                 levels = c("Yes","No")))


#remove variables if they have < 1% prevalence in the data set
tbl_summary(
  data = liver_cohort_clean,
  include = c(age_at_transplant_cat, donor_age_at_tx_cat, blood_type, donor_blood_type, albumin, artificial_liver, bmi_cat, cold_ischemic_time, donor_creat, creat, donor_ddavp, donor_death_circum, donor_death_mech, donor_hist_diab, diagnosis, dialysis_within_last_week, ebv_status, last_enceph, functional_status, sex, donor_sex, donor_hist_cancer, donor_hist_cig, donor_hist_hyperten, donor_hist_mi, donor_inotrop, inr, donor_intracranial_cancer, life_support, prev_malig, med_condition, donor_dcd, ventilator, pvt, prev_ab_surg, donor_sgot, bilirubin, tipss, donor_org_shared, exception, donor_cmv_status, prev_tx, donor_heavy_alc, donor_arginine, donor_blood_infect, donor_bun, donor_cdc_high_risk, donor_clinical_infect, donor_cod, donor_infect_other, donor_protein_urine, donor_diuretics, donor_steroids, donor_t3, donor_t4, donor_infect_lung, donor_recov_out_us, donor_bilirubin, donor_infect_urine, donor_vasodil, waitlist_time_cat, growth_failure, tx_year, donor_type, donor_bmi_cat, gs_outcome),
  by = gs_outcome
) %>%
  bold_labels() %>%
  add_overall()

liver_cohort_clean <- liver_cohort_clean %>%
  select(-c(artificial_liver, donor_hist_mi, donor_intracranial_cancer, donor_dcd, donor_hist_cancer, donor_recov_out_us, donor_t3))



#remove variables that are collected and applicable for adults but not necessarily adults and peds.
tbl_summary(
  data = liver_cohort_clean,
  include = c(age_at_transplant_cat, donor_age_at_tx_cat, blood_type, donor_blood_type, albumin, bmi_cat, cold_ischemic_time, donor_creat, creat, donor_ddavp, donor_death_circum, donor_death_mech, donor_hist_diab, diagnosis, dialysis_within_last_week, ebv_status, last_enceph, functional_status, sex, donor_sex, donor_hist_cig, donor_hist_hyperten, donor_inotrop, inr, life_support, prev_malig, med_condition, ventilator, pvt, prev_ab_surg, donor_sgot, bilirubin, tipss, donor_org_shared, exception, donor_cmv_status, prev_tx, donor_heavy_alc, donor_arginine, donor_blood_infect, donor_bun, donor_cdc_high_risk, donor_clinical_infect, donor_cod, donor_infect_other, donor_protein_urine, donor_diuretics, donor_steroids, donor_t4, donor_infect_lung, donor_bilirubin, donor_infect_urine, donor_vasodil, waitlist_time_cat, growth_failure, tx_year, donor_type, donor_bmi_cat),
  by = donor_age_at_tx_cat
) %>%
  bold_labels() %>%
  add_overall()

#donor_hist_cig, donor_heavy_alc

liver_cohort_clean <- liver_cohort_clean %>%
  select(-c(donor_hist_cig, donor_heavy_alc))




# filter the liver_cohort_clean data to just the variables that I want
liver_cohort_clean <- liver_cohort_clean %>%
  select(gs_time, gs_outcome, age_at_transplant_cat, age_at_transplant, donor_age_at_tx_cat, blood_type, donor_blood_type, albumin, bmi_cat, donor_bmi_cat, cold_ischemic_time, donor_creat, creat, donor_ddavp, donor_death_circum, donor_death_mech, donor_hist_diab, diagnosis, dialysis_within_last_week, ebv_status, last_enceph, functional_status, sex, donor_sex, donor_hist_hyperten, donor_inotrop, inr, life_support, prev_malig, med_condition, ventilator, pvt, prev_ab_surg, donor_sgot, bilirubin, tipss, donor_org_shared, exception, donor_cmv_status, prev_tx, donor_arginine, donor_blood_infect, donor_bun, donor_cdc_high_risk, donor_clinical_infect, donor_cod, donor_infect_other, donor_protein_urine, donor_diuretics, donor_steroids, donor_t4, donor_infect_lung, donor_bilirubin, donor_infect_urine, donor_vasodil, waitlist_time_cat, growth_failure, tx_year, donor_type)


#inspect variable levels
lapply(liver_cohort_clean, table)


#re-inspect missing data
liver_cohort_clean_miss <- miss_var_summary(liver_cohort_clean)





#############################
#perform multiple imputation#
#############################

pred <- quickpred(liver_cohort_clean)

pred[ ,"age_at_transplant"] <- 0
pred["age_at_transplant", ] <- 0


liver_cohort_imputed <- mice(liver_cohort_clean, m = 5, seed = 1145, printFlag = F,
                             predictorMatrix = pred)

#check to see if the algorithm converged
liver_cohort_imputed$loggedEvents


#inspection of imputed values
plot(liver_cohort_imputed, layout = c(5,5))

# - There is a lack of trends in the means which suggests that the algorithm converged. 


#extract the first imputed dataset
data_imputed_1 <- complete(liver_cohort_imputed, 1)
miss_data_imputed_1 <- miss_var_summary(data_imputed_1)

analysis_sens <- data_imputed_1


#write out the analysis dataset
save(analysis_sens, file = "L:/Projects/Investigator/Jonathan Merola/data/analysis_sens.RData")