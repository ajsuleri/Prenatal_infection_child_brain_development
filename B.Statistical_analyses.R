#------------------------------------------------------------------------#
#############PROJECT: infection & child brain development#################
#------------------------------------------------------------------------#

# Author: Anna Suleri

### Parts in this script:
#' Part 1: Transformations, mixed model, with interaction term, raw brain volumes, (with and without icv) + also main effect
#' Part 2: Additional analyses after peer review revisions 

# Of note, we adjust for fdr within each subset of analyses

# Clear environment
rm(list = ls()) 

# Opening libraries we need & setting wd 
setwd('Set_wd_path_to_data_folder')
wd <- getwd()

libraries <- c('foreign', 'haven', 'dplyr', 'tidyverse', 'tidyr', 'broom.mixed', 'xlsx', 'corrplot', 'stringi', 'miceadds', 'mitools', 'lme4', 'rlang', 'splines')

invisible(lapply(libraries, require, character.only = T))

#\

###############################################
#-------------------PART 1--------------------#
###############################################

# Load data of final sample 
load('no_outcome_imp_ipw_mids_final.RDS') 

# Set wd to results folder 
setwd('Set_wd_to_results_folder') 

# Standardizing vars in final sample
imp.test_long <- complete(ipw_mids_final, include = T, action = "long")

outcome_vars <- c("genr_tbv","Cerebellum_Cortex_vol_subcortical","Amygdala_vol_subcortical",
                  "Hippocampus_vol_subcortical","Caudate_vol_subcortical",
                  "Putamen_vol_subcortical", "Thalamus_Proper_vol_subcortical"
                  ,"Pallidum_vol_subcortical",
                  "bankssts_vol_cortical",
                  "caudalanteriorcingulate_vol_cortical",
                  "caudalmiddlefrontal_vol_cortical",
                  "cuneus_vol_cortical","entorhinal_vol_cortical",
                  "fusiform_vol_cortical",
                  "inferiorparietal_vol_cortical","inferiortemporal_vol_cortical",
                  "isthmuscingulate_vol_cortical","lateraloccipital_vol_cortical",
                  "lateralorbitofrontal_vol_cortical",
                  "lingual_vol_cortical","medialorbitofrontal_vol_cortical",
                  "middletemporal_vol_cortical",
                  "parahippocampal_vol_cortical","paracentral_vol_cortical",
                  "parsopercularis_vol_cortical",
                  "parsorbitalis_vol_cortical","parstriangularis_vol_cortical",
                  "pericalcarine_vol_cortical",
                  "postcentral_vol_cortical","posteriorcingulate_vol_cortical.x",
                  "precentral_vol_cortical.x",
                  "precuneus_vol_cortical.x","rostralanteriorcingulate_vol_cortical.x",
                  "rostralmiddlefrontal_vol_cortical.x","superiorfrontal_vol_cortical.x",
                  "superiorparietal_vol_cortical.x",
                  "superiortemporal_vol_cortical.x","supramarginal_vol_cortical.x",
                  "frontalpole_vol_cortical.x",
                  "temporalpole_vol_cortical.x", "transversetemporal_vol_cortical.x" ,"insula_vol_cortical.x") 

for (x in outcome_vars){ #making z-scores of raw brain volumes 
  t <- imp.test_long[x]
  colname <- paste0(colnames(t), "_standardized")
  imp.test_long[,colname] <- as.numeric(scale(t))
}

# Check structure of vars and recode if needed 
imp.test_long$timepoint <- as.factor(imp.test_long$timepoint) 
imp.test_long$timepoint <- as.numeric(imp.test_long$timepoint)

imp.test_long$sumscore_inf_tot <- as.numeric(imp.test_long$sumscore_inf_tot)
imp.test_long$sumscore_inf_tot_standardized <- as.numeric(scale(imp.test_long$sumscore_inf_tot))
imp.test_long$sumscore_inf_tri1 <- as.numeric(imp.test_long$sumscore_inf_tri1)
imp.test_long$sumscore_inf_tri2 <- as.numeric(imp.test_long$sumscore_inf_tri2)
imp.test_long$sumscore_inf_tri3 <- as.numeric(imp.test_long$sumscore_inf_tri3)

imp.test_long$sumscore_inf_tot_standardized <- as.numeric(scale(imp.test_long$sumscore_inf_tot))
imp.test_long$sumscore_inf_tri1_standardized <- as.numeric(scale(imp.test_long$sumscore_inf_tri1))
imp.test_long$sumscore_inf_tri2_standardized <- as.numeric(scale(imp.test_long$sumscore_inf_tri2))
imp.test_long$sumscore_inf_tri3_standardized <- as.numeric(scale(imp.test_long$sumscore_inf_tri3))

# Convert back to mids object 
imp.test_mids <- as.mids(imp.test_long) 

### Main analysis

# Select outcomes to loop over 
outcome_vars_ztrans <- paste0(outcome_vars, "_standardized")  

#Mixed model for total infection and every outcome (standardized coefficient)
results_mm_raw_brain <- data.frame() #create empty df 

for(x in outcome_vars_ztrans) {
  f <- paste0(x, "~ sumscore_inf_tot_standardized*timepoint + GENDER + AGE_M_v2 + ETHNMv2 + EDUCM_3groups+ GSI + SMOKE_ALL + (1|IDC)") #full model, with interaction, only with random intercept, without icv
  
  g <- paste0(x, "~ sumscore_inf_tot_standardized*timepoint + eTIV + GENDER + AGE_M_v2 + ETHNMv2 + EDUCM_3groups+ GSI + SMOKE_ALL + (1|IDC)") #full model, with interaction, only with random intercept, with icv
  
  h <- paste0(x, "~ sumscore_inf_tot_standardized + timepoint + GENDER + AGE_M_v2 + ETHNMv2 + EDUCM_3groups+ GSI + SMOKE_ALL + (1|IDC)") #full model, only with random intercept, main effect, without icv 
  
  #calculating beta, SE, confidence interval and pval for each model
  bval_int_no_icv <- summary(pool(with(imp.test_mids, lmer(as.formula(f), weights = combined_weights))),conf.int = T)[12,2]
  seval_int_no_icv <- summary(pool(with(imp.test_mids, lmer(as.formula(f), weights = combined_weights))),conf.int = T)[12,3]
  lowerCI_int_no_icv <- summary(pool(with(imp.test_mids, lmer(as.formula(f), weights = combined_weights))),conf.int = T)[12,7]
  upperCI_int_no_icv <- summary(pool(with(imp.test_mids, lmer(as.formula(f), weights = combined_weights))),conf.int = T)[12,8]
  pval_int_no_icv <- summary(pool(with(imp.test_mids, lmer(as.formula(f), weights = combined_weights))), conf.int = T)[12,6]
  
  #calculating beta, SE, confidence interval and pval for each model
  bval_int_icv <- summary(pool(with(imp.test_mids, lmer(as.formula(g), weights = combined_weights))),conf.int = T)[13,2]
  seval_int_icv <- summary(pool(with(imp.test_mids, lmer(as.formula(g), weights = combined_weights))),conf.int = T)[13,3]
  lowerCI_int_icv <- summary(pool(with(imp.test_mids, lmer(as.formula(g), weights = combined_weights))),conf.int = T)[13,7]
  upperCI_int_icv <- summary(pool(with(imp.test_mids, lmer(as.formula(g), weights = combined_weights))),conf.int = T)[13,8]
  pval_int_icv <- summary(pool(with(imp.test_mids, lmer(as.formula(g), weights = combined_weights))), conf.int = T)[13,6]
  
  #calculating beta, SE, confidence interval and pval for each model
  bval_main <- summary(pool(with(imp.test_mids, lmer(as.formula(h), weights = combined_weights))),conf.int = T)[2,2]
  seval_main <- summary(pool(with(imp.test_mids, lmer(as.formula(h), weights = combined_weights))),conf.int = T)[2,3]
  lowerCI_main <- summary(pool(with(imp.test_mids, lmer(as.formula(h), weights = combined_weights))),conf.int = T)[2,7]
  upperCI_main <- summary(pool(with(imp.test_mids, lmer(as.formula(h), weights = combined_weights))),conf.int = T)[2,8]
  pval_main <- summary(pool(with(imp.test_mids, lmer(as.formula(h), weights = combined_weights))), conf.int = T)[2,6]
  
  #storing values in empty datafile  
  results_mm_raw_brain[x,1] <- bval_int_no_icv
  results_mm_raw_brain[x,2] <- seval_int_no_icv
  results_mm_raw_brain[x,3] <- lowerCI_int_no_icv
  results_mm_raw_brain[x,4] <- upperCI_int_no_icv
  results_mm_raw_brain[x,5] <- pval_int_no_icv
  
  results_mm_raw_brain[x,6] <- bval_int_icv
  results_mm_raw_brain[x,7] <- seval_int_icv
  results_mm_raw_brain[x,8] <- lowerCI_int_icv
  results_mm_raw_brain[x,9] <- upperCI_int_icv
  results_mm_raw_brain[x,10] <- pval_int_icv
  
  results_mm_raw_brain[x,11] <- bval_main
  results_mm_raw_brain[x,12] <- seval_main
  results_mm_raw_brain[x,13] <- lowerCI_main
  results_mm_raw_brain[x,14] <- upperCI_main
  results_mm_raw_brain[x,15] <- pval_main
  
  #assigning names to columns 
  colnames(results_mm_raw_brain) <- c("bval_int_no_icv", "seval_int_no_icv", "lowerCI_int_no_icv", "upperCI_int_no_icv", "pval_int_no_icv", "bval_int_icv","seval_int_icv","lowerCI_int_icv","upperCI_int_icv", "pval_int_icv", "bval_main", "seval_main","lowerCI_main","upperCI_main","pval_main")
}

write.xlsx(results_mm_raw_brain, "all_results_raw_brain_ipw_standardized.xlsx") #saving results to excel 

# Mixed model for total infection and every outcome (unstandardized coefficient)

results_mm_raw_brain_unstandardized <- data.frame() #create empty df 

#Mixed model for total infection and every outcome (unstandardized coefficient)
for(x in outcome_vars) {
  f <- paste0(x, "~ sumscore_inf_tot*timepoint + GENDER + AGE_M_v2 + ETHNMv2 + EDUCM_3groups+ GSI + SMOKE_ALL + (1|IDC)") #full model, with interaction, only with random intercept, without icv
  
  g <- paste0(x, "~ sumscore_inf_tot*timepoint + eTIV + GENDER + AGE_M_v2 + ETHNMv2 + EDUCM_3groups+ GSI + SMOKE_ALL + (1|IDC)") #full model, with interaction, only with random intercept, with icv
  
  h <- paste0(x, "~ sumscore_inf_tot + timepoint + GENDER + AGE_M_v2 + ETHNMv2 + EDUCM_3groups+ GSI + SMOKE_ALL + (1|IDC)") #full model, only with random intercept, main effect, without icv 
  
  #calculating beta, SE, confidence interval and pval for each model
  bval_int_no_icv <- summary(pool(with(imp.test_mids, lmer(as.formula(f), weights = combined_weights))),conf.int = T)[12,2]
  seval_int_no_icv <- summary(pool(with(imp.test_mids, lmer(as.formula(f), weights = combined_weights))),conf.int = T)[12,3]
  lowerCI_int_no_icv <- summary(pool(with(imp.test_mids, lmer(as.formula(f), weights = combined_weights))),conf.int = T)[12,7]
  upperCI_int_no_icv <- summary(pool(with(imp.test_mids, lmer(as.formula(f), weights = combined_weights))),conf.int = T)[12,8]
  pval_int_no_icv <- summary(pool(with(imp.test_mids, lmer(as.formula(f), weights = combined_weights))), conf.int = T)[12,6]
  
  #calculating beta, SE, confidence interval and pval for each model
  bval_int_icv <- summary(pool(with(imp.test_mids, lmer(as.formula(g), weights = combined_weights))),conf.int = T)[13,2]
  seval_int_icv <- summary(pool(with(imp.test_mids, lmer(as.formula(g), weights = combined_weights))),conf.int = T)[13,3]
  lowerCI_int_icv <- summary(pool(with(imp.test_mids, lmer(as.formula(g), weights = combined_weights))),conf.int = T)[13,7]
  upperCI_int_icv <- summary(pool(with(imp.test_mids, lmer(as.formula(g), weights = combined_weights))),conf.int = T)[13,8]
  pval_int_icv <- summary(pool(with(imp.test_mids, lmer(as.formula(g), weights = combined_weights))), conf.int = T)[13,6]
  
  #calculating beta, SE, confidence interval and pval for each model
  bval_main <- summary(pool(with(imp.test_mids, lmer(as.formula(h), weights = combined_weights))),conf.int = T)[2,2]
  seval_main <- summary(pool(with(imp.test_mids, lmer(as.formula(h), weights = combined_weights))),conf.int = T)[2,3]
  lowerCI_main <- summary(pool(with(imp.test_mids, lmer(as.formula(h), weights = combined_weights))),conf.int = T)[2,7]
  upperCI_main <- summary(pool(with(imp.test_mids, lmer(as.formula(h), weights = combined_weights))),conf.int = T)[2,8]
  pval_main <- summary(pool(with(imp.test_mids, lmer(as.formula(h), weights = combined_weights))), conf.int = T)[2,6]
  
  #storing values in empty datafile  
  results_mm_raw_brain_unstandardized[x,1] <- bval_int_no_icv
  results_mm_raw_brain_unstandardized[x,2] <- seval_int_no_icv
  results_mm_raw_brain_unstandardized[x,3] <- lowerCI_int_no_icv
  results_mm_raw_brain_unstandardized[x,4] <- upperCI_int_no_icv
  results_mm_raw_brain_unstandardized[x,5] <- pval_int_no_icv
  
  results_mm_raw_brain_unstandardized[x,6] <- bval_int_icv
  results_mm_raw_brain_unstandardized[x,7] <- seval_int_icv
  results_mm_raw_brain_unstandardized[x,8] <- lowerCI_int_icv
  results_mm_raw_brain_unstandardized[x,9] <- upperCI_int_icv
  results_mm_raw_brain_unstandardized[x,10] <- pval_int_icv
  
  results_mm_raw_brain_unstandardized[x,11] <- bval_main
  results_mm_raw_brain_unstandardized[x,12] <- seval_main
  results_mm_raw_brain_unstandardized[x,13] <- lowerCI_main
  results_mm_raw_brain_unstandardized[x,14] <- upperCI_main
  results_mm_raw_brain_unstandardized[x,15] <- pval_main
  
  #assigning names to columns 
  colnames(results_mm_raw_brain_unstandardized) <- c("bval_int_no_icv", "seval_int_no_icv", "lowerCI_int_no_icv", "upperCI_int_no_icv", "pval_int_no_icv", "bval_int_icv","seval_int_icv","lowerCI_int_icv","upperCI_int_icv", "pval_int_icv", "bval_main", "seval_main","lowerCI_main","upperCI_main","pval_main")
}

write.xlsx(results_mm_raw_brain_unstandardized, "all_results_raw_brain_ipw_unstandardized.xlsx")

# Mixed model for trimester infection and every outcome (standardized coefficient)
results_trimesters_rawbrain <- data.frame() #create empty df 

for(x in outcome_vars_ztrans) {
  f <- paste0(x, "~ sumscore_inf_tri1_standardized*timepoint + AGE_M_v2 + ETHNMv2 + EDUCM_3groups+ GSI + SMOKE_ALL + (1|IDC)") #tri1, with interaction, no icv
  
  g <- paste0(x, "~ sumscore_inf_tri2_standardized*timepoint + AGE_M_v2 + ETHNMv2 + EDUCM_3groups+ GSI + SMOKE_ALL + (1|IDC)") #tri2, with interaction, and icv
  
  h <- paste0(x, "~ sumscore_inf_tri3_standardized*timepoint + AGE_M_v2 + ETHNMv2 + EDUCM_3groups+ GSI + SMOKE_ALL + (1|IDC)") #tri3, with interaction, no icv
  
  #calculating beta, SE, confidence interval and pval for each model
  bval_tri1 <- summary(pool(with(imp.test_mids, lmer(as.formula(f), weights = combined_weights))),conf.int = T)[11,2]
  seval_tri1 <- summary(pool(with(imp.test_mids, lmer(as.formula(f), weights = combined_weights))),conf.int = T)[11,3]
  lowerCI_tri1 <- summary(pool(with(imp.test_mids, lmer(as.formula(f), weights = combined_weights))),conf.int = T)[11,7]
  upperCI_tri1 <- summary(pool(with(imp.test_mids, lmer(as.formula(f), weights = combined_weights))),conf.int = T)[11,8]
  pval_tri1 <- summary(pool(with(imp.test_mids, lmer(as.formula(f), weights = combined_weights))), conf.int = T)[11,6]
  
  #calculating beta, SE, confidence interval and pval for each model
  bval_tri2 <- summary(pool(with(imp.test_mids, lmer(as.formula(g), weights = combined_weights))),conf.int = T)[11,2]
  seval_tri2 <- summary(pool(with(imp.test_mids, lmer(as.formula(g), weights = combined_weights))),conf.int = T)[11,3]
  lowerCI_tri2 <- summary(pool(with(imp.test_mids, lmer(as.formula(g), weights = combined_weights))),conf.int = T)[11,7]
  upperCI_tri2 <- summary(pool(with(imp.test_mids, lmer(as.formula(g), weights = combined_weights))),conf.int = T)[11,8]
  pval_tri2 <- summary(pool(with(imp.test_mids, lmer(as.formula(g), weights = combined_weights))), conf.int = T)[11,6]
  
  #calculating beta, SE, confidence interval and pval for each model
  bval_tri3 <- summary(pool(with(imp.test_mids, lmer(as.formula(h), weights = combined_weights))),conf.int = T)[11,2]
  seval_tri3 <- summary(pool(with(imp.test_mids, lmer(as.formula(h), weights = combined_weights))),conf.int = T)[11,3]
  lowerCI_tri3 <- summary(pool(with(imp.test_mids, lmer(as.formula(h), weights = combined_weights))),conf.int = T)[11,7]
  upperCI_tri3 <- summary(pool(with(imp.test_mids, lmer(as.formula(h), weights = combined_weights))),conf.int = T)[11,8]
  pval_tri3 <- summary(pool(with(imp.test_mids, lmer(as.formula(h), weights = combined_weights))), conf.int = T)[11,6]
  
  #storing values in empty datafile  
  results_trimesters_rawbrain[x,1] <- bval_tri1
  results_trimesters_rawbrain[x,2] <- seval_tri1
  results_trimesters_rawbrain[x,3] <- lowerCI_tri1
  results_trimesters_rawbrain[x,4] <- upperCI_tri1
  results_trimesters_rawbrain[x,5] <- pval_tri1
  
  results_trimesters_rawbrain[x,6] <- bval_tri2
  results_trimesters_rawbrain[x,7] <- seval_tri2
  results_trimesters_rawbrain[x,8] <- lowerCI_tri2
  results_trimesters_rawbrain[x,9] <- upperCI_tri2
  results_trimesters_rawbrain[x,10] <- pval_tri2
  
  results_trimesters_rawbrain[x,11] <- bval_tri3
  results_trimesters_rawbrain[x,12] <- seval_tri3
  results_trimesters_rawbrain[x,13] <- lowerCI_tri3
  results_trimesters_rawbrain[x,14] <- upperCI_tri3
  results_trimesters_rawbrain[x,15] <- pval_tri3
  
  #assigning names to columns 
  colnames(results_trimesters_rawbrain) <- c("bval_tri1", "seval_tri1", "lowerCI_tri1", "upperCI_tri1", "pval_tri1","bval_tri2", 'seval_tri2', "lowerCI_tri2", "upperCI_tri2", "pval_tri2", "bval_tri3", "seval_tri3", "lowerCI_tri3", "upperCI_tri3", "pval_tri3")
}

write.xlsx(results_trimesters_rawbrain, "results_trimesters_rawbrain_standardized.xlsx")

#\

# Mixed model for trimester infection and every outcome (unstandardized coefficient)
results_trimesters_rawbrain_unstandardized <- data.frame() #create empty df 

for(x in outcome_vars) {
  f <- paste0(x, "~ sumscore_inf_tri1*timepoint + AGE_M_v2 + ETHNMv2 + EDUCM_3groups+ GSI + SMOKE_ALL + (1|IDC)") #tri1, with interaction, no icv
  
  g <- paste0(x, "~ sumscore_inf_tri2*timepoint + AGE_M_v2 + ETHNMv2 + EDUCM_3groups+ GSI + SMOKE_ALL + (1|IDC)") #tri2, with interaction, and icv
  
  h <- paste0(x, "~ sumscore_inf_tri3*timepoint + AGE_M_v2 + ETHNMv2 + EDUCM_3groups+ GSI + SMOKE_ALL + (1|IDC)") #tri3, with interaction, no icv
  
  #calculating beta, SE, confidence interval and pval for each model
  bval_tri1 <- summary(pool(with(imp.test_mids, lmer(as.formula(f), weights = combined_weights))),conf.int = T)[11,2]
  seval_tri1 <- summary(pool(with(imp.test_mids, lmer(as.formula(f), weights = combined_weights))),conf.int = T)[11,3]
  lowerCI_tri1 <- summary(pool(with(imp.test_mids, lmer(as.formula(f), weights = combined_weights))),conf.int = T)[11,7]
  upperCI_tri1 <- summary(pool(with(imp.test_mids, lmer(as.formula(f), weights = combined_weights))),conf.int = T)[11,8]
  pval_tri1 <- summary(pool(with(imp.test_mids, lmer(as.formula(f), weights = combined_weights))), conf.int = T)[11,6]
  
  #calculating beta, SE, confidence interval and pval for each model
  bval_tri2 <- summary(pool(with(imp.test_mids, lmer(as.formula(g), weights = combined_weights))),conf.int = T)[11,2]
  seval_tri2 <- summary(pool(with(imp.test_mids, lmer(as.formula(g), weights = combined_weights))),conf.int = T)[11,3]
  lowerCI_tri2 <- summary(pool(with(imp.test_mids, lmer(as.formula(g), weights = combined_weights))),conf.int = T)[11,7]
  upperCI_tri2 <- summary(pool(with(imp.test_mids, lmer(as.formula(g), weights = combined_weights))),conf.int = T)[11,8]
  pval_tri2 <- summary(pool(with(imp.test_mids, lmer(as.formula(g), weights = combined_weights))), conf.int = T)[11,6]
  
  #calculating beta, SE, confidence interval and pval for each model
  bval_tri3 <- summary(pool(with(imp.test_mids, lmer(as.formula(h), weights = combined_weights))),conf.int = T)[11,2]
  seval_tri3 <- summary(pool(with(imp.test_mids, lmer(as.formula(h), weights = combined_weights))),conf.int = T)[11,3]
  lowerCI_tri3 <- summary(pool(with(imp.test_mids, lmer(as.formula(h), weights = combined_weights))),conf.int = T)[11,7]
  upperCI_tri3 <- summary(pool(with(imp.test_mids, lmer(as.formula(h), weights = combined_weights))),conf.int = T)[11,8]
  pval_tri3 <- summary(pool(with(imp.test_mids, lmer(as.formula(h), weights = combined_weights))), conf.int = T)[11,6]
  
  #storing values in empty datafile  
  results_trimesters_rawbrain_unstandardized[x,1] <- bval_tri1
  results_trimesters_rawbrain_unstandardized[x,2] <- seval_tri1
  results_trimesters_rawbrain_unstandardized[x,3] <- lowerCI_tri1
  results_trimesters_rawbrain_unstandardized[x,4] <- upperCI_tri1
  results_trimesters_rawbrain_unstandardized[x,5] <- pval_tri1
  
  results_trimesters_rawbrain_unstandardized[x,6] <- bval_tri2
  results_trimesters_rawbrain_unstandardized[x,7] <- seval_tri2
  results_trimesters_rawbrain_unstandardized[x,8] <- lowerCI_tri2
  results_trimesters_rawbrain_unstandardized[x,9] <- upperCI_tri2
  results_trimesters_rawbrain_unstandardized[x,10] <- pval_tri2
  
  results_trimesters_rawbrain_unstandardized[x,11] <- bval_tri3
  results_trimesters_rawbrain_unstandardized[x,12] <- seval_tri3
  results_trimesters_rawbrain_unstandardized[x,13] <- lowerCI_tri3
  results_trimesters_rawbrain_unstandardized[x,14] <- upperCI_tri3
  results_trimesters_rawbrain_unstandardized[x,15] <- pval_tri3
  
  #assigning names to columns 
  colnames(results_trimesters_rawbrain_unstandardized) <- c("bval_tri1", "seval_tri1", "lowerCI_tri1", "upperCI_tri1", "pval_tri1","bval_tri2", 'seval_tri2', "lowerCI_tri2", "upperCI_tri2", "pval_tri2", "bval_tri3", "seval_tri3", "lowerCI_tri3", "upperCI_tri3", "pval_tri3")
}

write.xlsx(results_trimesters_rawbrain_unstandardized, "results_trimesters_rawbrain_unstandardized.xlsx")

#\

# Sensitivity analysis for fdr sign results in trimester 3 after mutually adjusting for infections in trimester 1 and 2
fdr_outcomes <- c('middletemporal_vol_cortical_standardized', 'parsorbitalis_vol_cortical_standardized', 'rostralanteriorcingulate_vol_cortical_standardized', 'superiorfrontal_vol_cortical_standardized', 'temporalpole_vol_cortical_standardized')

results_sens_analysis_tri3 <- data.frame() #create empty df 

for(x in fdr_outcomes) {
  f <- paste0(x, "~ sumscore_inf_tri3_standardized*timepoint + AGE_M_v2 + ETHNMv2 + EDUCM_3groups+ GSI + SMOKE_ALL + sumscore_inf_tri1_standardized + sumscore_inf_tri2_standardized + (1|IDC)") 
  
  #calculating beta, SE, confidence interval and pval for each model
  bval_tri3 <- summary(pool(with(imp.test_mids, lmer(as.formula(f), weights = combined_weights))),conf.int = T)[13,2]
  seval_tri3 <- summary(pool(with(imp.test_mids, lmer(as.formula(f), weights = combined_weights))),conf.int = T)[13,3]
  lowerCI_tri3 <- summary(pool(with(imp.test_mids, lmer(as.formula(f), weights = combined_weights))),conf.int = T)[13,7]
  upperCI_tri3 <- summary(pool(with(imp.test_mids, lmer(as.formula(f), weights = combined_weights))),conf.int = T)[13,8]
  pval_tri3 <- summary(pool(with(imp.test_mids, lmer(as.formula(f), weights = combined_weights))), conf.int = T)[13,6]
  
  #storing values in empty datafile  
  results_sens_analysis_tri3[x,1] <- bval_tri3
  results_sens_analysis_tri3[x,2] <- seval_tri3
  results_sens_analysis_tri3[x,3] <- lowerCI_tri3
  results_sens_analysis_tri3[x,4] <- upperCI_tri3
  results_sens_analysis_tri3[x,5] <- pval_tri3
  
  #assigning names to columns 
  colnames(results_sens_analysis_tri3) <- c("bval_tri3", "seval_tri3", "lowerCI_tri3", "upperCI_tri3", "pval_tri3")
}

pval1 <- unlist(results_sens_analysis_tri3[, 5]) 
results_sens_analysis_tri3$pval_fdr_adjusted <- p.adjust(pval1, method = 'fdr')
is.num1 <- sapply(results_sens_analysis_tri3, is.numeric)
results_sens_analysis_tri3[is.num1] <- lapply(results_sens_analysis_tri3[is.num1], round, 3)

write.xlsx(results_sens_analysis_tri3, "results_sens_analysis_tri3.xlsx")
#\

###############################################
#-------------------PART 2--------------------#
###############################################

### Extra follow-up analyses or data checks after revisions
setwd('Set_wd_to_path_revisions_results_folder') 

## Add lmer analysis tbv ~ infections (interaction effect total and trimester infections + main effect, and for each effect show standardized and unstandardized beta)
df_tbv_results <- data.frame()

exposure_vars_z <- c('sumscore_inf_tot_standardized', 'sumscore_inf_tri1_standardized', 'sumscore_inf_tri2_standardized', 'sumscore_inf_tri3_standardized')

df_tbv_results_raw <- data.frame()

exposure_vars_raw <- c('sumscore_inf_tot', 'sumscore_inf_tri1', 'sumscore_inf_tri2', 'sumscore_inf_tri3')

for(x in exposure_vars_z) {
  f <- paste0("genr_tbv_standardized ~ ",x,"*timepoint + AGE_M_v2 + ETHNMv2 + EDUCM_3groups+ GSI + SMOKE_ALL  + (1|IDC)") 
  g <-  paste0("genr_tbv_standardized ~ ",x,"+ timepoint + AGE_M_v2 + ETHNMv2 + EDUCM_3groups+ GSI + SMOKE_ALL + (1|IDC)") 
  
  #calculating beta, SE, confidence interval and pval for each model
  bval_int <- summary(pool(with(imp.test_mids, lmer(as.formula(f), weights = combined_weights))),conf.int = T)[11,2]
  seval_int <- summary(pool(with(imp.test_mids, lmer(as.formula(f), weights = combined_weights))),conf.int = T)[11,3]
  lowerCI_int <- summary(pool(with(imp.test_mids, lmer(as.formula(f), weights = combined_weights))),conf.int = T)[11,7]
  upperCI_int <- summary(pool(with(imp.test_mids, lmer(as.formula(f), weights = combined_weights))),conf.int = T)[11,8]
  pval_int <- summary(pool(with(imp.test_mids, lmer(as.formula(f), weights = combined_weights))), conf.int = T)[11,6]
  
  bval_main <- summary(pool(with(imp.test_mids, lmer(as.formula(g), weights = combined_weights))),conf.int = T)[2,2]
  seval_main <- summary(pool(with(imp.test_mids, lmer(as.formula(g), weights = combined_weights))),conf.int = T)[2,3]
  lowerCI_main <- summary(pool(with(imp.test_mids, lmer(as.formula(g), weights = combined_weights))),conf.int = T)[2,7]
  upperCI_main <- summary(pool(with(imp.test_mids, lmer(as.formula(g), weights = combined_weights))),conf.int = T)[2,8]
  pval_main <- summary(pool(with(imp.test_mids, lmer(as.formula(g), weights = combined_weights))), conf.int = T)[2,6]
  
  #storing values in empty datafile  
  df_tbv_results[x,1] <- bval_int
  df_tbv_results[x,2] <- seval_int
  df_tbv_results[x,3] <- lowerCI_int
  df_tbv_results[x,4] <- upperCI_int
  df_tbv_results[x,5] <- pval_int
  df_tbv_results[x,6] <- bval_main
  df_tbv_results[x,7] <- seval_main
  df_tbv_results[x,8] <- lowerCI_main
  df_tbv_results[x,9] <- upperCI_main
  df_tbv_results[x,10] <- pval_main
  
  #assigning names to columns 
  colnames(df_tbv_results) <- c("bval_int", "seval_int", "lowerCI_int", "upperCI_int", "pval_int","bval_main", "seval_main", "lowerCI_main", "upperCI_main", "pval_main")
}

write.xlsx(df_tbv_results, 'tbv_results_z.xlsx')

for(x in exposure_vars_raw) {
  f <- paste0("genr_tbv ~ ",x,"*timepoint + AGE_M_v2 + ETHNMv2 + EDUCM_3groups+ GSI + SMOKE_ALL + (1|IDC)") 
  g <-  paste0("genr_tbv ~ ",x,"+ timepoint + AGE_M_v2 + ETHNMv2 + EDUCM_3groups+ GSI + SMOKE_ALL + (1|IDC)") 
  
  #calculating beta, SE, confidence interval and pval for each model
  bval_int <- summary(pool(with(imp.test_mids, lmer(as.formula(f), weights = combined_weights))),conf.int = T)[11,2]
  bval_main <- summary(pool(with(imp.test_mids, lmer(as.formula(g), weights = combined_weights))),conf.int = T)[2,2]

  #storing values in empty datafile  
  df_tbv_results_raw[x,1] <- bval_int
  df_tbv_results_raw[x,2] <- bval_main
  
  #assigning names to columns 
  colnames(df_tbv_results_raw) <- c("bval_int","bval_main")
}

write.xlsx(df_tbv_results_raw, 'tbv_results_raw_b.xlsx')

#\

## Rerun main analysis with interaction for child sex
results_sex_int <- data.frame()

for(x in outcome_vars_ztrans) {
  f <- paste0(x, "~ sumscore_inf_tot_standardized*timepoint*GENDER + AGE_M_v2 + ETHNMv2 + EDUCM_3groups+ GSI + SMOKE_ALL + (1|IDC)") 
  
  #calculating beta, SE, confidence interval and pval for each model
  bval_int_no_icv <- summary(pool(with(imp.test_mids, lmer(as.formula(f), weights = combined_weights))),conf.int = T)[15,2]
  seval_int_no_icv <- summary(pool(with(imp.test_mids, lmer(as.formula(f), weights = combined_weights))),conf.int = T)[15,3]
  lowerCI_int_no_icv <- summary(pool(with(imp.test_mids, lmer(as.formula(f), weights = combined_weights))),conf.int = T)[15,7]
  upperCI_int_no_icv <- summary(pool(with(imp.test_mids, lmer(as.formula(f), weights = combined_weights))),conf.int = T)[15,8]
  pval_int_no_icv <- summary(pool(with(imp.test_mids, lmer(as.formula(f), weights = combined_weights))), conf.int = T)[15,6]
  
  #storing values in empty datafile  
  results_sex_int[x,1] <- bval_int_no_icv
  results_sex_int[x,2] <- seval_int_no_icv
  results_sex_int[x,3] <- lowerCI_int_no_icv
  results_sex_int[x,4] <- upperCI_int_no_icv
  results_sex_int[x,5] <- pval_int_no_icv
  
  #assigning names to columns 
  colnames(results_sex_int) <- c("bval_int_no_icv", "seval_int_no_icv", "lowerCI_int_no_icv", "upperCI_int_no_icv", "pval_int_no_icv")
}

write.xlsx(results_sex_int, "results_sex_int.xlsx") 

## Rerun analysis for 5 sign regions after adding SES as interaction term
sign_outcomes <- c('middletemporal_vol_cortical_standardized', 'parsorbitalis_vol_cortical_standardized', 'rostralanteriorcingulate_vol_cortical.x_standardized', 'superiorfrontal_vol_cortical.x_standardized', 'temporalpole_vol_cortical.x_standardized')

ses_interaction <- data.frame()

for(x in sign_outcomes) {
  f <- paste0(x, "~ sumscore_inf_tri3_standardized*timepoint*EDUCM_3groups + GENDER + AGE_M_v2 + ETHNMv2 + GSI + SMOKE_ALL + (1|IDC)") 
  
  #calculating beta, SE, confidence interval and pval for each model
  bval <- summary(pool(with(imp.test_mids, lmer(as.formula(f), weights = combined_weights))),conf.int = T)[15,2]
  lowerCI <- summary(pool(with(imp.test_mids, lmer(as.formula(f), weights = combined_weights))),conf.int = T)[15,7]
  upperCI <- summary(pool(with(imp.test_mids, lmer(as.formula(f), weights = combined_weights))),conf.int = T)[15,8]
  pval <- summary(pool(with(imp.test_mids, lmer(as.formula(f), weights = combined_weights))), conf.int = T)[15,6]
  
  #storing values in empty datafile  
  ses_interaction[x,1] <- bval
  ses_interaction[x,2] <- lowerCI
  ses_interaction[x,3] <- upperCI
  ses_interaction[x,4] <- pval
  
  #assigning names to columns 
  colnames(ses_interaction) <- c("bval", 'lowerCI', 'upperCI', 'pval')
}

write.xlsx(ses_interaction, "results_ses_int.xlsx")

## Rerun analysis of 5 sign regions after adjusting for immune conditions and medication mother
immune_confounding_results <- data.frame()

for(x in sign_outcomes) {
  f <- paste0(x, "~ sumscore_inf_tri3_standardized*timepoint + EDUCM_3groups + GENDER + AGE_M_v2 + ETHNMv2 + GSI + SMOKE_ALL  + DIAB_GRA + PIH_v1 + preeclampsia + (1|IDC)") 
  
  #calculating beta, SE, confidence interval and pval for each model
  bval <- summary(pool(with(imp.test_mids, lmer(as.formula(f), weights = combined_weights))),conf.int = T)[15,2]
  lowerCI <- summary(pool(with(imp.test_mids, lmer(as.formula(f), weights = combined_weights))),conf.int = T)[15,7]
  upperCI <- summary(pool(with(imp.test_mids, lmer(as.formula(f), weights = combined_weights))),conf.int = T)[15,8]
  pval <- summary(pool(with(imp.test_mids, lmer(as.formula(f), weights = combined_weights))), conf.int = T)[15,6]
  
  #storing values in empty datafile  
  immune_confounding_results[x,1] <- bval
  immune_confounding_results[x,2] <- lowerCI
  immune_confounding_results[x,3] <- upperCI
  immune_confounding_results[x,4] <- pval
  
  #assigning names to columns 
  colnames(immune_confounding_results) <- c("bval", 'lowerCI', 'upperCI', 'pval')
}

write.xlsx(immune_confounding_results, "immune_confounding_results.xlsx")

## Rerun analysis of 5 sign regions with age instead of time
age_interaction <- data.frame()

for(x in sign_outcomes) {
  f <- paste0(x, "~ sumscore_inf_tri3_standardized*age_child_mri + EDUCM_3groups + GENDER + AGE_M_v2 + ETHNMv2 + GSI + SMOKE_ALL + (1|IDC)") 
  
  #calculating beta, SE, confidence interval and pval for each model
  bval <- summary(pool(with(imp.test_mids, lmer(as.formula(f), weights = combined_weights))),conf.int = T)[12,2]
  lowerCI <- summary(pool(with(imp.test_mids, lmer(as.formula(f), weights = combined_weights))),conf.int = T)[12,7]
  upperCI <- summary(pool(with(imp.test_mids, lmer(as.formula(f), weights = combined_weights))),conf.int = T)[12,8]
  pval <- summary(pool(with(imp.test_mids, lmer(as.formula(f), weights = combined_weights))), conf.int = T)[12,6]
  
  #storing values in empty datafile  
  age_interaction[x,1] <- bval
  age_interaction[x,2] <- lowerCI
  age_interaction[x,3] <- upperCI
  age_interaction[x,4] <- pval
  
  #assigning names to columns 
  colnames(age_interaction) <- c("bval", 'lowerCI', 'upperCI', 'pval')
}

pval1 <- unlist(age_interaction[, 4]) 
age_interaction$pval_fdr_adjusted <- p.adjust(pval1, method = 'fdr')
is.num1 <- sapply(age_interaction, is.numeric)
age_interaction[is.num1] <- lapply(age_interaction[is.num1], round, 3)

write.xlsx(age_interaction, "age_interaction.xlsx")

## Check if non linear spline improves model fit for 5 sign regions 
single_df <- complete(imp.test_mids, 30)

lm1a <- lmer(middletemporal_vol_cortical_standardized ~ sumscore_inf_tri3_standardized*timepoint + EDUCM_3groups + GENDER + AGE_M_v2 + ETHNMv2 + GSI + SMOKE_ALL + (1|IDC), weights = combined_weights, data = single_df)
lm1b <- lmer(middletemporal_vol_cortical_standardized ~ ns(sumscore_inf_tri3_standardized*timepoint,3) + EDUCM_3groups + GENDER + AGE_M_v2 + ETHNMv2 + GSI + SMOKE_ALL + (1|IDC), weights = combined_weights, data = single_df)
anova(lm1a, lm1b) #sign 

lm2a <- lmer(parsorbitalis_vol_cortical_standardized ~ sumscore_inf_tri3_standardized*timepoint + EDUCM_3groups + GENDER + AGE_M_v2 + ETHNMv2 + GSI + SMOKE_ALL + (1|IDC), weights = combined_weights, data = single_df)
lm2b <- lmer(parsorbitalis_vol_cortical_standardized ~ ns(sumscore_inf_tri3_standardized*timepoint,3) + EDUCM_3groups + GENDER + AGE_M_v2 + ETHNMv2 + GSI + SMOKE_ALL + (1|IDC), weights = combined_weights, data = single_df)
anova(lm2a, lm2b) #sign 

lm3a <- lmer(rostralanteriorcingulate_vol_cortical.x_standardized ~ sumscore_inf_tri3_standardized*timepoint + EDUCM_3groups + GENDER + AGE_M_v2 + ETHNMv2 + GSI + SMOKE_ALL + (1|IDC), weights = combined_weights, data = single_df)
lm3b <- lmer(rostralanteriorcingulate_vol_cortical.x_standardized ~ ns(sumscore_inf_tri3_standardized*timepoint,3) + EDUCM_3groups + GENDER + AGE_M_v2 + ETHNMv2 + GSI + SMOKE_ALL + (1|IDC), weights = combined_weights, data = single_df)
anova(lm3a, lm3b) #sign 

lm4a <- lmer(superiorfrontal_vol_cortical.x_standardized ~ sumscore_inf_tri3_standardized*timepoint + EDUCM_3groups + GENDER + AGE_M_v2 + ETHNMv2 + GSI + SMOKE_ALL + (1|IDC), weights = combined_weights, data = single_df)
lm4b <- lmer(superiorfrontal_vol_cortical.x_standardized ~ ns(sumscore_inf_tri3_standardized*timepoint,3) + EDUCM_3groups + GENDER + AGE_M_v2 + ETHNMv2 + GSI + SMOKE_ALL + (1|IDC), weights = combined_weights, data = single_df)
anova(lm4a, lm4b) #sign 

lm5a <- lmer(temporalpole_vol_cortical.x_standardized ~ sumscore_inf_tri3_standardized*timepoint + EDUCM_3groups + GENDER + AGE_M_v2 + ETHNMv2 + GSI + SMOKE_ALL + (1|IDC), weights = combined_weights, data = single_df)
lm5b <- lmer(temporalpole_vol_cortical.x_standardized ~ ns(sumscore_inf_tri3_standardized*timepoint,3) + EDUCM_3groups + GENDER + AGE_M_v2 + ETHNMv2 + GSI + SMOKE_ALL + (1|IDC), weights = combined_weights, data = single_df)
anova(lm5a, lm5b) #sign 

# run non linear models with splines for sum stats 
lm1b <- with(imp.test_mids, lmer(middletemporal_vol_cortical_standardized ~ ns(sumscore_inf_tri3_standardized*timepoint,3) + EDUCM_3groups + GENDER + AGE_M_v2 + ETHNMv2 + GSI + SMOKE_ALL + (1|IDC), weights = combined_weights), conf.int = T)
summary(pool(lm1b), conf.int = T)

lm2b <- with(imp.test_mids, lmer(parsorbitalis_vol_cortical_standardized ~ ns(sumscore_inf_tri3_standardized*timepoint,3) + EDUCM_3groups + GENDER + AGE_M_v2 + ETHNMv2 + GSI + SMOKE_ALL + (1|IDC), weights = combined_weights), conf.int = T)
summary(pool(lm2b), conf.int = T)

lm3b <- with(imp.test_mids, lmer(rostralanteriorcingulate_vol_cortical.x_standardized ~ ns(sumscore_inf_tri3_standardized*timepoint,3) + EDUCM_3groups + GENDER + AGE_M_v2 + ETHNMv2 + GSI + SMOKE_ALL + (1|IDC), weights = combined_weights), conf.int = T)
summary(pool(lm3b), conf.int = T)

lm4b <- with(imp.test_mids, lmer(superiorfrontal_vol_cortical.x_standardized ~ ns(sumscore_inf_tri3_standardized*timepoint,3) + EDUCM_3groups + GENDER + AGE_M_v2 + ETHNMv2 + GSI + SMOKE_ALL + (1|IDC), weights = combined_weights), conf.int = T)
summary(pool(lm4b), conf.int = T)

lm5b <- with(imp.test_mids, lmer(temporalpole_vol_cortical.x_standardized ~ ns(sumscore_inf_tri3_standardized*timepoint,3) + EDUCM_3groups + GENDER + AGE_M_v2 + ETHNMv2 + GSI + SMOKE_ALL + (1|IDC), weights = combined_weights), conf.int = T)
summary(pool(lm5b), conf.int = T)

# plot one of the nonlinear trajectories to better understand results 
library(sjPlot)

fit <- lmer(middletemporal_vol_cortical_standardized ~ sumscore_inf_tri3_standardized*timepoint + EDUCM_3groups + GENDER + AGE_M_v2 + ETHNMv2 + GSI + SMOKE_ALL + (1|IDC), weights = combined_weights, data = single_df)
ns_fit <- lmer(middletemporal_vol_cortical_standardized ~ ns(sumscore_inf_tri3_standardized*timepoint,3) + EDUCM_3groups + GENDER + AGE_M_v2 + ETHNMv2 + GSI + SMOKE_ALL + (1|IDC), weights = combined_weights, data = single_df)
plot_model(fit, type = "pred", terms = c("timepoint","sumscore_inf_tri3_standardized[-0.97198, 0.04738, 5.14415]"))

ns_fit <- lmer(middletemporal_vol_cortical_standardized ~ ns(sumscore_inf_tri3_standardized*timepoint,3) + EDUCM_3groups + GENDER + AGE_M_v2 + ETHNMv2 + GSI + SMOKE_ALL + (1|IDC), weights = combined_weights, data = single_df)
plot_model(ns_fit, type = "pred", terms = c("timepoint[all]","sumscore_inf_tri3_standardized[-0.97198, 0.04738, 5.14415]"))

## Rerun analysis of 5 sign regions with moderation of maternal psychosocial factors, child trauma, maternal alcohol use (second hits paper moderators)
secondhits_interaction <- data.frame()

for(x in sign_outcomes) {
  f <- paste0(x, "~ sumscore_inf_tri3_standardized*timepoint*post_life_events + EDUCM_3groups + GENDER + AGE_M_v2 + ETHNMv2 + GSI + SMOKE_ALL + (1|IDC)") 
  g <- paste0(x, "~ sumscore_inf_tri3_standardized*timepoint*post_direct_victimization + EDUCM_3groups + GENDER + AGE_M_v2 + ETHNMv2 + GSI + SMOKE_ALL + (1|IDC)") 
  h <- paste0(x, "~ sumscore_inf_tri3_standardized*timepoint*GSI + EDUCM_3groups + GENDER + AGE_M_v2 + ETHNMv2 + SMOKE_ALL + (1|IDC)") 
  i <- paste0(x, "~ sumscore_inf_tri3_standardized*timepoint*mdrink_updated + EDUCM_3groups + GENDER + AGE_M_v2 + ETHNMv2 + SMOKE_ALL + (1|IDC)") 
  
  #calculating beta, SE, confidence interval and pval for each model
  bval_ple <- summary(pool(with(imp.test_mids, lmer(as.formula(f), weights = combined_weights))),conf.int = T)[16,2]
  lowerCI_ple <- summary(pool(with(imp.test_mids, lmer(as.formula(f), weights = combined_weights))),conf.int = T)[16,7]
  upperCI_ple <- summary(pool(with(imp.test_mids, lmer(as.formula(f), weights = combined_weights))),conf.int = T)[16,8]
  pval_ple <- summary(pool(with(imp.test_mids, lmer(as.formula(f), weights = combined_weights))), conf.int = T)[16,6]
  
  #calculating beta, SE, confidence interval and pval for each model
  bval_pdv <- summary(pool(with(imp.test_mids, lmer(as.formula(g), weights = combined_weights))),conf.int = T)[16,2]
  lowerCI_pdv <- summary(pool(with(imp.test_mids, lmer(as.formula(g), weights = combined_weights))),conf.int = T)[16,7]
  upperCI_pdv <- summary(pool(with(imp.test_mids, lmer(as.formula(g), weights = combined_weights))),conf.int = T)[16,8]
  pval_pdv <- summary(pool(with(imp.test_mids, lmer(as.formula(g), weights = combined_weights))), conf.int = T)[16,6]
  
  #calculating beta, SE, confidence interval and pval for each model
  bval_gsi <- summary(pool(with(imp.test_mids, lmer(as.formula(h), weights = combined_weights))),conf.int = T)[15,2]
  lowerCI_gsi <- summary(pool(with(imp.test_mids, lmer(as.formula(h), weights = combined_weights))),conf.int = T)[15,7]
  upperCI_gsi <- summary(pool(with(imp.test_mids, lmer(as.formula(h), weights = combined_weights))),conf.int = T)[15,8]
  pval_gsi <- summary(pool(with(imp.test_mids, lmer(as.formula(h), weights = combined_weights))), conf.int = T)[15,6]
  
  #calculating beta, SE, confidence interval and pval for each model
  bval_alc <- summary(pool(with(imp.test_mids, lmer(as.formula(i), weights = combined_weights))),conf.int = T)[16,2]
  lowerCI_alc <- summary(pool(with(imp.test_mids, lmer(as.formula(i), weights = combined_weights))),conf.int = T)[16,7]
  upperCI_alc <- summary(pool(with(imp.test_mids, lmer(as.formula(i), weights = combined_weights))),conf.int = T)[16,8]
  pval_alc <- summary(pool(with(imp.test_mids, lmer(as.formula(i), weights = combined_weights))), conf.int = T)[16,6]
  
  #storing values in empty datafile  
  secondhits_interaction[x,1] <- bval_ple
  secondhits_interaction[x,2] <- lowerCI_ple
  secondhits_interaction[x,3] <- upperCI_ple
  secondhits_interaction[x,4] <- pval_ple
  
  secondhits_interaction[x,5] <- bval_pdv
  secondhits_interaction[x,6] <- lowerCI_pdv
  secondhits_interaction[x,7] <- upperCI_pdv
  secondhits_interaction[x,8] <- pval_pdv
  
  secondhits_interaction[x,9] <- bval_gsi
  secondhits_interaction[x,10] <- lowerCI_gsi
  secondhits_interaction[x,11] <- upperCI_gsi
  secondhits_interaction[x,12] <- pval_gsi
  
  secondhits_interaction[x,13] <- bval_alc
  secondhits_interaction[x,14] <- lowerCI_alc
  secondhits_interaction[x,15] <- upperCI_alc
  secondhits_interaction[x,16] <- pval_alc
  
  #assigning names to columns 
  colnames(secondhits_interaction) <- c("bval_ple", 'lowerCI_ple', 'upperCI_ple', 'pval_ple', "bval_pdv", 'lowerCI_pdv', 'upperCI_pdv', 'pval_pdv', "bval_gsi", 'lowerCI_gsi', 'upperCI_gsi', 'pval_gsi', "bval_alc", 'lowerCI_alc', 'upperCI_alc', 'pval_alc')
}

write.xlsx(secondhits_interaction, "secondhits_interaction.xlsx")

#\

## Sensitivity analysis for fdr sign results in trimester 3 after mutually adjusting for infections in trimester 1 and 2 (but now additionally interaction for time for also trimester 1 and 2)
fdr_outcomes <- c('middletemporal_vol_cortical_standardized', 'parsorbitalis_vol_cortical_standardized', 'rostralanteriorcingulate_vol_cortical.x_standardized', 'superiorfrontal_vol_cortical.x_standardized', 'temporalpole_vol_cortical.x_standardized')

results_sens_analysis_tri3_corrected <- data.frame() #create empty df 

for(x in fdr_outcomes) {
  f <- paste0(x, "~ timepoint*(sumscore_inf_tri1_standardized + sumscore_inf_tri2_standardized + sumscore_inf_tri3_standardized) + AGE_M_v2 + ETHNMv2 + EDUCM_3groups+ GSI + SMOKE_ALL +  (1|IDC)") 
  #calculating beta, SE, confidence interval and pval for each model
  bval_tri3 <- summary(pool(with(imp.test_mids, lmer(as.formula(f), weights = combined_weights))),conf.int = T)[15,2]
  seval_tri3 <- summary(pool(with(imp.test_mids, lmer(as.formula(f), weights = combined_weights))),conf.int = T)[15,3]
  lowerCI_tri3 <- summary(pool(with(imp.test_mids, lmer(as.formula(f), weights = combined_weights))),conf.int = T)[15,7]
  upperCI_tri3 <- summary(pool(with(imp.test_mids, lmer(as.formula(f), weights = combined_weights))),conf.int = T)[15,8]
  pval_tri3 <- summary(pool(with(imp.test_mids, lmer(as.formula(f), weights = combined_weights))), conf.int = T)[15,6]
  
  #storing values in empty datafile  
  results_sens_analysis_tri3_corrected[x,1] <- bval_tri3
  results_sens_analysis_tri3_corrected[x,2] <- seval_tri3
  results_sens_analysis_tri3_corrected[x,3] <- lowerCI_tri3
  results_sens_analysis_tri3_corrected[x,4] <- upperCI_tri3
  results_sens_analysis_tri3_corrected[x,5] <- pval_tri3
  
  #assigning names to columns 
  colnames(results_sens_analysis_tri3_corrected) <- c("bval_tri3", "seval_tri3", "lowerCI_tri3", "upperCI_tri3", "pval_tri3")
}

pval1 <- unlist(results_sens_analysis_tri3_corrected[, 5]) 
results_sens_analysis_tri3_corrected$pval_fdr_adjusted <- p.adjust(pval1, method = 'fdr')
is.num1 <- sapply(results_sens_analysis_tri3_corrected, is.numeric)
results_sens_analysis_tri3_corrected[is.num1] <- lapply(results_sens_analysis_tri3_corrected[is.num1], round, 3)

write.xlsx(results_sens_analysis_tri3_corrected, "results_sens_analysis_tri3_corrected.xlsx")

#\

## Run analysis with SES as exposure and brain outcoems 
results_ses <- data.frame() #create empty df 

for(x in fdr_outcomes) {
  f <- paste0(x, "~ EDUCM_3groups*timepoint + AGE_M_v2 + ETHNMv2 + GSI + SMOKE_ALL +  (1|IDC)")
  
  #calculating beta, SE, confidence interval and pval for each model
  bval <- summary(pool(with(imp.test_mids, lmer(as.formula(f), weights = combined_weights))),conf.int = T)[11,2]
  seval <- summary(pool(with(imp.test_mids, lmer(as.formula(f), weights = combined_weights))),conf.int = T)[11,3]
  lowerCI <- summary(pool(with(imp.test_mids, lmer(as.formula(f), weights = combined_weights))),conf.int = T)[11,7]
  upperCI <- summary(pool(with(imp.test_mids, lmer(as.formula(f), weights = combined_weights))),conf.int = T)[11,8]
  pval <- summary(pool(with(imp.test_mids, lmer(as.formula(f), weights = combined_weights))), conf.int = T)[11,6]
  
  #storing values in empty datafile  
  results_ses[x,1] <- bval
  results_ses[x,2] <- seval
  results_ses[x,3] <- lowerCI
  results_ses[x,4] <- upperCI
  results_ses[x,5] <- pval

  #assigning names to columns 
  colnames(results_ses) <- c("bval", "seval", "lowerCI", "upperCI", "pval")
}

pval1 <- unlist(results_ses[, 5]) 
results_ses$pval_fdr_adjusted <- p.adjust(pval1, method = 'fdr')
is.num1 <- sapply(results_ses, is.numeric)
results_ses[is.num1] <- lapply(results_ses[is.num1], round, 3)

#\

## Run individual analysis between infections & 5 fdr significant regions for each time point individually to understand the longitudinal patterns better (i.e., is it catch-up growth driven by changes at age 6?) 
df_age6 <- subset(imp.test_long, timepoint == 1)
df_age10 <- subset(imp.test_long, timepoint == 2)
df_age14 <- subset(imp.test_long, timepoint == 3)

df_age6_mids <- as.mids(df_age6)
df_age10_mids <- as.mids(df_age10)
df_age14_mids <- as.mids(df_age14)

fdr_outcomes <- c('middletemporal_vol_cortical_standardized', 'parsorbitalis_vol_cortical_standardized', 'rostralanteriorcingulate_vol_cortical.x_standardized', 'superiorfrontal_vol_cortical.x_standardized', 'temporalpole_vol_cortical.x_standardized')

# age 6 lm 
cross_sectional_results_f6 <- data.frame() #create empty df 

for(x in fdr_outcomes) {
  model <- paste0(x, "~ sumscore_inf_tri3_standardized + age_child_mri + EDUCM_3groups + GENDER + AGE_M_v2 + ETHNMv2 + GSI + SMOKE_ALL") 
  
  #calculating beta, SE, confidence interval and pval
  bval <- summary(pool(with(df_age6_mids, lm(as.formula(model), weights = combined_weights))),conf.int = T)[2,2]
  seval <- summary(pool(with(df_age6_mids, lm(as.formula(model), weights = combined_weights))),conf.int = T)[2,3]
  lowerCI <- summary(pool(with(df_age6_mids, lm(as.formula(model), weights = combined_weights))),conf.int = T)[2,7]
  upperCI <- summary(pool(with(df_age6_mids, lm(as.formula(model), weights = combined_weights))),conf.int = T)[2,8]
  pval <- summary(pool(with(df_age6_mids, lm(as.formula(model), weights = combined_weights))), conf.int = T)[2,6]
  
  #storing values in empty datafile  
  cross_sectional_results_f6[x,1] <- bval
  cross_sectional_results_f6[x,2] <- seval
  cross_sectional_results_f6[x,3] <- lowerCI
  cross_sectional_results_f6[x,4] <- upperCI
  cross_sectional_results_f6[x,5] <- pval
  
  #assigning names to columns 
  colnames(cross_sectional_results_f6) <- c("bval", "seval", "lowerCI", "upperCI", "pval")
}

pval1 <- unlist(cross_sectional_results_f6[, 5]) 
cross_sectional_results_f6$pval_fdr_adjusted <- p.adjust(pval1, method = 'fdr')
is.num1 <- sapply(cross_sectional_results_f6, is.numeric)
cross_sectional_results_f6[is.num1] <- lapply(cross_sectional_results_f6[is.num1], round, 3)

# age 10 lm 
cross_sectional_results_f10 <- data.frame() #create empty df 

for(x in fdr_outcomes) {
  model <- paste0(x, "~ sumscore_inf_tri3_standardized + age_child_mri + EDUCM_3groups + GENDER + AGE_M_v2 + ETHNMv2 + GSI + SMOKE_ALL") 
  
  #calculating beta, SE, confidence interval and pval
  bval <- summary(pool(with(df_age10_mids, lm(as.formula(model), weights = combined_weights))),conf.int = T)[2,2]
  seval <- summary(pool(with(df_age10_mids, lm(as.formula(model), weights = combined_weights))),conf.int = T)[2,3]
  lowerCI <- summary(pool(with(df_age10_mids, lm(as.formula(model), weights = combined_weights))),conf.int = T)[2,7]
  upperCI <- summary(pool(with(df_age10_mids, lm(as.formula(model), weights = combined_weights))),conf.int = T)[2,8]
  pval <- summary(pool(with(df_age10_mids, lm(as.formula(model), weights = combined_weights))), conf.int = T)[2,6]
  
  #storing values in empty datafile  
  cross_sectional_results_f10[x,1] <- bval
  cross_sectional_results_f10[x,2] <- seval
  cross_sectional_results_f10[x,3] <- lowerCI
  cross_sectional_results_f10[x,4] <- upperCI
  cross_sectional_results_f10[x,5] <- pval
  
  #assigning names to columns 
  colnames(cross_sectional_results_f10) <- c("bval", "seval", "lowerCI", "upperCI", "pval")
}

pval1 <- unlist(cross_sectional_results_f10[, 5]) 
cross_sectional_results_f10$pval_fdr_adjusted <- p.adjust(pval1, method = 'fdr')
is.num1 <- sapply(cross_sectional_results_f10, is.numeric)
cross_sectional_results_f10[is.num1] <- lapply(cross_sectional_results_f10[is.num1], round, 3)

# age 14 lm 
cross_sectional_results_f14 <- data.frame() #create empty df 

for(x in fdr_outcomes) {
  model <- paste0(x, "~ sumscore_inf_tri3_standardized + age_child_mri + EDUCM_3groups + GENDER + AGE_M_v2 + ETHNMv2 + GSI + SMOKE_ALL") 
  
  #calculating beta, SE, confidence interval and pval
  bval <- summary(pool(with(df_age14_mids, lm(as.formula(model), weights = combined_weights))),conf.int = T)[2,2]
  seval <- summary(pool(with(df_age14_mids, lm(as.formula(model), weights = combined_weights))),conf.int = T)[2,3]
  lowerCI <- summary(pool(with(df_age14_mids, lm(as.formula(model), weights = combined_weights))),conf.int = T)[2,7]
  upperCI <- summary(pool(with(df_age14_mids, lm(as.formula(model), weights = combined_weights))),conf.int = T)[2,8]
  pval <- summary(pool(with(df_age14_mids, lm(as.formula(model), weights = combined_weights))), conf.int = T)[2,6]
  
  #storing values in empty datafile  
  cross_sectional_results_f14[x,1] <- bval
  cross_sectional_results_f14[x,2] <- seval
  cross_sectional_results_f14[x,3] <- lowerCI
  cross_sectional_results_f14[x,4] <- upperCI
  cross_sectional_results_f14[x,5] <- pval
  
  #assigning names to columns 
  colnames(cross_sectional_results_f14) <- c("bval", "seval", "lowerCI", "upperCI", "pval")
}

pval1 <- unlist(cross_sectional_results_f14[, 5]) 
cross_sectional_results_f14$pval_fdr_adjusted <- p.adjust(pval1, method = 'fdr')
is.num1 <- sapply(cross_sectional_results_f14, is.numeric)
cross_sectional_results_f14[is.num1] <- lapply(cross_sectional_results_f14[is.num1], round, 3)

# create correlation plots for age 6, 10 and 14
df_age6_nonimp <- complete(df_age6_mids, 30)
df_age10_nonimp <- complete(df_age10_mids, 30)
df_age14_nonimp <- complete(df_age14_mids, 30)

library(corrplot)
library(RColorBrewer)

corr_vars_f6 <- dplyr::select(df_age6_nonimp, 'sumscore_inf_tri3', 'middletemporal_vol_cortical', 'parsorbitalis_vol_cortical', 'rostralanteriorcingulate_vol_cortical.x', 'superiorfrontal_vol_cortical.x', 'temporalpole_vol_cortical.x')
corr_vars_f6 <- rename(corr_vars_f6, 'Infection sum score (trimester 3)' = sumscore_inf_tri3, 'Middle temporal' = middletemporal_vol_cortical, 'Pars orbitalis' = parsorbitalis_vol_cortical, 'Rostral anterior cingulate' = rostralanteriorcingulate_vol_cortical.x, 'Superior frontal' = superiorfrontal_vol_cortical.x, 'Temporal pole' = temporalpole_vol_cortical.x)
corr_vars_f10 <- dplyr::select(df_age10_nonimp, 'sumscore_inf_tri3', 'middletemporal_vol_cortical', 'parsorbitalis_vol_cortical', 'rostralanteriorcingulate_vol_cortical.x', 'superiorfrontal_vol_cortical.x', 'temporalpole_vol_cortical.x')
corr_vars_f10 <- rename(corr_vars_f10, 'Infection sum score (trimester 3)' = sumscore_inf_tri3, 'Middle temporal' = middletemporal_vol_cortical, 'Pars orbitalis' = parsorbitalis_vol_cortical, 'Rostral anterior cingulate' = rostralanteriorcingulate_vol_cortical.x, 'Superior frontal' = superiorfrontal_vol_cortical.x, 'Temporal pole' = temporalpole_vol_cortical.x)
corr_vars_f14 <- dplyr::select(df_age14_nonimp, 'sumscore_inf_tri3', 'middletemporal_vol_cortical', 'parsorbitalis_vol_cortical', 'rostralanteriorcingulate_vol_cortical.x', 'superiorfrontal_vol_cortical.x', 'temporalpole_vol_cortical.x')
corr_vars_f14 <- rename(corr_vars_f14, 'Infection sum score (trimester 3)' = sumscore_inf_tri3, 'Middle temporal' = middletemporal_vol_cortical, 'Pars orbitalis' = parsorbitalis_vol_cortical, 'Rostral anterior cingulate' = rostralanteriorcingulate_vol_cortical.x, 'Superior frontal' = superiorfrontal_vol_cortical.x, 'Temporal pole' = temporalpole_vol_cortical.x)

correlation_f6 <- cor(corr_vars_f6, use="pairwise.complete.obs")
correlation_f10 <- cor(corr_vars_f10, use="pairwise.complete.obs")
correlation_f14 <- cor(corr_vars_f14, use="pairwise.complete.obs")

# Set up a 3x2 grid for plots
par(mfrow = c(3, 2))

# Plot 1: Neuroimaging visit 1
corrplot::corrplot(
  correlation_f6,
  method = 'number',
  order = 'FPC',
  type = 'lower',
  diag = FALSE,
  tl.col = 'black',
  tl.cex = 0.7,
  col = brewer.pal(n = 8, name = "RdBu"),
  tl.srt = 45
)
title("Neuroimaging visit 1")

# Plot 2: Neuroimaging visit 2
corrplot::corrplot(
  correlation_f10,
  method = 'number',
  order = 'FPC',
  type = 'lower',
  diag = FALSE,
  tl.col = 'black',
  tl.cex = 0.7,
  col = brewer.pal(n = 8, name = "RdBu"),
  tl.srt = 45
)
title("Neuroimaging visit 2")

# Plot 3: Neuroimaging visit 3
corrplot::corrplot(
  correlation_f14,
  method = 'number',
  order = 'FPC',
  type = 'lower',
  diag = FALSE,
  tl.col = 'black',
  tl.cex = 0.7,
  col = brewer.pal(n = 8, name = "RdBu"),
  tl.srt = 45
)
title("Neuroimaging visit 3")

# Compare infection vs no infection group at baseline ; so visit 1 (violin plots)
df_age6_nonimp$sumscore_inf_tri3_cat2 <- ifelse(df_age6_nonimp$sumscore_inf_tri3 > 1, 'Infection exposure', 'No infection exposure')

a <- ggplot(df_age6_nonimp, aes(x=sumscore_inf_tri3_cat2, y = middletemporal_vol_cortical)) + geom_violin(fill='thistle3') + theme_minimal() + stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange", color="maroon") + labs(x = 'Prenatal infection sum score in trimester 3', y = 'Middle temporal volume at T1')

b <- ggplot(df_age6_nonimp, aes(x=sumscore_inf_tri3_cat2, y = parsorbitalis_vol_cortical)) + geom_violin(fill='thistle3') + theme_minimal() + stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange", color="maroon") + labs(x = 'Prenatal infection sum score in trimester 3', y = 'Pars orbitalis volume at T1')

c <- ggplot(df_age6_nonimp, aes(x=sumscore_inf_tri3_cat2, y = rostralanteriorcingulate_vol_cortical.x)) + geom_violin(fill='thistle3') + theme_minimal() + stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange", color="maroon") + labs(x = 'Prenatal infection sum score in trimester 3', y = 'Rostral anterior cingulate volume at T1')

d <- ggplot(df_age6_nonimp, aes(x=sumscore_inf_tri3_cat2, y = superiorfrontal_vol_cortical.x)) + geom_violin(fill='thistle3') + theme_minimal() + stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange", color="maroon") + labs(x = 'Prenatal infection sum score in trimester 3', y = 'Superior frontal volume at T1')

e <- ggplot(df_age6_nonimp, aes(x=sumscore_inf_tri3_cat2, y = temporalpole_vol_cortical.x)) + geom_violin(fill='thistle3') + theme_minimal() + stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange", color="maroon") + labs(x = 'Prenatal infection sum score in trimester 3', y = 'Temporal pole volume at T1')

library(ggpubr)
ggarrange(a, b, c, d, e, labels = 'AUTO', ncol = 3, nrow=2)

# Compare infection vs no infection group at baseline ; so visit 1 (boxplots)
a <- ggplot(df_age6_nonimp, aes(x = sumscore_inf_tri3_cat2, y = middletemporal_vol_cortical)) +
  geom_boxplot(fill = 'lightcyan3', color = 'black') +  
  theme_minimal() +
  stat_summary(fun.data = mean_sdl, mult = 1, geom = "pointrange", color = "maroon") +
  labs(x = 'Prenatal infection sum score in trimester 3', y = 'Middle temporal volume at T1')

b <- ggplot(df_age6_nonimp, aes(x = sumscore_inf_tri3_cat2, y = parsorbitalis_vol_cortical)) +
  geom_boxplot(fill = 'lightcyan3', color = 'black') +  
  theme_minimal() +
  stat_summary(fun.data = mean_sdl, mult = 1, geom = "pointrange", color = "maroon") +
  labs(x = 'Prenatal infection sum score in trimester 3', y = 'Pars orbitalis volume at T1')

c <- ggplot(df_age6_nonimp, aes(x = sumscore_inf_tri3_cat2, y = rostralanteriorcingulate_vol_cortical.x)) +
  geom_boxplot(fill = 'lightcyan3', color = 'black') +  
  theme_minimal() +
  stat_summary(fun.data = mean_sdl, mult = 1, geom = "pointrange", color = "maroon") +
  labs(x = 'Prenatal infection sum score in trimester 3', y = 'Rostral anterior cingulate volume at T1')

d <- ggplot(df_age6_nonimp, aes(x = sumscore_inf_tri3_cat2, y = superiorfrontal_vol_cortical.x)) +
  geom_boxplot(fill = 'lightcyan3', color = 'black') + 
  theme_minimal() +
  stat_summary(fun.data = mean_sdl, mult = 1, geom = "pointrange", color = "maroon") +
  labs(x = 'Prenatal infection sum score in trimester 3', y = 'Superior frontal volume at T1')

e <- ggplot(df_age6_nonimp, aes(x = sumscore_inf_tri3_cat2, y = temporalpole_vol_cortical.x)) +
  geom_boxplot(fill = 'lightcyan3', color = 'black') +  
  theme_minimal() +
  stat_summary(fun.data = mean_sdl, mult = 1, geom = "pointrange", color = "maroon") +
  labs(x = 'Prenatal infection sum score in trimester 3', y = 'Temporal pole volume at T1')

ggarrange(a, b, c, d, e, labels = 'AUTO', ncol = 3, nrow = 2)

# Create plots comparing different quartiles of infections (violin plots)
quartiles <- quantile(df_age6_nonimp$sumscore_inf_tri3, probs = c(0.25, 0.5, 0.8))

df_age6_nonimp$sumscore_inf_tri3_quartiles <- cut(df_age6_nonimp$sumscore_inf_tri3, breaks = c(-Inf, quartiles, Inf), labels = c("Q1", "Q2", "Q3", "Q4"), include.lowest = TRUE)

a <- ggplot(df_age6_nonimp, aes(x=sumscore_inf_tri3_quartiles, y = middletemporal_vol_cortical)) + geom_violin(fill='thistle3') + theme_minimal() + stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange", color="maroon") + labs(x = 'Prenatal infection sum score in trimester 3', y = 'Middle temporal volume at T1')

b <- ggplot(df_age6_nonimp, aes(x=sumscore_inf_tri3_quartiles, y = parsorbitalis_vol_cortical)) + geom_violin(fill='thistle3') + theme_minimal() + stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange", color="maroon") + labs(x = 'Prenatal infection sum score in trimester 3', y = 'Pars orbitalis volume at T1')

c <- ggplot(df_age6_nonimp, aes(x=sumscore_inf_tri3_quartiles, y = rostralanteriorcingulate_vol_cortical.x)) + geom_violin(fill='thistle3') + theme_minimal() + stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange", color="maroon") + labs(x = 'Prenatal infection sum score in trimester 3', y = 'Rostral anterior cingulate volume at T1')

d <- ggplot(df_age6_nonimp, aes(x=sumscore_inf_tri3_quartiles, y = superiorfrontal_vol_cortical.x)) + geom_violin(fill='thistle3') + theme_minimal() + stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange", color="maroon") + labs(x = 'Prenatal infection sum score in trimester 3', y = 'Superior frontal volume at T1')

e <- ggplot(df_age6_nonimp, aes(x=sumscore_inf_tri3_quartiles, y = temporalpole_vol_cortical.x)) + geom_violin(fill='thistle3') + theme_minimal() + stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange", color="maroon") + labs(x = 'Prenatal infection sum score in trimester 3', y = 'Temporal pole volume at T1')

ggarrange(a, b, c, d, e, labels = 'AUTO', ncol = 3, nrow=2)

# Create plots comparing different quartiles of infections (boxplots)
a <- ggplot(df_age6_nonimp, aes(x = sumscore_inf_tri3_quartiles, y = middletemporal_vol_cortical)) +
  geom_boxplot(fill = 'lightcyan3', color = 'black') +  
  theme_minimal() +
  stat_summary(fun.data = mean_sdl, mult = 1, geom = "pointrange", color = "maroon") +
  labs(x = 'Prenatal infection sum score in trimester 3', y = 'Middle temporal volume at T1')

b <- ggplot(df_age6_nonimp, aes(x = sumscore_inf_tri3_quartiles, y = parsorbitalis_vol_cortical)) +
  geom_boxplot(fill = 'lightcyan3', color = 'black') +  
  theme_minimal() +
  stat_summary(fun.data = mean_sdl, mult = 1, geom = "pointrange", color = "maroon") +
  labs(x = 'Prenatal infection sum score in trimester 3', y = 'Pars orbitalis volume at T1')

c <- ggplot(df_age6_nonimp, aes(x = sumscore_inf_tri3_quartiles, y = rostralanteriorcingulate_vol_cortical.x)) +
  geom_boxplot(fill = 'lightcyan3', color = 'black') +  
  theme_minimal() +
  stat_summary(fun.data = mean_sdl, mult = 1, geom = "pointrange", color = "maroon") +
  labs(x = 'Prenatal infection sum score in trimester 3', y = 'Rostral anterior cingulate volume at T1')

d <- ggplot(df_age6_nonimp, aes(x = sumscore_inf_tri3_quartiles, y = superiorfrontal_vol_cortical.x)) +
  geom_boxplot(fill = 'lightcyan3', color = 'black') + 
  theme_minimal() +
  stat_summary(fun.data = mean_sdl, mult = 1, geom = "pointrange", color = "maroon") +
  labs(x = 'Prenatal infection sum score in trimester 3', y = 'Superior frontal volume at T1')

e <- ggplot(df_age6_nonimp, aes(x = sumscore_inf_tri3_quartiles, y = temporalpole_vol_cortical.x)) +
  geom_boxplot(fill = 'lightcyan3', color = 'black') +  
  theme_minimal() +
  stat_summary(fun.data = mean_sdl, mult = 1, geom = "pointrange", color = "maroon") +
  labs(x = 'Prenatal infection sum score in trimester 3', y = 'Temporal pole volume at T1')

ggarrange(a, b, c, d, e, labels = 'AUTO', ncol = 3, nrow = 2)

## Run analysis for FDR regions only for T2 and T2 visits
two_visit_df <- subset(imp.test_long, timepoint == 2 | timepoint == 3)
two_visit_df_mids <- as.mids(two_visit_df)

two_visits_results <- data.frame() #create empty df 

for(x in fdr_outcomes) {
  f <- paste0(x, "~ sumscore_inf_tot_standardized*timepoint + GENDER + AGE_M_v2 + ETHNMv2 + EDUCM_3groups+ GSI + SMOKE_ALL + (1|IDC)") #full model, with interaction, only with random intercept, without icv

  #calculating beta, SE, confidence interval and pval for each model
  bval <- summary(pool(with(two_visit_df_mids, lmer(as.formula(f), weights = combined_weights))),conf.int = T)[12,2]
  seval <- summary(pool(with(two_visit_df_mids, lmer(as.formula(f), weights = combined_weights))),conf.int = T)[12,3]
  lowerCI <- summary(pool(with(two_visit_df_mids, lmer(as.formula(f), weights = combined_weights))),conf.int = T)[12,7]
  upperCI <- summary(pool(with(two_visit_df_mids, lmer(as.formula(f), weights = combined_weights))),conf.int = T)[12,8]
  pval <- summary(pool(with(two_visit_df_mids, lmer(as.formula(f), weights = combined_weights))), conf.int = T)[12,6]
  
  #storing values in empty datafile  
  two_visits_results[x,1] <- bval
  two_visits_results[x,2] <- seval
  two_visits_results[x,3] <- lowerCI
  two_visits_results[x,4] <- upperCI
  two_visits_results[x,5] <- pval

  #assigning names to columns 
  colnames(two_visits_results) <- c("bval", "seval", "lowerCI", "upperCI", "pval")
}

pval1 <- unlist(two_visits_results[, 5]) 
two_visits_results$pval_fdr_adjusted <- p.adjust(pval1, method = 'fdr')
is.num1 <- sapply(two_visits_results, is.numeric)
two_visits_results[is.num1] <- lapply(two_visits_results[is.num1], round, 3)

## Run analysis for FDR regions only for T1 and T2 visits
t1_t2_visit_df <- subset(imp.test_long, timepoint == 1 | timepoint == 2)
t1_t2_visit_df_mids <- as.mids(t1_t2_visit_df)

t1_t2_visits_only_results <- data.frame() #create empty df 

for(x in fdr_outcomes) {
  f <- paste0(x, "~ sumscore_inf_tot_standardized*timepoint + GENDER + AGE_M_v2 + ETHNMv2 + EDUCM_3groups+ GSI + SMOKE_ALL + (1|IDC)") #full model, with interaction, only with random intercept, without icv
  
  #calculating beta, SE, confidence interval and pval for each model
  bval <- summary(pool(with(t1_t2_visit_df_mids, lmer(as.formula(f), weights = combined_weights))),conf.int = T)[12,2]
  seval <- summary(pool(with(t1_t2_visit_df_mids, lmer(as.formula(f), weights = combined_weights))),conf.int = T)[12,3]
  lowerCI <- summary(pool(with(t1_t2_visit_df_mids, lmer(as.formula(f), weights = combined_weights))),conf.int = T)[12,7]
  upperCI <- summary(pool(with(t1_t2_visit_df_mids, lmer(as.formula(f), weights = combined_weights))),conf.int = T)[12,8]
  pval <- summary(pool(with(t1_t2_visit_df_mids, lmer(as.formula(f), weights = combined_weights))), conf.int = T)[12,6]
  
  #storing values in empty datafile  
  t1_t2_visits_only_results[x,1] <- bval
  t1_t2_visits_only_results[x,2] <- seval
  t1_t2_visits_only_results[x,3] <- lowerCI
  t1_t2_visits_only_results[x,4] <- upperCI
  t1_t2_visits_only_results[x,5] <- pval
  
  #assigning names to columns 
  colnames(t1_t2_visits_only_results) <- c("bval", "seval", "lowerCI", "upperCI", "pval")
}

pval1 <- unlist(t1_t2_visits_only_results[, 5]) 
t1_t2_visits_only_results$pval_fdr_adjusted <- p.adjust(pval1, method = 'fdr')
is.num1 <- sapply(t1_t2_visits_only_results, is.numeric)
t1_t2_visits_only_results[is.num1] <- lapply(t1_t2_visits_only_results[is.num1], round, 3)

library(writexl)
write_xlsx(t1_t2_visits_only_results, 't1_t2_visits_only_results.xlsx')

#\ 

## Power analysis

#https://glimmpse.samplesizeshop.org --> calculate how much power we have 
#Calculate mean value of specific brain regions for each infection group at each time point 
calculate_mean_by_group <- function(data, variable_term) {
  # Get unique values of sumscore_inf_tot
  unique_values <- unique(data$sumscore_inf_tot)
  # Calculate mean for each group
  means <- sapply(unique_values, function(group) {
    mean(subset(data, sumscore_inf_tot == group)[[variable_term]], na.rm = TRUE)
  })
  # Combine results into a data frame
  result <- data.frame(sumscore_inf_tot = unique_values, mean_value = means)
  return(result)
}

calculate_mean_by_group(df_age6_nonimp, "temporalpole_vol_cortical.x")
calculate_mean_by_group(df_age10_nonimp, "temporalpole_vol_cortical.x")
calculate_mean_by_group(df_age14_nonimp, "temporalpole_vol_cortical.x")

calculate_mean_by_group(df_age6_nonimp, "middletemporal_vol_cortical")
calculate_mean_by_group(df_age10_nonimp, "middletemporal_vol_cortical")
calculate_mean_by_group(df_age14_nonimp, "middletemporal_vol_cortical")

calculate_mean_by_group(df_age6_nonimp, "rostralanteriorcingulate_vol_cortical.x")
calculate_mean_by_group(df_age10_nonimp, "rostralanteriorcingulate_vol_cortical.x")
calculate_mean_by_group(df_age14_nonimp, "rostralanteriorcingulate_vol_cortical.x")

calculate_mean_by_group(df_age6_nonimp, "superiorfrontal_vol_cortical.x")
calculate_mean_by_group(df_age10_nonimp, "superiorfrontal_vol_cortical.x")
calculate_mean_by_group(df_age14_nonimp, "superiorfrontal_vol_cortical.x")

calculate_mean_by_group(df_age6_nonimp, "parsorbitalis_vol_cortical")
calculate_mean_by_group(df_age10_nonimp, "parsorbitalis_vol_cortical")
calculate_mean_by_group(df_age14_nonimp, "parsorbitalis_vol_cortical")

#calculate raw SE
raw_fdr_outcomes <- c('middletemporal_vol_cortical', 'parsorbitalis_vol_cortical', 'rostralanteriorcingulate_vol_cortical.x', 'superiorfrontal_vol_cortical.x', 'temporalpole_vol_cortical.x')

results_mm_raw_brain_unstandardized <- data.frame() 

for(x in raw_fdr_outcomes) {
  #model specification 
  f <- paste0(x, "~ sumscore_inf_tot*timepoint + GENDER + AGE_M_v2 + ETHNMv2 + EDUCM_3groups+ GSI + SMOKE_ALL + (1|IDC)") 
  #calculating beta, SE, confidence interval and pval for each model
  seval_int_no_icv <- summary(pool(with(imp.test_mids, lmer(as.formula(f), weights = combined_weights))),conf.int = T)[12,3]
  #storing values in empty datafile  
  results_mm_raw_brain_unstandardized[x,1] <- seval_int_no_icv
  #assigning names to columns 
  colnames(results_mm_raw_brain_unstandardized) <- c("seval_int_no_icv")
}

#\

## Sample size calculation mixed-effects model (for each timepoint individually)
library(pwr)

# Set parameters
alpha <- 0.05
power_target <- 0.8
sample_size_f5 <- 582
sample_size_f10 <- 1705
sample_size_f14 <- 1543
sample_size_total <- 2155

# Perform power analysis to estimate effect size
power_result1 <- pwr::pwr.t.test(n = sample_size_f5, d = NULL, sig.level = alpha, power = power_target, type = "two.sample")
power_result2 <- pwr::pwr.t.test(n = sample_size_f10, d = NULL, sig.level = alpha, power = power_target, type = "two.sample")
power_result3 <- pwr::pwr.t.test(n = sample_size_f14, d = NULL, sig.level = alpha, power = power_target, type = "two.sample")
power_result4 <- pwr::pwr.t.test(n = sample_size_total, d = NULL, sig.level = alpha, power = power_target, type = "two.sample")

# Extract the estimated effect size
estimated_effect_size1 <- power_result1$d
estimated_effect_size2 <- power_result2$d
estimated_effect_size3 <- power_result3$d
estimated_effect_size4 <- power_result4$d

# Print the estimated effect size 
estimated_effect_size1
estimated_effect_size2
estimated_effect_size3
estimated_effect_size4

#\

## Sample size calculation mixed-effects model (for interaction term)
library(simr)
df <- complete(imp.test_mids, 30)

# Specify models 
model_unadjusted <- lmer(genr_tbv_standardized ~ sumscore_inf_tri3_standardized*timepoint + (1|IDC),  data = df)
model_with_ipw <- lmer(genr_tbv_standardized ~ sumscore_inf_tri3_standardized*timepoint + AGE_M_v2 + ETHNMv2 + EDUCM_3groups+ GSI + SMOKE_ALL + (1|IDC), weights = combined_weights, data = df)
model_without_ipw <- lmer(genr_tbv_standardized ~ sumscore_inf_tri3_standardized*timepoint + AGE_M_v2 + ETHNMv2 + EDUCM_3groups+ GSI + SMOKE_ALL + (1|IDC), data = df)

# Perform power analysis for the interaction term
# Here, you're interested in the power to detect the interaction effect
# 'test = F("fixedEffectTerm")' is a placeholder; replace 'fixedEffectTerm' with your actual term, like 'sumscore_inf_tri3_standardized:timepoint'

power_analysis_unadjusted <- powerSim(model_unadjusted, test = fixed("sumscore_inf_tri3_standardized:timepoint", "z"), nsim = 1000)
power_analysis_ipw_adjusted <- powerSim(model_with_ipw, test = fixed("sumscore_inf_tri3_standardized:timepoint", "z"), nsim = 1000)
power_analysis_no_ipw_adjusted <- powerSim(model_without_ipw, test = fixed("sumscore_inf_tri3_standardized:timepoint", "z"), nsim = 1000)

# Print the results
print(power_analysis_unadjusted)
print(power_analysis_ipw_adjusted)
print(power_analysis_no_ipw_adjusted)

#\

## Sensitivity analysis with fever
# Load data of final sample 
load('no_outcome_imp_ipw_mids_final_revisions.RDS') 

# Set wd to results folder 
setwd('V:\\medewerkers\\039088 Suleri, A\\Projecten\\05. Project_normative_development\\Results') 

# Standardizing vars in final sample
imp.test_long_rev <- complete(ipw_mids_final_revisions, include = T, action = "long")

outcome_vars <- c("genr_tbv","Cerebellum_Cortex_vol_subcortical"
                  ,"Amygdala_vol_subcortical",
                  "Hippocampus_vol_subcortical","Caudate_vol_subcortical",
                  "Putamen_vol_subcortical", "Thalamus_Proper_vol_subcortical"
                  ,"Pallidum_vol_subcortical",
                  "bankssts_vol_cortical",
                  "caudalanteriorcingulate_vol_cortical",
                  "caudalmiddlefrontal_vol_cortical",
                  "cuneus_vol_cortical","entorhinal_vol_cortical",
                  "fusiform_vol_cortical",
                  "inferiorparietal_vol_cortical","inferiortemporal_vol_cortical",
                  "isthmuscingulate_vol_cortical","lateraloccipital_vol_cortical",
                  "lateralorbitofrontal_vol_cortical",
                  "lingual_vol_cortical","medialorbitofrontal_vol_cortical",
                  "middletemporal_vol_cortical",
                  "parahippocampal_vol_cortical","paracentral_vol_cortical",
                  "parsopercularis_vol_cortical",
                  "parsorbitalis_vol_cortical","parstriangularis_vol_cortical",
                  "pericalcarine_vol_cortical",
                  "postcentral_vol_cortical","posteriorcingulate_vol_cortical",
                  "precentral_vol_cortical",
                  "precuneus_vol_cortical","rostralanteriorcingulate_vol_cortical",
                  "rostralmiddlefrontal_vol_cortical","superiorfrontal_vol_cortical",
                  "superiorparietal_vol_cortical",
                  "superiortemporal_vol_cortical","supramarginal_vol_cortical",
                  "frontalpole_vol_cortical",
                  "temporalpole_vol_cortical", "transversetemporal_vol_cortical" ,
                  "insula_vol_cortical.x") 

for (x in outcome_vars){ #making z-scores of raw brain volumes 
  t <- imp.test_long_rev[x]
  colname <- paste0(colnames(t), "_standardized")
  imp.test_long_rev[,colname] <- as.numeric(scale(t))
}

# Check structure of vars and recode if needed 
imp.test_long_rev$timepoint <- as.factor(imp.test_long_rev$timepoint) 
imp.test_long_rev$timepoint <- as.numeric(imp.test_long_rev$timepoint)

imp.test_long_rev$sumscore_inf_tot <- as.numeric(imp.test_long_rev$sumscore_inf_tot)
imp.test_long_rev$sumscore_inf_tot_standardized <- as.numeric(scale(imp.test_long_rev$sumscore_inf_tot))
imp.test_long_rev$sumscore_inf_tri1 <- as.numeric(imp.test_long_rev$sumscore_inf_tri1)
imp.test_long_rev$sumscore_inf_tri2 <- as.numeric(imp.test_long_rev$sumscore_inf_tri2)
imp.test_long_rev$sumscore_inf_tri3 <- as.numeric(imp.test_long_rev$sumscore_inf_tri3)

imp.test_long_rev$sumscore_inf_tot_standardized <- as.numeric(scale(imp.test_long_rev$sumscore_inf_tot))
imp.test_long_rev$sumscore_inf_tri1_standardized <- as.numeric(scale(imp.test_long_rev$sumscore_inf_tri1))
imp.test_long_rev$sumscore_inf_tri2_standardized <- as.numeric(scale(imp.test_long_rev$sumscore_inf_tri2))
imp.test_long_rev$sumscore_inf_tri3_standardized <- as.numeric(scale(imp.test_long_rev$sumscore_inf_tri3))

# Create infection score without fever
imp.test_long_rev$sumscore_inf_nofever_tri3 <- imp.test_long_rev$flu_tri3 + imp.test_long_rev$dermatitis_tri3 + imp.test_long_rev$herpeszoster_tri3 + imp.test_long_rev$jaundice_tri3 + imp.test_long_rev$clean_upper_resp_inf_tri3 + imp.test_long_rev$clean_uwi_tri3 + imp.test_long_rev$clean_GI_inf_tri3 + imp.test_long_rev$STD_tri3 + imp.test_long_rev$clean_lower_resp_inf_tri3

# Convert back to mids object 
imp.test_mids_revisions <- as.mids(imp.test_long_rev) 

# Run fdr sign regions with infection score without fever and with interaction of fever
outcomes <- c('middletemporal_vol_cortical_standardized', 'parsorbitalis_vol_cortical_standardized', 'rostralanteriorcingulate_vol_cortical_standardized', 'superiorfrontal_vol_cortical_standardized', 'temporalpole_vol_cortical_standardized')

results_fever_int <- data.frame()

for(x in outcomes) {
  f <- paste0(x, "~ sumscore_inf_nofever_tri3*timepoint + GENDER + AGE_M_v2 + ETHNMv2 + EDUCM_3groups+ GSI + SMOKE_ALL + (1|IDC)") #shortened infection score
  g <-  paste0(x, "~ sumscore_inf_nofever_tri3*fever_tri3*timepoint + GENDER + AGE_M_v2 + ETHNMv2 + EDUCM_3groups+ GSI + SMOKE_ALL + (1|IDC)") #interaction with fever
  
  #calculating beta, SE, confidence interval and pval for each model
  bval <- summary(pool(with(imp.test_mids_revisions, lmer(as.formula(f), weights = combined_weights))),conf.int = T)[12,2]
  seval <- summary(pool(with(imp.test_mids_revisions, lmer(as.formula(f), weights = combined_weights))),conf.int = T)[12,3]
  lowerCI <- summary(pool(with(imp.test_mids_revisions, lmer(as.formula(f), weights = combined_weights))),conf.int = T)[12,7]
  upperCI <- summary(pool(with(imp.test_mids_revisions, lmer(as.formula(f), weights = combined_weights))),conf.int = T)[12,8]
  pval <- summary(pool(with(imp.test_mids_revisions, lmer(as.formula(f), weights = combined_weights))), conf.int = T)[12,6]
  
  bval_int <- summary(pool(with(imp.test_mids_revisions, lmer(as.formula(g), weights = combined_weights))),conf.int = T)[16,2]
  seval_int <- summary(pool(with(imp.test_mids_revisions, lmer(as.formula(g), weights = combined_weights))),conf.int = T)[16,3]
  lowerCI_int <- summary(pool(with(imp.test_mids_revisions, lmer(as.formula(g), weights = combined_weights))),conf.int = T)[16,7]
  upperCI_int <- summary(pool(with(imp.test_mids_revisions, lmer(as.formula(g), weights = combined_weights))),conf.int = T)[16,8]
  pval_int <- summary(pool(with(imp.test_mids_revisions, lmer(as.formula(g), weights = combined_weights))), conf.int = T)[16,6]
  
  #storing values in empty datafile  
  results_fever_int[x,1] <- bval
  results_fever_int[x,2] <- seval
  results_fever_int[x,3] <- lowerCI
  results_fever_int[x,4] <- upperCI
  results_fever_int[x,5] <- pval
  results_fever_int[x,6] <- bval_int
  results_fever_int[x,7] <- seval_int
  results_fever_int[x,8] <- lowerCI_int
  results_fever_int[x,9] <- upperCI_int
  results_fever_int[x,10] <- pval_int
  
  #assigning names to columns 
  colnames(results_fever_int) <- c("bval", "seval", "lowerCI", "upperCI", "pval","bval_int", "seval_int", "lowerCI_int", "upperCI_int", "pval_int")
}

write.xlsx(results_fever_int, 'results_fever_int.xlsx')

#\ END OF SCRIPT, GO TO PART C FOR FIGURES 
