#------------------------------------------------------------------------#
#############PROJECT: infection & child brain development#################
#------------------------------------------------------------------------#

# Author: Anna Suleri

### Parts in this script:
#' Part 1: Loading, cleaning, renaming data and merging final data frame
#' Part 2: Applying inclusion and exclusion criteria & create an inclusion variable
#' Part 3: Adapt to long format

#########################################################
#-------------------Data preparation--------------------#
########################################################

rm(list = ls()) #clears the environment

### Opening libraries we need & setting wd 
setwd('Set_path_to_your_wd')
wd <- getwd()

libraries <- c('foreign', 'haven', 'dplyr', 'tidyverse' ,'mice', 'lme4', 'lattice', 'ggplot2', 'ggeffects', 'splines', 'reshape2', 'gridExtra', 'grid', 'ggpubr', 'tidyr', 'broom.mixed', 'xlsx', 'corrplot', 'stringi', 'miceadds', 'mitools', 'CBPS', 'survey', 'survival', 'Hmisc', 'ggpubr')

invisible(lapply(libraries, require, character.only = T))

### Loading in datasets and selecting variables we need for this project  

## Exposure + covariates
prenatal_stress_imp_df <- readRDS("imputation_list_full.rds")
prenatal_stress_df2 <- complete(prenatal_stress_imp_df, 30) 
prenatal_stress_df <- dplyr::select(prenatal_stress_df2, c('IDC', 'post_life_events', 'post_direct_victimization'))

# Create maternal inflammatory score and inflammatory medication variable  
df1 <- read.spss('Second_hits_DF.sav', to.data.frame = T)
df_med <- read.spss('MEDICATIONSELFREPORTPREGNANCY_30112017.sav', to.data.frame = T)
df_inflam<- read.spss('GR1001-D1-37_16072021.sav', to.data.frame = T)

inflam_vars <- merge(df1[, c( 'IDC','IDM', 'DIAB_GRA', 'PIH_v1', 'preeclampsia')], df_med[, c('IDM', 'NSAIDTOT', 'ABIOTOT', 'PMOLTOT', 'CORTTOT', 'MUCOTOT' ,'COUGHTOT')], by = 'IDM', all.x = T)
inflam_vars <- merge(inflam_vars, df_inflam[, c('IDM','d2300101', 'd2400101', 'd2500101', 'd1100101')], by = 'IDM', all.x = T)

inflam_vars$preg_ic <- as.factor(ifelse(inflam_vars$DIAB_GRA == 'Yes' | inflam_vars$preeclampsia == 'Yes' | inflam_vars$PIH_v == 'Yes', 1, ifelse(inflam_vars$DIAB_GRA == "No" & inflam_vars$preeclampsia == 'No' & inflam_vars$PIH_v1 == 'No', 0, NA)))

inflam_vars$ic <- as.factor(ifelse(inflam_vars$d2300101 == 'Yes' | inflam_vars$d2400101 == 'Yes' | inflam_vars$d2500101 == 'Yes' | inflam_vars$d1100101 == 'Yes', 1, ifelse(inflam_vars$d2300101 == 'No' & inflam_vars$d2400101 == 'No' & inflam_vars$d2500101 == 'No' & inflam_vars$d1100101 == 'No', 0, NA)))

inflam_vars$mat_inflam <- as.factor(ifelse(inflam_vars$preg_ic == 1 | inflam_vars$ic == 1, 1, ifelse(inflam_vars$preg_ic == 0 & inflam_vars$ic == 0, 0, NA)))

inflam_vars$im <- as.factor(ifelse(inflam_vars$NSAIDTOT == 'no use' & inflam_vars$ABIOTOT == 'no use' & inflam_vars$PMOLTOT == 'no use' & inflam_vars$CORTTOT == 'no use' & inflam_vars$MUCOTOT == 'no use' & inflam_vars$COUGHTOT == 'no use', 0, ifelse(is.na(inflam_vars$NSAIDTOT) | is.na(inflam_vars$ABIOTOT) | is.na(inflam_vars$PMOLTOT) | is.na(inflam_vars$CORTTOT) | is.na(inflam_vars$MUCOTOT) | is.na(inflam_vars$COUGHTOT), NA, 1)))

# Select exposure vars and covars of interest 
dd1 <- df1[, c('IDM', 'IDC', 'MOTHER', 'AGE_M_v2', 'ETHNMv2', 'SMOKE_ALL', 'GSI', 'sumscore_inf_tri1', 'sumscore_inf_tri2', 'sumscore_inf_tri3', 'sumscore_inf_tot', 'GENDER', 'INTAKEPERIOD', 'mdrink_updated', 'cbcl_sum_14', 'sum_int_14', 'sum_ext_14', 'f1100101', 'fever_tri1', 'fever_tri2', 'fever_tri3', 'flu_tri1', 'flu_tri2', 'flu_tri3', 'dermatitis_tri1', 'dermatitis_tri2', 'dermatitis_tri3', 'herpeszoster_tri1', 'herpeszoster_tri2', 'herpeszoster_tri3', 'jaundice_tri1', 'jaundice_tri2', 'jaundice_tri3', 'clean_upper_resp_inf_tri1', 'clean_upper_resp_inf_tri2', 'clean_upper_resp_inf_tri3', 'clean_uwi_tri1', 'clean_uwi_tri2', 'clean_uwi_tri3', 'clean_GI_inf_tri1', 'clean_GI_inf_tri2', 'clean_GI_inf_tri3', 'STD_tri1', 'STD_tri2', 'STD_tri3', 'clean_lower_resp_inf_tri1', 'clean_lower_resp_inf_tri2', 'clean_lower_resp_inf_tri3')]
df2 <- read.spss('Covariates_MRI_analyses.sav', to.data.frame = T)
dd2 <- df2[, c('IDC', 'EDUCM_3groups', 'INCOME')]
dd2$INCOME <- as.factor(ifelse(dd2$INCOME == 'less than 450' | dd2$INCOME == '450-600 euro' | dd2$INCOME == '600-700 euro' | dd2$INCOME == '700-800 euro' | dd2$INCOME == '800-900 euro' | dd2$INCOME == '900-1200 euro' | dd2$INCOME == '1200-1400 euro' | dd2$INCOME == '1400-1600 euro' | dd2$INCOME == '1600-1800 euro' | dd2$INCOME == '1800-2000 euro' | dd2$INCOME == '2000-2200 euro', '<2000', ifelse(is.na(dd2$INCOME), NA, '>2000')))

## Outcomes
# Select wave 1 vars 
df3 <- readRDS("f05_freesurfer_v6_24june2021_aparc_stats_pull18Aug2021.rds")
dd3 <- dplyr::select(df3, "idc" ,ends_with('vol_f05')) #cortical regions 
df4 <- readRDS("f05_freesurfer_v6_24june2021_aseg_stats_pull18Aug2021_v1.rds")
dd4 <- df4[, c('idc',
               'Left_Cerebellum_Cortex_vol_f05', 'Left_Amygdala_vol_f05', 'Left_Hippocampus_vol_f05', 'Left_Lateral_Ventricle_vol_f05', 'Left_Caudate_vol_f05', 'Left_Putamen_vol_f05', 'Left_Thalamus_Proper_vol_f05', 'Left_Pallidum_vol_f05', 
               'Right_Cerebellum_Cortex_vol_f05', 'Right_Amygdala_vol_f05', 'Right_Hippocampus_vol_f05', 'Right_Lateral_Ventricle_vol_f05', 'Brain_Stem_vol_f05', 'Right_Caudate_vol_f05', 'Right_Putamen_vol_f05', 'Right_Thalamus_Proper_vol_f05', 'Right_Pallidum_vol_f05',
               'CC_Posterior_vol_f05', 'CC_Mid_Posterior_vol_f05', 'CC_Central_vol_f05', 'CC_Mid_Anterior_vol_f05', 'CC_Anterior_vol_f05')] #subcortical regions 
df5 <- readRDS("f05_freesurfer_v6_24june2021_tbv_stats_pull18Aug2021_v2.rds")
dd5 <- df5[, c('idc', 'TotalGrayVol_f05', 'CerebralWhiteMatterVol_f05', 'genr_tbv_f05', 'eTIV_f05')]

# Select wave 2 vars 
df6 <- readRDS('f09_freesurfer_v6_09dec2016_aparc_stats_pull06june2017.rds')
dd6 <-  dplyr::select(df6, "idc" ,ends_with('vol_f09')) #cortical regions 
df7 <- readRDS("f09_freesurfer_v6_09dec2016_aseg_stats_pull06june2017_v1.rds")
dd7 <- df7[, c('idc',
               'Left_Cerebellum_Cortex_vol_f09', 'Left_Amygdala_vol_f09', 'Left_Hippocampus_vol_f09', 'Left_Lateral_Ventricle_vol_f09', 'Left_Caudate_vol_f09', 'Left_Putamen_vol_f09', 'Left_Thalamus_Proper_vol_f09', 'Left_Pallidum_vol_f09', 
               'Right_Cerebellum_Cortex_vol_f09', 'Right_Amygdala_vol_f09', 'Right_Hippocampus_vol_f09', 'Right_Lateral_Ventricle_vol_f09', 'Brain_Stem_vol_f09', 'Right_Caudate_vol_f09', 'Right_Putamen_vol_f09', 'Right_Thalamus_Proper_vol_f09', 'Right_Pallidum_vol_f09',
               'CC_Posterior_vol_f09', 'CC_Mid_Posterior_vol_f09', 'CC_Central_vol_f09', 'CC_Mid_Anterior_vol_f09', 'CC_Anterior_vol_f09')] #subcortical regions 
df8 <- readRDS("f09_freesurfer_v6_09dec2016_tbv_stats_pull20june2017_v2.rds")
dd8 <- df8[, c('idc', 'TotalGrayVol_f09', 'CerebralWhiteMatterVol_f09', 'genr_tbv_f09', 'eTIV_f09')]

# Select wave 3 vars 
df9 <- readRDS("f13_freesurfer_v6_14oct2020_aparc_stats_pull23Nov2020.rds")
dd9 <- dplyr::select(df9, "idc",ends_with('vol_f13')) #cortical regions
df10 <- readRDS("f13_freesurfer_v6_14oct2020_aseg_stats_pull23Nov2020_v1.rds")
dd10 <- df10[, c('idc',
                 'Left_Cerebellum_Cortex_vol_f13', 'Left_Amygdala_vol_f13', 'Left_Hippocampus_vol_f13', 'Left_Lateral_Ventricle_vol_f13', 'Brain_Stem_vol_f13', 'Left_Caudate_vol_f13', 'Left_Putamen_vol_f13', 'Left_Thalamus_Proper_vol_f13', 'Left_Pallidum_vol_f13', 
                 'Right_Cerebellum_Cortex_vol_f13', 'Right_Amygdala_vol_f13', 'Right_Hippocampus_vol_f13', 'Right_Lateral_Ventricle_vol_f13', 'Right_Caudate_vol_f13', 'Right_Putamen_vol_f13', 'Right_Thalamus_Proper_vol_f13', 'Right_Pallidum_vol_f13',
                 'CC_Posterior_vol_f13', 'CC_Mid_Posterior_vol_f13', 'CC_Central_vol_f13', 'CC_Mid_Anterior_vol_f13', 'CC_Anterior_vol_f13')] #subcortical regions 
df11 <- readRDS("f13_freesurfer_v6_14oct2020_tbv_stats_pull23Nov2020_v2.rds")
dd11 <- df11[, c('idc', 'TotalGrayVol_f13', 'CerebralWhiteMatterVol_f13', 'genr_tbv_f13', 'eTIV_f13')]

# Select core imaging file 
df12 <- readRDS("genr_mri_core_data_20221010.rds")
dd12 <- df12[, c('idc', 'mri_consent_f05', 'mri_consent_f09', 'mri_consent_f13', 'age_child_mri_f05', 'age_child_mri_f09', 'age_child_mri_f13', 't1_has_nii_f05', 't1_has_nii_f09', 't1_has_nii_f13', 't1_asset_has_nii_f09', 'has_braces_mri_f05', 'has_braces_mri_f09', 'has_braces_mri_f13', 'exclude_incidental_f05', 'exclude_incidental_f09', 'exclude_incidental_f13', 'freesurfer_qc_f05', 'freesurfer_qc_f09', 'freesurfer_qc_f13')]

# Renaming id vars for consistence across dataframes 
dd3 <- rename(dd3, IDC = idc)
dd4 <- rename(dd4, IDC = idc)
dd5 <- rename(dd5, IDC = idc)
dd6 <- rename(dd6, IDC = idc)
dd7 <- rename(dd7, IDC = idc)
dd8 <- rename(dd8, IDC = idc)
dd9 <- rename(dd9, IDC = idc)
dd10 <- rename(dd10, IDC = idc)
dd11 <- rename(dd11, IDC = idc)
dd12 <- rename(dd12, IDC = idc)

## Select variables for IPW
# Wave 1: CBCL tot, SRS
df13 <-  read.spss('CHILDCBCL_6_incl_Tscores_20201111.sav', to.data.frame = T)
dd13 <- df13[, c('IDC', 'cbcl_sum_5')]
df14 <- read.spss('Covariaten_Anna.sav', to.data.frame = T)
dd14 <- df14[, c('IDC', 'srs_weighted')]

# Baseline
df15 <- read.spss('MEDICATIONSELFREPORTPREGNANCY_30112017.sav', to.data.frame = T)
dd15 <- df15[, c('IDM', 'SSRITOT')]
df16 <- read.spss('MATERNALFOLICACID_23062010.sav', to.data.frame = T)
df17 <- read.spss('GR1001-F11-12_22112016.sav', to.data.frame = T)
dd17 <- df17[, c('IDM', 'f1200201', 'f1200101')]
dd17_use <- subset(dd17, f1200101 == 'Yes still' & f1200201 == 'Marihuana')

### Merging all dataframes 
merger1 <- merge(dd1, dd2, by = 'IDC')
merger2 <- merge(merger1, dd3, by = 'IDC', all.x = T)
merger3 <- merge(merger2, dd4, by = 'IDC', all.x = T)
merger4 <- merge(merger3, dd5, by = 'IDC', all.x = T)
merger5 <- merge(merger4, dd6, by = 'IDC', all.x = T)
merger6 <- merge(merger5, dd7, by = 'IDC', all.x = T)
merger7 <- merge(merger6, dd8, by = 'IDC', all.x = T)
merger8 <- merge(merger7, dd9, by = 'IDC', all.x = T)
merger9 <- merge(merger8, dd10, by = 'IDC', all.x = T)
merger10 <- merge(merger9, dd11, by = 'IDC', all.x = T)
merger11 <- merge(merger10, dd12, by = 'IDC', all.x = T)
merger12 <- merge(merger11, dd13, by = 'IDC', all.x = T)
merger13 <- merge(merger12, dd14, by = 'IDC', all.x = T)
merger14 <- merge(merger13, dd15, by = 'IDM', all.x = T)
merger15 <- merge(merger14, df16, by ='IDM', all.x = T)
merger16 <- merge(merger15, dd17_use, by = 'IDM', all.x = T)
merger17 <- merge(merger16, prenatal_stress_df, by = 'IDC', all.x=T)
merger18 <- merge(merger17, inflam_vars, by = 'IDC', all.x = T)

### Summing up both hemispheres (because no hypothesis on lateralized effects)
# Loop over all structures and average columns of the same structure in both hemispheres (at the same measurement occasion)
aseg_cols <- colnames(dplyr::select(merger18, starts_with("Left_"), starts_with("Right_")))
aparc_cols <- colnames(dplyr::select(merger18, starts_with("lh_"), starts_with('rh_')))
aseg_structs_temp1 <- stri_remove_empty(unlist(strsplit(aseg_cols, "Left_"))) #unlist to make vector
aseg_structs_temp2<- stri_remove_empty(unlist(strsplit(aseg_structs_temp1, "Right_")))
aparc_structs_temp1 <- stri_remove_empty(unlist(strsplit(aparc_cols, "rh_")))
aparc_structs_temp2 <- stri_remove_empty(unlist(strsplit(aparc_structs_temp1, "lh_")))
aseg_structs1 <- stri_remove_empty(unlist(strsplit(aseg_structs_temp2, "_f05")))
aparc_structs1 <- stri_remove_empty(unlist(strsplit(aparc_structs_temp2, "_f05")))
aseg_structs2 <- stri_remove_empty(unlist(strsplit(aseg_structs1, "_f09")))
aparc_structs2 <- stri_remove_empty(unlist(strsplit(aparc_structs1, "_f09")))
aseg_structs3 <- stri_remove_empty(unlist(strsplit(aseg_structs2, "_f13")))
aparc_structs3 <- stri_remove_empty(unlist(strsplit(aparc_structs2, "_f13")))

structs <- c(aseg_structs3, aparc_structs3)

# Go over all structures
for(x in structs){
  #Go over all columns to identify the first column (left hemisphere)
  for(y in 1:ncol(merger18)){
    #Store the column name to find the counterpart in the other (right) hemisphere
    a <- colnames(merger18[y])
    #Go over all columns to identify the second column (right hemisphere)
    for(z in 1:ncol(merger18)){
      #Store the column name to identify the same structures in each hemisphere
      b <- colnames(merger18[z])
      #If we have the same structure for both hemispheres, calculate the mean of those
      
      ##MRI T1
      #Subcortical volumes MRI wave 1
      if(a == paste0("Left_",x,"_f05") & b == paste0("Right_",x,"_f05")){
        print(x)
        newcolname <- paste0(x,"_subcortical_f05")
        merger18[,newcolname] <- (merger18[,y] + merger18[,z])/2
      }
      #Cortical volumes MRI wave 1
      if(a == paste0("lh_",x,"_f05") & b == paste0("rh_",x,"_f05")){
        print(x)
        newcolname <- paste0(x,"_cortical_f05")
        merger18[,newcolname] <- (merger18[,y] + merger18[,z])/2
      }
      
      ##MRI T2
      #Subcortical volumes MRI wave 2
      if(a == paste0("Left_",x,"_f09") & b == paste0("Right_",x,"_f09")){
        print(x)
        newcolname <- paste0(x,"_subcortical_f09")
        merger18[,newcolname] <- (merger18[,y] + merger18[,z])/2
      }
      #Cortical volumes MRI wave 2
      if(a == paste0("lh_",x,"_f09") & b == paste0("rh_",x,"_f09")){
        print(x)
        newcolname <- paste0(x,"_cortical_f09")
        merger18[,newcolname] <- (merger18[,y] + merger18[,z])/2
      }
      
      ##MRI T3
      #Subcortical volumes MRI wave 3
      if(a == paste0("Left_",x,"_f13") & b == paste0("Right_",x,"_f13")){
        print(x)
        newcolname <- paste0(x,"_subcortical_f13")
        merger18[,newcolname] <- (merger18[,y] + merger18[,z])/2
      }
      #Cortical volumes MRI wave 3
      if(a == paste0("lh_",x,"_f13") & b == paste0("rh_",x,"_f13")){
        print(x)
        newcolname <- paste0(x,"_cortical_f13")
        merger18[,newcolname] <- (merger18[,y] + merger18[,z])/2
      }
    }
  }
}

# Now tidy up the dataset and remove all original hemisphere specific variables
df_tidy <- dplyr::select(merger18, -c(starts_with("Left_"), starts_with("Right_"), starts_with("lh_"), starts_with("rh_"))) 

# Save the dataset in an .rds file, in the directory where the raw data are stored
save(df_tidy, file = 'df_tidy_updated.Rds') 

load("df_tidy_updated.Rds")

#\

############################################################
#-------------------Inclusion/Exclusion--------------------#
###########################################################

# Make one variable to indicate whether someone was included in our sample or not from the whole cohort. So for every IDC who is in df_tidy but not df we say 'non included' and for the ones that are in df and df_tidy we say 'included'. After which impute this df in wide format. 

# Inclusion criteria 
inclu1 <- subset(df_tidy, INTAKEPERIOD == '<18 weeks') #including only mothers who enrolled in tri1, n=7145
inclu2 <- subset(inclu1, complete.cases(sumscore_inf_tot)) #no missing information in prenatal infection, n=4259
inclu3 <- subset(inclu2, complete.cases(genr_tbv_f05) | complete.cases(genr_tbv_f09) | complete.cases(genr_tbv_f13)) #one of the time points should have complete cases

# Exclusion criteria
#'Excluding 1 twin/sibling based on most available information and if equal then at random 
exclu1 <- inclu3[sample(nrow(inclu3)),] #making a random order in the df
exclu1$na_count <- apply(exclu1, 1, function(x) sum(is.na(x))) #add new column to df with na count per row/IDC
exclu2 <- exclu1[order(exclu1$na_count),] #order dataset from least to most missing
exclu3 <- exclu2[!duplicated(exclu2$MOTHER, fromLast = T),] #delete 1 twin/sibling, deleting when a duplicate is found (so second one) which also has the most missings

#'Excluding everyone with unusable scan or poor quality or IF or braces scan
exclu3$inclu_f05 <- as.factor(ifelse(exclu3$mri_consent_f05 == 'yes' & exclu3$t1_has_nii_f05 == 'yes' & exclu3$has_braces_mri_f05 == 'no' & exclu3$exclude_incidental_f05 == 'include' & exclu3$freesurfer_qc_f05 == 'usable', 'include' , 'exclude'))
exclu3$inclu_f09 <- as.factor(ifelse(exclu3$mri_consent_f09 == 'yes' & exclu3$t1_has_nii_f09 == 'yes' & exclu3$t1_asset_has_nii_f09 != 'exclude' & exclu3$has_braces_mri_f09 == 'no' &  exclu3$exclude_incidental_f09 == 'include' & exclu3$freesurfer_qc_f09 == 'usable', 'include', 'exclude'))
exclu3$inclu_f13 <- as.factor(ifelse(exclu3$mri_consent_f13 == 'yes' & exclu3$t1_has_nii_f13 == 'yes' & exclu3$has_braces_mri_f13 == 'no' & exclu3$exclude_incidental_f13 == 'include' &exclu3$freesurfer_qc_f13 == 'usable', 'include', 'exclude'))

df_f05 <- subset(exclu3, complete.cases(genr_tbv_f05) & inclu_f05 == 'include')
df_f09 <- subset(exclu3, complete.cases(genr_tbv_f09) & inclu_f09 == 'include')
df_f13 <- subset(exclu3, complete.cases(genr_tbv_f13) & inclu_f13 == 'include')

dd <- merge(df_f05, df_f09, by = 'IDC', all = T)
dd2 <- merge(dd, df_f13, by = 'IDC', all = T) #ids of children with usable brain at each timepoint after excluding IF or braces or unusable scan [depending on complete cases time point]

df_final <- exclu3[(exclu3$IDC %in% dd2$IDC),] #final sample

# Create inclusion var 
df_tidy$include <- as.factor(ifelse(!df_tidy$IDC %in% df_final$IDC, 1, 0)) #1 = include in final sample

summary(df_tidy$include)

save(df_tidy, file = "df_tidy_include.Rds")

#\

####################################################
#-----------Figure infections over time------------#
####################################################

## Infection prevalence plot for each trimester
infection_columns <- c('clean_upper_resp_inf_tri1', 'clean_upper_resp_inf_tri2', 'clean_upper_resp_inf_tri3', 'clean_GI_inf_tri1', 'clean_GI_inf_tri2', 'clean_GI_inf_tri3', 'flu_tri1', 'flu_tri2', 'flu_tri3', 'clean_uwi_tri1', 'clean_uwi_tri2', 'clean_uwi_tri3', 'dermatitis_tri1', 'dermatitis_tri2', 'dermatitis_tri3', 'clean_lower_resp_inf_tri1', 'clean_lower_resp_inf_tri2', 'clean_lower_resp_inf_tri3', 'herpeszoster_tri1', 'herpeszoster_tri2', 'herpeszoster_tri3', 'jaundice_tri1', 'jaundice_tri2', 'jaundice_tri3', 'STD_tri1', 'STD_tri2', 'STD_tri3', 'fever_tri1', 'fever_tri2', 'fever_tri3')

# From numeric to factor, so we have a yes/no variable 
df_final[, infection_columns] <- lapply(df_final[, infection_columns], factor) 

# Calculate percentage of 'yes' for each infection type per trimester
for(i in infection_columns) {
  yes_count <- sum(df_final[[i]] == 1, na.rm = T)
  total_responses <- sum(!is.na(df_final[[i]]))
  percentage_yes <- (yes_count / total_responses) * 100
  message(i)
  print(percentage_yes)
}

# Create frequency plots per infection type per trimester 
tri1_df <- data.frame(
  infection = as.factor(c("Upper respiratory infection", "Gastro-intestinal infection", "Flu", "Urinary tract infection", "Dermatitis", "Lower respiratory infection", "Herpes Zoster", "Eye infection", "STD", "Fever")),
  prevalence = as.numeric(c(55.3, 18.3, 17.1, 3.0, 3.4, 0.4, 0.2, 0.1, 0.2, 9.3))
)
tri1_df$timing <- 'Trimester 1'

tri2_df <- data.frame(
  infection = as.factor(c("Upper respiratory infection", "Gastro-intestinal infection", "Flu", "Urinary tract infection", "Dermatitis", "Lower respiratory infection", "Herpes Zoster", "Eye infection", "STD","Fever")),
  prevalence = as.numeric(c(57.2, 15.7, 14.9, 2.7, 3.1, 0.3, 0.3, 0.2, 0.3, 6.3))
)
tri2_df$timing <- 'Trimester 2'

tri3_df <- data.frame(
  infection = as.factor(c("Upper respiratory infection", "Gastro-intestinal infection", "Flu", "Urinary tract infection", "Dermatitis", "Lower respiratory infection", "Herpes Zoster", "Eye infection", "STD","Fever")),
  prevalence = as.numeric(c(53.1, 15.9, 13.6, 3.7, 2.2, 0.4, 0.2, 0.1, 0.1, 6.1))
)
tri3_df$timing <- 'Trimester 3'

full_trimester_df <- rbind(tri1_df, tri2_df, tri3_df)
full_trimester_df$timing <- as.factor(full_trimester_df$timing)

infectionplot1 <- ggplot(full_trimester_df, aes(infection, prevalence, fill = timing)) + 
  geom_bar(stat = 'identity', position = 'dodge') + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, size = 10), legend.title = element_text(face = 'bold'), legend.position = c(0.009, 0.985), legend.justification = c(0,1), legend.box.background = element_rect(fill = 'darkgrey'),legend.box.margin=margin(1,1,1,1), legend.key = element_rect(fill = 'black'), legend.background = element_rect(fill = 'floralwhite', size = 0.5, linetype = 'solid')) + labs(x = "", y = 'Prevalence (%)') + 
  scale_fill_brewer(palette = 'Accent') + 
  guides(fill = guide_legend(title = 'Timing of infection')) + 
  geom_text(aes(label = paste0(round(prevalence, 2), "%")), position = position_dodge(width = 0.9), vjust = -0.5, size = 2)

## Visualize number and type of infections for each pregnant women over the course of gestation

# Create infections df 
df_infections <- dplyr::select(df_final, c('IDC','IDM.x', 'MOTHER', 'sumscore_inf_tot', 'sumscore_inf_tri1', 'sumscore_inf_tri2', 'sumscore_inf_tri3', 'fever_tri1', 'fever_tri2', 'fever_tri3', 'flu_tri1', 'flu_tri2', 'flu_tri3', 'dermatitis_tri1', 'dermatitis_tri2', 'dermatitis_tri3', 'herpeszoster_tri1', 'herpeszoster_tri2', 'herpeszoster_tri3', 'jaundice_tri1', 'jaundice_tri2', 'jaundice_tri3', 'clean_upper_resp_inf_tri1', 'clean_upper_resp_inf_tri2', 'clean_upper_resp_inf_tri3', 'clean_uwi_tri1', 'clean_uwi_tri2', 'clean_uwi_tri3', 'clean_GI_inf_tri1', 'clean_GI_inf_tri2', 'clean_GI_inf_tri3', 'STD_tri1', 'STD_tri2', 'STD_tri3', 'clean_lower_resp_inf_tri1', 'clean_lower_resp_inf_tri2', 'clean_lower_resp_inf_tri3'))

# Transform to long format 
df_infections_long <- reshape(df_infections, 
                              idvar = "MOTHER", 
                              varying = list(c("sumscore_inf_tri1", "sumscore_inf_tri2"
                                               , "sumscore_inf_tri3"),
                                             c("fever_tri1", "fever_tri2",
                                               "fever_tri3"),  
                                             c("flu_tri1", "flu_tri2", "flu_tri3"),
                                             c("dermatitis_tri1", "dermatitis_tri2", 
                                               "dermatitis_tri3"),
                                             c("herpeszoster_tri1", "herpeszoster_tri2"
                                               , "herpeszoster_tri3"),
                                             c("jaundice_tri1", "jaundice_tri2", 
                                               "jaundice_tri3"),
                                             c("clean_upper_resp_inf_tri1", 
                                               "clean_upper_resp_inf_tri2",
                                               "clean_upper_resp_inf_tri3"),  
                                             c("clean_uwi_tri1", "clean_uwi_tri2", 
                                               "clean_uwi_tri3"),
                                             c("clean_GI_inf_tri1", "clean_GI_inf_tri2"
                                               , "clean_GI_inf_tri3"),
                                             c("STD_tri1", "STD_tri2", "STD_tri3"),
                                             c("clean_lower_resp_inf_tri1", 
                                               "clean_lower_resp_inf_tri2",
                                               "clean_lower_resp_inf_tri3")), 
                              v.names = c("infection sum score", "fever", "flu", 
                                          "dermatitis",
                                          'herpes zoster', 'jaundice', 'upper 
                                          respiratory infections', 
                                          'urinary tract infections', 
                                          'gastrointestinal infections',
                                          'sextually transmitted diseases', 
                                          'lower respiratory infections'), 
                              direction = 'long')

df_infections_long <- df_infections_long[order(df_infections_long$IDC),] 

# Create histogram for infection sum score
hist_tot <- ggplot(df_final, aes(x = sumscore_inf_tot)) +
  geom_histogram(fill = 'mediumpurple4', color = 'black', bins = 30, alpha = 0.7) +
  labs(x = 'Total infection sum score', y = 'Count') +
  theme_minimal() + 
  geom_vline(xintercept = quantile(df_final$sumscore_inf_tot, 0.25), linetype = 'dashed', color = 'brown1') +
  geom_vline(xintercept = median(df_final$sumscore_inf_tot), linetype = 'dashed', color = 'seagreen2') +
  geom_vline(xintercept = quantile(df_final$sumscore_inf_tot, 0.75), linetype = 'dashed', color = 'orange')

hist_tri1 <- ggplot(df_final, aes(x = sumscore_inf_tri1)) +
  geom_histogram(fill = 'mediumpurple4', color = 'black', bins = 30, alpha = 0.7) +
  labs(x = 'Infection sum score (trimester 1)', y = 'Count') +
  theme_minimal() + 
  geom_vline(xintercept = quantile(df_final$sumscore_inf_tri1, 0.25), linetype = 'dashed', color = 'brown1') +
  geom_vline(xintercept = median(df_final$sumscore_inf_tri1), linetype = 'dashed', color = 'seagreen2') +
  geom_vline(xintercept = quantile(df_final$sumscore_inf_tri1, 0.75), linetype = 'dashed', color = 'orange')

hist_tri2 <- ggplot(df_final, aes(x = sumscore_inf_tri2)) +
  geom_histogram(fill = 'mediumpurple4', color = 'black', bins = 30, alpha = 0.7) +
  labs(x = 'Infection sum score (trimester 2)', y = 'Count') +
  theme_minimal() + 
  geom_vline(xintercept = quantile(df_final$sumscore_inf_tri2, 0.25), linetype = 'dashed', color = 'brown1') +
  geom_vline(xintercept = median(df_final$sumscore_inf_tri2), linetype = 'dashed', color = 'seagreen2') +
  geom_vline(xintercept = quantile(df_final$sumscore_inf_tri2, 0.75), linetype = 'dashed', color = 'orange')

hist_tri3 <- ggplot(df_final, aes(x = sumscore_inf_tri3,)) +
  geom_histogram(fill = 'mediumpurple4', color = 'black', bins = 30, alpha = 0.7) +
  labs(x = 'Infection sum score (trimester 3)', y = 'Count') +
  theme_minimal() + 
  geom_vline(xintercept = quantile(df_final$sumscore_inf_tri3, 0.25), linetype = 'dashed', color = 'brown1') +
  geom_vline(xintercept = median(df_final$sumscore_inf_tri3), linetype = 'dashed', color = 'seagreen2') +
  geom_vline(xintercept = quantile(df_final$sumscore_inf_tri3, 0.75), linetype = 'dashed', color = 'orange')

infectionplot2 <- ggarrange(hist_tot, hist_tri1, hist_tri2, hist_tri3, ncol = 2, nrow=2)

# Create spaghetti plot to show rate of infection over time for all mothers
infectionplot3 <- ggplot(df_infections_long, aes(x = time, y = `infection sum score`, group = IDC)) +
  geom_jitter(position = position_jitter(width = 0.5, height = 0.06), color = 'lightsalmon', alpha = 0.07) +
  geom_line(aes(group = IDC), color = 'mediumpurple4', alpha = 0.2) +
  labs(x = 'Pregnancy trimesters', y = 'Infection sum score') +
  theme_bw()

## Combine all infection plots into one plot for figure 1 in paper
ggarrange(infectionplot1,infectionplot2,infectionplot3, labels = 'AUTO', ncol = 2, nrow=2)

#\ 

####################################################
#-------------------Long format--------------------#
####################################################
# Removing columns i don't need for whole analysis 
df_final <- dplyr::select(df_tidy, -c("INTAKEPERIOD", "mri_consent_f05", "mri_consent_f09", "mri_consent_f13", "t1_has_nii_f05", "t1_has_nii_f09", "t1_has_nii_f13", "t1_asset_has_nii_f09", "has_braces_mri_f05", "has_braces_mri_f09", "has_braces_mri_f13", "exclude_incidental_f05", "exclude_incidental_f09", "exclude_incidental_f13", "freesurfer_qc_f05", "freesurfer_qc_f09", "freesurfer_qc_f13",'f1200201', 'fever_tri1', 'fever_tri2', 'fever_tri3', 'flu_tri1', 'flu_tri2', 'flu_tri3', 'dermatitis_tri1', 'dermatitis_tri2', 'dermatitis_tri3', 'herpeszoster_tri1', 'herpeszoster_tri2', 'herpeszoster_tri3', 'jaundice_tri1', 'jaundice_tri2', 'jaundice_tri3', 'clean_upper_resp_inf_tri1', 'clean_upper_resp_inf_tri2', 'clean_upper_resp_inf_tri3', 'clean_uwi_tri1', 'clean_uwi_tri2', 'clean_uwi_tri3', 'clean_GI_inf_tri1', 'clean_GI_inf_tri2', 'clean_GI_inf_tri3', 'STD_tri1', 'STD_tri2', 'STD_tri3', 'clean_lower_resp_inf_tri1', 'clean_lower_resp_inf_tri2', 'clean_lower_resp_inf_tri3', 'IDM.y', 'NSAIDTOT', 'ABIOTOT', 'PMOLTOT', 'CORTTOT', 'MUCOTOT', 'COUGHTOT', 'd2300101', 'd2400101', 'd2500101', 'd1100101', 'preg_ic', 'ic')) 
df_final <- rename(df_final, 'IDM' = IDM.x)

# Removing columns i don't need for additional analyses requested by reviewers (same as prior df_final but now we do include all the individual infection domains)
df_final_revisions <- dplyr::select(df_tidy, -c("INTAKEPERIOD", "mri_consent_f05", "mri_consent_f09", "mri_consent_f13", "t1_has_nii_f05", "t1_has_nii_f09", "t1_has_nii_f13", "t1_asset_has_nii_f09", "has_braces_mri_f05", "has_braces_mri_f09", "has_braces_mri_f13", "exclude_incidental_f05", "exclude_incidental_f09", "exclude_incidental_f13", "freesurfer_qc_f05", "freesurfer_qc_f09", "freesurfer_qc_f13",'f1200201', 'IDM.y', 'NSAIDTOT', 'ABIOTOT', 'PMOLTOT', 'CORTTOT', 'MUCOTOT', 'COUGHTOT', 'd2300101', 'd2400101', 'd2500101', 'd1100101', 'preg_ic', 'ic')) 
df_final_revisions <- rename(df_final_revisions, 'IDM' = IDM.x)

# Transforming from wide to long format for df final 
df_tidy_long <- pivot_longer(df_final, -c(IDM, IDC, MOTHER, AGE_M_v2,ETHNMv2,SMOKE_ALL,GSI,sumscore_inf_tri1,sumscore_inf_tri2, sumscore_inf_tri3, sumscore_inf_tot, GENDER,EDUCM_3groups,INCOME,srs_weighted,SSRITOT,FOLIUM_VALIDATED,f1200101, cbcl_sum_5, include, mdrink_updated, cbcl_sum_14, sum_int_14, sum_ext_14, f1100101, post_life_events,post_direct_victimization, DIAB_GRA,PIH_v1, mat_inflam, im,preeclampsia), names_to = c("variable", "timepoint"), values_to = c("value"), names_pattern = "(.*)_(.*)") #transform to long format 

long_df_ordered <- pivot_wider(df_tidy_long, names_from = c("variable"), values_from = c("value")) #transform back to wide, in order to not have 1 column for each outcome but different columns per outcome, but the time points are stacked under each other

# Save full df with inclusion var (which can now be loaded in next step for imputation)
save(long_df_ordered, file = 'full_df.Rds') #save full df with inclusion var
load("full_df.Rds")

# Transforming from wide to long format for df final for revisions
df_tidy_long_revisions <- pivot_longer(df_final_revisions, -c(IDM, IDC, MOTHER, AGE_M_v2,ETHNMv2,SMOKE_ALL,GSI,sumscore_inf_tri1,sumscore_inf_tri2, sumscore_inf_tri3, sumscore_inf_tot, GENDER,EDUCM_3groups,INCOME,srs_weighted,SSRITOT,FOLIUM_VALIDATED,f1200101, cbcl_sum_5, include, mdrink_updated, cbcl_sum_14, sum_int_14, sum_ext_14, f1100101, post_life_events,post_direct_victimization, DIAB_GRA,PIH_v1, mat_inflam, im,preeclampsia, fever_tri1, fever_tri2, fever_tri3, flu_tri1, flu_tri2, flu_tri3, dermatitis_tri1, dermatitis_tri2, dermatitis_tri3, herpeszoster_tri1, herpeszoster_tri2, herpeszoster_tri3, jaundice_tri1, jaundice_tri2, jaundice_tri3, clean_upper_resp_inf_tri1, clean_upper_resp_inf_tri2, clean_upper_resp_inf_tri3, clean_uwi_tri1, clean_uwi_tri2, clean_uwi_tri3, clean_GI_inf_tri1, clean_GI_inf_tri2, clean_GI_inf_tri3, STD_tri1, STD_tri2, STD_tri3, clean_lower_resp_inf_tri1, clean_lower_resp_inf_tri2, clean_lower_resp_inf_tri3), names_to = c("variable", "timepoint"), values_to = c("value"), names_pattern = "(.*)_(.*)") 

long_df_ordered_revisions <- pivot_wider(df_tidy_long_revisions, names_from = c("variable"), values_from = c("value"))

save(long_df_ordered_revisions, file = 'full_df_revisions.Rds')

# Write df to csv for bash
write.csv(long_df_ordered, "full_df.csv", row.names = F)
write.csv(long_df_ordered_revisions, "full_df_revisions.csv", row.names = F)

#\

####################################################
#-----------Imaging demographics-------------------#
####################################################

## Calculate how many people have 1, 2 or 3 scans for each exposure
one_scan <- subset(df_final, (complete.cases(genr_tbv_f05) & is.na(genr_tbv_f09) & is.na(genr_tbv_f13)) | 
                     (complete.cases(genr_tbv_f09) & is.na(genr_tbv_f05) & is.na(genr_tbv_f13)) |
                     (complete.cases(genr_tbv_f13) & is.na(genr_tbv_f09) & is.na(genr_tbv_f13)))

two_scans <- subset(df_final, (complete.cases(genr_tbv_f09) & complete.cases(genr_tbv_f13)) | 
                      (complete.cases(genr_tbv_f05) & complete.cases(genr_tbv_f09)) | 
                      (complete.cases(genr_tbv_f05) & complete.cases(genr_tbv_f13)))

three_scans <- subset(df_final, complete.cases(genr_tbv_f05) & 
                        complete.cases(genr_tbv_f09) & 
                        complete.cases(genr_tbv_f13)) 

## Visualize age of each participant at each study time point 
long_df_ordered2 <- long_df_ordered %>%
  group_by(IDC) %>%
  mutate(scan_count = sum(!is.na(genr_tbv))) %>%
  mutate(scan_dummy = case_when(
    scan_count == 1 ~ 1,
    scan_count == 2 ~ 2,
    scan_count == 3 ~ 3,
    TRUE ~ NA_real_
  )) %>%
  ungroup()

long_df_ordered3 <- subset(long_df_ordered2, scan_count == 1 | scan_count == 2 | scan_count == 3)

long_df_ordered4 <- long_df_ordered3 %>%
  arrange(IDC, scan_count)  

ggplot(long_df_ordered4, aes(x = age_child_mri, y = scan_count, group = IDC, color = as.factor(scan_count))) +
  geom_line(position = "stack", size = 0.02) +
  geom_point(position = position_stack(vjust = 0.5), size = 0.1) +
  labs(x = 'Child age (years)', y = "", color = 'Number of scans') +
  theme_minimal() +
  theme(axis.text.y = element_blank())

# Making histogram of age per time point
hist1 <- ggplot(df_final, aes(x = age_child_mri_f05)) + geom_histogram(bins = 50, color = 'white', fill = '#33CC99') + theme_bw() + xlab('Child age range at T1') + ylab('Count') + theme(axis.title=element_text(face = 'bold', family = 'serif', size = 12)) + ylim(c(0,250))
hist2 <- ggplot(df_final, aes(x = age_child_mri_f09)) + geom_histogram(bins = 50, color = 'white', fill = '#CC6666') + theme_bw() + xlab('Child age range at T2') + ylab('Count') + theme(axis.title=element_text(face = 'bold', family = 'serif', size = 12)) + ylim(c(0,250))
hist3 <- ggplot(df_final, aes(x = age_child_mri_f13)) + geom_histogram(bins = 50, color = 'white', fill = '#66CCFF') + theme_bw() + xlab('Child age range at T3') + ylab('Count') + theme(axis.title=element_text(face = 'bold', family = 'serif', size = 12)) + ylim(c(0,250))

ggarrange(hist1, hist2, hist3, labels = 'AUTO', ncol = 1, nrow=3)

#\

############################################################
#-------------------Multiple imputation--------------------#
############################################################

# running random forest multiple imputation 
# working in the server now so not on local r studio
# So from bash load R
# Type the following in bash
# Module load R/4.0.2

setwd("Set_wd_to_path_on_sever")

library(mice)
library(ranger)

## imputations for df_final 
df <- read.csv('full_df.csv') 

imp0 <- mice(df, 
             maxit = 0, 
             method = 'rf') 

meth <- imp0$method

meth[c(1:3, 32:86)] <- "" 

pred <- imp0$predictorMatrix

pred[, c(1:3, 32:86)] <- 0 

visSeq <- imp0$visitSequence

imp.test <- mice(df, 
                 method = meth, 
                 predictorMatrix = pred, 
                 visitSequence = visSeq, 
                 maxit = 50, 
                 m = 30, 
                 printFlag = T, 
                 seed = 2023) 

save(imp.test, file = 'set_name_file1.Rds')

## imputations for df_final_revisions
colnames(long_df_ordered_revisions)

imp0_revisions <- mice(long_df_ordered_revisions, 
             maxit = 0, 
             method = 'rf') 

meth_revisions <- imp0_revisions$method

meth_revisions[c(1:3, 18:47, 62:116)] <- "" 

pred_revisions <- imp0_revisions$predictorMatrix

pred_revisions[, c(1:3, 18:47, 62:116)] <- 0 

visSeq_revisions <- imp0_revisions$visitSequence

imp.test_revisions <- mice(long_df_ordered_revisions, 
                 method = meth_revisions, 
                 predictorMatrix = pred_revisions, 
                 visitSequence = visSeq_revisions, 
                 maxit = 50, 
                 m = 30, 
                 printFlag = T, 
                 seed = 2024) 

save(imp.test_revisions, file = 'set_name_file2.Rds')

#\ 

############################################
#-------------------IPW--------------------#
###########################################

### Opening libraries we need & setting wd 
# go back to local R studio                        
# Clears  environment

rm(list = ls()) 

# Set wd 
setwd('Set_path_to_wd')
wd <- getwd()

# Open libraries we need 
libraries <- c('foreign', 'haven', 'dplyr', 'tidyverse' ,'mice', 'lme4', 'lattice', 'ggplot2', 'ggeffects', 'splines', 'reshape2', 'gridExtra', 'grid', 'ggpubr', 'tidyr', 'broom.mixed', 'xlsx', 'corrplot', 'stringi', 'miceadds', 'mitools', 'CBPS', 'survey', 'survival')

invisible(lapply(libraries, require, character.only = T))

# Open imputed full df 
load('no_outcome_imp.Rds') 
load('no_outcome_imp_revisions.rds')

#-------------------PART 1--------------------#
# Part 1: calculating weights for baseline

### Creating weights for IPW
incl <- "include" #inclusion variable 
pred_baseline <- c("AGE_M_v2","ETHNMv2","SMOKE_ALL","GSI","sumscore_inf_tot","GENDER","EDUCM_3groups","INCOME","SSRITOT")  

# Use the reformulate function to create a formula
fit_eq <- reformulate(termlabels = pred_baseline, response = incl)
print(fit_eq)

# Create empty list for storing propensity scores after running model on each imputed set
ps <- list()
ps2 <- list()

# For loop to run CBPS model over each imputed set and get ps's
for (imp in 1:imp.test$m){
  imp_ipw.imp <- complete(imp.test, imp)
  #get weights for weighting back to baseline
  fit.imp <-  CBPS((fit_eq), data=imp_ipw.imp) 
  #propensity scores 
  ps[[imp]] <- fit.imp$fitted.values 
}

for (imp in 1:imp.test_revisions$m){
  imp_ipw.imp2 <- complete(imp.test_revisions, imp)
  fit.imp2 <-  CBPS((fit_eq), data=imp_ipw.imp2) 
  ps2[[imp]] <- fit.imp2$fitted.values 
}

# Use long fomat to paste ps's in imputed sets
long <- complete(imp.test, action = "long", include =TRUE)
len <- lengths(ps)

long2 <- complete(imp.test_revisions, action = 'long', include = TRUE)
len2 <- lengths(ps2)

#here we make a new df with a vector .imp en .id and each propensity score so we can merge it with our long df 
df <- data.frame(.imp = rep(seq_along(len), len), .id=sequence(len), ps=unlist(ps))
merge <- merge(long, df, by=c(".imp", ".id"), all.x=TRUE) 

df2 <- data.frame(.imp = rep(seq_along(len2), len2), .id=sequence(len2), ps2=unlist(ps2))
merge2 <- merge(long2, df2, by=c(".imp", ".id"), all.x=TRUE) 

# Calculate weights from ps's
merge2 <- mutate(merge, baseline_weights = ifelse(get(incl) == 1, (1 / ps), (1 / (1-ps)))) #you create weights because individuals with a lower probability of being included in the study would need a higher weights; hence we take the inverse of the probability 

merge3 <- mutate(merge2, baseline_weights = ifelse(get(incl) == 1, (1 / ps2), (1 / (1-ps2)))) 

#\ 

#-------------------PART 2--------------------#
# Part 2: calculating weights for 6

### Creating weights for IPW
incl <- "include" #inclusion variable 
pred_6 <- c("cbcl_sum_5","srs_weighted") #predictors of wave 6 

# Use the reformulate function to create a formula
fit_eq6 <- reformulate(termlabels = pred_6, response = incl)
print(fit_eq6)

# create empty list for storing propensity scores after running model on each imputed set
ps6 <- list()
ps6.2 <- list()

# For loop to run CBPS model over each imputed set and get ps's
for (imp in 1:imp.test$m){
  imp_ipw.imp6 <- complete(imp.test, imp)
  #get weights for weighting back to f6
  fit.imp6 <-  CBPS((fit_eq6), data=imp_ipw.imp6) 
  ps6[[imp]] <- fit.imp6$fitted.values #propensity scores (chance that ID is in the sample)
}

for (imp in 1:imp.test_revisions$m){
  imp_ipw.imp6.2 <- complete(imp.test_revisions, imp)
  fit.imp6.2 <-  CBPS((fit_eq6), data=imp_ipw.imp6.2) 
  ps6.2[[imp]] <- fit.imp6.2$fitted.values 
}

# Use long fomat to paste ps's in imputed sets
long6 <- complete(imp.test, action = "long", include =TRUE)
len6 <- lengths(ps6)

long6.2 <- complete(imp.test_revisions, action = 'long', include = TRUE)
len6.2 <- lengths(ps6.2)

#here we make a new df with a vector .imp en .id and each propensity score so we can merge it with our long df 
df6 <- data.frame(.imp = rep(seq_along(len6), len6), .id=sequence(len6), ps6=unlist(ps6))  
merge6 <- merge(long6, df6, by=c(".imp", ".id"), all.x=TRUE)

df6.2 <- data.frame(.imp = rep(seq_along(len6.2), len6.2), .id=sequence(len6.2), ps6.2=unlist(ps6.2))
merge6.2 <- merge(long6.2, df6.2, by = c('.imp', '.id'), all.x = TRUE)

# Calculate weights from ps's
merge4 <- mutate(merge6, f6_weights = ifelse(get(incl) == 1, (1 / ps6), (1 / (1-ps6)))) 
merge5 <- mutate(merge6.2, f6_weights = ifelse(get(incl) == 1, (1/ps6.2), (1/(1-ps6.2))))

#\

#-------------------PART 3--------------------#
# Part 3: creating combined weights by multiplying baseline and 6

# Merge df with weights baseline and 6
total_merge <- merge(merge2, merge4, by=c(1:75), all.x=TRUE)

total_merge2 <- merge(merge3, merge5, by=c(1:117), all.x = TRUE)

# Creating combined weights by multiplying baseline and f6
total_merge$combined_weights <- ifelse(is.na(total_merge$baseline_weights), NA, total_merge$baseline_weights * total_merge$f6_weights)

total_merge2$combined_weights <- ifelse(is.na(total_merge2$baseline_weights), NA, total_merge2$baseline_weights * total_merge2$f6_weights)

# Convert to mids
ipw_mids <- as.mids(total_merge)

ipw_mids2 <- as.mids(total_merge2)

# Create final analysis set with weights by selecting only in observations with inclusionvar = 1
ipw_mids_final <- filter(ipw_mids, get(incl) == 1)

ipw_mids_final_revisions <- filter(ipw_mids2, get(incl) == 1)

# save dataframes 
save(ipw_mids, file = 'no_outcome_imp_ipw_mids.RDS')
save(ipw_mids_final, file = 'no_outcome_imp_ipw_mids_final.RDS')

save(ipw_mids2, file = 'no_outcome_imp_ipw_mids_revisions.RDS')
save(ipw_mids_final_revisions, file = 'no_outcome_imp_ipw_mids_final_revisions.RDS')

#\

############################################
#----------Complete case analysis---------#
###########################################

### Conduct complete cases analysis across 3 timepoints for 5 fdr sign regions 

## select df
complete_df <- subset(df_final, complete.cases(genr_tbv_f05) & 
                        complete.cases(genr_tbv_f09) & 
                        complete.cases(genr_tbv_f13)) 

complete_df <- dplyr::select(complete_df, c('IDC', 'AGE_M_v2', 'ETHNMv2', 'SMOKE_ALL', 'GSI', 'sumscore_inf_tot', 'sumscore_inf_tri1', 'sumscore_inf_tri2', 'sumscore_inf_tri3', 'GENDER', 'EDUCM_3groups', 'middletemporal_vol_cortical_f05', 'middletemporal_vol_cortical_f09', 'middletemporal_vol_cortical_f13', 'parsorbitalis_vol_cortical_f05', 'parsorbitalis_vol_cortical_f09', 'parsorbitalis_vol_cortical_f13', 'rostralanteriorcingulate_vol_cortical_f05', 'rostralanteriorcingulate_vol_cortical_f09', 'rostralanteriorcingulate_vol_cortical_f13', 'superiorfrontal_vol_cortical_f05', 'superiorfrontal_vol_cortical_f09', 'superiorfrontal_vol_cortical_f13', 'temporalpole_vol_cortical_f05', 'temporalpole_vol_cortical_f09', 'temporalpole_vol_cortical_f13'))

## create long format 
df_tidy_long_complete <- pivot_longer(complete_df, -c(IDC, AGE_M_v2,ETHNMv2,SMOKE_ALL,GSI,sumscore_inf_tri1,sumscore_inf_tri2, sumscore_inf_tri3, sumscore_inf_tot, GENDER,EDUCM_3groups), names_to = c("variable", "timepoint"), values_to = c("value"), names_pattern = "(.*)_(.*)") 

long_df_ordered_complete <- pivot_wider(df_tidy_long_complete, names_from = c("variable"), values_from = c("value"))

## standardize vars 
fdr_outcomes <- c('middletemporal_vol_cortical', 'parsorbitalis_vol_cortical', 'rostralanteriorcingulate_vol_cortical', 'superiorfrontal_vol_cortical', 'temporalpole_vol_cortical')

for (x in fdr_outcomes){
  t <- long_df_ordered_complete[x]
  colname <- paste0(colnames(t), "_standardized")
  long_df_ordered_complete[,colname] <- as.numeric(scale(t))
}

fdr_outcomes2 <- paste0(fdr_outcomes, '_standardized')

## Repeat analysis for 5 sign regions for complete cases of 3 scans
library(lmerTest)

long_df_ordered_complete$timepoint <- factor(long_df_ordered_complete$timepoint, levels = c("f05", "f09", "f13"))

long_df_ordered_complete$timepoint <- as.numeric(long_df_ordered_complete$timepoint)

complete_cases_results <- data.frame()

for(x in fdr_outcomes2) {
  f <- paste0(x, "~ scale(sumscore_inf_tri1)*timepoint + EDUCM_3groups + GENDER + AGE_M_v2 + ETHNMv2 + GSI + SMOKE_ALL + (1|IDC)") 
  g <- paste0(x, "~ scale(sumscore_inf_tri2)*timepoint + EDUCM_3groups + GENDER + AGE_M_v2 + ETHNMv2 + GSI + SMOKE_ALL + (1|IDC)") 
  h <- paste0(x, "~ scale(sumscore_inf_tri3)*timepoint + EDUCM_3groups + GENDER + AGE_M_v2 + ETHNMv2 + GSI + SMOKE_ALL + (1|IDC)") 
  
  #calculating beta, SE, confidence interval and pval for each model
  bval1 <- summary(lmer(as.formula(f), data = long_df_ordered_complete))$coefficients[12,'Estimate']
  seval1 <- summary(lmer(as.formula(f), data = long_df_ordered_complete))$coefficients[12,'Std. Error']
  pval1 <- summary(lmer(as.formula(f), data = long_df_ordered_complete))$coefficients[12,'Pr(>|t|)']
  
  bval2 <- summary(lmer(as.formula(g), data = long_df_ordered_complete))$coefficients[12,'Estimate']
  seval2 <- summary(lmer(as.formula(g), data = long_df_ordered_complete))$coefficients[12,'Std. Error']
  pval2 <- summary(lmer(as.formula(g), data = long_df_ordered_complete))$coefficients[12,'Pr(>|t|)']
  
  bval3 <- summary(lmer(as.formula(h), data = long_df_ordered_complete))$coefficients[12,'Estimate']
  seval3 <- summary(lmer(as.formula(h), data = long_df_ordered_complete))$coefficients[12,'Std. Error']
  pval3 <- summary(lmer(as.formula(h), data = long_df_ordered_complete))$coefficients[12,'Pr(>|t|)']

  #storing values in empty datafile  
  complete_cases_results[x,1] <- bval1
  complete_cases_results[x,2] <- seval1
  complete_cases_results[x,3] <- pval1
  complete_cases_results[x,4] <- bval2
  complete_cases_results[x,5] <- seval2
  complete_cases_results[x,6] <- pval2
  complete_cases_results[x,7] <- bval3
  complete_cases_results[x,8] <- seval3
  complete_cases_results[x,9] <- pval3
  
  #assigning names to columns 
  colnames(complete_cases_results) <- c("bval_tri1", 'seval_tri1', 'pval_tri1',"bval_tri2", 'seval_tri2', 'pval_tri2', "bval_tri3", 'seval_tri3', 'pval_tri3')
}

write.xlsx(complete_cases_results, "complete_cases_results.xlsx")

#\

#\ END OF THIS SCRIPT, GO TO PART B TO RUN ANALYSES 
