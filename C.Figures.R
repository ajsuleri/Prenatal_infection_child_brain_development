#------------------------------------------------------------------------#
#############PROJECT: infection & child brain development#################
#------------------------------------------------------------------------#

# Author: Anna Suleri

### Goal of this script: make a flowchart and baseline table & figures to depict results
# Of note, other descripte figures for infection sum score e.g. are made in script A 

### Start script by setting wd & loading libraries 
setwd('V:\\medewerkers\\039088 Suleri, A\\Projecten\\05. Project_normative_development\\Data')
wd <- getwd()

libraries <- c('foreign', 'haven', 'dplyr', 'tidyverse' ,'mice', 'lme4', 'lattice', 'ggplot2', 'ggeffects', 'splines', 'reshape2', 'gridExtra', 'grid', 'ggpubr', 'tidyr', 'broom.mixed', 'corrplot', 'stringi', 'miceadds', 'mitools', 'survival')

invisible(lapply(libraries, require, character.only = T))

load('df_tidy.Rds')

###---Descriptive figure & table & non-response analysis---###

### Flowchart
#Inclusion criteria 
inclu1 <- subset(df_tidy, complete.cases(genr_tbv_f05) | complete.cases(genr_tbv_f09) | complete.cases(genr_tbv_f13)) #one of the time points should have complete cases

inclu2 <- subset(inclu1, INTAKEPERIOD == '<18 weeks') #including only mothers who enrolled in tri1, n=7145
inclu3 <- subset(inclu2, complete.cases(sumscore_inf_tot)) #no missing information in prenatal infection, n=4259

#Exclusion criteria
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

#\ 

### Baseline table
baselinevars <- c("AGE_M_v2", "sumscore_inf_tot","ETHNMv2","EDUCM_3groups", "INCOME", "SMOKE_ALL","GENDER")

for(i in baselinevars){ #for each baseline variable as specified in baselinevars decide whether cont or cat and give output accordingly 
  #x = i vars that are columns in the dataframe df
  x <- df_final[, i] 
  #show column name as heading per output 
  message(i) 
  #function for continuous variables
  summary_continuous <- function(x){
    standev <- sd(x, na.rm = T)
    meanvar <- mean(x, na.rm = T)
    print(paste0(round(meanvar, 1), '(', round(standev, 1), ')'))
  }
  #function for categorical variables 
  summary_categorical <- function(x){
    tab1 <- prop.table(table(x, useNA = 'always'))
    tab2 <- table(x, useNA = "always")
    print(paste(round(tab1 * 100, 1), '%', names(tab1), collapse = ','))
    print(paste(tab2, names(tab2)))
  }
  #if else to apply correct function for vars type 
  if (class(x) == 'numeric') {
    summary_continuous(x)
  } 
  else 
  {
    summary_categorical(x) 
  }
}

#\

### Non-response analysis
nr_vars <- c("AGE_M_v2","ETHNMv2","SMOKE_ALL","GSI","sumscore_inf_tot","GENDER","EDUCM_3groups","INCOME","SSRITOT", "cbcl_sum_5","srs_weighted")

for (i in nr_vars){
  #define df with excluded sample (so who are in genR but not in study sample)
  df_NR <- df_tidy[!(df_tidy$IDC %in% df_final$IDC),]
  #for every loop call name i 
  cat(i)
  #if cat vars then chi square test else t test 
  if (class(df_final[,i]) == 'factor') {
    print(chisq.test(cbind(table(df_NR[,i]), table(df_final[,i]))))
  } 
  else 
  {
    print(t.test(df_NR[,i], df_final[,i], paired = F))
  }
}

#\ 

###---Now we make plots for the results---###

### load libraries 
libraries <- c('foreign', 'haven', 'dplyr', 'tidyverse' ,'mice', 'lme4', 'lattice', 'ggplot2', 'ggeffects', 'splines', 'reshape2', 'gridExtra', 'grid', 'ggpubr', 'tidyr', 'broom.mixed',  'corrplot', 'stringi', 'miceadds', 'mitools', 'CBPS', 'survey', 'survival', 'ggseg', "patchwork", "openxlsx", 'jtools', 'interactions', 'corrplot', 'cowplot', 'ggcorrplot', 'RColorBrewer')

invisible(lapply(libraries, require, character.only = T))

### set wd 
setwd('V:\\medewerkers\\039088 Suleri, A\\Projecten\\05. Project_normative_development\\Data')

load('no_outcome_imp_ipw_mids_final.RDS')

setwd('V:\\medewerkers\\039088 Suleri, A\\Projecten\\05. Project_normative_development\\Figures')

### select single df out of imputed df for figures 
imp.test_long <- complete(ipw_mids_final, include = T, action = "long")

imp.test_long$timepoint <- as.factor(imp.test_long$timepoint) #making time continuous 
imp.test_long$timepoint <- as.numeric(imp.test_long$timepoint)
imp.test_long$sumscore_inf_tot <- as.numeric(imp.test_long$sumscore_inf_tot)
imp.test_long$sumscore_inf_tot_standardized <- as.numeric(scale(imp.test_long$sumscore_inf_tot))
imp.test_long <- as.mids(imp.test_long)

single_df <- complete(imp.test_long, 1)
single_df2<- data.frame(single_df)

### Correlation plot
corr_vars <- dplyr::select(single_df2, c(8:11, 22:29, 31, 33:74)) #removing covars
corr_vars2 <- dplyr::select(corr_vars, -c(6:10, 11:12, 17)) #removing brain vars we dont include
corr_vars3 <- rename(corr_vars2, 'Prenatal infection sum score (trimester 1)' = sumscore_inf_tri1, 'Prenatal infection sum score (trimester 2)' = sumscore_inf_tri2, 'Prenatal infection sum score (trimester 3)' = sumscore_inf_tri3, 'Prenatal infection sum score (whole pregnancy)' = sumscore_inf_tot, 'Brain stem' = Brain_Stem_vol, 'Intracranial volume' = eTIV, 'Cerebellum' = Cerebellum_Cortex_vol_subcortical, 'Amygdala' = Amygdala_vol_subcortical, 'Hippocampus' = Hippocampus_vol_subcortical, 'Caudate' = Caudate_vol_subcortical, 'Putamen' = Putamen_vol_subcortical, 'Thalamus' = Thalamus_Proper_vol_subcortical, 'Pallidum' = Pallidum_vol_subcortical, 'Banks of the superior temporal sulcus' = bankssts_vol_cortical, 'Caudal anterior cingulate' = caudalanteriorcingulate_vol_cortical, 'Caudal middle frontal' = caudalmiddlefrontal_vol_cortical, 'Cuneus' = cuneus_vol_cortical, 'Entorhinal' = entorhinal_vol_cortical, 'Fusiform' = fusiform_vol_cortical, 'Inferior parietal' = inferiorparietal_vol_cortical, 'Inferior temporal' = inferiortemporal_vol_cortical, 'Isthmus cingulate' = isthmuscingulate_vol_cortical, 'Lateral occipital' = lateraloccipital_vol_cortical, 'Lateral orbitofrontal' = lateralorbitofrontal_vol_cortical, 'Lingual' = lingual_vol_cortical, 'Medial orbitofrontal' = medialorbitofrontal_vol_cortical, 'Middle temporal' = middletemporal_vol_cortical, 'Parahippocampal' = parahippocampal_vol_cortical, 'Paracentral' = paracentral_vol_cortical, 'Pars opercularis' = parsopercularis_vol_cortical, 'Pars orbitalis' = parsorbitalis_vol_cortical, 'Pars triangularis' = parstriangularis_vol_cortical, 'Pericalcarine' = pericalcarine_vol_cortical, 'Postcentral' = postcentral_vol_cortical, 'POsterior cingulate' = posteriorcingulate_vol_cortical, 'Precentral' = precentral_vol_cortical, 'Precuneus' = precuneus_vol_cortical, 'Rostral anterior cingulate' = rostralanteriorcingulate_vol_cortical, 'Rostral middle frontal' = rostralmiddlefrontal_vol_cortical, "Superior frontal" = superiorfrontal_vol_cortical, 'Superior parietal' = superiorparietal_vol_cortical, 'Superior temporal' = superiortemporal_vol_cortical, 'Supramarginal' = supramarginal_vol_cortical, 'Frontal pole' = frontalpole_vol_cortical, 'Temporal pole' = temporalpole_vol_cortical, 'Transverse temporal' = transversetemporal_vol_cortical, 'Insula' = insula_vol_cortical.x)
corr_vars4 <- as.matrix(corr_vars3)

correlation <- cor(corr_vars4, use="pairwise.complete.obs")

corrplot::corrplot(correlation, method = 'color', order = 'FPC', type = 'lower',diag = F, tl.col = 'black', tl.cex = 0.7, col=brewer.pal(n=8, name="PuOr"), tl.srt = 45)

#\ 

### Show brain anatomy plot (all based on standardized beta's for interaction term)
#Create  relevant dataframes
sign_results_df_wholepreg_cortex <- data.frame(
  region = c('temporal pole'),
  beta = c(0.06))

sign_results_df_wholepreg_subcortex <- data.frame(
  region = c('cerebellum cortex', 'thalamus proper'),
  beta = c(0.038, -0.035))

sign_results_df_tri1_cortex <- data.frame(
  region = c('superior parietal', 'temporal pole'),
  beta = c(0.041,0.058))

sign_results_df_tri2_cortex <- data.frame(
  region = c('bankssts', 'caudal anterior cingulate', 'inferior parietal','pars opercularis','pars triangularis'),
  beta = c(-0.032, -0.05, -0.038, -0.037, -0.048))

sign_results_df_tri2_subcortex <- data.frame(
  region = c('cerebellum cortex', 'thalamus proper'),
  beta = c(0.037, -0.03))

sign_results_df_tri3_cortex <- data.frame(
  region = c('caudal middle frontal', 'fusiform', 'inferior temporal', 'lateral orbitofrontal', 'middle temporal', 'paracentral','pars orbitalis','rostral anterior cingulate', 'rostral middle frontal' ,'superior frontal', 'superior parietal' ,'temporal pole'),
  beta = c(0.045, 0.04, 0.04, 0.064, 0.064, 0.036,0.073, 0.073, 0.043, 0.059,0.051,0.076))

sign_results_df_tri3_subcortex <- data.frame(
  region = c('caudate'),
  beta = c(0.52))

#Create  4 panels for all sign cortex results
sign_results_df_wholepreg_cortex$Exposure <- 'Total pregnancy'
sign_results_df_tri1_cortex$Exposure <- 'Trimester 1'
sign_results_df_tri2_cortex$Exposure <- 'Trimester 2'
sign_results_df_tri3_cortex$Exposure <- 'Trimester 3'

sign_results_df_cortex <- rbind(sign_results_df_wholepreg_cortex, sign_results_df_tri1_cortex, sign_results_df_tri2_cortex, sign_results_df_tri3_cortex)

a <- sign_results_df_cortex %>%
  group_by(Exposure) %>%
  ggplot() +
  geom_brain(atlas = dk,
             position = position_brain(hemi ~ side),
             aes(fill = beta)) +
  facet_wrap(~Exposure) +
  theme_bw() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  theme(strip.background = element_rect(fill = 'grey83')) +
  theme(strip.text = element_text(colour = 'black', face = 'bold')) +
  theme(plot.title = element_text(size = 12, family = 'sans')) + 
  theme(legend.position = 'bottom',
        legend.text = element_text(size = 8, angle =45)) +
  scale_fill_gradient(low = '#330099', high = 'dark red') + 
  guides(fill = guide_colorbar(barheight = 1.5)) + 
  labs(fill = 'Beta coefficient')

#Create 3 panels for all sign subcortex results
sign_results_df_wholepreg_subcortex$Exposure <- 'Total pregnancy'
sign_results_df_tri2_subcortex$Exposure <- 'Trimester 2'
sign_results_df_tri3_subcortex$Exposure <- 'Trimester 3'

sign_results_df_subcortex <- rbind(sign_results_df_wholepreg_subcortex, sign_results_df_tri2_subcortex, sign_results_df_tri3_subcortex)

b <- sign_results_df_subcortex %>%
  group_by(Exposure) %>%
  ggplot() +
  geom_brain(atlas = aseg,
             aes(fill = beta)) +
  facet_wrap(~Exposure) +
  theme_bw() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  theme(strip.background = element_rect(fill = 'grey83')) +
  theme(strip.text = element_text(colour = 'black', face = 'bold')) +
  theme(plot.title = element_text(size = 12, family = 'sans')) + 
  theme(legend.position = 'bottom',
        legend.text = element_text(size = 8, angle =45)) +
  scale_fill_gradient(low = '#330099', high = 'dark red') + 
  guides(fill = guide_colorbar(barheight = 1.5)) + 
  labs(fill = 'Beta coefficient')

#Combing pval <0.05 sign plots for cortex and subcortex
plot_grid(a, b, labels = c('A', 'B'), ncol = 2)

### make longitudinal plots for nominal sign cortex and subcortex regions
single_df2$'Prenatal infection trimester 3' <- single_df2$sumscore_inf_tri3
single_df2$'Prenatal infection trimester 2' <- single_df2$sumscore_inf_tri2
single_df2$'Prenatal infection trimester 1' <- single_df2$sumscore_inf_tri1
single_df2$'Prenatal infection (total pregnancy)' <- single_df2$sumscore_inf_tot

# total preg 
totinf_plot1 <- lmer(temporalpole_vol_cortical.x ~ `Prenatal infection (total pregnancy)`*timepoint + AGE_M_v2 + GENDER + ETHNMv2 + EDUCM_3groups + GSI + SMOKE_ALL + (1|IDC), data = single_df2, weights = combined_weights) 
totinf_plot2 <- lmer(Cerebellum_Cortex_vol_subcortical ~ `Prenatal infection (total pregnancy)`*timepoint + AGE_M_v2 + GENDER + ETHNMv2 + EDUCM_3groups + GSI + SMOKE_ALL + (1|IDC), data = single_df2, weights = combined_weights) 
totinf_plot3 <- lmer(Thalamus_Proper_vol_subcortical ~ `Prenatal infection (total pregnancy)`*timepoint + AGE_M_v2 + GENDER + ETHNMv2 + EDUCM_3groups + GSI + SMOKE_ALL + (1|IDC), data = single_df2, weights = combined_weights) 

p1 <- interact_plot(totinf_plot1, pred = timepoint, modx = `Prenatal infection (total pregnancy)`,modx.values = 'plus-minus', x.label = 'Timepoint', y.label = 'Temporal pole', interval = T, plot.points = F,point.size = 0.3, point.alpha = 0.1, line.thickness = 1, jitter = 0.2, int.type = 'confidence', colors = 'Qual1', vary.lty = F) + theme(panel.background = element_rect(fill = 'white', colour = 'black', size = 0.5))
p2 <- interact_plot(totinf_plot2, pred = timepoint, modx = `Prenatal infection (total pregnancy)`,modx.values = 'plus-minus', x.label = 'Timepoint', y.label = 'Cerebellum cortex', interval = T, plot.points = F,point.size = 0.3, point.alpha = 0.1, line.thickness = 1, jitter = 0.2, int.type = 'confidence', colors = 'Qual1', vary.lty = F) + theme(panel.background = element_rect(fill = 'white', colour = 'black', size = 0.5))
p3 <- interact_plot(totinf_plot3, pred = timepoint, modx = `Prenatal infection (total pregnancy)`,modx.values = 'plus-minus', x.label = 'Timepoint', y.label = 'Thalamus', interval = T, plot.points = F,point.size = 0.3, point.alpha = 0.1, line.thickness = 1, jitter = 0.2, int.type = 'confidence', colors = 'Qual1', vary.lty = F) + theme(panel.background = element_rect(fill = 'white', colour = 'black', size = 0.5))

ggarrange(p1, p2, p3, common.legend = T)

# trimester 1
tri1_plot1 <- lmer(superiorparietal_vol_cortical.x ~ `Prenatal infection trimester 1`*timepoint + AGE_M_v2 + GENDER + ETHNMv2 + EDUCM_3groups + GSI + SMOKE_ALL + (1|IDC), data = single_df2, weights = combined_weights) 
tri1_plot2 <- lmer(temporalpole_vol_cortical.x ~ `Prenatal infection trimester 1`*timepoint + AGE_M_v2 + GENDER + ETHNMv2 + EDUCM_3groups + GSI + SMOKE_ALL + (1|IDC), data = single_df2, weights = combined_weights) 

p4 <- interact_plot(tri1_plot1, pred = timepoint, modx = `Prenatal infection trimester 1`,modx.values = 'plus-minus', x.label = 'Timepoint', y.label = 'Superior parietal', interval = T, plot.points = F,point.size = 0.3, point.alpha = 0.1, line.thickness = 1, jitter = 0.2, int.type = 'confidence', colors = 'Qual1', vary.lty = F) + theme(panel.background = element_rect(fill = 'white', colour = 'black', size = 0.5))
p5 <- interact_plot(tri1_plot2, pred = timepoint, modx = `Prenatal infection trimester 1`,modx.values = 'plus-minus', x.label = 'Timepoint', y.label = 'Temporal pole', interval = T, plot.points = F,point.size = 0.3, point.alpha = 0.1, line.thickness = 1, jitter = 0.2, int.type = 'confidence', colors = 'Qual1', vary.lty = F) + theme(panel.background = element_rect(fill = 'white', colour = 'black', size = 0.5))

ggarrange(p4, p5, common.legend = T)

# trimester 2
tri2_plot1 <- lmer(bankssts_vol_cortical ~ `Prenatal infection trimester 2`*timepoint + AGE_M_v2 + GENDER + ETHNMv2 + EDUCM_3groups + GSI + SMOKE_ALL + (1|IDC), data = single_df2, weights = combined_weights) 
tri2_plot2 <- lmer(caudalanteriorcingulate_vol_cortical ~ `Prenatal infection trimester 2`*timepoint + AGE_M_v2 + GENDER + ETHNMv2 + EDUCM_3groups + GSI + SMOKE_ALL + (1|IDC), data = single_df2, weights = combined_weights) 
tri2_plot3 <- lmer(inferiorparietal_vol_cortical ~ `Prenatal infection trimester 2`*timepoint + AGE_M_v2 + GENDER + ETHNMv2 + EDUCM_3groups + GSI + SMOKE_ALL + (1|IDC), data = single_df2, weights = combined_weights) 
tri2_plot4 <- lmer(parsopercularis_vol_cortical ~ `Prenatal infection trimester 2`*timepoint + AGE_M_v2 + GENDER + ETHNMv2 + EDUCM_3groups + GSI + SMOKE_ALL + (1|IDC), data = single_df2, weights = combined_weights) 
tri2_plot5 <- lmer(parstriangularis_vol_cortical ~ `Prenatal infection trimester 2`*timepoint + AGE_M_v2 + GENDER + ETHNMv2 + EDUCM_3groups + GSI + SMOKE_ALL + (1|IDC), data = single_df2, weights = combined_weights) 
tri2_plot6 <- lmer(Cerebellum_Cortex_vol_subcortical ~ `Prenatal infection trimester 2`*timepoint + AGE_M_v2 + GENDER + ETHNMv2 + EDUCM_3groups + GSI + SMOKE_ALL + (1|IDC), data = single_df2, weights = combined_weights) 
tri2_plot7 <- lmer(Thalamus_Proper_vol_subcortical ~ `Prenatal infection trimester 2`*timepoint + AGE_M_v2 + GENDER + ETHNMv2 + EDUCM_3groups + GSI + SMOKE_ALL + (1|IDC), data = single_df2, weights = combined_weights) 

p6 <- interact_plot(tri2_plot1, pred = timepoint, modx = `Prenatal infection trimester 2`,modx.values = 'plus-minus', x.label = 'Timepoint', y.label = 'Banks of the superior temporal sulcus', interval = T, plot.points = F,point.size = 0.3, point.alpha = 0.1, line.thickness = 1, jitter = 0.2, int.type = 'confidence', colors = 'Qual1', vary.lty = F) + theme(panel.background = element_rect(fill = 'white', colour = 'black', size = 0.5))
p7 <- interact_plot(tri2_plot2, pred = timepoint, modx = `Prenatal infection trimester 2`,modx.values = 'plus-minus', x.label = 'Timepoint', y.label = 'Caudal anterior cingulate', interval = T, plot.points = F,point.size = 0.3, point.alpha = 0.1, line.thickness = 1, jitter = 0.2, int.type = 'confidence', colors = 'Qual1', vary.lty = F) + theme(panel.background = element_rect(fill = 'white', colour = 'black', size = 0.5))
p8 <- interact_plot(tri2_plot3, pred = timepoint, modx = `Prenatal infection trimester 2`,modx.values = 'plus-minus', x.label = 'Timepoint', y.label = 'Inferior parietal', interval = T, plot.points = F,point.size = 0.3, point.alpha = 0.1, line.thickness = 1, jitter = 0.2, int.type = 'confidence', colors = 'Qual1', vary.lty = F) + theme(panel.background = element_rect(fill = 'white', colour = 'black', size = 0.5))
p9 <- interact_plot(tri2_plot4, pred = timepoint, modx = `Prenatal infection trimester 2`,modx.values = 'plus-minus', x.label = 'Timepoint', y.label = 'Pars opercularis', interval = T, plot.points = F,point.size = 0.3, point.alpha = 0.1, line.thickness = 1, jitter = 0.2, int.type = 'confidence', colors = 'Qual1', vary.lty = F) + theme(panel.background = element_rect(fill = 'white', colour = 'black', size = 0.5))
p10 <- interact_plot(tri2_plot5, pred = timepoint, modx = `Prenatal infection trimester 2`,modx.values = 'plus-minus', x.label = 'Timepoint', y.label = 'Pars triangularis', interval = T, plot.points = F,point.size = 0.3, point.alpha = 0.1, line.thickness = 1, jitter = 0.2, int.type = 'confidence', colors = 'Qual1', vary.lty = F) + theme(panel.background = element_rect(fill = 'white', colour = 'black', size = 0.5))
p11 <- interact_plot(tri2_plot6, pred = timepoint, modx = `Prenatal infection trimester 2`,modx.values = 'plus-minus', x.label = 'Timepoint', y.label = 'Cerebellum cortex', interval = T, plot.points = F,point.size = 0.3, point.alpha = 0.1, line.thickness = 1, jitter = 0.2, int.type = 'confidence', colors = 'Qual1', vary.lty = F) + theme(panel.background = element_rect(fill = 'white', colour = 'black', size = 0.5))
p12 <- interact_plot(tri2_plot7, pred = timepoint, modx = `Prenatal infection trimester 2`,modx.values = 'plus-minus', x.label = 'Timepoint', y.label = 'Thalamus', interval = T, plot.points = F,point.size = 0.3, point.alpha = 0.1, line.thickness = 1, jitter = 0.2, int.type = 'confidence', colors = 'Qual1', vary.lty = F) + theme(panel.background = element_rect(fill = 'white', colour = 'black', size = 0.5))

ggarrange(p6, p7,p8,p9,p10,p11,p12, common.legend = T)

# trimester 3
tri3_plot1 <- lmer(caudalmiddlefrontal_vol_cortical ~ `Prenatal infection trimester 3`*timepoint + AGE_M_v2 + GENDER + ETHNMv2 + EDUCM_3groups + GSI + SMOKE_ALL + (1|IDC), data = single_df2, weights = combined_weights) 
tri3_plot2 <- lmer(fusiform_vol_cortical ~ `Prenatal infection trimester 3`*timepoint + AGE_M_v2 + GENDER + ETHNMv2 + EDUCM_3groups + GSI + SMOKE_ALL + (1|IDC), data = single_df2, weights = combined_weights) 
tri3_plot3 <- lmer(inferiortemporal_vol_cortical ~ `Prenatal infection trimester 3`*timepoint + AGE_M_v2 + GENDER + ETHNMv2 + EDUCM_3groups + GSI + SMOKE_ALL + (1|IDC), data = single_df2, weights = combined_weights) 
tri3_plot4 <- lmer(lateralorbitofrontal_vol_cortical ~ `Prenatal infection trimester 3`*timepoint + AGE_M_v2 + GENDER + ETHNMv2 + EDUCM_3groups + GSI + SMOKE_ALL + (1|IDC), data = single_df2, weights = combined_weights) 
tri3_plot5 <- lmer(middletemporal_vol_cortical ~ `Prenatal infection trimester 3`*timepoint + AGE_M_v2 + GENDER + ETHNMv2 + EDUCM_3groups + GSI + SMOKE_ALL + (1|IDC), data = single_df2, weights = combined_weights) 
tri3_plot6 <- lmer(paracentral_vol_cortical ~ `Prenatal infection trimester 3`*timepoint + AGE_M_v2 + GENDER + ETHNMv2 + EDUCM_3groups + GSI + SMOKE_ALL + (1|IDC), data = single_df2, weights = combined_weights) 
tri3_plot7 <- lmer(parsorbitalis_vol_cortical ~ `Prenatal infection trimester 3`*timepoint + AGE_M_v2 + GENDER + ETHNMv2 + EDUCM_3groups + GSI + SMOKE_ALL + (1|IDC), data = single_df2, weights = combined_weights) 
tri3_plot8 <- lmer(rostralanteriorcingulate_vol_cortical.x ~ `Prenatal infection trimester 3`*timepoint + AGE_M_v2 + GENDER + ETHNMv2 + EDUCM_3groups + GSI + SMOKE_ALL + (1|IDC), data = single_df2, weights = combined_weights) 
tri3_plot9 <- lmer(rostralmiddlefrontal_vol_cortical.x ~ `Prenatal infection trimester 3`*timepoint + AGE_M_v2 + GENDER + ETHNMv2 + EDUCM_3groups + GSI + SMOKE_ALL + (1|IDC), data = single_df2, weights = combined_weights) 
tri3_plot10 <- lmer(superiorfrontal_vol_cortical.x ~ `Prenatal infection trimester 3`*timepoint + AGE_M_v2 + GENDER + ETHNMv2 + EDUCM_3groups + GSI + SMOKE_ALL + (1|IDC), data = single_df2, weights = combined_weights) 
tri3_plot11 <- lmer(superiorparietal_vol_cortical.x ~ `Prenatal infection trimester 3`*timepoint + AGE_M_v2 + GENDER + ETHNMv2 + EDUCM_3groups + GSI + SMOKE_ALL + (1|IDC), data = single_df2, weights = combined_weights) 
tri3_plot12 <- lmer(temporalpole_vol_cortical.x ~ `Prenatal infection trimester 3`*timepoint + AGE_M_v2 + GENDER + ETHNMv2 + EDUCM_3groups + GSI + SMOKE_ALL + (1|IDC), data = single_df2, weights = combined_weights) 
tri3_plot13 <- lmer(Caudate_vol_subcortical ~ `Prenatal infection trimester 3`*timepoint + AGE_M_v2 + GENDER + ETHNMv2 + EDUCM_3groups + GSI + SMOKE_ALL + (1|IDC), data = single_df2, weights = combined_weights) 

p13 <- interact_plot(tri3_plot1, pred = timepoint, modx = `Prenatal infection trimester 3`,modx.values = 'plus-minus', x.label = 'Timepoint', y.label = 'Caudal middle frontal', interval = T, plot.points = F,point.size = 0.3, point.alpha = 0.1, line.thickness = 1, jitter = 0.2, int.type = 'confidence', colors = 'Qual1', vary.lty = F) + theme(panel.background = element_rect(fill = 'white', colour = 'black', size = 0.5))
p14 <- interact_plot(tri3_plot2, pred = timepoint, modx = `Prenatal infection trimester 3`,modx.values = 'plus-minus', x.label = 'Timepoint', y.label = 'Fusiform', interval = T, plot.points = F,point.size = 0.3, point.alpha = 0.1, line.thickness = 1, jitter = 0.2, int.type = 'confidence', colors = 'Qual1', vary.lty = F) + theme(panel.background = element_rect(fill = 'white', colour = 'black', size = 0.5))
p15 <- interact_plot(tri3_plot3, pred = timepoint, modx = `Prenatal infection trimester 3`,modx.values = 'plus-minus', x.label = 'Timepoint', y.label = 'Inferior temporal', interval = T, plot.points = F,point.size = 0.3, point.alpha = 0.1, line.thickness = 1, jitter = 0.2, int.type = 'confidence', colors = 'Qual1', vary.lty = F) + theme(panel.background = element_rect(fill = 'white', colour = 'black', size = 0.5))
p16 <- interact_plot(tri3_plot4, pred = timepoint, modx = `Prenatal infection trimester 3`,modx.values = 'plus-minus', x.label = 'Timepoint', y.label = 'Lateral orbitofrontal', interval = T, plot.points = F,point.size = 0.3, point.alpha = 0.1, line.thickness = 1, jitter = 0.2, int.type = 'confidence', colors = 'Qual1', vary.lty = F) + theme(panel.background = element_rect(fill = 'white', colour = 'black', size = 0.5))
p17 <- interact_plot(tri3_plot5, pred = timepoint, modx = `Prenatal infection trimester 3`,modx.values = 'plus-minus', x.label = 'Timepoint', y.label = 'Middle temporal', interval = T, plot.points = F,point.size = 0.3, point.alpha = 0.1, line.thickness = 1, jitter = 0.2, int.type = 'confidence', colors = 'Qual1', vary.lty = F) + theme(panel.background = element_rect(fill = 'white', colour = 'black', size = 0.5))
p18 <- interact_plot(tri3_plot6, pred = timepoint, modx = `Prenatal infection trimester 3`,modx.values = 'plus-minus', x.label = 'Timepoint', y.label = 'Paracentral', interval = T, plot.points = F,point.size = 0.3, point.alpha = 0.1, line.thickness = 1, jitter = 0.2, int.type = 'confidence', colors = 'Qual1', vary.lty = F) + theme(panel.background = element_rect(fill = 'white', colour = 'black', size = 0.5))
p19 <- interact_plot(tri3_plot7, pred = timepoint, modx = `Prenatal infection trimester 3`,modx.values = 'plus-minus', x.label = 'Timepoint', y.label = 'Pars orbitalis', interval = T, plot.points = F,point.size = 0.3, point.alpha = 0.1, line.thickness = 1, jitter = 0.2, int.type = 'confidence', colors = 'Qual1', vary.lty = F) + theme(panel.background = element_rect(fill = 'white', colour = 'black', size = 0.5))
p20 <- interact_plot(tri3_plot8, pred = timepoint, modx = `Prenatal infection trimester 3`,modx.values = 'plus-minus', x.label = 'Timepoint', y.label = 'Rostral anterior cingulate', interval = T, plot.points = F,point.size = 0.3, point.alpha = 0.1, line.thickness = 1, jitter = 0.2, int.type = 'confidence', colors = 'Qual1', vary.lty = F) + theme(panel.background = element_rect(fill = 'white', colour = 'black', size = 0.5))
p21 <- interact_plot(tri3_plot9, pred = timepoint, modx = `Prenatal infection trimester 3`,modx.values = 'plus-minus', x.label = 'Timepoint', y.label = 'Rostral middle frontal', interval = T, plot.points = F,point.size = 0.3, point.alpha = 0.1, line.thickness = 1, jitter = 0.2, int.type = 'confidence', colors = 'Qual1', vary.lty = F) + theme(panel.background = element_rect(fill = 'white', colour = 'black', size = 0.5))
p22 <- interact_plot(tri3_plot10, pred = timepoint, modx = `Prenatal infection trimester 3`,modx.values = 'plus-minus', x.label = 'Timepoint', y.label = 'Superior frontal', interval = T, plot.points = F,point.size = 0.3, point.alpha = 0.1, line.thickness = 1, jitter = 0.2, int.type = 'confidence', colors = 'Qual1', vary.lty = F) + theme(panel.background = element_rect(fill = 'white', colour = 'black', size = 0.5))
p23 <- interact_plot(tri3_plot11, pred = timepoint, modx = `Prenatal infection trimester 3`,modx.values = 'plus-minus', x.label = 'Timepoint', y.label = 'Superior parietal', interval = T, plot.points = F,point.size = 0.3, point.alpha = 0.1, line.thickness = 1, jitter = 0.2, int.type = 'confidence', colors = 'Qual1', vary.lty = F) + theme(panel.background = element_rect(fill = 'white', colour = 'black', size = 0.5))
p24 <- interact_plot(tri3_plot12, pred = timepoint, modx = `Prenatal infection trimester 3`,modx.values = 'plus-minus', x.label = 'Timepoint', y.label = 'Temporal pole', interval = T, plot.points = F,point.size = 0.3, point.alpha = 0.1, line.thickness = 1, jitter = 0.2, int.type = 'confidence', colors = 'Qual1', vary.lty = F) + theme(panel.background = element_rect(fill = 'white', colour = 'black', size = 0.5))
p25 <- interact_plot(tri3_plot13, pred = timepoint, modx = `Prenatal infection trimester 3`,modx.values = 'plus-minus', x.label = 'Timepoint', y.label = 'Caudate', interval = T, plot.points = F,point.size = 0.3, point.alpha = 0.1, line.thickness = 1, jitter = 0.2, int.type = 'confidence', colors = 'Qual1', vary.lty = F) + theme(panel.background = element_rect(fill = 'white', colour = 'black', size = 0.5))

ggarrange(p13, p14, p15, p16, p17, p18, p19, p20, p21, p22, p23, p24, p25,common.legend = T)

#\ 

###---now we make plots for FDR results---###

### Make a brain anatomy plot only for FDR sign results and align that with longitudinal without jitter plot
tri3_fdr_results_df<- data.frame(
  region = c('middle temporal','pars orbitalis','rostral anterior cingulate', 'superior frontal', 'temporal pole'),
  beta = c(0.064, 0.073, 0.073, 0.059, 0.076))

fdr_brain <- tri3_fdr_results_df %>%
  ggplot() +
  geom_brain(atlas = dk,
             aes(fill = beta)) +
  theme_bw() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  theme(strip.text = element_text(colour = 'black', face = 'bold')) +
  theme(plot.title = element_text(size = 12)) + 
  theme(legend.position = 'bottom', legend.title = element_text(face = 'bold', family = 'serif')) +
  scale_fill_gradient(low = '#330099', high = 'dark red') + 
  guides(fill = guide_colorbar(barheight = 1.5)) + 
  labs(fill = 'Beta coefficient')

#\ 

### plots longitudinal raw brain volumes: trimester3, interaction, no icv, no jitter 
single_df2$'Prenatal infection trimester 3' <- single_df2$sumscore_inf_tri3

fig_cort0 <- lmer(middletemporal_vol_cortical ~ `Prenatal infection trimester 3`*timepoint + AGE_M_v2 + GENDER + ETHNMv2 + EDUCM_3groups + GSI + SMOKE_ALL + (1|IDC), data = single_df2, weights = combined_weights) 
fig_cort1 <- lmer(parsorbitalis_vol_cortical ~ `Prenatal infection trimester 3`*timepoint + AGE_M_v2 + GENDER + ETHNMv2 + EDUCM_3groups + GSI + SMOKE_ALL + (1|IDC), data = single_df2, weights = combined_weights) 
fig_cort2 <- lmer(rostralanteriorcingulate_vol_cortical.x ~ `Prenatal infection trimester 3`*timepoint + AGE_M_v2 + GENDER + ETHNMv2 + EDUCM_3groups + GSI + SMOKE_ALL + (1|IDC), data = single_df2, weights = combined_weights) 
fig_cort3<- lmer(superiorfrontal_vol_cortical.x ~ `Prenatal infection trimester 3`*timepoint + AGE_M_v2 + GENDER + ETHNMv2 + EDUCM_3groups + GSI + SMOKE_ALL + (1|IDC), data = single_df2, weights = combined_weights) 
fig_cort4 <- lmer(temporalpole_vol_cortical.x ~ `Prenatal infection trimester 3`*timepoint + AGE_M_v2 + GENDER + ETHNMv2 + EDUCM_3groups + GSI + SMOKE_ALL + (1|IDC), data = single_df2, weights = combined_weights) 

plot0 <- interact_plot(fig_cort0, 
                       pred = timepoint, 
                       modx = `Prenatal infection trimester 3`,
                       modx.values = 'plus-minus',
                       x.label = 'Timepoint', 
                       y.label = 'Middle temporal',
                       interval = T, 
                       plot.points = F, 
                       point.size = 0.3, 
                       point.alpha = 0.1, 
                       line.thickness = 1, 
                       jitter = 0.2,
                       int.type = 'confidence',
                       colors = 'Qual1',
                       vary.lty = F) + 
  theme(panel.background = element_rect(fill = 'white', 
                                        colour = 'black', 
                                        size = 0.5))
plot1 <- interact_plot(fig_cort1, 
                       pred = timepoint, 
                       modx = `Prenatal infection trimester 3`,
                       modx.values = 'plus-minus',
                       x.label = 'Timepoint', 
                       y.label = 'Pars orbitalis',
                       interval = T, 
                       plot.points = F, 
                       point.size = 0.3, 
                       point.alpha = 0.1, 
                       line.thickness = 1, 
                       jitter = 0.2,
                       int.type = 'confidence',
                       colors = 'Qual1',
                       vary.lty = F) + 
  theme(panel.background = element_rect(fill = 'white', 
                                        colour = 'black', 
                                        size = 0.5))

plot2 <- interact_plot(fig_cort2, 
                       pred = timepoint, 
                       modx = `Prenatal infection trimester 3`,
                       modx.values = 'plus-minus',
                       x.label = 'Timepoint', 
                       y.label = 'Rostral anterior cingulate',
                       interval = T, 
                       plot.points = F, 
                       point.size = 0.3, 
                       point.alpha = 0.1, 
                       line.thickness = 1, 
                       jitter = 0.2,
                       int.type = 'confidence',
                       colors = 'Qual1',
                       vary.lty = F) + 
  theme(panel.background = element_rect(fill = 'white', 
                                        colour = 'black', 
                                        size = 0.5))

plot3 <- interact_plot(fig_cort3, 
                       pred = timepoint, 
                       modx = `Prenatal infection trimester 3`,
                       modx.values = 'plus-minus',
                       x.label = 'Timepoint', 
                       y.label = 'Superior frontal',
                       interval = T, 
                       plot.points = F, 
                       point.size = 0.3, 
                       point.alpha = 0.1, 
                       line.thickness = 1, 
                       jitter = 0.2,
                       int.type = 'confidence',
                       colors = 'Qual1',
                       vary.lty = F) + 
  theme(panel.background = element_rect(fill = 'white', 
                                        colour = 'black', 
                                        size = 0.5))

plot4 <- interact_plot(fig_cort4, 
                       pred = timepoint, 
                       modx = `Prenatal infection trimester 3`,
                       modx.values = 'plus-minus',
                       x.label = 'Timepoint', 
                       y.label = 'Temporal pole',
                       interval = T, 
                       plot.points = F, 
                       point.size = 0.3, 
                       point.alpha = 0.1, 
                       line.thickness = 1, 
                       jitter = 0.2,
                       int.type = 'confidence',
                       colors = 'Qual1',
                       vary.lty = F) + 
  theme(panel.background = element_rect(fill = 'white', 
                                        colour = 'black', 
                                        size = 0.5))

#\

### forest plot of fdr significant trimester 3 results across all trimesters to show if confidence intervals are overlapping regardless of pvalue 
df_forestplot_totalpreg <- data.frame(
  outcome = c("Middle temporal", "Pars orbitalis","Rostral anterior cingulate", "Superior frontal", "Temporal pole"),
  index = 1:5,
  beta = c(0.023, 0.032,0.037, 0.021,0.060),
  lower = c(-0.014, -0.015,-0.010,-0.015,0.010),
  upper= c(0.060, 0.079,0.084,0.057,0.110)
)

df_forestplot_totalpreg$Timepoint <- 'Total pregnancy'

df_forestplot_tri1 <- data.frame(
  outcome = c("Middle temporal", "Pars orbitalis","Rostral anterior cingulate", "Superior frontal", "Temporal pole"),
  index = 1:5,
  beta = c(0.017, 0.023,0.039,0.012, 0.058),
  lower = c(-0.023, -0.025,-0.010,-0.025, 0.007),
  upper= c(0.057, 0.072,0.088,0.05, 0.109)
)

df_forestplot_tri1$Timepoint <- 'Trimester 1'

df_forestplot_tri2 <- data.frame(
  outcome = c("Middle temporal", "Pars orbitalis","Rostral anterior cingulate", "Superior frontal", "Temporal pole"),
  index = 1:5,
  beta = c(-0.002, -0.014, -0.016, -0.012, 0.006),
  lower = c(-0.039, -0.060, -0.062, -0.047,-0.043),
  upper= c(0.034, 0.032, 0.030,0.022,0.054)
)

df_forestplot_tri2$Timepoint <- 'Trimester 2'

df_forestplot_tri3 <- data.frame(
  outcome = c("Middle temporal", "Pars orbitalis","Rostral anterior cingulate", "Superior frontal", "Temporal pole"),
  index = 1:5,
  beta = c(0.064,0.073,0.073,0.059,0.076),
  lower = c(0.024,0.025,0.025,0.022,0.025),
  upper= c(0.104,0.120,0.121,0.096,0.127)
)

df_forestplot_tri3$Timepoint <- 'Trimester 3'

df_forestplot_merged <- rbind(df_forestplot_totalpreg, df_forestplot_tri1, df_forestplot_tri2, df_forestplot_tri3) #combine all dataframes for forest plot figure

model_order <- c('Total pregnancy', 
                 'Trimester 1', 
                 'Trimester 2', 
                 'Trimester 3') #Define the order of levels for the Model column

df_forestplot_merged$Timepoint <- factor(df_forestplot_merged$Timepoint, levels = model_order) #make a factor so it can be grouped by time point in figure later 

dotCOLS <- c('black', 'black','black','black','black')
barCOLS <- c('gold','#669900', 'blue4','#993300')

forestplot <- ggplot(df_forestplot_merged, 
                     aes(x=outcome, 
                         y=beta,
                         ymin=lower,
                         ymax=upper,
                         col=Timepoint,
                         fill=Timepoint)) +
  geom_linerange(size=4,
                 position=position_dodge(width= -0.5)) + 
  geom_point(size=3, 
             shape=22, 
             colour = 'black', 
             stroke = 0.5, 
             position=position_dodge(width= -0.5)) +
  geom_hline(yintercept = 0, 
             lty=2) + 
  coord_flip() + 
  theme_classic() + 
  scale_x_discrete(name = 'Brain morphology outcome') +
  scale_y_continuous(name='Beta coefficient', 
                     limits = c(-0.07, 0.15)) + 
  theme(legend.position = 'bottom',
        legend.title = element_text(face = 'bold', family = 'serif'),
        legend.text = element_text(size=11, family = 'serif'), 
        axis.text = element_text(size=11, family = 'serif'),
        axis.title.x = element_text(size=12,face = 'bold', family = 'serif'),
        axis.title.y = element_text(size=12,face = 'bold', family = 'serif'),
        axis.text.x = element_text(size=11, hjust =0.5), 
        axis.text.y = element_text(size=11),
        panel.background = element_rect(colour = 'black')) + 
  scale_fill_manual(values = dotCOLS) + 
  scale_color_manual(values=barCOLS) 

#\

### combine plots for main manuscript

#combine forest plot, brain anatomy plot and longitudinal plot 
combined_long_plot <- ggarrange(plot0, plot1, plot2, plot3, plot4, common.legend = T)
combined_brain <- ggarrange(fdr_brain, forestplot, labels = 'AUTO', ncol = 2)

com1 <- ggarrange(forestplot,combined_long_plot, ncol = 2)
com2 <- ggarrange(fdr_brain,com1, ncol = 1)
com2

#\

### longitudinal plots with jitter for supplements 
plot0b <- interact_plot(fig_cort0, 
                        pred = timepoint, 
                        modx = `Prenatal infection trimester 3`,
                        modx.values = 'plus-minus',
                        colors = 'seagreen', 
                        x.label = 'Timepoint', 
                        y.label = 'Middle temporal',
                        interval = T, 
                        plot.points = T, 
                        point.size = 0.3, 
                        point.alpha = 0.1, 
                        line.thickness = 1, 
                        jitter = 0.2,
                        int.type = 'confidence',
                        vary.lty = F) + 
  theme(panel.background = element_rect(fill = 'white', 
                                        colour = 'black', 
                                        size = 0.5))

plot1b <- interact_plot(fig_cort1, 
                        pred = timepoint, 
                        modx = `Prenatal infection trimester 3`,
                        modx.values = 'plus-minus',
                        colors = 'seagreen', 
                        x.label = 'Timepoint', 
                        y.label = 'Pars orbitalis',
                        interval = T, 
                        plot.points = T, 
                        point.size = 0.3, 
                        point.alpha = 0.1, 
                        line.thickness = 1, 
                        jitter = 0.2,
                        int.type = 'confidence',
                        vary.lty = F) + 
  theme(panel.background = element_rect(fill = 'white', 
                                        colour = 'black', 
                                        size = 0.5))

plot2b <- interact_plot(fig_cort2, 
                        pred = timepoint, 
                        modx = `Prenatal infection trimester 3`,
                        modx.values = 'plus-minus',
                        colors = 'seagreen', 
                        x.label = 'Timepoint', 
                        y.label = 'Rostral anterior cingulate',
                        interval = T, 
                        plot.points = T, 
                        point.size = 0.3, 
                        point.alpha = 0.1, 
                        line.thickness = 1, 
                        jitter = 0.2,
                        int.type = 'confidence',
                        vary.lty = F) + 
  theme(panel.background = element_rect(fill = 'white', 
                                        colour = 'black', 
                                        size = 0.5))

plot3b <- interact_plot(fig_cort3, 
                        pred = timepoint, 
                        modx = `Prenatal infection trimester 3`,
                        modx.values = 'plus-minus',
                        colors = 'seagreen', 
                        x.label = 'Timepoint', 
                        y.label = 'Superior frontal',
                        interval = T, 
                        plot.points = T, 
                        point.size = 0.3, 
                        point.alpha = 0.1, 
                        line.thickness = 1, 
                        jitter = 0.2,
                        int.type = 'confidence',
                        vary.lty = F) + 
  theme(panel.background = element_rect(fill = 'white', 
                                        colour = 'black', 
                                        size = 0.5))

plot4b <- interact_plot(fig_cort4, 
                        pred = timepoint, 
                        modx = `Prenatal infection trimester 3`,
                        modx.values = 'plus-minus',
                        colors = 'seagreen', 
                        x.label = 'Timepoint', 
                        y.label = 'Temporal pole',
                        interval = T, 
                        plot.points = T, 
                        point.size = 0.3, 
                        point.alpha = 0.1, 
                        line.thickness = 1, 
                        jitter = 0.2,
                        int.type = 'confidence',
                        vary.lty = F) + 
  theme(panel.background = element_rect(fill = 'white', 
                                        colour = 'black', 
                                        size = 0.5))

plot_grid(plot0b, plot1b, plot2b, plot3b, plot4b, ncol = 2, labels = 'AUTO')

#\

###---now we make plots for all results---###

### ggseg plot with all results for supplements; regardless of pvalue
# create dataframes with all results for each trimester 
all_results_df_total_preg_cortex<- data.frame(
  region = c('bankssts', 'caudal anterior cingulate', 'caudal middle frontal', 'cuneus', 'entorhinal', 'fusiform', 'inferior parietal', 'inferior temporal', 'isthmus cingulate', 'lateral occipital', 'lateral orbitofrontal', 'lingual', 'medial orbitofrontal', 'middle temporal', 'parahippocampal', 'paracentral', 'pars opercularis', 'pars orbitalis', 'pars triangularis', 'pericalcarine', 'postcentral', 'posterior cingulate', 'precentral', 'precuneus', 'rostral anterior cingulate', 'rostral middle frontal', 'superior frontal', 'superior parietal', 'superior temporal', 'supramarginal', 'frontal pole', 'temporal pole', 'transverse temporal', 'insula'),
  beta = c(-0.024, -0.003, 0.023, -0.008, 0.002, 0.017, -0.004, 0.023, -0.009, 0.008, 0.034, -0.013, 0.005, 0.032, -0.006, 0.003, -0.007, 0.032, -0.029, -0.008, 0.012, 0.003, 0.000, 0.007, 0.037, 0.011, 0.021, 0.023, 0.003, 0.004, 0.024, 0.060, -0.013, 0.024))
all_results_df_total_preg_cortex$timepoint <- 'Total pregnancy'

all_results_df_tri1_cortex<- data.frame(
  region = c('bankssts', 'caudal anterior cingulate', 'caudal middle frontal', 'cuneus', 'entorhinal', 'fusiform', 'inferior parietal', 'inferior temporal', 'isthmus cingulate', 'lateral occipital', 'lateral orbitofrontal', 'lingual', 'medial orbitofrontal', 'middle temporal', 'parahippocampal', 'paracentral', 'pars opercularis', 'pars orbitalis', 'pars triangularis', 'pericalcarine', 'postcentral', 'posterior cingulate', 'precentral', 'precuneus', 'rostral anterior cingulate', 'rostral middle frontal', 'superior frontal', 'superior parietal', 'superior temporal', 'supramarginal', 'frontal pole', 'temporal pole', 'transverse temporal', 'insula'),
  beta = c(-0.014, 0.025, 0.031, -0.010, 0.014, 0.018, 0.014, 0.013, -0.007, 0.018, 0.035, 0.003, 0.026, 0.017, -0.007, -0.005, 0.008, 0.023, -0.023, -0.011, 0.012, 0.014, 0.018, 0.018, 0.039, 0.018, 0.012, 0.041, 0.008, 0.020, 0.013, 0.058, -0.01, 0.019))
all_results_df_tri1_cortex$timepoint <-'Trimester 1'

all_results_df_tri2_cortex<- data.frame(
  region = c('bankssts', 'caudal anterior cingulate', 'caudal middle frontal', 'cuneus', 'entorhinal', 'fusiform', 'inferior parietal', 'inferior temporal', 'isthmus cingulate', 'lateral occipital', 'lateral orbitofrontal', 'lingual', 'medial orbitofrontal', 'middle temporal', 'parahippocampal', 'paracentral', 'pars opercularis', 'pars orbitalis', 'pars triangularis', 'pericalcarine', 'postcentral', 'posterior cingulate', 'precentral', 'precuneus', 'rostral anterior cingulate', 'rostral middle frontal', 'superior frontal', 'superior parietal', 'superior temporal', 'supramarginal', 'frontal pole', 'temporal pole', 'transverse temporal', 'insula'),
  beta = c(-0.032, -0.050, -0.018, -0.018, -0.012, -0.012, -0.038, 0.006, -0.013, -0.009, -0.01, -0.013, -0.034, -0.002, -0.03, -0.016, -0.037, -0.014, -0.048, -0.017, -0.005, -0.018, -0.028, -0.022, -0.016, -0.022, -0.012, -0.030, -0.023, -0.021, 0.005, 0.006, -0.019, 0.012))
all_results_df_tri2_cortex$timepoint <-'Trimester 2'

all_results_df_tri3_cortex<- data.frame(
  region = c('bankssts', 'caudal anterior cingulate', 'caudal middle frontal', 'cuneus', 'entorhinal', 'fusiform', 'inferior parietal', 'inferior temporal', 'isthmus cingulate', 'lateral occipital', 'lateral orbitofrontal', 'lingual', 'medial orbitofrontal', 'middle temporal', 'parahippocampal', 'paracentral', 'pars opercularis', 'pars orbitalis', 'pars triangularis', 'pericalcarine', 'postcentral', 'posterior cingulate', 'precentral', 'precuneus', 'rostral anterior cingulate', 'rostral middle frontal', 'superior frontal', 'superior parietal', 'superior temporal', 'supramarginal', 'frontal pole', 'temporal pole', 'transverse temporal', 'insula'),
  beta = c(0, 0.029, 0.045, 0.014, 0.012, 0.04, 0.028, 0.04, 0.005, 0.014, 0.064, -0.013, 0.036, 0.064, 0.029, 0.036, 0.023, 0.073, 0.017, 0.015, 0.026, 0.017, 0.02, 0.027, 0.073, 0.043, 0.059, 0.051, 0.031, 0.018, 0.042, 0.076, 0.004, 0.027))
all_results_df_tri3_cortex$timepoint <-'Trimester 3'

total_df_results_cortex <- rbind(all_results_df_total_preg_cortex, all_results_df_tri1_cortex, all_results_df_tri2_cortex, all_results_df_tri3_cortex)

all_results_df_total_preg_subcortex<- data.frame(
  region = c('cerebellum cortex', 'amygdala', 'hippocampus', 'caudate', 'putamen', 'thalamus proper', 'pallidum'),
  beta = c(0.038, 0.016, -0.001, 0.027, 0.013, -0.035, 0.025))
all_results_df_total_preg_subcortex$timepoint <-'Total pregnancy'

all_results_df_tri1_subcortex<- data.frame(
  region = c('cerebellum cortex', 'amygdala', 'hippocampus', 'caudate', 'putamen', 'thalamus proper', 'pallidum'),
  beta = c(0.023, 0.018, -0.008, 0.026, 0.018, -0.013, 0.019))
all_results_df_tri1_subcortex$timepoint <-'Trimester 1'

all_results_df_tri2_subcortex<- data.frame(
  region = c('cerebellum cortex', 'amygdala', 'hippocampus', 'caudate', 'putamen', 'thalamus proper', 'pallidum'),
  beta = c(0.037, -0.006, -0.006, -0.010, -0.006, -0.030, 0.009))
all_results_df_tri2_subcortex$timepoint <-'Trimester 2'

all_results_df_tri3_subcortex<- data.frame(
  region = c('cerebellum cortex', 'amygdala', 'hippocampus', 'caudate', 'putamen', 'thalamus proper', 'pallidum'),
  beta = c(0.025, 0.035, 0.017, 0.052, 0.022, -0.027, 0.034))
all_results_df_tri3_subcortex$timepoint <-'Trimester 3'

total_df_results_subcortex <- rbind(all_results_df_total_preg_subcortex, all_results_df_tri1_subcortex, all_results_df_tri2_subcortex, all_results_df_tri3_subcortex)

# create cortical and subcortical figure 
a2 <- total_df_results_cortex %>%
  group_by(timepoint) %>%
  ggplot() +
  geom_brain(atlas = dk,
             position = position_brain(hemi ~ side),
             aes(fill = beta)) +
  facet_wrap(~timepoint) +
  theme_bw() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  theme(strip.background = element_rect(fill = 'grey83')) +
  theme(strip.text = element_text(colour = 'black', face = 'bold')) +
  theme(plot.title = element_text(size = 12, family = 'sans')) + 
  theme(legend.position = 'bottom',
        legend.text = element_text(size = 8, angle =45))+ 
  guides(fill = guide_colorbar(barheight = 1.5)) + 
  labs(fill = 'Beta coefficient')

b2 <- total_df_results_subcortex %>%
  group_by(timepoint) %>%
  ggplot() +
  geom_brain(atlas = aseg,
             aes(fill = beta)) +
  facet_wrap(~timepoint) +
  theme_bw() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  theme(strip.background = element_rect(fill = 'grey83')) +
  theme(strip.text = element_text(colour = 'black', face = 'bold')) +
  theme(plot.title = element_text(size = 12, family = 'sans')) + 
  theme(legend.position = 'bottom',
        legend.text = element_text(size = 8, angle =45))+ 
  guides(fill = guide_colorbar(barheight = 1.5)) + 
  labs(fill = 'Beta coefficient')

#Combing plots for cortex and subcortex
plot_grid(a2, b2, labels = c('A', 'B'), ncol = 2)

#\ END OF ALL SCRIPTS FOR THIS PROJECT.
