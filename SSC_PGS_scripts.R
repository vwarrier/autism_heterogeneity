library(data.table)

Cases = fread("~/SFARI/cases.txt", header= T)
setnames(Cases, "id", "individual")

prs_scz = fread("~/ALSPAC/PRSice2results/SFARI_scz_eas_eur_lam2019.all.score", header = TRUE)
setnames(prs_scz, 3, "scz_0001")
setnames(prs_scz, 4, "scz_001")
setnames(prs_scz, 5, "scz_005")
setnames(prs_scz, 6, "scz_01")
setnames(prs_scz, 7, "scz_05")
setnames(prs_scz, 8, "scz_1")
setnames(prs_scz, 9, "scz_25")
setnames(prs_scz, 10, "scz_5")
setnames(prs_scz, 11, "scz_75")
setnames(prs_scz, 12, "scz_all")

prs_scz = prs_scz[,c(2:12)]



prs_edu = fread("~/ALSPAC/PRSice2results/SFARI_edu2018leeforPRSICE.all.score", header = TRUE)
setnames(prs_edu, 3, "edu_0001")
setnames(prs_edu, 4, "edu_001")
setnames(prs_edu, 5, "edu_005")
setnames(prs_edu, 6, "edu_01")
setnames(prs_edu, 7, "edu_05")
setnames(prs_edu, 8, "edu_1")
setnames(prs_edu, 9, "edu_25")
setnames(prs_edu, 10, "edu_5")
setnames(prs_edu, 11, "edu_75")
setnames(prs_edu, 12, "edu_all")

prs_edu = prs_edu[,c(2:12)]


prs_IQ = fread("~/ALSPAC/PRSice2results/SFARI_savageCP2018forprsice.all.score", header = TRUE)
setnames(prs_IQ, 3, "IQ_0001")
setnames(prs_IQ, 4, "IQ_001")
setnames(prs_IQ, 5, "IQ_005")
setnames(prs_IQ, 6, "IQ_01")
setnames(prs_IQ, 7, "IQ_05")
setnames(prs_IQ, 8, "IQ_1")
setnames(prs_IQ, 9, "IQ_25")
setnames(prs_IQ, 10, "IQ_5")
setnames(prs_IQ, 11, "IQ_75")
setnames(prs_IQ, 12, "IQ_all")

prs_merged = merge(prs_IQ, prs_edu, by = "IID")
prs_merged = merge(prs_merged, prs_scz, by = "IID")

pheno = fread("~/SFARI/ssccore.txt")

pca = fread("~/SFARI/liftOverPlink/files_imputed/SSC_cases_pca.eigenvec")
setnames(pca, 2, "IID")

merged = merge(prs_merged, pca, by = "IID")
merged = merge(Cases, merged, by = "IID")
merged = merge(merged, pheno, by = "individual")


#ADOS: soc affect
#ADOS: RRB
#ADI: communication total
#ADI: RRB
#ADI: Social
#RBS
#VABS: total score
#SRS_parent_raw
#Full_scale IQ
#Nonverbal IQ 
#Verbal IQ

ados_list = list("ados_restricted_repetitive")
other_list = list("adi_r_b_comm_verbal_total", "adi_r_comm_b_non_verbal_total", "adi_r_rrb_c_total","adi_r_soc_a_total", "rbs_r_overall_score", "srs_parent_raw_total", "vineland_ii_composite_standard_score")
IQ_list = list("ssc_diagnosis_full_scale_iq", "ssc_diagnosis_nonverbal_iq", "ssc_diagnosis_verbal_iq")


## Run regression on ADOS_social_affect

IQ_all = as.data.frame((summary(lm(scale(as.numeric(as.character(ados_social_affect))) ~ scale(IQ_all) + age_at_ados + ados_module + as.character(Sex) + as.character(dataset) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) + 
                                     scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = merged)))$coefficients[2,])

Edu_all = as.data.frame((summary(lm(scale(as.numeric(as.character(ados_social_affect))) ~ scale(edu_all) + age_at_ados + ados_module + as.character(Sex) + as.character(dataset) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) + 
                                     scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = merged)))$coefficients[2,])

SCZ_all = as.data.frame((summary(lm(scale(as.numeric(as.character(ados_social_affect))) ~ scale(scz_1) + age_at_ados + ados_module + as.character(Sex) + as.character(dataset) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) + 
                                     scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = merged)))$coefficients[2,])




for(i in ados_list){
  
  
  
  alpha = lm(paste0("scale(", i, ") ~ scale(IQ_all) + age_at_ados + ados_module + as.character(Sex) + as.character(dataset) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) + 
                                     scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17)"), data = merged)
  
  IQ = (summary(alpha))$coefficients[2,]
  
  beta = lm(paste0("scale(", i, ") ~ scale(edu_all) + age_at_ados + ados_module + as.character(Sex) + as.character(dataset) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) + 
                                     scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17)"), data = merged)
  
  
  Edu = (summary(beta))$coefficients[2,]
  
  gamma = lm(paste0("scale(", i, ") ~ scale(scz_1) + age_at_ados + ados_module + as.character(Sex) + as.character(dataset) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) + 
                                     scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17)"), data = merged)
  
  
  SCZ = (summary(gamma))$coefficients[2,]
  
  IQ_all = cbind(IQ_all, IQ)
  Edu_all = cbind(Edu_all, Edu)
  SCZ_all = cbind(SCZ_all, SCZ)
}


for(i in other_list){
  
  
  
  alpha = lm(paste0("scale(", i, ") ~ scale(IQ_all)  + as.character(Sex) + as.character(dataset) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) + 
                                     scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17)"), data = merged)
  
  IQ = (summary(alpha))$coefficients[2,]
  
  beta = lm(paste0("scale(", i, ") ~ scale(edu_all) +  as.character(Sex) + as.character(dataset) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) + 
                                     scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17)"), data = merged)
  
  
  Edu = (summary(beta))$coefficients[2,]
  
  gamma = lm(paste0("scale(", i, ") ~ scale(scz_1) +  as.character(Sex) + as.character(dataset) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) + 
                                     scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17)"), data = merged)
  
  
  SCZ = (summary(gamma))$coefficients[2,]
  
  IQ_all = cbind(IQ_all, IQ)
  Edu_all = cbind(Edu_all, Edu)
  SCZ_all = cbind(SCZ_all, SCZ)
}



for(i in IQ_list){
  
  
  
  alpha = lm(paste0("scale(", i, ") ~ scale(IQ_all)  + as.character(ssc_diagnosis_verbal_iq_type) + as.character(Sex) + as.character(dataset) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) + 
                                     scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17)"), data = merged)
  
  IQ = (summary(alpha))$coefficients[2,]
  
  beta = lm(paste0("scale(", i, ") ~ scale(edu_all) + as.character(ssc_diagnosis_verbal_iq_type) + as.character(Sex) + as.character(dataset) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) + 
                                     scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17)"), data = merged)
  
  
  Edu = (summary(beta))$coefficients[2,]
  
  gamma = lm(paste0("scale(", i, ") ~ scale(scz_1) + as.character(ssc_diagnosis_verbal_iq_type) + as.character(Sex) + as.character(dataset) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) + 
                                     scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17)"), data = merged)
  
  
  SCZ = (summary(gamma))$coefficients[2,]
  
  IQ_all = cbind(IQ_all, IQ)
  Edu_all = cbind(Edu_all, Edu)
  SCZ_all = cbind(SCZ_all, SCZ)
}




#############Run the analysis only in males###########################
male = subset(merged, sex == "male")


IQ_male = as.data.frame((summary(lm(scale(as.numeric(as.character(ados_social_affect))) ~ scale(IQ_all) + age_at_ados + ados_module + as.character(dataset) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) + 
                                     scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = male)))$coefficients[2,])

Edu_male = as.data.frame((summary(lm(scale(as.numeric(as.character(ados_social_affect))) ~ scale(edu_all) + age_at_ados + ados_module  + as.character(dataset) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) + 
                                      scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = male)))$coefficients[2,])

SCZ_male = as.data.frame((summary(lm(scale(as.numeric(as.character(ados_social_affect))) ~ scale(scz_1) + age_at_ados + ados_module  + as.character(dataset) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) + 
                                      scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = male)))$coefficients[2,])




for(i in ados_list){
  
  
  
  alpha = lm(paste0("scale(", i, ") ~ scale(IQ_all) + age_at_ados + ados_module  + as.character(dataset) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) + 
                    scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17)"), data = male)
  
  IQ = (summary(alpha))$coefficients[2,]
  
  beta = lm(paste0("scale(", i, ") ~ scale(edu_all) + age_at_ados + ados_module + as.character(dataset) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) + 
                   scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17)"), data = male)
  
  
  Edu = (summary(beta))$coefficients[2,]
  
  gamma = lm(paste0("scale(", i, ") ~ scale(scz_1) + age_at_ados + ados_module  + as.character(dataset) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) + 
                    scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17)"), data = male)
  
  
  SCZ = (summary(gamma))$coefficients[2,]
  
  IQ_male = cbind(IQ_male, IQ)
  Edu_male = cbind(Edu_male, Edu)
  SCZ_male = cbind(SCZ_male, SCZ)
}


for(i in other_list){
  
  
  
  alpha = lm(paste0("scale(", i, ") ~ scale(IQ_all)  + as.character(dataset) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) + 
                    scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17)"), data = male)
  
  IQ = (summary(alpha))$coefficients[2,]
  
  beta = lm(paste0("scale(", i, ") ~ scale(edu_all)  + as.character(dataset) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) + 
                   scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17)"), data = male)
  
  
  Edu = (summary(beta))$coefficients[2,]
  
  gamma = lm(paste0("scale(", i, ") ~ scale(scz_1) + as.character(dataset) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) + 
                    scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17)"), data = male)
  
  
  SCZ = (summary(gamma))$coefficients[2,]
  
  IQ_male = cbind(IQ_male, IQ)
  Edu_male = cbind(Edu_male, Edu)
  SCZ_male = cbind(SCZ_male, SCZ)
}



for(i in IQ_list){
  
  
  
  alpha = lm(paste0("scale(", i, ") ~ scale(IQ_all)  + as.character(ssc_diagnosis_verbal_iq_type)  + as.character(dataset) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) + 
                                     scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17)"), data = male)
  
  IQ = (summary(alpha))$coefficients[2,]
  
  beta = lm(paste0("scale(", i, ") ~ scale(edu_all) + as.character(ssc_diagnosis_verbal_iq_type)  + as.character(dataset) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) + 
                                     scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17)"), data = male)
  
  
  Edu = (summary(beta))$coefficients[2,]
  
  gamma = lm(paste0("scale(", i, ") ~ scale(scz_1) + as.character(ssc_diagnosis_verbal_iq_type)  + as.character(dataset) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) + 
                                     scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17)"), data = male)
  
  
  SCZ = (summary(gamma))$coefficients[2,]
  
  IQ_male = cbind(IQ_male, IQ)
  Edu_male = cbind(Edu_male, Edu)
  SCZ_male = cbind(SCZ_male, SCZ)
}





############################IQ > 70#########################


non_id = subset(merged, ssc_diagnosis_full_scale_iq > 70)

IQ_table_under70 = as.data.frame((summary(lm(scale(as.numeric(as.character(ados_social_affect))) ~ scale(IQ_all) + age_at_ados + ados_module + as.character(Sex) + as.character(dataset) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) + 
                                    scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = non_id)))$coefficients[2,])

Edu_table_under70 = as.data.frame((summary(lm(scale(as.numeric(as.character(ados_social_affect))) ~ scale(edu_all) + age_at_ados + ados_module + as.character(Sex) + as.character(dataset) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) + 
                                      scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = non_id)))$coefficients[2,])

SCZ_table_under70 = as.data.frame((summary(lm(scale(as.numeric(as.character(ados_social_affect))) ~ scale(scz_1) + age_at_ados + ados_module + as.character(Sex) + as.character(dataset) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) + 
                                      scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = non_id)))$coefficients[2,])




for(i in ados_list){
  
  
  
  alpha = lm(paste0("scale(", i, ") ~ scale(IQ_all) + age_at_ados + ados_module + as.character(Sex) + as.character(dataset) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) + 
                    scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17)"), data = non_id)
  
  IQ = (summary(alpha))$coefficients[2,]
  
  beta = lm(paste0("scale(", i, ") ~ scale(edu_all) + age_at_ados + ados_module + as.character(Sex) + as.character(dataset) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) + 
                   scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17)"), data = non_id)
  
  
  Edu = (summary(beta))$coefficients[2,]
  
  gamma = lm(paste0("scale(", i, ") ~ scale(scz_1) + age_at_ados + ados_module + as.character(Sex) + as.character(dataset) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) + 
                    scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17)"), data = non_id)
  
  
  SCZ = (summary(gamma))$coefficients[2,]
  
  IQ_table_under70 = cbind(IQ_table_under70, IQ)
  Edu_table_under70 = cbind(Edu_table_under70, Edu)
  SCZ_table_under70 = cbind(SCZ_table_under70, SCZ)
}


for(i in other_list){
  
  
  
  alpha = lm(paste0("scale(", i, ") ~ scale(IQ_all)  + as.character(Sex) + as.character(dataset) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) + 
                    scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17)"), data = non_id)
  
  IQ = (summary(alpha))$coefficients[2,]
  
  beta = lm(paste0("scale(", i, ") ~ scale(edu_all) +  as.character(Sex) + as.character(dataset) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) + 
                   scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17)"), data = non_id)
  
  
  Edu = (summary(beta))$coefficients[2,]
  
  gamma = lm(paste0("scale(", i, ") ~ scale(scz_1) +  as.character(Sex) + as.character(dataset) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) + 
                    scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17)"), data = non_id)
  
  
  SCZ = (summary(gamma))$coefficients[2,]
  
  IQ_table_under70 = cbind(IQ_table_under70, IQ)
  Edu_table_under70 = cbind(Edu_table_under70, Edu)
  SCZ_table_under70 = cbind(SCZ_table_under70, SCZ)
}



for(i in IQ_list){
  
  
  
  alpha = lm(paste0("scale(", i, ") ~ scale(IQ_all)  + as.character(ssc_diagnosis_verbal_iq_type) + as.character(Sex) + as.character(dataset) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) + 
                                     scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17)"), data = non_id)
  
  IQ = (summary(alpha))$coefficients[2,]
  
  beta = lm(paste0("scale(", i, ") ~ scale(edu_all) + as.character(ssc_diagnosis_verbal_iq_type) + as.character(Sex) + as.character(dataset) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) + 
                                     scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17)"), data = non_id)
  
  
  Edu = (summary(beta))$coefficients[2,]
  
  gamma = lm(paste0("scale(", i, ") ~ scale(scz_1) + as.character(ssc_diagnosis_verbal_iq_type) + as.character(Sex) + as.character(dataset) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) + 
                                     scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17)"), data = non_id)
  
  
  SCZ = (summary(gamma))$coefficients[2,]
  
  IQ_table_under70 = cbind(IQ_table_under70, IQ)
  Edu_table_under70 = cbind(Edu_table_under70, Edu)
  SCZ_table_under70 = cbind(SCZ_table_under70, SCZ)
}




############Cor_mat#################################################
data_for_matrix = merged[,c("ados_social_affect", "ados_restricted_repetitive", "adi_r_b_comm_verbal_total",
                            "adi_r_comm_b_non_verbal_total", "adi_r_rrb_c_total", "adi_r_soc_a_total", "rbs_r_overall_score", 
                            "vineland_ii_composite_standard_score", "srs_parent_raw_total", "ssc_diagnosis_full_scale_iq", 
                            "ssc_diagnosis_nonverbal_iq", "ssc_diagnosis_verbal_iq")]

data_for_matrix$ados_social_affect = as.numeric(as.character(data_for_matrix$ados_social_affect))
correlationmatrix = cor(data_for_matrix, use = "complete.obs", method = "pearson")


hM <- format(round(correlationmatrix, 2))
col<- colorRampPalette(c("blue", "white", "red"))(20)
heatmap.2(correlationmatrix, trace = "none", col=col, cellnote=hM, notecol = "black")


heatmap.2(correlationmatrix, trace = "none", col=col)

res2<-rcorr(as.matrix(data_for_matrix))
res2_results = flattenCorrMatrix(res2$r, res2$P)




