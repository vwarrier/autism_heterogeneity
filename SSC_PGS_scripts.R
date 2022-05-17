###Step 1: Create the dataset###

library(data.table)

Cases = fread("~/SFARI/cases.txt", header= T)
setnames(Cases, "id", "individual")

prs_scz = fread("./PGS/SSC_scz_finalscore.profile", header = TRUE)
setnames(prs_scz, "SCORE", "scz_prs")
prs_scz = prs_scz[,c("IID", "scz_prs")]


prs_edu = fread("./PGS/SSC_edu_finalscore.profile", header = TRUE)
setnames(prs_edu, "SCORE", "edu_prs")
prs_edu = prs_edu[,c("IID", "edu_prs")]

prs_IQ = fread("./PGS/SSC_IQ_finalscore.profile", header = TRUE)
setnames(prs_IQ, "SCORE", "IQ_prs")
prs_IQ = prs_IQ[,c("IID", "IQ_prs")]

prs_autism = fread("./PGS/SSC_autism_finalscore.profile", header = TRUE)
setnames(prs_autism, "SCORE", "autism_prs")
prs_autism = prs_autism[,c("IID", "autism_prs")]


prs_merged = merge(prs_IQ, prs_edu, by = "IID")
prs_merged = merge(prs_merged, prs_scz, by = "IID")
prs_merged = merge(prs_merged, prs_autism, by = "IID")

pheno = fread("~/SFARI/ssccore.txt")

pca = fread("~/SFARI/liftOverPlink/files_imputed/SSC_cases_pca.eigenvec")
setnames(pca, 2, "IID")

merged = merge(prs_merged, pca, by = "IID")
merged = merge(Cases, merged, by = "IID")
merged = merge(merged, pheno, by = "individual")

merged$ados_social_affect = as.numeric(as.character(merged$ados_social_affect))

merged$IQ_prs = scale(merged$IQ_prs)
merged$autism_prs = scale(merged$autism_prs)
merged$edu_prs = scale(merged$edu_prs)
merged$scz_prs = scale(merged$scz_prs)

dcdq = fread("~/SFARI/Phenotypes/dcdq.csv")
dcdq = dcdq[,c("individual", "total")]
setnames(dcdq, "total", "dcdq_total")

merged = merge(merged, dcdq, by = "individual", all.x = TRUE)

scq = fread("~/SFARI/Phenotypes/scq_life.csv")
scq = scq[,c("individual", "summary_score")]
setnames(scq, "summary_score", "scq_total")

merged = merge(merged, scq, by = "individual", all.x = TRUE)


###Read de novo###
data_in_sample = fread("./Sebat_data/all_samples_sebat.txt")

#subset to autistic probands
ssc = subset(data_in_sample, Cohort == "SSC")
ssc = subset(ssc, Phenotype == "1")

SNV = fread("./Sebat_data/denovo_SNVindel.txt")
SV = fread("./Sebat_data/denovo_SV.txt")

##loeuf_upper_decile_score = 0.37
SNV_upper_decile = subset(SNV, gnomad_loeuf < 0.37)
SV_upper_decile = subset(SV, gnomad_loeuf_minimum < 0.37)

##Count_data
count1 = as.data.table(table(SNV$fid))
count2 = as.data.table(table(SV$fid))
count_noconstraint = merge(count1, count2, by = "V1", all = TRUE)
count_noconstraint$total = rowSums(count_noconstraint[,c("N.x", "N.y")], na.rm=TRUE)
count_noconstraint = count_noconstraint[,c("V1", "total")]
setnames(count_noconstraint, old = c("V1", "total"), new = c("FID", "anydenovo_sum"))
count_noconstraint$FID = as.numeric(as.character(count_noconstraint$FID))

merged <- merge(merged, count_noconstraint, by = "FID", all.x = TRUE)
merged$anydenovo_sum[is.na(merged$anydenovo_sum)] <- 0


count1 = as.data.table(table(SNV_upper_decile$fid))
count2 = as.data.table(table(SV_upper_decile$fid))
count_constraint = merge(count1, count2, by = "V1", all = TRUE)
count_constraint$total = rowSums(count_constraint[,c("N.x", "N.y")], na.rm=TRUE)
count_constraint = count_constraint[,c("V1", "total")]
setnames(count_constraint, old = c("V1", "total"), new = c("FID", "anydenovo_constraint_sum"))
count_constraint$FID = as.numeric(as.character(count_constraint$FID))

merged <- merge(merged, count_constraint, by = "FID", all.x = TRUE)
merged$anydenovo_constraint_sum[is.na(merged$anydenovo_constraint_sum)] <- 0


merged$anydenovo = ifelse(merged$FID %in% SNV$fid | merged$FID %in% SV$fid, 1, 0)
merged$anydenovo_constraint = ifelse(merged$FID %in% SNV_upper_decile$fid | merged$FID %in% SV_upper_decile$fid, 1, 0)


##obtain the pTDT data###
autism = fread("./PGS/SSC_autism_midparent_ptdt.txt")
setnames(autism, old = c("midparent", "diff"), new = c("autism_midparent", "autism_deviation"))

scz = fread("./PGS/SSC_scz_midparent_ptdt.txt")
setnames(scz, old = c("midparent", "diff"), new = c("scz_midparent", "scz_deviation"))

IQ = fread("./PGS/SSC_IQ_midparent_ptdt.txt")
setnames(IQ, old = c("midparent", "diff"), new = c("IQ_midparent", "IQ_deviation"))

edu = fread("./PGS/SSC_edu_midparent_ptdt.txt")
setnames(edu, old = c("midparent", "diff"), new = c("edu_midparent", "edu_deviation"))

ptdt_merged = merge(IQ, edu, by = "IID")
ptdt_merged = merge(ptdt_merged, scz, by = "IID")
ptdt_merged = merge(ptdt_merged, autism, by = "IID")

merged = merge(merged, ptdt_merged, by = "IID")

save(merged, file = "SSC_dataforanalysis.Rdata")



######

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

autism_results = NULL
scz_results = NULL
IQ_results = NULL
edu_results = NULL
de_novo_results = NULL

ados_list = list("ados_restricted_repetitive", "ados_social_affect")
other_list = list("adi_r_b_comm_verbal_total", "adi_r_rrb_c_total","adi_r_soc_a_total", "rbs_r_overall_score", "srs_parent_raw_total", "vineland_ii_composite_standard_score")
IQ_list = list("ssc_diagnosis_full_scale_iq", "ssc_diagnosis_nonverbal_iq", "ssc_diagnosis_verbal_iq")


## Run regression in all

for(i in ados_list){
  
  results_all = summary(lm(paste0("scale(", i, ") ~ IQ_prs + autism_prs + edu_prs + scz_prs + age_at_ados + ados_module + as.character(Sex) + as.character(dataset) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) + 
                                     scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17)"), data = merged))
  
  IQ_results = rbind(IQ_results, cbind(i, t(results_all$coefficients[2,])))
  autism_results = rbind(autism_results, cbind(i, t(results_all$coefficients[3,])))
  edu_results = rbind(edu_results, cbind(i, t(results_all$coefficients[4,])))
  scz_results = rbind(scz_results, cbind(i, t(results_all$coefficients[5,])))
  
}


for(i in other_list){
  
  results_all = summary(lm(paste0("scale(", i, ") ~ IQ_prs + autism_prs + edu_prs + scz_prs + age_at_ados + as.character(Sex) + as.character(dataset) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) + 
                                     scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17)"), data = merged))
  
  IQ_results = rbind(IQ_results, cbind(i, t(results_all$coefficients[2,])))
  autism_results = rbind(autism_results, cbind(i, t(results_all$coefficients[3,])))
  edu_results = rbind(edu_results, cbind(i, t(results_all$coefficients[4,])))
  scz_results = rbind(scz_results, cbind(i, t(results_all$coefficients[5,])))
  
}



for(i in IQ_list){
  
  results_all = summary(lm(paste0("scale(", i, ") ~ IQ_prs + autism_prs + edu_prs + scz_prs + age_at_ados + as.character(Sex) + as.character(ssc_diagnosis_verbal_iq_type) + as.character(dataset) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) + 
                                     scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17)"), data = merged))
  
  IQ_results = rbind(IQ_results, cbind(i, t(results_all$coefficients[2,])))
  autism_results = rbind(autism_results, cbind(i, t(results_all$coefficients[3,])))
  edu_results = rbind(edu_results, cbind(i, t(results_all$coefficients[4,])))
  scz_results = rbind(scz_results, cbind(i, t(results_all$coefficients[5,])))
  
  }


###With denovo###

autism_results = NULL
scz_results = NULL
IQ_results = NULL
edu_results = NULL
de_novo_results = NULL

ados_list = list("ados_restricted_repetitive", "ados_social_affect")
other_list = list("adi_r_b_comm_verbal_total", "adi_r_rrb_c_total","adi_r_soc_a_total", "rbs_r_overall_score", "srs_parent_raw_total", "vineland_ii_composite_standard_score", "dcdq_total", "scq_total")
IQ_list = list("ssc_diagnosis_full_scale_iq", "ssc_diagnosis_nonverbal_iq", "ssc_diagnosis_verbal_iq")


## Run regression in all

for(i in ados_list){
  
  results_all = summary(lm(paste0("scale(", i, ") ~ IQ_prs + autism_prs + edu_prs + scz_prs + anydenovo_constraint_sum + age_at_ados + ados_module + as.character(Sex) + as.character(dataset) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) + 
                                     scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17)"), data = merged))
  
  IQ_results = rbind(IQ_results, cbind(i, t(results_all$coefficients[2,])))
  autism_results = rbind(autism_results, cbind(i, t(results_all$coefficients[3,])))
  edu_results = rbind(edu_results, cbind(i, t(results_all$coefficients[4,])))
  scz_results = rbind(scz_results, cbind(i, t(results_all$coefficients[5,])))
  de_novo_results = rbind(de_novo_results, cbind(i, t(results_all$coefficients[6,])))
  
}


for(i in other_list){
  
  results_all = summary(lm(paste0("scale(", i, ") ~ IQ_prs + autism_prs + edu_prs + scz_prs + anydenovo_constraint_sum + age_at_ados + as.character(Sex) + as.character(dataset) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) + 
                                     scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17)"), data = merged))
  
  IQ_results = rbind(IQ_results, cbind(i, t(results_all$coefficients[2,])))
  autism_results = rbind(autism_results, cbind(i, t(results_all$coefficients[3,])))
  edu_results = rbind(edu_results, cbind(i, t(results_all$coefficients[4,])))
  scz_results = rbind(scz_results, cbind(i, t(results_all$coefficients[5,])))
  de_novo_results = rbind(de_novo_results, cbind(i, t(results_all$coefficients[6,])))
}



for(i in IQ_list){
  
  results_all = summary(lm(paste0("scale(", i, ") ~ IQ_prs + autism_prs + edu_prs + scz_prs + anydenovo_constraint_sum + age_at_ados + as.character(Sex) + as.character(ssc_diagnosis_verbal_iq_type) + as.character(dataset) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) + 
                                     scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17)"), data = merged))
  
  IQ_results = rbind(IQ_results, cbind(i, t(results_all$coefficients[2,])))
  autism_results = rbind(autism_results, cbind(i, t(results_all$coefficients[3,])))
  edu_results = rbind(edu_results, cbind(i, t(results_all$coefficients[4,])))
  scz_results = rbind(scz_results, cbind(i, t(results_all$coefficients[5,])))
  de_novo_results = rbind(de_novo_results, cbind(i, t(results_all$coefficients[6,])))
}

dcdq = fread("~/SFARI/Phenotypes/dcdq.csv")
dcdq_merged = merge(dcdq, merged, by = "individual")


results_all = summary(lm(scale(total) ~ IQ_prs + autism_prs + edu_prs + scz_prs + anydenovo_constraint +  as.character(Sex) + age_at_ados + as.character(dataset) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) + scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = dcdq_merged))

i = "dcdq"

IQ_results = rbind(IQ_results, cbind(i, t(results_all$coefficients[2,])))
autism_results = rbind(autism_results, cbind(i, t(results_all$coefficients[3,])))
edu_results = rbind(edu_results, cbind(i, t(results_all$coefficients[4,])))
scz_results = rbind(scz_results, cbind(i, t(results_all$coefficients[5,])))
de_novo_results = rbind(de_novo_results, cbind(i, t(results_all$coefficients[6,])))


scq = fread("~/SFARI/Phenotypes/scq_life.csv")
scq_merged = merge(scq, merged, by = "individual")
results_all = summary(lm(scale(summary_score) ~ IQ_prs + autism_prs + edu_prs + scz_prs + anydenovo_constraint + as.character(Sex) + age_at_ados + as.character(dataset) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) + scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = scq_merged))

i = "scq"

IQ_results = rbind(IQ_results, cbind(i, t(results_all$coefficients[2,])))
autism_results = rbind(autism_results, cbind(i, t(results_all$coefficients[3,])))
edu_results = rbind(edu_results, cbind(i, t(results_all$coefficients[4,])))
scz_results = rbind(scz_results, cbind(i, t(results_all$coefficients[5,])))
de_novo_results = rbind(de_novo_results, cbind(i, t(results_all$coefficients[6,])))


IQ_results = as.data.frame(IQ_results)
edu_results = as.data.frame(edu_results)
scz_results = as.data.frame(scz_results)
autism_results = as.data.frame(autism_results)
de_novo_results = as.data.frame(de_novo_results)

IQ_results$category = "IQ"
autism_results$category = "autism"
edu_results$category = "edu"
scz_results$category = "scz"
de_novo_results$category = "de_novo"

results_table = cbind(IQ_results, autism_results, edu_results, scz_results, de_novo_results)


save(results_table, file = "./Results/SSC_denovo_results.RData")







#############Run the analysis only in males###########################
males = subset(merged, sex == "male")

autism_results = NULL
scz_results = NULL
IQ_results = NULL
edu_results = NULL

for(i in ados_list){
  
  results_all = summary(lm(paste0("scale(", i, ") ~ IQ_prs + autism_prs + edu_prs + scz_prs + age_at_ados + ados_module + as.character(dataset) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) + 
                                     scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17)"), data = merged))
  
  IQ_results = rbind(IQ_results, cbind(i, t(results_all$coefficients[2,])))
  autism_results = rbind(autism_results, cbind(i, t(results_all$coefficients[3,])))
  edu_results = rbind(edu_results, cbind(i, t(results_all$coefficients[4,])))
  scz_results = rbind(scz_results, cbind(i, t(results_all$coefficients[5,])))
  
}


for(i in other_list){
  
  results_all = summary(lm(paste0("scale(", i, ") ~ IQ_prs + autism_prs + edu_prs + scz_prs + age_at_ados  + as.character(dataset) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) + 
                                     scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17)"), data = merged))
  
  IQ_results = rbind(IQ_results, cbind(i, t(results_all$coefficients[2,])))
  autism_results = rbind(autism_results, cbind(i, t(results_all$coefficients[3,])))
  edu_results = rbind(edu_results, cbind(i, t(results_all$coefficients[4,])))
  scz_results = rbind(scz_results, cbind(i, t(results_all$coefficients[5,])))
  
}



for(i in IQ_list){
  
  results_all = summary(lm(paste0("scale(", i, ") ~ IQ_prs + autism_prs + edu_prs + scz_prs + age_at_ados  + as.character(ssc_diagnosis_verbal_iq_type) + as.character(dataset) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) + 
                                     scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17)"), data = merged))
  
  IQ_results = rbind(IQ_results, cbind(i, t(results_all$coefficients[2,])))
  autism_results = rbind(autism_results, cbind(i, t(results_all$coefficients[3,])))
  edu_results = rbind(edu_results, cbind(i, t(results_all$coefficients[4,])))
  scz_results = rbind(scz_results, cbind(i, t(results_all$coefficients[5,])))
  
}






############################IQ > 70#########################


non_id = subset(merged, ssc_diagnosis_full_scale_iq > 70)

autism_results = NULL
scz_results = NULL
IQ_results = NULL
edu_results = NULL

for(i in ados_list){
  
  results_all = summary(lm(paste0("scale(", i, ") ~ IQ_prs + autism_prs + edu_prs + scz_prs + age_at_ados + ados_module + as.character(Sex) + as.character(dataset) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) + 
                                     scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17)"), data = non_id))
  
  IQ_results = rbind(IQ_results, cbind(i, t(results_all$coefficients[2,])))
  autism_results = rbind(autism_results, cbind(i, t(results_all$coefficients[3,])))
  edu_results = rbind(edu_results, cbind(i, t(results_all$coefficients[4,])))
  scz_results = rbind(scz_results, cbind(i, t(results_all$coefficients[5,])))
  
}


for(i in other_list){
  
  results_all = summary(lm(paste0("scale(", i, ") ~ IQ_prs + autism_prs + edu_prs + scz_prs + age_at_ados + as.character(Sex) + as.character(dataset) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) + 
                                     scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17)"), data = non_id))
  
  IQ_results = rbind(IQ_results, cbind(i, t(results_all$coefficients[2,])))
  autism_results = rbind(autism_results, cbind(i, t(results_all$coefficients[3,])))
  edu_results = rbind(edu_results, cbind(i, t(results_all$coefficients[4,])))
  scz_results = rbind(scz_results, cbind(i, t(results_all$coefficients[5,])))
  
}



for(i in IQ_list){
  
  results_all = summary(lm(paste0("scale(", i, ") ~ IQ_prs + autism_prs + edu_prs + scz_prs + age_at_ados + as.character(Sex) + as.character(ssc_diagnosis_verbal_iq_type) + as.character(dataset) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) + 
                                     scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17)"), data = non_id))
  
  IQ_results = rbind(IQ_results, cbind(i, t(results_all$coefficients[2,])))
  autism_results = rbind(autism_results, cbind(i, t(results_all$coefficients[3,])))
  edu_results = rbind(edu_results, cbind(i, t(results_all$coefficients[4,])))
  scz_results = rbind(scz_results, cbind(i, t(results_all$coefficients[5,])))
  
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




#### Multivariate regression ###


alpha = lm(cbind(scale(edu_prs), scale(IQ_PRS), scale (scz_1)) ~  scale(ssc_diagnosis_full_scale_iq)  + as.character(Sex) + as.character(dataset) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) + scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = merged)




###dcdq###

dcdq = fread("~/SFARI/Phenotypes/dcdq.csv")
dcdq_merged = merge(dcdq, merged, by = "individual")


summary(lm(scale(total) ~ IQ_prs + autism_prs + edu_prs + scz_prs + as.character(Sex) + age_at_ados + as.character(dataset) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) + scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = dcdq_merged))


###scq##

scq = fread("~/SFARI/Phenotypes/scq_life.csv")
scq_merged = merge(scq, merged, by = "individual")
summary(lm(scale(summary_score) ~ IQ_prs + autism_prs + edu_prs + scz_prs + as.character(Sex) + age_at_ados + as.character(dataset) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) + scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = scq_merged))



merged2 = merge(merged, dcdq, by = "individual")
merged2 = merge(merged2, scq, by = "individual")

###

data_for_matrix = merged2[,c("ados_social_affect", "ados_restricted_repetitive", "adi_r_b_comm_verbal_total",
                            "adi_r_rrb_c_total", "adi_r_soc_a_total", "rbs_r_overall_score", 
                            "vineland_ii_composite_standard_score", "srs_parent_raw_total", "ssc_diagnosis_full_scale_iq", 
                            "ssc_diagnosis_nonverbal_iq", "ssc_diagnosis_verbal_iq", "total", "summary_score")]

setnames(data_for_matrix, "total", "ssc_dcdq")
setnames(data_for_matrix, "summary_score", "ssc_scq")

data_for_matrix$ados_social_affect = as.numeric(as.character(data_for_matrix$ados_social_affect))
correlationmatrix = cor(data_for_matrix, use = "complete.obs", method = "pearson")


hM <- format(round(correlationmatrix, 2))
col<- colorRampPalette(c("blue", "white", "red"))(20)
heatmap.2(correlationmatrix, trace = "none", col=col, cellnote=hM, notecol = "black")

pdf(filename="test.pdf",width=20,height=8)

heatmap.2(correlationmatrix, trace = "none", col=col)

res2<-rcorr(as.matrix(data_for_matrix))
res2_results = flattenCorrMatrix(res2$r, res2$P)



###SSC_pTDT based analyses

Cases = fread("~/SFARI/cases.txt", header= T)
setnames(Cases, "id", "individual")

prs_scz = fread("./PGS/SSC_scz_midparent_ptdt.txt", header = TRUE)
prs_scz$scz_deviation = prs_scz$SCORE - prs_scz$midparent
setnames(prs_scz,"midparent", "scz_midparent")
prs_scz = prs_scz[,c("IID", "scz_midparent", "scz_deviation")]


prs_edu = fread("./PGS/SSC_edu_midparent_ptdt.txt", header = TRUE)
prs_edu$edu_deviation = prs_edu$SCORE - prs_edu$midparent
setnames(prs_edu,"midparent", "edu_midparent")
prs_edu = prs_edu[,c("IID", "edu_midparent", "edu_deviation")]

prs_IQ = fread("./PGS/SSC_IQ_midparent_ptdt.txt", header = TRUE)
prs_IQ$IQ_deviation = prs_IQ$SCORE - prs_IQ$midparent
setnames(prs_IQ,"midparent", "IQ_midparent")
prs_IQ = prs_IQ[,c("IID", "IQ_midparent", "IQ_deviation")]

prs_autism = fread("./PGS/SSC_autism_midparent_ptdt.txt", header = TRUE)
prs_autism$autism_deviation = prs_autism$SCORE - prs_autism$midparent
setnames(prs_autism,"midparent", "autism_midparent")
prs_autism = prs_autism[,c("IID", "autism_midparent", "autism_deviation")]

prs_merged = merge(prs_IQ, prs_edu, by = "IID")
prs_merged = merge(prs_merged, prs_scz, by = "IID")
prs_merged = merge(prs_merged, prs_autism, by = "IID")

pheno = fread("~/SFARI/ssccore.txt")

pca = fread("~/SFARI/liftOverPlink/files_imputed/SSC_cases_pca.eigenvec")
setnames(pca, 2, "IID")

merged = merge(prs_merged, pca, by = "IID")
merged = merge(Cases, merged, by = "IID")
merged = merge(merged, pheno, by = "individual")

merged$ados_social_affect = as.numeric(as.character(merged$ados_social_affect))

merged$autism_midparent = scale(merged$autism_midparent)
merged$autism_deviation = scale(merged$autism_deviation)

merged$scz_midparent = scale(merged$scz_midparent)
merged$scz_deviation = scale(merged$scz_deviation)

merged$IQ_midparent = scale(merged$IQ_midparent)
merged$IQ_deviation = scale(merged$IQ_deviation)

merged$edu_midparent = scale(merged$edu_midparent)
merged$edu_deviation = scale(merged$edu_deviation)


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

autism_deviation_results = NULL
scz_deviation_results = NULL
IQ_deviation_results = NULL
edu_deviation_results = NULL

autism_midparent_results = NULL
scz_midparent_results = NULL
IQ_midparent_results = NULL
edu_midparent_results = NULL
de_novo_results = NULL

ados_list = list("ados_restricted_repetitive", "ados_social_affect")
other_list = list("adi_r_b_comm_verbal_total", "adi_r_rrb_c_total","adi_r_soc_a_total", "rbs_r_overall_score", "srs_parent_raw_total", "vineland_ii_composite_standard_score")
IQ_list = list("ssc_diagnosis_full_scale_iq", "ssc_diagnosis_nonverbal_iq", "ssc_diagnosis_verbal_iq")


## Run regression in all

for(i in ados_list){
  
  results_all = summary(lm(paste0("scale(", i, ") ~ IQ_midparent + autism_midparent + edu_midparent + scz_midparent + IQ_deviation + autism_deviation + edu_deviation + scz_deviation + age_at_ados + ados_module + as.character(Sex) + as.character(dataset) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) + 
                                     scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17)"), data = merged))
  
  IQ_midparent_results = rbind(IQ_midparent_results, cbind(i, t(results_all$coefficients[2,])))
  autism_midparent_results = rbind(autism_midparent_results, cbind(i, t(results_all$coefficients[3,])))
  edu_midparent_results = rbind(edu_midparent_results, cbind(i, t(results_all$coefficients[4,])))
  scz_midparent_results = rbind(scz_midparent_results, cbind(i, t(results_all$coefficients[5,])))
  
  IQ_deviation_results = rbind(IQ_deviation_results, cbind(i, t(results_all$coefficients[6,])))
  autism_deviation_results = rbind(autism_deviation_results, cbind(i, t(results_all$coefficients[7,])))
  edu_deviation_results = rbind(edu_deviation_results, cbind(i, t(results_all$coefficients[8,])))
  scz_deviation_results = rbind(scz_deviation_results, cbind(i, t(results_all$coefficients[9,])))
  
}


for(i in other_list){
  
  results_all = summary(lm(paste0("scale(", i, ") ~ IQ_midparent + autism_midparent + edu_midparent + scz_midparent + IQ_deviation + autism_deviation + edu_deviation + scz_deviation + age_at_ados + as.character(Sex) + as.character(dataset) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) + 
                                     scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17)"), data = merged))
  
  IQ_midparent_results = rbind(IQ_midparent_results, cbind(i, t(results_all$coefficients[2,])))
  autism_midparent_results = rbind(autism_midparent_results, cbind(i, t(results_all$coefficients[3,])))
  edu_midparent_results = rbind(edu_midparent_results, cbind(i, t(results_all$coefficients[4,])))
  scz_midparent_results = rbind(scz_midparent_results, cbind(i, t(results_all$coefficients[5,])))
  
  IQ_deviation_results = rbind(IQ_deviation_results, cbind(i, t(results_all$coefficients[6,])))
  autism_deviation_results = rbind(autism_deviation_results, cbind(i, t(results_all$coefficients[7,])))
  edu_deviation_results = rbind(edu_deviation_results, cbind(i, t(results_all$coefficients[8,])))
  scz_deviation_results = rbind(scz_deviation_results, cbind(i, t(results_all$coefficients[9,])))
  
}



for(i in IQ_list){
  
  results_all = summary(lm(paste0("scale(", i, ") ~ IQ_midparent + autism_midparent + edu_midparent + scz_midparent + IQ_deviation + autism_deviation + edu_deviation + scz_deviation + age_at_ados + as.character(Sex) + as.character(ssc_diagnosis_verbal_iq_type) + as.character(dataset) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) + 
                                     scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17)"), data = merged))
  
  IQ_midparent_results = rbind(IQ_midparent_results, cbind(i, t(results_all$coefficients[2,])))
  autism_midparent_results = rbind(autism_midparent_results, cbind(i, t(results_all$coefficients[3,])))
  edu_midparent_results = rbind(edu_midparent_results, cbind(i, t(results_all$coefficients[4,])))
  scz_midparent_results = rbind(scz_midparent_results, cbind(i, t(results_all$coefficients[5,])))
  
  IQ_deviation_results = rbind(IQ_deviation_results, cbind(i, t(results_all$coefficients[6,])))
  autism_deviation_results = rbind(autism_deviation_results, cbind(i, t(results_all$coefficients[7,])))
  edu_deviation_results = rbind(edu_deviation_results, cbind(i, t(results_all$coefficients[8,])))
  scz_deviation_results = rbind(scz_deviation_results, cbind(i, t(results_all$coefficients[9,])))
  
}




###With de novo (only IQ list)

autism_deviation_results = NULL
scz_deviation_results = NULL
IQ_deviation_results = NULL
edu_deviation_results = NULL

autism_midparent_results = NULL
scz_midparent_results = NULL
IQ_midparent_results = NULL
edu_midparent_results = NULL
de_novo_results = NULL

for(i in IQ_list){
  
  results_all = summary(lm(paste0("scale(", i, ") ~ IQ_midparent + autism_midparent + edu_midparent + scz_midparent + IQ_deviation + autism_deviation + edu_deviation + scz_deviation + anydenovo_constraint + age_at_ados + as.character(Sex) + as.character(ssc_diagnosis_verbal_iq_type) + as.character(dataset) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) + 
                                     scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17)"), data = merged))
  
  IQ_midparent_results = rbind(IQ_midparent_results, cbind(i, t(results_all$coefficients[2,])))
  autism_midparent_results = rbind(autism_midparent_results, cbind(i, t(results_all$coefficients[3,])))
  edu_midparent_results = rbind(edu_midparent_results, cbind(i, t(results_all$coefficients[4,])))
  scz_midparent_results = rbind(scz_midparent_results, cbind(i, t(results_all$coefficients[5,])))
  
  IQ_deviation_results = rbind(IQ_deviation_results, cbind(i, t(results_all$coefficients[6,])))
  autism_deviation_results = rbind(autism_deviation_results, cbind(i, t(results_all$coefficients[7,])))
  edu_deviation_results = rbind(edu_deviation_results, cbind(i, t(results_all$coefficients[8,])))
  scz_deviation_results = rbind(scz_deviation_results, cbind(i, t(results_all$coefficients[9,])))
  de_novo_results = rbind(de_novo_results, cbind(i, t(results_all$coefficients[10,])))
  
}


results_all = summary(lm(scale(ssc_diagnosis_nonverbal_iq) ~ IQ_midparent*anydenovo_constraint_sum + autism_midparent + edu_midparent + scz_midparent + IQ_deviation + autism_deviation + edu_deviation + scz_deviation  + age_at_ados + as.character(Sex) + as.character(ssc_diagnosis_verbal_iq_type) + as.character(dataset) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) + 
                                     scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = merged))

##Plot the data

midparent = ggplot(merged, aes(x= ssc_diagnosis_nonverbal_iq, y = IQ_midparent, colour = as.character(anydenovo_constraint))) +  geom_smooth(method=lm, se=TRUE) + theme_classic()

transmission = ggplot(merged, aes(x= ssc_diagnosis_nonverbal_iq, y = IQ_deviation, colour = as.character(anydenovo_constraint))) +  geom_smooth(method=lm, se=TRUE) + theme_classic()




####Interaction with sex###
interaction_table = NULL
list1 =  list("ssc_diagnosis_full_scale_iq", "ssc_diagnosis_nonverbal_iq", "ssc_diagnosis_verbal_iq", "vineland_ii_composite_standard_score")

for(i in list1){
  results_all = summary(lm(paste0("scale(", i, ") ~ IQ_prs + autism_prs + edu_prs + anydenovo_constraint_sum*Sex + scale(V3) + scale(V4) + scale(V5) + scale(V6)"), data = merged))
  
  interaction_table = rbind(interaction_table, cbind(i, t(results_all$coefficients[11,])))
  
}

