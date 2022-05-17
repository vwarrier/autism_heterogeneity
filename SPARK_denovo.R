##de novo analyses in SPARK

library(MASS)
library(data.table)
library(ggplot2)
library(dplyr)
library(lme4)
library(car)


##All genetic
prs_scz = fread("./PGS/SPARK_scz_finalscore.profile", header = TRUE)
setnames(prs_scz, "SCORE", "scz_prs")
prs_scz = prs_scz[,c("IID", "scz_prs")]


prs_edu = fread("./PGS/SPARK_edu_finalscore.profile", header = TRUE)
setnames(prs_edu, "SCORE", "edu_prs")
prs_edu = prs_edu[,c("IID", "edu_prs")]

prs_IQ = fread("./PGS/SPARK_IQ_finalscore.profile", header = TRUE)
setnames(prs_IQ, "SCORE", "IQ_prs")
prs_IQ = prs_IQ[,c("IID", "IQ_prs")]

prs_autism = fread("./PGS/SPARK_autism_finalscore.profile", header = TRUE)
setnames(prs_autism, "SCORE", "autism_prs")
prs_autism = prs_autism[,c("IID", "autism_prs", "FID")]


prs_merged = merge(prs_IQ, prs_edu, by = "IID")
prs_merged = merge(prs_merged, prs_scz, by = "IID")
prs_merged = merge(prs_merged, prs_autism, by = "IID")

PCS = fread("~/SPARK/SPARK_postimputation/CEU/SPARK_CEUPCsforGWAS.txt", header = TRUE)
setnames(PCS, "Sample_name", "IID")

merged = merge(prs_merged, PCS, by = "IID")

##All_phenotypic

pheno = fread("./GRMs/SPARK_phenos.txt")
pheno = pheno[,-1]

factor = fread("./GRMs/factor_scores.txt")
factor = factor[,-1]

pheno = merge(pheno, factor, by = "V2", all.x = TRUE)
setnames(pheno, "V2", "IID")

registration = fread("~/SPARK/Phenotypes/V5/individuals_registration.csv", fill = TRUE)
trios =  registration[!(registration$biofather_sp_id =="" | registration$biomother_sp_id==""), ] # 52936
asd_trios = subset(trios, asd == "TRUE") #35723

asd_trios = asd_trios[,c("subject_sp_id", "age_at_registration_months", "sex", "biomother_sp_id", "biofather_sp_id")]

merged_pgs_trio = merge(merged, asd_trios, by.x = "IID", by.y = "subject_sp_id") #4559

merged_pgs_trio = merged_pgs_trio[merged_pgs_trio$biomother_sp_id %in% merged$IID,] #3904
merged_pgs_trio = merged_pgs_trio[merged_pgs_trio$biofather_sp_id %in% merged$IID,] #3342

merged_final = merge(merged_pgs_trio, pheno, by = "IID", all.x  = TRUE)

##sibling_trios

sibling_trios = subset(trios, asd == "FALSE")
sibling_trios = sibling_trios[,c("subject_sp_id", "age_at_registration_months", "sex", "biomother_sp_id", "biofather_sp_id")]

merged_pgs_trio_sibling = merge(merged, sibling_trios, by.x = "IID", by.y = "subject_sp_id") #4559

merged_pgs_trio_sibling = merged_pgs_trio_sibling[merged_pgs_trio_sibling$biomother_sp_id %in% merged$IID,] #3904
merged_pgs_trio_sibling = merged_pgs_trio_sibling[merged_pgs_trio_sibling$biofather_sp_id %in% merged$IID,] #3342



###

SNV = fread("./Sebat_data/denovo_SNVindel.txt")
SNV = subset(SNV, cohort == "SPARK")
ddd = fread("~/Autism_heterogeneity/DDD_genes.txt", header = F)


##loeuf_upper_decile_score = 0.37
SNV_upper_decile = subset(SNV, gnomad_loeuf < 0.37)
SNV_ddd = SNV_upper_decile[SNV_upper_decile$gene %in% ddd$V1,]

##Count_data
count1 = as.data.table(table(SNV$iid))
setnames(count1, old = c("V1", "N"), new = c("IID", "anydenovo_sum"))
merged_final <- merge(merged_final, count1, by = "IID", all.x = TRUE)
merged_final$anydenovo_sum[is.na(merged$anydenovo_sum)] <- 0

merged_sibling_final = merge(merged_pgs_trio_sibling, count1, by = "IID", all.x = TRUE)
merged_sibling_final$anydenovo_sum[is.na(merged_sibling_final$anydenovo_sum)] <- 0

count1 = as.data.table(table(SNV_upper_decile$iid))
setnames(count1, old = c("V1", "N"), new = c("IID", "anydenovo_constraint_sum"))
merged_final <- merge(merged_final, count1, by = "IID", all.x = TRUE)
merged_final$anydenovo_constraint_sum[is.na(merged_final$anydenovo_constraint_sum)] <- 0

merged_sibling_final <- merge(merged_sibling_final, count1, by = "IID", all.x = TRUE)
merged_sibling_final$anydenovo_constraint_sum[is.na(merged_sibling_final$anydenovo_constraint_sum)] <- 0


count1 = as.data.table(table(SNV_ddd$iid))
setnames(count1, old = c("V1", "N"), new = c("IID", "ddd_constraint_sum"))
merged_final <- merge(merged_final, count1, by = "IID", all.x = TRUE)
merged_final$ddd_constraint_sum[is.na(merged_final$ddd_constraint_sum)] <- 0

merged_sibling_final <- merge(merged_sibling_final, count1, by = "IID", all.x = TRUE)
merged_sibling_final$ddd_constraint_sum[is.na(merged_sibling_final$ddd_constraint_sum)] <- 0


merged_final$nonddd_constraint_sum = merged_final$anydenovo_constraint_sum - merged_final$ddd_constraint_sum
merged_final$nonndd_constraint_sum_only = ifelse(merged_final$ddd_constraint_sum > 0, 0, merged_final$nonddd_constraint_sum)

merged_sibling_final$nonddd_constraint_sum = merged_sibling_final$anydenovo_constraint_sum - merged_sibling_final$ddd_constraint_sum
merged_sibling_final$nonndd_constraint_sum_only = ifelse(merged_sibling_final$ddd_constraint_sum > 0, 0, merged_sibling_final$nonddd_constraint_sum)

merged_final$anydenovo_constraint = ifelse(merged_final$anydenovo_constraint_sum > 0, 1, 0)
merged_sibling_final$anydenovo_constraint = ifelse(merged_sibling_final$anydenovo_constraint_sum > 0, 1, 0)


save(merged_final, file = "~/Autism_heterogeneity/SPARK_denovo_dataforanalysis.Rdata")
save(merged_sibling_final, file = "~/Autism_heterogeneity/SPARK_sibling_denovo_dataforanalysis.Rdata")

###diff_between siblings and high-impact de novo carriers
merged_final$category = "case"
merged_sibling_final$category = "sibling"

merged1 = merged_final[,c("IQ_prs", "autism_prs", "scz_prs", "edu_prs", "sex", "X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8", "X9", "X10", "anydenovo_constraint", "IID", "category")]
merged2 = merged_sibling_final[,c("IQ_prs", "autism_prs", "scz_prs", "edu_prs", "sex",  "X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8", "X9", "X10", "anydenovo_constraint", "IID", "category")]
merged1_denovo = subset(merged1, anydenovo_constraint == "1")
merged_sibling_denovo = rbind(merged2, merged1_denovo)
summary(lm(scale(autism_prs) ~ category + sex + X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10, data = merged_sibling_denovo))

###Association analyses##
load("~/Autism_heterogeneity/SPARK_denovo_dataforanalysis.Rdata")


autism_results = NULL
scz_results = NULL
IQ_results = NULL
edu_results = NULL
de_novo_results = NULL

list1 = list("rbs_r", "scq", "PA1", "PA2", "PA3", "PA4", "PA5", "PA6", "vabs", "dcdq", "IQ_score")



## Run regression in all

for(i in list1){
  
  results_all = summary(lm(paste0("scale(", i, ") ~ scale(`IQ_prs`) + scale(autism_prs) + scale(scz_prs) + scale(edu_prs) + anydenovo_constraint_sum  + sex  + scale(age_at_registration_months) + scale(X1) + scale(X2) + scale(X3) + scale(X4) + scale(X5) + scale(X6) + scale(X7) + scale(X8) + scale(X9) + scale(X10)"), data = merged_final))
  
  IQ_results = rbind(IQ_results, cbind(i, t(results_all$coefficients[2,])))
  autism_results = rbind(autism_results, cbind(i, t(results_all$coefficients[3,])))
  edu_results = rbind(edu_results, cbind(i, t(results_all$coefficients[4,])))
  scz_results = rbind(scz_results, cbind(i, t(results_all$coefficients[5,])))
  de_novo_results = rbind(de_novo_results, cbind(i, t(results_all$coefficients[6,])))
  
}

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


print(results_table)

save(results_table, file = "./Results/SPARK_denovo_results_autisticsiblingsallowed.RData")

##Without  de_novo

for(i in list1){
  
  results_all = summary(lm(paste0("scale(", i, ") ~ scale(`IQ_prs`) + scale(autism_prs) + scale(scz_prs) + scale(edu_prs)  + sex  + scale(age_at_registration_months) + scale(X1) + scale(X2) + scale(X3) + scale(X4) + scale(X5) + scale(X6) + scale(X7) + scale(X8) + scale(X9) + scale(X10)"), data = merged_final))
  
  IQ_results = rbind(IQ_results, cbind(i, t(results_all$coefficients[2,])))
  autism_results = rbind(autism_results, cbind(i, t(results_all$coefficients[3,])))
  edu_results = rbind(edu_results, cbind(i, t(results_all$coefficients[4,])))
  scz_results = rbind(scz_results, cbind(i, t(results_all$coefficients[5,])))
  
}

IQ_results = as.data.frame(IQ_results)
edu_results = as.data.frame(edu_results)
scz_results = as.data.frame(scz_results)
autism_results = as.data.frame(autism_results)

IQ_results$category = "IQ"
autism_results$category = "autism"
edu_results$category = "edu"
scz_results$category = "scz"

results_table = cbind(IQ_results, autism_results, edu_results, scz_results)

print(results_table)

save(results_table, file = "./Results/SPARK_denovo_results_autisticsiblingsallowed_withoutdenovo.RData")


###After accounting for IQ

autism_results = NULL
scz_results = NULL
IQ_results = NULL
edu_results = NULL
de_novo_results = NULL

list1 = list("rbs_r", "scq", "PA1", "PA2", "PA3", "PA4", "PA5", "PA6", "vabs", "dcdq")



## Run regression in all

for(i in list1){
  
  results_all = summary(lm(paste0("scale(", i, ") ~ scale(`IQ_prs`) + scale(autism_prs) + scale(scz_prs) + scale(edu_prs) + anydenovo_constraint_sum  + scale(IQ_score) + sex  + scale(age_at_registration_months) + scale(X1) + scale(X2) + scale(X3) + scale(X4) + scale(X5) + scale(X6) + scale(X7) + scale(X8) + scale(X9) + scale(X10)"), data = merged_final))
  
  IQ_results = rbind(IQ_results, cbind(i, t(results_all$coefficients[2,])))
  autism_results = rbind(autism_results, cbind(i, t(results_all$coefficients[3,])))
  edu_results = rbind(edu_results, cbind(i, t(results_all$coefficients[4,])))
  scz_results = rbind(scz_results, cbind(i, t(results_all$coefficients[5,])))
  de_novo_results = rbind(de_novo_results, cbind(i, t(results_all$coefficients[6,])))
  
}

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

results_table_noIQ = cbind(IQ_results, autism_results, edu_results, scz_results, de_novo_results)




###Dividing it into the two gene lists

autism_results = NULL
scz_results = NULL
IQ_results = NULL
edu_results = NULL
de_novoddd_results = NULL
de_novononddd_results = NULL

list1 = list("rbs_r", "scq", "PA1", "PA2", "PA3", "PA4", "PA5", "PA6", "vabs", "dcdq", "IQ_score")



## Run regression in all

for(i in list1){
  
  results_all = summary(lm(paste0("scale(", i, ") ~ scale(`IQ_prs`) + scale(autism_prs) + scale(scz_prs) + scale(edu_prs) + ddd_constraint_sum  + nonndd_constraint_sum_only +  sex  + scale(age_at_registration_months) + scale(X1) + scale(X2) + scale(X3) + scale(X4) + scale(X5) + scale(X6) + scale(X7) + scale(X8) + scale(X9) + scale(X10)"), data = merged_final))
  
  IQ_results = rbind(IQ_results, cbind(i, t(results_all$coefficients[2,])))
  autism_results = rbind(autism_results, cbind(i, t(results_all$coefficients[3,])))
  edu_results = rbind(edu_results, cbind(i, t(results_all$coefficients[4,])))
  scz_results = rbind(scz_results, cbind(i, t(results_all$coefficients[5,])))
  de_novoddd_results = rbind(de_novoddd_results, cbind(i, t(results_all$coefficients[6,])))
  de_novononddd_results = rbind(de_novononddd_results, cbind(i, t(results_all$coefficients[7,])))
  
}

IQ_results = as.data.frame(IQ_results)
edu_results = as.data.frame(edu_results)
scz_results = as.data.frame(scz_results)
autism_results = as.data.frame(autism_results)
de_novoddd_results = as.data.frame(de_novoddd_results)
de_novononddd_results = as.data.frame(de_novononddd_results)

IQ_results$category = "IQ"
autism_results$category = "autism"
edu_results$category = "edu"
scz_results$category = "scz"
de_novoddd_results$category = "de_novo_ddd"
de_novononddd_results$category = "de_novo_nonddd"


results_table = cbind(IQ_results, autism_results, edu_results, scz_results, de_novoddd_results, de_novononddd_results)


print(results_table)




# Are de novo mutations lower in individuals with high autism PGS?


summary(glm(anydenovo_constraint ~ scale(`IQ_prs`) + scale(autism_prs) + scale(scz_prs) + scale(edu_prs) + sex  + scale(age_at_registration_months) + scale(X1) + scale(X2) + scale(X3) + scale(X4) + scale(X5) + scale(X6) + scale(X7) + scale(X8) + scale(X9) + scale(X10), data = merged_final, family = "binomial"))


##By quantile

merged_final = merged_final %>% mutate(quantile = ntile(autism_prs, 10))


summary(glm(anydenovo_constraint ~ as.factor(quantile) + sex  + scale(age_at_registration_months) + scale(X1) + scale(X2) + scale(X3) + scale(X4) + scale(X5), data = merged_final, family = "binomial"))





###medical_diagnosis

###medical_history
med = fread("~/SPARK/Phenotypes/V5/basic_medical_screening.csv")
med[is.na(med)] <- 0 
setnames(med, 1, "IID")
med_asd = subset(med,asd == "TRUE")
merged1 = merge(med_asd, merged_final, by = "IID")
merged2 = merge(med, merged_sibling_final, by = "IID")

merged1$ddd_constraint = ifelse(merged1$ddd_constraint_sum > 0, 1, 0)
merged2$ddd_constraint = ifelse(merged2$ddd_constraint_sum > 0, 1, 0)

merged1$nonddd_constraint = ifelse(merged1$nonndd_constraint_sum_only > 0, 1, 0)
merged2$nonddd_constraint = ifelse(merged2$nonndd_constraint_sum_only > 0, 1, 0)



merged1$dev_diagnosis = merged1$dev_id.x + merged1$dev_lang_dis + merged1$dev_ld + merged1$dev_motor + merged1$dev_mutism + merged1$dev_soc_prag + merged1$dev_speech


summary(glm(dev_diagnosis ~ scale(`IQ_prs`) + scale(autism_prs) + scale(scz_prs) + scale(edu_prs) + anydenovo_constraint  + sex.x  + scale(age_at_eval_months) + scale(X1) + scale(X2) + scale(X3) + scale(X4) + scale(X5) + scale(X6) + scale(X7) + scale(X8) + scale(X9) + scale(X10), data = merged1, family = "quasipoisson"))
summary(glm(dev_diagnosis ~ scale(`IQ_prs`) + scale(autism_prs) + scale(scz_prs) + scale(edu_prs) + ddd_constraint  + nonddd_constraint + sex.x  + scale(age_at_eval_months) + scale(X1) + scale(X2) + scale(X3) + scale(X4) + scale(X5) + scale(X6) + scale(X7) + scale(X8) + scale(X9) + scale(X10), data = merged1, family = "quasipoisson"))

merged1$dev_diagnosis_2 = ifelse(merged1$dev_diagnosis == 0, "zero", NA)
merged1$dev_diagnosis_2 = ifelse(merged1$dev_diagnosis > 0 & merged1$dev_diagnosis < 3, "one_plus", merged1$dev_diagnosis_2)
merged1$dev_diagnosis_2 = ifelse(merged1$dev_diagnosis > 2 & merged1$dev_diagnosis < 5, "three_plus", merged1$dev_diagnosis_2)
merged1$dev_diagnosis_2 = ifelse(merged1$dev_diagnosis > 4 , "fiveplus", merged1$dev_diagnosis_2)

merged2$dev_diagnosis_2 = "a_sibling"

merged1_a = merged1[,c("IQ_prs", "autism_prs", "scz_prs", "edu_prs", "sex.x", "age_at_eval_months", "X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8", "X9", "X10", "anydenovo_constraint", "IID", "dev_diagnosis_2", "ddd_constraint", "nonddd_constraint")]
merged2_a = merged2[,c("IQ_prs", "autism_prs", "scz_prs", "edu_prs", "sex.x", "age_at_eval_months", "X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8", "X9", "X10", "anydenovo_constraint", "IID", "dev_diagnosis_2", "ddd_constraint", "nonddd_constraint")]

merged_aut_sibling = rbind(merged1_a, merged2_a)

summary(glm(anydenovo_constraint  ~ as.factor(dev_diagnosis_2) + scale(`IQ_prs`) + scale(autism_prs) + scale(scz_prs) + scale(edu_prs) + sex.x  + scale(age_at_eval_months) + scale(X1) + scale(X2) + scale(X3) + scale(X4) + scale(X5) + scale(X6) + scale(X7) + scale(X8) + scale(X9) + scale(X10), data = merged_aut_sibling, family = "binomial"))

summary(glm(ddd_constraint  ~ as.factor(dev_diagnosis_2) + scale(`IQ_prs`) + scale(autism_prs) + scale(scz_prs) + scale(edu_prs) + sex.x  + scale(age_at_eval_months) + scale(X1) + scale(X2) + scale(X3) + scale(X4) + scale(X5) + scale(X6) + scale(X7) + scale(X8) + scale(X9) + scale(X10), data = merged_aut_sibling, family = "binomial"))

summary(glm(nonddd_constraint  ~ as.factor(dev_diagnosis_2) + scale(`IQ_prs`) + scale(autism_prs) + scale(scz_prs) + scale(edu_prs) + sex.x  + scale(age_at_eval_months) + scale(X1) + scale(X2) + scale(X3) + scale(X4) + scale(X5) + scale(X6) + scale(X7) + scale(X8) + scale(X9) + scale(X10), data = merged_aut_sibling, family = "binomial"))






###Sibling
med = fread("~/SPARK/Phenotypes/V5/basic_medical_screening.csv")
med[is.na(med)] <- 0 
setnames(med, 1, "IID")
med_nonasd = subset(med, asd == "FALSE")
merged2 = merge(med_nonasd, merged_sibling_final, by = "IID")





########################## Age of walkinge etc

bgx_child = fread("~/SPARK/Phenotypes/V5/background_history_child.csv")
setnames(bgx_child, "subject_sp_id", "IID")

merged1 = merge(merged_final, bgx_child, by = "IID")
merged1_asd = subset(merged1, asd == "TRUE") #2579

merged1_asd$walked_age_mos = ifelse(merged1_asd$walked_age_mos == "888", NA, merged1_asd$walked_age_mos)
merged1_asd$smiled_age_mos = ifelse(merged1_asd$smiled_age_mos == "888", NA, merged1_asd$smiled_age_mos)
merged1_asd$used_words_age_mos = ifelse(merged1_asd$used_words_age_mos == "888", NA, merged1_asd$used_words_age_mos)
merged1_asd$crawled_age_mos = ifelse(merged1_asd$crawled_age_mos == "888", NA, merged1_asd$crawled_age_mos)
merged1_asd$sat_wo_support_age_mos = ifelse(merged1_asd$sat_wo_support_age_mos == "888", NA, merged1_asd$sat_wo_support_age_mos)


summary(glm.nb(used_words_age_mos ~ scale(anydenovo_constraint_sum )+ scale(`IQ_prs`) + scale(autism_prs) + scale(scz_prs) + scale(edu_prs)  + sex.x  + scale(age_at_eval_months) + scale(X1) + scale(X2) + scale(X3) + scale(X4) + scale(X5) + scale(X6) + scale(X7) + scale(X8) + scale(X9) + scale(X10) + scale(X11) + scale(X12) + scale(X13) + scale(X14), data = merged1_asd))
summary(glm.nb(walked_age_mos ~ scale(anydenovo_constraint_sum )+ scale(`IQ_prs`) + scale(autism_prs) + scale(scz_prs) + scale(edu_prs)  + sex.x  + scale(age_at_eval_months) + scale(X1) + scale(X2) + scale(X3) + scale(X4) + scale(X5) + scale(X6) + scale(X7) + scale(X8) + scale(X9) + scale(X10) + scale(X11) + scale(X12) + scale(X13) + scale(X14), data = merged1_asd))


summary(glm.nb(walked_age_mos ~ scale(anydenovo_constraint_sum )+ scale(IQ_score) + scale(`IQ_prs`) + scale(autism_prs) + scale(scz_prs) + scale(edu_prs)  + sex.x  + scale(age_at_eval_months) + scale(X1) + scale(X2) + scale(X3) + scale(X4) + scale(X5) + scale(X6) + scale(X7) + scale(X8) + scale(X9) + scale(X10) + scale(X11) + scale(X12) + scale(X13) + scale(X14), data = merged1_asd))

summary(glm.nb(used_words_age_mos ~ scale(anydenovo_constraint_sum )+ scale(IQ_score) + scale(`IQ_prs`) + scale(autism_prs) + scale(scz_prs) + scale(edu_prs)  + sex.x  + scale(age_at_eval_months) + scale(X1) + scale(X2) + scale(X3) + scale(X4) + scale(X5) + scale(X6) + scale(X7) + scale(X8) + scale(X9) + scale(X10) + scale(X11) + scale(X12) + scale(X13) + scale(X14), data = merged1_asd))



#in sibling
bgx_sibling = fread("~/SPARK/Phenotypes/V5/background_history_sibling.csv")
setnames(bgx_sibling, "subject_sp_id", "IID")
merged1 = merge(merged_sibling_final, bgx_sibling, by = "IID")
merged1_nonasd = subset(merged1, asd == "FALSE") #1745
merged1_nonasd$walked_age_mos = ifelse(merged1_nonasd$walked_age_mos == "888", NA, merged1_nonasd$walked_age_mos)
merged1_nonasd$used_words_age_mos = ifelse(merged1_nonasd$used_words_age_mos == "888", NA, merged1_nonasd$used_words_age_mos)

summary(lm(used_words_age_mos ~ scale(anydenovo_constraint_sum )+ scale(`IQ_prs`) + scale(autism_prs) + scale(scz_prs) + scale(edu_prs)  + sex.x  + scale(age_at_eval_months) + scale(X1) + scale(X2) + scale(X3) + scale(X4) + scale(X5) + scale(X6) + scale(X7) + scale(X8) + scale(X9) + scale(X10) + scale(X11) + scale(X12) + scale(X13) + scale(X14), data = merged1_nonasd))
summary(lm(walked_age_mos ~ scale(anydenovo_constraint_sum )+ scale(`IQ_prs`) + scale(autism_prs) + scale(scz_prs) + scale(edu_prs)  + sex.x  + scale(age_at_eval_months) + scale(X1) + scale(X2) + scale(X3) + scale(X4) + scale(X5) + scale(X6) + scale(X7) + scale(X8) + scale(X9) + scale(X10) + scale(X11) + scale(X12) + scale(X13) + scale(X14), data = merged1_nonasd))






##SCQ
merged1 = merge(scq, merged, by = "IID")
merged1 = subset(merged1, asd == "TRUE")

summary(lmer(scale(final_score) ~ scale(`IQ_prs`) + scale(autism_prs) + scale(scz_prs) + scale(edu_prs) + anydenovo_constraint_sum + (1|family_sf_id) + sex  + scale(age_at_eval_months) + scale(X1) + scale(X2) + scale(X3) + scale(X4) + scale(X5) + scale(X6) + scale(X7) + scale(X8) + scale(X9) + scale(X10) + scale(X11) + scale(X12) + scale(X13) + scale(X14), data = merged1))


##RBS
merged1 = merge(rbs, merged, by = "IID")

summary(lmer(scale(total_final_score) ~ scale(IQ_prs) + scale(autism_prs) + scale(scz_prs) + scale(edu_prs) + anydenovo_constraint_sum + (1|family_sf_id) + sex  + scale(age_at_eval_months) + scale(X1) + scale(X2) + scale(X3) + scale(X4) + scale(X5) + scale(X6) + scale(X7) + scale(X8) + scale(X9) + scale(X10) + scale(X11) + scale(X12) + scale(X13) + scale(X14), data = merged1))

summary(lmer(scale(i_stereotyped_behavior_score) ~ scale(`IQ_prs`) + scale(autism_prs) + scale(scz_prs) + scale(edu_prs) + (1|FID.x) + sex  + scale(age_at_eval_months) + scale(X1) + scale(X2) + scale(X3) + scale(X4) + scale(X5) + scale(X6) + scale(X7) + scale(X8) + scale(X9) + scale(X10) + scale(X11) + scale(X12) + scale(X13) + scale(X14), data = merged1))
summary(lmer(scale(ii_self_injurious_score) ~ scale(`IQ_prs`) + scale(autism_prs) + scale(scz_prs) + scale(edu_prs) + (1|FID) + sex  + scale(age_at_eval_months) + scale(X1) + scale(X2) + scale(X3) + scale(X4) + scale(X5) + scale(X6) + scale(X7) + scale(X8) + scale(X9) + scale(X10) + scale(X11) + scale(X12) + scale(X13) + scale(X14), data = merged1))
summary(lmer(scale(iii_compulsive_behavior_score) ~ scale(`IQ_prs`) + scale(autism_prs) + scale(scz_prs) + scale(edu_prs) + (1|FID) + sex  + scale(age_at_eval_months) + scale(X1) + scale(X2) + scale(X3) + scale(X4) + scale(X5) + scale(X6) + scale(X7) + scale(X8) + scale(X9) + scale(X10) + scale(X11) + scale(X12) + scale(X13) + scale(X14), data = merged1))
summary(lmer(scale(iv_ritualistic_behavior_score) ~ scale(`IQ_prs`) + scale(autism_prs) + scale(scz_prs) + scale(edu_prs) + (1|FID) + sex  + scale(age_at_eval_months) + scale(X1) + scale(X2) + scale(X3) + scale(X4) + scale(X5) + scale(X6) + scale(X7) + scale(X8) + scale(X9) + scale(X10) + scale(X11) + scale(X12) + scale(X13) + scale(X14), data = merged1))
summary(lmer(scale(v_sameness_behavior_score) ~ scale(`IQ_prs`) + scale(autism_prs) + scale(scz_prs) + scale(edu_prs) + (1|FID) + sex  + scale(age_at_eval_months) + scale(X1) + scale(X2) + scale(X3) + scale(X4) + scale(X5) + scale(X6) + scale(X7) + scale(X8) + scale(X9) + scale(X10) + scale(X11) + scale(X12) + scale(X13) + scale(X14), data = merged1))
summary(lmer(scale(vi_restricted_behavior_score) ~ scale(`IQ_prs`) + scale(autism_prs) + scale(scz_prs) + scale(edu_prs) + (1|FID) + sex  + scale(age_at_eval_months) + scale(X1) + scale(X2) + scale(X3) + scale(X4) + scale(X5) + scale(X6) + scale(X7) + scale(X8) + scale(X9) + scale(X10) + scale(X11) + scale(X12) + scale(X13) + scale(X14), data = merged1))


##DCDQ
merged1 = merge(dcdq, merged, by = "IID")
summary(lmer(scale(final_score) ~ scale(`IQ_prs`) + scale(autism_prs) + scale(scz_prs) + scale(edu_prs) + anydenovo_constraint_sum + (1|family_sf_id) + sex  + scale(age_at_eval_months) + scale(X1) + scale(X2) + scale(X3) + scale(X4) + scale(X5) + scale(X6) + scale(X7) + scale(X8) + scale(X9) + scale(X10) + scale(X11) + scale(X12) + scale(X13) + scale(X14), data = merged1))



###IQ
merged1 = merge(merged, bgx_child, by = "IID")
merged1 = subset(merged1, asd == "TRUE") #4862

#Age_walking
summary(lmer(scale(smiled_age_mos) ~ scale(`IQ_prs`) + scale(autism_prs) + scale(scz_prs) + scale(edu_prs) + anydenovo_constraint_sum + (1|family_sf_id) + sex  + scale(age_at_eval_months) + scale(X1) + scale(X2) + scale(X3) + scale(X4) + scale(X5) + scale(X6) + scale(X7) + scale(X8) + scale(X9) + scale(X10) + scale(X11) + scale(X12) + scale(X13) + scale(X14), data = merged1))
summary(lmer(scale(walked_age_mos) ~ scale(`IQ_prs`) + scale(autism_prs) + scale(scz_prs) + scale(edu_prs) + anydenovo_constraint_sum + (1|family_sf_id) + sex  + scale(age_at_eval_months) + scale(X1) + scale(X2) + scale(X3) + scale(X4) + scale(X5) + scale(X6) + scale(X7) + scale(X8) + scale(X9) + scale(X10) + scale(X11) + scale(X12) + scale(X13) + scale(X14), data = merged1))


merged1$IQ_score = ifelse(merged1$cog_test_score == "24_below", 1, NA)
merged1$IQ_score = ifelse(merged1$cog_test_score == "25_39", 2, merged1$IQ_score)
merged1$IQ_score = ifelse(merged1$cog_test_score == "40_54", 3, merged1$IQ_score)
merged1$IQ_score = ifelse(merged1$cog_test_score == "55_69", 4, merged1$IQ_score)
merged1$IQ_score = ifelse(merged1$cog_test_score == "70_79", 5, merged1$IQ_score)
merged1$IQ_score = ifelse(merged1$cog_test_score == "80_89", 6, merged1$IQ_score)
merged1$IQ_score = ifelse(merged1$cog_test_score == "90_109", 7, merged1$IQ_score)
merged1$IQ_score = ifelse(merged1$cog_test_score == "110_119", 8, merged1$IQ_score)
merged1$IQ_score = ifelse(merged1$cog_test_score == "120_129", 9, merged1$IQ_score)
merged1$IQ_score = ifelse(merged1$cog_test_score == "130_above", 10, merged1$IQ_score)



summary(lmer(scale(IQ_score) ~ scale(`IQ_prs`) + scale(autism_prs) + scale(scz_prs) + scale(edu_prs) + anydenovo_constraint_sum + (1|FID.x) + sex  + scale(age_at_eval_months) + scale(X1) + scale(X2) + scale(X3) + scale(X4) + scale(X5) + scale(X6) + scale(X7) + scale(X8) + scale(X9) + scale(X10) + scale(X11) + scale(X12) + scale(X13) + scale(X14), data = merged1))

##Factor analysis

factor = fread("~/Autism_heterogeneity/Factor_data/SPARK_factorscores_6.txt")
merged1 = merge(merged, factor, by.x = "IID", by.y = "individual")
bgx2 = bgx_child[,c("IID", "sex", "age_at_eval_months")]

merged1 = merge(merged1, bgx2, by = "IID")

summary(lmer(scale(PA3) ~ scale(`IQ_prs`) + scale(autism_prs) + scale(scz_prs) + scale(edu_prs) + anydenovo_constraint_sum + (1|FID) + sex + scale(age_at_eval_months) + scale(X1) + scale(X2) + scale(X3) + scale(X4) + scale(X5) + scale(X6) + scale(X7) + scale(X8) + scale(X9) + scale(X10) + scale(X11) + scale(X12) + scale(X13) + scale(X14), data = merged1))


###Conditioning on IQ
merged2 = merge(scq, merged1, by = "IID")
merged2 = subset(merged2, asd.x == "TRUE")

summary(lmer(scale(final_score) ~ scale(`1.000000`) +  scale(IQ_score) + (1|FID) + sex.x  + scale(age_at_eval_months.x) + scale(X1) + scale(X2) + scale(X3) + scale(X4) + scale(X5) + scale(X6) + scale(X7) + scale(X8) + scale(X9) + scale(X10) + scale(X11) + scale(X12) + scale(X13) + scale(X14), data = merged2))


##RBS
merged2 = merge(rbs, merged1, by = "IID")

summary(lmer(scale(total_final_score) ~ scale(`0.100000`)+  scale(IQ_score) + (1|FID) + sex.x  + scale(age_at_eval_months.x) + scale(X1) + scale(X2) + scale(X3) + scale(X4) + scale(X5) + scale(X6) + scale(X7) + scale(X8) + scale(X9) + scale(X10) + scale(X11) + scale(X12) + scale(X13) + scale(X14), data = merged2))

summary(lmer(scale(i_stereotyped_behavior_score) ~ scale(`1.000000`) + (1|FID) + sex  + scale(age_at_eval_months) + scale(X1) + scale(X2) + scale(X3) + scale(X4) + scale(X5) + scale(X6) + scale(X7) + scale(X8) + scale(X9) + scale(X10) + scale(X11) + scale(X12) + scale(X13) + scale(X14), data = merged1))
summary(lmer(scale(ii_self_injurious_score) ~ scale(`1.000000`) + (1|FID) + sex  + scale(age_at_eval_months) + scale(X1) + scale(X2) + scale(X3) + scale(X4) + scale(X5) + scale(X6) + scale(X7) + scale(X8) + scale(X9) + scale(X10) + scale(X11) + scale(X12) + scale(X13) + scale(X14), data = merged1))
summary(lmer(scale(iii_compulsive_behavior_score) ~ scale(`1.000000`) + (1|FID) + sex  + scale(age_at_eval_months) + scale(X1) + scale(X2) + scale(X3) + scale(X4) + scale(X5) + scale(X6) + scale(X7) + scale(X8) + scale(X9) + scale(X10) + scale(X11) + scale(X12) + scale(X13) + scale(X14), data = merged1))
summary(lmer(scale(iv_ritualistic_behavior_score) ~ scale(`1.000000`) + (1|FID) + sex  + scale(age_at_eval_months) + scale(X1) + scale(X2) + scale(X3) + scale(X4) + scale(X5) + scale(X6) + scale(X7) + scale(X8) + scale(X9) + scale(X10) + scale(X11) + scale(X12) + scale(X13) + scale(X14), data = merged1))
summary(lmer(scale(v_sameness_behavior_score) ~ scale(`1.000000`) + (1|FID) + sex  + scale(age_at_eval_months) + scale(X1) + scale(X2) + scale(X3) + scale(X4) + scale(X5) + scale(X6) + scale(X7) + scale(X8) + scale(X9) + scale(X10) + scale(X11) + scale(X12) + scale(X13) + scale(X14), data = merged1))
summary(lmer(scale(vi_restricted_behavior_score) ~ scale(`1.000000`) + (1|FID) + sex  + scale(age_at_eval_months) + scale(X1) + scale(X2) + scale(X3) + scale(X4) + scale(X5) + scale(X6) + scale(X7) + scale(X8) + scale(X9) + scale(X10) + scale(X11) + scale(X12) + scale(X13) + scale(X14), data = merged1))


##DCDQ
merged2 = merge(dcdq, merged2, by = "IID")
summary(lmer(scale(final_score) ~ scale(`0.100000`) +  scale(IQ_score) + (1|FID) + sex.x  + scale(age_at_eval_months.x) + scale(X1) + scale(X2) + scale(X3) + scale(X4) + scale(X5) + scale(X6) + scale(X7) + scale(X8) + scale(X9) + scale(X10) + scale(X11) + scale(X12) + scale(X13) + scale(X14), data = merged2))






####

all.length.vs.width <- ggplot(merged, aes(x = total_final_score, y = final_score))

ID_info = fread("~/SPARK/Phenotypes/basic_medical_screening.txt")
ID = subset(ID_info, dev_id == "1")

merged_nonID = merged[!merged$IID %in% ID$IID,]

all.length.vs.width_nonID <- ggplot(merged_nonID, aes(x = ii_self_injurious_score, y = final_score))
all.length.vs.width_nonID +   geom_density_2d(aes(color = sex.y)) + labs(subtitle = "All species")


###SPARK_denovo

denovo_list = fread("SPARK_pilot_WES.txt")
autism_denovo_list = subset(denovo_list, Phenotype == "2")
autism_denovo_list = subset(autism_denovo_list, Father != 0)
autism_denovo_list = subset(autism_denovo_list, Mother != 0) #461 individuals with parental info for denovo

denovo_mutations = fread("SPARK_pilot_denovo_SNVs.txt")
metrics = fread("Gnomad_Shet_metrics.txt")

denovo_with_metrics = merge(denovo_mutations, metrics, by.x = "ANN.GENE", by.y = "gene")

contributing = subset(denovo_with_metrics, ANN.EFFECT != "synonymous_variant" & ANN.EFFECT != "missense_variant")
missense_contributing = subset(denovo_with_metrics, ANN.EFFECT == "missense_variant" & MPC > 2)

contributing = rbind(contributing, missense_contributing)

autism_denovo_list = merge(autism_denovo_list, contributing, by = "SPID", all.x = TRUE)

autism_denovo_list$pli9 = ifelse(autism_denovo_list$pLI > 0.9, 1, 0)
autism_denovo_list$pli995 = ifelse(autism_denovo_list$pLI > 0.9, 1, 0)
autism_denovo_list$loeuf_decile = ifelse(autism_denovo_list$oe_lof_upper_bin == 0, 1, 0)
autism_denovo_list$loeuf35 = ifelse(autism_denovo_list$oe_lof_upper < 0.35,1, 0)


denovo_annotated_sum = autism_denovo_list %>%
  group_by(SPID) %>% 
  transmute(pli9sum=sum(na.omit(pli9)), pli995sum=sum(na.omit(pli995)), pli=sum(na.omit(pLI)), shetsum=sum(na.omit(s_het)), loef_sum = sum(na.omit(1.99 - oe_lof_upper)), loefdecile_sum = sum(na.omit(loeuf_decile)))


##DCDQ
dcdq_merged = merge(denovo_annotated_sum, dcdq, by.x = "SPID", by.y = "IID")
summary(lmer(scale(final_score) ~ scale(loefdecile_sum) + (1|family_sf_id)+ sex + scale(age_at_eval_years), data = dcdq_merged ))



### IQ
IQ_Merged = merge(denovo_annotated_sum, bgx, by.x = "SPID", by.y = "IID")

IQ_Merged$IQ_score = ifelse(IQ_Merged$cog_test_score == "24_below", 1, NA)
IQ_Merged$IQ_score = ifelse(IQ_Merged$cog_test_score == "25_39", 2, IQ_Merged$IQ_score)
IQ_Merged$IQ_score = ifelse(IQ_Merged$cog_test_score == "40_54", 3, IQ_Merged$IQ_score)
IQ_Merged$IQ_score = ifelse(IQ_Merged$cog_test_score == "55_69", 4, IQ_Merged$IQ_score)
IQ_Merged$IQ_score = ifelse(IQ_Merged$cog_test_score == "70_79", 5, IQ_Merged$IQ_score)
IQ_Merged$IQ_score = ifelse(IQ_Merged$cog_test_score == "80_89", 6, IQ_Merged$IQ_score)
IQ_Merged$IQ_score = ifelse(IQ_Merged$cog_test_score == "90_109", 7, IQ_Merged$IQ_score)
IQ_Merged$IQ_score = ifelse(IQ_Merged$cog_test_score == "110_119", 8, IQ_Merged$IQ_score)
IQ_Merged$IQ_score = ifelse(IQ_Merged$cog_test_score == "120_129", 9, IQ_Merged$IQ_score)
IQ_Merged$IQ_score = ifelse(IQ_Merged$cog_test_score == "130_above", 10, IQ_Merged$IQ_score)


summary(lmer(scale(IQ_score) ~ scale(loefdecile_sum) + (1|family_sf_id)+ sex + scale(age_at_eval_years), data = IQ_Merged ))

###
rbs_merged = merge(denovo_annotated_sum, rbs, by.x = "SPID", by.y = "IID")
summary(lmer(scale(total_final_score) ~ scale(loefdecile_sum) + (1|family_id)+ sex + scale(age_at_eval_years), data = rbs_merged ))


###
scq_merged = merge(denovo_annotated_sum, scq, by.x = "SPID", by.y = "IID")
summary(lmer(scale(final_score) ~ scale(loefdecile_sum) + (1|family_id)+ sex + scale(age_at_eval_years), data = scq_merged ))


IQ_prs_merged = merge(IQ_Merged, prs_IQ, by.x = "SPID", by.y = "IID")
IQ_prs_merged = merge(IQ_prs_merged, PCS, by.x = "SPID", by.y = "IID")

summary(lmer(scale(IQ_score) ~ scale(loefdecile_sum) + scale(IQ_all) + (1|family_sf_id)+ sex + scale(age_at_eval_years) + + scale(X1) + scale(X2) + scale(X3) + scale(X4) + scale(X5) + scale(X6) + scale(X7) + scale(X8) + scale(X9) + scale(X10) + scale(X11) + scale(X12) + scale(X13) + scale(X14), data = IQ_prs_merged ))
