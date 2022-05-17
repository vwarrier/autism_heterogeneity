##Combine both SPARK datasets and generate PCs

./plink --bfile ~/SPARK/SPARK_preimputation/QC6output-updated_CEU  --bmerge ~/SPARK/SPARK_v2/QC6output-updated --make-bed --out spark1and2_forPC --threads 20


~/ABCD/ABCDgenotype/king -b ~/Autism_heterogeneity/spark1and2_forPC.bed --kinship --prefix ~/Autism_heterogeneity/spark1and2_forPC_kinship


library(GWASTools)
library(GENESIS)
library(SNPRelate)

## 5.2: Read files and prune SNPs for LD

snpgdsBED2GDS(bed.fn = "spark1and2_forPC.bed", 
              bim.fn = "spark1and2_forPC.bim", 
              fam.fn = "spark1and2_forPC.fam", 
              out.gdsfn = "spark1and2_forPC.gds")


gds <- snpgdsOpen("spark1and2_forPC.gds")
snpset <- snpgdsLDpruning(gds, method="corr", slide.max.bp=10e6, 
                          ld.threshold=sqrt(0.1), verbose=FALSE)
pruned <- unlist(snpset, use.names=FALSE)
length(pruned)

snpgdsClose(gds)

## 5.3 Read the kinship matrix, and run PC
KINGmat <- kingToMatrix(c("spark1and2_forPC_kinship.kin0","spark1and2_forPC_kinship.kin"), estimator = "Kinship")
KINGmat[1:5,1:5]

data1 <- GdsGenotypeReader(filename = "spark1and2_forPC.gds")
data2 <- GenotypeData(data1)
data2

mypcair <- pcair(data2, kinobj = KINGmat, divobj = KINGmat,
                 snp.include = pruned)


PC = data.frame(mypcair$vectors)
PC$Sample_name = row.names(PC)

write.table(PC, file = "SPARK_all_PC.txt", row.names = F, col.names = T, quote = F)


####
library(MASS)
library(data.table)
library(ggplot2)
library(dplyr)
library(lme4)
library(car)

Mother = fread("~/SPARK/Phenotypes/SPARK_mother.txt", header= T)
setnames(Mother, c("subject_sp_id", "family_id"), c("IID", "FID"))

Father = read.table("~/SPARK/Phenotypes/SPARK_father.txt", header= T)
setnames(Father, c("subject_sp_id", "family_id"), c("IID", "FID"))

Cases = fread("~/SPARK/Phenotypes/Autistic_proband.txt", header= T)
setnames(Cases, c("subject_sp_id", "family_id"), c("IID", "FID"))

siblings = fread("~/SPARK/Phenotypes/Nonautistic_sibling.txt", header= T)
setnames(siblings, c("subject_sp_id", "family_id"), c("IID", "FID"))

#Next, read the prs scores

prs_scz = fread("./PGS/SPARK_scz_finalscore.profile", header = TRUE)
prs_scz2 = fread("./PGS/SPARK_v2_scz_finalscore.profile")
prs_scz$category = "SPARK1"
prs_scz2$category = "SPARK2"
prs_scz = rbind(prs_scz, prs_scz2)
setnames(prs_scz, old = c("SCORE", "category"), new = c("scz_prs", "category"))
prs_scz = prs_scz[,c("IID", "category", "scz_prs")]


prs_edu = fread("./PGS/SPARK_edu_finalscore.profile", header = TRUE)
prs_edu2 = fread("./PGS/SPARK_v2_edu_finalscore.profile")
prs_edu = rbind(prs_edu, prs_edu2)
setnames(prs_edu, "SCORE", "edu_prs")
prs_edu = prs_edu[,c("IID", "edu_prs")]

prs_IQ = fread("./PGS/SPARK_IQ_finalscore.profile", header = TRUE)
prs_IQ2 = fread("./PGS/SPARK_v2_IQ_finalscore.profile")
prs_IQ = rbind(prs_IQ, prs_IQ2)
setnames(prs_IQ, "SCORE", "IQ_prs")
prs_IQ = prs_IQ[,c("IID", "IQ_prs")]

prs_autism = fread("./PGS/SPARK_autism_finalscore.profile", header = TRUE)
prs_autism2 = fread("./PGS/SPARK_v2_autism_finalscore.profile")
prs_autism = rbind(prs_autism, prs_autism2)
setnames(prs_autism, "SCORE", "autism_prs")
prs_autism = prs_autism[,c("IID", "autism_prs", "FID")]

prs_haircolour = fread("./PGS/SPARK_haircolour_finalscore.profile")
prs_haircolour2 = fread("./PGS/SPARK_v2_haircolour_finalscore.profile")
prs_haircolour = rbind(prs_haircolour, prs_haircolour2)
setnames(prs_haircolour, "SCORE", "haircolour_prs")
prs_haircolour = prs_haircolour[,c("IID", "haircolour_prs", "FID")]

prs_adhd = fread("./PGS/SPARKv1_ADHD_finalscore.profile")
prs_adhd2 = fread("./PGS/SPARKv2_ADHD_finalscore.profile")
prs_adhd = rbind(prs_adhd, prs_adhd2)
setnames(prs_adhd, "SCORE", "adhd_prs")
prs_adhd = prs_adhd[,c("IID", "adhd_prs", "FID")]



prs_merged = merge(prs_IQ, prs_edu, by = "IID")
prs_merged = merge(prs_merged, prs_scz, by = "IID")
prs_merged = merge(prs_merged, prs_autism, by = "IID")
prs_merged = merge(prs_merged, prs_haircolour, by = "IID")
prs_merged = merge(prs_merged, prs_adhd, by = "IID")


PCS = fread("~/Autism_heterogeneity/SPARK_all_PC.txt", header = TRUE)
setnames(PCS, "Sample_name", "IID")

merged = merge(prs_merged, PCS, by = "IID")

merged = unique(merged)

scq = fread("~/SPARK/Phenotypes/V5/scq.csv")
setnames(scq, "subject_sp_id", "IID")

rbs = fread("~/SPARK/Phenotypes/V5/rbsr.csv")
setnames(rbs, "subject_sp_id", "IID")
rbs = subset(rbs, asd == "TRUE")

dcdq = fread("~/SPARK/Phenotypes/V5/dcdq.csv")
setnames(dcdq, "subject_sp_id", "IID")

bgx_child = fread("~/SPARK/Phenotypes/V5/background_history_child.csv")
setnames(bgx_child, "subject_sp_id", "IID")

bgx_adult = fread("~/SPARK/Phenotypes/V5/background_history_adult.csv")
setnames(bgx_adult, "subject_sp_id", "IID")

vabs = fread("~/SPARK/Phenotypes/V5/vineland.csv")
setnames(vabs, "subject_sp_id", "IID")

factor = fread("~/Autism_heterogeneity/Factor_data/SPARK_factorscores_6.txt")
setnames(factor, 1, "IID")

###Create an omnibus_phenotype_file
vabs1 = vabs[,c("IID", "abc_standard")]
setnames(vabs1, 2, "vabs_score")
bgx1 = bgx_child[,c("IID", "sex", "cog_test_score", "age_at_eval_years")]
bgx_adult$cog_test_score = NA
bgx2 = bgx_adult[,c("IID", "sex", "cog_test_score", "age_at_eval_years")]
bgx3 = rbind(bgx1, bgx2)
setnames(vabs1, 2, "VABS_score")

rbs1 = rbs[,c("IID", "total_final_score")]
setnames(rbs1, 2, "rbs_score")

scq1 = scq[,c("IID", "final_score")]
setnames(scq1, 2, "scq_score")

dcdq1 = dcdq[,c("IID", "final_score")]
setnames(dcdq1, 2, "dcdq_score")

pheno_merged = merge(bgx3, vabs1, by = "IID", all = TRUE)
pheno_merged = merge(pheno_merged, rbs1, by = "IID", all = TRUE)
pheno_merged = merge(pheno_merged, scq1, by = "IID", all = TRUE)
pheno_merged = merge(pheno_merged, dcdq1, by = "IID", all = TRUE)
pheno_merged = merge(pheno_merged, factor, by = "IID", all = TRUE)

pheno_merged$IQ_score = ifelse(pheno_merged$cog_test_score == "24_below", 1, NA)
pheno_merged$IQ_score = ifelse(pheno_merged$cog_test_score == "25_39", 2, pheno_merged$IQ_score)
pheno_merged$IQ_score = ifelse(pheno_merged$cog_test_score == "40_54", 3, pheno_merged$IQ_score)
pheno_merged$IQ_score = ifelse(pheno_merged$cog_test_score == "55_69", 4, pheno_merged$IQ_score)
pheno_merged$IQ_score = ifelse(pheno_merged$cog_test_score == "70_79", 5, pheno_merged$IQ_score)
pheno_merged$IQ_score = ifelse(pheno_merged$cog_test_score == "80_89", 6, pheno_merged$IQ_score)
pheno_merged$IQ_score = ifelse(pheno_merged$cog_test_score == "90_109", 7, pheno_merged$IQ_score)
pheno_merged$IQ_score = ifelse(pheno_merged$cog_test_score == "110_119", 8, pheno_merged$IQ_score)
pheno_merged$IQ_score = ifelse(pheno_merged$cog_test_score == "120_129", 9, pheno_merged$IQ_score)
pheno_merged$IQ_score = ifelse(pheno_merged$cog_test_score == "130_above", 10, pheno_merged$IQ_score)

#write.table(pheno_merged[,-c("cog_test_score")], file = "~/Autism_heterogeneity/SPARK_allphenotypes.txt", row.names = F, col.names = T, quote = F)





##SCQ
scq = subset(scq, asd == "TRUE")
scq = scq[,c("IID", "final_score", "family_sf_id", "sex", "age_at_eval_months")]
merged1 = merge(scq, merged, by = "IID")

summary(lmer(scale(final_score) ~ scale(adhd_prs) + scale(`IQ_prs`) + scale(autism_prs) + scale(scz_prs) + scale(edu_prs) + scale(haircolour_prs) + (1|family_sf_id) + sex  + scale(age_at_eval_months) + scale(X1) + scale(X2) + scale(X3) + scale(X4) + scale(X5) + scale(X6) + scale(X7) + scale(X8) + scale(X9) + scale(X10), data = merged1))


##RBS
merged1 = merge(rbs, merged, by = "IID")

summary(lmer(scale(total_final_score) ~ scale(adhd_prs) + scale(IQ_prs) + scale(autism_prs) + scale(scz_prs) + scale(edu_prs)+ scale(haircolour_prs)  + (1|family_sf_id) + sex  + scale(age_at_eval_months) + scale(X1) + scale(X2) + scale(X3) + scale(X4) + scale(X5) + scale(X6) + scale(X7) + scale(X8) + scale(X9) + scale(X10), data = merged1))



##DCDQ
merged1 = merge(dcdq, merged, by = "IID")
summary(lmer(scale(final_score) ~ scale(adhd_prs) + scale(`IQ_prs`) + scale(autism_prs) + scale(scz_prs) + scale(edu_prs) + scale(haircolour_prs) +  (1|family_sf_id) + sex  + scale(age_at_eval_months) + scale(X1) + scale(X2) + scale(X3) + scale(X4) + scale(X5) + scale(X6) + scale(X7) + scale(X8) + scale(X9) + scale(X10), data = merged1))



###IQ
merged1 = merge(merged, bgx_child, by = "IID")
merged1_asd = subset(merged1, asd == "TRUE") #4862
merged1_asd$walked_age_mos = ifelse(merged1_asd$walked_age_mos == "888", NA, merged1_asd$walked_age_mos)
merged1_asd$smiled_age_mos = ifelse(merged1_asd$smiled_age_mos == "888", NA, merged1_asd$smiled_age_mos)
merged1_asd$used_words_age_mos = ifelse(merged1_asd$used_words_age_mos == "888", NA, merged1_asd$used_words_age_mos)
merged1_asd$crawled_age_mos = ifelse(merged1_asd$crawled_age_mos == "888", NA, merged1_asd$crawled_age_mos)
merged1_asd$sat_wo_support_age_mos = ifelse(merged1_asd$sat_wo_support_age_mos == "888", NA, merged1_asd$sat_wo_support_age_mos)
merged1_asd$combined_words_age_mos = ifelse(merged1_asd$combined_words_age_mos == "888", NA, merged1_asd$combined_words_age_mos)


#Age_walking
summary(glm.nb(walked_age_mos ~ scale(adhd_prs) + scale(`IQ_prs`) + scale(autism_prs) + scale(scz_prs) + scale(edu_prs) + sex  + scale(age_at_eval_months) + scale(X1) + scale(X2) + scale(X3) + scale(X4) + scale(X5) + scale(X6) + scale(X7) + scale(X8) + scale(X9) + scale(X10) + scale(X11) + scale(X12) + scale(X13) + scale(X14), data = merged1_asd))


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



summary(lmer(scale(IQ_score) ~ scale(adhd_prs) + scale(`IQ_prs`) + scale(autism_prs) + scale(scz_prs) + scale(edu_prs) + scale(haircolour_prs) + (1|family_sf_id) + sex  + scale(age_at_eval_months) + scale(X1) + scale(X2) + scale(X3) + scale(X4) + scale(X5) + scale(X6) + scale(X7) + scale(X8) + scale(X9) + scale(X10), data = merged1))

merged1$walked_age_mos = ifelse(merged1$walked_age_mos == "888", NA, merged1$walked_age_mos)
merged1$used_words_age_mos = ifelse(merged1$used_words_age_mos == "888", NA, merged1$used_words_age_mos)

summary(glm.nb(walked_age_mos ~ scale(IQ_score) + scale(`IQ_prs`) + scale(autism_prs) + scale(scz_prs) + scale(edu_prs) + sex  + scale(age_at_eval_months) + scale(X1) + scale(X2) + scale(X3) + scale(X4) + scale(X5) + scale(X6) + scale(X7) + scale(X8) + scale(X9) + scale(X10) + scale(X11) + scale(X12) + scale(X13) + scale(X14), data = merged1))
summary(glm.nb(used_words_age_mos ~ scale(IQ_score) + scale(`IQ_prs`) + scale(autism_prs) + scale(scz_prs) + scale(edu_prs) + sex  + scale(age_at_eval_months) + scale(X1) + scale(X2) + scale(X3) + scale(X4) + scale(X5) + scale(X6) + scale(X7) + scale(X8) + scale(X9) + scale(X10) + scale(X11) + scale(X12) + scale(X13) + scale(X14), data = merged1))


###in sibling
bgx_sibling = fread("~/SPARK/Phenotypes/V5/background_history_sibling.csv")
setnames(bgx_sibling, "subject_sp_id", "IID")
merged1 = merge(merged, bgx_sibling, by = "IID")
merged1_nonasd = subset(merged1, asd == "FALSE") #1745
merged1_nonasd$walked_age_mos = ifelse(merged1_nonasd$walked_age_mos == "888", NA, merged1_nonasd$walked_age_mos)
merged1_nonasd$used_words_age_mos = ifelse(merged1_nonasd$used_words_age_mos == "888", NA, merged1_nonasd$used_words_age_mos)


##Factor analysis

factor = fread("~/Autism_heterogeneity/Factor_data/SPARK_factorscores_6.txt")
merged1 = merge(merged, factor, by.x = "IID", by.y = "individual")
bgx2 = bgx_child[,c("IID", "sex", "age_at_eval_months")]

merged1 = merge(merged1, bgx2, by = "IID")

summary(lmer(scale(PA6) ~scale(adhd_prs) + scale(`IQ_prs`) + scale(autism_prs) + scale(scz_prs) + scale(edu_prs) + scale(haircolour_prs) + (1|FID.y) + sex + scale(age_at_eval_months) + scale(X1) + scale(X2) + scale(X3) + scale(X4) + scale(X5) + scale(X6) + scale(X7) + scale(X8) + scale(X9) + scale(X10), data = merged1))


###VABS
merged1 = merge(vabs, merged, by = "IID")
summary(lmer(scale(abc_standard) ~ scale(`adhd_prs`) + scale(`IQ_prs`) + scale(autism_prs) + scale(scz_prs) + scale(edu_prs)+ scale(haircolour_prs) +  (1|family_sf_id) + sex  + scale(age_at_eval_months) + scale(X1) + scale(X2) + scale(X3) + scale(X4) + scale(X5) + scale(X6) + scale(X7) + scale(X8) + scale(X9) + scale(X10), data = merged1))



###Run with IQ_score

merged2 = merge(merged, pheno_merged, by ="IID")
#rm individuals without IQ score
merged2 = merged2[!is.na(merged2$IQ_score),]

summary(lmer(scale(VABS_score) ~ scale(`IQ_prs`) + scale(autism_prs) + scale(scz_prs) + scale(edu_prs) + (1|FID) + sex  +  IQ_score + scale(age_at_eval_years) + scale(X1) + scale(X2) + scale(X3) + scale(X4) + scale(X5) + scale(X6) + scale(X7) + scale(X8) + scale(X9) + scale(X10), data = merged2))
summary(lmer(scale(PA1) ~ scale(`IQ_prs`) + scale(autism_prs) + scale(scz_prs) + scale(edu_prs) + (1|FID) + sex  +  IQ_score + scale(age_at_eval_years) + scale(X1) + scale(X2) + scale(X3) + scale(X4) + scale(X5) + scale(X6) + scale(X7) + scale(X8) + scale(X9) + scale(X10), data = merged2))
summary(lmer(scale(PA2) ~ scale(`IQ_prs`) + scale(autism_prs) + scale(scz_prs) + scale(edu_prs) + (1|FID) + sex  +  IQ_score + scale(age_at_eval_years) + scale(X1) + scale(X2) + scale(X3) + scale(X4) + scale(X5) + scale(X6) + scale(X7) + scale(X8) + scale(X9) + scale(X10), data = merged2))
summary(lmer(scale(PA3) ~ scale(`IQ_prs`) + scale(autism_prs) + scale(scz_prs) + scale(edu_prs) + (1|FID) + sex  +  IQ_score + scale(age_at_eval_years) + scale(X1) + scale(X2) + scale(X3) + scale(X4) + scale(X5) + scale(X6) + scale(X7) + scale(X8) + scale(X9) + scale(X10), data = merged2))
summary(lmer(scale(PA4) ~ scale(`IQ_prs`) + scale(autism_prs) + scale(scz_prs) + scale(edu_prs) + (1|FID) + sex  +  IQ_score + scale(age_at_eval_years) + scale(X1) + scale(X2) + scale(X3) + scale(X4) + scale(X5) + scale(X6) + scale(X7) + scale(X8) + scale(X9) + scale(X10), data = merged2))
summary(lmer(scale(PA5) ~ scale(`IQ_prs`) + scale(autism_prs) + scale(scz_prs) + scale(edu_prs) + (1|FID) + sex  +  IQ_score + scale(age_at_eval_years) + scale(X1) + scale(X2) + scale(X3) + scale(X4) + scale(X5) + scale(X6) + scale(X7) + scale(X8) + scale(X9) + scale(X10), data = merged2))
summary(lmer(scale(scq_score) ~ scale(`IQ_prs`) + scale(autism_prs) + scale(scz_prs) + scale(edu_prs) + (1|FID) + sex  +  IQ_score + scale(age_at_eval_years) + scale(X1) + scale(X2) + scale(X3) + scale(X4) + scale(X5) + scale(X6) + scale(X7) + scale(X8) + scale(X9) + scale(X10), data = merged2))
summary(lmer(scale(rbs_score) ~ scale(`IQ_prs`) + scale(autism_prs) + scale(scz_prs) + scale(edu_prs) + (1|FID) + sex  +  IQ_score + scale(age_at_eval_years) + scale(X1) + scale(X2) + scale(X3) + scale(X4) + scale(X5) + scale(X6) + scale(X7) + scale(X8) + scale(X9) + scale(X10), data = merged2))






###medical_history
med = fread("~/SPARK/Phenotypes/V5/basic_medical_screening.csv")
med[is.na(med)] <- 0 
setnames(med, 1, "IID")
med_asd = subset(med,asd == "TRUE")
merged1 = merge(med_asd, merged, by = "IID")

merged1$dev_diagnosis = merged1$dev_id + merged1$dev_lang_dis + merged1$dev_ld + merged1$dev_motor + merged1$dev_mutism + merged1$dev_soc_prag + merged1$dev_speech

summary(glm(dev_diagnosis ~ scale(`IQ_prs`) + scale(autism_prs) + scale(scz_prs) + scale(edu_prs)  + sex  + scale(age_at_eval_months) + scale(X1) + scale(X2) + scale(X3) + scale(X4) + scale(X5) + scale(X6) + scale(X7) + scale(X8) + scale(X9) + scale(X10), data = merged1, family = "quasipoisson"))

#LOO
merged1$dev_diagnosis_2 = merged1$dev_id + merged1$dev_lang_dis + merged1$dev_ld + merged1$dev_motor + merged1$dev_mutism + merged1$dev_soc_prag

summary(glm(dev_diagnosis_2 ~ scale(`IQ_prs`) + scale(autism_prs) + scale(scz_prs) + scale(edu_prs)  + sex  + scale(age_at_eval_months) + scale(X1) + scale(X2) + scale(X3) + scale(X4) + scale(X5) + scale(X6) + scale(X7) + scale(X8) + scale(X9) + scale(X10), data = merged1, family = "quasipoisson"))





merged1$dev_diagnosis_2 = ifelse(merged1$asd == "FALSE", "a_sibling", merged1$dev_diagnosis)

merged1$dev_diagnosis_2 = ifelse(merged1$asd == "TRUE" & merged1$dev_diagnosis == 0, "zero", merged1$dev_diagnosis_2)
merged1$dev_diagnosis_2 = ifelse(merged1$asd == "TRUE" & merged1$dev_diagnosis > 0 & merged1$dev_diagnosis < 3, "one_two", merged1$dev_diagnosis_2)
merged1$dev_diagnosis_2 = ifelse(merged1$asd == "TRUE" & merged1$dev_diagnosis > 2 & merged1$dev_diagnosis < 5, "three_four", merged1$dev_diagnosis_2)
merged1$dev_diagnosis_2 = ifelse(merged1$asd == "TRUE" & merged1$dev_diagnosis > 4 , "five_plus", merged1$dev_diagnosis_2)

summary(lm(scale(autism_prs) ~ as.factor(dev_diagnosis_2) + sex  + scale(age_at_eval_months) + scale(X1) + scale(X2) + scale(X3) + scale(X4) + scale(X5) + scale(X6) + scale(X7) + scale(X8) + scale(X9) + scale(X10), data = merged1, family = "quasipoisson"))


merged1$dev_diagnosis_3 = ifelse(merged1$asd == "TRUE" & merged1$dev_diagnosis == 0, "zero", merged1$dev_diagnosis)
merged1$dev_diagnosis_3 = ifelse(merged1$asd == "TRUE" & merged1$dev_diagnosis > 0 & merged1$dev_diagnosis < 3, "one_two", merged1$dev_diagnosis_3)
merged1$dev_diagnosis_3 = ifelse(merged1$asd == "TRUE" & merged1$dev_diagnosis > 2 & merged1$dev_diagnosis < 5, "three_four", merged1$dev_diagnosis_3)
merged1$dev_diagnosis_3 = ifelse(merged1$asd == "TRUE" & merged1$dev_diagnosis > 4 , "five_plus", merged1$dev_diagnosis_3)

merged1$dev_diagnosis_3 = ifelse(merged1$asd == "FALSE" & merged1$dev_diagnosis == 0, "sibling_zero", merged1$dev_diagnosis_3)
merged1$dev_diagnosis_3 = ifelse(merged1$asd == "FALSE" & merged1$dev_diagnosis > 0 & merged1$dev_diagnosis < 3, "sibling_one_two", merged1$dev_diagnosis_3)
merged1$dev_diagnosis_3 = ifelse(merged1$asd == "FALSE" & merged1$dev_diagnosis > 2, "a_sibling_three_four", merged1$dev_diagnosis_3)


summary(lm(scale(autism_prs) ~ as.factor(dev_diagnosis_3) + sex  + scale(age_at_eval_months) + scale(X1) + scale(X2) + scale(X3) + scale(X4) + scale(X5) + scale(X6) + scale(X7) + scale(X8) + scale(X9) + scale(X10), data = merged1, family = "quasipoisson"))


###Mean in autistic individuals without any dev disorders


####

all.length.vs.width <- ggplot(merged, aes(x = total_final_score, y = final_score))

ID_info = fread("~/SPARK/Phenotypes/basic_medical_screening.txt")
ID = subset(ID_info, dev_id == "1")

merged_nonID = merged[!merged$IID %in% ID$IID,]

all.length.vs.width_nonID <- ggplot(merged_nonID, aes(x = ii_self_injurious_score, y = final_score))
all.length.vs.width_nonID +   geom_density_2d(aes(color = sex.y)) + labs(subtitle = "All species")



