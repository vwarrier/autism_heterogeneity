## LEAP v2 scripts

##Install packages

library(gplot2)
library(Hmisc)
library(ggplot2)
library(data.table)

##Background functions
# ++++++++++++++++++++++++++++
# flattenCorrMatrix
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

##Read all variable etc

prs_scz = fread("~/ALSPAC/PRSice2results/LEAP_v2_daner_natgen_pgc_eas_eur.all.score")

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

prs_edu = fread("~/ALSPAC/PRSice2results/LEAP_v2_edu2018leeforPRSICE.all.score", header = TRUE)
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

prs_IQ = fread("~/ALSPAC/PRSice2results/LEAP_v2_savageCP2018forprsice.all.score", header = TRUE)
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
prs_merged$IID = as.numeric(prs_merged$IID)

PCs = read.table("~/LEAP_PARIS/LEAP_v2/LEAP_ancestry.txt", header = T)

merged = merge(prs_merged, PCs, by = "IID")


phenotype1 = fread("~/LEAP_PARIS/LEAP_phenotype/LEAP_t1_Core_clinical_variables_withvalues.txt")
phenotype1$IID = as.numeric(phenotype1$IID)

merged_all = merge(merged, phenotype1, by = "IID") #530 observations

##keep participants
merged_all_2 = subset(merged_all, Relation == "participant") #509


##keep autistic indiviudals
merged_all_3 = subset(merged_all_2, t1_diagnosis == 2) #294

merged_all_3 = subset(merged_all_3, hdbscan == 9) #262


merged_all_3[merged_all_3==999]<-NA
merged_all_3[merged_all_3==777]<-NA

males = subset(merged_all_3, sex == "1")

non_id = subset(merged_all_3, t1_fsiq > 70 )

##Run regressions - verbal IQ

dim(as.table(na.omit(merged_all_3$t1_viq)))


summary(lm(scale(t1_viq) ~ scale(IQ_all) + as.character(sex) + scale(t1_age) + as.character(t1_iqtype) + scale(PC1) + scale(PC2) + scale(PC3) + scale(PC4) + scale(PC5) +
               scale(PC6) + scale(PC7) + scale(PC8) + scale(PC9) + scale(PC10) + scale(PC11) + scale(PC12) + scale(PC13) + scale(PC14) + scale(PC15), data = merged_all_3))


summary(lm(scale(t1_viq) ~ scale(edu_all) + as.character(sex) + scale(t1_age) + as.character(t1_iqtype) + scale(PC1) + scale(PC2) + scale(PC3) + scale(PC4) + scale(PC5) +
             scale(PC6) + scale(PC7) + scale(PC8) + scale(PC9) + scale(PC10) + scale(PC11) + scale(PC12) + scale(PC13) + scale(PC14) + scale(PC15), data = merged_all_3))


summary(lm(scale(t1_viq) ~ scale(scz_1) + as.character(sex) + scale(t1_age) + as.character(t1_iqtype) + scale(PC1) + scale(PC2) + scale(PC3) + scale(PC4) + scale(PC5) +
             scale(PC6) + scale(PC7) + scale(PC8) + scale(PC9) + scale(PC10) + scale(PC11) + scale(PC12) + scale(PC13) + scale(PC14) + scale(PC15), data = merged_all_3))



##Run regressions - FSIQ


dim(as.table(na.omit(merged_all_3$t1_fsiq)))


summary(lm(scale(t1_fsiq) ~ scale(IQ_all) + as.character(sex) + scale(t1_age) + as.character(t1_iqtype) + scale(PC1) + scale(PC2) + scale(PC3) + scale(PC4) + scale(PC5) +
             scale(PC6) + scale(PC7) + scale(PC8) + scale(PC9) + scale(PC10) + scale(PC11) + scale(PC12) + scale(PC13) + scale(PC14) + scale(PC15), data = merged_all_3))


summary(lm(scale(t1_fsiq) ~ scale(edu_all) + as.character(sex) + scale(t1_age) + as.character(t1_iqtype) + scale(PC1) + scale(PC2) + scale(PC3) + scale(PC4) + scale(PC5) +
             scale(PC6) + scale(PC7) + scale(PC8) + scale(PC9) + scale(PC10) + scale(PC11) + scale(PC12) + scale(PC13) + scale(PC14) + scale(PC15), data = merged_all_3))


summary(lm(scale(t1_fsiq) ~ scale(scz_1) + as.character(sex) + scale(t1_age) + as.character(t1_iqtype) + scale(PC1) + scale(PC2) + scale(PC3) + scale(PC4) + scale(PC5) +
             scale(PC6) + scale(PC7) + scale(PC8) + scale(PC9) + scale(PC10) + scale(PC11) + scale(PC12) + scale(PC13) + scale(PC14) + scale(PC15), data = merged_all_3))



##Run regressions - NVIQ


dim(as.table(na.omit(merged_all_3$t1_piq)))


summary(lm(scale(t1_piq) ~ scale(IQ_all) + as.character(sex) + scale(t1_age) + as.character(t1_iqtype) + scale(PC1) + scale(PC2) + scale(PC3) + scale(PC4) + scale(PC5) +
             scale(PC6) + scale(PC7) + scale(PC8) + scale(PC9) + scale(PC10) + scale(PC11) + scale(PC12) + scale(PC13) + scale(PC14) + scale(PC15), data = merged_all_3))


summary(lm(scale(t1_piq) ~ scale(edu_all) + as.character(sex) + scale(t1_age) + as.character(t1_iqtype) + scale(PC1) + scale(PC2) + scale(PC3) + scale(PC4) + scale(PC5) +
             scale(PC6) + scale(PC7) + scale(PC8) + scale(PC9) + scale(PC10) + scale(PC11) + scale(PC12) + scale(PC13) + scale(PC14) + scale(PC15), data = merged_all_3))


summary(lm(scale(t1_piq) ~ scale(scz_1) + as.character(sex) + scale(t1_age) + as.character(t1_iqtype) + scale(PC1) + scale(PC2) + scale(PC3) + scale(PC4) + scale(PC5) +
             scale(PC6) + scale(PC7) + scale(PC8) + scale(PC9) + scale(PC10) + scale(PC11) + scale(PC12) + scale(PC13) + scale(PC14) + scale(PC15), data = merged_all_3))



### Run regressions: ADOS


dim(as.table(na.omit(merged_all_3$t1_adi_communication_total)))

summary(lm(scale(t1_adi_communication_total) ~ scale(IQ_all) + as.character(sex) + scale(t1_age)  + scale(PC1) + scale(PC2) + scale(PC3) + scale(PC4) + scale(PC5) +
             scale(PC6) + scale(PC7) + scale(PC8) + scale(PC9) + scale(PC10) + scale(PC11) + scale(PC12) + scale(PC13) + scale(PC14) + scale(PC15), data = merged_all_3))


summary(lm(scale(t1_adi_communication_total) ~ scale(edu_all) + as.character(sex) + scale(t1_age)  + scale(PC1) + scale(PC2) + scale(PC3) + scale(PC4) + scale(PC5) +
             scale(PC6) + scale(PC7) + scale(PC8) + scale(PC9) + scale(PC10) + scale(PC11) + scale(PC12) + scale(PC13) + scale(PC14) + scale(PC15), data = merged_all_3))


summary(lm(scale(t1_adi_communication_total) ~ scale(scz_1) + as.character(sex) + scale(t1_age) + scale(PC1) + scale(PC2) + scale(PC3) + scale(PC4) + scale(PC5) +
             scale(PC6) + scale(PC7) + scale(PC8) + scale(PC9) + scale(PC10) + scale(PC11) + scale(PC12) + scale(PC13) + scale(PC14) + scale(PC15), data = merged_all_3))





## Run regressions ADI rrb
dim(as.table(na.omit(merged_all_3$t1_adi_rrb_total)))

summary(lm(scale(t1_adi_rrb_total) ~ scale(IQ_all) + as.character(sex) + scale(t1_age)  + scale(PC1) + scale(PC2) + scale(PC3) + scale(PC4) + scale(PC5) +
             scale(PC6) + scale(PC7) + scale(PC8) + scale(PC9) + scale(PC10) + scale(PC11) + scale(PC12) + scale(PC13) + scale(PC14) + scale(PC15), data = merged_all_3))


summary(lm(scale(t1_adi_rrb_total) ~ scale(edu_all) + as.character(sex) + scale(t1_age)  + scale(PC1) + scale(PC2) + scale(PC3) + scale(PC4) + scale(PC5) +
             scale(PC6) + scale(PC7) + scale(PC8) + scale(PC9) + scale(PC10) + scale(PC11) + scale(PC12) + scale(PC13) + scale(PC14) + scale(PC15), data = merged_all_3))


summary(lm(scale(t1_adi_rrb_total) ~ scale(scz_1) + as.character(sex) + scale(t1_age) + scale(PC1) + scale(PC2) + scale(PC3) + scale(PC4) + scale(PC5) +
             scale(PC6) + scale(PC7) + scale(PC8) + scale(PC9) + scale(PC10) + scale(PC11) + scale(PC12) + scale(PC13) + scale(PC14) + scale(PC15), data = merged_all_3))



## Run regressions ADI social
dim(as.table(na.omit(merged_all_3$t1_adi_social_total)))

summary(lm(scale(t1_adi_social_total) ~ scale(IQ_all) + as.character(sex) + scale(t1_age)  + scale(PC1) + scale(PC2) + scale(PC3) + scale(PC4) + scale(PC5) +
             scale(PC6) + scale(PC7) + scale(PC8) + scale(PC9) + scale(PC10) + scale(PC11) + scale(PC12) + scale(PC13) + scale(PC14) + scale(PC15), data = merged_all_3))


summary(lm(scale(t1_adi_social_total) ~ scale(edu_all) + as.character(sex) + scale(t1_age)  + scale(PC1) + scale(PC2) + scale(PC3) + scale(PC4) + scale(PC5) +
             scale(PC6) + scale(PC7) + scale(PC8) + scale(PC9) + scale(PC10) + scale(PC11) + scale(PC12) + scale(PC13) + scale(PC14) + scale(PC15), data = merged_all_3))


summary(lm(scale(t1_adi_social_total) ~ scale(scz_1) + as.character(sex) + scale(t1_age) + scale(PC1) + scale(PC2) + scale(PC3) + scale(PC4) + scale(PC5) +
             scale(PC6) + scale(PC7) + scale(PC8) + scale(PC9) + scale(PC10) + scale(PC11) + scale(PC12) + scale(PC13) + scale(PC14) + scale(PC15), data = merged_all_3))


## Run regressions on VABS
dim(as.table(na.omit(merged_all_3$t1_vabsabcabc_standard)))


summary(lm(scale(t1_vabsabcabc_standard) ~ scale(IQ_all) + as.character(sex)  + scale(t1_age) + scale(PC1) + scale(PC2) + scale(PC3) + scale(PC4) + scale(PC5) +
             scale(PC6) + scale(PC7) + scale(PC8) + scale(PC9) + scale(PC10) + scale(PC11) + scale(PC12) + scale(PC13) + scale(PC14) + scale(PC15), data = merged_all_3))


summary(lm(scale(t1_vabsabcabc_standard) ~ scale(edu_all) + as.character(sex) + scale(t1_age) + scale(PC1) + scale(PC2) + scale(PC3) + scale(PC4) + scale(PC5) +
             scale(PC6) + scale(PC7) + scale(PC8) + scale(PC9) + scale(PC10) + scale(PC11) + scale(PC12) + scale(PC13) + scale(PC14) + scale(PC15), data = merged_all_3))


summary(lm(scale(t1_vabsabcabc_standard) ~ scale(scz_1) + as.character(sex) + scale(t1_age) + scale(PC1) + scale(PC2) + scale(PC3) + scale(PC4) + scale(PC5) +
             scale(PC6) + scale(PC7) + scale(PC8) + scale(PC9) + scale(PC10) + scale(PC11) + scale(PC12) + scale(PC13) + scale(PC14) + scale(PC15), data = merged_all_3))

## Run regression on SRS

dim(as.table(na.omit(merged_all_3$t1_srs_rawscore)))


summary(lm(scale(t1_srs_rawscore) ~ scale(IQ_all) + as.character(sex)  + scale(t1_age) + scale(PC1) + scale(PC2) + scale(PC3) + scale(PC4) + scale(PC5) +
             scale(PC6) + scale(PC7) + scale(PC8) + scale(PC9) + scale(PC10) + scale(PC11) + scale(PC12) + scale(PC13) + scale(PC14) + scale(PC15), data = merged_all_3))


summary(lm(scale(t1_srs_rawscore) ~ scale(edu_all) + as.character(sex) + scale(t1_age) + scale(PC1) + scale(PC2) + scale(PC3) + scale(PC4) + scale(PC5) +
             scale(PC6) + scale(PC7) + scale(PC8) + scale(PC9) + scale(PC10) + scale(PC11) + scale(PC12) + scale(PC13) + scale(PC14) + scale(PC15), data = merged_all_3))


summary(lm(scale(t1_srs_rawscore) ~ scale(scz_1) + as.character(sex) + scale(t1_age) + scale(PC1) + scale(PC2) + scale(PC3) + scale(PC4) + scale(PC5) +
             scale(PC6) + scale(PC7) + scale(PC8) + scale(PC9) + scale(PC10) + scale(PC11) + scale(PC12) + scale(PC13) + scale(PC14) + scale(PC15), data = merged_all_3))


## Run regression on RBS-R

dim(as.table(na.omit(merged_all_3$t1_rbs_total)))


summary(lm(scale(t1_rbs_total) ~ scale(IQ_all) + as.character(sex)  + scale(t1_age) + scale(PC1) + scale(PC2) + scale(PC3) + scale(PC4) + scale(PC5) +
             scale(PC6) + scale(PC7) + scale(PC8) + scale(PC9) + scale(PC10) + scale(PC11) + scale(PC12) + scale(PC13) + scale(PC14) + scale(PC15), data = merged_all_3))


summary(lm(scale(t1_rbs_total) ~ scale(edu_all) + as.character(sex) + scale(PC1) + scale(PC2) + scale(PC3) + scale(PC4) + scale(PC5) +
             scale(PC6) + scale(PC7) + scale(PC8) + scale(PC9) + scale(PC10) + scale(PC11) + scale(PC12) + scale(PC13) + scale(PC14) + scale(PC15), data = merged_all_3))


summary(lm(scale(t1_rbs_total) ~ scale(scz_1) + as.character(sex) + scale(t1_age) + scale(PC1) + scale(PC2) + scale(PC3) + scale(PC4) + scale(PC5) +
             scale(PC6) + scale(PC7) + scale(PC8) + scale(PC9) + scale(PC10) + scale(PC11) + scale(PC12) + scale(PC13) + scale(PC14) + scale(PC15), data = merged_all_3))



## Run regression on ADOS_social_affect

dim(as.table(na.omit(merged_all_3$t1_sa_css_all)))


summary(lm(scale(t1_sa_css_all) ~ scale(IQ_all) + as.character(sex)  + scale(t1_age) + scale(PC1) + scale(PC2) + scale(PC3) + scale(PC4) + scale(PC5) +
             scale(PC6) + scale(PC7) + scale(PC8) + scale(PC9) + scale(PC10) + scale(PC11) + scale(PC12) + scale(PC13) + scale(PC14) + scale(PC15), data = merged_all_3))


summary(lm(scale(t1_sa_css_all) ~ scale(edu_all) + as.character(sex) + scale(PC1) + scale(PC2) + scale(PC3) + scale(PC4) + scale(PC5) +
             scale(PC6) + scale(PC7) + scale(PC8) + scale(PC9) + scale(PC10) + scale(PC11) + scale(PC12) + scale(PC13) + scale(PC14) + scale(PC15), data = merged_all_3))


summary(lm(scale(t1_sa_css_all) ~ scale(scz_1) + as.character(sex) + scale(t1_age) + scale(PC1) + scale(PC2) + scale(PC3) + scale(PC4) + scale(PC5) +
             scale(PC6) + scale(PC7) + scale(PC8) + scale(PC9) + scale(PC10) + scale(PC11) + scale(PC12) + scale(PC13) + scale(PC14) + scale(PC15), data = merged_all_3))




## Run regression on ADOS_rrb_affect

dim(as.table(na.omit(merged_all_3$t1_rrb_css_all)))


summary(lm(scale(t1_rrb_css_all) ~ scale(IQ_all) + as.character(sex)  + scale(t1_age) + scale(PC1) + scale(PC2) + scale(PC3) + scale(PC4) + scale(PC5) +
             scale(PC6) + scale(PC7) + scale(PC8) + scale(PC9) + scale(PC10) + scale(PC11) + scale(PC12) + scale(PC13) + scale(PC14) + scale(PC15), data = merged_all_3))


summary(lm(scale(t1_rrb_css_all) ~ scale(edu_all) + as.character(sex) + scale(PC1) + scale(PC2) + scale(PC3) + scale(PC4) + scale(PC5) +
             scale(PC6) + scale(PC7) + scale(PC8) + scale(PC9) + scale(PC10) + scale(PC11) + scale(PC12) + scale(PC13) + scale(PC14) + scale(PC15), data = merged_all_3))


summary(lm(scale(t1_rrb_css_all) ~ scale(scz_1) + as.character(sex) + scale(t1_age) + scale(PC1) + scale(PC2) + scale(PC3) + scale(PC4) + scale(PC5) +
             scale(PC6) + scale(PC7) + scale(PC8) + scale(PC9) + scale(PC10) + scale(PC11) + scale(PC12) + scale(PC13) + scale(PC14) + scale(PC15), data = merged_all_3))






##############################Run on males###########################################################
#####################################################################################################
#####################################################################################################

dim(as.table(na.omit(males$t1_viq)))


summary(lm(scale(t1_viq) ~ scale(IQ_all) + scale(t1_age) + as.character(t1_iqtype) + scale(PC1) + scale(PC2) + scale(PC3) + scale(PC4) + scale(PC5) +
             scale(PC6) + scale(PC7) + scale(PC8) + scale(PC9) + scale(PC10) + scale(PC11) + scale(PC12) + scale(PC13) + scale(PC14) + scale(PC15), data = males))


summary(lm(scale(t1_viq) ~ scale(edu_all) + scale(t1_age) + as.character(t1_iqtype) + scale(PC1) + scale(PC2) + scale(PC3) + scale(PC4) + scale(PC5) +
             scale(PC6) + scale(PC7) + scale(PC8) + scale(PC9) + scale(PC10) + scale(PC11) + scale(PC12) + scale(PC13) + scale(PC14) + scale(PC15), data = males))


summary(lm(scale(t1_viq) ~ scale(scz_1) + scale(t1_age) + as.character(t1_iqtype) + scale(PC1) + scale(PC2) + scale(PC3) + scale(PC4) + scale(PC5) +
             scale(PC6) + scale(PC7) + scale(PC8) + scale(PC9) + scale(PC10) + scale(PC11) + scale(PC12) + scale(PC13) + scale(PC14) + scale(PC15), data = males))



##Run regressions - FSIQ


dim(as.table(na.omit(males$t1_fsiq)))


summary(lm(scale(t1_fsiq) ~ scale(IQ_all) + scale(t1_age) + as.character(t1_iqtype) + scale(PC1) + scale(PC2) + scale(PC3) + scale(PC4) + scale(PC5) +
             scale(PC6) + scale(PC7) + scale(PC8) + scale(PC9) + scale(PC10) + scale(PC11) + scale(PC12) + scale(PC13) + scale(PC14) + scale(PC15), data = males))


summary(lm(scale(t1_fsiq) ~ scale(edu_all) + scale(t1_age) + as.character(t1_iqtype) + scale(PC1) + scale(PC2) + scale(PC3) + scale(PC4) + scale(PC5) +
             scale(PC6) + scale(PC7) + scale(PC8) + scale(PC9) + scale(PC10) + scale(PC11) + scale(PC12) + scale(PC13) + scale(PC14) + scale(PC15), data = males))


summary(lm(scale(t1_fsiq) ~ scale(scz_1) + scale(t1_age) + as.character(t1_iqtype) + scale(PC1) + scale(PC2) + scale(PC3) + scale(PC4) + scale(PC5) +
             scale(PC6) + scale(PC7) + scale(PC8) + scale(PC9) + scale(PC10) + scale(PC11) + scale(PC12) + scale(PC13) + scale(PC14) + scale(PC15), data = males))



##Run regressions - NVIQ


dim(as.table(na.omit(males$t1_piq)))


summary(lm(scale(t1_piq) ~ scale(IQ_all) + scale(t1_age) + as.character(t1_iqtype) + scale(PC1) + scale(PC2) + scale(PC3) + scale(PC4) + scale(PC5) +
             scale(PC6) + scale(PC7) + scale(PC8) + scale(PC9) + scale(PC10) + scale(PC11) + scale(PC12) + scale(PC13) + scale(PC14) + scale(PC15), data = males))


summary(lm(scale(t1_piq) ~ scale(edu_all) + scale(t1_age) + as.character(t1_iqtype) + scale(PC1) + scale(PC2) + scale(PC3) + scale(PC4) + scale(PC5) +
             scale(PC6) + scale(PC7) + scale(PC8) + scale(PC9) + scale(PC10) + scale(PC11) + scale(PC12) + scale(PC13) + scale(PC14) + scale(PC15), data = males))


summary(lm(scale(t1_piq) ~ scale(scz_1) + scale(t1_age) + as.character(t1_iqtype) + scale(PC1) + scale(PC2) + scale(PC3) + scale(PC4) + scale(PC5) +
             scale(PC6) + scale(PC7) + scale(PC8) + scale(PC9) + scale(PC10) + scale(PC11) + scale(PC12) + scale(PC13) + scale(PC14) + scale(PC15), data = males))



### Run regressions: ADOS


dim(as.table(na.omit(males$t1_adi_communication_total)))

summary(lm(scale(t1_adi_communication_total) ~ scale(IQ_all) + scale(t1_age)  + scale(PC1) + scale(PC2) + scale(PC3) + scale(PC4) + scale(PC5) +
             scale(PC6) + scale(PC7) + scale(PC8) + scale(PC9) + scale(PC10) + scale(PC11) + scale(PC12) + scale(PC13) + scale(PC14) + scale(PC15), data = males))


summary(lm(scale(t1_adi_communication_total) ~ scale(edu_all) + scale(t1_age)  + scale(PC1) + scale(PC2) + scale(PC3) + scale(PC4) + scale(PC5) +
             scale(PC6) + scale(PC7) + scale(PC8) + scale(PC9) + scale(PC10) + scale(PC11) + scale(PC12) + scale(PC13) + scale(PC14) + scale(PC15), data = males))


summary(lm(scale(t1_adi_communication_total) ~ scale(scz_1) + scale(t1_age) + scale(PC1) + scale(PC2) + scale(PC3) + scale(PC4) + scale(PC5) +
             scale(PC6) + scale(PC7) + scale(PC8) + scale(PC9) + scale(PC10) + scale(PC11) + scale(PC12) + scale(PC13) + scale(PC14) + scale(PC15), data = males))





## Run regressions ADI rrb
dim(as.table(na.omit(males$t1_adi_rrb_total)))

summary(lm(scale(t1_adi_rrb_total) ~ scale(IQ_all)  + scale(t1_age)  + scale(PC1) + scale(PC2) + scale(PC3) + scale(PC4) + scale(PC5) +
             scale(PC6) + scale(PC7) + scale(PC8) + scale(PC9) + scale(PC10) + scale(PC11) + scale(PC12) + scale(PC13) + scale(PC14) + scale(PC15), data = males))


summary(lm(scale(t1_adi_rrb_total) ~ scale(edu_all) + scale(t1_age)  + scale(PC1) + scale(PC2) + scale(PC3) + scale(PC4) + scale(PC5) +
             scale(PC6) + scale(PC7) + scale(PC8) + scale(PC9) + scale(PC10) + scale(PC11) + scale(PC12) + scale(PC13) + scale(PC14) + scale(PC15), data = males))


summary(lm(scale(t1_adi_rrb_total) ~ scale(scz_1) + scale(t1_age) + scale(PC1) + scale(PC2) + scale(PC3) + scale(PC4) + scale(PC5) +
             scale(PC6) + scale(PC7) + scale(PC8) + scale(PC9) + scale(PC10) + scale(PC11) + scale(PC12) + scale(PC13) + scale(PC14) + scale(PC15), data = males))



## Run regressions ADI social
dim(as.table(na.omit(males$t1_adi_social_total)))

summary(lm(scale(t1_adi_social_total) ~ scale(IQ_all) + scale(t1_age)  + scale(PC1) + scale(PC2) + scale(PC3) + scale(PC4) + scale(PC5) +
             scale(PC6) + scale(PC7) + scale(PC8) + scale(PC9) + scale(PC10) + scale(PC11) + scale(PC12) + scale(PC13) + scale(PC14) + scale(PC15), data = males))


summary(lm(scale(t1_adi_social_total) ~ scale(edu_all) + scale(t1_age)  + scale(PC1) + scale(PC2) + scale(PC3) + scale(PC4) + scale(PC5) +
             scale(PC6) + scale(PC7) + scale(PC8) + scale(PC9) + scale(PC10) + scale(PC11) + scale(PC12) + scale(PC13) + scale(PC14) + scale(PC15), data = males))


summary(lm(scale(t1_adi_social_total) ~ scale(scz_1) + scale(t1_age) + scale(PC1) + scale(PC2) + scale(PC3) + scale(PC4) + scale(PC5) +
             scale(PC6) + scale(PC7) + scale(PC8) + scale(PC9) + scale(PC10) + scale(PC11) + scale(PC12) + scale(PC13) + scale(PC14) + scale(PC15), data = males))


## Run regressions on VABS
dim(as.table(na.omit(males$t1_vabsabcabc_standard)))


summary(lm(scale(t1_vabsabcabc_standard) ~ scale(IQ_all)  + scale(t1_age) + scale(PC1) + scale(PC2) + scale(PC3) + scale(PC4) + scale(PC5) +
             scale(PC6) + scale(PC7) + scale(PC8) + scale(PC9) + scale(PC10) + scale(PC11) + scale(PC12) + scale(PC13) + scale(PC14) + scale(PC15), data = males))


summary(lm(scale(t1_vabsabcabc_standard) ~ scale(edu_all) + scale(t1_age) + scale(PC1) + scale(PC2) + scale(PC3) + scale(PC4) + scale(PC5) +
             scale(PC6) + scale(PC7) + scale(PC8) + scale(PC9) + scale(PC10) + scale(PC11) + scale(PC12) + scale(PC13) + scale(PC14) + scale(PC15), data = males))


summary(lm(scale(t1_vabsabcabc_standard) ~ scale(scz_1) + scale(t1_age) + scale(PC1) + scale(PC2) + scale(PC3) + scale(PC4) + scale(PC5) +
             scale(PC6) + scale(PC7) + scale(PC8) + scale(PC9) + scale(PC10) + scale(PC11) + scale(PC12) + scale(PC13) + scale(PC14) + scale(PC15), data = males))

## Run regression on SRS

dim(as.table(na.omit(males$t1_srs_rawscore)))


summary(lm(scale(t1_srs_rawscore) ~ scale(IQ_all)  + scale(t1_age) + scale(PC1) + scale(PC2) + scale(PC3) + scale(PC4) + scale(PC5) +
             scale(PC6) + scale(PC7) + scale(PC8) + scale(PC9) + scale(PC10) + scale(PC11) + scale(PC12) + scale(PC13) + scale(PC14) + scale(PC15), data = males))


summary(lm(scale(t1_srs_rawscore) ~ scale(edu_all)  + scale(t1_age) + scale(PC1) + scale(PC2) + scale(PC3) + scale(PC4) + scale(PC5) +
             scale(PC6) + scale(PC7) + scale(PC8) + scale(PC9) + scale(PC10) + scale(PC11) + scale(PC12) + scale(PC13) + scale(PC14) + scale(PC15), data = males))


summary(lm(scale(t1_srs_rawscore) ~ scale(scz_1)  + scale(t1_age) + scale(PC1) + scale(PC2) + scale(PC3) + scale(PC4) + scale(PC5) +
             scale(PC6) + scale(PC7) + scale(PC8) + scale(PC9) + scale(PC10) + scale(PC11) + scale(PC12) + scale(PC13) + scale(PC14) + scale(PC15), data = males))


## Run regression on RBS-R

dim(as.table(na.omit(males$t1_rbs_total)))


summary(lm(scale(t1_rbs_total) ~ scale(IQ_all)  + scale(t1_age) + scale(PC1) + scale(PC2) + scale(PC3) + scale(PC4) + scale(PC5) +
             scale(PC6) + scale(PC7) + scale(PC8) + scale(PC9) + scale(PC10) + scale(PC11) + scale(PC12) + scale(PC13) + scale(PC14) + scale(PC15), data = males))


summary(lm(scale(t1_rbs_total) ~ scale(edu_all) + scale(PC1) + scale(PC2) + scale(PC3) + scale(PC4) + scale(PC5) +
             scale(PC6) + scale(PC7) + scale(PC8) + scale(PC9) + scale(PC10) + scale(PC11) + scale(PC12) + scale(PC13) + scale(PC14) + scale(PC15), data = males))


summary(lm(scale(t1_rbs_total) ~ scale(scz_1) + scale(t1_age) + scale(PC1) + scale(PC2) + scale(PC3) + scale(PC4) + scale(PC5) +
             scale(PC6) + scale(PC7) + scale(PC8) + scale(PC9) + scale(PC10) + scale(PC11) + scale(PC12) + scale(PC13) + scale(PC14) + scale(PC15), data = males))



## Run regression on ADOS_social_affect

dim(as.table(na.omit(males$t1_sa_css_all)))


summary(lm(scale(t1_sa_css_all) ~ scale(IQ_all)  + scale(t1_age) + scale(PC1) + scale(PC2) + scale(PC3) + scale(PC4) + scale(PC5) +
             scale(PC6) + scale(PC7) + scale(PC8) + scale(PC9) + scale(PC10) + scale(PC11) + scale(PC12) + scale(PC13) + scale(PC14) + scale(PC15), data = males))


summary(lm(scale(t1_sa_css_all) ~ scale(edu_all) + scale(PC1) + scale(PC2) + scale(PC3) + scale(PC4) + scale(PC5) +
             scale(PC6) + scale(PC7) + scale(PC8) + scale(PC9) + scale(PC10) + scale(PC11) + scale(PC12) + scale(PC13) + scale(PC14) + scale(PC15), data = males))


summary(lm(scale(t1_sa_css_all) ~ scale(scz_1) + scale(t1_age) + scale(PC1) + scale(PC2) + scale(PC3) + scale(PC4) + scale(PC5) +
             scale(PC6) + scale(PC7) + scale(PC8) + scale(PC9) + scale(PC10) + scale(PC11) + scale(PC12) + scale(PC13) + scale(PC14) + scale(PC15), data = males))




## Run regression on ADOS_rrb_affect

dim(as.table(na.omit(males$t1_rrb_css_all)))


summary(lm(scale(t1_rrb_css_all) ~ scale(IQ_all)   + scale(t1_age) + scale(PC1) + scale(PC2) + scale(PC3) + scale(PC4) + scale(PC5) +
             scale(PC6) + scale(PC7) + scale(PC8) + scale(PC9) + scale(PC10) + scale(PC11) + scale(PC12) + scale(PC13) + scale(PC14) + scale(PC15), data = males))


summary(lm(scale(t1_rrb_css_all) ~ scale(edu_all)  + scale(PC1) + scale(PC2) + scale(PC3) + scale(PC4) + scale(PC5) +
             scale(PC6) + scale(PC7) + scale(PC8) + scale(PC9) + scale(PC10) + scale(PC11) + scale(PC12) + scale(PC13) + scale(PC14) + scale(PC15), data = males))


summary(lm(scale(t1_rrb_css_all) ~ scale(scz_1) + scale(t1_age) + scale(PC1) + scale(PC2) + scale(PC3) + scale(PC4) + scale(PC5) +
             scale(PC6) + scale(PC7) + scale(PC8) + scale(PC9) + scale(PC10) + scale(PC11) + scale(PC12) + scale(PC13) + scale(PC14) + scale(PC15), data = males))


#################################Run on individuals with IQ > 70 ##############################################
###############################################################################################################



###Order

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

variable_list = list("t1_sa_css_all", "t1_rrb_css_all", "t1_adi_communication_total", "t1_adi_rrb_total", "t1_adi_social_total", "t1_rbs_total", "t1_vabsabcabc_standard",
                     "t1_srs_rawscore")



###Individuals with IQ < 70

## Run regression on ADOS_social_affect

IQ_table_under70 = as.data.frame((summary(lm(scale(t1_sa_css_all) ~ scale(IQ_all) + as.character(sex) + scale(t1_age) + scale(PC1) + scale(PC2) + scale(PC3) + scale(PC4) + scale(PC5) +
                                        scale(PC6) + scale(PC7) + scale(PC8) + scale(PC9) + scale(PC10) + scale(PC11) + scale(PC12) + scale(PC13) + scale(PC14) + scale(PC15), data = non_id)))$coefficients[2,])


Edu_table_under70 = as.data.frame((summary(lm(scale(t1_sa_css_all) ~ scale(edu_all) + as.character(sex) + scale(PC1) + scale(t1_age) + scale(PC2) + scale(PC3) + scale(PC4) + scale(PC5) +
                                        scale(PC6) + scale(PC7) + scale(PC8) + scale(PC9) + scale(PC10) + scale(PC11) + scale(PC12) + scale(PC13) + scale(PC14) + scale(PC15), data = non_id)))$coefficients[2,])


SCZ_table_under70 = as.data.frame((summary(lm(scale(t1_sa_css_all) ~ scale(scz_1) + as.character(sex) + scale(t1_age)  + scale(PC1) + scale(PC2) + scale(PC3) + scale(PC4) + scale(PC5) +
                                        scale(PC6) + scale(PC7) + scale(PC8) + scale(PC9) + scale(PC10) + scale(PC11) + scale(PC12) + scale(PC13) + scale(PC14) + scale(PC15), data = non_id)))$coefficients[2,])

table_under70 = dim(as.table(na.omit(non_id$t1_sa_css_all)))


for(i in variable_list){
  
  
  
  alpha = lm(paste0("scale(", i, ") ~ scale(IQ_all) + as.character(sex) + scale(t1_age) + scale(PC1) + scale(PC2) + scale(PC3) + scale(PC4) + scale(PC5) +
                      scale(PC6) + scale(PC7) + scale(PC8) + scale(PC9) + scale(PC10) + scale(PC11) + scale(PC12) + scale(PC13) + scale(PC14) + scale(PC15)"), data = non_id)
  
  IQ = (summary(alpha))$coefficients[2,]
  
  beta = lm(paste0("scale(", i, ") ~ scale(edu_all) + as.character(sex) +  scale(t1_age) +scale(PC1) + scale(PC2) + scale(PC3) + scale(PC4) + scale(PC5) +
                      scale(PC6) + scale(PC7) + scale(PC8) + scale(PC9) + scale(PC10) + scale(PC11) + scale(PC12) + scale(PC13) + scale(PC14) + scale(PC15)"), data = non_id)
  
  
  Edu = (summary(beta))$coefficients[2,]
  
  gamma = lm(paste0("scale(", i, ") ~ scale(scz_1) + as.character(sex) +  scale(t1_age) + scale(PC1) + scale(PC2) + scale(PC3) + scale(PC4) + scale(PC5) +
                      scale(PC6) + scale(PC7) + scale(PC8) + scale(PC9) + scale(PC10) + scale(PC11) + scale(PC12) + scale(PC13) + scale(PC14) + scale(PC15)"), data = non_id)
  
  
  SCZ = (summary(gamma))$coefficients[2,]
  
  IQ_table_under70 = cbind(IQ_table_under70, IQ)
  Edu_table_under70 = cbind(Edu_table_under70, Edu)
  SCZ_table_under70 = cbind(SCZ_table_under70, SCZ)
}




##Run correlations and create a correlogram


dev.off() ##Or make the plot window larger if you get the invalid graphics state message

data_for_matrix = merged_all_3[,c("t1_viq", "t1_piq", "t1_fsiq", "t1_adi_social_total", 
                                  "t1_adi_communication_total", "t1_adi_communication_total", "t1_adi_rrb_total", "t1_sa_css_all", 
                                  "t1_rrb_css_all", "t1_srs_rawscore", "t1_vabsabcabc_standard", "t1_rbs_total")]

correlationmatrix = cor(data_for_matrix, use = "complete.obs", method = "pearson")



hM <- format(round(correlationmatrix, 2))
col<- colorRampPalette(c("blue", "white", "red"))(20)
heatmap.2(correlationmatrix, trace = "none", col=col, cellnote=hM, notecol = "black")

LEAP_heatmap = heatmap.2(correlationmatrix, trace = "none", col=col)

res2_LEAP<-rcorr(as.matrix(data_for_matrix))
res2_results = flattenCorrMatrix(res2$r, res2$P)



### Run the analysis only in males


