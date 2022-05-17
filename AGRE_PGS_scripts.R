### AGRE polygenic score analyses
library(data.table)
library(lme4)

setwd("~/Autism_heterogeneity")

prs_scz = fread("./PGS/CHOP_scz_finalscore.profile", header = TRUE)
setnames(prs_scz, "SCORE", "scz_prs")
prs_scz = prs_scz[,c("IID", "scz_prs")]

prs_edu = fread("./PGS/CHOP_edu_finalscore.profile", header = TRUE)
setnames(prs_edu, "SCORE", "edu_prs")
prs_edu = prs_edu[,c("IID", "edu_prs")]

prs_IQ = fread("./PGS/CHOP_IQ_finalscore.profile", header = TRUE)
setnames(prs_IQ, "SCORE", "IQ_prs")
prs_IQ = prs_IQ[,c("IID", "IQ_prs")]

prs_autism = fread("./PGS/CHOP_autism_finalscore.profile", header = TRUE)
setnames(prs_autism, "SCORE", "autism_prs")
prs_autism = prs_autism[,c("IID", "autism_prs")]

prs_haircolour = fread("./PGS/CHOP_haircolour_finalscore.profile", header = TRUE)
setnames(prs_haircolour, "SCORE", "haircolour_prs")
prs_haircolour = prs_haircolour[,c("IID", "haircolour_prs")]


prs_ADHD = fread("./PGS/CHOP_ADHD_finalscore.profile", header = TRUE)
setnames(prs_ADHD, "SCORE", "ADHD_prs")
prs_ADHD = prs_ADHD[,c("IID", "ADHD_prs")]


prs_merged = merge(prs_IQ, prs_edu, by = "IID")
prs_merged = merge(prs_merged, prs_scz, by = "IID")
prs_merged = merge(prs_merged, prs_autism, by = "IID")
prs_merged = merge(prs_merged, prs_haircolour, by = "IID")
prs_merged = merge(prs_merged, prs_ADHD, by = "IID")




foo <- data.frame(do.call('rbind', strsplit(as.character(prs_merged$IID),'_',fixed=TRUE)))
prs_merged = cbind(prs_merged, foo)
prs_merged = prs_merged[,-c("FID", "IID")]
setnames(prs_merged, "X1", "FID")
setnames(prs_merged, "X2", "IID")

##read PCs###
PC = fread("~/AGRE_data/CHOP/PCAall_chop.txt")
setnames(PC, "V2", "IID")
merged_all = merge(prs_merged, PC, by = "IID")




###Read_old_fam for sex

CHOPOldfam = fread("~/AGRE_data/CHOP/CHOPimputationfile.fam")
setnames(CHOPOldfam, "V2", "IID")
setnames(CHOPOldfam, "V5", "sex")

CHOPOldfam = CHOPOldfam[,c("IID", "sex")]

merged_all = merge(merged_all, CHOPOldfam, by = "IID")



### Read the phenotype data - Raven's#######################################

Ravens  = fread("~/AGRE_data/Phenotypes/Raven1.xls")

Ravens_1 = subset(Ravens, CSEst_Nonverbal_IQ != "Untestable")
Ravens_1 = Ravens_1[,c("Individual ID", "age", "Gender", "CSEst_Nonverbal_IQ")]
Ravens_1$CSEst_Nonverbal_IQ = as.numeric(as.character(Ravens_1$CSEst_Nonverbal_IQ))
setnames(Ravens_1, "Individual ID", "IID")

merged_all_ravens = merge(merged_all, Ravens_1, by = "IID") ##583 observations

males = subset(merged_all_ravens, sex == 1) ##450

nonID = subset(merged_all_ravens, CSEst_Nonverbal_IQ > 70)

##Run regressions
summary(lmer(scale(CSEst_Nonverbal_IQ) ~ scale(IQ_prs) + scale(edu_prs) +scale(scz_prs) + scale(autism_prs) + scale(ADHD_prs) + scale(haircolour_prs) + (1|FID) + as.character(sex) + scale(age) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
             scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13), data = merged_all_ravens))




### Read the phenotype data - ADOS
load("~/AGRE_data/Phenotypes/ADOS.RData")
setnames(ADOS, "Individual ID", "IID")
merged_all_ADOS = merge(merged_all, ADOS, by = "IID")

males = subset(merged_all_ADOS, sex == 1) ##450

ADOS_nonID = merged_all_ADOS[merged_all_ADOS$IID %in% nonID$IID,]

ADOS_with_IQ = merge(ADOS, merged_all_ravens, by = "IID")

##Run regressions
summary(lmer(scale(CSSocAff2010_total) ~ scale(IQ_prs) + scale(edu_prs) +scale(scz_prs) + scale(autism_prs) + scale(ADHD_prs) + scale(haircolour_prs) + as.character(Type) + (1|FID) + as.character(sex) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = merged_all_ADOS))





############# Run regressions - ADOS RRB##############################################################
summary(lmer(scale(CSRRB_2010_total) ~ scale(IQ_prs) + scale(edu_prs) +scale(scz_prs) + scale(autism_prs) + scale(ADHD_prs) + scale(haircolour_prs) + as.character(Type) + (1|FID) + as.character(sex) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = merged_all_ADOS))




############### Read the phenotype data - ADI-R#########################################

ADI_R = fread("~/AGRE_data/Phenotypes/ADIR2.xls")
setnames(ADI_R, "Individual ID", "IID")
merged_all_ADI = merge(merged_all, ADI_R, by = "IID")

males = subset(merged_all_ADI, sex == 1 )

ADI_nonID = merged_all_ADI[merged_all_ADI$IID %in% nonID$IID,]

ADI_with_IQ = merge(ADI_R, merged_all_ravens, by = "IID")

## Run regressions 1
summary(lmer(scale(COMVT_CS) ~ scale(IQ_prs) + scale(edu_prs) +scale(scz_prs) + scale(autism_prs) + scale(ADHD_prs) + scale(haircolour_prs) + (1|FID) + as.character(sex) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = merged_all_ADI))


### Run regressions 3 SOCT_CS


summary(lmer(scale(SOCT_CS) ~ scale(IQ_prs) + scale(edu_prs) +scale(scz_prs) + scale(autism_prs) + scale(ADHD_prs) + scale(haircolour_prs) + (1|FID) + as.character(sex) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = merged_all_ADI))



### Run regressions 4 BEHT_CS

summary(lmer(scale(BEHT_CS) ~ scale(IQ_prs) + scale(edu_prs) +scale(scz_prs) + scale(autism_prs)  + scale(ADHD_prs) + scale(haircolour_prs) + (1|FID) + as.character(sex) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = merged_all_ADI))



############### Read the phenotype data - RBS-R #########################################

RBS = fread("~/AGRE_data/Phenotypes/Repetitive_Behavior_Scales1.xls")
setnames(RBS, "Individual ID", "IID")
merged_all_RBS = merge(merged_all, RBS, by = "IID")

RBS_with_IQ = merge(RBS, merged_all_ravens, by = "IID")

males = subset(merged_all_RBS, sex == 1 )

RBS_nonID = merged_all_RBS[merged_all_RBS$IID %in% nonID$IID,]


### Run regressions 1 Overall_score

summary(lmer(scale(Overall_score) ~ scale(IQ_prs) + scale(edu_prs) +scale(scz_prs) + scale(autism_prs) + scale(ADHD_prs) + scale(haircolour_prs) + (1|FID) + as.character(sex) + scale(age) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = merged_all_RBS))



############### Read the phenotype data - VABS #########################################

VABS = fread("~/AGRE_data/Phenotypes/VINE1.xls")
setnames(VABS, "Individual ID", "IID")
merged_all_VABS = merge(merged_all, VABS, by = "IID")

dim(as.table(na.omit(merged_all_VABS$Composite_Standard_Score)))

VABS_with_IQ = merge(VABS, merged_all_ravens, by = "IID")

males = subset(merged_all_VABS, sex == 1 )

VABS_nonID = merged_all_VABS[merged_all_VABS$IID %in% nonID$IID,]



### Run regressions 1 Composite_Standard_Score

summary(lmer(scale(Composite_Standard_Score) ~ scale(IQ_prs) + scale(edu_prs) +scale(scz_prs) + scale(autism_prs) + scale(ADHD_prs) + scale(haircolour_prs) + (1|FID) + as.character(sex) + scale(age) +  scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = merged_all_VABS))




############### Read the phenotype data - SRS #########################################

SRS = fread("~/AGRE_data/Phenotypes/SRS2_SRS20021.xls")
setnames(SRS, "Individual ID", "IID")
merged_all_SRS = merge(merged_all, SRS, by = "IID")
merged_all_SRS_parent = subset(merged_all_SRS, SRS_Respond == "1")

dim(as.table(na.omit(merged_all_SRS_parent$Raw_total)))

SRS_with_IQ = merge(SRS, merged_all_ravens, by = "IID")


males = subset(merged_all_SRS_parent, sex == 1 )

SRS_nonID = merged_all_SRS_parent[merged_all_SRS_parent$IID %in% nonID$IID,]


### Run regressions 

summary(lmer(scale(Raw_total) ~ scale(IQ_prs) + scale(edu_prs) +scale(scz_prs) + scale(autism_prs) + scale(ADHD_prs) + scale(haircolour_prs) + (1|FID) + as.character(sex) + scale(age) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = merged_all_SRS_parent))




