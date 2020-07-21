### AGRE polygenic score analyses

prs_scz = fread("~/ALSPAC/PRSice2results/CHOP_scz_eas_eur_lam2019.all.score")

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

prs_edu = fread("~/ALSPAC/PRSice2results/CHOP_edu2018leeforPRSICE.all.score", header = TRUE)
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

prs_IQ = fread("~/ALSPAC/PRSice2results/CHOP_savageCP2018forprsice.all.score", header = TRUE)
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
summary(lmer(scale(CSEst_Nonverbal_IQ) ~ scale(IQ_all) + (1|FID) + as.character(sex) + scale(age) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
             scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = merged_all_ravens))


summary(lmer(scale(CSEst_Nonverbal_IQ) ~ scale(edu_all) + (1|FID) + as.character(sex) + scale(age) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = merged_all_ravens))


summary(lmer(scale(CSEst_Nonverbal_IQ) ~ scale(scz_1) + (1|FID) + as.character(sex) + scale(age) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = merged_all_ravens))



#Males

summary(lmer(scale(CSEst_Nonverbal_IQ) ~ scale(IQ_all) + (1|FID) + scale(age) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = males))


summary(lmer(scale(CSEst_Nonverbal_IQ) ~ scale(edu_all) + (1|FID)  + scale(age) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = males))


summary(lmer(scale(CSEst_Nonverbal_IQ) ~ scale(scz_1) + (1|FID)  + scale(age) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = males))



#Non-ID
summary(lmer(scale(CSEst_Nonverbal_IQ) ~ scale(IQ_all) + (1|FID) + as.character(sex) + scale(age) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = nonID))


summary(lmer(scale(CSEst_Nonverbal_IQ) ~ scale(edu_all) + (1|FID) + as.character(sex) + scale(age) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = nonID))


summary(lmer(scale(CSEst_Nonverbal_IQ) ~ scale(scz_1) + (1|FID) + as.character(sex) + scale(age) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = nonID))








### Read the phenotype data - ADOS
load("~/AGRE_data/Phenotypes/ADOS.RData")
setnames(ADOS, "Individual ID", "IID")
merged_all_ADOS = merge(merged_all, ADOS, by = "IID")

males = subset(merged_all_ADOS, sex == 1) ##450

ADOS_nonID = merged_all_ADOS[merged_all_ADOS$IID %in% nonID$IID,]

##Run regressions
summary(lmer(scale(CSSocAff2010_total) ~ scale(IQ_all) + as.character(Type) + (1|FID) + as.character(sex) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = merged_all_ADOS))


summary(lmer(scale(CSSocAff2010_total) ~ scale(edu_all) + as.character(Type) + (1|FID) + as.character(sex)  + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = merged_all_ADOS))


summary(lmer(scale(CSSocAff2010_total) ~ scale(scz_1) + as.character(Type) + (1|FID) + as.character(sex)  + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = merged_all_ADOS))



###Males-only

summary(lmer(scale(CSSocAff2010_total) ~ scale(IQ_all) + as.character(Type) + (1|FID)  + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = males))


summary(lmer(scale(CSSocAff2010_total) ~ scale(edu_all) + as.character(Type) + (1|FID)  + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = males))


summary(lmer(scale(CSSocAff2010_total) ~ scale(scz_1) + as.character(Type) + (1|FID)   + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = males))



###Non-ID
summary(lmer(scale(CSSocAff2010_total) ~ scale(IQ_all) + as.character(Type) + (1|FID) + as.character(sex) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = ADOS_nonID))


summary(lmer(scale(CSSocAff2010_total) ~ scale(edu_all) + as.character(Type) + (1|FID) + as.character(sex)  + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = ADOS_nonID))


summary(lmer(scale(CSSocAff2010_total) ~ scale(scz_1) + as.character(Type) + (1|FID) + as.character(sex)  + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = ADOS_nonID))





############# Run regressions - ADOS RRB##############################################################
summary(lmer(scale(CSRRB_2010_total) ~ scale(IQ_all) + as.character(Type) + (1|FID) + as.character(sex) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = merged_all_ADOS))


summary(lmer(scale(CSRRB_2010_total) ~ scale(edu_all) + as.character(Type) + (1|FID) + as.character(sex)  + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = merged_all_ADOS))


summary(lmer(scale(CSRRB_2010_total) ~ scale(scz_1) + as.character(Type) + (1|FID) + as.character(sex)  + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = merged_all_ADOS))



###Males-only

summary(lmer(scale(CSRRB_2010_total) ~ scale(IQ_all) + as.character(Type) + (1|FID)  + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = males))


summary(lmer(scale(CSRRB_2010_total) ~ scale(edu_all) + as.character(Type) + (1|FID)  + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = males))


summary(lmer(scale(CSRRB_2010_total) ~ scale(scz_1) + as.character(Type) + (1|FID)  + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = males))



###Non-ID
summary(lmer(scale(CSRRB_2010_total) ~ scale(IQ_all) + as.character(Type) + (1|FID) + as.character(sex) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = ADOS_nonID))


summary(lmer(scale(CSRRB_2010_total) ~ scale(edu_all) + as.character(Type) + (1|FID) + as.character(sex)  + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = ADOS_nonID))


summary(lmer(scale(CSRRB_2010_total) ~ scale(scz_1) + as.character(Type) + (1|FID) + as.character(sex)  + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = ADOS_nonID))





############### Read the phenotype data - ADI-R#########################################

ADI_R = fread("~/AGRE_data/Phenotypes/ADIR2.xls")
setnames(ADI_R, "Individual ID", "IID")
merged_all_ADI = merge(merged_all, ADI_R, by = "IID")

males = subset(merged_all_ADI, sex == 1 )

ADI_nonID = merged_all_ADI[merged_all_ADI$IID %in% nonID$IID,]

## Run regressions 1
summary(lmer(scale(COMVT_CS) ~ scale(IQ_all) + (1|FID) + as.character(sex) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = merged_all_ADI))


summary(lmer(scale(COMVT_CS) ~ scale(edu_all) + (1|FID) + as.character(sex)  + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = merged_all_ADI))


summary(lmer(scale(COMVT_CS) ~ scale(scz_1) + (1|FID)  + as.character(sex) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = merged_all_ADI))


##Males-only
summary(lmer(scale(COMVT_CS) ~ scale(IQ_all) + (1|FID)  + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = males))


summary(lmer(scale(COMVT_CS) ~ scale(edu_all) + (1|FID)   + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = males))


summary(lmer(scale(COMVT_CS) ~ scale(scz_1) + (1|FID)   + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = males))



##non-ID
summary(lmer(scale(COMVT_CS) ~ scale(IQ_all) + (1|FID) + as.character(sex) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = ADI_nonID))


summary(lmer(scale(COMVT_CS) ~ scale(edu_all) + (1|FID) + as.character(sex)  + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = ADI_nonID))


summary(lmer(scale(COMVT_CS) ~ scale(scz_1) + (1|FID)  + as.character(sex) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = ADI_nonID))






## Run regressions 2 #####
summary(lmer(scale(COMNVTCS) ~ scale(IQ_all) + (1|FID) + as.character(sex) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = merged_all_ADI))


summary(lmer(scale(COMNVTCS) ~ scale(edu_all) + (1|FID) + as.character(sex)  + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = merged_all_ADI))


summary(lmer(scale(COMNVTCS) ~ scale(scz_1) + (1|FID) + as.character(sex)  + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = merged_all_ADI))



##Males-only
summary(lmer(scale(COMNVTCS) ~ scale(IQ_all) + (1|FID)  + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = males))


summary(lmer(scale(COMNVTCS) ~ scale(edu_all) + (1|FID)   + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = males))


summary(lmer(scale(COMNVTCS) ~ scale(scz_1) + (1|FID)   + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = males))



##non-ID
summary(lmer(scale(COMNVTCS) ~ scale(IQ_all) + (1|FID) + as.character(sex) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = ADI_nonID))


summary(lmer(scale(COMNVTCS) ~ scale(edu_all) + (1|FID) + as.character(sex)  + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = ADI_nonID))


summary(lmer(scale(COMNVTCS) ~ scale(scz_1) + (1|FID)  + as.character(sex) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = ADI_nonID))




##Males-only
summary(lmer(scale(COMNVTCS) ~ scale(IQ_all) + (1|FID)  + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = males))


summary(lmer(scale(COMNVTCS) ~ scale(edu_all) + (1|FID)   + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = males))


summary(lmer(scale(COMNVTCS) ~ scale(scz_1) + (1|FID)   + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = males))



##non-ID
summary(lmer(scale(COMNVTCS) ~ scale(IQ_all) + (1|FID) + as.character(sex) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = ADI_nonID))


summary(lmer(scale(COMNVTCS) ~ scale(edu_all) + (1|FID) + as.character(sex)  + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = ADI_nonID))


summary(lmer(scale(COMNVTCS) ~ scale(scz_1) + (1|FID)  + as.character(sex) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = ADI_nonID))







### Run regressions 3 SOCT_CS


summary(lmer(scale(SOCT_CS) ~ scale(IQ_all) + (1|FID) + as.character(sex) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = merged_all_ADI))


summary(lmer(scale(SOCT_CS) ~ scale(edu_all) + (1|FID) + as.character(sex)  + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = merged_all_ADI))


summary(lmer(scale(SOCT_CS) ~ scale(scz_1) + (1|FID) + as.character(sex)  + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = merged_all_ADI))



##Males-only
summary(lmer(scale(SOCT_CS) ~ scale(IQ_all) + (1|FID)  + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = males))


summary(lmer(scale(SOCT_CS) ~ scale(edu_all) + (1|FID)   + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = males))


summary(lmer(scale(SOCT_CS) ~ scale(scz_1) + (1|FID)   + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = males))



##non-ID
summary(lmer(scale(SOCT_CS) ~ scale(IQ_all) + (1|FID) + as.character(sex) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = ADI_nonID))


summary(lmer(scale(SOCT_CS) ~ scale(edu_all) + (1|FID) + as.character(sex)  + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = ADI_nonID))


summary(lmer(scale(SOCT_CS) ~ scale(scz_1) + (1|FID)  + as.character(sex) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = ADI_nonID))




### Run regressions 4 BEHT_CS

summary(lmer(scale(BEHT_CS) ~ scale(IQ_all) + (1|FID) + as.character(sex) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = merged_all_ADI))


summary(lmer(scale(BEHT_CS) ~ scale(edu_all) + (1|FID) + as.character(sex)  + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = merged_all_ADI))


summary(lmer(scale(BEHT_CS) ~ scale(scz_1) + (1|FID) + as.character(sex)  + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = merged_all_ADI))




##Males-only
summary(lmer(scale(BEHT_CS) ~ scale(IQ_all) + (1|FID)  + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = males))


summary(lmer(scale(BEHT_CS) ~ scale(edu_all) + (1|FID)   + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = males))


summary(lmer(scale(BEHT_CS) ~ scale(scz_1) + (1|FID)   + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = males))



##non-ID
summary(lmer(scale(BEHT_CS) ~ scale(IQ_all) + (1|FID) + as.character(sex) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = ADI_nonID))


summary(lmer(scale(BEHT_CS) ~ scale(edu_all) + (1|FID) + as.character(sex)  + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = ADI_nonID))


summary(lmer(scale(BEHT_CS) ~ scale(scz_1) + (1|FID)  + as.character(sex) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = ADI_nonID))




############### Read the phenotype data - RBS-R #########################################

RBS = fread("~/AGRE_data/Phenotypes/Repetitive_Behavior_Scales1.xls")
setnames(RBS, "Individual ID", "IID")
merged_all_RBS = merge(merged_all, RBS, by = "IID")


males = subset(merged_all_RBS, sex == 1 )

RBS_nonID = merged_all_RBS[merged_all_RBS$IID %in% nonID$IID,]


### Run regressions 1 Overall_score

summary(lmer(scale(Overall_score) ~ scale(IQ_all) + (1|FID) + as.character(sex) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = merged_all_RBS))


summary(lmer(scale(Overall_score) ~ scale(edu_all) + (1|FID) + as.character(sex)  + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = merged_all_RBS))


summary(lmer(scale(Overall_score) ~ scale(scz_1) + (1|FID) + as.character(sex)  + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = merged_all_RBS))




###Males-only

summary(lmer(scale(Overall_score) ~ scale(IQ_all) + (1|FID)  + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = males))


summary(lmer(scale(Overall_score) ~ scale(edu_all) + (1|FID)  + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = males))


summary(lmer(scale(Overall_score) ~ scale(scz_1) + (1|FID)  + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = males))



##NonID

summary(lmer(scale(Overall_score) ~ scale(IQ_all) + (1|FID) + as.character(sex) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = RBS_nonID))


summary(lmer(scale(Overall_score) ~ scale(edu_all) + (1|FID) + as.character(sex)  + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = RBS_nonID))


summary(lmer(scale(Overall_score) ~ scale(scz_1) + (1|FID) + as.character(sex)  + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = RBS_nonID))



############### Read the phenotype data - VABS #########################################

VABS = fread("~/AGRE_data/Phenotypes/VINE1.xls")
setnames(VABS, "Individual ID", "IID")
merged_all_VABS = merge(merged_all, VABS, by = "IID")

dim(as.table(na.omit(merged_all_VABS$Composite_Standard_Score)))

males = subset(merged_all_VABS, sex == 1 )

VABS_nonID = merged_all_VABS[merged_all_VABS$IID %in% nonID$IID,]



### Run regressions 1 Composite_Standard_Score

summary(lmer(scale(Composite_Standard_Score) ~ scale(IQ_all) + (1|FID) + as.character(sex) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = merged_all_VABS))


summary(lmer(scale(Composite_Standard_Score) ~ scale(edu_all) + (1|FID) + as.character(sex)  + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = merged_all_VABS))


summary(lmer(scale(Composite_Standard_Score) ~ scale(scz_1) + (1|FID) + as.character(sex)  + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = merged_all_VABS))



###Males

summary(lmer(scale(Composite_Standard_Score) ~ scale(IQ_all) + (1|FID) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = males))


summary(lmer(scale(Composite_Standard_Score) ~ scale(edu_all) + (1|FID) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = males))


summary(lmer(scale(Composite_Standard_Score) ~ scale(scz_1) + (1|FID) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = males))




### Non_ID

summary(lmer(scale(Composite_Standard_Score) ~ scale(IQ_all) + (1|FID) + as.character(sex) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = VABS_nonID))


summary(lmer(scale(Composite_Standard_Score) ~ scale(edu_all) + (1|FID) + as.character(sex)  + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = VABS_nonID))


summary(lmer(scale(Composite_Standard_Score) ~ scale(scz_1) + (1|FID) + as.character(sex)  + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = VABS_nonID))









############### Read the phenotype data - SRS #########################################

SRS = fread("~/AGRE_data/Phenotypes/SRS2_SRS20021.xls")
setnames(SRS, "Individual ID", "IID")
merged_all_SRS = merge(merged_all, SRS, by = "IID")
merged_all_SRS_parent = subset(merged_all_SRS, SRS_Respond == "1")

dim(as.table(na.omit(merged_all_SRS_parent$Raw_total)))


males = subset(merged_all_SRS_parent, sex == 1 )

SRS_nonID = merged_all_SRS_parent[merged_all_SRS_parent$IID %in% nonID$IID,]


### Run regressions 1 Composite_Standard_Score

summary(lmer(scale(Raw_total) ~ scale(IQ_all) + (1|FID) + as.character(sex) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = merged_all_SRS_parent))


summary(lmer(scale(Raw_total) ~ scale(edu_all) + (1|FID) + as.character(sex)  + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = merged_all_SRS_parent))


summary(lmer(scale(Raw_total) ~ scale(scz_1) + (1|FID) + as.character(sex)  + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = merged_all_SRS_parent))

###males-only

summary(lmer(scale(Raw_total) ~ scale(IQ_all) + (1|FID) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = males))


summary(lmer(scale(Raw_total) ~ scale(edu_all) + (1|FID)  + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = males))


summary(lmer(scale(Raw_total) ~ scale(scz_1) + (1|FID)  + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = males))


###non-ID

summary(lmer(scale(Raw_total) ~ scale(IQ_all) + (1|FID) + as.character(sex) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = SRS_nonID))


summary(lmer(scale(Raw_total) ~ scale(edu_all) + (1|FID) + as.character(sex)  + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = SRS_nonID))


summary(lmer(scale(Raw_total) ~ scale(scz_1) + (1|FID) + as.character(sex)  + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = SRS_nonID))





############### Read the phenotype data - Stanford Binet #########################################

IQ = fread("~/AGRE_data/Phenotypes/Stanford_Binet1.xls")
setnames(IQ, "Individual ID", "IID")
merged_all_IQ = merge(merged_all, IQ, by = "IID")
merged_all_FSIQ = subset(merged_all_IQ, SSFSIQ > 0) 

dim(as.table(na.omit(merged_all_FSIQ$SSFSIQ)))

males = subset(merged_all_FSIQ, sex == 1)

merged_all_VIQ = subset(merged_all_IQ, SSVIQ > 0)

dim(as.table(na.omit(merged_all_VIQ$SSVIQ)))

### Run regressions 1 FSIQ

summary(lmer(scale(SSFSIQ) ~ scale(IQ_all) + (1|FID) + as.character(sex) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = merged_all_FSIQ))


summary(lmer(scale(SSFSIQ) ~ scale(edu_all) + (1|FID) + as.character(sex)  + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = merged_all_FSIQ))


summary(lmer(scale(SSFSIQ) ~ scale(scz_1) + (1|FID) + as.character(sex)  + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = merged_all_FSIQ))



### Run regressions 2 VIQ

summary(lmer(scale(SSVIQ) ~ scale(IQ_all) + (1|FID) + as.character(sex) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = merged_all_VIQ))


summary(lmer(scale(SSVIQ) ~ scale(edu_all) + (1|FID) + as.character(sex)  + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = merged_all_VIQ))


summary(lmer(scale(SSVIQ) ~ scale(scz_1) + (1|FID) + as.character(sex)  + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) +
               scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) + scale(V13) + scale(V14) + scale(V15) + scale(V16) + scale(V17), data = merged_all_VIQ))

### Scripts for correlations in ID
Ravens  = fread("~/AGRE_data/Phenotypes/Raven1.xls")
setnames(Ravens, "Individual ID", "IID")
merged_all_2 = merge(merged_all, Ravens, by = "IID", all.x = TRUE) ##583 observations
merged_all_2 = merged_all_2[,c("IID", "FID", "CSEst_Nonverbal_IQ")]
merged_all_2$CSEst_Nonverbal_IQ = ifelse(merged_all_2$CSEst_Nonverbal_IQ == "Untestable", NA, merged_all_2$CSEst_Nonverbal_IQ)
merged_all_2$CSEst_Nonverbal_IQ = ifelse(merged_all_2$CSEst_Nonverbal_IQ == "BTN", NA, merged_all_2$CSEst_Nonverbal_IQ)
merged_all_2$CSEst_Nonverbal_IQ = ifelse(merged_all_2$CSEst_Nonverbal_IQ == "ATN", NA, merged_all_2$CSEst_Nonverbal_IQ)


IQ = fread("~/AGRE_data/Phenotypes/Stanford_Binet1.xls")
setnames(IQ, "Individual ID", "IID")
IQ = IQ[,c("IID", "SSFSIQ", "SSVIQ")]
IQ$SSFSIQ = ifelse(IQ$SSFSIQ < 0, NA, IQ$SSFSIQ)
IQ$SSVIQ = ifelse(IQ$SSVIQ < 0, NA, IQ$SSVIQ)
merged_all_2 = merge(merged_all_2, IQ, by = "IID", all.x = TRUE)


SRS = fread("~/AGRE_data/Phenotypes/SRS2_SRS20021.xls")
setnames(SRS, "Individual ID", "IID")
SRS = SRS[,c("IID", "Raw_total")]
merged_all_2 = merge(merged_all_2, SRS, by = "IID", all.x = TRUE)


VABS = fread("~/AGRE_data/Phenotypes/VINE1.xls")
setnames(VABS, "Individual ID", "IID")
VABS = VABS[,c("IID", "Composite_Standard_Score")]
VABS$Composite_Standard_Score = ifelse(VABS$Composite_Standard_Score < 0, NA, VABS$Composite_Standard_Score)
merged_all_2 = merge(merged_all_2, VABS, by = "IID", all.x = TRUE)



RBS = fread("~/AGRE_data/Phenotypes/Repetitive_Behavior_Scales1.xls")
setnames(RBS, "Individual ID", "IID")
RBS = RBS[,c("IID", "Overall_score")]
RBS$Overall_Score = ifelse(RBS$Overall_Score < 0, NA, RBS$Overall_Score)
merged_all_2 = merge(merged_all_2, RBS, by = "IID", all.x = TRUE)


ADI_R = fread("~/AGRE_data/Phenotypes/ADIR2.xls")
setnames(ADI_R, "Individual ID", "IID")
ADI_R = ADI_R[,c("IID", "BEHT_CS", "SOCT_CS", "COMNVTCS", "COMVT_CS")]
ADI_R$BEHT_CS = ifelse(ADI_R$BEHT_CS < 0, NA, ADI_R$BEHT_CS)
ADI_R$SOCT_CS = ifelse(ADI_R$SOCT_CS  < 0, NA, ADI_R$SOCT_CS)
ADI_R$COMNVTCS = ifelse(ADI_R$COMNVTCS  < 0, NA, ADI_R$COMNVTCS)
ADI_R$COMVT_CS = ifelse(ADI_R$COMVT_CS  < 0, NA, ADI_R$COMVT_CS)
merged_all_2 = merge(merged_all_2, ADI_R, by = "IID", all.x = TRUE)

load("~/AGRE_data/Phenotypes/ADOS.RData")
setnames(ADOS, "Individual ID", "IID")
ADOS = ADOS[,c("IID", "CSSocAff2010_total", "CSRRB_2010_total")]
ADOS$CSSocAff2010_total = ifelse(ADOS$CSSocAff2010_total  < 0, NA, ADOS$CSSocAff2010_total)
ADOS$CSRRB_2010_total = ifelse(ADOS$CSRRB_2010_total  < 0, NA, ADOS$CSRRB_2010_total)
merged_all_2 = merge(merged_all_2, ADOS, by = "IID", all.x = TRUE)

merged_all_2$CSEst_Nonverbal_IQ = as.numeric(as.character(merged_all_2$CSEst_Nonverbal_IQ))

merged_all_3 = merged_all_2[!duplicated(merged_all_2$FID), ]


data_for_matrix = merged_all_3[,-c("IID", "FID")]

correlationmatrix = cor(data_for_matrix, method = "pearson", use = "pairwise.complete.obs")


flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}



### Scripts for correlations in ID
Ravens  = fread("~/AGRE_data/Phenotypes/Raven1.xls")
setnames(Ravens, "Individual ID", "IID")
setnames(Ravens, "AU/Family ID", "FID")
merged_all_2 = Ravens
merged_all_2 = merged_all_2[,c("IID", "FID", "CSEst_Nonverbal_IQ")]
merged_all_2$CSEst_Nonverbal_IQ = ifelse(merged_all_2$CSEst_Nonverbal_IQ == "Untestable", NA, merged_all_2$CSEst_Nonverbal_IQ)
merged_all_2$CSEst_Nonverbal_IQ = ifelse(merged_all_2$CSEst_Nonverbal_IQ == "BTN", NA, merged_all_2$CSEst_Nonverbal_IQ)
merged_all_2$CSEst_Nonverbal_IQ = ifelse(merged_all_2$CSEst_Nonverbal_IQ == "ATN", NA, merged_all_2$CSEst_Nonverbal_IQ)


IQ = fread("~/AGRE_data/Phenotypes/Stanford_Binet1.xls")
setnames(IQ, "Individual ID", "IID")
IQ = IQ[,c("IID", "SSFSIQ", "SSVIQ")]
IQ$SSFSIQ = ifelse(IQ$SSFSIQ < 0, NA, IQ$SSFSIQ)
IQ$SSVIQ = ifelse(IQ$SSVIQ < 0, NA, IQ$SSVIQ)
merged_all_2 = merge(merged_all_2, IQ, by = "IID", all.x = TRUE)


SRS = fread("~/AGRE_data/Phenotypes/SRS2_SRS20021.xls")
setnames(SRS, "Individual ID", "IID")
SRS = SRS[,c("IID", "Raw_total")]
merged_all_2 = merge(merged_all_2, SRS, by = "IID", all.x = TRUE)


VABS = fread("~/AGRE_data/Phenotypes/VINE1.xls")
setnames(VABS, "Individual ID", "IID")
VABS = VABS[,c("IID", "Composite_Standard_Score")]
VABS$Composite_Standard_Score = ifelse(VABS$Composite_Standard_Score < 0, NA, VABS$Composite_Standard_Score)
merged_all_2 = merge(merged_all_2, VABS, by = "IID", all.x = TRUE)



RBS = fread("~/AGRE_data/Phenotypes/Repetitive_Behavior_Scales1.xls")
setnames(RBS, "Individual ID", "IID")
RBS = RBS[,c("IID", "Overall_score")]
RBS$Overall_Score = ifelse(RBS$Overall_Score < 0, NA, RBS$Overall_Score)
merged_all_2 = merge(merged_all_2, RBS, by = "IID", all.x = TRUE)


ADI_R = fread("~/AGRE_data/Phenotypes/ADIR2.xls")
setnames(ADI_R, "Individual ID", "IID")
ADI_R = ADI_R[,c("IID", "BEHT_CS", "SOCT_CS", "COMNVTCS", "COMVT_CS")]
ADI_R$BEHT_CS = ifelse(ADI_R$BEHT_CS < 0, NA, ADI_R$BEHT_CS)
ADI_R$SOCT_CS = ifelse(ADI_R$SOCT_CS  < 0, NA, ADI_R$SOCT_CS)
ADI_R$COMNVTCS = ifelse(ADI_R$COMNVTCS  < 0, NA, ADI_R$COMNVTCS)
ADI_R$COMVT_CS = ifelse(ADI_R$COMVT_CS  < 0, NA, ADI_R$COMVT_CS)
merged_all_2 = merge(merged_all_2, ADI_R, by = "IID", all.x = TRUE)

load("~/AGRE_data/Phenotypes/ADOS.RData")
setnames(ADOS, "Individual ID", "IID")
ADOS = ADOS[,c("IID", "CSSocAff2010_total", "CSRRB_2010_total")]
ADOS$CSSocAff2010_total = ifelse(ADOS$CSSocAff2010_total  < 0, NA, ADOS$CSSocAff2010_total)
ADOS$CSRRB_2010_total = ifelse(ADOS$CSRRB_2010_total  < 0, NA, ADOS$CSRRB_2010_total)
merged_all_2 = merge(merged_all_2, ADOS, by = "IID", all.x = TRUE)

merged_all_2$CSEst_Nonverbal_IQ = as.numeric(as.character(merged_all_2$CSEst_Nonverbal_IQ))

merged_all_3 = merged_all_2[!duplicated(merged_all_2$FID), ]


data_for_matrix = merged_all_3[,-c("IID", "FID")]

correlationmatrix = cor(data_for_matrix, method = "pearson", use = "pairwise.complete.obs")


flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

res2<-rcorr(as.matrix(data_for_matrix))
res2_results = flattenCorrMatrix(res2$r, res2$P)


hM <- format(round(correlationmatrix, 2))
col<- colorRampPalette(c("blue", "white", "red"))(20)
heatmap.2(correlationmatrix, trace = "none", col=col, cellnote=hM, notecol = "black")

LEAP_heatmap = heatmap.2(correlationmatrix, trace = "none", col=col)




