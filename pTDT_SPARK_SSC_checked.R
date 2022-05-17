## PTDT script
### Author: Varun Warrier


### Paper: https://www.nature.com/articles/ng.3863#methods

###Formulas
###PGSmidparent = (PGSfather + PGSmother)/2 ; 
###pTDT deviation is (PGSchild - PGSmidparent)/SD(PGSmidparent);
###tPDT = mean (pTDT deviation)/ (SD(pTDT deviation)/sqrt(N))

#Checked on 7th June 2021

## STEP 1: Read files and merge

library(data.table)

#First, read the fam files
Mother = fread("~/SFARI/mothersdata.txt", header= T)
Father = read.table("~/SFARI/fathersdata.txt", header= T)
Cases = fread("~/SFARI/cases.txt", header= T)
siblings = fread("~/SFARI/siblingsdata.txt", header= T)

#Next, read the prs scores
PRS = fread("./PGS/SSC_autism_finalscore.profile", header = TRUE)

#Create merged files
fatherpgs = merge(Father, PRS, by = "IID")
motherpgs = merge(Mother, PRS, by = "IID")
casespgs = merge(Cases, PRS, by = "IID")

##Step 2: Calculate mid-parent PGS and SD of midparent PGS

parentpgs = merge(motherpgs, fatherpgs, by = "FID.x")
parentpgs$midparent = (parentpgs$`SCORE.x` + parentpgs$`SCORE.y`)/2


triopgs = merge(parentpgs, casespgs, by = "FID.x")

Sd_score = sd(triopgs$midparent)

## Step 3: Calculate pTDT deviation
triopgs$diff = (triopgs$SCORE - triopgs$midparent)/Sd_score


## Step 4: Calculate the T score for pTDT deviation
N = sqrt(nrow(triopgs))
Z = mean(triopgs$diff)/(sd(triopgs$diff)/N)
mean(triopgs$diff)
(sd(triopgs$diff)/N)
2*pnorm(-abs(Z))



###run it in de novo

SNV = fread("./Sebat_data/denovo_SNVindel.txt")
SV = fread("./Sebat_data/denovo_SV.txt")

##loeuf_upper_decile_score = 0.37
SNV_upper_decile = subset(SNV, gnomad_loeuf < 0.37)
SV_upper_decile = subset(SV, gnomad_loeuf_minimum < 0.37)

autism_SNV = subset(SNV_upper_decile, phenotype == "ASD")
autism_SV = subset(SV_upper_decile, phenotype == "ASD")

control_SNV = subset(SNV_upper_decile, phenotype == "Control")
control_SV = subset(SV_upper_decile, phenotype == "Control")

all = rbind(autism_SNV[,c("fid", "iid")], autism_SV[,c("fid", "iid")])
all_control = rbind(control_SNV[,c("fid", "iid")], control_SV[,c("fid", "iid")])

denovo_all = triopgs[triopgs$FID.x %in% all$fid,]
nondenovo_all = triopgs[!triopgs$FID.x %in% all$fid,]

N = sqrt(nrow(denovo_all))
mean(denovo_all$diff)
(sd(denovo_all$diff)/N)
Z = mean(denovo_all$diff)/(sd(denovo_all$diff)/N)
2*pnorm(-abs(Z))

N = sqrt(nrow(nondenovo_all))
mean(nondenovo_all$diff)
(sd(nondenovo_all$diff)/N)
Z = mean(nondenovo_all$diff)/(sd(nondenovo_all$diff)/N)
2*pnorm(-abs(Z))



### In siblings

siblingpgs = merge(siblings, PRS, by = "IID")
controltriopgs = merge(parentpgs, siblingpgs, by = "FID.x")

Sd_score = sd(controltriopgs$midparent)
controltriopgs$diff = (controltriopgs$SCORE - controltriopgs$midparent)/Sd_score
N = sqrt(nrow(controltriopgs))

Z = mean(controltriopgs$diff)/(sd(controltriopgs$diff)/N)
mean(controltriopgs$diff)
(sd(controltriopgs$diff)/N)

2*pnorm(-abs(Z))

rm(list = ls())



####SPARK_cohort_v2####

prs_autism = fread("./PGS/SPARK_autism_finalscore.profile", header = TRUE)
prs_autism2 = fread("./PGS/SPARK_v2_autism_finalscore.profile")
prs_autism = rbind(prs_autism, prs_autism2)
setnames(prs_autism, "SCORE", "autism_prs")
merged = prs_autism[,c("IID", "autism_prs", "FID")]

merged = unique(merged)

registration = fread("~/SPARK/Phenotypes/V5/individuals_registration.csv", fill = TRUE)
trios =  registration[!(registration$biofather_sp_id =="" | registration$biomother_sp_id==""), ] # 52936
asd_trios = subset(trios, asd == "TRUE") #35723

asd_trios = asd_trios[,c("subject_sp_id", "age_at_registration_months", "sex", "biomother_sp_id", "biofather_sp_id")]

merged_pgs_trio = merge(merged, asd_trios, by.x = "IID", by.y = "subject_sp_id") #7993

merged_pgs_trio = merged_pgs_trio[merged_pgs_trio$biomother_sp_id %in% merged$IID,] #5493
merged_pgs_trio = merged_pgs_trio[merged_pgs_trio$biofather_sp_id %in% merged$IID,] #4777

merged_final = merge(merged_pgs_trio, pheno, by = "IID", all.x  = TRUE)


fatherpgs = merged[merged$IID %in% asd_trios$biofather_sp_id,] # 5363
fatherpgs = fatherpgs[,c("IID", "FID", "autism_prs")]
motherpgs = merged[merged$IID %in% asd_trios$biomother_sp_id,]
motherpgs = motherpgs[,c("IID", "FID","autism_prs")]

casespgs = merged_pgs_trio


parentpgs = merge(motherpgs, fatherpgs, by = "FID")
parentpgs$midparent_autism = (parentpgs$autism_prs.x + parentpgs$autism_prs.y)/2

triopgs = merge(parentpgs, casespgs, by = "FID")
triopgs = triopgs[!duplicated(triopgs[,c("IID")]),] # 4747
Sd_autism = sd(triopgs$midparent_autism)




## Step 3: Calculate pTDT deviation
triopgs$autism_diff = (triopgs$autism_prs - triopgs$midparent_autism)/Sd_autism



## Step 4: Calculate the T score for pTDT deviation
N = sqrt(nrow(triopgs))
mean(triopgs$autism_diff)
sd(triopgs$autism_diff)/N

2*pnorm(-abs(mean(triopgs$autism_diff)/(sd(triopgs$autism_diff)/N)))





###For siblings


sibling_trios = subset(trios, asd == "FALSE") #17213

sibling_trios = sibling_trios[,c("subject_sp_id", "age_at_registration_months", "sex", "biomother_sp_id", "biofather_sp_id")]

merged_pgs_trio_sibling = merge(merged, sibling_trios, by.x = "IID", by.y = "subject_sp_id") #7993

merged_pgs_trio_sibling = merged_pgs_trio_sibling[merged_pgs_trio_sibling$biomother_sp_id %in% merged$IID,] #1964
merged_pgs_trio_sibling = merged_pgs_trio_sibling[merged_pgs_trio_sibling$biofather_sp_id %in% merged$IID,] #1715


fatherpgs = merged[merged$IID %in% sibling_trios$biofather_sp_id,] # 5363
fatherpgs = fatherpgs[,c("IID", "FID", "autism_prs")]
motherpgs = merged[merged$IID %in% sibling_trios$biomother_sp_id,]
motherpgs = motherpgs[,c("IID", "FID", "autism_prs")]

sibling = merged_pgs_trio_sibling


parentpgs = merge(motherpgs, fatherpgs, by = "FID")
parentpgs$midparent_autism = (parentpgs$autism_prs.x + parentpgs$autism_prs.y)/2

triopgs = merge(parentpgs, sibling, by = "FID")
triopgs = triopgs[!duplicated(triopgs[,c("IID")]),] # 4747
Sd_autism = sd(triopgs$midparent_autism)




## Step 3: Calculate pTDT deviation
triopgs$autism_diff = (triopgs$autism_prs - triopgs$midparent_autism)/Sd_autism

## Step 4: Calculate the T score for pTDT deviation
N = sqrt(nrow(triopgs))
mean(triopgs$autism_diff)
sd(triopgs$autism_diff)/N

2*pnorm(-abs(mean(triopgs$autism_diff)/(sd(triopgs$autism_diff)/N)))





###pTDT with denovo_only### 
###SPARK###
prs_autism = fread("./PGS/SPARK_autism_finalscore.profile", header = TRUE)
setnames(prs_autism, "SCORE", "autism_prs")
merged = prs_autism[,c("IID", "autism_prs", "FID")]
merged = unique(merged)

registration = fread("~/SPARK/Phenotypes/V5/individuals_registration.csv", fill = TRUE)
trios =  registration[!(registration$biofather_sp_id =="" | registration$biomother_sp_id==""), ] # 52936
asd_trios = subset(trios, asd == "TRUE") #35723

asd_trios = asd_trios[,c("subject_sp_id", "age_at_registration_months", "sex", "biomother_sp_id", "biofather_sp_id")]

merged_pgs_trio = merge(merged, asd_trios, by.x = "IID", by.y = "subject_sp_id") #7993

merged_pgs_trio = merged_pgs_trio[merged_pgs_trio$biomother_sp_id %in% merged$IID,] #5493
merged_pgs_trio = merged_pgs_trio[merged_pgs_trio$biofather_sp_id %in% merged$IID,] #4777


fatherpgs = merged[merged$IID %in% asd_trios$biofather_sp_id,] # 5363
fatherpgs = fatherpgs[,c("IID", "FID", "autism_prs")]
motherpgs = merged[merged$IID %in% asd_trios$biomother_sp_id,]
motherpgs = motherpgs[,c("IID", "FID","autism_prs")]

casespgs = merged_pgs_trio


parentpgs = merge(motherpgs, fatherpgs, by = "FID")
parentpgs$midparent_autism = (parentpgs$autism_prs.x + parentpgs$autism_prs.y)/2

triopgs = merge(parentpgs, casespgs, by = "FID")
triopgs = triopgs[!duplicated(triopgs[,c("IID")]),] # 4747


Sd_autism = sd(triopgs$midparent_autism)




## Step 3: Calculate pTDT deviation
triopgs$autism_diff = (triopgs$autism_prs - triopgs$midparent_autism)/Sd_autism


##de novo
SNV = fread("./Sebat_data/denovo_SNVindel.txt")
SNV = subset(SNV, cohort == "SPARK")


##loeuf_upper_decile_score = 0.37
SNV_upper_decile = subset(SNV, gnomad_loeuf < 0.37)


denovo = triopgs[triopgs$IID %in% SNV_upper_decile$iid,]
nondenovo = triopgs[!triopgs$IID %in% SNV_upper_decile$iid,]

N = sqrt(nrow(denovo))
mean(denovo$autism_diff)
sd(denovo$autism_diff)/N



