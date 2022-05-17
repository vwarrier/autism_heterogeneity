#SNPs from SPARK
setwd("~/SPARK/SPARK_v2/Imputed_files/")
library(data.table)
a = fread("chr1.info")
b = subset(a, Rsq > 0.90)
c = subset(b, ALT_Frq > 0.05 & ALT_Frq < 0.95)

for (i in 2:22){
  a = fread(paste0("chr", i, ".info"))
  b = subset(a, Rsq > 0.90)
  b = subset(b, ALT_Frq > 0.05 & ALT_Frq < 0.95)
  c = rbind(b,c)
}

spark_snps = c # 8044396

##SNPs from SSC

setwd("~/SFARI/liftOverPlink/SSC_Omni/Imputed/")
a = fread("chr1.info")
a$Rsq = as.numeric(as.character(a$Rsq))
b = subset(a, Rsq > 0.90)
c = subset(b, ALT_Frq > 0.05 & ALT_Frq < 0.95)

for (i in 2:22){
  a = fread(paste0("chr", i, ".info"))
  a$Rsq = as.numeric(as.character(a$Rsq))
  b = subset(a, Rsq > 0.90)
  b = subset(b, ALT_Frq > 0.05 & ALT_Frq < 0.95)
  c = rbind(b,c)
}

ssc_omni_snps = c #8367716

setwd("~/SFARI/liftOverPlink/SSC_1Mv3/Imputed_1000G/")

a = fread("chr1.info")
a$Rsq = as.numeric(as.character(a$Rsq))
b = subset(a, Rsq > 0.90)
c = subset(b, ALT_Frq > 0.05 & ALT_Frq < 0.95)

for (i in 2:22){
  a = fread(paste0("chr", i, ".info"))
  a$Rsq = as.numeric(as.character(a$Rsq))
  b = subset(a, Rsq > 0.90)
  b = subset(b, ALT_Frq > 0.05 & ALT_Frq < 0.95)
  c = rbind(b,c)
}

ssc_1mv3_snps = c # 7588247

setwd("/mnt/b2/home4/arc/vw260/SFARI/liftOverPlink/SSC_1Mv1/Imputed/")

a = fread("chr1.info")
a$Rsq = as.numeric(as.character(a$Rsq))
b = subset(a, Rsq > 0.90)
c = subset(b, ALT_Frq > 0.05 & ALT_Frq < 0.95)

for (i in 2:22){
  a = fread(paste0("chr", i, ".info"))
  a$Rsq = as.numeric(as.character(a$Rsq))
  b = subset(a, Rsq > 0.90)
  b = subset(b, ALT_Frq > 0.05 & ALT_Frq < 0.95)
  c = rbind(b,c)
}

ssc_1mv1_snps = c #7588247


### SNPs from ABCD

setwd("/mnt/b2/home4/arc/vw260/ABCD/ABCDgenotype/Genotype_postimputation/TOPMED/")

a = fread("chr1.info")
a$Rsq = as.numeric(as.character(a$Rsq))
b = subset(a, Rsq > 0.90)
c = subset(b, ALT_Frq > 0.05 & ALT_Frq < 0.95)

for (i in 2:22){
  a = fread(paste0("chr", i, ".info"))
  a$Rsq = as.numeric(as.character(a$Rsq))
  b = subset(a, Rsq > 0.90)
  b = subset(b, ALT_Frq > 0.05 & ALT_Frq < 0.95)
  c = rbind(b,c)
}

abcd_snps = c #9933596

#convert snps to correct rsids
spark_recode = fread("~/SPARK/SPARK_postimputation/CEU/plinkrecodingfile2.txt")
spark_snps2 = spark_recode[spark_recode$V2 %in% spark_snps$SNP,]


##ssc_all
merged1 = merge(ssc_1mv1_snps, ssc_omni_snps, by = "SNP")
merged1 = merge(merged1, ssc_1mv3_snps, by = "SNP") # 7425595
ssc_recode = fread("~/SFARI/liftOverPlink/plinkrecodingfile.txt")
ssc_snps2 = ssc_recode[ssc_recode$chrompos %in% merged1$SNP,] #7216951

spark_snps2 = spark_recode[spark_recode$V2 %in% spark_snps$SNP,]

spark_ssc = ssc_snps2[ssc_snps2$ID %in% spark_snps2$ID,] #6344423 snps vs 7216136

abcd_snps2 = spark_recode[spark_recode$V2 %in% abcd_snps$SNP,]

spark_abcd_ssc = spark_ssc[spark_ssc$ID %in% abcd_snps2$ID,] # 6,114,631 ##266411 all combined

write.table(spark_abcd_ssc[,c("ID")], file = "~/Autism_heterogeneity/Commonwellimputed_SNPs_ABCD_SSC_SPARK_highr2.txt", row.names = F, col.names = T, quote = F)


./plink --bfile ~/SPARK/SPARK_postimputation/CEU/Plink_files/SPARKCEU_hg19_allchrs --keep ~/Autism_heterogeneity/GRMs/Autism_europeans_SPARK1_keep.txt --extract ~/Autism_heterogeneity/Commonwellimputed_SNPs_ABCD_SSC_SPARK.txt --make-bed --out ~/Autism_heterogeneity/GRMs/GRMs_inc_UKB/SPARK_v1_forABCD --threads 20
./plink --bfile ~/SPARK/SPARK_v2/Imputed_files/Plink_files/SPARK_hg19_allchrs_v2 --keep ~/Autism_heterogeneity/GRMs/Autism_europeans_SPARK2_keep.txt --extract ~/Autism_heterogeneity/Commonwellimputed_SNPs_ABCD_SSC_SPARK.txt --make-bed --out ~/Autism_heterogeneity/GRMs/GRMs_inc_UKB/SPARK_v2_forABCD --threads 20
./plink --bfile ~/SFARI/liftOverPlink/files_imputed/SFARImergedallcasesonly  --extract ~/Autism_heterogeneity/Commonwellimputed_SNPs_ABCD_SSC_SPARK.txt --make-bed --out ~/Autism_heterogeneity/GRMs/GRMs_inc_UKB/SSC_forABCD --threads 20
./plink --bfile ~/ABCD/ABCDgenotype/Genotype_postimputation/TOPMED/Plink_files/ABCD_hg19_allchrs_europeanonly  --extract ~/Autism_heterogeneity/Commonwellimputed_SNPs_ABCD_SSC_SPARK.txt --make-bed --out ~/Autism_heterogeneity/GRMs/GRMs_inc_UKB/ABCD_euroonly_forABCD --threads 20

##Combine all into 1 big happy (?) file


./plink --bfile ABCD_euroonly_forABCD --merge-list allfilesmergelist_abcd.txt --make-bed --out allautism_ABCDmerge
./plink --bfile SPARK_v1_forABCD --exclude allautism_ABCDmerge-merge.missnp --make-bed --out SPARK_v1_forABCD_v2 --threads 20 
./plink --bfile SPARK_v2_forABCD --exclude allautism_ABCDmerge-merge.missnp --make-bed --out SPARK_v2_forABCD_v2 --threads 20 
./plink --bfile SSC_forABCD --exclude allautism_ABCDmerge-merge.missnp --make-bed --out SSC_forABCD_v2 --threads 20 
./plink --bfile ABCD_euroonly_forABCD --exclude allautism_ABCDmerge-merge.missnp --make-bed --out ABCD_euroonly_forABCD_v2 --threads 20 

./plink --bfile ABCD_euroonly_forABCD_v2 --merge-list allfilesmergelist.txt --make-bed --out allautism_ABCDmerge --threads 20

#Generate frequency files
./plink --bfile ABCD_euroonly_forABCD_v2 --freq --out ukb_freq
./plink --bfile SPARK_v1_forABCD_v2 --freq --out SPARK_v1_freq
./plink --bfile SPARK_v2_forABCD_v2 --freq --out SPARK_v2_freq
./plink --bfile SSC_forABCD_v2 --freq --out SSC_v1_freq --nonfounders


#In R
setwd("~/Autism_heterogeneity/GRMs/GRMs_inc_UKB")
SSC = read.table("SSC_v1_freq.frq", fill = TRUE)
SPARK_v1 = fread("SPARK_v1_freq.frq")
SPARK_v2 = fread("SPARK_v2_freq.frq")
UKB = fread("ukb_freq.frq")

SSC = SSC[,c("V2", "V5")]
setnames(SSC, old = c("V2", "V5"), new = c("SNP", "MAF"))
SPARK_v1 = SPARK_v1[,c("SNP", "MAF")]
SPARK_v2 = SPARK_v2[,c("SNP", "MAF")]
UKB = UKB[,c("SNP", "MAF")]

merged = merge(SSC, SPARK_v1, by = "SNP")
merged = merge(merged, UKB, by = "SNP")

merged$MAF.x = as.numeric(as.character(merged$MAF.x))
merged$MAF_diff1 = merged$MAF - merged$MAF.x
merged$MAF_diff2 = merged$MAF - merged$MAF.y
merged$MAF_diff3 = merged$MAF.x - merged$MAF.y


merged2 = subset(merged, abs(merged$MAF_diff1) < 0.05)
merged2 = subset(merged2, abs(merged2$MAF_diff2) < 0.05)
merged2 = subset(merged2, abs(merged2$MAF_diff3) < 0.05)

write.table(merged2[,c("SNP")], file = "frequency_OK_allSNPs.txt", row.names = F, col.names = F, quote = F)


#plink commands
./plink --bfile allautism_ABCDmerge --extract frequency_OK_allSNPs.txt --hwe 0.000001 --geno 0.01 --maf 0.01 --make-bed --out allautism_ABCDmerge_SNPQC --threads 20

#Remove related individuals through GRM
./plink --bfile allautism_ABCDmerge_SNPQC --indep-pairwise 1000 1000 0.1 --out SNPsforpruning
./plink --bfile allautism_ABCDmerge_SNPQC --extract SNPsforpruning.prune.in --make-bed --out allautism_ABCDmerge_forGRMandPC
./gcta64  --bfile allautism_ABCDBmerge_forGRMandPC  --autosome  --make-grm  --out ABCD_autism_greml --thread-num 20
./gcta64 --grm ABCD_autism_greml --grm-cutoff 0.025 --make-grm --out ABCD_autism_greml_unrelated --thread-num 20 #43173 individuals

./plink --bfile allautism_ABCDmerge_forGRMandPC --keep ABCD_autism_greml_unrelated.grm.id --make-bed --out allautism_ABCDmerge_forGRMandPC_unrelated

./plink --bfile allautism_ABCDmerge_forGRMandPC_unrelated --pca --out allABCD_autism_pcs --threads 15

#In R

pcs = fread("allABCD_autism_pcs.eigenvec")
abcd = fread("ABCD_euroonly_forABCD_v2.fam")
spark_v1 = fread("SPARK_v1_forABCD_v2.fam")
spark_v2 = fread("SPARK_v2_forABCD_v2.fam")
ssc = fread("SSC_forABCD_v2.fam")

pcs$cohort = ifelse(pcs$V2 %in% abcd$V2, "ABCD", NA)
pcs$cohort = ifelse(pcs$V2 %in% spark_v1$V2, "SPARK", pcs$cohort)
pcs$cohort = ifelse(pcs$V2 %in% spark_v2$V2, "SPARK", pcs$cohort)
pcs$cohort = ifelse(pcs$V2 %in% ssc$V2, "SSC", pcs$cohort)

cases = subset(pcs, cohort == "SPARK" | cohort == "SSC")
controls = subset(pcs, cohort == "ABCD")


a = ggplot(pcs,aes(x=V3,y=V4,color=cohort)) + geom_point() + xlab("PC1") + ylab("PC2")

ggplot(pcs,aes(x=V3,y=V4,color=cohort)) + geom_point(alpha = 1/2) + xlab("PC1") + ylab("PC2")
ggplot(controls,aes(x=V3,y=V5,color=cohort)) + geom_point() + xlab("PC1") + ylab("PC2")


a = ggplot(pcs,aes(x=V3,y=V4,color=cohort)) + geom_point(alpha = 1/5) + xlab("PC1") + ylab("PC2") + theme_classic()
b = ggplot(pcs,aes(x=V3,y=V5,color=cohort)) + geom_point(alpha = 1/5) + xlab("PC1") + ylab("PC3") + theme_classic()
c = ggplot(pcs,aes(x=V4,y=V5,color=cohort)) + geom_point(alpha = 1/5) + xlab("PC2") + ylab("PC3") + theme_classic()
d = ggplot(pcs,aes(x=V3,y=V6,color=cohort)) + geom_point(alpha = 1/5) + xlab("PC1") + ylab("PC4") + theme_classic()
e = ggplot(pcs,aes(x=V4,y=V6,color=cohort)) + geom_point(alpha = 1/5) + xlab("PC2") + ylab("PC4") + theme_classic()
f = ggplot(pcs,aes(x=V5,y=V6,color=cohort)) + geom_point(alpha = 1/5) + xlab("PC3") + ylab("PC4") + theme_classic()


multiplot(a,b,d,c,e,f, cols = 2)

# for ldms
./gcta64  --bfile allautism_ABCDmerge_SNPQC --keep ABCD_autism_greml_unrelated.grm.id --ld-score-region 200  --out ABCD_autism_ldscore_byregion --thread-num 15

lds_seg = read.table("ABCD_autism_ldscore_byregion.score.ld",header=T,colClasses=c("character",rep("numeric",8)))
quartiles=summary(lds_seg$ldscore_SNP)

lb1 = which(lds_seg$ldscore_SNP <= quartiles[2])
lb2 = which(lds_seg$ldscore_SNP > quartiles[2] & lds_seg$ldscore_SNP <= quartiles[3])
lb3 = which(lds_seg$ldscore_SNP > quartiles[3] & lds_seg$ldscore_SNP <= quartiles[5])
lb4 = which(lds_seg$ldscore_SNP > quartiles[5])

lb1_snp = lds_seg$SNP[lb1]
lb2_snp = lds_seg$SNP[lb2]
lb3_snp = lds_seg$SNP[lb3]
lb4_snp = lds_seg$SNP[lb4]

write.table(lb1_snp, "snp_group1_ABCD.txt", row.names=F, quote=F, col.names=F)
write.table(lb2_snp, "snp_group2_ABCD.txt", row.names=F, quote=F, col.names=F)
write.table(lb3_snp, "snp_group3_ABCD.txt", row.names=F, quote=F, col.names=F)
write.table(lb4_snp, "snp_group4_ABCD.txt", row.names=F, quote=F, col.names=F)

./gcta64 --bfile allautism_ABCDmerge_SNPQC --extract snp_group1_ABCD.txt --keep ABCD_autism_greml_unrelated.grm.id --make-grm --out ABCD_autism_greml_snpgroup1 --thread-num 20
./gcta64 --bfile allautism_ABCDmerge_SNPQC --extract snp_group2_ABCD.txt --keep ABCD_autism_greml_unrelated.grm.id --make-grm --out ABCD_autism_greml_snpgroup2 --thread-num 20
./gcta64 --bfile allautism_ABCDmerge_SNPQC --extract snp_group3_ABCD.txt --keep ABCD_autism_greml_unrelated.grm.id --make-grm --out ABCD_autism_greml_snpgroup3 --thread-num 20
./gcta64 --bfile allautism_ABCDmerge_SNPQC --extract snp_group4_ABCD.txt --keep ABCD_autism_greml_unrelated.grm.id --make-grm --out ABCD_autism_greml_snpgroup4 --thread-num 20
./gcta64 --bfile allautism_ABCDmerge_forGRMandPC_unrelated   --make-grm --out ABCD_autism_ldpruned --thread-num 20
./gcta64 --bfile allautism_ABCDmerge_SNPQC --extract ~/Autism_heterogeneity/Commonwellimputed_SNPs_ABCD_SSC_SPARK_highr2.txt  --make-grm --out ABCD_autism_ldpruned_highr2 --thread-num 20

#cases-only
./gcta64 --bfile allautism_ABCDmerge_SNPQC --keep cases_to_keep.txt --extract ~/Autism_heterogeneity/Commonwellimputed_SNPs_ABCD_SSC_SPARK_highr2.txt  --make-grm --out autismonly_ldpruned_highr2 --thread-num 20
./plink --bfile allautism_ABCDmerge_forGRMandPC_unrelated --keep cases_to_keep.txt --pca --out autismonly_pcs --threads 15


###

##########################################################################################
##Generate covariate and PC file
samples = fread("~/Autism_heterogeneity/GRMs/GRMs_inc_UKB/ABCD_autism_greml_unrelated.grm.id")

setwd("~/Autism_heterogeneity/GRMs/GRMs_inc_UKB/")

abcd = fread("ABCD_euroonly_forABCD_v2.fam")
spark_v1 = fread("SPARK_v1_forABCD_v2.fam")
spark_v2 = fread("SPARK_v2_forABCD_v2.fam")
ssc = fread("SSC_forABCD_v2.fam")

samples$cohort = ifelse(samples$V2 %in% abcd$V2, "abcd", NA)
samples$cohort = ifelse(samples$V2 %in% spark_v1$V2, "spark", samples$cohort)
samples$cohort = ifelse(samples$V2 %in% spark_v2$V2, "spark", samples$cohort)
samples$cohort = ifelse(samples$V2 %in% ssc$V2, "ssc", samples$cohort)

samples$case_control = ifelse(samples$cohort == "ssc" | samples$cohort =="spark", "cases", "controls")

cases = subset(samples, case_control == "cases")

write.table(cases[,c("V1", "V2")], file = "cases_to_keep.txt", row.names = F, col.names = T, quote = F)

controls = subset(samples, case_control == "controls")


###Read basicpheno_data
SNV = fread("~/Autism_heterogeneity/Sebat_data/denovo_SNVindel.txt")
SV = fread("~/Autism_heterogeneity/Sebat_data/denovo_SV.txt")

SNV_upper_decile = subset(SNV, gnomad_loeuf < 0.37)
SV_upper_decile = subset(SV, gnomad_loeuf_minimum < 0.37)
SNV_ssc = subset(SNV_upper_decile, cohort == "SSC")
SNV_spark = subset(SNV_upper_decile, cohort == "SPARK")

cases$denovo = ifelse(cases$V1 %in% SV_upper_decile$fid | cases$V1 %in% SNV_ssc$fid  | cases$V2 %in% SNV_spark$iid, 1, NA )
controls$denovo = 0
cases$nondenovo = ifelse(cases$case_control =="cases" & is.na(cases$denovo) == "TRUE", 1, 0)
controls$nondenovo = 0

###Read phenotype_data
load("~/Autism_heterogeneity/SSC_dataforanalysis.Rdata")
SSC = merged

SSC$IQ_score = ifelse(SSC$ssc_diagnosis_nonverbal_iq < 24, 1, NA)
SSC$IQ_score = ifelse(SSC$ssc_diagnosis_nonverbal_iq > 24 & SSC$ssc_diagnosis_nonverbal_iq  < 39,  2, SSC$IQ_score)
SSC$IQ_score = ifelse(SSC$ssc_diagnosis_nonverbal_iq > 38 & SSC$ssc_diagnosis_nonverbal_iq  < 54, 3, SSC$IQ_score)
SSC$IQ_score = ifelse(SSC$ssc_diagnosis_nonverbal_iq > 53 & SSC$ssc_diagnosis_nonverbal_iq  < 69, 4, SSC$IQ_score)
SSC$IQ_score = ifelse(SSC$ssc_diagnosis_nonverbal_iq > 69 & SSC$ssc_diagnosis_nonverbal_iq  < 80, 5, SSC$IQ_score)
SSC$IQ_score = ifelse(SSC$ssc_diagnosis_nonverbal_iq > 79 & SSC$ssc_diagnosis_nonverbal_iq  < 90, 6, SSC$IQ_score)
SSC$IQ_score = ifelse(SSC$ssc_diagnosis_nonverbal_iq > 89 & SSC$ssc_diagnosis_nonverbal_iq  < 110, 7, SSC$IQ_score)
SSC$IQ_score = ifelse(SSC$ssc_diagnosis_nonverbal_iq > 109 & SSC$ssc_diagnosis_nonverbal_iq  < 120, 8, SSC$IQ_score)
SSC$IQ_score = ifelse(SSC$ssc_diagnosis_nonverbal_iq > 119 & SSC$ssc_diagnosis_nonverbal_iq  < 130, 9, SSC$IQ_score)
SSC$IQ_score = ifelse(SSC$ssc_diagnosis_nonverbal_iq > 130, 10, SSC$IQ_score)

SSC$Sex = ifelse(SSC$Sex == 2, "Female", "Male")
SSC$age_at_ados = round((SSC$age_at_ados)/12)

setnames(SSC, old = c("Sex", "age_at_ados", "vineland_ii_composite_standard_score",
                      "rbs_r_overall_score", "scq_total", "dcdq_total"), new = c("sex", "age_at_eval_years",
                                                                                 "VABS_score", "rbs_score", "scq_score", "dcdq_score") )

ssc2 = SSC[,c("IID", "sex", "age_at_eval_years", "VABS_score", "rbs_score", "scq_score", "dcdq_score", "PA1", "PA2", "PA3", "PA4", "PA5", "PA6", "IQ_score")]


spark = fread("~/Autism_heterogeneity/SPARK_allphenotypes.txt", fill = TRUE)

pheno_all = rbind(ssc2, spark)
pheno_all$PA6 = -1*pheno_all$PA6

ID = subset(pheno_all, IQ_score < 5)
cases$ID = ifelse(cases$V2 %in% ID$IID, 1, NA)
controls$ID = 0


nonID = subset(pheno_all, IQ_score > 4)
cases$nonID = ifelse(cases$V2 %in% nonID$IID, 1, NA)
controls$nonID = 0

highIQ = subset(pheno_all, IQ_score > 6)
cases$highIQ = ifelse(cases$V2 %in% highIQ$IID, 1, NA)
controls$highIQ = 0

males = subset(pheno_all, sex == "Male")
cases$males = ifelse(cases$V2 %in% males$IID, 1, NA)
controls$males = 0

females = subset(pheno_all, sex == "Female")
cases$females = ifelse(cases$V2 %in% females$IID, 1, NA)
controls$females = 0

nonID_female = subset(nonID, sex == "Female")
cases$nonID_female = ifelse(cases$V2 %in% nonID_female$IID, 1, NA)
controls$nonID_female = 0

nonID_male = subset(nonID, sex == "Male")
cases$nonID_male = ifelse(cases$V2 %in% nonID_male$IID, 1, NA)
controls$nonID_male = 0

ID_female = subset(ID, sex == "Female")
cases$ID_female = ifelse(cases$V2 %in% ID_female$IID, 1, NA)
controls$ID_female = 0

ID_male = subset(ID, sex == "Male")
cases$ID_male = ifelse(cases$V2 %in% ID_male$IID, 1, NA)
controls$ID_male = 0

high_scq = subset(pheno_all, scq_score > 29) #mean + 1 SD
cases$high_scq = ifelse(cases$V2 %in% high_scq$IID, 1, NA)
controls$high_scq = 0

high_rbs_r = subset(pheno_all, rbs_score > 50) #median + 1 MAD
cases$high_rbs_r = ifelse(cases$V2 %in% high_rbs_r$IID, 1, NA)
controls$high_rbs_r = 0

high_dcdq = subset(pheno_all, dcdq_score > 50) #median + 1 MAD
cases$high_dcdq = ifelse(cases$V2 %in% high_dcdq$IID, 1, NA)
controls$high_dcdq = 0

high_PA1 = subset(pheno_all, PA1 > 1)
cases$high_PA1 = ifelse(cases$V2 %in% high_PA1$IID, 1, NA)
controls$high_PA1 = 0

high_PA2 = subset(pheno_all, PA2 > 1)
cases$high_PA2 = ifelse(cases$V2 %in% high_PA2$IID, 1, NA)
controls$high_PA2 = 0

high_PA3 = subset(pheno_all, PA3 > 1)
cases$high_PA3 = ifelse(cases$V2 %in% high_PA3$IID, 1, NA)
controls$high_PA3 = 0

high_PA4 = subset(pheno_all, PA4 > 1)
cases$high_PA4 = ifelse(cases$V2 %in% high_PA4$IID, 1, NA)
controls$high_PA4 = 0

high_PA5 = subset(pheno_all, PA5 > 1)
cases$high_PA5 = ifelse(cases$V2 %in% high_PA5$IID, 1, NA)
controls$high_PA5 = 0

high_PA6 = subset(pheno_all, PA6 > 1)
cases$high_PA6 = ifelse(cases$V2 %in% high_PA6$IID, 1, NA)
controls$high_PA6 = 0

higher_PA1 = subset(pheno_all, PA1 > PA2  )
cases$higher_PA1 = ifelse(cases$V2 %in% higher_PA1$IID, 1, NA)
controls$higher_PA1 = 0

higher_PA2 = subset(pheno_all, PA2>  PA1 )
cases$higher_PA2 = ifelse(cases$V2 %in% higher_PA2$IID, 1, NA)
controls$higher_PA2 = 0

highest_PA1 = subset(pheno_all, PA1 - PA2 > 1.5 )
cases$highest_PA1 = ifelse(cases$V2 %in% highest_PA1$IID, 1, NA)
controls$highest_PA1 = 0

highest_PA2 = subset(pheno_all, PA2 - PA1 > 1.5 )
cases$highest_PA2 = ifelse(cases$V2 %in% highest_PA2$IID, 1, NA)
controls$highest_PA2 = 0

highestIQ = subset(pheno_all, IQ_score > 8)
cases$highestIQ = ifelse(cases$V2 %in% highestIQ$IID, 1, NA)
controls$highestIQ = 0

high_PA1_PA2 = subset(pheno_all, PA1 > 1 & PA2 > 1  )
cases$high_PA1_PA2 = ifelse(cases$V2 %in% high_PA1_PA2$IID, 1, NA)
controls$high_PA1_PA2 = 0


### Generate covariates

pcs = fread("allABCD_autism_pcs.eigenvec")
pcs = pcs[,c(2:22)]

age_cases = pheno_all[,c("IID", "age_at_eval_years")]
setnames(age_cases, "age_at_eval_years", "Age")
#######
age_sex = fread("~/Autism_heterogeneity/GRMs/GRMs_inc_UKB/abcd_basic_demo.txt")
age_sex$Age = age_sex$`Age in months at the time of the interview/test/sampling/imaging.`/12
age_sex$sex = ifelse(age_sex$`Sex of the subject` == "M", "Male", "Female")
age_controls = age_sex[,c("ID1", "Age")]
setnames(age_controls, "ID1", "IID")

age_all = rbind(age_cases, age_controls)
qcovar_all = merge(age_all, pcs, by.x = "IID", by.y = "V2")

cases_discrete_covar = pheno_all[,c("IID", "sex")]
controls_discrete_covar = age_sex[,c("ID1", "sex")]
setnames(controls_discrete_covar, "ID1", "IID")

covar_discrete_all = rbind(cases_discrete_covar, controls_discrete_covar)


all = rbind(cases, controls)
all_data = all[,c("V1", "V2")]
qcovar_all_2 = merge(all_data, qcovar_all, by.x = "V2", by.y = "IID")
covar_all_2 = merge(qcovar_all_2, covar_discrete_all, by.x = "V2", by.y = "IID")
merged_all = merge(all_data, covar_all_2, by = "V2")
all = all[all$V2 %in% merged_all$V2,]
setnames(merged_all, "V2", "IID")
setnames(merged_all, "V1.x", "FID")

setnames(all, "V2", "IID")
setnames(all, "V1", "FID")

merged_all = merged_all[,-c("V1.y", "FID")]

merged_all2 = merge(all, merged_all, by = "IID")

merged_all2$males = ifelse(merged_all2$sex ==  "Male", merged_all2$males, NA)
merged_all2$nonID_male = ifelse(merged_all2$sex ==  "Male", merged_all2$nonID_male, NA)
merged_all2$ID_male = ifelse(merged_all2$sex ==  "Male", merged_all2$ID_male, NA)

merged_all2$females = ifelse(merged_all2$sex ==  "Female", merged_all2$females , NA)
merged_all2$nonID_female = ifelse(merged_all2$sex ==  "Female", merged_all2$nonID_female , NA)
merged_all2$ID_female = ifelse(merged_all2$sex ==  "Female", merged_all2$ID_female , NA)

merged_all2$case_control = ifelse(merged_all2$case_control == "cases", 1, 0)
merged_all2 = unique(merged_all2)

write.table(merged_all2[,c("FID", "IID", "Age", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10","V11")], file = "qcovar_ABCD.txt", row.names = F, col.names = T, quote = F)
write.table(merged_all2[,c("FID", "IID", "sex")], file = "covar_ABCD.txt", row.names = F, col.names = T, quote = F)


###Sample to keep the case:control ratio under 1

cases = subset(merged_all2, case_control == "1")
cases2 = sample_n(cases, 4481)
merged_all2$case_control_downsampled = ifelse(merged_all2$case_control == "1" & !merged_all2$IID %in% cases2$IID, NA, merged_all2$case_control)

males = subset(cases, males == "1")
males2 = sample_n(males, 2386)
merged_all2$males_downsampled = ifelse(merged_all2$males == "1" & !merged_all2$IID %in% males2$IID, NA, merged_all2$males)
males2 = sample_n(males, 2095)
merged_all2$males_femalesampled = ifelse(merged_all2$males == "1" & !merged_all2$IID %in% males2$IID, NA, merged_all2$males)

males_nonid = subset(males, nonID_male == "1")
males2 = sample_n(males_nonid, 2386)
merged_all2$nonID_male_downsampled = ifelse(merged_all2$nonID_male == "1" & !merged_all2$IID %in% males2$IID, NA, merged_all2$nonID_male)


females = subset(cases, females == "1")
females2 = sample_n(females, 2095)
merged_all2$females_downampled = ifelse(merged_all2$females == "1" & !merged_all2$IID %in% females2$IID, NA, merged_all2$females)


write.table(merged_all2[,c("FID", "IID", "case_control", "denovo", "nondenovo", "ID", "nonID",
                   "highIQ", "males", "females", "nonID_female", "nonID_male", "ID_female", "ID_male",
                   "high_scq", "high_rbs_r", "high_dcdq", "high_PA1", "high_PA2", "high_PA3", "high_PA4",
                   "high_PA5", "high_PA6", "higher_PA1", "higher_PA2", "highest_PA1", "highest_PA2", 
                   "highestIQ",  "high_PA1_PA2", "case_control_downsampled", "males_downsampled", "males_femalesampled", 
                   "nonID_male_downsampled", "females_downampled")], file = "pheno_ABCD.txt", row.names = F, col.names = T, quote = F)

IDs = merged_all2[,c("IID", "FID")]

qpheno = merge(IDs, pheno_all, by = "IID")
write.table(qpheno[,c("FID", "IID", "IQ_score", "VABS_score", "dcdq_score", "rbs_score", "scq_score", "PA1", "PA2", "PA3", "PA4", "PA5", "PA6")], file = "qpheno_ABCD.txt", row.names = F, col.names = T, quote = F)



##de novo vs autism
data1 = fread("pheno_ABCD.txt")
spark_v1 = fread("~/SPARK/SPARK_postimputation/CEU/Plink_files/SPARKCEU_hg19_allchrs.fam")
ssc = fread("~/SFARI/liftOverPlink/files_imputed/SFARImergedall.fam")
denov = subset(data1, denovo == 1)


data1$denovo_v_cases = ifelse(data1$IID %in% spark_v1$V2, 0, NA)
data1$denovo_v_cases = ifelse(data1$IID %in% ssc$V2, 0, data1$denovo_v_cases)
data1$denovo_v_cases = ifelse(data1$IID %in% denov$IID, 1, data1$denovo_v_cases)

write.table(data1, file = "pheno_ABCD.txt", row.names = F, col.names = T, quote = F)


### Run_the_analyses###

### SNP heritability ##

./gcta64 --reml --grm ABCD_autism_ldpruned_highr2 --pheno pheno_ABCD.txt --covar covar_ABCD.txt --qcovar qcovar_ABCD.txt --thread-num 20 --mpheno 1 --out ./Results_ABCD/Case_control --prevalence 0.01
./gcta64 --reml --grm ABCD_autism_ldpruned_highr2 --pheno pheno_ABCD.txt --covar covar_ABCD.txt --qcovar qcovar_ABCD.txt --thread-num 15 --mpheno 2 --out ./Results_ABCD/denovo --prevalence 0.002
./gcta64 --reml --grm ABCD_autism_ldpruned_highr2 --pheno pheno_ABCD.txt --covar covar_ABCD.txt --qcovar qcovar_ABCD.txt --thread-num 15 --mpheno 3 --out ./Results_ABCD/nondenovo_highr2 --prevalence 0.01
./gcta64 --reml --grm ABCD_autism_ldpruned_highr2 --pheno pheno_ABCD.txt --covar covar_ABCD.txt --qcovar qcovar_ABCD.txt --thread-num 15 --mpheno 4 --out ./Results_ABCD/ID_highr2 --prevalence 0.01
./gcta64 --reml --grm ABCD_autism_ldpruned_highr2 --pheno pheno_ABCD.txt --covar covar_ABCD.txt --qcovar qcovar_ABCD.txt --thread-num 15 --mpheno 5 --out ./Results_ABCD/nonID_highr2 --prevalence 0.01
./gcta64 --reml --grm ABCD_autism_ldpruned_highr2 --pheno pheno_ABCD.txt --covar covar_ABCD.txt --qcovar qcovar_ABCD.txt --thread-num 15 --mpheno 6 --out ./Results_ABCD/highIQ_highr2 --prevalence 0.01
./gcta64 --reml --grm ABCD_autism_ldpruned_highr2 --pheno pheno_ABCD.txt --qcovar qcovar_ABCD.txt --thread-num 15 --mpheno 7 --out ./Results_ABCD/males_highr2 --prevalence 0.0288
./gcta64 --reml --grm ABCD_autism_ldpruned_highr2 --pheno pheno_ABCD.txt --qcovar qcovar_ABCD.txt --thread-num 15 --mpheno 8 --out ./Results_ABCD/females_highr2 --prevalence 0.0072
./gcta64 --reml --grm ABCD_autism_ldpruned_highr2 --pheno pheno_ABCD.txt --qcovar qcovar_ABCD.txt --thread-num 15 --mpheno 9 --out ./Results_ABCD/nonID_females_highr2 --prevalence 0.004392
./gcta64 --reml --grm ABCD_autism_ldpruned_highr2 --pheno pheno_ABCD.txt --qcovar qcovar_ABCD.txt --thread-num 15 --mpheno 10 --out ./Results_ABCD/nonID_males_highr2 --prevalence 0.01188
./gcta64 --reml --grm ABCD_autism_ldpruned_highr2 --pheno pheno_ABCD.txt --qcovar qcovar_ABCD.txt --thread-num 15 --mpheno 11 --out ./Results_ABCD/ID_female_highr2 --prevalence 0.002808
./gcta64 --reml --grm ABCD_autism_ldpruned_highr2 --pheno pheno_ABCD.txt --qcovar qcovar_ABCD.txt --thread-num 15 --mpheno 12 --out ./Results_ABCD/ID_male_highr2 --prevalence 0.009216
./gcta64 --reml --grm ABCD_autism_ldpruned_highr2 --pheno pheno_ABCD.txt --covar covar_ABCD.txt --qcovar qcovar_ABCD.txt --thread-num 15 --mpheno 13 --out ./Results_ABCD/high_scq_highr2 --prevalence 0.01
./gcta64 --reml --grm ABCD_autism_ldpruned_highr2 --pheno pheno_ABCD.txt --covar covar_ABCD.txt --qcovar qcovar_ABCD.txt --thread-num 15 --mpheno 14 --out ./Results_ABCD/high_rbs_highr2 --prevalence 0.01
./gcta64 --reml --grm ABCD_autism_ldpruned_highr2 --pheno pheno_ABCD.txt --covar covar_ABCD.txt --qcovar qcovar_ABCD.txt --thread-num 15 --mpheno 15 --out ./Results_ABCD/high_dcdq_highr2 --prevalence 0.01
./gcta64 --reml --grm ABCD_autism_ldpruned_highr2 --pheno pheno_ABCD.txt --covar covar_ABCD.txt --qcovar qcovar_ABCD.txt --thread-num 15 --mpheno 16 --out ./Results_ABCD/high_PA1_highr2 --prevalence 0.01
./gcta64 --reml --grm ABCD_autism_ldpruned_highr2 --pheno pheno_ABCD.txt --covar covar_ABCD.txt --qcovar qcovar_ABCD.txt --thread-num 15 --mpheno 17 --out ./Results_ABCD/high_PA2_highr2 --prevalence 0.01
./gcta64 --reml --grm ABCD_autism_ldpruned_highr2 --pheno pheno_ABCD.txt --covar covar_ABCD.txt --qcovar qcovar_ABCD.txt --thread-num 15 --mpheno 18 --out ./Results_ABCD/high_PA3_highr2 --prevalence 0.01
./gcta64 --reml --grm ABCD_autism_ldpruned_highr2 --pheno pheno_ABCD.txt --covar covar_ABCD.txt --qcovar qcovar_ABCD.txt --thread-num 15 --mpheno 19 --out ./Results_ABCD/high_PA4_highr2 --prevalence 0.01
./gcta64 --reml --grm ABCD_autism_ldpruned_highr2 --pheno pheno_ABCD.txt --covar covar_ABCD.txt --qcovar qcovar_ABCD.txt --thread-num 15 --mpheno 20 --out ./Results_ABCD/high_PA5_highr2 --prevalence 0.01
./gcta64 --reml --grm ABCD_autism_ldpruned_highr2 --pheno pheno_ABCD.txt --covar covar_ABCD.txt --qcovar qcovar_ABCD.txt --thread-num 15 --mpheno 21 --out ./Results_ABCD/high_PA6_highr2 --prevalence 0.01
./gcta64 --reml --grm ABCD_autism_ldpruned_highr2 --pheno pheno_ABCD.txt --covar covar_ABCD.txt --qcovar qcovar_ABCD.txt --thread-num 15 --mpheno 22 --out ./Results_ABCD/higher_PA1_highr2 --prevalence 0.01
./gcta64 --reml --grm ABCD_autism_ldpruned_highr2 --pheno pheno_ABCD.txt --covar covar_ABCD.txt --qcovar qcovar_ABCD.txt --thread-num 15 --mpheno 23 --out ./Results_ABCD/higher_PA2_highr2 --prevalence 0.01
./gcta64 --reml --grm ABCD_autism_ldpruned_highr2 --pheno pheno_ABCD.txt --covar covar_ABCD.txt --qcovar qcovar_ABCD.txt --thread-num 15 --mpheno 26 --out ./Results_ABCD/highest_IQ --prevalence 0.00288
./gcta64 --reml --grm ABCD_autism_ldpruned_highr2 --pheno pheno_ABCD.txt --covar covar_ABCD.txt --qcovar qcovar_ABCD.txt --thread-num 15 --mpheno 27 --out ./Results_ABCD/highestPA1PA2 --prevalence 0.004

./gcta64 --reml --grm ABCD_autism_ldpruned_highr2 --pheno pheno_ABCD.txt --covar covar_ABCD.txt --qcovar qcovar_ABCD.txt --thread-num 15 --mpheno 28 --out ./Results_ABCD/casecontrol_2 --prevalence 0.01
./gcta64 --reml --grm ABCD_autism_ldpruned_highr2 --pheno pheno_ABCD.txt --qcovar qcovar_ABCD.txt --thread-num 15 --mpheno 29 --out ./Results_ABCD/males_downsampled --prevalence 0.0288
./gcta64 --reml --grm ABCD_autism_ldpruned_highr2 --pheno pheno_ABCD.txt --qcovar qcovar_ABCD.txt --thread-num 15 --mpheno 30 --out ./Results_ABCD/males_femalesampled --prevalence 0.0288
./gcta64 --reml --grm ABCD_autism_ldpruned_highr2 --pheno pheno_ABCD.txt --qcovar qcovar_ABCD.txt --thread-num 15 --mpheno 31 --out ./Results_ABCD/nonID_male_downsampled --prevalence 0.019584
./gcta64 --reml --grm ABCD_autism_ldpruned_highr2 --pheno pheno_ABCD.txt --qcovar qcovar_ABCD.txt --thread-num 15 --mpheno 32 --out ./Results_ABCD/females_downsampled --prevalence 0.0072
./gcta64 --reml --grm ABCD_autism_ldpruned_highr2 --pheno pheno_ABCD.txt --qcovar qcovar_ABCD.txt --thread-num 15 --mpheno 33 --out ./Results_ABCD/denovo_vs_cases --prevalence 0.12





./gcta64 --reml --grm autismonly_ldpruned_highr2 --pheno qpheno_ABCD.txt --covar covar_ABCD_withcohort.txt --qcovar qcovar_ABCD_casesonly.txt --thread-num 15 --mpheno 1 --out ./Results_ABCD/IQ
./gcta64 --reml --grm autismonly_ldpruned_highr2 --pheno qpheno_ABCD.txt --covar covar_ABCD_withcohort.txt --qcovar qcovar_ABCD_casesonly.txt --thread-num 15 --mpheno 2 --out ./Results_ABCD/vabs
./gcta64 --reml --grm autismonly_ldpruned_highr2 --pheno qpheno_ABCD.txt --covar covar_ABCD_withcohort.txt --qcovar qcovar_ABCD_casesonly.txt --thread-num 15 --mpheno 3 --out ./Results_ABCD/dcdq
./gcta64 --reml --grm autismonly_ldpruned_highr2 --pheno qpheno_ABCD.txt --covar covar_ABCD_withcohort.txt --qcovar qcovar_ABCD_casesonly.txt --thread-num 15 --mpheno 4 --out ./Results_ABCD/rbs
./gcta64 --reml --grm autismonly_ldpruned_highr2 --pheno qpheno_ABCD.txt --covar covar_ABCD_withcohort.txt --qcovar qcovar_ABCD_casesonly.txt --thread-num 15 --mpheno 5 --out ./Results_ABCD/scq
./gcta64 --reml --grm autismonly_ldpruned_highr2 --pheno qpheno_ABCD.txt --covar covar_ABCD_withcohort.txt --qcovar qcovar_ABCD_casesonly.txt --thread-num 15 --mpheno 6 --out ./Results_ABCD/PA1
./gcta64 --reml --grm autismonly_ldpruned_highr2 --pheno qpheno_ABCD.txt --covar covar_ABCD_withcohort.txt --qcovar qcovar_ABCD_casesonly.txt --thread-num 15 --mpheno 7 --out ./Results_ABCD/PA2
./gcta64 --reml --grm autismonly_ldpruned_highr2 --pheno qpheno_ABCD.txt --covar covar_ABCD_withcohort.txt --qcovar qcovar_ABCD_casesonly.txt --thread-num 15 --mpheno 8 --out ./Results_ABCD/PA3
./gcta64 --reml --grm autismonly_ldpruned_highr2 --pheno qpheno_ABCD.txt --covar covar_ABCD_withcohort.txt --qcovar qcovar_ABCD_casesonly.txt --thread-num 15 --mpheno 9 --out ./Results_ABCD/PA4
./gcta64 --reml --grm autismonly_ldpruned_highr2 --pheno qpheno_ABCD.txt --covar covar_ABCD_withcohort.txt --qcovar qcovar_ABCD_casesonly.txt --thread-num 15 --mpheno 10 --out ./Results_ABCD/PA5
./gcta64 --reml --grm autismonly_ldpruned_highr2 --pheno qpheno_ABCD.txt --covar covar_ABCD_withcohort.txt --qcovar qcovar_ABCD_casesonly.txt --thread-num 15 --mpheno 11 --out ./Results_ABCD/PA6



###qpheno gencor

./gcta64 --reml-bivar 1 2 --grm autismonly_ldpruned_highr2 --pheno qpheno_ABCD.txt --covar covar_ABCD_withcohort.txt --qcovar qcovar_ABCD_casesonly.txt --thread-num 15 --out ./Results_ABCD/IQ_vabs
./gcta64 --reml-bivar 1 3 --grm autismonly_ldpruned_highr2 --pheno qpheno_ABCD.txt --covar covar_ABCD_withcohort.txt --qcovar qcovar_ABCD_casesonly.txt --thread-num 15 --out ./Results_ABCD/IQ_dcdq
./gcta64 --reml-bivar 1 4 --grm autismonly_ldpruned_highr2 --pheno qpheno_ABCD.txt --covar covar_ABCD_withcohort.txt --qcovar qcovar_ABCD_casesonly.txt --thread-num 15 --out ./Results_ABCD/IQ_rbs
./gcta64 --reml-bivar 1 5 --grm autismonly_ldpruned_highr2 --pheno qpheno_ABCD.txt --covar covar_ABCD_withcohort.txt --qcovar qcovar_ABCD_casesonly.txt --thread-num 15 --out ./Results_ABCD/IQ_scq
./gcta64 --reml-bivar 1 6 --grm autismonly_ldpruned_highr2 --pheno qpheno_ABCD.txt --covar covar_ABCD_withcohort.txt --qcovar qcovar_ABCD_casesonly.txt --thread-num 15 --out ./Results_ABCD/IQ_PA1
./gcta64 --reml-bivar 1 7 --grm autismonly_ldpruned_highr2 --pheno qpheno_ABCD.txt --covar covar_ABCD_withcohort.txt --qcovar qcovar_ABCD_casesonly.txt --thread-num 15 --out ./Results_ABCD/IQ_PA2
./gcta64 --reml-bivar 1 8 --grm autismonly_ldpruned_highr2 --pheno qpheno_ABCD.txt --covar covar_ABCD_withcohort.txt --qcovar qcovar_ABCD_casesonly.txt --thread-num 15 --out ./Results_ABCD/IQ_PA3
./gcta64 --reml-bivar 1 9 --grm autismonly_ldpruned_highr2 --pheno qpheno_ABCD.txt --covar covar_ABCD_withcohort.txt --qcovar qcovar_ABCD_casesonly.txt --thread-num 15 --out ./Results_ABCD/IQ_Pa4
./gcta64 --reml-bivar 1 10 --grm autismonly_ldpruned_highr2 --pheno qpheno_ABCD.txt --covar covar_ABCD_withcohort.txt --qcovar qcovar_ABCD_casesonly.txt --thread-num 15 --out ./Results_ABCD/IQ_pa5
./gcta64 --reml-bivar 1 11 --grm autismonly_ldpruned_highr2 --pheno qpheno_ABCD.txt --covar covar_ABCD_withcohort.txt --qcovar qcovar_ABCD_casesonly.txt --thread-num 15 --out ./Results_ABCD/IQ_PA6


./gcta64 --reml-bivar 4 2 --grm autismonly_ldpruned_highr2 --pheno qpheno_ABCD.txt --covar covar_ABCD_withcohort.txt --qcovar qcovar_ABCD_casesonly.txt --thread-num 15 --out ./Results_ABCD/rbs_vabs
./gcta64 --reml-bivar 4 3 --grm autismonly_ldpruned_highr2 --pheno qpheno_ABCD.txt --covar covar_ABCD_withcohort.txt --qcovar qcovar_ABCD_casesonly.txt --thread-num 15 --out ./Results_ABCD/rbs_dcdq
./gcta64 --reml-bivar 4 5 --grm autismonly_ldpruned_highr2 --pheno qpheno_ABCD.txt --covar covar_ABCD_withcohort.txt --qcovar qcovar_ABCD_casesonly.txt --thread-num 15 --out ./Results_ABCD/rbs_scq
./gcta64 --reml-bivar 4 6 --grm autismonly_ldpruned_highr2 --pheno qpheno_ABCD.txt --covar covar_ABCD_withcohort.txt --qcovar qcovar_ABCD_casesonly.txt --thread-num 15 --out ./Results_ABCD/rbs_PA1
./gcta64 --reml-bivar 4 7 --grm autismonly_ldpruned_highr2 --pheno qpheno_ABCD.txt --covar covar_ABCD_withcohort.txt --qcovar qcovar_ABCD_casesonly.txt --thread-num 15 --out ./Results_ABCD/rbs_PA2
./gcta64 --reml-bivar 4 8 --grm autismonly_ldpruned_highr2 --pheno qpheno_ABCD.txt --covar covar_ABCD_withcohort.txt --qcovar qcovar_ABCD_casesonly.txt --thread-num 15 --out ./Results_ABCD/rbs_PA3
./gcta64 --reml-bivar 4 9 --grm autismonly_ldpruned_highr2 --pheno qpheno_ABCD.txt --covar covar_ABCD_withcohort.txt --qcovar qcovar_ABCD_casesonly.txt --thread-num 15 --out ./Results_ABCD/rbs_Pa4
./gcta64 --reml-bivar 4 10 --grm autismonly_ldpruned_highr2 --pheno qpheno_ABCD.txt --covar covar_ABCD_withcohort.txt --qcovar qcovar_ABCD_casesonly.txt --thread-num 15 --out ./Results_ABCD/rbs_pa5
./gcta64 --reml-bivar 4 11 --grm autismonly_ldpruned_highr2 --pheno qpheno_ABCD.txt --covar covar_ABCD_withcohort.txt --qcovar qcovar_ABCD_casesonly.txt --thread-num 15 --out ./Results_ABCD/rbs_PA6

./gcta64 --reml-bivar 5 2 --grm autismonly_ldpruned_highr2 --pheno qpheno_ABCD.txt --covar covar_ABCD_withcohort.txt --qcovar qcovar_ABCD_casesonly.txt --thread-num 15 --out ./Results_ABCD/scq_vabs
./gcta64 --reml-bivar 5 3 --grm autismonly_ldpruned_highr2 --pheno qpheno_ABCD.txt --covar covar_ABCD_withcohort.txt --qcovar qcovar_ABCD_casesonly.txt --thread-num 15 --out ./Results_ABCD/scq_dcdq
./gcta64 --reml-bivar 5 6 --grm autismonly_ldpruned_highr2 --pheno qpheno_ABCD.txt --covar covar_ABCD_withcohort.txt --qcovar qcovar_ABCD_casesonly.txt --thread-num 15 --out ./Results_ABCD/scq_PA1
./gcta64 --reml-bivar 5 7 --grm autismonly_ldpruned_highr2 --pheno qpheno_ABCD.txt --covar covar_ABCD_withcohort.txt --qcovar qcovar_ABCD_casesonly.txt --thread-num 15 --out ./Results_ABCD/scq_PA2
./gcta64 --reml-bivar 5 8 --grm autismonly_ldpruned_highr2 --pheno qpheno_ABCD.txt --covar covar_ABCD_withcohort.txt --qcovar qcovar_ABCD_casesonly.txt --thread-num 15 --out ./Results_ABCD/scq_PA3
./gcta64 --reml-bivar 5 9 --grm autismonly_ldpruned_highr2 --pheno qpheno_ABCD.txt --covar covar_ABCD_withcohort.txt --qcovar qcovar_ABCD_casesonly.txt --thread-num 15 --out ./Results_ABCD/scq_Pa4
./gcta64 --reml-bivar 5 10 --grm autismonly_ldpruned_highr2 --pheno qpheno_ABCD.txt --covar covar_ABCD_withcohort.txt --qcovar qcovar_ABCD_casesonly.txt --thread-num 15 --out ./Results_ABCD/scq_pa5
./gcta64 --reml-bivar 5 11 --grm autismonly_ldpruned_highr2 --pheno qpheno_ABCD.txt --covar covar_ABCD_withcohort.txt --qcovar qcovar_ABCD_casesonly.txt --thread-num 15 --out ./Results_ABCD/scq_PA6



./gcta64 --reml-bivar 6 7 --grm autismonly_ldpruned_highr2 --pheno qpheno_ABCD.txt --covar covar_ABCD_withcohort.txt --qcovar qcovar_ABCD_casesonly.txt --thread-num 15 --out ./Results_ABCD/PA1_PA2
./gcta64 --reml-bivar 6 8 --grm autismonly_ldpruned_highr2 --pheno qpheno_ABCD.txt --covar covar_ABCD_withcohort.txt --qcovar qcovar_ABCD_casesonly.txt --thread-num 15 --out ./Results_ABCD/PA1_PA3
./gcta64 --reml-bivar 6 9 --grm autismonly_ldpruned_highr2 --pheno qpheno_ABCD.txt --covar covar_ABCD_withcohort.txt --qcovar qcovar_ABCD_casesonly.txt --thread-num 15 --out ./Results_ABCD/PA1_PA4
./gcta64 --reml-bivar 6 10 --grm autismonly_ldpruned_highr2 --pheno qpheno_ABCD.txt --covar covar_ABCD_withcohort.txt --qcovar qcovar_ABCD_casesonly.txt --thread-num 15 --out ./Results_ABCD/PA1_PA5
./gcta64 --reml-bivar 6 11 --grm autismonly_ldpruned_highr2 --pheno qpheno_ABCD.txt --covar covar_ABCD_withcohort.txt --qcovar qcovar_ABCD_casesonly.txt --thread-num 15 --out ./Results_ABCD/PA1_PA6

./gcta64 --reml-bivar 7 8 --grm autismonly_ldpruned_highr2 --pheno qpheno_ABCD.txt --covar covar_ABCD_withcohort.txt --qcovar qcovar_ABCD_casesonly.txt --thread-num 15 --out ./Results_ABCD/PA2_PA3
./gcta64 --reml-bivar 7 9 --grm autismonly_ldpruned_highr2 --pheno qpheno_ABCD.txt --covar covar_ABCD_withcohort.txt --qcovar qcovar_ABCD_casesonly.txt --thread-num 15 --out ./Results_ABCD/PA2_PA4
./gcta64 --reml-bivar 7 10 --grm autismonly_ldpruned_highr2 --pheno qpheno_ABCD.txt --covar covar_ABCD_withcohort.txt --qcovar qcovar_ABCD_casesonly.txt --thread-num 15 --out ./Results_ABCD/PA2_PA5
./gcta64 --reml-bivar 7 11 --grm autismonly_ldpruned_highr2 --pheno qpheno_ABCD.txt --covar covar_ABCD_withcohort.txt --qcovar qcovar_ABCD_casesonly.txt --thread-num 15 --out ./Results_ABCD/PA2_PA6


./gcta64 --reml-bivar 8 9 --grm autismonly_ldpruned_highr2 --pheno qpheno_ABCD.txt --covar covar_ABCD_withcohort.txt --qcovar qcovar_ABCD_casesonly.txt --thread-num 15 --out ./Results_ABCD/PA3_PA4
./gcta64 --reml-bivar 8 10 --grm autismonly_ldpruned_highr2 --pheno qpheno_ABCD.txt --covar covar_ABCD_withcohort.txt --qcovar qcovar_ABCD_casesonly.txt --thread-num 15 --out ./Results_ABCD/PA3_PA5
./gcta64 --reml-bivar 8 11 --grm autismonly_ldpruned_highr2 --pheno qpheno_ABCD.txt --covar covar_ABCD_withcohort.txt --qcovar qcovar_ABCD_casesonly.txt --thread-num 15 --out ./Results_ABCD/PA3_PA6


./gcta64 --reml-bivar 9 10 --grm autismonly_ldpruned_highr2 --pheno qpheno_ABCD.txt --covar covar_ABCD_withcohort.txt --qcovar qcovar_ABCD_casesonly.txt --thread-num 15 --out ./Results_ABCD/PA4_PA5
./gcta64 --reml-bivar 9 11 --grm autismonly_ldpruned_highr2 --pheno qpheno_ABCD.txt --covar covar_ABCD_withcohort.txt --qcovar qcovar_ABCD_casesonly.txt --thread-num 15 --out ./Results_ABCD/PA4_PA6

./gcta64 --reml-bivar 10 11 --grm autismonly_ldpruned_highr2 --pheno qpheno_ABCD.txt --covar covar_ABCD_withcohort.txt --qcovar qcovar_ABCD_casesonly.txt --thread-num 15 --out ./Results_ABCD/PA5_PA6


###Genetic correlation###

./gcta64 --reml-bivar 1 2 --grm ABCD_autism_ldpruned --pheno pheno_ABCD.txt --thread-num 15 --covar covar_ABCD.txt --qcovar qcovar_ABCD.txt --out ./Results_ABCD/rg_all_nondenovo --reml-bivar-prevalence 0.018 0.001 --reml-bivar-nocove
./gcta64 --reml-bivar 1 3 --grm ABCD_autism_ldpruned --pheno pheno_ABCD.txt --thread-num 15 --covar covar_ABCD.txt --qcovar qcovar_ABCD.txt --out ./Results_ABCD/rg_all_denovo --reml-bivar-prevalence 0.018 0.001 --reml-bivar-nocove
./gcta64 --reml-bivar 1 4 --grm ABCD_autism_ldpruned --pheno pheno_ABCD.txt --thread-num 15 --covar covar_ABCD.txt --qcovar qcovar_ABCD.txt --out ./Results_ABCD/rg_denovo_nondenovo --reml-bivar-prevalence 0.018 0.001 --reml-bivar-nocove
./gcta64 --reml-bivar 1 5 --grm ABCD_autism_ldpruned --pheno pheno_ABCD.txt --thread-num 15 --covar covar_ABCD.txt --qcovar qcovar_ABCD.txt --out ./Results_ABCD/rg_denovo_nondenovo --reml-bivar-prevalence 0.018 0.001 --reml-bivar-nocove
./gcta64 --reml-bivar 1 6 --grm ABCD_autism_ldpruned --pheno pheno_ABCD.txt --thread-num 15 --covar covar_ABCD.txt --qcovar qcovar_ABCD.txt --out ./Results_ABCD/rg_denovo_nondenovo --reml-bivar-prevalence 0.018 0.001 --reml-bivar-nocove
./gcta64 --reml-bivar 1 7 --grm ABCD_autism_ldpruned --pheno pheno_ABCD.txt --thread-num 15 --covar covar_ABCD.txt --qcovar qcovar_ABCD.txt --out ./Results_ABCD/rg_denovo_nondenovo --reml-bivar-prevalence 0.018 0.001 --reml-bivar-nocove
./gcta64 --reml-bivar 1 8 --grm ABCD_autism_ldpruned --pheno pheno_ABCD.txt --thread-num 15 --covar covar_ABCD.txt --qcovar qcovar_ABCD.txt --out ./Results_ABCD/rg_denovo_nondenovo --reml-bivar-prevalence 0.018 0.001 --reml-bivar-nocove
./gcta64 --reml-bivar 1 9 --grm ABCD_autism_ldpruned --pheno pheno_ABCD.txt --thread-num 15 --covar covar_ABCD.txt --qcovar qcovar_ABCD.txt --out ./Results_ABCD/rg_denovo_nondenovo --reml-bivar-prevalence 0.018 0.001 --reml-bivar-nocove
./gcta64 --reml-bivar 1 10 --grm ABCD_autism_ldpruned --pheno pheno_ABCD.txt --thread-num 15 --covar covar_ABCD.txt --qcovar qcovar_ABCD.txt --out ./Results_ABCD/rg_denovo_nondenovo --reml-bivar-prevalence 0.018 0.001 --reml-bivar-nocove
./gcta64 --reml-bivar 1 11 --grm ABCD_autism_ldpruned --pheno pheno_ABCD.txt --thread-num 15 --covar covar_ABCD.txt --qcovar qcovar_ABCD.txt --out ./Results_ABCD/rg_denovo_nondenovo --reml-bivar-prevalence 0.018 0.001 --reml-bivar-nocove
./gcta64 --reml-bivar 1 12 --grm ABCD_autism_ldpruned --pheno pheno_ABCD.txt --thread-num 15 --covar covar_ABCD.txt --qcovar qcovar_ABCD.txt --out ./Results_ABCD/rg_denovo_nondenovo --reml-bivar-prevalence 0.018 0.001 --reml-bivar-nocove
./gcta64 --reml-bivar 1 13 --grm ABCD_autism_ldpruned --pheno pheno_ABCD.txt --thread-num 15 --covar covar_ABCD.txt --qcovar qcovar_ABCD.txt --out ./Results_ABCD/rg_denovo_nondenovo --reml-bivar-prevalence 0.018 0.001 --reml-bivar-nocove
./gcta64 --reml-bivar 1 14 --grm ABCD_autism_ldpruned --pheno pheno_ABCD.txt --thread-num 15 --covar covar_ABCD.txt --qcovar qcovar_ABCD.txt --out ./Results_ABCD/rg_denovo_nondenovo --reml-bivar-prevalence 0.018 0.001 --reml-bivar-nocove
./gcta64 --reml-bivar 1 15 --grm ABCD_autism_ldpruned --pheno pheno_ABCD.txt --thread-num 15 --covar covar_ABCD.txt --qcovar qcovar_ABCD.txt --out ./Results_ABCD/rg_denovo_nondenovo --reml-bivar-prevalence 0.018 0.001 --reml-bivar-nocove
./gcta64 --reml-bivar 1 16 --grm ABCD_autism_ldpruned --pheno pheno_ABCD.txt --thread-num 15 --covar covar_ABCD.txt --qcovar qcovar_ABCD.txt --out ./Results_ABCD/rg_denovo_nondenovo --reml-bivar-prevalence 0.018 0.001 --reml-bivar-nocove
./gcta64 --reml-bivar 1 17 --grm ABCD_autism_ldpruned --pheno pheno_ABCD.txt --thread-num 15 --covar covar_ABCD.txt --qcovar qcovar_ABCD.txt --out ./Results_ABCD/rg_denovo_nondenovo --reml-bivar-prevalence 0.018 0.001 --reml-bivar-nocove
./gcta64 --reml-bivar 1 18 --grm ABCD_autism_ldpruned --pheno pheno_ABCD.txt --thread-num 15 --covar covar_ABCD.txt --qcovar qcovar_ABCD.txt --out ./Results_ABCD/rg_denovo_nondenovo --reml-bivar-prevalence 0.018 0.001 --reml-bivar-nocove
./gcta64 --reml-bivar 1 20 --grm ABCD_autism_ldpruned --pheno pheno_ABCD.txt --thread-num 15 --covar covar_ABCD.txt --qcovar qcovar_ABCD.txt --out ./Results_ABCD/rg_denovo_nondenovo --reml-bivar-prevalence 0.018 0.001 --reml-bivar-nocove
./gcta64 --reml-bivar 1 21 --grm ABCD_autism_ldpruned --pheno pheno_ABCD.txt --thread-num 15 --covar covar_ABCD.txt --qcovar qcovar_ABCD.txt --out ./Results_ABCD/rg_denovo_nondenovo --reml-bivar-prevalence 0.018 0.001 --reml-bivar-nocove
./gcta64 --reml-bivar 1 22 --grm ABCD_autism_ldpruned --pheno pheno_ABCD.txt --thread-num 15 --covar covar_ABCD.txt --qcovar qcovar_ABCD.txt --out ./Results_ABCD/rg_denovo_nondenovo --reml-bivar-prevalence 0.018 0.001 --reml-bivar-nocove
./gcta64 --reml-bivar 1 19 --grm ABCD_autism_ldpruned --pheno pheno_ABCD.txt --thread-num 15 --covar covar_ABCD.txt --qcovar qcovar_ABCD.txt --out ./Results_ABCD/rg_denovo_nondenovo --reml-bivar-prevalence 0.018 0.001 --reml-bivar-nocove

./gcta64 --reml-bivar 2 3 --grm ABCD_autism_ldpruned --pheno pheno_ABCD.txt --thread-num 15 --covar covar_ABCD.txt --qcovar qcovar_ABCD.txt --out ./Results_ABCD/rg_denovo_nondenovo --reml-bivar-prevalence 0.001 0.001 --reml-bivar-nocove
./gcta64 --reml-bivar 9 10 --grm ABCD_autism_ldpruned --pheno pheno_ABCD.txt --thread-num 15 --covar covar_ABCD.txt --qcovar qcovar_ABCD.txt --out ./Results_ABCD/rg_nonidmales_nonidfemales --reml-bivar-prevalence 0.004392 0.019584 --reml-bivar-nocove










##PCGC analyses

#Create the files
./gcta64 --grm ABCD_autism_ldpruned_highr2 --keep covar_all_ABCD_LDAK_keepfile.txt --out ldakgrm1 --make-grm --thread-num 20
./ldak5.1.linux --adjust-grm ldakgrm2 --grm ldakgrm1 --covar covar_all_ABCD_LDAK.txt --kinship-details NO
./ldak5.1.linux --adjust-grm ldakgrm3 --grm ldakgrm1 --covar covar_all_withoutsex_ABCD_LDAK.txt --kinship-details NO

#Run REML
./ldak5.1.linux  --grm ldakgrm2 --covar covar_all_ABCD_LDAK.txt --pheno pheno_ABCD.txt --prevalence 0.018 --mpheno 1 --kinship-details NO --pcgc ./Results_ABCD/casecontrol_2_LDAK
./ldak5.1.linux  --grm ldakgrm2 --covar covar_all_ABCD_LDAK.txt --pheno pheno_ABCD.txt --prevalence 0.002 --mpheno 2 --kinship-details NO --pcgc ./Results_ABCD/denovo_2_LDAK
./ldak5.1.linux  --grm ldakgrm2 --covar covar_all_ABCD_LDAK.txt --pheno pheno_ABCD.txt --prevalence 0.00594 --mpheno 4 --kinship-details NO --pcgc ./Results_ABCD/ID_LDAK
./ldak5.1.linux  --grm ldakgrm2 --covar covar_all_ABCD_LDAK.txt --pheno pheno_ABCD.txt --prevalence 0.01188 --mpheno 5 --kinship-details NO --pcgc ./Results_ABCD/nonID_LDAK
./ldak5.1.linux  --grm ldakgrm2 --covar covar_all_ABCD_LDAK.txt --pheno pheno_ABCD.txt --prevalence 0.01188 --mpheno 6 --kinship-details NO --pcgc ./Results_ABCD/highIQ_LDAK
./ldak5.1.linux  --grm ldakgrm2 --covar covar_all_ABCD_LDAK.txt --pheno pheno_ABCD.txt --prevalence 0.00288 --mpheno 13 --kinship-details NO --pcgc ./Results_ABCD/high_scq_LDAK
./ldak5.1.linux  --grm ldakgrm2 --covar covar_all_ABCD_LDAK.txt --pheno pheno_ABCD.txt --prevalence 0.00288 --mpheno 14 --kinship-details NO --pcgc ./Results_ABCD/high_rbs_LDAK
./ldak5.1.linux  --grm ldakgrm2 --covar covar_all_ABCD_LDAK.txt --pheno pheno_ABCD.txt --prevalence 0.00288 --mpheno 15 --kinship-details NO --pcgc ./Results_ABCD/high_dcdq_LDAK
./ldak5.1.linux  --grm ldakgrm2 --covar covar_all_ABCD_LDAK.txt --pheno pheno_ABCD.txt --prevalence 0.00288 --mpheno 16 --kinship-details NO --pcgc ./Results_ABCD/high_PA1_LDAK
./ldak5.1.linux  --grm ldakgrm2 --covar covar_all_ABCD_LDAK.txt --pheno pheno_ABCD.txt --prevalence 0.00288 --mpheno 17 --kinship-details NO --pcgc ./Results_ABCD/high_PA2_LDAK
./ldak5.1.linux  --grm ldakgrm2 --covar covar_all_ABCD_LDAK.txt --pheno pheno_ABCD.txt --prevalence 0.00288 --mpheno 18 --kinship-details NO --pcgc ./Results_ABCD/high_PA3_LDAK
./ldak5.1.linux  --grm ldakgrm2 --covar covar_all_ABCD_LDAK.txt --pheno pheno_ABCD.txt --prevalence 0.00288 --mpheno 19 --kinship-details NO --pcgc ./Results_ABCD/high_PA4_LDAK
./ldak5.1.linux  --grm ldakgrm2 --covar covar_all_ABCD_LDAK.txt --pheno pheno_ABCD.txt --prevalence 0.00288 --mpheno 20 --kinship-details NO --pcgc ./Results_ABCD/high_PA5_LDAK
./ldak5.1.linux  --grm ldakgrm2 --covar covar_all_ABCD_LDAK.txt --pheno pheno_ABCD.txt --prevalence 0.00288 --mpheno 21 --kinship-details NO --pcgc ./Results_ABCD/high_PA6_LDAK
./ldak5.1.linux  --grm ldakgrm2 --covar covar_all_ABCD_LDAK.txt --pheno pheno_ABCD.txt --prevalence 0.009 --mpheno 22 --kinship-details NO --pcgc ./Results_ABCD/higher_PA1_LDAK
./ldak5.1.linux  --grm ldakgrm2 --covar covar_all_ABCD_LDAK.txt --pheno pheno_ABCD.txt --prevalence 0.009 --mpheno 23 --kinship-details NO --pcgc ./Results_ABCD/higher_PA2_LDAK
./ldak5.1.linux  --grm ldakgrm2 --covar covar_all_ABCD_LDAK.txt --pheno pheno_ABCD.txt --prevalence 0.007 --mpheno 26 --kinship-details NO --pcgc ./Results_ABCD/highest_IQ_LDAK

./ldak5.1.linux  --grm ldakgrm2 --covar covar_all_ABCD_LDAK.txt --pheno pheno_ABCD.txt --prevalence 0.12 --mpheno 33 --kinship-details NO --pcgc ./Results_ABCD/denovo_v_cases_LDAK
./ldak5.1.linux  --grm ldakgrm2 --covar covar_all_ABCD_LDAK.txt --pheno pheno_ABCD.txt --prevalence 0.0288 --mpheno 30 --kinship-details NO --pcgc ./Results_ABCD/


./ldak5.1.linux  --grm ldakgrm3 --covar covar_all_withoutsex_ABCD_LDAK.txt --pheno pheno_ABCD.txt --prevalence 0.0288 --mpheno 7 --kinship-details NO --adjusted NO --pcgc ./Results_ABCD/high_males_LDAK
./ldak5.1.linux  --grm ldakgrm3 --covar covar_all_withoutsex_ABCD_LDAK.txt --pheno pheno_ABCD.txt --prevalence 0.0072 --mpheno 8 --kinship-details NO --adjusted NO --pcgc ./Results_ABCD/high_females_LDAK
./ldak5.1.linux  --grm ldakgrm3 --covar covar_all_withoutsex_ABCD_LDAK.txt --pheno pheno_ABCD.txt --prevalence 0.00432 --mpheno 9 --kinship-details NO --adjusted NO --pcgc ./Results_ABCD/high_nonIDfemales_LDAK
./ldak5.1.linux  --grm ldakgrm3 --covar covar_all_withoutsex_ABCD_LDAK.txt --pheno pheno_ABCD.txt --prevalence 0.019584 --mpheno 10 --kinship-details NO --adjusted NO --pcgc ./Results_ABCD/high_nonIDmales_LDAK
./ldak5.1.linux  --grm ldakgrm3 --covar covar_all_withoutsex_ABCD_LDAK.txt --pheno pheno_ABCD.txt --prevalence 0.002808 --mpheno 11 --kinship-details NO --adjusted NO --pcgc ./Results_ABCD/high_IDfemales_LDAK
./ldak5.1.linux  --grm ldakgrm3 --covar covar_all_withoutsex_ABCD_LDAK.txt --pheno pheno_ABCD.txt --prevalence 0.009216 --mpheno 12 --kinship-details NO --adjusted NO --pcgc ./Results_ABCD/high_IDmales_LDAK
./ldak5.1.linux  --grm ldakgrm3 --covar covar_all_withoutsex_ABCD_LDAK.txt --pheno pheno_ABCD.txt --prevalence 0.0288 --mpheno 30 --kinship-details NO --adjusted NO --pcgc ./Results_ABCD/males_femalesdownsampled_LDAK










###For bivariate rg
pheno = fread("~/Autism_heterogeneity/GRMs/GRMs_inc_UKB/pheno_ABCD.txt")
cases = subset(pheno, case_control == 1)
controls = subset(pheno, case_control == 0)
controls1 = sample_n(controls, 2241)


pheno1 = pheno[,c("FID", "IID", "denovo", "ID")]
pheno1 = na.omit(pheno1) #146 samples
write.table(pheno1, file = "~/Autism_heterogeneity/GRMs/GRMs_inc_UKB/bivar_pheno.txt", row.names = F, col.names = T, quote = F)

./gcta64 --reml-bivar 1 2 --grm ABCD_autism_ldpruned --pheno bivar_pheno.txt --thread-num 15  --out ./Results_ABCD/rg_denovoID --reml-bivar-prevalence 0.01 0.01










./gcta64 --reml --grm ABCD_autism_ldpruned --pheno pheno_ABCD.txt --covar covar_ABCD.txt --qcovar qcovar_ABCD.txt --thread-num 20 --mpheno 1 --out ./Results_ABCD/Case_control_ldpruned --prevalence 0.01
./gcta64 --reml --grm ABCD_autism_ldpruned --pheno pheno_ABCD.txt --covar covar_ABCD.txt --qcovar qcovar_ABCD.txt --thread-num 15 --mpheno 4 --out ./Results_ABCD/ID_ldpruned --prevalence 0.01

./gcta64 --reml --grm ABCD_autism_ldpruned --pheno pheno_ABCD.txt --covar covar_ABCD.txt --qcovar qcovar_ABCD.txt --thread-num 20 --mpheno 1 --out ./Results_ABCD/Case_control_ldpruned --prevalence 0.01
./gcta64 --reml --grm ABCD_autism_ldpruned --pheno pheno_ABCD.txt --covar covar_ABCD.txt --qcovar qcovar_ABCD.txt --thread-num 15 --mpheno 2 --out ./Results_ABCD/denovo_ldpruned --prevalence 0.01
./gcta64 --reml --grm ABCD_autism_ldpruned --pheno pheno_ABCD.txt --covar covar_ABCD.txt --qcovar qcovar_ABCD.txt --thread-num 15 --mpheno 3 --out ./Results_ABCD/nondenovo_ldpruned --prevalence 0.01
./gcta64 --reml --grm ABCD_autism_ldpruned --pheno pheno_ABCD.txt --covar covar_ABCD.txt --qcovar qcovar_ABCD.txt --thread-num 15 --mpheno 4 --out ./Results_ABCD/ID_ldpruned --prevalence 0.01
./gcta64 --reml --grm ABCD_autism_ldpruned --pheno pheno_ABCD.txt --covar covar_ABCD.txt --qcovar qcovar_ABCD.txt --thread-num 15 --mpheno 5 --out ./Results_ABCD/nonID_ldpruned --prevalence 0.01
./gcta64 --reml --grm ABCD_autism_ldpruned --pheno pheno_ABCD.txt --covar covar_ABCD.txt --qcovar qcovar_ABCD.txt --thread-num 15 --mpheno 6 --out ./Results_ABCD/highIQ_ldpruned --prevalence 0.01
./gcta64 --reml --grm ABCD_autism_ldpruned --pheno pheno_ABCD.txt --qcovar qcovar_ABCD.txt --thread-num 15 --mpheno 7 --out ./Results_ABCD/males_ldpruned --prevalence 0.01
./gcta64 --reml --grm ABCD_autism_ldpruned --pheno pheno_ABCD.txt --qcovar qcovar_ABCD.txt --thread-num 15 --mpheno 8 --out ./Results_ABCD/females_ldpruned --prevalence 0.01
./gcta64 --reml --grm ABCD_autism_ldpruned --pheno pheno_ABCD.txt --qcovar qcovar_ABCD.txt --thread-num 15 --mpheno 9 --out ./Results_ABCD/nonID_females_ldpruned --prevalence 0.01
./gcta64 --reml --grm ABCD_autism_ldpruned --pheno pheno_ABCD.txt --qcovar qcovar_ABCD.txt --thread-num 15 --mpheno 10 --out ./Results_ABCD/nonID_males_ldpruned --prevalence 0.01
./gcta64 --reml --grm ABCD_autism_ldpruned --pheno pheno_ABCD.txt --qcovar qcovar_ABCD.txt --thread-num 15 --mpheno 11 --out ./Results_ABCD/ID_female_ldpruned --prevalence 0.01
./gcta64 --reml --grm ABCD_autism_ldpruned --pheno pheno_ABCD.txt --qcovar qcovar_ABCD.txt --thread-num 15 --mpheno 12 --out ./Results_ABCD/ID_male_ldpruned --prevalence 0.01
./gcta64 --reml --grm ABCD_autism_ldpruned --pheno pheno_ABCD.txt --covar covar_ABCD.txt --qcovar qcovar_ABCD.txt --thread-num 15 --mpheno 13 --out ./Results_ABCD/high_scq_ldpruned --prevalence 0.01
./gcta64 --reml --grm ABCD_autism_ldpruned --pheno pheno_ABCD.txt --covar covar_ABCD.txt --qcovar qcovar_ABCD.txt --thread-num 15 --mpheno 14 --out ./Results_ABCD/high_rbs_ldpruned --prevalence 0.01
./gcta64 --reml --grm ABCD_autism_ldpruned --pheno pheno_ABCD.txt --covar covar_ABCD.txt --qcovar qcovar_ABCD.txt --thread-num 15 --mpheno 15 --out ./Results_ABCD/high_dcdq_ldpruned --prevalence 0.01
./gcta64 --reml --grm ABCD_autism_ldpruned --pheno pheno_ABCD.txt --covar covar_ABCD.txt --qcovar qcovar_ABCD.txt --thread-num 15 --mpheno 16 --out ./Results_ABCD/high_PA1_ldpruned --prevalence 0.01
./gcta64 --reml --grm ABCD_autism_ldpruned --pheno pheno_ABCD.txt --covar covar_ABCD.txt --qcovar qcovar_ABCD.txt --thread-num 15 --mpheno 17 --out ./Results_ABCD/high_PA2_ldpruned --prevalence 0.01
./gcta64 --reml --grm ABCD_autism_ldpruned --pheno pheno_ABCD.txt --covar covar_ABCD.txt --qcovar qcovar_ABCD.txt --thread-num 15 --mpheno 18 --out ./Results_ABCD/high_PA3_ldpruned --prevalence 0.01
./gcta64 --reml --grm ABCD_autism_ldpruned --pheno pheno_ABCD.txt --covar covar_ABCD.txt --qcovar qcovar_ABCD.txt --thread-num 15 --mpheno 19 --out ./Results_ABCD/high_PA4_ldpruned --prevalence 0.01
./gcta64 --reml --grm ABCD_autism_ldpruned --pheno pheno_ABCD.txt --covar covar_ABCD.txt --qcovar qcovar_ABCD.txt --thread-num 15 --mpheno 20 --out ./Results_ABCD/high_PA5_ldpruned --prevalence 0.01
./gcta64 --reml --grm ABCD_autism_ldpruned --pheno pheno_ABCD.txt --covar covar_ABCD.txt --qcovar qcovar_ABCD.txt --thread-num 15 --mpheno 21 --out ./Results_ABCD/high_PA6_ldpruned --prevalence 0.01
./gcta64 --reml --grm ABCD_autism_ldpruned --pheno pheno_ABCD.txt --covar covar_ABCD.txt --qcovar qcovar_ABCD.txt --thread-num 15 --mpheno 22 --out ./Results_ABCD/higher_PA1_ldpruned --prevalence 0.01
./gcta64 --reml --grm ABCD_autism_ldpruned --pheno pheno_ABCD.txt --covar covar_ABCD.txt --qcovar qcovar_ABCD.txt --thread-num 15 --mpheno 23 --out ./Results_ABCD/higher_PA2_ldpruned --prevalence 0.01









####old lines
./gcta64 --reml --mgrm multi_GRMs.txt --pheno pheno_ABCD.txt --covar covar_ABCD.txt --qcovar qcovar_ABCD.txt --thread-num 20 --mpheno 1 --out ./Results_ABCD/Case_control --prevalence 0.01
./gcta64 --reml --mgrm multi_GRMs.txt --pheno pheno_ABCD.txt --covar covar_ABCD.txt --qcovar qcovar_ABCD.txt --thread-num 15 --mpheno 2 --out ./Results_ABCD/denovo --prevalence 0.01
./gcta64 --reml --mgrm multi_GRMs.txt --pheno pheno_ABCD.txt --covar covar_ABCD.txt --qcovar qcovar_ABCD.txt --thread-num 15 --mpheno 3 --out ./Results_ABCD/nondenovo --prevalence 0.01
./gcta64 --reml --mgrm multi_GRMs.txt --pheno pheno_ABCD.txt --covar covar_ABCD.txt --qcovar qcovar_ABCD.txt --thread-num 15 --mpheno 4 --out ./Results_ABCD/ID --prevalence 0.01
./gcta64 --reml --mgrm multi_GRMs.txt --pheno pheno_ABCD.txt --covar covar_ABCD.txt --qcovar qcovar_ABCD.txt --thread-num 15 --mpheno 5 --out ./Results_ABCD/nonID --prevalence 0.01
./gcta64 --reml --mgrm multi_GRMs.txt --pheno pheno_ABCD.txt --covar covar_ABCD.txt --qcovar qcovar_ABCD.txt --thread-num 15 --mpheno 6 --out ./Results_ABCD/highIQ --prevalence 0.01
./gcta64 --reml --mgrm multi_GRMs.txt --pheno pheno_ABCD.txt --qcovar qcovar_ABCD.txt --thread-num 15 --mpheno 7 --out ./Results_ABCD/males --prevalence 0.0288
./gcta64 --reml --mgrm multi_GRMs.txt --pheno pheno_ABCD.txt --qcovar qcovar_ABCD.txt --thread-num 15 --mpheno 8 --out ./Results_ABCD/females --prevalence 0.0072
./gcta64 --reml --mgrm multi_GRMs.txt --pheno pheno_ABCD.txt --qcovar qcovar_ABCD.txt --thread-num 15 --mpheno 9 --out ./Results_ABCD/nonID_females --prevalence 0.004392
./gcta64 --reml --mgrm multi_GRMs.txt --pheno pheno_ABCD.txt --qcovar qcovar_ABCD.txt --thread-num 15 --mpheno 10 --out ./Results_ABCD/nonID_males --prevalence 0.019584
./gcta64 --reml --mgrm multi_GRMs.txt --pheno pheno_ABCD.txt --qcovar qcovar_ABCD.txt --thread-num 15 --mpheno 11 --out ./Results_ABCD/ID_female --prevalence 0.002808
./gcta64 --reml --mgrm multi_GRMs.txt --pheno pheno_ABCD.txt --qcovar qcovar_ABCD.txt --thread-num 15 --mpheno 12 --out ./Results_ABCD/ID_male --prevalence 0.009216
./gcta64 --reml --mgrm multi_GRMs.txt --pheno pheno_ABCD.txt --covar covar_ABCD.txt --qcovar qcovar_ABCD.txt --thread-num 15 --mpheno 13 --out ./Results_ABCD/high_scq --prevalence 0.01
./gcta64 --reml --mgrm multi_GRMs.txt --pheno pheno_ABCD.txt --covar covar_ABCD.txt --qcovar qcovar_ABCD.txt --thread-num 15 --mpheno 14 --out ./Results_ABCD/high_rbs --prevalence 0.01
./gcta64 --reml --mgrm multi_GRMs.txt --pheno pheno_ABCD.txt --covar covar_ABCD.txt --qcovar qcovar_ABCD.txt --thread-num 15 --mpheno 15 --out ./Results_ABCD/high_dcdq --prevalence 0.01
./gcta64 --reml --mgrm multi_GRMs.txt --pheno pheno_ABCD.txt --covar covar_ABCD.txt --qcovar qcovar_ABCD.txt --thread-num 15 --mpheno 16 --out ./Results_ABCD/high_PA1 --prevalence 0.01
./gcta64 --reml --mgrm multi_GRMs.txt --pheno pheno_ABCD.txt --covar covar_ABCD.txt --qcovar qcovar_ABCD.txt --thread-num 15 --mpheno 17 --out ./Results_ABCD/high_PA2 --prevalence 0.01
./gcta64 --reml --mgrm multi_GRMs.txt --pheno pheno_ABCD.txt --covar covar_ABCD.txt --qcovar qcovar_ABCD.txt --thread-num 15 --mpheno 18 --out ./Results_ABCD/high_PA3 --prevalence 0.01
./gcta64 --reml --mgrm multi_GRMs.txt --pheno pheno_ABCD.txt --covar covar_ABCD.txt --qcovar qcovar_ABCD.txt --thread-num 15 --mpheno 19 --out ./Results_ABCD/high_PA4 --prevalence 0.01
./gcta64 --reml --mgrm multi_GRMs.txt --pheno pheno_ABCD.txt --covar covar_ABCD.txt --qcovar qcovar_ABCD.txt --thread-num 15 --mpheno 20 --out ./Results_ABCD/high_PA5 --prevalence 0.01
./gcta64 --reml --mgrm multi_GRMs.txt --pheno pheno_ABCD.txt --covar covar_ABCD.txt --qcovar qcovar_ABCD.txt --thread-num 15 --mpheno 21 --out ./Results_ABCD/high_PA6 --prevalence 0.01
./gcta64 --reml --mgrm multi_GRMs.txt --pheno pheno_ABCD.txt --covar covar_ABCD.txt --qcovar qcovar_ABCD.txt --thread-num 15 --mpheno 22 --out ./Results_ABCD/higher_PA1 --prevalence 0.01
./gcta64 --reml --mgrm multi_GRMs.txt --pheno pheno_ABCD.txt --covar covar_ABCD.txt --qcovar qcovar_ABCD.txt --thread-num 15 --mpheno 23 --out ./Results_ABCD/higher_PA2 --prevalence 0.01










for i in {1..30}; do ./gcta64 --reml --mgrm multi_GRMs.txt --pheno med_pheno_forGRMs_SPARK.txt --covar covar_ABCD.txt --qcovar qcovar_ABCD.txt --thread-num 20 --mpheno ${i} --out ./Results_medpheno/Case_control_${i} --prevalence 0.01; done
for i in {31..40}; do ./gcta64 --reml --mgrm multi_GRMs.txt --pheno med_pheno_forGRMs_SPARK.txt --covar covar_ABCD.txt --qcovar qcovar_ABCD.txt --thread-num 20 --mpheno ${i} --out ./Results_medpheno/quantitative_${i}; done


# for ldms
./gcta64  --bfile ~/SPARK/SPARK_v1v2/SPARK_v1v2_eur_allchrs_hg19 --keep autism_spark_all_unrelated.grm.id --ld-score-region 200  --out ABCD_autism_ldscore_byregion --thread-num 15

lds_seg = read.table("ABCD_autism_ldscore_byregion.score.ld",header=T,colClasses=c("character",rep("numeric",8)))
quartiles=summary(lds_seg$ldscore_SNP)

lb1 = which(lds_seg$ldscore_SNP <= quartiles[2])
lb2 = which(lds_seg$ldscore_SNP > quartiles[2] & lds_seg$ldscore_SNP <= quartiles[3])
lb3 = which(lds_seg$ldscore_SNP > quartiles[3] & lds_seg$ldscore_SNP <= quartiles[5])
lb4 = which(lds_seg$ldscore_SNP > quartiles[5])

lb1_snp = lds_seg$SNP[lb1]
lb2_snp = lds_seg$SNP[lb2]
lb3_snp = lds_seg$SNP[lb3]
lb4_snp = lds_seg$SNP[lb4]

write.table(lb1_snp, "snp_group1_SPARK.txt", row.names=F, quote=F, col.names=F)
write.table(lb2_snp, "snp_group2_SPARK.txt", row.names=F, quote=F, col.names=F)
write.table(lb3_snp, "snp_group3_SPARK.txt", row.names=F, quote=F, col.names=F)
write.table(lb4_snp, "snp_group4_SPARK.txt", row.names=F, quote=F, col.names=F)

./gcta64 --bfile ~/SPARK/SPARK_v1v2/SPARK_v1v2_eur_allchrs_hg19 --extract snp_group1_SPARK.txt --keep autism_spark_all_unrelated.grm.id --make-grm --out SPARK_autism_greml_snpgroup1_common --thread-num 20 --maf 0.05
./gcta64 --bfile ~/SPARK/SPARK_v1v2/SPARK_v1v2_eur_allchrs_hg19 --extract snp_group2_SPARK.txt --keep autism_spark_all_unrelated.grm.id --make-grm --out SPARK_autism_greml_snpgroup2_common --thread-num 20 --maf 0.05
./gcta64 --bfile ~/SPARK/SPARK_v1v2/SPARK_v1v2_eur_allchrs_hg19 --extract snp_group3_SPARK.txt --keep autism_spark_all_unrelated.grm.id --make-grm --out SPARK_autism_greml_snpgroup3_common --thread-num 20 --maf 0.05
./gcta64 --bfile ~/SPARK/SPARK_v1v2/SPARK_v1v2_eur_allchrs_hg19 --extract snp_group4_SPARK.txt --keep autism_spark_all_unrelated.grm.id --make-grm --out SPARK_autism_greml_snpgroup4_common --thread-num 20 --maf 0.05

./gcta64 --bfile ~/SPARK/SPARK_v1v2/SPARK_v1v2_eur_allchrs_hg19 --extract snp_group1_SPARK.txt --keep autism_spark_all_unrelated.grm.id --make-grm --out SPARK_autism_greml_snpgroup1_lowfreq --thread-num 20 --maf 0.005 --max-maf 0.0499
./gcta64 --bfile ~/SPARK/SPARK_v1v2/SPARK_v1v2_eur_allchrs_hg19 --extract snp_group2_SPARK.txt --keep autism_spark_all_unrelated.grm.id --make-grm --out SPARK_autism_greml_snpgroup2_lowfreq --thread-num 20 --maf 0.005 --max-maf 0.0499
./gcta64 --bfile ~/SPARK/SPARK_v1v2/SPARK_v1v2_eur_allchrs_hg19 --extract snp_group3_SPARK.txt --keep autism_spark_all_unrelated.grm.id --make-grm --out SPARK_autism_greml_snpgroup3_lowfreq --thread-num 20 --maf 0.005 --max-maf 0.0499
./gcta64 --bfile ~/SPARK/SPARK_v1v2/SPARK_v1v2_eur_allchrs_hg19 --extract snp_group4_SPARK.txt --keep autism_spark_all_unrelated.grm.id --make-grm --out SPARK_autism_greml_snpgroup4_lowfreq --thread-num 20 --maf 0.005 --max-maf 0.0499

./gcta64 --bfile ~/SPARK/SPARK_v1v2/SPARK_v1v2_eur_allchrs_hg19 --extract snp_group1_SPARK.txt --keep autism_spark_all_unrelated.grm.id --make-grm --out SPARK_autism_greml_snpgroup1_rare --thread-num 20 --maf 0.0005 --max-maf 0.00499
./gcta64 --bfile ~/SPARK/SPARK_v1v2/SPARK_v1v2_eur_allchrs_hg19 --extract snp_group2_SPARK.txt --keep autism_spark_all_unrelated.grm.id --make-grm --out SPARK_autism_greml_snpgroup2_rare --thread-num 20 --maf 0.0005 --max-maf 0.00499
./gcta64 --bfile ~/SPARK/SPARK_v1v2/SPARK_v1v2_eur_allchrs_hg19 --extract snp_group3_SPARK.txt --keep autism_spark_all_unrelated.grm.id --make-grm --out SPARK_autism_greml_snpgroup3_rare --thread-num 20 --maf 0.0005 --max-maf 0.00499
./gcta64 --bfile ~/SPARK/SPARK_v1v2/SPARK_v1v2_eur_allchrs_hg19 --extract snp_group4_SPARK.txt --keep autism_spark_all_unrelated.grm.id --make-grm --out SPARK_autism_greml_snpgroup4_rare --thread-num 20 --maf 0.0005 --max-maf 0.00499


for i in {1..30}; do ./gcta64 --reml --mgrm mgrm_list_all.txt --pheno ./GRMs_inc_UKB/med_pheno_forGRMs_SPARK.txt --covar ./GRMs_inc_UKB/covar_ABCD.txt --qcovar ./GRMs_inc_UKB/qcovar_ABCD.txt --thread-num 20 --mpheno ${i} --out ./GRMs_inc_UKB/Results_medpheno/Case_control_all${i} --prevalence 0.01; done
for i in {31..40}; do ./gcta64 --reml --mgrm mgrm_list_all.txt --pheno ./GRMs_inc_UKB/med_pheno_forGRMs_SPARK.txt --covar ./GRMs_inc_UKB/covar_ABCD.txt --qcovar ./GRMs_inc_UKB/qcovar_ABCD.txt --thread-num 20 --mpheno ${i} --out ./GRMs_inc_UKB/Results_medpheno/quantitative_${i}; done



./gcta64 --bfile allautism_ABCDmerge_forGRMandPC_unrelated   --make-grm --out ABCD_autism_ldpruned --thread-num 20

./gcta64 --bfile allautism_ABCDmerge_SNPQC --extract ~/Autism_heterogeneity/Commonwellimputed_SNPs_ABCD_SSC_SPARK_highr2.txt  --make-grm --out ABCD_autism_ldpruned_highr2 --thread-num 20

#cases-only
./gcta64 --bfile allautism_ABCDmerge_SNPQC --keep cases_to_keep.txt --extract ~/Autism_heterogeneity/Commonwellimputed_SNPs_ABCD_SSC_SPARK_highr2.txt  --make-grm --out autismonly_ldpruned_highr2 --thread-num 20
./plink --bfile allautism_ABCDmerge_forGRMandPC_unrelated --keep cases_to_keep.txt --pca --out autismonly_pcs --threads 15

