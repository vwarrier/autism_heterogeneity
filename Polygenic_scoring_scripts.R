#We are using PRScs to generate PGS. More information can be found here: https://github.com/getian107/PRScs
#PRScs requires the following columns
# SNP  A1   A2   BETA/OR      P

#Step 1: Generate the posterior SNP effects in bash

#SPARK
python PRScs.py --ref_dir=./ldblk_1kg_eur --bim_prefix=/mnt/b2/home4/arc/vw260/SPARK/SPARK_postimputation/CEU/Plink_files/SPARKCEU_hg19_allchrs --sst_file=/mnt/b2/home4/arc/vw260/Autism_heterogeneity/PGC_SCZ3.txt --n_gwas=161405 --out_dir=/mnt/b2/home4/arc/vw260/Autism_heterogeneity/SPARK_scz  --phi=1e-2 
python PRScs.py --ref_dir=./ldblk_1kg_eur --bim_prefix=/mnt/b2/home4/arc/vw260/SPARK/SPARK_postimputation/CEU/Plink_files/SPARKCEU_hg19_allchrs --sst_file=/mnt/b2/home4/arc/vw260/Autism_heterogeneity/Lee_edu.txt --n_gwas=766341 --out_dir=/mnt/b2/home4/arc/vw260/Autism_heterogeneity/SPARK_edu  --phi=1e-2 
python PRScs.py --ref_dir=./ldblk_1kg_eur --bim_prefix=/mnt/b2/home4/arc/vw260/SPARK/SPARK_postimputation/CEU/Plink_files/SPARKCEU_hg19_allchrs --sst_file=/mnt/b2/home4/arc/vw260/Autism_heterogeneity/Savage_cognition.txt --n_gwas=269867 --out_dir=/mnt/b2/home4/arc/vw260/Autism_heterogeneity/SPARK_IQ  --phi=1e-2 
python PRScs.py --ref_dir=./ldblk_1kg_eur --bim_prefix=/mnt/b2/home4/arc/vw260/SPARK/SPARK_postimputation/CEU/Plink_files/SPARKCEU_hg19_allchrs --sst_file=/mnt/b2/home4/arc/vw260/Autism_heterogeneity/autism_grove_ipsychonly2020.txt --n_gwas=58948 --out_dir=/mnt/b2/home4/arc/vw260/Autism_heterogeneity/SPARK_autism  --phi=1e-2 
python PRScs.py --ref_dir=./ldblk_1kg_eur --bim_prefix=/mnt/b2/home4/arc/vw260/SPARK/SPARK_postimputation/CEU/Plink_files/SPARKCEU_hg19_allchrs --sst_file=/mnt/b2/home4/arc/vw260/Autism_heterogeneity/haircolour_negcontrol.txt --n_gwas=385603 --out_dir=/mnt/b2/home4/arc/vw260/Autism_heterogeneity/SPARK_haircolour  --phi=1e-2 
python PRScs.py --ref_dir=./ldblk_1kg_eur --bim_prefix=/mnt/b2/home4/arc/vw260/SPARK/SPARK_postimputation/CEU/Plink_files/SPARKCEU_hg19_allchrs --sst_file=/mnt/b2/home4/arc/vw260/Autism_heterogeneity/ADHD.txt --n_gwas=55374 --out_dir=/mnt/b2/home4/arc/vw260/Autism_heterogeneity/temp_pgs_files/SPARK_ADHD  --phi=1e-2 

#SFARI
python PRScs.py --ref_dir=./ldblk_1kg_eur --bim_prefix=/mnt/b2/home4/arc/vw260/SFARI/liftOverPlink/files_imputed/SFARImergedall --sst_file=/mnt/b2/home4/arc/vw260/Autism_heterogeneity/PGC_SCZ3.txt --n_gwas=161405 --out_dir=/mnt/b2/home4/arc/vw260/Autism_heterogeneity/SSC_scz  --phi=1e-2 
python PRScs.py --ref_dir=./ldblk_1kg_eur --bim_prefix=/mnt/b2/home4/arc/vw260/SFARI/liftOverPlink/files_imputed/SFARImergedall --sst_file=/mnt/b2/home4/arc/vw260/Autism_heterogeneity/Lee_edu.txt --n_gwas=766341 --out_dir=/mnt/b2/home4/arc/vw260/Autism_heterogeneity/SSC_edu  --phi=1e-2 
python PRScs.py --ref_dir=./ldblk_1kg_eur --bim_prefix=/mnt/b2/home4/arc/vw260/SFARI/liftOverPlink/files_imputed/SFARImergedall --sst_file=/mnt/b2/home4/arc/vw260/Autism_heterogeneity/Savage_cognition.txt --n_gwas=269867 --out_dir=/mnt/b2/home4/arc/vw260/Autism_heterogeneity/SSC_IQ  --phi=1e-2 
python PRScs.py --ref_dir=./ldblk_1kg_eur --bim_prefix=/mnt/b2/home4/arc/vw260/SFARI/liftOverPlink/files_imputed/SFARImergedall --sst_file=/mnt/b2/home4/arc/vw260/Autism_heterogeneity/autism_grove_ipsychonly2020.txt --n_gwas=58948 --out_dir=/mnt/b2/home4/arc/vw260/Autism_heterogeneity/SSC_autism  --phi=1e-2 
python PRScs.py --ref_dir=./ldblk_1kg_eur --bim_prefix=/mnt/b2/home4/arc/vw260/SFARI/liftOverPlink/files_imputed/SFARImergedall --sst_file=/mnt/b2/home4/arc/vw260/Autism_heterogeneity/haircolour_negcontrol.txt --n_gwas=385603 --out_dir=/mnt/b2/home4/arc/vw260/Autism_heterogeneity/SSC_haircolour  --phi=1e-2 
python PRScs.py --ref_dir=./ldblk_1kg_eur --bim_prefix=/mnt/b2/home4/arc/vw260/SFARI/liftOverPlink/files_imputed/SFARImergedall --sst_file=/mnt/b2/home4/arc/vw260/Autism_heterogeneity/ADHD.txt --n_gwas=55374 --out_dir=/mnt/b2/home4/arc/vw260/Autism_heterogeneity/temp_pgs_files/SSC_ADHD  --phi=1e-2 

#LEAP
python PRScs.py --ref_dir=./ldblk_1kg_eur --bim_prefix=/mnt/b2/home4/arc/vw260/LEAP_PARIS/LEAP_v2/LEAP-ImputedData-Freeze-V2-June-2019 --sst_file=/mnt/b2/home4/arc/vw260/Autism_heterogeneity/PGC_SCZ3.txt --n_gwas=161405 --out_dir=/mnt/b2/home4/arc/vw260/Autism_heterogeneity/LEAP_scz  --phi=1e-2 
python PRScs.py --ref_dir=./ldblk_1kg_eur --bim_prefix=/mnt/b2/home4/arc/vw260/LEAP_PARIS/LEAP_v2/LEAP-ImputedData-Freeze-V2-June-2019 --sst_file=/mnt/b2/home4/arc/vw260/Autism_heterogeneity/Lee_edu.txt --n_gwas=766341 --out_dir=/mnt/b2/home4/arc/vw260/Autism_heterogeneity/LEAP_edu  --phi=1e-2 
python PRScs.py --ref_dir=./ldblk_1kg_eur --bim_prefix=/mnt/b2/home4/arc/vw260/LEAP_PARIS/LEAP_v2/LEAP-ImputedData-Freeze-V2-June-2019 --sst_file=/mnt/b2/home4/arc/vw260/Autism_heterogeneity/Savage_cognition.txt --n_gwas=269867 --out_dir=/mnt/b2/home4/arc/vw260/Autism_heterogeneity/LEAP_IQ  --phi=1e-2 
python PRScs.py --ref_dir=./ldblk_1kg_eur --bim_prefix=/mnt/b2/home4/arc/vw260/LEAP_PARIS/LEAP_v2/LEAP-ImputedData-Freeze-V2-June-2019 --sst_file=/mnt/b2/home4/arc/vw260/Autism_heterogeneity/autism_grove_ipsychonly2020.txt --n_gwas=58948 --out_dir=/mnt/b2/home4/arc/vw260/Autism_heterogeneity/LEAP_autism  --phi=1e-2 
python PRScs.py --ref_dir=./ldblk_1kg_eur --bim_prefix=/mnt/b2/home4/arc/vw260/LEAP_PARIS/LEAP_v2/LEAP-ImputedData-Freeze-V2-June-2019 --sst_file=/mnt/b2/home4/arc/vw260/Autism_heterogeneity/haircolour_negcontrol.txt --n_gwas=385603 --out_dir=/mnt/b2/home4/arc/vw260/Autism_heterogeneity/LEAP_haircolour  --phi=1e-2 
python PRScs.py --ref_dir=./ldblk_1kg_eur --bim_prefix=/mnt/b2/home4/arc/vw260/LEAP_PARIS/LEAP_v2/LEAP-ImputedData-Freeze-V2-June-2019 --sst_file=/mnt/b2/home4/arc/vw260/Autism_heterogeneity/ADHD.txt --n_gwas=55374 --out_dir=/mnt/b2/home4/arc/vw260/Autism_heterogeneity/temp_pgs_files/LEAP_ADHD  --phi=1e-2 

#AGRE

python PRScs.py --ref_dir=./ldblk_1kg_eur --bim_prefix=/mnt/b2/home4/arc/vw260/AGRE_data/CHOP/imputed_plinkfile/CHOP_mergedQC --sst_file=/mnt/b2/home4/arc/vw260/Autism_heterogeneity/PGC_SCZ3.txt --n_gwas=161405 --out_dir=/mnt/b2/home4/arc/vw260/Autism_heterogeneity/CHOP_scz  --phi=1e-2 
python PRScs.py --ref_dir=./ldblk_1kg_eur --bim_prefix=/mnt/b2/home4/arc/vw260/AGRE_data/CHOP/imputed_plinkfile/CHOP_mergedQC --sst_file=/mnt/b2/home4/arc/vw260/Autism_heterogeneity/Lee_edu.txt --n_gwas=766341 --out_dir=/mnt/b2/home4/arc/vw260/Autism_heterogeneity/CHOP_edu  --phi=1e-2 
python PRScs.py --ref_dir=./ldblk_1kg_eur --bim_prefix=/mnt/b2/home4/arc/vw260/AGRE_data/CHOP/imputed_plinkfile/CHOP_mergedQC --sst_file=/mnt/b2/home4/arc/vw260/Autism_heterogeneity/Savage_cognition.txt --n_gwas=269867 --out_dir=/mnt/b2/home4/arc/vw260/Autism_heterogeneity/CHOP_IQ  --phi=1e-2 
python PRScs.py --ref_dir=./ldblk_1kg_eur --bim_prefix=/mnt/b2/home4/arc/vw260/AGRE_data/CHOP/imputed_plinkfile/CHOP_mergedQC --sst_file=/mnt/b2/home4/arc/vw260/Autism_heterogeneity/autism_grove_ipsychonly2020.txt --n_gwas=58948 --out_dir=/mnt/b2/home4/arc/vw260/Autism_heterogeneity/CHOP_autism  --phi=1e-2 
python PRScs.py --ref_dir=./ldblk_1kg_eur --bim_prefix=/mnt/b2/home4/arc/vw260/AGRE_data/CHOP/imputed_plinkfile/CHOP_mergedQC --sst_file=/mnt/b2/home4/arc/vw260/Autism_heterogeneity/haircolour_negcontrol.txt --n_gwas=385603 --out_dir=/mnt/b2/home4/arc/vw260/Autism_heterogeneity/CHOP_haircolour --phi=1e-2 
python PRScs.py --ref_dir=./ldblk_1kg_eur --bim_prefix=/mnt/b2/home4/arc/vw260/AGRE_data/CHOP/imputed_plinkfile/CHOP_mergedQC --sst_file=/mnt/b2/home4/arc/vw260/Autism_heterogeneity/ADHD.txt --n_gwas=55374 --out_dir=/mnt/b2/home4/arc/vw260/Autism_heterogeneity/temp_pgs_files/CHOP_ADHD  --phi=1e-2 


##SPARK_v2
python PRScs.py --ref_dir=./ldblk_1kg_eur --bim_prefix=/mnt/b2/home4/arc/vw260/SPARK/SPARK_v2/Imputed_files/Plink_files/SPARK_hg19_allchrs_v2 --sst_file=/mnt/b2/home4/arc/vw260/Autism_heterogeneity/PGC_SCZ3.txt --n_gwas=161405 --out_dir=/mnt/b2/home4/arc/vw260/Autism_heterogeneity/temp_pgs_files/SPARK_v2_scz  --phi=1e-2 
python PRScs.py --ref_dir=./ldblk_1kg_eur --bim_prefix=/mnt/b2/home4/arc/vw260/SPARK/SPARK_v2/Imputed_files/Plink_files/SPARK_hg19_allchrs_v2 --sst_file=/mnt/b2/home4/arc/vw260/Autism_heterogeneity/Lee_edu.txt --n_gwas=766341 --out_dir=/mnt/b2/home4/arc/vw260/Autism_heterogeneity/temp_pgs_files/SPARK_v2_edu  --phi=1e-2 
python PRScs.py --ref_dir=./ldblk_1kg_eur --bim_prefix=/mnt/b2/home4/arc/vw260/SPARK/SPARK_v2/Imputed_files/Plink_files/SPARK_hg19_allchrs_v2 --sst_file=/mnt/b2/home4/arc/vw260/Autism_heterogeneity/Savage_cognition.txt --n_gwas=269867 --out_dir=/mnt/b2/home4/arc/vw260/Autism_heterogeneity/temp_pgs_files/SPARK_v2_IQ  --phi=1e-2 
python PRScs.py --ref_dir=./ldblk_1kg_eur --bim_prefix=/mnt/b2/home4/arc/vw260/SPARK/SPARK_v2/Imputed_files/Plink_files/SPARK_hg19_allchrs_v2 --sst_file=/mnt/b2/home4/arc/vw260/Autism_heterogeneity/autism_grove_ipsychonly2020.txt --n_gwas=58948 --out_dir=/mnt/b2/home4/arc/vw260/Autism_heterogeneity/temp_pgs_files/SPARK_v2_autism  --phi=1e-2 
python PRScs.py --ref_dir=./ldblk_1kg_eur --bim_prefix=/mnt/b2/home4/arc/vw260/SPARK/SPARK_v2/Imputed_files/Plink_files/SPARK_hg19_allchrs_v2 --sst_file=/mnt/b2/home4/arc/vw260/Autism_heterogeneity/haircolour_negcontrol.txt --n_gwas=385603 --out_dir=/mnt/b2/home4/arc/vw260/Autism_heterogeneity/temp_pgs_files/SPARK_v2_haircolour  --phi=1e-2 
python PRScs.py --ref_dir=./ldblk_1kg_eur --bim_prefix=/mnt/b2/home4/arc/vw260/SPARK/SPARK_v2/Imputed_files/Plink_files/SPARK_hg19_allchrs_v2 --sst_file=/mnt/b2/home4/arc/vw260/Autism_heterogeneity/ADHD.txt --n_gwas=55374 --out_dir=/mnt/b2/home4/arc/vw260/Autism_heterogeneity/temp_pgs_files/SPARK_v2_ADHD  --phi=1e-2 


###Step 2, generate a single SNP file, in R

##SPARK ###SSC ###LEAP ###CHOP

data1 = fread("CHOP_ADHD_pst_eff_a1_b0.5_phi1e-02_chr1.txt")

for (i in 2:22){
  data2 = fread(paste0("CHOP_ADHD_pst_eff_a1_b0.5_phi1e-02_chr", i, ".txt"))
  data1 = rbind(data1, data2)
}

write.table(data1, file = "~/Autism_heterogeneity/final_data_forPGS/ADHD_CHOP_score.txt", row.names = F, col.names = F, quote = F)

rm(list = ls())


### Step 3: Run score in bash

./plink --bfile ~/SPARK/SPARK_postimputation/CEU/Plink_files/SPARKCEU_hg19_allchrs --score ./final_data_forPGS/scz_spark_score.txt 2 4 6 center --out ./PGS/SPARK_scz_finalscore
./plink --bfile ~/SFARI/liftOverPlink/files_imputed/SFARImergedall --score ./final_data_forPGS/scz_ssc_score.txt 2 4 6 center --out ./PGS/SSC_scz_finalscore
./plink --bfile ~/LEAP_PARIS/LEAP_v2/LEAP-ImputedData-Freeze-V2-June-2019 --score ./final_data_forPGS/scz_LEAP_score.txt 2 4 6 center --out ./PGS/LEAP_scz_finalscore
./plink --bfile ~/AGRE_data/CHOP/imputed_plinkfile/CHOP_mergedQC --score ./final_data_forPGS/scz_CHOP_score.txt 2 4 6 center --out ./PGS/CHOP_scz_finalscore
./plink --bfile ~/SPARK/SPARK_v2/Imputed_files/Plink_files/SPARK_hg19_allchrs_v2 --score ./final_data_forPGS/scz_SPARK_v2_score.txt 2 4 6 center --out ./PGS/SPARK_v2_scz_finalscore


./plink --bfile ~/SPARK/SPARK_postimputation/CEU/Plink_files/SPARKCEU_hg19_allchrs --score ./final_data_forPGS/edu_SPARK_score.txt 2 4 6 center --out ./PGS/SPARK_edu_finalscore
./plink --bfile ~/SFARI/liftOverPlink/files_imputed/SFARImergedall --score ./final_data_forPGS/edu_SSC_score.txt 2 4 6 center --out ./PGS/SSC_edu_finalscore
./plink --bfile ~/LEAP_PARIS/LEAP_v2/LEAP-ImputedData-Freeze-V2-June-2019 --score ./final_data_forPGS/edu_LEAP_score.txt 2 4 6 center --out ./PGS/LEAP_edu_finalscore
./plink --bfile ~/AGRE_data/CHOP/imputed_plinkfile/CHOP_mergedQC --score ./final_data_forPGS/edu_CHOP_score.txt 2 4 6 center --out ./PGS/CHOP_edu_finalscore
./plink --bfile ~/SPARK/SPARK_v2/Imputed_files/Plink_files/SPARK_hg19_allchrs_v2 --score ./final_data_forPGS/edu_SPARK_v2_score.txt 2 4 6 center --out ./PGS/SPARK_v2_edu_finalscore


./plink --bfile ~/SPARK/SPARK_postimputation/CEU/Plink_files/SPARKCEU_hg19_allchrs --score ./final_data_forPGS/IQ_SPARK_score.txt 2 4 6 center --out ./PGS/SPARK_IQ_finalscore
./plink --bfile ~/SFARI/liftOverPlink/files_imputed/SFARImergedall --score ./final_data_forPGS/IQ_SSC_score.txt 2 4 6 center --out ./PGS/SSC_IQ_finalscore
./plink --bfile ~/LEAP_PARIS/LEAP_v2/LEAP-ImputedData-Freeze-V2-June-2019 --score ./final_data_forPGS/IQ_LEAP_score.txt 2 4 6 center --out ./PGS/LEAP_IQ_finalscore
./plink --bfile ~/AGRE_data/CHOP/imputed_plinkfile/CHOP_mergedQC --score ./final_data_forPGS/IQ_CHOP_score.txt 2 4 6 center --out ./PGS/CHOP_IQ_finalscore
./plink --bfile ~/SPARK/SPARK_v2/Imputed_files/Plink_files/SPARK_hg19_allchrs_v2 --score ./final_data_forPGS/IQ_SPARK_v2_score.txt 2 4 6 center --out ./PGS/SPARK_v2_IQ_finalscore


./plink --bfile ~/SPARK/SPARK_postimputation/CEU/Plink_files/SPARKCEU_hg19_allchrs --score ./final_data_forPGS/autism_SPARK_score.txt 2 4 6 center --out ./PGS/SPARK_autism_finalscore
./plink --bfile ~/SFARI/liftOverPlink/files_imputed/SFARImergedall --score ./final_data_forPGS/autism_SSC_score.txt 2 4 6 center --out ./PGS/SSC_autism_finalscore
./plink --bfile ~/LEAP_PARIS/LEAP_v2/LEAP-ImputedData-Freeze-V2-June-2019 --score ./final_data_forPGS/autism_LEAP_score.txt 2 4 6 center --out ./PGS/LEAP_autism_finalscore
./plink --bfile ~/AGRE_data/CHOP/imputed_plinkfile/CHOP_mergedQC --score ./final_data_forPGS/autism_CHOP_score.txt 2 4 6 center --out ./PGS/CHOP_autism_finalscore
./plink --bfile ~/SPARK/SPARK_v2/Imputed_files/Plink_files/SPARK_hg19_allchrs_v2 --score ./final_data_forPGS/autism_SPARK_v2_score.txt 2 4 6 center --out ./PGS/SPARK_v2_autism_finalscore


./plink --bfile ~/SPARK/SPARK_postimputation/CEU/Plink_files/SPARKCEU_hg19_allchrs --score ./final_data_forPGS/haircolour_SPARK_v1_score.txt 2 4 6 center --out ./PGS/SPARK_haircolour_finalscore
./plink --bfile ~/SFARI/liftOverPlink/files_imputed/SFARImergedall --score ./final_data_forPGS/haircolour_SSC_score.txt 2 4 6 center --out ./PGS/SSC_haircolour_finalscore
./plink --bfile ~/LEAP_PARIS/LEAP_v2/LEAP-ImputedData-Freeze-V2-June-2019 --score ./final_data_forPGS/haircolour_LEAP_score.txt 2 4 6 center --out ./PGS/LEAP_haircolour_finalscore
./plink --bfile ~/AGRE_data/CHOP/imputed_plinkfile/CHOP_mergedQC --score ./final_data_forPGS/haircolour_CHOP_score.txt 2 4 6 center --out ./PGS/CHOP_haircolour_finalscore
./plink --bfile ~/SPARK/SPARK_v2/Imputed_files/Plink_files/SPARK_hg19_allchrs_v2 --score ./final_data_forPGS/haircolour_SPARK_v2_score.txt 2 4 6 center --out ./PGS/SPARK_v2_haircolour_finalscore

./plink --bfile ~/SPARK/SPARK_postimputation/CEU/Plink_files/SPARKCEU_hg19_allchrs --score ./final_data_forPGS/ADHD_SPARK_v1_score.txt 2 4 6 center --out ./PGS/SPARK_ADHD_finalscore
./plink --bfile ~/SFARI/liftOverPlink/files_imputed/SFARImergedall --score ./final_data_forPGS/ADHD_SSC_score.txt 2 4 6 center --out ./PGS/SSC_ADHD_finalscore
./plink --bfile ~/LEAP_PARIS/LEAP_v2/LEAP-ImputedData-Freeze-V2-June-2019 --score ./final_data_forPGS/ADHD_LEAP_score.txt 2 4 6 center --out ./PGS/LEAP_ADHD_finalscore
./plink --bfile ~/AGRE_data/CHOP/imputed_plinkfile/CHOP_mergedQC --score ./final_data_forPGS/ADHD_CHOP_score.txt 2 4 6 center --out ./PGS/CHOP_ADHD_finalscore
./plink --bfile ~/SPARK/SPARK_v2/Imputed_files/Plink_files/SPARK_hg19_allchrs_v2 --score ./final_data_forPGS/ADHD_SPARK_v2_score.txt 2 4 6 center --out ./PGS/SPARK_v2_ADHD_finalscore
