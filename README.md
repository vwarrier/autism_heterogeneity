# autism_heterogeneity
This git has the code for the heterogeneity in autism project


Polygenic scoring: 


LEAP-sample
for i in savageCP2018forprsice edu2018leeforPRSICE; do Rscript PRSice2.R --dir . --prsice ./PRSice_linux --base ./pgssumstats/${i}.txt --target ~/LEAP_PARIS/LEAP_v2/LEAP-ImputedData-Freeze-V2-June-2019 --base-maf 0.01 --thread 20 --stat BETA --binary-target T --out ./PRSice2results/LEAP_v2_${i} --bar-levels 0.0001,0.001,0.005,0.01,0.05,0.1,0.25,0.5,0.75,1 --no-regress --fastscore; done
