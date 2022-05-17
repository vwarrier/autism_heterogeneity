###Factor analysis and invariance analyses
###The scripts for comprehensive set of factor analyses are provided elsewhere###


##Step 1: we will generate factor for ssc and spark combined###
library(lavaan)
library(semTools)
library(semPlot)


factor = fread("~/Autism_heterogeneity/data_for_factorall.txt")
factor_ssc = fread("~/Autism_heterogeneity/ssc_fa.txt")
setnames(factor_ssc, "FID", "individual")


factor_all = rbind(factor, factor_ssc)

model_gen_6 <-'
PA1 =~ r38+r29+r39+r33+r25+r26+r37+r34+r35+r27+r15+r24+r31+r32+r21+r30+r41+r16+s08+r23+r19+r20+r28+r36+r18+r40+r17
PA2 =~ s30+s28+s29+s21+s20+s22+s36+s34+s27+s37+s40+s31+s33+s39+s32+s26+s38+s35+s23+s09+0*s19
PA3 =~ r03+s15+r04+r05+s16+s12+r01+r42+r06+r43+r02+s14+r22
PA4 =~ r11+r07+s17+r10+r08+r09+r12+r14+r13
PA5 =~ s03+s07+s05+s06+s11+s10+0*s18
PA6 =~ s04+s02+s24+s13+s25
'
fit_gen_6 <- cfa(model_gen_6, data = factor_all, std.lv = TRUE,ordered =  c('s02', 's03', 's04', 's05', 's06', 's07', 's08', 's09', 's10', 's11', 's12', 's13', 's14', 's15', 's16', 's17', 's18', 's19', 's20', 's21', 's22', 's23', 's24', 's25', 's26', 's27', 's28', 's29', 's30', 's31', 's32', 's33', 's34', 's35', 's36', 's37', 's38', 's39', 's40', 'r01', 'r02', 'r03', 'r04', 'r05', 'r06', 'r07', 'r08', 'r09', 'r10', 'r11', 'r12', 'r13', 'r14', 'r15', 'r16', 'r17', 'r18', 'r19', 'r20', 'r21', 'r22', 'r23', 'r24', 'r25', 'r26', 'r27', 'r28', 'r29', 'r30', 'r31', 'r32', 'r33', 'r34', 'r35', 'r36', 'r37', 'r38', 'r39','r40', 'r41', 'r42', 'r43'), std.lv = TRUE, auto.fix.first = FALSE, estimator = "DWLS",test =  "Satorra.Bentler")

summary(fit_gen_6, fit.measures = TRUE)

pred = predict(fit_gen_6)

pred_2 = as.data.frame(cbind(factor_all$individual, pred))
setnames(pred_2, 1, "IID")

write.table(pred_2, file = "~/Autism_heterogeneity/Factor_data/Combined_factor_predictive.txt", row.names = F, col.names = T, quote = F)



###step 2: run secondary analyses###
##Factor
factor = fread("~/Autism_heterogeneity/Factor_data/Combined_factor_predictive.txt")
factor2 = as.data.frame(scale(factor[,2:7]))
factor3 = as.data.frame(cbind(factor[,1], factor2))
factor = factor3

#sex_age
spark_basic = fread("~/SPARK/Phenotypes/V5/background_history_child.csv")
spark_basic_2 = spark_basic[,c("subject_sp_id", "sex", "age_at_eval_months", "cog_test_score")]

spark_basic_2$IQ_score = ifelse(spark_basic_2$cog_test_score == "24_below", 1, NA)
spark_basic_2$IQ_score = ifelse(spark_basic_2$cog_test_score == "25_39", 2, spark_basic_2$IQ_score)
spark_basic_2$IQ_score = ifelse(spark_basic_2$cog_test_score == "40_54", 3, spark_basic_2$IQ_score)
spark_basic_2$IQ_score = ifelse(spark_basic_2$cog_test_score == "55_69", 4, spark_basic_2$IQ_score)
spark_basic_2$IQ_score = ifelse(spark_basic_2$cog_test_score == "70_79", 5, spark_basic_2$IQ_score)
spark_basic_2$IQ_score = ifelse(spark_basic_2$cog_test_score == "80_89", 6, spark_basic_2$IQ_score)
spark_basic_2$IQ_score = ifelse(spark_basic_2$cog_test_score == "90_109", 7, spark_basic_2$IQ_score)
spark_basic_2$IQ_score = ifelse(spark_basic_2$cog_test_score == "110_119", 8, spark_basic_2$IQ_score)
spark_basic_2$IQ_score = ifelse(spark_basic_2$cog_test_score == "120_129", 9, spark_basic_2$IQ_score)
spark_basic_2$IQ_score = ifelse(spark_basic_2$cog_test_score == "130_above", 10, spark_basic_2$IQ_score)

setnames(spark_basic_2, old = c("subject_sp_id", "age_at_eval_months"), new = c("individual", "age"))
spark_basic_3 = spark_basic_2[,c("individual", "sex", "age", "IQ_score")]

ssc_basic = fread("~/SFARI/ssccore.txt")
ssc_basic_2 = ssc_basic[,c("individual", "sex", "age_at_ados", "ssc_diagnosis_full_scale_iq")]

ssc_basic_2$individual <- gsub('.p1', '', ssc_basic_2$individual)
ssc_basic_2$IQ_score = ifelse(ssc_basic_2$ssc_diagnosis_full_scale_iq < 24, 1, NA)
ssc_basic_2$IQ_score = ifelse(ssc_basic_2$ssc_diagnosis_full_scale_iq > 24 & ssc_basic_2$ssc_diagnosis_full_scale_iq  < 39,  2, ssc_basic_2$IQ_score)
ssc_basic_2$IQ_score = ifelse(ssc_basic_2$ssc_diagnosis_full_scale_iq > 38 & ssc_basic_2$ssc_diagnosis_full_scale_iq  < 54, 3, ssc_basic_2$IQ_score)
ssc_basic_2$IQ_score = ifelse(ssc_basic_2$ssc_diagnosis_full_scale_iq > 53 & ssc_basic_2$ssc_diagnosis_full_scale_iq  < 69, 4, ssc_basic_2$IQ_score)
ssc_basic_2$IQ_score = ifelse(ssc_basic_2$ssc_diagnosis_full_scale_iq > 69 & ssc_basic_2$ssc_diagnosis_full_scale_iq  < 80, 5, ssc_basic_2$IQ_score)
ssc_basic_2$IQ_score = ifelse(ssc_basic_2$ssc_diagnosis_full_scale_iq > 79 & ssc_basic_2$ssc_diagnosis_full_scale_iq  < 90, 6, ssc_basic_2$IQ_score)
ssc_basic_2$IQ_score = ifelse(ssc_basic_2$ssc_diagnosis_full_scale_iq > 89 & ssc_basic_2$ssc_diagnosis_full_scale_iq  < 110, 7, ssc_basic_2$IQ_score)
ssc_basic_2$IQ_score = ifelse(ssc_basic_2$ssc_diagnosis_full_scale_iq > 109 & ssc_basic_2$ssc_diagnosis_full_scale_iq  < 120, 8, ssc_basic_2$IQ_score)
ssc_basic_2$IQ_score = ifelse(ssc_basic_2$ssc_diagnosis_full_scale_iq > 119 & ssc_basic_2$ssc_diagnosis_full_scale_iq  < 130, 9, ssc_basic_2$IQ_score)
ssc_basic_2$IQ_score = ifelse(ssc_basic_2$ssc_diagnosis_full_scale_iq > 130, 10, ssc_basic_2$IQ_score)

ssc_basic_2$sex = ifelse(ssc_basic_2$sex == "male", "Male", "Female")

setnames(ssc_basic_2, old = c("age_at_ados"), new = c("age"))

ssc_basic_3 = ssc_basic_2[,c("individual", "sex", "age", "IQ_score")]

basic_all = rbind(ssc_basic_3, spark_basic_3)

all = merge(basic_all, factor, by.x = "individual", by.y = "IID")
all$PA6 = -1*all$PA6

library(stats)
library(ggplot2)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}


##Sex-based data

sex_all = all[,c(2, 5:10)]
d = melt(sex_all, id.vars="sex")
gamma <- summarySE(d, measurevar= "value", groupvars=c("variable", "sex"))

pd <- position_dodge(width = 0.3)
a = ggplot(gamma, aes(x = variable, y = value, fill = sex, colour = sex)) + 
  geom_point(size = 5, position=pd) + geom_errorbar(width = 0, aes(ymin = value-1.96*se, ymax = value+1.96*se), position = pd) +
  theme_classic() + ylab("Mean scores") + geom_hline(yintercept = 0) + theme(axis.text.x=element_text(size=rel(1.5), angle=90)) +  scale_colour_manual(values=c("#00798c", "#edae49"))

## IQ based data

IQ_all = all[,c(4:10)]
d = melt(IQ_all, id.vars="IQ_score")
setnames(d, "variable", "factor")
gamma <- summarySE(d, measurevar= "value", groupvars=c("factor", "IQ_score"))
gamma = na.omit(gamma)


pd <- position_dodge(width = 0.5)

b = ggplot(gamma, aes(x = as.factor(IQ_score), y =value, fill = factor, colour = factor)) + 
  geom_point(size = 5, position=pd) + geom_errorbar(width = 0, aes(ymin = value-1.96*se, ymax = value+1.96*se), position = pd) +
  theme_classic() + ylab("Mean scores") + geom_hline(yintercept = 0) + theme(axis.text.x=element_text(size=rel(1.5), angle=90)) +  scale_colour_manual(values=c("#00798c", "#edae49", "#2e4057", "#d1495b", "#8d96a3", "#66a182"))

##Age based data analysis

age_all = all[,c(3,5:10)]
d = melt(age_all, id.vars="age")
setnames(d, "variable", "factor")


c = ggplot(d, aes(x = age, y = value, fill = factor, colour = factor)) + geom_smooth() +
  theme_classic() + ylab("Mean scores") + scale_colour_manual(values=c("#00798c", "#edae49", "#2e4057", "#d1495b", "#8d96a3", "#66a182")) + scale_fill_manual(values=c("#00798c", "#edae49", "#2e4057", "#d1495b", "#8d96a3", "#66a182"))


##Age and sex

males = subset(all, sex == "Male")

age_all_male = males[,c(3,5:10)]
d = melt(age_all_male, id.vars="age")
setnames(d, "variable", "factor")


d_plot = ggplot(d, aes(x = age, y = value, fill = factor, colour = factor)) + geom_smooth() +
  theme_classic() + ylab("Mean scores") + scale_colour_manual(values=c("#00798c", "#edae49", "#2e4057", "#d1495b", "#8d96a3", "#66a182")) + scale_fill_manual(values=c("#00798c", "#edae49", "#2e4057", "#d1495b", "#8d96a3", "#66a182"))


females = subset(all, sex == "Female")

age_all_female = females[,c(3,5:10)]
d = melt(age_all_female, id.vars="age")
setnames(d, "variable", "factor")

e = ggplot(d, aes(x = age, y = value, fill = factor, colour = factor)) + geom_smooth() +
  theme_classic() + ylab("Mean scores") + scale_colour_manual(values=c("#00798c", "#edae49", "#2e4057", "#d1495b", "#8d96a3", "#66a182")) + scale_fill_manual(values=c("#00798c", "#edae49", "#2e4057", "#d1495b", "#8d96a3", "#66a182"))

##IQ and sex

IQ_sex_all = all[,c(2, 4:10)]
d = melt(IQ_sex_all, id.vars=c("IQ_score", "sex"))
setnames(d, "variable", "factor")
gamma <- summarySE(d, measurevar= "value", groupvars=c("IQ_score", "factor", "sex"))
gamma = na.omit(gamma)

ggplot(gamma, aes(x = as.factor(IQ_score), y =value, fill = sex, colour = sex)) + 
  geom_point(size = 5, position=pd) + geom_errorbar(width = 0, aes(ymin = value-1.96*se, ymax = value+1.96*se), position = pd) +
  theme_classic() + ylab("Mean scores") + geom_hline(yintercept = 0) + theme(axis.text.x=element_text(size=rel(1.5), angle=90)) +  scale_colour_manual(values=c("#00798c", "#edae49", "#2e4057", "#d1495b", "#8d96a3", "#66a182")) + facet_wrap(~ factor)

###Correlation plot
#Corplots from here: http://www.sthda.com/english/wiki/ggplot2-quick-correlation-matrix-heatmap-r-software-and-data-visualization


reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}



cormat <- round(cor(all[,c(5:10)]),2)
cormat <- reorder_cormat(cormat)

# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

upper_tri <- get_upper_tri(cormat)
melted_cormat <- melt(upper_tri, na.rm = TRUE)
corplot = ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "#00798c", high = "#d1495b", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed() +
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.6, 0.7),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))

#Arrange stuff from here: http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/81-ggplot2-easy-way-to-mix-multiple-graphs-on-the-same-page/
library(ggpubr)
ggarrange(ggarrange(corplot, a, labels = c("A", "B"), ncol = 2), b,c, nrow = 3, labels = c("C", "D"))

ggarrange(d_plot, e, labels = c("A", "B"), nrow = 2)

