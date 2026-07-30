
# Copy data from:
#https://github.com/StevisonLab/Combination-effects-of-hypoxia-and-low-temperature-in-Drosophila/blob/main/rawdata/merged.csv

#Note: merged data from previous study
merged=read.csv(file="rawdata/merged.csv")

# Anonymize observer columm
merged$observer <- as.numeric(as.factor(merged$observer))

#define haplotype classes
haplotypes=c("wt.wt.wt.wt","y.cv.v.f","y.wt.wt.wt","wt.cv.wt.wt","wt.wt.v.wt","wt.wt.wt.f",      "y.cv.wt.wt","y.cv.v.wt","y.wt.v.wt","y.wt.wt.f","wt.wt.v.f","wt.cv.v.f","wt.cv.v.wt" ,     "wt.cv.wt.f","y.wt.v.f","y.cv.wt.f")
haps_nco=c("wt.wt.wt.wt","mt.mt.mt.mt")
haps_sco=c("y.wt.wt.wt","wt.wt.wt.f","y.cv.wt.wt","y.cv.v.wt","wt.wt.v.f","wt.cv.v.f")
haps_dco=c("wt.cv.wt.wt","wt.wt.v.wt","y.wt.wt.f","wt.cv.v.wt","y.wt.v.f","y.cv.wt.f")
haps_tco=c("y.wt.v.wt","wt.cv.wt.f")

#Clean-up based on potential observer bias as in Novak et al. 2026:
#Check observer triple crossover rates
obs_check=aggregate(merged[, c(haplotypes, "num_males")], by= list(merged$observer), sum)
obs_check$TCO_rate=(obs_check$y.wt.v.wt+obs_check$wt.cv.wt.f)/obs_check$num_males
hist(obs_check$TCO_rate)

#remove small sets
##obs_check=subset(obs_check,obs_check$num_males>50)
#this will happen later

# Levene's Test for Homogeneity of Variance across dataset:
leveneTest(num_offspring~interaction(stock, treatment,day),merged)
# p=1e-06
# Significant heterogeneity of variance detected!!!

#Similar but this time across observers:
leveneTest(num_males~interaction(as.factor(observer),as.factor(F1_cross_vial_num)),merged)
#p=9.469e-10

#TCO rate estimates
mean_tco=mean(obs_check$TCO_rate)
sd_tco=sd(obs_check$TCO_rate)

md=mad(obs_check$TCO_rate)
md

iqr=IQR(obs_check$TCO_rate,na.rm = TRUE)
iqr

# Expected TCO rate across dataset using Comeron Map Distances
exp_tco=(0.12 * 0.21 * 0.24)
exp_tco

#Get list above md/iqr:
obs_check$Group.1[obs_check$TCO_rate>mean_tco+md]
obs_check$Group.1[obs_check$TCO_rate>mean_tco+iqr]

# Observers removed from dataset were:
## Individual #6 did not observe >50 males in the total dataset
## Individual #8 had an excess of TCO events as compared to the rest of the dataset and far beyond expectation

# difference from exp TCO
obs_check$TCO_rate[obs_check$Group.1=="8"]/exp_tco
#15.35577, more than an order of magnitude above expected TCO rate

#Obs #8 was over an order of magnitude above the Expected TCO rate

# Mean TCO without removed observers:
exp_tco
mean_tco
adj_mean=mean(obs_check$TCO_rate[obs_check$Group.1!="6" & obs_check$Group.1!="8"])
adj_mean

# Make plot for supplement: 
#Note: this relies on running the adj_mean code below!
tco_check=ggplot(data=obs_check,aes(x=Group.1,y=TCO_rate))+geom_point()+geom_hline(yintercept = mean_tco ,col="red")+geom_hline(yintercept = exp_tco,col="purple")+geom_hline(yintercept = mean_tco+md,col="grey50")+geom_hline(yintercept = mean_tco+iqr,col="grey30")+geom_hline(yintercept = adj_mean,col="orange")+xlab("Anonymized observer ID")+ylab("Triple Crossover Rate")+geom_text(aes(label=ifelse(TCO_rate>mean_tco+md,as.character(Group.1),ifelse(num_males<50,as.character(Group.1),''))),hjust=0,vjust=0)
tco_check
ggsave("images/Observer_Bias_TCO.png")

#Now remove the observer from the merged dataset
merged=subset(merged,merged$observer!="8")
merged=subset(merged,merged$observer!="6")

# Repeat Levene's Test for Homogeneity of Variance across subsetted dataset:
leveneTest(num_males~interaction(as.factor(observer),as.factor(F1_cross_vial_num)),merged)
#p=0.4531
#No longer detectable heterogeneity across the dataset

#quick summaries
haplotypes=c("wt.wt.wt.wt","y.cv.v.f","y.wt.wt.wt","wt.cv.wt.wt","wt.wt.v.wt","wt.wt.wt.f",      "y.cv.wt.wt","y.cv.v.wt","y.wt.v.wt","y.wt.wt.f","wt.wt.v.f","wt.cv.v.f","wt.cv.v.wt" ,     "wt.cv.wt.f","y.wt.v.f","y.cv.wt.f")
by_treatment=aggregate(merged[, c(haplotypes,"num_offspring")], by =list(merged$treatment,merged$day,merged$stock,merged$F1_cross_vial_num),sum)

by_vialday=aggregate(merged[, c(haplotypes,"num_offspring", "num_females")], by =list(merged$vial,merged$day,merged$stock,merged$treatment,merged$F1_cross_vial_num),sum)
names(by_vialday)[1:5] <- c("vial_num", "vial_letter","Stock","Treatment","F1_vial")

write_csv(by_vialday, "rawdata/byvialday.csv")

by_vial=aggregate(merged[, c(haplotypes, "num_offspring")], by= list(merged$vial_f1col_num), sum)

count_unique <- function(x) length(unique(x))
bc_summary = summaryBy(F1_cross_vial_num~oxygen+temperature+stock,data=merged,FUN = count_unique)
crosses_summary = summaryBy(vial~oxygen+temperature+stock,data=merged,FUN = count_unique)
