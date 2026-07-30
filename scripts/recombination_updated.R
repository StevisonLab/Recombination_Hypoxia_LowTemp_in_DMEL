#Kosambi correction function to convert recombination frequency to map distance
kosambi_correction <- function(input) {
  ifelse(input < 0.5, (100/4)*(log((1+(2*input))/(1-(2*input)))), NA)
}

#reset treatment order for plotting
treat_order=c("Control","ColdTemp","Hypoxia","Combined")

#define haplotype classes
haplotypes=c("wt.wt.wt.wt","y.cv.v.f","y.wt.wt.wt","wt.cv.wt.wt","wt.wt.v.wt","wt.wt.wt.f",      "y.cv.wt.wt","y.cv.v.wt","y.wt.v.wt","y.wt.wt.f","wt.wt.v.f","wt.cv.v.f","wt.cv.v.wt" ,     "wt.cv.wt.f","y.wt.v.f","y.cv.wt.f")
haps_nco=c("wt.wt.wt.wt","mt.mt.mt.mt")
haps_sco=c("y.wt.wt.wt","wt.wt.wt.f","y.cv.wt.wt","y.cv.v.wt","wt.wt.v.f","wt.cv.v.f")
haps_dco=c("wt.cv.wt.wt","wt.wt.v.wt","y.wt.wt.f","wt.cv.v.wt","y.wt.v.f","y.cv.wt.f")
haps_tco=c("y.wt.v.wt","wt.cv.wt.f")

ycv_haps=c("y.wt.wt.wt","wt.cv.wt.wt","y.wt.v.wt","y.wt.wt.f","wt.cv.v.f","wt.cv.v.wt","wt.cv.wt.f","y.wt.v.f")
cvv_haps=c("wt.cv.wt.wt","wt.wt.v.wt","y.cv.wt.wt","y.wt.v.wt","wt.wt.v.f","wt.cv.wt.f","y.wt.v.f","y.cv.wt.f")
vf_haps=c("wt.wt.v.wt","wt.wt.wt.f","y.cv.v.wt","y.wt.v.wt","y.wt.wt.f","wt.cv.v.wt","wt.cv.wt.f","y.cv.wt.f")

#store by_vialday from merge script as new object called recomb
recomb <- by_vialday
recomb$Stock=as.factor(recomb$Stock)

recomb$nco_inds <- rowSums(recomb[names(recomb) %in% haps_nco], )
recomb$sco_inds <- rowSums(recomb[names(recomb) %in% haps_sco], )
recomb$dco_inds <- rowSums(recomb[names(recomb) %in% haps_dco], )
recomb$tco_inds <- rowSums(recomb[names(recomb) %in% haps_tco], )
recomb$co_inds <- recomb$sco_inds + 2*recomb$dco_inds + 3*recomb$tco_inds
recomb$num_M <- rowSums(recomb[names(recomb) %in% haplotypes], )

#totals:
sum(recomb$num_M)
sum(recomb$num_offspring)
length(unique(sort(recomb$F1_vial)))

recomb$recomb_rate <- recomb$co_inds / recomb$num_M
#recomb$kosambi <- kosambi_correction(recomb$recomb_rate)

recomb <- na.omit(recomb)

## Rough recombination rate by interval for validation purposes
recomb$ycv_count <-rowSums(recomb[names(recomb) %in% ycv_haps], )
recomb$cvv_count <-rowSums(recomb[names(recomb) %in% cvv_haps], )
recomb$vf_count <-rowSums(recomb[names(recomb) %in% vf_haps], )
recomb$mapl_count <-recomb$ycv_count+recomb$cvv_count+recomb$vf_count

recomb$ycv_rr <- recomb$ycv_count / recomb$num_M
recomb$cvv_rr <-recomb$cvv_count / recomb$num_M
recomb$vf_rr <-recomb$vf_count / recomb$num_M
recomb$ml_rr <- recomb$ycv_rr+recomb$cvv_rr+recomb$vf_rr

total_ycv_count <- sum(recomb[names(recomb) %in% ycv_haps])
total_cvv_count <- sum(recomb[names(recomb) %in% cvv_haps])
total_vf_count <- sum(recomb[names(recomb) %in% vf_haps])
total_map_count <- total_ycv_count + total_cvv_count + total_vf_count
total_num_M <- sum(recomb["num_M"])

total_ycv_rr <- total_ycv_count / total_num_M
total_cvv_rr <-total_cvv_count / total_num_M
total_vf_rr <-total_vf_count / total_num_M
total_map_length <- total_map_count / total_num_M

recomb_intervals <- data.frame(interval = c("y-cv", "cv-v", "v-f", "map_length"),
                               observed = c(total_ycv_rr, total_cvv_rr, total_vf_rr,total_map_length) * 100,
                               expected = c(13.7, 19.3, 23.7,56.7))

# Run COI script
source("scripts/COI.R")

## Remove haplotype columns
recomb <- recomb[,!names(recomb) %in% haplotypes]
recomb$oxygen=ifelse(recomb$Treatment=="Combined","Low",ifelse(recomb$Treatment=="Hypoxia","Low","Control"))
recomb$temperature=ifelse(recomb$Treatment=="Control","Control",ifelse(recomb$Treatment=="Hypoxia","Control","Low"))

# Create RR summaries with Kosambi Correction for Table in Paper:
recomb_rate=summaryBy(nco_inds+co_inds+ycv_count+cvv_count+vf_count+num_M~vial_num+vial_letter+Stock+Treatment+F1_vial,data=recomb,FUN=sum,na.rm=T)

recomb_rate=recomb_rate[recomb_rate$num_M.sum>=5,]
recomb_rate$ycv_rr=kosambi_correction(recomb_rate$ycv_count.sum/recomb_rate$num_M.sum)
recomb_rate$cvv_rr=kosambi_correction(recomb_rate$cvv_count.sum/recomb_rate$num_M.sum)
recomb_rate$vf_rr=kosambi_correction(recomb_rate$vf_count.sum/recomb_rate$num_M.sum)
recomb_rate$recomb_rate=rowSums(recomb_rate[,c("ycv_rr","cvv_rr","vf_rr")],na.rm=TRUE)

recomb_rate$Stock=as.factor(recomb_rate$Stock)
recomb_rate$oxygen=ifelse(recomb_rate$Treatment=="Combined","Low",ifelse(recomb_rate$Treatment=="Hypoxia","Low","Control"))
recomb_rate$temperature=ifelse(recomb_rate$Treatment=="Control","Control",ifelse(recomb_rate$Treatment=="Hypoxia","Control","Low"))

#Get summary:
recomb_summary <- aggregate(recomb_rate[, c("recomb_rate", "ycv_rr", "cvv_rr", "vf_rr")], by = list(recomb_rate$Treatment, recomb_rate$Stock), mean, na.rm=TRUE)
names(recomb_summary)[1:2] <- c("Treatment", "Stock")

recomb_summary$Stock=as.factor(recomb_summary$Stock)
recomb_summary$oxygen=ifelse(recomb_summary$Treatment=="Combined","Low",ifelse(recomb_summary$Treatment=="Hypoxia","Low","Control"))
recomb_summary$temperature=ifelse(recomb_summary$Treatment=="Control","Control",ifelse(recomb_summary$Treatment=="Hypoxia","Control","Low"))

#Write output
write.csv(recomb_summary, "output/recombination_summary_Kosambi.csv", row.names = FALSE)

#Repeat without Kosambi for models:
recomb_rate=summaryBy(nco_inds+co_inds+ycv_count+cvv_count+vf_count+num_M~vial_num+vial_letter+Stock+Treatment+F1_vial,data=recomb,FUN=sum,na.rm=T)

recomb_rate=recomb_rate[recomb_rate$num_M.sum>=5,]

recomb_rate$recomb_rate=recomb_rate$co_inds.sum/recomb_rate$num_M.sum
recomb_rate$ycv_rr=recomb_rate$ycv_count.sum/recomb_rate$num_M.sum
recomb_rate$cvv_rr=recomb_rate$cvv_count.sum/recomb_rate$num_M.sum
recomb_rate$vf_rr=recomb_rate$vf_count.sum/recomb_rate$num_M.sum
recomb_rate$recomb_rate=rowSums(recomb_rate[,c("ycv_rr","cvv_rr","vf_rr")],na.rm=TRUE)

recomb_rate$Stock=as.factor(recomb_rate$Stock)
recomb_rate$oxygen=ifelse(recomb_rate$Treatment=="Combined","Low",ifelse(recomb_rate$Treatment=="Hypoxia","Low","Control"))
recomb_rate$temperature=ifelse(recomb_rate$Treatment=="Control","Control",ifelse(recomb_rate$Treatment=="Hypoxia","Control","Low"))

####recomb_rate$ml_rr=recomb_rate$ycv_rr+recomb_rate$cvv_rr+recomb_rate$vf_rr
##recomb_summary <- aggregate(recomb_rate ~ Stock + vial_num+vial_letter+Treatment + F1_vial, data = recomb, sum)

#Get summary:
recomb_summary <- aggregate(recomb_rate[, c("recomb_rate", "ycv_rr", "cvv_rr", "vf_rr")], by = list(recomb_rate$Treatment, recomb_rate$Stock), mean, na.rm=TRUE)
names(recomb_summary)[1:2] <- c("Treatment", "Stock")

recomb_summary$Stock=as.factor(recomb_summary$Stock)
recomb_summary$oxygen=ifelse(recomb_summary$Treatment=="Combined","Low",ifelse(recomb_summary$Treatment=="Hypoxia","Low","Control"))
recomb_summary$temperature=ifelse(recomb_summary$Treatment=="Control","Control",ifelse(recomb_summary$Treatment=="Hypoxia","Control","Low"))

write.csv(recomb_summary, "output/recombination_summary.csv", row.names = FALSE)

#Go back and rerun so recomb_summary is raw counts and not Kosambi corrected!

## Stats and figures (NOTE: all below is fully copied and pasted from Novak et al 2026. Needs to be edited for hypoxia project!)

## FULL interval recomb model

#recomb_rate$nco_inds=recomb_rate$num_M.sum-recomb_rate$co_inds.sum
#NOTE: this gives negative numbers!
fit_full <- glmer(cbind(co_inds.sum,nco_inds.sum) ~ (1|F1_vial) + Treatment * Stock * vial_letter, data=recomb_rate,family=binomial(link="logit"),control = glmerControl(
  optimizer ='optimx', optCtrl=list(method='nlminb')))
anova_full <- Anova(fit_full, test = "Chisq")
anova_full

write.csv(anova_full, "output/recomb-stats_anova-full_comboTreatment.csv")

## Post hoc test for full
fit_contrast_full <- emmeans::emmeans(fit_full, specs = pairwise ~ Treatment, by=c("Stock","vial_letter"), mode="kenward-roger")
fit_contr_full <- contrast(fit_contrast_full, method="pairwise")

pheno_contr_full <- as.data.frame(summary(fit_contr_full))
pheno_contr_full$full_sig=ifelse(pheno_contr_full$p.value<0.001,"***",ifelse(pheno_contr_full$p.value<0.01,"**",ifelse(pheno_contr_full$p.value<0.05,"*","")))
pheno_contr_full

write.csv(pheno_contr_full, "output/recomb-stats_anova-full_posthoc.csv")

# OR Plots for Full Model
source("scripts/OR_Plots_Full_Model.R")

#Repeat with oxygen and temperature as two different parameters
fit_full_split <- glmer(cbind(co_inds.sum,nco_inds.sum) ~ (1|F1_vial) + oxygen * temperature * Stock * vial_letter, data=recomb_rate,family=binomial(link="logit"),control = glmerControl(
  optimizer ='optimx', optCtrl=list(method='nlminb')))
anova_full_split <- Anova(fit_full_split, test = "Chisq")
anova_full_split

write.csv(anova_full_split, "output/recomb-stats_anova-full_splitTreatment.csv")

## Post hoc test for full
fit_contrast_full <- emmeans::emmeans(fit_full_split, specs = pairwise ~ oxygen, by=c("Stock","vial_letter","temperature"), mode="kenward-roger")
fit_contr_full <- contrast(fit_contrast_full, method="pairwise")

pheno_contr_full_oxy <- as.data.frame(summary(fit_contr_full))
pheno_contr_full_oxy$full_sig=ifelse(pheno_contr_full_oxy$p.value<0.001,"***",ifelse(pheno_contr_full_oxy$p.value<0.01,"**",ifelse(pheno_contr_full_oxy$p.value<0.05,"*","")))
pheno_contr_full_oxy

write.csv(pheno_contr_full_oxy, "output/recomb-stats_anova-full_posthoc_oxygen.csv")

#Repeat for temperature
fit_contrast_full <- emmeans::emmeans(fit_full_split, specs = pairwise ~ temperature, by=c("Stock","vial_letter","oxygen"), mode="kenward-roger")
fit_contr_full <- contrast(fit_contrast_full, method="pairwise")

pheno_contr_full_temp <- as.data.frame(summary(fit_contr_full))
pheno_contr_full_temp$full_sig=ifelse(pheno_contr_full_temp$p.value<0.001,"***",ifelse(pheno_contr_full_temp$p.value<0.01,"**",ifelse(pheno_contr_full_temp$p.value<0.05,"*","")))
pheno_contr_full_temp

write.csv(pheno_contr_full_temp, "output/recomb-stats_anova-full_posthoc_temperature.csv")

# OR Plots for Full Factorial Split Model
source("scripts/GxE_Plots_Split_Model.R")
source("scripts/Treatment_Plots.R")

#Compare models:
AIC(fit_full, fit_full_split)  # Lower AIC is better
#Identical AIC, suggestion is to pick simpler model

# Interval Models ---------------------------------------------------------

## y-cv interval recomb model
fit_ycv <- glmer(cbind(ycv_count,nco_inds) ~ (1|F1_vial) + Treatment * Stock * vial_letter, data=recomb,family=binomial(link="logit"),control = glmerControl(
  optimizer ='optimx', optCtrl=list(method='nlminb')))
anova_ycv <- Anova(fit_ycv, test = "Chisq")
anova_ycv
write.csv(anova_ycv, "output/recomb-stats_anova-ycv.csv")

## Post hoc test for y-cv
fit_contrast_ycv <- emmeans::emmeans(fit_ycv, specs = pairwise ~ Treatment, by=c("Stock","vial_letter"), mode="kenward-roger")
fit_contr_ycv <- contrast(fit_contrast_ycv, method="pairwise")

pheno_contr_ycv <- as.data.frame(summary(fit_contr_ycv))
pheno_contr_ycv

write.csv(pheno_contr_ycv, "output/recomb-stats_anova-ycv_posthoc.csv")

## cv-v interval recomb model
fit_cvv <- glmer(cbind(cvv_count,nco_inds) ~ (1|F1_vial) + Treatment * Stock * vial_letter, data=recomb,family=binomial(link="logit"),control = glmerControl(
  optimizer ='optimx', optCtrl=list(method='nlminb')))
anova_cvv <- Anova(fit_cvv, test = "Chisq")
anova_cvv
write.csv(anova_cvv, "output/recomb-stats_anova-cvv.csv")

## Post hoc test for cv-v
fit_contrast_cvv <- emmeans::emmeans(fit_cvv, specs = pairwise ~ Treatment, by=c("Stock","vial_letter"), mode="kenward-roger")
fit_contr_cvv <- contrast(fit_contrast_cvv, method="pairwise")

pheno_contr_cvv <- as.data.frame(summary(fit_contr_cvv))
pheno_contr_cvv

write.csv(pheno_contr_cvv, "output/recomb-stats_anova-cvv_posthoc.csv")

## v-f interval recomb model
fit_vf <- glmer(cbind(vf_count,nco_inds) ~ (1|F1_vial) + Treatment * Stock * vial_letter, data=recomb,family=binomial(link="logit"),control = glmerControl(
  optimizer ='optimx', optCtrl=list(method='nlminb')))
anova_vf <- Anova(fit_vf, test = "Chisq")
anova_vf

write.csv(anova_vf, "output/recomb-stats_anova-vf.csv")

## Post hoc test for v-f
fit_contrast_vf <- emmeans::emmeans(fit_vf, specs = pairwise ~ Treatment, by=c("Stock","vial_letter"), mode="kenward-roger")
fit_contr_vf <- contrast(fit_contrast_vf, method="pairwise")

pheno_contr_vf <- as.data.frame(summary(fit_contr_vf))
pheno_contr_vf

write.csv(pheno_contr_vf, "output/recomb-stats_anova-vf_posthoc.csv")

source("scripts/Per_Interval_Map_Plot.R")

# Summary Values for Paper (OLD): -----------------------------------------------

recomb_rate=summaryBy(co_inds+ycv_count+cvv_count+vf_count+num_M~vial_letter+Stock+Treatment,data=recomb,FUN=sum,na.rm=T)
recomb_rate=recomb_rate[recomb_rate$num_M.sum>=5,]
recomb_rate$ycv_rr=kosambi_correction(recomb_rate$ycv_count.sum/recomb_rate$num_M.sum)
recomb_rate$cvv_rr=kosambi_correction(recomb_rate$cvv_count.sum/recomb_rate$num_M.sum)
recomb_rate$vf_rr=kosambi_correction(recomb_rate$vf_count.sum/recomb_rate$num_M.sum)
recomb_rate$recomb_rate=rowSums(recomb_rate[,c("ycv_rr","cvv_rr","vf_rr")],na.rm=TRUE)
recomb_summary <- aggregate(recomb_rate[, c("recomb_rate", "ycv_rr", "cvv_rr", "vf_rr")], by = list(recomb_rate$Treatment, recomb_rate$Stock), mean)
names(recomb_summary)[1:2] <- c("Treatment", "Stock")

#Average Treatment Differences
mean(recomb_summary$recomb_rate[recomb_summary$Treatment=="Hypoxia"]-recomb_summary$recomb_rate[recomb_summary$Treatment=="Control"])
#-0.916243 Control vs. Hypoxia
mean(recomb_summary$recomb_rate[recomb_summary$Treatment=="ColdTemp"]-recomb_summary$recomb_rate[recomb_summary$Treatment=="Control"])
#4.315604 Control vs. ColdTemp
mean(recomb_summary$recomb_rate[recomb_summary$Treatment=="Combined"]-recomb_summary$recomb_rate[recomb_summary$Treatment=="Control"])
#6.13872 Control vs. Combined

# Differences in recomb rate between diets for strain 57
sum(recomb_summary$ycv_rr[recomb_summary$Stock=="57" & recomb_summary$Treatment=="Control"]+recomb_summary$cvv_rr[recomb_summary$Stock=="57" & recomb_summary$Treatment=="Control"]+recomb_summary$vf_rr[recomb_summary$Stock=="57" & recomb_summary$Treatment=="Control"])
#51.09812 #control
sum(recomb_summary$ycv_rr[recomb_summary$Stock=="57" & recomb_summary$Treatment=="Hypoxia"]+recomb_summary$cvv_rr[recomb_summary$Stock=="57" & recomb_summary$Treatment=="Hypoxia"]+recomb_summary$vf_rr[recomb_summary$Stock=="57" & recomb_summary$Treatment=="Hypoxia"])
#47.22197 Hypoxia
sum(recomb_summary$ycv_rr[recomb_summary$Stock=="57" & recomb_summary$Treatment=="ColdTemp"]+recomb_summary$cvv_rr[recomb_summary$Stock=="57" & recomb_summary$Treatment=="ColdTemp"]+recomb_summary$vf_rr[recomb_summary$Stock=="57" & recomb_summary$Treatment=="ColdTemp"])
#59.08704 ColdTemp
sum(recomb_summary$ycv_rr[recomb_summary$Stock=="57" & recomb_summary$Treatment=="Combined"]+recomb_summary$cvv_rr[recomb_summary$Stock=="57" & recomb_summary$Treatment=="Combined"]+recomb_summary$vf_rr[recomb_summary$Stock=="57" & recomb_summary$Treatment=="Combined"])
#60.81918 Combined

#Treatment Differences from Control
recomb_summary$recomb_rate[recomb_summary$Stock=="57" & recomb_summary$Treatment=="Hypoxia"]-recomb_summary$recomb_rate[recomb_summary$Stock=="57" & recomb_summary$Treatment=="Control"]
#-3.876155 Control vs. Hypoxia
recomb_summary$recomb_rate[recomb_summary$Stock=="57" & recomb_summary$Treatment=="ColdTemp"]-recomb_summary$recomb_rate[recomb_summary$Stock=="57" & recomb_summary$Treatment=="Control"]
#7.988914 Control vs. ColdTemp
recomb_summary$recomb_rate[recomb_summary$Stock=="57" & recomb_summary$Treatment=="Combined"]-recomb_summary$recomb_rate[recomb_summary$Stock=="57" & recomb_summary$Treatment=="Control"]
#9.721051 Control vs. Combined


# Differences in recomb rate between diets for strain 508
sum(recomb_summary$ycv_rr[recomb_summary$Stock=="508" & recomb_summary$Treatment=="Control"]+recomb_summary$cvv_rr[recomb_summary$Stock=="508" & recomb_summary$Treatment=="Control"]+recomb_summary$vf_rr[recomb_summary$Stock=="508" & recomb_summary$Treatment=="Control"])
#48.53728 #control
sum(recomb_summary$ycv_rr[recomb_summary$Stock=="508" & recomb_summary$Treatment=="Hypoxia"]+recomb_summary$cvv_rr[recomb_summary$Stock=="508" & recomb_summary$Treatment=="Hypoxia"]+recomb_summary$vf_rr[recomb_summary$Stock=="508" & recomb_summary$Treatment=="Hypoxia"])
#54.18651 Hypoxia
sum(recomb_summary$ycv_rr[recomb_summary$Stock=="508" & recomb_summary$Treatment=="ColdTemp"]+recomb_summary$cvv_rr[recomb_summary$Stock=="508" & recomb_summary$Treatment=="ColdTemp"]+recomb_summary$vf_rr[recomb_summary$Stock=="508" & recomb_summary$Treatment=="ColdTemp"])
#53.6045 ColdTemp
sum(recomb_summary$ycv_rr[recomb_summary$Stock=="508" & recomb_summary$Treatment=="Combined"]+recomb_summary$cvv_rr[recomb_summary$Stock=="508" & recomb_summary$Treatment=="Combined"]+recomb_summary$vf_rr[recomb_summary$Stock=="508" & recomb_summary$Treatment=="Combined"])
#58.19794 Combined

#Treatment Differences from Control
recomb_summary$recomb_rate[recomb_summary$Stock=="508" & recomb_summary$Treatment=="Hypoxia"]-recomb_summary$recomb_rate[recomb_summary$Stock=="508" & recomb_summary$Treatment=="Control"]
#5.649229 Control vs. Hypoxia
recomb_summary$recomb_rate[recomb_summary$Stock=="508" & recomb_summary$Treatment=="ColdTemp"]-recomb_summary$recomb_rate[recomb_summary$Stock=="508" & recomb_summary$Treatment=="Control"]
#5.067214 Control vs. ColdTemp
recomb_summary$recomb_rate[recomb_summary$Stock=="508" & recomb_summary$Treatment=="Combined"]-recomb_summary$recomb_rate[recomb_summary$Stock=="508" & recomb_summary$Treatment=="Control"]
#9.660652 Control vs. Combined



# Differences in recomb rate between diets for strain 391
sum(recomb_summary$ycv_rr[recomb_summary$Stock=="391" & recomb_summary$Treatment=="Control"]+recomb_summary$cvv_rr[recomb_summary$Stock=="391" & recomb_summary$Treatment=="Control"]+recomb_summary$vf_rr[recomb_summary$Stock=="391" & recomb_summary$Treatment=="Control"])
#55.77273 #control
sum(recomb_summary$ycv_rr[recomb_summary$Stock=="391" & recomb_summary$Treatment=="Hypoxia"]+recomb_summary$cvv_rr[recomb_summary$Stock=="391" & recomb_summary$Treatment=="Hypoxia"]+recomb_summary$vf_rr[recomb_summary$Stock=="391" & recomb_summary$Treatment=="Hypoxia"])
#55.73224 Hypoxia
sum(recomb_summary$ycv_rr[recomb_summary$Stock=="391" & recomb_summary$Treatment=="ColdTemp"]+recomb_summary$cvv_rr[recomb_summary$Stock=="391" & recomb_summary$Treatment=="ColdTemp"]+recomb_summary$vf_rr[recomb_summary$Stock=="391" & recomb_summary$Treatment=="ColdTemp"])
#61.64082 ColdTemp
sum(recomb_summary$ycv_rr[recomb_summary$Stock=="391" & recomb_summary$Treatment=="Combined"]+recomb_summary$cvv_rr[recomb_summary$Stock=="391" & recomb_summary$Treatment=="Combined"]+recomb_summary$vf_rr[recomb_summary$Stock=="391" & recomb_summary$Treatment=="Combined"])
#60.27328 Combined

#Treatment Differences from Control
recomb_summary$recomb_rate[recomb_summary$Stock=="391" & recomb_summary$Treatment=="Hypoxia"]-recomb_summary$recomb_rate[recomb_summary$Stock=="391" & recomb_summary$Treatment=="Control"]
#-0.04048637 Control vs. Hypoxia
recomb_summary$recomb_rate[recomb_summary$Stock=="391" & recomb_summary$Treatment=="ColdTemp"]-recomb_summary$recomb_rate[recomb_summary$Stock=="391" & recomb_summary$Treatment=="Control"]
#5.868091 Control vs. ColdTemp
recomb_summary$recomb_rate[recomb_summary$Stock=="391" & recomb_summary$Treatment=="Combined"]-recomb_summary$recomb_rate[recomb_summary$Stock=="391" & recomb_summary$Treatment=="Control"]
#4.500559 Control vs. Combined


# Differences in recomb rate between diets for strain 491
sum(recomb_summary$ycv_rr[recomb_summary$Stock=="491" & recomb_summary$Treatment=="Control"]+recomb_summary$cvv_rr[recomb_summary$Stock=="491" & recomb_summary$Treatment=="Control"]+recomb_summary$vf_rr[recomb_summary$Stock=="491" & recomb_summary$Treatment=="Control"])
#56.8071 #control
sum(recomb_summary$ycv_rr[recomb_summary$Stock=="491" & recomb_summary$Treatment=="Hypoxia"]+recomb_summary$cvv_rr[recomb_summary$Stock=="491" & recomb_summary$Treatment=="Hypoxia"]+recomb_summary$vf_rr[recomb_summary$Stock=="491" & recomb_summary$Treatment=="Hypoxia"])
#50.57746 Hypoxia
sum(recomb_summary$ycv_rr[recomb_summary$Stock=="491" & recomb_summary$Treatment=="ColdTemp"]+recomb_summary$cvv_rr[recomb_summary$Stock=="491" & recomb_summary$Treatment=="ColdTemp"]+recomb_summary$vf_rr[recomb_summary$Stock=="491" & recomb_summary$Treatment=="ColdTemp"])
#61.1726 ColdTemp
sum(recomb_summary$ycv_rr[recomb_summary$Stock=="491" & recomb_summary$Treatment=="Combined"]+recomb_summary$cvv_rr[recomb_summary$Stock=="491" & recomb_summary$Treatment=="Combined"]+recomb_summary$vf_rr[recomb_summary$Stock=="491" & recomb_summary$Treatment=="Combined"])
#63.1126 Combined

#Treatment Differences from Control
recomb_summary$recomb_rate[recomb_summary$Stock=="491" & recomb_summary$Treatment=="Hypoxia"]-recomb_summary$recomb_rate[recomb_summary$Stock=="491" & recomb_summary$Treatment=="Control"]
#-6.229642 Control vs. Hypoxia
recomb_summary$recomb_rate[recomb_summary$Stock=="491" & recomb_summary$Treatment=="ColdTemp"]-recomb_summary$recomb_rate[recomb_summary$Stock=="491" & recomb_summary$Treatment=="Control"]
#4.365505 Control vs. ColdTemp
recomb_summary$recomb_rate[recomb_summary$Stock=="491" & recomb_summary$Treatment=="Combined"]-recomb_summary$recomb_rate[recomb_summary$Stock=="491" & recomb_summary$Treatment=="Control"]
#6.305502 Control vs. Combined


# Differences in recomb rate between diets for strain 42
sum(recomb_summary$ycv_rr[recomb_summary$Stock=="42" & recomb_summary$Treatment=="Control"]+recomb_summary$cvv_rr[recomb_summary$Stock=="42" & recomb_summary$Treatment=="Control"]+recomb_summary$vf_rr[recomb_summary$Stock=="42" & recomb_summary$Treatment=="Control"])
#51.11953 #control
sum(recomb_summary$ycv_rr[recomb_summary$Stock=="42" & recomb_summary$Treatment=="Hypoxia"]+recomb_summary$cvv_rr[recomb_summary$Stock=="42" & recomb_summary$Treatment=="Hypoxia"]+recomb_summary$vf_rr[recomb_summary$Stock=="42" & recomb_summary$Treatment=="Hypoxia"])
#51.03537 Hypoxia
sum(recomb_summary$ycv_rr[recomb_summary$Stock=="42" & recomb_summary$Treatment=="ColdTemp"]+recomb_summary$cvv_rr[recomb_summary$Stock=="42" & recomb_summary$Treatment=="ColdTemp"]+recomb_summary$vf_rr[recomb_summary$Stock=="42" & recomb_summary$Treatment=="ColdTemp"])
#49.40783 ColdTemp
sum(recomb_summary$ycv_rr[recomb_summary$Stock=="42" & recomb_summary$Treatment=="Combined"]+recomb_summary$cvv_rr[recomb_summary$Stock=="42" & recomb_summary$Treatment=="Combined"]+recomb_summary$vf_rr[recomb_summary$Stock=="42" & recomb_summary$Treatment=="Combined"])
#51.62537 Combined

#Treatment Differences from Control
recomb_summary$recomb_rate[recomb_summary$Stock=="42" & recomb_summary$Treatment=="Hypoxia"]-recomb_summary$recomb_rate[recomb_summary$Stock=="42" & recomb_summary$Treatment=="Control"]
#-0.08416133 Control vs. Hypoxia
recomb_summary$recomb_rate[recomb_summary$Stock=="42" & recomb_summary$Treatment=="ColdTemp"]-recomb_summary$recomb_rate[recomb_summary$Stock=="42" & recomb_summary$Treatment=="Control"]
#-1.711707 Control vs. ColdTemp
recomb_summary$recomb_rate[recomb_summary$Stock=="42" & recomb_summary$Treatment=="Combined"]-recomb_summary$recomb_rate[recomb_summary$Stock=="42" & recomb_summary$Treatment=="Control"]
#0.5058342 Control vs. Combined
