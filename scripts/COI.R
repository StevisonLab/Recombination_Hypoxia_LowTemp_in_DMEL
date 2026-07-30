
## Crossover Interference between the y-cv Interval and the cv-v Interval

### Crossover interference calculation

COI <- recomb
COI=summaryBy(co_inds+ycv_count+cvv_count+vf_count+num_M+wt.cv.wt.wt+wt.wt.v.wt+y.wt.wt.f+wt.cv.v.wt+y.wt.v.f+y.cv.wt.f+y.wt.v.wt+wt.cv.wt.f~Stock+Treatment+vial_letter,data=recomb,FUN=sum,na.rm=T)
COI$COI=COI$co_inds.sum/COI$num_M.sum
COI$ycv_rr=COI$ycv_count.sum/COI$num_M.sum
COI$cvv_rr=COI$cvv_count.sum/COI$num_M.sum
COI$vf_rr=COI$vf_count.sum/COI$num_M.sum
COI=COI[COI$num_M.sum>=5,]

# expected vs observed
COI$Exp_DCO_ycv_cvv <-  round((COI$ycv_rr * COI$cvv_rr)*COI$num_M,0)
COI$Obs_DCO_ycv_cvv <-  COI$wt.cv.wt.wt + COI$y.wt.v.f

COI$Exp_DCO_cvv_vf <-  round((COI$cvv_rr * COI$vf_rr)*COI$num_M,0)
COI$Obs_DCO_cvv_vf <-  COI$wt.wt.v.wt + COI$y.cv.wt.f

COI$Exp_TCO_ycv_vf <-  round((COI$ycv_rr * COI$cvv_rr * COI$vf_rr)*COI$num_M,0)
COI$Obs_TCO_ycv_vf <- COI$y.wt.v.wt + COI$wt.cv.wt.f

# coefficient of coincidence and interference
COI$COC_ycv_cvv <- COI$Obs_DCO_ycv_cvv / COI$Exp_DCO_ycv_cvv
COI$Interference_ycv_cvv <- 1 - COI$COC_ycv_cvv

COI$COC_cvv_vf <- COI$Obs_DCO_cvv_vf / COI$Exp_DCO_cvv_vf
COI$Interference_cvv_vf <- 1 - COI$COC_cvv_vf

COI$COC_ycv_vf <- COI$Obs_TCO_ycv_vf / COI$Exp_TCO_ycv_vf
COI$Interference_ycv_vf <- 1 - COI$COC_ycv_vf

hist(COI$Interference_ycv_cvv)
hist(COI$Interference_ycv_vf)
hist(COI$Interference_cvv_vf)

### Statistical models
COI$Stock=as.factor(COI$Stock)
COI$Treatment=as.factor(COI$Treatment)
#ycv_cvv
#remove expected zero
COI_test <- COI %>%
  filter(Exp_DCO_ycv_cvv > 0 )
fit_ycv_cvv <- glm(cbind(Obs_DCO_ycv_cvv) ~ Treatment * Stock + offset(log(Exp_DCO_ycv_cvv)), data=COI_test,family=quasipoisson(link = "log"))
#summary(fit_ycv_cvv)

anova_coi_ycv_cvv <- anova(fit_ycv_cvv, test="F")
anova_coi_ycv_cvv
#Significant interaction detected
write.csv(anova_coi_ycv_cvv,"output/interference_ycv_cvv_model_table.csv")

fit_contrast_ycv_cvv <- emmeans::emmeans(fit_ycv_cvv, "Treatment", by="Stock", mode="kenward-roger")
fit_contr_ycv_cvv <- contrast(fit_contrast_ycv_cvv, method="trt.vs.ctrl")

pheno_contr_ycv_cvv <- as.data.frame(summary(fit_contr_ycv_cvv))
pheno_contr_ycv_cvv
write.csv(pheno_contr_ycv_cvv,"output/interference_ycv_cvv_posthoc_table.csv")

#ycv_vf
#remove expected zero
COI_test <- COI %>%
  filter(Exp_TCO_ycv_vf > 0 )
fit_ycv_vf <- glm(Obs_TCO_ycv_vf ~ Treatment * Stock + offset(log(Exp_TCO_ycv_vf)), data=COI_test,family=quasipoisson(link = "log"))
#summary(fit_ycv_vf)

anova_coi_ycv_vf <- anova(fit_ycv_vf, test="F")
anova_coi_ycv_vf
#No significant interaction 
write.csv(anova_coi_ycv_vf,"output/interference_ycv_vf_model_table.csv")

#fit_contrast_ycv_vf <- emmeans::emmeans(fit_ycv_vf, "Treatment", by="Stock", mode="kenward-roger")
#fit_contr_ycv_vf <- contrast(fit_contrast_ycv_vf, method="trt.vs.ctrl")

#pheno_contr_ycv_vf <- as.data.frame(summary(fit_contr_ycv_vf))
#pheno_contr_ycv_vf
#write.csv(pheno_contr_ycv_vf,"output/interference_ycv_vf_posthoc_table.csv")

#cvv_vf
#remove expected zero
COI_test <- COI %>%
  filter(Exp_DCO_cvv_vf > 0 )
fit_cvv_vf <- glm(Obs_DCO_cvv_vf ~ Treatment * Stock + offset(log(Exp_DCO_cvv_vf)), data=COI_test,family=quasipoisson(link = "log"))
#summary(fit_cvv_vf)

anova_coi_cvv_vf <- anova(fit_cvv_vf, test="F")
anova_coi_cvv_vf
#No significant interaction
write.csv(anova_coi_cvv_vf,"output/interference_cvv_vf_model_table.csv")

# fit_contrast_cvv_vf <- emmeans::emmeans(fit_cvv_vf, "Treatment", by="Stock", mode="kenward-roger")
# fit_contr_cvv_vf <- contrast(fit_contrast_cvv_vf, method="trt.vs.ctrl")
# 
# pheno_contr_cvv_vf <- as.data.frame(summary(fit_contr_cvv_vf))
# pheno_contr_cvv_vf
# write.csv(pheno_contr_cvv_vf,"output/interference_cvv_vf_posthoc_table.csv")

### Figures 
stock_order <- c("508","42","57","391","491")
col_scheme2 <- c("#1b9e77","#d95f02","#7570b3","#66a61e","#e7298a")

#Get numeric treatment scale:
COI$treat_order=ifelse(COI$Treatment=="Control",1,ifelse(COI$Treatment=="Hypoxia",2,ifelse(COI$Treatment=="ColdTemp",3,4)))

COI_figure_ycv_cvv <- ggplot(aes(y=Interference_ycv_cvv, x=as.numeric(treat_order), col=factor(Stock,levels = stock_order)), data=COI) + 
  xlab("Treatment") + ylab("Crossover Interference") +
  scale_x_continuous(limits=c(0.9,4.1),breaks=c(1,2,3,4),labels=c("Control","Hypoxia","ColdTemp","Combined")) + ylim(-1.5,1.05) + 
  ggtitle("Interference of y-cv-v interval") + 
  theme_base() + scale_color_manual(values = col_scheme2)+theme(legend.title = element_blank(),plot.title = element_text(size = 8)) +
    #geom_point(size=1.5,na.rm=TRUE) + 
  stat_summary(fun.data = mean_se, geom = "errorbar", width=0.1) +
  stat_summary(fun = median, geom="line", aes(group=Stock),na.rm=TRUE, linewidth=1) 

COI_figure_ycv_cvv

COI_figure_ycv_vf <- ggplot(aes(y=Interference_ycv_vf, x=as.numeric(treat_order), col=factor(Stock,levels = stock_order)), data=COI) + 
  xlab("Treatment") + ylab("Crossover Interference") +
  scale_x_continuous(limits=c(0.9,4.1),breaks=c(1,2,3,4),labels=c("Control","Hypoxia","ColdTemp","Combined")) + ylim(-1.5,1.05) + 
  ggtitle("Interference of y-cv-v-f interval") + 
  theme_base() + scale_color_manual(values = col_scheme2)+
  #geom_point(size=1.5,na.rm=TRUE) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width=0.1) +
  stat_summary(fun = median, geom="line", aes(group=Stock),na.rm=TRUE, linewidth=1) +
  theme(legend.title = element_blank(),plot.title = element_text(size = 8)) 

COI_figure_ycv_vf

COI_figure_cvv_vf <- ggplot(aes(y=Interference_cvv_vf, x=as.numeric(treat_order), col=factor(Stock,levels = stock_order)), data=COI) + 
  xlab("Treatment") + ylab("Crossover Interference") +
  scale_x_continuous(limits=c(0.9,4.1),breaks=c(1,2,3,4),labels=c("Control","Hypoxia","ColdTemp","Combined")) + ylim(-1.5,1.05) + 
  ggtitle("Interference of cv-v-f interval") + 
  theme_base() + scale_color_manual(values = col_scheme2)+
  #geom_point(size=1.5,na.rm=TRUE) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width=0.1) +
  stat_summary(fun = median, geom="line", aes(group=Stock),na.rm=TRUE, linewidth=1) +
  theme(legend.title = element_blank(),plot.title = element_text(size = 8)) 

COI_figure_cvv_vf

### Save figures
ggsave("images/COI-ycv-cvv.png", plot=COI_figure_ycv_cvv, height=5)
ggsave("images/COI-cvv-vf.png", plot=COI_figure_cvv_vf, height=5)
ggsave("images/COI-ycv-vf.png", plot=COI_figure_ycv_vf, height=5)

# Plotting interference against genetic map distance

# Need to reorganize the data frame: Interference measure and Genetic Distance Measure
wide_df=COI[,c("Stock","Treatment","vial_letter","Interference_cvv_vf","Interference_ycv_cvv","Interference_ycv_vf")]
colnames(wide_df)=c("Stock","Treatment","vial_letter","cvv_vf","ycv_cvv","ycv_vf")
long_df=melt(wide_df,id.vars=c("Stock","Treatment","vial_letter"),variable.name = "Interval",value.name = "Interference")

wide_df2=COI[,c("Stock","Treatment","vial_letter","ycv_count.sum","cvv_count.sum","vf_count.sum","co_inds.sum")]
wide_df2$ycv_cvv=rowSums(wide_df2[,c("ycv_count.sum","cvv_count.sum")],na.rm=TRUE)
wide_df2$cvv_vf=rowSums(wide_df2[,c("cvv_count.sum","vf_count.sum")],na.rm=TRUE)
wide_df2$ycv_vf=wide_df2$co_inds.sum
wide_df2=wide_df2[,c("Stock","Treatment","vial_letter","cvv_vf","ycv_cvv","ycv_vf")]
long_df2=melt(wide_df2,id.vars=c("Stock","Treatment","vial_letter"),variable.name = "Interval",value.name = "Num_COs")

wide_df3=COI[,c("Stock","Treatment","vial_letter","num_M.sum")]
wide_df3$sub1=wide_df3$num_M.sum
wide_df3$sub2=wide_df3$num_M.sum
long_df3=melt(wide_df3,id.vars=c("Stock","Treatment","vial_letter"),variable.name = "Interval",value.name = "Total")

wide_df5=COI[,c("Stock","Treatment","vial_letter","Exp_DCO_cvv_vf","Exp_DCO_ycv_cvv","Exp_TCO_ycv_vf")]
colnames(wide_df5)=c("Stock","Treatment","vial_letter","cvv_vf","ycv_cvv","ycv_vf")
long_df5=melt(wide_df5,id.vars=c("Stock","Treatment","vial_letter"),variable.name = "Interval",value.name = "Exp_DCO")

wide_df6=COI[,c("Stock","Treatment","vial_letter","Obs_DCO_cvv_vf","Obs_DCO_ycv_cvv","Obs_TCO_ycv_vf")]
colnames(wide_df6)=c("Stock","Treatment","vial_letter","cvv_vf","ycv_cvv","ycv_vf")
long_df6=melt(wide_df6,id.vars=c("Stock","Treatment","vial_letter"),variable.name = "Interval",value.name = "Obs_DCO")

long_df4=cbind(long_df,long_df2[,5],long_df3[,5],long_df5[,5],long_df6[,5])
colnames(long_df4)=c("Stock","Treatment","vial_letter","Interval","Interference","Num_COs","Total","Exp_DCO","Obs_DCO")
long_df4$Gendist=long_df4$Num_COs/long_df4$Total
long_df4$Gendist=ifelse(long_df4$Num_COs/long_df4$Total<0.5,long_df4$Num_COs/long_df4$Total,0.499+(((long_df4$Num_COs/long_df4$Total)-0.5)/100))
long_df4$kosambi_distances <- (100/4)*(log((1+(2*long_df4$Gendist))/(1-(2*long_df4$Gendist))))
long_df4$Stock=as.factor(long_df4$Stock)

ggplot(data=long_df4,aes(x=kosambi_distances,y=Interference,col=factor(Stock,levels = stock_order)))+geom_line()+facet_wrap(~Treatment)
ggsave("images/COI-vs-RR_byDay.png", height=5)

#Summarize for further plotting:
COI_final=summaryBy(Interference+kosambi_distances~Stock+Treatment+Interval,data=long_df4,FUN=mean,na.rm=T)

#Add CO type to distinguish TCO and DCO values
COI_final$crossover_type=ifelse(COI_final$Interval=="ycv_vf","TCO","DCO")
COI_final$Treatment=as.factor(COI_final$Treatment)
COI_final$treat_order=ifelse(COI_final$Treatment=="Control",1,ifelse(COI_final$Treatment=="Hypoxia",2,ifelse(COI_final$Treatment=="ColdTemp",3,4)))

#Plot per treatment of COI:
gxe_all <- ggplot(data = COI_final,
                       aes(x = kosambi_distances.mean, y = Interference.mean,
                           col = factor(Stock,levels = stock_order),
                           group = interaction(factor(Stock,levels = stock_order), treat_order, crossover_type),
                           linetype = crossover_type,
                           shape = crossover_type)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 5) +
  facet_wrap(~factor(Treatment,levels=c("Control","Hypoxia","ColdTemp","Combined"))) +
  xlab("Genetic Map Length (cM)") + ylab("Crossover Interference") +
  ggtitle("GxE in Rate and Distribution of Recombination") +
  theme_base() +
  scale_color_manual(values = col_scheme2, name = "Stock") +
  scale_linetype_manual(values = c("DCO" = "solid", "TCO" = "dashed"), name = "Crossover\nType") +
  scale_shape_discrete(name = "Crossover\nType")

gxe_all
ggsave("images/COI-vs-RR_byTreatment.png", height=11,width=8.5,units = "in")

#Make recombination Plot to match

#extract total map length
ycv_cvv=subset(long_df4,long_df4$Interval=="ycv_cvv")
ycv_cvv$treat_order=ifelse(ycv_cvv$Treatment=="Control",1,ifelse(ycv_cvv$Treatment=="Hypoxia",2,ifelse(ycv_cvv$Treatment=="ColdTemp",3,4)))

#plot
gxe_map <- ggplot(ycv_cvv, aes(x=treat_order, y=kosambi_distances,col=factor(Stock,levels = stock_order)))  + 
  xlab("Treatment") + ylab("Genetic Map Length (cM)") + 
  theme_base()+ scale_x_continuous(limits=c(0.9,4.1),breaks=c(1,2,3,4),labels=c("Control","Hypoxia","ColdTemp","Combined")) + 
  scale_color_manual(values = col_scheme2)+
  theme(legend.title=element_blank())+ theme(plot.title = element_text(size = 8)) +
  stat_summary(fun = median, geom="line", aes(group=Stock),na.rm=TRUE, linewidth=1) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width=0.1) +
  ggtitle("Recombination Rate of y-f interval") 
gxe_map

ggsave("images/RR-total_ycv_cvv.png", plot=gxe_map, height=5)

ggarrange(gxe_map,COI_figure_ycv_cvv,ncol=2,nrow=1,legend="bottom",common.legend = TRUE,labels=c("A","B"),font.label = list(size = 16))
ggsave("images/ycv-cvv_twoPanel.png",width=8,height=4)


#Get values for Table in Paper:
COI_summary=summaryBy(Interference~Stock+Treatment+Interval,data=long_df4,FUN=mean,na.rm=T)
COI_summary2=summaryBy(Interference~Stock+Treatment+Interval,data=long_df4,FUN=sd,na.rm=T)
COI_summary$Interference.sd=COI_summary2$Interference.sd
write.csv(COI_summary,"output/interference_summary.csv")
