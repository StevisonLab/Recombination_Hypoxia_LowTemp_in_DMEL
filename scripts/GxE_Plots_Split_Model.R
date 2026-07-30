
stock_order=c("508","42","57","391","491")
col_scheme2=c("#1b9e77","#d95f02","#7570b3","#66a61e","#e7298a")


# Oxygen Reaction Norm Plots -----------------------------------------

#Get recomb_rates from split model:
fit_contrast_means <- emmeans(fit_full_split, ~ Stock+temperature+oxygen, type="response")
fit_contr_yf <- contrast(fit_contrast_means)
#fit_contr_yf
pheno_contr_yf <- as.data.frame(summary(fit_contr_yf))
#pheno_contr_yf

cld_results <- cld(fit_contrast_means, adjust = "sidak", Letters = letters)
cld_results$temp_mod=ifelse(cld_results$temperature=="Low",18,25)
#cld_results$temperature=ifelse(cld_results$temperature=="Low","18°C","25°C")
cld_results$oxy_mod=ifelse(cld_results$oxygen=="Low",8,21)

#Add significance (manually based on posthoc contrast table):
#Note: here if ANY of the timepoints are significant, the strain is indicated as signficant here. The ONLY signficant reaction norm here is the temperature slop of 57 at low oxygen!
low_stocks=unique(sort(pheno_contr_full_oxy$Stock[pheno_contr_full_oxy$full_sig!="" & pheno_contr_full_oxy$temperature=="Low"]))
con_stocks=unique(sort(pheno_contr_full_oxy$Stock[pheno_contr_full_oxy$full_sig!="" & pheno_contr_full_oxy$temperature=="Control"]))
cld_results$oxy_highlights = ifelse((cld_results$temp_mod == "18" & cld_results$Stock %in% low_stocks) | 
                                       (cld_results$temp_mod == "25" & cld_results$Stock %in% con_stocks),1, 0)

recom_merge=merge(cld_results,recomb_summary,by.x=c("Stock","oxygen","temperature"),by.y=c("Stock","oxygen","temperature"))
recom_merge$recomb_rate=recom_merge$recomb_rate/100

full_figure_GxE <- ggplot(recom_merge, aes(x=oxy_mod, y=100*prob,col=factor(Stock,levels = stock_order))) +
  facet_wrap(~temp_mod, labeller = labeller(temp_mod =c("18" = "18°C","25" = "25°C"))) +
  xlab("Oxygen") + theme(plot.title = element_blank()) +
  ylab(bquote(Map~Length~(cM))) +
  theme_base()+ylim(50,70)+
  scale_color_manual( values = col_scheme2)+
  scale_x_continuous(limits=c(8,21),breaks=c(8,21),labels=c("8%","21%")) +
  geom_line(linewidth = 2, lineend = "round")+ theme(legend.title=element_blank()) +
  geom_errorbar(data = recom_merge, 
                   aes(x = oxy_mod, ymin = 100*asymp.LCL, ymax = 100*asymp.UCL, 
                       group = Stock, color = factor(Stock, levels = stock_order)), 
                   width = 0.1, size = 2) +
  gghighlight(oxy_highlights > 0, calculate_per_facet = TRUE,unhighlighted_params = list(colour = NULL, alpha = 0.3))
full_figure_GxE

ggsave("images/GxE_oxy_full.png",width = 8,height=6,units="in")


# Temperature Reaction Norm Overall ---------------------------------------

# Repeat to have a similar plot by Temperatures instead of Oxygen:

#Add significance (based on posthoc contrast table):
#Note: here if ANY of the timepoints are significant, the strain is indicated as signficant here. The ONLY signficant reaction norm here is the temperature slop of 57 at low oxygen!
low_stocks=unique(sort(pheno_contr_full_temp$Stock[pheno_contr_full_temp$full_sig!="" & pheno_contr_full_temp$oxygen=="Low"]))
con_stocks=unique(sort(pheno_contr_full_temp$Stock[pheno_contr_full_temp$full_sig!="" & pheno_contr_full_temp$oxygen=="Control"]))
cld_results$temp_highlights = ifelse((cld_results$oxygen == "Low" & cld_results$Stock %in% low_stocks) | 
                                      (cld_results$oxygen == "Control" & cld_results$Stock %in% con_stocks),1, 0)

recom_merge=merge(cld_results,recomb_summary,by.x=c("Stock","oxygen","temperature"),by.y=c("Stock","oxygen","temperature"))
recom_merge$recomb_rate=recom_merge$recomb_rate/100

full_figure_GxE_temp <- ggplot(recom_merge, aes(x=temp_mod, y=prob*100,col=factor(Stock,levels = stock_order))) +
  facet_wrap(~oxy_mod, labeller = labeller(oxy_mod =c("8" = "8%","21" = "21%"))) +
  xlab("Temperature") + theme(plot.title = element_blank()) +
  ylab(bquote(Map~Length~(cM))) +
  theme_base()+ylim(50,70)+
  scale_color_manual( values = col_scheme2)+
  scale_x_continuous(limits=c(18,25),breaks=c(18,25),labels=c("18°C","25°C")) +
  geom_line(linewidth = 2, lineend = "round")+ theme(legend.title=element_blank()) +
  geom_errorbar(data = recom_merge, 
                aes(x = temp_mod, ymin = 100*asymp.LCL, ymax = 100*asymp.UCL, 
                    group = Stock, color = factor(Stock, levels = stock_order)), 
                width = 0.1, size = 2) +
  gghighlight(temp_highlights > 0, calculate_per_facet = TRUE,unhighlighted_params = list(colour = NULL, alpha = 0.3))
full_figure_GxE_temp

ggsave("images/GxE_temp_full.png",width = 8,height=6,units="in")



# Per Day Oxygen Reaction Norm Plots --------------------------------------

# Repeat to get plots PER timepoint:
#Get recomb_rates from split model:
fit_contrast_means <- emmeans(fit_full_split, ~ Stock+temperature+oxygen+vial_letter, type="response")

cld_results2 <- cld(fit_contrast_means, adjust = "sidak", Letters = letters)
cld_results2$temp_mod=ifelse(cld_results2$temperature=="Low",18,25)
#cld_results2$temperature=ifelse(cld_results2$temperature=="Low","18°C","25°C")
cld_results2$oxy_mod=ifelse(cld_results2$oxygen=="Low",8,21)

#Regenerate recomb summary broken down by day!
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
recomb_summary <- aggregate(recomb_rate[, c("recomb_rate", "ycv_rr", "cvv_rr", "vf_rr")], by = list(recomb_rate$Treatment, recomb_rate$Stock,recomb_rate$vial_letter), mean, na.rm=TRUE)
names(recomb_summary)[1:2] <- c("Treatment", "Stock")

recomb_summary$Stock=as.factor(recomb_summary$Stock)
recomb_summary$oxygen=ifelse(recomb_summary$Treatment=="Combined","Low",ifelse(recomb_summary$Treatment=="Hypoxia","Low","Control"))
recomb_summary$temperature=ifelse(recomb_summary$Treatment=="Control","Control",ifelse(recomb_summary$Treatment=="Hypoxia","Control","Low"))

recom_merge=merge(cld_results2,recomb_summary,by.x=c("Stock","oxygen","temperature","vial_letter"),by.y=c("Stock","oxygen","temperature","Group.3"))
recom_merge$recomb_rate=recom_merge$recomb_rate/100

#Add significance (manually based on posthoc contrast table):
#Note: here if ANY of the timepoints are significant, the strain is indicated as signficant here. The ONLY signficant reaction norm here is the temperature slop of 57 at low oxygen!

A_plot=subset(recom_merge,recom_merge$vial_letter=="A")

low_stocks=unique(sort(pheno_contr_full_oxy$Stock[pheno_contr_full_oxy$full_sig!="" & pheno_contr_full_oxy$temperature=="Low" & pheno_contr_full_oxy$vial_letter=="A"]))

con_stocks=unique(sort(pheno_contr_full_oxy$Stock[pheno_contr_full_oxy$full_sig!="" & pheno_contr_full_oxy$temperature=="Control" & pheno_contr_full_oxy$vial_letter=="A"]))

A_plot$oxy_highlights = ifelse((A_plot$temp_mod == "18" & A_plot$Stock %in% low_stocks) | 
                                 (A_plot$temp_mod == "25" & A_plot$Stock %in% con_stocks),1, 0)

full_figure_GxE_A <- ggplot(A_plot, aes(x=oxy_mod, y=prob*100,col=factor(Stock,levels = stock_order))) +
  facet_wrap(~temp_mod, labeller = labeller(temp_mod =c("18" = "18°C","25" = "25°C"))) +
  xlab("Oxygen") + theme(plot.title = element_blank()) +
  ylab(bquote(Map~Length~(cM))) +
  theme_base()+ylim(27.5,88)+
  scale_color_manual( values = col_scheme2)+
  scale_x_continuous(limits=c(8,21),breaks=c(8,21),labels=c("8%","21%")) +
  geom_line(linewidth = 2, lineend = "round")+ theme(legend.title=element_blank()) +
  geom_errorbar(data = A_plot, 
                aes(x = oxy_mod, ymin = 100*asymp.LCL, ymax = 100*asymp.UCL, 
                    group = Stock, color = factor(Stock, levels = stock_order)), 
                width = 0.1, size = 2) +
  gghighlight(oxy_highlights > 0, calculate_per_facet = TRUE,unhighlighted_params = list(colour = NULL, alpha = 0.3))
full_figure_GxE_A

ggsave("images/GxE_oxy_full_A.png",width = 8,height=6,units="in")

#Repeat for B:
B_plot=subset(recom_merge,recom_merge$vial_letter=="B")

low_stocks=unique(sort(pheno_contr_full_oxy$Stock[pheno_contr_full_oxy$full_sig!="" & pheno_contr_full_oxy$temperature=="Low" & pheno_contr_full_oxy$vial_letter=="B"]))

con_stocks=unique(sort(pheno_contr_full_oxy$Stock[pheno_contr_full_oxy$full_sig!="" & pheno_contr_full_oxy$temperature=="Control" & pheno_contr_full_oxy$vial_letter=="B"]))

B_plot$oxy_highlights = ifelse((B_plot$temp_mod == "18" & B_plot$Stock %in% low_stocks) | 
                                 (B_plot$temp_mod == "25" & B_plot$Stock %in% con_stocks),1, 0)

full_figure_GxE_B <- ggplot(B_plot, aes(x=oxy_mod, y=prob*100,col=factor(Stock,levels = stock_order))) +
  facet_wrap(~temp_mod, labeller = labeller(temp_mod =c("18" = "18°C","25" = "25°C"))) +
  xlab("Oxygen") + theme(plot.title = element_blank()) +
  ylab(bquote(Map~Length~(cM))) +
  theme_base()+ylim(27.5,88)+
  scale_color_manual( values = col_scheme2)+
  scale_x_continuous(limits=c(8,21),breaks=c(8,21),labels=c("8%","21%")) +
  geom_line(linewidth = 2, lineend = "round")+ theme(legend.title=element_blank()) +
  geom_errorbar(data = B_plot, 
                aes(x = oxy_mod, ymin = 100*asymp.LCL, ymax = 100*asymp.UCL, 
                    group = Stock, color = factor(Stock, levels = stock_order)), 
                width = 0.1, size = 2) +
  gghighlight(oxy_highlights > 0, calculate_per_facet = TRUE,unhighlighted_params = list(colour = NULL, alpha = 0.3))
full_figure_GxE_B

ggsave("images/GxE_oxy_full_B.png",width = 8,height=6,units="in")

#Repeat for C:
C_plot=subset(recom_merge,recom_merge$vial_letter=="C")

low_stocks=unique(sort(pheno_contr_full_oxy$Stock[pheno_contr_full_oxy$full_sig!="" & pheno_contr_full_oxy$temperature=="Low" & pheno_contr_full_oxy$vial_letter=="C"]))

con_stocks=unique(sort(pheno_contr_full_oxy$Stock[pheno_contr_full_oxy$full_sig!="" & pheno_contr_full_oxy$temperature=="Control" & pheno_contr_full_oxy$vial_letter=="C"]))

C_plot$oxy_highlights = ifelse((C_plot$temp_mod == "18" & C_plot$Stock %in% low_stocks) | 
                                 (C_plot$temp_mod == "25" & C_plot$Stock %in% con_stocks),1, 0)

full_figure_GxE_C <- ggplot(C_plot, aes(x=oxy_mod, y=prob*100,col=factor(Stock,levels = stock_order))) +
  facet_wrap(~temp_mod, labeller = labeller(temp_mod =c("18" = "18°C","25" = "25°C"))) +
  xlab("Oxygen") + theme(plot.title = element_blank()) +
  ylab(bquote(Map~Length~(cM))) +
  theme_base()+ylim(27.5,88)+
  scale_color_manual( values = col_scheme2)+
  scale_x_continuous(limits=c(8,21),breaks=c(8,21),labels=c("8%","21%")) +
  geom_line(linewidth = 2, lineend = "round")+ theme(legend.title=element_blank()) +
  geom_errorbar(data = C_plot, 
                aes(x = oxy_mod, ymin = 100*asymp.LCL, ymax = 100*asymp.UCL, 
                    group = Stock, color = factor(Stock, levels = stock_order)), 
                width = 0.1, size = 2) +
  gghighlight(oxy_highlights > 0, calculate_per_facet = TRUE,unhighlighted_params = list(colour = NULL, alpha = 0.3))
full_figure_GxE_C

ggsave("images/GxE_oxy_full_C.png",width = 8,height=6,units="in")

#Repeat for D:
D_plot=subset(recom_merge,recom_merge$vial_letter=="D")

low_stocks=unique(sort(pheno_contr_full_oxy$Stock[pheno_contr_full_oxy$full_sig!="" & pheno_contr_full_oxy$temperature=="Low" & pheno_contr_full_oxy$vial_letter=="D"]))

con_stocks=unique(sort(pheno_contr_full_oxy$Stock[pheno_contr_full_oxy$full_sig!="" & pheno_contr_full_oxy$temperature=="Control" & pheno_contr_full_oxy$vial_letter=="D"]))

D_plot$oxy_highlights = ifelse((D_plot$temp_mod == "18" & D_plot$Stock %in% low_stocks) | 
                                 (D_plot$temp_mod == "25" & D_plot$Stock %in% con_stocks),1, 0)

full_figure_GxE_D <- ggplot(D_plot, aes(x=oxy_mod, y=prob*100,col=factor(Stock,levels = stock_order))) +
  facet_wrap(~temp_mod, labeller = labeller(temp_mod =c("18" = "18°C","25" = "25°C"))) +
  xlab("Oxygen") + theme(plot.title = element_blank()) +
  ylab(bquote(Map~Length~(cM))) +
  theme_base()+ylim(27.5,88)+
  scale_color_manual( values = col_scheme2)+
  scale_x_continuous(limits=c(8,21),breaks=c(8,21),labels=c("8%","21%")) +
  geom_line(linewidth = 2, lineend = "round")+ theme(legend.title=element_blank()) +
  geom_errorbar(data = D_plot, 
                aes(x = oxy_mod, ymin = 100*asymp.LCL, ymax = 100*asymp.UCL, 
                    group = Stock, color = factor(Stock, levels = stock_order)), 
                width = 0.1, size = 2) +
  gghighlight(oxy_highlights > 0, calculate_per_facet = TRUE,unhighlighted_params = list(colour = NULL, alpha = 0.3))
full_figure_GxE_D

ggsave("images/GxE_oxy_full_D.png",width = 8,height=6,units="in")


# Temperature Reaction Per Day Norm Plots -----------------------------------------


# Repeat to have a similar plot by Temperatures instead of Oxygen:

#Add significance (based on posthoc contrast table):
#Note: here if ANY of the timepoints are significant, the strain is indicated as signficant here. The ONLY signficant reaction norm here is the temperature slop of 57 at low oxygen!
low_stocks=unique(sort(pheno_contr_full_temp$Stock[pheno_contr_full_temp$full_sig!="" & pheno_contr_full_temp$oxygen=="Low" & pheno_contr_full_oxy$vial_letter=="A"]))
con_stocks=unique(sort(pheno_contr_full_temp$Stock[pheno_contr_full_temp$full_sig!="" & pheno_contr_full_temp$oxygen=="Control" & pheno_contr_full_oxy$vial_letter=="A"]))

A_plot$temp_highlights = ifelse((A_plot$oxygen == "Low" & A_plot$Stock %in% low_stocks) | 
                                 (A_plot$oxygen == "Control" & A_plot$Stock %in% con_stocks),1, 0)

full_figure_GxE_temp_A <- ggplot(A_plot, aes(x=temp_mod, y=prob*100,col=factor(Stock,levels = stock_order))) +
  facet_wrap(~oxy_mod, labeller = labeller(oxy_mod =c("8" = "8%","21" = "21%"))) +
  xlab("Temperature") + theme(plot.title = element_blank()) +
  ylab(bquote(Map~Length~(cM))) +
  theme_base()+ylim(27.5,88)+
  scale_color_manual( values = col_scheme2)+
  scale_x_continuous(limits=c(18,25),breaks=c(18,25),labels=c("18°C","25°C")) +
  geom_line(linewidth = 2, lineend = "round")+ theme(legend.title=element_blank()) +
  geom_errorbar(data = A_plot, 
                aes(x = temp_mod, ymin = 100*asymp.LCL, ymax = 100*asymp.UCL, 
                    group = Stock, color = factor(Stock, levels = stock_order)), 
                width = 0.1, size = 2) +
  gghighlight(temp_highlights > 0, calculate_per_facet = TRUE,unhighlighted_params = list(colour = NULL, alpha = 0.3))
full_figure_GxE_temp_A

ggsave("images/GxE_temp_full_A.png",width = 8,height=6,units="in")

#Repeat for B:
low_stocks=unique(sort(pheno_contr_full_temp$Stock[pheno_contr_full_temp$full_sig!="" & pheno_contr_full_temp$oxygen=="Low" & pheno_contr_full_oxy$vial_letter=="B"]))
con_stocks=unique(sort(pheno_contr_full_temp$Stock[pheno_contr_full_temp$full_sig!="" & pheno_contr_full_temp$oxygen=="Control" & pheno_contr_full_oxy$vial_letter=="B"]))

B_plot$temp_highlights = ifelse((B_plot$oxygen == "Low" & B_plot$Stock %in% low_stocks) | 
                                  (B_plot$oxygen == "Control" & B_plot$Stock %in% con_stocks),1, 0)

full_figure_GxE_temp_B <- ggplot(B_plot, aes(x=temp_mod, y=prob*100,col=factor(Stock,levels = stock_order))) +
  facet_wrap(~oxy_mod, labeller = labeller(oxy_mod =c("8" = "8%","21" = "21%"))) +
  xlab("Temperature") + theme(plot.title = element_blank()) +
  ylab(bquote(Map~Length~(cM))) +
  theme_base()+ylim(27.5,88)+
  scale_color_manual( values = col_scheme2)+
  scale_x_continuous(limits=c(18,25),breaks=c(18,25),labels=c("18°C","25°C")) +
  geom_line(linewidth = 2, lineend = "round")+ theme(legend.title=element_blank()) +
  geom_errorbar(data = B_plot, 
                aes(x = temp_mod, ymin = 100*asymp.LCL, ymax = 100*asymp.UCL, 
                    group = Stock, color = factor(Stock, levels = stock_order)), 
                width = 0.1, size = 2) +
  gghighlight(temp_highlights > 0, calculate_per_facet = TRUE,unhighlighted_params = list(colour = NULL, alpha = 0.3))
full_figure_GxE_temp_B

ggsave("images/GxE_temp_full_B.png",width = 8,height=6,units="in")


# Repeat for C:
low_stocks=unique(sort(pheno_contr_full_temp$Stock[pheno_contr_full_temp$full_sig!="" & pheno_contr_full_temp$oxygen=="Low" & pheno_contr_full_oxy$vial_letter=="C"]))
con_stocks=unique(sort(pheno_contr_full_temp$Stock[pheno_contr_full_temp$full_sig!="" & pheno_contr_full_temp$oxygen=="Control" & pheno_contr_full_oxy$vial_letter=="C"]))

C_plot$temp_highlights = ifelse((C_plot$oxygen == "Low" & C_plot$Stock %in% low_stocks) | 
                                  (C_plot$oxygen == "Control" & C_plot$Stock %in% con_stocks),1, 0)

full_figure_GxE_temp_C <- ggplot(C_plot, aes(x=temp_mod, y=prob*100,col=factor(Stock,levels = stock_order))) +
  facet_wrap(~oxy_mod, labeller = labeller(oxy_mod =c("8" = "8%","21" = "21%"))) +
  xlab("Temperature") + theme(plot.title = element_blank()) +
  ylab(bquote(Map~Length~(cM))) +
  theme_base()+ylim(27.5,88)+
  scale_color_manual( values = col_scheme2)+
  scale_x_continuous(limits=c(18,25),breaks=c(18,25),labels=c("18°C","25°C")) +
  geom_line(linewidth = 2, lineend = "round")+ theme(legend.title=element_blank()) +
  geom_errorbar(data = C_plot, 
                aes(x = temp_mod, ymin = 100*asymp.LCL, ymax = 100*asymp.UCL, 
                    group = Stock, color = factor(Stock, levels = stock_order)), 
                width = 0.1, size = 2) +
  gghighlight(temp_highlights > 0, calculate_per_facet = TRUE,unhighlighted_params = list(colour = NULL, alpha = 0.3))
full_figure_GxE_temp_C

ggsave("images/GxE_temp_full_C.png",width = 8,height=6,units="in")

# Repeat for D:
low_stocks=unique(sort(pheno_contr_full_temp$Stock[pheno_contr_full_temp$full_sig!="" & pheno_contr_full_temp$oxygen=="Low" & pheno_contr_full_oxy$vial_letter=="D"]))
con_stocks=unique(sort(pheno_contr_full_temp$Stock[pheno_contr_full_temp$full_sig!="" & pheno_contr_full_temp$oxygen=="Control" & pheno_contr_full_oxy$vial_letter=="D"]))

D_plot$temp_highlights = ifelse((D_plot$oxygen == "Low" & D_plot$Stock %in% low_stocks) | 
                                  (D_plot$oxygen == "Control" & D_plot$Stock %in% con_stocks),1, 0)

full_figure_GxE_temp_D <- ggplot(D_plot, aes(x=temp_mod, y=prob*100,col=factor(Stock,levels = stock_order))) +
  facet_wrap(~oxy_mod, labeller = labeller(oxy_mod =c("8" = "8%","21" = "21%"))) +
  xlab("Temperature") + theme(plot.title = element_blank()) +
  ylab(bquote(Map~Length~(cM))) +
  theme_base()+ylim(27.5,88)+
  scale_color_manual( values = col_scheme2)+
  scale_x_continuous(limits=c(18,25),breaks=c(18,25),labels=c("18°C","25°C")) +
  geom_line(linewidth = 2, lineend = "round")+ theme(legend.title=element_blank()) +
  geom_errorbar(data = D_plot, 
                aes(x = temp_mod, ymin = 100*asymp.LCL, ymax = 100*asymp.UCL, 
                    group = Stock, color = factor(Stock, levels = stock_order)), 
                width = 0.1, size = 2) +
  gghighlight(temp_highlights > 0, calculate_per_facet = TRUE,unhighlighted_params = list(colour = NULL, alpha = 0.3))
full_figure_GxE_temp_D

ggsave("images/GxE_temp_full_D.png",width = 8,height=6,units="in")

