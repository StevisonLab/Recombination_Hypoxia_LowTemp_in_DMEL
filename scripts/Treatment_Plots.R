#inherit from GxE_Plots_Split_Model.R script:
recom_merge$Treatment=factor(recom_merge$Treatment,levels=c("Control","Hypoxia","ColdTemp","Combined"))
recom_merge$treat_order=ifelse(recom_merge$Treatment=="Control",1,ifelse(recom_merge$Treatment=="Hypoxia",2,ifelse(recom_merge$Treatment=="ColdTemp",3,4)))

sig_stocks_A=unique(sort(pheno_contr_full_temp$Stock[pheno_contr_full_temp$full_sig!="" & pheno_contr_full_oxy$vial_letter=="A"]))
sig_stocks_B=unique(sort(pheno_contr_full_temp$Stock[pheno_contr_full_temp$full_sig!="" & pheno_contr_full_oxy$vial_letter=="B"]))
sig_stocks_C=unique(sort(pheno_contr_full_temp$Stock[pheno_contr_full_temp$full_sig!="" & pheno_contr_full_oxy$vial_letter=="C"]))
sig_stocks_D=unique(sort(pheno_contr_full_temp$Stock[pheno_contr_full_temp$full_sig!="" & pheno_contr_full_oxy$vial_letter=="D"]))

recom_merge$highlights = ifelse((recom_merge$Stock %in% sig_stocks_A & recom_merge$vial_letter=="A") | (recom_merge$Stock %in% sig_stocks_B & recom_merge$vial_letter=="B") | (recom_merge$Stock %in% sig_stocks_C & recom_merge$vial_letter=="C") | (recom_merge$Stock %in% sig_stocks_D & recom_merge$vial_letter=="D"),1, 0)

full_figure <- ggplot(recom_merge, aes(x=treat_order, y=recomb_rate*100,col=factor(Stock,levels = stock_order))) +
  facet_wrap(~vial_letter, labeller = labeller(vial_letter =c("A" = "1-2 Days Post-mating","B" = "3-4 Days Post-mating","C" = "5-6 Days Post-mating","D" = "7-8 Days Post-mating"))) +
  xlab("Treatment") + theme(plot.title = element_blank()) +
  ylab(bquote(Genetic~map~length~(cM))) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    strip.text = element_text(size = 8),
    legend.position = "top"
  )+ylim(27.5,88)+
  scale_color_manual( values = col_scheme2)+
  scale_x_continuous(labels=c("Control","8% Oxygen","18°C","18°C + 8% Oxygen")) +
  geom_line(linewidth = 2, lineend = "round")+ theme(legend.title=element_blank()) +
   geom_errorbar(data = recom_merge, 
                 aes(x = treat_order, ymin = 100*asymp.LCL, ymax = 100*asymp.UCL, 
                     group = Stock, color = factor(Stock, levels = stock_order)), 
                 width = 0.1, size = 2) +
   gghighlight(highlights > 0, calculate_per_facet = TRUE,unhighlighted_params = list(colour = NULL, alpha = 0.3))

full_figure

ggsave("images/Line_plot_byDay.png",width = 8,height=6,units="in")


# Repeat but using emmeans values instead:
full_figure_emmeans <- ggplot(recom_merge, aes(x=treat_order, y=prob*100,col=factor(Stock,levels = stock_order))) +
  facet_wrap(~vial_letter, labeller = labeller(vial_letter =c("A" = "1-2 Days Post-mating","B" = "3-4 Days Post-mating","C" = "5-6 Days Post-mating","D" = "7-8 Days Post-mating"))) +
  xlab("Treatment") + theme(plot.title = element_blank()) +
  ylab(bquote(Genetic~map~length~(cM))) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    strip.text = element_text(size = 8),
    legend.position = "top"
  )+ylim(27.5,88)+
  scale_color_manual( values = col_scheme2)+
  scale_x_continuous(labels=c("Control","8% Oxygen","18°C","18°C + 8% Oxygen")) +
  geom_line(linewidth = 2, lineend = "round")+ theme(legend.title=element_blank()) +
  geom_errorbar(data = recom_merge, 
                aes(x = treat_order, ymin = 100*asymp.LCL, ymax = 100*asymp.UCL, 
                    group = Stock, color = factor(Stock, levels = stock_order)), 
                width = 0.1, size = 2) +
  gghighlight(highlights > 0, calculate_per_facet = TRUE,unhighlighted_params = list(colour = NULL, alpha = 0.3))

full_figure_emmeans

ggsave("images/Line_plot_byDay_emmeans.png",width = 8,height=6,units="in")

