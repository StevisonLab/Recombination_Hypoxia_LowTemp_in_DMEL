# Odds Ratio Plots --------------------------------------------------

col_scheme=c("#e41a1c",  "#377eb8",  "#4daf4a",  "#984ea3",  "#ff7f00", "#ffff33")

odds_ratio=subset(pheno_contr_full,pheno_contr_full$Stock=="42")
odds_ratio$full_or=exp(odds_ratio$estimate)
odds_ratio$full_sig=ifelse(odds_ratio$p.value<0.001,"***",ifelse(odds_ratio$p.value<0.01,"**",ifelse(odds_ratio$p.value<0.05,"*","")))
#odds_ratio

pdf("images/full_odds_42.pdf")

odds_figure_full=ggplot(aes(y=full_or,x=vial_letter,group=contrast,col=contrast,shape=contrast),data=odds_ratio)+ggtitle("Full Model Results DGRP_42")+scale_colour_manual(values=col_scheme)+ 
  geom_point(size=3)+ylab("Odds Ratios")+theme_base()+geom_hline(yintercept = 1,linetype="dashed",color="grey")+ylim(0.5,2.75)+theme(axis.text.x = element_text(angle = 45))+
  scale_x_discrete(name="Days post-mating",labels=c("1-2","3-4","5-6","7-8"))+geom_line()+geom_errorbar(aes(ymin=full_or-SE,ymax=full_or+SE))+
  annotate(geom="text", x=1, y=2, label=odds_ratio$full_sig[1],color="#e41a1c",size=10)+annotate(geom="text", x=2, y=2, label=odds_ratio$full_sig[7],color="#e41a1c",size=10)+annotate(geom="text", x=3, y=2, label=odds_ratio$full_sig[13],color="#e41a1c",size=10)+annotate(geom="text", x=4, y=2, label=odds_ratio$full_sig[19],color="#e41a1c",size=10)+
  annotate(geom="text", x=1, y=2.1, label=odds_ratio$full_sig[2],color="#377eb8",size=10)+annotate(geom="text", x=2, y=2.1, label=odds_ratio$full_sig[8],color="#377eb8",size=10)+annotate(geom="text", x=3, y=2.1, label=odds_ratio$full_sig[14],color="#377eb8",size=10)+annotate(geom="text", x=4, y=2.1, label=odds_ratio$full_sig[20],color="#377eb8",size=10)+
  annotate(geom="text", x=1, y=2.2, label=odds_ratio$full_sig[3],color="#4daf4a",size=10)+annotate(geom="text", x=2, y=2.2, label=odds_ratio$full_sig[9],color="#4daf4a",size=10)+annotate(geom="text", x=3, y=2.2, label=odds_ratio$full_sig[15],color="#4daf4a",size=10)+annotate(geom="text", x=4, y=2.2, label=odds_ratio$full_sig[21],color="#4daf4a",size=10)+
  annotate(geom="text", x=1, y=2.3, label=odds_ratio$full_sig[4],color="#984ea3",size=10)+annotate(geom="text", x=2, y=2.3, label=odds_ratio$full_sig[10],color="#984ea3",size=10)+annotate(geom="text", x=3, y=2.3, label=odds_ratio$full_sig[16],color="#984ea3",size=10)+annotate(geom="text", x=4, y=2.3, label=odds_ratio$full_sig[22],color="#984ea3",size=10)+
  annotate(geom="text", x=1, y=2.4, label=odds_ratio$full_sig[5],color="#ff7f00",size=10)+annotate(geom="text", x=2, y=2.4, label=odds_ratio$full_sig[11],color="#ff7f00",size=10)+annotate(geom="text", x=3, y=2.4, label=odds_ratio$full_sig[17],color="#ff7f00",size=10)+annotate(geom="text", x=4, y=2.4, label=odds_ratio$full_sig[23],color="#ff7f00",size=10)+
  annotate(geom="text", x=1, y=2.5, label=odds_ratio$full_sig[6],color="#ffff33",size=10)+annotate(geom="text", x=2, y=2.5, label=odds_ratio$full_sig[12],color="#ffff33",size=10)+annotate(geom="text", x=3, y=2.5, label=odds_ratio$full_sig[18],color="#ffff33",size=10)+annotate(geom="text", x=4, y=2.5, label=odds_ratio$full_sig[24],color="#ffff33",size=10)
odds_figure_full
dev.off()

#repeat for 57
odds_ratio=subset(pheno_contr_full,pheno_contr_full$Stock=="57")
odds_ratio$full_or=exp(odds_ratio$estimate)
odds_ratio$full_sig=ifelse(odds_ratio$p.value<0.001,"***",ifelse(odds_ratio$p.value<0.01,"**",ifelse(odds_ratio$p.value<0.05,"*","")))
#odds_ratio

pdf("images/full_odds_57.pdf")
odds_figure_full=ggplot(aes(y=full_or,x=vial_letter,group=contrast,col=contrast,shape=contrast),data=odds_ratio)+ggtitle("Full Model Results DGRP_57")+scale_colour_manual(values=col_scheme)+ 
  geom_point(size=3)+ylab("Odds Ratios")+theme_base()+geom_hline(yintercept = 1,linetype="dashed",color="grey")+ylim(0.5,2.75)+theme(axis.text.x = element_text(angle = 45))+
  scale_x_discrete(name="Days post-mating",labels=c("1-2","3-4","5-6","7-8"))+geom_line()+geom_errorbar(aes(ymin=full_or-SE,ymax=full_or+SE))+
  annotate(geom="text", x=1, y=2, label=odds_ratio$full_sig[1],color="#e41a1c",size=10)+annotate(geom="text", x=2, y=2, label=odds_ratio$full_sig[7],color="#e41a1c",size=10)+annotate(geom="text", x=3, y=2, label=odds_ratio$full_sig[13],color="#e41a1c",size=10)+annotate(geom="text", x=4, y=2, label=odds_ratio$full_sig[19],color="#e41a1c",size=10)+
  annotate(geom="text", x=1, y=2.1, label=odds_ratio$full_sig[2],color="#377eb8",size=10)+annotate(geom="text", x=2, y=2.1, label=odds_ratio$full_sig[8],color="#377eb8",size=10)+annotate(geom="text", x=3, y=2.1, label=odds_ratio$full_sig[14],color="#377eb8",size=10)+annotate(geom="text", x=4, y=2.1, label=odds_ratio$full_sig[20],color="#377eb8",size=10)+
  annotate(geom="text", x=1, y=2.2, label=odds_ratio$full_sig[3],color="#4daf4a",size=10)+annotate(geom="text", x=2, y=2.2, label=odds_ratio$full_sig[9],color="#4daf4a",size=10)+annotate(geom="text", x=3, y=2.2, label=odds_ratio$full_sig[15],color="#4daf4a",size=10)+annotate(geom="text", x=4, y=2.2, label=odds_ratio$full_sig[21],color="#4daf4a",size=10)+
  annotate(geom="text", x=1, y=2.3, label=odds_ratio$full_sig[4],color="#984ea3",size=10)+annotate(geom="text", x=2, y=2.3, label=odds_ratio$full_sig[10],color="#984ea3",size=10)+annotate(geom="text", x=3, y=2.3, label=odds_ratio$full_sig[16],color="#984ea3",size=10)+annotate(geom="text", x=4, y=2.3, label=odds_ratio$full_sig[22],color="#984ea3",size=10)+
  annotate(geom="text", x=1, y=2.4, label=odds_ratio$full_sig[5],color="#ff7f00",size=10)+annotate(geom="text", x=2, y=2.4, label=odds_ratio$full_sig[11],color="#ff7f00",size=10)+annotate(geom="text", x=3, y=2.4, label=odds_ratio$full_sig[17],color="#ff7f00",size=10)+annotate(geom="text", x=4, y=2.4, label=odds_ratio$full_sig[23],color="#ff7f00",size=10)+
  annotate(geom="text", x=1, y=2.5, label=odds_ratio$full_sig[6],color="#ffff33",size=10)+annotate(geom="text", x=2, y=2.5, label=odds_ratio$full_sig[12],color="#ffff33",size=10)+annotate(geom="text", x=3, y=2.5, label=odds_ratio$full_sig[18],color="#ffff33",size=10)+annotate(geom="text", x=4, y=2.5, label=odds_ratio$full_sig[24],color="#ffff33",size=10)
odds_figure_full
dev.off()

#repeat for 491
odds_ratio=subset(pheno_contr_full,pheno_contr_full$Stock=="491")
odds_ratio$full_or=exp(odds_ratio$estimate)
odds_ratio$full_sig=ifelse(odds_ratio$p.value<0.001,"***",ifelse(odds_ratio$p.value<0.01,"**",ifelse(odds_ratio$p.value<0.05,"*","")))
#odds_ratio

pdf("images/full_odds_491.pdf")
odds_figure_full=ggplot(aes(y=full_or,x=vial_letter,group=contrast,col=contrast,shape=contrast),data=odds_ratio)+ggtitle("Full Model Results DGRP_491")+scale_colour_manual(values=col_scheme)+ 
  geom_point(size=3)+ylab("Odds Ratios")+theme_base()+geom_hline(yintercept = 1,linetype="dashed",color="grey")+ylim(0.5,2.75)+theme(axis.text.x = element_text(angle = 45))+
  scale_x_discrete(name="Days post-mating",labels=c("1-2","3-4","5-6","7-8"))+geom_line()+geom_errorbar(aes(ymin=full_or-SE,ymax=full_or+SE))+
  annotate(geom="text", x=1, y=2, label=odds_ratio$full_sig[1],color="#e41a1c",size=10)+annotate(geom="text", x=2, y=2, label=odds_ratio$full_sig[7],color="#e41a1c",size=10)+annotate(geom="text", x=3, y=2, label=odds_ratio$full_sig[13],color="#e41a1c",size=10)+annotate(geom="text", x=4, y=2, label=odds_ratio$full_sig[19],color="#e41a1c",size=10)+
  annotate(geom="text", x=1, y=2.1, label=odds_ratio$full_sig[2],color="#377eb8",size=10)+annotate(geom="text", x=2, y=2.1, label=odds_ratio$full_sig[8],color="#377eb8",size=10)+annotate(geom="text", x=3, y=2.1, label=odds_ratio$full_sig[14],color="#377eb8",size=10)+annotate(geom="text", x=4, y=2.1, label=odds_ratio$full_sig[20],color="#377eb8",size=10)+
  annotate(geom="text", x=1, y=2.2, label=odds_ratio$full_sig[3],color="#4daf4a",size=10)+annotate(geom="text", x=2, y=2.2, label=odds_ratio$full_sig[9],color="#4daf4a",size=10)+annotate(geom="text", x=3, y=2.2, label=odds_ratio$full_sig[15],color="#4daf4a",size=10)+annotate(geom="text", x=4, y=2.2, label=odds_ratio$full_sig[21],color="#4daf4a",size=10)+
  annotate(geom="text", x=1, y=2.3, label=odds_ratio$full_sig[4],color="#984ea3",size=10)+annotate(geom="text", x=2, y=2.3, label=odds_ratio$full_sig[10],color="#984ea3",size=10)+annotate(geom="text", x=3, y=2.3, label=odds_ratio$full_sig[16],color="#984ea3",size=10)+annotate(geom="text", x=4, y=2.3, label=odds_ratio$full_sig[22],color="#984ea3",size=10)+
  annotate(geom="text", x=1, y=2.4, label=odds_ratio$full_sig[5],color="#ff7f00",size=10)+annotate(geom="text", x=2, y=2.4, label=odds_ratio$full_sig[11],color="#ff7f00",size=10)+annotate(geom="text", x=3, y=2.4, label=odds_ratio$full_sig[17],color="#ff7f00",size=10)+annotate(geom="text", x=4, y=2.4, label=odds_ratio$full_sig[23],color="#ff7f00",size=10)+
  annotate(geom="text", x=1, y=2.5, label=odds_ratio$full_sig[6],color="#ffff33",size=10)+annotate(geom="text", x=2, y=2.5, label=odds_ratio$full_sig[12],color="#ffff33",size=10)+annotate(geom="text", x=3, y=2.5, label=odds_ratio$full_sig[18],color="#ffff33",size=10)+annotate(geom="text", x=4, y=2.5, label=odds_ratio$full_sig[24],color="#ffff33",size=10)
odds_figure_full
dev.off()

#repeat for 391
odds_ratio=subset(pheno_contr_full,pheno_contr_full$Stock=="391")
odds_ratio$full_or=exp(odds_ratio$estimate)
odds_ratio$full_sig=ifelse(odds_ratio$p.value<0.001,"***",ifelse(odds_ratio$p.value<0.01,"**",ifelse(odds_ratio$p.value<0.05,"*","")))
#odds_ratio

pdf("images/full_odds_391.pdf")
odds_figure_full=ggplot(aes(y=full_or,x=vial_letter,group=contrast,col=contrast,shape=contrast),data=odds_ratio)+ggtitle("Full Model Results DGRP_391")+scale_colour_manual(values=col_scheme)+ 
  geom_point(size=3)+ylab("Odds Ratios")+theme_base()+geom_hline(yintercept = 1,linetype="dashed",color="grey")+ylim(0.5,2.75)+theme(axis.text.x = element_text(angle = 45))+
  scale_x_discrete(name="Days post-mating",labels=c("1-2","3-4","5-6","7-8"))+geom_line()+geom_errorbar(aes(ymin=full_or-SE,ymax=full_or+SE))+
  annotate(geom="text", x=1, y=2, label=odds_ratio$full_sig[1],color="#e41a1c",size=10)+annotate(geom="text", x=2, y=2, label=odds_ratio$full_sig[7],color="#e41a1c",size=10)+annotate(geom="text", x=3, y=2, label=odds_ratio$full_sig[13],color="#e41a1c",size=10)+annotate(geom="text", x=4, y=2, label=odds_ratio$full_sig[19],color="#e41a1c",size=10)+
  annotate(geom="text", x=1, y=2.1, label=odds_ratio$full_sig[2],color="#377eb8",size=10)+annotate(geom="text", x=2, y=2.1, label=odds_ratio$full_sig[8],color="#377eb8",size=10)+annotate(geom="text", x=3, y=2.1, label=odds_ratio$full_sig[14],color="#377eb8",size=10)+annotate(geom="text", x=4, y=2.1, label=odds_ratio$full_sig[20],color="#377eb8",size=10)+
  annotate(geom="text", x=1, y=2.2, label=odds_ratio$full_sig[3],color="#4daf4a",size=10)+annotate(geom="text", x=2, y=2.2, label=odds_ratio$full_sig[9],color="#4daf4a",size=10)+annotate(geom="text", x=3, y=2.2, label=odds_ratio$full_sig[15],color="#4daf4a",size=10)+annotate(geom="text", x=4, y=2.2, label=odds_ratio$full_sig[21],color="#4daf4a",size=10)+
  annotate(geom="text", x=1, y=2.3, label=odds_ratio$full_sig[4],color="#984ea3",size=10)+annotate(geom="text", x=2, y=2.3, label=odds_ratio$full_sig[10],color="#984ea3",size=10)+annotate(geom="text", x=3, y=2.3, label=odds_ratio$full_sig[16],color="#984ea3",size=10)+annotate(geom="text", x=4, y=2.3, label=odds_ratio$full_sig[22],color="#984ea3",size=10)+
  annotate(geom="text", x=1, y=2.4, label=odds_ratio$full_sig[5],color="#ff7f00",size=10)+annotate(geom="text", x=2, y=2.4, label=odds_ratio$full_sig[11],color="#ff7f00",size=10)+annotate(geom="text", x=3, y=2.4, label=odds_ratio$full_sig[17],color="#ff7f00",size=10)+annotate(geom="text", x=4, y=2.4, label=odds_ratio$full_sig[23],color="#ff7f00",size=10)+
  annotate(geom="text", x=1, y=2.5, label=odds_ratio$full_sig[6],color="#ffff33",size=10)+annotate(geom="text", x=2, y=2.5, label=odds_ratio$full_sig[12],color="#ffff33",size=10)+annotate(geom="text", x=3, y=2.5, label=odds_ratio$full_sig[18],color="#ffff33",size=10)+annotate(geom="text", x=4, y=2.5, label=odds_ratio$full_sig[24],color="#ffff33",size=10)
odds_figure_full
dev.off()

#repeat for 508
odds_ratio=subset(pheno_contr_full,pheno_contr_full$Stock=="508")
odds_ratio$full_or=exp(odds_ratio$estimate)
odds_ratio$full_sig=ifelse(odds_ratio$p.value<0.001,"***",ifelse(odds_ratio$p.value<0.01,"**",ifelse(odds_ratio$p.value<0.05,"*","")))
#odds_ratio

pdf("images/full_odds_508.pdf")
odds_figure_full=ggplot(aes(y=full_or,x=vial_letter,group=contrast,col=contrast,shape=contrast),data=odds_ratio)+ggtitle("Full Model Results DGRP_508")+scale_colour_manual(values=col_scheme)+ 
  geom_point(size=3)+ylab("Odds Ratios")+theme_base()+geom_hline(yintercept = 1,linetype="dashed",color="grey")+ylim(0.5,2.75)+theme(axis.text.x = element_text(angle = 45))+
  scale_x_discrete(name="Days post-mating",labels=c("1-2","3-4","5-6","7-8"))+geom_line()+geom_errorbar(aes(ymin=full_or-SE,ymax=full_or+SE))+
  annotate(geom="text", x=1, y=2, label=odds_ratio$full_sig[1],color="#e41a1c",size=10)+annotate(geom="text", x=2, y=2, label=odds_ratio$full_sig[7],color="#e41a1c",size=10)+annotate(geom="text", x=3, y=2, label=odds_ratio$full_sig[13],color="#e41a1c",size=10)+annotate(geom="text", x=4, y=2, label=odds_ratio$full_sig[19],color="#e41a1c",size=10)+
  annotate(geom="text", x=1, y=2.1, label=odds_ratio$full_sig[2],color="#377eb8",size=10)+annotate(geom="text", x=2, y=2.1, label=odds_ratio$full_sig[8],color="#377eb8",size=10)+annotate(geom="text", x=3, y=2.1, label=odds_ratio$full_sig[14],color="#377eb8",size=10)+annotate(geom="text", x=4, y=2.1, label=odds_ratio$full_sig[20],color="#377eb8",size=10)+
  annotate(geom="text", x=1, y=2.2, label=odds_ratio$full_sig[3],color="#4daf4a",size=10)+annotate(geom="text", x=2, y=2.2, label=odds_ratio$full_sig[9],color="#4daf4a",size=10)+annotate(geom="text", x=3, y=2.2, label=odds_ratio$full_sig[15],color="#4daf4a",size=10)+annotate(geom="text", x=4, y=2.2, label=odds_ratio$full_sig[21],color="#4daf4a",size=10)+
  annotate(geom="text", x=1, y=2.3, label=odds_ratio$full_sig[4],color="#984ea3",size=10)+annotate(geom="text", x=2, y=2.3, label=odds_ratio$full_sig[10],color="#984ea3",size=10)+annotate(geom="text", x=3, y=2.3, label=odds_ratio$full_sig[16],color="#984ea3",size=10)+annotate(geom="text", x=4, y=2.3, label=odds_ratio$full_sig[22],color="#984ea3",size=10)+
  annotate(geom="text", x=1, y=2.4, label=odds_ratio$full_sig[5],color="#ff7f00",size=10)+annotate(geom="text", x=2, y=2.4, label=odds_ratio$full_sig[11],color="#ff7f00",size=10)+annotate(geom="text", x=3, y=2.4, label=odds_ratio$full_sig[17],color="#ff7f00",size=10)+annotate(geom="text", x=4, y=2.4, label=odds_ratio$full_sig[23],color="#ff7f00",size=10)+
  annotate(geom="text", x=1, y=2.5, label=odds_ratio$full_sig[6],color="#ffff33",size=10)+annotate(geom="text", x=2, y=2.5, label=odds_ratio$full_sig[12],color="#ffff33",size=10)+annotate(geom="text", x=3, y=2.5, label=odds_ratio$full_sig[18],color="#ffff33",size=10)+annotate(geom="text", x=4, y=2.5, label=odds_ratio$full_sig[24],color="#ffff33",size=10)
odds_figure_full
dev.off()