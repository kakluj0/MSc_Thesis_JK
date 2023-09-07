# plot ATE
library(ggplot2)
library(ggforce)
library(lemon)
library(patchwork)

setwd("/exports/igmm/datastore/UK-Biobank-proj76173/ava")

#names nmr
names_nmr <- read_csv('all_variables_names.csv')
names_nmr <- names_nmr[,-1]
colnames(names_nmr) <- c('outcome','name')
names_nmr$name <- paste(names_nmr$outcome,names_nmr$name)

#names full blood count
outcome <- c('f.21022.0.0','f.22001.0.0','f.30000.0.0',
          'f.30240.0.0','f.30250.0.0','f.30070.0.0',
          'f.30010.0.0','f.30110.0.0','f.30090.0.0','f.30080.0.0',
          'f.30230.0.0','f.30170.0.0','f.30200.0.0','f.30140.0.0',
          'f.30190.0.0','f.30130.0.0','f.30270.0.0','f.30260.0.0',
          'f.30100.0.0','f.30040.0.0','f.30050.0.0','f.30060.0.0',
          'f.30180.0.0','f.30120.0.0','f.30280.0.0','f.30290.0.0',
          'f.30300.0.0','f.30020.0.0','f.30030.0.0','f.30210.0.0',
          'f.30150.0.0','f.30220.0.0','f.30160.0.0')

name <- c('21022 Age at recruitment',
             '22001 Genetic sex',
             '30000 White blood cell (leukocyte) count',
             '30240 Reticulocyte percentage',
             '30250 Reticulocyte count',
             '30070 Red blood cell (erythrocyte) distribution width',
             '30010 Red blood cell (erythrocyte) count',
             '30110 Platelet distribution width',
             '30090 Platelet crit',
             '30080 Platelet count',
             '30230 Nucleated red blood cell percentage',
             '30170 Nucleated red blood cell count',
             '30200 Neutrophil percentage',
             '30140 Neutrophil count',
             '30190 Monocte percentage',
             '30130 Monocyte count',
             '30270 Mean sphered cell volume',
             '30260 Mean reticulocyte volume',
             '30100 Mean platelet (thrombocyte) volume',
             '30040 Mean corpuscular volume',
             '30050 Mean corpuscular haemoglobin',
             '30060 Mean corpuscular haemoglobin concentration',
             '30180 Lymphocyte percentage',
             '30120 Lymphocyte count',
             '30280 Immature reticulocyte fraction',
             '30290 High light scatter reticulocyte percentage',
             '30300 High light scatter reticulocyte count',
             '30020 Haemoglobin concentration',
             '30030 Haematocrit percentage',
             '30210 Eosinophill percentage',
             '30150 Eosinophill count',
             '30220 Basophill percentage',
             '30160 Basophill count')

fbc <- cbind(outcome,name)

outcome <- c( "f.30600.0.0","f.30610.0.0",  "f.30620.0.0" , "f.30630.0.0", 
           "f.30640.0.0" , "f.30650.0.0",  "f.30660.0.0" , "f.30670.0.0",  "f.30680.0.0", 
           "f.30690.0.0",  "f.30700.0.0",  "f.30710.0.0",  "f.30720.0.0" , "f.30730.0.0" ,
           "f.30740.0.0",  "f.30750.0.0",  "f.30760.0.0" , "f.30770.0.0",  "f.30780.0.0" ,
           "f.30790.0.0" , "f.30800.0.0",  "f.30810.0.0",  "f.30820.0.0" , "f.30830.0.0" ,
           "f.30840.0.0" , "f.30850.0.0",  "f.30860.0.0" , "f.30870.0.0",  "f.30880.0.0" ,
           "f.30890.0.0")

name <- c(
  '30600 Albumin',
  '30610 Alkaline phosphatase',
  '30620 Alanine aminotransferase',
  '30630 Apolipoprotein A',
  '30640 Apolipoprotein B',
  '30650 Aspartate aminotransferase',
  '30660 Direct bilirubin',
  '30670 Urea',
  '30680 Calcium',
  '30690 Cholesterol',
  '30700 Creatinine',
  '30710 C-reactive protein',
  '30720 Cystatin C',
  '30730 Gamma glutamyltransferase',
  '30740 Glucose',
  '30750 Glycated haemoglobin (HbA1c)',
  '30760 HDL cholesterol',
  '30770 IGF-1',
  '30780 LDL direct',
  '30790 Lipoprotein A',
  '30800 Oestradiol',
  '30810 Phosphate',
  '30820 Rheumatoid factor',
  '30830 SHBG',
  '30840 Total bilirubin',
  '30850 Testosterone',
  '30860 Total protein',
  '30870 Triglycerides',
  '30880 Urate',
  '30890 Vitamin D'
)

bio <- cbind(outcome,name)
names_all <- rbind(fbc,bio)

setwd("/exports/igmm/datastore/UK-Biobank-proj76173/julia/walks")
# combined metabolites 
a1 <- read_csv('combined.csv')
a1 <- merge(a1, names_all, by='outcome', all.x=TRUE)
#a1 <- read_csv('nmr_combined.csv')
#a1 <- merge(a1, names_nmr, by='outcome', all.x=TRUE)
a1$group <- rep('Combined',nrow(a1))
a1 <- a1[which(a1$outcome!='f.eid'),]
a1 <- a1[which(a1$parameter=='E{Y(1)-Y(0)}'),]

# females
a2 <- read_csv('females.csv')
a2 <- merge(a2, names_all, by='outcome', all.x=TRUE)
#a2 <- read_csv('nmr_females.csv')
#a2 <- merge(a2, names_nmr, by='outcome', all.x=TRUE)
a2$group <- rep('Females',nrow(a2))
a2 <- a2[which(a2$outcome!='f.eid'),]
a2 <- a2[which(a2$parameter=='E{Y(1)-Y(0)}'),]

# males
a3 <- read_csv('males.csv')
a3 <- merge(a3, names_all, by='outcome', all.x=TRUE)
#a3 <- read_csv('nmr_males.csv')
#a3 <- merge(a3, names_nmr, by='outcome', all.x=TRUE)
a3$group <- rep('Males',nrow(a3))
a3 <- a3[which(a3$outcome!='f.eid'),]
a3 <- a3[which(a3$parameter=='E{Y(1)-Y(0)}'),]

metabolites <- rbind(a1,a2,a3)

# palette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# p-values
a1$pvalue <- 2*pnorm(-abs(a1$est/a1$se))
a2$pvalue <- 2*pnorm(-abs(a2$est/a2$se))
a3$pvalue <- 2*pnorm(-abs(a3$est/a3$se))

# BH p-vals
a1$BH_pvalue <- p.adjust(a1$pvalue,method="BH")
a2$BH_pvalue <- p.adjust(a2$pvalue,method="BH")
a3$BH_pvalue <- p.adjust(a3$pvalue,method="BH")

# significant BH p-vals
a1$asterisk_f <- rep('',nrow(a1))
a1$asterisk_m <- rep('',nrow(a1))
a1$asterisk_c <- rep('',nrow(a1))
a2$asterisk_f <- rep('',nrow(a2))
a2$asterisk_m <- rep('',nrow(a2))
a2$asterisk_c <- rep('',nrow(a2))
a3$asterisk_f <- rep('',nrow(a3))
a3$asterisk_m <- rep('',nrow(a3))
a3$asterisk_c <- rep('',nrow(a3))

a1$asterisk_f[which(a1$BH_pvalue<=0.05)] <- rep('*',length(which(a1$BH_pvalue<=0.05)))
a2$asterisk_m[which(a2$BH_pvalue<=0.05)] <- rep('*',length(which(a2$BH_pvalue<=0.05)))
a3$asterisk_c[which(a3$BH_pvalue<=0.05)] <- rep('*',length(which(a3$BH_pvalue<=0.05)))

metabolites <- rbind(a1,a2,a3)

# names of significant metabolites
# females:
f_plot <- metabolites[which(metabolites$asterisk_f=='*'),'name']
# combined:
c_plot <- metabolites[which(metabolites$asterisk_c=='*'),'name']
# males:
m_plot <- metabolites[which(metabolites$asterisk_m=='*'),'name']

# significant metabolites in groups
metabolites$sig <- rep(NA,nrow(metabolites))

# females + males:
fm_plot <- f_plot[which(f_plot%in%m_plot)]

# females + males + combined:
fmc_plot <- fm_plot[which(fm_plot%in%c_plot)]
metabolites$sig[which(metabolites$name%in%fmc_plot)] <- rep('fmc',length(which(metabolites$name%in%fmc_plot)))

# females + males:
fm_plot <- fm_plot[which(!(fm_plot%in%fmc_plot))]
metabolites$sig[which(metabolites$name%in%fm_plot)] <- rep('fm',length(which(metabolites$name%in%fm_plot)))

# females + combined:
fc_plot <- c_plot[which(c_plot%in%f_plot)]
fc_plot <- fc_plot[which(!(fc_plot%in%fmc_plot))]
metabolites$sig[which(metabolites$name%in%fc_plot)] <- rep('fc',length(which(metabolites$name%in%fc_plot)))

# males + combined:
#mc_plot <- c_plot[which(c_plot%in%m_plot)]
#mc_plot <- mc_plot[which(!(mc_plot%in%fmc_plot))]

# female only:
f_only_plot <- f_plot[which(!(f_plot%in%c(c_plot,m_plot)))]
metabolites$sig[which(metabolites$name%in%f_only_plot)] <- rep('f',length(which(metabolites$name%in%f_only_plot)))

# male only: 
m_only_plot <- m_plot[which(!(m_plot%in%c(f_plot,c_plot)))]
metabolites$sig[which(metabolites$name%in%m_only_plot)] <- rep('m',length(which(metabolites$name%in%m_only_plot)))

# combined only:
c_only_plot <- c_plot[which(!(c_plot%in%c(f_plot,m_plot)))]
metabolites$sig[which(metabolites$name%in%c_only_plot)] <- rep('c',length(which(metabolites$name%in%c_only_plot)))

metab <- metabolites[which(complete.cases(metabolites)),]

# zeros
scaleFUN <- function(x) x

# plot with intervals
p1 <- ggplot(metab) + 
  geom_vline(xintercept = 0,alpha=0.1)+
  #geom_vline(xintercept = 1.26,alpha=0.1)+
  #geom_vline(xintercept = 1.44,alpha=0.1)+
  geom_errorbarh(aes(y=name,xmax = ci.ul, xmin = ci.ll,colour=group),position=position_dodge(width = 1.1),height = 0.7,lwd=1)+
  geom_point(aes(x=est,y=name,colour=group),position=position_dodge(width = 1.1),size=2.2)+
  ##geom_text_repel(aes(x=1,y=name,label=asterisk,color=group),min.segment.length = Inf, box.padding = 0.1,nudge_x =0.1,max.overlaps = 20,direction='x',size=10,hjust=0,vjust=0)+
  #annotate("rect", xmin = 1.26, xmax = 1.44,ymin=-Inf,ymax=Inf,fill = "white")+
  #geom_text(aes(x=ci.ul+0.05,y=name,label=asterisk,colour=group),position=position_dodge(width = 0.8),size=8)+
  #geom_text(aes(x=1.3,y=name,label=asterisk_f,colour=group),size=10)+
  #geom_text(aes(x=1.35,y=name,label=asterisk_m,colour=group),size=10)+
  #geom_text(aes(x=1.4,y=name,label=asterisk_c,colour=group),size=10)+
  #geom_text(aes(x=ci.ul,y=name,label=asterisk,colour=group),size=7,vjust=0.8,hjust=-0.8,position=position_dodge(width = 0.8))+
  labs(x="", y = "",title = '',caption='')+
  facet_col(factor(sig, levels=c('fmc','fm','fc','f','m','c')) ~ name ,scales = "free",space = "free")+
  #coord_cartesian(xlim = c(-1.5, 1.5))+
  theme_minimal()+
  scale_colour_manual(values=cbPalette,name='Group')+
  guides( color = guide_legend(override.aes = aes(label = "")))+
  theme(strip.text.x = element_blank(),axis.text.y=element_text(size=15),legend.text = element_text(size=15),legend.title = element_text(size=17))+
  scale_x_symmetric(mid=0,labels = scaleFUN)

# plot with asterisks
p2 <- ggplot(metab) + 
  geom_vline(xintercept = 1.26,alpha=0.1)+
  geom_vline(xintercept = 1.44,alpha=0.1)+
  ##geom_text_repel(aes(x=1,y=name,label=asterisk,color=group),min.segment.length = Inf, box.padding = 0.1,nudge_x =0.1,max.overlaps = 20,direction='x',size=10,hjust=0,vjust=0)+
  annotate("rect", xmin = 1.26, xmax = 1.44,ymin=-Inf,ymax=Inf,fill = "white")+
  #geom_text(aes(x=ci.ul+0.05,y=name,label=asterisk,colour=group),position=position_dodge(width = 0.8),size=8)+
  geom_text(aes(x=1.3,y=name,label=asterisk_f,colour=group),size=10,vjust=0.8)+
  geom_text(aes(x=1.35,y=name,label=asterisk_m,colour=group),size=10,vjust=0.8)+
  geom_text(aes(x=1.4,y=name,label=asterisk_c,colour=group),size=10,vjust=0.8)+
  #geom_text(aes(x=ci.ul,y=name,label=asterisk,colour=group),size=7,vjust=0.8,hjust=-0.8,position=position_dodge(width = 0.8))+
  labs(x="", y = "",title = '',caption='')+
  facet_col(factor(sig, levels=c('fmc','fm','fc','f','m','c')) ~ name ,scales = "free",space = "free")+
  #coord_cartesian(xlim = c(-1.5, 1.5))+
  theme_minimal()+
  scale_colour_manual(values=cbPalette,name='Group')+
  guides( color = guide_legend(override.aes = aes(label = "")))+
  theme(strip.text.x = element_blank(),legend.text = element_blank(),legend.title = element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank())
  #scale_x_symmetric(mid=0,labels = scaleFUN)

# plot
p1+p2 + plot_layout(widths = c(25,2),guides='collect')
