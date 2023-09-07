# plot shift interventions
setwd("/exports/igmm/datastore/UK-Biobank-proj76173/julia/shift_onestep/")


name <- c('f.21022.0.0','f.22001.0.0','f.30000.0.0',
             'f.30240.0.0','f.30250.0.0','f.30070.0.0',
             'f.30010.0.0','f.30110.0.0','f.30090.0.0','f.30080.0.0',
             'f.30230.0.0','f.30170.0.0','f.30200.0.0','f.30140.0.0',
             'f.30190.0.0','f.30130.0.0','f.30270.0.0','f.30260.0.0',
             'f.30100.0.0','f.30040.0.0','f.30050.0.0','f.30060.0.0',
             'f.30180.0.0','f.30120.0.0','f.30280.0.0','f.30290.0.0',
             'f.30300.0.0','f.30020.0.0','f.30030.0.0','f.30210.0.0',
             'f.30150.0.0','f.30220.0.0','f.30160.0.0')

outcome <- c('21022 Age at recruitment',
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

name <- c( "f.30600.0.0","f.30610.0.0",  "f.30620.0.0" , "f.30630.0.0", 
              "f.30640.0.0" , "f.30650.0.0",  "f.30660.0.0" , "f.30670.0.0",  "f.30680.0.0", 
              "f.30690.0.0",  "f.30700.0.0",  "f.30710.0.0",  "f.30720.0.0" , "f.30730.0.0" ,
              "f.30740.0.0",  "f.30750.0.0",  "f.30760.0.0" , "f.30770.0.0",  "f.30780.0.0" ,
              "f.30790.0.0" , "f.30800.0.0",  "f.30810.0.0",  "f.30820.0.0" , "f.30830.0.0" ,
              "f.30840.0.0" , "f.30850.0.0",  "f.30860.0.0" , "f.30870.0.0",  "f.30880.0.0" ,
              "f.30890.0.0")

outcome <- c(
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

# combined metabolites 
m1 <- read_csv('msm_combined.csv')
p1 <- read_csv('points_combined.csv')
m1 <- merge(m1, names_all, by='name', all.x=TRUE)
p1 <- merge(p1, names_all, by='name', all.x=TRUE)

m1$group <- rep('Combined',nrow(m1))
m1 <- m1[which(m1$outcome!='f.eid'),]
p1$group <- rep('Combined',nrow(p1))
p1 <- p1[which(p1$outcome!='f.eid'),]

# female metabolites 
m2 <- read_csv('msm_females.csv')
p2 <- read_csv('points_females.csv')
m2 <- merge(m2, names_all, by='name', all.x=TRUE)
p2 <- merge(p2, names_all, by='name', all.x=TRUE)

m2$group <- rep('Females',nrow(m2))
m2 <- m2[which(m2$outcome!='f.eid'),]
p2$group <- rep('Females',nrow(p2))
p2 <- p2[which(p2$outcome!='f.eid'),]

# male metabolites 
m3 <- read_csv('msm_females.csv')
p3 <- read_csv('points_females.csv')
m3 <- merge(m3, names_all, by='name', all.x=TRUE)
p3 <- merge(p3, names_all, by='name', all.x=TRUE)

m3$group <- rep('Males',nrow(m3))
m3 <- m3[which(m3$outcome!='f.eid'),]
p3$group <- rep('Males',nrow(p3))
p3 <- p3[which(p3$outcome!='f.eid'),]

msm <- rbind(m1,m2,m3)
points <- rbind(p1,p2,p3)

all <- unique(points$name)
all_n <- unique(points$outcome)


for (k in 1:length(all)){
b <- all[k]
b1 <- all_n[k]
p <- points[which(points$name==b),]
m <- msm[which(msm$name==b),]

i <- m[which(m$param=='(Intercept)'),c('param_est','group','param_se','p_value')]
colnames(i) <- c('intercept','group','sei','pi')
s <- m[which(m$param=='delta'),c('param_est','group','param_se','p_value')]
colnames(s) <- c('slope','group','ses','ps')
m <- merge(i,s,by='group')

supp.labs <- c(paste('Combined\nIntercept:',as.character(round(as.numeric(m$intercept[which(m$group=='Combined')]),4)),
                     '+-',as.character(round(as.numeric(m$sei[which(m$group=='Combined')]),5)),
                     '\np-value=',as.character(round(as.numeric(m$pi[which(m$group=='Combined')]),5)),
                     '\nSlope:', as.character(round(as.numeric(m$slope[which(m$group=='Combined')]),5)),
                     '+-', as.character(round(as.numeric(m$ses[which(m$group=='Combined')]),5)),
                     '\np-value=',as.character(round(as.numeric(m$ps[which(m$group=='Combined')]),5))),
               paste('Males\nIntercept:',as.character(round(as.numeric(m$intercept[which(m$group=='Males')]),5)),
                     '+-',as.character(round(as.numeric(m$sei[which(m$group=='Males')]),5)),
                     '\np-value=',as.character(round(as.numeric(m$pi[which(m$group=='Males')]),5)),
                     '\nSlope:', as.character(round(as.numeric(m$slope[which(m$group=='Males')]),5)),
                     '+-', as.character(round(as.numeric(m$ses[which(m$group=='Males')]),5)),
                     '\np-value=',as.character(round(as.numeric(m$ps[which(m$group=='Males')]),5))),
               paste('Females\nIntercept:',as.character(round(as.numeric(m$intercept[which(m$group=='Females')]),5)),
                     '+-',as.character(round(as.numeric(m$sei[which(m$group=='Females')]),5)),
                     '\np-value=',as.character(round(as.numeric(m$pi[which(m$group=='Females')]),5)),
                     '\nSlope:', as.character(round(as.numeric(m$slope[which(m$group=='Females')]),5)),
                     '+-', as.character(round(as.numeric(m$ses[which(m$group=='Females')]),5)),
                     '\np-value=',as.character(round(as.numeric(m$ps[which(m$group=='Females')]),5)))
)
names(supp.labs) <- c('Combined', "Males",'Females')


ggplot(p) +
  geom_abline(data=m, aes(intercept = intercept, slope = slope, alpha=0.5)) +
  geom_point(aes(x=delta,y=psi))+
  geom_errorbar(aes(x=delta,ymin=ci_lwr, ymax=ci_upr))+
  labs(title=b1, x="delta", y = "TSM") +
  facet_grid(~group,labeller = labeller(group = supp.labs),scale='free_x')+
  #facet_grid(~group)+
  theme_minimal()+ scale_alpha(guide = 'none')

dir <- paste('/exports/igmm/datastore/UK-Biobank-proj76173/julia/shift_onestep/plots/',b,'.pdf',sep='')
ggsave(dir,width = 8, height = 6)
}
