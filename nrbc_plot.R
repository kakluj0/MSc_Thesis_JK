# Nucleated Red Blood Count/Percentage plots

library(ggplot2)

setwd("/exports/igmm/datastore/UK-Biobank-proj76173/ava/walks_f874/")
cases <- read_csv('walks_metabolites_CFS_gold_cases.csv')
controls <- read_csv('walks_metabolites_CFS_gold_controls_v2.csv')
cases$f.120010.0.0 <- 1
controls$f.120010.0.0 <- 0
full <- rbind(cases,controls)

full_red <- full[,c('f.30230.0.0','f.120010.0.0','f.22001.0.0')]
full_red$metabolite <- rep('percentage',nrow(full_red))
colnames(full_red) <- c('value','mecfs','sex','metabolite')
full_red2 <- full[,c('f.30170.0.0','f.120010.0.0','f.22001.0.0')]
full_red2$metabolite <- rep('count',nrow(full_red2))
colnames(full_red2) <- c('value','mecfs','sex','metabolite')
full_red <- rbind(full_red,full_red2)
full_red <- full_red[complete.cases(full_red),]
full_red3 <- full_red
full_red3$sex <- rep('combined',nrow(full_red3))
full_red <- rbind(full_red,full_red3)
full_red[which(full_red$sex==0),'sex'] <- rep('female',length(which(full_red$sex==0)))
full_red[which(full_red$sex==1),'sex'] <- rep('male',length(which(full_red$sex==1)))

full_red$mecfs[which(full_red$mecfs==1)] <- rep('Cases',length(which(full_red$mecfs==1)))
p2 <- ggplot(full_red[which(full_red$metabolite=='count'),],aes(x=value,fill=as.factor(mecfs)))+
  geom_histogram(binwidth = 0.01, position="identity")+
  coord_cartesian(ylim = c(0, 100))+
  labs(y='Count')+
  facet_wrap(mecfs~as.factor(sex))+
  theme_minimal()
