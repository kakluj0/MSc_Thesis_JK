# Point Bi-serial Correaltion

library(tidyverse)
library(broom)

setwd("/exports/igmm/datastore/UK-Biobank-proj76173/ava/")
cases <- read_csv('metabolites_CFS_gold_cases.csv')
controls <- read_csv('metabolites_CFS_gold_controls_v2.csv')
cases$f.120010.0.0 <- 1
controls$f.120010.0.0 <- 0
full <- rbind(cases,controls)

cols <- colnames(full)[!(colnames(full)%in%c("f.eid","f.120010.0.0","f.21022.0.0","f.22001.0.0"))]

full <- as.matrix(full)

est <- rep(NA,length(cols))
ci.l <- rep(NA,length(cols))
ci.u <- rep(NA,length(cols))
name <- rep(NA,length(cols))
pval <- rep(NA,length(cols))
# Point-Biserial Correlation
for (i in 1:(length(cols))){
  a <- cols[i]
  correlation <- cor.test(full[,a],full[,'f.120010.0.0'])
  est[i] <- correlation$estimate
  ci.l[i] <- correlation$conf.int[1]
  ci.u[i] <- correlation$conf.int[2]
  name[i] <- substr(a,1,7)
  pval[i] <- correlation$p.value
}

data <- cbind(est,ci.l,ci.u,name,pval)
data <- as.data.frame(data)

data$est <- as.numeric(data$est)
data$ci.l <- as.numeric(data$ci.l)
data$ci.u <- as.numeric(data$ci.u)

data$name <- factor(data$name, levels = data$name[order(data$est)])

ggplot(data)+
  geom_point(aes(x=name,y=est))+
  geom_errorbar(aes(x=name,ymin=ci.l,ymax=ci.u))+
  labs(x=NULL, y='Point-Biserial Correlation')+
  theme_minimal()+
  theme(axis.text.x = element_text(angle=90))

p_val_combined <- pval
BH_combined <- p.adjust(pval,method="BH")
p_val_combined <- cbind(name,p_val_combined,BH_combined)
p_val_combined <- merge(p_val_combined,names_all,by='name',allx=TRUE)

########################################### females
cases <- read_csv('metabolites_CFS_gold_cases.csv')
controls <- read_csv('metabolites_CFS_gold_controls_v2.csv')
cases$f.120010.0.0 <- 1
controls$f.120010.0.0 <- 0
full <- rbind(cases,controls)

cols <- colnames(full)[!(colnames(full)%in%c("f.eid","f.120010.0.0","f.21022.0.0","f.22001.0.0"))]
full <- full[which(full$f.22001.0.0==0),]

full <- as.matrix(full)

est <- rep(NA,length(cols))
ci.l <- rep(NA,length(cols))
ci.u <- rep(NA,length(cols))
name <- rep(NA,length(cols))
pval <- rep(NA,length(cols))
# Point-Biserial Correlation
for (i in 1:(length(cols))){
  a <- cols[i]
  correlation <- cor.test(full[,a],full[,'f.120010.0.0'])
  est[i] <- correlation$estimate
  ci.l[i] <- correlation$conf.int[1]
  ci.u[i] <- correlation$conf.int[2]
  name[i] <- substr(a,1,7)
  pval[i] <- correlation$p.value
}

data <- cbind(est,ci.l,ci.u,name,pval)
data <- as.data.frame(data)

data$est <- as.numeric(data$est)
data$ci.l <- as.numeric(data$ci.l)
data$ci.u <- as.numeric(data$ci.u)

data$name <- factor(data$name, levels = data$name[order(data$est)])

ggplot(data)+
  geom_point(aes(x=name,y=est))+
  geom_errorbar(aes(x=name,ymin=ci.l,ymax=ci.u))+
  labs(x=NULL, y='Point-Biserial Correlation',title='Males')+
  theme_minimal()+
  theme(axis.text.x = element_text(angle=90))

p_val_females <- pval
BH_females <- p.adjust(pval,method="BH")
p_val_females <- cbind(name,p_val_females,BH_females)
p_val_females <- merge(p_val_females,names_all,by='name',all.x=TRUE)

pval_corr <- merge(p_val_combined,p_val_females,by='name')
pval_corr <- merge(pval_corr, p_val_males,by='name',all.x=TRUE)
write.csv(pval_corr,'/exports/igmm/datastore/UK-Biobank-proj76173/julia/reverseonestepcate/pvals_corr.csv')
