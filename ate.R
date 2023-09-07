library(data.table)
library(haldensify)
library(sl3)
library(tmle3)
library(tmle3shift)
library(devtools)
library(npcausal)
library("dplyr")

#cases <- data.frame(read.csv('/exports/igmm/eddie/UK-Biobank-proj76173/julia/walks_data/walks_metabolites_CFS_gold_cases.csv'))
#controls <- data.frame(read.csv('/exports/igmm/eddie/UK-Biobank-proj76173/julia/walks_data/walks_metabolites_CFS_gold_controls_v2.csv'))

cases <- data.frame(read.csv('/exports/igmm/eddie/UK-Biobank-proj76173/sjoerd/metabolites_CFS_gold_cases.csv'))
controls <- data.frame(read.csv('/exports/igmm/eddie/UK-Biobank-proj76173/sjoerd/metabolites_CFS_gold_controls_v2.csv'))
head(cases)
head(controls)

#setwd("/exports/igmm/datastore/UK-Biobank-proj76173/ava/")
#cases <- read_csv('metabolites_CFS_gold_cases.csv')
#controls <- read_csv('metabolites_CFS_gold_controls_v2.csv')

############################################################
#NMR
#cases <- data.frame(read.csv('/exports/igmm/eddie/UK-Biobank-proj76173/julia/walks_data/walks_NMR_metabolomics_metabolites_CFS_gold_cases.csv'))
#controls <- data.frame(read.csv('/exports/igmm/eddie/UK-Biobank-proj76173/julia/walks_data/walks_NMR_metabolomics_metabolites_CFS_gold_controls_v2.csv'))
#cases <- select(cases,-starts_with('f.3'))
#controls <- select(controls,-starts_with('f.3'))
#############################################################

cases$f.120010.0.0 <- 1
controls$f.120010.0.0 <- 0

# cate
#sex_cases <- which(cases$f.22001.0.0==1)
#sex_controls <- which(controls$f.22001.0.0==1)
#cases <- cases[sex_cases,]
#controls <- controls[sex_controls,]
##cases <- select(cases,-c(f.22001.0.0))
##controls <- select(controls,-c(f.22001.0.0))

full <- rbind(cases,controls)
# permute ME column
#index <- sample(1:nrow(full),nrow(full))
#full$f.120010.0.0[index] <- full$f.120010.0.0

var <- colnames(full)
var <- var[-(var=='f.eid')]
WT <- c("f.21022.0.0", "f.22001.0.0","f.120010.0.0",'f.874.0.0')
var <- var[!(var%in%WT)]

i=1
outcome <- var[i]
# y is the feature (blood trait), a is ME status, and x is sex and age

##mini <- full[,c("f.21022.0.0", "f.22001.0.0","f.120010.0.0",'f.874.0.0',outcome)]
#mini <- full[,c("f.21022.0.0", "f.120010.0.0",'f.874.0.0',outcome)] #cate

mini <- full[,c('f.21022.0.0','f.120010.0.0','f.22001.0.0',outcome)] # no walks

mini <- mini[complete.cases(mini),]
a <- mini[,"f.120010.0.0"]

##x <- mini[,c("f.21022.0.0", "f.22001.0.0",'f.874.0.0')]
#x <- mini[,c("f.21022.0.0",'f.874.0.0')] #cate
x <- mini[,c('f.21022.0.0','f.22001.0.0')] # no walk

y <- mini[,outcome]
data <- ate(y,a,x, sl.lib=c("SL.mean","SL.gam","SL.glmnet", "SL.glm.interaction", "SL.ranger"))
data <- data$res
data$outcome <- outcome
data$cases <- sum(mini$f.120010.0.0==1)
data$controls <- sum(mini$f.120010.0.0==0)

for (i in 2:length(var)){
  outcome <- var[i]

  ##mini <- full[,c("f.21022.0.0", "f.22001.0.0","f.120010.0.0",outcome,'f.874.0.0')]
  #mini <- full[,c("f.21022.0.0", "f.120010.0.0",'f.874.0.0',outcome)] #cate

  mini <- full[,c('f.21022.0.0','f.120010.0.0','f.22001.0.0',outcome)] # no walks

  mini <- mini[complete.cases(mini),]
  a <- mini[,"f.120010.0.0"]

  ##x <- mini[,c("f.21022.0.0",'f.874.0.0')] #cate
  #x <- mini[,c("f.21022.0.0", "f.22001.0.0",'f.874.0.0')]

  x <- mini[,c('f.21022.0.0','f.22001.0.0')] # no walk

#  if (i!=18){
  y <- mini[,outcome]
  cat('iter',i)
  try(onestep <- ate(y,a,x, sl.lib=c("SL.mean","SL.gam","SL.glmnet", "SL.glm.interaction", "SL.ranger")))
  onestep <- onestep$res
  onestep$outcome <- outcome
  onestep$cases <- sum(mini$f.120010.0.0==1)
  onestep$controls <- sum(mini$f.120010.0.0==0)

  # Save data
  data <- rbind(data,onestep)
#  }
}

#dir <- '/exports/igmm/eddie/UK-Biobank-proj76173/julia/walks/p_combined.csv'
dir <- '/exports/igmm/eddie/UK-Biobank-proj76173/julia/nobatches/all_no_controls.csv'
write.csv(data,dir)
