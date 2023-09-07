# K-NN

library(readr)
library(mice)
library(VIM)
library(naniar)
library(visdat)
library(ggcorrplot)
library(finalfit)
library(GGally)
library(misty)
library(plot.matrix)
library(caret)
library(RANN)
library(philentropy)
library(miceadds)

#setwd("/exports/igmm/eddie/UK-Biobank-proj76173/julia/amputed/")
a_1 <- data.frame(read.csv('/exports/igmm/eddie/UK-Biobank-proj76173/julia/amputed/a_1.csv'))
a_1 <- a_1[,-1]

set.seed(1)
#colnames(a_1) <- full_blood_count
knn.model <- preProcess(a_1, "knnImpute")
knn_12 <- predict(knn.model,a_1,na.action = na.pass)
procNames <- data.frame(col = names(knn.model$mean), mean = knn.model$mean, sd = knn.model$std)
for(i in procNames$col){
  knn_12[i] <- knn_12[i]*knn.model$std[i]+knn.model$mean[i]
}
write.csv(knn_12, '/exports/igmm/eddie/UK-Biobank-proj76173/julia/imputed/knn/knn_1.csv')
