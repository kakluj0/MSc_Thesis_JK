# Random Forest

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
library(missRanger)

setwd("/exports/igmm/eddie/UK-Biobank-proj76173/julia/amputed/")
a_1 <- data.frame(read.csv('a_1.csv'))
a_1 <- a_1[,-1]
non_miss <- rowSums(!is.na(a_1))
a_1 <- missRanger(a_1, num.trees = 100, pmm.k = 3, seed = 5, verbose = 0, case.weights = non_miss,splitrule = "extratrees")
write.csv(a_1, '/exports/igmm/eddie/UK-Biobank-proj76173/julia/imputed/forest/a_1.csv')

# for larger datasets (simpulated datasets 7-15)
non_miss <- rowSums(!is.na(a_7))
a_1 <- missRanger(a_7, num.trees = 50, pmm.k = 3, seed = 5, verbose = 0, case.weights = non_miss,splitrule = "extratrees",sample.fraction=0.1)
write.csv(a_1, '/exports/igmm/eddie/UK-Biobank-proj76173/julia/imputed/forest/a_7.csv')
