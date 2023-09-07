# imputation with mice
# mean: method="mean"
# regression: method="norm.predict"
# stochastic regression: method="norm.nob"
# bayesian regression: method="norm"

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

a <- mice(a_1, method="mean", m=1, maxit=1)
a <- complete(a)
write.csv(a, '/exports/igmm/eddie/UK-Biobank-proj76173/julia/imputed/mean/a_1.csv')
