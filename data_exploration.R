# load libraries
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

# read data
setwd("/exports/igmm/datastore/UK-Biobank-proj76173/ava/")
cases <- read_csv('metabolites_CFS_gold_cases.csv')
controls <- read_csv('metabolites_CFS_gold_controls_v2.csv')

# full blood count
full_blood_count <- c('f.21022.0.0','f.22001.0.0','f.30000.0.0',
                      'f.30240.0.0','f.30250.0.0','f.30070.0.0',
                      'f.30010.0.0','f.30110.0.0','f.30090.0.0','f.30080.0.0',
                      'f.30230.0.0','f.30170.0.0','f.30200.0.0','f.30140.0.0',
                      'f.30190.0.0','f.30130.0.0','f.30270.0.0','f.30260.0.0',
                      'f.30100.0.0','f.30040.0.0','f.30050.0.0','f.30060.0.0',
                      'f.30180.0.0','f.30120.0.0','f.30280.0.0','f.30290.0.0',
                      'f.30300.0.0','f.30020.0.0','f.30030.0.0','f.30210.0.0',
                      'f.30150.0.0','f.30220.0.0','f.30160.0.0')
colnames <- colnames(cases)
index <- !(colnames%in%full_blood_count)
index <- colnames[index]
index <- index[-1]
biochem <- append(c('f.21022.0.0','f.22001.0.0'),index)
biochem <- biochem[-33]

blood_cases <- cases[full_blood_count]
blood_cases <- blood_cases[rowSums(is.na(blood_cases)) != ncol(blood_cases),] #delete empty rows
blood_controls <- controls[full_blood_count]
blood_controls <- blood_controls[rowSums(is.na(blood_controls)) != ncol(blood_controls),]
complete_blood_cases <- na.omit(blood_cases)
complete_blood_controls <- na.omit(blood_controls)

biochem_cases <- cases[biochem]
biochem_cases <- biochem_cases[rowSums(is.na(biochem_cases)) != ncol(biochem_cases),]
biochem_controls <- controls[biochem]
biochem_controls <- biochem_controls[rowSums(is.na(biochem_controls)) != ncol(biochem_controls),]

#summary
ff_glimpse(blood_cases)
ff_glimpse(blood_controls)
ff_glimpse(biochem_cases)
ff_glimpse(biochem_controls)

# missingness in a row
narow_b <- rowSums(is.na(blood_cases))/33*100
narow_b <- c(sum(narow_b==0),sum(0<narow_b&narow_b<=10),sum(10<narow_b&narow_b<=20),sum(20<narow_b&narow_b<=30),sum(30<narow_b&narow_b<=40),sum(40<narow_b&narow_b<=50),sum(50<narow_b&narow_b<=60),sum(60<narow_b&narow_b<=70),sum(70<narow_b&narow_b<=80),sum(80<narow_b&narow_b<=90),sum(90<narow_b&narow_b<=100))/length(narow_b)
narow_bc <- rowSums(is.na(blood_controls))/33*100
narow_bc <- c(sum(narow_bc==0),sum(0<narow_bc&narow_bc<=10),sum(10<narow_bc&narow_bc<=20),sum(20<narow_bc&narow_bc<=30),sum(30<narow_bc&narow_bc<=40),sum(40<narow_bc&narow_bc<=50),sum(50<narow_bc&narow_bc<=60),sum(60<narow_bc&narow_bc<=70),sum(70<narow_bc&narow_bc<=80),sum(80<narow_bc&narow_bc<=90),sum(90<narow_bc&narow_bc<=100))/length(narow_bc)
narow_c <- rowSums(is.na(biochem_cases))/32*100
narow_c <- c(sum(narow_c==0),sum(0<narow_c&narow_c<=10),sum(10<narow_c&narow_c<=20),sum(20<narow_c&narow_c<=30),sum(30<narow_c&narow_c<=40),sum(40<narow_c&narow_c<=50),sum(50<narow_c&narow_c<=60),sum(60<narow_c&narow_c<=70),sum(70<narow_c&narow_c<=80),sum(80<narow_c&narow_c<=90),sum(90<narow_c&narow_c<=100))/length(narow_c)
narow_cc <- rowSums(is.na(biochem_cases))/32*100
narow_cc <- c(sum(narow_cc==0),sum(0<narow_cc&narow_cc<=10),sum(10<narow_cc&narow_cc<=20),sum(20<narow_cc&narow_cc<=30),sum(30<narow_cc&narow_cc<=40),sum(40<narow_cc&narow_cc<=50),sum(50<narow_cc&narow_cc<=60),sum(60<narow_cc&narow_cc<=70),sum(70<narow_cc&narow_cc<=80),sum(80<narow_cc&narow_cc<=90),sum(90<narow_cc&narow_cc<=100))/length(narow_cc)
narow <- data.frame(narow_b,narow_bc,narow_c,narow_cc)
rownames(narow) <- c('0%','0-10%','10-20%','20-30%','30-40%','40-50%','50-60%','60-70%','70-80%','80-90%','90-100%')
colnames(narow) <- c('Blood cases','Blood contr', 'Biochem cases','Biochem contr')
round(narow,4)

# missingness summary
na.descript(blood_cases, table=TRUE)
na.descript(blood_controls, table=TRUE)
na.descript(biochem_cases,table=TRUE)
na.descript(biochem_controls,table=TRUE)

# visualisation
md.pattern(blood_cases)
md.pattern(blood_controls)
#md.pattern(biochem_cases)
#md.pattern(biochem_controls)

#delete almost empty rows
blood_cases <- cases[full_blood_count]
blood_cases <- blood_cases[which(rowMeans(is.na(blood_cases)) < 0.9),which(colMeans(is.na(blood_cases)) < 0.9) ]
#blood_cases <- blood_cases[rowSums(is.na(blood_cases['f.30000.0.0'])) != ncol(blood_cases['f.30000.0.0']),]
complete_blood_cases <- na.omit(blood_cases)

blood_controls <- controls[full_blood_count]
blood_controls <- blood_controls[which(rowMeans(is.na(blood_controls)) < 0.9), which(colMeans(is.na(blood_controls)) < 0.9)]
#blood_controls <- blood_controls[rowSums(is.na(blood_controls[,c('f.30010.0.0','f.30020.0.0')])) != ncol(blood_controls[,c('f.30010.0.0','f.30020.0.0')]),]
complete_blood_controls <- na.omit(blood_controls)

# mcar test
mcar_test(blood_cases)
mcar_test(blood_controls)

# nullity correlation matrix
## Create data frame indicating missingness by 1
x <- as.data.frame(abs(is.na(blood_cases)))
## Select columns with some (but not all) missing values
y <- x[,sapply(x, sd) > 0]
## Create a correlation matrix: Variables missing together have high correlation
ggcorrplot::ggcorrplot(cor(y), type = "lower")

## Create data frame indicating missingness by 1
x <- as.data.frame(abs(is.na(blood_controls)))
## Select columns with some (but not all) missing values
y <- x[,sapply(x, sd) > 0]
## Create a correlation matrix: Variables missing together have high correlation
ggcorrplot::ggcorrplot(cor(y), type = "lower")

# biochemistry

biochem_controls <- controls[biochem]
biochem_controls <- biochem_controls[which(rowMeans(is.na(biochem_controls)) < 0.9), which(colMeans(is.na(biochem_controls)) < 0.9)]
## Create data frame indicating missingness by 1
x <- as.data.frame(abs(is.na(biochem_controls)))
## Select columns with some (but not all) missing values
y <- x[,sapply(x, sd) > 0]
## Create a correlation matrix: Variables missing together have high correlation
ggcorrplot::ggcorrplot(cor(y), type = "lower")

biochem_cases <- cases[biochem]
biochem_cases <- biochem_cases[which(rowMeans(is.na(biochem_cases)) < 0.9), which(colMeans(is.na(biochem_cases)) < 0.9)]
## Create data frame indicating missingness by 1
x <- as.data.frame(abs(is.na(biochem_cases)))
## Select columns with some (but not all) missing values
y <- x[,sapply(x, sd) > 0]
## Create a correlation matrix: Variables missing together have high correlation
ggcorrplot::ggcorrplot(cor(y), type = "lower")

# mcar test
mcar_test(biochem_cases)
mcar_test(biochem_controls)
