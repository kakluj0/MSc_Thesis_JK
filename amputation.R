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

#setwd("/exports/igmm/datastore/UK-Biobank-proj76173/ava/")
#cases <- read_csv('metabolites_CFS_gold_cases.csv')
#controls <- read_csv('metabolites_CFS_gold_controls_v2.csv')

cases <- data.frame(read.csv('../../sjoerd/metabolites_CFS_gold_cases.csv'))
controls <- data.frame(read.csv('../../sjoerd/metabolites_CFS_gold_controls_v2.csv'))

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

#delete almost empty rows
blood_cases <- cases[full_blood_count]
blood_cases <- blood_cases[which(rowMeans(is.na(blood_cases)) < 0.9),which(colMeans(is.na(blood_cases)) < 0.9) ]
#blood_cases <- blood_cases[rowSums(is.na(blood_cases['f.30000.0.0'])) != ncol(blood_cases['f.30000.0.0']),]
complete_blood_cases <- na.omit(blood_cases)

#delete almost empty rows
blood_controls <- controls[full_blood_count]
blood_controls <- blood_controls[which(rowMeans(is.na(blood_controls)) < 0.9), which(colMeans(is.na(blood_controls)) < 0.9)]
#blood_controls <- blood_controls[rowSums(is.na(blood_controls[,c('f.30010.0.0','f.30020.0.0')])) != ncol(blood_controls[,c('f.30010.0.0','f.30020.0.0')]),]
complete_blood_controls <- na.omit(blood_controls)


####### AMPUTATION ###################
####################################################################################
prop <- 0.05
p_cases <- md.pattern(blood_cases)
freq <- as.numeric(row.names(p_cases)[2:5])/sum(as.numeric(row.names(p_cases)[2:5]))

#MAR
set.seed(1)
amputed_cases_1_MAR <- ampute(complete_blood_cases, prop=prop, bycases=FALSE, patterns=p_cases[2:5,1:33],
                              freq=freq)
a_blood_cases_MAR <- amputed_cases_1_MAR$amp
write.csv(a_blood_cases_MAR, '/exports/igmm/eddie/UK-Biobank-proj76173/julia/amputed/a_1.csv')

#MNAR
set.seed(1)
amputed_cases_1_MNAR <- ampute(complete_blood_cases, prop=prop, bycases=FALSE, patterns=p_cases[2:5,1:33],
                               freq=freq, mech = "MNAR")
a_blood_cases_MNAR <- amputed_cases_1_MNAR$amp
write.csv(a_blood_cases_MNAR, '/exports/igmm/eddie/UK-Biobank-proj76173/julia/amputed/a_2.csv')

#MCAR&MAR
# ampute the complete data once for every mechanism
set.seed(1)
ampdata1 <- ampute(complete_blood_cases, patterns = p_cases[2:5,1:33], prop = 0.6, mech = "MAR")$amp
ampdata2 <- ampute(complete_blood_cases, patterns = p_cases[2:5,1:33], prop = 0.6, mech = "MCAR")$amp

# create a random allocation vector
# use the prob argument to specify how much of each mechanism should be created
# here, 0.5 of the missingness should be MAR and 0.5 should be MCAR
indices <- sample(x = c(1, 2), size = nrow(complete_blood_cases), 
                  replace = TRUE, prob = c(0.5, 0.5))

# create an empty data matrix
# fill this matrix with values from either of the two amputed datasets
a_blood_cases_MCAR <- matrix(NA, nrow = nrow(complete_blood_cases), ncol = ncol(complete_blood_cases))
a_blood_cases_MCAR[indices == 1, ] <- as.matrix(ampdata1[indices == 1, ])
a_blood_cases_MCAR[indices == 2, ] <- as.matrix(ampdata2[indices == 2, ])
write.csv(a_blood_cases_MCAR, '/exports/igmm/eddie/UK-Biobank-proj76173/julia/amputed/a_3.csv')

###############################################################################
sam <- sample(1:nrow(complete_blood_controls), nrow(blood_cases))
control_small <- complete_blood_controls[ sam,]
write.csv(control_small, '/exports/igmm/eddie/UK-Biobank-proj76173/julia/amputed/c_1.csv')

prop <- 0.05
p_cases <- md.pattern(blood_cases)
freq <- as.numeric(row.names(p_cases)[2:5])/sum(as.numeric(row.names(p_cases)[2:5]))

#MAR
set.seed(1)
amputed_controls_1_MAR <- ampute(control_small, prop=prop, bycases=FALSE, patterns=p_cases[2:5,1:33],
                              freq=freq)
a_blood_controls_MAR <- amputed_controls_1_MAR$amp
write.csv(a_blood_controls_MAR, '/exports/igmm/eddie/UK-Biobank-proj76173/julia/amputed/a_4.csv')

#MNAR
set.seed(1)
amputed_controls_1_MNAR <- ampute(control_small, prop=prop, bycases=FALSE, patterns=p_cases[2:5,1:33],
                               freq=freq, mech = "MNAR")
a_blood_controls_MNAR <- amputed_controls_1_MNAR$amp
write.csv(a_blood_controls_MNAR, '/exports/igmm/eddie/UK-Biobank-proj76173/julia/amputed/a_5.csv')

#MCAR&MAR
# ampute the complete data once for every mechanism
set.seed(1)
ampdata1 <- ampute(control_small, patterns = p_cases[2:5,1:33], prop = 0.6, mech = "MAR")$amp
ampdata2 <- ampute(control_small, patterns = p_cases[2:5,1:33], prop = 0.6, mech = "MCAR")$amp

# create a random allocation vector
# use the prob argument to specify how much of each mechanism should be created
# here, 0.5 of the missingness should be MAR and 0.5 should be MCAR
indices <- sample(x = c(1, 2), size = nrow(control_small), 
                  replace = TRUE, prob = c(0.5, 0.5))

# create an empty data matrix
# fill this matrix with values from either of the two amputed datasets
a_blood_controls_MCAR <- matrix(NA, nrow = nrow(control_small), ncol = ncol(control_small))
a_blood_controls_MCAR[indices == 1, ] <- as.matrix(ampdata1[indices == 1, ])
a_blood_controls_MCAR[indices == 2, ] <- as.matrix(ampdata2[indices == 2, ])
write.csv(a_blood_controls_MCAR, '/exports/igmm/eddie/UK-Biobank-proj76173/julia/amputed/a_6.csv')

################################################################################
p <- md.pattern(biochem_cases)
prop <- 0.07
freq <- as.numeric(row.names(p)[2:159])/sum(as.numeric(row.names(p)[2:159]))

#MAR
set.seed(1)
amputed_controls_2_MAR <- ampute(complete_blood_controls[,-2], prop=prop, bycases=FALSE, 
                              patterns=p[2:159,1:32],freq=freq)
a_controls_2_MAR <- amputed_controls_2_MAR$amp
write.csv(a_controls_2_MAR, '/exports/igmm/eddie/UK-Biobank-proj76173/julia/amputed/a_7.csv')

#MNAR
set.seed(1)
amputed_controls_2_MNAR <- ampute(complete_blood_controls[,-2], prop=prop, bycases=FALSE, patterns=p[2:159,1:32],
                               freq=freq, mech = "MNAR")
a_controls_2_MNAR <- amputed_controls_2_MNAR$amp
write.csv(a_controls_2_MNAR, '/exports/igmm/eddie/UK-Biobank-proj76173/julia/amputed/a_8.csv')

#MCAR&MAR
# ampute the complete data once for every mechanism
ampdata1 <- ampute(complete_blood_controls[,-2], patterns = p[2:159,1:32], prop = 0.85, mech = "MAR")$amp
ampdata2 <- ampute(complete_blood_controls[,-2], patterns = p[2:159,1:32], prop = 0.85, mech = "MCAR")$amp

# create a random allocation vector
# use the prob argument to specify how much of each mechanism should be created
# here, 0.5 of the missingness should be MAR and 0.5 should be MCAR
indices <- sample(x = c(1, 2), size = nrow(complete_blood_controls[,-2]), 
                  replace = TRUE, prob = c(0.5, 0.5))

# create an empty data matrix
# fill this matrix with values from either of the two amputed datasets
a_controls_2_MCAR <- matrix(NA, nrow = nrow(complete_blood_controls[,-2]), ncol = ncol(complete_blood_controls[,-2]))
a_controls_2_MCAR[indices == 1, ] <- as.matrix(ampdata1[indices == 1, ])
a_controls_2_MCAR[indices == 2, ] <- as.matrix(ampdata2[indices == 2, ])
write.csv(a_controls_2_MCAR, '/exports/igmm/eddie/UK-Biobank-proj76173/julia/amputed/a_9.csv')

################################################################################
prop <- 0.13
p_cases <- md.pattern(blood_controls)
set.seed(1)
a <- sample(14)
freq <- a/sum(a)

#MAR
set.seed(1)
amputed_controls_1_MAR <- ampute(complete_blood_controls, prop=prop, bycases=FALSE, patterns=p_cases[2:15,1:33],
                                 freq=freq)
a_blood_controls_2_MAR <- amputed_controls_1_MAR$amp
write.csv(a_blood_controls_2_MAR, '/exports/igmm/eddie/UK-Biobank-proj76173/julia/amputed/a_10.csv')

#MNAR
set.seed(1)
amputed_controls_1_MNAR <- ampute(complete_blood_controls, prop=prop, bycases=FALSE, patterns=p_cases[2:15,1:33],
                                  freq=freq, mech = "MNAR")
a_blood_controls_2_MNAR <- amputed_controls_1_MNAR$amp
write.csv(a_blood_controls_2_MNAR, '/exports/igmm/eddie/UK-Biobank-proj76173/julia/amputed/a_11.csv')

#MCAR&MAR
# ampute the complete data once for every mechanism
set.seed(1)
ampdata1 <- ampute(complete_blood_controls, patterns = p_cases[2:15,1:33], prop = 0.99, mech = "MAR")$amp
ampdata2 <- ampute(complete_blood_controls, patterns = p_cases[2:15,1:33], prop = 0.99, mech = "MCAR")$amp

# create a random allocation vector
# use the prob argument to specify how much of each mechanism should be created
# here, 0.5 of the missingness should be MAR and 0.5 should be MCAR
indices <- sample(x = c(1, 2), size = nrow(complete_blood_controls), 
                  replace = TRUE, prob = c(0.5, 0.5))

# create an empty data matrix
# fill this matrix with values from either of the two amputed datasets
a_blood_controls_2_MCAR <- matrix(NA, nrow = nrow(complete_blood_controls), ncol = ncol(complete_blood_controls))
a_blood_controls_2_MCAR[indices == 1, ] <- as.matrix(ampdata1[indices == 1, ])
a_blood_controls_2_MCAR[indices == 2, ] <- as.matrix(ampdata2[indices == 2, ])
write.csv(a_blood_controls_2_MCAR, '/exports/igmm/eddie/UK-Biobank-proj76173/julia/amputed/a_12.csv')

################################################################################
prop <- 0.05
p_cases <- md.pattern(blood_controls)
set.seed(1)
a <- sample(14)
freq <- a/sum(a)

#MAR
set.seed(1)
amputed_controls_1_MAR <- ampute(complete_blood_controls, prop=prop, bycases=FALSE, patterns=p_cases[2:15,1:33],
                                 freq=freq)
a_blood_controls_3_MAR <- amputed_controls_1_MAR$amp
write.csv(a_blood_controls_3_MAR, '/exports/igmm/eddie/UK-Biobank-proj76173/julia/amputed/a_13.csv')

#MNAR
set.seed(1)
amputed_controls_1_MNAR <- ampute(complete_blood_controls, prop=prop, bycases=FALSE, patterns=p_cases[2:15,1:33],
                                  freq=freq, mech = "MNAR")
a_blood_controls_3_MNAR <- amputed_controls_1_MNAR$amp
write.csv(a_blood_controls_3_MNAR, '/exports/igmm/eddie/UK-Biobank-proj76173/julia/amputed/a_14.csv')

#MCAR&MAR
# ampute the complete data once for every mechanism
set.seed(1)
ampdata1 <- ampute(complete_blood_controls, patterns = p_cases[2:15,1:33], prop = 0.5, mech = "MAR")$amp
ampdata2 <- ampute(complete_blood_controls, patterns = p_cases[2:15,1:33], prop = 0.5, mech = "MCAR")$amp

# create a random allocation vector
# use the prob argument to specify how much of each mechanism should be created
# here, 0.5 of the missingness should be MAR and 0.5 should be MCAR
indices <- sample(x = c(1, 2), size = nrow(complete_blood_controls), 
                  replace = TRUE, prob = c(0.5, 0.5))

# create an empty data matrix
# fill this matrix with values from either of the two amputed datasets
a_blood_controls_3_MCAR <- matrix(NA, nrow = nrow(complete_blood_controls), ncol = ncol(complete_blood_controls))
a_blood_controls_3_MCAR[indices == 1, ] <- as.matrix(ampdata1[indices == 1, ])
a_blood_controls_3_MCAR[indices == 2, ] <- as.matrix(ampdata2[indices == 2, ])
write.csv(a_blood_controls_3_MCAR, '/exports/igmm/eddie/UK-Biobank-proj76173/julia/amputed/a_15.csv')
################################################################################

sam <- sample(1:nrow(complete_blood_controls), 2*nrow(blood_cases))
control_small <- complete_blood_controls[ sam,]
write.csv(control_small, '/exports/igmm/eddie/UK-Biobank-proj76173/julia/amputed/c_2.csv')

prop <- 0.05
p_cases <- md.pattern(blood_cases)
freq <- as.numeric(row.names(p_cases)[2:5])/sum(as.numeric(row.names(p_cases)[2:5]))

#MAR
set.seed(1)
amputed_controls_1_MAR <- ampute(control_small, prop=prop, bycases=FALSE, patterns=p_cases[2:5,1:33],
                              freq=freq)
a_blood_controls_MAR <- amputed_controls_1_MAR$amp
write.csv(a_blood_controls_MAR, '/exports/igmm/eddie/UK-Biobank-proj76173/julia/amputed/a_16.csv')

#MNAR
set.seed(1)
amputed_controls_1_MNAR <- ampute(control_small, prop=prop, bycases=FALSE, patterns=p_cases[2:5,1:33],
                               freq=freq, mech = "MNAR")
a_blood_controls_MNAR <- amputed_controls_1_MNAR$amp
write.csv(a_blood_controls_MNAR, '/exports/igmm/eddie/UK-Biobank-proj76173/julia/amputed/a_17.csv')

#MCAR&MAR
# ampute the complete data once for every mechanism
set.seed(1)
ampdata1 <- ampute(control_small, patterns = p_cases[2:5,1:33], prop = 0.6, mech = "MAR")$amp
ampdata2 <- ampute(control_small, patterns = p_cases[2:5,1:33], prop = 0.6, mech = "MCAR")$amp

# create a random allocation vector
# use the prob argument to specify how much of each mechanism should be created
# here, 0.5 of the missingness should be MAR and 0.5 should be MCAR
indices <- sample(x = c(1, 2), size = nrow(control_small), 
                  replace = TRUE, prob = c(0.5, 0.5))

# create an empty data matrix
# fill this matrix with values from either of the two amputed datasets
a_blood_controls_MCAR <- matrix(NA, nrow = nrow(control_small), ncol = ncol(control_small))
a_blood_controls_MCAR[indices == 1, ] <- as.matrix(ampdata1[indices == 1, ])
a_blood_controls_MCAR[indices == 2, ] <- as.matrix(ampdata2[indices == 2, ])
write.csv(a_blood_controls_MCAR, '/exports/igmm/eddie/UK-Biobank-proj76173/julia/amputed/a_18.csv')

####################################################################################################

sam <- sample(1:nrow(complete_blood_controls), 4*nrow(blood_cases))

control_small <- complete_blood_controls[ sam,]
write.csv(control_small, '/exports/igmm/eddie/UK-Biobank-proj76173/julia/amputed/c_4.csv')

prop <- 0.05
p_cases <- md.pattern(blood_cases)
freq <- as.numeric(row.names(p_cases)[2:5])/sum(as.numeric(row.names(p_cases)[2:5]))

#MAR
set.seed(1)
amputed_controls_1_MAR <- ampute(control_small, prop=prop, bycases=FALSE, patterns=p_cases[2:5,1:33],
                              freq=freq)
a_blood_controls_MAR <- amputed_controls_1_MAR$amp
write.csv(a_blood_controls_MAR, '/exports/igmm/eddie/UK-Biobank-proj76173/julia/amputed/a_19.csv')

#MNAR
set.seed(1)
amputed_controls_1_MNAR <- ampute(control_small, prop=prop, bycases=FALSE, patterns=p_cases[2:5,1:33],
                               freq=freq, mech = "MNAR")
a_blood_controls_MNAR <- amputed_controls_1_MNAR$amp
write.csv(a_blood_controls_MNAR, '/exports/igmm/eddie/UK-Biobank-proj76173/julia/amputed/a_20.csv')

#MCAR&MAR
# ampute the complete data once for every mechanism
set.seed(1)
ampdata1 <- ampute(control_small, patterns = p_cases[2:5,1:33], prop = 0.6, mech = "MAR")$amp
ampdata2 <- ampute(control_small, patterns = p_cases[2:5,1:33], prop = 0.6, mech = "MCAR")$amp

# create a random allocation vector
# use the prob argument to specify how much of each mechanism should be created
# here, 0.5 of the missingness should be MAR and 0.5 should be MCAR
indices <- sample(x = c(1, 2), size = nrow(control_small), 
                  replace = TRUE, prob = c(0.5, 0.5))

# create an empty data matrix
# fill this matrix with values from either of the two amputed datasets
a_blood_controls_MCAR <- matrix(NA, nrow = nrow(control_small), ncol = ncol(control_small))
a_blood_controls_MCAR[indices == 1, ] <- as.matrix(ampdata1[indices == 1, ])
a_blood_controls_MCAR[indices == 2, ] <- as.matrix(ampdata2[indices == 2, ])
write.csv(a_blood_controls_MCAR, '/exports/igmm/eddie/UK-Biobank-proj76173/julia/amputed/a_21.csv')

#################################################################################################

sam <- sample(1:nrow(complete_blood_controls), 8*nrow(blood_cases))
control_small <- complete_blood_controls[ sam,]
write.csv(control_small, '/exports/igmm/eddie/UK-Biobank-proj76173/julia/amputed/c_8.csv')

prop <- 0.05
p_cases <- md.pattern(blood_cases)
freq <- as.numeric(row.names(p_cases)[2:5])/sum(as.numeric(row.names(p_cases)[2:5]))

#MAR
set.seed(1)
amputed_controls_1_MAR <- ampute(control_small, prop=prop, bycases=FALSE, patterns=p_cases[2:5,1:33],
                              freq=freq)
a_blood_controls_MAR <- amputed_controls_1_MAR$amp
write.csv(a_blood_controls_MAR, '/exports/igmm/eddie/UK-Biobank-proj76173/julia/amputed/a_22.csv')

#MNAR
set.seed(1)
amputed_controls_1_MNAR <- ampute(control_small, prop=prop, bycases=FALSE, patterns=p_cases[2:5,1:33],
                               freq=freq, mech = "MNAR")
a_blood_controls_MNAR <- amputed_controls_1_MNAR$amp
write.csv(a_blood_controls_MNAR, '/exports/igmm/eddie/UK-Biobank-proj76173/julia/amputed/a_23.csv')

#MCAR&MAR
# ampute the complete data once for every mechanism
set.seed(1)
ampdata1 <- ampute(control_small, patterns = p_cases[2:5,1:33], prop = 0.6, mech = "MAR")$amp
ampdata2 <- ampute(control_small, patterns = p_cases[2:5,1:33], prop = 0.6, mech = "MCAR")$amp

# create a random allocation vector
# use the prob argument to specify how much of each mechanism should be created
# here, 0.5 of the missingness should be MAR and 0.5 should be MCAR
indices <- sample(x = c(1, 2), size = nrow(control_small), 
                  replace = TRUE, prob = c(0.5, 0.5))

# create an empty data matrix
# fill this matrix with values from either of the two amputed datasets
a_blood_controls_MCAR <- matrix(NA, nrow = nrow(control_small), ncol = ncol(control_small))
a_blood_controls_MCAR[indices == 1, ] <- as.matrix(ampdata1[indices == 1, ])
a_blood_controls_MCAR[indices == 2, ] <- as.matrix(ampdata2[indices == 2, ])
write.csv(a_blood_controls_MCAR, '/exports/igmm/eddie/UK-Biobank-proj76173/julia/amputed/a_24.csv')

####################################################################################################

