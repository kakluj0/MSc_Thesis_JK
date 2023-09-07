# per variable RMSE/NRMSE/MAE

j=24
#h='MNAR'
h='MCAR+MAR'

dir <- paste('/exports/igmm/eddie/UK-Biobank-proj76173/julia/amputed/a_',j,'.csv',sep='')
data <- data.frame(read.csv(dir))
data <- data[,-1]


cases <- data.frame(read.csv('/exports/igmm/eddie/UK-Biobank-proj76173/sjoerd/metabolites_CFS_gold_cases.csv'))
controls <- data.frame(read.csv('/exports/igmm/eddie/UK-Biobank-proj76173/sjoerd/metabolites_CFS_gold_controls_v2.csv'))


full_blood_count <- c('f.21022.0.0','f.22001.0.0','f.30000.0.0','f.30240.0.0','f.30250.0.0','f.30070.0.0',
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




mean <- data.frame(read.csv(paste('/exports/igmm/eddie/UK-Biobank-proj76173/julia/imputed/mean/a_',j,'.csv',sep='')))
mean <- mean[,-1]
mean <- mean[,-1]
meanreg <- data.frame(read.csv(paste('/exports/igmm/eddie/UK-Biobank-proj76173/julia/imputed/mean_reg/a_',j,'.csv',sep='')))
meanreg <- meanreg[,-1]
s_meanreg <- data.frame(read.csv(paste('/exports/igmm/eddie/UK-Biobank-proj76173/julia/imputed/s_mean_reg/a_',j,'.csv',sep='')))
s_meanreg <- s_meanreg[,-1]
bays_mean <- data.frame(read.csv(paste('/exports/igmm/eddie/UK-Biobank-proj76173/julia/imputed/bays_mean/a_',j,'.csv',sep='')))
bays_mean <- bays_mean[,-1]
bays_mean <- bays_mean[,-1]
forest <- data.frame(read.csv(paste('/exports/igmm/eddie/UK-Biobank-proj76173/julia/imputed/forest/a_',j,'.csv',sep='')))
forest <- forest[,-1]
knn <- data.frame(read.csv(paste('/exports/igmm/eddie/UK-Biobank-proj76173/julia/imputed/knn/knn_',j,'.csv',sep='')))
knn <- knn[,-1]
glrm <- data.frame(read.csv(paste('/exports/igmm/eddie/UK-Biobank-proj76173/julia/imputed/glrm/a',j,'.csv',sep='')))
glrm <- glrm[,-1]



library(Metrics)

cols <- colnames(data)[colSums(is.na(data)) > 0]

# if colnames are V1...
cols <- as.numeric(substr(cols,2,4))


rmse <- matrix(NA,length(cols),9)
rmse <- as.data.frame(rmse)
colnames(rmse) <- c('Mean','Mean Regression','Stochatics Mean Regression','Bayesian Regression','RF','KNN','GLRM','Mis','Var')
for (i in 1:length(cols)){
rmse[i,1] <- rmse(complete_blood_cases[is.na(data[,cols[i]]),cols[i]],mean[is.na(data[,cols[i]]),cols[i]])
rmse[i,2] <- rmse(complete_blood_cases[is.na(data[,cols[i]]),cols[i]],meanreg[is.na(data[,cols[i]]),cols[i]])
rmse[i,3] <- rmse(complete_blood_cases[is.na(data[,cols[i]]),cols[i]],s_meanreg[is.na(data[,cols[i]]),cols[i]])
rmse[i,4] <- rmse(complete_blood_cases[is.na(data[,cols[i]]),cols[i]],bays_mean[is.na(data[,cols[i]]),cols[i]])
rmse[i,5] <- rmse(complete_blood_cases[is.na(data[,cols[i]]),cols[i]],forest[is.na(data[,cols[i]]),cols[i]])
rmse[i,6] <- rmse(complete_blood_cases[is.na(data[,cols[i]]),cols[i]],knn[is.na(data[,cols[i]]),cols[i]])
rmse[i,7] <- rmse(complete_blood_cases[is.na(data[,cols[i]]),cols[i]],glrm[is.na(data[,cols[i]]),cols[i]])
rmse[i,8] <- h
#rmse[i,9] <- substr(cols[i],3,7)
rmse[i,9] <- substr(colnames(complete_blood_cases)[cols[i]],3,7)
}

dir <- paste('/exports/igmm/eddie/UK-Biobank-proj76173/julia/imp_res/rmse/a',j,'.csv',sep='')
write.csv(rmse,dir)

mae <- matrix(NA,length(cols),9)
mae <- as.data.frame(mae)
colnames(mae) <- c('Mean','Mean Regression','Stochatics Mean Regression','Bayesian Regression','RF','KNN','GLRM','Mis','Var')
for (i in 1:length(cols)){
mae[i,1] <- mae(complete_blood_cases[is.na(data[,cols[i]]),cols[i]],mean[is.na(data[,cols[i]]),cols[i]])
mae[i,2] <- mae(complete_blood_cases[is.na(data[,cols[i]]),cols[i]],meanreg[is.na(data[,cols[i]]),cols[i]])
mae[i,3] <- mae(complete_blood_cases[is.na(data[,cols[i]]),cols[i]],s_meanreg[is.na(data[,cols[i]]),cols[i]])
mae[i,4] <- mae(complete_blood_cases[is.na(data[,cols[i]]),cols[i]],bays_mean[is.na(data[,cols[i]]),cols[i]])
mae[i,5] <- mae(complete_blood_cases[is.na(data[,cols[i]]),cols[i]],forest[is.na(data[,cols[i]]),cols[i]])
mae[i,6] <- mae(complete_blood_cases[is.na(data[,cols[i]]),cols[i]],knn[is.na(data[,cols[i]]),cols[i]])
mae[i,7] <- mae(complete_blood_cases[is.na(data[,cols[i]]),cols[i]],glrm[is.na(data[,cols[i]]),cols[i]])
mae[i,8] <- h
#mae[i,9] <- substr(cols[i],3,7)
mae[i,9] <- substr(colnames(complete_blood_cases)[cols[i]],3,7)
}

dir <- paste('/exports/igmm/eddie/UK-Biobank-proj76173/julia/imp_res/mae/a',j,'.csv',sep='')
write.csv(mae,dir)

library(missForest)

nrmse <- matrix(NA,length(cols),9)
nrmse <- as.data.frame(nrmse)
colnames(nrmse) <- c('Mean','Mean Regression','Stochatics Mean Regression','Bayesian Regression','RF','KNN','GLRM','Mis','Var')
for (i in 1:length(cols)){
nrmse[i,1] <- nrmse(mean[,cols[i]], data[,cols[i]], complete_blood_cases[,cols[i]])
nrmse[i,2] <- nrmse(meanreg[,cols[i]], data[,cols[i]], complete_blood_cases[,cols[i]])
nrmse[i,3] <- nrmse(s_meanreg[,cols[i]], data[,cols[i]], complete_blood_cases[,cols[i]])
nrmse[i,4] <- nrmse(bays_mean[,cols[i]], data[,cols[i]], complete_blood_cases[,cols[i]])
nrmse[i,5] <- nrmse(forest[,cols[i]], data[,cols[i]], complete_blood_cases[,cols[i]])
nrmse[i,6] <- nrmse(knn[,cols[i]], data[,cols[i]], complete_blood_cases[,cols[i]])
nrmse[i,7] <- nrmse(glrm[,cols[i]], data[,cols[i]], complete_blood_cases[,cols[i]])
nrmse[i,8] <- h
#nrmse[i,9] <- substr(cols[i],3,7)
nrmse[i,9] <- substr(colnames(complete_blood_cases)[cols[i]],3,7)
}

dir <- paste('/exports/igmm/eddie/UK-Biobank-proj76173/julia/imp_res/nrmse/a',j,'.csv',sep='')
write.csv(nrmse,dir)
