# Overall RMSE, NRMSE, MAE
j=9

dir <- paste('/exports/igmm/eddie/UK-Biobank-proj76173/julia/amputed/a_',j,'.csv',sep='')
data <- data.frame(read.csv(dir))
data <- data[,-1]

cases <- data.frame(read.csv('/exports/igmm/eddie/UK-Biobank-proj76173/sjoerd/metabolites_CFS_gold_cases.csv'))
controls <- data.frame(read.csv('/exports/igmm/eddie/UK-Biobank-proj76173/sjoerd/metabolites_CFS_gold_controls_v2.csv'))

full_blood_count <- c('f.21022.0.0','f.22001.0.0','f.30000.0.0',                                        'f.30240.0.0','f.30250.0.0','f.30070.0.0',
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

complete_blood_cases <- as.data.frame(read.csv('/exports/igmm/eddie/UK-Biobank-proj76173/julia/amputed/c_1.csv'))
complete_blood_cases <- complete_blood_cases[,-1]
complete_blood_cases <- complete_blood_controls[,-2]

sqrt(sum((mean[is.na(data)]-complete_blood_cases[is.na(data)])^2)/sum(is.na(data)))
sqrt(sum((meanreg[is.na(data)]-complete_blood_cases[is.na(data)])^2)/sum(is.na(data)))
sqrt(sum((s_meanreg[is.na(data)]-complete_blood_cases[is.na(data)])^2)/sum(is.na(data)))
sqrt(sum((bays_mean[is.na(data)]-complete_blood_cases[is.na(data)])^2)/sum(is.na(data)))
sqrt(sum((forest[is.na(data)]-complete_blood_cases[is.na(data)])^2)/sum(is.na(data)))
sqrt(sum((knn[is.na(data)]-complete_blood_cases[is.na(data)])^2)/sum(is.na(data)))
sqrt(sum((glrm[is.na(data)]-complete_blood_cases[is.na(data)])^2)/sum(is.na(data)))

library(Metrics)
mae(complete_blood_cases[is.na(data)],mean[is.na(data)])
mae(complete_blood_cases[is.na(data)],meanreg[is.na(data)])
mae(complete_blood_cases[is.na(data)],s_meanreg[is.na(data)])
mae(complete_blood_cases[is.na(data)],bays_mean[is.na(data)])
mae(complete_blood_cases[is.na(data)],forest[is.na(data)])
mae(complete_blood_cases[is.na(data)],knn[is.na(data)])
mae(complete_blood_cases[is.na(data)],glrm[is.na(data)])

library(missForest)
nrmse(mean, data, complete_blood_cases)
nrmse(meanreg, data, complete_blood_cases)
nrmse(s_meanreg, data, complete_blood_cases)
nrmse(bays_mean, data, complete_blood_cases)
nrmse(forest, data, complete_blood_cases)
nrmse(knn, data, complete_blood_cases)
nrmse(glrm, data, complete_blood_cases)
