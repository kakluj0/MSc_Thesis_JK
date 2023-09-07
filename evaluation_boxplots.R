# per variable densities/boxplots

library(ggplot2)

# dataset index
j=1

# load data
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

# load amputed data
dir <- paste('/exports/igmm/eddie/UK-Biobank-proj76173/julia/amputed/a_',j,'.csv',sep='')
data <- data.frame(read.csv(dir))
data <- data[,-1]

# load imputed data
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

# columns with missingness
cols <- colnames(data)[colSums(is.na(data)) > 0]

# if colnames are V1...
#cols <- as.numeric(substr(cols,2,4))

# first iteration
i=1
pdf(file = paste('/exports/igmm/eddie/UK-Biobank-proj76173/julia/dens/boxplot',j,'.pdf',sep=''),width = 8,height = 8)
par(mfrow = c(4, 4)) 

# variable
name <- cols[i]

# load data by variable
com <- complete_blood_cases[,name]
mean_com <- mean[,name]
meanreg_com <- meanreg[,name]
meanreg_com <- meanreg[,name]
s_meanreg_com <- s_meanreg[,name]
bayes_com <- bays_mean[,name]
forest_com <- forest[,name]
knn_com <- knn[,name]
glrm_com <- glrm[,name]
data_com <- data[,name]

# create dataframe
com7 <- cbind(com, mean_com,meanreg_com,s_meanreg_com,bayes_com,forest_com,knn_com,glrm_com,data_com)
com7 <- com7[which(!complete.cases(com7)),]
com7 <- as.data.frame(com7)
#com7$com <- as.numeric(com7$com)
#com7$mean_com <- as.numeric(com7$mean_com)
#com7$meanreg_com <- as.numeric(com7$meanreg_com)
#com7$s_meanreg_com <- as.numeric(com7$s_meanreg_com)
#com7$bayes_com <- as.numeric(com7$bayes_com)
#com7$forest_com <- as.numeric(com7$forest_com)
#com7$knn_com <- as.numeric(com7$knn_com)
#com7$glrm_com <- as.numeric(com7$glrm_com)


# check how much missigness is in the column
percentagena <- round(sum(is.na(data_com))/length(data_com)*100,1)

# id of the variable
#name <- substr(colnames(complete_blood_cases)[name],3,7)
name <- substr(name,3,7)

# plot
#p <- ggplot(com7, aes(x=com)) +
#  geom_density()+
#  geom_density(aes(x=mean_com,col='Mean'))+
#  geom_density(aes(x=meanreg_com,col='Mean Regression'))+
#  geom_density(aes(x=s_meanreg_com,col='Stochastics MR'))+
#  geom_density(aes(x=bayes_com,col='Bayes Regression'))+
#  geom_density(aes(x=forest_com,col='Random Forest'))+
#  geom_density(aes(x=knn_com,col='K-NN'))+
#  geom_density(aes(x=glrm_com,col='GLRM'))+
#  labs(title=paste(name,':',percentagena,'%'))+
#  #coord_cartesian(xlim = c(-0.25, 0.25))+
#  #facet_grid(~Mis)+
#  theme_minimal()

# diffrent plot
q1 <- summary(com7$com)[2]
q3 <- summary(com7$com)[5]
iqr <- IQR(com7$com)

#plot(density(com7$com,bw = bw.bcv(com7$com)),xlim=c(q1-3*iqr,q3+3*iqr),ylim=c(0,2*(max(density(com7$com)$y))),main=paste(name,':',percentagena,'% NA'),ylab='',xlab='')
#lines(density(com7$mean_com),col=2)
#lines(density(com7$meanreg_com), col = 3, lwd = 1)
#lines(density(com7$s_meanreg_com), col = 4, lwd = 1)
#lines(density(com7$bayes_com), col = 5, lwd = 1)
#lines(density(com7$forest_com), col = 6, lwd = 1)
#lines(density(com7$knn_com), col = 7, lwd = 1)
#lines(density(com7$glrm_com), col = 8, lwd = 1)

com7 <- com7[,-9]
boxplot(com7,xaxt = "n",outline=FALSE, col = c('red','black','blue','green','purple','grey','brown','orange'),main=paste(name,':',percentagena,'% NA'))

# loop over next iteration
for (i in 2:length(cols)){

# name of the variable
name <- cols[i]

# load data by variable
com <- complete_blood_cases[,name]
mean_com <- mean[,name]
meanreg_com <- meanreg[,name]
meanreg_com <- meanreg[,name]
s_meanreg_com <- s_meanreg[,name]
bayes_com <- bays_mean[,name]
forest_com <- forest[,name]
knn_com <- knn[,name]
glrm_com <- glrm[,name]
data_com <- data[,name]

# create dataframe
com7 <- cbind(com, mean_com,meanreg_com,s_meanreg_com,bayes_com,forest_com,knn_com,glrm_com,data_com)
com7 <- com7[which(!complete.cases(com7)),]
com7 <- as.data.frame(com7)

# check how much missigness is in the column
percentagena <- round(sum(is.na(data_com))/length(data_com)*100,1)

# id of the variable
#name <- substr(colnames(complete_blood_cases)[name],3,7)
name <- substr(name,3,7)

# plot
#q <- ggplot(com7, aes(x=com)) +
#  geom_density()+
#  geom_density(aes(x=mean_com,col='Mean'))+
#  geom_density(aes(x=meanreg_com,col='Mean Regression'))+
#  geom_density(aes(x=s_meanreg_com,col='Stochastics MR'))+
#  geom_density(aes(x=bayes_com,col='Bayes Regression'))+
#  geom_density(aes(x=forest_com,col='Random Forest'))+
#  geom_density(aes(x=knn_com,col='K-NN'))+
#  geom_density(aes(x=glrm_com,col='GLRM'))+
#  labs(title=paste(name,':',percentagena,'%'))+
#  #coord_cartesian(xlim = c(-0.25, 0.25))+
#  #facet_grid(~Mis)+
#  theme_minimal()
#
#p <- p+q
#

# diffrent plot
#q1 <- summary(com7$com)[2]
#q3 <- summary(com7$com)[5]
#iqr <- IQR(com7$com)

#plot(density(com7$com,bw = bw.bcv(com7$com)),xlim=c(q1-3*iqr,q3+3*iqr),ylim=c(0,2*(max(density(com7$com)$y))),main=paste(name,':',percentagena,'% NA'),ylab='',xlab='')
#lines(density(com7$mean_com),col=2)
#lines(density(com7$meanreg_com), col = 3, lwd = 1)
#lines(density(com7$s_meanreg_com), col = 4, lwd = 1)
#lines(density(com7$bayes_com), col = 5, lwd = 1)
#lines(density(com7$forest_com), col = 6, lwd = 1)
#lines(density(com7$knn_com), col = 7, lwd = 1)
#lines(density(com7$glrm_com), col = 8, lwd = 1)

com7 <- com7[,-9]
boxplot(com7,xaxt = "n",outline=FALSE, col = c('red','black','blue','green','purple','grey','brown','orange'),main=paste(name,':',percentagena,'% NA'))

}

dev.off()
#dir <- paste('/exports/igmm/eddie/UK-Biobank-proj76173/julia/dens/boxplot_',j,'.pdf',sep='')
#write.csv(rmse,dir)
#ggsave(p,dir)
