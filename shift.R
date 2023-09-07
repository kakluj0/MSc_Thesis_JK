library(txshift)
library(sl3)

cases <- data.frame(read.csv('/exports/igmm/eddie/UK-Biobank-proj76173/sjoerd/metabolites_CFS_gold_cases.csv'))
controls <- data.frame(read.csv('/exports/igmm/eddie/UK-Biobank-proj76173/sjoerd/metabolites_CFS_gold_controls_v2.csv'))
#setwd("/exports/igmm/datastore/UK-Biobank-proj76173/ava/")
#cases <- read_csv('metabolites_CFS_gold_cases.csv')
#controls <- read_csv('metabolites_CFS_gold_controls_v2.csv')
cases$f.120010.0.0 <- 1
controls$f.120010.0.0 <- 0

############################################################################################################################
#CATE
sex_cases <- which(cases$f.22001.0.0==0)
sex_controls <- which(controls$f.22001.0.0==0)
cases <- cases[sex_cases,]
controls <- controls[sex_controls,]
#############################################################################################################################

full <- rbind(cases,controls)

# permute ME column
set.seed(1)
#index <- sample(1:nrow(full),nrow(full))
#full$f.120010.0.0[index] <- full$f.120010.0.0

blood <- colnames(full)
blood <- blood[-(blood=='f.eid')]
WY <- c("f.21022.0.0", "f.22001.0.0","f.120010.0.0")
blood <- blood[!(blood%in%WY)]

#res <- rep(NA,5*length(blood))
#resse <- rep(NA,5*length(blood))
#name <- rep(NA,5*length(blood))
#delta <- rep(NA,5*length(blood))
#cases <- rep(NA,5*length(blood))
#controls <- rep(NA,5*length(blood))

msm <- as.data.frame(matrix(NA,2*length(blood),7))
colnames(msm) <- c('param','ci_lwr','param_est','ci_upr','param_se','p_value','name')
points <- as.data.frame(matrix(NA,5*length(blood),5))
colnames(points) <- c('delta','ci_lwr','psi','ci_upr','name')
k=1
l=1
for (j in 1:length(blood)){
  var <- blood[j]
  new <- full[,c('f.21022.0.0','f.22001.0.0','f.120010.0.0',var)]
  new <- new[complete.cases(new),]
  cas <- sum(which(new$f.120010.0.0==1))
  contr <- sum(which(new$f.120010.0.0==0))
  new <- as.matrix(new)
  W <- new[,c('f.21022.0.0','f.22001.0.0')]
  Y <- new[,'f.120010.0.0']
  A <- new[,var]

  # permutation
  index <- sample(1:length(A),length(A))
  A[index] <- A

  # grid of shifts
  std <- sd(A)
  delta_grid <- c(-std,-std/2,0,std/2,std)
  # one-step estimation
  #for (i in 1:5){
    #del <- delta_grid[i]
    os <- msm_vimshift(
      W = W, A = A, Y = Y,
      estimator = "onestep",
      g_exp_fit_args = list(
        fit_type = "sl",
        sl_learners_density = Lrnr_density_hse$new(Lrnr_hal9001$new())
      ),
      Q_fit_args = list(fit_type = "glm", glm_formula = "Y ~ ."),
delta_grid=delta_grid
)
    #res[k]<- os$psi
    #resse[k] <- sqrt(os$var)
    #delta[k] <- os$.delta
    #name[k] <- var
    #cases[k] <- cas
    #controls[k] <- contr
    #k <- k+1
  #}

os$msm_est$name <- rep(var,2)
os$param_est$name <- rep(var,5)

msm[k,] <- os$msm_est[1,]
points[l,] <- os$param_est[1,]

k <- k+1
l <- l+1
msm[k,] <- os$msm_est[2,]
points[l,] <- os$param_est[2,]

k <- k+1
l <- l+1
points[l,] <- os$param_est[3,]

l <- l+1
points[l,] <- os$param_est[4,]

l <- l+1
points[l,] <- os$param_est[5,]

l <- l+1

}

#group <- rep('Combined',length(res))
#data <- cbind(res,resse,name,delta,group)

#dir <- '/exports/igmm/eddie/UK-Biobank-proj76173/julia/shift_nobatches/females.csv'
#write.csv(data,dir)

dir <- '/exports/igmm/eddie/UK-Biobank-proj76173/julia/shift_nobatches/p2_msm_females.csv'
write.csv(msm,dir)
dir <- '/exports/igmm/eddie/UK-Biobank-proj76173/julia/shift_nobatches/p2_points_females.csv'
write.csv(points,dir)
