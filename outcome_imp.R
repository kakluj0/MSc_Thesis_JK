library(readr)
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
library(h2o)

# upload data from Eddie
cases <- data.frame(read.csv('/exports/igmm/eddie/UK-Biobank-proj76173/sjoerd/metabolites_CFS_gold_cases.csv'))
controls <- data.frame(read.csv('/exports/igmm/eddie/UK-Biobank-proj76173/sjoerd/metabolites_CFS_gold_controls_v2.csv'))

# upload data from Datastore
#setwd("/exports/igmm/datastore/UK-Biobank-proj76173/ava/")
#cases <- read_csv('metabolites_CFS_gold_cases.csv')
#controls <- read_csv('metabolites_CFS_gold_controls_v2.csv')

# sample 2086 controls for the 2nd run
#set.seed(10)
#sam <- sample(1:nrow(controls),nrow(cases))
#controls <- controls[sam,]

# join controls + cases
full <- rbind(cases,controls)
full <- full[,colnames(full)!='f.eid']

#################################################
# Random Forest imputation
non_miss <- rowSums(!is.na(full))
imp <- missRanger(full, num.trees = 50, pmm.k = 3, seed = 5, verbose = 0, case.weights = non_miss,splitrule = "extratrees",sample.fraction=0.1)
write.csv(imp, '/exports/igmm/eddie/UK-Biobank-proj76173/julia/MEimputationforestsmall.csv')
#################################################

#################################################
# K-NN imputation
set.seed(1)
knn.model <- preProcess(full, "knnImpute")
knn_12 <- predict(knn.model,full,na.action = na.pass)
procNames <- data.frame(col = names(knn.model$mean), mean = knn.model$mean, sd = knn.model$std)
for(i in procNames$col){
  knn_12[i] <- knn_12[i]*knn.model$std[i]+knn.model$mean[i]
}

write.csv(knn_12, '/exports/igmm/eddie/UK-Biobank-proj76173/julia/MEimputationknn.csv')
################################################

################################################
# imputation with mice
# mean: method="mean"
# regression: method="norm.predict"
# stochastic regression: method="norm.nob"
# bayesian regression: method="norm"

m <- mice(full, method="norm.nob", m=1, maxit=1)
m <- complete(m)
write.csv(m, '/exports/igmm/eddie/UK-Biobank-proj76173/julia/MEimputationmice.csv')
###############################################

###############################################
# Imputaion with GLRM

impute_df <- full

# Convert binary variables to logical
impute_df[['f.22001.0.0']] = as.logical(impute_df[['f.22001.0.0']])
impute_df[['f.120010.0.0']] = as.logical(impute_df[['f.120010.0.0']])

table(sapply(impute_df, class))

# Create a dataframe describing the loss function by variable; the first variable must have index = 0

losses = data.frame("index" = seq(ncol(impute_df)) - 1,
                    "feature" = colnames(impute_df),
                    "class" = sapply(impute_df, class),
                    "type" = sapply(impute_df, typeof),
                    stringsAsFactors = FALSE)

# loss function for numeric/logical class
losses$loss[losses$class == "numeric"] = "Huber"
losses$loss[losses$class == "logical"] = "Hinge"

# Initialize h2o
h2o::h2o.no_progress()  # Turn off progress bars
analyst_name = "jk"
h2o::h2o.init(max_mem_size = "20g",
              name = paste0("h2o-", analyst_name),
              # Default port is 54321, but other analysts may be using that.
              port = 54320,
              # This can reduce accidental sharing of h2o processes on a shared server.
              username = analyst_name,
              password = paste0("pw-", analyst_name),
              # Use half of available cores for h2o.
              #nthreads = get_cores()
              )

# Convert data to h2o object
h2o_df = h2o::as.h2o(impute_df)

# Split data into train & validation
split = h2o::h2o.splitFrame(h2o_df, ratios = 0.75, seed = 1)
train = split[[1]]
valid = split[[2]]

# h2o required validation frame to have the same number of rows as train data
val_df = as.data.frame(valid)
tra_df <- as.data.frame(train)
tra_df <- tra_df[1:nrow(val_df),]
train <- h2o::as.h2o(tra_df)
valid <- h2o::as.h2o(val_df)

# hyperparameter search grid
params = expand.grid(
  k = 24,
  regularization_x = c("None", "Quadratic", "L1"),
  regularization_y = c("None", "Quadratic",'L1'),
  gamma_x = c(0, 1, 4),
  gamma_y = c(0, 1, 4),
  error_num = NA,
  error_cat = NA,
  objective = NA,
  stringsAsFactors = FALSE)

# Remove combinations
params = subset(params, regularization_x != "None" | gamma_x == 0)
params = subset(params, regularization_x == "None" | gamma_x != 0)
params = subset(params, regularization_y != "None" | gamma_y == 0)
params = subset(params, regularization_y == "None" | gamma_y != 0)

# Randomly order the params
set.seed(1)
params = params[sample(nrow(params)), ]

glrm_metrics = list()
glrm_sum = list()

# Grid search
for (i in seq_len(nrow(params))) {
    # Create model
    glrm_model = h2o::h2o.glrm(
      training_frame = train,
      validation_frame = valid,
      k = params$k[i],
      loss = "Quadratic",
      regularization_x = params$regularization_x[i],
      regularization_y = params$regularization_y[i],
      gamma_x = params$gamma_x[i],
      gamma_y = params$gamma_y[i],
      transform = "STANDARDIZE",
      max_iterations = 2000,
      max_runtime_secs = 1000,
      seed = 1,
      loss_by_col_idx = losses$index,
      loss_by_col = losses$loss)

 summ_text = capture.output({ h2o::summary(glrm_model) })
    glrm_sum[[i]] = summ_text

    params$objective[i] = glrm_model@model$objective

    try({
      validate = h2o::h2o.performance(glrm_model, valid)
      glrm_metrics[[i]] = validate@metrics

      params$error_num[i] = validate@metrics$numerr
      params$error_cat[i] = validate@metrics$caterr

    })

    # Remove model
    h2o::h2o.rm(glrm_model);rm(glrm_model);gc()
  }

params$error = params$error_num + params$error_cat

params <- as.data.frame(arrange(params,error))

best_params <- head(as.data.frame(arrange(params,error)),1)


# Run on full data
  glrm_result =
    h2o::h2o.glrm(training_frame = h2o_df, cols = colnames(h2o_df),
                  loss = "Quadratic",
                  model_id = "impute_glrm",
                  seed = 1,
                  k = best_params$k,
                  max_iterations = 2000,
                  transform = "standardize",
                  regularization_x = best_params$regularization_x,
                  regularization_y = best_params$regularization_y,
                  gamma_x = best_params$gamma_x,
                  gamma_y = best_params$gamma_y,
                  loss_by_col_idx = losses$index,
                  loss_by_col = losses$loss)

# Extract compressed dataset.
new_data = as.data.frame(h2o::h2o.getFrame(glrm_result@model$representation_name))

# Reconstructed data from GLRM.
recon_df = h2o::h2o.reconstruct(glrm_result, h2o_df,
                                reverse_transform = TRUE)
# Column names
names(recon_df) = names(impute_df)

# Data frame
recon_df = as.data.frame(recon_df)

write.csv(recon_df,'/exports/igmm/eddie/UK-Biobank-proj76173/julia/MEimputeGLRM_full.csv')

