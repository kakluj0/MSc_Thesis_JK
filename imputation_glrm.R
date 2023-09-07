# GLRM
library(dplyr)

# load data
j = 11
dir <- paste('/exports/igmm/eddie/UK-Biobank-proj76173/julia/amputed/a_',j,'.csv',sep='')
data <- data.frame(read.csv(dir))
data <- data[,-1]
impute_df <- data

# change data type for binary variables
impute_df[['f.22001.0.0']] = as.logical(impute_df[['f.22001.0.0']])
impute_df[['f.120010.0.0']] = as.logical(impute_df[['f.120010.0.0']])

# loss function by variable
losses = data.frame("index" = seq(ncol(impute_df)) - 1,
                    "feature" = colnames(impute_df),
                    "class" = sapply(impute_df, class),
                    "type" = sapply(impute_df, typeof),
                    stringsAsFactors = FALSE)

losses$loss[losses$class == "numeric"] = "Huber"
losses$loss[losses$class == "logical"] = "Hinge"
losses$loss[losses$class == "integer"] = "Huber"

# Initialize h2o
h2o::h2o.no_progress()  # Turn off progress bars
analyst_name = "jk"
h2o::h2o.init(max_mem_size = "12g",
              name = paste0("h2o-", analyst_name),
              username = analyst_name,
              password = paste0("pw-", analyst_name),
              )

# data to h2o object
h2o_df = h2o::as.h2o(impute_df)

# train and validation split
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
  # Try 3 values on the exponential scale up to the maximum number of predictors.
  k = 24,
  regularization_x = c("None", "Quadratic", "L1"),
  regularization_y = c("None", "Quadratic", 'L1'),
  gamma_x = c(0, 1, 4),
  gamma_y = c(0, 1, 4),
  error_num = NA,
  error_cat = NA,
  objective = NA,
  stringsAsFactors = FALSE)

# remove unnecessary combinations
params = subset(params, regularization_x != "None" | gamma_x == 0)
params = subset(params, regularization_x == "None" | gamma_x != 0)
params = subset(params, regularization_y != "None" | gamma_y == 0)
params = subset(params, regularization_y == "None" | gamma_y != 0)

set.seed(1)
params = params[sample(nrow(params)), ]

glrm_metrics = list()
glrm_sum = list()

# grid search
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

    # Predict on validation set and extract error
    try({
      validate = h2o::h2o.performance(glrm_model, valid)
      #print(validate@metrics)
      glrm_metrics[[i]] = validate@metrics

      params$error_num[i] = validate@metrics$numerr
      params$error_cat[i] = validate@metrics$caterr

    })

    # Removing the model
    h2o::h2o.rm(glrm_model)
  }


params$error = params$error_num + params$error_cat
params <- as.data.frame(arrange(params,error))
best_params <- head(as.data.frame(arrange(params,error)),1)

# parameters from grid search
dir <- paste('/exports/igmm/eddie/UK-Biobank-proj76173/julia/imputed/glrm/params/a',j,'.csv',sep='')
write.csv(best_params,dir)

best_params <- as.data.frame(read.csv('/exports/igmm/eddie/UK-Biobank-proj76173/julia/imputed/glrm/params/a11.csv')[,-1])

# run on entire data
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

# complessed dataset
new_data = as.data.frame(h2o::h2o.getFrame(glrm_result@model$representation_name))

# reconstructed data
recon_df = h2o::h2o.reconstruct(glrm_result, h2o_df,
                                reverse_transform = TRUE)
# column names
names(recon_df) = names(impute_df)

# Convert from h2o object back to an R df.
recon_df = as.data.frame(recon_df)

# Impute data
impute <- impute_df
for (i in 1:ncol(impute)){
  impute[which(is.na(impute[,i])),i] <- recon_df[which(is.na(impute[,i])),i]
}

dir <- paste('/exports/igmm/eddie/UK-Biobank-proj76173/julia/imputed/glrm/a',j,'.csv',sep='')
write.csv(impute,dir)

#Shut down h2o
h2o::h2o.shutdown(prompt = FALSE)
