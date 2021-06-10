
library(gbm)
library(glmnet)
library(lmtest)
library(lubridate)
library(rjson)
library(sandwich)
library(tidyverse)

# Load settings, data, and utility functions
settings <- fromJSON(file = "config.json")
load("data/model_data_completeObs.RData")
source("lib/utilities.R")

# Create settings variables
set.seed(settings$random_seed)
num_folds <- settings$num_folds
shrinkage_vals <- seq(from = settings$boost_tuning$shrinkage_vals$from, 
                      by = settings$boost_tuning$shrinkage_vals$by, 
                      length.out = settings$boost_tuning$shrinkage_vals$length.out)
depth_vals <- (settings$boost_tuning$depth_vals$min):(settings$boost_tuning$depth_vals$max)
num_trees <- settings$boost_tuning$num_trees

########### BEGIN: TUNING OF BOOSTED REGRESSION TREE #############
if(settings$boost_tuning$run_boost_tuning) {
  tuning_result <- model_data %>% tuneBoostedTree(shrinkage_vals, depth_vals, num_trees)
  save(tuning_result, file = "output/tuning_result.RData")
}

########### END: TUNING OF BOOSTED REGRESSION TREE #############

tuning_result %>% ggplot(aes(x = Shrinkage, y = MSE, color = as.factor(Depth))) +
  geom_line()


min(tuning_result$MSE)
########### BEGIN: BACKWARD STEPWISE BOOSTED ALGORITHM ###########

min.depth <- 9
min.shrinkage <- 0.007
covariates <- colnames(model_data)
covariates <- covariates[!(covariates %in% c("Date","CC6630","CC6620","CC66200","Year"))]
dep_var <- "log(CC6620)"
max_iter <- length(covariates)
threshold <- 0.0005

buildFormula <- function(dep_var, covariates) {
  f <- paste0(covariates, collapse = " + ")
  f <- paste0(dep_var, " ~ ", f)
  return(as.formula(f))
}

# Function to fit cross-fit GBM
cf_gbm <- function(model_f, data, folds, min.depth, min.shrinkage) {
  
  K <- length(folds)
  y.hat.full <- c()

  for (k in 1:K) {
    cat(paste0("  Fold ",k,"\n"))
    test_ix <- folds[[k]]
    boost.uplift <- gbm(model_f, 
                      data = data[-test_ix,], 
                      distribution = "gaussian", n.trees = 5000, 
                      interaction.depth = min.depth, 
                      shrinkage = min.shrinkage)
    y.hat.boost <- predict(boost.uplift, newdata = data[test_ix,], n.trees = 5000)
    y.hat.full <- c(y.hat.full, y.hat.boost)
  }
  
  y.hat.ordered <- rep(NA, nrow(data))
  y.hat.ordered[unlist(folds)] <- y.hat.full
  return(y.hat.ordered)
}

# Fit full model
cat(paste0("Running:\n"))
model_f <- buildFormula(dep_var, covariates)
y_hat <- cf_gbm(model_f, model_data, test_ids, min.depth, min.shrinkage)
mse_min <- mean((y_hat - log(model_data$CC6620))^2)
  
result_table <- tibble(
  Vars.Dropped = 0,
  Var.Name = "Base",
  MSE = mse_min
)

# Start backward stepwise selection loop
iter <- 1
while(iter <= max_iter) {
  mse_c <- c()
  for (c in covariates) {
    cat(paste0("Drop var: ",c,"\n"))
    covariates_i <- covariates[covariates != c]
    model_f <- buildFormula(dep_var, covariates_i)
    y_hat <- cf_gbm(model_f, model_data, test_ids, min.depth, min.shrinkage)
    mse_c <- c(mse_c, mean((y_hat - log(model_data$CC6620))^2))
  }
  
  if (min(mse_c - mse_min) > threshold) { break() }
  
  drop_var <- covariates[which.min(mse_c)]
  covariates <- covariates[covariates != drop_var]
  
  result_table_i <- tibble(
    Vars.Dropped = iter,
    Var.Name = drop_var,
    MSE = min(mse_c)
  )
  
  result_table <- result_table %>% bind_rows(result_table_i)
  mse_min <- min(mse_c)
  iter <- iter + 1
}

################ END: BACKWARD STEPWISE BOOSTED ALGORITHM ###########

########### BEGIN: FIT CROSS-FIT REDUCED MODEL ###########

# Fit reduced model
cat(paste0("Running:\n"))
model_f <- buildFormula(dep_var, covariates)
y_hat <- cf_gbm(model_f, model_data, test_ids, min.depth, min.shrinkage)

tibble(
  Type = c(rep("Predicted", 727), rep("Actual", 727)) %>% as.factor,
  Index = c(1:727, 1:727),
  Value = c(y_hat, log(model_data$CC6620))
  ) %>%
    ggplot(aes(x = Index, y = Value, color = Type)) +
      geom_line()

save(y_hat, file = "boosting_results/yhat_reduced_CF_boost.RData")

########### END: FIT CROSS-FIT REDUCED MODEL ###########

########### BEGIN: FIT FULL MODEL ###########
## Fit full dataset, to view variable importance
boost.uplift.full <- gbm(model_f, 
                      data = model_data, 
                      distribution = "gaussian", n.trees = 5000, 
                      interaction.depth = min.depth, 
                      shrinkage = min.shrinkage)
save(boost.uplift.full, file="boosting_results/boost_uplift_noCF.RData")

plot(boost.uplift.full, i="Load.Mileage.MW")

########### END: FIT FULL MODEL ###########

