
library(tidyverse)
library(lubridate)
library(glmnet)
library(lmtest)
library(sandwich)
library(RColorBrewer)
library(gbm)

load("data/model_data_completeObs.RData")

## SETTINGS ##
num_folds <- 5
shrinkage_vals <- seq(from = 0.001, by = 0.001, length.out = 20)
depth_vals <- 5:12
num_trees <- 1000

source("lib/utilities.R")

set.seed(20210504)
test_ids <- nrow(model_data) %>% drawFolds(K = num_folds)

########### BEGIN: TUNING OF BOOSTED REGRESSION TREE #############

grid <- shrinkage_vals
result <- tibble()

for (d in depth_vals) {
  cat(paste0("Depth: ",d,"\n"))
  for (g in grid) {
    cat(paste0("  Shrinkage: ",g,"\n"))
    sq.error <- c()
    for (k in 1:5) {
      cat(paste0("    Fold: ",k,"\n"))
      test_ix <- test_ids[[k]]
      boost.uplift <- gbm(log(CC6620) ~ . - Date - Year - CC6630 - CC66200, 
                        data = model_data[-test_ix,], 
                        distribution = "gaussian", n.trees = num_trees, interaction.depth = d, 
                        shrinkage = g)
      y.hat.boost <- predict(boost.uplift, newdata = model_data[test_ix,], n.trees = 5000)
      uplift.test <- model_data$CC6620[test_ix] %>% log
      sq.error <- c(sq.error, (y.hat.boost - uplift.test)^2)
    }
    result_i <- tibble(
      Depth = d,
      Shrinkage = g,
      MSE = mean(sq.error)
    )
    result <- result %>% bind_rows(result_i)
  }
}

########### END: TUNING OF BOOSTED REGRESSION TREE #############

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

