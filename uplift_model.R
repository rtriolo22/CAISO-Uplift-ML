
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
n <- nrow(model_data)
num_folds <- settings$num_folds
shrinkage_vals <- seq(from = settings$boost_tuning$shrinkage_vals$from, 
                      by = settings$boost_tuning$shrinkage_vals$by, 
                      length.out = settings$boost_tuning$shrinkage_vals$length.out)
depth_vals <- (settings$boost_tuning$depth_vals$min):(settings$boost_tuning$depth_vals$max)
num_trees <- settings$num_trees
dep_var <- settings$y_var
covariates <- colnames(model_data)
covariates <- covariates[!(covariates %in% settings$exclude_variables)]
min_depth <- settings$boost_min_depth
min_shrinkage <- settings$boost_min_shrinkage
selection_threshold <- settings$boost_backward_selection$exit_threshold
cf_market <- settings$boost_reduced_cf$market

# Cap minimum values at 1 (log values at 0)
model_data[,dep_var] <- model_data[,dep_var] %>% unlist %>% map_dbl(function(x){ifelse(x < 1, 1, x)})

# Draw folds for cross-fitting
test_folds <- drawFolds(n, K = num_folds)

########### TUNING OF BOOSTED REGRESSION TREE #############
if (settings$boost_tuning$run_boost_tuning) {
  tuning_result <- model_data %>% 
    tuneBoostedTree(dep_var, covariates, shrinkage_vals, depth_vals, num_trees, test_folds)
  save(tuning_result, settings, file = "output/tuning_result.RData")
}

########### BACKWARD STEPWISE BOOSTED REGRESSION TREE ###########
if (settings$boost_backward_selection$run_backward_stepwise) {
  selection_result <- model_data %>% 
    backwardStepwiseSelection(dep_var, covariates, test_folds, min_depth, min_shrinkage, selection_threshold,
                              num_trees)
  save(selection_result, settings, file = "output/selection_result.RData")
}

########### FIT CROSS-FIT REDUCED MODEL ###########
if (settings$boost_reduced_cf$run_reduced_cf) {
  y_hat_cf <- model_data %>%
    crossFitReducedGBM(dep_var, covariates, test_folds, min_depth, min_shrinkage, num_trees, market)
}
# 
# # Fit reduced model
# cat(paste0("Running:\n"))
# model_f <- buildFormula(dep_var, covariates)
# y_hat <- cf_gbm(model_f, model_data, test_ids, min.depth, min.shrinkage)
# 
# tibble(
#   Type = c(rep("Predicted", 727), rep("Actual", 727)) %>% as.factor,
#   Index = c(1:727, 1:727),
#   Value = c(y_hat, log(model_data$CC6620))
#   ) %>%
#     ggplot(aes(x = Index, y = Value, color = Type)) +
#       geom_line()
# 
# save(y_hat, file = "boosting_results/yhat_reduced_CF_boost.RData")
# 
# ########### END: FIT CROSS-FIT REDUCED MODEL ###########
# 
# ########### BEGIN: FIT FULL MODEL ###########
# ## Fit full dataset, to view variable importance
# boost.uplift.full <- gbm(model_f, 
#                       data = model_data, 
#                       distribution = "gaussian", n.trees = 5000, 
#                       interaction.depth = min.depth, 
#                       shrinkage = min.shrinkage)
# save(boost.uplift.full, file="boosting_results/boost_uplift_noCF.RData")
# 
# plot(boost.uplift.full, i="Load.Mileage.MW")
# 
# ########### END: FIT FULL MODEL ###########
# 
