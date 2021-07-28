
#####################################################################################
#                                                                                   #
#  Triolo, R.C. and Wolak, F.A. (Forthcoming)                                       #
#  Market Conditions and Uplift Costs in the California Electricity Market          #
#                                                                                   #
#  Open source code for reproduction of results                                     #
#                                                                                   #
#####################################################################################

# Load necessary packages
suppressPackageStartupMessages(library(gbm))
suppressPackageStartupMessages(library(glmnet))
suppressPackageStartupMessages(library(lmtest))
suppressPackageStartupMessages(library(lubridate))
suppressPackageStartupMessages(library(rjson))
suppressPackageStartupMessages(library(sandwich))
suppressPackageStartupMessages(library(tidyverse))

# Set local path, if provided
args <- commandArgs(trailingOnly = TRUE)
local_path <- ifelse(length(args) > 0, args[1], "")

# Load settings, data, and utility functions
config_f <- paste0(local_path, "config.json")
settings <- fromJSON(file = config_f)
model_data_f <- paste0(local_path, "data/model_data_completeObs.RData")
load(model_data_f)
utilities_f <- paste0(local_path, "lib/utilities.R")
source(utilities_f)
messages_f <- paste0(local_path, "lib/messages.R")
source(messages_f)

# Print begin message
beginMessage()

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

# Print status messages
dependentVariableMessage(dep_var)
printProcessStatus(settings)

# Cap minimum values at 1 (log values at 0)
if (min(model_data[,dep_var]) < 1) {
  printZeroValueAdjustment(dep_var)
  model_data[,dep_var] <- model_data[,dep_var] + 1
}

# Draw folds for cross-fitting
test_folds <- drawFolds(n, K = num_folds)

########### TUNING OF BOOSTED REGRESSION TREE #############
if (settings$boost_tuning$run_boost_tuning) {
  tuning_result <- model_data %>% 
    tuneBoostedTree(dep_var, covariates, shrinkage_vals, depth_vals, num_trees, test_folds)
  result_f <- paste0(local_path, "output/tuning_result_",dep_var,".RData")
  save(tuning_result, settings, test_folds, file = result_f)
}

########### BACKWARD STEPWISE BOOSTED REGRESSION TREE ###########
if (settings$boost_backward_selection$run_backward_stepwise) {
  selection_result <- model_data %>% 
    backwardStepwiseSelection(dep_var, covariates, test_folds, min_depth, min_shrinkage, selection_threshold,
                              num_trees)
  result_f <- paste0(local_path, "output/selection_result_",dep_var,".RData")
  save(selection_result, settings, file = result_f)
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
