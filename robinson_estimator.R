
setwd("~/Documents/SU/Research/uplift_study/data/CAISO")

load("cleaned_data/model_data_completeObs.RData")
#load("cleaned_data/model_data_allObs.RData")
load("boosting_results/backward_selection.RData")

library(dbarts)
library(tidyverse)
library(gbm)
library(glmnet)
library(lmtest)
library(sandwich)

set.seed(20210504)
fold_counts <- c(146,146,145,145,145)
sample_ix <- sample(1:nrow(model_data))
test <- list()
for (i in 1:5) {
  test[[i]] <- sample_ix[(sum(fold_counts[0:(i-1)])+1):sum(fold_counts[0:i])]
}

# Tibble for results
result <- tibble()

# Specify focal covariate and desired units
focal_covariate <- "Load.SCE"
dep_var <- "log(CC6620)"
# cov_reduce_mag <- c("Load.Mileage.MW", "DA.LoadError.Over", "Convergence.Max.Supply")
# Reduce order of magnitude
# if (focal_covariate %in% cov_reduce_mag) {
#   model_data[,focal_covariate] <- model_data[,focal_covariate] / 1000
# }
# Normalize
model_data[,focal_covariate] <- (model_data[,focal_covariate] - mean(unlist(model_data[,focal_covariate]))) / sd(unlist(model_data[,focal_covariate]))

# Function to build model formulas
buildFormula <- function(dep_var, covariates, poly.2 = FALSE, cov_fct = c()) {
  if (poly.2) {
    covariates_cont <- covariates[!(covariates %in% cov_fct)]
    covariates_fct <- covariates[(covariates %in% cov_fct)]
    f <- paste0(covariates_cont, collapse = ", degree = 2, raw = TRUE) + poly(")
    f <- paste0("poly(",f,", degree = 2, raw = TRUE)")
    f <- paste0(f, " + ", paste0(covariates_fct, collapse = " + "))
  } else {
    f <- paste0(covariates, collapse = " + ")
  }
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

# Fit full OLS
ols_covariates <- colnames(model_data)
ols_covariates <- ols_covariates[!(ols_covariates %in% c("Date","CC6630","CC6620","CC66200","Year"))]
model_f_ols <- buildFormula(dep_var, ols_covariates)
model_ols <- lm(model_f_ols, data = model_data)
coef_ols <- coeftest(model_ols, vcov=vcovHC(model_ols, type='HC2')) #[2,]
coef_ols <- coef_ols[focal_covariate,] 
result <- result %>%
  bind_rows(tibble(Coefficient = focal_covariate, 
                   Method = "OLS",
                   Estimate = coef_ols[1],
                   `Std. Error` = coef_ols[2],
                   `t value` = coef_ols[3],
                   `Pr(>|t|)` = coef_ols[4],
                   `MSE E[Y|Z]` = NA,
                   `MSE E[X|Z]` = NA))

# Fit boosted regression tree Robinson estimator
min.depth <- 9
min.shrinkage <- 0.007
covariates <- colnames(model_data)
covariates <- covariates[!(covariates %in% c("Date","CC6630","CC6620","CC66200","Year"))]
# Remove variables omitted from backward stepwise selection
covariates <- covariates[!(covariates %in% result_table$Var.Name[2:9])]

X.cov <- focal_covariate
Y.cov <- dep_var
Z.cov <- covariates[covariates != X.cov]

# Find E(Y|Z) with a cross-fit boosted tree
model_f <- buildFormula(Y.cov, Z.cov)
EY.Z <- cf_gbm(model_f, model_data, test, min.depth, min.shrinkage)

# Find E(X|Z) with a boosted tree
model_f <- buildFormula(X.cov, Z.cov)
EX.Z <- cf_gbm(model_f, model_data, test, min.depth, min.shrinkage)

resid_data <- tibble(
  Y.star = log(model_data$CC6620) - EY.Z,
  X.star = unlist(model_data[,X.cov] - EX.Z)
)

model_gbm <- lm(Y.star ~ X.star, data = resid_data)
coef_gbm <- coeftest(model_gbm, vcov=vcovHC(model_gbm, type='HC2')) #[2,]
coef_gbm <- coef_gbm["X.star",] 

result <- result %>%
  bind_rows(tibble(Coefficient = focal_covariate, 
                   Method = "GBM",
                   Estimate = coef_gbm[1],
                   `Std. Error` = coef_gbm[2],
                   `t value` = coef_gbm[3],
                   `Pr(>|t|)` = coef_gbm[4],
                   `MSE E[Y|Z]` = mean(resid_data$Y.star^2),
                   `MSE E[X|Z]` = mean(resid_data$X.star^2)))

### BART

# Function to fit cross-fit BART
cf_bart <- function(dep_var, covariates, data, folds) {
  
  K <- length(folds)
  y.hat.full <- c()
  
  for (k in 1:K) {
    cat(paste0("  Fold ",k,"\n"))
    test_ix <- folds[[k]]
    
    X.train <- data[-test_ix, covariates] 
    X.train <- model.matrix(~ ., X.train)
    Y.train <- data[-test_ix, dep_var] %>% as.matrix
    X.test <- data[test_ix, covariates] 
    X.test <- model.matrix(~ ., X.test)
    
    bart.model <- bart(X.train, Y.train, keeptrees = TRUE)
    y.hat.bart <- predict(bart.model, newdata = X.test) %>% colMeans
    y.hat.full <- c(y.hat.full, y.hat.bart)
  }
  
  y.hat.ordered <- rep(NA, nrow(data))
  y.hat.ordered[unlist(folds)] <- y.hat.full
  return(y.hat.ordered)
}

model_data_bart <- model_data %>% 
  mutate(CC6620.log = log(CC6620))
bart_covariates_Z <- ols_covariates[ols_covariates != focal_covariate]

# Find E(Y|Z) with a cross-fit boosted tree
bart_Y <- "CC6620.log"
EY.Z <- cf_bart(dep_var = bart_Y, covariates = bart_covariates_Z, data = model_data_bart, folds = test)

# Find E(X|Z) with a boosted tree
bart_X <- focal_covariate
EX.Z <- cf_bart(dep_var = bart_X, covariates = bart_covariates_Z, data = model_data_bart, folds = test)

resid_data <- tibble(
  Y.star = log(model_data$CC6620) - EY.Z,
  X.star = unlist(model_data[,bart_X] - EX.Z)
)

model_bart <- lm(Y.star ~ X.star, data = resid_data)
coef_bart <- coeftest(model_bart, vcov=vcovHC(model_bart, type='HC2')) 
coef_bart <- coef_bart["X.star",] 
result <- result %>%
  bind_rows(tibble(Coefficient = focal_covariate, 
                   Method = "BART",
                   Estimate = coef_bart[1],
                   `Std. Error` = coef_bart[2],
                   `t value` = coef_bart[3],
                   `Pr(>|t|)` = coef_bart[4],
                   `MSE E[Y|Z]` = mean(resid_data$Y.star^2),
                   `MSE E[X|Z]` = mean(resid_data$X.star^2)))

####### LASSO #######
model_data_lasso <- model_data_bart
dep_var_lasso <- "CC6620.log"
lasso_covariates_Z <- ols_covariates[ols_covariates != focal_covariate]

# Function to fit cross-fit LASSO
cf_lasso <- function(dep_var, covariates, data, folds) {
  
  K <- length(folds)
  y.hat.full <- c()
  
  # data matrices
  X <- data[, covariates] 
  model_f_lasso <- buildFormula("", covariates, poly.2 = TRUE, cov_fct = c("Month","Weekend.Ind"))
  X <- model.matrix(model_f_lasso, X)
  Y <- data[, dep_var] %>% as.matrix
  
  # CV for finding lamda value
  cv.lasso <- cv.glmnet(X, Y, alpha=1)
  best.lam <- cv.lasso$lambda.min
  
  for (k in 1:K) {
    cat(paste0("  Fold ",k,"\n"))
    test_ix <- folds[[k]]
    
    lasso.model <- glmnet(X[-test_ix,], Y[-test_ix], alpha=1, lambda=best.lam)
    y.hat.lasso <- predict(lasso.model, s = best.lam, newx = X[test_ix,])
    y.hat.full <- c(y.hat.full, y.hat.lasso)
  }
  
  y.hat.ordered <- rep(NA, nrow(data))
  y.hat.ordered[unlist(folds)] <- y.hat.full
  return(y.hat.ordered)
}

# Find E(Y|Z) with a cross-fit boosted tree
lasso_Y <- "CC6620.log"
EY.Z <- cf_lasso(dep_var = lasso_Y, covariates = lasso_covariates_Z, data = model_data_lasso, folds = test)

# Find E(X|Z) with a boosted tree
lasso_X <- focal_covariate
EX.Z <- cf_lasso(dep_var = lasso_X, covariates = lasso_covariates_Z, data = model_data_lasso, folds = test)

resid_data <- tibble(
  Y.star = log(model_data$CC6620) - EY.Z,
  X.star = unlist(model_data[,lasso_X] - EX.Z)
)

model_lasso <- lm(Y.star ~ X.star, data = resid_data)
coef_lasso <- coeftest(model_lasso, vcov=vcovHC(model_lasso, type='HC2')) 
coef_lasso <- coef_lasso["X.star",] 
result <- result %>%
  bind_rows(tibble(Coefficient = focal_covariate, 
                   Method = "LASSO",
                   Estimate = coef_lasso[1],
                   `Std. Error` = coef_lasso[2],
                   `t value` = coef_lasso[3],
                   `Pr(>|t|)` = coef_lasso[4],
                   `MSE E[Y|Z]` = mean(resid_data$Y.star^2),
                   `MSE E[X|Z]` = mean(resid_data$X.star^2)))





