
# Compute number of elements per fold
computeFoldCounts <- function(n, K=5) {
  folds <- (n/K) %>% floor %>% rep(K)
  mod <- n%%K
  folds <- folds + c(rep(1, mod), rep(0,(K-mod)))
  return(folds)
}

# Sample from index list to create list of fold indices
drawFolds <- function(n, K=5) {
  fold_counts <- computeFoldCounts(n, K)
  sample_ix <- sample(1:n)
  test_ids <- list()
  for (i in 1:5) {
    test_ids[[i]] <- sample_ix[(sum(fold_counts[0:(i-1)])+1):sum(fold_counts[0:i])]
  }
  return(test_ids)
}

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

# Run boosted tree for all shrinkage and depth values
tuneBoostedTree <- function(data, dep_var, covariates, shrinkage_vals, depth_vals, num_trees, test_ids) {
  
  cat("\nRunning tuning of boosted gradient tree\n")
  cat(paste0("Dependent variable: log(", dep_var, ")\n\n"))
  result <- tibble()
  f <- paste0("log(", dep_var, ")") %>% buildFormula(covariates)
  
  for (d in depth_vals) {
    cat(paste0("Depth: ",d,"\n"))
    for (g in shrinkage_vals) {
      cat(paste0("  Shrinkage: ",g,"\n"))
      sq.error <- c()
      cat("    Fold: ")
      for (k in 1:(length(test_folds))) {
        cat(paste0(" ",k," "))
        test_ix <- test_ids[[k]]
        boost.uplift <- gbm(f, 
                            data = data[-test_ix,], 
                            distribution = "gaussian", n.trees = num_trees, interaction.depth = d, 
                            shrinkage = g)
        y.hat.boost <- boost.uplift %>% predict(newdata = data[test_ix,], n.trees = num_trees)
        uplift.test <- data[test_ix, dep_var] %>% log %>% unlist
        sq.error <- c(sq.error, (y.hat.boost - uplift.test)^2)
      }
      cat("\n")
      result_i <- tibble(
        Depth = d,
        Shrinkage = g,
        MSE = mean(sq.error)
      )
      result <- result %>% bind_rows(result_i)
    }
  }
  
  return(result)
  
}

# Function to fit cross-fit GBM
cf_gbm <- function(data, model_f, folds, min.depth, min.shrinkage, num_trees) {
  
  K <- length(folds)
  y.hat.full <- c()
  
  for (k in 1:K) {
    cat(paste0("  Fold ",k,"\n"))
    test_ix <- folds[[k]]
    boost.uplift <- gbm(model_f,
                        data = data[-test_ix,],
                        distribution = "gaussian", n.trees = num_trees,
                        interaction.depth = min.depth,
                        shrinkage = min.shrinkage)
    y.hat.boost <- predict(boost.uplift, newdata = data[test_ix,], n.trees = num_trees)
    y.hat.full <- c(y.hat.full, y.hat.boost)
  }
  
  y.hat.ordered <- rep(NA, nrow(data))
  y.hat.ordered[unlist(folds)] <- y.hat.full
  return(y.hat.ordered)
}

# Function to perform backward stepwise selection algorithm
backwardStepwiseSelection <- function(data, dep_var, covariates, test_ids, min_depth, min_shrinkage, 
                                      selection_threshold, num_trees) {
  max_iter <- length(covariates)
  
  # Fit full model
  model_f <- paste0("log(", dep_var, ")") %>% buildFormula(covariates)
  y_hat <- data %>% cf_gbm(model_f, test_ids, min_depth, min_shrinkage, num_trees)
  mse_min <- mean((y_hat - unlist(log(data[dep_var])))^2)
  
  # Save result to tibble
  result_table <- tibble(
    Vars.Dropped = 0,
    Var.Name = "Base",
    MSE = mse_min
  )
  
  # Start backward stepwise selection loop
  iter <- 1
  while(iter < max_iter) {
    mse_c <- c()
    for (c in covariates) {
      cat(paste0("Drop var: ",c,"\n"))
      covariates_i <- covariates[covariates != c]
      model_f <- paste0("log(", dep_var, ")") %>% buildFormula(covariates_i)
      y_hat <- data %>% cf_gbm(model_f, test_ids, min_depth, min_shrinkage, num_trees)
      mse_c <- c(mse_c, mean((y_hat - unlist(log(data[dep_var])))^2))
    }

    if (min(mse_c - mse_min) > selection_threshold) { break() }

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
  
  result <- list(result_table = result_table, covariates = covariates)
  
  return(result)
}

# Function to estimate a cross-fit full set of data
crossFitReducedGBM <- function(data, dep_var, covariates, test_ids, min_depth, min_shrinkage, num_trees, result_file) {
  load(result_file)
  drop_vars <- selection_result$result_table$Var.Name[2:(which.min(selection_result$result_table$MSE))]
  covariates <- covariates[!(covariates %in% drop_vars)]
  model_f <- paste0("log(", dep_var, ")") %>% buildFormula(covariates)
  y_hat <- data %>% cf_gbm(model_f, test_ids, min_depth, min_shrinkage, num_trees)
  return(y_hat)
}

# Function to fit reduced model, use all data, not cross-fit
fitReducedGBM <- function(data, dep_var, covariates, min_depth, min_shrinkage, num_trees, result_file) {
  load(result_file)
  drop_vars <- selection_result$result_table$Var.Name[2:(which.min(selection_result$result_table$MSE))]
  covariates <- covariates[!(covariates %in% drop_vars)]
  model_f <- paste0("log(", dep_var, ")") %>% buildFormula(covariates)
  boost_uplift_model <- gbm(model_f, data = data, distribution = "gaussian", n.trees = num_trees, 
                           interaction.depth = min_depth, 
                           shrinkage = min_shrinkage)
  return(boost_uplift_model)
}


