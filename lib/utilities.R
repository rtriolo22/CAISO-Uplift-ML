
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

# Run boosted tree for all shrinkage and depth values
tuneBoostedTree <- function(data, shrinkage_vals, depth_vals, num_trees) {
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
                            data = data[-test_ix,], 
                            distribution = "gaussian", n.trees = num_trees, interaction.depth = d, 
                            shrinkage = g)
        y.hat.boost <- predict(boost.uplift, newdata = data[test_ix,], n.trees = 5000)
        uplift.test <- data$CC6620[test_ix] %>% log
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
  
  return(result)
  
}
