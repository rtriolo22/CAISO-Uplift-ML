
# Print message upon running the 'uplift_model.R' script
beginMessage <- function() {
  cat("\n")
  cat("-------------------------------------------------------------------\n\n")
  cat("  Beginning running of script for reproduction of results for:\n")
  cat("  Triolo, R.C. and Wolak, F.A. (Forthcoming)\n")                                  
  cat("  Market Conditions and Uplift Costs in the California Electricity Market\n\n")
  cat("-------------------------------------------------------------------\n\n")
}

# Print message to indicate the dependent variable being run
dependentVariableMessage <- function(dep_var) {
  cat("Running process for dependent variable: ")
  cat(dep_var)
  cat("\n\n")
}

# Print message to screen to indicate the process being run
printProcessStatus <- function(settings) {
  if (settings$boost_tuning$run_boost_tuning) {
    cat("Running tuning procedure for boosted gradient tree...\n\n")
  }
  if (settings$boost_backward_selection$run_backward_stepwise) {
    cat("Running backward stepwise selection of covariates...\n\n")
  }
  if (settings$boost_reduced_cf$run_reduced_cf) {
    cat("Fitting boosted gradient tree on reduced model...\n\n")
  }
}

# Print message if adding one to dependent variable to avoid -Inf log values
printZeroValueAdjustment <- function(dep_var) {
  cat("Minimum value of dependent variable (")
  cat(dep_var)
  cat(") less than 1.\n")
  cat("Add 1 to this variable.\n")
}