# CAISO-Uplift-ML
ML approach to assessment of uplift costs in the CAISO market

Reproduction code for working paper: "Market Conditions and Uplift Costs in the California Electricity Market"

**Usage**

Two files are provided to recreate results in the above paper:

- `uplift_model_boost.R`: Recreates boosted tree tuning, cross-fitting, backward stepwise selection, and full model fit
- `robinson_estimator.R`: Recreates results for marginal effect estimation

**`uplift_model_boost.R`**:

This file begins with a header that must be run to load necessary packages, load data, and select data folds

After this there are four labeled sections in the code:

1. Model tuning (*NOTE: this code fits many models and takes very long to be run. User can skip running this section and still run the remaining code.*)
2. Backward stepwise selection
3. Cross-fitting of the reduced model selected in part 2 (*Must run section 2 to run this section*)
4. Fit full model for variable influence (*Must run sections 2 and 3 to run this*)

**`robinson_estimator.R`**:

This file computes a Robinsom marginal effect estimate for BOOST, BART, and LASSO methods, and a conventional OLS for comparison.

To run the code specify the variable for which you would like to compute the estimate as `focal_covariate`.

Valid options are:

- `"Load.Mileage.MW"`
- `"Load.PGE"`
- `"DA.LoadError.Over"`
- `"Load.SDGE"`
- `"Load.SCE"`
- `"Convergence.Max.Supply"`
- `"Load.EIM"`
- `"Load.CAISO.Other"`
- `"Thermal_NP15"`
- `"ZP26_Solar.Actual"`
- `"Convergence.NetSupply"`
