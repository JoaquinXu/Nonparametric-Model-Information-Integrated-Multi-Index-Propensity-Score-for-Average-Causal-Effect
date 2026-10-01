# npMiPS: educational R code for the revised workflow

This repository provides a compact, runnable example of the revised
**nonparametric model information-integrated multi-index propensity score (npMiPS)**
workflow.

The code is designed for **method learning and reproducibility of the workflow**.
It is intentionally simpler than the simulation settings used in the manuscript:

- the example uses 6 baseline covariates instead of the larger manuscript setting;
- variable types, coefficients, nonlinear terms, and treatment prevalence are different;
- the example uses a small candidate architecture set and a small number of bootstrap samples so that it can run quickly;
- the example data-generating mechanism is therefore **not intended to reproduce any manuscript table**.

The main goal is to show how the revised procedure is implemented.

## Main revision reflected in this code

The integration ANN is selected by **3-fold out-of-fold mean absolute standardized mean difference (OOF-MASMD)** rather than by in-sample prediction accuracy for treatment (PAT).

For each candidate integration-ANN architecture:

1. split the data into treatment-stratified folds;
2. fit all nuisance models using the training folds only;
3. construct the model indexes in the training and held-out folds;
4. fit the integration ANN in the training folds;
5. predict propensity scores in the held-out fold;
6. combine the held-out predictions across folds;
7. calculate MASMD using ordinary ATE inverse-probability weights;
8. select the architecture with the smallest OOF-MASMD.

If a numerical tie occurs, the code prefers larger effective sample size (ESS), then the simpler ANN architecture.

The selected architecture is then refitted on the full observed dataset and held fixed during bootstrap resampling.

## Files

### `npMiPS_functions.R`
Core functions for:

- generating a small demonstration dataset;
- fitting parametric and ANN-based PS models;
- fitting parametric and ANN-based outcome regression models;
- selecting ANN.PS by MASMD;
- selecting ANN.OcR by mean observed prediction absolute error (MOPAE);
- selecting the integration ANN by 3-fold OOF-MASMD;
- calculating overlap, balance, calibration, ESS, and weight diagnostics;
- estimating ACE by normalized IPW;
- bootstrap inference with the selected ANN structures fixed.

### `01_run_npMiPS_OOF_demo.R`
A complete educational example of the revised npMiPS workflow.

The example:

1. simulates one dataset;
2. selects ANN.PS and ANN.OcR structures;
3. selects the integration ANN using OOF-MASMD;
4. estimates the ACE;
5. performs a small bootstrap;
6. saves the selected structures and diagnostic results.

### `02_run_crossfitted_SL_AIPW_demo.R`
A compact example of the cross-fitted Super Learner AIPW comparator used in the revision.

This script uses:

- outer 3-fold cross-fitting;
- Super Learner nuisance estimation for both PS and outcome regression;
- AIPW estimation from out-of-fold nuisance predictions;
- influence-function standard errors.

The candidate learner library is deliberately small in this educational example.

### `AMORE_0.2-15.tar.gz`
A local source archive retained from the previous code package for convenience when AMORE is not available through the user's usual package repository.

## Required R packages

### AMORE package

The npMiPS demonstration uses the `AMORE` package for fitting artificial neural networks.  
Because `AMORE` is currently available through the CRAN Archive rather than the main CRAN repository, users who do not already have it installed can obtain the archived source package from:

https://cran.r-project.org/src/contrib/Archive/AMORE/

For example:

```r
install.packages(
  "https://cran.r-project.org/src/contrib/Archive/AMORE/AMORE_0.2-15.tar.gz",
  repos = NULL,
  type = "source"
)

For the Super Learner AIPW demo:

```r
install.packages(c("SuperLearner", "nnet"))
```


### Compatibility note for `SuperLearner`

The Super Learner demo explicitly evaluates learner and screening functions in the
`SuperLearner` namespace. This avoids an `object 'All' not found` error that can occur
in some R/SuperLearner installations when the package is used through `::` without
being attached to the search path.

## How to run

The two demonstration scripts automatically try to locate their own folder when run with
RStudio, `source()`, or `Rscript`. Keep `npMiPS_functions.R` in the same folder as the two
demo scripts.

In RStudio, you can open either script and click **Source**, or run:

```r
source("/path/to/npMiPS_GitHub_code_v3_2/01_run_npMiPS_OOF_demo.R")
source("/path/to/npMiPS_GitHub_code_v3_2/02_run_crossfitted_SL_AIPW_demo.R")
```

If your R/RStudio environment cannot detect the script path automatically, first set the
working directory to the repository folder:

```r
setwd("/path/to/npMiPS_GitHub_code_v3_2")
source("01_run_npMiPS_OOF_demo.R")
source("02_run_crossfitted_SL_AIPW_demo.R")
```

Or from a terminal:

```bash
Rscript 01_run_npMiPS_OOF_demo.R
Rscript 02_run_crossfitted_SL_AIPW_demo.R
```

Outputs are written to the `results/` folder.

## Important implementation notes

- The OOF folds are stratified by treatment and are fixed across candidate integration-ANN architectures.
- ANN.PS and ANN.OcR structures are selected once in the observed dataset in this demonstration.
- During integration-ANN tuning, nuisance models are re-estimated within each outer training fold to avoid leakage from held-out observations.
- The propensity score is clipped only at `1e-7` and `1 - 1e-7` for numerical stability; this is not substantive truncation.
- OOF-MASMD uses ordinary ATE-IPW weights.
- Weight-stability diagnostics use stabilized weights.
- The bootstrap keeps the selected ANN structures fixed and therefore does not propagate architecture-selection uncertainty.
- The demonstration is **not a fully cross-fitted npMiPS estimator**. OOF is used for architecture tuning; the final npMiPS propensity score is refitted using the full observed dataset.

## Full candidate integration-ANN set

The demonstration uses only a few architectures for speed. The function
`full_integration_candidates()` returns the larger candidate set used by the revised manuscript workflow.

To use it, replace:

```r
integration_candidates <- demo_integration_candidates()
```

with:

```r
integration_candidates <- full_integration_candidates()
```

## Scope

The scripts are intended to make the workflow transparent and easy to learn. They should not be interpreted as exact code for reproducing every numerical result in the manuscript, because the demonstration data-generating mechanism and computational settings are intentionally simplified.
