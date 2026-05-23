# Nonparametric-Model-Information-Integrated-Multi-Index-Propensity-Score-for-Average-Causal-Effect
npMiPS GitHub code archive
==========================

This archive provides a minimal, runnable R implementation of the proposed
nonparametric model information-integrated multi-index propensity score (npMiPS)
workflow described in the manuscript.

Files
-----
1. npMiPS_functions.R
   Core functions for:
   - generating simulated data under the manuscript data-generating mechanism;
   - fitting parametric and ANN-based PS models;
   - fitting parametric and ANN-based outcome regression (OcR) models;
   - selecting ANN.PS, ANN.OcR, and integration-ANN hidden-layer structures in a
     single observed dataset;
   - constructing npMiPS from multiple PS/OcR model indexes;
   - estimating ACE by IPW using npMiPS;
   - bootstrap inference with ANN structures fixed after selection.

2. run_npMiPS_example_once.R
   A minimal example script. It generates one simulated dataset, selects ANN
   structures within that dataset, estimates npMiPS-111111, performs a small
   number of bootstrap resamples, and saves example outputs in the results folder.

Required R packages
-------------------
- AMORE
- MASS

You can install them with:
install.packages(c("AMORE", "MASS"))

How to run
----------
From the folder containing the R files, run:

Rscript run_npMiPS_example_once.R

Relationship to the manuscript
------------------------------
The simulated data-generating mechanism follows the manuscript setting:
- 12 baseline covariates are generated.
- X1-X4 are associated with both treatment and outcome.
- X5-X7 are associated only with treatment.
- X8-X10 are associated only with outcome.
- X11-X12 are independent noise covariates.
- alpha0 = -0.75 gives approximately 25% treated subjects.
- alpha0 = 0 gives approximately 50% treated subjects.
- The true ACE is 1.

The argument model_set controls the parametric model set used in npMiPS:
- model_set = "with_correct" corresponds to model sets A and B in the manuscript.
  The two parametric PS/OcR models include one correctly specified model and one
  misspecified model.
- model_set = "without_correct" corresponds to model sets P and M in the manuscript.
  The two parametric PS/OcR models are both misspecified.

The six-digit npMiPS code follows the manuscript notation. For example,
npMiPS-111111 includes ANN.PS, parametric PS Model 1, parametric PS Model 2,
ANN.OcR, parametric OcR Model 1, and parametric OcR Model 2.

Important note
--------------
The example script is intended to demonstrate the complete npMiPS workflow on one
simulated dataset. It is not intended to reproduce the full Monte Carlo simulation
tables in the manuscript. To run manuscript-scale simulations, increase the number
of Monte Carlo replicates and bootstrap resamples, and use:

candidate_structures <- default_hidden_candidates()

instead of the reduced candidate set used in the example script.
