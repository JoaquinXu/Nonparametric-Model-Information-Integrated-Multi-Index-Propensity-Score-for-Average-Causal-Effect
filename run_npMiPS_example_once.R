################################################################################
# Example code for npMiPS
#
# This script demonstrates one complete, runnable npMiPS workflow:
#   1. generate one simulated dataset under the manuscript data-generating mechanism;
#   2. select ANN.PS, ANN.OcR, and integration-ANN structures within this dataset;
#   3. estimate npMiPS-111111;
#   4. run a small bootstrap only to verify the code structure;
#   5. save example outputs.
#
# This example is intentionally small and is not intended to reproduce the full
# 1000-replicate Monte Carlo simulation tables in the manuscript.
################################################################################

# Set the working directory to the folder containing this script when called by Rscript.
cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- sub(file_arg, "", cmd_args[grep(file_arg, cmd_args)])
if (length(script_path) > 0) {
  setwd(dirname(normalizePath(script_path)))
}

source("npMiPS_functions.R")

# -----------------------------------------------------------------------------
# User-adjustable settings
# -----------------------------------------------------------------------------

# model_set = "with_correct" corresponds to model sets A and B in the manuscript.
# model_set = "without_correct" corresponds to model sets P and M in the manuscript.
model_set <- "with_correct"

# alpha0 = -0.75 gives approximately 25% treatment prevalence.
# alpha0 = 0 gives approximately 50% treatment prevalence.
alpha0 <- -0.75

# The manuscript-scale analyses used the full candidate set in Appendix Table 1:
# candidate_structures <- default_hidden_candidates()
# To make this example fast, we use a reduced candidate set.
candidate_structures <- list(4, c(4, 4), c(4, 4, 4))

# The manuscript analyses used many bootstrap resamples.
# Here boot_num is deliberately small for a quick code check.
boot_num <- 5

# -----------------------------------------------------------------------------
# 1. Generate one simulated dataset
# -----------------------------------------------------------------------------

set.seed(123)
dat <- generate_npMiPS_data(seed = 123, n = 300, alpha0 = alpha0, true_ate = 1)
Y <- dat$Y
A <- dat$A
X <- dat$X

cat("Treatment prevalence:", mean(A), "\n")
cat("True ACE:", dat$true_ate, "\n")
cat("Model set:", model_set, "\n")

# -----------------------------------------------------------------------------
# 2. Select ANN structures within this single dataset
# -----------------------------------------------------------------------------

structure_selection <- select_single_dataset_structures(
  Y = Y,
  A = A,
  X = X,
  code = "111111",
  model_set = model_set,
  candidate_structures = candidate_structures
)

cat("Selected ANN.PS structure:", structure_to_string(structure_selection$h_ps_ann), "\n")
cat("Selected ANN.OcR structure:", structure_to_string(structure_selection$h_or_ann), "\n")
cat("Selected integration ANN structure:", structure_to_string(structure_selection$h_npMiPS), "\n")

# -----------------------------------------------------------------------------
# 3. Estimate ACE using npMiPS and bootstrap with fixed ANN structures
# -----------------------------------------------------------------------------
# In a real application, structures are selected once in the observed dataset
# and then fixed during bootstrap resampling.

result_npMiPS <- estimate_npMiPS_with_bootstrap(
  Y = Y,
  A = A,
  X = X,
  code = "111111",
  model_set = model_set,
  h_npMiPS = structure_selection$h_npMiPS,
  h_ps_ann = structure_selection$h_ps_ann,
  h_or_ann = structure_selection$h_or_ann,
  boot_num = boot_num,
  seed = 123,
  save_each_bootstrap = TRUE,
  output_dir = "results",
  file_prefix = "example_npMiPS"
)

cat("\nnpMiPS-111111 point estimate:\n")
print(result_npMiPS$point)

cat("\nBootstrap summary based on a small number of bootstrap resamples, for code checking only:\n")
print(result_npMiPS$summary)

# -----------------------------------------------------------------------------
# 4. Optional benchmark estimators on the same dataset
# -----------------------------------------------------------------------------

benchmark <- data.frame(
  Estimator = c("MPIPW.ANN", "MPIPW.correct", "G-comp.ANN", "G-comp.correct", "npMiPS-111111"),
  Estimate = c(
    estimate_mpipw(Y, A, X, ps_type = "ANN", hidden.neurons = structure_selection$h_ps_ann),
    estimate_mpipw(Y, A, X, ps_type = "correct"),
    estimate_gcomp(Y, A, X, or_type = "ANN", hidden.neurons = structure_selection$h_or_ann),
    estimate_gcomp(Y, A, X, or_type = "correct"),
    result_npMiPS$point
  )
)

if (!dir.exists("results")) dir.create("results", recursive = TRUE)
write.csv(benchmark, file = file.path("results", "example_estimates.csv"), row.names = FALSE)
save(result_npMiPS, structure_selection, benchmark,
     file = file.path("results", "example_npMiPS_result.RData"))

cat("\nBenchmark estimates:\n")
print(benchmark)
cat("\nFiles saved in the 'results' folder.\n")
