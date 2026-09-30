################################################################################
# 01_run_npMiPS_OOF_demo.R
#
# Educational demonstration of the revised npMiPS workflow.
# The DGM and computational settings are intentionally simpler than the manuscript.
################################################################################

# Automatically move to the script folder when run with Rscript.
# Locate this script folder robustly in RStudio, source(), and Rscript.
get_script_dir <- function() {
  # 1) When called via source(), R may expose the current file in a source frame.
  frs <- sys.frames()
  for (i in rev(seq_along(frs))) {
    of <- frs[[i]]$ofile
    if (!is.null(of) && nzchar(of)) {
      p <- tryCatch(normalizePath(of, mustWork = TRUE), error = function(e) "")
      if (nzchar(p)) return(dirname(p))
    }
  }

  # 2) When called with Rscript, use the --file argument.
  args_all <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args_all, value = TRUE)
  if (length(file_arg) > 0L) {
    p <- tryCatch(
      normalizePath(sub("^--file=", "", file_arg[1]), mustWork = TRUE),
      error = function(e) ""
    )
    if (nzchar(p)) return(dirname(p))
  }

  # 3) When sourced interactively in RStudio, use the active editor path.
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    p <- tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) "")
    if (nzchar(p)) {
      p <- tryCatch(normalizePath(p, mustWork = TRUE), error = function(e) "")
      if (nzchar(p)) return(dirname(p))
    }
  }

  # 4) Final fallback: current working directory.
  getwd()
}

script_dir <- get_script_dir()
helper_file <- file.path(script_dir, "npMiPS_functions.R")
if (!file.exists(helper_file)) {
  stop(
    "Cannot locate 'npMiPS_functions.R'. Please keep all repository files in the same folder, ",
    "or set the R working directory to that folder before running the script."
  )
}
setwd(script_dir)
source(helper_file)


if (!dir.exists("results")) dir.create("results", recursive = TRUE)

# -----------------------------------------------------------------------------
# User-adjustable settings
# -----------------------------------------------------------------------------

n <- 600
true_ate <- 1.5
code <- "111111"
K <- 3

# Small candidate sets for a quick educational run.
ann_candidates <- demo_ann_candidates()
integration_candidates <- demo_integration_candidates()

# For a larger integration-ANN search, use:
# integration_candidates <- full_integration_candidates()

# Kept small so the example finishes quickly.
bootstrap_B <- 20

# -----------------------------------------------------------------------------
# 1. Generate demonstration data
# -----------------------------------------------------------------------------

dat <- generate_demo_data(n = n, seed = 2026, true_ate = true_ate)
Y <- dat$Y
A <- dat$A
X <- dat$X

cat("Sample size:", n, "\n")
cat("Treatment prevalence:", round(mean(A), 3), "\n")
cat("True ACE:", true_ate, "\n\n")

# -----------------------------------------------------------------------------
# 2. Select ANN.PS and ANN.OcR structures in the observed dataset
# -----------------------------------------------------------------------------

ps_sel <- select_ann_ps_structure(A, X, ann_candidates, seed = 1100)
or_sel <- select_ann_or_structure(Y, A, X, ann_candidates, seed = 1200)

cat("Selected ANN.PS structure:", structure_label(ps_sel$hidden), "\n")
cat("Selected ANN.OcR structure:", structure_label(or_sel$hidden), "\n\n")

# -----------------------------------------------------------------------------
# 3. Select integration ANN by 3-fold OOF-MASMD
# -----------------------------------------------------------------------------

int_sel <- select_integration_oof_masmd(
  Y = Y,
  A = A,
  X = X,
  code = code,
  h_ps_ann = ps_sel$hidden,
  h_or_ann = or_sel$hidden,
  candidate_structures = integration_candidates,
  K = K,
  seed = 1300
)

cat("Selected integration ANN structure:", structure_label(int_sel$hidden), "\n")
cat("OOF-MASMD:", round(int_sel$oof_diagnostics[["MASMD"]], 4), "\n")
cat("OOF-MaxASMD:", round(int_sel$oof_diagnostics[["MaxASMD"]], 4), "\n")
cat("OOF-ESS:", round(int_sel$oof_diagnostics[["ESS_total"]], 1), "\n\n")

# -----------------------------------------------------------------------------
# 4. Refit selected npMiPS model on the full dataset
# -----------------------------------------------------------------------------

final_fit <- fit_npMiPS_full(
  Y = Y,
  A = A,
  X = X,
  code = code,
  h_ps_ann = ps_sel$hidden,
  h_or_ann = or_sel$hidden,
  h_integration = int_sel$hidden,
  seed = 1400
)

cat("npMiPS point estimate:", round(final_fit$ate, 4), "\n")
cat("True ACE:", true_ate, "\n\n")

# -----------------------------------------------------------------------------
# 5. Small bootstrap with selected structures fixed
# -----------------------------------------------------------------------------

boots <- bootstrap_npMiPS(
  Y = Y,
  A = A,
  X = X,
  code = code,
  h_ps_ann = ps_sel$hidden,
  h_or_ann = or_sel$hidden,
  h_integration = int_sel$hidden,
  B = bootstrap_B,
  seed = 1500
)

boot_sum <- bootstrap_summary(final_fit$ate, boots)
print(boot_sum)

# -----------------------------------------------------------------------------
# 6. Save outputs
# -----------------------------------------------------------------------------

write.csv(ps_sel$table,
          file.path("results", "demo_ANN_PS_selection.csv"), row.names = FALSE)
write.csv(or_sel$table,
          file.path("results", "demo_ANN_OcR_selection.csv"), row.names = FALSE)
write.csv(int_sel$candidate_table,
          file.path("results", "demo_integration_OOF_selection.csv"), row.names = FALSE)
write.csv(data.frame(t(final_fit$diagnostics)),
          file.path("results", "demo_full_data_PS_diagnostics.csv"), row.names = FALSE)
write.csv(boot_sum,
          file.path("results", "demo_npMiPS_summary.csv"), row.names = FALSE)
write.csv(data.frame(bootstrap_estimate = boots),
          file.path("results", "demo_npMiPS_bootstrap_estimates.csv"), row.names = FALSE)

selected_structures <- data.frame(
  Component = c("ANN.PS", "ANN.OcR", "Integration ANN"),
  Structure = c(
    structure_label(ps_sel$hidden),
    structure_label(or_sel$hidden),
    structure_label(int_sel$hidden)
  )
)
write.csv(selected_structures,
          file.path("results", "demo_selected_structures.csv"), row.names = FALSE)

cat("\nDone. Results were saved in the 'results' folder.\n")
