################################################################################
# 02_run_crossfitted_SL_AIPW_demo.R
#
# Educational example of 3-fold cross-fitted Super Learner AIPW.
# This script is intentionally compact and is not intended to reproduce the
# manuscript simulation table exactly.
################################################################################

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


if (!requireNamespace("SuperLearner", quietly = TRUE)) {
  stop("Package 'SuperLearner' is required. Please install it first.")
}
if (!requireNamespace("nnet", quietly = TRUE)) {
  stop("Package 'nnet' is required for SL.nnet. Please install it first.")
}

if (!dir.exists("results")) dir.create("results", recursive = TRUE)

# -----------------------------------------------------------------------------
# 1. Generate the same small educational dataset used in the npMiPS demo
# -----------------------------------------------------------------------------

dat <- generate_demo_data(n = 600, seed = 2026, true_ate = 1.5)
Y <- dat$Y
A <- dat$A
X <- dat$X
n <- length(Y)

# -----------------------------------------------------------------------------
# 2. Outer cross-fitting setup
# -----------------------------------------------------------------------------

K <- 3
folds <- make_stratified_folds(A, K = K, seed = 2100)

# Small learner library for readability and speed.
# SL.glm provides a parametric learner and SL.nnet provides a flexible learner.
SL_library <- c("SL.glm", "SL.nnet")

ps_oof <- rep(NA_real_, n)
m1_oof <- rep(NA_real_, n)
m0_oof <- rep(NA_real_, n)
weight_rows <- list()

# -----------------------------------------------------------------------------
# 3. Cross-fit nuisance functions
# -----------------------------------------------------------------------------

for (k in seq_len(K)) {
  tr <- which(folds != k)
  te <- which(folds == k)

  # --- Propensity score Super Learner ---
  set.seed(seed_hash(2200, "ps", k))
  sl_ps <- SuperLearner::SuperLearner(
    Y = A[tr],
    X = X[tr, , drop = FALSE],
    family = binomial(),
    SL.library = SL_library,
    method = "method.NNLS",
    cvControl = list(V = 5),
    env = asNamespace("SuperLearner")
  )

  ps_oof[te] <- clip_ps(
    predict(sl_ps, newdata = X[te, , drop = FALSE])$pred,
    eps = 1e-6
  )

  # --- Outcome regression Super Learner ---
  Xout_tr <- data.frame(A = A[tr], X[tr, , drop = FALSE])
  Xout_1 <- data.frame(A = 1, X[te, , drop = FALSE])
  Xout_0 <- data.frame(A = 0, X[te, , drop = FALSE])

  set.seed(seed_hash(2200, "or", k))
  sl_or <- SuperLearner::SuperLearner(
    Y = Y[tr],
    X = Xout_tr,
    family = gaussian(),
    SL.library = SL_library,
    method = "method.NNLS",
    cvControl = list(V = 5),
    env = asNamespace("SuperLearner")
  )

  m1_oof[te] <- as.numeric(predict(sl_or, newdata = Xout_1)$pred)
  m0_oof[te] <- as.numeric(predict(sl_or, newdata = Xout_0)$pred)

  # Save ensemble weights for transparency.
  weight_rows[[length(weight_rows) + 1]] <- data.frame(
    Fold = k,
    Nuisance = "PS",
    Learner = names(sl_ps$coef),
    Weight = as.numeric(sl_ps$coef)
  )
  weight_rows[[length(weight_rows) + 1]] <- data.frame(
    Fold = k,
    Nuisance = "Outcome",
    Learner = names(sl_or$coef),
    Weight = as.numeric(sl_or$coef)
  )
}

if (anyNA(ps_oof) || anyNA(m1_oof) || anyNA(m0_oof)) {
  stop("Cross-fitted nuisance predictions are incomplete.")
}

# -----------------------------------------------------------------------------
# 4. AIPW estimate and influence-function standard error
# -----------------------------------------------------------------------------

psi <- m1_oof - m0_oof +
  A * (Y - m1_oof) / ps_oof -
  (1 - A) * (Y - m0_oof) / (1 - ps_oof)

ate <- mean(psi)
if_se <- sd(psi - ate) / sqrt(n)
lcl <- ate - 1.96 * if_se
ucl <- ate + 1.96 * if_se

summary_tab <- data.frame(
  Estimator = "Cross-fitted SL-AIPW",
  Estimate = ate,
  IF_SE = if_se,
  LCL = lcl,
  UCL = ucl,
  True_ATE = dat$true_ate,
  Bias = ate - dat$true_ate
)

print(summary_tab)

# Propensity-score diagnostics based on the OOF predictions.
ps_diag <- data.frame(t(ps_diagnostics(A, X, ps_oof, eps = 1e-6)))

weights_tab <- do.call(rbind, weight_rows)

write.csv(summary_tab,
          file.path("results", "demo_SL_AIPW_summary.csv"), row.names = FALSE)
write.csv(ps_diag,
          file.path("results", "demo_SL_AIPW_PS_diagnostics.csv"), row.names = FALSE)
write.csv(weights_tab,
          file.path("results", "demo_SL_AIPW_ensemble_weights.csv"), row.names = FALSE)

cat("\nDone. Cross-fitted SL-AIPW results were saved in the 'results' folder.\n")
