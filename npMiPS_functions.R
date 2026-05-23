################################################################################
# npMiPS functions
#
# This file contains the core functions for the nonparametric model
# information-integrated multi-index propensity score (npMiPS) estimator.
# It is written as a self-contained function library for the GitHub code archive.
#
# Required packages: AMORE, MASS
################################################################################

require_package <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required. Please install it before running this code.", pkg),
         call. = FALSE)
  }
}

require_package("AMORE")
require_package("MASS")

################################################################################
# Utility functions
################################################################################

expit <- function(x) 1 / (1 + exp(-x))

clip_probability <- function(p, eps = 1e-7) {
  p <- as.numeric(p)
  p[p < eps] <- eps
  p[p > 1 - eps] <- 1 - eps
  p
}

logit_safe <- function(p, eps = 1e-7) {
  p <- clip_probability(p, eps = eps)
  log(p / (1 - p))
}

weighted_sd <- function(x, w) {
  m <- sum(w * x) / sum(w)
  sqrt(sum(w * (x - m)^2) / sum(w))
}

absolute_mean_difference <- function(x, A, w) {
  abs(sum(x * A * w) / sum(A * w) -
        sum(x * (1 - A) * w) / sum((1 - A) * w))
}

estimate_ipw <- function(Y, A, ps, eps = 1e-7) {
  ps <- clip_probability(ps, eps = eps)
  mu1 <- sum(Y * A / ps) / sum(A / ps)
  mu0 <- sum(Y * (1 - A) / (1 - ps)) / sum((1 - A) / (1 - ps))
  as.numeric(mu1 - mu0)
}

bootstrap_summary <- function(point, boots, null_value = 0, digits = 3) {
  se <- stats::sd(boots)
  lcl <- point - 1.96 * se
  ucl <- point + 1.96 * se
  p_value <- 2 * stats::pt(abs((point - null_value) / se),
                           df = length(boots) - 1,
                           lower.tail = FALSE)
  data.frame(
    Estimate = round(point, digits),
    BSSE = round(se, digits),
    LCL = round(lcl, digits),
    UCL = round(ucl, digits),
    P_value = ifelse(p_value < 0.001, "<0.001", as.character(round(p_value, digits))),
    stringsAsFactors = FALSE
  )
}

# Candidate ANN hidden-layer structures used in the manuscript.
default_hidden_candidates <- function() {
  list(
    2, 3, 4, 5, 6, 7, 8, 9,
    c(2, 2), c(3, 3), c(4, 4), c(5, 5), c(6, 6), c(7, 7), c(8, 8), c(9, 9),
    c(2, 2, 2), c(3, 3, 3), c(4, 4, 4), c(5, 5, 5), c(6, 6, 6),
    c(2, 2, 2, 2), c(3, 3, 3, 3), c(4, 4, 4, 4), c(5, 5, 5, 5), c(6, 6, 6, 6),
    c(2, 2, 2, 2, 2), c(3, 3, 3, 3, 3), c(4, 4, 4, 4, 4),
    c(5, 5, 5, 5, 5), c(6, 6, 6, 6, 6)
  )
}

structure_to_string <- function(h) paste(h, collapse = "-")

################################################################################
# Data-generating mechanism used in the simulation example
#
# This function follows the data-generating setting described in the manuscript:
# - 12 baseline covariates are generated.
# - X1-X4 are associated with both treatment and outcome.
# - X5-X7 are associated only with treatment.
# - X8-X10 are associated only with outcome.
# - X11-X12 are independent noise covariates.
# - alpha0 = -0.75 gives approximately 25% treated subjects.
# - alpha0 = 0 gives approximately 50% treated subjects.
# - true_ate = 1 is the true average causal effect.
################################################################################

generate_npMiPS_data <- function(seed = 1,
                                n = 300,
                                alpha0 = -0.75,
                                true_ate = 1) {
  set.seed(seed)

  alpha <- c(0.18, -0.17, -0.08, 0.19, -0.16, 0.20, -0.11, rep(0, 3))
  beta0 <- -1.2
  beta <- c(0.58, -0.29, -0.58, 0.65, rep(0, 3), 0.67, -0.31, 0.62)

  X1.5 <- MASS::mvrnorm(n, mu = c(0, 0), Sigma = matrix(c(1, 0.2, 0.2, 1), 2))
  X1 <- X1.5[, 1]
  X5 <- X1.5[, 2]

  X2.6 <- MASS::mvrnorm(n, mu = c(0, 0), Sigma = matrix(c(1, 0.9, 0.9, 1), 2))
  X2 <- X2.6[, 1]
  X6 <- X2.6[, 2]

  X3.8 <- MASS::mvrnorm(n, mu = c(0, 0), Sigma = matrix(c(1, 0.2, 0.2, 1), 2))
  X3 <- X3.8[, 1]
  X8 <- X3.8[, 2]

  X4.9 <- MASS::mvrnorm(n, mu = c(0, 0), Sigma = matrix(c(1, 0.9, 0.9, 1), 2))
  X4 <- X4.9[, 1]
  X9 <- X4.9[, 2]

  X7 <- stats::rnorm(n)
  X10 <- stats::rnorm(n)

  X2 <- ifelse(X2 > mean(X2), 1, 0)
  X5 <- ifelse(X5 > mean(X5), 1, 0)
  X8 <- ifelse(X8 > mean(X8), 1, 0)

  X_for_dgm <- cbind(X1, X2, X3, X4, X5, X6, X7, X8, X9, X10)
  ps_true <- expit(as.numeric(alpha0 + X_for_dgm %*% alpha))
  A <- stats::rbinom(n, size = 1, prob = ps_true)

  Y1 <- true_ate + beta0 + as.numeric(X_for_dgm %*% beta)
  Y0 <- beta0 + as.numeric(X_for_dgm %*% beta)
  Y <- A * Y1 + (1 - A) * Y0 + stats::rnorm(n)

  X11 <- stats::rnorm(n)
  X12 <- stats::rbinom(n, size = 1, prob = 0.5)
  X <- cbind(X1, X2, X3, X4, X5, X6, X7, X8, X9, X10, X11, X12)
  colnames(X) <- paste0("X", seq_len(ncol(X)))

  list(Y = as.numeric(Y), A = as.numeric(A), X = X, true_ps = ps_true, true_ate = true_ate)
}

################################################################################
# Single-model PS and OcR estimators
################################################################################

# Parametric PS covariate sets used in the manuscript simulation.
# correct    : true PS covariates, X1-X7.
# incorrect1 : wrong covariate selection, X1, X2, X5, X11.
# incorrect2 : functional-form misspecification, squared X1-X7.
get_ps_covariates <- function(X, model = c("correct", "incorrect1", "incorrect2")) {
  model <- match.arg(model)
  if (model == "correct") {
    return(as.data.frame(X[, 1:7, drop = FALSE]))
  }
  if (model == "incorrect1") {
    return(as.data.frame(X[, c(1, 2, 5, 11), drop = FALSE]))
  }
  as.data.frame((X^2)[, 1:7, drop = FALSE])
}

# Parametric OcR covariate sets used in the manuscript simulation.
# correct    : true OcR covariates, X1-X4 and X8-X10.
# incorrect1 : wrong covariate selection, X1, X2, X8, X11.
# incorrect2 : functional-form misspecification, squared true OcR covariates.
get_or_covariates <- function(X, model = c("correct", "incorrect1", "incorrect2")) {
  model <- match.arg(model)
  if (model == "correct") {
    return(as.data.frame(X[, c(1:4, 8:10), drop = FALSE]))
  }
  if (model == "incorrect1") {
    return(as.data.frame(X[, c(1, 2, 8, 11), drop = FALSE]))
  }
  as.data.frame((X^2)[, c(1:4, 8:10), drop = FALSE])
}

fit_parametric_ps <- function(A, X, model = c("correct", "incorrect1", "incorrect2")) {
  Z <- get_ps_covariates(X, model = model)
  dat <- data.frame(A = A, Z)
  fit <- stats::glm(A ~ ., data = dat, family = stats::binomial())
  ps <- as.numeric(stats::predict(fit, type = "response"))
  index <- as.numeric(as.matrix(Z) %*% stats::coef(fit)[-1])
  list(fit = fit, ps = clip_probability(ps), index = index)
}

fit_parametric_or <- function(Y, A, X, model = c("correct", "incorrect1", "incorrect2")) {
  Z <- get_or_covariates(X, model = model)
  dat <- data.frame(Y = Y, A = A, Z)
  fit <- stats::lm(Y ~ ., data = dat)

  dat1 <- dat
  dat1$A <- 1
  dat0 <- dat
  dat0$A <- 0
  mu1 <- stats::predict(fit, newdata = dat1)
  mu0 <- stats::predict(fit, newdata = dat0)

  coefs <- stats::coef(fit)
  z_names <- colnames(Z)
  index <- as.numeric(as.matrix(Z) %*% coefs[z_names])

  list(fit = fit, mu1 = as.numeric(mu1), mu0 = as.numeric(mu0),
       ate = mean(mu1 - mu0), index = index)
}

fit_ann_ps <- function(A, X,
                       hidden.neurons = 5,
                       learning.rate.global = 0.001,
                       momentum.global = 0.5,
                       error.criterium = "LMS",
                       hidden.layer = "tansig",
                       output.layer = "sigmoid",
                       method = "ADAPTgdwm") {
  dat <- data.frame(A = A, X)
  n.neurons <- c(ncol(dat) - 1, hidden.neurons, 1)
  net <- AMORE::newff(n.neurons,
                      learning.rate.global = learning.rate.global,
                      momentum.global = momentum.global,
                      error.criterium = error.criterium,
                      Stao = NA,
                      hidden.layer = hidden.layer,
                      output.layer = output.layer,
                      method = method)
  fit <- AMORE::train(net, dat[, -1, drop = FALSE], as.numeric(dat[, 1]),
                      error.criterium = error.criterium,
                      report = FALSE,
                      show.step = 100,
                      n.shows = 5)
  ps <- clip_probability(AMORE::sim(fit$net, dat[, -1, drop = FALSE]))
  list(net = fit$net, ps = ps, index = logit_safe(ps))
}

fit_ann_or <- function(Y, A, X,
                       hidden.neurons = c(9, 9),
                       learning.rate.global = 0.001,
                       momentum.global = 0.5,
                       error.criterium = "LMS",
                       hidden.layer = "tansig",
                       output.layer = "purelin",
                       method = "ADAPTgdwm") {
  dat <- data.frame(Y = Y, A = A, X)
  n.neurons <- c(ncol(dat) - 1, hidden.neurons, 1)
  net <- AMORE::newff(n.neurons,
                      learning.rate.global = learning.rate.global,
                      momentum.global = momentum.global,
                      error.criterium = error.criterium,
                      Stao = NA,
                      hidden.layer = hidden.layer,
                      output.layer = output.layer,
                      method = method)
  fit <- AMORE::train(net, dat[, -1, drop = FALSE], as.numeric(dat[, 1]),
                      error.criterium = error.criterium,
                      report = FALSE,
                      show.step = 100,
                      n.shows = 5)

  dat1 <- dat
  dat1$A <- 1
  dat0 <- dat
  dat0$A <- 0
  mu1 <- as.numeric(AMORE::sim(fit$net, dat1[, -1, drop = FALSE]))
  mu0 <- as.numeric(AMORE::sim(fit$net, dat0[, -1, drop = FALSE]))
  fitted_y <- as.numeric(AMORE::sim(fit$net, dat[, -1, drop = FALSE]))

  list(net = fit$net, mu1 = mu1, mu0 = mu0, ate = mean(mu1 - mu0),
       index = mu0, fitted_y = fitted_y)
}

estimate_mpipw <- function(Y, A, X,
                           ps_type = c("ANN", "correct", "incorrect1", "incorrect2"),
                           hidden.neurons = 5) {
  ps_type <- match.arg(ps_type)
  if (ps_type == "ANN") {
    ps <- fit_ann_ps(A, X, hidden.neurons = hidden.neurons)$ps
  } else {
    ps <- fit_parametric_ps(A, X, model = ps_type)$ps
  }
  estimate_ipw(Y, A, ps)
}

estimate_gcomp <- function(Y, A, X,
                           or_type = c("ANN", "correct", "incorrect1", "incorrect2"),
                           hidden.neurons = c(9, 9)) {
  or_type <- match.arg(or_type)
  if (or_type == "ANN") {
    return(fit_ann_or(Y, A, X, hidden.neurons = hidden.neurons)$ate)
  }
  fit_parametric_or(Y, A, X, model = or_type)$ate
}

################################################################################
# npMiPS construction
################################################################################

# Mapping between example code and manuscript model sets.
#
# model_set = "with_correct" corresponds to model sets A and B in the manuscript:
#   pi1(X) = ANN.PS, pi2(X) = correctly specified parametric PS,
#   pi3(X) = misspecified parametric PS,
#   mA1(X) = ANN.OcR, mA2(X) = correctly specified parametric OcR,
#   mA3(X) = misspecified parametric OcR.
#
# model_set = "without_correct" corresponds to model sets P and M in the manuscript:
#   pi1(X) = ANN.PS, pi2(X) = misspecified PS by wrong covariate selection,
#   pi3(X) = misspecified PS by wrong functional form,
#   mA1(X) = ANN.OcR, mA2(X) = misspecified OcR by wrong covariate selection,
#   mA3(X) = misspecified OcR by wrong functional form.
model_set_map <- function(model_set = c("with_correct", "without_correct")) {
  model_set <- match.arg(model_set)
  if (model_set == "with_correct") {
    return(list(ps1 = "correct", ps2 = "incorrect2",
                or1 = "correct", or2 = "incorrect2"))
  }
  list(ps1 = "incorrect1", ps2 = "incorrect2",
       or1 = "incorrect1", or2 = "incorrect2")
}

validate_code <- function(code) {
  code <- as.character(code)
  if (!grepl("^[01]{6}$", code)) {
    stop("The model code must be a six-digit string containing only 0 and 1, e.g., '100100' or '111111'.",
         call. = FALSE)
  }
  as.integer(strsplit(code, "")[[1]])
}

build_npMiPS_indexes <- function(Y, A, X,
                                 code = "111111",
                                 model_set = c("with_correct", "without_correct"),
                                 h_ps_ann = 5,
                                 h_or_ann = c(9, 9)) {
  digits <- validate_code(code)
  models <- model_set_map(model_set)
  index_list <- list(A = A)

  if (digits[1] == 1) {
    index_list$PS_ANN <- fit_ann_ps(A, X, hidden.neurons = h_ps_ann)$index
  }
  if (digits[2] == 1) {
    index_list$PS_1 <- fit_parametric_ps(A, X, model = models$ps1)$index
  }
  if (digits[3] == 1) {
    index_list$PS_2 <- fit_parametric_ps(A, X, model = models$ps2)$index
  }
  if (digits[4] == 1) {
    index_list$OcR_ANN <- fit_ann_or(Y, A, X, hidden.neurons = h_or_ann)$index
  }
  if (digits[5] == 1) {
    index_list$OcR_1 <- fit_parametric_or(Y, A, X, model = models$or1)$index
  }
  if (digits[6] == 1) {
    index_list$OcR_2 <- fit_parametric_or(Y, A, X, model = models$or2)$index
  }

  if (length(index_list) == 1) {
    stop("No PS or OcR model was selected. At least one digit in the model code must be 1.",
         call. = FALSE)
  }

  as.data.frame(index_list)
}

fit_npMiPS <- function(Y, A, X,
                       code = "111111",
                       model_set = c("with_correct", "without_correct"),
                       h_npMiPS = 4,
                       h_ps_ann = 5,
                       h_or_ann = c(9, 9),
                       learning.rate.global = 0.001,
                       momentum.global = 0.5,
                       error.criterium = "LMS",
                       hidden.layer = "tansig",
                       output.layer = "sigmoid",
                       method = "ADAPTgdwm") {
  indexes <- build_npMiPS_indexes(Y, A, X,
                                  code = code,
                                  model_set = model_set,
                                  h_ps_ann = h_ps_ann,
                                  h_or_ann = h_or_ann)
  x_train <- indexes[, -1, drop = FALSE]
  n.neurons <- c(ncol(x_train), h_npMiPS, 1)
  net <- AMORE::newff(n.neurons,
                      learning.rate.global = learning.rate.global,
                      momentum.global = momentum.global,
                      error.criterium = error.criterium,
                      Stao = NA,
                      hidden.layer = hidden.layer,
                      output.layer = output.layer,
                      method = method)
  fit <- AMORE::train(net, x_train, as.numeric(indexes[, 1]),
                      error.criterium = error.criterium,
                      report = FALSE,
                      show.step = 100,
                      n.shows = 5)
  ps_npMiPS <- clip_probability(AMORE::sim(fit$net, x_train))
  list(net = fit$net, ps = ps_npMiPS, indexes = indexes)
}

estimate_npMiPS <- function(Y, A, X,
                            code = "111111",
                            model_set = c("with_correct", "without_correct"),
                            h_npMiPS = 4,
                            h_ps_ann = 5,
                            h_or_ann = c(9, 9)) {
  fit <- fit_npMiPS(Y, A, X,
                    code = code,
                    model_set = model_set,
                    h_npMiPS = h_npMiPS,
                    h_ps_ann = h_ps_ann,
                    h_or_ann = h_or_ann)
  estimate_ipw(Y, A, fit$ps)
}

################################################################################
# ANN structure selection on one observed dataset
################################################################################

select_ann_ps_structure <- function(A, X,
                                    candidate_structures = default_hidden_candidates()) {
  sd_treated <- apply(X[A == 1, , drop = FALSE], 2, stats::sd)
  scores <- numeric(length(candidate_structures))

  for (i in seq_along(candidate_structures)) {
    h <- candidate_structures[[i]]
    ps <- fit_ann_ps(A, X, hidden.neurons = h)$ps
    sw <- ifelse(A == 1, mean(A) / ps, (1 - mean(A)) / (1 - ps))
    amd <- apply(X, 2, absolute_mean_difference, A = A, w = sw)
    asmd <- amd / sd_treated
    scores[i] <- mean(asmd, na.rm = TRUE)
  }

  best <- which.min(scores)
  list(best_structure = candidate_structures[[best]],
       criterion = "MASMD",
       scores = data.frame(structure = sapply(candidate_structures, structure_to_string),
                           MASMD = scores,
                           stringsAsFactors = FALSE))
}

select_ann_or_structure <- function(Y, A, X,
                                    candidate_structures = default_hidden_candidates()) {
  scores <- numeric(length(candidate_structures))

  for (i in seq_along(candidate_structures)) {
    h <- candidate_structures[[i]]
    fit <- fit_ann_or(Y, A, X, hidden.neurons = h)
    scores[i] <- mean(abs(Y - fit$fitted_y))
  }

  best <- which.min(scores)
  list(best_structure = candidate_structures[[best]],
       criterion = "MOPAE",
       scores = data.frame(structure = sapply(candidate_structures, structure_to_string),
                           MOPAE = scores,
                           stringsAsFactors = FALSE))
}

select_npMiPS_structure <- function(Y, A, X,
                                    code = "111111",
                                    model_set = c("with_correct", "without_correct"),
                                    h_ps_ann = 5,
                                    h_or_ann = c(9, 9),
                                    candidate_structures = default_hidden_candidates()) {
  scores <- numeric(length(candidate_structures))

  for (i in seq_along(candidate_structures)) {
    h <- candidate_structures[[i]]
    ps <- fit_npMiPS(Y, A, X,
                     code = code,
                     model_set = model_set,
                     h_npMiPS = h,
                     h_ps_ann = h_ps_ann,
                     h_or_ann = h_or_ann)$ps
    pred_A <- ifelse(ps > 0.5, 1, 0)
    scores[i] <- mean(pred_A == A)
  }

  best <- which.max(scores)
  list(best_structure = candidate_structures[[best]],
       criterion = "PAT",
       scores = data.frame(structure = sapply(candidate_structures, structure_to_string),
                           PAT = scores,
                           stringsAsFactors = FALSE))
}

select_single_dataset_structures <- function(Y, A, X,
                                             code = "111111",
                                             model_set = c("with_correct", "without_correct"),
                                             candidate_structures = default_hidden_candidates(),
                                             default_ps_structure = 5,
                                             default_or_structure = c(9, 9)) {
  digits <- validate_code(code)

  ps_sel <- NULL
  or_sel <- NULL
  h_ps <- default_ps_structure
  h_or <- default_or_structure

  if (digits[1] == 1) {
    ps_sel <- select_ann_ps_structure(A, X, candidate_structures)
    h_ps <- ps_sel$best_structure
  }

  if (digits[4] == 1) {
    or_sel <- select_ann_or_structure(Y, A, X, candidate_structures)
    h_or <- or_sel$best_structure
  }

  np_sel <- select_npMiPS_structure(Y, A, X,
                                    code = code,
                                    model_set = model_set,
                                    h_ps_ann = h_ps,
                                    h_or_ann = h_or,
                                    candidate_structures = candidate_structures)

  list(h_ps_ann = h_ps,
       h_or_ann = h_or,
       h_npMiPS = np_sel$best_structure,
       ps_selection = ps_sel,
       or_selection = or_sel,
       npMiPS_selection = np_sel)
}

################################################################################
# Point estimate and bootstrap using fixed ANN structures
################################################################################

estimate_npMiPS_with_bootstrap <- function(Y, A, X,
                                           code = "111111",
                                           model_set = c("with_correct", "without_correct"),
                                           h_npMiPS = 4,
                                           h_ps_ann = 5,
                                           h_or_ann = c(9, 9),
                                           boot_num = 100,
                                           seed = 2026,
                                           save_each_bootstrap = FALSE,
                                           output_dir = NULL,
                                           file_prefix = "npMiPS") {
  point <- estimate_npMiPS(Y, A, X,
                           code = code,
                           model_set = model_set,
                           h_npMiPS = h_npMiPS,
                           h_ps_ann = h_ps_ann,
                           h_or_ann = h_or_ann)

  boots <- rep(NA_real_, boot_num)
  dat <- data.frame(Y = Y, A = A, X)

  if (save_each_bootstrap) {
    if (is.null(output_dir)) output_dir <- getwd()
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  }

  for (b in seq_len(boot_num)) {
    set.seed(seed + b)
    id <- sample(seq_len(nrow(dat)), size = nrow(dat), replace = TRUE)
    boot_dat <- dat[id, , drop = FALSE]
    Yb <- boot_dat$Y
    Ab <- boot_dat$A
    Xb <- as.matrix(boot_dat[, -(1:2), drop = FALSE])

    boots[b] <- estimate_npMiPS(Yb, Ab, Xb,
                                code = code,
                                model_set = model_set,
                                h_npMiPS = h_npMiPS,
                                h_ps_ann = h_ps_ann,
                                h_or_ann = h_or_ann)

    if (save_each_bootstrap) {
      save(point, boots,
           file = file.path(output_dir,
                            sprintf("%s_%s_bootstrap_progress.RData", file_prefix, code)))
    }
  }

  list(point = point,
       boots = boots,
       summary = bootstrap_summary(point, boots),
       code = code,
       model_set = match.arg(model_set),
       h_npMiPS = h_npMiPS,
       h_ps_ann = h_ps_ann,
       h_or_ann = h_or_ann)
}
