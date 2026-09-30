################################################################################
# npMiPS educational functions
#
# Purpose
# -------
# A compact implementation of the revised npMiPS workflow for teaching and
# code transparency. The demonstration data-generating mechanism is intentionally
# simpler and different from that used in the manuscript.
#
# Required package: AMORE
################################################################################

require_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required but is not installed.", pkg), call. = FALSE)
  }
}

require_pkg("AMORE")

################################################################################
# 1. General utilities
################################################################################

expit <- function(x) 1 / (1 + exp(-x))

clip_ps <- function(p, eps = 1e-7) {
  pmin(pmax(as.numeric(p), eps), 1 - eps)
}

structure_label <- function(h) paste(as.integer(h), collapse = "-")

structure_complexity <- function(h) {
  # Fewer hidden layers first, then fewer total hidden nodes.
  c(n_layers = length(h), total_nodes = sum(h))
}

seed_hash <- function(...) {
  vals <- unlist(list(...), use.names = FALSE)
  vals <- vapply(vals, function(z) {
    if (is.character(z)) sum(utf8ToInt(z)) else as.numeric(z)
  }, numeric(1))
  h <- 104729
  for (z in vals) h <- (h * 1009 + z * 9176 + 12345) %% 2147483000
  as.integer(h + 1)
}

make_stratified_folds <- function(A, K = 3, seed = 1001) {
  fold <- integer(length(A))
  for (g in c(0, 1)) {
    id <- which(A == g)
    if (length(id) < K) stop("A treatment group has fewer observations than folds.")
    set.seed(seed_hash(seed, g))
    id <- sample(id, length(id), replace = FALSE)
    fold[id] <- rep(seq_len(K), length.out = length(id))
  }
  fold
}

weighted_mean_safe <- function(x, w) {
  sw <- sum(w)
  if (!is.finite(sw) || sw <= 0) return(NA_real_)
  sum(x * w) / sw
}

effective_sample_size <- function(w) {
  if (length(w) == 0 || any(!is.finite(w)) || sum(w^2) <= 0) return(NA_real_)
  (sum(w)^2) / sum(w^2)
}

ate_weights <- function(A, ps, eps = 1e-7) {
  ps <- clip_ps(ps, eps)
  ifelse(A == 1, 1 / ps, 1 / (1 - ps))
}

stabilized_weights <- function(A, ps, eps = 1e-7) {
  ps <- clip_ps(ps, eps)
  pa <- mean(A)
  ifelse(A == 1, pa / ps, (1 - pa) / (1 - ps))
}

estimate_ipw_ate <- function(Y, A, ps, eps = 1e-7) {
  ps <- clip_ps(ps, eps)
  mu1 <- sum(A * Y / ps) / sum(A / ps)
  mu0 <- sum((1 - A) * Y / (1 - ps)) / sum((1 - A) / (1 - ps))
  as.numeric(mu1 - mu0)
}

################################################################################
# 2. A small demonstration data-generating mechanism
################################################################################

# The DGM below is intentionally different from the manuscript.
# It uses only six baseline covariates and modest nonlinear terms.

generate_demo_data <- function(n = 600, seed = 2026, true_ate = 1.5) {
  set.seed(seed)

  X1 <- rnorm(n)
  X2 <- rnorm(n)
  X3 <- rbinom(n, 1, 0.50)
  X4 <- rnorm(n)
  X5 <- rbinom(n, 1, 0.40)
  X6 <- rnorm(n)

  X <- data.frame(X1, X2, X3, X4, X5, X6)

  # Treatment mechanism: nonlinear, but intentionally simple.
  lp_a <- -0.35 +
    0.55 * X1 - 0.45 * X2 + 0.60 * X3 +
    0.25 * X1 * X2 - 0.20 * X4^2
  ps_true <- expit(lp_a)
  A <- rbinom(n, 1, ps_true)

  # Outcome mechanism: constant treatment effect plus nonlinear prognostic terms.
  mu0 <- 0.40 +
    0.75 * X1 - 0.55 * X2 + 0.45 * X4 + 0.50 * X5 +
    0.25 * X1^2 - 0.30 * X2 * X3
  Y <- mu0 + true_ate * A + rnorm(n, sd = 1)

  list(
    Y = as.numeric(Y),
    A = as.numeric(A),
    X = X,
    true_ate = true_ate,
    true_ps = ps_true
  )
}

################################################################################
# 3. Candidate ANN structures
################################################################################

demo_ann_candidates <- function() {
  list(3, 5, c(4, 4))
}

demo_integration_candidates <- function() {
  list(3, 5, c(4, 4), c(5, 5))
}

# Larger candidate set used by the revised workflow.
full_integration_candidates <- function() {
  c(
    lapply(2:9, function(x) x),
    lapply(2:9, function(x) c(x, x)),
    lapply(2:6, function(x) c(x, x, x)),
    lapply(2:6, function(x) c(x, x, x, x)),
    lapply(2:6, function(x) c(x, x, x, x, x))
  )
}

################################################################################
# 4. AMORE wrappers
################################################################################

fit_amore <- function(x, y, hidden.neurons, output.layer, seed = 1) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  y <- as.numeric(y)

  set.seed(seed)
  net <- AMORE::newff(
    c(ncol(x), as.integer(hidden.neurons), 1),
    learning.rate.global = 0.001,
    momentum.global = 0.5,
    error.criterium = "LMS",
    Stao = NA,
    hidden.layer = "tansig",
    output.layer = output.layer,
    method = "ADAPTgdwm"
  )

  fit <- AMORE::train(
    net, x, y,
    error.criterium = "LMS",
    report = FALSE,
    show.step = 100,
    n.shows = 5
  )
  fit$net
}

predict_amore <- function(net, x) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  as.numeric(AMORE::sim(net, x))
}

################################################################################
# 5. Parametric nuisance models used in the educational example
################################################################################

# These are deliberately simple and not identical to manuscript specifications.

fit_parametric_ps_models <- function(A_train, X_train, X_new) {
  d1 <- data.frame(A = A_train, X_train)
  fit1 <- glm(A ~ X1 + X2 + X3 + X4 + X5 + X6,
              family = binomial(), data = d1)

  fit2 <- glm(A ~ X1 + X2 + X3,
              family = binomial(), data = d1)

  p1_tr <- clip_ps(predict(fit1, type = "response"))
  p1_ne <- clip_ps(predict(fit1, newdata = X_new, type = "response"))
  p2_tr <- clip_ps(predict(fit2, type = "response"))
  p2_ne <- clip_ps(predict(fit2, newdata = X_new, type = "response"))

  list(
    train = data.frame(ps1 = qlogis(p1_tr), ps2 = qlogis(p2_tr)),
    new   = data.frame(ps1 = qlogis(p1_ne), ps2 = qlogis(p2_ne))
  )
}

fit_parametric_or_models <- function(Y_train, A_train, X_train, X_new) {
  d1 <- data.frame(Y = Y_train, A = A_train, X_train)

  fit1 <- lm(Y ~ A + X1 + X2 + X3 + X4 + X5 + X6, data = d1)
  fit2 <- lm(Y ~ A + X1 + X2 + X4 + X5, data = d1)

  tr0 <- data.frame(A = 0, X_train)
  ne0 <- data.frame(A = 0, X_new)

  list(
    train = data.frame(
      or1 = as.numeric(predict(fit1, newdata = tr0)),
      or2 = as.numeric(predict(fit2, newdata = tr0))
    ),
    new = data.frame(
      or1 = as.numeric(predict(fit1, newdata = ne0)),
      or2 = as.numeric(predict(fit2, newdata = ne0))
    )
  )
}

################################################################################
# 6. ANN nuisance models
################################################################################

fit_ann_ps_index <- function(A_train, X_train, X_new, hidden, seed = 1) {
  net <- fit_amore(X_train, A_train, hidden, output.layer = "sigmoid", seed = seed)
  p_tr <- clip_ps(predict_amore(net, X_train))
  p_ne <- clip_ps(predict_amore(net, X_new))
  list(train = qlogis(p_tr), new = qlogis(p_ne), ps_train = p_tr)
}

fit_ann_or_index <- function(Y_train, A_train, X_train, X_new, hidden, seed = 1) {
  xtr <- data.frame(A = A_train, X_train)
  net <- fit_amore(xtr, Y_train, hidden, output.layer = "purelin", seed = seed)

  tr0 <- data.frame(A = 0, X_train)
  ne0 <- data.frame(A = 0, X_new)

  fitted_obs <- predict_amore(net, xtr)
  list(
    train = predict_amore(net, tr0),
    new = predict_amore(net, ne0),
    fitted_obs = fitted_obs
  )
}

################################################################################
# 7. Balance and PS diagnostics
################################################################################

balance_diagnostics <- function(A, X, ps, eps = 1e-7) {
  ps <- clip_ps(ps, eps)
  w <- ate_weights(A, ps, eps)
  X <- as.data.frame(X)

  sds <- vapply(X, sd, numeric(1))
  smd <- vapply(names(X), function(nm) {
    den <- sds[[nm]]
    if (!is.finite(den) || den <= 0) return(NA_real_)
    m1 <- weighted_mean_safe(X[[nm]][A == 1], w[A == 1])
    m0 <- weighted_mean_safe(X[[nm]][A == 0], w[A == 0])
    abs(m1 - m0) / den
  }, numeric(1))

  c(MASMD = mean(smd, na.rm = TRUE), MaxASMD = max(smd, na.rm = TRUE))
}

ps_diagnostics <- function(A, X, ps, eps = 1e-7) {
  ps <- clip_ps(ps, eps)
  bal <- balance_diagnostics(A, X, ps, eps)
  sw <- stabilized_weights(A, ps, eps)
  lp <- qlogis(ps)

  cal_fit <- suppressWarnings(try(glm(A ~ lp, family = binomial()), silent = TRUE))
  if (inherits(cal_fit, "try-error")) {
    cal_int <- NA_real_
    cal_slope <- NA_real_
  } else {
    cf <- coef(cal_fit)
    cal_int <- unname(cf[1])
    cal_slope <- unname(cf[2])
  }

  c(
    Brier = mean((A - ps)^2),
    MASMD = bal[["MASMD"]],
    MaxASMD = bal[["MaxASMD"]],
    CalIntercept = cal_int,
    CalSlope = cal_slope,
    ESS_total = effective_sample_size(sw),
    ESS_treated = effective_sample_size(sw[A == 1]),
    ESS_control = effective_sample_size(sw[A == 0]),
    P99_stabilized_weight = unname(quantile(sw, 0.99)),
    Max_stabilized_weight = max(sw),
    Extreme_PS_fraction = mean(ps < 0.05 | ps > 0.95),
    PS_P01 = unname(quantile(ps, 0.01)),
    PS_P99 = unname(quantile(ps, 0.99))
  )
}

################################################################################
# 8. Select ANN.PS and ANN.OcR structures in the observed dataset
################################################################################

select_ann_ps_structure <- function(A, X, candidates = demo_ann_candidates(), seed = 100) {
  ans <- lapply(seq_along(candidates), function(i) {
    h <- candidates[[i]]
    fit <- fit_ann_ps_index(A, X, X, h, seed_hash(seed, "ps", i))
    bal <- balance_diagnostics(A, X, fit$ps_train)
    data.frame(
      candidate = i,
      structure = structure_label(h),
      MASMD = bal[["MASMD"]],
      stringsAsFactors = FALSE
    )
  })
  tab <- do.call(rbind, ans)
  best <- which.min(tab$MASMD)
  list(hidden = candidates[[best]], table = tab)
}

select_ann_or_structure <- function(Y, A, X, candidates = demo_ann_candidates(), seed = 200) {
  ans <- lapply(seq_along(candidates), function(i) {
    h <- candidates[[i]]
    fit <- fit_ann_or_index(Y, A, X, X, h, seed_hash(seed, "or", i))
    mopae <- mean(abs(Y - fit$fitted_obs))
    data.frame(
      candidate = i,
      structure = structure_label(h),
      MOPAE = mopae,
      stringsAsFactors = FALSE
    )
  })
  tab <- do.call(rbind, ans)
  best <- which.min(tab$MOPAE)
  list(hidden = candidates[[best]], table = tab)
}

################################################################################
# 9. Build the six npMiPS model indexes
################################################################################

build_all_indexes <- function(Y_train, A_train, X_train, X_new,
                              h_ps_ann, h_or_ann, seed = 1) {
  pmod <- fit_parametric_ps_models(A_train, X_train, X_new)
  omod <- fit_parametric_or_models(Y_train, A_train, X_train, X_new)
  pann <- fit_ann_ps_index(A_train, X_train, X_new, h_ps_ann,
                           seed_hash(seed, "ANNPS"))
  oann <- fit_ann_or_index(Y_train, A_train, X_train, X_new, h_or_ann,
                           seed_hash(seed, "ANNOR"))

  train <- data.frame(
    ANN_PS = pann$train,
    PS_1 = pmod$train$ps1,
    PS_2 = pmod$train$ps2,
    ANN_OcR = oann$train,
    OcR_1 = omod$train$or1,
    OcR_2 = omod$train$or2
  )

  new <- data.frame(
    ANN_PS = pann$new,
    PS_1 = pmod$new$ps1,
    PS_2 = pmod$new$ps2,
    ANN_OcR = oann$new,
    OcR_1 = omod$new$or1,
    OcR_2 = omod$new$or2
  )

  list(train = train, new = new)
}

subset_indexes <- function(index_df, code = "111111") {
  if (!grepl("^[01]{6}$", code)) stop("code must be a six-digit 0/1 string")
  keep <- as.integer(strsplit(code, "", fixed = TRUE)[[1]]) == 1
  if (!any(keep)) stop("At least one index must be included.")
  index_df[, keep, drop = FALSE]
}

################################################################################
# 10. Integration ANN and OOF-MASMD selection
################################################################################

fit_integration_ann <- function(index_train, A_train, hidden, seed = 1) {
  fit_amore(index_train, A_train, hidden, output.layer = "sigmoid", seed = seed)
}

predict_integration_ann <- function(net, index_new, eps = 1e-7) {
  clip_ps(predict_amore(net, index_new), eps)
}

select_integration_oof_masmd <- function(Y, A, X,
                                         code = "111111",
                                         h_ps_ann,
                                         h_or_ann,
                                         candidate_structures = demo_integration_candidates(),
                                         K = 3,
                                         seed = 300,
                                         eps = 1e-7,
                                         tie_tol = 1e-12) {
  folds <- make_stratified_folds(A, K = K, seed = seed)

  # Cache nuisance indexes once per outer fold.
  fold_cache <- vector("list", K)
  for (k in seq_len(K)) {
    tr <- which(folds != k)
    te <- which(folds == k)
    idx <- build_all_indexes(
      Y_train = Y[tr], A_train = A[tr], X_train = X[tr, , drop = FALSE],
      X_new = X[te, , drop = FALSE],
      h_ps_ann = h_ps_ann, h_or_ann = h_or_ann,
      seed = seed_hash(seed, "nuisance", k)
    )
    fold_cache[[k]] <- list(tr = tr, te = te, index_train = idx$train, index_test = idx$new)
  }

  rows <- vector("list", length(candidate_structures))
  oof_store <- vector("list", length(candidate_structures))

  for (j in seq_along(candidate_structures)) {
    h <- candidate_structures[[j]]
    oof_ps <- rep(NA_real_, length(A))

    for (k in seq_len(K)) {
      fc <- fold_cache[[k]]
      ztr <- subset_indexes(fc$index_train, code)
      zte <- subset_indexes(fc$index_test, code)
      net <- fit_integration_ann(
        ztr, A[fc$tr], h,
        seed = seed_hash(seed, "integration", j, k)
      )
      oof_ps[fc$te] <- predict_integration_ann(net, zte, eps)
    }

    if (anyNA(oof_ps)) stop("OOF propensity predictions are incomplete.")
    diag <- ps_diagnostics(A, X, oof_ps, eps)

    rows[[j]] <- data.frame(
      candidate = j,
      structure = structure_label(h),
      MASMD = diag[["MASMD"]],
      MaxASMD = diag[["MaxASMD"]],
      Brier = diag[["Brier"]],
      ESS_total = diag[["ESS_total"]],
      Max_stabilized_weight = diag[["Max_stabilized_weight"]],
      stringsAsFactors = FALSE
    )
    oof_store[[j]] <- oof_ps
  }

  tab <- do.call(rbind, rows)
  min_m <- min(tab$MASMD)
  tied <- which(abs(tab$MASMD - min_m) <= tie_tol)

  if (length(tied) > 1) {
    max_ess <- max(tab$ESS_total[tied])
    tied <- tied[abs(tab$ESS_total[tied] - max_ess) <= tie_tol]
  }

  if (length(tied) > 1) {
    comp <- t(vapply(candidate_structures[tied], structure_complexity, numeric(2)))
    ord <- order(comp[, 1], comp[, 2])
    best <- tied[ord[1]]
  } else {
    best <- tied[1]
  }

  list(
    hidden = candidate_structures[[best]],
    selected_row = tab[best, , drop = FALSE],
    candidate_table = tab,
    oof_ps = oof_store[[best]],
    folds = folds,
    oof_diagnostics = ps_diagnostics(A, X, oof_store[[best]], eps)
  )
}

################################################################################
# 11. Final full-data npMiPS fit
################################################################################

fit_npMiPS_full <- function(Y, A, X, code,
                            h_ps_ann, h_or_ann, h_integration,
                            seed = 400, eps = 1e-7) {
  idx <- build_all_indexes(
    Y_train = Y, A_train = A, X_train = X,
    X_new = X,
    h_ps_ann = h_ps_ann, h_or_ann = h_or_ann,
    seed = seed_hash(seed, "full_nuisance")
  )

  z <- subset_indexes(idx$train, code)
  net <- fit_integration_ann(z, A, h_integration,
                             seed = seed_hash(seed, "full_integration"))
  ps <- predict_integration_ann(net, z, eps)

  list(
    ate = estimate_ipw_ate(Y, A, ps, eps),
    ps = ps,
    diagnostics = ps_diagnostics(A, X, ps, eps)
  )
}

################################################################################
# 12. Bootstrap with structures fixed after selection
################################################################################

bootstrap_npMiPS <- function(Y, A, X, code,
                             h_ps_ann, h_or_ann, h_integration,
                             B = 20, seed = 500, eps = 1e-7) {
  n <- length(Y)
  boots <- numeric(B)

  for (b in seq_len(B)) {
    set.seed(seed_hash(seed, "bootstrap", b))
    id <- sample.int(n, size = n, replace = TRUE)

    fit_b <- fit_npMiPS_full(
      Y = Y[id], A = A[id], X = X[id, , drop = FALSE],
      code = code,
      h_ps_ann = h_ps_ann,
      h_or_ann = h_or_ann,
      h_integration = h_integration,
      seed = seed_hash(seed, "fit", b),
      eps = eps
    )
    boots[b] <- fit_b$ate
  }

  boots
}

bootstrap_summary <- function(point, boots) {
  se <- sd(boots)
  z <- point / se
  data.frame(
    Estimate = point,
    BSSE = se,
    LCL = point - 1.96 * se,
    UCL = point + 1.96 * se,
    P_value = 2 * pnorm(abs(z), lower.tail = FALSE)
  )
}
