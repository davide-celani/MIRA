# Native longitudinal QTE: regression and mathematical sanity checks.
# Synthetic CmdStanMCMC-compatible mocks require no CmdStan installation,
# Stan compilation, or posterior sampling. Helpers are prefixed to avoid
# collisions with the other package test fixtures. Monte Carlo tolerances
# exceed the expected sampling error; structural checks use small simulations.

qte_expect_true <- function(condition, message) {
  expect_true(isTRUE(condition), info = message)
}

qte_expect_error <- function(expr, pattern = NULL) {
  expect_error(force(expr), regexp = pattern)
}

qte_expect_near <- function(actual, expected, tol, label) {
  expect_equal(length(actual), length(expected), info = paste(label, "length mismatch"))
  error <- max(abs(actual - expected))
  expect_true(is.finite(error) && error <= tol,
    info = sprintf("%s: max error %.6g > tolerance %.6g", label, error, tol))
}

qte_mock_fit_fixture <- function(K = 3L, G = 3L, D = 4L, likelihood = 2L,
                    X = NULL, map = NULL, L = matrix(c(.6, 0, .1, .2), 2, 2, byrow = TRUE),
                    sigma = .4, nu = 7, identical_arms = FALSE,
                    mu = NULL, offsets = NULL, treatment = NULL, cov_effect = NULL,
                    bounds = c(NA_real_, NA_real_), posterior_variation = 0) {
  S <- 12L
  time_value <- seq_len(K) * (seq_len(K) + 1) / 2 - 1
  if (is.null(X)) X <- matrix(numeric(), S, 0)
  S <- nrow(X); P <- ncol(X)
  covariate_names <- if (P) colnames(X) else character()
  if (P && is.null(covariate_names)) covariate_names <- paste0("x", seq_len(P))
  colnames(X) <- covariate_names
  if (is.null(map)) map <- data.frame(
    index = seq_len(P), name = covariate_names, label = covariate_names,
    original_name = covariate_names, original_label = covariate_names,
    type = rep("binary_numeric", P), encoding = rep("identity_binary", P),
    level = rep("1", P), reference_level = rep("0", P),
    center = rep(0, P), scale = rep(1, P), unit = rep("1 vs 0", P),
    stringsAsFactors = FALSE)
  if (is.null(mu)) mu <- .2 + .08 * time_value
  if (is.null(offsets)) offsets <- if (identical_arms) rep(0, G - 1L) else .1 * seq_len(G - 1L)
  if (is.null(treatment)) treatment <- if (identical_arms) matrix(0, G - 1L, K) else outer(seq_len(G - 1L), .07 * time_value)
  if (is.null(cov_effect)) cov_effect <- matrix(.2, P, K)
  columns <- list()
  add <- function(name, value) columns[[name]] <<- rep_len(as.numeric(value), D)
  add("baseline_mean", mu[1] + posterior_variation * seq_len(D))
  add("beta_time", .08); add("tau_common", .02); add("arm_baseline_sd", .1)
  add("sigma", sigma); add("nu_value", nu); add("nu[1]", nu)
  add("sigma_intercept", sqrt(sum(L[1, ]^2)))
  add("sigma_slope", sqrt(sum(L[2, ]^2)))
  add("rho_subject", if (prod(sqrt(rowSums(L^2))) > 0) sum(L[1, ] * L[2, ]) / prod(sqrt(rowSums(L^2))) else 0)
  add("mcid", .1); add("residual_sd", sigma); add("residual_cv", sqrt(expm1(sigma^2)))
  for (i in 1:2) for (j in 1:2) add(sprintf("subject_cholesky[%d,%d]", i, j), L[i, j])
  for (k in seq_len(K)) add(sprintf("mu_reference[%d]", k), mu[k] + posterior_variation * seq_len(D))
  for (g in seq_len(G - 1L)) {
    add(sprintf("arm_baseline_offset[%d]", g), offsets[g])
    add(sprintf("beta_treatment[%d]", g), .07 * g)
    add(sprintf("tau_treatment[%d]", g), .02)
    for (k in seq_len(K)) add(sprintf("treatment_change[%d,%d]", g, k), treatment[g, k])
  }
  if (P) for (p in seq_len(P)) {
    add(sprintf("covariate_baseline_effect[%d]", p), cov_effect[p, 1])
    add(sprintf("beta_covariate[%d]", p), .01)
    add(sprintf("tau_covariate[%d]", p), .02)
    for (k in seq_len(K)) {
      add(sprintf("covariate_effect[%d,%d]", p, k), cov_effect[p, k])
      for (prefix in c("covariate_level_difference", "covariate_change_difference", "covariate_directional_change_difference"))
        add(sprintf("%s[%d,%d]", prefix, p, k), cov_effect[p, k])
    }
  }
  clamp <- function(z) { if (is.finite(bounds[1])) z <- pmax(bounds[1], z); if (is.finite(bounds[2])) z <- pmin(bounds[2], z); z }
  location <- function(g, k) {
    eta <- mu[k] + posterior_variation * seq_len(D) + if (g == 1L) 0 else offsets[g - 1L] + treatment[g - 1L, k]
    clamp(if (likelihood == 3L) exp(eta) else eta)
  }
  for (g in seq_len(G)) for (k in seq_len(K)) {
    level <- location(g, k); change <- level - location(g, 1L)
    for (prefix in c("population_mean", "population_median")) add(sprintf("%s[%d,%d]", prefix, g, k), level)
    add(sprintf("population_model_location[%d,%d]", g, k), mu[k])
    for (prefix in c("population_change_from_baseline", "population_directional_change", "population_standardized_change", "new_subject_latent_change_draw", "new_subject_predictive_change_draw"))
      add(sprintf("%s[%d,%d]", prefix, g, k), change)
    for (prefix in c("new_subject_latent_any_improvement_draw", "new_subject_latent_responder_draw", "new_subject_predictive_responder_draw"))
      add(sprintf("%s[%d,%d]", prefix, g, k), as.numeric(change > .1))
    add(sprintf("population_ratio_from_baseline[%d,%d]", g, k), level / location(g, 1L))
    add(sprintf("population_percent_change_from_baseline[%d,%d]", g, k), 100 * (level / location(g, 1L) - 1))
  }
  for (g in seq_len(G - 1L)) for (k in seq_len(K)) {
    z <- location(g + 1L, k) - location(g + 1L, 1L) - (location(1L, k) - location(1L, 1L))
    for (prefix in c("treatment_change_difference", "directional_treatment_benefit", "treatment_responder_probability_difference")) add(sprintf("%s[%d,%d]", prefix, g, k), z)
    for (prefix in c("treatment_benefit_positive_draw", "treatment_benefit_meaningful_draw")) add(sprintf("%s[%d,%d]", prefix, g, k), as.numeric(z > .1))
    add(sprintf("treatment_ratio_of_ratios[%d,%d]", g, k), exp(z))
  }
  arm <- rep_len(seq_len(G), S)
  for (s in seq_len(S)) for (k in seq_len(K)) {
    z <- location(arm[s], k) - location(arm[s], 1L)
    for (prefix in c("individual_change_from_baseline", "individual_directional_change", "individual_change_minus_mcid")) add(sprintf("%s[%d,%d]", prefix, k, s), z)
    for (prefix in c("individual_any_improvement_draw", "individual_meaningful_responder_draw")) add(sprintf("%s[%d,%d]", prefix, k, s), as.numeric(z > .1))
  }
  N <- S * K
  y_range <- if (likelihood == 3L) c(.8, 1.8) else c(.2, .8)
  if (all(is.finite(bounds))) y_range <- bounds[1] + c(.1, .9) * diff(bounds)
  else if (is.finite(bounds[1])) y_range <- bounds[1] + c(.1, .7)
  else if (is.finite(bounds[2])) y_range <- bounds[2] - c(.7, .1)
  if (likelihood == 3L) y_range <- pmax(y_range, .01)
  y <- seq(y_range[1], y_range[2], length.out = N)
  for (n in seq_len(N)) { add(sprintf("y_rep[%d]", n), y[n]); add(sprintf("log_lik[%d]", n), -1) }
  draws <- do.call(cbind, columns)
  fit <- structure(list(
    draws = function(format = "draws_matrix") draws,
    summary = function() data.frame(variable = colnames(draws), rhat = 1, ess_bulk = 500, ess_tail = 500)
  ), class = "CmdStanMCMC")
  stan_data <- list(K = K, G = G, S = S, P = P, N = N, X = X, y = y,
    time_value = time_value, time = rep(seq_len(K), each = S), subject = rep(seq_len(S), K),
    arm = arm, arm_labels = paste0("arm_", seq_len(G)), subject_labels = as.character(seq_len(S)),
    direction = 1L, likelihood_id = likelihood, outcome_name = "Synthetic outcome", outcome = "custom",
    has_lower_bound = as.integer(is.finite(bounds[1])), outcome_lower_bound = if (is.finite(bounds[1])) bounds[1] else 0,
    has_upper_bound = as.integer(is.finite(bounds[2])), outcome_upper_bound = if (is.finite(bounds[2])) bounds[2] else 0,
    meaningful_between_arm_difference = .1, covariate_names = covariate_names,
    covariate_original_names = map$original_name, covariate_labels = map$label,
    covariate_types = map$type, covariate_reference_levels = map$reference_level,
    covariate_centers = map$center, covariate_scales = map$scale, covariate_map = map,
    covariate_metadata = list(columns = map))
  list(fit = fit, data = stan_data, draws = draws, L = L, mu = mu, offsets = offsets, treatment = treatment)
}
qte_mock_summary <- function(fx, ..., .warn_tail = FALSE) {
  defaults <- list(fit = fx$fit, stan_data = fx$data, verbose = FALSE,
    qte_probs = c(.1, .5, .9), qte_distribution = "both", qte_standardization = "reference",
    qte_covariates = FALSE, qte_n_sim = 300L, qte_max_draws = 4L, qte_seed = 917)
  additions <- list(...)
  defaults[names(additions)] <- additions
  withCallingHandlers(do.call(MIRA::mira_summary_long, defaults), warning = function(w) {
    if (!.warn_tail && grepl("^QTE extreme quantiles", conditionMessage(w))) {
      invokeRestart("muffleWarning")
    }
  })
}
qte_estimand_tables <- function(qte) list(outcome = qte$outcome, time = qte$time$temporal_quantile_shift,
  change = qte$change$quantiles, level = qte$treatment$level,
  change_qte = qte$treatment$change, did = qte$treatment$did)
qte_select_quantiles <- function(table, distribution, time = NULL, arm = 1L, standardization = "reference") {
  rows <- table$distribution == distribution & table$standardization == standardization & table$arm == arm
  if (!is.null(time)) rows <- rows & table$time == time
  selected <- table[rows, , drop = FALSE]
  selected <- selected[order(selected$tau), , drop = FALSE]
  selected$mean
}

test_that("199 default probabilities; exact sorted custom grid; argument validation", {
  fx <- qte_mock_fit_fixture(K = 2L, G = 2L)
  warnings <- character()
  q <- withCallingHandlers(qte_mock_summary(fx, qte_probs = NULL, qte_n_sim = 50L, qte_distribution = "latent", .warn_tail = TRUE),
    warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") })$qte
  qte_expect_true(any(grepl("QTE extreme quantiles", warnings)), "An undersampled tail should produce an explicit Monte Carlo diagnostic warning.")
  qte_expect_true(length(q$probs) == 199L, "Default grid must contain exactly 199 probabilities.")
  qte_expect_near(q$probs, seq(.005, .995, by = .005), 1e-12, "Default 0.5%-99.5% grid")
  custom <- c(.9, .1, .5)
  q <- qte_mock_summary(fx, qte_probs = custom)$qte
  qte_expect_true(identical(q$probs, sort(custom)), "Custom probabilities must be sorted without adding probabilities.")
  for (bad in list(numeric(), c(0, .5), c(.5, 1), c(.5, .5), c(NA, .5), Inf, "0.5", TRUE))
    qte_expect_error(qte_mock_summary(fx, qte_probs = bad))
  for (arg in c("qte_n_sim", "qte_max_draws")) for (bad in list(0, -1, 1.5, NA_real_, Inf, TRUE)) {
    options <- setNames(list(bad), arg)
    qte_expect_error(do.call(summarize, c(list(fx), options)))
  }
  for (arg in c("qte_covariates", "qte_keep_draws")) for (bad in list(NA, 1, c(TRUE, FALSE)))
    qte_expect_error(do.call(summarize, c(list(fx), setNames(list(bad), arg))))
  for (arg in c("qte_distribution", "qte_standardization", "qte_arm_contrasts"))
    qte_expect_error(do.call(summarize, c(list(fx), setNames(list("invalid"), arg))))
  for (bad in list(-1, 1.5, NA_real_, Inf, TRUE, c(1, 2))) qte_expect_error(qte_mock_summary(fx, qte_seed = bad))
  qte_expect_error(qte_mock_summary(fx, qte_probs = matrix(c(.1, .5), 1, 2)))
})

test_that("all K=2..5 and G=2..4 pairs, estimands, metadata, posterior columns", {
  mandatory <- c("estimand", "tau", "distribution", "standardization", "mean", "median", "sd", "mad", "lower", "upper", "CrI_width", "P_gt_0", "P_lt_0", "mean_directional_effect", "P_benefit")
  pair_columns <- c("from", "from_label", "from_time_value", "to", "to_label", "to_time_value", "elapsed")
  for (K in 2:5) for (G in 2:4) {
    fx <- qte_mock_fit_fixture(K = K, G = G)
    q <- qte_mock_summary(fx, qte_probs = c(.25, .5), qte_standardization = "both", qte_n_sim = 40L)$qte
    tables <- qte_estimand_tables(q)
    expected <- c(outcome = G * K, time = G * choose(K, 2), change = G * choose(K, 2),
      level = choose(G, 2) * K, change_qte = choose(G, 2) * choose(K, 2), did = choose(G, 2) * choose(K, 2)) * 2 * 2 * 2
    for (name in names(tables)) {
      table <- tables[[name]]
      qte_expect_true(nrow(table) == expected[[name]], paste("Unexpected rows", name, "K", K, "G", G, nrow(table), expected[[name]]))
      qte_expect_true(all(mandatory %in% names(table)), paste("Missing common columns in", name))
      qte_expect_true(length(unique(table$estimand)) == 1L, paste("Ambiguous estimand in", name))
      qte_expect_true(setequal(unique(table$distribution), c("latent", "predictive")), "Both distribution labels are needed.")
      qte_expect_near(table$percentile, 100 * table$tau, 0, "percentile metadata")
      if (name %in% c("time", "change", "change_qte", "did")) {
        qte_expect_true(all(pair_columns %in% names(table)), paste("Missing pair metadata in", name))
        pairs <- unique(table[c("from", "to")])
        qte_expect_true(nrow(pairs) == choose(K, 2) && all(pairs$to > pairs$from), "All unique forward temporal pairs are required.")
        qte_expect_near(table$elapsed, fx$data$time_value[table$to] - fx$data$time_value[table$from], 0, "elapsed")
      }
    }
    ref <- qte_mock_summary(fx, qte_probs = .5, qte_arm_contrasts = "reference", qte_distribution = "latent", qte_n_sim = 30L)$qte
    qte_expect_true(nrow(ref$treatment$level) == (G - 1L) * K, "Reference-only arm comparison count is wrong.")
  }
})

test_that("Gaussian and lognormal analytical outcome quantiles, including median", {
  probs <- c(.1, .5, .9)
  for (family in c(2L, 3L)) {
    fx <- qte_mock_fit_fixture(K = 2L, G = 2L, D = 4L, likelihood = family)
    q <- qte_mock_summary(fx, qte_probs = probs, qte_n_sim = 40000L)$qte
    covariance <- fx$L %*% t(fx$L)
    for (distribution in c("latent", "predictive")) for (k in 1:2) {
      t <- fx$data$time_value[k] - fx$data$time_value[1]
      variance <- covariance[1, 1] + 2 * t * covariance[1, 2] + t^2 * covariance[2, 2] + if (distribution == "predictive") .4^2 else 0
      expected <- fx$mu[k] + qnorm(probs) * sqrt(variance)
      if (family == 3L) expected <- exp(expected)
      actual <- qte_select_quantiles(q$outcome, distribution, time = k)
      qte_expect_near(actual, expected, if (family == 3L) .045 * max(expected) else .04, paste("Analytic outcome", family, distribution, k))
      median <- actual[2L]
      qte_expect_near(median, if (family == 3L) exp(fx$mu[k]) else fx$mu[k], .025, "Reference population median")
    }
  }
})

test_that("Gaussian within-subject changes reuse random effects and distinguish C from B", {
  fx <- qte_mock_fit_fixture(K = 3L, G = 2L, L = matrix(c(.8, 0, .12, .18), 2, 2, byrow = TRUE))
  q <- qte_mock_summary(fx, qte_n_sim = 40000L)$qte
  pairs <- unique(q$change$quantiles[c("from", "to")])
  for (distribution in c("latent", "predictive")) for (i in seq_len(nrow(pairs))) {
    from <- pairs$from[i]; to <- pairs$to[i]
    dt <- fx$data$time_value[to] - fx$data$time_value[from]
    sd_change <- sqrt(dt^2 * sum(fx$L[2, ]^2) + if (distribution == "predictive") 2 * .4^2 else 0)
    expected <- fx$mu[to] - fx$mu[from] + qnorm(c(.1, .5, .9)) * sd_change
    rows <- q$change$quantiles$from == from & q$change$quantiles$to == to
    actual <- qte_select_quantiles(q$change$quantiles[rows, ], distribution)
    qte_expect_near(actual, expected, .04, paste("Analytic individual changes", distribution, from, to))
  }
  rows <- function(table) table$arm == 1L & table$from == 1L & table$to == 2L & table$distribution == "predictive" & table$tau == .9
  B <- q$time$temporal_quantile_shift$mean[rows(q$time$temporal_quantile_shift)]
  C <- q$change$quantiles$mean[rows(q$change$quantiles)]
  qte_expect_true(abs(C - B) > .1, "The quantile of paired changes C must differ from the shift of marginal quantiles B.")
})

test_that("Student-t likelihood and fallback subject Cholesky", {
  fx <- qte_mock_fit_fixture(K = 2L, G = 2L, likelihood = 1L, L = matrix(0, 2, 2), sigma = .7, nu = 5)
  q <- qte_mock_summary(fx, qte_distribution = "predictive", qte_n_sim = 50000L)$qte
  expected <- fx$mu[1] + .7 * qt(c(.1, .5, .9), 5)
  qte_expect_near(qte_select_quantiles(q$outcome, "predictive", time = 1L), expected, .035, "Student-t sigma parameterization")
  fx <- qte_mock_fit_fixture(K = 2L, G = 2L)
  direct <- qte_mock_summary(fx, qte_n_sim = 500L)$qte
  covariance <- fx$L %*% t(fx$L)
  sds <- sqrt(diag(covariance)); corr_L <- diag(1 / sds) %*% fx$L
  alternate <- fx$draws[, !grepl("^subject_cholesky\\[", colnames(fx$draws)), drop = FALSE]
  for (i in 1:2) alternate <- cbind(alternate, setNames(data.frame(rep(sds[i], nrow(alternate))), sprintf("sigma_subject[%d]", i)))
  for (i in 1:2) for (j in 1:2) alternate <- cbind(alternate, setNames(data.frame(rep(corr_L[i, j], nrow(alternate))), sprintf("L_subject[%d,%d]", i, j)))
  alternate <- as.matrix(alternate)
  fx$fit$draws <- function(format = "draws_matrix") alternate
  reconstructed <- qte_mock_summary(fx, qte_n_sim = 500L)$qte
  qte_expect_near(reconstructed$outcome$mean, direct$outcome$mean, 1e-12, "Cholesky fallback reconstruction")
})

test_that("four bound configurations preserve floor/ceiling point masses", {
  for (family in c(1L, 2L, 3L)) for (bounds in list(c(NA, NA), c(.5, NA), c(NA, .6), c(.5, .6))) {
    fx <- qte_mock_fit_fixture(K = 2L, G = 2L, likelihood = family, bounds = bounds)
    q <- qte_mock_summary(fx, qte_probs = c(.005, .01, .99, .995), qte_n_sim = 2000L)$qte
    values <- q$outcome$mean
    qte_expect_true(all(is.finite(values)), "Bounded quantiles must remain finite.")
    if (is.finite(bounds[1])) {
      qte_expect_true(all(values >= bounds[1]), "Lower-bound violation.")
      qte_expect_true(any(values == bounds[1]), "The lower endpoint point mass was lost.")
    }
    if (is.finite(bounds[2])) {
      qte_expect_true(all(values <= bounds[2]), "Upper-bound violation.")
      qte_expect_true(any(values == bounds[2]), "The upper endpoint point mass was lost.")
    }
    if (family == 3L && !is.finite(bounds[1])) qte_expect_true(all(values > 0), "Lognormal natural-scale quantiles must be positive.")
  }
})

test_that("identical arms give exact null effects with common random numbers", {
  for (family in 1:3) {
    fx <- qte_mock_fit_fixture(K = 3L, G = 4L, likelihood = family, identical_arms = TRUE, bounds = c(.1, 3))
    q <- qte_mock_summary(fx, qte_n_sim = 600L, qte_standardization = "both")$qte
    for (table in q$treatment) if (is.data.frame(table) && nrow(table))
      qte_expect_true(all(table$mean == 0 & table$lower == 0 & table$upper == 0), "Identical arms must produce exactly null QTE, Change-QTE and DiD.")
  }
})

test_that("lognormal Change-QTE and quantile DiD remain distinct", {
  fx <- qte_mock_fit_fixture(K = 2L, G = 2L, likelihood = 3L,
    L = matrix(c(.25, 0, 0, .15), 2, 2, byrow = TRUE), sigma = .3,
    mu = c(0, .3), offsets = .5, treatment = matrix(c(0, .2), 1, 2))
  q <- qte_mock_summary(fx, qte_n_sim = 50000L, qte_probs = c(.1, .5, .9))$qte
  change <- q$treatment$change; did <- q$treatment$did
  select <- function(table) table$mean[table$distribution == "predictive" & table$tau == .9]
  qte_expect_true(abs(select(change) - select(did)) > .05, "Natural-scale Change-QTE and quantile DiD were conflated.")
  qte_expect_true(!identical(unique(change$estimand), unique(did$estimand)), "Different estimands need different labels.")
})

test_that("empirical standardization samples whole rows and forms quantiles after mixing", {
  X <- cbind(x1 = rep(c(0, 1), 20), x2 = rep(c(0, 1), 20))
  fx <- qte_mock_fit_fixture(K = 2L, G = 2L, X = X, L = matrix(0, 2, 2), sigma = 1e-8,
    cov_effect = matrix(c(3, -3), 2, 2))
  q <- qte_mock_summary(fx, qte_standardization = "both", qte_distribution = "latent", qte_n_sim = 4000L)$qte
  reference <- q$outcome[q$outcome$standardization == "reference", ]
  empirical <- q$outcome[q$outcome$standardization == "empirical", ]
  qte_expect_near(empirical$mean, reference$mean, 1e-12, "Complete empirical rows preserve perfect dependence")
  X <- matrix(rep(c(0, 1), 20), ncol = 1, dimnames = list(NULL, "binary"))
  fx <- qte_mock_fit_fixture(K = 2L, G = 2L, X = X, L = matrix(0, 2, 2), sigma = 1e-8,
    mu = c(0, 0), cov_effect = matrix(10, 1, 2))
  q <- qte_mock_summary(fx, qte_probs = c(.25, .75), qte_distribution = "latent",
    qte_standardization = "empirical", qte_n_sim = 6000L)$qte
  qte_expect_near(qte_select_quantiles(q$outcome, "latent", time = 1L, standardization = "empirical"), c(0, 10), 1e-12,
    "Mixture quantiles cannot be averages of subgroup quantiles")
})

test_that("binary, continuous and categorical profiles reuse original encoding", {
  X <- cbind(binary = rep(c(0, 1), 12), age = seq(-2, 2, length.out = 24),
    categoryB = rep(c(0, 1, 0), 8), categoryC = rep(c(0, 0, 1), 8))
  map <- data.frame(index = 1:4, name = colnames(X), label = colnames(X),
    original_name = c("binary", "age", "category", "category"), original_label = c("binary", "age", "category", "category"),
    type = c("binary_numeric", "numeric", "factor", "factor"),
    encoding = c("identity_binary", "center_scale", "treatment", "treatment"),
    level = c("1", NA, "B", "C"), reference_level = c("0", NA, "A", "A"),
    center = c(0, 50, 0, 0), scale = c(1, 10, 1, 1), unit = "unit", stringsAsFactors = FALSE)
  fx <- qte_mock_fit_fixture(K = 3L, G = 3L, X = X, map = map)
  q <- qte_mock_summary(fx, qte_covariates = TRUE, qte_n_sim = 100L)$qte
  profiles <- q$covariates$profiles
  original_column <- intersect(c("original_covariate", "covariate", "original_name"), names(profiles))[1]
  qte_expect_true(!is.na(original_column), "Profiles need original covariate names.")
  counts <- table(profiles[[original_column]])
  qte_expect_true(all(counts[c("binary", "age", "category")] == c(2, 3, 3)), "Profiles must group categorical dummy columns and analyze one covariate at a time.")
  age <- profiles[profiles[[original_column]] == "age", , drop = FALSE]
  qte_expect_true(all(c("original_value", "encoded_value", "center", "scale") %in% names(age)), "Continuous profile metadata must retain original-scale values and encoding.")
  qte_expect_near(sort(age$original_value), 50 + 10 * quantile(X[, "age"], c(.25, .5, .75), names = FALSE), 1e-12, "Original-scale age quartiles")
  qte_expect_near(age$encoded_value, (age$original_value - age$center) / age$scale, 1e-12, "Continuous profile encoding")
  qte_expect_true(isFALSE(q$metadata$explicit_treatment_covariate_interaction), "QTE metadata must state absence of explicit treatment-covariate interactions.")
  encoded_profiles <- q$covariates$encoded_profiles
  qte_expect_true(is.matrix(encoded_profiles) && identical(dim(encoded_profiles), c(8L, 4L)), "Encoded profile matrix is missing or has unexpected dimensions.")
  qte_expect_true(all(rowSums(encoded_profiles[, c("categoryB", "categoryC"), drop = FALSE]) <= 1), "Categorical profiles must never activate two levels of the same variable.")
  qte_expect_true(all(rowSums(encoded_profiles != 0) <= 1), "Profiles must vary only one original covariate at a time.")
  for (table in qte_estimand_tables(q$covariates)) qte_expect_true(is.data.frame(table) && nrow(table) > 0, "All six profile-specific estimand families are required.")
  qte_expect_true(is.data.frame(q$covariates$covariate_quantile_contrast) && nrow(q$covariates$covariate_quantile_contrast) > 0, "Covariate level contrasts are missing.")
  qte_expect_true(is.data.frame(q$covariates$covariate_change_quantile_contrast) && nrow(q$covariates$covariate_change_quantile_contrast) > 0, "Covariate change contrasts are missing.")
  contrasts <- q$covariates$covariate_quantile_contrast
  qte_expect_near(contrasts$mean, .2 * (contrasts$profile_b_encoded_value - contrasts$profile_a_encoded_value), 1e-12,
    "Identity-link covariate contrasts reuse the fitted encoding")
  qte_expect_near(q$covariates$covariate_change_quantile_contrast$mean, rep(0, nrow(q$covariates$covariate_change_quantile_contrast)), 1e-12,
    "Time-constant identity covariate effects cancel in changes")
  profile_treatment <- q$covariates$treatment$level
  key <- function(table) paste(table$distribution, table$tau, table$arm_a, table$arm_b, table$time)
  qte_expect_near(profile_treatment$mean, q$treatment$level$mean[match(key(profile_treatment), key(q$treatment$level))], 1e-12,
    "Additive identity model has no invented treatment-covariate interaction")
  # Legacy serialized metadata can retain original variables but omit the
  # columns table. The summary's synthetic fallback map must not split dummies.
  legacy <- fx
  legacy$data[c("covariate_map", "covariate_original_names", "covariate_labels",
    "covariate_types", "covariate_reference_levels", "covariate_centers", "covariate_scales")] <- NULL
  legacy$data$covariate_metadata <- list(variables = data.frame(
    name = c("binary", "age", "category"), label = c("binary", "age", "category"),
    type = c("binary_numeric", "numeric", "factor"),
    encoding = c("identity_binary", "center_scale", "treatment"),
    reference_level = c("0", NA, "A"), center = c(0, 50, 0), scale = c(1, 10, 1),
    column_start = c(1L, 2L, 3L), column_end = c(1L, 2L, 4L),
    levels = I(list(c("0", "1"), character(), c("A", "B", "C")))))
  # Actual X must take precedence over stale marginal quartile metadata.
  legacy$data$covariate_distributions <- data.frame(name = colnames(X),
    q25 = 100, median = 100, q75 = 100)
  recovered <- qte_mock_summary(legacy, qte_covariates = TRUE, qte_n_sim = 40L)$qte$covariates$profiles
  qte_expect_true(identical(table(recovered$original_covariate), table(profiles$original_covariate)),
    "Variables-only metadata must recover original categorical grouping through the public wrapper.")
  qte_expect_near(recovered$original_value[recovered$original_covariate == "age"], age$original_value, 1e-12,
    "Variables-only center/scale and observed X take precedence over stale marginals")
  qte_expect_true(all(recovered$metadata_source == "covariate_metadata_variables"),
    "Recovered profiles must identify their real metadata source.")
  no_X <- fx; no_X$data$X <- NULL
  missing_X <- suppressWarnings(qte_mock_summary(no_X, qte_standardization = "empirical"))$qte
  qte_expect_true("empirical" %in% missing_X$metadata$skipped_standardization && nrow(missing_X$outcome) == 0L,
    "Without joint X, empirical standardization must be explicitly unavailable rather than fabricated.")
  qte_expect_true(missing_X$metadata$status == "unavailable" &&
    missing_X$monte_carlo$posterior_draws_used == 0L &&
    length(missing_X$monte_carlo$draw_indices) == 0L,
    "When no requested estimand can be computed, metadata must not claim posterior draws were used.")
})

test_that("reproducibility, seed restoration, thinning, posterior draws and direction", {
  fx <- qte_mock_fit_fixture(K = 2L, G = 2L, D = 9L, posterior_variation = .02)
  original_rng_kind <- RNGkind()
  RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
  set.seed(871)
  before <- .Random.seed
  a <- qte_mock_summary(fx, qte_keep_draws = TRUE, qte_max_draws = 4L)$qte
  qte_expect_true(identical(before, .Random.seed), "QTE generation must restore the caller's RNG state.")
  qte_expect_true(identical(RNGkind(), c("L'Ecuyer-CMRG", "Inversion", "Rejection")), "QTE generation must restore the caller's RNG algorithm.")
  b <- qte_mock_summary(fx, qte_keep_draws = TRUE, qte_max_draws = 4L)$qte
  qte_expect_true(identical(a, b), "The same seed and options must reproduce QTE output exactly.")
  qte_expect_true(length(a$monte_carlo$draw_indices) == 4L && !anyDuplicated(a$monte_carlo$draw_indices), "Thinning must retain four unique posterior indices.")
  qte_expect_true(!is.null(a$draws), "Retained QTE draws are missing.")
  draw_values <- a$draws$strata$reference$outcome[, 2L, 1L, 1L, 1L]
  one_row <- a$outcome[a$outcome$distribution == "latent" & a$outcome$tau == .5 & a$outcome$arm == 1L & a$outcome$time == 1L, ]
  qte_expect_near(one_row$mean, mean(draw_values), 0, "Posterior mean summarizes draw-specific quantiles")
  qte_expect_near(c(one_row$lower, one_row$upper), quantile(draw_values, c(.05, .95), names = FALSE), 1e-12,
    "Credible intervals summarize posterior quantile draws")
  qte_expect_near(c(one_row$sd, one_row$mad), c(sd(draw_values), mad(draw_values)), 1e-12, "Posterior SD and MAD")
  qte_expect_true(all(a$outcome$lower <= a$outcome$median & a$outcome$median <= a$outcome$upper), "Posterior interval ordering is invalid.")
  c <- qte_mock_summary(fx, qte_seed = 918)$qte
  qte_expect_true(!identical(a$outcome$mean, c$outcome$mean), "Different seeds should alter finite-Monte-Carlo quantiles.")
  fx$data$direction <- -1L
  negative <- qte_mock_summary(fx)$qte
  fx$data$direction <- 1L
  positive <- qte_mock_summary(fx)$qte
  qte_expect_near(negative$treatment$level$mean, positive$treatment$level$mean, 0, "Raw effect cannot change with clinical direction")
  qte_expect_near(negative$treatment$level$mean_directional_effect, -positive$treatment$level$mean_directional_effect, 0, "Directional effects")
  qte_expect_true(all(is.na(positive$treatment$level$P_benefit_meaningful)) && all(is.na(positive$treatment$did$P_benefit_meaningful)),
    "A change threshold must not be silently applied to level QTE or quantile DiD.")
  qte_expect_true(all(c("P_benefit_meaningful", "meaningful_between_arm_difference", "threshold_estimand") %in% names(positive$treatment$change)),
    "Change-QTE threshold application must be explicit.")
  qte_expect_true(all(positive$treatment$change$threshold_estimand == "treatment_change_qte"), "Threshold estimand is mislabeled.")
  # When no RNG existed initially, do not leave a newly created RNG behind.
  rm(".Random.seed", envir = .GlobalEnv)
  invisible(qte_mock_summary(fx))
  qte_expect_true(!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE), "The helper introduced global RNG state when none existed before.")
  do.call(RNGkind, as.list(original_rng_kind))
})

test_that("singleton posterior/probability dimensions and exact clinical threshold equality", {
  fx <- qte_mock_fit_fixture(K = 2L, G = 2L, D = 1L, L = matrix(0, 2, 2), sigma = 0,
    mu = c(0, 0), offsets = 0, treatment = matrix(c(0, .1), 1, 2))
  q <- qte_mock_summary(fx, qte_probs = .5, qte_max_draws = 1L, qte_keep_draws = TRUE)$qte
  for (table in qte_estimand_tables(q)) qte_expect_true(nrow(table) > 0 && all(is.finite(table$mean)), "Singleton array dimensions collapsed unexpectedly.")
  qte_expect_true(identical(dim(q$draws$strata$reference$outcome), c(1L, 1L, 2L, 2L, 2L)), "Singleton retained outcome dimensions are wrong.")
  qte_expect_true(identical(dim(q$draws$strata$reference$change), c(1L, 1L, 2L, 1L, 2L)), "Singleton retained change dimensions are wrong.")
  qte_expect_true(all(q$treatment$change$mean == fx$data$meaningful_between_arm_difference), "Fixture must produce exact threshold equality.")
  qte_expect_true(all(q$treatment$change$P_benefit_meaningful == 1), "Clinically meaningful threshold equality must use the Stan >= convention.")
})

test_that("integrated summary retains legacy components and positional APIs; compact print", {
  fx <- qte_mock_fit_fixture(K = 3L, G = 3L)
  result <- qte_mock_summary(fx, qte_probs = NULL, qte_n_sim = 40L)
  qte_expect_true(inherits(result, "mira_summary_long"), "Summary class changed.")
  essential <- c("population", "change", "treatment", "covariates", "individual_change", "responders", "heterogeneity", "diagnostics", "ppc", "loo", "draws", "qte")
  qte_expect_true(all(essential %in% names(result)), "A legacy component was removed.")
  qte_expect_true(identical(colnames(result$draws), colnames(fx$draws)), "Raw posterior draw columns changed.")
  qte_expect_near(as.numeric(as.matrix(result$draws)), as.numeric(fx$draws), 0, "Raw posterior values")
  positional <- suppressWarnings(MIRA::mira_summary_long(fx$fit, fx$data, NULL, NULL, .9, c(.5, .8, .95), FALSE,
    qte_n_sim = 40L, qte_max_draws = 4L))
  qte_expect_true(inherits(positional, "mira_summary_long") && length(positional$qte$probs) == 199L, "The original seven-position signature is broken.")
  # Invoke the sourced print method explicitly, because source() uses a private environment.
  text <- capture.output(getS3method("print", "mira_summary_long", envir = asNamespace("MIRA"))(result))
  qte_expect_true(any(grepl("QTE|QUANTILE", text)), "Print report does not advertise native QTE.")
  qte_expect_true(sum(grepl("0\\.005|0\\.995", text)) < 30L, "Default print appears to enumerate the complete 199-quantile grid.")
  qte_expect_true(length(text) < 500L, "Default QTE print must remain compact.")
  print_positions <- c("x", "digits", "max_rows", "trajectories", "treatment", "covariates", "responders", "heterogeneity", "ppc", "diagnostics")
  qte_expect_true(identical(names(formals(getS3method("print", "mira_summary_long", envir = asNamespace("MIRA"))))[seq_along(print_positions)], print_positions), "The original positional print switches changed.")
  legacy_switches <- capture.output(getS3method("print", "mira_summary_long", envir = asNamespace("MIRA"))(result, 3, 4, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, qte = FALSE))
  qte_expect_true(!any(grepl("TREATMENT EFFECTS VS REFERENCE ARM|SELECTED COVARIATE EFFECTS|CLINICAL AND RESPONDER SUMMARY|HETEROGENEITY / VARIANCE COMPONENTS|POSTERIOR PREDICTIVE CHECKS", legacy_switches)), "Legacy positional print switches no longer disable their sections.")
})

test_that("fit-only metadata fallback, unavailable empirical standardization", {
  X <- matrix(rep(c(0, 1), 6), ncol = 1, dimnames = list(NULL, "binary"))
  fx <- qte_mock_fit_fixture(K = 3L, G = 3L, X = X)
  fit_info <- fx$data
  fit_info$X <- fit_info$y <- NULL
  fit_info$time_value <- c(0, 2, 5)
  fit_info$time_labels <- c("Baseline", "Week 2", "Week 5")
  fit_info$arm_labels <- c("Control", "Treatment A", "Treatment B")
  fit_info$direction <- -1L
  fit_info$meaningful_between_arm_difference <- .2
  attr(fx$fit, "mira_fit_info") <- fit_info
  warnings <- character()
  result <- withCallingHandlers(MIRA::mira_summary_long(fx$fit, verbose = FALSE,
    qte_probs = .5, qte_standardization = "both", qte_covariates = TRUE,
    qte_n_sim = 100L, qte_max_draws = 4L), warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
    })
  q <- result$qte
  qte_expect_true(q$metadata$status == "partial" && "empirical" %in% q$metadata$skipped_standardization,
    "Fit-only metadata must explicitly skip empirical standardization when joint X is unavailable.")
  qte_expect_true(any(grepl("stan_data\\$X", warnings)), "Missing joint X must produce an informative warning.")
  qte_expect_true(identical(q$metadata$time_value, fit_info$time_value), "Actual irregular time metadata were lost.")
  qte_expect_true(identical(q$metadata$time_labels, fit_info$time_labels), "Fit-only time labels were lost.")
  qte_expect_true(identical(q$metadata$arm_labels, fit_info$arm_labels), "Fit-only arm labels were lost.")
  qte_expect_true(identical(q$metadata$direction, -1L), "Fit-only clinical direction was lost.")
  qte_expect_true(all(q$treatment$change$meaningful_between_arm_difference == .2), "Fit-only threshold metadata were lost.")
  qte_expect_true(setequal(q$outcome$time_value, fit_info$time_value), "QTE tables use invented equally spaced times.")
  qte_expect_true(setequal(q$outcome$time_label, fit_info$time_labels), "QTE tables must retain time labels.")
  qte_expect_true(setequal(q$outcome$arm_label, fit_info$arm_labels), "QTE tables must retain arm labels.")
  qte_expect_true(identical(unique(q$outcome$standardization), "reference"), "An empirical distribution was fabricated without X.")
  qte_expect_near(q$treatment$level$mean_directional_effect, -q$treatment$level$mean, 0, "Direction from fit-only metadata")
  qte_expect_true(nrow(q$covariates$profiles) == 2L, "Fit-only map should still support binary profiles.")
})

test_that("missing generative posterior variables preserve legacy summary and RNG", {
  fx <- qte_mock_fit_fixture(K = 3L, G = 3L)
  incomplete <- fx$draws[, !grepl("^mu_reference\\[", colnames(fx$draws)), drop = FALSE]
  fx$fit$draws <- function(format = "draws_matrix") incomplete
  set.seed(335)
  rng_before <- .Random.seed
  warnings <- character()
  result <- withCallingHandlers(qte_mock_summary(fx), warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
  })
  qte_expect_true(result$qte$metadata$status == "unavailable", "An incomplete older fit must explicitly report unavailable QTE.")
  qte_expect_true(any(grepl("Missing QTE posterior variable.*mu_reference", warnings)), "Missing source parameter warning is not informative.")
  qte_expect_true(identical(.Random.seed, rng_before), "Unavailable QTE must restore caller RNG state.")
  qte_expect_true(inherits(result, "mira_summary_long") && nrow(result$population_time_means) == 9L && nrow(result$change) == 9L,
    "QTE unavailability must preserve existing population/change components.")
})

