# The complete-profile oracles use the multivariate-normal REML solution:
# cell means are arithmetic means and Sigma is the pooled subject covariance.
# They do not call mmrm, emmeans or a MIRA helper for the expected numbers.
mmrm_oracle_matrix <- function() {
  rbind(c(2, 5, 7), c(4, 6, 12), c(7, 8, 9), c(3, 7, 13),
    c(9, 11, 10), c(6, 9, 16), c(8, 14, 18), c(12, 13, 15))
}

mmrm_oracle_long <- function(y = mmrm_oracle_matrix(), arms = NULL) {
  n <- nrow(y)
  data <- data.frame(patient_factor = factor(rep(sprintf("S%02d", seq_len(n)), each = 3L)),
    time_factor = factor(rep(c("Zbase", "Amid", "Mfinal"), n),
      levels = c("Zbase", "Amid", "Mfinal")),
    time_index = rep(1:3, n), value = as.vector(t(y)))
  if (!is.null(arms)) data$arm_factor <- factor(rep(arms, each = 3L),
    levels = c("Zcontrol", "Aactive"))
  data
}

mmrm_oracle_keys <- function(data) {
  paste(as.character(data$patient_factor), as.character(data$time_factor), sep = "::")
}

test_that("unstructured MMRM matches subject-level covariance and exact paired inference", {
  withr::local_preserve_seed()
  set.seed(2026102601)
  seed_before <- .Random.seed
  kind_before <- RNGkind()
  y <- mmrm_oracle_matrix()
  data <- mmrm_oracle_long(y)
  before <- data
  result <- .mira_run_mmrm_long(data, value ~ time_factor, FALSE, 0.9, "holm")
  expect_true(result$performed)
  expect_true(result$converged)
  expect_s3_class(result$model, "mmrm")
  expect_null(result$error)
  expect_length(result$errors, 0L)
  expect_null(result$reason_skipped)
  expect_identical(result$ic_penalty_parameters, 6L)
  expect_identical(result$bic_sample_size, 8L)
  expect_equal(result$AIC, -2 * result$logLik + 2 * 6, tolerance = 1e-12)
  expect_equal(result$BIC, -2 * result$logLik + log(8) * 6, tolerance = 1e-12)
  expect_match(result$criteria_note, "not harmonized", fixed = TRUE)
  expect_identical(result$method, "REML")
  expect_identical(result$df_method, "Satterthwaite")
  expect_identical(result$start_method, "empirical")
  expect_identical(result$n, 24L)
  expect_identical(result$n_subjects, 8L)
  expect_identical(result$covariance_parameters, 6L)
  expect_identical(result$covariance_rank, 3L)
  expect_identical(mmrm_oracle_keys(result$data), mmrm_oracle_keys(data))
  expect_identical(data, before)
  expect_identical(.Random.seed, seed_before)
  expect_identical(RNGkind(), kind_before)
  means <- colSums(y) / 8
  centered <- sweep(y, 2L, means, "-")
  covariance <- crossprod(centered) / 7
  dimnames(covariance) <- list(levels(data$time_factor), levels(data$time_factor))
  # Empirical initialization starts at the analytical US solution here.
  # 2e-6 allows only optimizer/Hessian arithmetic at this fixture's scale.
  expect_equal(result$covariance, covariance, tolerance = 2e-6)
  basis <- rbind(c(1, 0, 0), c(-1, 1, 0), c(-1, 0, 1))
  expected <- as.vector(basis %*% means)
  se <- sqrt(diag(basis %*% (covariance / 8) %*% t(basis)))
  coefficients <- result$coefficients
  expect_identical(coefficients$term,
    c("(Intercept)", "time_factorAmid", "time_factorMfinal"))
  expect_equal(coefficients$estimate, expected, tolerance = 1e-10)
  expect_equal(coefficients$SE, se, tolerance = 2e-6)
  expect_equal(coefficients$df, rep(7, 3L), tolerance = 2e-5)
  expect_equal(coefficients$t.ratio, expected / se, tolerance = 2e-6)
  expect_equal(coefficients$p_raw, 2 * stats::pt(-abs(expected / se), 7), tolerance = 2e-6)
  critical <- stats::qt(0.95, 7)
  expect_equal(coefficients$lower.CL, expected - critical * se, tolerance = 2e-6)
  expect_equal(coefficients$upper.CL, expected + critical * se, tolerance = 2e-6)
  expect_true(all(coefficients$estimable))
  # Counting 24 repeated rows as independent would give df=21 and p<0.001.
  # The independent subject-level contrast has df=7 and p about 0.002735.
  final <- result$emmeans$baseline_final$tidy
  delta <- y[, 3L] - y[, 1L]
  center <- sum(delta) / 8
  delta_se <- sqrt(sum((delta - center)^2) / 7) / sqrt(8)
  expect_identical(nrow(final), 1L)
  expect_identical(final$contrast, "final - baseline")
  expect_equal(final$estimate, 6.125, tolerance = 1e-10)
  expect_equal(final$SE, delta_se, tolerance = 2e-6)
  expect_equal(final$df, 7, tolerance = 2e-5)
  expect_equal(final$p_raw, 2 * stats::pt(-abs(center / delta_se), 7), tolerance = 2e-6)
  expect_gt(final$p_raw, 0.002)
  expect_equal(result$AIC, stats::AIC(result$model), tolerance = 1e-12)
  expect_equal(result$BIC, stats::BIC(result$model), tolerance = 1e-12)
  expect_equal(result$logLik, as.numeric(stats::logLik(result$model)), tolerance = 1e-12)
})

test_that("unequal MMRM arms retain cell identity and pooled subject uncertainty", {
  y <- mmrm_oracle_matrix()
  arms <- c(rep("Zcontrol", 3L), rep("Aactive", 5L))
  data <- mmrm_oracle_long(y, arms)
  result <- .mira_run_mmrm_long(data, value ~ time_factor * arm_factor, TRUE, 0.9, "holm")
  expect_true(result$performed)
  groups <- list(Zcontrol = 1:3, Aactive = 4:8)
  means <- vapply(groups, function(index) colSums(y[index, , drop = FALSE]) / length(index), numeric(3L))
  residuals <- y
  for (arm in names(groups)) residuals[groups[[arm]], ] <-
    sweep(y[groups[[arm]], , drop = FALSE], 2L, means[, arm], "-")
  covariance <- crossprod(residuals) / 6
  dimnames(covariance) <- list(levels(data$time_factor), levels(data$time_factor))
  expect_equal(result$covariance, covariance, tolerance = 2e-6)
  cells <- result$emmeans$arm_time$tidy
  keys <- paste(as.character(cells$arm_factor), as.character(cells$time_factor), sep = "::")
  expected <- expand.grid(arm_factor = names(groups),
    time_factor = levels(data$time_factor), stringsAsFactors = FALSE)
  expected_keys <- paste(expected$arm_factor, expected$time_factor, sep = "::")
  expect_identical(nrow(cells), 6L)
  expect_identical(anyDuplicated(keys), 0L)
  expect_identical(sort(keys), sort(expected_keys))
  index <- match(expected_keys, keys)
  positions <- cbind(match(expected$time_factor, levels(data$time_factor)),
    match(expected$arm_factor, names(groups)))
  counts <- c(Zcontrol = 3, Aactive = 5)[expected$arm_factor]
  expected_mean <- means[positions]
  expected_se <- sqrt(diag(covariance)[positions[, 1L]] / counts)
  expect_identical(keys[index], expected_keys)
  expect_equal(cells$emmean[index], unname(expected_mean), tolerance = 1e-10)
  expect_equal(cells$SE[index], unname(expected_se), tolerance = 2e-6)
  expect_equal(cells$df[index], rep(6, 6L), tolerance = 3e-5)
  expect_equal(cells$lower.CL[index], unname(expected_mean - stats::qt(0.95, 6) * expected_se), tolerance = 2e-6)
  expect_equal(cells$upper.CL[index], unname(expected_mean + stats::qt(0.95, 6) * expected_se), tolerance = 2e-6)
  simple <- result$emmeans$arm_time$simple_arm_tidy
  expect_identical(nrow(simple), 3L)
  expect_identical(as.character(simple$time_factor), levels(data$time_factor))
  expect_identical(as.character(simple$contrast), rep("Zcontrol - Aactive", 3L))
  effect <- means[, "Zcontrol"] - means[, "Aactive"]
  se <- sqrt(diag(covariance) * (1 / 3 + 1 / 5))
  p <- 2 * stats::pt(-abs(effect / se), 6)
  expect_equal(simple$estimate, unname(effect), tolerance = 1e-10)
  expect_equal(simple$SE, unname(se), tolerance = 2e-6)
  expect_equal(simple$df, rep(6, 3L), tolerance = 3e-5)
  expect_equal(simple$p_raw, unname(p), tolerance = 2e-6)
  expect_equal(simple$p_primary, unname(stats::p.adjust(p, "holm")), tolerance = 2e-6)
})

mmrm_incomplete_prepared <- function() {
  y <- mmrm_oracle_matrix()
  y <- rbind(y, y + cbind(c(1, -2, 3, 0, 2, -1, 4, -3),
    c(3, 1, -2, 4, -1, 2, 0, -3), c(-1, 4, 2, -3, 1, 0, -2, 3)))
  raw <- data.frame(patient = rep(sprintf("S%02d", 1:16), each = 3L),
    time_label = rep(c("Zbase", "Amid", "Mfinal"), 16L),
    time_index = rep(1:3, 16L), value = as.vector(t(y)),
    arm = rep(c(rep("Zcontrol", 7L), rep("Aactive", 9L)), each = 3L),
    time_arm_nuisance = rep(c(-2, 1, 3, -1, 0, 2, -3, 4, 1, -2, 5, 0, 3, -1, 2, 4), each = 3L))
  raw$time_arm_nuisance[raw$patient == "S03"] <- NA_real_
  raw$arm[raw$patient == "S04"] <- NA_character_
  raw$value[paste(raw$patient, raw$time_label, sep = "::") %in%
    c("S01::Zbase", "S02::Amid", "S09::Mfinal", "S13::Amid")] <- NA_real_
  .mira_prepare_advanced_data_long(raw, c("Zbase", "Amid", "Mfinal"), TRUE,
    c("Zcontrol", "Aactive"), "time_arm_nuisance", character(0))
}

test_that("incomplete MMRM samples and term-aware joint tests match a direct engine fit", {
  data <- mmrm_incomplete_prepared()
  formula <- value ~ time_factor * arm_factor + time_arm_nuisance
  result <- .mira_run_mmrm_long(data, formula, TRUE, 0.9, "BH")
  expect_true(result$performed)
  visits <- c("Zbase", "Amid", "Mfinal")
  expected <- unlist(lapply(c("S01", "S02", sprintf("S%02d", 5:16)), function(id) {
    selected <- if (id == "S01") visits[2:3] else {
      if (id %in% c("S02", "S13")) visits[c(1, 3)] else {
        if (id == "S09") visits[1:2] else visits
      }
    }
    paste(id, selected, sep = "::")
  }), use.names = FALSE)
  expect_identical(length(expected), 38L)
  expect_identical(result$n, 38L)
  expect_identical(result$n_subjects, 14L)
  expect_identical(mmrm_oracle_keys(result$data), expected)
  reference_formula <- value ~ time_factor * arm_factor + time_arm_nuisance +
    us(time_factor | patient_factor)
  environment(reference_formula) <- asNamespace("mmrm")
  reference <- mmrm::mmrm(reference_formula, data = data, reml = TRUE,
    control = mmrm::mmrm_control(method = "Satterthwaite", start = mmrm::emp_start,
      accept_singular = FALSE))
  expect_identical(mmrm_oracle_keys(mmrm::component(reference, "full_frame")), expected)
  expect_equal(result$coefficients$estimate, unname(stats::coef(reference)), tolerance = 1e-8)
  expect_equal(result$covariance, mmrm::component(reference, "varcor"), tolerance = 2e-6)
  expected_terms <- list(TIME = c("time_factorAmid", "time_factorMfinal",
      "time_factorAmid:arm_factorAactive", "time_factorMfinal:arm_factorAactive"),
    ARM = c("arm_factorAactive", "time_factorAmid:arm_factorAactive",
      "time_factorMfinal:arm_factorAactive"),
    TIME_X_ARM = c("time_factorAmid:arm_factorAactive", "time_factorMfinal:arm_factorAactive"))
  expect_identical(result$effect_tests$question, names(expected_terms))
  expect_identical(anyDuplicated(result$effect_tests$question), 0L)
  beta <- stats::coef(reference)
  for (question in names(expected_terms)) {
    diagnostic <- result$effect_diagnostics[[question]]
    expect_identical(diagnostic$requested_coefficients, expected_terms[[question]])
    expect_identical(diagnostic$tested_coefficients, expected_terms[[question]])
    expect_false("time_arm_nuisance" %in% diagnostic$tested_coefficients)
    L <- diag(length(beta))[match(expected_terms[[question]], names(beta)), , drop = FALSE]
    engine <- mmrm::df_md(reference, L)
    row <- result$effect_tests[result$effect_tests$question == question, , drop = FALSE]
    expect_identical(nrow(row), 1L)
    expect_equal(c(row$statistic, row$df_num, row$df_den, row$p_raw),
      c(engine$f_stat, engine$num_df, engine$denom_df, engine$p_val), tolerance = 2e-6)
  }
})

test_that("MMRM skips structurally unsupported samples without fabricating inference", {
  data <- mmrm_oracle_long()
  duplicate <- rbind(data, data[1L, , drop = FALSE])
  result <- .mira_run_mmrm_long(duplicate, value ~ time_factor, FALSE, 0.9, "holm")
  expect_false(result$performed)
  expect_match(result$reason_skipped, "duplicate subject-visit records")
  expect_true("error" %in% names(result))
  expect_null(result$error)
  expect_null(result$model)
  expect_null(result$coefficients)
  expect_null(result$effect_tests)
  empty <- data
  empty$time_factor <- factor(empty$time_factor,
    levels = c(levels(empty$time_factor), "Unused"))
  result <- .mira_run_mmrm_long(empty, value ~ time_factor, FALSE, 0.9, "holm")
  expect_false(result$performed)
  expect_match(result$reason_skipped, "empty modeled visit")
  disjoint <- data[c(1L, 4L, 8L, 11L), , drop = FALSE]
  disjoint$time_factor <- droplevels(disjoint$time_factor)
  result <- .mira_run_mmrm_long(disjoint, value ~ time_factor, FALSE, 0.9, "holm")
  expect_false(result$performed)
  expect_match(result$reason_skipped, "no subject observed at both")
  expect_identical(result$pair_coverage[1L, 2L], 0L)
  deficient <- data
  deficient$constant_covariate <- 1
  result <- .mira_run_mmrm_long(deficient, value ~ time_factor + constant_covariate, FALSE, 0.9, "holm")
  expect_false(result$performed)
  expect_match(result$reason_skipped, "full-rank fixed-effect design")
  expect_null(result$coefficients)
  unavailable <- data
  unavailable$value[[1L]] <- NaN
  result <- .mira_run_mmrm_long(unavailable, value ~ time_factor, FALSE, 0.9, "holm")
  expect_false(result$performed)
  expect_match(result$reason_skipped, "unavailable subject, visit or response")
})

test_that("MMRM engine failures remain explicit and package absence is isolated", {
  data <- mmrm_oracle_long()
  failure <- testthat::with_mocked_bindings(
    .mira_run_mmrm_long(data, value ~ time_factor, FALSE, 0.9, "holm"),
    mmrm = function(...) stop("intentional MMRM engine failure"), .package = "mmrm")
  expect_false(failure$performed)
  expect_identical(failure$error, "intentional MMRM engine failure")
  expect_match(failure$reason_skipped, "could not be estimated")
  expect_null(failure$model)
  expect_null(failure$coefficients)
  missing <- testthat::with_mocked_bindings(
    .mira_run_mmrm_long(data, value ~ time_factor, FALSE, 0.9, "holm"),
    requireNamespace = function(...) FALSE, .package = "base")
  expect_false(missing$performed)
  expect_identical(missing$package, "mmrm")
  expect_identical(missing$reason_skipped, "Package 'mmrm' is not installed.")
  expect_null(missing$model)
})

test_that("MMRM preserves RNG and row permutations preserve the scientific model", {
  withr::local_preserve_seed()
  set.seed(20261026)
  before <- .Random.seed
  data <- mmrm_oracle_long()
  original <- data
  fit <- .mira_run_mmrm_long(data, value ~ time_factor, FALSE, 0.9, "holm")
  expect_identical(.Random.seed, before)
  expect_identical(data, original)
  permuted <- .mira_run_mmrm_long(data[rev(seq_len(nrow(data))), ],
    value ~ time_factor, FALSE, 0.9, "holm")
  expect_true(permuted$performed)
  expect_identical(.Random.seed, before)
  expect_identical(permuted$coefficients$term, fit$coefficients$term)
  expect_equal(permuted$coefficients, fit$coefficients, tolerance = 2e-6)
  expect_equal(permuted$covariance, fit$covariance, tolerance = 2e-6)
})

test_that("MMRM initialization fallback is explicit without changing the model target", {
  data <- mmrm_oracle_long()
  reference <- .mira_run_mmrm_long(data, value ~ time_factor, FALSE, 0.9, "holm")
  result <- testthat::with_mocked_bindings(
    .mira_run_mmrm_long(data, value ~ time_factor, FALSE, 0.9, "holm"),
    emp_start = function(...) stop("intentional empirical initialization failure"),
    .package = "mmrm")
  expect_true(result$performed)
  expect_true(result$converged)
  expect_identical(result$start_method, "standard fallback")
  expect_identical(result$start_fallback_reason, "intentional empirical initialization failure")
  expect_identical(result$method, "REML")
  expect_identical(result$df_method, "Satterthwaite")
  expect_equal(result$emmeans$baseline_final$tidy$estimate, 6.125, tolerance = 1e-10)
  expect_identical(result$n, 24L)
  expect_identical(result$data, reference$data)
  expect_identical(deparse(result$formula), deparse(reference$formula))
  expect_identical(mmrm::component(result$model, "cov_type"), "us")
  expect_true(mmrm::component(result$model, "reml"))
  # Different starting values can terminate at slightly different points
  # within the native optimizer criterion. Compare the fallback against the
  # actual standard-start engine, while checking its target/frame above.
  standard <- mmrm::mmrm(result$formula, data = data, reml = TRUE,
    control = mmrm::mmrm_control(method = "Satterthwaite", accept_singular = FALSE,
      start = mmrm::std_start))
  expect_identical(result$data, mmrm::component(standard, "full_frame"))
  expect_equal(result$covariance, mmrm::component(standard, "varcor"), tolerance = 1e-12)
  expect_equal(result$logLik, as.numeric(stats::logLik(standard)), tolerance = 1e-12)
  for (i in seq_len(nrow(result$coefficients))) {
    L <- numeric(nrow(result$coefficients))
    L[[i]] <- 1
    native <- mmrm::df_1d(standard, L)
    expect_equal(unname(unlist(result$coefficients[i, c("estimate", "SE", "df", "t.ratio", "p_raw")])),
      unname(unlist(native[c("est", "se", "df", "t_stat", "p_val")])), tolerance = 1e-12)
  }
})

test_that("MMRM rejects nonconverged and singular covariance inference while retaining evidence", {
  data <- mmrm_oracle_long()
  valid <- .mira_run_mmrm_long(data, value ~ time_factor, FALSE, 0.9, "holm")
  original_component <- mmrm::component
  for (corruption in c("convergence", "singular", "nonfinite")) {
    changed <- function(object, name, ...) {
      if (corruption == "convergence" && identical(name, "convergence")) return(1L)
      value <- original_component(object, name, ...)
      if (identical(name, "varcor") && corruption != "convergence") {
        if (corruption == "singular") value[] <- 1 else value[1L, 1L] <- NaN
      }
      value
    }
    result <- testthat::with_mocked_bindings(
      .mira_run_mmrm_long(data, value ~ time_factor, FALSE, 0.9, "holm"),
      mmrm = function(...) valid$model, component = changed, .package = "mmrm")
    expect_false(result$performed, info = corruption)
    expect_identical(result$model, valid$model, info = corruption)
    expect_null(result$coefficients, info = corruption)
    expect_null(result$effect_tests, info = corruption)
    expect_match(result$reason_skipped, "fit did not converge or its coefficient/residual covariance", fixed = TRUE)
  }
  scaled <- data
  scaled$value <- scaled$value * 1e-6
  scaled_result <- .mira_run_mmrm_long(scaled, value ~ time_factor, FALSE, 0.9, "holm")
  expect_true(scaled_result$performed)
  expect_identical(scaled_result$covariance_rank, valid$covariance_rank)
  expect_equal(scaled_result$coefficients$estimate / 1e-6,
    valid$coefficients$estimate, tolerance = 1e-8)
  expect_equal(scaled_result$coefficients$p_raw, valid$coefficients$p_raw, tolerance = 2e-6)
})

test_that("MMRM diagnoses actual degenerate repeated profiles without valid-looking tests", {
  data <- mmrm_oracle_long()
  data$value <- rep(c(3, 5, 4, 7, 6, 8, 5, 9), each = 3L)
  result <- .mira_run_mmrm_long(data, value ~ time_factor, FALSE, 0.9, "holm")
  expect_false(result$performed)
  expect_null(result$coefficients)
  expect_null(result$effect_tests)
  expect_true(is.character(result$reason_skipped) && length(result$reason_skipped) == 1L &&
    nzchar(result$reason_skipped))
  if (!is.null(result$model)) {
    expect_s3_class(result$model, "mmrm")
    expect_true(!result$converged || result$covariance_rank < 3L)
  } else {
    expect_true(is.character(result$error) && nzchar(result$error))
  }
})

# Boundary faults below alter only MIRA's local evaluated-result binding.
# Every fit and every unaltered df/emmeans calculation still uses real mmrm.
test_that("MMRM rejects impossible finite coefficient and joint-test probabilities locally", {
  data <- mmrm_oracle_long()
  reference <- .mira_run_mmrm_long(data, value ~ time_factor, FALSE, 0.9, "holm")
  original_eval <- MIRA:::.mira_eval_long
  for (fault in c("p_low", "p_high", "se_zero", "df_zero")) {
    changed <- function(expr) {
      call <- paste(deparse(substitute(expr)), collapse = " ")
      evaluated <- original_eval(expr)
      if (startsWith(call, "mmrm::df_1d(") &&
          isTRUE(all.equal(evaluated$value$est, reference$coefficients$estimate[[2L]]))) {
        field <- switch(fault, p_low = "p_val", p_high = "p_val", se_zero = "se", df_zero = "df")
        evaluated$value[[field]] <- switch(fault, p_low = -0.01, p_high = 1.01, se_zero = 0, df_zero = 0)
      }
      evaluated
    }
    actual <- testthat::with_mocked_bindings(
      .mira_run_mmrm_long(data, value ~ time_factor, FALSE, 0.9, "holm"),
      .mira_eval_long = changed, .package = "MIRA")
    expect_true(actual$performed, info = fault)
    expect_false(actual$coefficients$estimable[[2L]], info = fault)
    expect_true(all(is.na(actual$coefficients[2L, c("t.ratio", "p_raw", "lower.CL", "upper.CL")])), info = fault)
    expect_true(nzchar(actual$errors$coefficient_time_factorAmid), info = fault)
    expect_equal(actual$coefficients[c(1L, 3L), ], reference$coefficients[c(1L, 3L), ], tolerance = 0, info = fault)
    expect_equal(actual$effect_tests, reference$effect_tests, tolerance = 0, info = fault)
    expect_true(actual$emmeans$performed, info = fault)
  }
  for (fault in c("negative_F", "p_low", "p_high", "df_zero")) {
    changed <- function(expr) {
      call <- paste(deparse(substitute(expr)), collapse = " ")
      evaluated <- original_eval(expr)
      if (startsWith(call, "mmrm::df_md(")) {
        field <- switch(fault, negative_F = "f_stat", p_low = "p_val", p_high = "p_val", df_zero = "denom_df")
        evaluated$value[[field]] <- switch(fault, negative_F = -0.1, p_low = -0.01, p_high = 1.01, df_zero = 0)
      }
      evaluated
    }
    actual <- testthat::with_mocked_bindings(
      .mira_run_mmrm_long(data, value ~ time_factor, FALSE, 0.9, "holm"),
      .mira_eval_long = changed, .package = "MIRA")
    expect_true(actual$performed, info = fault)
    expect_identical(nrow(actual$effect_tests), 0L, info = fault)
    expect_true(is.character(actual$effect_diagnostics$TIME$reason_skipped) && grepl("inference", actual$effect_diagnostics$TIME$reason_skipped), info = fault)
    expect_equal(actual$coefficients, reference$coefficients, tolerance = 0, info = fault)
    expect_true(actual$emmeans$performed, info = fault)
  }
})

test_that("MMRM refuses covariance axes or coefficient dimensions that contradict the fitted design", {
  data <- mmrm_oracle_long()
  original_eval <- MIRA:::.mira_eval_long
  for (fault in c("visit_names", "coefficient_dimension", "coefficient_names")) {
    changed <- function(expr) {
      evaluated <- original_eval(expr)
      value <- evaluated$value
      if (is.list(value) && all(c("beta", "covariance", "beta_vcov", "data", "convergence") %in% names(value))) {
        if (fault == "visit_names") dimnames(value$covariance) <- list(rev(levels(data$time_factor)), rev(levels(data$time_factor)))
        if (fault == "coefficient_dimension") value$beta_vcov <- value$beta_vcov[1:2, 1:2, drop = FALSE]
        if (fault == "coefficient_names") dimnames(value$beta_vcov) <- list(rev(names(value$beta)), rev(names(value$beta)))
        evaluated$value <- value
      }
      evaluated
    }
    actual <- testthat::with_mocked_bindings(
      .mira_run_mmrm_long(data, value ~ time_factor, FALSE, 0.9, "holm"),
      .mira_eval_long = changed, .package = "MIRA")
    expect_false(actual$performed, info = fault)
    expect_s3_class(actual$model, "mmrm")
    expect_null(actual$coefficients)
    expect_null(actual$effect_tests)
    expect_true(is.character(actual$reason_skipped) && grepl("covariance", actual$reason_skipped), info = fault)
    expect_identical(mmrm_oracle_keys(actual$data), mmrm_oracle_keys(data), info = fault)
  }
})



test_that("MMRM keeps independent inference when a summary coefficient joint test or emmeans step fails", {
  data <- mmrm_oracle_long(arms = c(rep("Zcontrol", 3L), rep("Aactive", 5L)))
  reference <- .mira_run_mmrm_long(data, value ~ time_factor * arm_factor, TRUE, 0.9, "holm")
  original_eval <- MIRA:::.mira_eval_long
  for (fault in c("summary", "coefficient", "joint", "emmeans")) {
    coefficient_calls <- 0L
    changed <- function(expr) {
      call <- paste(deparse(substitute(expr)), collapse = " ")
      evaluated <- original_eval(expr)
      selected <- identical(fault, "summary") && identical(call, "summary(model)")
      if (startsWith(call, "mmrm::df_1d(")) {
        coefficient_calls <<- coefficient_calls + 1L
        selected <- selected || (fault == "coefficient" && coefficient_calls == 2L)
      }
      selected <- selected || (fault == "joint" && startsWith(call, "mmrm::df_md(") &&
        isTRUE(evaluated$value$num_df == 3)) ||
        (fault == "emmeans" && inherits(evaluated$value, "emmGrid") &&
          identical(get0("name", envir = parent.frame(), inherits = FALSE), "time"))
      if (selected) return(list(value = NULL, warnings = paste("controlled", fault, "warning"),
        error = paste("controlled", fault, "error")))
      evaluated
    }
    actual <- testthat::with_mocked_bindings(
      .mira_run_mmrm_long(data, value ~ time_factor * arm_factor, TRUE, 0.9, "holm"),
      .mira_eval_long = changed, .package = "MIRA")
    expect_true(actual$performed, info = fault)
    expect_identical(actual$data, reference$data)
    expect_equal(actual$covariance, reference$covariance, tolerance = 0, info = fault)
    if (fault == "coefficient") {
      expect_identical(which(!actual$coefficients$estimable), 2L)
      expect_identical(actual$errors$coefficient_time_factorAmid, "controlled coefficient error")
      expect_identical(actual$coefficients[-2L, ], reference$coefficients[-2L, ])
    } else expect_identical(actual$coefficients, reference$coefficients)
    if (fault == "joint") {
      expect_identical(actual$effect_tests$question, c("TIME", "TIME_X_ARM"))
      expect_identical(actual$errors$ARM, "controlled joint error")
      expected <- reference$effect_tests[c(1L, 3L), ]
      rownames(expected) <- NULL
      expect_identical(actual$effect_tests, expected)
    } else expect_identical(actual$effect_tests, reference$effect_tests)
    if (fault == "summary") {
      expect_null(actual$summary)
      expect_identical(actual$errors$summary, "controlled summary error")
    } else expect_s3_class(actual$summary, "summary.mmrm")
    if (fault == "emmeans") {
      expect_false(actual$emmeans$performed)
      expect_identical(actual$emmeans$error, "controlled emmeans error")
      expect_identical(actual$emmeans$warnings, "controlled emmeans warning")
    } else expect_true(actual$emmeans$performed, info = fault)
  }
})

test_that("MMRM temporal contrast families follow independent paired arithmetic at two confidence levels", {
  y <- mmrm_oracle_matrix()
  data <- mmrm_oracle_long(y)
  means <- colMeans(y)
  covariance <- crossprod(sweep(y, 2L, means)) / 7 / 8
  families <- list(
    pairwise = rbind(c(1,-1,0), c(1,0,-1), c(0,1,-1)),
    baseline_followup = rbind(c(-1,1,0), c(-1,0,1)),
    baseline_final = rbind(c(-1,0,1)),
    consecutive = rbind(c(-1,1,0), c(0,-1,1)),
    trends = rbind(c(-1,0,1), c(1,-2,1)))
  labels <- list(pairwise = c("Zbase - Amid", "Zbase - Mfinal", "Amid - Mfinal"),
    baseline_followup = c("Amid - Zbase", "Mfinal - Zbase"),
    baseline_final = "final - baseline", consecutive = c("Amid - Zbase", "Mfinal - Amid"),
    trends = c("linear", "quadratic"))
  previous <- NULL
  for (level in c(0.9, 0.95)) {
    actual <- .mira_run_mmrm_long(data, value ~ time_factor, FALSE, level, "holm")
    tables <- list(pairwise = actual$emmeans$time$pairwise_tidy,
      baseline_followup = actual$emmeans$baseline_followup$tidy,
      baseline_final = actual$emmeans$baseline_final$tidy,
      consecutive = actual$emmeans$consecutive$tidy, trends = actual$emmeans$trends$tidy)
    for (family in names(families)) {
      L <- families[[family]]
      expected <- as.vector(L %*% means)
      se <- sqrt(diag(L %*% covariance %*% t(L)))
      p <- 2 * stats::pt(-abs(expected / se), 7)
      observed <- tables[[family]]
      expect_identical(as.character(observed$contrast), labels[[family]])
      expect_equal(observed$estimate, expected, tolerance = 1e-10, info = family)
      expect_equal(observed$SE, se, tolerance = 2e-6, info = family)
      expect_equal(observed$df, rep(7, length(expected)), tolerance = 2e-5, info = family)
      expect_equal(observed$t.ratio, expected / se, tolerance = 2e-6, info = family)
      expect_equal(observed$p_raw, p, tolerance = 2e-6, info = family)
      expect_equal(observed$p_primary, stats::p.adjust(p, "holm"), tolerance = 2e-6, info = family)
      expect_equal(observed$lower.CL, expected - stats::qt((1+level)/2,7) * se, tolerance = 2e-6, info = family)
      expect_equal(observed$upper.CL, expected + stats::qt((1+level)/2,7) * se, tolerance = 2e-6, info = family)
      if (!is.null(previous)) expect_identical(observed[c("contrast","estimate","SE","df","t.ratio","p_raw","p_primary")],
        previous[[family]][c("contrast","estimate","SE","df","t.ratio","p_raw","p_primary")])
    }
    previous <- tables
  }
})

test_that("KR and MMRM obey positive affine response and bijective subject-ID transformations", {
  data <- mmrm_oracle_long()
  kr_model <- lme4::lmer(value ~ time_factor + (1 | patient_factor), data = data, REML = TRUE)
  kr <- MIRA:::.mira_run_kenward_roger_long(kr_model, FALSE)
  mmrm <- .mira_run_mmrm_long(data, value ~ time_factor, FALSE, 0.9, "holm")
  transformed <- data
  transformed$value <- 13 + 2.5 * data$value
  id_map <- stats::setNames(paste0("id ", rev(seq_len(8L))), levels(data$patient_factor))
  transformed$patient_factor <- factor(unname(id_map[as.character(data$patient_factor)]))
  transformed <- transformed[rev(seq_len(nrow(transformed))), ]
  kr_new_model <- lme4::lmer(value ~ time_factor + (1 | patient_factor), data = transformed, REML = TRUE)
  kr_new <- MIRA:::.mira_run_kenward_roger_long(kr_new_model, FALSE)
  mmrm_new <- .mira_run_mmrm_long(transformed, value ~ time_factor, FALSE, 0.9, "holm")
  expect_true(kr$performed && kr_new$performed && mmrm$performed && mmrm_new$performed)
  expect_identical(kr_new$tidy[c("question", "n", "n_subjects")], kr$tidy[c("question", "n", "n_subjects")])
  # Relative optimizer tolerance: moderate affine changes retain the same
  # Gaussian model and hypothesis; likelihood criteria themselves need not match.
  expect_equal(kr_new$tidy, kr$tidy, tolerance = 2e-6)
  expect_identical(mmrm_new$coefficients$term, mmrm$coefficients$term)
  offset <- c(13,0,0)
  expect_equal(mmrm_new$coefficients$estimate, offset + 2.5 * mmrm$coefficients$estimate, tolerance = 1e-9)
  expect_equal(mmrm_new$coefficients$SE, 2.5 * mmrm$coefficients$SE, tolerance = 2e-6)
  expect_equal(mmrm_new$coefficients$lower.CL, offset + 2.5 * mmrm$coefficients$lower.CL, tolerance = 2e-6)
  expect_equal(mmrm_new$coefficients$upper.CL, offset + 2.5 * mmrm$coefficients$upper.CL, tolerance = 2e-6)
  expect_equal(mmrm_new$coefficients$df, mmrm$coefficients$df, tolerance = 2e-5)
  # The intercept tests a changed numerical origin and is not invariant.
  expect_equal(mmrm_new$coefficients$p_raw[-1L], mmrm$coefficients$p_raw[-1L], tolerance = 2e-6)
  expect_equal(mmrm_new$effect_tests, mmrm$effect_tests, tolerance = 2e-6)
  expect_identical(dimnames(mmrm_new$covariance), list(levels(data$time_factor), levels(data$time_factor)))
  expect_equal(mmrm_new$covariance, 2.5^2 * mmrm$covariance, tolerance = 2e-6)
  expected_keys <- paste(transformed$patient_factor, transformed$time_factor, sep = "::")
  actual_keys <- mmrm_oracle_keys(mmrm_new$data)
  expect_identical(anyDuplicated(actual_keys), 0L)
  expect_identical(sort(actual_keys), sort(expected_keys))
  expect_identical(mmrm_new$data$value, transformed$value[match(actual_keys, expected_keys)])
})

test_that("three-arm MMRM joint tests retain nuisance separation and explicit restrictions", {
  y <- mmrm_oracle_matrix()
  data <- mmrm_oracle_long(rbind(y, sweep(y, 2L, c(4,-2,1), "+"),
    sweep(y, 2L, c(-3,2,5), "+")))
  data$arm_factor <- factor(rep(rep(c("Zcontrol", "Aactive", "Bthird"), each = 8L), each = 3L),
    levels = c("Zcontrol", "Aactive", "Bthird"))
  data$time_arm_nuisance <- rep(rep(c(-2,1,3,-1,0,2,-3,4), 3L), each = 3L)
  original_eval <- MIRA:::.mira_eval_long
  captured <- list()
  observed <- function(expr) {
    call <- paste(deparse(substitute(expr)), collapse = " ")
    if (startsWith(call, "mmrm::df_md(")) {
      question <- get("question", envir = parent.frame())
      captured[[question]] <<- get("contrasts", envir = parent.frame())
    }
    original_eval(expr)
  }
  actual <- testthat::with_mocked_bindings(
    .mira_run_mmrm_long(data, value ~ time_factor * arm_factor + time_arm_nuisance,
      TRUE, 0.9, "holm"), .mira_eval_long = observed, .package = "MIRA")
  expect_true(actual$performed)
  expect_identical(actual$n, 72L)
  expect_identical(actual$n_subjects, 24L)
  expect_identical(mmrm_oracle_keys(actual$data), mmrm_oracle_keys(data))
  interactions <- c("time_factorAmid:arm_factorAactive", "time_factorMfinal:arm_factorAactive",
    "time_factorAmid:arm_factorBthird", "time_factorMfinal:arm_factorBthird")
  targets <- list(TIME = c("time_factorAmid", "time_factorMfinal", interactions),
    ARM = c("arm_factorAactive", "arm_factorBthird", interactions), TIME_X_ARM = interactions)
  beta <- stats::coef(actual$model)
  expect_equal(unname(beta[c("arm_factorAactive", "arm_factorBthird")]), c(4,-3), tolerance = 1e-9)
  expect_identical(actual$effect_tests$question, names(targets))
  for (question in names(targets)) {
    expect_identical(actual$effect_diagnostics[[question]]$requested_coefficients, targets[[question]])
    expect_identical(actual$effect_diagnostics[[question]]$tested_coefficients, targets[[question]])
    L <- diag(length(beta))[match(targets[[question]], names(beta)), , drop = FALSE]
    expect_identical(captured[[question]], L)
    native <- mmrm::df_md(actual$model, L)
    row <- actual$effect_tests[actual$effect_tests$question == question, ]
    expect_equal(unname(unlist(row[c("statistic", "df_num", "df_den", "p_raw")])),
      unname(unlist(native[c("f_stat", "num_df", "denom_df", "p_val")])), tolerance = 1e-12)
  }
})
