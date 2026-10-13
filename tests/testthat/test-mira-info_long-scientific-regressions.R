scientific_regression_long <- function() {
  y <- rbind(c(3, 5, 8), c(7, 6, 10), c(4, 8, 7), c(9, 12, 11),
             c(6, 9, 12), c(12, 11, 15), c(8, 13, 16), c(15, 14, 18))
  data.frame(
    patient_factor = factor(rep(seq_len(nrow(y)), each = ncol(y))),
    time_factor = factor(rep(c("Base", "Mid", "Final"), nrow(y)),
                         levels = c("Base", "Mid", "Final")),
    time_index = rep(1:3, nrow(y)), value = as.vector(t(y)),
    z1 = rep(c(4, 6, 2, 5, 3, 9, 8, 7), each = 3L),
    z2 = rep(c(5, 3, 8, 6, 9, 7, 2, 4), each = 3L))
}

test_that("response-derived and ordinal model columns cannot be adjustment covariates", {
  data <- data.frame(subject_id = 1:30, score_t0 = 1:30,
    score_t1 = (1:30) + rep(c(1, 3, -1), 10L),
    patient_mean = 51:80, time_ordinal = 31:60, age = 21:50)
  # Both omitted names are overwritten downstream; accepting either silently
  # changes the adjustment variable. Every internal column must be reserved.
  reserved <- c("patient", "outcome", "time", "time_index", "time_label", "arm",
    "value", "patient_factor", "time_factor", "arm_factor", "patient_mean", "time_ordinal")
  for (name in reserved) {
    input <- data
    input[[name]] <- seq_len(nrow(input)) + 100
    expect_error(mira_info_long(input, id = "subject_id", outcomes = "score",
      covariates = name, inspect_only = TRUE, analyses = "none", verbose = FALSE),
      "reserved for internal", info = name)
  }
  automatic <- mira_info_long(data, id = "subject_id", outcomes = "score",
    inspect_only = TRUE, analyses = "none", verbose = FALSE)
  expect_identical(automatic$config$covariates, "age")
  expect_false(any(c("patient_mean", "time_ordinal") %in%
    automatic$detected_variables$covariates_numeric))
})

test_that("Wald rank and inference are invariant to the units of measurement", {
  # Q=beta'V^-1 beta and rank(V) are invariant under beta -> c beta,
  # V -> c^2 V for nonzero c. The rank threshold must be relative.
  for (scale in c(1, 1e-6, 1e6)) {
    full <- MIRA:::.mira_wald_chisq_long(c(1, 2) * scale,
      diag(c(1, 4)) * scale^2, 1:2)
    expect_type(full, "list")
    expect_identical(full$df, 2L)
    expect_equal(full$statistic, 2, tolerance = 1e-12)
    expect_equal(full$p_value, exp(-1), tolerance = 1e-12)
    singular <- MIRA:::.mira_wald_chisq_long(c(1, 1) * scale,
      matrix(1, 2L, 2L) * scale^2, 1:2)
    expect_identical(singular$df, 1L)
    expect_equal(singular$statistic, 1, tolerance = 1e-12)
  }
  expect_null(MIRA:::.mira_wald_chisq_long(c(0, 0), matrix(0, 2L, 2L), 1:2))
})

test_that("effect selection uses model terms rather than covariate substrings", {
  data <- scientific_regression_long()
  data$age_time_factor <- data$z1
  data$pre_arm_factor <- data$z2
  fitted <- stats::lm(value ~ time_factor + age_time_factor + pre_arm_factor, data)
  terms <- names(stats::coef(fitted))
  expect_identical(MIRA:::.mira_effect_indices_long(terms, "TIME", model = fitted), 2:3)
  expect_identical(MIRA:::.mira_effect_indices_long(terms, "ARM", model = fitted), integer(0))
  expect_identical(MIRA:::.mira_effect_indices_long(terms, "TIME_X_ARM", model = fitted), integer(0))
  data$arm_factor <- factor(rep(c("A", "B"), each = 12L))
  fitted_arm <- stats::lm(value ~ time_factor * arm_factor + age_time_factor + pre_arm_factor, data)
  arm_terms <- names(stats::coef(fitted_arm))
  expect_identical(MIRA:::.mira_effect_indices_long(arm_terms, "TIME", model = fitted_arm),
    c(2L, 3L, 7L, 8L))
  expect_identical(MIRA:::.mira_effect_indices_long(arm_terms, "ARM", model = fitted_arm),
    c(4L, 7L, 8L))
  expect_identical(MIRA:::.mira_effect_indices_long(arm_terms, "TIME_X_ARM", model = fitted_arm),
    c(7L, 8L))
  # A factor level can also contain a misleading variable name. Its text
  # cannot create an arm effect when the formula has no arm variable.
  levels(data$time_factor)[2L] <- "arm_factor_shadow"
  data$time_factor_shadow <- data$z1
  deceptive_levels <- stats::lm(value ~ time_factor + time_factor_shadow, data)
  deceptive_terms <- names(stats::coef(deceptive_levels))
  expect_identical(MIRA:::.mira_effect_indices_long(deceptive_terms, "TIME",
    model = deceptive_levels), 2:3)
  expect_identical(MIRA:::.mira_effect_indices_long(deceptive_terms, "ARM",
    model = deceptive_levels), integer(0))
  expect_identical(MIRA:::.mira_effect_indices_long(deceptive_terms, "TIME_X_ARM",
    model = deceptive_levels), integer(0))
})

test_that("RM ANOVA sensitivity selects only the requested canonical effects", {
  table <- data.frame(effect = c("lifetime", "farm", "time_factor", "arm_factor",
    "time_factor:arm_factor", "age_time_factor", "pre_arm_factor",
    "age_time_factor:pre_arm_factor", "Time_factor", "Arm_factor"))
  expect_identical(MIRA:::.mira_match_effect_long(table, "TIME"), 3L)
  expect_identical(MIRA:::.mira_match_effect_long(table, "ARM"), 4L)
  expect_identical(MIRA:::.mira_match_effect_long(table, "TIME_X_ARM"), 5L)
})

test_that("real RM ANOVA sensitivity preserves case-sensitive covariate identities", {
  skip_if_not_installed("afex")
  data <- scientific_regression_long()
  original <- MIRA:::.mira_run_rm_anova_long(data, FALSE, c("z1", "z2"), character(0))
  names(data)[names(data) == "z1"] <- "Time_factor"
  names(data)[names(data) == "z2"] <- "Arm_factor"
  renamed <- MIRA:::.mira_run_rm_anova_long(data, FALSE,
    c("Time_factor", "Arm_factor"), character(0))
  expect_true(original$performed)
  expect_true(renamed$performed)
  # R names are case-sensitive: Time_factor is a distinct baseline covariate,
  # and Arm_factor cannot create an arm hypothesis in a model without arms.
  # This is a real afex wrapper/selection check, not an independent afex oracle.
  expect_identical(MIRA:::.mira_match_effect_long(renamed$tidy, "TIME"), 3L)
  expect_identical(MIRA:::.mira_match_effect_long(renamed$tidy, "ARM"), integer(0))
  expect_identical(MIRA:::.mira_match_effect_long(renamed$tidy, "TIME_X_ARM"), integer(0))
  sensitivity <- function(rm) MIRA:::.mira_build_sensitivity_long(
    mixed_tests = list(), rm_anova = rm, friedman = list(performed = FALSE),
    gee = list(), random_slope = list(), robust = list(performed = FALSE),
    n_default = nrow(data))
  expected <- sensitivity(original)
  observed <- sensitivity(renamed)
  expect_identical(observed$question, rep("TIME", 3L))
  expect_equal(observed, expected, tolerance = 1e-10)
  row <- renamed$tidy[renamed$tidy$effect == "time_factor", ]
  expect_equal(observed$statistic, rep(row$statistic, 3L), tolerance = 1e-10)
  expect_equal(observed$p_raw, c(row$p_value, row$p_gg, row$p_hf), tolerance = 1e-10)
})

test_that("GEE scientific hypotheses are unchanged by consistent covariate renaming", {
  skip_if_not_installed("geepack")
  data <- scientific_regression_long()
  original <- MIRA:::.mira_run_gee_long(data, value ~ time_factor + z1 + z2, 0.95)
  names(data)[names(data) == "z1"] <- "age_time_factor"
  names(data)[names(data) == "z2"] <- "pre_arm_factor"
  renamed <- MIRA:::.mira_run_gee_long(data,
    value ~ time_factor + age_time_factor + pre_arm_factor, 0.95)
  for (structure in c("independence", "exchangeable", "ar1")) {
    expect_true(original[[structure]]$performed)
    expect_true(renamed[[structure]]$performed)
    expect_identical(renamed[[structure]]$effect_tests$question, "TIME")
    expect_identical(renamed[[structure]]$effect_tests$df, 2L)
    expect_equal(renamed[[structure]]$effect_tests,
      original[[structure]]$effect_tests, tolerance = 1e-10)
  }
  # Independently form the OLS cluster sandwich for the independence GEE.
  reference <- stats::lm(value ~ time_factor + age_time_factor + pre_arm_factor, data)
  design <- stats::model.matrix(reference)
  bread <- solve(crossprod(design))
  scores <- rowsum(design * stats::residuals(reference), data$patient_factor)
  covariance <- bread %*% crossprod(scores) %*% bread
  beta <- stats::coef(reference)[2:3]
  statistic <- as.numeric(crossprod(beta, solve(covariance[2:3, 2:3], beta)))
  observed <- renamed$independence$effect_tests
  expect_equal(observed$statistic, statistic, tolerance = 1e-10)
  expect_equal(observed$p_raw, exp(-statistic / 2), tolerance = 1e-10)
})

test_that("CR2 joint hypotheses are unchanged by consistent covariate renaming", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("clubSandwich")
  data <- scientific_regression_long()
  fitted <- lme4::lmer(value ~ time_factor + z1 + z2 + (1 | patient_factor), data)
  original <- MIRA:::.mira_run_robust_inference_long(fitted, 0.95)
  names(data)[names(data) == "z1"] <- "age_time_factor"
  names(data)[names(data) == "z2"] <- "pre_arm_factor"
  fitted_renamed <- lme4::lmer(value ~ time_factor + age_time_factor + pre_arm_factor +
    (1 | patient_factor), data)
  renamed <- MIRA:::.mira_run_robust_inference_long(fitted_renamed, 0.95)
  expect_true(original$performed)
  expect_true(renamed$performed)
  expect_identical(names(renamed$effect_tests), "TIME")
  expect_equal(as.data.frame(renamed$effect_tests$TIME),
    as.data.frame(original$effect_tests$TIME), tolerance = 1e-10)
})

test_that("rescaling an observed outcome preserves every GEE Wald hypothesis", {
  skip_if_not_installed("geepack")
  data <- scientific_regression_long()
  original <- MIRA:::.mira_run_gee_long(data, value ~ time_factor, 0.95)
  data$value <- data$value * 1e-6
  scaled <- MIRA:::.mira_run_gee_long(data, value ~ time_factor, 0.95)
  # The Gaussian identity-link fit rescales coefficients and covariance;
  # Wald hypotheses, df, and chi-squared inference must remain unchanged.
  for (structure in c("independence", "exchangeable", "ar1")) {
    expect_true(original[[structure]]$performed)
    expect_true(scaled[[structure]]$performed)
    expect_identical(scaled[[structure]]$effect_tests$question, "TIME")
    expect_identical(scaled[[structure]]$effect_tests$df, 2L)
    expect_equal(scaled[[structure]]$effect_tests,
      original[[structure]]$effect_tests, tolerance = 1e-10)
  }
})

test_that("Friedman statistic and Kendall W have a hand-ranked tie-correction oracle", {
  data <- data.frame(subject_id = 1:4, score_t0 = c(1, 2, 6, 1),
    score_t1 = c(1, 4, 5, 2), score_t2 = c(3, 3, 5, 3))
  result <- mira_info_long(data, id = "subject_id", outcomes = "score",
    covariates = character(0), analyses = "none", verbose = FALSE)
  # Ranks: (1.5,1.5,3),(1,3,2),(3,1.5,1.5),(1,2,3).
  # Rank sums=(6.5,8,9.5), Q_uncorrected=9/8, tie correction=7/8,
  # hence Q=9/7 and W=9/56. Chi-squared_2 survival is exp(-Q/2).
  global <- result$advanced_tests$friedman
  expect_true(global$performed)
  expect_identical(global$tidy$n, 4L)
  expect_equal(global$tidy$df, 2)
  expect_equal(global$tidy$statistic, 9 / 7, tolerance = 1e-12)
  expect_equal(global$tidy$kendalls_w, 9 / 56, tolerance = 1e-12)
  expect_equal(global$tidy$p_raw, exp(-9 / 14), tolerance = 1e-12)
  expect_equal(result$effect_sizes$friedman$kendalls_w, 9 / 56, tolerance = 1e-12)
})

