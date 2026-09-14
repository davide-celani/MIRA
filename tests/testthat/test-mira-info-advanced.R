test_that("mira_info adjusts finite p-values within a family", {
  p <- c(0, 0.01, NA_real_, 0.04, 1)
  table <- data.frame(p_raw = p)
  finite <- is.finite(p)

  for (method in c("none", "holm", "bonferroni", "BH", "BY")) {
    adjusted <- MIRA:::.mira_adjust_family(table, primary = method)
    expect_equal(adjusted$p_primary[finite],
                 stats::p.adjust(p[finite], method = method))
    expect_true(is.na(adjusted$p_primary[!finite]))
    expect_identical(adjusted$primary_method, rep(method, length(p)))
    expect_equal(adjusted$p_bonferroni[finite],
                 stats::p.adjust(p[finite], method = "bonferroni"))
    expect_equal(adjusted$p_holm[finite],
                 stats::p.adjust(p[finite], method = "holm"))
    expect_equal(adjusted$p_bh[finite],
                 stats::p.adjust(p[finite], method = "BH"))
    expect_equal(adjusted$p_by[finite],
                 stats::p.adjust(p[finite], method = "BY"))
  }

  detected <- MIRA:::.mira_adjust_family(
    data.frame(effect = c("time", "arm"), p.value = c("<0.01", "0.5"))
  )
  expect_equal(detected$p_raw, c(0.01, 0.5))
  expect_true(all(is.na(MIRA:::.mira_adjust_family(
    data.frame(effect = c("time", "arm"))
  )$p_raw)))
  expect_null(MIRA:::.mira_adjust_family(NULL))
})

test_that("mira_info builds model formulas with arm and non-syntactic covariates", {
  with_arm <- MIRA:::.mira_make_fixed_formula(
    TRUE, covariates = c("age years", "baseline-score")
  )
  expect_setequal(all.vars(with_arm),
                  c("value", "time_factor", "arm_factor", "age years", "baseline-score"))
  expect_match(MIRA:::.mira_formula_text(with_arm),
               "time_factor \\* arm_factor")

  without_arm <- MIRA:::.mira_make_fixed_formula(FALSE)
  expect_setequal(all.vars(without_arm), c("value", "time_factor"))
  expect_false("arm_factor" %in% all.vars(without_arm))

  random_slope <- MIRA:::.mira_make_lmer_formula(
    "time_factor", "age years", "(1 + time_ordinal | patient_factor)"
  )
  expect_setequal(all.vars(random_slope),
                  c("value", "time_factor", "age years", "time_ordinal", "patient_factor"))
  expect_match(MIRA:::.mira_formula_text(random_slope), "time_ordinal")
  expect_true(is.na(MIRA:::.mira_formula_text(NULL)))
})

test_that("mira_info's protected evaluator records distinct warnings and errors", {
  warned <- MIRA:::.mira_eval({
    warning("imprecise estimate")
    warning("imprecise estimate")
    7L
  })
  expect_identical(warned$value, 7L)
  expect_identical(warned$warnings, "imprecise estimate")
  expect_null(warned$error)

  failed <- MIRA:::.mira_eval(stop("model did not converge"))
  expect_null(failed$value)
  expect_identical(failed$error, "model did not converge")
  expect_length(failed$warnings, 0L)
})

test_that("mira_info advanced data preparation filters unusable values and orders visits", {
  long <- data.frame(
    patient = c(2, 1, 1, 2, 3),
    time_index = c(2, 1, 2, 1, 1),
    time_label = c("Follow", "Base", "Follow", "Base", "Base"),
    value = c(6, 2, Inf, 4, 5),
    arm = c("drug", "placebo", "placebo", "drug", "drug"),
    age = c(30, 20, 20, Inf, 40),
    group = c("A", "B", "B", "A", NA_character_),
    flag = c(TRUE, FALSE, FALSE, TRUE, TRUE)
  )
  prepared <- MIRA:::.mira_prepare_advanced_data(
    long, time_levels = c("Base", "Follow"), arm_enabled = TRUE,
    arm_levels = c("placebo", "drug"),
    covariates = c("age", "group", "flag"),
    categorical_covariates = c("group", "flag")
  )
  expect_equal(prepared$patient, c(1, 2))
  expect_equal(prepared$value, c(2, 6))
  expect_identical(levels(prepared$time_factor), c("Base", "Follow"))
  expect_identical(levels(prepared$arm_factor), c("placebo", "drug"))
  expect_true(is.factor(prepared$group))
  expect_true(is.factor(prepared$flag))
  expect_true(all(is.finite(prepared$value)))
  expect_false(anyNA(prepared[c("age", "group", "flag")]))

})

test_that("mira_info standardizes optional analysis tables across package formats", {
  afex_table <- data.frame(
    Effect = "time", F = 4.5, `num Df` = 2, `den Df` = 10,
    `Pr(>F)` = "<0.02", pes = 0.3, ges = 0.2,
    check.names = FALSE
  )
  tidy <- MIRA:::.mira_standardize_afex(afex_table)
  expect_identical(tidy$effect, "time")
  expect_equal(tidy$statistic, 4.5)
  expect_equal(tidy$df_num, 2)
  expect_equal(tidy$df_den, 10)
  expect_equal(tidy$p_value, 0.02)
  expect_equal(tidy$partial_eta_squared, 0.3)
  expect_equal(tidy$generalized_eta_squared, 0.2)
  expect_null(MIRA:::.mira_standardize_afex(NULL))

  expect_equal(MIRA:::.mira_qic_value(c(QIC = 12, CIC = 10), "qic"), 12)
  expect_equal(MIRA:::.mira_qic_value(
    matrix(c(12, 10), nrow = 1L, dimnames = list(NULL, c("QIC", "CIC"))), "cic"
  ), 10)
  expect_true(is.na(MIRA:::.mira_qic_value(NULL, "QIC")))
  expect_true(is.na(MIRA:::.mira_qic_value(c(QIC = 12), "CIC")))
})

test_that("mira_info Wald helper respects selected effects and covariance rank", {
  terms <- c("(Intercept)", "time_factorFollow", "arm_factorDrug",
             "time_factorFollow:arm_factorDrug")
  expect_identical(MIRA:::.mira_effect_indices(terms, "TIME_X_ARM"), 4L)
  expect_identical(MIRA:::.mira_effect_indices(terms, "unknown"), integer(0))
  effect_table <- data.frame(effect = c("time_factor", "arm_factor",
                                       "time_factor:arm_factor", "Residual"))
  expect_identical(MIRA:::.mira_match_effect(effect_table, "TIME"), 1L)
  expect_identical(MIRA:::.mira_match_effect(effect_table, "ARM"), 2L)
  expect_identical(MIRA:::.mira_match_effect(effect_table, "TIME_X_ARM"), 3L)

  wald <- MIRA:::.mira_wald_chisq(
    beta = c(1, 2, -1), covariance = diag(c(1, 4, 1)), indices = c(2L, 3L)
  )
  expect_equal(wald$statistic, 2)
  expect_equal(wald$df, 2L)
  expect_equal(wald$p_value, stats::pchisq(2, df = 2, lower.tail = FALSE))
  expect_equal(dim(wald$constraints), c(2L, 3L))
  expect_null(MIRA:::.mira_wald_chisq(c(1, 2), diag(2), integer(0)))
  expect_null(MIRA:::.mira_wald_chisq(c(1, 2), diag(c(0, 1)), 1L))
})

test_that("mira_info sensitivity extraction selects the tested comparison", {
  lrt <- data.frame(
    Chisq = c(NA_real_, 6.2), Df = c(NA_real_, 2),
    `Pr(>Chisq)` = c(NA_real_, 0.045), check.names = FALSE
  )
  extracted <- MIRA:::.mira_lrt_sensitivity(
    lrt, question = "TIME", method = "mixed-model LRT",
    object_path = "$model$global_time_test", n = 24L
  )
  expect_equal(extracted$statistic, 6.2)
  expect_equal(extracted$df, "2")
  expect_equal(extracted$p_raw, 0.045)
  expect_equal(extracted$n, 24L)
  expect_identical(extracted$question, "TIME")
  expect_null(MIRA:::.mira_lrt_sensitivity(
    NULL, "TIME", "LRT", "$model$global_time_test", 24L
  ))
})

test_that("mira_info Friedman statistics agree with independent base-R calculations", {
  analysis <- data.frame(
    score_t0 = c(1, 3, 2, 6, 5, 4),
    score_t1 = c(2, 4, 3, 5, 7, 6),
    score_t2 = c(4, 5, 6, 7, 8, 9)
  )
  variables <- names(analysis)
  labels <- stats::setNames(c("Base", "Middle", "Final"), variables)
  result <- MIRA:::.mira_run_friedman(analysis, variables, labels, "BH")
  reference <- stats::friedman.test(as.matrix(analysis))

  expect_true(result$performed)
  expect_equal(result$tidy$statistic, unname(reference$statistic))
  expect_equal(result$tidy$p_raw, reference$p.value)
  expect_equal(result$tidy$kendalls_w,
               unname(reference$statistic) / (nrow(analysis) * (ncol(analysis) - 1L)))
  expect_equal(result$tidy$n, nrow(analysis))
  expect_equal(nrow(result$posthoc$tidy), choose(ncol(analysis), 2L))

  pair <- result$posthoc$tidy
  expected <- vapply(seq_len(nrow(pair)), function(i) {
    stats::wilcox.test(analysis[[pair$to[[i]]]], analysis[[pair$from[[i]]]],
                       paired = TRUE, exact = FALSE)$p.value
  }, numeric(1L))
  expect_equal(pair$p_raw, expected)
  expect_equal(pair$p_primary, stats::p.adjust(expected, method = "BH"))
  expect_true(all(pair$rank_biserial >= -1 & pair$rank_biserial <= 1))
})

test_that("mira_info paired rank-biserial effects handle signs and unusable changes", {
  expect_equal(MIRA:::.mira_rank_biserial_paired(c(1, 2, 3)), 1)
  expect_equal(MIRA:::.mira_rank_biserial_paired(c(-1, -2, -3)), -1)
  expect_equal(MIRA:::.mira_rank_biserial_paired(c(1, 2, -3)), 0)
  expect_true(is.na(MIRA:::.mira_rank_biserial_paired(c(0, NA, Inf))))
})

test_that("mira_info skip and protected-failure objects preserve diagnostics", {
  skipped <- MIRA:::.mira_skipped_module("emmeans", "model unavailable")
  expect_false(skipped$performed)
  expect_identical(skipped$package, "emmeans")
  expect_identical(skipped$reason_skipped, "model unavailable")
  expect_null(skipped$object)
  expect_null(skipped$error)

  no_model <- MIRA:::.mira_run_emmeans(
    NULL, arm_enabled = FALSE, primary_method = "holm", confidence_level = 0.95
  )
  expect_false(no_model$performed)
  expect_match(no_model$reason_skipped, "mixed model")

  failed <- MIRA:::.mira_failed_advanced("upstream failure", "fit warning")
  expect_identical(failed$error, "upstream failure")
  expect_identical(failed$warnings, "fit warning")
  expect_false(failed$advanced_tests$friedman$performed)
  expect_false(failed$advanced_models$random_slope$performed)
  expect_identical(failed$advanced_tests$emmeans$error, "upstream failure")
  expect_s3_class(failed$sensitivity, "data.frame")
  expect_equal(nrow(failed$sensitivity), 0L)
})

test_that("mira_info keeps Friedman inference while model-based modules are disabled", {
  data <- data.frame(
    patient = 1:6,
    score_t0 = c(1, 3, 2, 6, 5, 4),
    score_t1 = c(2, 4, 3, 5, 7, 6),
    score_t2 = c(4, 5, 6, 7, 8, 9)
  )
  result <- mira_info(
    data, id = "patient", time_vars = c("score_t0", "score_t1", "score_t2"),
    outcomes = "score", analyses = "none", verbose = FALSE
  )
  expect_true(result$advanced_tests$friedman$performed)
  expect_equal(result$advanced_tests$friedman$tidy$n, nrow(data))
  for (module in c("emmeans", "rm_anova")) {
    expect_false(result$advanced_tests[[module]]$performed)
    expect_match(result$advanced_tests[[module]]$reason_skipped, "model = FALSE")
  }
  expect_false(result$advanced_models$random_slope$performed)
  expect_false(result$advanced_models$gee$performed)
  expect_false(result$advanced_models$nlme$performed)
  expect_false(result$robustness$club_sandwich$performed)
  expect_false(result$config$analyses$model)
  expect_true(is.null(result$model$fitted_model))
})

test_that("mira_info fits the main and advanced models on sufficient repeated data", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("lmerTest")
  skip_if_not_installed("emmeans")
  skip_if_not_installed("afex")
  skip_if_not_installed("geepack")

  data <- data.frame(
    patient = seq_len(12),
    score_t0 = c(3, 5, 4, 7, 6, 8, 5, 9, 7, 10, 8, 11),
    score_t1 = c(4, 6, 5, 6, 8, 9, 7, 10, 8, 11, 9, 12),
    score_t2 = c(6, 7, 7, 8, 9, 11, 8, 12, 10, 13, 11, 14)
  )
  result <- mira_info(
    data, id = "patient", time_vars = c("score_t0", "score_t1", "score_t2"),
    outcomes = "score", analyses = "model", verbose = FALSE
  )

  expect_true(inherits(result$model$fitted_model, "merMod"))
  expect_equal(stats::nobs(result$model$fitted_model), nrow(data) * 3L)
  expect_null(result$model$error)
  expect_true("time_factor" %in% all.vars(stats::formula(result$model$fitted_model)))
  expect_true(result$advanced_tests$emmeans$performed)
  expect_true(result$advanced_tests$rm_anova$performed)
  expect_true(result$advanced_tests$friedman$performed)
  expect_s3_class(result$advanced_tests$rm_anova$tidy, "data.frame")
  expect_s3_class(result$advanced_models$model_comparison, "data.frame")
  expect_setequal(result$advanced_models$model_comparison$model,
                  c("mixed_random_intercept", "mixed_random_slope",
                    "nlme_compound_symmetry", "nlme_ar1", "gee_independence",
                    "gee_exchangeable", "gee_ar1"))
  expect_true(all(is.na(result$advanced_models$model_comparison$AIC[
    grepl("^gee_", result$advanced_models$model_comparison$model)
  ])))
  expect_s3_class(result$sensitivity, "data.frame")
  expect_identical(result$multiplicity$primary_method, "holm")
  expect_equal(result$multiplicity$friedman_posthoc$table$p_primary,
               stats::p.adjust(result$advanced_tests$friedman$posthoc$tidy$p_raw,
                               method = "holm"))
  # On these complete repeated data, base geepack::geeglm() can estimate the
  # independence model; the public pipeline should expose that branch too.
  expect_true(result$advanced_models$gee$independence$performed)
})

test_that("blank categorical covariates use the same subjects across model families", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("lmerTest")
  skip_if_not_installed("geepack")

  data <- data.frame(
    patient = seq_len(12),
    score_t0 = c(3, 5, 4, 7, 6, 8, 5, 9, 7, 10, 8, 11),
    score_t1 = c(4, 6, 5, 6, 8, 9, 7, 10, 8, 11, 9, 12),
    score_t2 = c(6, 7, 7, 8, 9, 11, 8, 12, 10, 13, 11, 14),
    sex = rep(c("F", "M"), each = 6L),
    stringsAsFactors = FALSE
  )
  data$sex[[1L]] <- " "

  for (as_factor in c(FALSE, TRUE)) {
    input <- data
    if (as_factor) input$sex <- factor(input$sex)
    result <- mira_info(input, id = "patient", outcomes = "score",
                        covariates = "sex", analyses = "model", verbose = FALSE)

    expect_equal(result$overview$n_patients, 12L)
    expect_equal(nrow(result$long_data), 36L)
    expect_equal(stats::nobs(result$model$fitted_model), 33L)
    expect_true(result$advanced_models$gee$independence$performed)
    expect_equal(stats::nobs(result$advanced_models$gee$independence$model), 33L)
    comparison <- result$advanced_models$model_comparison
    relevant <- comparison$model %in% c("mixed_random_intercept", "gee_independence")
    expect_equal(comparison$n[relevant], c(33L, 33L))
    expect_equal(comparison$n_subjects[relevant], c(11L, 11L))
  }
})
