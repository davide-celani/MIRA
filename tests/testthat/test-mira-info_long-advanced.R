mira_advanced_asymmetric_fixture <- function() {
  subject <- sprintf("S%02d", 1:16)
  x <- c(-2, .4, 1.3, -1.1, 2, .8, -.7, 1.6, .2, -1.8, 1, -.3, 2.2, -.9, 1.4, -.5)
  intercept <- c(-1.5, .8, 2.1, -.7, 1.3, -2, .4, 1.7, -.3, -1.2, 2.4, .6, -1.8, 1, -.5, 1.5)
  slope <- c(.8, -.6, .2, -1, .5, -.3, 1.2, -.8, .4, 1, -.5, .7, -.9, .1, -.2, .6)
  residual <- matrix(c(
    .2, -.7, .5, .9, -.4, .8, -.6, .3, .6, -.2, .9, -.8, -.9, .4, .1, -.5,
    .3, .7, -.8, .2, -.6, .1, .4, -.7, .8, -.5, .3, .6, -.2, .9, -.4, .1,
    .7, -.3, -.9, .5, -.5, .6, .2, -.4, .4, .3, -.7, .8, -.8, -.1, .6, .4,
    .1, -.9, .7, -.2, -.7, .5, -.3, .9, .9, .2, -.5, -.6, -.3, -.4, .8, .7
  ), ncol = 4L, byrow = TRUE)
  active <- as.numeric(seq_len(16L) > 8L)
  values <- outer(10 + .8 * x + intercept + 1.2 * active, rep(1, 4)) +
    outer(slope, 0:3) + matrix(rep(c(0, 1.1, 2.5, 2.9), each = 16), nrow = 16) +
    outer(active, c(0, .4, 1.3, 2.1)) + residual
  time_vars <- c("score_t0", "score_t3", "score_t14", "score_t28")
  labels <- c("Start", "Visit10", "Visit2", "End")
  data <- data.frame(subject_id = subject,
    arm = rep(c("Z_control", "A_active"), each = 8), x = x, values)
  names(data)[4:7] <- time_vars
  data$score_t3[[3L]] <- NA_real_
  data[8L, c("score_t0", "score_t14", "score_t28")] <- NA_real_
  data$x[[2L]] <- Inf

  # Independent sample specification: exclude S02, preserve S03's gap and
  # S08's single visit. No MIRA preparation or fitted frame supplies this oracle.
  subject_index <- rep(seq_len(16L)[-2L], each = 4L)
  visit <- rep(1:4, 15L)
  keep <- !(subject_index == 3L & visit == 2L) &
    !(subject_index == 8L & visit != 2L)
  subject_index <- subject_index[keep]
  visit <- visit[keep]
  reference <- data.frame(
    patient_factor = factor(subject[subject_index]),
    time_factor = factor(labels[visit], levels = labels),
    time_index = visit, time_ordinal = visit - 1,
    arm_factor = factor(data$arm[subject_index], levels = c("Z_control", "A_active")),
    x = x[subject_index], value = values[cbind(subject_index, visit)]
  )
  # Neither input order, lexical visit labels, nor elapsed day numbers encode
  # the ordered-visit scale specified by the public model contract.
  data <- data[c(16, 2, 9, 1, 14, 4, 11, 6, 3, 15, 8, 12, 5, 10, 7, 13), ]
  list(data = data, reference = reference, labels = labels, time_vars = time_vars)
}

test_that("mira_info_long adjusts finite p-values within a family", {
  p <- c(0, 0.01, NA_real_, 0.04, 1)
  table <- data.frame(p_raw = p)
  finite <- is.finite(p)

  for (method in c("none", "holm", "bonferroni", "BH", "BY")) {
    adjusted <- MIRA:::.mira_adjust_family_long(table, primary = method)
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

  detected <- MIRA:::.mira_adjust_family_long(
    data.frame(effect = c("time", "arm"), p.value = c("<0.01", "0.5"))
  )
  expect_equal(detected$p_raw, c(0.01, 0.5))
  expect_true(all(is.na(MIRA:::.mira_adjust_family_long(
    data.frame(effect = c("time", "arm"))
  )$p_raw)))
  expect_null(MIRA:::.mira_adjust_family_long(NULL))
})

test_that("mira_info_long builds model formulas with arm and non-syntactic covariates", {
  with_arm <- MIRA:::.mira_make_fixed_formula_long(
    TRUE, covariates = c("age years", "baseline-score")
  )
  expect_setequal(all.vars(with_arm),
                  c("value", "time_factor", "arm_factor", "age years", "baseline-score"))
  expect_match(MIRA:::.mira_formula_text_long(with_arm),
               "time_factor \\* arm_factor")

  without_arm <- MIRA:::.mira_make_fixed_formula_long(FALSE)
  expect_setequal(all.vars(without_arm), c("value", "time_factor"))
  expect_false("arm_factor" %in% all.vars(without_arm))

  random_slope <- MIRA:::.mira_make_lmer_formula_long(
    "time_factor", "age years", "(1 + time_ordinal | patient_factor)"
  )
  expect_setequal(all.vars(random_slope),
                  c("value", "time_factor", "age years", "time_ordinal", "patient_factor"))
  expect_match(MIRA:::.mira_formula_text_long(random_slope), "time_ordinal")
  expect_true(is.na(MIRA:::.mira_formula_text_long(NULL)))
})

test_that("mira_info_long's protected evaluator records distinct warnings and errors", {
  warned <- MIRA:::.mira_eval_long({
    warning("imprecise estimate")
    warning("imprecise estimate")
    7L
  })
  expect_identical(warned$value, 7L)
  expect_identical(warned$warnings, "imprecise estimate")
  expect_null(warned$error)

  failed <- MIRA:::.mira_eval_long(stop("model did not converge"))
  expect_null(failed$value)
  expect_identical(failed$error, "model did not converge")
  expect_length(failed$warnings, 0L)
})

test_that("mira_info_long advanced data preparation filters unusable values and orders visits", {
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
  prepared <- MIRA:::.mira_prepare_advanced_data_long(
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

test_that("mira_info_long standardizes optional analysis tables across package formats", {
  afex_table <- data.frame(
    Effect = "time", F = 4.5, `num Df` = 2, `den Df` = 10,
    `Pr(>F)` = "<0.02", pes = 0.3, ges = 0.2,
    check.names = FALSE
  )
  tidy <- MIRA:::.mira_standardize_afex_long(afex_table)
  expect_identical(tidy$effect, "time")
  expect_equal(tidy$statistic, 4.5)
  expect_equal(tidy$df_num, 2)
  expect_equal(tidy$df_den, 10)
  expect_equal(tidy$p_value, 0.02)
  expect_equal(tidy$partial_eta_squared, 0.3)
  expect_equal(tidy$generalized_eta_squared, 0.2)
  expect_null(MIRA:::.mira_standardize_afex_long(NULL))

  expect_equal(MIRA:::.mira_qic_value_long(c(QIC = 12, CIC = 10), "qic"), 12)
  expect_equal(MIRA:::.mira_qic_value_long(
    matrix(c(12, 10), nrow = 1L, dimnames = list(NULL, c("QIC", "CIC"))), "cic"
  ), 10)
  expect_true(is.na(MIRA:::.mira_qic_value_long(NULL, "QIC")))
  expect_true(is.na(MIRA:::.mira_qic_value_long(c(QIC = 12), "CIC")))
})

test_that("mira_info_long Wald helper respects selected effects and covariance rank", {
  terms <- c("(Intercept)", "time_factorFollow", "arm_factorDrug",
             "time_factorFollow:arm_factorDrug")
  expect_identical(MIRA:::.mira_effect_indices_long(terms, "TIME_X_ARM"), 4L)
  expect_identical(MIRA:::.mira_effect_indices_long(terms, "unknown"), integer(0))
  effect_table <- data.frame(effect = c("time_factor", "arm_factor",
                                       "time_factor:arm_factor", "Residual"))
  expect_identical(MIRA:::.mira_match_effect_long(effect_table, "TIME"), 1L)
  expect_identical(MIRA:::.mira_match_effect_long(effect_table, "ARM"), 2L)
  expect_identical(MIRA:::.mira_match_effect_long(effect_table, "TIME_X_ARM"), 3L)

  wald <- MIRA:::.mira_wald_chisq_long(
    beta = c(1, 2, -1), covariance = diag(c(1, 4, 1)), indices = c(2L, 3L)
  )
  expect_equal(wald$statistic, 2)
  expect_equal(wald$df, 2L)
  expect_equal(wald$p_value, stats::pchisq(2, df = 2, lower.tail = FALSE))
  expect_equal(dim(wald$constraints), c(2L, 3L))
  expect_null(MIRA:::.mira_wald_chisq_long(c(1, 2), diag(2), integer(0)))
  expect_null(MIRA:::.mira_wald_chisq_long(c(1, 2), diag(c(0, 1)), 1L))
})

test_that("mira_info_long sensitivity extraction selects the tested comparison", {
  lrt <- data.frame(
    Chisq = c(NA_real_, 6.2), Df = c(NA_real_, 2),
    `Pr(>Chisq)` = c(NA_real_, 0.045), check.names = FALSE
  )
  extracted <- MIRA:::.mira_lrt_sensitivity_long(
    lrt, question = "TIME", method = "mixed-model LRT",
    object_path = "$model$global_time_test", n = 24L
  )
  expect_equal(extracted$statistic, 6.2)
  expect_equal(extracted$df, "2")
  expect_equal(extracted$p_raw, 0.045)
  expect_equal(extracted$n, 24L)
  expect_identical(extracted$question, "TIME")
  expect_null(MIRA:::.mira_lrt_sensitivity_long(
    NULL, "TIME", "LRT", "$model$global_time_test", 24L
  ))
})

test_that("mira_info_long Friedman statistics agree with independent base-R calculations", {
  analysis <- data.frame(
    score_t0 = c(1, 3, 2, 6, 5, 4),
    score_t1 = c(2, 4, 3, 5, 7, 6),
    score_t2 = c(4, 5, 6, 7, 8, 9)
  )
  variables <- names(analysis)
  labels <- stats::setNames(c("Base", "Middle", "Final"), variables)
  result <- MIRA:::.mira_run_friedman_long(analysis, variables, labels, "BH")
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

test_that("mira_info_long paired rank-biserial effects handle signs and unusable changes", {
  expect_equal(MIRA:::.mira_rank_biserial_paired_long(c(1, 2, 3)), 1)
  expect_equal(MIRA:::.mira_rank_biserial_paired_long(c(-1, -2, -3)), -1)
  expect_equal(MIRA:::.mira_rank_biserial_paired_long(c(1, 2, -3)), 0)
  expect_true(is.na(MIRA:::.mira_rank_biserial_paired_long(c(0, NA, Inf))))
})

test_that("mira_info_long skip and protected-failure objects preserve diagnostics", {
  skipped <- MIRA:::.mira_skipped_module_long("emmeans", "model unavailable")
  expect_false(skipped$performed)
  expect_identical(skipped$package, "emmeans")
  expect_identical(skipped$reason_skipped, "model unavailable")
  expect_null(skipped$object)
  expect_null(skipped$error)

  no_model <- MIRA:::.mira_run_emmeans_long(
    NULL, arm_enabled = FALSE, primary_method = "holm", confidence_level = 0.95
  )
  expect_false(no_model$performed)
  expect_match(no_model$reason_skipped, "mixed model")

  failed <- MIRA:::.mira_failed_advanced_long("upstream failure", "fit warning")
  expect_identical(failed$error, "upstream failure")
  expect_identical(failed$warnings, "fit warning")
  expect_false(failed$advanced_tests$friedman$performed)
  expect_false(failed$advanced_models$random_slope$performed)
  expect_identical(failed$advanced_tests$emmeans$error, "upstream failure")
  expect_s3_class(failed$sensitivity, "data.frame")
  expect_equal(nrow(failed$sensitivity), 0L)
})

test_that("mira_info_long keeps Friedman inference while model-based modules are disabled", {
  data <- data.frame(
    patient = 1:6,
    score_t0 = c(1, 3, 2, 6, 5, 4),
    score_t1 = c(2, 4, 3, 5, 7, 6),
    score_t2 = c(4, 5, 6, 7, 8, 9)
  )
  result <- mira_info_long(
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

test_that("mira_info_long fits the main and advanced models on sufficient repeated data", {
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
  result <- mira_info_long(
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

  fixture <- mira_advanced_asymmetric_fixture()
  reference <- fixture$reference
  asymmetric <- mira_info_long(
    fixture$data, id = "subject_id", outcomes = "score",
    time_vars = fixture$time_vars, time_labels = fixture$labels,
    arm = "arm", reference_arm = "Z_control", covariates = "x",
    alpha = .10, analyses = c("model", "arm_tests"), verbose = FALSE
  )
  fixed <- value ~ time_factor * arm_factor + x
  expected_design <- stats::model.matrix(fixed, reference)
  expected_keys <- paste(reference$patient_factor, reference$time_factor, sep = "::")
  expect_length(expected_keys, 56L)
  expect_identical(anyDuplicated(expected_keys), 0L)

  check_frame <- function(frame, id) {
    keys <- paste(id, frame$time_factor, sep = "::")
    expect_identical(keys, expected_keys)
    expect_identical(anyDuplicated(keys), 0L)
    expect_identical(levels(frame$time_factor), fixture$labels)
    expect_identical(levels(frame$arm_factor), c("Z_control", "A_active"))
    expect_equal(frame$value, reference$value, tolerance = 0)
    expect_equal(frame$x, reference$x, tolerance = 0)
    expect_identical(as.character(frame$arm_factor), as.character(reference$arm_factor))
  }
  check_formula <- function(formula) {
    expect_identical(all.vars(formula), all.vars(fixed))
    expect_identical(attr(stats::terms(formula), "term.labels"),
                     attr(stats::terms(fixed), "term.labels"))
    expect_identical(attr(stats::terms(formula), "intercept"), 1L)
  }
  # Same-engine oracles isolate wrapper specification and extraction. 1e-8 is
  # numerical agreement for the same ordered data/engine, not a statistical
  # uncertainty band or an equivalence requirement across model families.
  for (structure in c("exchangeable", "ar1")) {
    item <- asymmetric$advanced_models$gee[[structure]]
    expect_true(item$performed)
    expect_true(item$converged)
    expect_null(item$error)
    expect_length(item$warnings, 0L)
    expect_identical(item$model$corstr, structure)
    expect_identical(item$model$std.err, "san.se")
    expect_identical(item$model$family$family, "gaussian")
    expect_identical(item$model$family$link, "identity")
    frame <- stats::model.frame(item$model)
    check_frame(frame, as.character(item$model$id))
    check_formula(stats::formula(item$model))
    expect_equal(unname(stats::model.matrix(item$model)), unname(expected_design), tolerance = 0)
    oracle <- geepack::geeglm(fixed, data = reference,
      id = patient_factor, waves = time_index, family = stats::gaussian(),
      corstr = structure, std.err = "san.se")
    expect_equal(stats::coef(item$model), stats::coef(oracle), tolerance = 1e-8)
    expect_equal(stats::vcov(item$model), stats::vcov(oracle), tolerance = 1e-8)
    expect_equal(item$model$geese$alpha, oracle$geese$alpha, tolerance = 1e-8)
    robust_se <- sqrt(diag(stats::vcov(oracle)))
    expect_identical(item$coefficients$term, names(stats::coef(oracle)))
    expect_equal(item$coefficients$Estimate, unname(stats::coef(oracle)), tolerance = 1e-8)
    expect_equal(item$coefficients$robust_se, unname(robust_se), tolerance = 1e-8)
    expect_equal(item$coefficients$conf_low,
      unname(stats::coef(oracle) - stats::qnorm(.95) * robust_se), tolerance = 1e-8)
    expect_equal(item$coefficients$conf_high,
      unname(stats::coef(oracle) + stats::qnorm(.95) * robust_se), tolerance = 1e-8)
  }

  for (structure in c("compound_symmetry", "ar1")) {
    item <- asymmetric$advanced_models$nlme[[structure]]
    expect_true(item$performed)
    expect_true(item$converged)
    expect_null(item$error)
    expect_length(item$warnings, 0L)
    expect_identical(item$model$method, "REML")
    frame <- nlme::getData(item$model)
    check_frame(frame, as.character(frame$patient_factor))
    expect_equal(stats::nobs(item$model), 56L, tolerance = 0)
    expect_identical(as.character(item$model$groups$patient_factor),
      as.character(reference$patient_factor))
    check_formula(stats::formula(item$model))
    expect_identical(frame$time_index, reference$time_index)
    correlation <- if (structure == "ar1") {
      nlme::corAR1(form = ~time_index | patient_factor)
    } else nlme::corCompSymm(form = ~1 | patient_factor)
    oracle <- nlme::lme(fixed, random = ~1 | patient_factor,
      correlation = correlation, data = reference, method = "REML",
      na.action = stats::na.omit, control = nlme::lmeControl(returnObject = TRUE))
    # nlme initializes AR1 as corARMA when a subject has a visit gap.
    expect_identical(class(item$model$modelStruct$corStruct),
      class(oracle$modelStruct$corStruct))
    expect_identical(deparse(stats::formula(item$model$modelStruct$corStruct)),
      deparse(stats::formula(oracle$modelStruct$corStruct)))
    expect_equal(nlme::fixef(item$model), nlme::fixef(oracle), tolerance = 1e-8)
    expect_equal(stats::vcov(item$model), stats::vcov(oracle), tolerance = 1e-8)
    expect_equal(item$coefficients, summary(oracle)$tTable, tolerance = 1e-8)
    expect_equal(item$logLik, as.numeric(stats::logLik(oracle)), tolerance = 1e-8)
    # A random intercept plus CS residuals does not separately identify both
    # constant covariance components. Compare the estimable total covariance,
    # not a supposedly unique split into intercept variance and CS rho.
    expect_equal(nlme::getVarCov(item$model, individual = "S01", type = "marginal"),
      nlme::getVarCov(oracle, individual = "S01", type = "marginal"), tolerance = 1e-8)
    if (structure == "ar1") {
      expect_equal(nlme::getCovariate(item$model$modelStruct$corStruct),
        split(reference$time_index, reference$patient_factor), tolerance = 0)
      expect_equal(item$correlation_parameter,
        unname(stats::coef(oracle$modelStruct$corStruct, unconstrained = FALSE)),
        tolerance = 1e-8)
    }
  }

  item <- asymmetric$advanced_models$random_slope
  expect_true(item$performed)
  expect_true(item$converged)
  expect_false(item$singular)
  expect_null(item$error)
  expect_length(item$warnings, 0L)
  expect_true(lme4::isREML(item$model))
  frame <- stats::model.frame(item$model)
  check_frame(frame, as.character(frame$patient_factor))
  expect_equal(frame$time_ordinal, reference$time_ordinal, tolerance = 0)
  slope_formula <- value ~ time_factor * arm_factor + x +
    (1 + time_ordinal | patient_factor)
  expect_identical(all.vars(stats::formula(item$model)), all.vars(slope_formula))
  slope_design <- lme4::getME(item$model, "X")
  expect_identical(dim(slope_design), dim(expected_design))
  expect_identical(colnames(slope_design), colnames(expected_design))
  expect_equal(as.vector(slope_design), as.vector(expected_design), tolerance = 0)
  oracle <- lme4::lmer(slope_formula, data = reference, REML = TRUE,
    na.action = stats::na.omit)
  expect_false(lme4::isSingular(oracle))
  expect_equal(lme4::fixef(item$model), lme4::fixef(oracle), tolerance = 1e-8)
  expect_equal(as.matrix(stats::vcov(item$model)), as.matrix(stats::vcov(oracle)), tolerance = 1e-8)
  expect_equal(item$variance_components, as.data.frame(lme4::VarCorr(oracle)), tolerance = 1e-8)
  expect_equal(item$logLik, as.numeric(stats::logLik(oracle)), tolerance = 1e-8)
  expect_identical(item$fixed_effects$term, names(lme4::fixef(oracle)))
  expect_equal(item$fixed_effects$Estimate, unname(lme4::fixef(oracle)), tolerance = 1e-8)

  comparison <- asymmetric$advanced_models$model_comparison
  selected <- comparison$model %in% c("gee_exchangeable", "gee_ar1",
    "nlme_compound_symmetry", "nlme_ar1", "mixed_random_slope")
  expect_equal(comparison$n[selected], rep(56L, 5L), tolerance = 0)
  expect_equal(comparison$n_subjects[selected], rep(15L, 5L), tolerance = 0)
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
    result <- mira_info_long(input, id = "patient", outcomes = "score",
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
