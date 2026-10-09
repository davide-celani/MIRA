model_oracle_measurements <- function() {
  rbind(c(3, 5, 8), c(7, 6, 10), c(4, 8, 7), c(9, 12, 11),
        c(6, 9, 12), c(12, 11, 15), c(8, 13, 16), c(15, 14, 18))
}

model_oracle_long <- function() {
  y <- model_oracle_measurements()
  data.frame(
    patient_factor = factor(rep(seq_len(nrow(y)), each = ncol(y))),
    time_factor = factor(rep(c("Base", "Month 4", "Month 12"), nrow(y)),
                         levels = c("Base", "Month 4", "Month 12")),
    time_index = rep(1:3, nrow(y)),
    time_ordinal = rep(0:2, nrow(y)),
    value = as.vector(t(y))
  )
}

test_that("repeated-measures ANOVA F and effect sizes match a balanced linear-model oracle", {
  skip_if_not_installed("afex")
  data <- model_oracle_long()
  result <- MIRA:::.mira_run_rm_anova_long(data, FALSE, character(0), character(0))
  # Subject and time sums of squares are orthogonal in this balanced design.
  reference <- stats::anova(stats::lm(value ~ patient_factor + time_factor, data))
  ss_time <- reference["time_factor", "Sum Sq"]
  ss_error <- reference["Residuals", "Sum Sq"]
  ss_subject <- reference["patient_factor", "Sum Sq"]
  row <- result$tidy[result$tidy$effect == "time_factor", ]
  expect_true(result$performed)
  expect_identical(nrow(row), 1L)
  expect_identical(result$n_subjects, 8L)
  expect_equal(row$df_num, 2)
  expect_equal(row$df_den, 14)
  expect_equal(row$statistic, reference["time_factor", "F value"], tolerance = 1e-10)
  expect_equal(row$p_value, stats::pf(row$statistic, 2, 14, lower.tail = FALSE),
               tolerance = 1e-12)
  expect_equal(row$partial_eta_squared, ss_time / (ss_time + ss_error), tolerance = 1e-10)
  expect_equal(row$generalized_eta_squared,
               ss_time / (ss_time + ss_error + ss_subject), tolerance = 1e-10)
  expect_true(row$df_num_gg > 0 && row$df_num_gg <= row$df_num)
  expect_true(row$df_den_gg > 0 && row$df_den_gg <= row$df_den)
  expect_equal(row$df_num_gg / row$df_num, row$df_den_gg / row$df_den, tolerance = 1e-10)
  expect_equal(row$p_gg, stats::pf(row$statistic, row$df_num_gg, row$df_den_gg,
                                  lower.tail = FALSE), tolerance = 1e-10)
  expect_null(result$error)
})

test_that("independence GEE agrees with OLS and an independently formed cluster sandwich", {
  skip_if_not_installed("geepack")
  data <- model_oracle_long()
  result <- MIRA:::.mira_run_gee_long(data, value ~ time_factor, confidence_level = 0.9)
  model <- result$independence
  reference <- stats::lm(value ~ time_factor, data)
  design <- stats::model.matrix(reference)
  bread <- solve(crossprod(design))
  cluster_scores <- rowsum(design * stats::residuals(reference), data$patient_factor)
  covariance <- bread %*% crossprod(cluster_scores) %*% bread
  se <- sqrt(diag(covariance))
  beta <- stats::coef(reference)
  critical <- stats::qnorm(0.95)

  expect_true(model$performed)
  expect_identical(model$coefficients$term, names(beta))
  expect_equal(model$coefficients$Estimate, unname(beta), tolerance = 1e-10)
  expect_equal(model$coefficients$robust_se, unname(se), tolerance = 1e-10)
  expect_equal(model$coefficients$conf_low, unname(beta - critical * se), tolerance = 1e-10)
  expect_equal(model$coefficients$conf_high, unname(beta + critical * se), tolerance = 1e-10)
  time_indices <- 2:3
  time_statistic <- as.numeric(crossprod(beta[time_indices],
    solve(covariance[time_indices, time_indices], beta[time_indices])))
  test <- model$effect_tests[model$effect_tests$question == "TIME", ]
  expect_identical(nrow(test), 1L)
  expect_identical(test$df, 2L)
  expect_equal(test$statistic, time_statistic, tolerance = 1e-10)
  expect_equal(test$p_raw, stats::pchisq(time_statistic, 2, lower.tail = FALSE), tolerance = 1e-10)
  expect_false(any(model$effect_tests$question %in% c("ARM", "TIME_X_ARM")))
  expect_identical(stats::nobs(model$model), 24L)
  expect_null(model$error)
})

test_that("emmeans contrast families preserve arithmetic estimands and ordinal trends", {
  skip_if_not_installed("emmeans")
  data <- model_oracle_long()
  fitted <- stats::lm(value ~ time_factor, data)
  result <- MIRA:::.mira_run_emmeans_long(fitted, FALSE, "BH", confidence_level = 0.9)
  means <- colMeans(model_oracle_measurements())
  expect_true(result$performed)
  expect_equal(result$time$tidy$emmean, means, tolerance = 1e-12)
  expect_equal(result$time$pairwise_tidy$estimate,
               c(means[[1L]] - means[[2L]], means[[1L]] - means[[3L]],
                 means[[2L]] - means[[3L]]), tolerance = 1e-12)
  expect_equal(result$baseline_followup$tidy$estimate, means[-1L] - means[[1L]],
               tolerance = 1e-12)
  expect_equal(result$baseline_final$tidy$estimate, means[[3L]] - means[[1L]],
               tolerance = 1e-12)
  expect_equal(result$consecutive$tidy$estimate, diff(means), tolerance = 1e-12)
  expect_identical(as.character(result$trends$tidy$contrast), c("linear", "quadratic"))
  expect_equal(result$trends$tidy$estimate,
               c(means[[3L]] - means[[1L]], means[[1L]] - 2 * means[[2L]] + means[[3L]]),
               tolerance = 1e-12)
  expect_identical(result$trends$scale, "ordered timepoint levels")
  # Irregular Month 4 / Month 12 labels still use ordinal polynomial contrasts.
  expect_match(result$trends$note, "do not assume")
  for (family in list(result$time$pairwise_tidy, result$baseline_followup$tidy,
                      result$baseline_final$tidy, result$consecutive$tidy,
                      result$trends$tidy)) {
    expect_equal(family$p_primary, stats::p.adjust(family$p_raw, "BH"), tolerance = 1e-12)
    expect_equal(family$lower.CL, family$estimate - stats::qt(0.95, family$df) * family$SE,
                 tolerance = 1e-10)
    expect_equal(family$upper.CL, family$estimate + stats::qt(0.95, family$df) * family$SE,
                 tolerance = 1e-10)
    expect_equal(family$p_raw, 2 * stats::pt(-abs(family$estimate / family$SE), family$df),
                 tolerance = 1e-10)
  }
  expect_null(result$arm$tidy)
  expect_null(result$arm_time$tidy)
})

test_that("public mixed-model estimates match complete balanced cell means and ML likelihood tests", {
  skip_if_not_installed("lme4")
  y <- model_oracle_measurements()
  data <- data.frame(subject_id = seq_len(nrow(y)), score_t0 = y[, 1L],
                     score_t1 = y[, 2L], score_t2 = y[, 3L])
  result <- suppressMessages(mira_info_long(data, id = "subject_id", outcomes = "score",
    covariates = character(0), analyses = "model", verbose = FALSE))
  fitted <- result$model$fitted_model
  expect_true(inherits(fitted, "merMod"))
  expect_null(result$model$error)
  expect_identical(stats::nobs(fitted), 24L)
  expect_equal(unname(lme4::fixef(fitted)),
               c(mean(y[, 1L]), mean(y[, 2L] - y[, 1L]), mean(y[, 3L] - y[, 1L])),
               tolerance = 1e-8)
  # Independently fit the nested models by ML; the stored primary model uses REML.
  model_data <- stats::model.frame(fitted)
  null_ml <- lme4::lmer(value ~ 1 + (1 | patient_factor), data = model_data, REML = FALSE)
  full_ml <- lme4::lmer(value ~ time_factor + (1 | patient_factor),
                        data = model_data, REML = FALSE)
  statistic <- 2 * (as.numeric(stats::logLik(full_ml)) - as.numeric(stats::logLik(null_ml)))
  lrt <- as.data.frame(result$model$global_time_test)
  expect_equal(tail(lrt$Chisq, 1L), statistic, tolerance = 1e-8)
  expect_equal(tail(lrt$`Pr(>Chisq)`, 1L), stats::pchisq(statistic, 2, lower.tail = FALSE),
               tolerance = 1e-8)
  variances <- as.data.frame(lme4::VarCorr(fitted))
  subject <- variances$vcov[variances$grp == "patient_factor"]
  residual <- variances$vcov[variances$grp == "Residual"]
  expect_equal(result$variability$ICC_model, subject / (subject + residual), tolerance = 1e-12)
  expect_equal(result$variability$ICC, result$variability$ICC_model)
  comparison <- result$advanced_models$model_comparison
  primary <- comparison[comparison$model == "mixed_random_intercept", ]
  expect_equal(primary$AIC, stats::AIC(fitted), tolerance = 1e-12)
  expect_equal(primary$BIC, stats::BIC(fitted), tolerance = 1e-12)
  expect_equal(primary$logLik, as.numeric(stats::logLik(fitted)), tolerance = 1e-12)
  expect_identical(primary$n_subjects, 8L)
})

test_that("available emmeans records incompatible-model errors without a fabricated estimate", {
  skip_if_not_installed("emmeans")
  data <- model_oracle_long()
  result <- MIRA:::.mira_run_emmeans_long(stats::lm(value ~ 1, data), FALSE, "holm", 0.95)
  expect_false(result$performed)
  expect_identical(result$package, "emmeans")
  expect_match(result$error, "time_factor")
  expect_identical(result$errors$time, result$error)
  expect_match(result$reason_skipped, "could not be estimated")
  expect_null(result$object)
  expect_null(result$tidy)
})

test_that("available GEE engines record estimation failures separately for every working correlation", {
  skip_if_not_installed("geepack")
  result <- MIRA:::.mira_run_gee_long(model_oracle_long(), value ~ missing_covariate, 0.95)
  expect_false(result$performed)
  for (structure in c("independence", "exchangeable", "ar1")) {
    item <- result[[structure]]
    expect_false(item$performed)
    expect_false(item$converged)
    expect_identical(item$correlation_structure, structure)
    expect_match(item$error, "missing_covariate")
    expect_match(item$reason_skipped, "estimation failed")
    expect_null(item$model)
    expect_null(item$coefficients)
    expect_null(item$qic)
  }
})
