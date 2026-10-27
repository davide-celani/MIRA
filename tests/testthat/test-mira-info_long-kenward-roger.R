# Kenward-Roger complements asymptotic LRT and cluster-robust CR2 inference.
# The arithmetic oracles below concern balanced Gaussian random-intercept
# models. They do not assert equality across different covariance models.
kr_oracle_time_data <- function() {
  y <- rbind(c(3, 5, 8), c(7, 6, 10), c(4, 8, 7), c(9, 12, 11),
    c(6, 9, 12), c(12, 11, 15), c(8, 13, 16), c(15, 14, 18))
  data.frame(patient_factor = factor(rep(paste0("P", 1:8), each = 3L)),
    time_factor = factor(rep(c("Base", "Mid", "Final"), 8L),
      levels = c("Base", "Mid", "Final")), value = as.vector(t(y)))
}

kr_oracle_arm_data <- function() {
  baseline <- c(3, 5, 4, 7, 6, 8, 5, 9, 7, 10, 8, 11)
  amplitude <- c(0, 1, -1, 0, -1, 1, 0, 2, 1, 0, -2, -1)
  arm <- rep(c("control", "active"), each = 6L)
  data <- data.frame(patient_factor = factor(rep(sprintf("P%02d", 1:12), each = 3L)),
    arm_factor = factor(rep(arm, each = 3L), levels = c("control", "active")),
    time_factor = factor(rep(c("Base", "Mid", "Final"), 12L),
      levels = c("Base", "Mid", "Final")),
    value = rep(baseline, each = 3L) + rep(c(0, 1, 3), 12L) +
      rep(amplitude, each = 3L) * rep(c(0, 1, -1), 12L) +
      rep(as.numeric(arm == "active"), each = 3L) * rep(c(0, 2, 1), 12L))
  data
}

kr_oracle_fit <- function(data = kr_oracle_time_data(), arm = FALSE, covariate = FALSE) {
  formula <- if (arm) {
    if (covariate) value ~ time_factor * arm_factor + age_time_factor + (1 | patient_factor)
    else value ~ time_factor * arm_factor + (1 | patient_factor)
  } else value ~ time_factor + (1 | patient_factor)
  lme4::lmer(formula, data = data, REML = TRUE)
}

test_that("balanced time KR inference matches independent within-subject sums of squares", {
  data <- kr_oracle_time_data()
  model <- kr_oracle_fit(data)
  result <- MIRA:::.mira_run_kenward_roger_long(model, FALSE)
  raw <- matrix(data$value, nrow = 8L, byrow = TRUE)
  # Orthogonal subject and visit effects: SS_time = n*sum((mean_t-grand)^2),
  # SS_error from removing both effects; df = 2 and (8-1)*(3-1)=14.
  residual <- sweep(sweep(raw, 1L, rowMeans(raw)), 2L, colMeans(raw)) + mean(raw)
  ss_time <- 8 * sum((colMeans(raw) - mean(raw))^2)
  ss_error <- sum(residual^2)
  expected_f <- (ss_time / 2) / (ss_error / 14)
  expect_true(result$performed)
  expect_identical(result$tidy$question, "TIME")
  expect_identical(result$tidy$n, 24L)
  expect_identical(result$tidy$n_subjects, 8L)
  expect_equal(result$tidy$statistic, expected_f, tolerance = 1e-7)
  expect_equal(result$tidy$df_num, 2, tolerance = 1e-12)
  expect_equal(result$tidy$df_den, 14, tolerance = 1e-7)
  expect_equal(result$tidy$F_scaling, 1, tolerance = 1e-12)
  expect_equal(result$tidy$p_raw, stats::pf(expected_f, 2, 14, lower.tail = FALSE),
    tolerance = 1e-7)
  expect_s3_class(result$tests$TIME, "KRmodcomp")
  expect_identical(result$effect_diagnostics$TIME$tested_coefficients,
    c("time_factorMid", "time_factorFinal"))
  expect_false(result$effect_diagnostics$TIME$hypothesis_reduced)
  expect_true(result$effect_diagnostics$TIME$estimable)
  expect_null(result$reason_skipped)
  expect_null(result$error)
  expect_length(result$errors, 0L)
  expect_identical(result$warnings, character(0))
})

test_that("balanced joint time and interaction KR hypotheses match factorial arithmetic", {
  data <- kr_oracle_arm_data()
  model <- kr_oracle_fit(data, arm = TRUE)
  result <- MIRA:::.mira_run_kenward_roger_long(model, TRUE)
  expect_true(result$performed)
  expect_identical(result$tidy$question, c("TIME", "ARM", "TIME_X_ARM"))
  expect_identical(anyDuplicated(result$tidy$question), 0L)
  # Twelve subjects, six per arm. After subject means and cell means are
  # removed, RSS=28 on (12-2)*(3-1)=20 df. The null of no change in either arm
  # has RSS=108; common time with no interaction has RSS=34. Thus joint TIME
  # is (108-28)/4/(28/20)=100/7 and interaction is (34-28)/2/(28/20)=15/7.
  time <- result$tidy[result$tidy$question == "TIME", , drop = FALSE]
  interaction <- result$tidy[result$tidy$question == "TIME_X_ARM", , drop = FALSE]
  expect_identical(nrow(time), 1L)
  expect_identical(nrow(interaction), 1L)
  expect_equal(time$statistic, 100 / 7, tolerance = 1e-7)
  expect_equal(c(time$df_num, time$df_den), c(4, 20), tolerance = 1e-7)
  expect_equal(time$p_raw, stats::pf(100 / 7, 4, 20, lower.tail = FALSE), tolerance = 1e-7)
  expect_equal(interaction$statistic, 15 / 7, tolerance = 1e-7)
  expect_equal(c(interaction$df_num, interaction$df_den), c(2, 20), tolerance = 1e-7)
  expect_equal(interaction$p_raw, stats::pf(15 / 7, 2, 20, lower.tail = FALSE), tolerance = 1e-7)
  expect_equal(time$F_scaling, 1, tolerance = 1e-12)
  expect_equal(interaction$F_scaling, 1, tolerance = 1e-12)
  expect_identical(result$n, 36L)
  expect_identical(result$n_subjects, 12L)
})

test_that("KR keeps the fitted sample and explicit term-aware hypotheses with missing cells", {
  data <- kr_oracle_arm_data()
  data$age_time_factor <- rep(c(31, 44, 29, 58, 36, 47, 33, 52, 40, 61, 37, 49), each = 3L)
  data$age_time_factor[data$patient_factor == "P02"] <- NA_real_
  data$value[data$patient_factor == "P05" & data$time_factor == "Mid"] <- NA_real_
  data$value[data$patient_factor == "P09" & data$time_factor == "Final"] <- NA_real_
  model <- kr_oracle_fit(data, arm = TRUE, covariate = TRUE)
  result <- MIRA:::.mira_run_kenward_roger_long(model, TRUE)
  expected <- data[data$patient_factor != "P02" &
    !(data$patient_factor == "P05" & data$time_factor == "Mid") &
    !(data$patient_factor == "P09" & data$time_factor == "Final"), , drop = FALSE]
  key <- function(x) paste(x$patient_factor, x$time_factor, sep = "::")
  expect_identical(result$n, 31L)
  expect_identical(result$n_subjects, 11L)
  expect_identical(key(result$data), key(expected))
  expect_identical(anyDuplicated(key(result$data)), 0L)
  expect_identical(result$data$value, expected$value)
  expect_identical(result$data$age_time_factor, expected$age_time_factor)
  expect_identical(result$data, stats::model.frame(model))
  beta <- lme4::fixef(model)
  # Explicit coefficient identities are an engine-wrapper oracle. In
  # particular the nuisance name contains time_factor but is never TIME.
  targets <- list(TIME = c("time_factorMid", "time_factorFinal",
      "time_factorMid:arm_factoractive", "time_factorFinal:arm_factoractive"),
    ARM = c("arm_factoractive", "time_factorMid:arm_factoractive",
      "time_factorFinal:arm_factoractive"),
    TIME_X_ARM = c("time_factorMid:arm_factoractive", "time_factorFinal:arm_factoractive"))
  for (question in names(targets)) {
    row <- result$tidy[result$tidy$question == question, , drop = FALSE]
    expect_identical(nrow(row), 1L)
    restriction <- diag(length(beta))[match(targets[[question]], names(beta)), , drop = FALSE]
    dimnames(restriction) <- list(targets[[question]], names(beta))
    expect_identical(result$constraints[[question]], restriction)
    expect_identical(result$effect_diagnostics[[question]]$tested_coefficients,
      targets[[question]])
    reference <- pbkrtest::getKR(pbkrtest::KRmodcomp(model, restriction))
    expect_equal(c(row$statistic, row$df_num, row$df_den, row$F_scaling, row$p_raw),
      unname(unlist(reference[c("Fstat", "ndf", "ddf", "F.scaling", "p.value")])),
      tolerance = 1e-12, info = question)
  }
})

test_that("rank-deficient KR hypotheses expose reduced and unavailable coefficient sets", {
  for (all_missing in c(FALSE, TRUE)) {
    data <- kr_oracle_arm_data()
    missing <- data$arm_factor == "active" &
      if (all_missing) data$time_factor != "Base" else data$time_factor == "Final"
    data$value[missing] <- NA_real_
    expect_message(model <- kr_oracle_fit(data, arm = TRUE), "rank deficient")
    result <- MIRA:::.mira_run_kenward_roger_long(model, TRUE)
    expect_true(result$performed)
    diagnostic <- result$effect_diagnostics$TIME_X_ARM
    requested <- c("time_factorMid:arm_factoractive", "time_factorFinal:arm_factoractive")
    expect_identical(diagnostic$requested_coefficients, requested)
    expect_true(diagnostic$hypothesis_reduced)
    expect_match(diagnostic$note, "Dropped coefficients", fixed = TRUE)
    expect_identical(diagnostic$dropped_coefficients,
      if (all_missing) requested else requested[[2L]])
    row <- result$tidy[result$tidy$question == "TIME_X_ARM", , drop = FALSE]
    expect_identical(nrow(row), 1L)
    if (all_missing) {
      expect_false(diagnostic$estimable)
      expect_identical(diagnostic$tested_coefficients, character(0))
      expect_match(diagnostic$reason_skipped, "No retained coefficient", fixed = TRUE)
      expect_null(result$tests$TIME_X_ARM)
      expect_null(result$constraints$TIME_X_ARM)
      expect_true(all(is.na(row[c("statistic", "df_num", "df_den", "F_scaling", "p_raw")])))
    } else {
      expect_true(diagnostic$estimable)
      expect_identical(diagnostic$tested_coefficients, requested[[1L]])
      expect_null(diagnostic$reason_skipped)
      beta <- lme4::fixef(model)
      restriction <- diag(length(beta))[match(requested[[1L]], names(beta)), , drop = FALSE]
      reference <- pbkrtest::getKR(pbkrtest::KRmodcomp(model, restriction))
      expect_equal(row$statistic, reference$Fstat, tolerance = 1e-12)
      expect_equal(row$p_raw, reference$p.value, tolerance = 1e-12)
      expect_equal(row$df_num, 1, tolerance = 0)
    }
  }
})

test_that("KR prerequisite guards retain sample evidence without invoking the engine", {
  model <- kr_oracle_fit()
  called <- 0L
  results <- testthat::with_mocked_bindings(
    list(absent = MIRA:::.mira_run_kenward_roger_long(NULL, FALSE),
      wrong = MIRA:::.mira_run_kenward_roger_long(stats::lm(mpg ~ wt, mtcars), FALSE),
      not_converged = MIRA:::.mira_run_kenward_roger_long(model, FALSE, converged = FALSE),
      singular = MIRA:::.mira_run_kenward_roger_long(model, FALSE, singular = TRUE)),
    KRmodcomp = function(...) { called <<- called + 1L; stop("Unexpected KR engine invocation") },
    .package = "pbkrtest")
  expect_identical(called, 0L)
  for (item in results) {
    expect_false(item$performed)
    expect_identical(nrow(item$tidy), 0L)
    expect_null(item$error)
  }
  expect_match(results$absent$reason_skipped, "not available")
  expect_match(results$wrong$reason_skipped, "Gaussian lmer")
  expect_match(results$not_converged$reason_skipped, "convergence")
  expect_match(results$singular$reason_skipped, "singular")
  expect_identical(results$not_converged$data, stats::model.frame(model))
  expect_identical(results$singular$n_subjects, 8L)
  original_require <- base::requireNamespace
  absent <- testthat::with_mocked_bindings(
    MIRA:::.mira_run_kenward_roger_long(model, FALSE),
    requireNamespace = function(package, ...) {
      if (identical(package, "pbkrtest")) FALSE else original_require(package, ...)
    }, .package = "base")
  expect_false(absent$performed)
  expect_match(absent$reason_skipped, "pbkrtest.*not installed")
  expect_identical(nrow(absent$tidy), 0L)
})

test_that("KR engine failure, nonfinite output and a failed question remain distinguishable", {
  model <- kr_oracle_fit()
  failed <- testthat::with_mocked_bindings(
    MIRA:::.mira_run_kenward_roger_long(model, FALSE),
    KRmodcomp = function(...) {
      warning("controlled KR engine warning")
      stop("controlled KR engine error")
    }, .package = "pbkrtest")
  expect_false(failed$performed)
  expect_identical(failed$error, "controlled KR engine error")
  expect_identical(failed$errors$TIME, "controlled KR engine error")
  expect_identical(failed$warnings, "controlled KR engine warning")
  expect_identical(failed$tidy$question, "TIME")
  expect_true(is.na(failed$tidy$p_raw))
  expect_match(failed$effect_diagnostics$TIME$reason_skipped, "estimation failed")

  nonfinite <- testthat::with_mocked_bindings(
    MIRA:::.mira_run_kenward_roger_long(model, FALSE),
    KRmodcomp = function(...) structure(list(stats = list(Fstat = Inf, ndf = 2,
      ddf = 0, F.scaling = 1, p.value = NaN)), class = "KRmodcomp"),
    .package = "pbkrtest")
  expect_false(nonfinite$performed)
  expect_s3_class(nonfinite$tests$TIME, "KRmodcomp")
  expect_identical(nonfinite$tests$TIME$stats$Fstat, Inf)
  expect_true(all(is.na(nonfinite$tidy[c("statistic", "df_num", "df_den", "F_scaling", "p_raw")])))
  expect_null(nonfinite$error)
  expect_length(nonfinite$errors, 0L)
  expect_match(nonfinite$effect_diagnostics$TIME$note, "no finite estimable", fixed = TRUE)

  model <- kr_oracle_fit(kr_oracle_arm_data(), arm = TRUE)
  original <- pbkrtest::KRmodcomp
  partial <- testthat::with_mocked_bindings(
    MIRA:::.mira_run_kenward_roger_long(model, TRUE),
    KRmodcomp = function(largeModel, smallModel, ...) {
      if (nrow(smallModel) == 3L) stop("controlled ARM-only failure")
      original(largeModel, smallModel, ...)
    }, .package = "pbkrtest")
  expect_true(partial$performed)
  expect_identical(partial$tidy$question, c("TIME", "ARM", "TIME_X_ARM"))
  expect_identical(partial$errors$ARM, "controlled ARM-only failure")
  expect_null(partial$error)
  expect_false(partial$effect_diagnostics$ARM$estimable)
  expect_true(is.na(partial$tidy$p_raw[partial$tidy$question == "ARM"]))
  expect_true(all(is.finite(partial$tidy$p_raw[partial$tidy$question != "ARM"])))
  expect_s3_class(partial$tests$TIME, "KRmodcomp")
  expect_s3_class(partial$tests$TIME_X_ARM, "KRmodcomp")
})

test_that("a nuisance name does not invent a KR scientific time hypothesis", {
  data <- kr_oracle_time_data()
  data$age_time_factor <- rep(c(31, 43, 29, 52, 37, 46, 33, 58), each = 3L)
  model <- lme4::lmer(value ~ age_time_factor + (1 | patient_factor), data = data)
  result <- MIRA:::.mira_run_kenward_roger_long(model, FALSE)
  expect_false(result$performed)
  expect_identical(nrow(result$tidy), 0L)
  expect_length(result$constraints, 0L)
  expect_length(result$tests, 0L)
  expect_match(result$reason_skipped, "No time or arm hypothesis")
  expect_identical(result$n, 24L)
})

test_that("KR canonical metric names preserve native evidence and structural sample keys", {
  result <- MIRA:::.mira_run_kenward_roger_long(kr_oracle_fit(), FALSE)
  raw <- list(advanced_tests = list(kenward_roger = result))
  bases <- c("statistic", "df_num", "df_den", "F_scaling", "p_raw")
  for (declared in list(NULL, "Quantitative score")) {
    output <- MIRA:::.mira_finalize_output_long(raw, declared)
    module <- output$advanced_tests$kenward_roger
    expected <- if (is.null(declared)) bases else paste0(declared, "_", bases)
    expect_identical(unname(attr(module$tidy, "mira_metric_names")), expected)
    expect_true(all(expected %in% names(module$tidy)))
    expect_identical(module$tidy$question, "TIME")
    expect_identical(module$tidy$n, 24L)
    expect_identical(module$tidy$n_subjects, 8L)
    expect_identical(module$tests, result$tests)
    expect_identical(module$constraints, result$constraints)
    expect_identical(module$data, result$data)
    for (i in seq_along(bases)) {
      expect_identical(module$tidy[[expected[i]]], result$tidy[[bases[i]]])
    }
  }
})
