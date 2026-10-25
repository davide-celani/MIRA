# Independent arithmetic supplements the existing reference-engine comparisons.
# No t.test(), p.adjust(), production helper or model fit supplies these expected
# Welch quantities. The 1e-12 tolerance covers only floating-point arithmetic.
audit_welch_long <- function(a, b, alpha = 0.1) {
  n_a <- length(a)
  n_b <- length(b)
  if (min(n_a, n_b) < 2L) {
    return(c(estimate = NA_real_, lower = NA_real_, upper = NA_real_,
             p = NA_real_, g = NA_real_))
  }
  center_a <- sum(a) / n_a
  center_b <- sum(b) / n_b
  variance_a <- sum((a - center_a)^2) / (n_a - 1L)
  variance_b <- sum((b - center_b)^2) / (n_b - 1L)
  variance_mean_a <- variance_a / n_a
  variance_mean_b <- variance_b / n_b
  se <- sqrt(variance_mean_a + variance_mean_b)
  df <- (variance_mean_a + variance_mean_b)^2 /
    (variance_mean_a^2 / (n_a - 1L) + variance_mean_b^2 / (n_b - 1L))
  estimate <- center_b - center_a
  half_width <- stats::qt(1 - alpha / 2, df) * se
  pooled_df <- n_a + n_b - 2L
  pooled_variance <- ((n_a - 1L) * variance_a +
    (n_b - 1L) * variance_b) / pooled_df
  correction <- 1 - 3 / (4 * pooled_df - 1)
  c(estimate = estimate, lower = estimate - half_width,
    upper = estimate + half_width,
    p = 2 * stats::pt(-abs(estimate / se), df),
    g = correction * estimate / sqrt(pooled_variance))
}

audit_holm_long <- function(p) {
  result <- rep(NA_real_, length(p))
  finite <- which(is.finite(p))
  order_p <- finite[order(p[finite])]
  result[order_p] <- pmin(1, cummax(p[order_p] * rev(seq_along(order_p))))
  result
}

test_that("asymmetric arm samples have independent Welch, Hedges and finite-family Holm oracles", {
  # C4 has baseline only, A5 lacks baseline, X has an unknown arm. Their
  # extreme values distinguish time-specific, paired-change and core samples.
  # Mid has two active observations but only one active baseline pair, so its
  # Welch change test is unavailable and must not enlarge the finite family.
  data <- data.frame(
    subject_id = c(paste0("C", 1:4), paste0("A", 1:5), "X"),
    arm = c(rep("control", 4L), rep("active", 5L), NA_character_),
    score_t0 = c(1, 2, 5, 80, 3, 7, 9, 13, NA, 1000),
    score_t1 = c(2, 4, 9, NA, 5, NA, NA, NA, -40, 1500),
    score_t2 = c(4, 1, 10, NA, 5, 12, 10, 19, -40, 2000),
    stringsAsFactors = FALSE
  )
  expect_warning(result <- mira_test_info_base(data, id = "subject_id", outcomes = "score",
    time_vars = c("score_t0", "score_t1", "score_t2"), arm = "arm",
    reference_arm = "control", covariates = character(0), analyses = "arm_tests",
    alpha = 0.1, p_adjust_method = "holm", verbose = FALSE),
    "1 subjects have a missing treatment arm", fixed = TRUE)

  selected_time <- list(
    list(control = paste0("C", 1:4), active = paste0("A", 1:4)),
    list(control = paste0("C", 1:3), active = c("A1", "A5")),
    list(control = paste0("C", 1:3), active = paste0("A", 1:5))
  )
  selected_change <- list(
    list(control = paste0("C", 1:3), active = "A1"),
    list(control = paste0("C", 1:3), active = paste0("A", 1:4))
  )
  expected_time_values <- list(
    list(a = c(1, 2, 5, 80), b = c(3, 7, 9, 13)),
    list(a = c(2, 4, 9), b = c(5, -40)),
    list(a = c(4, 1, 10), b = c(5, 12, 10, 19, -40))
  )
  expected_change_values <- list(
    list(a = c(1, 2, 4), b = 2),
    list(a = c(3, -1, 5), b = c(2, 5, 1, 6))
  )

  # The public long data preserve the subjects/visits from which each table
  # is constructed; these explicit key sets distinguish identical N with
  # different observations, in addition to the independent numeric expected.
  long <- result$long_data
  for (i in seq_along(selected_time)) {
    for (group in c("control", "active")) {
      keep <- which(long$time == result$time_vars[[i]] &
        long$arm == group & is.finite(long$value))
      expect_setequal(as.character(long$patient[keep]), selected_time[[i]][[group]])
    }
  }
  for (i in seq_along(selected_change)) {
    for (group in c("control", "active")) {
      baseline <- long[which(long$time == "score_t0" & long$arm == group &
        is.finite(long$value)), "patient"]
      followup <- long[which(long$time == result$time_vars[[i + 1L]] &
        long$arm == group & is.finite(long$value)), "patient"]
      expect_setequal(intersect(as.character(baseline), as.character(followup)),
        selected_change[[i]][[group]])
    }
  }
  expect_identical(result$descriptives$n, c(9L, 6L, 9L))
  expect_identical(result$arm_analysis$time_omnibus$n, c(8L, 5L, 8L))

  for (kind in c("time", "change")) {
    values <- if (kind == "time") expected_time_values else expected_change_values
    table <- result$arm_analysis[[paste0(kind, "_pairwise")]]
    expected <- t(vapply(values, function(pair) audit_welch_long(pair$a, pair$b),
                         numeric(5L)))
    estimate_column <- if (kind == "time") "mean_difference_b_minus_a" else {
      "difference_in_change_b_minus_a"
    }
    expect_identical(table$arm_a, rep("control", length(values)))
    expect_identical(table$arm_b, rep("active", length(values)))
    expect_equal(table$n_a, vapply(values, function(pair) length(pair$a), integer(1L)))
    expect_equal(table$n_b, vapply(values, function(pair) length(pair$b), integer(1L)))
    observed <- as.matrix(table[, c(estimate_column, "ci_lower", "ci_upper",
                                    "welch_t_p", "hedges_g")])
    expect_equal(unname(observed), unname(expected), tolerance = 1e-12)
    expect_equal(table$welch_t_p_adj, audit_holm_long(expected[, "p"]), tolerance = 1e-12)
  }

  # Demonstrate that the fixture rejects plausible selection/denominator errors.
  correct <- audit_welch_long(c(3, -1, 5), c(2, 5, 1, 6))
  wrong_available_rows <- audit_welch_long(c(3, -1, 5, -80), c(2, 5, 1, 6, -40))
  wrong_unknown_arm <- audit_welch_long(c(3, -1, 5), c(2, 5, 1, 6, 1000))
  expect_gt(abs(correct[["estimate"]] - wrong_available_rows[["estimate"]]), 1)
  expect_gt(abs(correct[["estimate"]] - wrong_unknown_arm[["estimate"]]), 1)
  time_p <- result$arm_analysis$time_pairwise$welch_t_p
  expect_gt(max(abs(audit_holm_long(time_p) - time_p)), 0.01)
  change_p <- result$arm_analysis$change_pairwise$welch_t_p
  expect_identical(is.finite(change_p), c(FALSE, TRUE))
  expect_equal(result$arm_analysis$change_pairwise$welch_t_p_adj[[2L]], change_p[[2L]])
})

test_that("mixed and independence GEE models retain the exact complete ID and visit keys", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("geepack")
  # Independent sample oracle: S02 has blank site, S09 infinite age, S11 blank
  # arm; S03 misses Mid and S08 has only Mid. Keeping a single-visit subject
  # and an incomplete profile distinguishes long-row selection from complete
  # profiles. Outcome values and covariates are asymmetric across subjects.
  data <- data.frame(subject_id = sprintf("S%02d", 1:12),
    arm = rep(c("control", "active"), each = 6L),
    age = c(21, 34, 27, 45, 38, 52, 23, 36, Inf, 47, 40, 55),
    site = rep(c("N", "S"), 6L),
    score_t0 = c(3, 5, 4, 7, 6, 8, 5, NA, 7, 10, 8, 11),
    score_t1 = c(4, 6, NA, 6, 8, 9, 7, 10, 8, 11, 9, 12),
    score_t2 = c(6, 7, 7, 8, 9, 11, 8, NA, 10, 13, 11, 14),
    stringsAsFactors = FALSE)
  data$site[[2L]] <- " "
  data$arm[[11L]] <- " "
  expect_warning(result <- mira_test_info_base(data, id = "subject_id", outcomes = "score",
    time_vars = c("score_t0", "score_t1", "score_t2"),
    time_labels = c("Base", "Mid", "Final"), arm = "arm", reference_arm = "control",
    covariates = c("age", "site"), analyses = c("model", "arm_tests"), verbose = FALSE),
    "1 subjects have a missing treatment arm", fixed = TRUE)
  expected <- c(
    "S01::Base", "S01::Mid", "S01::Final", "S03::Base", "S03::Final",
    "S04::Base", "S04::Mid", "S04::Final", "S05::Base", "S05::Mid", "S05::Final",
    "S06::Base", "S06::Mid", "S06::Final", "S07::Base", "S07::Mid", "S07::Final",
    "S08::Mid", "S10::Base", "S10::Mid", "S10::Final", "S12::Base", "S12::Mid",
    "S12::Final")
  fitted <- result$model$fitted_model
  expect_s4_class(fitted, "merMod")
  frame <- stats::model.frame(fitted)
  observed <- paste(frame$patient_factor, frame$time_factor, sep = "::")
  expect_setequal(observed, expected)
  expect_identical(anyDuplicated(observed), 0L)
  expect_identical(stats::nobs(fitted), 24L)
  expect_identical(result$model$covariates_used, c("age", "site"))
  expect_identical(result$overview$n_patients, 12L)
  expect_identical(nrow(result$long_data), 36L)
  gee <- result$advanced_models$gee$independence
  expect_true(gee$performed)
  expect_true(gee$converged)
  expect_identical(stats::nobs(gee$model), 24L)
  gee_frame <- stats::model.frame(gee$model)
  # geeglm retains the supplied subject IDs in its id vector; the fixed model
  # frame retains the corresponding time factor and response in the same order.
  gee_keys <- paste(as.character(gee$model$id), gee_frame$time_factor, sep = "::")
  expect_setequal(gee_keys, expected)
  expect_identical(anyDuplicated(gee_keys), 0L)
  # Count-only checks would miss swapping S03's unavailable Mid for Final.
  wrong_keys <- sub("S03::Final", "S03::Mid", expected, fixed = TRUE)
  expect_length(setdiff(expected, wrong_keys), 1L)
  expect_length(setdiff(wrong_keys, expected), 1L)
  comparison <- result$advanced_models$model_comparison
  rows <- comparison$model %in% c("mixed_random_intercept", "gee_independence")
  expect_identical(comparison$n[rows], c(24L, 24L))
  expect_identical(comparison$n_subjects[rows], c(9L, 9L))
})

test_that("formula quoting preserves literal backticks and backslashes in selected covariates", {
  # Legal R names must resolve to the selected columns, including two names
  # that would collide if a literal backtick were removed. The expected design
  # is explicit numeric columns, independent of the production formula builder.
  data <- data.frame(value = c(3, 5, 7, 4, 8, 10),
    time_factor = factor(rep(c("Base", "Final"), 3L), levels = c("Base", "Final")),
    patient_factor = factor(rep(c("P1", "P2", "P3"), each = 2L)),
    check.names = FALSE)
  data[["age`years"]] <- c(21, 21, 35, 35, 48, 48)
  data[["ageyears"]] <- c(5, 5, 8, 8, 2, 2)
  data[["dose\\level"]] <- c(3, 3, 1, 1, 9, 9)
  selections <- list(c("age`years", "ageyears"), "dose\\level",
                     c("age`years", "ageyears", "dose\\level"))
  for (variables in selections) {
    expected <- cbind(1, rep(c(0, 1), 3L), as.matrix(data[variables]))
    for (builder in c("fixed", "mixed")) {
      formula <- tryCatch(
        if (builder == "fixed") MIRA:::.mira_make_fixed_formula_long(FALSE, variables) else {
          MIRA:::.mira_make_lmer_formula_long("time_factor", variables, "(1 | patient_factor)")
        }, error = identity)
      expect_false(inherits(formula, "error"),
        info = paste("literal column names must parse:", builder, paste(variables, collapse = ", ")))
      if (inherits(formula, "error")) next
      expect_setequal(all.vars(formula), c("value", "time_factor", variables,
        if (builder == "mixed") "patient_factor"))
      if (builder == "fixed") {
        observed <- stats::model.matrix(formula, data)
        expect_identical(dim(observed), dim(expected))
        expect_equal(as.numeric(observed), as.numeric(expected), tolerance = 0)
      }
    }
  }
  expect_identical(MIRA:::.mira_quote_formula_name_long("age"), "`age`")
  expect_identical(MIRA:::.mira_quote_formula_name_long("age years"), "`age years`")
})

test_that("literal covariate names retain both adjustment columns in primary and advanced models", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("geepack")
  # A consistent column rename must leave this fixed Gaussian specification
  # and its design unchanged. No engine serves as an independent inferential
  # oracle here: explicit model columns prove specification, and coefficients
  # compare the wrapper before/after renaming at fitting precision (1e-10).
  y <- rbind(c(3, 5, 8), c(7, 6, 10), c(4, 8, 7), c(9, 12, 11),
    c(6, 9, 12), c(12, 11, 15), c(8, 13, 16), c(15, 14, 18))
  data <- data.frame(subject_id = paste0("P", 1:8), score_t0 = y[, 1L],
    score_t1 = y[, 2L], score_t2 = y[, 3L],
    first_adjustment = c(4, 6, 2, 5, 3, 9, 8, 7),
    second_adjustment = c(5, 3, 8, 6, 9, 7, 2, 4), check.names = FALSE)
  run <- function(input, variables) mira_test_info_base(input, id = "subject_id", outcomes = "score",
    time_vars = c("score_t0", "score_t1", "score_t2"), covariates = variables,
    analyses = "model", verbose = FALSE)
  original <- run(data, c("first_adjustment", "second_adjustment"))
  names(data)[names(data) == "first_adjustment"] <- "age`years"
  names(data)[names(data) == "second_adjustment"] <- "ageyears"
  renamed <- run(data, c("age`years", "ageyears"))
  expect_identical(renamed$model$covariates_used, c("age`years", "ageyears"))
  fitted <- renamed$model$fitted_model
  expect_s4_class(fitted, "merMod")
  frame <- stats::model.frame(fitted)
  expect_true(all(c("age`years", "ageyears") %in% names(frame)))
  expect_setequal(all.vars(stats::formula(fitted)),
    c("value", "time_factor", "age`years", "ageyears", "patient_factor"))
  expect_equal(unname(lme4::fixef(fitted)),
    unname(lme4::fixef(original$model$fitted_model)), tolerance = 1e-10)
  expected_design <- unname(stats::model.matrix(original$model$fitted_model))
  expect_equal(unname(stats::model.matrix(fitted)), expected_design, tolerance = 0)
  gee <- renamed$advanced_models$gee$independence
  expect_true(gee$performed)
  expect_setequal(all.vars(gee$formula), c("value", "time_factor", "age`years", "ageyears"))
  expect_equal(unname(stats::coef(gee$model)),
    unname(stats::coef(original$advanced_models$gee$independence$model)), tolerance = 1e-10)
  expect_equal(unname(stats::model.matrix(gee$model)),
    unname(stats::model.matrix(original$advanced_models$gee$independence$model)), tolerance = 0)
  expect_equal(renamed$descriptives, original$descriptives)
  expect_equal(renamed$change, original$change)
})
