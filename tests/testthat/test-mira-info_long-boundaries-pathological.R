mira_boundary_fit <- function(data, analyses = "none") {
  mira_info_long(
    data,
    id = "subject_id",
    time_vars = c("score_t0", "score_t1"),
    outcomes = "score",
    covariates = character(0),
    analyses = analyses,
    verbose = FALSE
  )
}

mira_boundary_arm_fit <- function(data) {
  mira_info_long(
    data,
    id = "subject_id",
    time_vars = c("score_t0", "score_t1"),
    outcomes = "score",
    arm = "arm",
    covariates = character(0),
    analyses = "arm_tests",
    verbose = FALSE
  )
}

test_that("automatic arm detection includes its cardinality cutoff", {
  # With 44 rows the structural ceiling is max(10, floor(44 / 4)) = 11.
  make_data <- function(groups) {
    data.frame(
      arm = sprintf("G%02d", rep(seq_len(groups), length.out = 44L)),
      stringsAsFactors = FALSE
    )
  }

  below <- MIRA:::.mira_detect_arm_long(make_data(10L))
  at <- MIRA:::.mira_detect_arm_long(make_data(11L))
  above <- MIRA:::.mira_detect_arm_long(make_data(12L))

  expect_identical(below$selected, "arm")
  expect_identical(below$candidates$levels_n, 10L)
  expect_identical(below$warnings, character(0))

  expect_identical(at$selected, "arm")
  expect_identical(at$source, "auto")
  expect_identical(at$candidates$levels_n, 11L)
  expect_identical(at$warnings, character(0))

  expect_null(above$selected)
  expect_identical(above$source, "none")
  expect_identical(nrow(above$candidates), 0L)
  expect_identical(
    above$warnings,
    paste0(
      "Treatment-arm candidates were excluded because they lacked at least ",
      "two usable groups or had implausible cardinality: arm."
    )
  )
})

test_that("automatic covariates include the degrees-of-freedom limit", {
  # Thirty rows permit exactly three automatically selected adjustment df.
  make_data <- function(number) {
    values <- lapply(seq_len(number), function(index) {
      seq_len(30L) + index * 100L
    })
    as.data.frame(setNames(values, paste0("x", seq_len(number))))
  }

  below <- MIRA:::.mira_detect_covariates_long(make_data(2L))
  at <- MIRA:::.mira_detect_covariates_long(make_data(3L))
  above <- MIRA:::.mira_detect_covariates_long(make_data(4L))

  expect_identical(below$selected, c("x1", "x2"))
  expect_identical(below$numeric, c("x1", "x2"))
  expect_identical(below$warnings, character(0))

  expect_identical(at$selected, c("x1", "x2", "x3"))
  expect_identical(at$numeric, c("x1", "x2", "x3"))
  expect_identical(at$warnings, character(0))

  expect_identical(above$selected, character(0))
  expect_identical(above$numeric, character(0))
  expect_identical(above$detected, c("x1", "x2", "x3", "x4"))
  expect_identical(
    above$warnings,
    paste0(
      "Detected 4 candidate covariates (x1, x2, x3, x4), but their ",
      "complexity exceeds the conservative data-driven limit (3 df): they ",
      "will not be included automatically in the model. Specify ",
      "covariates= to choose."
    )
  )
})

test_that("correlations require two variable timepoints with two observations", {
  exact_data <- data.frame(
    subject_id = 1:3,
    score_t0 = c(1, 2, NA_real_),
    score_t1 = c(2, 4, NA_real_)
  )
  expect_no_warning(exact <- mira_boundary_fit(exact_data, "correlations"))

  expected_names <- c("score_t0", "score_t1")
  expected_correlation <- matrix(
    1,
    nrow = 2L,
    ncol = 2L,
    dimnames = list(expected_names, expected_names)
  )
  expected_n <- matrix(
    2L,
    nrow = 2L,
    ncol = 2L,
    dimnames = list(expected_names, expected_names)
  )
  adaptation <- exact$diagnostics$adaptation$score
  expect_identical(adaptation$available_by_time,
                   c(score_t0 = 2L, score_t1 = 2L))
  expect_identical(adaptation$variable_timepoints, 2L)
  expect_identical(adaptation$disabled, character(0))
  expect_identical(exact$correlations$pearson, expected_correlation)
  expect_equal(exact$correlations$spearman, expected_correlation,
               tolerance = 1e-15)
  expect_identical(exact$correlations$pairwise_n, expected_n)

  constant_data <- exact_data
  constant_data$score_t1 <- c(5, 5, NA_real_)
  expect_no_warning(constant <- mira_boundary_fit(constant_data, "correlations"))
  constant_adaptation <- constant$diagnostics$adaptation$score
  expect_identical(constant_adaptation$available_by_time,
                   c(score_t0 = 2L, score_t1 = 2L))
  expect_identical(constant_adaptation$variable_timepoints, 1L)
  expect_identical(
    constant_adaptation$disabled,
    "correlations: fewer than two estimable timepoints"
  )
  expect_identical(
    constant$correlations,
    list(pearson = NULL, spearman = NULL, pairwise_n = NULL)
  )

  insufficient_data <- exact_data
  insufficient_data$score_t1 <- c(2, NA_real_, NA_real_)
  expect_no_warning(insufficient <-
                      mira_boundary_fit(insufficient_data, "correlations"))
  insufficient_adaptation <- insufficient$diagnostics$adaptation$score
  expect_identical(insufficient_adaptation$available_by_time,
                   c(score_t0 = 2L, score_t1 = 1L))
  expect_identical(insufficient_adaptation$variable_timepoints, 1L)
  expect_identical(
    insufficient_adaptation$disabled,
    "correlations: fewer than two estimable timepoints"
  )
  expect_identical(
    insufficient$correlations,
    list(pearson = NULL, spearman = NULL, pairwise_n = NULL)
  )
})

test_that("IQR flags start at four values and use strict outer fences", {
  n3 <- mira_boundary_fit(data.frame(
    subject_id = 1:3,
    score_t0 = rep(0, 3L),
    score_t1 = c(0, 1, 2)
  ), "outliers")
  expect_identical(nrow(n3$outliers$by_time$score_t1), 0L)
  expect_identical(
    names(n3$outliers$by_time$score_t1),
    c("row", "patient", "value", "lower_bound", "upper_bound")
  )
  expect_identical(nrow(n3$outliers$change$score_t0_to_score_t1), 0L)
  expect_identical(
    names(n3$outliers$change$score_t0_to_score_t1),
    c("patient", "change", "lower_bound", "upper_bound")
  )

  at_fence <- c(0, 1, 2, 7)
  quartiles <- stats::quantile(at_fence, c(0.25, 0.75), names = FALSE)
  upper_fence <- quartiles[[2L]] + 1.5 * diff(quartiles)
  expect_equal(at_fence[[4L]], upper_fence)
  n4 <- mira_boundary_fit(data.frame(
    subject_id = 1:4,
    score_t0 = rep(0, 4L),
    score_t1 = at_fence
  ), "outliers")
  expect_identical(nrow(n4$outliers$by_time$score_t1), 0L)
  expect_identical(nrow(n4$outliers$change$score_t0_to_score_t1), 0L)

  just_outside <- c(0, 1, 2, 7 + 1e-6)
  outside_quartiles <- stats::quantile(
    just_outside, c(0.25, 0.75), names = FALSE
  )
  outside_lower <- outside_quartiles[[1L]] - 1.5 * diff(outside_quartiles)
  outside_upper <- outside_quartiles[[2L]] + 1.5 * diff(outside_quartiles)
  flagged <- mira_boundary_fit(data.frame(
    subject_id = 1:4,
    score_t0 = rep(0, 4L),
    score_t1 = just_outside
  ), "outliers")

  time_flag <- flagged$outliers$by_time$score_t1
  expect_identical(time_flag$row, 4L)
  expect_identical(time_flag$patient, 4L)
  expect_equal(time_flag$value, just_outside[[4L]])
  expect_equal(time_flag$lower_bound, outside_lower)
  expect_equal(time_flag$upper_bound, outside_upper)
  expect_gt(time_flag$value, time_flag$upper_bound)

  change_flag <- flagged$outliers$change$score_t0_to_score_t1
  expect_identical(change_flag$patient, 4L)
  expect_equal(change_flag$change, just_outside[[4L]])
  expect_equal(change_flag$lower_bound, outside_lower)
  expect_equal(change_flag$upper_bound, outside_upper)
})

test_that("core summaries handle two rows, all missing values, and constants", {
  expect_no_warning(two <- mira_boundary_fit(data.frame(
    subject_id = 1:2,
    score_t0 = c(1, 2),
    score_t1 = c(2, 3)
  )))
  expect_identical(class(two), c("mira_info_long", "list"))
  expect_identical(two$overview$n_patients, 2L)
  expect_identical(two$overview$complete_profiles, 2L)
  expect_identical(two$descriptives$n, c(2L, 2L))
  expect_equal(two$descriptives$mean, c(1.5, 2.5))
  expect_identical(two$change$n, 2L)
  expect_equal(two$change$mean_change, 1)
  expect_equal(two$change$sd_change, 0)
  expect_true(is.na(two$change$cohens_dz))
  expect_identical(two$trajectories$direction,
                   c("Increase", "Increase"))

  expect_no_warning(all_missing <- mira_boundary_fit(data.frame(
    subject_id = 1:3,
    score_t0 = rep(NA_real_, 3L),
    score_t1 = rep(NA_real_, 3L)
  )))
  expect_identical(class(all_missing), c("mira_info_long", "list"))
  expect_identical(all_missing$overview$complete_profiles, 0L)
  expect_identical(all_missing$descriptives$n, c(0L, 0L))
  expect_identical(all_missing$descriptives$missing, c(3L, 3L))
  expect_identical(all_missing$descriptives$unavailable, c(3L, 3L))
  expect_true(all(is.na(all_missing$descriptives$mean)))
  expect_true(all(is.na(all_missing$descriptives$sd)))
  expect_identical(nrow(all_missing$change), 0L)
  expect_identical(
    names(all_missing$change),
    c(
      "from", "to", "from_label", "to_label", "n", "mean_from",
      "mean_to", "mean_change", "sd_change", "se_change", "ci_lower",
      "ci_upper", "cohens_dz", "increased_n", "increased_pct",
      "decreased_n", "decreased_pct", "stable_n", "stable_pct",
      "improved_n", "improved_pct", "worsened_n", "worsened_pct",
      "paired_t_p", "paired_t_p_adj", "wilcoxon_p", "wilcoxon_p_adj"
    )
  )
  expect_true(all(is.na(all_missing$trajectories$absolute_change)))
  expect_true(is.na(all_missing$variability$grand_mean))

  expect_no_warning(constant <- mira_boundary_fit(data.frame(
    subject_id = 1:4,
    score_t0 = rep(5, 4L),
    score_t1 = rep(5, 4L)
  )))
  expect_equal(constant$descriptives$mean, c(5, 5))
  expect_equal(constant$descriptives$sd, c(0, 0))
  expect_equal(constant$descriptives$variance, c(0, 0))
  expect_equal(constant$change$mean_change, 0)
  expect_equal(constant$change$sd_change, 0)
  expect_true(is.na(constant$change$cohens_dz))
  expect_identical(constant$change$stable_n, 4L)
  expect_equal(constant$change$stable_pct, 100)
  expect_true(is.nan(constant$change$paired_t_p))
  expect_true(is.nan(constant$change$wilcoxon_p))
  expect_true(is.na(constant$change$paired_t_p_adj))
  expect_true(is.na(constant$change$wilcoxon_p_adj))
  expect_equal(constant$variability$grand_mean, 5)
  expect_equal(constant$variability$between_subject_sd, 0)
  expect_equal(constant$variability$within_subject_sd, 0)
  expect_true(is.na(constant$variability$ICC))
  expect_identical(constant$trajectories$direction, rep("Stable", 4L))
})

test_that("factor arms drop unused levels and isolate blank or single groups", {
  unused_data <- data.frame(
    subject_id = 1:6,
    arm = factor(
      rep(c("control", "active"), each = 3L),
      levels = c("control", "active", "unused")
    ),
    score_t0 = 1:6,
    score_t1 = 2:7
  )
  expect_no_warning(unused <- mira_boundary_arm_fit(unused_data))
  expect_true(unused$arm_analysis$enabled)
  expect_identical(unused$arm_analysis$levels, c("control", "active"))
  expect_identical(names(unused$arm_analysis$counts), c("control", "active"))
  expect_identical(as.integer(unused$arm_analysis$counts), c(3L, 3L))
  expect_false("unused" %in% unused$arm_analysis$levels)

  blank_data <- data.frame(
    subject_id = 1:6,
    arm = factor(
      c("control", "active", " ", "", "control", "active"),
      levels = c("control", "active", " ", "", "unused")
    ),
    score_t0 = 1:6,
    score_t1 = 2:7
  )
  expect_warning(
    blank <- mira_boundary_arm_fit(blank_data),
    "2 subjects have a missing treatment arm"
  )
  expect_identical(blank$overview$missing_arm, 2L)
  expect_true(blank$arm_analysis$enabled)
  expect_identical(blank$arm_analysis$levels, c("control", "active"))
  expect_identical(as.integer(blank$arm_analysis$counts), c(2L, 2L))

  single_data <- data.frame(
    subject_id = 1:4,
    arm = factor(rep("control", 4L), levels = c("control", "unused")),
    score_t0 = 1:4,
    score_t1 = 2:5
  )
  expect_warning(
    single <- mira_boundary_arm_fit(single_data),
    "contains fewer than two observed groups"
  )
  expect_identical(single$overview$arm_levels, "control")
  expect_identical(as.integer(single$overview$arm_counts), 4L)
  expect_false(single$overview$arm_analysis)
  expect_false(single$arm_analysis$enabled)
  expect_true(single$config$analyses$arm_tests)
  expect_false(single$settings$arm_tests)
  expect_identical(
    single$diagnostics$adaptation$score$disabled,
    "arm_tests: no usable arm"
  )
})

test_that("explicit configuration preserves non-syntactic column names", {
  data <- data.frame(
    "patient id" = paste0("P", 1:6),
    "treatment group" = rep(c("Control", "Active"), each = 3L),
    "Outcome baseline" = c(1, 2, 3, 4, 5, 6),
    "Outcome - week 4" = c(2, 3, 5, 5, 7, 9),
    "age (years)" = c(21, 30, 40, 25, 35, 45),
    "covariate+1" = rep(c("A", "B"), 3L),
    check.names = FALSE
  )
  variables <- c("Outcome baseline", "Outcome - week 4")
  labels <- stats::setNames(c("Baseline", "Week 4"), variables)

  expect_no_warning(result <- mira_info_long(
    data,
    id = "patient id",
    time_vars = list("Outcome score" = variables),
    time_labels = list("Outcome score" = unname(labels)),
    arm = "treatment group",
    reference_arm = "Control",
    covariates = c("age (years)", "covariate+1"),
    analyses = "none",
    verbose = FALSE
  ))

  expect_identical(class(result), c("mira_info_long", "list"))
  expect_identical(result$outcome, "Outcome score")
  expect_identical(result$config$id, "patient id")
  expect_identical(result$config$arm, "treatment group")
  expect_identical(result$config$reference_arm, "Control")
  expect_identical(result$config$covariates,
                   c("age (years)", "covariate+1"))
  expect_identical(result$config$covariates_numeric, "age (years)")
  expect_identical(result$config$covariates_categorical, "covariate+1")
  expect_identical(result$config$time_vars[["Outcome score"]], variables)
  expect_identical(result$config$time_labels[["Outcome score"]], labels)
  expect_identical(result$time_vars, variables)
  expect_identical(result$time_labels, labels)
  expect_equal(result$descriptives$mean,
               c(mean(data[[variables[[1L]]]]), mean(data[[variables[[2L]]]])))
  expect_identical(
    names(result$long_data),
    c(
      "patient", "outcome", "time", "time_index", "time_label", "arm",
      "value", "age (years)", "covariate+1", "patient_mean"
    )
  )
  expect_identical(nrow(result$long_data), 12L)
})
