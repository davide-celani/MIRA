mira_schema_print_fixture_20260918 <- function() {
  data.frame(
    participant_id = sprintf("S%02d", seq_len(6L)),
    score_t0 = c(1, 2, 3, 4, 5, 6),
    score_t1 = c(3, 1, 5, 4.5, 7, 8),
    score_t2 = c(4, 3, 6, 8, 6, 10),
    stringsAsFactors = FALSE
  )
}

mira_schema_print_single_20260918 <- function() {
  mira_info_long(
    mira_schema_print_fixture_20260918(),
    id = "participant_id",
    outcomes = "score",
    time_vars = c("score_t0", "score_t1", "score_t2"),
    time_labels = c("Baseline", "Middle", "Final"),
    covariates = character(0),
    analyses = "none",
    improvement_direction = "higher",
    verbose = FALSE
  )
}

mira_schema_print_multi_20260918 <- function() {
  data <- mira_schema_print_fixture_20260918()
  data$fatigue_t0 <- c(10, 9, 8, 7, 6, 5)
  data$fatigue_t1 <- c(9, 8, 8, 6, 5, 3)
  mira_info_long(
    data,
    id = "participant_id",
    outcomes = c("score", "fatigue"),
    time_vars = list(
      score = c("score_t0", "score_t1", "score_t2"),
      fatigue = c("fatigue_t0", "fatigue_t1")
    ),
    time_labels = list(
      score = c("Baseline", "Middle", "Final"),
      fatigue = c("Baseline", "Final")
    ),
    covariates = character(0),
    analyses = "none",
    improvement_direction = c(score = "higher", fatigue = "lower"),
    verbose = FALSE
  )
}

mira_schema_print_capture_20260918 <- function(expr) {
  visible <- NULL
  output <- capture.output(visible <- withVisible(force(expr)))
  list(output = output, value = visible$value, visible = visible$visible)
}

mira_schema_print_types_20260918 <- function(x) {
  vapply(x, typeof, character(1L))
}

test_that("single-outcome public schema preserves exact order, classes, and types", {
  x <- mira_schema_print_single_20260918()

  expect_identical(class(x), c("mira_info_long", "list"))
  expected_top_types <- c(
    call = "language", config = "list", settings = "list",
    outcome = "character", outcome_display = "character",
    time_vars = "character", time_labels = "character",
    data_overview = "list", detected_variables = "list", overview = "list",
    missing = "list", descriptives = "list", change = "list",
    trajectories = "list", arm_analysis = "list", correlations = "list",
    variability = "list", model = "list", advanced_tests = "list",
    advanced_models = "list", effect_sizes = "list", multiplicity = "list",
    robustness = "list", sensitivity = "list", outliers = "list",
    diagnostics = "list", plots = "list", plot_error = "NULL", long_data = "list"
  )
  expect_identical(names(x), names(expected_top_types))
  expect_identical(mira_schema_print_types_20260918(x), expected_top_types)

  expect_identical(names(x$settings), c(
    "alpha", "confidence_level", "confidence_percent", "p_adjust_method",
    "improvement_direction", "stable_threshold", "covariates", "arm_variable",
    "reference_arm", "arm_tests", "strict_id", "non_finite_handling"
  ))
  expect_identical(mira_schema_print_types_20260918(x$settings), c(
    alpha = "double", confidence_level = "double", confidence_percent = "double",
    p_adjust_method = "character", improvement_direction = "character",
    stable_threshold = "double", covariates = "character", arm_variable = "NULL",
    reference_arm = "NULL", arm_tests = "logical", strict_id = "logical",
    non_finite_handling = "character"
  ))
  expect_identical(names(x$overview), c(
    "outcome", "outcome_display", "id", "n_rows", "n_patients", "n_timepoints",
    "timepoints", "time_labels", "supplied_time_vars", "missing_ids",
    "duplicated_ids", "complete_profiles", "complete_profiles_pct",
    "non_finite_values", "covariates", "arm_variable", "reference_arm",
    "arm_levels", "arm_counts", "missing_arm", "arm_analysis"
  ))
  expect_identical(mira_schema_print_types_20260918(x$overview), c(
    outcome = "character", outcome_display = "character", id = "character",
    n_rows = "integer", n_patients = "integer", n_timepoints = "integer",
    timepoints = "character", time_labels = "character",
    supplied_time_vars = "logical", missing_ids = "integer",
    duplicated_ids = "integer", complete_profiles = "integer",
    complete_profiles_pct = "double", non_finite_values = "integer",
    covariates = "character", arm_variable = "NULL", reference_arm = "NULL",
    arm_levels = "character", arm_counts = "NULL", missing_arm = "integer",
    arm_analysis = "logical"
  ))

  descriptive_types <- c(
    time = "character", label = "character", n = "integer", missing = "integer",
    missing_pct = "double", non_finite = "integer", non_finite_pct = "double",
    unavailable = "integer", unavailable_pct = "double",
    score_mean = "double", score_sd = "double", score_variance = "double",
    score_se = "double", score_ci_lower = "double", score_ci_upper = "double",
    score_median = "double", score_q1 = "double", score_q3 = "double",
    score_iqr = "double", score_min = "double", score_max = "double",
    score_cv_percent = "double"
  )
  expect_identical(names(x$descriptives), names(descriptive_types))
  expect_identical(mira_schema_print_types_20260918(x$descriptives), descriptive_types)

  expect_identical(names(x$missing), c("by_time", "by_patient"))
  expect_identical(mira_schema_print_types_20260918(x$missing$by_time), c(
    time = "character", label = "character", missing_n = "integer",
    missing_pct = "double", non_finite_n = "integer",
    non_finite_pct = "double", unavailable_n = "integer",
    unavailable_pct = "double"
  ))
  expect_identical(mira_schema_print_types_20260918(x$missing$by_patient), c(
    patient = "character", missing_n = "double", non_finite_n = "double",
    unavailable_n = "double", unavailable_pct = "double"
  ))
  change_types <- c(
    from = "character", from_label = "character", to = "character",
    to_label = "character", n = "integer",
    score_mean_from = "double", score_mean_to = "double",
    score_mean_change = "double", score_sd_change = "double",
    score_se_change = "double", score_ci_lower = "double", score_ci_upper = "double",
    score_cohens_dz = "double", score_increased_n = "integer", score_increased_pct = "double",
    score_decreased_n = "integer", score_decreased_pct = "double", score_stable_n = "integer",
    score_stable_pct = "double", score_improved_n = "integer", score_improved_pct = "double",
    score_worsened_n = "integer", score_worsened_pct = "double",
    score_paired_t_p = "double", score_paired_t_p_adj = "double",
    score_wilcoxon_p = "double", score_wilcoxon_p_adj = "double"
  )
  expect_identical(names(x$change), names(change_types))
  expect_identical(mira_schema_print_types_20260918(x$change), change_types)

  expect_identical(names(x$arm_analysis), c(
    "enabled", "arm_variable", "reference_arm", "levels", "counts",
    "descriptives", "baseline_balance", "missingness", "time_omnibus",
    "time_pairwise", "change_descriptives", "change_omnibus", "change_pairwise"
  ))
  expect_identical(names(x$correlations), c("pearson", "spearman", "pairwise_n"))
  expect_identical(names(x$model), c(
    "fitted_model", "summary", "anova", "global_time_test", "global_arm_test",
    "arm_time_interaction_test", "singular", "converged", "warnings", "error",
    "covariates_requested", "covariates_used", "covariates_skipped",
    "fixed_parameters"
  ))
  expect_identical(names(x$advanced_tests), c("emmeans", "rm_anova", "friedman",
    "kenward_roger", "rank_longitudinal"))
  expect_identical(names(x$advanced_models),
                   c("random_slope", "nlme", "gee", "model_comparison", "mmrm"))
  expect_identical(names(x$robustness), "club_sandwich")
  expect_identical(names(x$effect_sizes), c(
    "paired_cohens_dz", "paired_rank_biserial", "rm_anova", "friedman",
    "mixed_models"
  ))
  expect_identical(names(x$multiplicity), c(
    "primary_method", "note", "time_pairwise", "baseline_vs_followup",
    "consecutive_time", "ordinal_trends", "arm_pairwise", "simple_effects",
    "interaction_contrasts", "friedman_posthoc"
  ))
  expect_identical(names(x$outliers), c("note", "by_time", "change"))
  expect_identical(names(x$config), c(
    "id", "id_generated", "outcomes", "outcome_names", "time_vars", "time_labels", "arm",
    "reference_arm", "covariates", "covariates_detected", "covariates_numeric",
    "covariates_categorical", "improvement_direction", "stable_threshold",
    "variable_pattern", "auto_detected", "specified_manually", "sources",
    "alternatives", "analyses"
  ))
  expect_identical(names(x$data_overview), c("n_rows", "n_columns", "column_profile"))
  expect_identical(names(x$detected_variables), c(
    "id_candidates", "longitudinal", "longitudinal_map",
    "ambiguous_longitudinal", "outcomes", "arm_candidates",
    "covariates_numeric", "covariates_categorical"
  ))
  expect_identical(names(x$diagnostics), c("warnings", "adaptation"))
  expect_false("outcomes" %in% names(x))
  expect_identical(levels(x$long_data$time_label), c("Baseline", "Middle", "Final"))
  expect_identical(typeof(x$long_data$time_label), "integer")
})

test_that("zero-row change preserves its distinct typed schema", {
  data <- data.frame(
    participant_id = sprintf("S%02d", seq_len(4L)),
    score_t0 = c(1, 2, NA, NA),
    score_t1 = c(NA, NA, 3, 4),
    stringsAsFactors = FALSE
  )
  empty <- mira_info_long(
    data, id = "participant_id", outcomes = "score",
    time_vars = c("score_t0", "score_t1"), covariates = character(0),
    analyses = "none", verbose = FALSE
  )
  nonempty <- mira_schema_print_single_20260918()

  expect_identical(nrow(empty$change), 0L)
  expect_identical(names(empty$change), names(nonempty$change))
  expect_identical(mira_schema_print_types_20260918(empty$change),
                   mira_schema_print_types_20260918(nonempty$change))
  expect_identical(tail(names(empty$change), 4L),
    c("score_paired_t_p", "score_paired_t_p_adj", "score_wilcoxon_p", "score_wilcoxon_p_adj"))
  expect_false(any(c("mean_change", "score_change") %in% names(empty$change)))
  expect_false(any(c("mean_change", "score_change") %in% names(nonempty$change)))
})

test_that("inspect-only, multi-outcome, and failed-outcome schemas are exact", {
  inspect <- mira_info_long(
    mira_schema_print_fixture_20260918(),
    id = "participant_id", outcomes = "score",
    time_vars = c("score_t0", "score_t1", "score_t2"),
    covariates = character(0), analyses = "none", inspect_only = TRUE,
    verbose = FALSE
  )
  expect_identical(class(inspect), c("mira_detect_long", "list"))
  expect_identical(names(inspect), c(
    "call", "config", "data_overview", "detected_variables",
    "diagnostics"
  ))
  expect_identical(mira_schema_print_types_20260918(inspect), c(
    call = "language", config = "list",
    data_overview = "list", detected_variables = "list", diagnostics = "list"
  ))
  expect_identical(names(inspect$diagnostics), "warnings")
  expect_false(any(c("settings", "descriptives", "change", "model") %in%
                     names(inspect)))

  multi <- mira_schema_print_multi_20260918()
  expect_identical(class(multi), c("mira_info_multi_long", "mira_info_long", "list"))
  expect_identical(names(multi), c(
    "call", "config", "data_overview", "detected_variables",
    "diagnostics", "outcomes"
  ))
  expect_identical(mira_schema_print_types_20260918(multi), c(
    call = "language", config = "list",
    data_overview = "list", detected_variables = "list", diagnostics = "list",
    outcomes = "list"
  ))
  expect_identical(names(multi$outcomes), c("score", "fatigue"))
  expect_true(all(vapply(
    multi$outcomes,
    function(item) identical(class(item), c("mira_info_long", "list")),
    logical(1L)
  )))
  expect_identical(names(multi$diagnostics),
                   c("warnings", "adaptation", "failed_outcomes"))
  expect_identical(multi$diagnostics$failed_outcomes, character(0))

  failing <- data.frame(
    participant_id = sprintf("S%02d", seq_len(6L)),
    good_t0 = seq_len(6L),
    good_t1 = seq_len(6L) + 1L,
    stringsAsFactors = FALSE
  )
  failing$bad_t0 <- I(matrix(seq_len(12L), nrow = 6L))
  failing$bad_t1 <- seq_len(6L) + 2L
  expect_warning(
    failed_multi <- mira_info_long(
      failing, id = "participant_id", covariates = character(0),
      analyses = "none", verbose = FALSE
    ),
    "Analyses were not completed for: bad"
  )
  failure <- failed_multi$outcomes$bad
  expect_identical(class(failure), c("mira_info_error_long", "list"))
  expect_identical(names(failure), c("outcome", "error"))
  expect_identical(mira_schema_print_types_20260918(failure),
                   c(outcome = "character", error = "character"))
  expect_identical(failure$outcome, "bad")
  expect_true(length(failure$error) == 1L && nzchar(failure$error))
  expect_identical(failed_multi$diagnostics$failed_outcomes, "bad")

  failed_summary <- summary(failed_multi)
  expect_identical(class(failed_summary), c("summary.mira_info_multi_long", "list"))
  expect_identical(names(failed_summary),
                   c("config", "outcomes", "failed_outcomes"))
  expect_identical(failed_summary$outcomes$bad, failure)
  expect_identical(failed_summary$failed_outcomes, "bad")
})

test_that("core public objects survive RDS serialization exactly", {
  bundle <- list(
    single = mira_schema_print_single_20260918(),
    inspect = mira_info_long(
      mira_schema_print_fixture_20260918(),
      id = "participant_id", outcomes = "score",
      time_vars = c("score_t0", "score_t1", "score_t2"),
      covariates = character(0), analyses = "none", inspect_only = TRUE,
      verbose = FALSE
    ),
    multi = mira_schema_print_multi_20260918()
  )
  path <- tempfile(fileext = ".rds")
  withr::defer(unlink(path), envir = parent.frame())
  saveRDS(bundle, path)
  restored <- readRDS(path)

  expect_identical(restored, bundle)
  expect_identical(class(restored$single), c("mira_info_long", "list"))
  expect_identical(class(restored$inspect), c("mira_detect_long", "list"))
  expect_identical(class(restored$multi), c("mira_info_multi_long", "mira_info_long", "list"))
})

test_that("detect and single-outcome print methods are semantic and invisible", {
  inspect <- mira_info_long(
    mira_schema_print_fixture_20260918(),
    id = "participant_id", outcomes = "score",
    time_vars = c("score_t0", "score_t1", "score_t2"),
    covariates = character(0), analyses = "none", inspect_only = TRUE,
    verbose = FALSE
  )
  detect_print <- mira_schema_print_capture_20260918(print(inspect))
  detect_text <- paste(detect_print$output, collapse = "\n")
  expect_false(detect_print$visible)
  expect_identical(detect_print$value, inspect)
  expect_true(grepl("DETECTION OVERVIEW", detect_text, fixed = TRUE))
  expect_true(grepl("ID: participant_id", detect_text, fixed = TRUE))
  expect_true(grepl("Outcomes: score", detect_text, fixed = TRUE))
  expect_true(grepl("AVAILABLE OUTPUT ELEMENTS", detect_text, fixed = TRUE))
  expect_true(grepl("$config", detect_text, fixed = TRUE))

  x <- mira_schema_print_single_20260918()
  info_print <- mira_schema_print_capture_20260918(print(
    x, digits = 2, max_rows = Inf, correlations = FALSE, model = FALSE,
    outliers = FALSE, trajectories = FALSE, missingness = FALSE, plots = FALSE
  ))
  info_text <- paste(info_print$output, collapse = "\n")
  expect_false(info_print$visible)
  expect_identical(info_print$value, x)
  expect_true(grepl("KEY RESULTS", info_text, fixed = TRUE))
  expect_true(grepl("Outcome: SCORE", info_text, fixed = TRUE))
  expect_true(grepl("DESCRIPTIVE STATISTICS", info_text, fixed = TRUE))
  expect_true(grepl("LONGITUDINAL CHANGE", info_text, fixed = TRUE))
  expect_true(grepl("WITHIN- / BETWEEN-SUBJECT VARIABILITY", info_text, fixed = TRUE))
  expect_true(grepl("COMPLETE OUTPUT GUIDE", info_text, fixed = TRUE))
  expect_false(grepl("\nCORRELATIONS\n", info_text, fixed = TRUE))
  expect_false(grepl("\nMIXED-EFFECTS MODEL\n", info_text, fixed = TRUE))
})

test_that("print row limits and numeric formatting preserve current validation", {
  x <- mira_schema_print_single_20260918()
  print_args <- list(
    correlations = FALSE, model = FALSE, outliers = FALSE,
    trajectories = FALSE, missingness = FALSE, plots = FALSE
  )

  truncated <- mira_schema_print_capture_20260918(do.call(
    print, c(list(x = x, digits = 2, max_rows = 1), print_args)
  ))
  expect_true(any(grepl(
    "Showing the first 1 of 3 pairwise comparisons", truncated$output, fixed = TRUE
  )))
  expect_false(truncated$visible)
  expect_identical(truncated$value, x)

  complete <- mira_schema_print_capture_20260918(do.call(
    print, c(list(x = x, digits = 2, max_rows = Inf), print_args)
  ))
  expect_false(any(grepl("Showing the first", complete$output, fixed = TRUE)))
  expect_false(complete$visible)
  expect_identical(complete$value, x)

  for (value in list(-1, NA_real_, "2", c(1, 2), NULL)) {
    expect_error(print(x, digits = value), "digits must be an integer >= 0")
  }
  for (value in list(0, -1, NA_real_, "2", c(1, 2), NULL)) {
    expect_error(print(x, max_rows = value), "max_rows must be > 0 or Inf")
  }

  # Fractional digits are currently accepted and truncated by as.integer().
  integer_digits <- mira_schema_print_capture_20260918(do.call(
    print, c(list(x = x, digits = 1, max_rows = Inf), print_args)
  ))
  fractional_digits <- mira_schema_print_capture_20260918(do.call(
    print, c(list(x = x, digits = 1.9, max_rows = Inf), print_args)
  ))
  expect_identical(fractional_digits$output, integer_digits$output)
})

test_that("multi, error, and summary print methods return invisibly", {
  multi <- mira_schema_print_multi_20260918()
  multi_print <- mira_schema_print_capture_20260918(print(multi))
  multi_text <- paste(multi_print$output, collapse = "\n")
  expect_false(multi_print$visible)
  expect_identical(multi_print$value, multi)
  expect_true(grepl("MULTI-OUTCOME OVERVIEW", multi_text, fixed = TRUE))
  expect_true(grepl("Outcomes: 2", multi_text, fixed = TRUE))
  expect_true(grepl("score", multi_text, fixed = TRUE))
  expect_true(grepl("fatigue", multi_text, fixed = TRUE))
  expect_true(grepl("Access an outcome-specific result", multi_text, fixed = TRUE))

  failure <- structure(
    list(outcome = "bad", error = "deterministic analysis failure"),
    class = c("mira_info_error_long", "list")
  )
  error_print <- mira_schema_print_capture_20260918(print(failure))
  error_text <- paste(error_print$output, collapse = "\n")
  expect_false(error_print$visible)
  expect_identical(error_print$value, failure)
  expect_true(grepl("ANALYSIS ERROR", error_text, fixed = TRUE))
  expect_true(grepl(
    "Outcome bad: analysis not completed - deterministic analysis failure",
    error_text, fixed = TRUE
  ))

  x <- mira_schema_print_single_20260918()
  summary_x <- summary(x)
  expect_identical(class(summary_x), c("summary.mira_info_long", "list"))
  expect_identical(names(summary_x), c(
    "outcome", "n_patients", "n_timepoints", "complete_profiles",
    "complete_profiles_pct", "descriptives", "baseline_final", "variability",
    "arm_analysis", "global_time_test", "global_arm_test",
    "arm_time_interaction_test", "model_singular", "model_converged", "advanced"
  ))
  summary_print <- mira_schema_print_capture_20260918(print(summary_x, digits = 2))
  summary_text <- paste(summary_print$output, collapse = "\n")
  expect_false(summary_print$visible)
  expect_identical(summary_print$value, summary_x)
  expect_true(grepl("SUMMARY", summary_text, fixed = TRUE))
  expect_true(grepl("Outcome: score", summary_text, fixed = TRUE))
  expect_true(grepl("Subjects: 6 | Timepoints: 3", summary_text, fixed = TRUE))
  expect_true(grepl("Baseline-final mean change", summary_text, fixed = TRUE))
  expect_true(grepl("ICC:", summary_text, fixed = TRUE))
  expect_true(grepl("Friedman: available", summary_text, fixed = TRUE))

  summary_multi <- summary(multi)
  expect_identical(class(summary_multi), c("summary.mira_info_multi_long", "list"))
  expect_identical(names(summary_multi),
                   c("config", "outcomes", "failed_outcomes"))
  summary_multi_print <- mira_schema_print_capture_20260918(
    print(summary_multi, digits = 2)
  )
  summary_multi_text <- paste(summary_multi_print$output, collapse = "\n")
  expect_false(summary_multi_print$visible)
  expect_identical(summary_multi_print$value, summary_multi)
  expect_true(grepl(
    "mira_info_long multi-outcome summary - 2 outcomes",
    summary_multi_text, fixed = TRUE
  ))
  expect_true(grepl("Outcome: score", summary_multi_text, fixed = TRUE))
  expect_true(grepl("Outcome: fatigue", summary_multi_text, fixed = TRUE))
})
