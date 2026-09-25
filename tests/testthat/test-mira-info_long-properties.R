mira_info_property_data <- function() {
  data.frame(
    subject_id = sprintf("S%02d", seq_len(12L)),
    pain_t0 = c(12, 8, 15, 10, 18, 11, 14, 9, 16, 13, 17, 7),
    pain_t1 = c(13, 9, 14, NA, 19, 10, 16, 9, 17, 15, 16, 8),
    pain_t2 = c(15, 10, 16, 13, 21, 12, 17, 11, NA, 16, 19, 9),
    stringsAsFactors = FALSE
  )
}

mira_info_property_run <- function(data, ...) {
  args <- list(
    data = data,
    analyses = "none",
    verbose = FALSE,
    covariates = character(0),
    improvement_direction = "higher",
    stable_threshold = 0.5
  )
  dots <- list(...)
  args[names(dots)] <- dots
  do.call(mira_info_long, args)
}

mira_info_property_order <- function(x, by) {
  stopifnot(is.data.frame(x), all(by %in% names(x)))
  ordered <- x[do.call(order, x[by]), , drop = FALSE]
  rownames(ordered) <- NULL
  ordered
}

mira_info_property_select <- function(x, columns) {
  stopifnot(is.data.frame(x), all(columns %in% names(x)))
  selected <- x[, columns, drop = FALSE]
  rownames(selected) <- NULL
  selected
}

mira_info_property_scientific_view <- function(x) {
  descriptive_columns <- c(
    "time", "n", "missing", "missing_pct", "non_finite",
    "non_finite_pct", "unavailable", "unavailable_pct", "mean", "sd",
    "variance", "se", "ci_lower", "ci_upper", "median", "q1", "q3",
    "iqr", "min", "max", "cv_percent"
  )
  missing_time_columns <- c(
    "time", "missing_n", "missing_pct", "non_finite_n",
    "non_finite_pct", "unavailable_n", "unavailable_pct"
  )
  change_columns <- c(
    "from", "to", "n", "mean_from", "mean_to", "mean_change",
    "sd_change", "se_change", "ci_lower", "ci_upper", "cohens_dz",
    "increased_n", "increased_pct", "decreased_n", "decreased_pct",
    "stable_n", "stable_pct", "improved_n", "improved_pct", "worsened_n",
    "worsened_pct", "paired_t_p", "wilcoxon_p", "paired_t_p_adj",
    "wilcoxon_p_adj"
  )
  variability_columns <- c(
    "outcome", "grand_mean", "between_subject_sd", "within_subject_sd",
    "ICC_anova_complete_profiles", "ICC_model", "ms_between", "ms_within",
    "df_between", "df_within", "ICC"
  )
  trajectory_columns <- c(
    "patient", "outcome", "arm", "baseline", "final", "absolute_change",
    "relative_change_percent", "direction", "clinical_direction"
  )

  descriptives <- mira_info_property_select(x$descriptives, descriptive_columns)
  descriptives <- mira_info_property_order(descriptives, "time")
  missing_by_time <- mira_info_property_select(x$missing$by_time,
                                                missing_time_columns)
  missing_by_time <- mira_info_property_order(missing_by_time, "time")
  missing_by_patient <- mira_info_property_order(x$missing$by_patient, "patient")
  change <- mira_info_property_select(x$change, change_columns)
  change <- mira_info_property_order(change, c("from", "to"))
  variability <- mira_info_property_select(x$variability, variability_columns)
  trajectories <- mira_info_property_select(x$trajectories, trajectory_columns)
  trajectories <- mira_info_property_order(trajectories, "patient")

  list(
    overview = x$overview[c(
      "outcome", "n_rows", "n_patients", "n_timepoints",
      "complete_profiles", "complete_profiles_pct", "non_finite_values"
    )],
    descriptives = descriptives,
    missing_by_time = missing_by_time,
    missing_by_patient = missing_by_patient,
    change = change,
    variability = variability,
    trajectories = trajectories,
    correlations = x$correlations
  )
}

test_that("scientific results are invariant to input row order", {
  data <- mira_info_property_data()
  permutation <- c(7L, 2L, 11L, 4L, 9L, 1L, 12L, 5L, 3L, 10L, 6L, 8L)
  permuted <- data[permutation, , drop = FALSE]
  rownames(permuted) <- NULL

  original <- mira_info_property_run(data, analyses = "correlations")
  shuffled <- mira_info_property_run(permuted, analyses = "correlations")

  expect_equal(
    mira_info_property_scientific_view(shuffled),
    mira_info_property_scientific_view(original),
    tolerance = 1e-12
  )
})

test_that("automatic and explicit configuration agree on scientific results", {
  data <- mira_info_property_data()
  automatic <- mira_info_property_run(data)
  explicit <- mira_info_property_run(
    data,
    id = "subject_id",
    outcomes = "pain",
    time_vars = c("pain_t0", "pain_t1", "pain_t2"),
    time_labels = c("t0", "t1", "t2")
  )

  expect_false(automatic$config$specified_manually$id)
  expect_true(explicit$config$specified_manually$id)
  expect_equal(
    mira_info_property_scientific_view(automatic),
    mira_info_property_scientific_view(explicit)
  )
})

test_that("display labels do not change numerical results", {
  data <- mira_info_property_data()
  short_labels <- mira_info_property_run(
    data,
    id = "subject_id",
    outcomes = "pain",
    time_vars = c("pain_t0", "pain_t1", "pain_t2"),
    time_labels = c("T0", "T1", "T2")
  )
  display_labels <- mira_info_property_run(
    data,
    id = "subject_id",
    outcomes = "pain",
    time_vars = c("pain_t0", "pain_t1", "pain_t2"),
    time_labels = c("Baseline visit", "Middle visit", "Final visit")
  )

  expect_false(identical(short_labels$descriptives$label,
                         display_labels$descriptives$label))
  expect_equal(
    mira_info_property_scientific_view(short_labels),
    mira_info_property_scientific_view(display_labels)
  )
})

test_that("irrelevant columns do not affect explicitly selected analysis", {
  data <- mira_info_property_data()
  augmented <- data
  augmented$unused_metric <- seq(100, 111)
  augmented$audit_note <- paste0("not analysed ", seq_len(nrow(augmented)))

  arguments <- list(
    id = "subject_id",
    outcomes = "pain",
    time_vars = c("pain_t0", "pain_t1", "pain_t2"),
    time_labels = c("t0", "t1", "t2")
  )
  base_result <- do.call(mira_info_property_run, c(list(data = data), arguments))
  augmented_result <- do.call(
    mira_info_property_run,
    c(list(data = augmented), arguments)
  )

  expect_equal(
    mira_info_property_scientific_view(augmented_result),
    mira_info_property_scientific_view(base_result)
  )
})

test_that("each multi-outcome result agrees with its standalone analysis", {
  data <- mira_info_property_data()
  data$fatigue_t0 <- c(6, 8, 5, 9, 7, 10, 4, 8, 6, 9, 5, 7)
  data$fatigue_t1 <- c(5, 8, 4, 8, 6, 9, 5, 7, 5, 8, 4, 6)
  data$fatigue_t2 <- c(4, 7, 3, 7, 5, 8, 4, 6, 4, 7, 3, 5)

  multi <- mira_info_property_run(
    data,
    id = "subject_id",
    outcomes = c("pain", "fatigue")
  )
  pain_only <- mira_info_property_run(data, id = "subject_id", outcomes = "pain")
  fatigue_only <- mira_info_property_run(
    data,
    id = "subject_id",
    outcomes = "fatigue"
  )

  expect_equal(
    mira_info_property_scientific_view(multi$outcomes$pain),
    mira_info_property_scientific_view(pain_only)
  )
  expect_equal(
    mira_info_property_scientific_view(multi$outcomes$fatigue),
    mira_info_property_scientific_view(fatigue_only)
  )
})

test_that("analysis is deterministic and does not consume random numbers", {
  checked <- withr::with_seed(20260918, {
    seed_before <- .Random.seed
    first <- mira_info_property_run(mira_info_property_data())
    seed_after_first <- .Random.seed
    second <- mira_info_property_run(mira_info_property_data())
    seed_after_second <- .Random.seed
    list(
      seed_before = seed_before,
      seed_after_first = seed_after_first,
      seed_after_second = seed_after_second,
      first = mira_info_property_scientific_view(first),
      second = mira_info_property_scientific_view(second)
    )
  })

  expect_identical(checked$seed_after_first, checked$seed_before)
  expect_identical(checked$seed_after_second, checked$seed_before)
  expect_equal(checked$second, checked$first)
})

test_that("mira_info_long leaves its input data frame unchanged", {
  data <- mira_info_property_data()
  data$site <- factor(rep(c("north", "south"), 6L))
  before <- unserialize(serialize(data, NULL))

  invisible(mira_info_property_run(
    data,
    id = "subject_id",
    outcomes = "pain",
    time_vars = c("pain_t0", "pain_t1", "pain_t2")
  ))

  expect_identical(data, before)
})

test_that("moderate deterministic data complete with structural and numerical sanity", {
  n_subjects <- 250L
  time_index <- 0:5
  time_vars <- paste0("score_t", time_index)
  measurements <- stats::setNames(
    lapply(time_index, function(index) seq_len(n_subjects) / 10 + index),
    time_vars
  )
  data <- data.frame(
    patient_id = sprintf("P%03d", seq_len(n_subjects)),
    measurements,
    check.names = FALSE
  )

  result <- mira_info_long(
    data,
    id = "patient_id",
    outcomes = "score",
    time_vars = time_vars,
    covariates = character(0),
    analyses = "none",
    verbose = FALSE
  )

  expect_identical(class(result), c("mira_info_long", "list"))
  expect_identical(result$overview$n_patients, n_subjects)
  expect_identical(result$overview$n_timepoints, length(time_index))
  expect_identical(nrow(result$long_data), n_subjects * length(time_index))
  expect_identical(result$descriptives$n, rep(n_subjects, length(time_index)))
  expect_equal(
    result$descriptives$mean,
    vapply(measurements, mean, numeric(1L)),
    tolerance = 1e-12,
    ignore_attr = TRUE
  )
  expect_identical(
    nrow(result$change),
    as.integer(choose(length(time_index), 2L))
  )
  expect_identical(result$change$n, rep(n_subjects, nrow(result$change)))

  from_index <- as.numeric(sub("score_t", "", result$change$from, fixed = TRUE))
  to_index <- as.numeric(sub("score_t", "", result$change$to, fixed = TRUE))
  expect_equal(result$change$mean_change, to_index - from_index,
               tolerance = 1e-12)
})
