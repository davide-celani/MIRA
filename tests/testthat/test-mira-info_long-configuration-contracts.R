# The resolved chronology defines baseline/final and all change signs. These
# model-free cases complement numerical-oracle and schema/printing tests.
config_contract_data <- function() {
  data.frame(subject_id = sprintf("S%02d", 1:6),
             score_t0 = c(11, 12, 14, 15, 19, 20),
             score_t1 = c(4, 6, 5, 9, 7, 10), stringsAsFactors = FALSE)
}

config_contract_info <- function(data = config_contract_data(), ...) {
  args <- list(data = data, id = "subject_id", covariates = character(0),
               analyses = "none", verbose = FALSE)
  dots <- list(...)
  args[names(dots)] <- dots
  do.call(mira_info_long, args)
}

config_contract_capture <- function(expr) {
  warnings <- character(0)
  output <- capture.output(value <- withCallingHandlers(
    withVisible(force(expr)), warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }))
  list(value = value$value, visible = value$visible,
       warnings = warnings, output = output)
}

test_that("direct detection and inspection agree without attempting analysis", {
  data <- config_contract_data()
  args <- list(data = data, id = "subject_id", covariates = character(0),
               time_vars = c("score_t1", "score_t0"),
               time_labels = c("Follow-up", "Baseline"),
               improvement_direction = "lower", stable_threshold = 2,
               verbose = FALSE)
  results <- testthat::with_mocked_bindings(
    list(direct = config_contract_capture(do.call(mira_detect_long, args)),
         inspect = config_contract_capture(do.call(
           mira_info_long, c(args, list(inspect_only = TRUE, analyses = "all"))))),
    .mira_analyse_single_outcome_long = function(...) {
      stop("Detection unexpectedly attempted statistical analysis")
    }, .package = "MIRA")
  direct <- results$direct$value
  inspect <- results$inspect$value
  expect_identical(class(direct), c("mira_detect_long", "list"))
  expect_identical(inspect$config[names(direct$config)], direct$config)
  expect_identical(inspect$data_overview, direct$data_overview)
  expect_identical(inspect$detected_variables, direct$detected_variables)
  expect_identical(inspect$diagnostics, direct$diagnostics)
  expect_null(direct$config$analyses)
  expect_true(all(unlist(inspect$config$analyses[-1L])))
  for (captured in results) {
    expect_false(captured$visible)
    expect_identical(captured$output, character(0))
    expect_identical(captured$warnings, character(0))
  }
  for (bad in list(NULL, NA, 1, "FALSE", c(TRUE, FALSE))) {
    expect_error(mira_detect_long(data, verbose = bad), "verbose must be TRUE or FALSE")
  }
})

test_that("analyses NULL honors flags and matches the equivalent optional subset", {
  by_flags <- config_contract_info(
    analyses = NULL, model = FALSE, plots = FALSE, arm_tests = FALSE,
    correlations = TRUE, outliers = TRUE)
  by_subset <- config_contract_info(analyses = c("correlations", "outliers"))
  expect_identical(by_flags$config$analyses, by_subset$config$analyses)
  for (component in c("overview", "descriptives", "change", "missing",
                       "trajectories", "variability", "correlations", "outliers")) {
    expect_equal(by_flags[[component]], by_subset[[component]], info = component)
  }
  expect_null(by_flags$model$fitted_model)
  expect_length(by_flags$plots, 0L)
  expect_false(by_flags$advanced_tests$rm_anova$performed)
  expect_match(by_flags$advanced_tests$rm_anova$reason_skipped, "model = FALSE")
  expect_true(by_flags$advanced_tests$friedman$performed)
  expect_identical(by_flags$correlations$pairwise_n,
                   matrix(6L, 2L, 2L, dimnames = rep(list(by_flags$time_vars), 2L)))
})

test_that("flat manual columns split families and match heterogeneous labels by source", {
  data <- config_contract_data()
  data$stress_t2 <- c(5, 8, 4, 10, 6, 11)
  data$stress_t0 <- c(7, 9, 6, 12, 8, 13)
  data$stress_t1 <- c(6, 8, 5, 11, 7, 12)
  columns <- c("stress_t2", "score_t1", "stress_t0", "score_t0", "stress_t1")
  labels <- c(score_t1 = "Score final", stress_t1 = "Stress middle",
              score_t0 = "Score baseline", stress_t2 = "Stress final",
              stress_t0 = "Stress baseline")
  result <- config_contract_info(data, time_vars = columns, time_labels = labels)
  expect_identical(names(result$outcomes), c("stress", "score"))
  expect_identical(result$config$time_vars$stress, c("stress_t0", "stress_t1", "stress_t2"))
  expect_identical(result$config$time_vars$score, c("score_t0", "score_t1"))
  for (name in names(result$outcomes)) {
    one <- result$outcomes[[name]]
    expect_identical(one$time_labels, labels[one$time_vars])
    expect_identical(one$descriptives$label, unname(labels[one$time_vars]))
    expect_identical(levels(one$long_data$time_label), unname(labels[one$time_vars]))
    expect_identical(one$config$time_vars, one$time_vars)
    expect_identical(one$config$outcomes, name)
  }
  expect_true(result$config$specified_manually$time_vars)
  expect_true(result$config$specified_manually$time_labels)
})

test_that("unnamed groups use requested outcome names with independent time labels", {
  data <- config_contract_data()
  data$first <- c(7, 9, 6, 12, 8, 13)
  data$last <- c(5, 8, 4, 10, 6, 11)
  groups <- list(c("score_t1", "score_t0"), c("first", "last"))
  labels <- list(`Pain score` = c("After", "Before"), `Stress score` = c("Start", "End"))
  result <- config_contract_info(data, time_vars = groups,
    outcomes = c("Pain score", "Stress score"), time_labels = labels)
  expect_identical(result$config$outcomes, c("Pain score", "Stress score"))
  expect_identical(unname(result$config$time_labels[["Pain score"]]), c("Before", "After"))
  expect_identical(unname(result$config$time_labels[["Stress score"]]), c("Start", "End"))
  expect_equal(result$outcomes[["Stress score"]]$trajectories$absolute_change,
               data$last - data$first)
  expect_error(config_contract_info(data, time_vars = groups, outcomes = "one"),
               "unnamed time_vars list requires one outcome per element")
  expect_error(config_contract_info(data, time_vars = list()), "time_vars cannot be an empty list")
  expect_error(config_contract_info(data, time_vars = groups,
    time_labels = list(c("A", "B"), c("C", "D"))),
    "For multiple outcomes, time_labels must be a named list")
})

test_that("manual labels repair duplicate temporal suffixes in the supplied order", {
  data <- config_contract_data()
  data$score_time1 <- data$score_t0
  data$score_t0 <- NULL
  result <- config_contract_info(data, time_vars = c("score_t1", "score_time1"),
    time_labels = c("Before procedure", "After procedure"))
  expect_identical(result$time_vars, c("score_t1", "score_time1"))
  expect_identical(unname(result$time_labels), c("Before procedure", "After procedure"))
  expect_equal(result$change$mean_change, mean(data$score_time1 - data$score_t1))
  expect_identical(result$diagnostics$warnings, character(0))
  expect_identical(result$config$alternatives$ambiguous_longitudinal, "score")
  expect_identical(result$detected_variables$ambiguous_longitudinal$score,
                   c("score_t1", "score_time1"))
})

test_that("custom chronology determines baseline and final throughout public results", {
  data <- config_contract_data()
  # Conventional suffixes deliberately oppose the custom study chronology.
  # Re-sorting after configuration reverses the scientific direction of change.
  parser <- function(variable) {
    if (!variable %in% c("score_t0", "score_t1")) return(NULL)
    first <- variable == "score_t1"
    list(outcome = "score", time_label = if (first) "Start" else "End",
         time_order = if (first) 0 else 1)
  }
  for (explicit in c(FALSE, TRUE)) {
    args <- list(variable_pattern = parser)
    if (explicit) args$time_vars <- c("score_t0", "score_t1")
    result <- do.call(config_contract_info, c(list(data = data), args))
    expected_vars <- c("score_t1", "score_t0")
    expected_labels <- c(score_t1 = "Start", score_t0 = "End")
    expect_identical(result$config$time_vars$score, expected_vars)
    expect_identical(result$time_vars, expected_vars, info = paste("explicit", explicit))
    expect_identical(result$time_labels, expected_labels)
    expect_identical(result$overview$timepoints, expected_vars)
    expect_identical(result$descriptives$time, expected_vars)
    expect_equal(result$descriptives$mean, c(mean(data$score_t1), mean(data$score_t0)))
    expect_identical(result$change$from, "score_t1")
    expect_identical(result$change$to, "score_t0")
    expect_equal(result$change$mean_change, mean(data$score_t0 - data$score_t1))
    expect_identical(levels(result$long_data$time_label), c("Start", "End"))
    expect_identical(result$long_data$time, rep(expected_vars, each = nrow(data)))
    expect_equal(result$trajectories$baseline, data$score_t1)
    expect_equal(result$trajectories$final, data$score_t0)
    expect_equal(result$trajectories$absolute_change, data$score_t0 - data$score_t1)
    expect_identical(summary(result)$baseline_final, result$change)
  }
})

test_that("unselected longitudinal groups and retained ID-arm alternatives cannot become covariates", {
  n <- 40L
  data <- data.frame(patient_id = sprintf("P%02d", seq_len(n)),
    record_id = sprintf("R%02d", seq_len(n)),
    arm = rep(c("control", "active"), each = 20L),
    treatment_group = rep(c("A", "B"), length.out = n),
    score_t0 = seq_len(n), score_t1 = seq_len(n) + 2L,
    other_t0 = seq_len(n) * 2L, other_t1 = seq_len(n) * 3L,
    orphan_t0 = seq_len(n) * 4L,
    conflict_t1 = seq_len(n), conflict_time1 = seq_len(n) + 1L,
    age = rep(30:39, 4L), site = rep(c("north", "south"), 20L))
  captured <- config_contract_capture(mira_info_long(data, outcomes = "score",
    analyses = "none", inspect_only = TRUE, verbose = FALSE))
  x <- captured$value
  expect_identical(x$config$id, "patient_id")
  expect_identical(x$config$arm, "arm")
  expect_setequal(x$config$covariates, c("age", "site"))
  expect_identical(x$config$alternatives$id, "record_id")
  expect_identical(x$config$alternatives$arm, "treatment_group")
  expect_true(any(grepl("ambiguous timepoints", captured$warnings)))
  expect_identical(sort(captured$warnings), sort(x$diagnostics$warnings))
  expect_false("orphan_t0" %in% unlist(x$config$time_vars, use.names = FALSE))
  for (excluded in c("record_id", "treatment_group", "other_t0", "orphan_t0", "conflict_t1")) {
    expect_error(mira_info_long(data, outcomes = "score", covariates = excluded,
      inspect_only = TRUE, verbose = FALSE),
      "already used as the ID, treatment arm, or outcome", info = excluded)
  }
})

test_that("a structurally invalid high-priority ID does not mask a usable lower-priority ID", {
  data <- config_contract_data()
  data$patient_id <- c("duplicate", "duplicate", NA, "P4", "P5", "P6")
  data$subject_id <- NULL
  data$case_id <- sprintf("C%02d", seq_len(nrow(data)))
  result <- mira_detect_long(data, covariates = character(0), verbose = FALSE)
  expect_identical(result$config$id, "case_id")
  expect_false(result$config$id_generated)
  candidates <- result$detected_variables$id_candidates
  invalid <- candidates[candidates$variable == "patient_id", , drop = FALSE]
  expect_identical(invalid$missing_n, 1L)
  expect_identical(invalid$duplicated_n, 1L)
  expect_false(invalid$structurally_valid)
  expect_identical(result$diagnostics$warnings, character(0))
})

test_that("strict_id FALSE audits IDs while retaining all core records", {
  data <- config_contract_data()
  data$subject_id <- c("A", "A", "B", NA, "C", "D")
  before <- data
  captured <- config_contract_capture(config_contract_info(data, strict_id = FALSE))
  result <- captured$value
  expect_identical(data, before)
  expect_true(any(grepl("missing", captured$warnings)))
  expect_true(any(grepl("duplicat", captured$warnings)))
  expect_identical(result$overview$n_rows, 6L)
  expect_identical(result$overview$n_patients, 4L)
  expect_identical(result$overview$missing_ids, 1L)
  expect_identical(result$overview$duplicated_ids, 1L)
  expect_identical(result$descriptives$n, c(6L, 6L))
  expect_identical(result$change$n, 6L)
  expect_identical(result$trajectories$patient, data$subject_id)
  expect_identical(result$long_data$patient, rep(data$subject_id, 2L))
  expect_equal(result$trajectories$absolute_change, data$score_t1 - data$score_t0)
  expect_true(any(grepl("missing", result$diagnostics$warnings)))
  expect_true(any(grepl("duplicates", result$diagnostics$warnings)))
})

test_that("equal semantic arm candidates require an explicit override", {
  data <- config_contract_data()
  data[["treatment group"]] <- rep(c("control", "active"), 3L)
  data[["treatment-group"]] <- rep(c("A", "B"), 3L)
  captured <- config_contract_capture(config_contract_info(data))
  expect_null(captured$value$config$arm)
  expect_setequal(captured$value$config$alternatives$arm,
                  c("treatment group", "treatment-group"))
  expect_true(any(grepl("Ambiguous treatment-arm variable", captured$warnings)))
  expect_false(captured$value$arm_analysis$enabled)
  manual <- config_contract_info(data, arm = "treatment-group", reference_arm = "B")
  expect_identical(manual$config$arm, "treatment-group")
  expect_identical(manual$config$sources$arm, "manual")
  expect_identical(manual$overview$arm_levels, c("B", "A"))
  expect_identical(as.character(manual$long_data$arm), rep(data[["treatment-group"]], 2L))
})

test_that("numeric arm codes retain reference order even with arm tests disabled", {
  data <- config_contract_data()
  data$arm <- c(2L, 1L, 2L, 1L, 2L, 1L)
  result <- config_contract_info(data, reference_arm = "1")
  expect_identical(result$config$arm, "arm")
  expect_identical(result$config$sources$arm, "auto")
  expect_identical(result$overview$arm_levels, c("1", "2"))
  expect_identical(levels(result$long_data$arm), c("1", "2"))
  expect_identical(as.character(result$trajectories$arm), as.character(data$arm))
  expect_false(result$arm_analysis$enabled)
  expect_identical(as.integer(result$overview$arm_counts), c(3L, 3L))
  # With 44 rows eleven categorical groups are plausible. The safeguard for
  # numerical treatment codes still limits numerical arms to ten values.
  codes <- rep(1:11, length.out = 44L)
  numeric <- MIRA:::.mira_detect_arm_long(data.frame(arm = codes))
  categorical <- MIRA:::.mira_detect_arm_long(data.frame(arm = factor(codes)))
  expect_null(numeric$selected)
  expect_match(numeric$warnings, "implausible cardinality")
  expect_identical(categorical$selected, "arm")
})

test_that("reference aliases and ties depend on observed groups rather than factor levels", {
  for (label in c("CTRL", "Placebo", "standard", "usual care", "usual_care",
                  "usual-care", "comparator", "sham", "untreated")) {
    data <- config_contract_data()
    data$arm <- factor(c("drug", "drug", "drug", label, "", NA),
                       levels = c("unused", label, "drug", ""))
    captured <- config_contract_capture(config_contract_info(data))
    expect_identical(captured$value$config$reference_arm, label, info = label)
    expect_identical(captured$value$overview$arm_levels, c(label, "drug"))
    expect_identical(captured$value$overview$missing_arm, 2L)
    expect_true(any(grepl("2 subjects have a missing treatment arm", captured$warnings)))
  }
  data <- config_contract_data()
  data$arm <- factor(c("B", "A", "B", "A", "B", "A"), levels = c("A", "B"))
  expect_identical(config_contract_info(data)$config$reference_arm, "B")
  data$arm <- c("Placebo", "Control", "Placebo", "Control", "drug", "drug")
  captured <- config_contract_capture(config_contract_info(data))
  expect_identical(captured$value$config$reference_arm, "Placebo")
  expect_true(any(grepl("Multiple plausible reference arms", captured$warnings)))
})

test_that("observed factor levels set covariate cost and five-valued numeric inputs are categorical", {
  n <- 40L
  data <- data.frame(subject_id = seq_len(n), score_t0 = seq_len(n),
    score_t1 = seq_len(n) + 1L, age = rep(30:39, 4L),
    site = factor(rep(c("A", "B", "C", "D"), 10L), levels = LETTERS[1:8]))
  at_limit <- mira_detect_long(data, id = "subject_id", verbose = FALSE)
  expect_identical(at_limit$config$covariates, c("age", "site"))
  expect_identical(at_limit$config$covariates_categorical, "site")
  expect_identical(at_limit$diagnostics$warnings, character(0))
  data$site[1L] <- "E"
  above <- config_contract_capture(mira_detect_long(data, id = "subject_id", verbose = FALSE))
  expect_identical(above$value$config$covariates, character(0))
  expect_identical(above$value$config$covariates_detected, c("age", "site"))
  expect_true(any(grepl("complexity exceeds", above$warnings)))
  data$five_values <- rep(1:5, length.out = n)
  data$six_values <- rep(1:6, length.out = n)
  data$all_missing <- rep(NA_real_, n)
  data$constant <- 5L
  selected <- mira_detect_long(data, id = "subject_id",
    covariates = c("five_values", "six_values", "all_missing", "constant"), verbose = FALSE)
  expect_identical(selected$config$covariates_categorical, "five_values")
  expect_identical(selected$config$covariates_numeric, "six_values")
  expect_identical(selected$config$covariates,
                   c("five_values", "six_values", "all_missing", "constant"))
})

test_that("per-outcome directions and thresholds preserve provenance and response labels", {
  data <- data.frame(subject_id = 1:6, BCVA_t0 = rep(10, 6L),
    BCVA_t1 = c(8, 9, 10, 11, 12, NA),
    CMT_t0 = rep(100, 6L), CMT_t1 = c(97, 98, 100, 102, 103, 100))
  result <- config_contract_info(data, outcomes = c("cmt", "BcVa"),
    improvement_direction = c(CMT = "auto", bcva = "higher"),
    stable_threshold = c(cmt = 2, BCVA = 1))
  expect_identical(names(result$outcomes), c("cmt", "BcVa"))
  expect_identical(result$config$improvement_direction, c(cmt = "lower", BcVa = "higher"))
  expect_identical(result$config$stable_threshold, c(cmt = 2, BcVa = 1))
  expect_identical(result$config$sources$improvement_direction,
                   c(cmt = "internal_outcome_map", BcVa = "user"))
  expect_identical(result$config$auto_detected$improvement_direction, c(cmt = TRUE, BcVa = FALSE))
  expect_identical(result$config$auto_detected$stable_threshold, c(cmt = FALSE, BcVa = FALSE))
  expect_identical(result$outcomes$cmt$trajectories$clinical_direction,
                   c("Improved", "Stable", "Stable", "Stable", "Worsened", "Stable"))
  expect_identical(result$outcomes$BcVa$trajectories$clinical_direction,
                   c("Worsened", "Stable", "Stable", "Stable", "Improved", NA_character_))
  for (outcome in names(result$outcomes)) {
    one <- result$outcomes[[outcome]]
    expect_identical(one$config$outcomes, outcome)
    expect_identical(one$settings$improvement_direction,
                     result$config$improvement_direction[[outcome]])
    expect_identical(one$settings$stable_threshold, result$config$stable_threshold[[outcome]])
  }
})

test_that("normalized duplicates and incomplete outcome settings produce useful errors", {
  data <- config_contract_data()
  for (directions in list(c(score = "higher", SCORE = "lower"), c(other = "higher"))) {
    expect_error(config_contract_info(data, improvement_direction = directions),
      "improvement_direction does not contain exactly one value for 'score'")
  }
  expect_error(config_contract_info(data, improvement_direction = c(score = "higher", "lower")),
               "improvement_direction must be named by outcome")
  for (thresholds in list(c(score = 1, SCORE = 2), c(other = 1))) {
    expect_error(config_contract_info(data, stable_threshold = thresholds),
                 "stable_threshold does not contain exactly one value for 'score'")
  }
  expect_error(config_contract_info(data, stable_threshold = c(score = 1, 2)),
               "stable_threshold must be named by outcome")
  for (bad in list(NA_character_, "", " ", character(0), c("subject_id", "score_t0"))) {
    expect_error(config_contract_info(data, id = bad), "id must name exactly one variable")
    expect_error(config_contract_info(data, arm = bad), "arm must name exactly one variable")
  }
})

test_that("multi-outcome summaries preserve order and select baseline-final rows by name", {
  data <- config_contract_data()
  data$score_t2 <- c(13, 15, 16, 18, 20, 24)
  data$stress_t0 <- c(8, 9, 10, 11, 12, 13)
  data$stress_t1 <- c(7, 8, 8, 10, 11, 11)
  result <- config_contract_info(data, outcomes = c("stress", "score"))
  summarized <- summary(result)
  expect_identical(names(summarized$outcomes), c("stress", "score"))
  for (outcome in names(result$outcomes)) {
    one <- result$outcomes[[outcome]]
    expected <- one$change[one$change$from == one$time_vars[[1L]] &
      one$change$to == tail(one$time_vars, 1L), , drop = FALSE]
    expect_identical(summarized$outcomes[[outcome]]$baseline_final, expected)
    expect_identical(summarized$outcomes[[outcome]], summary(one))
  }
  empty <- config_contract_info(data, outcomes = character(0))
  empty_summary <- summary(empty)
  expect_identical(empty_summary$outcomes, list())
  expect_length(empty_summary$failed_outcomes, 0L)
  printed <- config_contract_capture(print(empty_summary))
  expect_false(printed$visible)
  expect_identical(printed$value, empty_summary)
  expect_true(any(grepl("0 outcomes", printed$output)))
  expect_identical(printed$warnings, character(0))
})

test_that("multi-outcome plot routing accepts normalized names and excludes failures", {
  data <- config_contract_data()
  data$stress_t0 <- c(8, 9, 10, 11, 12, 13)
  data$stress_t1 <- c(7, 8, 8, 10, 11, 11)
  result <- config_contract_info(data, outcomes = c("stress", "score"))
  expect_error(plot(result), "Specify outcome=.*stress, score")
  expect_error(plot(result, outcome = "absent"), "Outcome 'absent' is not available")
  expect_error(plot(result, outcome = "score", which = "mean_ci"), "Plot 'mean_ci' is not available")
  expect_error(plot(result$outcomes$score, which = "nonexistent"), "boxplot.*mean_ci")
  routed <- testthat::with_mocked_bindings(
    withVisible(plot(result, outcome = "SCORE", which = "mean_ci", marker = "forwarded")),
    plot.mira_info_long = function(x, which, ...) {
      invisible(list(outcome = x$outcome, which = which, dots = list(...)))
    }, .package = "MIRA")
  expect_false(routed$visible)
  expect_identical(routed$value,
    list(outcome = "score", which = "mean_ci", dots = list(marker = "forwarded")))
  failure <- structure(list(outcome = "stress", error = "unavailable measurement"),
                       class = c("mira_info_error_long", "list"))
  result$outcomes$stress <- failure
  fallback <- testthat::with_mocked_bindings(plot(result, which = "mean_ci"),
    plot.mira_info_long = function(x, which, ...) invisible(x$outcome), .package = "MIRA")
  expect_identical(fallback, "score")
  expect_error(plot(result, outcome = "stress"), "Outcome 'stress' is not available")
  result$outcomes$score <- failure
  expect_error(plot(result), "No outcome has available plots")
})

test_that("custom parsing without temporal metadata preserves authoritative manual order", {
  data <- config_contract_data()
  result <- config_contract_info(data, time_vars = c("score_t1", "score_t0"),
    outcomes = "score", time_labels = c("Start", "End"),
    variable_pattern = function(variable) NULL)
  expect_identical(result$config$time_vars$score, c("score_t1", "score_t0"))
  expect_identical(result$time_vars, result$config$time_vars$score)
  expect_identical(unname(result$time_labels), c("Start", "End"))
  expect_equal(result$change$mean_change, mean(data$score_t0 - data$score_t1))
  expect_equal(result$trajectories$baseline, data$score_t1)
  expect_equal(result$trajectories$final, data$score_t0)
  expect_identical(result$diagnostics$warnings, character(0))
})

test_that("unrelated measurements and incomplete suffixes do not invent longitudinal outcomes", {
  data <- data.frame(subject_id = 1:6, age = 31:36, pressure = 101:106,
    score_t = 1:6, scoretime1 = 2:7, score_1x = 3:8,
    orphan_t0 = 4:9, text_t0 = letters[1:6], text_t1 = letters[2:7])
  expect_error(mira_info_long(data, inspect_only = TRUE, analyses = "none", verbose = FALSE),
               "No outcome with at least two numeric longitudinal columns was detected")
  expect_error(mira_detect_long(data, verbose = FALSE),
               "Specify time_vars or variable_pattern")
})
