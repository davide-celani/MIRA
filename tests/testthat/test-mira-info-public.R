public_info_data <- function(n = 30L) {
  data.frame(
    subject_id = seq_len(n),
    pain_t0 = seq_len(n),
    pain_t1 = seq_len(n) + rep(c(-1, 0, 2), length.out = n),
    pain_t2 = seq_len(n) + rep(c(0, 1, 3), length.out = n),
    arm = rep(c("control", "active"), length.out = n),
    age = rep(30:39, length.out = n),
    sex = factor(rep(c("F", "M"), length.out = n)),
    stringsAsFactors = FALSE
  )
}

public_inspect <- function(data = public_info_data(), ...) {
  args <- list(data = data, inspect_only = TRUE, verbose = FALSE,
               analyses = "none")
  dots <- list(...)
  args[names(dots)] <- dots
  do.call(mira_info, args)
}

test_that("public scalar arguments reject missing, wrong-type, and nonscalar values", {
  data <- public_info_data()
  for (argument in c("arm_tests", "plots", "model", "outliers",
                     "correlations", "verbose", "strict_id", "inspect_only")) {
    for (value in list(NA, NULL, 1, 0, "TRUE", c(TRUE, FALSE), logical(0))) {
      args <- list(data = data, inspect_only = TRUE, verbose = FALSE,
                   analyses = "none")
      args[argument] <- list(value)
      expect_error(do.call(mira_info, args),
                   paste0(argument, " must be TRUE or FALSE"), info = argument)
    }
  }
  for (value in list(0, 1, -0.1, 1.1, NA_real_, NaN, Inf, -Inf,
                     c(0.05, 0.1), "0.05", NULL)) {
    expect_error(public_inspect(data, alpha = value),
                 "alpha must be a number strictly between 0 and 1")
  }
  for (value in list("invalid", NA_character_, c("holm", "BH"), 1, NULL)) {
    expect_error(public_inspect(data, p_adjust_method = value),
                 "p_adjust_method must be one of")
  }
  for (method in stats::p.adjust.methods) {
    expect_identical(public_inspect(data, p_adjust_method = method)$config$analyses$core,
                     c("descriptives", "missingness", "change", "variability", "trajectories"))
  }
  expect_error(mira_info(matrix(1:4, 2), inspect_only = TRUE, verbose = FALSE),
               "data must be a data.frame")
  expect_error(mira_info(data[FALSE, ], inspect_only = TRUE, verbose = FALSE),
               "no observations")
  names(data)[2:3] <- "duplicate"
  expect_error(mira_info(data, inspect_only = TRUE, verbose = FALSE),
               "duplicate column names")
})

test_that("analyses replaces module switches but never core summaries", {
  data <- public_info_data()
  modes <- list(
    none = character(0),
    one = "model",
    several = c("model", "correlations"),
    mixed_case = c("model", "plots"),
    all = c("plots", "model", "outliers", "correlations", "arm_tests")
  )
  inputs <- list(none = "none", one = "model",
                 several = c("model", "correlations"),
                 mixed_case = c("MODEL", "model", "Plots"), all = "all")
  for (name in names(inputs)) {
    actual <- public_inspect(data, analyses = inputs[[name]])$config$analyses
    for (module in c("plots", "model", "outliers", "correlations", "arm_tests")) {
      expect_identical(actual[[module]], module %in% modes[[name]], info = name)
    }
    expect_length(actual$core, 5L)
  }
  overridden <- public_inspect(data, analyses = "model", plots = TRUE,
                               model = FALSE, outliers = TRUE,
                               correlations = TRUE, arm_tests = TRUE)
  expect_true(overridden$config$analyses$model)
  expect_false(overridden$config$analyses$plots)
  expect_false(overridden$config$analyses$outliers)
  expect_false(overridden$config$analyses$correlations)
  expect_false(overridden$config$analyses$arm_tests)
  expect_true(public_inspect(data, analyses = c("all", "none"))$config$analyses$model)
  expect_false(public_inspect(data, analyses = c("none", "model"))$config$analyses$model)
  expect_false(public_inspect(data, analyses = character(0))$config$analyses$model)
  for (value in list("unknown", NA_character_, c("model", "unknown"), 1, TRUE)) {
    expect_error(public_inspect(data, analyses = value),
                 "analyses must be|Unrecognized optional analyses")
  }
})

test_that("inspection reports resolved configuration and exits before analysis", {
  data <- public_info_data()
  visible <- withVisible(public_inspect(data))
  expect_false(visible$visible)
  x <- visible$value
  expect_s3_class(x, "mira_detect")
  expect_identical(class(x), c("mira_detect", "list"))
  expect_named(x, c("call", "version", "config", "data_overview",
                    "detected_variables", "diagnostics"))
  expect_identical(x$version, "4.0.0")
  expect_identical(x$config$id, "subject_id")
  expect_false(x$config$id_generated)
  expect_identical(x$config$outcomes, "pain")
  expect_identical(x$config$time_vars$pain, c("pain_t0", "pain_t1", "pain_t2"))
  expect_identical(unname(x$config$time_labels$pain), c("t0", "t1", "t2"))
  expect_identical(x$config$arm, "arm")
  expect_identical(x$config$reference_arm, "control")
  expect_identical(x$data_overview$n_rows, nrow(data))
  expect_identical(x$data_overview$n_columns, ncol(data))
  expect_equal(nrow(x$data_overview$column_profile), ncol(data))
  expect_true(all(c("variable", "class", "numeric", "unique_n", "missing_n",
                    "missing_pct") %in% names(x$data_overview$column_profile)))
  expect_true(all(c("id_candidates", "longitudinal", "longitudinal_map",
                    "ambiguous_longitudinal", "outcomes", "arm_candidates",
                    "covariates_numeric", "covariates_categorical") %in%
                  names(x$detected_variables)))
  expect_false("descriptives" %in% names(x))
  expect_false("model" %in% names(x))
  expect_identical(x$config$sources$id, "auto")
  expect_true(x$config$auto_detected$id)
  expect_false(x$config$specified_manually$id)
})

test_that("manual and automatic ID rules are distinguished", {
  data <- public_info_data()
  data$record_code <- paste0("rec", seq_len(nrow(data)))
  x <- public_inspect(data, id = "record_code")
  expect_identical(x$config$id, "record_code")
  expect_identical(x$config$sources$id, "manual")
  expect_true(x$config$specified_manually$id)
  expect_error(public_inspect(data, id = "absent"), "does not exist")
  expect_error(public_inspect(data, id = 1), "id must name exactly one variable")

  duplicate <- data
  duplicate$record_code[2] <- duplicate$record_code[1]
  expect_error(public_inspect(duplicate, id = "record_code"), "duplicate values")
  expect_warning(mira_info(duplicate, id = "record_code", strict_id = FALSE,
                           inspect_only = TRUE, analyses = "none", verbose = FALSE,
                           covariates = character(0)),
                 "duplicates")
  missing <- data
  missing$record_code[2] <- NA_character_
  expect_error(public_inspect(missing, id = "record_code"), "missing values")
  expect_warning(mira_info(missing, id = "record_code", strict_id = FALSE,
                           inspect_only = TRUE, analyses = "none", verbose = FALSE,
                           covariates = character(0)),
                 "missing values")

  no_id <- data[c("pain_t0", "pain_t1", "pain_t2")]
  expect_warning(generated <- public_inspect(no_id), "No unique ID")
  expect_true(generated$config$id_generated)
  expect_identical(generated$config$sources$id, "generated_row_id")
  expect_identical(generated$config$id, ".mira_subject_id")
  no_id$.mira_subject_id <- rep(1L, nrow(no_id))
  expect_warning(generated_collision <- public_inspect(no_id), "No unique ID")
  expect_identical(generated_collision$config$id, ".mira_subject_id_")

  ambiguous <- data
  ambiguous$patient_id <- seq_len(nrow(data))
  expect_warning(detected <- public_inspect(ambiguous), "Ambiguous ID")
  expect_true(detected$config$id_generated)
  expect_setequal(detected$config$alternatives$id[1:2],
                  c("patient_id", "subject_id"))
  expect_true(any(grepl("Ambiguous ID", detected$diagnostics$warnings)))
})

test_that("time grouping, ordering, labels and outcome selection honor explicit input", {
  data <- public_info_data()
  x <- public_inspect(data, time_vars = c("pain_t2", "pain_t0", "pain_t1"),
                      time_labels = c("Two", "Zero", "One"))
  expect_identical(x$config$time_vars$pain, c("pain_t0", "pain_t1", "pain_t2"))
  expect_identical(unname(x$config$time_labels$pain), c("Zero", "One", "Two"))
  expect_true(x$config$specified_manually$time_vars)
  expect_true(x$config$specified_manually$time_labels)
  named <- public_inspect(data, time_vars = c("pain_t2", "pain_t0", "pain_t1"),
                          time_labels = c(pain_t0 = "A", pain_t1 = "B", pain_t2 = "C"))
  expect_identical(unname(named$config$time_labels$pain), c("A", "B", "C"))
  for (value in list("pain_t0", c("pain_t0", "absent"),
                     c("pain_t0", "pain_t0"), c("pain_t0", NA_character_),
                     c("pain_t0", "arm"))) {
    expect_error(public_inspect(data, time_vars = value),
                 "time_vars|Longitudinal variables")
  }
  for (value in list(c("only"), c("", "follow"), c("same", "same"),
                     c("start", NA_character_))) {
    expect_error(public_inspect(data, time_vars = c("pain_t0", "pain_t1"),
                                time_labels = value), "time_labels")
  }

  multi <- data
  multi$stress_t0 <- seq_len(nrow(data)) * 2
  multi$stress_t1 <- seq_len(nrow(data)) * 3
  selected <- public_inspect(multi, outcomes = "stress")
  expect_identical(selected$config$outcomes, "stress")
  expect_identical(public_inspect(multi, outcomes = c("stress", "pain"))$config$outcomes,
                   c("stress", "pain"))
  expect_error(public_inspect(multi, outcomes = "absent"), "Outcomes not detected")
  expect_error(public_inspect(multi, outcomes = c("pain", "PAIN")), "unique names")
  expect_error(public_inspect(multi, outcomes = NA_character_), "outcomes must")
  expect_length(public_inspect(multi, outcomes = character(0))$config$outcomes, 0L)
  by_column <- public_inspect(data, outcomes = c("pain_t0", "pain_t1"))
  expect_identical(by_column$config$time_vars$pain, c("pain_t0", "pain_t1"))
  by_list <- public_inspect(multi,
                            time_vars = list(stress = c("stress_t0", "stress_t1"),
                                             pain = c("pain_t0", "pain_t1")),
                            time_labels = list(stress = c("S0", "S1"),
                                               pain = c("P0", "P1")))
  expect_identical(by_list$config$outcomes, c("stress", "pain"))
  expect_identical(unname(by_list$config$time_labels$stress), c("S0", "S1"))
})

test_that("custom parsers and ambiguous names affect only valid outcome groups", {
  data <- data.frame(subject_id = 1:6, score_v1 = 1:6,
                     score_v0 = 2:7)
  regex <- public_inspect(data, variable_pattern = "^(.+)_v([0-9]+)$")
  expect_identical(regex$config$time_vars$score, c("score_v0", "score_v1"))
  parser <- function(name) {
    if (!grepl("^score_v[0-9]+$", name)) return(NULL)
    list(outcome = "score", time_label = sub("score_", "", name),
         time_order = as.numeric(sub("score_v", "", name)))
  }
  custom <- public_inspect(data, variable_pattern = parser)
  expect_identical(custom$config$time_vars$score, c("score_v0", "score_v1"))
  expect_identical(custom$config$variable_pattern, "<custom function>")
  ambiguous <- data.frame(subject_id = 1:6, pain_t1 = 1:6,
                          pain_time1 = 2:7, score_t0 = 1:6, score_t1 = 2:7)
  expect_warning(detected <- public_inspect(ambiguous), "ambiguous timepoints")
  expect_identical(detected$config$outcomes, "score")
  expect_identical(detected$config$alternatives$ambiguous_longitudinal, "pain")
  expect_error(public_inspect(ambiguous, outcomes = "pain"), "ambiguous timepoints")
  expect_error(public_inspect(data, variable_pattern = 1), "variable_pattern")
  expect_error(public_inspect(data, variable_pattern = function(x) list(time_label = x)),
               "must return NULL or a list")
})

test_that("arm, reference, and covariate selection retain provenance", {
  data <- public_info_data(30L)
  x <- public_inspect(data)
  expect_identical(x$config$arm, "arm")
  expect_identical(x$config$reference_arm, "control")
  expect_setequal(x$config$covariates, c("age", "sex"))
  expect_identical(x$config$covariates_numeric, "age")
  expect_identical(x$config$covariates_categorical, "sex")
  explicit <- public_inspect(data, covariates = c("sex", "age"),
                             reference_arm = "active")
  expect_identical(explicit$config$covariates, c("sex", "age"))
  expect_identical(explicit$config$reference_arm, "active")
  expect_true(explicit$config$specified_manually$covariates)
  expect_true(explicit$config$specified_manually$reference_arm)
  expect_error(public_inspect(data, arm = "missing_arm"), "does not exist")
  expect_error(public_inspect(data, reference_arm = "missing"), "not present")
  expect_error(public_inspect(data[c("subject_id", "pain_t0", "pain_t1")],
                              reference_arm = "control"), "requires an arm")
  for (value in list("absent", c("age", "age"), "subject_id", "arm",
                     "pain_t0", 1, NA_character_)) {
    expect_error(public_inspect(data, covariates = value),
                 "covariates|Covariates|already used")
  }
  expect_warning(small <- public_inspect(public_info_data(12L)), "complexity exceeds")
  expect_length(small$config$covariates, 0L)
  expect_true(any(grepl("complexity exceeds", small$diagnostics$warnings)))
  expect_identical(public_inspect(data, covariates = character(0))$config$covariates,
                   character(0))
})

test_that("single and multi-outcome public results have stable shape", {
  data <- public_info_data(8L)
  single <- mira_info(data, id = "subject_id", analyses = "none", verbose = FALSE,
                      covariates = character(0))
  expect_s3_class(single, "mira_info")
  expect_identical(single$version, "4.0.0")
  expect_named(single$outcomes, "pain")
  expect_true(all(c("descriptives", "missing", "change", "variability",
                    "trajectories", "long_data", "config", "diagnostics") %in% names(single)))
  expect_identical(single$overview$n_patients, 8L)
  expect_identical(single$overview$n_timepoints, 3L)
  expect_length(single$diagnostics$adaptation, 1L)

  data$stress_t0 <- seq_len(nrow(data)) * 2
  data$stress_t1 <- seq_len(nrow(data)) * 3
  multi <- mira_info(data, id = "subject_id", analyses = "none", verbose = FALSE,
                     covariates = character(0))
  expect_identical(class(multi), c("mira_info_multi", "mira_info", "list"))
  expect_setequal(names(multi$outcomes), c("pain", "stress"))
  expect_true(all(vapply(multi$outcomes, inherits, logical(1), "mira_info")))
  expect_identical(multi$diagnostics$failed_outcomes, character(0))
  expect_setequal(names(multi$diagnostics$adaptation), c("pain", "stress"))
  empty <- mira_info(data, id = "subject_id", outcomes = character(0),
                     analyses = "none", verbose = FALSE,
                     covariates = character(0))
  expect_s3_class(empty, "mira_info_multi")
  expect_length(empty$outcomes, 0L)
})

test_that("requested modules record per-outcome data adaptations", {
  empty <- data.frame(subject_id = 1:6, empty_t0 = rep(NA_real_, 6),
                      empty_t1 = rep(NA_real_, 6))
  x <- mira_info(empty, id = "subject_id", covariates = character(0),
                 analyses = "all", verbose = FALSE)
  adaptation <- x$diagnostics$adaptation$empty
  expect_identical(unname(adaptation$available_by_time), c(0L, 0L))
  expect_identical(adaptation$variable_timepoints, 0L)
  expect_setequal(adaptation$disabled, c(
    "model: insufficient data/variability",
    "correlations: fewer than two estimable timepoints",
    "plots: no finite outcome values",
    "arm_tests: no usable arm"
  ))
  expect_true(x$config$analyses$model)
  expect_null(x$model$fitted_model)
  expect_null(x$correlations$pearson)
  expect_length(x$plots, 0L)
  expect_length(x$outliers$by_time, 2L)
  expect_true(all(vapply(x$outliers$by_time, nrow, integer(1)) == 0L))

  mixed <- data.frame(subject_id = 1:6, good_t0 = 1:6, good_t1 = 2:7,
                      empty_t0 = rep(NA_real_, 6), empty_t1 = rep(NA_real_, 6))
  y <- mira_info(mixed, id = "subject_id", covariates = character(0),
                 analyses = "correlations", verbose = FALSE)
  expect_s3_class(y, "mira_info_multi")
  expect_length(y$diagnostics$adaptation$good$disabled, 0L)
  expect_match(y$diagnostics$adaptation$empty$disabled,
               "correlations: fewer than two estimable timepoints")
  expect_true(is.matrix(y$outcomes$good$correlations$pearson))
  expect_null(y$outcomes$empty$correlations$pearson)
})

test_that("one outcome failure is isolated in a multi-outcome call", {
  data <- data.frame(subject_id = 1:6, good_t0 = 1:6, good_t1 = 2:7)
  # A numeric matrix column is accepted by data.frame and initial numeric
  # detection, but cannot represent one longitudinal value per subject.
  data$bad_t0 <- I(matrix(1:12, nrow = 6))
  data$bad_t1 <- 2:7
  expect_error(mira_info(data, id = "subject_id", outcomes = "bad",
                         covariates = character(0), analyses = "none",
                         verbose = FALSE))
  expect_warning(
    multi <- mira_info(data, id = "subject_id", covariates = character(0),
                       analyses = "none", verbose = FALSE),
    "Analyses were not completed for: bad"
  )
  expect_s3_class(multi, "mira_info_multi")
  expect_s3_class(multi$outcomes$good, "mira_info")
  expect_s3_class(multi$outcomes$bad, "mira_info_error")
  expect_identical(multi$diagnostics$failed_outcomes, "bad")
  expect_identical(multi$outcomes$bad$outcome, "bad")
  expect_true(is.character(multi$outcomes$bad$error) &&
                nzchar(multi$outcomes$bad$error))
})
