name_provenance_data <- function(name = "FluX_mean_sd") {
  data <- data.frame(
    subject_id = paste0("P", seq_len(8)),
    first = c(3, 5, 4, 8, 7, 11, 9, 13),
    last = c(4, 4, 6, 7, 10, 12, 8, 15),
    check.names = FALSE
  )
  names(data)[2:3] <- paste0(name, c("_t0", "_t1"))
  data
}

name_provenance_info <- function(data, ...) {
  mira_info_long(data, id = "subject_id", covariates = character(0),
                 analyses = "none", improvement_direction = "unknown",
                 verbose = FALSE, ...)
}

test_that("metric naming follows public declaration rather than resolved labels", {
  for (name in c("FluX_mean_sd", "deltaChange", "q_INDEX", "Measure space")) {
    data <- name_provenance_data(name)
    variables <- names(data)[2:3]
    inferred <- list(
      automatic = name_provenance_info(data),
      auto_selector = name_provenance_info(data, outcomes = "auto"),
      time_columns = name_provenance_info(data, time_vars = rev(variables)),
      outcome_columns = name_provenance_info(data, outcomes = rev(variables)),
      unnamed_list = name_provenance_info(data, time_vars = list(rev(variables)))
    )
    for (result in inferred) {
      identity <- result$config$outcome_names[[name]]
      expect_identical(identity$key, name)
      expect_null(identity$declared)
      expect_identical(identity$source, "inferred")
      expect_identical(identity$label, toupper(name))
      expect_true(all(c("mean", "sd", "se") %in% names(result$descriptives)))
      expect_false(any(paste0(name, c("_mean", "_sd", "_se")) %in%
                         names(result$descriptives)))
      expect_identical(result$time_vars, variables)
      expect_null(result$outcomes)
    }
    declared <- list(
      outcomes = name_provenance_info(data, outcomes = name),
      time_vars_names = name_provenance_info(
        data, time_vars = stats::setNames(list(rev(variables)), name)),
      unnamed_labels = name_provenance_info(
        data, time_vars = list(rev(variables)), outcomes = name)
    )
    for (mode in names(declared)) {
      result <- declared[[mode]]
      identity <- result$config$outcome_names[[name]]
      expect_identical(identity$declared, name)
      expect_identical(identity$source,
                       if (mode == "time_vars_names") mode else "outcomes")
      expect_false(any(c("mean", "sd", "se") %in% names(result$descriptives)))
      expect_true(all(paste0(name, c("_mean", "_sd", "_se")) %in%
                        names(result$descriptives)))
      expect_equal(result$descriptives[[paste0(name, "_mean")]],
                   inferred$automatic$descriptives$mean)
      expect_equal(result$change[[paste0(name, "_mean_change")]],
                   inferred$automatic$change$mean_change)
      expect_identical(result$time_vars, variables)
    }
  }
})

test_that("custom detection rules infer names and preserve declared capitalization", {
  data <- name_provenance_data()
  names(data)[2:3] <- c("first", "last")
  parser <- function(variable) {
    if (!variable %in% c("first", "last")) return(NULL)
    list(outcome = "Bio_index", time_label = variable,
         time_order = match(variable, c("first", "last")))
  }
  inferred <- name_provenance_info(data, variable_pattern = parser)
  explicit <- name_provenance_info(data, variable_pattern = parser,
                                  outcomes = "BIO_index")
  expect_null(inferred$config$outcome_names$Bio_index$declared)
  expect_identical(explicit$config$outcome_names$BIO_index$declared, "BIO_index")
  expect_true("mean" %in% names(inferred$descriptives))
  expect_true("BIO_index_mean" %in% names(explicit$descriptives))
  expect_identical(explicit$time_vars, c("first", "last"))
  expect_equal(explicit$descriptives$BIO_index_mean, inferred$descriptives$mean)
})

test_that("partially named lists retain independent declaration provenance", {
  data <- name_provenance_data("First")
  second <- name_provenance_data("Second")
  data[names(second)[2:3]] <- second[2:3]
  variables <- list(c("First_t1", "First_t0"), c("Second_t1", "Second_t0"))
  names(variables) <- c("Mean_change_Mix", "")
  result <- name_provenance_info(data, time_vars = variables)
  expect_identical(names(result$outcomes), c("Mean_change_Mix", "Second"))
  expect_identical(names(result$config$outcome_names), names(result$outcomes))
  first <- result$outcomes[[1L]]
  second <- result$outcomes[[2L]]
  expect_identical(first$config$outcome_names$Mean_change_Mix$source, "time_vars_names")
  expect_identical(second$config$outcome_names$Second$source, "inferred")
  expect_identical(first$config$outcome_names$Mean_change_Mix$declared, "Mean_change_Mix")
  expect_null(second$config$outcome_names$Second$declared)
  expect_true("Mean_change_Mix_mean" %in% names(first$descriptives))
  expect_true("mean" %in% names(second$descriptives))
  expect_identical(first$trajectories$patient, data$subject_id)
  expect_identical(second$trajectories$patient, data$subject_id)
  expect_identical(first$time_vars, c("First_t0", "First_t1"))
  expect_identical(second$time_vars, c("Second_t0", "Second_t1"))
})

test_that("declaration collisions fail explicitly without overwriting groups", {
  data <- name_provenance_data()
  variables <- names(data)[2:3]
  for (labels in list(c("same", "same"), c("Mean sd", "mean_sd"))) {
    expect_error(name_provenance_info(
      data, time_vars = stats::setNames(list(variables, variables), labels)),
      "Declared outcome names must be unique")
  }
  expect_error(name_provenance_info(
    data, outcomes = c("FluX_mean_sd", "flux_mean_sd")), "unique names")
})

test_that("inspection and direct detection expose the same name provenance", {
  data <- name_provenance_data()
  variables <- names(data)[2:3]
  for (args in list(list(), list(outcomes = rev(variables)),
                   list(time_vars = list(Declared = variables)))) {
    common <- c(list(data = data, id = "subject_id", covariates = character(0),
                     verbose = FALSE), args)
    direct <- do.call(mira_detect_long, common)
    inspect <- do.call(mira_info_long,
                       c(common, list(inspect_only = TRUE, analyses = "none")))
    expect_identical(direct$config$outcome_names, inspect$config$outcome_names)
    expect_false("version" %in% names(direct))
    expect_false("version" %in% names(inspect))
    expect_identical(names(direct), names(inspect))
  }
})


test_that("public metric collisions preserve source data and explain the conflicting path", {
  data <- name_provenance_data()
  data$Declared_patient_mean <- seq_len(nrow(data)) + 20
  original <- data
  expect_error(mira_info_long(
    data, id = "subject_id", time_vars = names(data)[2:3], outcomes = "Declared",
    covariates = "Declared_patient_mean", analyses = "none", verbose = FALSE),
    "collide.*long_data")
  expect_identical(data, original)
})

test_that("source names containing dollar signs remain intact in metric path references", {
  data <- name_provenance_data("Measure$mean")
  result <- mira_info_long(
    data, id = "subject_id", outcomes = "Measure$mean", covariates = character(0),
    analyses = "outliers", verbose = FALSE)
  table <- result$outliers$by_time[["Measure$mean_t0"]]
  expect_true("Measure$mean_lower_bound" %in% names(table))
  base <- MIRA:::.mira_metric_base_view_long(result)
  expect_true("lower_bound" %in% names(base$outliers$by_time[["Measure$mean_t0"]]))
  expect_identical(result$time_vars, names(data)[2:3])
  expect_identical(result$config$outcome_names[["Measure$mean"]]$declared, "Measure$mean")
  printed <- capture.output(returned <- withVisible(print(result)))
  expect_false(returned$visible)
  expect_true(length(printed) > 0L)
})
