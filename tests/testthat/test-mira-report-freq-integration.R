report_integration_data <- function() {
  data.frame(
    patient = seq_len(6L),
    score_t0 = c(1, 3, 2, 6, 5, 4),
    score_t1 = c(2, 4, 3, 5, 7, 6),
    score_t2 = c(4, 5, 6, 7, 8, 9)
  )
}

report_integration_multi_data <- function() {
  data.frame(
    patient = seq_len(6L),
    score_t0 = c(1, 3, 2, 6, 5, 4),
    score_t1 = c(2, 4, 3, 5, 7, 6),
    fatigue_t0 = c(9, 8, 7, 6, 5, 4),
    fatigue_t1 = c(8, 8, 6, 5, 4, 3)
  )
}

report_integration_fit <- function(data = report_integration_data(), ...) {
  args <- list(
    data = data,
    id = "patient",
    outcomes = "score",
    time_vars = c("score_t0", "score_t1", "score_t2"),
    covariates = character(0),
    analyses = "none",
    verbose = FALSE
  )
  dots <- list(...)
  args[names(dots)] <- dots
  do.call(mira_info, args)
}

report_integration_with_dir <- function(code) {
  path <- tempfile("mira-report-integration-")
  dir.create(path, recursive = TRUE)
  on.exit(unlink(path, recursive = TRUE, force = TRUE), add = TRUE)
  code(path)
}

report_integration_manifest <- function(report) {
  utils::read.csv(
    report$results_manifest,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

report_integration_read_export <- function(report, path, outcome) {
  manifest <- report_integration_manifest(report)
  index <- which(manifest$path == path & manifest$outcome == outcome)
  testthat::expect_length(index, 1L)
  utils::read.csv(
    file.path(report$output_dir, manifest$file[[index]]),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

test_that("a real single-outcome result is preserved in payload and scientific CSVs", {
  data <- report_integration_data()
  fit <- report_integration_fit(data)

  report_integration_with_dir(function(output_dir) {
    report <- report_test_run(
      fit,
      output_dir = output_dir,
      output_file = "single-integration",
      format = "html",
      render = FALSE,
      open = FALSE
    )

    expect_s3_class(report, "mira_report_freq")
    expect_true(all(file.exists(c(
      report$source,
      report$payload,
      report$runtime,
      report$results_manifest,
      report$result_files
    ))))

    payload <- readRDS(report$payload)
    expect_s3_class(payload$result, "mira_info")
    expect_equal(payload$result$descriptives, fit$descriptives)
    expect_equal(payload$result$missing$by_time, fit$missing$by_time)
    expect_equal(payload$result$change, fit$change)
    expect_equal(
      payload$result$advanced_tests$friedman$tidy,
      fit$advanced_tests$friedman$tidy
    )
    expect_identical(payload$result$config$analyses, fit$config$analyses)

    descriptives <- report_integration_read_export(
      report, "outcome$descriptives", "score"
    )
    expect_equal(descriptives$n, fit$descriptives$n)
    expect_equal(descriptives$mean, fit$descriptives$mean, tolerance = 1e-12)
    expect_equal(descriptives$sd, fit$descriptives$sd, tolerance = 1e-12)

    missing <- report_integration_read_export(
      report, "outcome$missing$by_time", "score"
    )
    expect_equal(missing$unavailable_n, fit$missing$by_time$unavailable_n)
    expect_equal(missing$unavailable_pct,
                 fit$missing$by_time$unavailable_pct, tolerance = 1e-12)

    change <- report_integration_read_export(report, "outcome$change", "score")
    expect_true(is.numeric(change$paired_t_p))
    expect_equal(change$n, fit$change$n)
    expect_equal(change$mean_change, fit$change$mean_change, tolerance = 1e-12)
    expect_equal(change$paired_t_p, fit$change$paired_t_p, tolerance = 1e-12)
    expect_false(any(grepl("<", as.character(change$paired_t_p), fixed = TRUE)))
  })
})

test_that("data and data-frame x inputs forward mira_info arguments and verbosity", {
  data <- report_integration_data()

  report_integration_with_dir(function(output_dir) {
    report <- report_test_run(
      data = data,
      id = "patient",
      outcomes = "score",
      time_vars = c("score_t0", "score_t1", "score_t2"),
      covariates = character(0),
      analyses = "none",
      p_adjust_method = "BH",
      improvement_direction = "higher",
      stable_threshold = 1,
      output_dir = output_dir,
      output_file = "data-argument",
      format = "html",
      render = FALSE,
      open = FALSE
    )

    expect_s3_class(report$mira_result, "mira_info")
    expect_identical(report$mira_result$settings$p_adjust_method, "BH")
    expect_identical(report$mira_result$settings$improvement_direction, "higher")
    expect_identical(report$mira_result$settings$stable_threshold, 1)
    expect_false(report$mira_result$config$analyses$model)
    expect_false(any(grepl(
      "^MIRA INFO$",
      attr(report, "test_printed_output")
    )))
  })

  report_integration_with_dir(function(output_dir) {
    report <- report_test_run(
      x = data,
      id = "patient",
      outcomes = "score",
      time_vars = c("score_t0", "score_t1", "score_t2"),
      time_labels = c("Baseline", "Middle", "Final"),
      covariates = character(0),
      analyses = "none",
      verbose = TRUE,
      output_dir = output_dir,
      output_file = "x-data-frame",
      format = "html",
      render = FALSE,
      open = FALSE
    )

    expect_s3_class(report$mira_result, "mira_info")
    expect_identical(
      unname(report$mira_result$time_labels),
      c("Baseline", "Middle", "Final")
    )
    expect_true(report$mira_result$config$specified_manually$time_labels)
    expect_true(any(grepl(
      "^MIRA INFO$",
      attr(report, "test_printed_output")
    )))
  })
})

test_that("a real multi-outcome result exports each outcome to its own directory", {
  data <- report_integration_multi_data()
  fit <- mira_info(
    data,
    id = "patient",
    outcomes = c("score", "fatigue"),
    time_vars = list(
      score = c("score_t0", "score_t1"),
      fatigue = c("fatigue_t0", "fatigue_t1")
    ),
    covariates = character(0),
    analyses = "none",
    verbose = FALSE
  )
  expect_identical(class(fit), c("mira_info_multi", "mira_info", "list"))

  report_integration_with_dir(function(output_dir) {
    report <- report_test_run(
      fit,
      output_dir = output_dir,
      output_file = "multi-integration",
      format = "html",
      render = FALSE,
      open = FALSE
    )
    manifest <- report_integration_manifest(report)
    outcome_rows <- !is.na(manifest$outcome) & nzchar(manifest$outcome)

    expect_setequal(unique(manifest$outcome[outcome_rows]), names(fit$outcomes))
    for (outcome in names(fit$outcomes)) {
      rows <- manifest$outcome == outcome
      prefix <- paste0("results/", .mira_report_file_slug(outcome), "/")
      expect_true(any(rows))
      expect_true(all(startsWith(manifest$file[rows], prefix)))
      expect_true(dir.exists(file.path(
        report$output_dir, "results", .mira_report_file_slug(outcome)
      )))
    }
    expect_true(all(file.exists(report$result_files)))
    expect_s3_class(readRDS(report$payload)$result, "mira_info_multi")
  })
})

test_that("mira_detect and inspect-only inputs create context without outcomes", {
  data <- report_integration_data()
  inputs <- list(
    direct = mira_detect(
      data,
      id = "patient",
      outcomes = "score",
      time_vars = c("score_t0", "score_t1", "score_t2"),
      covariates = character(0),
      verbose = FALSE
    ),
    inspect_only = mira_info(
      data,
      id = "patient",
      outcomes = "score",
      time_vars = c("score_t0", "score_t1", "score_t2"),
      covariates = character(0),
      analyses = "none",
      inspect_only = TRUE,
      verbose = FALSE
    )
  )

  for (name in names(inputs)) {
    report_integration_with_dir(function(output_dir) {
      report <- report_test_run(
        inputs[[name]],
        output_dir = output_dir,
        output_file = paste0("detect-", name),
        format = "html",
        render = FALSE,
        open = FALSE
      )
      payload_result <- readRDS(report$payload)$result
      manifest <- report_integration_manifest(report)

      expect_s3_class(payload_result, "mira_detect")
      expect_length(.mira_report_extract_outcomes(payload_result), 0L)
      expect_false(any(grepl("^outcome\\$", manifest$path)))
      expect_true(any(grepl("^result\\$", manifest$path)))
      expect_false(any(!is.na(manifest$outcome) & nzchar(manifest$outcome)))
    })
  }
})

test_that("a nested mira_info_error stays in the payload and only surfaces as diagnostics", {
  data <- data.frame(
    patient = seq_len(6L),
    good_t0 = c(1, 3, 2, 6, 5, 4),
    good_t1 = c(2, 4, 3, 5, 7, 6)
  )
  data$bad_t0 <- I(matrix(seq_len(12L), nrow = 6L))
  data$bad_t1 <- 2:7

  expect_warning(
    fit <- mira_info(
      data,
      id = "patient",
      outcomes = c("good", "bad"),
      covariates = character(0),
      analyses = "none",
      verbose = FALSE
    ),
    "Analyses were not completed for: bad"
  )
  expect_s3_class(fit$outcomes$bad, "mira_info_error")

  report_integration_with_dir(function(output_dir) {
    report <- report_test_run(
      fit,
      output_dir = output_dir,
      output_file = "partial-outcome-failure",
      format = "html",
      render = FALSE,
      open = FALSE
    )
    saved <- readRDS(report$payload)$result
    manifest <- report_integration_manifest(report)
    context <- .mira_report_context_tree(saved, report_test_cfg())

    expect_s3_class(saved$outcomes$bad, "mira_info_error")
    expect_identical(saved$outcomes$bad$error, fit$outcomes$bad$error)
    expect_false("bad" %in% manifest$outcome)
    expect_false(any(grepl("/bad/", report_test_relpath(manifest$file),
                           fixed = TRUE)))
    expect_identical(
      context$diagnostics$children$outcome_errors$children$bad$value,
      fit$outcomes$bad$error
    )
    expect_length(.mira_report_outcome_tree(saved$outcomes$bad,
                                            report_test_cfg()), 0L)
  })
})

test_that("Quarto renders a minimal HTML report when its package and CLI are available", {
  skip_if_not_installed("knitr")
  skip_if_not_installed("quarto")
  quarto_cli <- tryCatch(
    suppressWarnings(quarto::quarto_path()),
    error = function(e) NULL
  )
  skip_if(
    is.null(quarto_cli) || length(quarto_cli) != 1L ||
      is.na(quarto_cli) || !nzchar(quarto_cli),
    "Quarto CLI is not available"
  )

  fit <- report_integration_fit()
  report_integration_with_dir(function(output_dir) {
    report <- report_test_run(
      fit,
      output_dir = output_dir,
      output_file = "quarto-html-integration",
      format = "html",
      render = TRUE,
      open = FALSE,
      quiet = TRUE
    )

    expect_true(unname(report$rendered[["html"]]))
    expect_true(file.exists(report$expected_files[["html"]]))
    expect_identical(unname(report$files),
                     unname(report$expected_files[["html"]]))
    expect_true(is.na(report$render_errors[["html"]]))
  })
})
