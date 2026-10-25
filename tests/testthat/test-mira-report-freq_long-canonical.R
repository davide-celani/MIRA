test_that("report metric resolution uses exact canonical references", {
  for (outcome in c("Measure_X", "mean_sd_change", "p_value_Signal")) {
    tab <- data.frame(time = c("t0", "t1"), n = c(6L, 5L),
                      mean = c(12.34567, 14.56789), p_raw = c(0.00001, 0.025),
                      p_holm = c(0.00002, 0.05), check.names = FALSE)
    map <- stats::setNames(paste0(outcome, "_", c("mean", "p_raw", "p_holm")),
                          c("mean", "p_raw", "p_holm"))
    names(tab)[3:5] <- unname(map)
    attr(tab, "mira_metric_names") <- map
    expect_identical(.mira_report_base_names_long(tab),
                     c("time", "n", "mean", "p_raw", "p_holm"))
    expect_identical(.mira_report_metric_name_long(tab, "mean"), unname(map[["mean"]]))
    expect_identical(.mira_report_metric_name_long(tab, "n"), "n")
    expect_null(.mira_report_metric_name_long(tab, "median"))
    formatted <- .mira_report_prepare_table_long(tab)
    expect_identical(formatted[[map[["mean"]]]], c(12.34567, 14.56789))
    expect_identical(formatted[[map[["p_raw"]]]], c("<0.001", "0.025"))
    expect_identical(formatted[[map[["p_holm"]]]], c("<0.001", "0.050"))
    expect_identical(names(formatted), names(tab))
    expect_identical(tab[[map[["p_raw"]]]], c(0.00001, 0.025))
    tab$p_value_subject <- c(12345, 23456)
    expect_identical(.mira_report_prepare_table_long(tab)$p_value_subject,
                     c(12345, 23456))
  }
})

test_that("bundles preserve each naming provenance through x and data", {
  root <- withr::local_tempdir(pattern = "mira-canonical-report-")
  for (outcome in c("Measure_X", "mean_sd_change", "p_value_Signal")) {
    data <- data.frame(patient = seq_len(6L), first = c(1, 3, 2, 6, 5, 4),
                       middle = c(2, 4, NA, 5, 7, 6),
                       last = c(4, 5, 6, 7, 8, 9))
    columns <- paste0(outcome, c("_t0", "_t1", "_t2"))
    names(data)[2:4] <- columns
    modes <- list(
      automatic = list(),
      time_columns = list(time_vars = columns),
      outcome_columns = list(outcomes = columns),
      pattern = list(variable_pattern = "^(.+)_t([0-9]+)$"),
      declared = list(outcomes = outcome),
      named_list = list(time_vars = stats::setNames(list(columns), outcome))
    )
    for (mode in names(modes)) {
      arguments <- c(list(id = "patient", covariates = character(),
                          analyses = "none", verbose = FALSE), modes[[mode]])
      fit <- do.call(mira_info_long, c(list(data = data), arguments))
      expected <- if (mode %in% c("declared", "named_list")) {
        paste0(outcome, "_mean")
      } else "mean"
      expect_true(expected %in% names(fit$descriptives), info = paste(outcome, mode))
      expect_identical(.mira_report_metric_name_long(fit$descriptives, "mean"), expected)
      for (input in c("x", "data")) {
        directory <- file.path(root, outcome, mode, input)
        supplied <- if (input == "x") list(x = fit) else c(list(data = data), arguments)
        capture.output(report <- do.call(mira_report_freq_long, c(supplied, list(
          output_dir = directory, format = "html", render = FALSE, open = FALSE,
          extra_objects = list(user_metadata = list(version = "user-release", value = 7))
        ))))
        payload <- readRDS(report$payload)
        expect_identical(payload$result, report$mira_result)
        expect_identical(payload$result$descriptives, fit$descriptives)
        expect_identical(payload$result$change, fit$change)
        expect_identical(payload$extra_objects$user_metadata$version, "user-release")
        expect_false("version" %in% names(payload$result))
        expect_true(all(file.exists(c(report$source, report$runtime, report$payload,
                                       report$results_manifest, report$result_files))))
        manifest <- utils::read.csv(report$results_manifest, stringsAsFactors = FALSE)
        expect_false(anyDuplicated(manifest$file) > 0L)
        expect_false(any(grepl("version", manifest$path, fixed = TRUE)))
        expect_false(any(grepl("multiplicity.*\\$table$", manifest$path)))
        for (component in c("descriptives", "change")) {
          source <- fit[[component]]
          row <- which(manifest$path == paste0("outcome$", component))
          expect_length(row, 1L)
          # CSV has no type schema: all-NA numeric columns and numeric-looking
          # time labels need their known source types when reading the archive.
          column_types <- vapply(source, function(column) {
            if (is.integer(column)) "integer" else if (is.numeric(column)) "numeric" else
              if (is.logical(column)) "logical" else "character"
          }, character(1L))
          csv <- utils::read.csv(file.path(directory, manifest$file[[row]]),
                                 check.names = FALSE, stringsAsFactors = FALSE,
                                 na.strings = "", colClasses = column_types)
          expect_identical(names(csv), names(source))
          expect_identical(nrow(csv), nrow(source))
          expect_equal(manifest$columns[[row]], ncol(source))
          expect_equal(manifest$rows[[row]], nrow(source))
          for (column in names(source)) {
            # Decimal CSV serialization may round the final binary digits.
            expect_equal(csv[[column]], source[[column]], tolerance = 1e-12,
                         info = paste(outcome, mode, input, component, column))
          }
        }
      }
    }
  }
})

test_that("mixed multi-outcome bundles retain independent naming and order", {
  data <- data.frame(patient = seq_len(6L),
                     Signal_t0 = c(1, 3, 2, 6, 5, 4),
                     Signal_t1 = c(2, 4, 3, 5, 7, 6),
                     inferred_t0 = c(9, 8, 7, 6, 5, 4),
                     inferred_t1 = c(8, 8, 6, 5, 4, 3))
  arguments <- list(id = "patient", time_vars = stats::setNames(
    list(c("Signal_t0", "Signal_t1"), c("inferred_t0", "inferred_t1")),
    c("Declared_mean", "")), covariates = character(), analyses = "none", verbose = FALSE)
  fit <- do.call(mira_info_long, c(list(data = data), arguments))
  expect_s3_class(fit, "mira_info_multi_long")
  expect_identical(names(fit$outcomes), c("Declared_mean", "inferred"))
  expect_true("Declared_mean_mean" %in% names(fit$outcomes[[1L]]$descriptives))
  expect_true("mean" %in% names(fit$outcomes[[2L]]$descriptives))
  root <- withr::local_tempdir(pattern = "mira-mixed-canonical-")
  for (input in c("x", "data")) {
    supplied <- if (input == "x") list(x = fit) else c(list(data = data), arguments)
    capture.output(report <- do.call(mira_report_freq_long, c(supplied, list(
      output_dir = file.path(root, input), render = FALSE, open = FALSE, format = "html"
    ))))
    saved <- readRDS(report$payload)$result
    expect_identical(saved$outcomes, fit$outcomes)
    manifest <- utils::read.csv(report$results_manifest, stringsAsFactors = FALSE)
    expect_identical(unique(manifest$outcome[nzchar(manifest$outcome)]), names(fit$outcomes))
    for (i in seq_along(fit$outcomes)) {
      row <- which(manifest$outcome == names(fit$outcomes)[[i]] &
                     manifest$path == "outcome$descriptives")
      expect_length(row, 1L)
      csv <- utils::read.csv(file.path(report$output_dir, manifest$file[[row]]),
                             check.names = FALSE)
      expect_identical(names(csv), names(fit$outcomes[[i]]$descriptives))
    }
  }
})

test_that("isolated saved runtimes execute prefixed conclusions and diagnostics", {
  skip_if_not_installed("knitr")
  root <- withr::local_tempdir(pattern = "mira-canonical-runtime-")
  old_options <- knitr::opts_knit$get()
  old_chunks <- knitr::opts_chunk$get()
  withr::defer(knitr::opts_knit$restore(old_options))
  withr::defer(knitr::opts_chunk$restore(old_chunks))
  grDevices::pdf(file.path(root, "knit-device.pdf"))
  device <- grDevices::dev.cur()
  withr::defer(grDevices::dev.off(device))
  knitr::opts_knit$set(global.device = TRUE)
  data <- data.frame(patient = seq_len(6L),
                     p_value_Signal_t0 = c(1, 3, 2, 6, 5, 4),
                     p_value_Signal_t1 = c(2, 4, 3, 5, 7, 6))
  for (explicit in c(FALSE, TRUE)) {
    fit <- mira_info_long(data, id = "patient",
      outcomes = if (explicit) "p_value_Signal" else NULL,
      covariates = character(), analyses = "none", verbose = FALSE)
    capture.output(report <- mira_report_freq_long(fit,
      output_dir = file.path(root, as.character(explicit)), format = "html",
      render = FALSE, open = FALSE, sections = c("descriptives", "outliers", "conclusions")))
    isolated <- new.env(parent = baseenv())
    sys.source(report$runtime, envir = isolated)
    payload <- readRDS(report$payload)
    expect_true(exists(".mira_report_metric_name_long", isolated, inherits = FALSE))
    expect_true(exists(".mira_report_base_names_long", isolated, inherits = FALSE))
    output <- capture.output(isolated$.mira_report_knit_long(payload))
    markdown_path <- file.path(report$output_dir, "executed.md")
    expect_no_error(knitr::knit(report$source, output = markdown_path,
                               envir = new.env(parent = baseenv()), quiet = TRUE))
    markdown <- paste(readLines(markdown_path, warn = FALSE), collapse = "\n")
    expect_match(markdown, "# Descriptive summary", fixed = TRUE)
    expect_match(markdown, "Differences between observed means", fixed = TRUE)
    expect_false(grepl("Error in|could not find function|cannot open", markdown))
    expect_true(any(grepl("# Descriptive summary", output, fixed = TRUE)))
    expect_true(any(grepl("Differences between observed means", output, fixed = TRUE)))
    expected <- if (explicit) "p value Signal mean" else "mean"
    expect_true(any(grepl(expected, output, fixed = TRUE)))
    expect_false(any(grepl("Error in|could not find function", output)))
    nodes <- isolated$.mira_report_outcome_tree_long(payload$result, payload$report)
    expect_true("descriptives" %in% names(nodes))
    expect_identical(names(nodes$descriptives$value), names(fit$descriptives))
  }
})

test_that("multiplicity metadata points to the canonical table without restoring aliases", {
  fit <- mira_info_long(data.frame(patient = 1:6,
    Signal_t0 = c(1, 3, 2, 6, 5, 4), Signal_t1 = c(2, 4, 3, 5, 7, 6)),
    id = "patient", outcomes = "Signal", covariates = character(),
    analyses = "none", verbose = FALSE)
  family <- fit$multiplicity$friedman_posthoc
  expect_false("table" %in% names(family))
  expect_identical(family$object_path, "$advanced_tests$friedman$posthoc$tidy")
  node <- .mira_report_outcome_tree_long(fit, report_test_cfg())$multiplicity
  expect_identical(node$children$friedman_posthoc$children$object_path$value,
                   family$object_path)
  paths <- vapply(.mira_report_collect_tables_long(fit, report_test_cfg()),
                  `[[`, character(1L), "path")
  expect_true("outcome$advanced_tests$friedman$posthoc$tidy" %in% paths)
  expect_false(any(grepl("multiplicity.*\\$table$", paths)))
})
