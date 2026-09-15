test_that("render = FALSE creates a coherent reproducibility bundle", {
  output_dir <- tempfile("mira-report-files-")
  on.exit(unlink(output_dir, recursive = TRUE, force = TRUE), add = TRUE)

  input <- report_test_mira()
  extras <- list(
    audit = data.frame(term = "registered", value = 1.25,
                       stringsAsFactors = FALSE)
  )
  visible <- NULL
  printed <- utils::capture.output(
    visible <- withVisible(mira_report_freq(
      input,
      output_dir = output_dir,
      output_file = "Study Bundle",
      format = "all",
      title = "O'Brien \"trial\"\nphase 2",
      subtitle = "Path C:\\analysis and 'notes'",
      author = c("Ana O'Brien", "Bao \"B\""),
      date = as.Date("2025-01-15"),
      sections = c("Overview", "MODELS", "overview"),
      include_complete_output = FALSE,
      extra_objects = extras,
      max_table_rows = 5L,
      max_table_columns = 4L,
      digits = 4L,
      max_depth = 7L,
      render = FALSE,
      open = FALSE
    ))
  )
  report <- visible$value

  expect_false(visible$visible)
  expect_s3_class(report, "mira_report_freq")
  expect_identical(report$mira_result, input)
  expect_named(
    report,
    c("call", "mira_result", "formats", "rendered", "render_errors",
      "files", "expected_files", "source", "payload", "runtime",
      "results_manifest", "result_files", "output_dir", "report_config"),
    ignore.order = FALSE
  )
  expect_identical(report$formats, c("html", "pdf", "docx"))
  expect_identical(names(report$rendered), report$formats)
  expect_identical(names(report$expected_files), report$formats)
  expect_true(all(!report$rendered))
  expect_true(all(is.na(report$render_errors)))
  expect_length(report$files, 0L)
  expect_identical(
    basename(unname(report$expected_files)),
    c("study_bundle.html", "study_bundle.pdf", "study_bundle.docx")
  )
  expect_true(all(!file.exists(report$expected_files)))
  expect_true(all(file.exists(c(
    report$source, report$payload, report$runtime,
    report$results_manifest, report$result_files
  ))))
  expect_match(paste(printed, collapse = "\n"), "MIRA Quarto report")
  expect_match(paste(printed, collapse = "\n"),
               "Not rendered: html, pdf, docx")

  print_visible <- NULL
  print_output <- utils::capture.output(
    print_visible <- withVisible(print(report))
  )
  expect_false(print_visible$visible)
  expect_identical(print_visible$value, report)
  expect_match(paste(print_output, collapse = "\n"),
               "Quarto source: study_bundle\\.qmd")
  expect_match(paste(print_output, collapse = "\n"),
               "Not rendered: html, pdf, docx")

  qmd <- readLines(report$source, warn = FALSE, encoding = "UTF-8")
  expect_identical(qmd[[1L]], "---")
  expect_equal(sum(qmd == "---"), 2L)
  expect_true("title: 'O''Brien \"trial\" phase 2'" %in% qmd)
  expect_true("subtitle: 'Path C:\\analysis and ''notes'''" %in% qmd)
  expect_true("author:" %in% qmd)
  expect_true("  - 'Ana O''Brien'" %in% qmd)
  expect_true("  - 'Bao \"B\"'" %in% qmd)
  expect_true("date: '2025-01-15'" %in% qmd)
  expect_true("lang: en" %in% qmd)
  expect_equal(sum(qmd == "format:"), 1L)
  expect_equal(sum(qmd == "  html:"), 1L)
  expect_equal(sum(qmd == "  pdf:"), 1L)
  expect_equal(sum(qmd == "  docx:"), 1L)
  expect_true(all(c(
    "    theme: cosmo", "    documentclass: article",
    "    papersize: a4", "    geometry: margin=25mm",
    "    number-sections: true", "execute:", "  echo: false",
    "  warning: false", "  message: false", "  error: false",
    "#| label: setup", "#| label: mira-report", "#| results: asis",
    ".mira_report_knit(payload)"
  ) %in% qmd))
  expect_true(
    "source(\"study_bundle-runtime.R\", local = knitr::knit_global())" %in% qmd
  )
  expect_true(
    "payload <- readRDS(\"study_bundle-payload.rds\")" %in% qmd
  )
  if (requireNamespace("yaml", quietly = TRUE)) {
    delimiter <- which(qmd == "---")
    metadata <- yaml::yaml.load(paste(
      qmd[seq.int(delimiter[[1L]] + 1L, delimiter[[2L]] - 1L)],
      collapse = "\n"
    ))
    expect_identical(metadata$title, "O'Brien \"trial\" phase 2")
    expect_identical(metadata$subtitle, "Path C:\\analysis and 'notes'")
    expect_identical(unlist(metadata$author, use.names = FALSE),
                     c("Ana O'Brien", "Bao \"B\""))
    expect_named(metadata$format, c("html", "pdf", "docx"))
    expect_false(metadata$execute$error)
  }

  parsed <- parse(file = report$runtime)
  expect_gt(length(parsed), 0L)
  runtime_env <- new.env(parent = baseenv())
  expect_no_error(sys.source(report$runtime, envir = runtime_env))
  runtime_names <- .mira_report_runtime_names()
  expect_setequal(ls(runtime_env, all.names = TRUE), runtime_names)
  expect_true(all(vapply(runtime_names, function(name) {
    exists(name, envir = runtime_env, inherits = FALSE) &&
      is.function(get(name, envir = runtime_env, inherits = FALSE))
  }, logical(1L))))

  payload <- readRDS(report$payload)
  expect_named(payload, c("result", "extra_objects", "report"),
               ignore.order = FALSE)
  expect_identical(payload$result, input)
  expect_identical(payload$extra_objects, extras)
  expect_identical(payload$report, report$report_config)
  expect_identical(payload$report$sections, c("overview", "models"))
  expect_false(payload$report$include_complete_output)
  expect_identical(payload$report$max_table_rows, 5L)
  expect_identical(payload$report$max_table_columns, 4L)
  expect_identical(payload$report$digits, 4L)
  expect_identical(payload$report$max_depth, 7L)
  expect_match(
    payload$report$created_at,
    "^[0-9]{4}-[0-9]{2}-[0-9]{2} [0-9]{2}:[0-9]{2}:[0-9]{2} UTC$"
  )
})

test_that("CSV exports and the manifest preserve raw scientific values", {
  output_dir <- tempfile("mira-report-csv-")
  on.exit(unlink(output_dir, recursive = TRUE, force = TRUE), add = TRUE)

  input <- report_test_mira()
  report <- report_test_run(
    x = input,
    output_dir = output_dir,
    output_file = "scientific-values",
    format = "html",
    render = FALSE,
    open = FALSE
  )
  manifest <- utils::read.csv(
    report$results_manifest, stringsAsFactors = FALSE, check.names = FALSE
  )

  expect_named(
    manifest,
    c("outcome", "path", "file", "original_class", "rows", "columns"),
    ignore.order = FALSE
  )
  expect_equal(nrow(manifest), length(report$result_files))
  expect_identical(anyDuplicated(tolower(manifest$file)), 0L)
  expect_true(all(grepl(
    "^results/([A-Za-z0-9_-]+/)*[A-Za-z0-9_-]+\\.csv$",
    manifest$file
  )))
  manifest_files <- file.path(output_dir, manifest$file)
  expect_true(all(file.exists(manifest_files)))
  expect_setequal(
    normalizePath(manifest_files, winslash = "/", mustWork = TRUE),
    normalizePath(report$result_files, winslash = "/", mustWork = TRUE)
  )

  for (i in seq_len(nrow(manifest))) {
    exported <- utils::read.csv(
      manifest_files[[i]], stringsAsFactors = FALSE, check.names = FALSE
    )
    expect_equal(nrow(exported), manifest$rows[[i]], info = manifest$path[[i]])
    expect_equal(ncol(exported), manifest$columns[[i]], info = manifest$path[[i]])
    expect_true(report_test_atomic_table(exported), info = manifest$path[[i]])
  }

  desc_index <- match("outcome$descriptives", manifest$path)
  expect_false(is.na(desc_index))
  expect_identical(manifest$original_class[[desc_index]], "data.frame")
  descriptives <- utils::read.csv(
    manifest_files[[desc_index]], stringsAsFactors = FALSE, check.names = FALSE
  )
  expect_type(descriptives[["p.value"]], "double")
  expect_equal(descriptives[["p.value"]], input$descriptives[["p.value"]])
  expect_equal(descriptives$mean, input$descriptives$mean)
  expect_false(any(startsWith(as.character(descriptives[["p.value"]]), "<")))

  change_index <- match("outcome$change", manifest$path)
  expect_false(is.na(change_index))
  change <- utils::read.csv(
    manifest_files[[change_index]], stringsAsFactors = FALSE, check.names = FALSE
  )
  expect_type(change$p_value, "double")
  expect_equal(change$p_value, input$change$p_value)

  correlation_index <- match("outcome$correlations$pearson", manifest$path)
  expect_false(is.na(correlation_index))
  correlation <- utils::read.csv(
    manifest_files[[correlation_index]], stringsAsFactors = FALSE,
    check.names = FALSE
  )
  expect_identical(names(correlation), c("row", "t0", "t1"))
  expect_identical(correlation$row, c("t0", "t1"))
  expect_equal(unname(as.matrix(correlation[c("t0", "t1")])),
               unname(input$correlations$pearson))
})

test_that("a display row limit never truncates the scientific CSV archive", {
  output_dir <- tempfile("mira-report-row-limit-")
  on.exit(unlink(output_dir, recursive = TRUE, force = TRUE), add = TRUE)

  input <- report_test_mira()
  input$large_table <- data.frame(
    id = seq_len(10L), estimate = seq(0.1, 1, length.out = 10L),
    p.value = seq(0.001, 0.01, length.out = 10L), check.names = FALSE
  )
  report <- report_test_run(
    input, output_dir = output_dir, format = "html",
    max_table_rows = 5L, render = FALSE, open = FALSE
  )
  manifest <- utils::read.csv(report$results_manifest, stringsAsFactors = FALSE)
  index <- match("outcome$large_table", manifest$path)
  expect_false(is.na(index))
  expect_identical(manifest$rows[[index]], 10L)
  archived <- utils::read.csv(file.path(output_dir, manifest$file[[index]]))
  expect_identical(nrow(archived), 10L)
  expect_equal(archived$estimate, input$large_table$estimate)
  expect_type(archived$p.value, "double")

  node <- .mira_report_prune(
    input$large_table, "large_table", "outcome$large_table",
    report$report_config
  )
  display <- capture.output(
    displayed <- .mira_report_render_node(node, 2L, report$report_config)
  )
  expect_null(displayed)
  expect_match(paste(display, collapse = "\n"), "Showing 5 of 10 rows",
               fixed = TRUE)
})

test_that("case-insensitive export collisions remain distinct in outcome directories", {
  output_dir <- tempfile("mira-report-collisions-")
  on.exit(unlink(output_dir, recursive = TRUE, force = TRUE), add = TRUE)

  input <- report_test_multi()
  input$outcomes$score[["Case Key"]] <- data.frame(value = 1)
  input$outcomes$score[["case_key"]] <- data.frame(value = 2)
  report <- report_test_run(
    x = input,
    output_dir = output_dir,
    format = "html",
    render = FALSE,
    open = FALSE
  )
  manifest <- utils::read.csv(
    report$results_manifest, stringsAsFactors = FALSE, check.names = FALSE
  )

  expect_identical(anyDuplicated(tolower(manifest$file)), 0L)
  for (outcome in c("score", "fatigue")) {
    outcome_files <- manifest$file[manifest$outcome == outcome]
    expect_gt(length(outcome_files), 0L)
    expect_true(all(startsWith(outcome_files, paste0("results/", outcome, "/"))))
  }
  collisions <- manifest[
    manifest$outcome == "score" &
      manifest$path %in% c("outcome$Case Key", "outcome$case_key"),
    , drop = FALSE
  ]
  expect_equal(nrow(collisions), 2L)
  expect_setequal(
    basename(collisions$file),
    c("case_key.csv", "case_key_2.csv")
  )
  expect_true(all(file.exists(file.path(output_dir, manifest$file))))
})

test_that("overwrite protects and then replaces every reproducibility artifact", {
  output_dir <- tempfile("mira-report-overwrite-")
  on.exit(unlink(output_dir, recursive = TRUE, force = TRUE), add = TRUE)

  first_input <- report_test_mira()
  first <- report_test_run(
    x = first_input,
    output_dir = output_dir,
    output_file = "replaceable",
    title = "First report",
    format = "html",
    render = FALSE,
    open = FALSE
  )
  qmd_before <- readLines(first$source, warn = FALSE, encoding = "UTF-8")

  second_input <- report_test_mira()
  second_input$descriptives$mean[[1L]] <- 999
  expect_error(
    report_test_run(
      x = second_input,
      output_dir = output_dir,
      output_file = "replaceable",
      title = "Blocked replacement",
      format = "html",
      render = FALSE,
      open = FALSE,
      overwrite = FALSE
    ),
    "Output files already exist.*overwrite=TRUE"
  )
  expect_identical(
    readLines(first$source, warn = FALSE, encoding = "UTF-8"),
    qmd_before
  )

  second <- report_test_run(
    x = second_input,
    output_dir = output_dir,
    output_file = "replaceable",
    title = "Replacement report",
    format = "html",
    render = FALSE,
    open = FALSE,
    overwrite = TRUE
  )
  expect_true("title: 'Replacement report'" %in%
                readLines(second$source, warn = FALSE, encoding = "UTF-8"))
  expect_identical(readRDS(second$payload)$result$descriptives$mean[[1L]], 999)
  expect_gt(length(parse(file = second$runtime)), 0L)

  manifest <- utils::read.csv(
    second$results_manifest, stringsAsFactors = FALSE, check.names = FALSE
  )
  desc_index <- match("outcome$descriptives", manifest$path)
  descriptives <- utils::read.csv(
    file.path(output_dir, manifest$file[[desc_index]]),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  expect_identical(descriptives$mean[[1L]], 999)
  expect_true(all(file.exists(c(
    second$source, second$payload, second$runtime,
    second$results_manifest, second$result_files
  ))))
})

test_that("stale cleanup removes only safe files recorded by the old manifest", {
  output_dir <- tempfile("mira-report-stale-")
  on.exit(unlink(output_dir, recursive = TRUE, force = TRUE), add = TRUE)

  original <- report_test_mira()
  original$stale_results <- data.frame(term = "obsolete", estimate = 1)
  first <- report_test_run(
    x = original,
    output_dir = output_dir,
    format = "html",
    render = FALSE,
    open = FALSE
  )
  previous <- utils::read.csv(
    first$results_manifest, stringsAsFactors = FALSE, check.names = FALSE
  )
  stale_index <- match("outcome$stale_results", previous$path)
  expect_false(is.na(stale_index))
  stale_file <- file.path(output_dir, previous$file[[stale_index]])
  expect_true(file.exists(stale_file))

  manual_file <- file.path(output_dir, "results", "manual.csv")
  outside_results <- file.path(output_dir, "outside.csv")
  writeLines("must remain", manual_file, useBytes = TRUE)
  writeLines("must remain", outside_results, useBytes = TRUE)

  unsafe <- previous[rep(1L, 2L), , drop = FALSE]
  unsafe$outcome <- ""
  unsafe$path <- c("tampered_parent", "tampered_absolute")
  unsafe$file <- c(
    "results/../outside.csv",
    report_test_relpath(normalizePath(outside_results, winslash = "/",
                                     mustWork = TRUE))
  )
  unsafe$original_class <- "data.frame"
  unsafe$rows <- 1L
  unsafe$columns <- 1L
  utils::write.csv(
    rbind(previous, unsafe), first$results_manifest,
    row.names = FALSE, fileEncoding = "UTF-8", na = ""
  )

  second <- report_test_run(
    x = report_test_mira(),
    output_dir = output_dir,
    format = "html",
    render = FALSE,
    open = FALSE,
    overwrite = TRUE
  )
  expect_false(file.exists(stale_file))
  expect_true(file.exists(manual_file))
  expect_true(file.exists(outside_results))

  current <- utils::read.csv(
    second$results_manifest, stringsAsFactors = FALSE, check.names = FALSE
  )
  expect_false(any(current$path %in% c(
    "outcome$stale_results", "tampered_parent", "tampered_absolute"
  )))
  expect_true(all(file.exists(file.path(output_dir, current$file))))
})

test_that("stale cleanup cannot follow a results link outside the report", {
  output_dir <- tempfile("mira-report-link-safety-")
  outside_dir <- tempfile("mira-report-outside-")
  on.exit(unlink(output_dir, recursive = TRUE, force = TRUE), add = TRUE)
  on.exit(unlink(outside_dir, recursive = TRUE, force = TRUE), add = TRUE)
  dir.create(file.path(output_dir, "results"), recursive = TRUE)
  dir.create(outside_dir, recursive = TRUE)

  link <- file.path(output_dir, "results", "external")
  linked <- suppressWarnings(tryCatch({
    if (.Platform$OS.type == "windows" && exists("Sys.junction", mode = "function")) {
      Sys.junction(outside_dir, link)
    } else {
      file.symlink(outside_dir, link)
    }
  }, error = function(e) FALSE))
  skip_if_not(isTRUE(linked) && dir.exists(link),
              "This platform cannot create a temporary directory link")

  victim <- file.path(outside_dir, "victim.csv")
  writeLines("scientific data", victim, useBytes = TRUE)
  manifest_path <- file.path(output_dir, "results_manifest.csv")
  utils::write.csv(
    data.frame(file = "results/external/victim.csv"), manifest_path,
    row.names = FALSE
  )

  expect_false(.mira_report_path_is_within(
    file.path(link, "victim.csv"), file.path(output_dir, "results")
  ))
  report_test_run(
    report_test_mira(), output_dir = output_dir, format = "html",
    render = FALSE, open = FALSE, overwrite = TRUE
  )
  expect_true(file.exists(victim))
  expect_identical(readLines(victim, warn = FALSE), "scientific data")
})

test_that("a result without tables produces a typed empty manifest", {
  output_dir <- tempfile("mira-report-empty-manifest-")
  on.exit(unlink(output_dir, recursive = TRUE, force = TRUE), add = TRUE)

  input <- structure(
    list(outcome = "score", note = "No tabular results are available."),
    class = c("mira_info", "list")
  )
  report <- report_test_run(
    x = input,
    output_dir = output_dir,
    format = "html",
    render = FALSE,
    open = FALSE
  )
  manifest <- utils::read.csv(
    report$results_manifest, stringsAsFactors = FALSE, check.names = FALSE
  )

  expect_equal(nrow(manifest), 0L)
  expect_named(
    manifest,
    c("outcome", "path", "file", "original_class", "rows", "columns"),
    ignore.order = FALSE
  )
  expect_identical(report$result_files, character())
  expect_true(all(file.exists(c(report$source, report$payload, report$runtime,
                                report$results_manifest))))
})

test_that("partial rendering records failures and opens only the first success", {
  if (!exists("local_mocked_bindings", envir = asNamespace("testthat"),
              inherits = FALSE)) {
    skip("testthat::local_mocked_bindings() is unavailable")
  }
  output_dir <- tempfile("mira-report-render-mock-")
  on.exit(unlink(output_dir, recursive = TRUE, force = TRUE), add = TRUE)

  render_calls <- list()
  opened <- character()
  testthat::local_mocked_bindings(
    quarto_render = function(input, output_format, output_file, execute_dir,
                             quiet, as_job, ...) {
      render_calls[[length(render_calls) + 1L]] <<- list(
        input = input, format = output_format, output_file = output_file,
        execute_dir = execute_dir, quiet = quiet, as_job = as_job
      )
      if (identical(output_format, "pdf")) {
        stop("controlled PDF failure", call. = FALSE)
      }
      target <- file.path(execute_dir, output_file)
      if (!file.create(target)) stop("could not create mocked output")
      invisible(target)
    },
    .package = "quarto"
  )
  testthat::local_mocked_bindings(
    browseURL = function(url, ...) {
      opened <<- c(opened, url)
      invisible(TRUE)
    },
    .package = "utils"
  )

  expect_warning(
    report <- report_test_run(
      x = report_test_mira(),
      output_dir = output_dir,
      output_file = "partial",
      format = c("html", "pdf"),
      render = TRUE,
      open = TRUE,
      quiet = FALSE
    ),
    "Rendering did not complete for pdf.*controlled PDF failure"
  )

  expect_identical(vapply(render_calls, `[[`, character(1L), "format"),
                   c("html", "pdf"))
  expect_true(all(vapply(render_calls, `[[`, logical(1L), "as_job") == FALSE))
  expect_true(all(vapply(render_calls, `[[`, logical(1L), "quiet") == FALSE))
  expect_identical(report$rendered, c(html = TRUE, pdf = FALSE))
  expect_true(is.na(report$render_errors[["html"]]))
  expect_match(report$render_errors[["pdf"]], "controlled PDF failure")
  expect_identical(report$files, report$expected_files["html"])
  expect_true(file.exists(report$expected_files[["html"]]))
  expect_false(file.exists(report$expected_files[["pdf"]]))
  expect_identical(opened, unname(report$expected_files[["html"]]))
})

test_that("a renderer that produces no file records an actionable error", {
  output_dir <- tempfile("mira-report-render-no-file-")
  on.exit(unlink(output_dir, recursive = TRUE, force = TRUE), add = TRUE)
  opened <- character()
  testthat::local_mocked_bindings(
    quarto_render = function(...) invisible(NULL),
    .package = "quarto"
  )
  testthat::local_mocked_bindings(
    browseURL = function(url, ...) {
      opened <<- c(opened, url)
      invisible(TRUE)
    },
    .package = "utils"
  )

  expect_warning(
    report <- report_test_run(
      report_test_mira(), output_dir = output_dir, format = "html",
      render = TRUE, open = TRUE
    ),
    "without creating the expected html file",
    fixed = TRUE
  )
  expect_identical(report$rendered, c(html = FALSE))
  expect_match(report$render_errors[["html"]],
               "without creating the expected html file", fixed = TRUE)
  expect_length(report$files, 0L)
  expect_length(opened, 0L)
})

test_that("open = TRUE does not browse when rendering is disabled", {
  if (!exists("local_mocked_bindings", envir = asNamespace("testthat"),
              inherits = FALSE)) {
    skip("testthat::local_mocked_bindings() is unavailable")
  }
  output_dir <- tempfile("mira-report-no-open-")
  on.exit(unlink(output_dir, recursive = TRUE, force = TRUE), add = TRUE)

  opened <- character()
  testthat::local_mocked_bindings(
    browseURL = function(url, ...) {
      opened <<- c(opened, url)
      invisible(TRUE)
    },
    .package = "utils"
  )
  report <- report_test_run(
    x = report_test_mira(),
    output_dir = output_dir,
    format = "html",
    render = FALSE,
    open = TRUE
  )

  expect_false(any(report$rendered))
  expect_length(opened, 0L)
})
