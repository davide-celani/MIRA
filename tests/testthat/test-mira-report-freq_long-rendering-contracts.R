test_that("missing rendering packages retain an inspectable source bundle", {
  root <- withr::local_tempdir(pattern = "mira-render-dependencies-")
  for (missing_package in c("knitr", "quarto")) {
    isolated <- new.env(parent = environment(mira_report_freq_long))
    isolated$requireNamespace <- local({
      unavailable <- missing_package
      function(package, quietly = FALSE, ...) !identical(package, unavailable)
    })
    report_function <- mira_report_freq_long
    environment(report_function) <- isolated
    directory <- file.path(root, missing_package)
    expect_error(
      capture.output(report_function(
        report_test_mira(), output_dir = directory, format = "html",
        render = TRUE, open = FALSE
      )),
      paste0("The R package '", missing_package,
             "' is required. Source files were created in:"), fixed = TRUE
    )
    expect_true(all(file.exists(file.path(directory, c(
      "mira_report.qmd", "mira_report-runtime.R", "mira_report-payload.rds",
      "results_manifest.csv"
    )))))
    expect_gt(length(parse(file = file.path(directory, "mira_report-runtime.R"))), 0L)
    expect_identical(readRDS(file.path(directory, "mira_report-payload.rds"))$result,
                     report_test_mira())
    expect_false(file.exists(file.path(directory, "mira_report.html")))
  }
})

test_that("complete render failures retain sources and never open old rendered files", {
  skip_if_not_installed("knitr")
  skip_if_not_installed("quarto")
  root <- withr::local_tempdir(pattern = "mira-render-failed-")
  formats <- c("docx", "html", "pdf")
  old_paths <- file.path(root, paste0("mira_report.", formats))
  for (path in old_paths) writeLines("old report", path)
  calls <- character()
  opened <- character()
  testthat::local_mocked_bindings(
    quarto_render = function(output_format, ...) {
      calls <<- c(calls, output_format)
      stop(paste("controlled failure for", output_format), call. = FALSE)
    }, .package = "quarto"
  )
  testthat::local_mocked_bindings(
    browseURL = function(url, ...) { opened <<- c(opened, url) }, .package = "utils"
  )
  expect_warning(
    report <- report_test_call(list(
      format = formats, render = TRUE, open = TRUE, overwrite = TRUE
    ), root),
    "Rendering did not complete for docx, html, pdf. Source files were retained.",
    fixed = TRUE
  )
  expect_identical(calls, formats)
  expect_identical(report$rendered, stats::setNames(rep(FALSE, 3L), formats))
  expect_identical(unname(report$render_errors), paste("controlled failure for", formats))
  expect_length(report$files, 0L)
  expect_length(opened, 0L)
  expect_true(all(file.exists(c(report$source, report$runtime, report$payload))))
  for (path in old_paths) expect_identical(readLines(path), "old report")
})

test_that("successful formats follow request order and browser failure preserves the descriptor", {
  skip_if_not_installed("knitr")
  skip_if_not_installed("quarto")
  root <- withr::local_tempdir(pattern = "mira-render-order-")
  opened <- character()
  testthat::local_mocked_bindings(
    quarto_render = function(output_file, execute_dir, ...) {
      writeLines("successful rendered report", file.path(execute_dir, output_file))
    }, .package = "quarto"
  )
  testthat::local_mocked_bindings(
    browseURL = function(url, ...) {
      opened <<- c(opened, url)
      stop("browser unavailable", call. = FALSE)
    }, .package = "utils"
  )
  report <- report_test_call(list(
    format = c("DOCX", "html", "docx"), render = TRUE, open = TRUE
  ), root)
  expect_identical(report$formats, c("docx", "html"))
  expect_identical(report$rendered, c(docx = TRUE, html = TRUE))
  expect_identical(report$files, report$expected_files)
  expect_true(all(is.na(report$render_errors)))
  expect_identical(opened, unname(report$expected_files[["docx"]]))
  expect_true(all(file.info(report$files)$size > 0))
})

test_that("print returns the descriptor invisibly for complete and partial rendering", {
  root <- withr::local_tempdir(pattern = "mira-render-print-")
  descriptor <- report_test_call(list(format = c("html", "pdf")), root)
  # Printing is intentionally a descriptor operation; these files need not
  # exist and print must not inspect, open, or attempt to render them.
  descriptor$rendered[] <- TRUE
  descriptor$files <- descriptor$expected_files
  complete <- capture.output(value <- withVisible(print(descriptor)))
  expect_false(value$visible)
  expect_identical(value$value, descriptor)
  expect_match(paste(complete, collapse = "\n"), "Rendered:", fixed = TRUE)
  expect_false(any(grepl("Not rendered", complete, fixed = TRUE)))
  descriptor$rendered[["pdf"]] <- FALSE
  descriptor$files <- descriptor$expected_files["html"]
  partial <- capture.output(value <- withVisible(print(descriptor)))
  expect_false(value$visible)
  expect_identical(value$value, descriptor)
  expect_match(paste(partial, collapse = "\n"), "Not rendered successfully: pdf", fixed = TRUE)
  expect_match(paste(partial, collapse = "\n"), "mira_report.html", fixed = TRUE)
})

test_that("Quarto PDF integration requires an existing CLI and LaTeX toolchain", {
  skip_if_not_installed("knitr")
  skip_if_not_installed("quarto")
  cli <- tryCatch(suppressWarnings(quarto::quarto_path()), error = function(e) NULL)
  skip_if(is.null(cli) || length(cli) != 1L || is.na(cli) || !nzchar(cli),
          "Quarto CLI is not available")
  latex_engines <- Sys.which(c("pdflatex", "xelatex", "lualatex"))
  skip_if(!any(nzchar(latex_engines)),
          "No existing LaTeX executable is available; tests never install one")
  selected_engine <- names(latex_engines)[which(nzchar(latex_engines))[[1L]]]
  root <- withr::local_tempdir(pattern = "mira-quarto-pdf-")
  # Add the supported Quarto PDF setting to this test's source so that a
  # missing TeX package fails locally instead of triggering an installation.
  original_qmd <- .mira_report_qmd_lines_long
  testthat::local_mocked_bindings(
    .mira_report_qmd_lines_long = function(...) {
      lines <- original_qmd(...)
      append(lines, c("    latex-auto-install: false",
                      paste0("    pdf-engine: ", selected_engine)),
             after = which(lines == "  pdf:"))
    }, .package = "MIRA"
  )
  report <- report_test_call(list(
    x = report_test_mira(), format = "pdf", sections = "descriptives",
    render = TRUE, date = "2025-01-15"
  ), root)
  expect_identical(report$rendered, c(pdf = TRUE))
  expect_true(is.na(report$render_errors[["pdf"]]))
  expect_gt(file.info(report$files[["pdf"]])$size, 0)
  header <- readBin(report$files[["pdf"]], what = "raw", n = 5L)
  expect_identical(rawToChar(header), "%PDF-")
})
