test_that("old rendered files cannot make a no-op renderer look successful", {
  root <- withr::local_tempdir(pattern = "mira-stale-render-")
  old_path <- file.path(root, "mira_report.html")
  writeLines("old report sentinel", old_path)
  opened <- character()
  testthat::local_mocked_bindings(
    quarto_render = function(...) invisible(NULL), .package = "quarto"
  )
  testthat::local_mocked_bindings(
    browseURL = function(url, ...) opened <<- c(opened, url), .package = "utils"
  )
  expect_warning(
    report <- report_test_call(list(render = TRUE, open = TRUE, overwrite = TRUE), root),
    "without creating the expected html file", fixed = TRUE
  )
  expect_identical(report$rendered, c(html = FALSE))
  expect_length(report$files, 0L)
  expect_length(opened, 0L)
  expect_identical(readLines(old_path), "old report sentinel")
})

test_that("a partial renderer failure preserves a previous successful document", {
  root <- withr::local_tempdir(pattern = "mira-partial-old-render-")
  old_path <- file.path(root, "mira_report.html")
  writeLines("old report sentinel", old_path)
  testthat::local_mocked_bindings(
    quarto_render = function(output_file, execute_dir, ...) {
      writeLines("incomplete document", file.path(execute_dir, output_file))
      stop("controlled failure after output", call. = FALSE)
    }, .package = "quarto"
  )
  expect_warning(
    report <- report_test_call(list(render = TRUE, overwrite = TRUE), root),
    "controlled failure after output", fixed = TRUE
  )
  expect_identical(report$rendered, c(html = FALSE))
  expect_length(report$files, 0L)
  expect_identical(readLines(old_path), "old report sentinel")
  expect_length(list.files(root, pattern = "render-"), 0L)
})

test_that("scientific headings and metadata are emitted as literal labels", {
  label <- "A <control> *dose* [site](x) & $change_1"
  heading <- paste(capture.output(.mira_report_heading_long(label)), collapse = "\n")
  expect_match(heading, "A \\<control\\> \\*dose\\* \\[site\\]\\(x\\) \\& \\$change\\_1",
               fixed = TRUE)
  node <- list(kind = "value", key = "arm_reference", value = label)
  rendered <- paste(capture.output(.mira_report_render_node_long(
    node, 2L, report_test_cfg()
  )), collapse = "\n")
  expect_match(rendered, "**Arm reference:** A \\<control\\> \\*dose\\*", fixed = TRUE)
})

test_that("failed advanced modules surface error diagnostics without stale results", {
  input <- report_test_mira()
  input$advanced_models <- list(gee = list(ar1 = list(
    performed = FALSE, error = "controlled AR1 estimation failure",
    reason_skipped = "Too few eligible clusters.",
    warnings = "controlled sparse design warning",
    fixed_effects = data.frame(term = "stale result", estimate = 999)
  )))
  nodes <- .mira_report_outcome_tree_long(input, report_test_cfg())
  expect_null(nodes$advanced_models)
  diagnostic <- if (is.null(nodes$diagnostics)) "" else paste(capture.output(
    .mira_report_render_node_long(nodes$diagnostics, 2L, report_test_cfg())
  ), collapse = "\n")
  expect_match(diagnostic, "controlled AR1 estimation failure", fixed = TRUE)
  expect_match(diagnostic, "Too few eligible clusters.", fixed = TRUE)
  expect_match(diagnostic, "controlled sparse design warning", fixed = TRUE)
  expect_false(grepl("stale result|999", diagnostic))
})

test_that("wrong output paths and promotion failures never validate an old report", {
  root <- withr::local_tempdir(pattern = "mira-render-output-faults-")
  for (fault in c("wrong-path", "promotion")) {
    directory <- file.path(root, fault)
    dir.create(directory)
    old_path <- file.path(directory, "mira_report.html")
    writeLines("old report sentinel", old_path)
    testthat::with_mocked_bindings({
      if (identical(fault, "wrong-path")) {
        expect_warning(report <- report_test_call(
          list(render = TRUE, overwrite = TRUE), directory
        ), "without creating the expected html file", fixed = TRUE)
      } else {
        testthat::with_mocked_bindings(
          expect_warning(report <- report_test_call(
            list(render = TRUE, overwrite = TRUE), directory
          ), "Could not place the rendered html file", fixed = TRUE),
          file.copy = function(...) FALSE, .package = "base"
        )
      }
    }, quarto_render = function(output_file, execute_dir, ...) {
      target <- if (identical(fault, "wrong-path")) "elsewhere.html" else output_file
      writeLines("new document", file.path(execute_dir, target))
    }, .package = "quarto")
    expect_identical(report$rendered, c(html = FALSE), info = fault)
    expect_length(report$files, 0L)
    expect_identical(readLines(old_path), "old report sentinel", info = fault)
    expect_length(list.files(directory, pattern = "mira_report-render-"), 0L)
  }
})

test_that("freshly rendered identical content succeeds and leaves no staging file", {
  root <- withr::local_tempdir(pattern = "mira-render-identical-")
  path <- file.path(root, "mira_report.html")
  writeLines("same valid document", path)
  testthat::local_mocked_bindings(
    quarto_render = function(output_file, execute_dir, ...) {
      writeLines("same valid document", file.path(execute_dir, output_file))
    }, .package = "quarto"
  )
  report <- report_test_call(list(render = TRUE, overwrite = TRUE), root)
  expect_identical(report$rendered, c(html = TRUE))
  expect_identical(report$files, report$expected_files)
  expect_true(is.na(report$render_errors[["html"]]))
  expect_identical(readLines(path), "same valid document")
  expect_length(list.files(root, pattern = "mira_report-render-"), 0L)
})

test_that("a directory conflict cannot destroy an unrelated file", {
  root <- withr::local_tempdir(pattern = "mira-directory-conflict-")
  sentinel <- file.path(root, "blocked")
  writeLines("unrelated data", sentinel)
  expect_error(
    report_test_call(output_dir = file.path(sentinel, "nested")),
    "Cannot create directory:", fixed = TRUE
  )
  expect_identical(readLines(sentinel), "unrelated data")
  expect_identical(list.files(root), "blocked")
})

test_that("serialization and CSV faults report their stage without deleting foreign files", {
  root <- withr::local_tempdir(pattern = "mira-export-faults-")
  directory <- file.path(root, "serialization")
  dir.create(directory)
  sentinel <- file.path(directory, "manual.csv")
  writeLines("unrelated observations", sentinel)
  isolated <- new.env(parent = environment(mira_report_freq_long))
  isolated$saveRDS <- function(...) stop("controlled serialization failure", call. = FALSE)
  local_report <- mira_report_freq_long
  environment(local_report) <- isolated
  expect_error(capture.output(local_report(
    report_test_mira(), output_dir = directory, render = FALSE, open = FALSE
  )), "controlled serialization failure", fixed = TRUE)
  expect_identical(readLines(sentinel), "unrelated observations")
  expect_identical(list.files(directory), "manual.csv")

  directory <- file.path(root, "csv")
  dir.create(directory)
  sentinel <- file.path(directory, "manual.csv")
  writeLines("unrelated observations", sentinel)
  original_write <- utils::write.csv
  testthat::local_mocked_bindings(
    write.csv = function(x, file, ...) {
      if (identical(basename(file), "descriptives.csv")) {
        stop("controlled scientific CSV failure", call. = FALSE)
      }
      original_write(x, file, ...)
    }, .package = "utils"
  )
  expect_error(report_test_call(output_dir = directory),
               "Could not export table at outcome$descriptives to results/descriptives.csv: controlled scientific CSV failure",
               fixed = TRUE)
  expect_identical(readLines(sentinel), "unrelated observations")
  # This is an explicitly partial bundle: no transactional guarantee is
  # inferred. Source/payload remain inspectable and no complete manifest exists.
  expect_true(all(file.exists(file.path(directory, c(
    "mira_report.qmd", "mira_report-runtime.R", "mira_report-payload.rds"
  )))))
  expect_identical(readRDS(file.path(directory, "mira_report-payload.rds"))$result,
                   report_test_mira())
  expect_false(file.exists(file.path(directory, "results_manifest.csv")))
})

test_that("explicit real-render profiles cannot silently skip a missing CLI", {
  withr::local_envvar(c(MIRA_RUN_QUARTO_TESTS = "true"))
  testthat::local_mocked_bindings(quarto_path = function() NULL, .package = "quarto")
  expect_error(report_test_require_quarto(),
               "MIRA_RUN_QUARTO_TESTS requires an existing Quarto CLI", fixed = TRUE)
})

test_that("published HTML cells retain independent values and literal scientific labels", {
  input <- data.frame(time = c("A < B", "B & C"), n = c(2L, 3L),
                      mean = c(1.234567, 2.345678), p.value = c(0.00001, 0.123456))
  rendered <- testthat::with_mocked_bindings(
    capture.output(.mira_report_table_long(input, caption = "Scientific results",
                                         cfg = report_test_cfg())),
    is_html_output = function() TRUE, is_latex_output = function() FALSE,
    .package = "knitr"
  )
  published <- report_test_html_table(paste(rendered, collapse = "\n"), "Scientific results")
  expect_identical(published$time, input$time)
  expect_equal(as.numeric(published$n), c(2, 3), tolerance = 0)
  expect_equal(as.numeric(published$mean), c(1.235, 2.346), tolerance = 0)
  expect_identical(published[["p value"]], c("<0.001", "0.123"))
})
test_that("missingness headings describe the displayed grouping", {
  input <- report_test_mira()
  input$plots <- list(missingness = grid::rectGrob())
  cfg <- report_test_cfg(c("figures", "arms"))
  nodes <- .mira_report_outcome_tree_long(input, cfg)
  expect_identical(nodes$plots$children$missingness$kind, "plot")
  # Rendering the heading does not require drawing a figure in this regression.
  isolated <- new.env(parent = environment(.mira_report_render_node_long))
  isolated$.mira_report_plot_long <- function(...) invisible(TRUE)
  render_node <- .mira_report_render_node_long
  environment(render_node) <- isolated
  plot_output <- capture.output(render_node(nodes$plots$children$missingness, 3L, cfg))
  expect_true("### Missing data" %in% plot_output)
  expect_false(any(grepl("by arm", plot_output, fixed = TRUE)))

  input$arm_analysis <- list(enabled = TRUE, missingness = data.frame(
    arm = "Control", time = "t0", n_missing = 0L
  ))
  input$plots$arm_missingness <- grid::rectGrob()
  nodes <- .mira_report_outcome_tree_long(input, cfg)
  arm_output <- capture.output(.mira_report_render_node_long(nodes$arm_analysis, 2L, cfg))
  expect_true("## Arm-specific analyses" %in% arm_output)
  expect_true("### Missing data" %in% arm_output)
  expect_match(paste(arm_output, collapse = "\n"), "Control", fixed = TRUE)
  arm_plot <- capture.output(render_node(nodes$plots$children$arm_missingness, 3L, cfg))
  expect_true("### Missing data by arm" %in% arm_plot)
})

test_that("DOCX pipe output preserves separate table lines and scientific cells", {
  input <- data.frame(time = c("t0", "t1"), n = c(6L, 5L), mean = c(3.5, 4.5))
  rendered <- testthat::with_mocked_bindings(
    capture.output(.mira_report_table_long(input, caption = "Scientific results",
                                         cfg = report_test_cfg())),
    is_html_output = function() FALSE, is_latex_output = function() FALSE,
    .package = "knitr"
  )
  # Pandoc requires a header, delimiter, and each data row on distinct lines.
  # Inspect real knitr output here; no mocked formatter can hide concatenation.
  lines <- rendered[grepl("^\\|", rendered)]
  cells <- lapply(strsplit(lines, "|", fixed = TRUE), function(row) {
    trimws(row[-1L])
  })
  expect_length(lines, 4L)
  expect_identical(cells[c(1L, 3L, 4L)], list(
    c("time", "n", "mean"), c("t0", "6", "3.5"), c("t1", "5", "4.5")
  ))
  expect_true(any(grepl("^\\|:?-+", lines)))
})

