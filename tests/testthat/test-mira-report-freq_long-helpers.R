test_that("small scalar and presentation helpers handle empty and unusual values", {
  expect_identical(.mira_report_or_long(NULL, "fallback"), "fallback")
  expect_identical(.mira_report_or_long(character(), "fallback"), "fallback")
  expect_identical(.mira_report_or_long(FALSE, TRUE), FALSE)
  expect_identical(.mira_report_or_long(0, 1), 0)

  empty_table <- table(factor(character(), levels = character()))
  expect_true(.mira_report_is_empty_long(NULL))
  expect_true(.mira_report_is_empty_long(empty_table))
  expect_true(.mira_report_is_empty_long(data.frame()))
  expect_true(.mira_report_is_empty_long(matrix(numeric(), nrow = 0L)))
  expect_true(.mira_report_is_empty_long(list()))
  expect_false(.mira_report_is_empty_long(new.env()))
  expect_false(.mira_report_is_empty_long(data.frame(x = NA_real_)))

  heading <- capture.output(.mira_report_heading_long("A\nheading", 99L))
  expect_match(paste(heading, collapse = "\n"), "###### A heading", fixed = TRUE)
  expect_identical(.mira_report_human_name_long("data_overview"), "Dataset overview")
  expect_identical(.mira_report_human_name_long("new.result_name"), "New result name")
  expect_identical(.mira_report_human_name_long("___"), "Unnamed item")

  expect_identical(.mira_report_inline_long(NULL), "not available")
  expect_identical(.mira_report_inline_long(as.Date(c("2025-01-01", "2025-01-02"))),
                   "2025-01-01, 2025-01-02")
  expect_identical(.mira_report_inline_long(c(a = 1, b = NA)), "1, NA")
  expect_match(.mira_report_inline_long(list(alpha = 1)), "alpha")
})

test_that("p-value formatting is exact at the reporting threshold", {
  values <- c(NA, NaN, Inf, -Inf, 0, 0.0001, 0.0009, 0.001,
              0.0011, 0.05, 1, -0.1)
  expect_identical(
    .mira_report_format_p_long(values, digits = 3L),
    c(NA_character_, NA_character_, "Inf", "-Inf", "<0.001", "<0.001",
      "<0.001", "0.001", "0.001", "0.050", "1.000", "-0.100")
  )

  for (digits in c(1L, 2L, 4L, 10L)) {
    threshold <- 10^(-digits)
    formatted <- .mira_report_format_p_long(
      c(threshold / 2, threshold, threshold * 1.1), digits = digits
    )
    label <- formatC(threshold, format = "f", digits = digits)
    expect_identical(formatted[[1L]], paste0("<", label))
    expect_identical(formatted[[2L]], label)
    expect_false(startsWith(formatted[[3L]], "<"))
  }

  expect_identical(.mira_report_format_number_long(NULL), "not estimable")
  expect_identical(.mira_report_format_number_long(NA_real_), "not estimable")
  expect_identical(.mira_report_format_number_long(1.23456, 3L), "1.23")
  expect_identical(.mira_report_format_number_long(Inf), "Inf")
  expect_identical(.mira_report_format_number_long(-Inf), "-Inf")
  expect_true(is.na(.mira_report_format_number_long("not numeric")))
})

test_that("prepare_table normalizes core tabular input classes", {
  numeric_frame <- data.frame(estimate = c(1.2345, NA), p.value = c(0.0004, 0.05),
                              check.names = FALSE)
  prepared <- .mira_report_prepare_table_long(numeric_frame, digits = 3L)
  expect_equal(prepared$estimate, numeric_frame$estimate)
  expect_identical(prepared$p.value, c("<0.001", "0.050"))
  expect_true(report_test_atomic_table(prepared))

  mixed <- data.frame(
    group = factor(c("A", "B")),
    day = as.Date(c("2025-01-01", "2025-01-02")),
    moment = as.POSIXct(c("2025-01-01 01:00:00", "2025-01-02 02:00:00"),
                       tz = "UTC"),
    stringsAsFactors = FALSE
  )
  mixed_out <- .mira_report_prepare_table_long(mixed)
  expect_type(mixed_out$group, "character")
  expect_type(mixed_out$day, "character")
  expect_type(mixed_out$moment, "character")
  expect_true(report_test_atomic_table(mixed_out))

  matrix_value <- matrix(
    1:4, nrow = 2L,
    dimnames = list(c("row-a", "row-b"), c("left", "right"))
  )
  class(matrix_value) <- c("anova", "matrix", "array")
  matrix_out <- .mira_report_prepare_table_long(matrix_value)
  expect_identical(matrix_out$.row, c("row-a", "row-b"))
  expect_equal(unname(as.matrix(matrix_out[c("left", "right")])),
               matrix(1:4, nrow = 2L))
  expect_true(report_test_atomic_table(matrix_out))

  tab <- table(group = c("A", "A", "B"), visit = c("one", "two", "one"))
  table_out <- .mira_report_prepare_table_long(tab)
  expect_setequal(names(table_out), c("group", "visit", "Freq"))
  expect_equal(sum(table_out$Freq), 3)
  expect_true(report_test_atomic_table(table_out))
})

test_that("prepare_table handles vectors, lists, and complex data-frame columns", {
  named_vector <- .mira_report_prepare_table_long(c(alpha = 1, beta = 2))
  expect_identical(named_vector$item, c("alpha", "beta"))
  expect_identical(named_vector$value, c("1", "2"))

  unnamed_vector <- .mira_report_prepare_table_long(c(4, 5))
  expect_identical(unnamed_vector$item, c("1", "2"))
  named_list <- .mira_report_prepare_table_long(list(alpha = 1:2, beta = NULL))
  expect_identical(named_list$parameter, c("alpha", "beta"))
  expect_identical(named_list$value, c("1, 2", "not available"))
  unnamed_list <- .mira_report_prepare_table_long(list("x", 2))
  expect_identical(unnamed_list$parameter, c("item_1", "item_2"))

  complex <- data.frame(id = 1:2)
  complex$list_value <- I(list(c(1, 2), list(label = "x")))
  complex$interval <- I(matrix(
    c(1, 2, 3, 4), nrow = 2L,
    dimnames = list(NULL, c("lower", "upper"))
  ))
  complex_out <- .mira_report_prepare_table_long(complex)
  expect_true(report_test_atomic_table(complex_out))
  expect_identical(complex_out$list_value[[1L]], "1, 2")
  expect_identical(complex_out$interval, c("1, 3", "2, 4"))

  raw_matrix_column <- structure(
    list(id = 1:2, ci = matrix(c(1, 2, 3, 4), nrow = 2L)),
    row.names = c(NA_integer_, -2L), class = "data.frame"
  )
  raw_out <- .mira_report_prepare_table_long(raw_matrix_column)
  expect_true(report_test_atomic_table(raw_out))
  expect_identical(raw_out$ci, c("1, 3", "2, 4"))

  zero_rows <- .mira_report_prepare_table_long(data.frame(a = numeric(), b = character()))
  expect_identical(dim(zero_rows), c(0L, 2L))
  expect_true(report_test_atomic_table(zero_rows))
  zero_columns <- .mira_report_prepare_table_long(data.frame(row.names = 1:3))
  expect_identical(dim(zero_columns), c(3L, 0L))
  expect_true(report_test_atomic_table(zero_columns))

  expect_error(
    .mira_report_prepare_table_long(new.env()),
    "cannot be converted to a table",
    fixed = TRUE
  )
  expect_error(
    .mira_report_prepare_table_long(structure(list(alpha = 1), class = "custom")),
    "cannot be converted to a table",
    fixed = TRUE
  )
})

test_that("p-column detection does not reformat unrelated numeric columns", {
  input <- data.frame(
    topvalue = c(0, 1), p = c(0, 0.05), p_value = c(0.001, 0.2),
    adjusted.p.value = c(0.0001, 0.3), check.names = FALSE
  )
  output <- .mira_report_prepare_table_long(input, digits = 3L)
  expect_identical(output$topvalue, input$topvalue)
  expect_identical(output$p, c("<0.001", "0.050"))
  expect_identical(output$p_value, c("0.001", "0.200"))
  expect_identical(output$adjusted.p.value, c("<0.001", "0.300"))
})

test_that("numeric-column recognition uses values and complete name tokens", {
  recognized <- c("n", "count", "percent", "pct", "mean", "median", "sd",
                  "se", "ci", "estimate", "statistic", "df", "p", "p.value",
                  "p_value", "mean_change", "ci_lower")
  for (name in recognized) {
    expect_true(
      .mira_report_numeric_column_long(c("not", "numeric"), name),
      info = name
    )
  }

  false_positives <- c("sex", "city", "country", "meaning", "topvalue")
  for (name in false_positives) {
    expect_false(
      .mira_report_numeric_column_long(c("low", "high"), name),
      info = name
    )
  }

  expect_true(.mira_report_numeric_column_long(c(rep("1", 9), "x"), "label"))
  expect_false(.mira_report_numeric_column_long(c(rep("1", 8), "x", "y"), "label"))
  expect_true(.mira_report_numeric_column_long(c("<0.001", "0.2", "--"), "label"))
  expect_false(.mira_report_numeric_column_long(c(NA_character_, "", "--"), "label"))
  expect_true(.mira_report_numeric_column_long(factor(c("1", "2")), "label"))
})

test_that("column widths and partitions preserve every non-key column in order", {
  expect_identical(.mira_report_column_widths_long(data.frame()), numeric())
  expect_identical(.mira_report_partition_columns_long(data.frame(), 4L),
                   list(integer()))

  small <- data.frame(id = 1:3, label = c("short", "a much longer value", "mid"),
                      estimate = c(1.1, 2.2, 3.3), stringsAsFactors = FALSE)
  widths <- .mira_report_column_widths_long(small)
  expect_length(widths, ncol(small))
  expect_true(all(is.finite(widths) & widths >= 4))
  expect_lte(widths[[1L]], 14)
  expect_lte(widths[[2L]], 34)
  expect_identical(.mira_report_partition_columns_long(small, 4L), list(1:3))

  large <- data.frame(
    id = seq_len(250L), key = paste0("K", seq_len(250L)),
    a = 1, b = 2, c = 3, d = 4, e = 5,
    stringsAsFactors = FALSE
  )
  blocks <- .mira_report_partition_columns_long(large, max_columns = 4L,
                                            width_budget = 22)
  expect_gt(length(blocks), 1L)
  expect_true(all(vapply(blocks, function(block) all(1:2 %in% block), logical(1L))))
  non_key <- unlist(lapply(blocks, function(block) setdiff(block, 1:2)),
                    use.names = FALSE)
  expect_identical(non_key, 3:ncol(large))
  expect_false(anyDuplicated(non_key) > 0L)
  expect_true(all(vapply(blocks, function(block) identical(block, sort(block)),
                         logical(1L))))

  unlimited <- .mira_report_partition_columns_long(large, max_columns = ncol(large),
                                                width_budget = Inf)
  expect_identical(unlimited, list(seq_len(ncol(large))))
})

test_that("LaTeX helpers create valid alignment and repeat headers idempotently", {
  data <- data.frame(label = c("A", "B"), estimate = c(1.2, 3.4),
                     p.value = c("<0.001", "0.20"), check.names = FALSE)
  alignment <- .mira_report_latex_alignment_long(data)
  expect_length(alignment, ncol(data))
  expect_match(alignment[[1L]], "raggedright", fixed = TRUE)
  expect_match(alignment[[2L]], "raggedleft", fixed = TRUE)
  expect_match(alignment[[3L]], "raggedleft", fixed = TRUE)
  expect_true(all(grepl("p\\{[0-9.]+\\\\linewidth\\}", alignment)))

  escaped <- "alpha\\_beta and price\\$value"
  broken <- .mira_report_latex_breaks_long(escaped)
  expect_match(broken, "\\_\\allowbreak{}", fixed = TRUE)
  expect_match(broken, "\\$\\allowbreak{}", fixed = TRUE)
  expect_identical(gsub("\\\\allowbreak\\{\\}", "", broken), escaped)
  expect_identical(.mira_report_latex_breaks_long("plain text"), "plain text")

  longtable <- paste(
    "\\begin{longtable}{lr}", "\\toprule", "term & estimate\\\\",
    "\\midrule", "a & 1\\\\", "\\bottomrule", "\\end{longtable}",
    sep = "\n"
  )
  repeated <- .mira_report_repeat_longtable_header_long(longtable)
  expect_match(repeated, "\\endfirsthead", fixed = TRUE)
  expect_match(repeated, "\\endhead", fixed = TRUE)
  expect_identical(.mira_report_repeat_longtable_header_long(repeated), repeated)
  expect_identical(.mira_report_repeat_longtable_header_long("minimal"), "minimal")
  no_midrule <- "\\toprule\nheader"
  expect_identical(.mira_report_repeat_longtable_header_long(no_midrule), no_midrule)
  empty_header <- "\\toprule\n\\midrule"
  expect_identical(.mira_report_repeat_longtable_header_long(empty_header), empty_header)
})

test_that("table rendering supports pipe output, truncation, and fallback", {
  cfg <- report_test_cfg(max_table_rows = 2L, max_table_columns = 4L)
  data <- data.frame(id = 1:4, estimate = c(1.2345, 2, 3, 4),
                     p.value = c(0.0001, 0.01, 0.2, 0.5), check.names = FALSE)
  old <- getOption("knitr.kable.NA")
  output <- capture.output(value <- .mira_report_table_long(data, "Results", cfg))
  expect_identical(nrow(value), 2L)
  expect_match(paste(output, collapse = "\n"), "Showing 2 of 4 rows", fixed = TRUE)
  expect_match(paste(output, collapse = "\n"), "Results", fixed = TRUE)
  expect_identical(getOption("knitr.kable.NA"), old)

  fractional <- capture.output(
    fractional_value <- .mira_report_table_long(
      data, cfg = report_test_cfg(max_table_rows = 0.5)
    )
  )
  expect_identical(nrow(fractional_value), 1L)
  expect_match(paste(fractional, collapse = "\n"), "Showing 1 of 4 rows",
               fixed = TRUE)
  expect_no_warning(capture.output(
    huge_value <- .mira_report_table_long(
      data, cfg = report_test_cfg(max_table_rows = .Machine$integer.max + 1)
    )
  ))
  expect_identical(nrow(huge_value), nrow(data))

  empty <- capture.output(value <- .mira_report_table_long(data[0, ], cfg = cfg))
  expect_null(value)
  expect_match(paste(empty, collapse = "\n"), "No rows are available", fixed = TRUE)

  fallback <- capture.output(
    returned <- testthat::with_mocked_bindings(
      .mira_report_table_long(data[1, ], cfg = cfg),
      is_latex_output = function() FALSE,
      is_html_output = function() FALSE,
      kable = function(...) stop("controlled kable failure"),
      .package = "knitr"
    )
  )
  expect_s3_class(returned, "data.frame")
  expect_identical(nrow(returned), 1L)
  expect_match(paste(fallback, collapse = "\n"), "Table formatting fallback",
               fixed = TRUE)
  expect_match(paste(fallback, collapse = "\n"), "controlled kable failure",
               fixed = TRUE)
  expect_match(paste(fallback, collapse = "\n"), "```text", fixed = TRUE)
})

test_that("HTML and LaTeX table branches add their format-specific safeguards", {
  cfg <- report_test_cfg(max_table_columns = 8L)
  data <- data.frame(
    "(term)" = c("(Intercept)", "a_b$c"), estimate = c(1, 2),
    check.names = FALSE, stringsAsFactors = FALSE
  )

  html <- testthat::with_mocked_bindings(
    capture.output(.mira_report_table_long(data, "Model", cfg)),
    is_latex_output = function() FALSE,
    is_html_output = function() TRUE,
    .package = "knitr"
  )
  html <- paste(html, collapse = "\n")
  expect_match(html, "<div class=\"mira-table-wrap\">", fixed = TRUE)
  expect_match(html, "<table class=\"mira-table\">", fixed = TRUE)

  latex <- testthat::with_mocked_bindings(
    capture.output(.mira_report_table_long(data, "Model", cfg)),
    is_latex_output = function() TRUE,
    is_html_output = function() FALSE,
    .package = "knitr"
  )
  latex <- paste(latex, collapse = "\n")
  expect_match(latex, "\\begin{longtable}", fixed = TRUE)
  expect_match(latex, "\\endfirsthead", fixed = TRUE)
  expect_match(latex, "\\toprule\\relax", fixed = TRUE)
  expect_match(latex, "\\_\\allowbreak{}", fixed = TRUE)
  expect_match(latex, "\\$\\allowbreak{}", fixed = TRUE)
  expect_match(latex, "\\begingroup", fixed = TRUE)
  expect_match(latex, "\\footnotesize", fixed = TRUE)
})

test_that("CSV normalization preserves raw scientific values and atomic columns", {
  frame <- data.frame(id = 1:2, p.value = c(0.0001, 0.05), check.names = FALSE)
  frame$interval <- I(matrix(
    c(1, 2, 3, 4), nrow = 2L,
    dimnames = list(NULL, c("lower", "upper"))
  ))
  frame$details <- I(list(c("a", "b"), list(code = 2)))
  csv <- .mira_report_csv_data_long(frame)
  expect_true(report_test_atomic_table(csv))
  expect_type(csv$p.value, "double")
  expect_equal(csv$p.value, frame$p.value)
  expect_equal(csv$interval_lower, c(1, 2))
  expect_equal(csv$interval_upper, c(3, 4))
  expect_identical(csv$details[[1L]], "a, b")
  expect_match(csv$details[[2L]], "code: num 2", fixed = TRUE)

  matrix_value <- matrix(
    1:4, nrow = 2L,
    dimnames = list(c("r1", "r2"), c("same", "same"))
  )
  matrix_csv <- .mira_report_csv_data_long(matrix_value)
  expect_identical(matrix_csv$row, c("r1", "r2"))
  expect_identical(names(matrix_csv), c("row", "same", "same_1"))
  expect_true(report_test_atomic_table(matrix_csv))

  tab <- table(group = c("A", "A", "B"))
  expect_equal(sum(.mira_report_csv_data_long(tab)$Freq), 3)
  zero <- .mira_report_csv_data_long(data.frame(a = numeric(), p = numeric()))
  expect_identical(dim(zero), c(0L, 2L))
  expect_true(report_test_atomic_table(zero))
  expect_error(.mira_report_csv_data_long(new.env()))
})

test_that("filename and quoting helpers are portable and reversible", {
  expect_identical(.mira_report_file_slug_long("My analysis"), "my_analysis")
  expect_identical(.mira_report_file_slug_long("a/b:c"), "a_b_c")
  expect_identical(.mira_report_file_slug_long("résumé"), "resume")
  expect_identical(.mira_report_file_slug_long("___"), "item")
  expect_identical(.mira_report_file_slug_long(""), "item")
  expect_identical(.mira_report_slug_long("___"), "mira_report")

  reserved <- c("CON", "PRN", "AUX", "NUL", "COM1", "LPT9")
  expect_identical(
    unname(vapply(reserved, .mira_report_file_slug_long, character(1L))),
    paste0("item_", tolower(reserved))
  )
  expect_identical(
    unname(vapply(reserved, .mira_report_slug_long, character(1L))),
    paste0("mira_", tolower(reserved))
  )
  expect_lte(nchar(.mira_report_file_slug_long(strrep("a", 200L))), 90L)
  expect_lte(nchar(.mira_report_slug_long(strrep("a", 200L))), 90L)

  greek <- .mira_report_file_slug_long("α β γ")
  expect_true(nzchar(greek))
  expect_match(greek, "^[a-z0-9_]+$")

  expect_identical(.mira_report_yaml_quote_long("O'Brien\nanalysis"),
                   "'O''Brien analysis'")
  r_value <- "C:\\analysis\\a\"b"
  expect_identical(eval(parse(text = .mira_report_r_quote_long(r_value))), r_value)

  format_lines <- .mira_report_format_lines_long(c("html", "pdf", "docx"))
  expect_identical(sum(format_lines == "format:"), 1L)
  expect_identical(sum(format_lines == "  html:"), 1L)
  expect_identical(sum(format_lines == "  pdf:"), 1L)
  expect_identical(sum(format_lines == "  docx:"), 1L)
  expect_true(all(c("    theme: cosmo", "    documentclass: article",
                    "    papersize: a4", "    geometry: margin=25mm") %in%
                  format_lines))
})

test_that("plot and model recognition covers base and optional display classes", {
  fit <- stats::lm(mpg ~ wt, data = mtcars)
  test <- stats::t.test(1:5, 2:6)
  expect_true(.mira_report_is_model_long(fit))
  expect_true(.mira_report_is_model_long(test))
  expect_true(.mira_report_is_model_long(structure(
    list(), class = c("lmerModLmerTest", "lmerMod", "merMod")
  )))
  expect_false(.mira_report_is_model_long(data.frame(x = 1)))
  expect_true(.mira_report_is_plain_list_long(list(x = 1)))
  expect_false(.mira_report_is_plain_list_long(data.frame(x = 1)))
  expect_false(.mira_report_is_plain_list_long(fit))

  expect_true(.mira_report_is_plot_long(structure(list(), class = "recordedplot")))
  expect_true(.mira_report_is_plot_long(grid::rectGrob()))
  expect_false(.mira_report_is_plot_long(list(x = 1)))

  skip_if_not_installed("ggplot2")
  plot <- ggplot2::ggplot(mtcars, ggplot2::aes(wt, mpg)) + ggplot2::geom_point()
  expect_true(.mira_report_is_plot_long(plot))
  device <- withr::local_tempfile(fileext = ".pdf")
  grDevices::pdf(device)
  on.exit(grDevices::dev.off(), add = TRUE)
  output <- capture.output(printed <- .mira_report_plot_long(plot, "mean_ci"))
  expect_true(printed)
  expect_match(paste(output, collapse = "\n"), "confidence intervals",
               ignore.case = TRUE)

  broken <- structure(list(), class = "recordedplot")
  output <- capture.output(
    printed <- testthat::with_mocked_bindings(
      .mira_report_plot_long(broken, "broken"),
      replayPlot = function(...) stop("controlled plot failure"),
      .package = "grDevices"
    )
  )
  expect_false(printed)
  expect_match(paste(output, collapse = "\n"), "could not be rendered",
               ignore.case = TRUE)
})
