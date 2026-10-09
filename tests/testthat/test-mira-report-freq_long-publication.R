# Scientific document tests control format detection but retain real kable
# formatting. They require neither Quarto nor a TeX installation.
report_publication_table_capture <- function(x, cfg = report_test_cfg(),
                                             target = "html", caption = NULL,
                                             fail_block = NULL) {
  skip_if_not_installed("knitr")
  original_kable <- knitr::kable
  calls <- list()
  lines <- testthat::with_mocked_bindings(
    capture.output(returned <- withVisible(
      .mira_report_table_long(x, caption = caption, cfg = cfg)
    )),
    is_latex_output = function() identical(target, "latex"),
    is_html_output = function() identical(target, "html"),
    kable = function(...) {
      args <- list(...)
      calls[[length(calls) + 1L]] <<- args
      if (identical(length(calls), fail_block)) {
        stop("controlled block formatting failure", call. = FALSE)
      }
      do.call(original_kable, args)
    },
    .package = "knitr"
  )
  list(lines = lines, text = paste(lines, collapse = "\n"), calls = calls,
       value = returned$value, visible = returned$visible)
}

test_that("LaTeX wrapping preserves escaped delimiters at the end of identifiers", {
  # strsplit() drops trailing empty pieces: reconstructing identifiers from its
  # pieces used to remove the last escaped underscore or dollar sign.
  identifiers <- c("arm\\_visit\\_", "cost\\$unit\\$", "\\_\\_", "\\$\\$",
                   "arm\\_\\$", "plain\\_", "plain\\$", "plain text")
  for (identifier in identifiers) {
    wrapped <- .mira_report_latex_breaks_long(identifier)
    expect_identical(
      gsub("\\allowbreak{}", "", wrapped, fixed = TRUE), identifier,
      info = identifier
    )
  }
  expect_match(.mira_report_latex_breaks_long("arm\\_visit\\_"),
               "visit\\_", fixed = TRUE)
  expect_match(.mira_report_latex_breaks_long("cost\\$unit\\$"),
               "unit\\$", fixed = TRUE)
})

test_that("multidimensional table cells retain row identity and missing values", {
  array_column <- array(seq_len(12L), dim = c(3L, 2L, 2L))
  array_column[3L, 2L, 2L] <- NA_integer_
  input <- data.frame(
    id = c("P3", "P1", "P2"),
    group = factor(c("B", NA, "A"), levels = c("A", "B", "unused")),
    date = as.Date(c("2025-01-03", NA, "2025-01-01")),
    moment = as.POSIXct(c("2025-01-03 12:00:00", NA, "2025-01-01 12:00:00"),
                       tz = "UTC"),
    stringsAsFactors = FALSE
  )
  input$intervals <- I(array_column)
  input$details <- I(list(NULL, c(NA_real_, 2), c("first", "last")))
  before <- unserialize(serialize(input, NULL))
  prepared <- .mira_report_prepare_table_long(input)
  expect_true(report_test_atomic_table(prepared))
  expect_identical(dim(prepared), c(3L, 6L))
  expect_identical(prepared$id, input$id)
  expect_identical(prepared$group, c("B", NA_character_, "A"))
  expect_identical(prepared$date, c("2025-01-03", NA_character_, "2025-01-01"))
  expect_identical(prepared$moment, as.character(input$moment))
  expect_identical(prepared$intervals, c("1, 4, 7, 10", "2, 5, 8, 11", "3, 6, 9, NA"))
  expect_identical(prepared$details, c("not available", "NA, 2", "first, last"))
  expect_identical(input, before)
})

test_that("publication p columns include model headers and preserve preformatted labels", {
  input <- data.frame(
    "Pr(>|t|)" = c(0.00004, NA_real_, 0.05),
    "P.VALUE" = c(0, 0.001, 1),
    "paired_t_p_adj" = c(0.0004, 0.2, NaN),
    "p_value_label" = c("<0.0001", "not estimable", "0.0500"),
    "probability" = c(0.00004, NA_real_, 0.05),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  before <- input
  prepared <- .mira_report_prepare_table_long(input, digits = 4L)
  expect_identical(prepared[["Pr(>|t|)"]], c("<0.0001", NA_character_, "0.0500"))
  expect_identical(prepared[["P.VALUE"]], c("<0.0001", "0.0010", "1.0000"))
  expect_identical(prepared$paired_t_p_adj, c("0.0004", "0.2000", NA_character_))
  expect_identical(prepared$p_value_label, input$p_value_label)
  expect_identical(prepared$probability, input$probability)
  expect_identical(.mira_report_prepare_table_long(input, format_p = FALSE), input)
  expect_identical(input, before)
})

test_that("rendered table blocks preserve every statistic and repeat participant keys", {
  input <- data.frame(participant_id = c("P2", "P1", "P3"),
                      time_label = c("Baseline", "Final", "Middle"))
  for (j in seq_len(8L)) input[[paste0("estimate_", j)]] <- j + seq_len(3L) / 10
  rendered <- report_publication_table_capture(
    input, report_test_cfg(max_table_rows = 2L, max_table_columns = 4L),
    caption = "Observed estimates"
  )
  expect_false(rendered$visible)
  expect_identical(rendered$value, input[1:2, , drop = FALSE])
  expect_length(rendered$calls, 4L)
  statistic_columns <- unlist(lapply(rendered$calls, function(call) {
    expect_identical(call$x$participant_id, input$participant_id[1:2])
    expect_identical(call$x$time_label, input$time_label[1:2])
    expect_lte(ncol(call$x), 4L)
    expect_identical(call$col.names, gsub("[._]+", " ", names(call$x)))
    expect_identical(call$escape, TRUE)
    expect_identical(call$row.names, FALSE)
    names(call$x)[-(1:2)]
  }), use.names = FALSE)
  expect_identical(statistic_columns, names(input)[-(1:2)])
  expect_identical(vapply(rendered$calls, `[[`, character(1L), "caption"),
                   sprintf("Observed estimates (block %d of 4)", 1:4))
  expect_identical(sum(grepl("Showing 2 of 3 rows", rendered$lines, fixed = TRUE)), 1L)
  expect_identical(sum(grepl("<div class=\"mira-table-wrap\">", rendered$lines,
                             fixed = TRUE)), 4L)
  expect_identical(sum(rendered$lines == "</div>"), 4L)
})

test_that("width-aware PDF partitions retain long identifiers and all columns", {
  input <- data.frame(
    participant_identifier = rep(strrep("subject_", 8L), 2L),
    analysis_object_path = rep(strrep("model$", 8L), 2L),
    stringsAsFactors = FALSE
  )
  for (j in seq_len(10L)) input[[paste0("statistic_", j)]] <- c(j, j + 0.5)
  blocks <- .mira_report_partition_columns_long(input, max_columns = 8L,
                                               width_budget = 76)
  # Two long identifier columns consume too much page width; the second must
  # still appear once as content when only the first remains a repeated key.
  expect_gt(length(blocks), 1L)
  expect_true(all(vapply(blocks, function(block) block[[1L]] == 1L, logical(1L))))
  expect_identical(unlist(lapply(blocks, function(block) block[-1L]), use.names = FALSE),
                   seq.int(2L, ncol(input)))
  expect_true(all(vapply(blocks, length, integer(1L)) <= 8L))
  rendered <- report_publication_table_capture(input, target = "latex",
                                                caption = "Model diagnostics")
  expect_true(all(vapply(rendered$calls, function(call) {
    isTRUE(call$longtable) && isTRUE(call$booktabs) && ncol(call$x) <= 8L
  }, logical(1L))))
  expect_match(rendered$text, "subject\\_\\allowbreak{}", fixed = TRUE)
  expect_match(rendered$text, "model\\$\\allowbreak{}", fixed = TRUE)
  expect_identical(sum(grepl("\\endfirsthead", rendered$lines, fixed = TRUE)),
                   length(rendered$calls))
  expect_identical(sum(grepl("\\endhead", rendered$lines, fixed = TRUE)),
                   length(rendered$calls))
})

test_that("HTML tables escape scientific labels and apply numeric precision", {
  input <- data.frame(
    "label<&>" = c("A < B & C > D", "quoted \"name\""),
    estimate = c(1.234567, NA_real_), p.value = c(0.00001, NA_real_),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  rendered <- report_publication_table_capture(input, report_test_cfg(digits = 4L),
                                                caption = "Statistical results")
  expect_match(rendered$text, "label&lt;&amp;&gt;", fixed = TRUE)
  expect_match(rendered$text, "A &lt; B &amp; C &gt; D", fixed = TRUE)
  expect_false(grepl("A < B & C > D", rendered$text, fixed = TRUE))
  expect_match(rendered$text, "1.2346", fixed = TRUE)
  expect_match(rendered$text, "&lt;0.0001", fixed = TRUE)
  expect_match(rendered$text, "--", fixed = TRUE)
  expect_match(rendered$text, "Statistical results", fixed = TRUE)
  expect_false(grepl("not estimable", rendered$text, fixed = TRUE))
  expect_false(grepl("NULL", rendered$text, fixed = TRUE))
})

test_that("PDF table guards cover bracketed first cells and repeated headers", {
  for (label in c("(Intercept)", "[reference]")) {
    header <- if (startsWith(label, "(")) "(term)" else "[term]"
    input <- data.frame(term = c(label, "arm_visit_$"), estimate = c(1, NA_real_))
    names(input)[1L] <- header
    rendered <- report_publication_table_capture(input, target = "latex",
                                                  caption = "Coefficients")
    # knitr versions may prefix brackets with '{}'; the report also has a
    # relax guard. Neither header may remain a bare optional booktabs argument.
    expect_false(grepl("\\\\toprule\\s*[([]", rendered$text, perl = TRUE))
    expect_false(grepl("\\\\midrule\\s*[([]", rendered$text, perl = TRUE))
    expect_match(rendered$text, "\\_\\allowbreak{}", fixed = TRUE)
    expect_match(rendered$text, "\\$\\allowbreak{}", fixed = TRUE)
    expect_match(rendered$text, "--", fixed = TRUE)
    expect_match(rendered$text, "\\end{longtable}", fixed = TRUE)
    expect_match(rendered$text, "\\endgroup", fixed = TRUE)
    expect_identical(sum(grepl(header, rendered$lines, fixed = TRUE)), 2L)
  }
})

test_that("longtable continuation headers retain every original header line", {
  original <- paste(
    "\\begin{longtable}{lr}", "  \\toprule",
    "\\multicolumn{2}{c}{Repeated observations}\\\\",
    "Visit & Estimate\\\\", "  \\midrule",
    "Baseline & 2\\\\", "\\bottomrule", "\\end{longtable}", sep = "\n"
  )
  continued <- .mira_report_repeat_longtable_header_long(original)
  lines <- strsplit(continued, "\n", fixed = TRUE)[[1L]]
  expect_identical(sum(lines == "Visit & Estimate\\\\"), 2L)
  expect_identical(sum(lines == "\\multicolumn{2}{c}{Repeated observations}\\\\"), 2L)
  expect_identical(sum(lines == "Baseline & 2\\\\"), 1L)
  expect_lt(match("\\endfirsthead", lines), match("\\endhead", lines))
  expect_lt(match("\\endhead", lines), match("Baseline & 2\\\\", lines))
  expect_identical(.mira_report_repeat_longtable_header_long(continued), continued)
})

test_that("one unformattable table block does not suppress later scientific results", {
  input <- data.frame(participant = c("P1", "P2"), visit = c("First", "Last"),
                      early_estimate = c(1, 2), early_sd = c(3, 4),
                      later_estimate = c(5, 6), later_sd = c(7, 8))
  withr::local_options(knitr.kable.NA = "pre-existing missing label")
  rendered <- report_publication_table_capture(
    input, report_test_cfg(max_table_columns = 4L),
    caption = "Arm estimates", fail_block = 1L
  )
  expect_identical(getOption("knitr.kable.NA"), "pre-existing missing label")
  expect_length(rendered$calls, 2L)
  expect_identical(rendered$value, input)
  expect_false(rendered$visible)
  expect_match(rendered$text, "Table formatting fallback", fixed = TRUE)
  expect_match(rendered$text, "controlled block formatting failure", fixed = TRUE)
  expect_match(rendered$text, "early_estimate", fixed = TRUE)
  expect_match(rendered$text, "```text", fixed = TRUE)
  expect_match(rendered$text, "Arm estimates (block 2 of 2)", fixed = TRUE)
  expect_match(rendered$text, "later estimate", fixed = TRUE)
  expect_match(rendered$text, "<table", fixed = TRUE)
})

test_that("report text fallback handles broken printers and embedded Markdown fences", {
  unreadable <- testthat::with_mocked_bindings(
    capture.output(.mira_report_text_block_long(structure(1, class = "unprintable"))),
    print = function(...) stop("controlled printing failure", call. = FALSE),
    .package = "base"
  )
  expect_match(paste(unreadable, collapse = "\n"),
               "The object could not be printed: controlled printing failure", fixed = TRUE)
  expect_identical(sum(unreadable == "```text"), 1L)
  expect_identical(sum(unreadable == "```"), 1L)
  fenced <- testthat::with_mocked_bindings(
    capture.output(.mira_report_text_block_long("source notes")),
    print = function(...) cat("stored text with ```embedded fence```\n"),
    .package = "base"
  )
  expect_match(paste(fenced, collapse = "\n"), "'''embedded fence'''", fixed = TRUE)
  expect_identical(sum(fenced == "```text"), 1L)
  expect_identical(sum(fenced == "```"), 1L)
})

test_that("rendered statistical tables retain explanatory attributes including empty results", {
  skip_if_not_installed("knitr")
  result <- data.frame(term = "Time", estimate = 1.25)
  attr(result, "note") <- "Confidence intervals use the prespecified confidence level."
  attr(result, "interpretation") <- "A positive estimate denotes an increase."
  node <- .mira_report_prune_long(result, "fixed_effects", "outcome$model$fixed_effects",
                                  report_test_cfg())
  rendered <- capture.output(.mira_report_render_node_long(node, 3L, report_test_cfg()))
  text <- paste(rendered, collapse = "\n")
  expect_match(text, "### Fixed effects", fixed = TRUE)
  expect_match(text, attr(result, "note"), fixed = TRUE)
  expect_match(text, attr(result, "interpretation"), fixed = TRUE)
  empty <- result[0L, , drop = FALSE]
  attr(empty, "note") <- "No eligible contrasts were observed."
  empty_node <- .mira_report_prune_long(empty, "contrasts", "outcome$advanced_tests$contrasts",
                                        report_test_cfg(), parent_performed = TRUE)
  empty_output <- capture.output(.mira_report_render_node_long(empty_node, 3L, report_test_cfg()))
  expect_match(paste(empty_output, collapse = "\n"), "Observed results:** 0 rows", fixed = TRUE)
  expect_match(paste(empty_output, collapse = "\n"), attr(empty, "note"), fixed = TRUE)
  expect_false(any(grepl("No rows are available for display", empty_output, fixed = TRUE)))
})

test_that("assembled report sections follow scientific order and omit unselected content", {
  skip_if_not_installed("knitr")
  input <- report_test_mira()
  input$diagnostics <- list(warnings = "Prespecified adaptation was used.")
  input$plot_error <- "controlled figure failure"
  payload <- list(result = input, extra_objects = list(audit = data.frame(result = 42)),
                  report = report_test_cfg(c("appendix", "diagnostics", "change",
                                              "descriptives", "overview", "methods")))
  before <- unserialize(serialize(payload, NULL))
  output <- capture.output(visible <- withVisible(.mira_report_knit_long(payload)))
  headings <- output[grepl("^#{1,2} ", output)]
  expect_identical(headings, c(
    "# Analysis context", "## Resolved configuration", "# Outcome: Outcome score",
    "## Analysis settings", "## Overview", "## Descriptive statistics",
    "## Within-participant changes", "## Diagnostics and adaptations",
    "# Supplementary objects", "## Audit"
  ))
  text <- paste(output, collapse = "\n")
  expect_match(text, "Prespecified adaptation was used.", fixed = TRUE)
  expect_match(text, "controlled figure failure", fixed = TRUE)
  # The methods legitimately describe enabled analyses. Selection omits their
  # result sections, while retaining the prespecified analysis configuration.
  expect_false("## Correlations across time points" %in% headings)
  expect_false("## Missing data" %in% headings)
  expect_false(any(grepl("^# Abstract$|^# Descriptive summary$", output)))
  expect_false(grepl("NULL", text, fixed = TRUE))
  expect_false(visible$visible)
  expect_null(visible$value)
  expect_identical(payload, before)
})

test_that("descriptive conclusions compare endpoint means only for eligible outcomes", {
  skip_if_not_installed("knitr")
  outcomes <- list(
    improving = list(descriptives = data.frame(time = c("First", "Middle", "Last"),
                                                mean = c(10, 99, 12.5))),
    declining = list(descriptives = data.frame(mean = c(25, 19))),
    one_visit = list(descriptives = data.frame(mean = 1)),
    missing_first = list(descriptives = data.frame(mean = c(NA_real_, 2))),
    infinite_last = list(descriptives = data.frame(mean = c(2, Inf))),
    missing_column = list(descriptives = data.frame(n = c(2L, 2L)))
  )
  observed <- NULL
  output <- testthat::with_mocked_bindings(
    capture.output(.mira_report_render_conclusions_long(outcomes, report_test_cfg("conclusions"))),
    .mira_report_table_long = function(x, caption = NULL, cfg) {
      observed <<- list(table = x, caption = caption)
      invisible(x)
    },
    .package = "MIRA"
  )
  expect_identical(output[grepl("^# ", output)], "# Descriptive summary")
  expect_identical(observed$table$outcome, c("improving", "declining"))
  expect_equal(observed$table$mean_first_visit, c(10, 25), tolerance = 1e-12)
  expect_equal(observed$table$mean_last_visit, c(12.5, 19), tolerance = 1e-12)
  expect_equal(observed$table$descriptive_difference, c(2.5, -6), tolerance = 1e-12)
  expect_identical(observed$caption, "Differences between observed means")
  expect_false(any(grepl("significant|causal|improved", output, ignore.case = TRUE)))
  expect_length(capture.output(.mira_report_render_conclusions_long(
    outcomes[-(1:2)], report_test_cfg("conclusions")
  )), 0L)
  expect_length(capture.output(.mira_report_render_conclusions_long(
    outcomes, report_test_cfg("overview")
  )), 0L)
})

test_that("supplementary models, statistical tests and nested objects render semantically", {
  skip_if_not_installed("knitr")
  data <- data.frame(x = 1:6, y = c(1, 2, 4, 3, 6, 8))
  fit <- stats::lm(y ~ x, data = data)
  test <- stats::t.test(data$y)
  payload <- list(
    result = report_test_detect(),
    extra_objects = list(
      fitted_model = fit, one_sample_test = test, registered_formula = y ~ x,
      nested = list(table = data.frame(term = "Time", estimate = 2),
                    missing = NULL, deeper = list(value = "registered result")),
      missing = NULL, unsupported_function = identity,
      unsupported_environment = new.env(parent = emptyenv())
    ),
    report = report_test_cfg("appendix")
  )
  output <- capture.output(.mira_report_knit_long(payload))
  text <- paste(output, collapse = "\n")
  expect_identical(output[grepl("^# ", output)], "# Supplementary objects")
  expect_match(text, "Coefficients:", fixed = TRUE)
  expect_match(text, "One Sample t-test", fixed = TRUE)
  expect_match(text, "y ~ x", fixed = TRUE)
  expect_match(text, "registered result", fixed = TRUE)
  expect_false(grepl("$coefficients", text, fixed = TRUE))
  expect_false(grepl("unsupported", text, fixed = TRUE))
  expect_false(grepl("NULL", text, fixed = TRUE))
  payload$report <- report_test_cfg("overview")
  expect_length(capture.output(.mira_report_knit_long(payload)), 0L)
})

test_that("figure captions choose the most specific scientific description", {
  device_path <- withr::local_tempfile(fileext = ".pdf")
  grDevices::pdf(device_path)
  device_id <- grDevices::dev.cur()
  withr::defer(grDevices::dev.off(device_id))
  figure <- grid::rectGrob()
  cases <- c(
    arm_mean_ci = "Longitudinal mean profile by arm with confidence intervals.",
    arm_change_ci = "Mean change by arm with confidence intervals.",
    arm_difference_ci = "Estimated arm differences with confidence intervals.",
    custom_diagnostic = "Plot `custom_diagnostic` produced by MIRA."
  )
  for (name in names(cases)) {
    output <- capture.output(visible <- withVisible(.mira_report_plot_long(figure, name)))
    expect_false(visible$visible)
    expect_true(visible$value)
    expect_match(paste(output, collapse = "\n"), paste0("**Figure.** ", cases[[name]]),
                 fixed = TRUE)
    expect_false(any(grepl("callout-warning", output, fixed = TRUE)))
  }
  failure <- testthat::with_mocked_bindings(
    capture.output(visible <- withVisible(.mira_report_plot_long(figure, "arm_mean_ci"))),
    grid.draw = function(...) stop("controlled drawing failure", call. = FALSE),
    .package = "grid"
  )
  expect_false(visible$visible)
  expect_false(visible$value)
  expect_match(paste(failure, collapse = "\n"), "Figure could not be rendered: arm_mean_ci",
               fixed = TRUE)
  expect_match(paste(failure, collapse = "\n"), "controlled drawing failure", fixed = TRUE)
  expect_false(any(grepl("**Figure.**", failure, fixed = TRUE)))
})
