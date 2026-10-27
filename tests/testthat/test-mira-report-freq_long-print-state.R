report_print_state_gee <- function() {
  y <- rbind(c(2, 5, 7), c(4, 6, 12), c(7, 8, 9), c(3, 7, 13),
    c(9, 11, 10), c(6, 9, 16), c(8, 14, 18), c(12, 13, 15))
  data <- data.frame(patient_factor = factor(rep(1:8, each = 3)),
    time_factor = factor(rep(c("Zbase", "Amid", "Mfinal"), 8),
      levels = c("Zbase", "Amid", "Mfinal")), value = as.vector(t(y)))
  geepack::geeglm(value ~ time_factor, id = patient_factor, data = data,
    family = stats::gaussian(), corstr = "exchangeable")
}

test_that("native GEE report printing preserves caller digits and later scientific precision", {
  skip_if_not_installed("geepack")
  withr::local_options(digits = 8L)
  model <- report_print_state_gee()
  before <- serialize(model, NULL)
  for (object in list(model, summary(model))) {
    text <- paste(capture.output(.mira_report_text_block_long(object)), collapse = "\n")
    expect_match(text, "Coefficients", fixed = TRUE)
    expect_false(grepl("The object could not be printed", text, fixed = TRUE))
    expect_identical(getOption("digits"), 8L)
  }
  expect_identical(serialize(model, NULL), before)
  text <- testthat::with_mocked_bindings(
    paste(capture.output(.mira_report_table_long(
      data.frame(time = c("Zbase", "Amid"), estimate = c(1.234567, 2.345678)),
      caption = "Following scientific table", cfg = report_test_cfg(digits = 4L))),
      collapse = "\n"),
    is_html_output = function(...) TRUE, is_latex_output = function(...) FALSE,
    .package = "knitr")
  table <- report_test_html_table(text, "Following scientific table")
  expect_identical(table$time, c("Zbase", "Amid"))
  expect_identical(table$estimate, c("1.2346", "2.3457"))
  expect_identical(getOption("digits"), 8L)
})

test_that("dumped runtime isolates native GEE printing from later formatted tables", {
  skip_if_not_installed("geepack")
  withr::local_options(digits = 8L)
  model <- report_print_state_gee()
  outcome <- report_test_mira()
  outcome$advanced_models$gee <- list(exchangeable = list(performed = TRUE,
    model = model, summary = summary(model)))
  report <- report_test_call(list(x = outcome, include_complete_output = FALSE),
    withr::local_tempdir(pattern = "mira-native-print-state-"))
  runtime <- new.env(parent = baseenv())
  sys.source(report$runtime, envir = runtime)
  payload <- readRDS(report$payload)
  text <- testthat::with_mocked_bindings(
    paste(capture.output(runtime$.mira_report_knit_long(payload)), collapse = "\n"),
    is_html_output = function(...) TRUE, is_latex_output = function(...) FALSE,
    .package = "knitr")
  expect_match(text, "Coefficients", fixed = TRUE)
  expect_false(grepl("The object could not be printed", text, fixed = TRUE))
  expect_identical(getOption("digits"), 8L)
  formatted <- testthat::with_mocked_bindings(
    paste(capture.output(runtime$.mira_report_table_long(
      data.frame(time = c("A < B", "B & C"), mean = c(1.234567, 2.345678)),
      caption = "Independent means", cfg = report_test_cfg())), collapse = "\n"),
    is_html_output = function(...) TRUE, is_latex_output = function(...) FALSE,
    .package = "knitr")
  table <- report_test_html_table(formatted, "Independent means")
  expect_identical(table$time, c("A < B", "B & C"))
  expect_identical(table$mean, c("1.235", "2.346"))
  expect_identical(getOption("digits"), 8L)
})

test_that("failed native printers restore digits in package and saved runtime helpers", {
  skip_if_not_installed("geepack")
  withr::local_options(digits = 8L)
  # geepack's printer changes global digits before it encounters the malformed
  # coefficient table. This exercises the actual external error path.
  malformed <- structure(list(), class = "summary.geeglm")
  report <- report_test_call(output_dir =
    withr::local_tempdir(pattern = "mira-failed-native-print-state-"))
  runtime <- new.env(parent = baseenv())
  sys.source(report$runtime, envir = runtime)
  for (printer in list(.mira_report_text_block_long,
                       runtime$.mira_report_text_block_long)) {
    text <- paste(capture.output(printer(malformed)), collapse = "\n")
    expect_match(text, "The object could not be printed:", fixed = TRUE)
    expect_identical(getOption("digits"), 8L)
  }
})
