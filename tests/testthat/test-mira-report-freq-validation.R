test_that("mira_report_freq is part of the public package API", {
  expect_true("mira_report_freq" %in% getNamespaceExports("MIRA"))
  expect_s3_class(report_test_mira(), "mira_info")
})

test_that("logical controls accept only one non-missing logical value", {
  invalid <- list(NA, NULL, 1, 0, "TRUE", logical(), c(TRUE, FALSE))
  flags <- c("include_complete_output", "render", "open", "overwrite", "quiet")

  for (flag in flags) {
    for (value in invalid) {
      args <- list()
      args[flag] <- list(value)
      expect_error(
        report_test_call(args),
        paste0(flag, " must be TRUE or FALSE"),
        fixed = TRUE,
        info = paste(flag, "with", paste(value, collapse = ","))
      )
    }
  }

  root <- withr::local_tempdir(pattern = "mira-flags-")
  report <- report_test_call(
    list(
      include_complete_output = FALSE,
      render = FALSE,
      open = FALSE,
      overwrite = TRUE,
      quiet = FALSE
    ),
    file.path(root, "valid")
  )
  expect_false(report$report_config$include_complete_output)
})

test_that("output_file is validated and normalized to a portable stem", {
  invalid <- list("", "   ", NA_character_, c("one", "two"), 42, character())
  for (value in invalid) {
    expect_error(
      report_test_call(list(output_file = value)),
      "output_file must be a non-empty string",
      fixed = TRUE
    )
  }

  root <- withr::local_tempdir(pattern = "mira-output-file-")
  values <- c("My analysis", "a/b:c", "résumé", "α β γ", "___")
  for (i in seq_along(values)) {
    report <- report_test_call(
      list(output_file = values[[i]]), file.path(root, paste0("case-", i))
    )
    stem <- tools::file_path_sans_ext(basename(report$source))
    expect_match(stem, "^[a-z0-9_]+$", info = values[[i]])
    expect_true(nzchar(stem), info = values[[i]])
  }

  reserved <- c("CON", "PRN", "AUX", "NUL", "COM1", "LPT9")
  for (i in seq_along(reserved)) {
    report <- report_test_call(
      list(output_file = reserved[[i]]), file.path(root, paste0("reserved-", i))
    )
    stem <- tolower(tools::file_path_sans_ext(basename(report$source)))
    expect_false(stem %in% tolower(reserved), info = reserved[[i]])
  }

  long_name <- paste(rep("publication", 30L), collapse = "-")
  report <- report_test_call(
    list(output_file = long_name), file.path(root, "long-name")
  )
  expect_lte(nchar(tools::file_path_sans_ext(basename(report$source))), 90L)
})

test_that("title, subtitle, author, and date follow their documented contracts", {
  bad_title <- list("", "  ", NA_character_, c("a", "b"), 1, character())
  for (value in bad_title) {
    expect_error(
      report_test_call(list(title = value)),
      "title must be a non-empty string",
      fixed = TRUE
    )
  }

  bad_subtitle <- list(NA_character_, c("a", "b"), 1, character())
  for (value in bad_subtitle) {
    expect_error(
      report_test_call(list(subtitle = value)),
      "subtitle must be NULL or a single string",
      fixed = TRUE
    )
  }

  bad_author <- list("", "  ", c("Ada", NA_character_), 1, list("Ada"))
  for (value in bad_author) {
    expect_error(
      report_test_call(list(author = value)),
      "author must be NULL or a vector of non-empty strings",
      fixed = TRUE
    )
  }

  bad_date <- list(character(), c("2025-01-01", "2025-01-02"),
                   NA_character_, as.Date(NA))
  for (value in bad_date) {
    expect_error(
      report_test_call(list(date = value)),
      "date must have length one and must not be missing",
      fixed = TRUE
    )
  }

  root <- withr::local_tempdir(pattern = "mira-metadata-")
  valid <- list(
    list(title = "Quoted 'title' \"value\"", subtitle = NULL,
         author = NULL, date = as.Date("2025-01-15")),
    list(title = "Unicode résumé", subtitle = "", author = character(),
         date = "2025-01-15"),
    list(title = "Study", subtitle = "Line one", author = "Ada Lovelace",
         date = as.POSIXct("2025-01-15 12:00:00", tz = "UTC")),
    list(title = "Study", subtitle = "Sub", author = c("Ada", "R. Fisher"),
         date = "fixed date")
  )
  for (i in seq_along(valid)) {
    report <- report_test_call(valid[[i]], file.path(root, paste0("valid-", i)))
    expect_true(file.exists(report$source))
  }
})

test_that("output_dir handles existing, nested, and invalid paths safely", {
  root <- withr::local_tempdir(pattern = "mira-output-dir-")
  existing <- file.path(root, "existing")
  dir.create(existing)
  report_existing <- report_test_call(output_dir = existing)
  expect_identical(
    report_existing$output_dir,
    normalizePath(existing, winslash = "/", mustWork = TRUE)
  )

  nested <- file.path(root, "one", "two", "three")
  report_nested <- report_test_call(output_dir = nested)
  expect_true(dir.exists(nested))
  expect_identical(
    report_nested$output_dir,
    normalizePath(nested, winslash = "/", mustWork = TRUE)
  )

  invalid <- list("", "   ", NA_character_, c("a", "b"), 1, character())
  for (value in invalid) {
    expect_error(
      report_test_call(output_dir = value),
      "output_dir must be a non-empty path",
      fixed = TRUE
    )
  }

  blocker <- withr::local_tempfile(lines = "ordinary file")
  expect_error(
    report_test_call(output_dir = file.path(blocker, "child")),
    "Cannot create directory",
    fixed = TRUE
  )
})

test_that("all documented sections validate case-insensitively", {
  allowed <- c(
    "all", "abstract", "methods", "data_quality", "overview",
    "descriptives", "frequencies", "missingness", "change", "arms",
    "correlations", "variability", "models", "robustness", "outliers",
    "trajectories", "figures", "diagnostics", "conclusions", "appendix"
  )
  root <- withr::local_tempdir(pattern = "mira-sections-")
  for (i in seq_along(allowed)) {
    report <- report_test_call(
      list(sections = toupper(allowed[[i]])),
      file.path(root, paste0("section-", i))
    )
    expect_identical(report$report_config$sections, allowed[[i]])
  }

  report <- report_test_call(
    list(sections = c("OVERVIEW", "models", "overview")),
    file.path(root, "combined")
  )
  expect_identical(report$report_config$sections, c("overview", "models"))

  invalid <- list(character(), NA_character_, 1, "unknown",
                  c("overview", "unknown"))
  for (value in invalid) {
    expect_error(
      report_test_call(list(sections = value)),
      "sections|Unrecognized sections"
    )
  }
})

test_that("formats normalize without hiding invalid values", {
  root <- withr::local_tempdir(pattern = "mira-formats-")
  report <- report_test_call(
    list(format = c("HTML", "pdf", "html")), file.path(root, "combined")
  )
  expect_identical(report$formats, c("html", "pdf"))
  expect_identical(names(report$expected_files), report$formats)
  expect_identical(tools::file_ext(unname(report$expected_files)), report$formats)

  all_report <- report_test_call(
    list(format = "ALL"), file.path(root, "all")
  )
  expect_identical(all_report$formats, c("html", "pdf", "docx"))
  expect_identical(names(all_report$expected_files), all_report$formats)

  invalid <- list(character(), NA_character_, 1, "epub", c("html", "epub"),
                  c("all", "epub"))
  for (value in invalid) {
    expect_error(
      report_test_call(list(format = value)),
      "format|Unrecognized formats"
    )
  }
})

test_that("numeric controls validate bounds and record intentional coercion", {
  invalid <- list(
    max_table_rows = list(0, -1, NA_real_, NaN, "5", numeric(), c(1, 2)),
    max_table_columns = list(1, 0, -1, Inf, NA_real_, "5", numeric(), c(2, 3),
                             .Machine$integer.max + 1),
    digits = list(0, 11, Inf, NA_real_, "3", numeric(), c(2, 3)),
    max_depth = list(0, -1, Inf, NA_real_, "3", numeric(), c(2, 3),
                     .Machine$integer.max + 1)
  )
  patterns <- c(
    max_table_rows = "max_table_rows must be a positive number or Inf",
    max_table_columns = "max_table_columns must be a finite integer >= 2",
    digits = "digits must be an integer between 1 and 10",
    max_depth = "max_depth must be a finite integer >= 1"
  )
  for (name in names(invalid)) {
    for (value in invalid[[name]]) {
      args <- list()
      args[name] <- list(value)
      expect_error(report_test_call(args), patterns[[name]], fixed = TRUE)
    }
  }

  root <- withr::local_tempdir(pattern = "mira-numeric-")
  report <- report_test_call(
    list(max_table_rows = 1.5, max_table_columns = 5.9,
         digits = 4.9, max_depth = 7.9),
    file.path(root, "coercion")
  )
  expect_identical(report$report_config$max_table_rows, 1.5)
  expect_identical(report$report_config$max_table_columns, 5L)
  expect_identical(report$report_config$digits, 4L)
  expect_identical(report$report_config$max_depth, 7L)

  infinite <- report_test_call(
    list(max_table_rows = Inf, digits = 10), file.path(root, "infinite")
  )
  expect_identical(infinite$report_config$max_table_rows, Inf)
  expect_identical(infinite$report_config$digits, 10L)
})

test_that("extra_objects are named uniquely and preserved in the payload", {
  expect_error(
    report_test_call(list(extra_objects = 1)),
    "extra_objects must be a list",
    fixed = TRUE
  )

  root <- withr::local_tempdir(pattern = "mira-extra-")
  objects <- list(1:3, data.frame(x = 1), letters[1:2], pi)
  names(objects) <- c("", "named", "named", NA_character_)
  report <- report_test_call(
    list(extra_objects = objects, sections = "appendix"), root
  )
  saved <- readRDS(report$payload)$extra_objects
  expect_identical(names(saved), c("extra_1", "named", "named.1", "extra_4"))
  expect_identical(unname(saved), unname(objects))

  empty <- report_test_call(
    list(extra_objects = list()), file.path(root, "empty")
  )
  expect_identical(readRDS(empty$payload)$extra_objects, list())
})

test_that("input dispatch accepts documented classes and rejects conflicts", {
  root <- withr::local_tempdir(pattern = "mira-inputs-")
  accepted <- list(
    info = report_test_mira(),
    multi = report_test_multi(),
    detect = report_test_detect()
  )
  for (name in names(accepted)) {
    report <- report_test_call(
      list(x = accepted[[name]]), file.path(root, name)
    )
    expect_identical(report$mira_result, accepted[[name]])
  }

  expect_error(
    mira_report_freq(),
    "Supply x (the output of mira_info) or data",
    fixed = TRUE
  )
  expect_error(
    report_test_call(list(x = structure(list(), class = "mira_info_error"))),
    "x must inherit from mira_info, mira_info_multi, or mira_detect",
    fixed = TRUE
  )
  expect_error(
    report_test_call(list(x = list())),
    "x must inherit from mira_info, mira_info_multi, or mira_detect",
    fixed = TRUE
  )
  expect_error(
    report_test_call(list(data = 1, x = NULL)),
    "data must be a data frame",
    fixed = TRUE
  )

  data <- data.frame(patient = 1:2, score_t0 = 1:2, score_t1 = 2:3)
  expect_error(
    report_test_call(list(x = data, data = data)),
    "Supply data only once, through x or data",
    fixed = TRUE
  )
  expect_error(
    report_test_call(list(x = report_test_mira(), data = data)),
    "Supply a MIRA object in x or a data frame in data, but not both",
    fixed = TRUE
  )
  expect_error(
    report_test_call(list(x = report_test_mira(), unsupported_dot = TRUE)),
    "Arguments in ... are allowed only when data is supplied",
    fixed = TRUE
  )

  duplicate <- list(
    data, data, file.path(root, "duplicate"), FALSE, FALSE
  )
  names(duplicate) <- c("data", "data", "output_dir", "render", "open")
  expect_error(
    do.call(mira_report_freq, duplicate),
    "matched by multiple actual arguments|formal argument"
  )
})

test_that("data dispatch forwards dots and supplies verbose FALSE by default", {
  root <- withr::local_tempdir(pattern = "mira-forwarding-")
  data <- data.frame(patient = 1:3, score_t0 = 1:3, score_t1 = 2:4)
  seen <- new.env(parent = emptyenv())
  fake_mira_info <- function(data, ...) {
    seen$data <- data
    seen$dots <- list(...)
    report_test_mira()
  }
  testthat::local_mocked_bindings(mira_info = fake_mira_info, .package = "MIRA")

  report <- report_test_call(
    list(x = NULL, data = data, forwarded_token = "kept"),
    file.path(root, "data")
  )
  expect_identical(seen$data, data)
  expect_identical(seen$dots$forwarded_token, "kept")
  expect_identical(seen$dots$verbose, FALSE)
  expect_s3_class(report$mira_result, "mira_info")

  report_test_call(
    list(x = data, verbose = TRUE), file.path(root, "x-data-frame")
  )
  expect_identical(seen$data, data)
  expect_identical(seen$dots$verbose, TRUE)
})

test_that("the missing mira_info dependency fails without changing global bindings", {
  expect_true(exists("mira_info", mode = "function", inherits = TRUE))
  isolated <- new.env(parent = baseenv())
  isolated$.mira_report_validate_scalar_flag <- .mira_report_validate_scalar_flag
  isolated_function <- mira_report_freq
  environment(isolated_function) <- isolated

  expect_error(
    isolated_function(
      data = data.frame(x = 1), render = FALSE, open = FALSE
    ),
    "mira_info() is not available",
    fixed = TRUE
  )
  expect_true(exists("mira_info", mode = "function", inherits = TRUE))
})
