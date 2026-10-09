test_that("a moved source-only bundle knits using only its saved runtime and payload", {
  skip_if_not_installed("knitr")
  root <- withr::local_tempdir(pattern = "mira-portable-bundle-")
  original <- file.path(root, "original")
  relocated <- file.path(root, "relocated bundle")
  input <- report_test_mira()
  report <- report_test_call(list(
    x = input, sections = c("descriptives", "change"),
    title = "Portable scientific report", date = "2025-01-15"
  ), original)
  dir.create(relocated)
  companions <- c(report$source, report$runtime, report$payload)
  expect_true(all(file.copy(companions, relocated)))

  # No attached package or test helper can rescue a missing dumped dependency.
  isolated <- new.env(parent = baseenv())
  old_options <- knitr::opts_knit$get()
  old_chunks <- knitr::opts_chunk$get()
  withr::defer(knitr::opts_knit$restore(old_options))
  withr::defer(knitr::opts_chunk$restore(old_chunks))
  # Use a fixture-owned device: some knitr versions unlink their temporary PNG
  # before closing its device, which then recreates an orphaned scratch file.
  grDevices::pdf(file.path(root, "knit-device.pdf"))
  device_id <- grDevices::dev.cur()
  withr::defer(grDevices::dev.off(device_id))
  knitr::opts_knit$set(global.device = TRUE)
  withr::local_dir(relocated)
  markdown_file <- file.path(relocated, "deferred.md")
  expect_no_error(knitr::knit(
    input = basename(report$source), output = markdown_file,
    envir = isolated, quiet = TRUE
  ))
  markdown <- paste(readLines(markdown_file, warn = FALSE), collapse = "\n")
  expect_match(markdown, "Portable scientific report", fixed = TRUE)
  expect_match(markdown, "Outcome: Outcome score", fixed = TRUE)
  expect_match(markdown, "Descriptive statistics", fixed = TRUE)
  expect_match(markdown, "(&lt;|<)0\\.001")
  expect_false(grepl("Error in|could not find function|cannot open", markdown))
  expect_false(grepl(report$output_dir, markdown, fixed = TRUE))
  expect_identical(isolated$payload$result, input)
  expect_true(exists(".mira_report_knit_long", isolated, inherits = FALSE))
})

test_that("display settings and selected sections preserve every archived scientific column", {
  root <- withr::local_tempdir(pattern = "mira-complete-archive-")
  input <- report_test_mira()
  input$wide_results <- data.frame(subject = c("A_1", "B$2", "C[3]"))
  for (i in seq_len(13L)) {
    input$wide_results[[paste0("estimate_", i)]] <- c(pi, NA_real_, -pi) * i
  }
  input$wide_results$p.value <- c(0.000012345, NA_real_, 0.123456789)
  before <- input
  reports <- lapply(seq_len(2L), function(i) {
    report_test_call(list(
      x = input, sections = if (i == 1L) "overview" else "all",
      include_complete_output = i == 2L, digits = if (i == 1L) 1L else 10L,
      max_table_rows = if (i == 1L) 1L else Inf,
      max_table_columns = if (i == 1L) 2L else 12L
    ), file.path(root, paste0("case-", i)))
  })
  manifests <- lapply(reports, function(report) {
    utils::read.csv(report$results_manifest, stringsAsFactors = FALSE)
  })
  expect_identical(manifests[[1L]], manifests[[2L]])
  for (i in seq_along(reports)) {
    manifest <- manifests[[i]]
    row <- match("outcome$wide_results", manifest$path)
    expect_false(is.na(row))
    exported <- utils::read.csv(
      file.path(reports[[i]]$output_dir, manifest$file[[row]]),
      stringsAsFactors = FALSE, check.names = FALSE
    )
    # CSV text round-tripping can round the last few binary digits.
    expect_equal(exported, input$wide_results, tolerance = 1e-12)
    expect_identical(dim(exported), c(3L, 15L))
    expect_identical(readRDS(reports[[i]]$payload)$result, before)
  }
  expect_identical(input, before)
})

test_that("matrix-valued CSV columns retain raw values and disambiguate expanded names", {
  root <- withr::local_tempdir(pattern = "mira-matrix-archive-")
  input <- report_test_mira()
  tab <- data.frame(id = c("subject_1", "subject_2"), ci_lower = c(99, 100))
  tab$ci <- I(matrix(c(1.25, NA_real_, 2.75, 4.5), nrow = 2L,
                     dimnames = list(NULL, c("lower", "upper"))))
  tab$labels <- I(list(c("A", "B"), NULL))
  tab$p.value <- c(0.000000123, NA_real_)
  input$intervals <- tab
  report <- report_test_call(list(x = input, digits = 1L), root)
  manifest <- utils::read.csv(report$results_manifest, stringsAsFactors = FALSE)
  row <- match("outcome$intervals", manifest$path)
  archived <- utils::read.csv(
    file.path(root, manifest$file[[row]]), stringsAsFactors = FALSE,
    check.names = FALSE, na.strings = ""
  )
  expect_identical(names(archived),
                   c("id", "ci_lower", "ci_lower_1", "ci_upper", "labels", "p.value"))
  expect_equal(archived$ci_lower, c(99, 100))
  expect_equal(archived$ci_lower_1, c(1.25, NA_real_))
  expect_equal(archived$ci_upper, c(2.75, 4.5))
  expect_identical(archived$labels, c("A, B", "not available"))
  expect_equal(archived$p.value, tab$p.value, tolerance = 1e-12)
  expect_identical(readRDS(report$payload)$result$intervals, tab)
})

test_that("max_depth limits archive traversal while preserving the complete payload", {
  root <- withr::local_tempdir(pattern = "mira-depth-archive-")
  input <- report_test_mira()
  input$nested <- list(level_one = list(result = data.frame(estimate = pi)))
  paths <- vector("list", 2L)
  for (depth in 1:2) {
    report <- report_test_call(list(x = input, max_depth = depth),
                               file.path(root, paste0("depth-", depth)))
    paths[[depth]] <- utils::read.csv(report$results_manifest)$path
    expect_identical(readRDS(report$payload)$result$nested, input$nested)
  }
  target <- "outcome$nested$level_one$result"
  expect_false(target %in% paths[[1L]])
  expect_true(target %in% paths[[2L]])
  expect_true("outcome$descriptives" %in% paths[[1L]])
})

test_that("each protected artifact independently blocks overwrite before any bundle writes", {
  root <- withr::local_tempdir(pattern = "mira-artifact-protection-")
  probe <- report_test_call(output_dir = file.path(root, "probe"))
  protected <- c(
    basename(probe$source), basename(probe$payload), basename(probe$runtime),
    basename(probe$expected_files), "results_manifest.csv",
    report_test_relpath(sub(paste0(probe$output_dir, "/"), "",
                           probe$result_files[[1L]], fixed = TRUE))
  )
  for (i in seq_along(protected)) {
    directory <- file.path(root, paste0("protected-", i))
    sentinel <- file.path(directory, protected[[i]])
    dir.create(dirname(sentinel), recursive = TRUE)
    writeLines("preserve this file", sentinel)
    before <- list.files(directory, recursive = TRUE, all.files = TRUE)
    expect_error(
      report_test_call(output_dir = directory),
      "Output files already exist. Use overwrite=TRUE", fixed = TRUE,
      info = protected[[i]]
    )
    expect_identical(readLines(sentinel), "preserve this file", info = protected[[i]])
    expect_identical(list.files(directory, recursive = TRUE, all.files = TRUE), before)
  }
})

test_that("different stems in one directory still protect the shared analytical archive", {
  root <- withr::local_tempdir(pattern = "mira-shared-archive-")
  first <- report_test_call(list(output_file = "first"), root)
  paths <- c(first$source, first$payload, first$runtime,
             first$results_manifest, first$result_files)
  checksums <- tools::md5sum(paths)
  expect_error(
    report_test_call(list(output_file = "second"), root),
    "Output files already exist.*results_manifest\\.csv"
  )
  expect_identical(tools::md5sum(paths), checksums)
  expect_false(file.exists(file.path(root, "second.qmd")))
})

test_that("default output placement is relative to the caller's working directory", {
  root <- withr::local_tempdir(pattern = "mira-default-placement-")
  withr::local_dir(root)
  report <- report_test_call(output_dir = NULL)
  expected <- file.path(root, "mira_analyses", "mira_report_flong")
  expect_identical(report$output_dir,
                   normalizePath(expected, winslash = "/", mustWork = TRUE))
  expect_true(all(file.exists(c(report$source, report$runtime, report$payload))))
})

test_that("an unreadable old manifest does not authorize deletion of unlisted CSV files", {
  root <- withr::local_tempdir(pattern = "mira-malformed-manifest-")
  dir.create(file.path(root, "results"))
  manual <- file.path(root, "results", "manual.csv")
  writeLines("irreplaceable observations", manual)
  writeLines('file\n"unterminated', file.path(root, "results_manifest.csv"))
  # A malformed CSV can warn before read.csv returns or errors; the deletion
  # invariant applies in either case and does not depend on R's parser wording.
  report <- suppressWarnings(report_test_call(list(overwrite = TRUE), root))
  expect_identical(readLines(manual), "irreplaceable observations")
  expect_true(all(file.exists(report$result_files)))
  manifest <- utils::read.csv(report$results_manifest, stringsAsFactors = FALSE)
  expect_false("results/manual.csv" %in% manifest$file)
})
