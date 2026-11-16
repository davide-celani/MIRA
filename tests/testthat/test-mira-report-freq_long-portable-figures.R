test_that("a relocated real report knits its saved plots and scientific results", {
  skip_if_not_installed("knitr")
  # Drawing saved plots can consume RNG state (for example through jitter).
  withr::local_preserve_seed()
  root <- withr::local_tempdir(pattern = "mira-portable-figures-")
  original <- file.path(root, "original")
  relocated <- file.path(root, "relocated bundle")
  y <- rbind(c(3, 5, 8), c(7, 6, 10), c(4, 8, 7), c(9, 12, 11),
    c(6, 9, 12), c(12, 11, 15), c(8, 13, 16), c(15, 14, 18))
  data <- data.frame(patient = paste0("S", 1:8),
    treatment = rep(c("Control", "Treatment"), each = 4L),
    score_t0 = y[, 1], score_t1 = y[, 2], score_t2 = y[, 3])
  fit <- mira_info_long(data, id = "patient", outcomes = "score",
    time_vars = c("score_t0", "score_t1", "score_t2"),
    arm = "treatment", reference_arm = "Control", covariates = character(),
    analyses = c("correlations", "plots", "arm_tests"), verbose = FALSE)
  report <- report_test_call(list(x = fit, output_file = "scientific",
    digits = 3L, max_table_columns = 100L, date = "2026-10-12"), original)

  # Visit sums are 64, 78 and 97; final-minus-baseline differences sum to 33
  # with squared differences summing to 165, so their SE is sqrt(33) / 8.
  expect_identical(fit$descriptives$time, paste0("score_t", 0:2))
  expect_equal(fit$descriptives$score_mean, c(8, 9.75, 12.125), tolerance = 1e-12)
  final <- fit$change[fit$change$from == "score_t0" &
    fit$change$to == "score_t2", , drop = FALSE]
  expect_equal(nrow(final), 1L)
  expect_identical(final$n, 8L)
  expect_equal(final$score_mean_change, 33 / 8, tolerance = 1e-12)
  expected_ci <- 33 / 8 + c(-1, 1) * stats::qt(0.975, 7) * sqrt(33) / 8
  expect_equal(c(final$score_ci_lower, final$score_ci_upper), expected_ci,
    tolerance = 1e-12)
  expect_equal(final$score_paired_t_p, 2 * stats::pt(-sqrt(33), 7),
    tolerance = 1e-12)

  manifest <- utils::read.csv(report$results_manifest, stringsAsFactors = FALSE)
  expect_true(all(file.exists(file.path(original, manifest$file))))
  expect_false(any(grepl("^([A-Za-z]:|/|\\\\\\\\)", manifest$file)))
  entry <- which(manifest$path == "outcome$change")
  expect_length(entry, 1L)
  archived <- utils::read.csv(file.path(original, manifest$file[[entry]]),
    check.names = FALSE)
  expect_identical(names(archived), names(fit$change))
  expect_equal(nrow(archived), 3L)
  expect_true(is.numeric(archived$score_paired_t_p))
  expect_equal(archived$score_paired_t_p, fit$change$score_paired_t_p,
    tolerance = 1e-12)

  dir.create(relocated)
  expect_true(all(file.copy(list.files(original, full.names = TRUE,
    all.files = TRUE, no.. = TRUE), relocated, recursive = TRUE)))
  # A hidden absolute reference must not be rescued by the original bundle.
  original_path <- normalizePath(original, winslash = "/", mustWork = TRUE)
  root_path <- normalizePath(root, winslash = "/", mustWork = TRUE)
  expect_identical(dirname(original_path), root_path)
  expect_identical(unlink(original_path, recursive = TRUE), 0L)
  expect_false(dir.exists(original))
  isolated <- new.env(parent = baseenv())
  old_options <- knitr::opts_knit$get()
  old_chunks <- knitr::opts_chunk$get()
  withr::defer(knitr::opts_knit$restore(old_options))
  withr::defer(knitr::opts_chunk$restore(old_chunks))
  grDevices::pdf(file.path(root, "knit-device.pdf"))
  device <- grDevices::dev.cur()
  withr::defer(grDevices::dev.off(device))
  knitr::opts_knit$set(global.device = FALSE, rmarkdown.pandoc.to = "html")
  withr::local_dir(relocated)
  expect_no_error(knitr::knit("scientific.qmd", output = "scientific.md",
    envir = isolated, quiet = TRUE))
  markdown <- paste(readLines("scientific.md", warn = FALSE, encoding = "UTF-8"),
    collapse = "\n")
  expect_match(markdown, "Descriptive statistics", fixed = TRUE)
  expect_match(markdown, "Within-participant changes", fixed = TRUE)
  expect_match(markdown, "12.125", fixed = TRUE)
  expect_match(markdown, "4.125", fixed = TRUE)
  expect_false(grepl("Error in|could not find function|cannot open", markdown))
  expect_false(grepl(original, markdown, fixed = TRUE))
  expect_identical(isolated$payload$result$descriptives, fit$descriptives)
  expect_identical(isolated$payload$result$change, fit$change)
  figures <- list.files(relocated, pattern = "\\.(png|pdf|svg)$", recursive = TRUE,
    full.names = TRUE)
  expect_gt(length(figures), 0L)
  expect_true(all(file.info(figures)$size > 0))
})
