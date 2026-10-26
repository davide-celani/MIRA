# Real analysis -> saved payload/runtime/QMD/CSV -> relocated, isolated knitr.
# This exercises available rendering through knitr, NOT the Quarto CLI.
# Rscript --vanilla dev/validate-bundle.R <existing-output-directory>
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L || !dir.exists(args[[1L]])) stop("Supply one existing output directory.")
output <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
library(MIRA)
validate <- function() {
  y <- rbind(c(3, 5, 8), c(7, 6, 10), c(4, 8, 7), c(9, 12, 11),
             c(6, 9, 12), c(12, 11, 15), c(8, 13, 16), c(15, 14, 18))
  data <- data.frame(patient = paste0("S", 1:8),
    treatment = rep(c("Control", "Treatment"), each = 4L),
    score_t0 = y[, 1], score_t1 = y[, 2], score_t2 = y[, 3])
  fit <- mira_info_long(data, id = "patient", outcomes = "score",
    time_vars = c("score_t0", "score_t1", "score_t2"),
    arm = "treatment", reference_arm = "Control", covariates = character(),
    analyses = c("correlations", "plots", "arm_tests"), verbose = FALSE)
  source_dir <- file.path(output, "source")
  if (file.exists(source_dir)) stop("Use a fresh output directory.")
  printed <- capture.output(report <- mira_report_freq_long(fit,
    output_dir = source_dir, output_file = "scientific", format = "html",
    render = FALSE, open = FALSE, digits = 3L, max_table_columns = 100L,
    title = "Real longitudinal scientific validation", date = "2026-10-12"))
  saved <- readRDS(report$payload)$result
  stopifnot(identical(saved$descriptives, fit$descriptives),
            identical(saved$change, fit$change))
  # Hand totals: visit sums=64,78,97; final-baseline differences sum=33,
  # sum of squared differences=165 -> variance=33/8, SE=sqrt(33)/8.
  # The explicit outcome declaration prefixes these metrics. Check shape and
  # finiteness first so absent columns cannot satisfy stopifnot(logical(0)).
  stopifnot(identical(fit$descriptives$time, c("score_t0", "score_t1", "score_t2")),
    length(fit$descriptives$score_mean) == 3L,
    all(is.finite(fit$descriptives$score_mean)),
    max(abs(fit$descriptives$score_mean - c(8, 9.75, 12.125))) < 1e-12,
    nrow(fit$change) == 3L, !anyDuplicated(fit$change[c("from", "to")]))
  final <- fit$change[fit$change$from == "score_t0" & fit$change$to == "score_t2", ]
  expected_ci <- 33 / 8 + c(-1, 1) * stats::qt(0.975, 7) * sqrt(33) / 8
  observed <- c(final$score_mean_change, final$score_ci_lower,
    final$score_ci_upper, final$score_paired_t_p)
  stopifnot(nrow(final) == 1L, final$n == 8L, length(observed) == 4L,
    all(is.finite(observed)), abs(final$score_mean_change - 33 / 8) < 1e-12,
    max(abs(c(final$score_ci_lower, final$score_ci_upper) - expected_ci)) < 1e-12,
    abs(final$score_paired_t_p - 2 * stats::pt(-sqrt(33), 7)) < 1e-12)
  manifest <- utils::read.csv(report$results_manifest, stringsAsFactors = FALSE)
  stopifnot(all(file.exists(file.path(report$output_dir, manifest$file))),
            !any(grepl("^([A-Za-z]:|/|\\\\)", manifest$file)))
  change_file <- manifest$file[manifest$path == "outcome$change"]
  stopifnot(length(change_file) == 1L)
  csv <- utils::read.csv(file.path(source_dir, change_file), check.names = FALSE)
  stopifnot(is.numeric(csv$score_paired_t_p),
            length(csv$score_paired_t_p) == 3L,
            max(abs(csv$score_paired_t_p - fit$change$score_paired_t_p)) < 1e-12,
            identical(names(csv), names(fit$change)), nrow(csv) == 3L)
  relocated <- file.path(output, "relocated bundle")
  dir.create(relocated)
  stopifnot(all(file.copy(list.files(source_dir, full.names = TRUE, all.files = TRUE, no.. = TRUE),
                         relocated, recursive = TRUE)))
  isolated <- new.env(parent = baseenv())
  options_before <- knitr::opts_knit$get()
  chunks_before <- knitr::opts_chunk$get()
  withr::defer(knitr::opts_knit$restore(options_before))
  withr::defer(knitr::opts_chunk$restore(chunks_before))
  grDevices::pdf(file.path(output, "knit-device.pdf"))
  device <- grDevices::dev.cur()
  withr::defer(grDevices::dev.off(device))
  knitr::opts_knit$set(global.device = FALSE, rmarkdown.pandoc.to = "html")
  withr::local_dir(relocated)
  md <- knitr::knit("scientific.qmd", output = "scientific.md", envir = isolated, quiet = TRUE)
  text <- paste(readLines(md, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  stopifnot(grepl("Descriptive statistics", text, fixed = TRUE),
    grepl("Within-participant changes", text, fixed = TRUE),
    grepl("12.125", text, fixed = TRUE), grepl("4.125", text, fixed = TRUE),
    !grepl("Error in|could not find function|cannot open", text),
    !grepl(source_dir, text, fixed = TRUE),
    identical(isolated$payload$result$change, fit$change))
  figures <- list.files(relocated, pattern = "\\.(png|pdf|svg)$", recursive = TRUE, full.names = TRUE)
  stopifnot(length(figures) > 0L, all(file.info(figures)$size > 0))
  writeLines(capture.output(sessionInfo()), file.path(output, "bundle-session.txt"))
  cat("PASS: real analysis, hand-derived final contrast, complete raw CSV, relative manifest, relocated isolated runtime and plots.\n")
  cat("Figures:", length(figures), "CSV files:", nrow(manifest), "\n")
  cat("Quarto HTML/PDF/DOCX was NOT executed by this script.\n")
}
validate()
