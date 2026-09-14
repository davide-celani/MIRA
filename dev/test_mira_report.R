# Reproducible examples for mira_info() and mira_report_freq().
# Run this file from the directory in which you want mira_example_checks/.
# Set render_reports <- TRUE to generate HTML files when Quarto is installed.

render_reports <- FALSE
output_root <- file.path(getwd(), "mira_example_checks")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

check <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

has_result <- function(manifest, pattern) {
  any(grepl(pattern, manifest$path, perl = TRUE))
}

run_example <- function(name, fit, expected = character(), absent = character()) {
  output_dir <- file.path(output_root, name)
  result <- mira_report_freq(
    x = fit,
    output_dir = output_dir,
    format = "html",
    render = render_reports,
    open = FALSE,
    overwrite = TRUE
  )
  manifest <- utils::read.csv(result$results_manifest, stringsAsFactors = FALSE)
  csv_files <- list.files(file.path(output_dir, "results"), pattern = "\\.csv$",
                          recursive = TRUE, full.names = TRUE)
  check(nrow(manifest) == length(csv_files),
        paste(name, "manifest and CSV file counts differ"))
  check(all(file.exists(file.path(output_dir, manifest$file))),
        paste(name, "manifest contains a missing CSV file"))
  check(file.exists(result$source) && file.exists(result$payload),
        paste(name, "report source or data payload is missing"))
  for (i in seq_len(nrow(manifest))) {
    table <- utils::read.csv(file.path(output_dir, manifest$file[[i]]),
                             check.names = FALSE)
    check(nrow(table) == manifest$rows[[i]] &&
            ncol(table) == manifest$columns[[i]],
          paste(name, "CSV dimensions differ from the manifest:",
                manifest$path[[i]]))
  }
  for (pattern in expected) {
    check(has_result(manifest, pattern),
          paste(name, "is missing the expected result:", pattern))
  }
  for (pattern in absent) {
    check(!has_result(manifest, pattern),
          paste(name, "contains an unexpected result:", pattern))
  }
  if (render_reports) {
    check(isTRUE(result$rendered[["html"]]) &&
            file.exists(result$expected_files[["html"]]),
          paste(name, "HTML rendering did not complete"))
  }
  cat(sprintf("PASS %-24s %3d CSV tables  %s\n",
              name, nrow(manifest), output_dir))
  invisible(list(fit = fit, report = result, manifest = manifest))
}

# 1. Core analysis only: no optional arm, correlation, model, outlier, or plot module.
set.seed(101)
n <- 30L
baseline <- rnorm(n, 50, 7)
core_data <- data.frame(
  subject_id = sprintf("C%03d", seq_len(n)),
  score_t0 = baseline,
  score_t1 = baseline - 2 + rnorm(n, sd = 2),
  score_t2 = baseline - 4 + rnorm(n, sd = 2)
)
core_fit <- mira_info(core_data, id = "subject_id", outcomes = "score",
                      time_vars = c("score_t0", "score_t1", "score_t2"),
                      analyses = "none", verbose = FALSE)
core_run <- run_example("01_core_only", core_fit,
                        expected = c("outcome\\$descriptives$", "outcome\\$missing\\$by_time$",
                                     "outcome\\$change$"),
                        absent = c("outcome\\$correlations\\$", "outcome\\$model\\$",
                                   "outcome\\$arm_analysis\\$", "outcome\\$outliers\\$"))

# 2. Correlations with genuine missing values.
set.seed(102)
n <- 36L
baseline <- rnorm(n, 45, 6)
correlation_data <- data.frame(
  subject_id = sprintf("R%03d", seq_len(n)),
  score_t0 = baseline,
  score_t1 = baseline + rnorm(n, sd = 2),
  score_t2 = baseline + rnorm(n, sd = 3),
  score_t3 = baseline + rnorm(n, sd = 4)
)
correlation_data$score_t1[c(3, 9, 16)] <- NA_real_
correlation_data$score_t2[c(5, 11, 25)] <- NA_real_
correlation_fit <- mira_info(correlation_data, id = "subject_id",
                             outcomes = "score", time_vars = paste0("score_t", 0:3),
                             analyses = "correlations", verbose = FALSE)
correlation_run <- run_example("02_correlations", correlation_fit,
                               expected = c("outcome\\$correlations\\$pearson$",
                                            "outcome\\$correlations\\$spearman$",
                                            "outcome\\$correlations\\$pairwise_n$"),
                               absent = c("outcome\\$arm_analysis\\$", "outcome\\$outliers\\$"))

# 3. Two arms and a categorical gender covariate, with all available analyses.
# This case exercises the classed ANOVA and covariance matrices that previously
# caused a CSV export error in R 4.4.
set.seed(202601)
n <- 48L
arm <- factor(rep(c("Control", "Treatment"), each = n / 2),
              levels = c("Control", "Treatment"))
gender <- factor(rep(c("Female", "Male"), length.out = n))
baseline <- rnorm(n, mean = 50, sd = 7) + ifelse(gender == "Male", 2, 0)
arm_data <- data.frame(
  subject_id = sprintf("S%03d", seq_len(n)),
  arm = arm,
  gender = gender,
  score_t0 = baseline,
  score_t1 = baseline - 1 - ifelse(arm == "Treatment", 2, 0) +
    rnorm(n, sd = 1.8),
  score_t2 = baseline - 2 - ifelse(arm == "Treatment", 5, 0) +
    rnorm(n, sd = 2.0)
)
arm_fit <- mira_info(arm_data, id = "subject_id", outcomes = "score",
                     time_vars = paste0("score_t", 0:2), arm = "arm",
                     reference_arm = "Control", covariates = "gender",
                     analyses = "all", p_adjust_method = "holm", verbose = FALSE)
arm_run <- run_example("03_arms_gender_all", arm_fit,
                       expected = c("outcome\\$arm_analysis\\$descriptives$"))

# 4. Two outcomes with different availability. Each outcome is processed
# independently, and its CSV files are stored in a separate subdirectory.
set.seed(104)
n <- 32L
score_base <- rnorm(n, 55, 6)
fatigue_base <- rnorm(n, 20, 4)
multi_data <- data.frame(
  subject_id = sprintf("M%03d", seq_len(n)),
  score_t0 = score_base,
  score_t1 = score_base - 2 + rnorm(n),
  score_t2 = score_base - 4 + rnorm(n),
  fatigue_t0 = fatigue_base,
  fatigue_t1 = fatigue_base + 1 + rnorm(n),
  fatigue_t2 = fatigue_base + 2 + rnorm(n)
)
multi_data$fatigue_t1[1:5] <- NA_real_
multi_fit <- mira_info(multi_data, id = "subject_id",
                       outcomes = c("score", "fatigue"), analyses = "correlations",
                       verbose = FALSE)
multi_run <- run_example("04_two_outcomes", multi_fit,
                         expected = c("outcome\\$descriptives$"))
check(length(unique(multi_run$manifest$outcome[multi_run$manifest$outcome != ""])) == 2L,
      "The multi-outcome manifest does not contain both outcomes")
check(all(grepl("^results/[^/]+/", multi_run$manifest$file[
  multi_run$manifest$outcome != ""])),
  "Multi-outcome CSV files are not separated by outcome")

# 5. Sparse measurements: correlations are requested but cannot be estimated.
# The report should keep core results without creating a correlation CSV.
set.seed(105)
n <- 16L
sparse_data <- data.frame(
  subject_id = sprintf("P%03d", seq_len(n)),
  score_t0 = rnorm(n, 40, 5),
  score_t1 = rep(NA_real_, n),
  score_t2 = rep(NA_real_, n)
)
sparse_data$score_t1[1] <- 39
sparse_data$score_t2[2] <- 37
sparse_fit <- mira_info(sparse_data, id = "subject_id", outcomes = "score",
                        time_vars = paste0("score_t", 0:2), analyses = "correlations",
                        verbose = FALSE)
sparse_run <- run_example("05_sparse_data", sparse_fit,
                          expected = "outcome\\$descriptives$",
                          absent = "outcome\\$correlations\\$")

# 6. A new result not known to the report's label map. The adaptive traversal
# should still display and export its table without changing report code.
extended_fit <- core_fit
extended_fit$outcomes[[1L]]$new_analysis <- list(
  results = data.frame(metric = c("A", "B"), estimate = c(0, 1.25))
)
extended_run <- run_example("06_new_component", extended_fit,
                            expected = "outcome\\$new_analysis\\$results$")

cat("All examples passed. Set render_reports <- TRUE for HTML reports.\n")
