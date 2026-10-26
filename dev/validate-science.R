# Extended offline Monte Carlo study; run from the repository root:
# Rscript --vanilla dev/validate-science.R <existing-output-directory>
# All required R packages must already be installed. No downloads/installers.
# This validates paired mean-change inference under its assumptions only.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L || !dir.exists(args[[1L]])) {
  stop("Supply one existing output directory.")
}
output <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
library(MIRA)

study <- function() {
  withr::local_preserve_seed()
  withr::local_seed(2026101201)
  # Preregistered before execution: n=24, B=500 in each of delta=0,2;
  # independent subjects with baseline ~ N(10, 3^2) and change ~ N(delta,2^2),
  # independent of baseline, no missingness. Target E[Y1-Y0]=delta.
  # Bias MCSE = empirical SD(estimates)/sqrt(B); binomial MCSE for coverage
  # and rejection. Interpret deviations >4 nominal MCSE as an investigation
  # signal (not a fast-suite pass threshold). Normal-theory t intervals are
  # exact here. No conclusion extends to informative missingness or engines.
  n <- 24L
  B <- 500L
  alpha <- 0.05
  rows <- vector("list", 2L * B)
  index <- 0L
  for (target in c(0, 2)) {
    for (replicate in seq_len(B)) {
      baseline <- stats::rnorm(n, 10, 3)
      change <- stats::rnorm(n, target, 2)
      data <- data.frame(patient = seq_len(n), score_t0 = baseline,
                         score_t1 = baseline + change)
      fit <- mira_info_long(data, id = "patient", outcomes = "score",
        time_vars = c("score_t0", "score_t1"), covariates = character(),
        analyses = "none", alpha = alpha, verbose = FALSE)
      stopifnot(nrow(fit$change) == 1L,
        identical(fit$change$from, "score_t0"),
        identical(fit$change$to, "score_t1"), fit$change$n == n)
      row <- fit$change[1L, ]
      # Independent per-replicate arithmetic confirms the package estimate
      # and CI before aggregating statistical performance.
      center <- sum(change) / n
      se <- sqrt(sum((change - center)^2) / (n - 1L)) / sqrt(n)
      ci <- center + c(-1, 1) * stats::qt(0.975, n - 1L) * se
      # outcomes = "score" declares a metric prefix. Require all expected
      # cells before arithmetic: numeric(0) must never pass a verifier.
      observed <- c(row$score_mean_change, row$score_ci_lower, row$score_ci_upper)
      expected_p <- 2 * stats::pt(-abs(center / se), n - 1L)
      if (length(observed) != 3L || any(!is.finite(observed)) ||
          length(row$score_paired_t_p) != 1L || !is.finite(row$score_paired_t_p) ||
          max(abs(observed - c(center, ci))) > 1e-10 ||
          abs(row$score_paired_t_p - expected_p) > 1e-10) {
        stop(sprintf("Arithmetic mismatch: seed=2026101201 target=%g replicate=%d", target, replicate))
      }
      index <- index + 1L
      rows[[index]] <- data.frame(target, replicate, n, estimate = row$score_mean_change,
        lower = row$score_ci_lower, upper = row$score_ci_upper, p = row$score_paired_t_p)
    }
  }
  draws <- do.call(rbind, rows)
  summary <- do.call(rbind, lapply(split(draws, draws$target), function(d) {
    coverage <- mean(d$lower <= d$target & d$upper >= d$target)
    rejection <- mean(d$p < alpha)
    bias <- mean(d$estimate - d$target)
    bias_mcse <- stats::sd(d$estimate) / sqrt(B)
    nominal_mcse <- sqrt(0.95 * 0.05 / B)
    data.frame(target = d$target[1L], n, replications = B, seed = 2026101201,
      bias, bias_mcse, coverage, coverage_mcse = sqrt(coverage * (1 - coverage) / B),
      rejection, rejection_mcse = sqrt(rejection * (1 - rejection) / B),
      bias_investigate = abs(bias) > 4 * bias_mcse,
      coverage_investigate = abs(coverage - 0.95) > 4 * nominal_mcse,
      type1_investigate = if (d$target[1L] == 0) abs(rejection - alpha) > 4 * nominal_mcse else NA)
  }))
  write.csv(draws, file.path(output, "simulation-replicates.csv"), row.names = FALSE)
  write.csv(summary, file.path(output, "simulation-summary.csv"), row.names = FALSE)
  writeLines(capture.output(sessionInfo()), file.path(output, "simulation-session.txt"))
  print(summary, row.names = FALSE)
  invisible(summary)
}
started <- proc.time()
study()
print(proc.time() - started)
