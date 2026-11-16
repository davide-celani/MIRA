# Enable with MIRA_RUN_SCIENCE_TESTS=true. Evidence goes to a fresh science-*
# directory under MIRA_VALIDATION_OUTPUT_DIR (an existing external directory),
# or the R session temporary directory when no output base is supplied.
test_that("extended paired-change Monte Carlo study matches independent arithmetic", {
  skip_if_not(mira_validation_campaign_enabled("MIRA_RUN_SCIENCE_TESTS"),
              "Set MIRA_RUN_SCIENCE_TESTS=true for the 1000-replicate study")
  output <- mira_validation_campaign_output("science", mira_validation_campaign_root())
  withr::local_preserve_seed()
  withr::local_seed(2026101201)

  # Preregistered design: 24 independent subjects, normal baseline N(10, 3^2)
  # and independent normal change N(delta, 2^2), delta=0,2; no missingness.
  # Target E[Y1-Y0]=delta. No conclusion extends to informative missingness,
  # Bayesian inference, clinical thresholds, or other model engines.
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
      label <- sprintf("seed=2026101201 target=%g replicate=%d", target, replicate)
      expect_identical(nrow(fit$change), 1L, info = label)
      expect_identical(fit$change$from, "score_t0", info = label)
      expect_identical(fit$change$to, "score_t1", info = label)
      expect_equal(fit$change$n, n, info = label)
      row <- fit$change[1L, ]
      # Check shape and finiteness first so numeric(0) cannot pass arithmetic.
      observed <- c(row$score_mean_change, row$score_ci_lower, row$score_ci_upper,
                    row$score_paired_t_p)
      expect_length(observed, 4L)
      expect_true(all(is.finite(observed)), info = label)
      if (length(observed) != 4L || any(!is.finite(observed))) next
      center <- sum(change) / n
      se <- sqrt(sum((change - center)^2) / (n - 1L)) / sqrt(n)
      ci <- center + c(-1, 1) * stats::qt(0.975, n - 1L) * se
      expected_p <- 2 * stats::pt(-abs(center / se), n - 1L)
      expect_true(max(abs(observed - c(center, ci, expected_p))) <= 1e-10,
                  info = label)
      index <- index + 1L
      rows[[index]] <- data.frame(target, replicate, n,
        estimate = row$score_mean_change, lower = row$score_ci_lower,
        upper = row$score_ci_upper, p = row$score_paired_t_p)
    }
  }
  expect_identical(index, 2L * B)
  if (index != 2L * B) stop("Incomplete Monte Carlo evidence after arithmetic failures.")
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
  # MCSE flags remain investigation signals, not stochastic pass thresholds.
  utils::write.csv(draws, file.path(output, "simulation-replicates.csv"), row.names = FALSE)
  utils::write.csv(summary, file.path(output, "simulation-summary.csv"), row.names = FALSE)
  writeLines(capture.output(sessionInfo()), file.path(output, "simulation-session.txt"))
  message("Monte Carlo evidence: ", output)
  print(summary, row.names = FALSE)
})
