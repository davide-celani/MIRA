oracle_info <- function(data, ...) {
  args <- list(
    data = data,
    analyses = "none",
    covariates = character(0),
    verbose = FALSE
  )
  dots <- list(...)
  args[names(dots)] <- dots
  do.call(mira_info, args)
}

oracle_descriptive_row <- function(x, alpha) {
  x <- x[is.finite(x)]
  n <- length(x)
  mean_x <- mean(x)
  sd_x <- stats::sd(x)
  se_x <- sd_x / sqrt(n)
  critical <- stats::qt(1 - alpha / 2, df = n - 1L)
  quartiles <- stats::quantile(
    x,
    probs = c(0.25, 0.75),
    names = FALSE,
    type = 7
  )
  c(
    n = n,
    mean = mean_x,
    sd = sd_x,
    variance = stats::var(x),
    se = se_x,
    ci_lower = mean_x - critical * se_x,
    ci_upper = mean_x + critical * se_x,
    median = stats::median(x),
    q1 = quartiles[[1L]],
    q3 = quartiles[[2L]],
    iqr = diff(quartiles),
    min = min(x),
    max = max(x),
    cv_percent = sd_x / abs(mean_x) * 100
  )
}

oracle_hedges_g <- function(a, b) {
  n_a <- length(a)
  n_b <- length(b)
  df <- n_a + n_b - 2L
  pooled_variance <- (
    (n_a - 1L) * stats::var(a) + (n_b - 1L) * stats::var(b)
  ) / df
  correction <- 1 - 3 / (4 * df - 1)
  correction * (mean(b) - mean(a)) / sqrt(pooled_variance)
}

test_that("descriptives and confidence intervals agree with independent base-R calculations", {
  data <- data.frame(
    patient_id = seq_len(6),
    score_t0 = c(1, 2, 4, 8, NA, 16),
    score_t1 = c(3, 5, 7, 9, 11, 13)
  )
  alpha <- 0.10
  result <- oracle_info(data, id = "patient_id", alpha = alpha)

  expected <- rbind(
    oracle_descriptive_row(data$score_t0, alpha),
    oracle_descriptive_row(data$score_t1, alpha)
  )
  observed <- as.matrix(result$descriptives[, colnames(expected), drop = FALSE])

  expect_equal(observed, expected, tolerance = 1e-12)
  expect_identical(result$descriptives$missing, c(1L, 0L))
  expect_identical(result$descriptives$unavailable, c(1L, 0L))
  expect_equal(result$descriptives$missing_pct, c(100 / 6, 0))
  expect_equal(result$settings$confidence_level, 0.90)
  expect_equal(result$settings$confidence_percent, 90)
})

test_that("pairwise change estimates, tests, intervals, and response counts have independent oracles", {
  data <- data.frame(
    patient_id = seq_len(5),
    score_t0 = c(1, 2, 4, 7, 11),
    score_t1 = c(2, 4, 3, 10, 15)
  )
  alpha <- 0.05
  result <- oracle_info(
    data,
    id = "patient_id",
    alpha = alpha,
    improvement_direction = "higher",
    stable_threshold = 1
  )
  row <- result$change[1L, ]
  delta <- data$score_t1 - data$score_t0
  delta_se <- stats::sd(delta) / sqrt(length(delta))
  critical <- stats::qt(1 - alpha / 2, df = length(delta) - 1L)
  paired_t <- stats::t.test(
    data$score_t1,
    data$score_t0,
    paired = TRUE,
    conf.level = 1 - alpha
  )
  paired_wilcox <- suppressWarnings(stats::wilcox.test(
    data$score_t1,
    data$score_t0,
    paired = TRUE,
    exact = FALSE
  ))

  expect_identical(row$n, 5L)
  expect_equal(row$mean_from, mean(data$score_t0))
  expect_equal(row$mean_to, mean(data$score_t1))
  expect_equal(row$mean_change, mean(delta))
  expect_equal(row$sd_change, stats::sd(delta))
  expect_equal(row$se_change, delta_se)
  expect_equal(row$ci_lower, mean(delta) - critical * delta_se)
  expect_equal(row$ci_upper, mean(delta) + critical * delta_se)
  expect_equal(row$cohens_dz, mean(delta) / stats::sd(delta))
  expect_equal(row$paired_t_p, unname(paired_t$p.value))
  expect_equal(row$wilcoxon_p, unname(paired_wilcox$p.value))
  expect_equal(row$paired_t_p_adj, row$paired_t_p)
  expect_equal(row$wilcoxon_p_adj, row$wilcoxon_p)

  # Changes of exactly +/-1 are stable because the implementation uses strict
  # inequalities at both response thresholds.
  expect_identical(
    unname(unlist(row[c("increased_n", "decreased_n", "stable_n")], use.names = FALSE)),
    c(3L, 0L, 2L)
  )
  expect_identical(
    unname(unlist(row[c("improved_n", "worsened_n")], use.names = FALSE)),
    c(3L, 0L)
  )
  expect_equal(row$increased_pct + row$decreased_pct + row$stable_pct, 100)
  expect_equal(row$improved_pct + row$worsened_pct + row$stable_pct, 100)
  expect_equal(
    result$advanced_tests$friedman$posthoc$tidy$rank_biserial,
    0.8
  )
})

test_that("perfect positive and negative correlations and pair counts are exact", {
  data <- data.frame(
    patient_id = seq_len(6),
    score_t0 = 1:6,
    score_t1 = 2 * (1:6) + 3,
    score_t2 = 7:2
  )
  result <- mira_info(
    data,
    id = "patient_id",
    analyses = "correlations",
    covariates = character(0),
    verbose = FALSE
  )
  expected <- matrix(
    c(1, 1, -1, 1, 1, -1, -1, -1, 1),
    nrow = 3L,
    dimnames = list(names(data)[-1L], names(data)[-1L])
  )

  expect_equal(result$correlations$pearson, expected, tolerance = 1e-15)
  expect_equal(result$correlations$spearman, expected, tolerance = 1e-15)
  expect_identical(
    result$correlations$pairwise_n,
    matrix(
      6L,
      nrow = 3L,
      ncol = 3L,
      dimnames = list(names(data)[-1L], names(data)[-1L])
    )
  )
})

test_that("relative trajectory change uses the absolute baseline denominator", {
  data <- data.frame(
    patient_id = seq_len(4),
    score_t0 = c(-10, 0, 5, NA_real_),
    score_t1 = c(-5, 2, 0, 3)
  )
  result <- oracle_info(data, id = "patient_id")

  expect_equal(
    result$trajectories$absolute_change,
    c(5, 2, -5, NA_real_)
  )
  expect_equal(
    result$trajectories$relative_change_percent,
    c(50, NA_real_, -100, NA_real_)
  )
  expect_identical(
    result$trajectories$direction,
    c("Increase", "Increase", "Decrease", NA_character_)
  )
})

test_that("arm Welch tests, confidence intervals, Hedges g, and adjustments match base R", {
  arm_a_t0 <- c(1, 2, 3, 5, 8)
  arm_b_t0 <- c(2, 4, 7, 8, 12)
  change_a <- c(1, 2, 1, 3, 2)
  change_b <- c(3, 5, 4, 6, 5)
  data <- data.frame(
    patient_id = seq_len(10),
    arm = rep(c("control", "active"), each = 5L),
    score_t0 = c(arm_a_t0, arm_b_t0),
    score_t1 = c(arm_a_t0 + change_a, arm_b_t0 + change_b)
  )
  result <- mira_info(
    data,
    id = "patient_id",
    arm = "arm",
    reference_arm = "control",
    analyses = "arm_tests",
    covariates = character(0),
    p_adjust_method = "BH",
    verbose = FALSE
  )

  time_values <- list(
    list(a = arm_a_t0, b = arm_b_t0),
    list(a = arm_a_t0 + change_a, b = arm_b_t0 + change_b)
  )
  for (i in seq_along(time_values)) {
    expected <- stats::t.test(
      time_values[[i]]$b,
      time_values[[i]]$a,
      var.equal = FALSE
    )
    observed <- result$arm_analysis$time_pairwise[i, ]
    expect_equal(observed$mean_difference_b_minus_a,
                 mean(time_values[[i]]$b) - mean(time_values[[i]]$a))
    expect_equal(observed$ci_lower, unname(expected$conf.int[[1L]]))
    expect_equal(observed$ci_upper, unname(expected$conf.int[[2L]]))
    expect_equal(observed$welch_t_p, unname(expected$p.value))
    expect_equal(observed$hedges_g,
                 oracle_hedges_g(time_values[[i]]$a, time_values[[i]]$b))
  }
  expect_equal(
    result$arm_analysis$time_pairwise$welch_t_p_adj,
    stats::p.adjust(result$arm_analysis$time_pairwise$welch_t_p, method = "BH")
  )
  expect_equal(
    result$arm_analysis$time_pairwise$wilcoxon_p_adj,
    stats::p.adjust(result$arm_analysis$time_pairwise$wilcoxon_p, method = "BH")
  )

  change_test <- stats::t.test(change_b, change_a, var.equal = FALSE)
  change_row <- result$arm_analysis$change_pairwise[1L, ]
  expect_equal(change_row$difference_in_change_b_minus_a,
               mean(change_b) - mean(change_a))
  expect_equal(change_row$ci_lower, unname(change_test$conf.int[[1L]]))
  expect_equal(change_row$ci_upper, unname(change_test$conf.int[[2L]]))
  expect_equal(change_row$welch_t_p, unname(change_test$p.value))
  expect_equal(change_row$hedges_g, oracle_hedges_g(change_a, change_b))

  for (i in seq_along(time_values)) {
    combined <- c(time_values[[i]]$a, time_values[[i]]$b)
    group <- factor(rep(c("control", "active"), each = 5L),
                    levels = c("control", "active"))
    expect_equal(
      result$arm_analysis$time_omnibus$welch_anova_p[[i]],
      unname(stats::oneway.test(combined ~ group, var.equal = FALSE)$p.value)
    )
    expect_equal(
      result$arm_analysis$time_omnibus$kruskal_p[[i]],
      unname(stats::kruskal.test(combined ~ group)$p.value)
    )
  }
})

test_that("arm missingness switches from Fisher to chi-square exactly at expected count five", {
  make_data <- function(per_arm) {
    n <- 2L * per_arm
    baseline <- seq_len(n)
    missing_per_arm <- per_arm %/% 2L
    baseline[c(
      seq_len(missing_per_arm),
      per_arm + seq_len(missing_per_arm)
    )] <- NA_real_
    data.frame(
      patient_id = seq_len(n),
      arm = rep(c("control", "active"), each = per_arm),
      score_t0 = baseline,
      score_t1 = seq_len(n) + 1
    )
  }

  below <- mira_info(
    make_data(8L),
    id = "patient_id",
    arm = "arm",
    analyses = "arm_tests",
    covariates = character(0),
    verbose = FALSE
  )
  at <- mira_info(
    make_data(10L),
    id = "patient_id",
    arm = "arm",
    analyses = "arm_tests",
    covariates = character(0),
    verbose = FALSE
  )

  expect_identical(below$arm_analysis$missingness$test[[1L]], "Fisher")
  expect_identical(at$arm_analysis$missingness$test[[1L]], "Chi-square")
  expect_equal(below$arm_analysis$missingness$p_value[[1L]], 1)
  expect_equal(at$arm_analysis$missingness$p_value[[1L]], 1)
})

test_that("p-value family adjustment handles empty and wholly non-finite families", {
  empty <- MIRA:::.mira_adjust_family(data.frame(p_raw = numeric(0)), primary = "holm")
  expect_identical(nrow(empty), 0L)
  expect_identical(
    names(empty),
    c("p_raw", "p_primary", "primary_method", "p_bonferroni", "p_holm", "p_bh", "p_by")
  )
  expect_identical(vapply(empty, typeof, character(1L)), c(
    p_raw = "double", p_primary = "double", primary_method = "character",
    p_bonferroni = "double", p_holm = "double", p_bh = "double", p_by = "double"
  ))

  raw <- c(NA_real_, NaN, Inf, -Inf)
  non_finite <- MIRA:::.mira_adjust_family(
    data.frame(p_raw = raw),
    primary = "BH"
  )
  expect_identical(non_finite$p_raw, raw)
  for (column in c("p_primary", "p_bonferroni", "p_holm", "p_bh", "p_by")) {
    expect_true(all(is.na(non_finite[[column]])))
  }

  one_finite <- MIRA:::.mira_adjust_family(
    data.frame(p_raw = c(NA_real_, 0.03, Inf)),
    primary = "hochberg"
  )
  expect_equal(one_finite$p_primary, c(NA_real_, 0.03, NA_real_))
  expect_identical(one_finite$primary_method, rep("hochberg", 3L))
})

test_that("model formulas evaluate non-syntactic additive covariates and the arm interaction", {
  formula <- MIRA:::.mira_make_fixed_formula(
    arm_enabled = TRUE,
    covariates = c("age years", "site+code")
  )
  model_data <- data.frame(
    value = seq_len(8),
    time_factor = factor(rep(c("Base", "Follow"), 4L)),
    arm_factor = factor(rep(c("Control", "Active"), each = 4L)),
    check.names = FALSE
  )
  model_data[["age years"]] <- seq(40, 47)
  model_data[["site+code"]] <- factor(rep(c("N", "S"), 4L))
  matrix <- stats::model.matrix(formula, data = model_data)

  expect_identical(
    colnames(matrix),
    c(
      "(Intercept)", "time_factorFollow", "arm_factorControl",
      "`age years`", "`site+code`S",
      "time_factorFollow:arm_factorControl"
    )
  )
  expect_equal(nrow(matrix), nrow(model_data))
  expect_equal(unname(matrix[, "`age years`"]), model_data[["age years"]])
})

test_that("rank-deficient Wald inference uses only estimable covariance directions", {
  covariance <- matrix(c(1, 1, 1, 1), nrow = 2L)
  result <- MIRA:::.mira_wald_chisq(
    beta = c(1, 1),
    covariance = covariance,
    indices = 1:2
  )

  expect_identical(result$df, 1L)
  expect_equal(result$statistic, 1, tolerance = 1e-12)
  expect_equal(
    result$p_value,
    stats::pchisq(1, df = 1, lower.tail = FALSE),
    tolerance = 1e-12
  )
})
