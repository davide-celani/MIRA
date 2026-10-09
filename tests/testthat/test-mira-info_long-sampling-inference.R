# The rows include complete profiles, monotone dropout, intermittent missingness,
# late entry, a single observed visit, and a wholly unobserved subject.
sampling_info_data <- function() {
  data.frame(subject_id = paste0("S", 1:8),
    score_t0 = c(2, 4, 10, 20, 3, NA, NA, NA),
    score_t1 = c(3, 6, 12, NA, NA, 5, NA, NA),
    score_t2 = c(5, 7, NA, NA, 7, 8, 30, NA),
    score_t3 = c(8, 11, NA, NA, 9, 12, NA, NA),
    stringsAsFactors = FALSE)
}

sampling_info_run <- function(data, ...) {
  args <- list(data = data, id = "subject_id", outcomes = "score",
               covariates = character(0), analyses = "none", verbose = FALSE)
  dots <- list(...)
  args[names(dots)] <- dots
  do.call(mira_info_long, args)
}

test_that("mixed missing profiles preserve observation accounting", {
  data <- sampling_info_data()
  values <- as.matrix(data[-1L])
  available <- !is.na(values)
  result <- sampling_info_run(data, analyses = "correlations")
  expect_identical(result$overview$n_patients, 8L)
  expect_identical(result$overview$complete_profiles, 2L)
  expect_equal(result$overview$complete_profiles_pct, 25)
  expect_equal(result$descriptives$n, c(5, 4, 5, 4))
  expect_equal(result$missing$by_patient$unavailable_n, c(0, 0, 2, 3, 1, 1, 3, 4))
  expect_equal(result$missing$by_patient$unavailable_pct,
               rowSums(!available) / ncol(values) * 100)
  expect_equal(result$missing$by_time$unavailable_n, unname(colSums(!available)))
  expect_equal(sum(result$missing$by_patient$unavailable_n),
               sum(result$missing$by_time$unavailable_n))
  expect_identical(nrow(result$long_data), 32L)
  expect_identical(sum(is.finite(result$long_data$value)), 18L)
  expect_true(all(is.na(result$long_data$value[result$long_data$patient == "S8"])))
  expect_equal(result$descriptives$mean, unname(colMeans(values, na.rm = TRUE)))
  expect_equal(result$correlations$pairwise_n, crossprod(available))
  expect_equal(unname(diag(result$correlations$pairwise_n)), result$descriptives$n)
})

test_that("contrasts use pairwise subjects while Friedman uses complete profiles", {
  data <- sampling_info_data()
  result <- sampling_info_run(data, analyses = "correlations", p_adjust_method = "BY")
  pairs <- utils::combn(result$time_vars, 2L)
  expect_identical(result$change$from, pairs[1L, ])
  expect_identical(result$change$to, pairs[2L, ])
  for (i in seq_len(ncol(pairs))) {
    from <- pairs[1L, i]
    to <- pairs[2L, i]
    keep <- is.finite(data[[from]]) & is.finite(data[[to]])
    x <- data[[from]][keep]
    y <- data[[to]][keep]
    delta <- y - x
    row <- result$change[i, ]
    ttest <- stats::t.test(delta)
    wtest <- stats::wilcox.test(delta, exact = FALSE)
    expect_equal(row$n, sum(keep), info = paste(from, to))
    expect_equal(row$n, result$correlations$pairwise_n[from, to])
    expect_equal(row$mean_from, mean(x))
    expect_equal(row$mean_to, mean(y))
    expect_equal(row$mean_change, mean(delta))
    expect_equal(c(row$ci_lower, row$ci_upper), as.numeric(ttest$conf.int), tolerance = 1e-12)
    expect_equal(row$cohens_dz, mean(delta) / stats::sd(delta), tolerance = 1e-12)
    expect_equal(row$paired_t_p, ttest$p.value, tolerance = 1e-12)
    expect_equal(row$wilcoxon_p, wtest$p.value, tolerance = 1e-12)
    expect_equal(result$correlations$pearson[from, to],
                 stats::cov(x, y) / (stats::sd(x) * stats::sd(y)), tolerance = 1e-12)
    expect_equal(result$correlations$spearman[from, to],
                 stats::cor(rank(x), rank(y)), tolerance = 1e-12)
  }
  expect_equal(result$change$paired_t_p_adj, stats::p.adjust(result$change$paired_t_p, "BY"))
  complete <- as.matrix(data[stats::complete.cases(data[result$time_vars]), result$time_vars])
  friedman <- result$advanced_tests$friedman
  expected_global <- stats::friedman.test(complete)
  expect_true(friedman$performed)
  expect_identical(friedman$tidy$n, 2L)
  expect_equal(friedman$tidy$statistic, unname(expected_global$statistic))
  expect_equal(friedman$tidy$p_raw, expected_global$p.value)
  expect_equal(friedman$posthoc$tidy$n, result$change$n)
  expect_equal(friedman$posthoc$tidy$p_raw, result$change$wilcoxon_p)
  expect_equal(friedman$posthoc$tidy$p_primary, stats::p.adjust(friedman$posthoc$tidy$p_raw, "BY"))
  expect_equal(result$multiplicity$friedman_posthoc$table, friedman$posthoc$tidy)
})

test_that("incomplete-profile variability and complete-profile ICC target different samples", {
  data <- sampling_info_data()
  result <- sampling_info_run(data)
  values <- as.matrix(data[-1L])
  subject_means <- rowMeans(values, na.rm = TRUE)
  residuals <- values - subject_means
  expect_equal(result$variability$grand_mean, mean(values, na.rm = TRUE))
  expect_equal(result$variability$between_subject_sd,
               stats::sd(subject_means[is.finite(subject_means)]))
  expect_equal(result$variability$within_subject_sd, stats::sd(residuals[is.finite(residuals)]))
  expect_equal(as.numeric(result$long_data$patient_mean),
               unname(rep(replace(subject_means, is.nan(subject_means), NA_real_), 4L)))
  complete <- values[stats::complete.cases(values), , drop = FALSE]
  k <- ncol(complete)
  # A one-way ANOVA is an independent mean-square oracle for the balanced ICC.
  balanced <- data.frame(value = as.vector(t(complete)),
                        subject = factor(rep(seq_len(nrow(complete)), each = k)))
  anova <- stats::anova(stats::lm(value ~ subject, data = balanced))
  ms_between <- anova$`Mean Sq`[[1L]]
  ms_within <- anova$`Mean Sq`[[2L]]
  expected_icc <- (ms_between - ms_within) / (ms_between + (k - 1L) * ms_within)
  expect_equal(result$variability$ms_between, ms_between, tolerance = 1e-12)
  expect_equal(result$variability$ms_within, ms_within, tolerance = 1e-12)
  expect_equal(result$variability$df_between, anova$Df[[1L]])
  expect_equal(result$variability$df_within, anova$Df[[2L]])
  expect_equal(result$variability$ICC_anova_complete_profiles, expected_icc, tolerance = 1e-12)
  expect_equal(result$variability$ICC, expected_icc, tolerance = 1e-12)
})

sampling_three_arm_data <- function() {
  data.frame(subject_id = seq_len(14L),
    arm = c(rep("control", 4L), rep("A", 4L), rep("B", 5L), NA_character_),
    score_t0 = c(2, 4, 7, 9, 3, 6, 10, 12, 4, 8, 13, 15, 17, 100),
    score_t1 = c(4, NA, 6, 12, 6, 11, NA, 15, NA, NA, 20, NA, NA, 200),
    score_t2 = c(5, 9, 8, 14, 7, 12, 13, 18, 6, 10, 16, 18, 22, 300),
    stringsAsFactors = FALSE)
}

test_that("unbalanced three-arm contrasts adjust only estimable tests and retain sparse arms", {
  data <- sampling_three_arm_data()
  expect_warning(result <- sampling_info_run(data, arm = "arm", reference_arm = "control",
    analyses = "arm_tests", p_adjust_method = "holm", alpha = 0.1,
    improvement_direction = "higher"), "1 subjects have a missing treatment arm")
  arm <- result$arm_analysis
  expect_identical(arm$levels, c("control", "A", "B"))
  expect_identical(nrow(arm$time_pairwise), 9L)
  expect_identical(nrow(arm$change_pairwise), 6L)
  expect_equal(result$descriptives$n, c(14, 8, 14))
  expect_equal(arm$time_omnibus$n, c(13, 7, 13))
  expect_equal(arm$descriptives$group_n, rep(c(4, 4, 5), each = 3L))
  expect_equal(arm$descriptives$n + arm$descriptives$unavailable_n, arm$descriptives$group_n)
  for (kind in c("time", "change")) {
    table <- arm[[paste0(kind, "_pairwise")]]
    expected_t <- expected_w <- rep(NA_real_, nrow(table))
    for (i in seq_len(nrow(table))) {
      row <- table[i, ]
      value <- if (kind == "time") data[[row$time]] else data[[row$to]] - data[[row$from]]
      a <- value[which(data$arm == row$arm_a & is.finite(value))]
      b <- value[which(data$arm == row$arm_b & is.finite(value))]
      expect_equal(c(row$n_a, row$n_b), c(length(a), length(b)))
      expected_w[[i]] <- stats::wilcox.test(b, a, exact = FALSE)$p.value
      estimate <- if (kind == "time") row$mean_difference_b_minus_a else {
        row$difference_in_change_b_minus_a
      }
      if (min(length(a), length(b)) < 2L) {
        expect_true(all(is.na(c(estimate, row$ci_lower, row$ci_upper, row$hedges_g, row$welch_t_p))))
      } else {
        test <- stats::t.test(b, a, var.equal = FALSE, conf.level = 0.9)
        expected_t[[i]] <- test$p.value
        expect_equal(estimate, mean(b) - mean(a), tolerance = 1e-12)
        expect_equal(c(row$ci_lower, row$ci_upper), as.numeric(test$conf.int), tolerance = 1e-12)
      }
    }
    expect_equal(table$welch_t_p, expected_t, tolerance = 1e-12)
    expect_equal(table$wilcoxon_p, expected_w, tolerance = 1e-12)
    finite <- is.finite(expected_t)
    expect_equal(table$welch_t_p_adj[finite], stats::p.adjust(expected_t[finite], "holm"))
    expect_true(all(is.na(table$welch_t_p_adj[!finite])))
    expect_equal(table$wilcoxon_p_adj, stats::p.adjust(expected_w, "holm"))
  }
  for (i in seq_along(result$time_vars)) {
    variable <- result$time_vars[[i]]
    keep <- !is.na(data$arm)
    tab <- table(factor(data$arm[keep], levels = arm$levels), is.na(data[[variable]][keep]))
    if (ncol(tab) == 2L) {
      expect_identical(arm$missingness$test[[i]], "Fisher")
      expect_equal(arm$missingness$p_value[[i]], stats::fisher.test(tab)$p.value)
    } else {
      expect_true(is.na(arm$missingness$test[[i]]))
      expect_true(is.na(arm$missingness$p_value[[i]]))
    }
  }
  expect_equal(arm$baseline_balance, arm$time_pairwise[1:3, , drop = FALSE])
})

test_that("reversing the arm reference reverses estimates and intervals but preserves inference", {
  data <- sampling_three_arm_data()
  data <- data[which(data$arm %in% c("control", "A")), ]
  run <- function(reference) sampling_info_run(data, arm = "arm",
    reference_arm = reference, analyses = "arm_tests")
  control <- run("control")
  active <- run("A")
  for (kind in c("time", "change")) {
    a <- control$arm_analysis[[paste0(kind, "_pairwise")]]
    b <- active$arm_analysis[[paste0(kind, "_pairwise")]]
    estimate <- if (kind == "time") "mean_difference_b_minus_a" else "difference_in_change_b_minus_a"
    expect_identical(a$arm_a, b$arm_b)
    expect_identical(a$arm_b, b$arm_a)
    expect_equal(a$n_a, b$n_b)
    expect_equal(a$n_b, b$n_a)
    expect_equal(a[[estimate]], -b[[estimate]], tolerance = 1e-12)
    expect_equal(a$ci_lower, -b$ci_upper, tolerance = 1e-12)
    expect_equal(a$ci_upper, -b$ci_lower, tolerance = 1e-12)
    expect_equal(a$hedges_g, -b$hedges_g, tolerance = 1e-12)
    expect_equal(a$welch_t_p, b$welch_t_p, tolerance = 1e-12)
    expect_equal(a$wilcoxon_p, b$wilcoxon_p, tolerance = 1e-12)
    expect_equal(a$welch_t_p_adj, b$welch_t_p_adj, tolerance = 1e-12)
  }
  expect_equal(active$descriptives, control$descriptives)
  expect_equal(active$change, control$change)
  expect_equal(active$correlations, control$correlations)
})

test_that("a singleton subject has estimable means but unavailable sampling uncertainty", {
  data <- data.frame(subject_id = "only", score_t0 = 10, score_t1 = 11, score_t2 = 12)
  result <- sampling_info_run(data, analyses = c("model", "correlations", "outliers"))
  expect_identical(result$overview$n_patients, 1L)
  expect_identical(result$overview$complete_profiles, 1L)
  expect_equal(result$descriptives$mean, c(10, 11, 12))
  expect_identical(result$descriptives$n, rep(1L, 3L))
  expect_true(all(is.na(result$descriptives$sd)))
  expect_true(all(is.na(result$descriptives$ci_lower)))
  expect_true(all(is.na(result$descriptives$ci_upper)))
  expect_equal(result$missing$by_patient$non_finite_n, 0)
  expect_equal(result$missing$by_patient$unavailable_n, 0)
  expect_identical(result$change$n, rep(1L, 3L))
  expect_equal(result$change$mean_change, c(1, 2, 1))
  expect_true(all(is.na(result$change$sd_change)))
  expect_true(all(is.na(result$change$paired_t_p)))
  expect_true(all(is.na(result$change$cohens_dz)))
  expect_null(result$model$fitted_model)
  expect_null(result$correlations$pearson)
  expect_true(all(vapply(result$outliers$by_time, nrow, integer(1L)) == 0L))
  expect_true(is.na(result$variability$between_subject_sd))
  expect_true(is.na(result$variability$ICC))
  expect_equal(result$trajectories$absolute_change, 2)
})

test_that("singleton missing and nonfinite visits retain separate audit counts", {
  data <- data.frame(subject_id = "only", score_t0 = 10, score_t1 = NA_real_,
                     score_t2 = Inf, score_t3 = 12)
  expect_warning(result <- sampling_info_run(data), "1 Inf/-Inf values")
  expect_identical(result$overview$complete_profiles, 0L)
  expect_identical(result$overview$non_finite_values, 1L)
  expect_identical(result$descriptives$n, c(1L, 0L, 0L, 1L))
  expect_identical(result$descriptives$missing, c(0L, 1L, 0L, 0L))
  expect_identical(result$descriptives$non_finite, c(0L, 0L, 1L, 0L))
  expect_identical(result$descriptives$unavailable, c(0L, 1L, 1L, 0L))
  expect_equal(result$missing$by_patient$missing_n, 1)
  expect_equal(result$missing$by_patient$non_finite_n, 1)
  expect_equal(result$missing$by_patient$unavailable_n, 2)
  expect_equal(result$missing$by_patient$unavailable_pct, 50)
  expect_equal(result$long_data$value, c(10, NA_real_, NA_real_, 12))
  expect_identical(result$change$from, "score_t0")
  expect_identical(result$change$to, "score_t3")
  expect_identical(result$change$n, 1L)
  expect_equal(result$change$mean_change, 2)
})

test_that("an entirely missing middle visit preserves estimable endpoint contrasts", {
  data <- data.frame(subject_id = 1:6, score_t0 = c(2, 4, 6, 7, 10, 13),
                     score_t1 = rep(NA_real_, 6L), score_t2 = c(3, 7, 8, 6, 12, 16))
  result <- sampling_info_run(data, analyses = "correlations")
  expect_identical(result$descriptives$n, c(6L, 0L, 6L))
  expect_equal(result$descriptives$unavailable_pct, c(0, 100, 0))
  expect_true(is.na(result$descriptives$mean[[2L]]))
  expect_identical(result$overview$complete_profiles, 0L)
  expect_identical(result$change$from, "score_t0")
  expect_identical(result$change$to, "score_t2")
  expect_identical(result$change$n, 6L)
  expect_equal(result$change$mean_change, mean(data$score_t2 - data$score_t0))
  expect_equal(result$change$paired_t_p, stats::t.test(data$score_t2 - data$score_t0)$p.value)
  expect_true(all(result$correlations$pairwise_n["score_t1", ] == 0L))
  expect_true(all(is.na(result$correlations$pearson["score_t1", ])))
  expect_equal(result$correlations$pearson["score_t0", "score_t2"],
               stats::cor(data$score_t0, data$score_t2))
  friedman <- result$advanced_tests$friedman
  expect_false(friedman$performed)
  expect_identical(friedman$tidy$n, 0L)
  expect_match(friedman$reason_skipped, "complete repeated profiles")
  expect_true(is.character(friedman$error) && nzchar(friedman$error))
  endpoint <- friedman$posthoc$tidy[ friedman$posthoc$tidy$from == "score_t0" &
                                     friedman$posthoc$tidy$to == "score_t2", ]
  expect_identical(endpoint$n, 6L)
  expect_equal(endpoint$p_raw, result$change$wilcoxon_p)
  expect_equal(endpoint$p_primary, endpoint$p_raw)
})
