core_mira_data <- function() {
  data.frame(
    patient_id = paste0("P", seq_len(6)),
    score_t0 = c(1, 2, 3, 4, 5, 6),
    score_t1 = c(2, 4, 4, 7, 6, 9),
    score_t2 = c(3, 5, 6, 8, 8, 11)
  )
}

test_that("core summaries remain numerically correct with every optional analysis off", {
  d <- core_mira_data()
  x <- mira_info(d, id = "patient_id", analyses = "none", verbose = FALSE)

  expect_s3_class(x, "mira_info")
  expect_equal(x$overview$n_patients, nrow(d))
  expect_equal(x$overview$n_timepoints, 3L)
  expect_equal(x$overview$complete_profiles, nrow(d))
  expect_equal(x$overview$complete_profiles_pct, 100)
  expect_identical(x$time_vars, c("score_t0", "score_t1", "score_t2"))
  expect_equal(x$descriptives$mean, unname(vapply(d[-1], mean, numeric(1))))
  expect_equal(x$descriptives$median,
               unname(vapply(d[-1], stats::median, numeric(1))))
  expect_equal(x$descriptives$sd, unname(vapply(d[-1], stats::sd, numeric(1))))
  expect_equal(x$descriptives$n, rep(nrow(d), 3))
  expect_equal(nrow(x$change), choose(3, 2))
  bf <- x$change[x$change$from == "score_t0" & x$change$to == "score_t2", ]
  expect_equal(nrow(bf), 1L)
  expect_equal(bf$n, nrow(d))
  expect_equal(bf$mean_change, mean(d$score_t2 - d$score_t0))
  expect_equal(bf$sd_change, stats::sd(d$score_t2 - d$score_t0))
  expect_equal(bf$cohens_dz, mean(d$score_t2 - d$score_t0) /
                 stats::sd(d$score_t2 - d$score_t0))
  expect_equal(x$trajectories$absolute_change, d$score_t2 - d$score_t0)
  expect_equal(x$variability$grand_mean, mean(as.matrix(d[-1])))
  expect_null(x$correlations$pearson)
  expect_null(x$outliers$by_time)
  expect_length(x$plots, 0)
  expect_null(x$model$fitted_model)
})

test_that("missing and nonfinite measurements retain distinct audit counts", {
  d <- core_mira_data()
  d$score_t0[1] <- NA_real_
  d$score_t1[2] <- NaN
  d$score_t1[3] <- Inf
  d$score_t2[4] <- -Inf
  expect_warning(
    x <- mira_info(d, id = "patient_id", analyses = "none", verbose = FALSE),
    "Inf/-Inf values"
  )

  expected_observed <- vapply(d[-1], function(z) sum(is.finite(z)), integer(1))
  expect_equal(x$descriptives$n, unname(expected_observed))
  expect_equal(x$descriptives$missing, c(1, 1, 0))
  expect_equal(x$descriptives$non_finite, c(0, 1, 1))
  expect_equal(x$descriptives$unavailable, c(1, 2, 1))
  expect_equal(x$descriptives$n + x$descriptives$unavailable, rep(nrow(d), 3))
  expect_equal(x$overview$non_finite_values, 2)
  expect_equal(x$overview$complete_profiles, 2)
  expect_equal(x$missing$by_patient$unavailable_n, c(1, 1, 1, 1, 0, 0))
  expect_true(all(is.finite(x$long_data$value[!is.na(x$long_data$value)])))
})

test_that("no pairwise overlap produces a typed empty change table", {
  d <- data.frame(patient_id = 1:4,
                  score_t0 = c(1, 2, NA, NA),
                  score_t1 = c(NA, NA, 3, 4))
  x <- mira_info(d, id = "patient_id", analyses = "none", verbose = FALSE)
  expect_equal(nrow(x$change), 0L)
  expect_true(all(c("from", "to", "paired_t_p", "paired_t_p_adj") %in%
                  names(x$change)))
  expect_true(all(is.na(x$trajectories$absolute_change)))
  expect_equal(x$overview$complete_profiles, 0L)
  expect_equal(x$descriptives$n, c(2L, 2L))
})

test_that("response classification includes the exact stable threshold", {
  d <- data.frame(patient_id = 1:5,
                  score_t0 = rep(10, 5),
                  score_t1 = c(11, 9, 10, 12, 8))
  high <- mira_info(d, id = "patient_id", analyses = "none", verbose = FALSE,
                    improvement_direction = "higher", stable_threshold = 1)
  low <- mira_info(d, id = "patient_id", analyses = "none", verbose = FALSE,
                   improvement_direction = "lower", stable_threshold = 1)
  expect_identical(high$trajectories$direction,
                   c("Stable", "Stable", "Stable", "Increase", "Decrease"))
  expect_identical(high$trajectories$clinical_direction,
                   c("Stable", "Stable", "Stable", "Improved", "Worsened"))
  expect_identical(low$trajectories$clinical_direction,
                   c("Stable", "Stable", "Stable", "Worsened", "Improved"))
  expect_equal(high$change$stable_n, 3L)
  expect_equal(high$change$improved_n, 1L)
  expect_equal(high$change$worsened_n, 1L)
  expect_equal(high$change$stable_pct, 60)
})

test_that("pairwise p-value adjustment matches stats p.adjust within each family", {
  for (method in c("holm", "bonferroni", "BH", "BY", "none")) {
    x <- mira_info(core_mira_data(), id = "patient_id", analyses = "none",
                   p_adjust_method = method, verbose = FALSE)
    for (family in c("paired_t", "wilcoxon")) {
      raw <- x$change[[paste0(family, "_p")]]
      adjusted <- x$change[[paste0(family, "_p_adj")]]
      finite <- is.finite(raw)
      expect_equal(adjusted[finite], stats::p.adjust(raw[finite], method = method))
      expect_true(all(is.na(adjusted[!finite])))
      expect_true(all(adjusted[finite] >= 0 & adjusted[finite] <= 1))
    }
  }
})

test_that("correlations use pairwise complete samples and standard estimators", {
  d <- core_mira_data()
  d$score_t1[2] <- NA_real_
  d$score_t2[5] <- NA_real_
  x <- mira_info(d, id = "patient_id", analyses = "correlations", verbose = FALSE)
  expect_equal(x$correlations$pearson,
               stats::cor(d[-1], use = "pairwise.complete.obs", method = "pearson"))
  expect_equal(x$correlations$spearman,
               stats::cor(d[-1], use = "pairwise.complete.obs", method = "spearman"))
  expect_equal(x$correlations$pairwise_n["score_t1", "score_t2"], 4L)
  expect_equal(x$change$n[x$change$from == "score_t1" &
                          x$change$to == "score_t2"], 4L)
})

test_that("IQR flags identify a known extreme observation and change", {
  d <- data.frame(patient_id = seq_len(10), score_t0 = seq_len(10),
                  score_t1 = c(2:10, 100))
  x <- mira_info(d, id = "patient_id", analyses = "outliers", verbose = FALSE)
  expect_equal(x$outliers$by_time$score_t1$patient, 10L)
  expect_equal(x$outliers$by_time$score_t1$value, 100)
  expect_equal(x$outliers$change$score_t0_to_score_t1$patient, 10L)
  expect_equal(x$outliers$change$score_t0_to_score_t1$change, 90)
  expect_equal(nrow(x$outliers$by_time$score_t0), 0L)
  expect_equal(x$overview$n_patients, 10L)
})

test_that("summary invariants hold over complete and incomplete data", {
  datasets <- list(core_mira_data(), transform(core_mira_data(),
                                               score_t0 = replace(score_t0, 1, NA_real_),
                                               score_t2 = replace(score_t2, 6, NA_real_)))
  for (d in datasets) {
    x <- mira_info(d, id = "patient_id", analyses = "none", verbose = FALSE)
    expect_equal(x$overview$n_timepoints, length(x$time_vars))
    expect_equal(x$overview$n_patients, length(unique(d$patient_id)))
    expect_true(x$overview$complete_profiles >= 0L)
    expect_true(x$overview$complete_profiles <= nrow(d))
    expect_true(all(x$descriptives$n + x$descriptives$unavailable == nrow(d)))
    expect_true(all(x$descriptives$unavailable_pct >= 0 &
                    x$descriptives$unavailable_pct <= 100))
    expect_true(all(x$change$n >= 0 & x$change$n <= nrow(d)))
    expect_equal(nrow(x$long_data), nrow(d) * length(x$time_vars))
    expect_identical(levels(x$long_data$time_label), unname(x$time_labels))
  }
})

test_that("complete-profile ANOVA ICC and variability match hand calculations", {
  d <- data.frame(patient_id = 1:3,
                  score_t0 = c(1, 2, 3),
                  score_t1 = c(2, 4, 6))
  x <- mira_info(d, id = "patient_id", analyses = "none", verbose = FALSE)
  expect_equal(x$variability$between_subject_sd, stats::sd(c(1.5, 3, 4.5)))
  expect_equal(x$variability$within_subject_sd, sqrt(7 / 5))
  expect_equal(x$variability$ms_between, 9 / 2)
  expect_equal(x$variability$ms_within, 7 / 3)
  expect_equal(x$variability$ICC_anova_complete_profiles,
               ((9 / 2) - (7 / 3)) / ((9 / 2) + (7 / 3)))
  expect_equal(x$variability$ICC,
               x$variability$ICC_anova_complete_profiles)
})

test_that("arm summaries use the chosen reference and estimate known differences", {
  d <- data.frame(
    patient_id = seq_len(10),
    arm = rep(c("control", "active"), each = 5),
    score_t0 = seq_len(10),
    score_t1 = seq_len(10) + c(1, 2, 1, 2, 1, 3, 4, 3, 4, 3)
  )
  x <- mira_info(d, id = "patient_id", arm = "arm", reference_arm = "control",
                 analyses = "arm_tests", improvement_direction = "higher",
                 verbose = FALSE)
  expect_true(x$arm_analysis$enabled)
  expect_identical(x$arm_analysis$reference_arm, "control")
  expect_identical(x$arm_analysis$levels, c("control", "active"))
  expect_equal(x$arm_analysis$counts, table(factor(d$arm,
                                                   levels = c("control", "active"))))
  expect_equal(nrow(x$arm_analysis$descriptives), 4L)
  expect_equal(x$arm_analysis$time_pairwise$mean_difference_b_minus_a,
               c(5, 7))
  expect_equal(x$arm_analysis$time_pairwise$n_a, c(5L, 5L))
  expect_equal(x$arm_analysis$time_pairwise$n_b, c(5L, 5L))
  expect_equal(x$arm_analysis$change_descriptives$mean_change, c(1.4, 3.4))
  expect_equal(nrow(x$arm_analysis$baseline_balance), 1L)

  off <- mira_info(d, id = "patient_id", arm = "arm", analyses = "none",
                   verbose = FALSE)
  expect_false(off$arm_analysis$enabled)
  expect_null(off$arm_analysis$descriptives)
  expect_false(off$config$analyses$arm_tests)
})

test_that("plot requests produce ggplot objects for the estimable core views", {
  skip_if_not_installed("ggplot2")
  d <- core_mira_data()
  x <- mira_info(d, id = "patient_id",
                 analyses = c("plots", "correlations"), verbose = FALSE)
  expected <- c("boxplot", "spaghetti", "mean_ci", "change",
                "change_from_baseline", "change_ci", "missingness",
                "correlation_heatmap", "response")
  expect_true(all(expected %in% names(x$plots)))
  expect_true(all(vapply(x$plots[expected], inherits, logical(1L), what = "ggplot")))
  expect_null(x$plot_error)
  expect_null(x$plots$arm_mean_ci)

  off <- mira_info(d, id = "patient_id", analyses = "none", verbose = FALSE)
  expect_error(plot(off, which = "mean_ci"), "not available")
})
