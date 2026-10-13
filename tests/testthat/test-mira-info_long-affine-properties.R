affine_info_data <- function() {
  data.frame(subject_id = paste0("S", 1:8),
    arm = rep(c("control", "active"), each = 4L),
    score_t0 = c(3, 7, 4, 9, 6, 12, 8, 15),
    score_t1 = c(5, NA, 8, 12, 9, 11, 13, 14),
    score_t2 = c(8, 10, 7, 11, 12, 15, NA, 18))
}

affine_info_run <- function(data, direction = "higher", threshold = 0.5,
                            reference = "control") {
  mira_info_long(data, id = "subject_id", outcomes = "score",
    time_vars = c("score_t0", "score_t1", "score_t2"),
    time_labels = c("Base", "Mid", "Final"), arm = "arm", reference_arm = reference,
    covariates = character(0), improvement_direction = direction,
    stable_threshold = threshold, analyses = c("arm_tests", "correlations", "outliers"),
    verbose = FALSE)
}

test_that("positive affine changes preserve inference and transform estimates in their units", {
  data <- affine_info_data()
  original <- affine_info_run(data)
  # For y*=a+b*y with b>0, location estimates and CIs transform affinely;
  # SD/SE scale by b, variance by b^2, paired changes by b. Standardized
  # effects, ranks, correlations, and two-sided tests are invariant. Scaling
  # the practical threshold by b also preserves response classifications.
  for (configuration in list(c(offset = 100, scale = 1), c(offset = 0, scale = 4))) {
    a <- unname(configuration[["offset"]])
    b <- unname(configuration[["scale"]])
    transformed <- data
    transformed[original$time_vars] <- lapply(data[original$time_vars], function(y) a + b * y)
    result <- affine_info_run(transformed, threshold = 0.5 * b)
    for (column in c("mean", "ci_lower", "ci_upper", "median", "q1", "q3", "min", "max")) {
      expect_equal(result$descriptives[[column]], a + b * original$descriptives[[column]],
                   tolerance = 1e-12)
    }
    for (column in c("sd", "se", "iqr")) {
      expect_equal(result$descriptives[[column]], b * original$descriptives[[column]], tolerance = 1e-12)
    }
    expect_equal(result$descriptives$variance, b^2 * original$descriptives$variance, tolerance = 1e-12)
    for (column in c("mean_change", "sd_change", "se_change", "ci_lower", "ci_upper")) {
      expect_equal(result$change[[column]], b * original$change[[column]], tolerance = 1e-12)
    }
    for (column in c("cohens_dz", "paired_t_p", "wilcoxon_p", "paired_t_p_adj", "wilcoxon_p_adj",
                     "increased_n", "decreased_n", "stable_n", "improved_n", "worsened_n")) {
      expect_equal(result$change[[column]], original$change[[column]], tolerance = 1e-12)
    }
    expect_equal(result$variability$grand_mean, a + b * original$variability$grand_mean, tolerance = 1e-12)
    for (column in c("within_subject_sd", "between_subject_sd")) {
      expect_equal(result$variability[[column]], b * original$variability[[column]], tolerance = 1e-12)
    }
    expect_equal(result$variability$ICC, original$variability$ICC, tolerance = 1e-12)
    expect_equal(result$correlations, original$correlations, tolerance = 1e-12)
    expect_equal(result$advanced_tests$friedman$tidy, original$advanced_tests$friedman$tidy,
                 tolerance = 1e-12)
    for (kind in c("time", "change")) {
      table_name <- paste0(kind, "_pairwise")
      observed <- result$arm_analysis[[table_name]]
      expected <- original$arm_analysis[[table_name]]
      estimate <- if (kind == "time") "mean_difference_b_minus_a" else "difference_in_change_b_minus_a"
      for (column in c(estimate, "ci_lower", "ci_upper")) {
        expect_equal(observed[[column]], b * expected[[column]], tolerance = 1e-12)
      }
      for (column in c("hedges_g", "welch_t_p", "wilcoxon_p", "welch_t_p_adj", "wilcoxon_p_adj")) {
        expect_equal(observed[[column]], expected[[column]], tolerance = 1e-12)
      }
    }
  }
})

test_that("sign reversal with reversed improvement direction preserves clinical responses", {
  data <- affine_info_data()
  original <- affine_info_run(data)
  reversed <- data
  reversed[original$time_vars] <- lapply(data[original$time_vars], function(y) -y)
  result <- affine_info_run(reversed, direction = "lower")
  # Reflection reverses signed estimates and swaps CI endpoints. Clinical
  # labels remain the same when the direction of improvement is also reversed.
  expect_equal(result$change$mean_change, -original$change$mean_change, tolerance = 1e-12)
  expect_equal(result$change$ci_lower, -original$change$ci_upper, tolerance = 1e-12)
  expect_equal(result$change$ci_upper, -original$change$ci_lower, tolerance = 1e-12)
  expect_equal(result$change$cohens_dz, -original$change$cohens_dz, tolerance = 1e-12)
  expect_equal(result$change$paired_t_p, original$change$paired_t_p, tolerance = 1e-12)
  expect_equal(result$change$wilcoxon_p, original$change$wilcoxon_p, tolerance = 1e-12)
  expect_equal(result$change$increased_n, original$change$decreased_n)
  expect_equal(result$change$decreased_n, original$change$increased_n)
  expect_equal(result$change$improved_n, original$change$improved_n)
  expect_equal(result$change$worsened_n, original$change$worsened_n)
  expect_identical(result$trajectories$clinical_direction, original$trajectories$clinical_direction)
  expect_equal(result$advanced_tests$friedman$posthoc$tidy$rank_biserial,
               -original$advanced_tests$friedman$posthoc$tidy$rank_biserial, tolerance = 1e-12)
  expect_equal(result$correlations, original$correlations, tolerance = 1e-12)
})

test_that("an entirely unobserved subject changes availability denominators only", {
  data <- affine_info_data()
  original <- affine_info_run(data)
  augmented <- rbind(data, data.frame(subject_id = "S9", arm = "control",
    score_t0 = NA_real_, score_t1 = NA_real_, score_t2 = NA_real_))
  result <- affine_info_run(augmented)
  # With explicit configuration, a subject without outcomes contributes to
  # availability denominators but cannot contribute to an estimate or test.
  expect_identical(result$overview$n_patients, original$overview$n_patients + 1L)
  expect_identical(result$overview$complete_profiles, original$overview$complete_profiles)
  expect_equal(result$overview$complete_profiles_pct,
               original$overview$complete_profiles / 9 * 100)
  expect_identical(result$descriptives$n, original$descriptives$n)
  expect_identical(result$descriptives$missing, original$descriptives$missing + 1L)
  expect_equal(result$descriptives$unavailable_pct,
               (original$descriptives$unavailable + 1) / 9 * 100)
  for (column in c("mean", "sd", "se", "ci_lower", "ci_upper", "variance")) {
    expect_equal(result$descriptives[[column]], original$descriptives[[column]], tolerance = 1e-12)
  }
  expect_equal(result$change, original$change, tolerance = 1e-12)
  expect_equal(result$correlations, original$correlations, tolerance = 1e-12)
  expect_equal(result$variability, original$variability, tolerance = 1e-12)
  expect_equal(result$advanced_tests$friedman$tidy, original$advanced_tests$friedman$tidy,
               tolerance = 1e-12)
  expect_equal(result$arm_analysis$time_pairwise, original$arm_analysis$time_pairwise,
               tolerance = 1e-12)
  expect_equal(result$arm_analysis$change_pairwise, original$arm_analysis$change_pairwise,
               tolerance = 1e-12)
  subject <- result$missing$by_patient[result$missing$by_patient$patient == "S9", ]
  expect_equal(subject$unavailable_n, 3)
  expect_equal(subject$unavailable_pct, 100)
  expect_true(all(is.na(result$trajectories$absolute_change[result$trajectories$patient == "S9"])))
  for (group in c("control", "active")) {
    before <- original$arm_analysis$descriptives[original$arm_analysis$descriptives$arm == group, ]
    after <- result$arm_analysis$descriptives[result$arm_analysis$descriptives$arm == group, ]
    increment <- as.integer(group == "control")
    expect_identical(after$group_n, before$group_n + increment)
    expect_identical(after$n, before$n)
    expect_identical(after$unavailable_n, before$unavailable_n + increment)
  }
})
