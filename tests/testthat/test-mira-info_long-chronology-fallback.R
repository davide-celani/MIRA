test_that("the direct analyzer sorts parsed visit chronology without a resolved outcome", {
  data <- data.frame(
    patient = paste0("S", 1:5),
    score_t0 = c(2, 4, 7, 9, 13),
    score_t4 = c(3, 7, 8, 11, 14),
    score_t12 = c(6, 8, 12, 15, 19)
  )
  # Without a resolved outcome, the analyzer's own fallback must infer numeric
  # visit order. Labels are supplied in the same unsorted order as the columns.
  result <- MIRA:::.mira_analyse_single_outcome_long(
    data, id = "patient", time_vars = c("score_t12", "score_t0", "score_t4"),
    time_labels = c("Final", "Baseline", "Interim"),
    arm_tests = FALSE, plots = FALSE, model = FALSE, outliers = FALSE,
    correlations = FALSE, verbose = FALSE
  )
  chronological <- c("score_t0", "score_t4", "score_t12")
  expect_identical(result$time_vars, chronological)
  expect_identical(result$time_labels,
                   c(score_t0 = "Baseline", score_t4 = "Interim", score_t12 = "Final"))
  expect_identical(result$change$from, c("score_t0", "score_t0", "score_t4"))
  expect_identical(result$change$to, c("score_t4", "score_t12", "score_t12"))
  expect_equal(result$change$mean_change, c(
    mean(data$score_t4 - data$score_t0), mean(data$score_t12 - data$score_t0),
    mean(data$score_t12 - data$score_t4)
  ))
  expect_equal(result$trajectories$baseline, data$score_t0)
  expect_equal(result$trajectories$final, data$score_t12)
  expect_equal(result$trajectories$absolute_change, data$score_t12 - data$score_t0)
})
