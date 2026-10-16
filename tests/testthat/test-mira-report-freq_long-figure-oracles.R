# Fixture totals are listed independently below; incomplete visits differ by
# arm, unused factor levels cannot become groups, and supplied visits are out of
# order. These figures must show available-visit means, Welch differences and
# missingness percentages with the full arm size as denominator.
test_that("priority longitudinal figures preserve independent arm values and displayed coordinates", {
  data <- data.frame(
    patient = c("C1", "C2", "C3", "C4", "T1", "T2", "T3", "T4"),
    arm = factor(rep(c("Control", "Treatment"), each = 4L),
                 levels = c("Treatment", "unused", "Control")),
    score_t0 = c(1, 4, 7, 10, 2, 6, 12, 16),
    score_t1 = c(3, NA, 9, 12, 5, 8, NA, 20),
    score_t2 = c(4, 6, NA, 16, 10, 13, 18, 22)
  )
  fit <- mira_info_long(data, id = "patient", outcomes = "score",
    time_vars = c("score_t2", "score_t0", "score_t1"),
    time_labels = c(score_t0 = "Baseline", score_t1 = "Middle", score_t2 = "Final"),
    arm = "arm", reference_arm = "Control", covariates = character(),
    analyses = c("plots", "arm_tests"), verbose = FALSE)
  root <- withr::local_tempdir(pattern = "mira-figure-oracles-")
  report <- report_test_call(list(x = fit, sections = "figures"), root)
  saved <- readRDS(report$payload)$result
  nodes <- .mira_report_outcome_tree_long(saved, report_test_cfg("figures"))$plots$children
  visits <- c("Baseline", "Middle", "Final")
  # Control then Treatment: hand sums/squared sums, excluding only the named
  # missing observations. CI arithmetic uses the standard t distribution.
  n <- c(4, 3, 3, 4, 3, 4)
  sums <- c(22, 24, 26, 36, 33, 63)
  squares <- c(166, 234, 308, 440, 489, 1077)
  expected_mean <- sums / n
  variance <- (squares - sums^2 / n) / (n - 1)
  margin <- stats::qt(0.975, n - 1) * sqrt(variance / n)
  expected_lower <- expected_mean - margin
  expected_upper <- expected_mean + margin
  p <- nodes$arm_mean_ci$value
  expect_identical(as.character(p$data$arm), rep(c("Control", "Treatment"), each = 3L))
  expect_identical(as.character(p$data$time_label), rep(visits, 2L))
  expect_identical(levels(p$data$time_label), visits)
  expect_equal(p$data$n, n, tolerance = 0)
  expect_equal(p$data$mean, expected_mean, tolerance = 1e-12)
  expect_equal(p$data$ci_lower, expected_lower, tolerance = 1e-12)
  expect_equal(p$data$ci_upper, expected_upper, tolerance = 1e-12)
  built <- ggplot2::ggplot_build(p)
  # Check actual graphical mappings too: correct stored data cannot rescue a
  # wrong y/group/errorbar aesthetic.
  points <- built$data[[2L]]
  intervals <- built$data[[3L]]
  expect_equal(as.numeric(points$x), rep(1:3, 2L), tolerance = 0)
  expect_identical(as.integer(points$group), rep(1:2, each = 3L))
  expect_equal(points$y, expected_mean, tolerance = 1e-12)
  expect_equal(intervals$ymin, expected_lower, tolerance = 1e-12)
  expect_equal(intervals$ymax, expected_upper, tolerance = 1e-12)
  expect_identical(p$labels$x, "Time")
  expect_identical(p$labels$y, "Mean SCORE")
  expect_identical(p$labels$shape, "Arm")

  difference <- expected_mean[4:6] - expected_mean[1:3]
  a <- variance[1:3] / n[1:3]
  b <- variance[4:6] / n[4:6]
  se <- sqrt(a + b)
  df <- (a + b)^2 / (a^2 / (n[1:3] - 1) + b^2 / (n[4:6] - 1))
  difference_margin <- stats::qt(0.975, df) * se
  p <- nodes$arm_difference_ci$value
  expect_identical(p$data$comparison, rep("Treatment - Control", 3L))
  expect_identical(levels(p$data$time_label), rev(visits))
  expect_equal(p$data$mean_difference_b_minus_a, difference, tolerance = 1e-12)
  built <- ggplot2::ggplot_build(p)
  expect_identical(as.character(built$layout$layout$comparison), "Treatment - Control")
  expect_equal(as.numeric(built$data[[3L]]$x), 3:1, tolerance = 0)
  expect_equal(built$data[[3L]]$y, difference, tolerance = 1e-12)
  expect_equal(built$data[[2L]]$ymin, difference - difference_margin, tolerance = 1e-12)
  expect_equal(built$data[[2L]]$ymax, difference + difference_margin, tolerance = 1e-12)
  expect_s3_class(p$coordinates, "CoordFlip")
  expect_match(p$labels$y, "Mean difference", fixed = TRUE)
  expect_match(p$labels$y, "95% CI", fixed = TRUE)

  p <- nodes$arm_missingness$value
  expect_identical(levels(p$data$arm), c("Control", "Treatment"))
  expect_identical(as.character(p$data$time_label), rep(visits, 2L))
  built <- ggplot2::ggplot_build(p)
  expect_equal(built$data[[2L]]$y, c(0, 25, 25, 0, 25, 0), tolerance = 0)
  expect_equal(as.numeric(built$data[[2L]]$x), rep(1:3, 2L), tolerance = 0)
  expect_identical(as.integer(built$data[[2L]]$group), rep(1:2, each = 3L))
  expect_identical(p$labels$y, "Unavailable observations (%)")
  # Carry the saved plots through the report renderer to their captions; drawing
  # is real, but this ordinary test does not claim a Quarto document render.
  grDevices::pdf(file.path(root, "priority-figures.pdf"))
  device <- grDevices::dev.cur()
  withr::defer(grDevices::dev.off(device))
  captions <- c(
    arm_mean_ci = "Longitudinal mean profile by arm with confidence intervals.",
    arm_difference_ci = "Estimated arm differences with confidence intervals.",
    arm_missingness = "Percentage of missing observations by arm and time point."
  )
  for (name in names(captions)) {
    rendered <- paste(capture.output(.mira_report_render_node_long(
      nodes[[name]], 3L, report_test_cfg("figures")
    )), collapse = "\n")
    expect_match(rendered, captions[[name]], fixed = TRUE)
    expect_false(grepl("could not be rendered", rendered, fixed = TRUE))
  }
})
