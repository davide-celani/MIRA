# Cell means have no fitted-model oracle: the saturated Gaussian model estimates
# every cell mean, and Var(mean_cell) = sigma^2 / n_cell. The residual sums of
# squares are 3*(2 + 5 + 10) = 51 on 36 - 9 = 27 df, hence sigma^2 = 17/9.
# Unsymmetrical means, unequal counts and nonlexical levels expose key/contrast
# swaps, incorrect equal-vs-observed weighting, signs and uncertainty errors.
estimability_cell_fixture_long <- function() {
  means <- matrix(c(4, 9, 15, 8.5, 3, 21, 12, 23.25, 5.5), 3L, 3L,
    dimnames = list(c("Zbase", "Amid", "Mfinal"), c("Zcontrol", "Bactive", "Aother")))
  counts <- matrix(c(3L, 4L, 5L, 4L, 5L, 3L, 5L, 3L, 4L), 3L, 3L)
  cells <- list()
  for (arm in seq_len(3L)) {
    for (time in seq_len(3L)) {
      n <- counts[time, arm]
      cells[[length(cells) + 1L]] <- data.frame(
        arm_factor = colnames(means)[arm], time_factor = rownames(means)[time],
        value = means[time, arm] + seq_len(n) - (n + 1) / 2)
    }
  }
  data <- do.call(rbind, cells)
  data$arm_factor <- factor(data$arm_factor, levels = colnames(means))
  data$time_factor <- factor(data$time_factor, levels = rownames(means))
  list(data = data, means = means, variances = (17 / 9) / counts)
}

# Validate cardinality and complete scientific keys before aligning numeric
# columns. match() is one-to-one after these checks; it cannot silently accept
# missing or duplicate results, unlike a permissive join or separate sorts.
estimability_check_contrasts_long <- function(observed, expected, keys,
                                              estimate = "estimate") {
  key <- function(x) do.call(paste, c(lapply(x[keys], as.character), sep = "::"))
  observed_key <- key(observed)
  expected_key <- key(expected)
  expect_identical(nrow(observed), nrow(expected))
  expect_identical(anyDuplicated(observed_key), 0L)
  expect_identical(anyDuplicated(expected_key), 0L)
  expect_identical(sort(observed_key), sort(expected_key))
  index <- match(expected_key, observed_key)
  observed <- observed[index, , drop = FALSE]
  expect_identical(key(observed), expected_key)
  # Same saturated fit: 1e-10 covers floating-point linear algebra at unit
  # scale. This is not an inferential equivalence margin or display rounding.
  expect_equal(observed[[estimate]], expected$estimate, tolerance = 1e-10)
  expect_equal(observed$SE, expected$SE, tolerance = 1e-10)
  expect_equal(observed$df, rep(27, nrow(expected)), tolerance = 1e-12)
  critical <- stats::qt(0.95, 27)
  expect_equal(observed$lower.CL, expected$estimate - critical * expected$SE,
    tolerance = 1e-10)
  expect_equal(observed$upper.CL, expected$estimate + critical * expected$SE,
    tolerance = 1e-10)
  if (identical(estimate, "estimate")) {
    ratio <- expected$estimate / expected$SE
    expect_equal(observed$t.ratio, ratio, tolerance = 1e-10)
    expect_equal(observed$p_raw, 2 * stats::pt(-abs(ratio), 27), tolerance = 1e-10)
  }
  invisible(NULL)
}

test_that("emmeans arm and cell means use their complete keys and equal marginal weights", {
  fixture <- estimability_cell_fixture_long()
  means <- fixture$means
  variance <- fixture$variances
  fitted <- stats::lm(value ~ time_factor * arm_factor, fixture$data)
  result <- MIRA:::.mira_run_emmeans_long(fitted, TRUE, "holm", 0.9)
  expect_true(result$performed)
  expect_null(result$error)
  expect_length(result$errors, 0L)
  expected_arm <- data.frame(arm_factor = colnames(means),
    estimate = colMeans(means), SE = sqrt(colSums(variance) / 9))
  estimability_check_contrasts_long(result$arm$tidy, expected_arm,
    "arm_factor", estimate = "emmean")
  expected_cells <- expand.grid(arm_factor = colnames(means),
    time_factor = rownames(means), stringsAsFactors = FALSE)
  positions <- cbind(match(expected_cells$time_factor, rownames(means)),
    match(expected_cells$arm_factor, colnames(means)))
  expected_cells$estimate <- means[positions]
  expected_cells$SE <- sqrt(variance[positions])
  estimability_check_contrasts_long(result$arm_time$tidy, expected_cells,
    c("arm_factor", "time_factor"), estimate = "emmean")
})

test_that("emmeans simple and interaction contrasts match independent cell arithmetic", {
  fixture <- estimability_cell_fixture_long()
  means <- fixture$means
  variance <- fixture$variances
  fitted <- stats::lm(value ~ time_factor * arm_factor, fixture$data)
  result <- MIRA:::.mira_run_emmeans_long(fitted, TRUE, "holm", 0.9)
  arm_pairs <- utils::combn(seq_len(3L), 2L)
  time_pairs <- utils::combn(seq_len(3L), 2L)
  arm_names <- apply(arm_pairs, 2L, function(i) paste(colnames(means)[i], collapse = " - "))
  time_names <- apply(time_pairs, 2L, function(i) paste(rownames(means)[i], collapse = " - "))
  marginal <- colMeans(means)
  marginal_variance <- colSums(variance) / 9
  expected_arm <- data.frame(contrast = arm_names,
    estimate = marginal[arm_pairs[1L, ]] - marginal[arm_pairs[2L, ]],
    SE = sqrt(marginal_variance[arm_pairs[1L, ]] + marginal_variance[arm_pairs[2L, ]]))
  estimability_check_contrasts_long(result$arm$pairwise_tidy, expected_arm, "contrast")

  expected_simple_arm <- expected_simple_time <- expected_interaction <- list()
  for (i in seq_len(3L)) {
    a <- arm_pairs[1L, i]
    b <- arm_pairs[2L, i]
    for (time in seq_len(3L)) {
      expected_simple_arm[[length(expected_simple_arm) + 1L]] <- data.frame(
        contrast = arm_names[i], time_factor = rownames(means)[time],
        estimate = means[time, a] - means[time, b],
        SE = sqrt(variance[time, a] + variance[time, b]))
    }
  }
  for (i in seq_len(3L)) {
    from <- time_pairs[1L, i]
    to <- time_pairs[2L, i]
    for (arm in seq_len(3L)) {
      expected_simple_time[[length(expected_simple_time) + 1L]] <- data.frame(
        contrast = time_names[i], arm_factor = colnames(means)[arm],
        estimate = means[from, arm] - means[to, arm],
        SE = sqrt(variance[from, arm] + variance[to, arm]))
    }
    for (j in seq_len(3L)) {
      a <- arm_pairs[1L, j]
      b <- arm_pairs[2L, j]
      expected_interaction[[length(expected_interaction) + 1L]] <- data.frame(
        arm_factor_pairwise = arm_names[j], time_factor_pairwise = time_names[i],
        estimate = (means[from, a] - means[from, b]) - (means[to, a] - means[to, b]),
        SE = sqrt(variance[from, a] + variance[from, b] + variance[to, a] + variance[to, b]))
    }
  }
  estimability_check_contrasts_long(result$arm_time$simple_arm_tidy,
    do.call(rbind, expected_simple_arm), c("time_factor", "contrast"))
  estimability_check_contrasts_long(result$arm_time$simple_time_tidy,
    do.call(rbind, expected_simple_time), c("arm_factor", "contrast"))
  estimability_check_contrasts_long(result$arm_time$interaction_tidy,
    do.call(rbind, expected_interaction), c("arm_factor_pairwise", "time_factor_pairwise"))
})

estimability_static_profiles_long <- function(with_arm = FALSE) {
  data <- data.frame(patient_factor = factor(rep(seq_len(8L), each = 3L)),
    time_factor = factor(rep(c("Base", "Mid", "Final"), 8L),
      levels = c("Base", "Mid", "Final")),
    value = rep(c(3, 5, 4, 7, 6, 8, 5, 9), each = 3L))
  if (with_arm) data$arm_factor <- factor(rep(c("control", "active"), each = 12L))
  data
}

test_that("RM ANOVA does not label wholly nonestimable inference as performed", {
  # The real engine returns an afex object with NaN F/p: subject values vary,
  # but there is exactly zero change and zero within-subject residual error.
  result <- MIRA:::.mira_run_rm_anova_long(estimability_static_profiles_long(),
    FALSE, character(0), character(0))
  expect_false(result$performed)
  expect_match(paste(result$reason_skipped, collapse = ""), "no finite estimable ANOVA effect", fixed = TRUE)
  expect_s3_class(result$object, "afex_aov")
  expect_identical(result$n_subjects, 8L)
  expect_identical(result$tidy$effect, "time_factor")
  expect_true(all(is.nan(result$tidy$statistic)))
  expect_true(all(is.nan(result$tidy$p_value)))
  expect_true(all(is.nan(result$tables$uncorrected_partial_eta[["Pr(>F)"]])))
  expect_null(result$error)
  expect_true(any(grepl("error SSP matrix", result$warnings, fixed = TRUE)))
  sensitivity <- MIRA:::.mira_build_sensitivity_long(mixed_tests = list(),
    rm_anova = result, friedman = list(performed = FALSE), gee = list(),
    random_slope = list(), robust = list(performed = FALSE), n_default = 24L)
  expect_identical(nrow(sensitivity), 0L)
})

test_that("a valid between-arm ANOVA effect survives unavailable within-subject inference", {
  result <- MIRA:::.mira_run_rm_anova_long(estimability_static_profiles_long(TRUE),
    TRUE, character(0), character(0))
  expect_true(result$performed)
  expect_null(result$reason_skipped)
  expect_s3_class(result$object, "afex_aov")
  # Subject-level values give SS_arm=81/8, SS_error=75/4 and F=81/25
  # with df=1/6. Repeating each subject three times does not increase the
  # between-subject denominator degrees of freedom.
  arm <- result$tidy[result$tidy$effect == "arm_factor", , drop = FALSE]
  expect_identical(nrow(arm), 1L)
  expect_equal(arm$statistic, 81 / 25, tolerance = 1e-12)
  expect_equal(c(arm$df_num, arm$df_den), c(1, 6), tolerance = 0)
  expect_equal(arm$p_value, stats::pf(81 / 25, 1, 6, lower.tail = FALSE), tolerance = 1e-12)
  within <- result$tidy[result$tidy$effect != "arm_factor", , drop = FALSE]
  expect_identical(nrow(within), 2L)
  expect_true(all(is.nan(within$statistic)))
  expect_true(all(is.nan(within$p_value)))
})



test_that("nonestimable real RM inference survives the source-report archive and diagnostics", {
  baseline <- c(3, 5, 4, 7, 6, 8, 5, 9)
  data <- data.frame(subject_id = paste0("S", seq_len(8L)),
    score_t0 = baseline, score_t1 = baseline, score_t2 = baseline)
  result <- mira_info_long(data, id = "subject_id", outcomes = "score",
    covariates = character(0), analyses = "model", verbose = FALSE)
  rm <- result$advanced_tests$rm_anova
  expect_false(rm$performed)
  expect_s3_class(rm$object, "afex_aov")
  expect_false(any(grepl("RM-ANOVA", result$sensitivity$method, fixed = TRUE)))
  directory <- withr::local_tempdir(pattern = "mira-nonestimable-rm-")
  report <- mira_report_freq_long(result, output_dir = directory, format = "html",
    sections = c("models", "diagnostics"), render = FALSE, open = FALSE, quiet = TRUE)
  payload <- readRDS(report$payload)
  retained <- payload$result$advanced_tests$rm_anova
  expect_false(retained$performed)
  expect_identical(retained$reason_skipped, rm$reason_skipped)
  expect_s3_class(retained$object, "afex_aov")
  expect_identical(retained$tidy, rm$tidy)
  manifest <- utils::read.csv(report$results_manifest, stringsAsFactors = FALSE)
  # The archive follows reportability for CSV; the complete RDS above retains
  # the raw failed ANOVA object/tables. No valid-looking RM CSV is published.
  expect_length(which(manifest$outcome == "score" &
    manifest$path == "outcome$advanced_tests$rm_anova$tidy"), 0L)
  entry <- which(manifest$outcome == "score" & manifest$path == "outcome$change")
  expect_length(entry, 1L)
  archived <- utils::read.csv(file.path(directory, manifest$file[[entry]]),
    stringsAsFactors = FALSE)
  expect_identical(nrow(archived), 3L)
  expect_identical(archived$from, c("score_t0", "score_t0", "score_t1"))
  expect_identical(archived$to, c("score_t1", "score_t2", "score_t2"))
  # Observed zero changes differ from unavailable constant-difference t tests.
  expect_equal(archived$score_mean_change, c(0, 0, 0), tolerance = 0)
  expect_true(all(is.na(archived$score_paired_t_p)))
  nodes <- MIRA:::.mira_report_outcome_tree_long(payload$result, payload$report)
  expect_null(nodes$advanced_tests$children$rm_anova)
  diagnostic <- nodes$diagnostics$children[["advanced_tests$rm_anova"]]
  expect_type(diagnostic, "list")
  visible <- paste(capture.output(MIRA:::.mira_report_render_node_long(
    diagnostic, 2L, payload$report)), collapse = "\n")
  expect_match(visible, rm$reason_skipped, fixed = TRUE)
  # A manual corruption removes engine warnings but keeps the actual failed
  # afex object and nonestimability reason. A valid reason must remain visible
  # without relying on a particular afex warning message/version.
  reason_only <- payload$result
  reason_only$advanced_tests$rm_anova$warnings <- character(0)
  reason_only$advanced_tests$rm_anova$error <- NULL
  reason_nodes <- MIRA:::.mira_report_outcome_tree_long(reason_only, payload$report)
  reason_node <- reason_nodes$diagnostics$children[["advanced_tests$rm_anova"]]
  expect_type(reason_node, "list")
  reason_visible <- if (is.null(reason_node)) "" else {
    paste(capture.output(MIRA:::.mira_report_render_node_long(
      reason_node, 2L, payload$report)), collapse = "\n")
  }
  expect_match(reason_visible, rm$reason_skipped, fixed = TRUE)
})

