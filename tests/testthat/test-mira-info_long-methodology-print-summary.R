methodology_print_fit <- function(analyses = "model", inspect_only = FALSE) {
  y <- rbind(c(2, 5, 7), c(4, 6, 12), c(7, 8, 9), c(3, 7, 13),
    c(9, 11, 10), c(6, 9, 16), c(8, 14, 18), c(12, 13, 15))
  data <- data.frame(patient = sprintf("P%02d", 1:8),
    signal_t0 = y[, 1L], signal_t1 = y[, 2L], signal_t2 = y[, 3L])
  mira_info_long(data, id = "patient", outcomes = "Signal score",
    time_vars = c("signal_t0", "signal_t1", "signal_t2"),
    time_labels = c("Zbase", "Amid", "Mfinal"), covariates = character(),
    analyses = analyses, inspect_only = inspect_only, verbose = FALSE)
}

methodology_print_capture <- function(expr) {
  value <- NULL
  text <- capture.output(value <- withVisible(force(expr)))
  list(text = paste(text, collapse = "\n"),
    value = value$value, visible = value$visible)
}

test_that("summary and print retain all added inferential tables and declared prefixes without mutation", {
  withr::local_preserve_seed()
  fit <- methodology_print_fit()
  expect_true(fit$advanced_tests$kenward_roger$performed)
  expect_true(fit$advanced_tests$rank_longitudinal$performed)
  expect_true(fit$advanced_models$mmrm$performed)
  snapshot <- serialize(fit, NULL)
  set.seed(2026102607)
  seed_before <- .Random.seed
  kind_before <- RNGkind()
  summarized <- summary(fit)
  expect_s3_class(summarized, "summary.mira_info_long")
  pairs <- list(
    kenward_roger = fit$advanced_tests$kenward_roger$tidy,
    rank_effects = fit$advanced_tests$rank_longitudinal$effects,
    rank_tests = fit$advanced_tests$rank_longitudinal$tests)
  for (name in names(pairs)) {
    expect_identical(summarized$advanced[[name]], pairs[[name]], info = name)
  }
  expect_true(summarized$advanced$mmrm$performed)
  expect_identical(summarized$advanced$mmrm$coefficients,
    fit$advanced_models$mmrm$coefficients)
  expect_identical(summarized$advanced$mmrm$effect_tests,
    fit$advanced_models$mmrm$effect_tests)
  expect_identical(summarized$advanced$sensitivity, fit$sensitivity)
  tables <- c(pairs, list(mmrm_coefficients = summarized$advanced$mmrm$coefficients,
    mmrm_effect_tests = summarized$advanced$mmrm$effect_tests))
  for (name in names(tables)) {
    table <- tables[[name]]
    map <- attr(table, "mira_metric_names", exact = TRUE)
    expect_true(length(map) > 0L, info = name)
    expect_true(all(startsWith(unname(map), "Signal score_")), info = name)
    expect_true(all(unname(map) %in% names(table)), info = name)
  }
  expected_paths <- list(c("advanced", "kenward_roger"),
    c("advanced", "rank_effects"), c("advanced", "rank_tests"),
    c("advanced", "mmrm", "coefficients"),
    c("advanced", "mmrm", "effect_tests"))
  for (path in expected_paths) {
    expect_true(any(vapply(attr(summarized, "mira_metric_paths", exact = TRUE),
      identical, logical(1L), path)))
  }
  printed <- methodology_print_capture(print(fit, digits = 8, max_rows = Inf,
    correlations = FALSE, outliers = FALSE, trajectories = FALSE,
    missingness = FALSE, plots = FALSE))
  expect_false(printed$visible)
  expect_identical(printed$value, fit)
  expect_match(printed$text, "ADDITIONAL LONGITUDINAL PERSPECTIVES", fixed = TRUE)
  expect_match(printed$text, "Kenward-Roger joint F tests", fixed = TRUE)
  expect_match(printed$text, "Nonparametric relative marginal effects", fixed = TRUE)
  expect_match(printed$text, "MMRM: unstructured / Satterthwaite", fixed = TRUE)
  added_section <- strsplit(printed$text,
    "ADDITIONAL LONGITUDINAL PERSPECTIVES", fixed = TRUE)[[1L]][2L]
  added_section <- strsplit(added_section, "COMPLETE OUTPUT GUIDE", fixed = TRUE)[[1L]][1L]
  for (metric in c("F_scaling", "relative_effect", "estimable", "df_num", "df_den")) {
    expect_match(added_section, metric, fixed = TRUE)
  }
  printed_summary <- methodology_print_capture(print(summarized))
  expect_false(printed_summary$visible)
  expect_identical(printed_summary$value, summarized)
  expect_match(printed_summary$text,
    "Kenward-Roger: available | Relative marginal effects: available | MMRM: available",
    fixed = TRUE)
  hidden <- methodology_print_capture(print(fit, model = FALSE,
    correlations = FALSE, outliers = FALSE, trajectories = FALSE,
    missingness = FALSE, plots = FALSE))
  expect_false(grepl("ADDITIONAL LONGITUDINAL PERSPECTIVES", hidden$text, fixed = TRUE))
  expect_identical(hidden$value, fit)
  expect_identical(serialize(fit, NULL), snapshot)
  expect_identical(.Random.seed, seed_before)
  expect_identical(RNGkind(), kind_before)
})

test_that("disabled and inspect-only public consumers never fit the added methods", {
  withr::local_preserve_seed()
  calls <- character()
  unexpected <- function(...) {
    calls <<- c(calls, "unexpected fit")
    stop("A disabled or inspect-only analysis invoked a fitting adapter.")
  }
  disabled <- with_mocked_bindings(methodology_print_fit(analyses = "none"),
    .mira_run_kenward_roger_long = unexpected,
    .mira_run_rank_longitudinal_long = unexpected,
    .mira_run_mmrm_long = unexpected, .package = "MIRA")
  for (module in list(disabled$advanced_tests$kenward_roger,
                     disabled$advanced_tests$rank_longitudinal,
                     disabled$advanced_models$mmrm)) {
    expect_false(module$performed)
    expect_false(module$requested)
  }
  set.seed(2026102608)
  seed_before <- .Random.seed
  kind_before <- RNGkind()
  printed <- methodology_print_capture(print(disabled))
  expect_false(grepl("ADDITIONAL LONGITUDINAL PERSPECTIVES", printed$text, fixed = TRUE))
  summarized <- summary(disabled)
  expect_null(summarized$advanced$kenward_roger)
  expect_null(summarized$advanced$rank_effects)
  expect_null(summarized$advanced$rank_tests)
  expect_false(summarized$advanced$mmrm$performed)
  expect_null(summarized$advanced$mmrm$coefficients)
  expect_null(summarized$advanced$mmrm$effect_tests)
  detected <- with_mocked_bindings(methodology_print_fit(inspect_only = TRUE),
    .mira_analyse_single_outcome_long = unexpected,
    .mira_run_kenward_roger_long = unexpected,
    .mira_run_rank_longitudinal_long = unexpected,
    .mira_run_mmrm_long = unexpected, .package = "MIRA")
  expect_s3_class(detected, "mira_detect_long")
  expect_null(detected$advanced_tests)
  expect_null(detected$advanced_models)
  detected_print <- methodology_print_capture(print(detected))
  expect_match(detected_print$text, "DETECTION OVERVIEW", fixed = TRUE)
  expect_false(grepl("ADDITIONAL LONGITUDINAL PERSPECTIVES", detected_print$text, fixed = TRUE))
  expect_identical(detected_print$value, detected)
  expect_length(calls, 0L)
  expect_identical(.Random.seed, seed_before)
  expect_identical(RNGkind(), kind_before)
})
