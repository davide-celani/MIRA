# Public orchestration must keep the unadjusted rank schedule separate from
# the outcome/covariate-complete Gaussian model sample. The expected samples
# below are explicit subject-by-visit keys, independent of preparation helpers.
methodology_sample_fixture <- function() {
  y <- rbind(c(2, 5, 7), c(4, 6, 12), c(7, 8, 9), c(3, 7, 13),
    c(9, 11, 10), c(6, 9, 16), c(8, 14, 18), c(12, 13, 15),
    c(3, 8, 12), c(7, 5, 11), c(6, 11, 15), c(5, 9, 14),
    c(11, 14, 17), c(8, 12, 13), c(10, 13, 20), c(13, 16, 19))
  data.frame(patient = sprintf("P%02d", 1:16),
    signal_t0 = y[, 1L], signal_t1 = y[, 2L], signal_t2 = y[, 3L],
    arm = factor(c(rep("Zcontrol", 7L), rep("Aactive", 9L)),
      levels = c("Zcontrol", "Aactive")),
    age = c(31, 42, 38, 55, 47, 36, 61, 35, 44, 58, 49, 33, 53, 63, 40, 46))
}

methodology_sample_public <- function(data, covariates = character(),
                                       analyses = c("model", "arm_tests")) {
  mira_info_long(data, id = "patient", outcomes = "Signal",
    time_vars = c("signal_t0", "signal_t1", "signal_t2"),
    time_labels = c("Zbase", "Amid", "Mfinal"), arm = "arm",
    reference_arm = "Zcontrol", covariates = covariates,
    analyses = analyses, verbose = FALSE)
}

methodology_sample_rank_keys <- function(data) {
  paste(as.character(data$patient), as.character(data$time_label), sep = "::")
}

methodology_sample_mmrm_keys <- function(data) {
  paste(as.character(data$patient_factor), as.character(data$time_factor), sep = "::")
}

test_that("public rank inference retains explicit missing visits outside the adjusted MMRM sample", {
  data <- methodology_sample_fixture()
  data$signal_t1[1L] <- NA_real_
  data$age[3L] <- NA_real_
  fit <- methodology_sample_public(data, "age")
  rank <- fit$advanced_tests$rank_longitudinal
  mmrm <- fit$advanced_models$mmrm
  expect_true(fit$advanced_tests$kenward_roger$requested)
  expect_true(rank$requested)
  expect_true(mmrm$requested)
  expect_true(rank$performed)
  expect_null(rank$error)
  expect_null(rank$reason_skipped)
  expect_identical(rank$n_subjects, 16L)
  expect_identical(rank$n_observations, 47L)
  expect_identical(nrow(rank$data), 48L)
  expect_setequal(as.character(unique(rank$data$patient)), sprintf("P%02d", 1:16))
  full_keys <- paste(rep(sprintf("P%02d", 1:16), each = 3L),
    rep(c("Zbase", "Amid", "Mfinal"), 16L), sep = "::")
  expect_setequal(methodology_sample_rank_keys(rank$data), full_keys)
  unavailable <- rank$data[!is.finite(rank$data$value), , drop = FALSE]
  expect_identical(methodology_sample_rank_keys(unavailable), "P01::Amid")
  expect_false("age" %in% names(rank$data))
  expect_true(rank$unadjusted)
  expect_identical(rank$missingness_assumption, "MCAR")
  expect_true(mmrm$performed)
  expect_true(mmrm$converged)
  expect_null(mmrm$error)
  expect_identical(mmrm$n, 44L)
  expect_identical(mmrm$n_subjects, 15L)
  expected_mmrm <- setdiff(full_keys,
    c("P01::Amid", "P03::Zbase", "P03::Amid", "P03::Mfinal"))
  expect_setequal(methodology_sample_mmrm_keys(mmrm$data), expected_mmrm)
  expect_true(all(is.finite(mmrm$data$age)))
  expect_true("age" %in% all.vars(mmrm$mean_formula))
  comparison <- fit$advanced_models$model_comparison
  mmrm_row <- comparison[comparison$model == "mmrm_unstructured", , drop = FALSE]
  expect_identical(nrow(mmrm_row), 1L)
  expect_identical(mmrm_row$n, 44L)
  expect_identical(mmrm_row$n_subjects, 15L)
  expect_identical(mmrm_row$object_path, "$advanced_models$mmrm$model")
})

test_that("public rank sample does not change when a selected adjustment covariate is missing", {
  data <- methodology_sample_fixture()
  data$age[3L] <- NA_real_
  fit <- methodology_sample_public(data, "age")
  rank <- fit$advanced_tests$rank_longitudinal
  mmrm <- fit$advanced_models$mmrm
  expect_true(rank$performed)
  expect_identical(rank$n_subjects, 16L)
  expect_identical(rank$n_observations, 48L)
  expect_identical(nrow(rank$data), 48L)
  expect_setequal(as.character(unique(rank$data$patient)), sprintf("P%02d", 1:16))
  expect_true(all(is.finite(rank$data$value)))
  expect_true(mmrm$performed)
  expect_identical(mmrm$n_subjects, 15L)
  expect_identical(mmrm$n, 45L)
  expect_false("P03" %in% as.character(mmrm$data$patient_factor))
  expect_true(all(c("TIME", "ARM", "TIME_X_ARM") %in% rank$tests$question))
  effects <- MIRA:::.mira_base_fields_long(rank$effects)
  expect_equal(effects$n[effects$arm == "Zcontrol"], rep(7L, 3L))
  expect_equal(effects$n[effects$arm == "Aactive"], rep(9L, 3L))
})

test_that("public adaptive switches disable all three added model modules", {
  fit <- methodology_sample_public(methodology_sample_fixture(), analyses = "none")
  for (module in list(fit$advanced_tests$kenward_roger,
                     fit$advanced_tests$rank_longitudinal,
                     fit$advanced_models$mmrm)) {
    expect_false(module$requested)
    expect_false(module$performed)
  }
})

test_that("MMRM preserves an entirely unavailable requested baseline interior or final visit", {
  for (unavailable_visit in c("signal_t0", "signal_t1", "signal_t2")) {
    data <- methodology_sample_fixture()
    data[[unavailable_visit]] <- NA_real_
    fit <- methodology_sample_public(data)
    mmrm <- fit$advanced_models$mmrm
    expect_true(mmrm$requested, info = unavailable_visit)
    expect_false(mmrm$performed, info = unavailable_visit)
    expect_identical(mmrm$n_timepoints, 3L, info = unavailable_visit)
    expect_identical(mmrm$covariance_parameters, 6L, info = unavailable_visit)
    expect_match(mmrm$reason_skipped, "empty modeled visit", fixed = TRUE,
      info = unavailable_visit)
    expect_null(mmrm$model, info = unavailable_visit)
    expect_null(mmrm$coefficients, info = unavailable_visit)
    expect_null(mmrm$effect_tests, info = unavailable_visit)
    expect_false(mmrm$emmeans$performed, info = unavailable_visit)
    expect_null(mmrm$emmeans$baseline_final, info = unavailable_visit)
    expect_null(fit$multiplicity$mmrm, info = unavailable_visit)
    expect_false(any(fit$advanced_models$model_comparison$model == "mmrm_unstructured"),
      info = unavailable_visit)
    expect_false(any(grepl("^MMRM:", fit$sensitivity$method)), info = unavailable_visit)
    # This gate concerns only the newly added marginal model: the existing
    # estimable two-visit mixed-model branch remains available for inspection.
    expect_true(inherits(fit$model$fitted_model, "lmerMod"), info = unavailable_visit)
    expect_true(fit$advanced_tests$emmeans$performed, info = unavailable_visit)
  }
})

test_that("covariate exclusion cannot silently remove an interior MMRM visit", {
  data <- methodology_sample_fixture()
  observed_middle <- c(1L, 2L, 8L, 9L)
  data$signal_t1[-observed_middle] <- NA_real_
  data$age[observed_middle] <- NA_real_
  fit <- methodology_sample_public(data, "age")
  module <- fit$advanced_models$mmrm
  expect_true(module$requested)
  expect_false(module$performed)
  expect_identical(module$n_timepoints, 3L)
  expect_identical(module$covariance_parameters, 6L)
  expect_identical(module$n, 24L)
  expect_identical(module$n_subjects, 12L)
  expect_match(module$reason_skipped, "empty modeled visit", fixed = TRUE)
  expect_null(module$model)
  expect_null(module$coefficients)
  expect_null(module$effect_tests)
  expect_null(module$emmeans$baseline_final)
  expect_null(fit$multiplicity$mmrm)
  expect_false(any(fit$advanced_models$model_comparison$model == "mmrm_unstructured"))
  expect_false(any(grepl("^MMRM:", fit$sensitivity$method)))
  rank <- fit$advanced_tests$rank_longitudinal
  expected <- paste(rep(sprintf("P%02d", 1:16), each = 3L),
    rep(c("Zbase", "Amid", "Mfinal"), 16L), sep = "::")
  expect_true(rank$requested)
  expect_identical(rank$n_subjects, 16L)
  expect_identical(rank$n_observations, 36L)
  expect_identical(nrow(rank$data), 48L)
  expect_identical(anyDuplicated(methodology_sample_rank_keys(rank$data)), 0L)
  expect_identical(sort(methodology_sample_rank_keys(rank$data)), sort(expected))
  observed <- rank$data$patient[rank$data$time_label == "Amid" & is.finite(rank$data$value)]
  expect_identical(sort(as.character(observed)), sprintf("P%02d", observed_middle))
  expect_false("age" %in% names(rank$data))
})

test_that("public outcomes keep their own samples and survive a single rank adapter failure", {
  data <- methodology_sample_fixture()
  variables <- paste0("signal_t", 0:2)
  second <- paste0("other_t", 0:2)
  data[second] <- lapply(data[variables], function(value) 30 - 2 * value)
  data[4L, second] <- NA_real_
  data$signal_t1[1L] <- NA_real_
  input <- data
  run <- function() mira_info_long(data, id = "patient",
    time_vars = list(Signal = variables, Other = second),
    time_labels = list(Signal = c("Zbase", "Amid", "Mfinal"),
      Other = c("Zbase", "Amid", "Mfinal")),
    arm = "arm", reference_arm = "Zcontrol", covariates = character(),
    analyses = c("model", "arm_tests"), verbose = FALSE)
  reference <- run()
  expect_identical(names(reference$outcomes), c("Signal", "Other"))
  visits <- c("Zbase", "Amid", "Mfinal")
  for (name in names(reference$outcomes)) {
    one <- reference$outcomes[[name]]
    vars <- if (name == "Signal") variables else second
    ids <- if (name == "Signal") data$patient else data$patient[-4L]
    expected_n <- if (name == "Signal") 47L else 45L
    expect_true(one$advanced_tests$kenward_roger$performed, info = name)
    expect_true(one$advanced_tests$rank_longitudinal$performed, info = name)
    expect_true(one$advanced_models$mmrm$performed, info = name)
    rank <- one$advanced_tests$rank_longitudinal
    expect_identical(rank$n_subjects, length(ids), info = name)
    expect_identical(rank$n_observations, expected_n, info = name)
    calendar <- paste(rep(ids, each = 3L), rep(visits, length(ids)), sep = "::")
    expect_setequal(methodology_sample_rank_keys(rank$data), calendar)
    finite <- unlist(lapply(seq_along(vars), function(k) {
      paste(data$patient[is.finite(data[[vars[[k]]]])], visits[[k]], sep = "::")
    }), use.names = FALSE)
    expect_setequal(methodology_sample_mmrm_keys(one$advanced_models$mmrm$data),
      finite)
    frame <- stats::model.frame(one$model$fitted_model)
    expect_setequal(paste(frame$patient_factor, frame$time_factor, sep = "::"),
      finite)
    # Independent empirical-distribution arithmetic, in the same explicit
    # arm/visit keys as the output. Each of the six reference cells has weight
    # 1/6 even with unequal group sizes or missing observations.
    cells <- lapply(c("Zcontrol", "Aactive"), function(arm) {
      lapply(vars, function(variable) {
        value <- data[[variable]][data$arm == arm]
        value[is.finite(value)]
      })
    })
    cells <- unlist(cells, recursive = FALSE)
    expected <- vapply(cells, function(target) {
      mean(vapply(cells, function(reference_cell) {
        mean(outer(reference_cell, target, function(a, b) (a < b) + .5 * (a == b)))
      }, numeric(1L)))
    }, numeric(1L))
    effects <- MIRA:::.mira_base_fields_long(rank$effects)
    expect_identical(effects$arm, rep(c("Zcontrol", "Aactive"), each = 3L))
    expect_identical(effects$time_label, rep(visits, 2L))
    expect_identical(effects$time, rep(vars, 2L))
    expect_equal(effects$relative_effect, expected, tolerance = 1e-12, info = name)
  }
  original_rank <- getFromNamespace(".mira_run_rank_longitudinal_long", "MIRA")
  called <- character()
  failed <- with_mocked_bindings(run(), .mira_run_rank_longitudinal_long = function(long_data, ...) {
    called <<- c(called, unique(long_data$time)[[1L]])
    if (identical(unique(long_data$time)[[1L]], "signal_t0")) {
      warning("controlled rank diagnostic")
      stop("controlled rank failure for Signal")
    }
    original_rank(long_data, ...)
  }, .package = "MIRA")
  expect_identical(called, c("signal_t0", "other_t0"))
  broken <- failed$outcomes$Signal$advanced_tests$rank_longitudinal
  expect_true(broken$requested)
  expect_false(broken$performed)
  expect_identical(broken$error, "controlled rank failure for Signal")
  expect_identical(broken$warnings, "controlled rank diagnostic")
  expect_match(broken$reason_skipped, "requested extended method could not be completed")
  expect_identical(failed$diagnostics$failed_outcomes, character())
  for (name in names(reference$outcomes)) {
    before <- reference$outcomes[[name]]
    after <- failed$outcomes[[name]]
    expect_identical(after$descriptives, before$descriptives, info = name)
    expect_identical(after$change, before$change, info = name)
    expect_equal(lme4::fixef(after$model$fitted_model),
      lme4::fixef(before$model$fitted_model), tolerance = 1e-10, info = name)
    expect_equal(after$advanced_tests$kenward_roger$tidy,
      before$advanced_tests$kenward_roger$tidy, tolerance = 1e-10, info = name)
    expect_equal(after$advanced_models$mmrm$coefficients,
      before$advanced_models$mmrm$coefficients, tolerance = 1e-10, info = name)
    expect_equal(after$advanced_models$mmrm$effect_tests,
      before$advanced_models$mmrm$effect_tests, tolerance = 1e-10, info = name)
    expect_equal(after$advanced_tests$emmeans$time$tidy,
      before$advanced_tests$emmeans$time$tidy, tolerance = 1e-10, info = name)
    expect_equal(after$advanced_tests$rm_anova$tidy,
      before$advanced_tests$rm_anova$tidy, tolerance = 1e-10, info = name)
    expect_identical(after$advanced_tests$friedman$tidy,
      before$advanced_tests$friedman$tidy, info = name)
    expect_equal(after$advanced_models$model_comparison,
      before$advanced_models$model_comparison, tolerance = 1e-10, info = name)
    old_rows <- !grepl("^Nonparametric relative", before$sensitivity$method)
    retained <- !grepl("^Nonparametric relative", after$sensitivity$method)
    old <- before$sensitivity[old_rows, , drop = FALSE]
    kept <- after$sensitivity[retained, , drop = FALSE]
    rownames(old) <- rownames(kept) <- NULL
    expect_equal(kept, old, tolerance = 1e-10, info = name)
  }
  expect_equal(failed$outcomes$Other$advanced_tests$rank_longitudinal$effects,
    reference$outcomes$Other$advanced_tests$rank_longitudinal$effects, tolerance = 1e-12)
  expect_equal(failed$outcomes$Other$advanced_tests$rank_longitudinal$tests,
    reference$outcomes$Other$advanced_tests$rank_longitudinal$tests, tolerance = 1e-12)
  expect_false(any(grepl("^Nonparametric relative", failed$outcomes$Signal$sensitivity$method)))
  expect_identical(data, input)
})

test_that("clinical direction and practical thresholds do not change model inference", {
  data <- methodology_sample_fixture()
  run <- function(direction, threshold) mira_info_long(data, id = "patient", outcomes = "Signal",
    time_vars = paste0("signal_t", 0:2), time_labels = c("Zbase", "Amid", "Mfinal"),
    arm = "arm", reference_arm = "Zcontrol", covariates = character(),
    analyses = c("model", "arm_tests"), improvement_direction = direction,
    stable_threshold = threshold, verbose = FALSE)
  original <- run("higher", 0)
  interpreted <- run("lower", 100)
  for (name in c("kenward_roger", "rank_longitudinal")) {
    before <- original$advanced_tests[[name]]
    after <- interpreted$advanced_tests[[name]]
    expect_true(before$performed, info = name)
    expect_true(after$performed, info = name)
    for (table in if (name == "kenward_roger") "tidy" else c("effects", "tests", "covariance")) {
      expect_equal(after[[table]], before[[table]], tolerance = 1e-10, info = paste(name, table))
    }
  }
  expect_true(original$advanced_models$mmrm$performed)
  expect_true(interpreted$advanced_models$mmrm$performed)
  for (table in c("coefficients", "effect_tests", "covariance")) {
    expect_equal(interpreted$advanced_models$mmrm[[table]],
      original$advanced_models$mmrm[[table]], tolerance = 1e-10, info = table)
  }
  expect_equal(interpreted$advanced_models$mmrm$emmeans$baseline_final$tidy,
    original$advanced_models$mmrm$emmeans$baseline_final$tidy, tolerance = 1e-10)
  expect_equal(lme4::fixef(interpreted$model$fitted_model),
    lme4::fixef(original$model$fitted_model), tolerance = 1e-10)
  expect_equal(interpreted$sensitivity, original$sensitivity, tolerance = 1e-10)
  first <- MIRA:::.mira_base_fields_long(original$trajectories)
  second <- MIRA:::.mira_base_fields_long(interpreted$trajectories)
  expect_true(any(first$clinical_direction == "Improved"))
  expect_identical(second$clinical_direction, rep("Stable", nrow(data)))
  expect_identical(second$absolute_change, first$absolute_change)
})
