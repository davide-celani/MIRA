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

test_that("MMRM never relabels a follow-up as an entirely unavailable requested baseline or final visit", {
  for (unavailable_visit in c("signal_t0", "signal_t2")) {
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
