canonical_contract_data <- function() {
  withr::local_preserve_seed()
  set.seed(78321)
  n <- 24L
  subject <- rnorm(n, 12, 2)
  data.frame(patient = sprintf("P%02d", seq_len(n)),
    arm = factor(rep(c("control", "active"), each = n / 2), levels = c("control", "active")),
    marker_t0 = subject + rnorm(n), marker_t1 = subject + rnorm(n, 1),
    marker_t2 = subject + rnorm(n, 2))
}

canonical_contract_fit <- function(data, declared = NULL, analyses = "none") {
  args <- list(data = data, id = "patient", time_vars = grep("^marker_t", names(data), value = TRUE),
    covariates = character(0), analyses = analyses, verbose = FALSE,
    improvement_direction = "higher", stable_threshold = 0)
  if (!is.null(declared)) args$outcomes <- declared
  withr::local_preserve_seed()
  set.seed(9432)
  MIRA:::.mira_eval_long(do.call(mira_info_long, args))
}

test_that("canonical fixtures preserve the caller RNG kind and seed presence", {
  original_kind <- RNGkind()
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) original_seed <- get(".Random.seed", envir = .GlobalEnv)
  withr::defer({
    do.call(RNGkind, as.list(original_kind))
    if (had_seed) assign(".Random.seed", original_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  })
  RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
  set.seed(20261026)
  seed <- .Random.seed
  kind <- RNGkind()
  data <- canonical_contract_data()
  expect_identical(.Random.seed, seed)
  expect_identical(RNGkind(), kind)
  expect_identical(canonical_contract_data(), data)
  captured <- canonical_contract_fit(data)
  expect_null(captured$error)
  expect_identical(.Random.seed, seed)
  expect_identical(RNGkind(), kind)
  invalid <- data
  invalid$patient[2L] <- invalid$patient[1L]
  expect_match(canonical_contract_fit(invalid)$error, "[Dd]uplicat")
  expect_identical(.Random.seed, seed)

  rm(".Random.seed", envir = .GlobalEnv)
  expect_identical(canonical_contract_data(), data)
  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
  expect_null(canonical_contract_fit(data)$error)
  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
  expect_match(canonical_contract_fit(invalid)$error, "[Dd]uplicat")
  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
  expect_identical(RNGkind(), kind)
})

test_that("canonical fields are unique while distinct equal and unavailable metrics survive", {
  data <- data.frame(patient = 1:4, marker_t0 = rep(7, 4), marker_t1 = rep(7, 4))
  for (name in c("Vital_Index", "mean_sd_change", "MiXeD")) {
    captured <- canonical_contract_fit(data, name)
    expect_null(captured$error)
    x <- captured$value
    expect_false("outcomes" %in% names(x))
    expected <- paste0(name, "_", c("mean", "sd", "variance", "se", "ci_lower", "ci_upper",
      "median", "q1", "q3", "iqr", "min", "max", "cv_percent"))
    expect_identical(tail(names(x$descriptives), length(expected)), expected)
    expect_false(any(c("mean", "sd", "variance", "se") %in% names(x$descriptives)))
    expect_identical(x$descriptives[[paste0(name, "_sd")]], c(0, 0))
    expect_identical(x$descriptives[[paste0(name, "_variance")]], c(0, 0))
    expect_identical(x$descriptives[[paste0(name, "_se")]], c(0, 0))
    expect_true(all(c(paste0(name, "_paired_t_p"), paste0(name, "_paired_t_p_adj"),
      paste0(name, "_wilcoxon_p"), paste0(name, "_wilcoxon_p_adj")) %in% names(x$change)))
    expect_false(paste0(name, "_change") %in% names(x$change))
    expect_true(paste0(name, "_mean_change") %in% names(x$change))
    expect_true(paste0(name, "_absolute_change") %in% names(x$trajectories))
    expect_true(all(c("direction", "clinical_direction") %in% names(x$trajectories)))
    expect_true(all(paste0(name, "_", c("ICC", "ICC_model", "ICC_anova_complete_profiles")) %in%
      names(x$variability)))
    for (path in attr(x, "mira_metric_paths")) {
      table <- MIRA:::.mira_result_path_long(x, path)
      expect_false(anyDuplicated(names(table)) > 0L, info = paste(path, collapse = "$"))
    }
  }
})

test_that("all registered metric families preserve values and semantic order across declared names", {
  data <- canonical_contract_data()
  generic <- canonical_contract_fit(data, analyses = "all")
  expect_null(generic$error)
  original <- generic$value
  for (declared in c("Measure_A", "sd_mean_MiX")) {
    captured <- canonical_contract_fit(data, declared, analyses = "all")
    expect_null(captured$error)
    x <- captured$value
    expect_identical(captured$warnings, generic$warnings)
    expect_identical(names(x), names(original))
    expect_identical(attr(x, "mira_metric_paths"), attr(original, "mira_metric_paths"))
    for (path in attr(original, "mira_metric_paths")) {
      before <- MIRA:::.mira_result_path_long(original, path)
      after <- MIRA:::.mira_result_path_long(x, path)
      old_map <- attr(before, "mira_metric_names")
      new_map <- attr(after, "mira_metric_names")
      label <- paste(path, collapse = "$")
      expect_identical(names(new_map), names(old_map), info = label)
      expect_identical(unname(old_map), names(old_map), info = label)
      expected_metrics <- if (length(old_map)) paste0(declared, "_", names(old_map)) else character(0)
      expect_identical(unname(new_map), expected_metrics, info = label)
      expected_names <- names(before)
      expected_names[match(unname(old_map), expected_names)] <- unname(new_map)
      expect_identical(names(after), expected_names, info = label)
      for (base in names(old_map)) {
        expect_identical(after[[new_map[[base]]]], before[[old_map[[base]]]],
          info = paste(label, base))
      }
      if (is.data.frame(before)) {
        expect_identical(rownames(after), rownames(before), info = label)
        for (id in intersect(c("patient", "time", "time_index", "time_label", "from", "to",
          "from_label", "to_label", "arm", "arm_a", "arm_b", "n", "n_a", "n_b"), names(before))) {
          expect_identical(after[[id]], before[[id]], info = paste(label, id))
        }
      }
    }
    expect_identical(stats::coef(x$model$fitted_model), stats::coef(original$model$fitted_model))
    expect_identical(attributes(x$model$summary), attributes(original$model$summary))
    expect_identical(names(x$advanced_tests$emmeans$time$tidy),
      names(original$advanced_tests$emmeans$time$tidy))
    expect_identical(x$plots$mean_ci$data, original$plots$mean_ci$data)
  }
})

test_that("empty change tables have the same canonical columns as populated tables", {
  populated <- canonical_contract_fit(canonical_contract_data(), "Total_Mean")$value
  empty_data <- data.frame(patient = 1:4, marker_t0 = c(1, 2, NA, NA), marker_t1 = c(NA, NA, 3, 4))
  empty <- canonical_contract_fit(empty_data, "Total_Mean")$value
  expect_identical(nrow(empty$change), 0L)
  expect_identical(names(empty$change), names(populated$change))
  expect_identical(vapply(empty$change, typeof, character(1)),
    vapply(populated$change, typeof, character(1)))
  expect_identical(attr(empty$change, "mira_metric_names"), attr(populated$change, "mira_metric_names"))
  expect_null(empty$model$fitted_model)
  expect_null(empty$correlations$pearson)
  expect_length(empty$plots, 0L)
  expect_null(empty$arm_analysis$descriptives)
})

test_that("only proven nested aliases disappear and references resolve to retained results", {
  x <- canonical_contract_fit(canonical_contract_data(), "Metric_X", "all")$value
  expect_false(any(c("object", "summary", "tidy") %in% names(x$advanced_tests$emmeans)))
  expect_false(any(c("object", "summary") %in% names(x$advanced_tests$friedman)))
  expect_false("n" %in% names(x$advanced_tests$rm_anova))
  expect_false(any(c("object", "tidy") %in% names(x$advanced_models$random_slope)))
  expect_false(any(c("object", "summary", "tidy") %in% names(x$robustness$club_sandwich)))
  for (name in c("compound_symmetry", "ar1")) {
    expect_false("object" %in% names(x$advanced_models$nlme[[name]]))
  }
  for (name in c("independence", "exchangeable", "ar1")) {
    gee <- x$advanced_models$gee[[name]]
    expect_false(any(c("object", "robust_standard_errors") %in% names(gee)))
    expect_false("Std.err" %in% names(gee$coefficients))
    if (isTRUE(gee$performed)) expect_true("Metric_X_robust_se" %in% names(gee$coefficients))
  }
  for (name in setdiff(names(x$multiplicity), c("primary_method", "note"))) {
    family <- x$multiplicity[[name]]
    expect_false("table" %in% names(family))
    if (!is.null(family$object_path) && !is.na(family$object_path)) {
      expect_true(is.data.frame(MIRA:::.mira_result_path_long(x, family$object_path)))
    }
  }
  expect_s3_class(x$advanced_tests$friedman$test, "htest")
  expect_true(is.data.frame(x$advanced_tests$friedman$tidy))
  expect_true(all(c("Metric_X_ICC", "Metric_X_ICC_model", "Metric_X_ICC_anova_complete_profiles") %in%
    names(x$variability)))
})

test_that("canonical naming rejects collisions explicitly without overwriting distinct fields", {
  object <- list(descriptives = data.frame(mean = 1, Named_mean = 2))
  expect_error(MIRA:::.mira_finalize_output_long(object, "Named"), "collide.*descriptives")
  expect_identical(object$descriptives$mean, 1)
  expect_identical(object$descriptives$Named_mean, 2)
})

test_that("column ordering preserves automatic and explicit row-name semantics", {
  for (row_names in list(NULL, c("subject-b", "subject-a"))) {
    before <- data.frame(mean = c(2, 3), sd = c(1, 1), row.names = row_names)
    attr(before, "note") <- "scientific context"
    after <- MIRA:::.mira_order_fields_long(before, c("sd", "mean"))
    expect_identical(.row_names_info(after, 0L), .row_names_info(before, 0L))
    expect_identical(rownames(as.matrix(after)), rownames(as.matrix(before)))
    expect_identical(attr(after, "note"), attr(before, "note"))
  }
})

test_that("analysis, detection, summary and native metadata retain their scoped contracts", {
  data <- canonical_contract_data()
  x <- canonical_contract_fit(data, "Measure_Q")$value
  inspect <- mira_info_long(data, id = "patient", analyses = "none", inspect_only = TRUE, verbose = FALSE)
  detected <- mira_detect_long(data, id = "patient", verbose = FALSE)
  data$other_t0 <- data$marker_t0 * 2
  data$other_t1 <- data$marker_t1 * 2
  multi <- mira_info_long(data, id = "patient", analyses = "none", verbose = FALSE)
  for (item in list(x, inspect, detected, multi, summary(multi))) {
    expect_false("version" %in% names(item))
    expect_null(attr(item, "version", exact = TRUE))
  }
  expect_identical(names(multi$outcomes), c("marker", "other"))
  for (item in multi$outcomes) expect_false("version" %in% names(item))
  native <- structure(data.frame(value = 1), version = "external-library-metadata")
  descriptive <- structure(data.frame(sd = 2, mean = 1), note = "scientific note", source = "study")
  manual <- list(descriptives = descriptive, model = list(summary = native))
  finalized <- MIRA:::.mira_finalize_output_long(manual, "Signal")
  expect_identical(finalized$model$summary, native)
  expect_identical(attr(finalized$descriptives, "note"), "scientific note")
  expect_identical(attr(finalized$descriptives, "source"), "study")
  expect_identical(names(finalized$descriptives), c("Signal_mean", "Signal_sd"))
})
