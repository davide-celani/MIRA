test_that("section routing filters outcome and context trees", {
  routes <- c(
    settings = "methods", time_vars = "methods", time_labels = "methods",
    overview = "overview", descriptives = "descriptives",
    missing = "missingness", change = "change", arm_analysis = "arms",
    correlations = "correlations", variability = "variability",
    model = "models", advanced_tests = "models", advanced_models = "models",
    robustness = "robustness", effect_sizes = "robustness",
    multiplicity = "robustness", sensitivity = "robustness",
    outliers = "outliers", trajectories = "trajectories", plots = "figures",
    diagnostics = "diagnostics"
  )
  expect_identical(vapply(names(routes), .mira_report_section_long, character(1L)), routes)
  expect_identical(.mira_report_section_long("frequencies"), "frequencies")
  expect_true(.mira_report_selected_long(report_test_cfg("all"), "models"))
  expect_true(.mira_report_selected_long(report_test_cfg(c("overview", "models")), "models"))
  expect_false(.mira_report_selected_long(report_test_cfg("overview"), "models"))

  outcome <- report_test_mira()
  fitted <- stats::lm(y ~ x, data = data.frame(x = 1:4, y = c(1, 2, 4, 8)))
  outcome$model <- list(fitted_model = fitted)
  outcome$advanced_tests <- list(
    friedman = list(
      performed = TRUE,
      tidy = data.frame(statistic = 3.5, p.value = 0.04, check.names = FALSE)
    )
  )
  nodes <- .mira_report_outcome_tree_long(outcome, report_test_cfg(c("overview", "models")))
  expect_setequal(names(nodes), c("overview", "model", "advanced_tests"))
  expect_null(nodes$descriptives)
  expect_null(nodes$correlations)

  frequency_outcome <- structure(
    list(frequencies = data.frame(level = "yes", n = 4L)),
    class = c("mira_info_long", "list")
  )
  expect_named(
    .mira_report_outcome_tree_long(frequency_outcome, report_test_cfg("frequencies")),
    "frequencies"
  )

  detected <- report_test_detect()
  expect_named(.mira_report_context_tree_long(detected, report_test_cfg("methods")), "config")
  expect_setequal(
    names(.mira_report_context_tree_long(detected, report_test_cfg("data_quality"))),
    c("data_overview", "detected_variables")
  )
})

test_that("abstract availability is independent of display-section filtering", {
  output <- capture.output(
    .mira_report_render_abstract_long(
      list(score = report_test_mira()), report_test_cfg("abstract")
    )
  )
  expect_true(any(grepl("^# Abstract$", output)))
  expect_true(any(grepl("score", output, fixed = TRUE)))
})

test_that("pruning classifies reportable leaves and rejects machinery", {
  cfg <- report_test_cfg()
  tabular <- list(
    data_frame = data.frame(x = 1:2),
    matrix = matrix(1:4, nrow = 2L),
    table = table(c("a", "b", "a"))
  )
  for (name in names(tabular)) {
    node <- .mira_report_prune_long(tabular[[name]], name, paste0("outcome$", name), cfg)
    expect_identical(node$kind, "table", info = name)
    expect_identical(node$path, paste0("outcome$", name), info = name)
  }

  plot <- structure(list(), class = "ggplot")
  expect_identical(
    .mira_report_prune_long(plot, "plot", "outcome$plots$plot", cfg)$kind, "plot"
  )
  fitted <- stats::lm(y ~ x, data = data.frame(x = 1:3, y = c(1, 3, 6)))
  expect_identical(
    .mira_report_prune_long(fitted, "fit", "outcome$model$fit", cfg)$kind, "model"
  )
  expect_identical(
    .mira_report_prune_long(y ~ x, "formula", "outcome$formula", cfg)$kind, "formula"
  )
  expect_identical(
    .mira_report_prune_long(quote(mean(x)), "call_expression", "outcome$call_expression", cfg)$kind,
    "formula"
  )
  expect_identical(
    .mira_report_prune_long(expression(x + 1), "expression", "outcome$expression", cfg)$kind,
    "formula"
  )
  expect_identical(
    .mira_report_prune_long(as.name("x"), "symbol", "outcome$symbol", cfg)$kind, "formula"
  )
  expect_identical(
    .mira_report_prune_long(c(alpha = 1), "estimate", "outcome$estimate", cfg)$kind, "value"
  )

  rejected <- list(
    null = NULL,
    function_value = identity,
    environment = new.env(parent = emptyenv()),
    outcome_error = structure(
      list(error = "controlled"), class = c("mira_info_error_long", "list")
    ),
    all_missing = c(NA_real_, NA_real_)
  )
  for (name in names(rejected)) {
    expect_null(
      .mira_report_prune_long(rejected[[name]], name, paste0("outcome$", name), cfg),
      info = name
    )
  }

  hidden_keys <- c(
    "call", "outcomes", "long_data", "reason_skipped", "performed", "enabled",
    "package", "status", "object_path", "plot_error", "disabled", "failed_outcomes"
  )
  for (key in hidden_keys) {
    expect_null(
      .mira_report_prune_long("internal", key, paste0("outcome$", key), cfg), info = key
    )
  }
  expect_null(.mira_report_prune_long(
    data.frame(x = 1), "data", "outcome$advanced_models$random_slope$data", cfg
  ))
  expect_null(.mira_report_prune_long(
    TRUE, "arm_analysis", "outcome$overview$arm_analysis", cfg
  ))

  table_node <- .mira_report_prune_long(data.frame(x = 1), "table", "outcome$table", cfg)
  value_node <- .mira_report_prune_long(1, "value", "outcome$value", cfg)
  expect_true(.mira_report_node_has_result_long(table_node))
  expect_false(.mira_report_node_has_result_long(value_node))
  expect_true(.mira_report_node_has_result_long(list(
    kind = "group", children = list(table_node)
  )))
})

test_that("performed FALSE is authoritative while performed TRUE preserves empty results", {
  cfg <- report_test_cfg()
  stale_fit <- stats::lm(y ~ x, data = data.frame(x = 1:3, y = 1:3))
  failed <- list(
    performed = FALSE,
    summary = data.frame(term = "stale", estimate = 99),
    model = stale_fit,
    error = "estimation failed",
    reason_skipped = "not estimable"
  )
  expect_null(.mira_report_prune_long(
    failed, "failed_module", "outcome$advanced_models$failed_module", cfg
  ))

  completed <- list(
    performed = TRUE,
    tidy = data.frame(term = character(), estimate = numeric())
  )
  node <- .mira_report_prune_long(
    completed, "completed_module", "outcome$advanced_tests$completed_module", cfg
  )
  expect_identical(node$kind, "group")
  expect_identical(node$children$tidy$kind, "table")
  expect_identical(nrow(node$children$tidy$value), 0L)
})

test_that("equivalent aliases and duplicate model objects appear only once", {
  cfg <- report_test_cfg()
  shared <- data.frame(term = c("time", "arm"), estimate = c(1.5, -0.25))
  different <- data.frame(term = "interaction", estimate = 0.75)
  module <- list(
    performed = TRUE,
    summary = shared,
    tidy = shared,
    fixed_effects = different,
    coefficient_tests = different
  )
  node <- .mira_report_prune_long(module, "module", "outcome$advanced_tests$module", cfg)
  expect_identical(names(node$children), c("summary", "fixed_effects"))
  expect_identical(node$children$summary$value, shared)
  expect_identical(node$children$fixed_effects$value, different)

  distinct <- module
  distinct$tidy <- transform(shared, estimate = estimate + 0.01)
  distinct_node <- .mira_report_prune_long(
    distinct, "module", "outcome$advanced_tests$module", cfg
  )
  expect_true(all(c("summary", "tidy", "fixed_effects") %in% names(distinct_node$children)))

  fitted <- stats::lm(y ~ x, data = data.frame(x = 1:4, y = c(1, 2, 5, 9)))
  model_node <- .mira_report_prune_long(
    list(performed = TRUE, model = fitted, object = fitted),
    "fit", "outcome$advanced_models$fit", cfg
  )
  expect_named(model_node$children, "model")
})

test_that("model comparison requires a nonblank fitted-object path", {
  cfg <- report_test_cfg()
  comparison <- data.frame(
    model = c("fitted", "empty", "missing", "whitespace"),
    object_path = c("$advanced_models$fitted$model", "", NA_character_, "  \t"),
    converged = c(TRUE, FALSE, FALSE, FALSE), stringsAsFactors = FALSE
  )
  node <- .mira_report_prune_long(
    comparison, "model_comparison", "outcome$advanced_models$model_comparison", cfg
  )
  expect_identical(node$value, comparison[1L, , drop = FALSE])

  comparison$object_path <- c("", NA_character_, " ", "\t")
  expect_null(.mira_report_prune_long(
    comparison, "model_comparison", "outcome$advanced_models$model_comparison", cfg
  ))

  descriptive_only <- comparison[c("model", "converged")]
  expect_identical(
    .mira_report_prune_long(
      descriptive_only, "model_comparison", "outcome$advanced_models$model_comparison", cfg
    )$value,
    descriptive_only
  )
})

test_that("NLME keeps only genuinely fitted correlation structures", {
  cfg <- report_test_cfg()
  fitted <- stats::lm(y ~ x, data = data.frame(x = 1:4, y = c(1, 2, 4, 7)))
  nlme <- list(
    performed = TRUE,
    package = "nlme",
    compound_symmetry = list(
      performed = TRUE,
      model = fitted,
      object = fitted,
      summary = data.frame(term = "x", estimate = stats::coef(fitted)[[2L]])
    ),
    ar1 = list(
      performed = FALSE, model = NULL, object = NULL,
      error = "controlled AR(1) failure", reason_skipped = "not estimable"
    ),
    declared_without_fit = list(performed = TRUE, model = NULL, object = NULL),
    comparison = data.frame(
      model = c("compound_symmetry", "ar1", "declared_without_fit"),
      AIC = c(10, NA, NA), stringsAsFactors = FALSE
    )
  )
  node <- .mira_report_prune_long(nlme, "nlme", "outcome$advanced_models$nlme", cfg)
  expect_true("compound_symmetry" %in% names(node$children))
  expect_false("ar1" %in% names(node$children))
  expect_false("declared_without_fit" %in% names(node$children))
  expect_identical(node$children$comparison$value$model, "compound_symmetry")
  expect_named(node$children$compound_symmetry$children, c("model", "summary"))

  unavailable <- nlme
  unavailable$performed <- FALSE
  expect_null(.mira_report_prune_long(
    unavailable, "nlme", "outcome$advanced_models$nlme", cfg
  ))
})

test_that("arm analysis requires an explicit enabled TRUE flag", {
  cfg <- report_test_cfg()
  arm_results <- list(
    enabled = TRUE,
    levels = c("control", "active"),
    descriptives = data.frame(
      arm = c("control", "active"), mean = c(10, 12), stringsAsFactors = FALSE
    )
  )
  enabled <- .mira_report_prune_long(
    arm_results, "arm_analysis", "outcome$arm_analysis", cfg
  )
  expect_identical(enabled$kind, "group")
  expect_true(all(c("levels", "descriptives") %in% names(enabled$children)))
  expect_false("enabled" %in% names(enabled$children))

  disabled <- arm_results
  disabled$enabled <- FALSE
  expect_null(.mira_report_prune_long(
    disabled, "arm_analysis", "outcome$arm_analysis", cfg
  ))
  expect_null(.mira_report_prune_long(
    arm_results[names(arm_results) != "enabled"],
    "arm_analysis", "outcome$arm_analysis", cfg
  ))
})

test_that("all-missing correlations are omitted without hiding estimable matrices", {
  all_missing <- matrix(
    NA_real_, 2L, 2L, dimnames = list(c("t0", "t1"), c("t0", "t1"))
  )
  partial <- all_missing
  diag(partial) <- 1
  outcome <- structure(
    list(correlations = list(
      pearson = all_missing,
      spearman = partial,
      pairwise_n = matrix(c(6L, 5L, 5L, 6L), 2L)
    )),
    class = c("mira_info_long", "list")
  )
  node <- .mira_report_outcome_tree_long(outcome, report_test_cfg("correlations"))$correlations
  expect_false("pearson" %in% names(node$children))
  expect_true(all(c("spearman", "pairwise_n") %in% names(node$children)))
  expect_identical(node$children$spearman$value, partial)

  outcome$correlations$spearman <- all_missing
  outcome$correlations$pairwise_n <- NULL
  expect_null(
    .mira_report_outcome_tree_long(outcome, report_test_cfg("correlations"))$correlations
  )
})

test_that("outlier trees distinguish insufficient samples from zero outliers", {
  empty_flags <- data.frame(patient = character(), value = numeric())
  observed_flag <- data.frame(patient = "P1", value = 100)
  outcome <- structure(
    list(
      descriptives = data.frame(
        time = c("small", "adequate", "observed"), n = c(3L, 4L, 2L),
        stringsAsFactors = FALSE
      ),
      change = data.frame(
        from = c("t0", "t0"), to = c("t1", "t2"), n = c(3L, 4L),
        stringsAsFactors = FALSE
      ),
      outliers = list(
        note = "IQR diagnostic flags",
        by_time = list(
          small = empty_flags, adequate = empty_flags, observed = observed_flag
        ),
        change = list(t0_to_t1 = empty_flags, t0_to_t2 = empty_flags)
      )
    ),
    class = c("mira_info_long", "list")
  )
  node <- .mira_report_outcome_tree_long(outcome, report_test_cfg("outliers"))$outliers
  expect_false("small" %in% names(node$children$by_time$children))
  expect_true(all(c("adequate", "observed") %in% names(node$children$by_time$children)))
  expect_identical(nrow(node$children$by_time$children$adequate$value), 0L)
  expect_identical(nrow(node$children$by_time$children$observed$value), 1L)
  expect_false("t0_to_t1" %in% names(node$children$change$children))
  expect_named(node$children$change$children, "t0_to_t2")

  rendered_empty <- capture.output(.mira_report_render_node_long(
    node$children$by_time$children$adequate, 3L, report_test_cfg()
  ))
  expect_true(any(grepl("Outliers detected:** 0", rendered_empty, fixed = TRUE)))
})

test_that("max_depth is enforced at the recursive boundary", {
  nested <- list(level_one = list(result = data.frame(estimate = 1)))
  shallow <- .mira_report_prune_long(
    nested, "root", "outcome$root", report_test_cfg(max_depth = 1L)
  )
  expect_null(shallow)

  boundary <- .mira_report_prune_long(
    nested, "root", "outcome$root", report_test_cfg(max_depth = 2L)
  )
  expect_identical(boundary$children$level_one$children$result$kind, "table")
  expect_identical(
    boundary$children$level_one$children$result$path,
    "outcome$root$level_one$result"
  )
})

test_that("an empty mira_info_multi_long has no synthetic outcome", {
  empty_multi <- structure(
    list(
      version = "4.0.0",
      config = list(id = "patient", outcomes = character()),
      data_overview = data.frame(records = 4L, variables = 3L),
      detected_variables = list(outcomes = character()),
      outcomes = list(),
      diagnostics = list(warnings = "configuration-only warning")
    ),
    class = c("mira_info_multi_long", "mira_info_long", "list")
  )
  expect_identical(.mira_report_extract_outcomes_long(empty_multi), list())

  context <- .mira_report_context_tree_long(empty_multi, report_test_cfg())
  expect_true(all(c("config", "data_overview", "diagnostics") %in% names(context)))
  archived <- .mira_report_collect_tables_long(empty_multi, report_test_cfg())
  paths <- vapply(archived, `[[`, character(1L), "path")
  expect_false(any(startsWith(paths, "outcome$")))

  rendered <- capture.output(.mira_report_knit_long(list(
    result = empty_multi, extra_objects = list(), report = report_test_cfg()
  )))
  expect_false(any(grepl("^# Outcome:", rendered)))
})

test_that("mira_info_error_long outcomes are surfaced only in diagnostics", {
  multi <- report_test_multi(include_error = TRUE)
  outcomes <- .mira_report_extract_outcomes_long(multi)
  expect_named(outcomes, c("score", "fatigue", "failed"))
  expect_identical(.mira_report_outcome_tree_long(outcomes$failed, report_test_cfg()), list())
  expect_null(.mira_report_prune_long(
    outcomes$failed, "failed", "result$outcomes$failed", report_test_cfg()
  ))

  context <- .mira_report_context_tree_long(multi, report_test_cfg("diagnostics"))
  error_node <- context$diagnostics$children$outcome_errors$children$failed
  expect_identical(error_node$kind, "value")
  expect_identical(error_node$value, "controlled outcome failure")
  expect_identical(error_node$path, "result$outcomes$failed$error")

  without_diagnostics <- .mira_report_context_tree_long(multi, report_test_cfg("overview"))
  expect_null(without_diagnostics$diagnostics)
})
