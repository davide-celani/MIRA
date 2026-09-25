test_that("the report omits model sections without fitted results", {
  cfg <- list(sections = "all", max_depth = 50L)
  gee <- list(
    performed = FALSE, package = "geepack", model = NULL,
    error = "GEE estimation failed", warnings = "numerical warning",
    reason_skipped = "GEE estimation failed for this working correlation."
  )
  node <- .mira_report_prune_long(gee, "independence",
                             "outcome$advanced_models$gee$independence", cfg)
  expect_null(node)

  primary <- list(fitted_model = NULL, error = "lmer failed",
                  warnings = "singular fit", covariates_requested = "age",
                  covariates_used = character(), covariates_skipped = "age")
  outcome <- structure(list(model = primary), class = c("mira_info_long", "list"))
  nodes <- .mira_report_outcome_tree_long(outcome, cfg)
  expect_null(nodes$model)
  expect_identical(nodes$diagnostics$children$model_error$value, primary$error)
  expect_identical(nodes$diagnostics$children$model_warnings$value,
                   primary$warnings)
})

test_that("model comparison omits unavailable fits", {
  cfg <- list(sections = "all", max_depth = 50L)
  comparison <- data.frame(
    model = c("gee_independence", "gee_ar1"),
    object_path = c("$advanced_models$gee$independence$model", NA_character_),
    converged = c(TRUE, FALSE)
  )
  node <- .mira_report_prune_long(comparison, "model_comparison",
                             "outcome$advanced_models$model_comparison", cfg)
  expect_equal(node$value, comparison[1L, , drop = FALSE])
  comparison$object_path <- NA_character_
  expect_null(.mira_report_prune_long(comparison, "model_comparison",
                                 "outcome$advanced_models$model_comparison", cfg))
})

test_that("the report shows the analyses selected in mira_info_long configuration", {
  cfg <- list(sections = "all", max_depth = 50L)
  detected <- structure(
    list(config = list(analyses = list(core = TRUE, model = FALSE,
                                         correlations = TRUE))),
    class = c("mira_detect_long", "list")
  )
  context <- .mira_report_context_tree_long(detected, cfg)
  expect_identical(context$config$children$analyses$children$model$value,
                   FALSE)
  expect_identical(context$config$children$analyses$children$correlations$value,
                   TRUE)
})

test_that("GEE summaries are printed as summaries rather than internal lists", {
  cfg <- list(sections = "all", max_depth = 50L)
  summary <- structure(list(coefficients = matrix(1, 1, 1),
                            geese = list(control = list(maxit = 25))),
                       class = "summary.geeglm")
  node <- .mira_report_prune_long(summary, "summary",
                             "outcome$advanced_models$gee$ar1$summary", cfg)
  expect_identical(node$kind, "model")
  expect_identical(node$value, summary)
})

test_that("a report built from mira_info_long preserves the result and exports its tables", {
  data <- data.frame(
    patient = 1:6,
    score_t0 = c(1, 3, 2, 6, 5, 4),
    score_t1 = c(2, 4, 3, 5, 7, 6),
    score_t2 = c(4, 5, 6, 7, 8, 9)
  )
  fit <- mira_info_long(data, id = "patient", outcomes = "score",
                   time_vars = c("score_t0", "score_t1", "score_t2"),
                   analyses = "none", verbose = FALSE)
  output <- tempfile("mira-report-")
  capture.output(report <- mira_report_freq_long(
    fit, output_dir = output, format = "html", render = FALSE, open = FALSE
  ))
  saved <- readRDS(report$payload)$result
  expect_equal(saved$descriptives, fit$descriptives)
  expect_equal(saved$advanced_tests$friedman$tidy,
               fit$advanced_tests$friedman$tidy)
  expect_identical(saved$advanced_tests$emmeans$reason_skipped,
                   fit$advanced_tests$emmeans$reason_skipped)
  nodes <- .mira_report_outcome_tree_long(fit, report$report_config)
  expect_null(nodes$model)
  expect_null(nodes$advanced_models)
  expect_true("friedman" %in% names(nodes$advanced_tests$children))
  manifest <- utils::read.csv(report$results_manifest)
  expect_true("outcome$descriptives" %in% manifest$path)
  expect_true("outcome$advanced_tests$friedman$summary" %in% manifest$path)
  expect_true(all(file.exists(report$result_files)))
})
