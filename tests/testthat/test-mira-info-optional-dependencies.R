test_that("missing optional engines yield structured module skips", {
  dummy_model <- structure(list(), class = "mira_dummy_model")
  formula <- stats::as.formula("value ~ time_factor")

  modules <- testthat::with_mocked_bindings(
    code = list(
      emmeans = MIRA:::.mira_run_emmeans(
        dummy_model,
        arm_enabled = FALSE,
        primary_method = "holm",
        confidence_level = 0.95
      ),
      rm_anova = MIRA:::.mira_run_rm_anova(
        data.frame(),
        arm_enabled = FALSE,
        covariates = character(0),
        categorical_covariates = character(0)
      ),
      gee = MIRA:::.mira_run_gee(data.frame(), formula, 0.95),
      random_slope = MIRA:::.mira_run_random_slope(
        data.frame(),
        arm_enabled = FALSE,
        covariates = character(0),
        mixed_intercept_model = NULL
      ),
      nlme = MIRA:::.mira_run_nlme(data.frame(), formula),
      robust = MIRA:::.mira_run_robust_inference(dummy_model, 0.95),
      r2 = MIRA:::.mira_run_r2(dummy_model)
    ),
    requireNamespace = function(package, quietly = TRUE) FALSE,
    .package = "base"
  )

  expect_false(modules$emmeans$performed)
  expect_identical(modules$emmeans$package, "emmeans")
  expect_match(modules$emmeans$reason_skipped, "not installed")

  expect_false(modules$rm_anova$performed)
  expect_identical(modules$rm_anova$package, "afex")
  expect_match(modules$rm_anova$reason_skipped, "not installed")

  expect_false(modules$gee$performed)
  expect_identical(modules$gee$package, "geepack")
  expect_match(modules$gee$reason_skipped, "not installed")
  for (structure in c("independence", "exchangeable", "ar1")) {
    expect_false(modules$gee[[structure]]$performed)
    expect_identical(modules$gee[[structure]]$package, "geepack")
  }

  expect_false(modules$random_slope$performed)
  expect_identical(modules$random_slope$package, "lme4")
  expect_match(modules$random_slope$reason_skipped, "not installed")

  expect_false(modules$nlme$performed)
  expect_identical(modules$nlme$package, "nlme")
  expect_match(modules$nlme$reason_skipped, "not installed")
  expect_false(modules$nlme$compound_symmetry$performed)
  expect_false(modules$nlme$ar1$performed)

  expect_false(modules$robust$performed)
  expect_identical(modules$robust$package, "clubSandwich")
  expect_match(modules$robust$reason_skipped, "not installed")

  expect_false(modules$r2$performed)
  expect_identical(modules$r2$package, "performance")
  expect_match(modules$r2$reason_skipped, "not installed")
})

test_that("missing model and plotting packages do not destroy core public results", {
  data <- data.frame(
    patient_id = seq_len(8),
    score_t0 = c(1, 3, 2, 5, 4, 7, 6, 8),
    score_t1 = c(2, 5, 4, 7, 7, 9, 8, 11),
    score_t2 = c(4, 6, 5, 9, 8, 11, 10, 13)
  )

  result <- testthat::with_mocked_bindings(
    code = mira_info(
      data,
      id = "patient_id",
      analyses = c("model", "plots"),
      covariates = character(0),
      verbose = FALSE
    ),
    requireNamespace = function(package, quietly = TRUE) FALSE,
    .package = "base"
  )

  expect_s3_class(result, "mira_info")
  expect_equal(result$descriptives$mean, unname(vapply(data[-1L], mean, numeric(1L))))
  expect_equal(result$change$mean_change, c(2.125, 3.75, 1.625))
  expect_true(result$advanced_tests$friedman$performed)

  expect_null(result$model$fitted_model)
  expect_match(result$model$error, "lme4.*not installed")
  expect_false(result$advanced_tests$emmeans$performed)
  expect_false(result$advanced_tests$rm_anova$performed)
  expect_false(result$advanced_models$random_slope$performed)
  expect_false(result$advanced_models$gee$performed)
  expect_false(result$advanced_models$nlme$performed)
  expect_false(result$robustness$club_sandwich$performed)

  expect_length(result$plots, 0L)
  expect_match(result$plot_error, "ggplot2.*not installed")
  expect_identical(result$config$analyses$model, TRUE)
  expect_identical(result$config$analyses$plots, TRUE)
})
