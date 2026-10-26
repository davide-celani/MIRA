# Explicit reference choices must be validated even if every arm is unavailable.
# Otherwise config claims a manual choice while silently recording NULL.
empty_reference_data <- function() {
  data.frame(subject_id = paste0("P", 1:4),
    arm = c(NA_character_, "", " ", NA_character_),
    score_t0 = c(3, 7, 4, 11), score_t1 = c(5, 6, 8, 10))
}

empty_reference_call <- function(data, reference, mode) {
  args <- list(data = data, id = "subject_id", arm = "arm",
    reference_arm = reference, covariates = character(0), verbose = FALSE)
  if (mode == "detect") return(do.call(mira_detect_long, args))
  do.call(mira_info_long, c(args,
    list(analyses = "none", inspect_only = mode == "inspect")))
}

test_that("unobserved explicit references fail in every public configuration workflow", {
  data <- empty_reference_data()
  before <- data
  for (mode in c("detect", "inspect", "analyse")) {
    expect_error(empty_reference_call(data, "control", mode),
      "reference_arm='control' is not present in 'arm'", fixed = TRUE, info = mode)
    for (bad in list(character(0), NA_character_, "", c("control", "active"), 1)) {
      expect_error(empty_reference_call(data, bad, mode),
        "reference_arm must identify exactly one non-empty group", info = mode)
    }
  }
  expect_identical(data, before)
})

test_that("an automatic reference remains unavailable with auditable empty arm data", {
  data <- empty_reference_data()
  # The automatic NULL case stays supported: no reference is invented and the
  # missing-arm warning and provenance survive both detection workflows.
  for (mode in c("detect", "inspect")) {
    expect_warning(result <- empty_reference_call(data, NULL, mode),
      "contains no observed groups")
    expect_null(result$config$reference_arm)
    expect_true(result$config$auto_detected$reference_arm)
    expect_false(result$config$specified_manually$reference_arm)
    expect_identical(result$diagnostics$warnings,
      "Treatment-arm variable 'arm' contains no observed groups.")
  }
})

test_that("empty outcome settings identify the invalid public argument", {
  data <- empty_reference_data()[c("subject_id", "score_t0", "score_t1")]
  for (mode in c("detect", "inspect", "analyse")) {
    args <- list(data = data, id = "subject_id", covariates = character(0),
      verbose = FALSE)
    fun <- if (mode == "detect") mira_detect_long else mira_info_long
    if (mode != "detect") args <- c(args,
      list(analyses = "none", inspect_only = mode == "inspect"))
    expect_error(do.call(fun, c(args, list(improvement_direction = character(0)))),
      "improvement_direction must be", info = mode)
    expect_error(do.call(fun, c(args, list(stable_threshold = numeric(0)))),
      "stable_threshold must be", info = mode)
  }
})
