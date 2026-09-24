test_that("outcome registry lists the supported outcomes", {
  registry <- mira_outcome_registry()

  expect_s3_class(registry, "data.frame")
  expect_identical(names(registry), c("outcome", "description"))
  expect_identical(registry$outcome, c("BCVA", "CMT", "generic"))
  expect_identical(nrow(registry), 3L)
  expect_false(anyNA(registry))
  expect_identical(anyDuplicated(registry$outcome), 0L)
  expect_true(all(nzchar(trimws(registry$description))))
})
