compilation_test_data_long <- function() {
  mira_data_long(
    data = data.frame(
      patient = paste0("P", seq_len(4L)),
      arm = rep(c("Control", "Treatment"), 2L),
      score_t0 = c(1, 2, 3, 4),
      score_t1 = c(2, 3, 3, 5)
    ),
    time_value = c(0, 1),
    outcome = "generic",
    likelihood = "gaussian",
    meaningful_change = 1,
    meaningful_change_sd = 0.2,
    direction = "higher",
    reference_arm = "Control",
    covariates = character(0)
  )
}

test_that("packaged Stan source compiles outside its installation directory", {
  calls <- list()
  testthat::local_mocked_bindings(
    cmdstan_model = function(stan_file, quiet, dir) {
      calls[[length(calls) + 1L]] <<- list(source = stan_file, dir = dir)
      stop("compilation intercepted before sampling", call. = FALSE)
    },
    .package = "cmdstanr"
  )

  stan_data <- compilation_test_data_long()
  expect_error(
    mira_fit_long(stan_data, verbose = FALSE),
    "compilation intercepted before sampling"
  )
  expect_length(calls, 1L)
  call <- calls[[1L]]
  expect_identical(
    normalizePath(call$source, winslash = "/"),
    normalizePath(system.file("stan", "mira_longitudinal.stan", package = "MIRA"),
                  winslash = "/")
  )
  compile_dir <- normalizePath(call$dir, winslash = "/")
  source_dir <- normalizePath(dirname(call$source), winslash = "/")
  session_dir <- normalizePath(tempdir(), winslash = "/")
  expect_false(identical(compile_dir, source_dir))
  expect_true(startsWith(compile_dir, paste0(session_dir, "/")))

  # Verify the directory can actually hold a compiler output, rather than
  # merely checking that a path was passed to the mocked compiler.
  probe <- file.path(call$dir, "mira-write-probe")
  expect_true(file.create(probe))
  unlink(probe)
})

test_that("custom Stan source is preserved and changed content gets a fresh cache", {
  calls <- list()
  testthat::local_mocked_bindings(
    cmdstan_model = function(stan_file, quiet, dir) {
      calls[[length(calls) + 1L]] <<- list(source = stan_file, dir = dir)
      stop("compilation intercepted before sampling", call. = FALSE)
    },
    .package = "cmdstanr"
  )

  custom_dir <- tempfile("mira-custom-source-")
  dir.create(custom_dir)
  on.exit(unlink(custom_dir, recursive = TRUE), add = TRUE)
  custom_source <- file.path(custom_dir, "custom_model.stan")
  original <- c("parameters { real x; }", "model { x ~ normal(0, 1); }")
  writeLines(original, custom_source)
  stan_data <- compilation_test_data_long()

  for (i in seq_len(2L)) {
    expect_error(
      mira_fit_long(stan_data, stan_file = custom_source, verbose = FALSE),
      "compilation intercepted before sampling"
    )
  }
  expect_identical(readLines(custom_source), original)
  expect_identical(calls[[1L]]$source, custom_source)
  expect_identical(calls[[1L]]$dir, calls[[2L]]$dir)
  expect_false(identical(
    normalizePath(calls[[1L]]$dir, winslash = "/"),
    normalizePath(custom_dir, winslash = "/")
  ))

  updated <- c(original, "// changed source")
  writeLines(updated, custom_source)
  expect_error(
    mira_fit_long(stan_data, stan_file = custom_source, verbose = FALSE),
    "compilation intercepted before sampling"
  )
  expect_length(calls, 3L)
  expect_identical(calls[[3L]]$source, custom_source)
  expect_false(identical(calls[[1L]]$dir, calls[[3L]]$dir))
  expect_identical(readLines(custom_source), updated)
})
