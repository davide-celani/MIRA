# Source integration has no access to MIRA helpers in the search path. Supply
# only the three imports declared in NAMESPACE, then source the main file.
consolidation_source_file <- function(filename) {
  candidates <- c(
    testthat::test_path("..", "..", "R", filename),
    file.path("R", filename),
    file.path(find.package("MIRA"), "R", filename),
    file.path(dirname(find.package("MIRA")), "00_pkg_src", "MIRA", "R", filename)
  )
  available <- candidates[file.exists(candidates)]
  if (!length(available)) {
    skip(paste("Original R source is unavailable in this installed-package test:", filename))
  }
  normalizePath(available[[1L]], winslash = "/", mustWork = TRUE)
}

consolidation_source <- function(report = FALSE) {
  imports <- new.env(parent = baseenv())
  imports$coef <- stats::coef
  imports$p.adjust.methods <- stats::p.adjust.methods
  imports$tail <- utils::tail
  isolated <- new.env(parent = imports)
  sys.source(consolidation_source_file("mira_info_long.R"), envir = isolated)
  if (report) {
    sys.source(consolidation_source_file("mira_report_freq_long.R"), envir = isolated)
  }
  isolated
}

consolidation_data <- function() {
  i <- seq_len(24L)
  active <- as.integer(i > 12L)
  baseline <- 10 + i / 3 + 2 * sin(i)
  data.frame(
    patient = sprintf("S%02d", i),
    treatment = factor(ifelse(active == 1L, "active", "control"),
                       levels = c("control", "active")),
    age = 30 + (i * 7) %% 25,
    sex = factor(rep(c("F", "F", "M", "M"), 6L)),
    score_t0 = baseline,
    score_t1 = baseline + 1 + 2 * cos(0.7 * i) + 0.5 * active,
    score_t2 = baseline + 3 + 3 * sin(1.7 * i) + 1.3 * active,
    fatigue_t0 = 35 - baseline / 2 + cos(0.8 * i),
    fatigue_t1 = 34 - baseline / 2 + 2 * sin(0.6 * i) - active,
    fatigue_t2 = 32 - baseline / 2 + 2 * cos(1.3 * i) - 2 * active
  )
}

test_that("the main source owns every extended adapter and has no unresolved internal helper", {
  source_file <- consolidation_source_file("mira_info_long.R")
  expressions <- parse(source_file)
  function_names <- vapply(expressions, function(expr) {
    if (is.call(expr) && identical(expr[[1L]], as.name("<-")) &&
        is.symbol(expr[[2L]]) && is.call(expr[[3L]]) &&
        identical(expr[[3L]][[1L]], as.name("function"))) {
      as.character(expr[[2L]])
    } else ""
  }, character(1L))
  function_names <- function_names[nzchar(function_names)]
  expect_identical(anyDuplicated(function_names), 0L)
  adapters <- c(".mira_protected_method_long", ".mira_extend_sensitivity_long",
    ".mira_extend_inference_long", ".mira_run_kenward_roger_long",
    ".mira_run_mmrm_long", ".mira_run_rank_longitudinal_long")
  expect_true(all(adapters %in% function_names))

  isolated <- consolidation_source()
  functions <- Filter(is.function,
    mget(ls(isolated, all.names = TRUE), envir = isolated, inherits = FALSE))
  referenced <- unique(unlist(lapply(functions, function(fun) {
    all.names(body(fun), functions = TRUE)
  }), use.names = FALSE))
  internal <- grep("^\\.mira_", referenced, value = TRUE)
  expect_setequal(setdiff(internal, ls(isolated, all.names = TRUE)), character())
  for (name in adapters) {
    expect_identical(body(isolated[[name]]), body(getFromNamespace(name, "MIRA")))
    expect_identical(formals(isolated[[name]]), formals(getFromNamespace(name, "MIRA")))
  }
  expect_identical(formals(isolated$mira_info_long), formals(mira_info_long))
})

test_that("sourcing only the main analysis file preserves complete and incomplete core contracts", {
  isolated <- consolidation_source()
  complete <- consolidation_data()[c("patient", "score_t0", "score_t1", "score_t2")]
  incomplete <- complete
  incomplete$score_t1[c(2L, 7L)] <- NA_real_
  incomplete$score_t2[4L] <- NA_real_
  for (data in list(complete, incomplete)) {
    args <- list(data = data, id = "patient", outcomes = "score",
      covariates = character(), analyses = "none", verbose = FALSE)
    sourced <- do.call(isolated$mira_info_long, args)
    packaged <- do.call(mira_info_long, args)
    expect_identical(class(sourced), class(packaged))
    expect_identical(names(sourced), names(packaged))
    expect_identical(attributes(sourced), attributes(packaged))
    for (component in c("descriptives", "missing", "change", "arm_analysis",
        "correlations", "variability", "multiplicity", "sensitivity",
        "outliers", "diagnostics", "plots")) {
      expect_equal(sourced[[component]], packaged[[component]], info = component)
    }
    for (module in list(sourced$advanced_tests$kenward_roger,
        sourced$advanced_tests$rank_longitudinal, sourced$advanced_models$mmrm)) {
      expect_false(module$requested)
      expect_false(module$performed)
    }
  }
})

test_that("the two standalone sources support result and raw-data report modes with real exports", {
  withr::local_preserve_seed()
  isolated <- consolidation_source(report = TRUE)
  data <- consolidation_data()
  data$score_t1[c(2L, 15L)] <- NA_real_
  data$fatigue_t2[c(4L, 19L)] <- NA_real_
  data$age[6L] <- NA_real_
  root <- withr::local_tempdir(pattern = "mira-consolidated-pipeline-")

  # Mode A: an existing single-outcome analysis without treatment/covariates.
  single <- isolated$mira_info_long(
    data[c("patient", "score_t0", "score_t1", "score_t2")],
    id = "patient", outcomes = "score", covariates = character(),
    analyses = "none", verbose = FALSE)
  capture.output(report_a <- isolated$mira_report_freq_long(single,
    output_dir = file.path(root, "existing"), format = c("html", "pdf"),
    render = FALSE, open = FALSE))
  expect_s3_class(report_a, "mira_report_freq_long")
  expect_identical(readRDS(report_a$payload)$result, single)
  expect_named(report_a$expected_files, c("html", "pdf"))

  # Mode B: all analyses, multiple outcomes, treatment, covariates and missing
  # profiles. This must work without ever sourcing the former extra modules.
  capture.output(report_b <- isolated$mira_report_freq_long(data = data,
    id = "patient", outcomes = c("score", "fatigue"), arm = "treatment",
    reference_arm = "control", covariates = c("age", "sex"),
    analyses = "all", verbose = FALSE,
    output_dir = file.path(root, "raw-data"), output_file = "consolidated",
    format = "html", render = FALSE, open = FALSE))
  expect_s3_class(report_b$mira_result, "mira_info_multi_long")
  expect_named(report_b$mira_result$outcomes, c("score", "fatigue"))
  for (report in list(report_a, report_b)) {
    expect_true(all(file.exists(c(report$source, report$payload, report$runtime,
      report$results_manifest, report$result_files))))
    payload <- readRDS(report$payload)
    runtime <- new.env(parent = baseenv())
    expect_no_error(sys.source(report$runtime, envir = runtime))
    manifest <- utils::read.csv(report$results_manifest, check.names = FALSE)
    expect_gt(nrow(manifest), 0L)
    expect_false(anyDuplicated(manifest$file) > 0L)
    expect_true(all(file.exists(file.path(report$output_dir, manifest$file))))
    outcomes <- if (inherits(payload$result, "mira_info_multi_long")) {
      payload$result$outcomes
    } else stats::setNames(list(payload$result), payload$result$outcome)
    for (name in names(outcomes)) {
      fit <- outcomes[[name]]
      expect_s3_class(fit, "mira_info_long")
      expect_false(inherits(fit, "mira_info_error_long"))
      rows <- manifest[manifest$outcome == name &
        manifest$path == "outcome$descriptives", , drop = FALSE]
      expect_equal(nrow(rows), 1L)
      csv <- utils::read.csv(file.path(report$output_dir, rows$file),
                            check.names = FALSE)
      expect_equal(csv[[paste0(name, "_mean")]],
                   fit$descriptives[[paste0(name, "_mean")]], tolerance = 1e-12)
      expect_equal(csv$n, fit$descriptives$n)
      tree <- runtime$.mira_report_outcome_tree_long(fit, payload$report)
      expect_true(all(c("descriptives", "missing", "change") %in% names(tree)))
      if (isTRUE(fit$config$analyses$model)) {
        expect_true(fit$arm_analysis$enabled)
        expect_equal(sum(fit$missing$by_time$unavailable_n), 2L)
        for (method in c("kenward_roger", "rank_longitudinal", "mmrm")) {
          module <- if (method == "mmrm") fit$advanced_models[[method]] else {
            fit$advanced_tests[[method]]
          }
          expect_true(module$requested)
          if (isTRUE(module$performed)) {
            prefix <- paste0("outcome$",
              if (method == "mmrm") "advanced_models$" else "advanced_tests$", method, "$")
            expect_true(any(startsWith(manifest$path[manifest$outcome == name], prefix)))
          } else {
            reason <- c(module$reason_skipped, module$error)
            expect_true(any(!is.na(reason) & nzchar(reason)))
            path <- paste0(if (method == "mmrm") "advanced_models$" else "advanced_tests$", method)
            expect_false(is.null(tree$diagnostics$children[[path]]))
          }
        }
        expect_gt(length(fit$plots), 0L)
        expect_false(is.null(tree$plots))
      }
    }
  }
})
