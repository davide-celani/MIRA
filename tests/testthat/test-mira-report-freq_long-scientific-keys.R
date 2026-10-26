# A visit does not identify a between-arm contrast: with three arms there are
# three different comparisons per visit. Unequal groups and visit-specific
# missingness make row swaps and wrong denominators visible.
report_arm_identity_data <- function() {
  data.frame(patient = sprintf("P%02d", seq_len(12L)),
    arm = factor(rep(c("control", "zeta", "alpha"), c(3L, 4L, 5L)),
      levels = c("zeta", "alpha", "control")),
    score_t0 = c(1, 4, 9, 2, 7, 12, 18, 6, 8, 17, 20, 27),
    score_t1 = c(3, NA, 11, 3, 9, 15, 19, 8, 12, 18, NA, 31))
}

report_arm_identity_oracle <- function(data) {
  pairs <- list(c("control", "zeta"), c("control", "alpha"), c("zeta", "alpha"))
  rows <- list()
  for (visit in c("score_t0", "score_t1")) {
    for (pair in pairs) {
      a <- data[[visit]][as.character(data$arm) == pair[[1L]]]
      b <- data[[visit]][as.character(data$arm) == pair[[2L]]]
      a <- a[is.finite(a)]
      b <- b[is.finite(b)]
      # A direct engine oracle verifies sample, orientation and Welch interval,
      # without using a MIRA analytical or report helper for expected values.
      test <- stats::t.test(b, a, paired = FALSE, var.equal = FALSE,
        alternative = "two.sided", conf.level = 0.90)
      rows[[length(rows) + 1L]] <- data.frame(time = visit,
        arm_a = pair[[1L]], arm_b = pair[[2L]], n_a = length(a), n_b = length(b),
        score_mean_difference_b_minus_a = mean(b) - mean(a),
        score_ci_lower = test$conf.int[[1L]], score_ci_upper = test$conf.int[[2L]],
        score_welch_t_p = test$p.value, stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, rows)
}

report_arm_identity_key <- function(table) {
  paste(table$time, table$arm_a, table$arm_b, sep = "::")
}

test_that("arm comparisons retain complete scientific keys in every published table block", {
  data <- report_arm_identity_data()
  fit <- mira_info_long(data, id = "patient", outcomes = "score", arm = "arm",
    reference_arm = "control", time_vars = c("score_t1", "score_t0"),
    time_labels = c(score_t0 = "Baseline", score_t1 = "Final"),
    covariates = character(), analyses = "arm_tests", alpha = 0.10, verbose = FALSE)
  expected <- report_arm_identity_oracle(data)
  expect_identical(anyDuplicated(report_arm_identity_key(expected)), 0L)
  expect_identical(expected$n_a, c(3L, 3L, 4L, 2L, 2L, 4L))
  expect_identical(expected$n_b, c(4L, 5L, 5L, 4L, 4L, 4L))
  root <- withr::local_tempdir(pattern = "mira-arm-identity-")
  report <- report_test_call(list(x = fit, sections = "arms"), root)
  manifest <- utils::read.csv(report$results_manifest, stringsAsFactors = FALSE)
  expect_identical(readRDS(report$payload)$result, fit)
  keys <- c("time", "arm_a", "arm_b")
  original_kable <- knitr::kable

  for (component in c("time_pairwise", "baseline_balance")) {
    source <- fit$arm_analysis[[component]]
    oracle <- if (component == "baseline_balance") {
      expected[expected$time == "score_t0", , drop = FALSE]
    } else expected
    # Establish exact cardinality and unique full keys before values.
    expect_identical(nrow(source), nrow(oracle), info = component)
    expect_identical(anyDuplicated(report_arm_identity_key(source)), 0L, info = component)
    expect_identical(report_arm_identity_key(source), report_arm_identity_key(oracle), info = component)
    for (column in setdiff(names(oracle), keys)) {
      # Deterministic engine-wrapper comparison; no optimiser or printed
      # rounding. 1e-12 allows only floating arithmetic differences.
      expect_equal(source[[column]], oracle[[column]], tolerance = 1e-12,
        info = paste(component, column))
    }
    row <- which(manifest$outcome == "score" &
      manifest$path == paste0("outcome$arm_analysis$", component))
    expect_length(row, 1L)
    csv <- utils::read.csv(file.path(root, manifest$file[[row]]),
      stringsAsFactors = FALSE, check.names = FALSE)
    expect_identical(nrow(csv), nrow(source))
    expect_identical(report_arm_identity_key(csv), report_arm_identity_key(oracle))
    for (column in setdiff(names(oracle), keys)) {
      expect_equal(csv[[column]], oracle[[column]], tolerance = 1e-12,
        info = paste("CSV", component, column))
    }
    for (target in c("html", "latex", "pipe")) {
      blocks <- list()
      rendered <- testthat::with_mocked_bindings(
        capture.output(.mira_report_table_long(source, caption = "Arm comparisons",
          cfg = report_test_cfg(max_table_columns = 7L))),
        is_latex_output = function() identical(target, "latex"),
        is_html_output = function() identical(target, "html"),
        kable = function(x, ...) {
          blocks[[length(blocks) + 1L]] <<- x
          original_kable(x, ...)
        }, .package = "knitr")
      expect_true(length(blocks) > 1L, info = paste(component, target))
      complete_keys <- vapply(blocks, function(block) all(keys %in% names(block)), logical(1L))
      expect_true(all(complete_keys), info = paste(component, target,
        "must repeat time, arm_a and arm_b in every block"))
      if (all(complete_keys)) {
        for (block in blocks) {
          expect_identical(report_arm_identity_key(block), report_arm_identity_key(oracle),
            info = paste(component, target))
          expect_identical(anyDuplicated(report_arm_identity_key(block)), 0L)
          for (column in intersect(names(block), setdiff(names(oracle), keys))) {
            # p-values have already received their publication formatting.
            if (column == "score_welch_t_p") {
              expect_identical(block[[column]], sprintf("%.3f", oracle[[column]]),
                info = paste(component, target, column))
            } else {
              expect_equal(block[[column]], oracle[[column]], tolerance = 1e-12,
                info = paste(component, target, column))
            }
          }
        }
      }
      if (target == "html" && all(complete_keys)) {
        html <- paste(rendered, collapse = "\n")
        for (i in seq_along(blocks)) {
          published <- report_test_html_table(html,
            sprintf("Arm comparisons (block %d of %d)", i, length(blocks)))
          expect_identical(nrow(published), nrow(oracle))
          published_keys <- paste(published$time, published[["arm a"]],
            published[["arm b"]], sep = "::")
          expect_identical(published_keys, report_arm_identity_key(oracle))
          expect_identical(anyDuplicated(published_keys), 0L)
        }
      }
      expect_true(all(lengths(lapply(blocks, names)) <= 7L))
      non_keys <- unlist(lapply(blocks, function(block) setdiff(names(block), keys)),
        use.names = FALSE)
      # Missing repeated keys would also duplicate incomplete time identifiers;
      # cardinality of non-keys is meaningful once the identity check passed.
      if (all(complete_keys)) {
        expect_identical(non_keys, setdiff(names(source), keys), info = paste(component, target))
      }
    }
    # At four columns a full visit/arm key and one statistic fit, even under a
    # tiny width budget. At three columns the documented fallback still keeps
    # every column and respects the requested limit.
    for (limit in c(3L, 4L)) {
      narrow <- .mira_report_partition_columns_long(source, limit, width_budget = 1)
      expect_true(all(lengths(narrow) <= limit))
      expect_identical(sort(unique(unlist(narrow))), seq_len(ncol(source)))
      if (limit == 4L) {
        key_columns <- match(keys, names(source))
        expect_true(all(vapply(narrow, function(block) all(key_columns %in% block), logical(1L))))
      }
    }
  }
  expect_identical(readRDS(report$payload)$result, fit)
})


