# These tests exercise the public analysis -> saved payload -> CSV -> exported
# runtime -> real HTML kable chain. Expected scientific quantities come from
# subject-level arithmetic, not another MIRA adapter or its formatted output.
report_methodology_bundle <- function(fit, output_dir) {
  report <- report_test_run(fit, output_dir = output_dir, format = "html",
    render = FALSE, open = FALSE, include_complete_output = FALSE,
    max_table_columns = 30L, digits = 3L)
  runtime <- new.env(parent = baseenv())
  sys.source(report$runtime, envir = runtime)
  payload <- readRDS(report$payload)
  list(report = report, runtime = runtime, payload = payload,
    manifest = utils::read.csv(report$results_manifest, check.names = FALSE),
    nodes = runtime$.mira_report_outcome_tree_long(payload$result, payload$report))
}

report_methodology_csv <- function(bundle, path) {
  selected <- bundle$manifest[bundle$manifest$path == path, , drop = FALSE]
  expect_identical(nrow(selected), 1L, info = path)
  if (nrow(selected) != 1L) return(data.frame())
  file <- file.path(bundle$report$output_dir, selected$file)
  expect_true(file.exists(file), info = path)
  utils::read.csv(file, check.names = FALSE, stringsAsFactors = FALSE)
}

report_methodology_html <- function(bundle, node = NULL) {
  testthat::with_mocked_bindings(
    paste(capture.output(if (is.null(node)) {
      bundle$runtime$.mira_report_knit_long(bundle$payload)
    } else {
      bundle$runtime$.mira_report_render_node_long(node, 2L, bundle$payload$report)
    }), collapse = "\n"),
    is_html_output = function(...) TRUE,
    is_latex_output = function(...) FALSE, .package = "knitr")
}

report_methodology_numbers <- function(actual, expected) {
  # kable's three decimal places introduce at most half of one unit in the
  # last displayed place. A missing/shifted/duplicated cell must still fail.
  parsed <- suppressWarnings(as.numeric(actual))
  expect_length(parsed, length(expected))
  expect_true(all(is.finite(parsed)))
  if (length(parsed) == length(expected)) {
    expect_true(all(abs(parsed - expected) <= 0.0005000001))
  }
}

test_that("KR and MMRM subject-level oracles survive the real report artifacts", {
  skip_if_not_installed("pbkrtest")
  skip_if_not_installed("mmrm")
  skip_if_not_installed("nparLD", "2.3.3")
  y <- report_methodology_matrix()
  fit <- report_methodology_fit(y)
  expect_true(fit$advanced_tests$kenward_roger$performed)
  expect_true(fit$advanced_tests$rank_longitudinal$performed)
  expect_true(fit$advanced_models$mmrm$performed)
  before <- fit$advanced_models$mmrm$coefficients
  output <- withr::local_tempdir(pattern = "mira-method-report-")
  # A supplied analysis must be serialized, never silently fitted again.
  bundle <- testthat::with_mocked_bindings(
    report_methodology_bundle(fit, output),
    mira_info_long = function(...) stop("unexpected refit"), .package = "MIRA")
  saved <- bundle$payload$result
  expect_identical(fit$advanced_models$mmrm$coefficients, before)
  expect_identical(saved$outcome, "Quantitative score")
  expect_identical(saved$outcome_display, "QUANTITATIVE SCORE")
  expect_identical(unname(saved$time_labels), c("Zbase", "Amid", "Mfinal"))
  expect_identical(saved$time_vars, c("signal_t0", "signal_t1", "signal_t2"))
  expect_identical(saved$advanced_models$mmrm$coefficients, before)

  # Balanced random-intercept inference reduces to the repeated-measures
  # time/error sums of squares, with 2 and (8-1)*(3-1)=14 degrees of freedom.
  means <- colSums(y) / 8
  grand <- sum(y) / 24
  residual <- y - rowSums(y) / 3 - rep(means, each = 8) + grand
  f <- (8 * sum((means - grand)^2) / 2) / (sum(residual^2) / 14)
  p <- stats::pf(f, 2, 14, lower.tail = FALSE)
  kr <- saved$advanced_tests$kenward_roger$tidy
  expect_identical(kr$question, "TIME")
  expect_equal(kr[["Quantitative score_statistic"]], f, tolerance = 1e-8)
  expect_equal(kr[["Quantitative score_df_num"]], 2, tolerance = 1e-8)
  expect_equal(kr[["Quantitative score_df_den"]], 14, tolerance = 1e-8)
  # Tail probabilities amplify the sub-1e-8 relative optimizer tolerance in F.
  expect_lte(abs(kr[["Quantitative score_p_raw"]] - p), 1e-10)
  expect_identical(kr$n, 24L)
  expect_identical(kr$n_subjects, 8L)
  kr_csv <- report_methodology_csv(bundle, "outcome$advanced_tests$kenward_roger$tidy")
  expect_identical(kr_csv$question, "TIME")
  expect_equal(kr_csv[["Quantitative score_statistic"]], f, tolerance = 1e-8)
  expect_equal(kr_csv[["Quantitative score_df_den"]], 14, tolerance = 1e-8)
  expect_lte(abs(kr_csv[["Quantitative score_p_raw"]] - p), 1e-10)
  kr_node <- bundle$nodes$advanced_tests$children$kenward_roger$children$tidy
  kr_html <- report_test_html_table(report_methodology_html(bundle, kr_node), "Tidy")
  expect_identical(kr_html$question, "TIME")
  report_methodology_numbers(kr_html[["Quantitative score statistic"]], f)
  report_methodology_numbers(kr_html[["Quantitative score df den"]], 14)
  expect_identical(kr_html[["Quantitative score p raw"]], "<0.001")

  # Complete-profile REML has Sigma = subject residual crossproduct/(8-1).
  # The final-minus-baseline contrast has subject df=7, not repeated-row df.
  centered <- sweep(y, 2L, means, "-")
  sigma <- crossprod(centered) / 7
  basis <- rbind(c(1, 0, 0), c(-1, 1, 0), c(-1, 0, 1))
  estimate <- as.vector(basis %*% means)
  se <- sqrt(diag(basis %*% (sigma / 8) %*% t(basis)))
  terms <- c("(Intercept)", "time_factorAmid", "time_factorMfinal")
  mmrm <- saved$advanced_models$mmrm
  expect_identical(mmrm$coefficients$term, terms)
  expect_equal(mmrm$coefficients[["Quantitative score_estimate"]], estimate,
    tolerance = 1e-10)
  expect_equal(mmrm$coefficients[["Quantitative score_SE"]], se, tolerance = 2e-6)
  expect_equal(mmrm$coefficients[["Quantitative score_df"]], rep(7, 3),
    tolerance = 2e-5)
  expect_equal(unname(mmrm$covariance), unname(sigma), tolerance = 2e-6)
  coefficients_csv <- report_methodology_csv(bundle,
    "outcome$advanced_models$mmrm$coefficients")
  expect_identical(coefficients_csv$term, terms)
  expect_equal(coefficients_csv[["Quantitative score_estimate"]], estimate,
    tolerance = 1e-10)
  expect_equal(coefficients_csv[["Quantitative score_SE"]], se, tolerance = 2e-6)
  coefficients_node <- bundle$nodes$advanced_models$children$mmrm$children$coefficients
  coefficients_html <- report_test_html_table(
    report_methodology_html(bundle, coefficients_node), "Coefficients")
  expect_identical(coefficients_html$term, terms)
  report_methodology_numbers(coefficients_html[["Quantitative score estimate"]], estimate)
  report_methodology_numbers(coefficients_html[["Quantitative score SE"]], se)
  report_methodology_numbers(coefficients_html[["Quantitative score df"]], rep(7, 3))

  delta <- y[, 3L] - y[, 1L]
  delta_mean <- sum(delta) / 8
  delta_se <- sqrt(sum((delta - delta_mean)^2) / 7) / sqrt(8)
  delta_p <- 2 * stats::pt(-abs(delta_mean / delta_se), 7)
  contrast_csv <- report_methodology_csv(bundle,
    "outcome$advanced_models$mmrm$emmeans$baseline_final$tidy")
  expect_identical(contrast_csv$contrast, "final - baseline")
  expect_equal(contrast_csv[["Quantitative score_estimate"]], 6.125, tolerance = 1e-10)
  expect_equal(contrast_csv[["Quantitative score_SE"]], delta_se, tolerance = 2e-6)
  expect_equal(contrast_csv[["Quantitative score_p_raw"]], delta_p, tolerance = 2e-6)
  contrast_node <- bundle$nodes$advanced_models$children$mmrm$children$emmeans$children$
    baseline_final$children$tidy
  contrast_html <- report_test_html_table(report_methodology_html(bundle, contrast_node), "Tidy")
  expect_identical(contrast_html$contrast, "final - baseline")
  report_methodology_numbers(contrast_html[["Quantitative score estimate"]], 6.125)
  report_methodology_numbers(contrast_html[["Quantitative score SE"]], delta_se)
  expect_identical(contrast_html[["Quantitative score p raw"]], "0.003")

  # Prefixing applies only to declared metric owners. Native package objects
  # are saved intact and form terminal print nodes rather than export subtrees.
  expect_s3_class(saved$advanced_tests$kenward_roger$tests$TIME, "KRmodcomp")
  expect_s3_class(saved$advanced_tests$rank_longitudinal$object, "nparld_fit")
  expect_s3_class(mmrm$model, "mmrm")
  native_nodes <- list(
    bundle$nodes$advanced_tests$children$kenward_roger$children$tests$children$TIME,
    bundle$nodes$advanced_tests$children$rank_longitudinal$children$object,
    bundle$nodes$advanced_models$children$mmrm$children$model,
    bundle$nodes$advanced_models$children$mmrm$children$summary)
  for (node in native_nodes) {
    expect_identical(node$kind, "model")
    expect_null(node$children)
  }
  expect_false(any(grepl(
    "^outcome\\$(advanced_tests\\$(kenward_roger\\$tests|rank_longitudinal\\$object)|advanced_models\\$mmrm\\$(model|summary))(\\$|$)",
    bundle$manifest$path)))
  metric_paths <- vapply(attr(saved, "mira_metric_paths", exact = TRUE),
    paste, character(1L), collapse = "$")
  expect_true(all(c("advanced_tests$kenward_roger$tidy",
    "advanced_tests$rank_longitudinal$effects", "advanced_tests$rank_longitudinal$tests",
    "advanced_models$mmrm$coefficients", "advanced_models$mmrm$effect_tests",
    "advanced_models$mmrm$emmeans$baseline_final$tidy") %in% metric_paths))
  expect_identical(unname(attr(kr, "mira_metric_names")[["statistic"]]),
    "Quantitative score_statistic")
  expect_identical(unname(attr(mmrm$coefficients, "mira_metric_names")[["SE"]]),
    "Quantitative score_SE")
  references <- saved$sensitivity$object_path
  expect_true(all(c("$advanced_tests$kenward_roger$tidy",
    "$advanced_tests$rank_longitudinal$tests", "$advanced_models$mmrm$effect_tests") %in%
    references))
  expect_true("$advanced_models$mmrm$model" %in%
    saved$advanced_models$model_comparison$object_path)
  expect_identical(saved$multiplicity$mmrm$baseline_vs_followup$object_path,
    "$advanced_models$mmrm$emmeans$baseline_followup$tidy")
  full_html <- report_methodology_html(bundle)
  expect_match(full_html, "Kenward-Roger joint F tests", fixed = TRUE)
  expect_match(full_html, "Nonparametric relative marginal effects", fixed = TRUE)
  expect_match(full_html, "Marginal repeated-measures model", fixed = TRUE)
  expect_match(full_html, "Outcome: QUANTITATIVE SCORE", fixed = TRUE)
  expect_match(full_html, "Zbase", fixed = TRUE)
  expect_match(full_html, "Mfinal", fixed = TRUE)
})

test_that("rank effects and tests retain unequal observed-cell weights in the report", {
  skip_if_not_installed("nparLD", "2.3.3")
  # Three subjects, five observations. Equal visit weights yield 5/12 and
  # 7/12; pooling the five raw ranks would answer a different estimand.
  y <- cbind(c(1, 3, 3), c(2, 4, NA_real_))
  fit <- report_methodology_fit(y, outcome = "Rank panel",
    labels = c("Zbase", "Afinal"))
  expect_true(fit$advanced_tests$rank_longitudinal$performed)
  bundle <- report_methodology_bundle(fit,
    withr::local_tempdir(pattern = "mira-rank-report-"))
  rank <- bundle$payload$result$advanced_tests$rank_longitudinal
  effects <- rank$effects
  expected <- c(5, 7) / 12
  se <- rep(1 / 12, 2)
  expect_identical(effects$time, c("signal_t0", "signal_t1"))
  expect_identical(effects$time_label, c("Zbase", "Afinal"))
  expect_identical(effects$n, c(3L, 2L))
  expect_identical(effects$missing_n, c(0L, 1L))
  expect_equal(effects[["Rank panel_relative_effect"]], expected, tolerance = 1e-12)
  expect_equal(effects[["Rank panel_SE"]], se, tolerance = 1e-12)
  critical <- stats::qnorm(0.95)
  logit_se <- se / (expected * (1 - expected))
  expect_equal(effects[["Rank panel_ci_lower"]],
    stats::plogis(stats::qlogis(expected) - critical * logit_se), tolerance = 1e-12)
  expect_equal(effects[["Rank panel_ci_upper"]],
    stats::plogis(stats::qlogis(expected) + critical * logit_se), tolerance = 1e-12)
  effects_csv <- report_methodology_csv(bundle,
    "outcome$advanced_tests$rank_longitudinal$effects")
  expect_identical(effects_csv$time_label, c("Zbase", "Afinal"))
  expect_equal(effects_csv[["Rank panel_relative_effect"]], expected, tolerance = 1e-12)
  expect_equal(effects_csv[["Rank panel_SE"]], se, tolerance = 1e-12)
  node <- bundle$nodes$advanced_tests$children$rank_longitudinal
  effects_html <- report_test_html_table(
    report_methodology_html(bundle, node$children$effects), "Effects")
  expect_identical(effects_html[["time label"]], c("Zbase", "Afinal"))
  report_methodology_numbers(effects_html[["Rank panel relative effect"]], expected)
  report_methodology_numbers(effects_html[["Rank panel SE"]], se)

  expected_p <- c(stats::pf(1, 1, 1e10, lower.tail = FALSE),
    stats::pchisq(1, 1, lower.tail = FALSE))
  expect_identical(rank$tests$question, c("TIME", "TIME"))
  expect_identical(rank$tests$method, c("ATS", "WTS"))
  expect_equal(rank$tests[["Rank panel_statistic"]], c(1, 1), tolerance = 1e-12)
  expect_equal(rank$tests[["Rank panel_p_raw"]], expected_p, tolerance = 1e-12)
  tests_csv <- report_methodology_csv(bundle,
    "outcome$advanced_tests$rank_longitudinal$tests")
  expect_identical(tests_csv$method, c("ATS", "WTS"))
  expect_equal(tests_csv[["Rank panel_statistic"]], c(1, 1), tolerance = 1e-12)
  expect_equal(tests_csv[["Rank panel_p_raw"]], expected_p, tolerance = 1e-12)
  expect_true(is.na(tests_csv[["Rank panel_df_den"]][2L]))
  tests_html <- report_test_html_table(
    report_methodology_html(bundle, node$children$tests), "Tests")
  expect_identical(tests_html$method, c("ATS", "WTS"))
  report_methodology_numbers(tests_html[["Rank panel statistic"]], c(1, 1))
  expect_identical(tests_html[["Rank panel p raw"]], c("0.317", "0.317"))
  expect_identical(tests_html[["Rank panel df den"]][2L], "--")
  expect_identical(node$children$object$kind, "model")
  expect_null(node$children$object$children)
  expect_false(any(startsWith(bundle$manifest$path,
    "outcome$advanced_tests$rank_longitudinal$object$")))
})

test_that("disabled and inspection-only analyses never run or advertise the new methods", {
  stop_method <- function(...) stop("a disabled method was evaluated")
  fits <- testthat::with_mocked_bindings({
    list(disabled = report_methodology_fit(analyses = "none"),
      inspect = report_methodology_fit(inspect_only = TRUE))
  }, .mira_run_kenward_roger_long = stop_method,
    .mira_run_rank_longitudinal_long = stop_method,
    .mira_run_mmrm_long = stop_method, .package = "MIRA")
  expect_s3_class(fits$inspect, "mira_detect_long")
  modules <- list(fits$disabled$advanced_tests$kenward_roger,
    fits$disabled$advanced_tests$rank_longitudinal,
    fits$disabled$advanced_models$mmrm)
  for (module in modules) {
    expect_false(module$requested)
    expect_false(module$performed)
    expect_null(module$error)
  }
  for (fit in fits) {
    bundle <- report_methodology_bundle(fit,
      withr::local_tempdir(pattern = "mira-disabled-method-report-"))
    expect_false(any(grepl("kenward_roger|rank_longitudinal|\\$mmrm(\\$|$)",
      bundle$manifest$path)))
    html <- report_methodology_html(bundle)
    expect_false(grepl("Kenward-Roger joint F tests", html, fixed = TRUE))
    expect_false(grepl("Nonparametric relative marginal effects", html, fixed = TRUE))
    expect_false(grepl("Marginal repeated-measures model", html, fixed = TRUE))
    expect_false(grepl("a disabled method was evaluated", html, fixed = TRUE))
  }
})

test_that("requested reason-only failures remain visible without inventing method tables", {
  reasons <- c(kenward_roger = "KR unavailable: covariance is not identified.",
    rank_longitudinal = "Rank unavailable: no within-subject variation.",
    mmrm = "MMRM unavailable: a visit pair has no joint observations.")
  unavailable <- function(reason, package) {
    force(reason); force(package)
    function(...) list(performed = FALSE, package = package,
      reason_skipped = reason, error = NULL, warnings = character())
  }
  # Exercise the public isolation boundary while deliberately supplying no
  # error or warning: requested + reason is sufficient diagnostic evidence.
  fit <- testthat::with_mocked_bindings(report_methodology_fit(),
    .mira_run_kenward_roger_long = unavailable(reasons[["kenward_roger"]], "pbkrtest"),
    .mira_run_rank_longitudinal_long = unavailable(reasons[["rank_longitudinal"]], "nparLD"),
    .mira_run_mmrm_long = unavailable(reasons[["mmrm"]], "mmrm"), .package = "MIRA")
  expect_s4_class(fit$model$fitted_model, "merMod")
  expect_true(fit$advanced_tests$friedman$performed)
  bundle <- report_methodology_bundle(fit,
    withr::local_tempdir(pattern = "mira-unavailable-method-report-"))
  saved <- bundle$payload$result
  modules <- list(saved$advanced_tests$kenward_roger,
    saved$advanced_tests$rank_longitudinal, saved$advanced_models$mmrm)
  for (module in modules) {
    expect_true(module$requested)
    expect_false(module$performed)
    expect_null(module$error)
    expect_length(module$warnings, 0L)
  }
  keys <- c("advanced_tests$kenward_roger", "advanced_tests$rank_longitudinal",
    "advanced_models$mmrm")
  for (i in seq_along(keys)) {
    diagnostic <- bundle$nodes$diagnostics$children[[keys[[i]]]]
    expect_identical(diagnostic$children$reason$value, unname(reasons[[i]]))
  }
  expect_null(bundle$nodes$advanced_tests$children$kenward_roger)
  expect_null(bundle$nodes$advanced_tests$children$rank_longitudinal)
  expect_null(bundle$nodes$advanced_models$children$mmrm)
  expect_false(any(grepl("^outcome\\$(advanced_tests\\$(kenward_roger|rank_longitudinal)|advanced_models\\$mmrm)(\\$|$)",
    bundle$manifest$path)))
  html <- report_methodology_html(bundle)
  for (reason in reasons) expect_match(html, reason, fixed = TRUE)
  expect_true("outcome$advanced_tests$friedman$tidy" %in% bundle$manifest$path)
})

test_that("split rank inference tables repeat the method beside every scientific value", {
  fit <- report_methodology_fit()
  table <- fit$advanced_tests$rank_longitudinal$tests
  expect_true(fit$advanced_tests$rank_longitudinal$performed)
  keys <- c('question', 'term', 'method')
  old <- knitr::opts_knit$get()
  withr::defer(knitr::opts_knit$restore(old))
  for (target in c('html', 'docx', 'latex')) {
    # Real format detection and kable; no external namespace replacement.
    knitr::opts_knit$set(rmarkdown.pandoc.to=target)
    width <- switch(target,html=Inf,docx=68,latex=76)
    blocks <- .mira_report_partition_columns_long(table,7L,width_budget=width)
    expect_gt(length(blocks),1L)
    for (indices in blocks) {
      block <- table[,indices,drop=FALSE]
      expect_true(all(keys %in% names(block)),info=paste(target,paste(names(block),collapse=', ')))
      for (key in intersect(keys,names(block)))expect_identical(block[[key]],table[[key]])
    }
    text <- paste(capture.output(.mira_report_table_long(table,'Rank tests',
      report_test_cfg(max_table_columns=7L))),collapse='\n')
    if(target=='html') for(i in seq_along(blocks)) {
      observed <- report_test_html_table(text,sprintf('Rank tests (block %d of %d)',i,length(blocks)))
      expect_identical(observed$question,table$question)
      expect_identical(observed$term,table$term)
      expect_identical(observed$method,c('ATS','WTS'))
    }
    if(target=='docx') {
      lines<-strsplit(text,'\n',fixed=TRUE)[[1L]]
      headers<-lines[grepl('^\\|[[:space:]]*question[[:space:]]*\\|',lines)]
      expect_length(headers,length(blocks))
      expect_true(all(grepl('\\|[[:space:]]*term[[:space:]]*\\|[[:space:]]*method[[:space:]]*\\|',headers)))
    }
    if(target=='latex') {
      headers<-regmatches(text,gregexpr('question & term & method',text,fixed=TRUE))[[1L]]
      expect_length(headers,2L*length(blocks))
    }
    narrow <- .mira_report_partition_columns_long(table,3L,width_budget=1)
    expect_true(all(lengths(narrow)<=3L))
    expect_identical(sort(unique(unlist(narrow))),seq_len(ncol(table)))
  }
})
