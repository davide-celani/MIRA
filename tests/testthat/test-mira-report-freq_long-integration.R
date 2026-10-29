report_integration_data <- function() {
  data.frame(
    patient = seq_len(6L),
    score_t0 = c(1, 3, 2, 6, 5, 4),
    score_t1 = c(2, 4, 3, 5, 7, 6),
    score_t2 = c(4, 5, 6, 7, 8, 9)
  )
}

report_integration_multi_data <- function() {
  data.frame(
    patient = seq_len(6L),
    score_t0 = c(1, 3, 2, 6, 5, 4),
    score_t1 = c(2, 4, 3, 5, 7, 6),
    fatigue_t0 = c(9, 8, 7, 6, 5, 4),
    fatigue_t1 = c(8, 8, 6, 5, 4, 3)
  )
}

report_integration_fit <- function(data = report_integration_data(), ...) {
  args <- list(
    data = data,
    id = "patient",
    outcomes = "score",
    time_vars = c("score_t0", "score_t1", "score_t2"),
    covariates = character(0),
    analyses = "none",
    verbose = FALSE
  )
  dots <- list(...)
  args[names(dots)] <- dots
  do.call(mira_info_long, args)
}

report_integration_with_dir <- function(code) {
  path <- tempfile("mira-report-integration-")
  dir.create(path, recursive = TRUE)
  on.exit(unlink(path, recursive = TRUE, force = TRUE), add = TRUE)
  code(path)
}

report_integration_manifest <- function(report) {
  utils::read.csv(
    report$results_manifest,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

report_integration_read_export <- function(report, path, outcome) {
  manifest <- report_integration_manifest(report)
  index <- which(manifest$path == path & manifest$outcome == outcome)
  testthat::expect_length(index, 1L)
  utils::read.csv(
    file.path(report$output_dir, manifest$file[[index]]),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

test_that("a real single-outcome result is preserved in payload and scientific CSVs", {
  data <- report_integration_data()
  fit <- report_integration_fit(data)

  report_integration_with_dir(function(output_dir) {
    report <- report_test_run(
      fit,
      output_dir = output_dir,
      output_file = "single-integration",
      format = "html",
      render = FALSE,
      open = FALSE
    )

    expect_s3_class(report, "mira_report_freq_long")
    expect_true(all(file.exists(c(
      report$source,
      report$payload,
      report$runtime,
      report$results_manifest,
      report$result_files
    ))))

    payload <- readRDS(report$payload)
    expect_s3_class(payload$result, "mira_info_long")
    expect_equal(payload$result$descriptives, fit$descriptives)
    expect_equal(payload$result$missing$by_time, fit$missing$by_time)
    expect_equal(payload$result$change, fit$change)
    expect_equal(
      payload$result$advanced_tests$friedman$tidy,
      fit$advanced_tests$friedman$tidy
    )
    expect_identical(payload$result$config$analyses, fit$config$analyses)

    descriptives <- report_integration_read_export(
      report, "outcome$descriptives", "score"
    )
    expect_equal(descriptives$n, fit$descriptives$n)
    expect_equal(descriptives$score_mean, fit$descriptives$score_mean, tolerance = 1e-12)
    expect_equal(descriptives$score_sd, fit$descriptives$score_sd, tolerance = 1e-12)

    missing <- report_integration_read_export(
      report, "outcome$missing$by_time", "score"
    )
    expect_equal(missing$unavailable_n, fit$missing$by_time$unavailable_n)
    expect_equal(missing$unavailable_pct,
                 fit$missing$by_time$unavailable_pct, tolerance = 1e-12)

    change <- report_integration_read_export(report, "outcome$change", "score")
    expect_true(is.numeric(change$score_paired_t_p))
    expect_equal(change$n, fit$change$n)
    expect_equal(change$score_mean_change, fit$change$score_mean_change, tolerance = 1e-12)
    expect_equal(change$score_paired_t_p, fit$change$score_paired_t_p, tolerance = 1e-12)
    expect_false(any(grepl("<", as.character(change$score_paired_t_p), fixed = TRUE)))
  })
})

test_that("data and data-frame x inputs forward mira_info_long arguments and verbosity", {
  data <- report_integration_data()

  report_integration_with_dir(function(output_dir) {
    report <- report_test_run(
      data = data,
      id = "patient",
      outcomes = "score",
      time_vars = c("score_t0", "score_t1", "score_t2"),
      covariates = character(0),
      analyses = "none",
      p_adjust_method = "BH",
      improvement_direction = "higher",
      stable_threshold = 1,
      output_dir = output_dir,
      output_file = "data-argument",
      format = "html",
      render = FALSE,
      open = FALSE
    )

    expect_s3_class(report$mira_result, "mira_info_long")
    expect_identical(report$mira_result$settings$p_adjust_method, "BH")
    expect_identical(report$mira_result$settings$improvement_direction, "higher")
    expect_identical(report$mira_result$settings$stable_threshold, 1)
    expect_false(report$mira_result$config$analyses$model)
    expect_false(any(grepl(
      "^MIRA INFO$",
      attr(report, "test_printed_output")
    )))
  })

  report_integration_with_dir(function(output_dir) {
    report <- report_test_run(
      x = data,
      id = "patient",
      outcomes = "score",
      time_vars = c("score_t0", "score_t1", "score_t2"),
      time_labels = c("Baseline", "Middle", "Final"),
      covariates = character(0),
      analyses = "none",
      verbose = TRUE,
      output_dir = output_dir,
      output_file = "x-data-frame",
      format = "html",
      render = FALSE,
      open = FALSE
    )

    expect_s3_class(report$mira_result, "mira_info_long")
    expect_identical(
      unname(report$mira_result$time_labels),
      c("Baseline", "Middle", "Final")
    )
    expect_true(report$mira_result$config$specified_manually$time_labels)
    expect_true(any(grepl(
      "^MIRA INFO$",
      attr(report, "test_printed_output")
    )))
  })
})

test_that("a real multi-outcome result exports each outcome to its own directory", {
  data <- report_integration_multi_data()
  fit <- mira_info_long(
    data,
    id = "patient",
    outcomes = c("score", "fatigue"),
    time_vars = list(
      score = c("score_t0", "score_t1"),
      fatigue = c("fatigue_t0", "fatigue_t1")
    ),
    covariates = character(0),
    analyses = "none",
    verbose = FALSE
  )
  expect_identical(class(fit), c("mira_info_multi_long", "mira_info_long", "list"))

  report_integration_with_dir(function(output_dir) {
    report <- report_test_run(
      fit,
      output_dir = output_dir,
      output_file = "multi-integration",
      format = "html",
      render = FALSE,
      open = FALSE
    )
    manifest <- report_integration_manifest(report)
    outcome_rows <- !is.na(manifest$outcome) & nzchar(manifest$outcome)

    expect_setequal(unique(manifest$outcome[outcome_rows]), names(fit$outcomes))
    for (outcome in names(fit$outcomes)) {
      rows <- manifest$outcome == outcome
      prefix <- paste0("results/", .mira_report_file_slug_long(outcome), "/")
      expect_true(any(rows))
      expect_true(all(startsWith(manifest$file[rows], prefix)))
      expect_true(dir.exists(file.path(
        report$output_dir, "results", .mira_report_file_slug_long(outcome)
      )))
    }
    expect_true(all(file.exists(report$result_files)))
    expect_s3_class(readRDS(report$payload)$result, "mira_info_multi_long")
  })
})

test_that("mira_detect_long and inspect-only inputs create context without outcomes", {
  data <- report_integration_data()
  inputs <- list(
    direct = mira_detect_long(
      data,
      id = "patient",
      outcomes = "score",
      time_vars = c("score_t0", "score_t1", "score_t2"),
      covariates = character(0),
      verbose = FALSE
    ),
    inspect_only = mira_info_long(
      data,
      id = "patient",
      outcomes = "score",
      time_vars = c("score_t0", "score_t1", "score_t2"),
      covariates = character(0),
      analyses = "none",
      inspect_only = TRUE,
      verbose = FALSE
    )
  )

  for (name in names(inputs)) {
    report_integration_with_dir(function(output_dir) {
      report <- report_test_run(
        inputs[[name]],
        output_dir = output_dir,
        output_file = paste0("detect-", name),
        format = "html",
        render = FALSE,
        open = FALSE
      )
      payload_result <- readRDS(report$payload)$result
      manifest <- report_integration_manifest(report)

      expect_s3_class(payload_result, "mira_detect_long")
      expect_length(.mira_report_extract_outcomes_long(payload_result), 0L)
      expect_false(any(grepl("^outcome\\$", manifest$path)))
      expect_true(any(grepl("^result\\$", manifest$path)))
      expect_false(any(!is.na(manifest$outcome) & nzchar(manifest$outcome)))
    })
  }
})

test_that("a nested mira_info_error_long stays in the payload and only surfaces as diagnostics", {
  data <- data.frame(
    patient = seq_len(6L),
    good_t0 = c(1, 3, 2, 6, 5, 4),
    good_t1 = c(2, 4, 3, 5, 7, 6)
  )
  data$bad_t0 <- I(matrix(seq_len(12L), nrow = 6L))
  data$bad_t1 <- 2:7

  expect_warning(
    fit <- mira_info_long(
      data,
      id = "patient",
      outcomes = c("good", "bad"),
      covariates = character(0),
      analyses = "none",
      verbose = FALSE
    ),
    "Analyses were not completed for: bad"
  )
  expect_s3_class(fit$outcomes$bad, "mira_info_error_long")

  report_integration_with_dir(function(output_dir) {
    report <- report_test_run(
      fit,
      output_dir = output_dir,
      output_file = "partial-outcome-failure",
      format = "html",
      render = FALSE,
      open = FALSE
    )
    saved <- readRDS(report$payload)$result
    manifest <- report_integration_manifest(report)
    context <- .mira_report_context_tree_long(saved, report_test_cfg())

    expect_s3_class(saved$outcomes$bad, "mira_info_error_long")
    expect_identical(saved$outcomes$bad$error, fit$outcomes$bad$error)
    expect_false("bad" %in% manifest$outcome)
    expect_false(any(grepl("/bad/", report_test_relpath(manifest$file),
                           fixed = TRUE)))
    expect_identical(
      context$diagnostics$children$outcome_errors$children$bad$value,
      fit$outcomes$bad$error
    )
    expect_length(.mira_report_outcome_tree_long(saved$outcomes$bad,
                                            report_test_cfg()), 0L)
  })
})

test_that("Quarto renders a minimal HTML report when its package and CLI are available", {
  report_test_require_quarto()

  fit <- report_integration_fit()
  report_integration_with_dir(function(output_dir) {
    report <- report_test_run(
      fit,
      output_dir = output_dir,
      output_file = "quarto-html-integration",
      format = "html",
      render = TRUE,
      open = FALSE,
      quiet = TRUE,
      max_table_columns = 100L
    )

    expect_true(unname(report$rendered[["html"]]))
    expect_true(file.exists(report$expected_files[["html"]]))
    expect_gt(file.info(report$expected_files[["html"]])$size, 0)
    html <- paste(readLines(report$expected_files[["html"]], warn = FALSE),
                  collapse = "\n")
    expect_match(html, "<!DOCTYPE html>|<html", ignore.case = TRUE)
    expect_match(html, "MIRA: Longitudinal Analysis", fixed = TRUE)
    expect_identical(unname(report$files),
                     unname(report$expected_files[["html"]]))
    expect_true(is.na(report$render_errors[["html"]]))
    descriptives <- report_test_html_table(html, "Descriptive statistics")
    raw <- report_integration_data()[c("score_t0", "score_t1", "score_t2")]
    # Published cells have three decimal places: a half-unit in the last
    # displayed place is the rounding bound, distinct from numerical tolerance.
    expect_identical(descriptives$time, names(raw))
    expect_equal(as.numeric(descriptives$n), rep(6, 3), tolerance = 0)
    expect_true(all(abs(as.numeric(descriptives[["score mean"]]) - vapply(raw, mean, numeric(1L))) <=
                      0.0005000001))
    expect_true(all(abs(as.numeric(descriptives[["score sd"]]) - vapply(raw, stats::sd, numeric(1L))) <=
                      0.0005000001))
    change <- report_test_html_table(html, "Within-participant changes")
    report_test_contrast_keys(change, names(raw))
    for (i in seq_len(nrow(change))) {
      differences <- raw[[change$to[[i]]]] - raw[[change$from[[i]]]]
      expect_equal(as.numeric(change$n[[i]]), length(differences), tolerance = 0)
      statistic <- mean(differences) / (stats::sd(differences) / sqrt(length(differences)))
      expected_p <- 2 * stats::pt(-abs(statistic), df = length(differences) - 1L)
      expected_label <- if (expected_p < 0.001) "<0.001" else sprintf("%.3f", expected_p)
      expect_identical(change[["score paired t p"]][[i]], expected_label)
      expect_lte(abs(as.numeric(change[["score mean change"]][[i]]) - mean(differences)),
                  0.0005000001)
    }
    references <- regmatches(html, gregexpr("(?:src|href)=[\"'][^\"']*[\"']", html,
                                             perl = TRUE))[[1L]]
    expect_false(any(grepl("quarto-html-integration-render-", references, fixed = TRUE)))
    expect_length(list.files(output_dir, pattern = "quarto-html-integration-render-"), 0L)
  })
})

test_that("Quarto DOCX output contains scientific table cells", {
  report_test_require_quarto()
  data <- report_test_inference_data()
  fit <- report_integration_fit(data, outcomes = NULL, alpha = 0.10)
  report_integration_with_dir(function(output_dir) {
    report <- report_test_run(
      fit, output_dir = output_dir, format = "docx", sections = c("descriptives", "change"),
      render = TRUE, open = FALSE, date = "2025-01-15"
    )
    expect_identical(report$rendered, c(docx = TRUE))
    entries <- utils::unzip(report$files[["docx"]], list = TRUE)$Name
    expect_true("word/document.xml" %in% entries)
    expanded <- file.path(output_dir, "docx-inspection")
    utils::unzip(report$files[["docx"]], files = "word/document.xml", exdir = expanded)
    xml <- paste(readLines(file.path(expanded, "word/document.xml"), warn = FALSE),
                  collapse = "\n")
    raw <- data[c("score_t0", "score_t1", "score_t2")]
    expect_no_error(report_test_docx_inference(xml, raw, alpha = 0.10))
    expect_length(list.files(output_dir, pattern = "mira_report-render-"), 0L)
  })
})
# The fixtures are encoded independently of MIRA's formatter. Mutations change
# only plausible cells/keys and keep headings, rectangular tables and files.
test_that('new-method document verifiers reject intact but scientifically corrupted tables', {
  for (format in c('html','docx')) {
    document <- report_methodology_fixture_document(format)
    expect_no_error(report_methodology_document(document,format))
    # Pandoc's smart reader turns the missing-value marker into an en dash,
    # and its HTML writer may put cell text on a separate line.
    rendered <- gsub('>--<', '>&#x2013;<', document, fixed=TRUE)
    rendered <- gsub('>signal_t1<', '>\n signal_t1 \n<', rendered, fixed=TRUE)
    expect_no_error(report_methodology_document(rendered,format))
    broken_denominator <- gsub('>&#x2013;<', '>7.000<', rendered, fixed=TRUE)
    expect_error(report_methodology_document(broken_denominator,format),
      'WTS absent denominator',fixed=TRUE)
    broken_rendered_key <- sub('(>[[:space:]]*)signal_t1([[:space:]]*<)',
      '\\1signal_t0\\2',rendered,perl=TRUE)
    expect_false(identical(broken_rendered_key,rendered))
    expect_error(report_methodology_document(broken_rendered_key,format),
      'rank visits',fixed=TRUE)
    # KR numerator df=2 is unique here as a three-decimal scientific cell.
    broken_number <- sub(if(format=='html') '>2.000<' else '>2.000<', '>9.000<',document,fixed=TRUE)
    expect_false(identical(broken_number,document))
    expect_error(report_methodology_document(broken_number,format),'Incorrect scientific cells: df_num',fixed=TRUE)
    broken_key <- sub('>signal_t1<','>signal_t0<',document,fixed=TRUE)
    expect_false(identical(broken_key,document))
    expect_error(report_methodology_document(broken_key,format),'Incorrect methodology keys/cells: rank visits',fixed=TRUE)
    swapped <- sub('>ATS<','>WTS<',document,fixed=TRUE)
    expect_error(report_methodology_document(swapped,format),'rank method',fixed=TRUE)
  }
})

test_that('numbered DOCX headings retain exact methodology scopes', {
  document <- report_methodology_fixture_document('docx')
  for (variant in c('tab', 'break', 'period')) {
    numbered <- document
    for (level in c(1L, 3L)) {
      start <- paste0('<w:p><w:pPr><w:pStyle w:val="Heading', level,
        '"/></w:pPr>')
      number <- if (level == 1L) '1' else '1.0.1'
      if (variant == 'period' && level == 1L) number <- paste0(number, '.')
      separator <- switch(variant, tab = '<w:tab />', 'break' = '<w:br/>',
        period = '<w:t xml:space="preserve"> </w:t>')
      runs <- paste0('<w:r><w:rPr><w:rStyle w:val="SectionNumber"/></w:rPr>',
        '<w:t>', number, '</w:t></w:r><w:r>', separator, '</w:r>')
      numbered <- gsub(start, paste0(start, runs), numbered, fixed = TRUE)
    }
    expect_no_error(report_methodology_document(numbered, 'docx'))
    # Normalising labels must not change offsets used to slice the original
    # XML, admit another outcome, or select a body paragraph as a heading.
    sibling <- paste0('<w:p><w:pPr><w:pStyle w:val="Heading1"/></w:pPr>',
      '<w:r><w:t>Other outcome</w:t></w:r></w:p>',
      '<w:p><w:r><w:t>OTHER_OUTCOME_SENTINEL</w:t></w:r></w:p>')
    body <- '<w:p><w:r><w:t>Outcome: QUANTITATIVE SCORE</w:t></w:r></w:p>'
    expect_identical(report_methodology_scope(paste0(body, numbered, sibling),
      'Outcome: QUANTITATIVE SCORE', 'docx'), numbered)
    wrong <- sub('Outcome: QUANTITATIVE SCORE', 'Outcome: OTHER SCORE',
      numbered, fixed = TRUE)
    expect_error(report_methodology_scope(paste0(body, wrong),
      'Outcome: QUANTITATIVE SCORE', 'docx'),
      'Expected exactly one document heading', fixed = TRUE)
    expect_error(report_methodology_scope(paste0(numbered, numbered),
      'Outcome: QUANTITATIVE SCORE', 'docx'),
      'Expected exactly one document heading', fixed = TRUE)
    broken <- sub('>signal_t1<', '>signal_t0<', numbered, fixed = TRUE)
    expect_error(report_methodology_document(broken, 'docx'),
      'Incorrect methodology keys/cells: rank visits', fixed = TRUE)
  }
})

test_that('all new engines survive relocation and knitting in a fresh R process', {
  fit <- report_methodology_fit()
  expect_true(fit$advanced_tests$kenward_roger$performed)
  expect_true(fit$advanced_tests$rank_longitudinal$performed)
  expect_true(fit$advanced_models$mmrm$performed)
  root <- withr::local_tempdir(pattern='mira-cold-methodology-')
  report <- report_test_call(list(x=fit,sections='models',max_table_columns=30L),file.path(root,'original'))
  moved <- file.path(root,'moved bundle'); dir.create(moved)
  expect_true(all(file.copy(c(report$source,report$runtime,report$payload),moved)))
  driver <- file.path(root,'cold.R'); log <- file.path(root,'cold.log')
  writeLines(c('main <-',deparse(report_methodology_cold_main,width.cutoff=100L),
    'main()'),driver)
  executable <- file.path(R.home('bin'),if(.Platform$OS.type=='windows')'Rscript.exe' else 'Rscript')
  status <- system2(executable,c('--vanilla',shQuote(driver),shQuote(moved)),stdout=log,stderr=log)
  expect_identical(status,0L,info=paste(readLines(log,warn=FALSE),collapse='\n'))
  expect_true(file.exists(file.path(moved,'cold-status.txt')))
  if(status==0L) {
    state <- readRDS(file.path(moved,'cold-status.rds'))
    expect_false(state$initially_loaded)
    expect_identical(state$refit_calls,0L)
    # knitr executes the real QMD setup/runtime; scoped table checks below
    # parse its actual HTML cells independently of native models and MIRA.
    # The runtime emits Markdown headings; retain their levels when turning
    # only those heading lines into a small HTML inspection envelope.
    lines <- readLines(file.path(moved,'cold.md'),warn=FALSE)
    lines <- vapply(lines,function(line) {
      if(!grepl('^#{1,6} ',line)) return(line)
      level <- nchar(sub(' .*','',line)); label <- sub('^#{1,6} ','',line)
      paste0('<h',level,'>',gsub('\\','',label,fixed=TRUE),'</h',level,'>')
    },character(1))
    expect_no_error(report_methodology_document(paste(lines,collapse='\n'),'html'))
    expect_false(grepl(report$output_dir,paste(lines,collapse='\n'),fixed=TRUE))
  }
})

test_that('Quarto HTML and DOCX publish keyed KR MMRM and rank scientific cells', {
  report_test_require_quarto()
  fit <- report_methodology_fit()
  expect_true(fit$advanced_tests$kenward_roger$performed)
  expect_true(fit$advanced_tests$rank_longitudinal$performed)
  expect_true(fit$advanced_models$mmrm$performed)
  root <- withr::local_tempdir(pattern='mira-quarto-methodology-')
  report <- report_test_run(fit,output_dir=root,format=c('html','docx'),sections='models',
    max_table_columns=30L,render=TRUE,open=FALSE,date='2026-10-27')
  expect_identical(report$rendered,c(html=TRUE,docx=TRUE))
  if(isTRUE(report$rendered[['html']])) {
    document <- paste(readLines(report$files[['html']],warn=FALSE,encoding='UTF-8'),collapse='\n')
    expect_no_error(report_methodology_document(document,'html'))
    scope <- report_methodology_scope(document, 'Nonparametric relative marginal effects', 'html')
    broken <- sub(scope, sub('(>[[:space:]]*)signal_t1([[:space:]]*<)',
      '\\1signal_t0\\2', scope, perl=TRUE), document, fixed=TRUE)
    path <- file.path(root,'corrupted.html'); writeLines(enc2utf8(broken),path,useBytes=TRUE)
    expect_true(file.exists(path)); expect_false(identical(broken,document))
    expect_error(report_methodology_document(paste(readLines(path,warn=FALSE,encoding='UTF-8'),collapse='\n'),'html'),
      '^Incorrect methodology keys/cells: rank (count )?visits$')
  }
  if(isTRUE(report$rendered[['docx']])) {
    expanded <- file.path(root,'docx'); utils::unzip(report$files[['docx']],exdir=expanded)
    xml <- file.path(expanded,'word','document.xml')
    document <- paste(readLines(xml,warn=FALSE,encoding='UTF-8'),collapse='\n')
    expect_no_error(report_methodology_document(document,'docx'))
    scope <- report_methodology_scope(document, 'Nonparametric relative marginal effects', 'docx')
    broken <- sub(scope, sub('(>[[:space:]]*)signal_t1([[:space:]]*<)',
      '\\1signal_t0\\2', scope, perl=TRUE), document, fixed=TRUE)
    expect_false(identical(broken,document)); writeLines(enc2utf8(broken),xml,useBytes=TRUE)
    # Split DOCX tables can put the visit/count block before effect estimates.
    expect_error(report_methodology_document(paste(readLines(xml,warn=FALSE,encoding='UTF-8'),collapse='\n'),'docx'),
      '^Incorrect methodology keys/cells: rank (count )?visits$')
    # All package entries are retained; only a scientific XML key was changed.
    # Zip creation is unnecessary for cell verification: XML is read directly
    # from the actual DOCX archive, just as in the established parser tests.
  }
})


test_that('PDF methodology inspection rejects wrong cells and visit attribution within their sections', {
  text <- report_methodology_pdf_fixture()
  expect_no_error(report_methodology_pdf(text))
  wrong <- sub('TIME  24  2.000','TIME  24  9.000',text,fixed=TRUE)
  expect_false(identical(wrong,text))
  expect_error(report_methodology_pdf(wrong),'Incorrect scientific cells: PDF kr df_num',fixed=TRUE)
  wrong <- sub('signal_t1','signal_t0',text,fixed=TRUE)
  expect_error(report_methodology_pdf(wrong),'Incorrect PDF methodology keys: rank time',fixed=TRUE)
  wrong <- sub('ATS','WTS',text,fixed=TRUE)
  expect_error(report_methodology_pdf(wrong),'Incorrect PDF methodology keys: rank method',fixed=TRUE)
})
