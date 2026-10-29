report_test_cfg <- function(sections = "all", ...) {
  defaults <- list(
    sections = sections,
    include_complete_output = TRUE,
    max_table_rows = Inf,
    max_table_columns = 12L,
    digits = 3L,
    max_depth = 50L
  )
  utils::modifyList(defaults, list(...))
}

report_test_mira <- function(outcome = "score") {
  pearson <- matrix(
    c(1, 0.5, 0.5, 1), nrow = 2L,
    dimnames = list(c("t0", "t1"), c("t0", "t1"))
  )
  result <- list(
    config = list(
      id = "patient",
      outcomes = outcome,
      analyses = list(
        core = c("descriptives", "missingness", "change", "variability",
                 "trajectories"),
        plots = FALSE,
        model = FALSE,
        outliers = FALSE,
        correlations = TRUE,
        arm_tests = FALSE
      )
    ),
    data_overview = data.frame(
      records = 6L, variables = 3L, stringsAsFactors = FALSE
    ),
    detected_variables = list(id = "patient", outcomes = outcome),
    outcome = outcome,
    outcome_display = paste("Outcome", outcome),
    settings = list(alpha = 0.05, p_adjust_method = "holm"),
    overview = data.frame(
      metric = c("participants", "timepoints"), value = c(6, 2),
      stringsAsFactors = FALSE
    ),
    descriptives = data.frame(
      time = c("t0", "t1"), n = c(6L, 6L), mean = c(10.25, 11.75),
      sd = c(1.5, 1.25), p.value = c(0.0004, 0.05),
      check.names = FALSE, stringsAsFactors = FALSE
    ),
    missing = list(
      by_time = data.frame(
        time = c("t0", "t1"), n_missing = c(0L, 1L),
        percent_missing = c(0, 100 / 6), stringsAsFactors = FALSE
      )
    ),
    change = data.frame(
      from = "t0", to = "t1", n = 5L, mean_change = 1.5,
      p_value = 0.012345, check.names = FALSE, stringsAsFactors = FALSE
    ),
    correlations = list(
      pearson = pearson,
      spearman = pearson - matrix(c(0, 0.1, 0.1, 0), 2L),
      pairwise_n = matrix(c(6L, 5L, 5L, 5L), 2L,
                          dimnames = dimnames(pearson))
    ),
    variability = list(
      summary = data.frame(component = c("within", "between"),
                           variance = c(1.25, 2.5), stringsAsFactors = FALSE)
    ),
    trajectories = data.frame(
      patient = paste0("P", 1:3), baseline = c(9, 10, 11),
      final = c(10, 12, 12), change = c(1, 2, 1),
      stringsAsFactors = FALSE
    )
  )
  structure(result, class = c("mira_info_long", "list"))
}

report_test_multi <- function(include_error = FALSE) {
  score <- report_test_mira("score")
  fatigue <- report_test_mira("fatigue")
  fatigue$descriptives$mean <- c(20.5, 18.25)
  outcomes <- list(score = score, fatigue = fatigue)
  if (include_error) {
    outcomes$failed <- structure(
      list(outcome = "failed", error = "controlled outcome failure"),
      class = c("mira_info_error_long", "list")
    )
  }
  structure(
    list(
      config = list(outcomes = names(outcomes)),
      data_overview = data.frame(records = 6L, variables = 5L),
      detected_variables = list(outcomes = names(outcomes)),
      outcomes = outcomes,
      diagnostics = list(
        warnings = character(),
        failed_outcomes = if (include_error) "failed" else character()
      )
    ),
    class = c("mira_info_multi_long", "mira_info_long", "list")
  )
}

report_test_detect <- function() {
  structure(
    list(
      config = list(
        id = "patient",
        outcomes = "score",
        time_vars = c("score_t0", "score_t1"),
        analyses = list(core = TRUE, model = FALSE)
      ),
      data_overview = data.frame(records = 4L, variables = 3L),
      detected_variables = list(id = "patient", outcomes = "score"),
      diagnostics = list(warnings = character())
    ),
    class = c("mira_detect_long", "list")
  )
}

report_test_run <- function(..., .capture = TRUE) {
  if (!.capture) return(mira_report_freq_long(...))
  printed <- utils::capture.output(value <- mira_report_freq_long(...))
  attr(value, "test_printed_output") <- printed
  value
}

report_test_call <- function(args = list(), output_dir = tempfile("mira-report-test-")) {
  call_args <- list(
    x = report_test_mira(),
    output_dir = output_dir,
    format = "html",
    render = FALSE,
    open = FALSE
  )
  call_args[names(args)] <- args
  printed <- utils::capture.output(value <- do.call(mira_report_freq_long, call_args))
  attr(value, "test_printed_output") <- printed
  value
}

report_test_atomic_table <- function(x) {
  is.data.frame(x) && all(vapply(x, function(column) {
    is.atomic(column) && is.null(dim(column)) && length(column) == nrow(x)
  }, logical(1L)))
}

report_test_relpath <- function(path) {
  gsub("\\\\", "/", path)
}

# Real integrations have an explicit profile. Enabling one turns missing
# prerequisites into an error, so an integration job cannot pass by skipping it.
report_test_require_quarto <- function(profile = "MIRA_RUN_QUARTO_TESTS") {
  if (!identical(tolower(Sys.getenv(profile)), "true")) {
    testthat::skip(paste0("Real rendering requires ", profile, "=true"))
  }
  for (package in c("knitr", "quarto")) {
    if (!requireNamespace(package, quietly = TRUE)) {
      stop(profile, " requires the installed R package '", package, "'.", call. = FALSE)
    }
  }
  cli <- tryCatch(quarto::quarto_path(), error = function(e) NULL)
  if (is.null(cli) || length(cli) != 1L || is.na(cli) || !nzchar(cli) ||
      !file.exists(cli)) {
    stop(profile, " requires an existing Quarto CLI; tests never install it.",
         call. = FALSE)
  }
  invisible(cli)
}

# These parsers read published table cells, independently of the report's table
# preparation helpers. They intentionally support the simple tables used here.
report_test_html_text <- function(html) {
  html <- gsub("(?is)<br\\b[^>]*>|</(?:p|div|li)>", " ", html, perl = TRUE)
  text <- gsub("(?s)<[^>]+>", "", html, perl = TRUE)
  entities <- gregexpr("&(?:amp|lt|gt|quot|apos|nbsp|#[0-9]+|#[xX][0-9A-Fa-f]+);",
                      text, perl = TRUE)
  matched <- regmatches(text, entities)
  replacements <- lapply(matched, function(values) {
    vapply(values, function(entity) {
      known <- c("&amp;" = "&", "&lt;" = "<", "&gt;" = ">", "&quot;" = '"',
                 "&apos;" = "'", "&nbsp;" = " ")
      if (entity %in% names(known)) return(unname(known[[entity]]))
      code <- substr(entity, 3L, nchar(entity) - 1L)
      value <- if (grepl("^[xX]", code)) {
        strtoi(substring(code, 2L), base = 16L)
      } else {
        suppressWarnings(as.integer(code))
      }
      if (is.na(value) || value <= 0L || value > 1114111L ||
          (value >= 55296L && value <= 57343L)) return(entity)
      intToUtf8(value)
    }, character(1L), USE.NAMES = FALSE)
  })
  # Replace matches from the original string once: &amp;lt; means literal &lt;.
  regmatches(text, entities) <- replacements
  text <- gsub("\u00a0", " ", text, fixed = TRUE)
  trimws(gsub("[[:space:]]+", " ", text))
}

report_test_html_table <- function(html, caption) {
  matches <- function(text, pattern) regmatches(text, gregexpr(pattern, text, perl = TRUE))[[1L]]
  tables <- matches(html, "(?s)<table\\b[^>]*>.*?</table>")
  if (length(matches(html, "<table\\b")) != length(tables) ||
      length(matches(html, "</table>")) != length(tables)) {
    stop("HTML has malformed or missing table boundaries.", call. = FALSE)
  }
  captions <- vapply(tables, function(table) {
    element <- regmatches(table, regexpr("(?s)<caption\\b[^>]*>.*?</caption>",
                                         table, perl = TRUE))
    if (!length(element)) return(NA_character_)
    report_test_html_text(element)
  }, character(1L))
  selected <- tables[!is.na(captions) & captions == caption]
  if (length(selected) != 1L) {
    stop(sprintf("Expected exactly one HTML table with caption '%s'; found %d.",
                 caption, length(selected)), call. = FALSE)
  }
  rows <- regmatches(selected, gregexpr("(?s)<tr\\b[^>]*>.*?</tr>", selected,
                                       perl = TRUE))[[1L]]
  fail <- function(message) stop(sprintf("HTML table '%s': %s.", caption, message), call. = FALSE)
  if (length(matches(selected, "<tr\\b")) != length(rows) ||
      length(matches(selected, "</tr>")) != length(rows)) fail("malformed or missing rows")
  if (!length(rows)) fail("no header row")
  cells <- lapply(seq_along(rows), function(i) {
    row <- rows[[i]]
    if (grepl("(?i)\\b(rowspan|colspan)\\s*=", row, perl = TRUE)) {
      fail("rowspan and colspan are unsupported")
    }
    # Match each tag to its own closing tag; a td closed by th is malformed.
    values <- matches(row, "(?s)<(th|td)\\b[^>]*>.*?</\\1>")
    if (length(matches(row, "<t[hd]\\b")) != length(values) ||
        length(matches(row, "</t[hd]>")) != length(values)) {
      fail(sprintf("malformed or mismatched cells in row %d", i))
    }
    if (i == 1L && (!length(values) || any(!grepl("^<th\\b", values)))) {
      fail("first row must contain only header cells")
    }
    if (i > 1L && any(!grepl("^<td\\b", values))) {
      fail(sprintf("data row %d must contain only data cells", i - 1L))
    }
    report_test_html_text(values)
  })
  headers <- cells[[1L]]
  if (any(!nzchar(headers)) || anyDuplicated(headers)) fail("header names must be nonempty and unique")
  widths <- lengths(cells[-1L])
  if (any(widths != length(headers))) {
    i <- which(widths != length(headers))[[1L]]
    fail(sprintf("data row %d has %d cells; expected %d", i, widths[[i]], length(headers)))
  }
  result <- as.data.frame(matrix(as.character(unlist(cells[-1L], use.names = FALSE)),
    ncol = length(headers), byrow = TRUE), stringsAsFactors = FALSE)
  names(result) <- headers
  result
}

# DOCX has no th tags. Identify the scientific table by a unique required-header
# signature and validate rectangular cells before selecting or indexing it.
report_test_docx_table <- function(xml, required_headers) {
  matches <- function(text, pattern) regmatches(text, gregexpr(pattern, text, perl = TRUE))[[1L]]
  fail <- function(message) stop(paste0("DOCX table: ", message, "."), call. = FALSE)
  tables <- matches(xml, "(?s)<w:tbl\\b[^>]*>.*?</w:tbl>")
  if (length(matches(xml, "<w:tbl\\b")) != length(tables) ||
      length(matches(xml, "</w:tbl>")) != length(tables)) fail("malformed or missing table boundaries")
  parsed <- lapply(seq_along(tables), function(index) {
    table <- tables[[index]]
    rows <- matches(table, "(?s)<w:tr\\b[^>]*>.*?</w:tr>")
    if (length(matches(table, "<w:tr\\b")) != length(rows) ||
        length(matches(table, "</w:tr>")) != length(rows)) fail("malformed or missing rows")
    if (!length(rows)) fail("no header row")
    cells <- lapply(seq_along(rows), function(i) {
      row <- rows[[i]]
      if (grepl("<w:(gridSpan|vMerge)\\b", row, perl = TRUE)) fail("merged cells are unsupported")
      pieces <- matches(row, "(?s)<w:tc\\b[^>]*>.*?</w:tc>")
      if (length(matches(row, "<w:tc\\b")) != length(pieces) ||
          length(matches(row, "</w:tc>")) != length(pieces)) fail("malformed or missing cells")
      vapply(pieces, function(cell) {
        cell <- gsub("<w:(br|tab)\\b[^>]*/>", "<w:t> </w:t>", cell, perl = TRUE)
        cell <- gsub("<w:t\\b[^>]*/>", "<w:t></w:t>", cell, perl = TRUE)
        text <- matches(cell, "(?s)<w:t\\b[^>]*>.*?</w:t>")
        if (length(matches(cell, "<w:t\\b")) != length(text) ||
            length(matches(cell, "</w:t>")) != length(text)) fail("malformed text nodes")
        report_test_html_text(paste0(text, collapse = ""))
      }, character(1L))
    })
    headers <- cells[[1L]]
    if (!length(headers) || any(!nzchar(headers)) || anyDuplicated(headers)) {
      fail("header names must be nonempty and unique")
    }
    widths <- lengths(cells[-1L])
    if (any(widths != length(headers))) fail("data row width differs from the header")
    result <- as.data.frame(matrix(as.character(unlist(cells[-1L], use.names = FALSE)),
      ncol = length(headers), byrow = TRUE), stringsAsFactors = FALSE)
    names(result) <- headers
    result
  })
  selected <- Filter(function(table) all(required_headers %in% names(table)), parsed)
  if (length(selected) != 1L) {
    fail(sprintf("expected exactly one table with headers [%s]; found %d",
      paste(required_headers, collapse = ", "), length(selected)))
  }
  selected[[1L]]
}
# Shared small real-engine fixture; the methodology and real-document profiles
# deliberately use the same individual observations.
report_methodology_matrix <- function() {
  rbind(c(2, 5, 7), c(4, 6, 12), c(7, 8, 9), c(3, 7, 13),
    c(9, 11, 10), c(6, 9, 16), c(8, 14, 18), c(12, 13, 15))
}

report_methodology_fit <- function(y = report_methodology_matrix(),
    outcome = "Quantitative score", labels = c("Zbase", "Amid", "Mfinal"),
    analyses = "model", inspect_only = FALSE) {
  data <- data.frame(patient = sprintf("P%02d", seq_len(nrow(y))),
    y, check.names = FALSE)
  names(data)[-1L] <- paste0("signal_t", seq_len(ncol(y)) - 1L)
  mira_info_long(data, id = "patient", outcomes = outcome,
    time_vars = names(data)[-1L], time_labels = labels,
    covariates = character(), analyses = analyses, alpha = 0.1,
    inspect_only = inspect_only, verbose = FALSE)
}

report_methodology_oracle <- function(y = report_methodology_matrix()) {
  n <- nrow(y); d <- ncol(y)
  means <- colMeans(y)
  centered <- sweep(y, 2L, means, "-")
  sigma <- crossprod(centered) / (n - 1L)
  basis <- rbind(c(1, 0, 0), c(-1, 1, 0), c(-1, 0, 1))
  estimate <- as.vector(basis %*% means)
  se <- sqrt(diag(basis %*% (sigma / n) %*% t(basis)))
  residual <- y - rowMeans(y) - rep(means, each = n) + mean(y)
  kr_df <- (n - 1L) * (d - 1L)
  kr_f <- (n * sum((means - mean(y))^2) / (d - 1L)) /
    (sum(residual^2) / kr_df)
  # Complete one-group H0p: F_t is the empirical mid-CDF of visit t.
  # p_j = mean_k mean_t F_t(Y_kj). The subject influence vector is
  # Z_kj = mean_t {F_t(Y_kj) - F_j(Y_kt)}; Cov(p) = cov(Z)/N.
  # Compute pair comparisons directly, independently of nparLD or MIRA.
  mid <- function(x, sample) mean((x > sample) + 0.5 * (x == sample))
  p <- vapply(seq_len(d), function(j) mean(vapply(y[, j], function(x) {
    mean(vapply(seq_len(d), function(t) mid(x, y[, t]), numeric(1L)))
  }, numeric(1L))), numeric(1L))
  z <- matrix(0, n, d)
  for (k in seq_len(n)) for (j in seq_len(d)) {
    z[k, j] <- mean(vapply(seq_len(d), function(t) {
      mid(y[k, j], y[, t]) - mid(y[k, t], y[, j])
    }, numeric(1L)))
  }
  rank_cov <- stats::cov(z) / n
  rank_se <- sqrt(diag(rank_cov))
  projection <- diag(d) - matrix(1 / d, d, d)
  projected <- projection %*% rank_cov
  trace <- sum(diag(projected))
  ats <- as.numeric(t(p) %*% projection %*% p / trace)
  ats_df <- trace^2 / sum(diag(projected %*% projected))
  contrast <- cbind(-1, diag(d - 1L))
  difference <- contrast %*% p
  wts <- as.numeric(t(difference) %*%
    solve(contrast %*% rank_cov %*% t(contrast), difference))
  rank_logit_se <- rank_se / (p * (1 - p))
  list(kr = c(statistic = kr_f, df_num = d - 1L, df_den = kr_df,
      p_raw = stats::pf(kr_f, d - 1L, kr_df, lower.tail = FALSE)),
    coefficients = cbind(estimate = estimate, SE = se, df = n - 1L,
      lower.CL = estimate - stats::qt(.95, n - 1L) * se,
      upper.CL = estimate + stats::qt(.95, n - 1L) * se),
    rank_effects = cbind(relative_effect = p, SE = rank_se,
      ci_lower = stats::plogis(stats::qlogis(p) - stats::qnorm(.95) * rank_logit_se),
      ci_upper = stats::plogis(stats::qlogis(p) + stats::qnorm(.95) * rank_logit_se)),
    rank_tests = rbind(c(statistic = ats, df_num = ats_df, df_den = 1e10,
        p_raw = stats::pf(ats, ats_df, 1e10, lower.tail = FALSE)),
      c(statistic = wts, df_num = d - 1L, df_den = NA_real_,
        p_raw = stats::pchisq(wts, d - 1L, lower.tail = FALSE))))
}

# Scope by the actual document heading, then reuse the established cell parsers.
# This prevents an identically printed number in another method/outcome from
# satisfying a scientific assertion. DOCX is inspected as XML table cells.
report_methodology_scope <- function(document, heading, format) {
  pattern <- if (format == "html") "(?s)<h([1-6])\\b[^>]*>.*?</h\\1>" else
    "(?s)<w:p\\b[^>]*>.*?</w:p>"
  positions <- gregexpr(pattern, document, perl = TRUE)[[1L]]
  pieces <- regmatches(document, list(positions))[[1L]]
  if (format == "docx") {
    keep <- grepl('<w:pStyle[^>]*w:val="Heading[1-6]"', pieces)
    positions <- positions[keep]; pieces <- pieces[keep]
    # Keep Word's separators in heading text; offsets still refer to document.
    pieces <- gsub("<w:(tab|br)\\b[^>]*/>", " ", pieces, perl = TRUE)
  }
  levels <- as.integer(if (format == "html") {
    sub('^<h([1-6]).*', '\\1', pieces)
  } else sub('(?s).*<w:pStyle[^>]*w:val="Heading([1-6])".*', '\\1', pieces, perl = TRUE))
  labels <- report_test_html_text(pieces)
  labels <- sub('^[0-9]+(?:\\.[0-9]+)*\\.?[[:space:]]+', '', labels, perl = TRUE)
  selected <- which(labels == heading)
  if (length(selected) != 1L) stop("Expected exactly one document heading: ", heading,
    call. = FALSE)
  after <- which(seq_along(levels) > selected & levels <= levels[[selected]])
  end <- if (length(after)) positions[[after[[1L]]]] - 1L else nchar(document)
  substr(document, positions[[selected]], end)
}

report_methodology_table <- function(document, format, headers, key = NULL, value = NULL) {
  if (format == "docx" && is.null(key)) return(report_test_docx_table(document, headers))
  pattern <- if (format == "html") "(?s)<table\\b[^>]*>.*?</table>" else
    "(?s)<w:tbl\\b[^>]*>.*?</w:tbl>"
  tables <- regmatches(document, gregexpr(pattern, document, perl = TRUE))[[1L]]
  parsed <- lapply(tables, function(table) {
    if (format == "docx") {
      return(report_test_docx_table(table, character()))
    }
    caption <- regmatches(table, regexpr('(?s)<caption\\b[^>]*>.*?</caption>', table, perl = TRUE))
    report_test_html_table(table, report_test_html_text(caption))
  })
  selected <- Filter(function(table) all(headers %in% names(table)) &&
    (is.null(key) || (key %in% names(table) && identical(table[[key]], value))), parsed)
  if (length(selected) != 1L) stop("Expected exactly one methodology table: ",
    paste(headers, collapse = ", "), call. = FALSE)
  selected[[1L]]
}

report_methodology_document <- function(document, format, expected = report_methodology_oracle()) {
  document <- report_methodology_scope(document, "Outcome: QUANTITATIVE SCORE", format)
  scope <- function(label) report_methodology_scope(document, label, format)
  numeric <- function(table, metric, value) report_test_numeric_cells(
    table[[paste("Quantitative score", gsub('[._]+', ' ', metric))]], value, metric)
  p_label <- function(x) unname(ifelse(x < .001, '<0.001', sprintf('%.3f', x)))
  same <- function(actual, expected, label) {
    if (!identical(actual, expected)) stop("Incorrect methodology keys/cells: ", label, call. = FALSE)
  }
  kr <- scope("Kenward-Roger joint F tests")
  for (metric in c('statistic', 'df_num', 'df_den', 'p_raw')) {
    table <- report_methodology_table(kr, format, c('question',
      paste('Quantitative score', gsub('[._]+', ' ', metric))))
    same(table$question, 'TIME', 'KR question')
    if (metric == 'p_raw') same(table[['Quantitative score p raw']], p_label(expected$kr[[metric]]), 'KR p') else
      numeric(table, metric, expected$kr[[metric]])
  }
  rank <- scope("Nonparametric relative marginal effects")
  for (metric in colnames(expected$rank_effects)) {
    table <- report_methodology_table(rank, format, c('time',
      paste('Quantitative score', gsub('[._]+', ' ', metric))))
    same(table$time, paste0('signal_t', 0:2), 'rank visits')
    numeric(table, metric, expected$rank_effects[, metric])
  }
  for (metric in c('statistic', 'df_num', 'df_den', 'p_raw')) {
    table <- report_methodology_table(rank, format, c('question', 'method',
      paste('Quantitative score', gsub('[._]+', ' ', metric))))
    same(table$question, c('TIME', 'TIME'), 'rank question')
    same(table$method, c('ATS', 'WTS'), 'rank method')
    if (metric == 'p_raw') same(table[['Quantitative score p raw']],
      p_label(expected$rank_tests[, metric]), 'rank p') else if (metric == 'df_den') {
      numeric(table[1L, , drop = FALSE], metric, 1e10)
      # Pandoc's smart reader renders the missing-value marker as an en dash.
      # Accept only the two spellings; a numeric WTS denominator is incorrect.
      same(table[['Quantitative score df den']][[2L]] %in% c('--', '\u2013'),
        TRUE, 'WTS absent denominator')
    } else numeric(table, metric, expected$rank_tests[, metric])
  }
  mmrm <- scope('Marginal repeated-measures model (MMRM)')
  for (metric in colnames(expected$coefficients)) {
    headers <- c('term', paste('Quantitative score', gsub('[._]+', ' ', metric)))
    # DOCX repeats the estimate beside term in each split coefficient block.
    if (metric == 'estimate') headers <- c(headers, 'Quantitative score SE')
    table <- report_methodology_table(mmrm, format, headers)
    same(table$term, c('(Intercept)', 'time_factorAmid', 'time_factorMfinal'), 'MMRM coefficient terms')
    numeric(table, metric, expected$coefficients[, metric])
  }
  for (metric in c('estimate', 'SE', 'df', 'lower.CL', 'upper.CL', 'p_raw')) {
    headers <- c('contrast', paste('Quantitative score', gsub('[._]+', ' ', metric)))
    if (metric == 'estimate') headers <- c(headers, 'Quantitative score SE')
    table <- report_methodology_table(mmrm, format, headers,
      key = 'contrast', value = 'final - baseline')
    value <- if (metric == 'p_raw') 2 * stats::pt(-abs(expected$coefficients[3L, 'estimate'] /
      expected$coefficients[3L, 'SE']), 7) else expected$coefficients[3L, metric]
    if (metric == 'p_raw') same(table[['Quantitative score p raw']], p_label(value), 'MMRM contrast p') else
      numeric(table, metric, value)
  }
  # Sample units remain counts, not display-rounded inferential quantities.
  table <- report_methodology_table(kr, format, c('question', 'n', 'n subjects'))
  same(table$question, 'TIME', 'KR count question')
  report_test_numeric_cells(table$n, 24, 'KR observations', rounding = 0)
  report_test_numeric_cells(table[['n subjects']], 8, 'KR subjects', rounding = 0)
  table <- report_methodology_table(rank, format, c('time', 'n subjects', 'n'))
  same(table$time, paste0('signal_t', 0:2), 'rank count visits')
  report_test_numeric_cells(table[['n subjects']], rep(8, 3), 'rank subjects', rounding = 0)
  report_test_numeric_cells(table$n, rep(8, 3), 'rank observations', rounding = 0)
  invisible(TRUE)
}

report_methodology_fixture_document <- function(format) {
  oracle <- report_methodology_oracle()
  header <- function(label, level) if (format == 'html') paste0('<h',level,'>',label,'</h',level,'>') else
    paste0('<w:p><w:pPr><w:pStyle w:val="Heading',level,'"/></w:pPr><w:r><w:t>',label,'</w:t></w:r></w:p>')
  table <- function(keys, metrics) {
    metric <- as.data.frame(metrics, check.names = FALSE)
    for (nm in names(metric)) metric[[nm]] <- if (nm == 'p_raw') ifelse(metric[[nm]] < .001,
      '<0.001',sprintf('%.3f',metric[[nm]])) else ifelse(is.na(metric[[nm]]),'--',sprintf('%.3f',metric[[nm]]))
    names(metric) <- paste('Quantitative score',gsub('[._]+',' ',names(metric)))
    tab <- cbind(keys,metric)
    names(tab) <- gsub('_', ' ', names(tab), fixed = TRUE)
    escape <- function(x) gsub('<','&lt;',gsub('&','&amp;',x,fixed=TRUE),fixed=TRUE)
    row <- function(x, head = FALSE) {
      tag <- if (head) 'th' else 'td'
      if (format == 'html') paste0('<tr>',paste0('<',tag,'>',escape(x),'</',tag,'>',collapse=''),'</tr>') else
        paste0('<w:tr>',paste0('<w:tc><w:p><w:r><w:t>',escape(x),'</w:t></w:r></w:p></w:tc>',collapse=''),'</w:tr>')
    }
    if(format=='html') paste0('<table><caption>Independent cells</caption>',row(names(tab),TRUE),
      paste(vapply(seq_len(nrow(tab)),function(i)row(as.character(tab[i,])),character(1)),collapse=''),'</table>') else
      paste0('<w:tbl>',row(names(tab),TRUE),paste(vapply(seq_len(nrow(tab)),function(i)row(as.character(tab[i,])),character(1)),collapse=''),'</w:tbl>')
  }
  paste0(header('Outcome: QUANTITATIVE SCORE',1), header('Kenward-Roger joint F tests',3),
    table(data.frame(question='TIME',n=24L,n_subjects=8L),as.data.frame(as.list(oracle$kr))),
    header('Nonparametric relative marginal effects',3),
    table(data.frame(arm='--',time=paste0('signal_t',0:2),n_subjects=8L,n=8L),oracle$rank_effects),
    table(data.frame(question=c('TIME','TIME'),term='time_factor',method=c('ATS','WTS')),oracle$rank_tests),
    header('Marginal repeated-measures model (MMRM)',3),
    table(data.frame(term=c('(Intercept)','time_factorAmid','time_factorMfinal')),oracle$coefficients),
    table(data.frame(contrast='final - baseline'),
      cbind(oracle$coefficients[3L,,drop=FALSE],p_raw=2*stats::pt(-abs(oracle$coefficients[3L,'estimate']/
        oracle$coefficients[3L,'SE']),7))))
}



# Narrow PDF verifier for this explicit four-column integration layout. It reads
# keyed rows under the actual method heading and named statistic header; a
# matching number elsewhere in the PDF cannot satisfy it. Unsupported wrapping
# fails inspection rather than silently falling back to an unscoped search.
report_methodology_pdf <- function(pages, expected = report_methodology_oracle()) {
  if (!length(pages) || anyNA(pages)) stop('PDF contains no methodology pages.',call.=FALSE)
  lines <- unlist(strsplit(paste(pages,collapse='\n'),'\n',fixed=TRUE),use.names=FALSE)
  labels <- sub('^[0-9]+(?:\\.[0-9]+)*[[:space:]]+','',trimws(lines),perl=TRUE)
  if (sum(labels == 'Outcome: Q') != 1L) stop('PDF requires one methodology outcome: Q.',call.=FALSE)
  scope <- function(start, end) {
    at <- which(labels==start)
    if(length(at)!=1L) stop('PDF requires one methodology heading: ',start,call.=FALSE)
    finish <- which(seq_along(labels)>at & labels %in% end)
    lines[seq.int(at,if(length(finish))finish[[1L]]-1L else length(lines))]
  }
  scopes <- list(
    kr=scope('Kenward-Roger joint F tests',c('Nonparametric relative marginal effects','Advanced longitudinal models')),
    rank=scope('Nonparametric relative marginal effects',c('Advanced longitudinal models','Marginal repeated-measures model (MMRM)')),
    mmrm=scope('Marginal repeated-measures model (MMRM)',c('Random-slope mixed-effects model','NLME models','GEE models','Descriptive model comparison')))
  cell <- function(module,metric,keys,expected_values,p=FALSE,extra_header=NULL) {
    section <- scopes[[module]]
    header <- paste('Q',gsub('[._]+',' ',metric))
    candidates <- which(vapply(section, function(line) {
      fields <- strsplit(trimws(line), '[[:space:]]{2,}', perl=TRUE)[[1L]]
      all(c(header, names(keys)) %in% fields)
    },logical(1L)))
    if(!is.null(extra_header)) candidates<-candidates[grepl(extra_header,section[candidates],fixed=TRUE)]
    # Estimate is repeated as a coefficient key in subsequent PDF blocks;
    # select its SE-containing block explicitly, just as for DOCX.
    if(length(candidates)!=1L) stop('PDF requires one keyed metric header: ',module,' ',metric,call.=FALSE)
    index<-candidates[[1L]]
    headers<-strsplit(trimws(section[[index]]),'[[:space:]]{2,}',perl=TRUE)[[1L]]
    column<-match(header,headers)
    if(is.na(column)||!all(names(keys)%in%headers)) stop('PDF methodology header wrapping is unsupported: ',metric,call.=FALSE)
    next_headers<-which(seq_along(section)>index &
      grepl('(^|[[:space:]])(question|arm|term)[[:space:]]{2,}',section,perl=TRUE))
    end<-if(length(next_headers))next_headers[[1L]]-1L else length(section)
    segment<-section[seq.int(index+1L,end)]
    first<-names(keys)[[1L]]
    parsed<-lapply(segment,function(line)strsplit(trimws(line),'[[:space:]]{2,}',perl=TRUE)[[1L]])
    parsed<-Filter(function(row)length(row)==length(headers) && row[[match(first,headers)]] %in% keys[[first]],parsed)
    if(length(parsed)!=length(expected_values)) stop('Incorrect PDF methodology row cardinality: ',module,' ',metric,call.=FALSE)
    table<-as.data.frame(do.call(rbind,parsed),stringsAsFactors=FALSE);names(table)<-headers
    for(key in names(keys)) if(!identical(table[[key]],keys[[key]])) stop('Incorrect PDF methodology keys: ',module,' ',key,call.=FALSE)
    if(p) {
      values<-unname(ifelse(expected_values<.001,'<0.001',sprintf('%.3f',expected_values)))
      if(!identical(table[[header]],values)) stop('Incorrect PDF methodology p: ',module,call.=FALSE)
    } else report_test_numeric_cells(table[[header]],expected_values,paste('PDF',module,metric))
  }
  kr_keys<-list(question='TIME',n='24')
  for(metric in names(expected$kr)) cell('kr',metric,kr_keys,unname(expected$kr[[metric]]),p=metric=='p_raw')
  rank_keys<-list(arm=rep('--',3),time=paste0('signal_t',0:2))
  for(metric in colnames(expected$rank_effects)) cell('rank',metric,rank_keys,expected$rank_effects[,metric])
  test_keys<-list(question=c('TIME','TIME'),term=rep('time_factor',2),method=c('ATS','WTS'))
  for(metric in c('statistic','df_num','p_raw')) cell('rank',metric,test_keys,expected$rank_tests[,metric],p=metric=='p_raw')
  coefficient_keys<-list(term=c('(Intercept)','time_factorAmid','time_factorMfinal'))
  for(metric in colnames(expected$coefficients)) cell('mmrm',metric,coefficient_keys,expected$coefficients[,metric],
    extra_header=if(metric=='estimate')'Q SE' else NULL)
  invisible(TRUE)
}

report_methodology_pdf_fixture <- function() {
  o <- report_methodology_oracle()
  block <- function(keys, metrics) {
    metrics<-as.data.frame(metrics,check.names=FALSE)
    for(nm in names(metrics))metrics[[nm]]<-if(nm=='p_raw')ifelse(metrics[[nm]]<.001,'<0.001',sprintf('%.3f',metrics[[nm]])) else sprintf('%.3f',metrics[[nm]])
    names(metrics)<-paste('Q',gsub('[._]+',' ',names(metrics)))
    table<-cbind(as.data.frame(keys),metrics)
    c(paste(names(table),collapse='  '),apply(table,1,paste,collapse='  '))
  }
  lines<-c('Outcome: Q','Kenward-Roger joint F tests')
  for(metric in names(o$kr))lines<-c(lines,block(list(question='TIME',n='24'),as.data.frame(as.list(o$kr[metric]))))
  lines<-c(lines,'Nonparametric relative marginal effects')
  for(metric in colnames(o$rank_effects))lines<-c(lines,block(list(arm=rep('--',3),time=paste0('signal_t',0:2)),o$rank_effects[,metric,drop=FALSE]))
  for(metric in c('statistic','df_num','p_raw'))lines<-c(lines,block(list(question=rep('TIME',2),term=rep('time_factor',2),method=c('ATS','WTS')),o$rank_tests[,metric,drop=FALSE]))
  lines<-c(lines,'Advanced longitudinal models','Marginal repeated-measures model (MMRM)')
  for(metrics in list(c('estimate','SE','df'),c('estimate','lower.CL','upper.CL')))lines<-c(lines,block(list(term=c('(Intercept)','time_factorAmid','time_factorMfinal')),o$coefficients[,metrics,drop=FALSE]))
  paste(lines,collapse='\n')
}

# Namespace environments are serialized by name (R Internals, Serialization
# Formats). Restoring an installed native fit can autoload MIRA without running
# its analyses. Start with MIRA absent, record that restoration, then prove no
# refit while the actual QMD reads its own payload and executes its runtime.
report_methodology_cold_main <- function() {
  stopifnot(!'MIRA' %in% loadedNamespaces())
  args <- commandArgs(TRUE)
  setwd(args[[1L]])
  restored <- readRDS('mira_report-payload.rds')
  namespace_restored <- 'MIRA' %in% loadedNamespaces()
  targets <- c('mira_info_long', '.mira_analyse_single_outcome_long',
    '.mira_run_advanced_inference_long', '.mira_run_kenward_roger_long',
    '.mira_run_rank_longitudinal_long', '.mira_run_mmrm_long')
  calls <- 0L
  blocker <- function(...) {
    calls <<- calls + 1L
    stop('Unexpected MIRA analytical refit during cold reporting.',call.=FALSE)
  }
  isolated <- new.env(parent=baseenv())
  for(name in targets) assign(name,blocker,isolated)
  knitr::opts_knit$set(rmarkdown.pandoc.to='html')
  knit <- function() knitr::knit('mira_report.qmd',output='cold.md',envir=isolated,quiet=TRUE)
  if(namespace_restored) {
    # Only MIRA bindings, restored by testthat even if knitting fails. This
    # development guard is confined to the fresh child process.
    testthat::with_mocked_bindings(knit(),mira_info_long=blocker,
      .mira_analyse_single_outcome_long=blocker,
      .mira_run_advanced_inference_long=blocker,
      .mira_run_kenward_roger_long=blocker,
      .mira_run_rank_longitudinal_long=blocker,
      .mira_run_mmrm_long=blocker,.package='MIRA')
  } else knit()
  stopifnot(identical(calls,0L))
  saveRDS(list(initially_loaded=FALSE,namespace_restored=namespace_restored,
    loaded_after_knit='MIRA' %in% loadedNamespaces(),refit_calls=calls), 'cold-status.rds')
  writeLines('COLD_RUNTIME_OK','cold-status.txt')
}
