.mira_report_or <-
function(x, y) {
  if (is.null(x) || length(x) == 0L) y else x
}
.mira_report_is_empty <-
function(x) {
  if (is.null(x)) return(TRUE)
  if (is.table(x)) return(length(x) == 0L)
  if (is.data.frame(x) || is.matrix(x)) return(nrow(x) == 0L)
  if (is.atomic(x) || is.list(x)) return(length(x) == 0L)
  FALSE
}
.mira_report_heading <-
function(text, level = 1L) {
  text <- gsub("[\r\n]+", " ", as.character(text)[1L])
  level <- max(1L, min(6L, as.integer(level)[1L]))
  cat("\n", strrep("#", level), " ", text, "\n\n", sep = "")
}
.mira_report_human_name <-
function(x) {
  labels <- c(
    call = "Analysis call",
    version = "MIRA version",
    settings = "Analysis settings",
    overview = "Overview",
    config = "Resolved configuration",
    data_overview = "Dataset overview",
    detected_variables = "Detected variables",
    diagnostics = "Diagnostics and adaptations",
    descriptives = "Descriptive statistics",
    missing = "Missing data",
    by_time = "By time point",
    by_patient = "By participant",
    change = "Within-participant changes",
    arm_analysis = "Arm-specific analyses",
    baseline_balance = "Baseline balance",
    missingness = "Missing data by arm",
    time_omnibus = "Omnibus arm comparisons by time point",
    time_pairwise = "Pairwise arm comparisons by time point",
    change_descriptives = "Change by arm",
    change_omnibus = "Omnibus change comparisons",
    change_pairwise = "Pairwise change comparisons",
    correlations = "Correlations across time points",
    pearson = "Pearson correlation",
    spearman = "Spearman correlation",
    pairwise_n = "Pairwise sample sizes",
    variability = "Variability and ICC",
    trajectories = "Individual trajectories",
    model = "Primary mixed-effects model",
    fitted_model = "Fitted model object",
    summary = "Model summary",
    anova = "Model ANOVA table",
    global_time_test = "Global time test",
    global_arm_test = "Global arm test",
    arm_time_interaction_test = "Global arm-by-time interaction test",
    fixed_parameters = "Fixed-effect parameters",
    advanced_tests = "Advanced longitudinal tests",
    advanced_models = "Advanced longitudinal models",
    robustness = "Robust inference",
    effect_sizes = "Effect sizes",
    multiplicity = "Multiplicity families",
    sensitivity = "Sensitivity analyses",
    outliers = "Diagnostic outliers",
    plots = "Figures",
    plot_error = "Figure diagnostics",
    long_data = "Longitudinal analysis dataset",
    outcomes = "Outcomes",
    warnings = "Warnings",
    adaptation = "Automatic adaptations",
    failed_outcomes = "Failed outcomes",
    emmeans = "Estimated marginal means and contrasts",
    rm_anova = "Repeated-measures ANOVA",
    friedman = "Friedman test",
    random_slope = "Random-slope mixed-effects model",
    nlme = "NLME models",
    gee = "GEE models",
    model_comparison = "Descriptive model comparison",
    club_sandwich = "Cluster-robust CR2 inference",
    compound_symmetry = "Compound-symmetry correlation",
    ar1 = "AR(1) correlation",
    independence = "Independent working correlation",
    exchangeable = "Exchangeable working correlation",
    data = "Module analysis data",
    formula = "Model formula",
    error = "Module error",
    note = "Note",
    boxplot = "Distribution by time point",
    mean_ci = "Mean profile with confidence intervals",
    change_from_baseline = "Change from baseline",
    change_ci = "Mean change with confidence intervals",
    correlation_heatmap = "Correlation matrix",
    arm_mean_ci = "Mean profile by arm with confidence intervals",
    arm_boxplot = "Distribution by arm and time point",
    arm_change = "Change by arm",
    arm_change_ci = "Mean change by arm with confidence intervals",
    arm_difference_ci = "Arm differences with confidence intervals",
    arm_missingness = "Missing data by arm",
    response = "Individual clinical response"
  )
  key <- as.character(x)[1L]
  if (key %in% names(labels)) return(unname(labels[[key]]))
  text <- gsub("[._]+", " ", key)
  text <- trimws(text)
  if (!nzchar(text)) return("Unnamed item")
  paste0(toupper(substr(text, 1L, 1L)), substr(text, 2L, nchar(text)))
}
.mira_report_format_p <-
function(x, digits = 3L) {
  x <- suppressWarnings(as.numeric(x))
  vapply(x, function(value) {
    if (is.na(value)) return(NA_character_)
    if (!is.finite(value)) return(as.character(value))
    threshold <- 10^(-digits)
    if (value < threshold) {
      paste0("<", formatC(threshold, format = "f", digits = digits))
    } else {
      formatC(value, format = "f", digits = digits)
    }
  }, character(1L))
}
.mira_report_format_number <-
function(x, digits = 3L) {
  if (length(x) == 0L || is.na(x[1L])) return("not estimable")
  value <- suppressWarnings(as.numeric(x[1L]))
  if (!is.finite(value)) return(as.character(value))
  formatC(value, format = "fg", digits = digits, flag = "#")
}
.mira_report_inline <-
function(x) {
  if (is.null(x) || length(x) == 0L) return("not available")
  if (inherits(x, c("Date", "POSIXct", "POSIXlt"))) {
    return(paste(as.character(x), collapse = ", "))
  }
  if (is.atomic(x)) {
    values <- as.character(x)
    values[is.na(values)] <- "NA"
    return(paste(values, collapse = ", "))
  }
  paste(capture.output(str(x, give.attr = FALSE, vec.len = 8L)), collapse = " ")
}
.mira_report_prepare_table <-
function(x, digits = 3L, format_p = TRUE) {
  if (is.table(x)) x <- as.data.frame(x, stringsAsFactors = FALSE)

  if (is.matrix(x)) {
    row_labels <- rownames(x)
    # Some analytical matrices carry an S3 class such as "anova". Calling
    # as.data.frame() would dispatch that class's method and may change the
    # row count; convert by matrix structure instead.
    matrix_value <- x
    class(matrix_value) <- NULL
    x <- as.data.frame.matrix(matrix_value, stringsAsFactors = FALSE,
                              check.names = FALSE)
    if (!is.null(row_labels)) {
      x <- data.frame(.row = row_labels, x, check.names = FALSE,
                      stringsAsFactors = FALSE)
    }
  }

  if (!is.data.frame(x)) {
    if (is.list(x) && !is.object(x)) {
      nms <- names(x)
      if (is.null(nms)) nms <- paste0("item_", seq_along(x))
      x <- data.frame(
        parameter = nms,
        value = vapply(x, .mira_report_inline, character(1L)),
        stringsAsFactors = FALSE
      )
    } else if (is.atomic(x)) {
      nms <- names(x)
      if (is.null(nms)) nms <- seq_along(x)
      x <- data.frame(
        item = as.character(nms),
        value = as.character(x),
        stringsAsFactors = FALSE
      )
    } else {
      stop("The object cannot be converted to a table.", call. = FALSE)
    }
  }

  x <- as.data.frame(x, stringsAsFactors = FALSE, check.names = FALSE)
  for (j in seq_along(x)) {
    column <- x[[j]]

    # Some model summaries contain matrix/array columns (for example a
    # two-column confidence interval stored as one data-frame column). knitr's
    # numeric formatter cannot safely assign those values back to a data.frame.
    # Collapse every row of such columns explicitly before calling kable().
    if (!is.null(dim(column))) {
      row_n <- nrow(x)
      if (row_n > 0L && length(column) %% row_n == 0L) {
        flattened <- matrix(column, nrow = row_n)
        x[[j]] <- vapply(seq_len(row_n), function(i) {
          .mira_report_inline(flattened[i, , drop = TRUE])
        }, character(1L))
      } else {
        x[[j]] <- rep(.mira_report_inline(column), row_n)
      }
    } else if (inherits(column, c("POSIXct", "POSIXlt", "Date"))) {
      x[[j]] <- as.character(column)
    } else if (is.factor(column)) {
      x[[j]] <- as.character(column)
    } else if (is.list(column)) {
      row_n <- nrow(x)
      if (length(column) == row_n) {
        x[[j]] <- vapply(column, .mira_report_inline, character(1L))
      } else {
        x[[j]] <- rep(.mira_report_inline(column), row_n)
      }
    }
  }

  # Final invariant required by knitr::kable(): every data-frame column must be
  # a one-dimensional atomic vector with exactly nrow(x) elements.
  row_n <- nrow(x)
  for (j in seq_along(x)) {
    column <- x[[j]]
    valid_column <- is.atomic(column) && is.null(dim(column)) &&
      length(column) == row_n
    if (!valid_column) {
      if (row_n > 0L && length(column) == row_n) {
        x[[j]] <- vapply(seq_len(row_n), function(i) {
          .mira_report_inline(column[[i]])
        }, character(1L))
      } else {
        x[[j]] <- rep(.mira_report_inline(column), row_n)
      }
    }
  }

  if (format_p) {
    p_columns <- grepl(
      "(^p$|(^|[._])p([._]|$)|p.value|p_value|pvalue|^pr\\()",
      names(x), ignore.case = TRUE, perl = TRUE
    )
    for (j in which(p_columns)) {
      if (is.numeric(x[[j]]) || is.integer(x[[j]])) {
        x[[j]] <- .mira_report_format_p(x[[j]], digits = digits)
      }
    }
  }
  x
}
.mira_report_numeric_column <-
function(x, name = "") {
  if (is.numeric(x) || is.integer(x)) return(TRUE)
  if (grepl(
    "(^n$|count|percent|pct|mean|median|sd|se|ci|estimate|statistic|df|(^|[._])p([._]|$)|p.value|p_value|pvalue)",
    name, ignore.case = TRUE, perl = TRUE
  )) return(TRUE)
  values <- trimws(as.character(x))
  values <- values[!is.na(values) & nzchar(values) & values != "--"]
  if (length(values) == 0L) return(FALSE)
  values <- sub("^[<>]=?", "", values)
  mean(!is.na(suppressWarnings(as.numeric(values)))) >= 0.9
}
.mira_report_column_widths <-
function(x) {
  if (!is.data.frame(x) || ncol(x) == 0L) return(numeric())
  row_n <- nrow(x)
  sample_index <- if (row_n <= 200L) {
    seq_len(row_n)
  } else {
    unique(as.integer(round(seq(1, row_n, length.out = 200L))))
  }

  vapply(seq_along(x), function(j) {
    column_name <- names(x)[[j]]
    header <- nchar(gsub("[._]+", " ", column_name), type = "width")
    values <- if (length(sample_index) == 0L) character() else {
      as.character(x[[j]][sample_index])
    }
    values[is.na(values)] <- "--"
    values <- gsub("[[:space:]]+", " ", values)
    value_widths <- nchar(values, type = "width", allowNA = FALSE)
    typical <- if (length(value_widths) == 0L) 0 else {
      as.numeric(stats::quantile(
        value_widths, probs = 0.85, names = FALSE, type = 7, na.rm = TRUE
      ))
    }
    numeric_column <- .mira_report_numeric_column(x[[j]], column_name)
    cap <- if (numeric_column) 14 else 34
    max(4, min(cap, max(header, typical, na.rm = TRUE)))
  }, numeric(1L))
}
.mira_report_partition_columns <-
function(x, max_columns, width_budget = Inf) {
  column_n <- ncol(x)
  if (column_n == 0L) return(list(integer()))
  widths <- .mira_report_column_widths(x)
  if (column_n <= max_columns && sum(widths) <= width_budget) {
    return(list(seq_len(column_n)))
  }

  key_n <- min(2L, max_columns - 1L, column_n - 1L)
  if (key_n > 1L && is.finite(width_budget) &&
      sum(widths[seq_len(key_n)]) > width_budget * 0.4) {
    key_n <- 1L
  }
  key_n <- max(1L, key_n)
  key_columns <- seq_len(key_n)
  other_columns <- setdiff(seq_len(column_n), key_columns)
  if (length(other_columns) == 0L) return(list(key_columns))

  available_width <- if (is.finite(width_budget)) {
    max(10, width_budget - sum(widths[key_columns]))
  } else {
    Inf
  }
  available_columns <- max(1L, max_columns - length(key_columns))
  chunks <- list()
  current <- integer()
  for (column in other_columns) {
    exceeds_columns <- length(current) >= available_columns
    exceeds_width <- length(current) > 0L &&
      sum(widths[c(current, column)]) > available_width
    if (exceeds_columns || exceeds_width) {
      chunks[[length(chunks) + 1L]] <- current
      current <- integer()
    }
    current <- c(current, column)
  }
  if (length(current) > 0L) chunks[[length(chunks) + 1L]] <- current
  lapply(chunks, function(indices) unique(c(key_columns, indices)))
}
.mira_report_latex_alignment <-
function(x) {
  widths <- .mira_report_column_widths(x)
  if (length(widths) == 0L) return(NULL)
  # Reserve the exact inter-column padding (2 * tabcolsep per column) while
  # using nearly all of the available text width. The constants correspond to
  # A4 with 25 mm margins and the 3 pt tabcolsep set around every table.
  usable_fraction <- max(0.82, 0.965 - 0.0135 * length(widths))
  fractions <- usable_fraction * widths / sum(widths)
  vapply(seq_along(x), function(j) {
    direction <- if (.mira_report_numeric_column(x[[j]], names(x)[[j]])) {
      "\\raggedleft"
    } else {
      "\\raggedright"
    }
    sprintf(
      ">%s\\arraybackslash}p{%.4f\\linewidth}",
      paste0("{", direction), fractions[[j]]
    )
  }, character(1L))
}
.mira_report_latex_breaks <-
function(x) {
  for (token in c("\\_", "\\$")) {
    pieces <- strsplit(x, token, fixed = TRUE)[[1L]]
    if (length(pieces) > 1L) {
      x <- paste(pieces, collapse = paste0(token, "\\allowbreak{}"))
    }
  }
  x
}
.mira_report_repeat_longtable_header <-
function(x) {
  lines <- strsplit(x, "\n", fixed = TRUE)[[1L]]
  trimmed <- trimws(lines)
  if (any(startsWith(trimmed, "\\endfirsthead"))) return(x)
  top_index <- which(startsWith(trimmed, "\\toprule"))
  if (length(top_index) == 0L) return(x)
  top_index <- top_index[[1L]]
  mid_index <- which(
    seq_along(lines) > top_index & startsWith(trimmed, "\\midrule")
  )
  if (length(mid_index) == 0L || mid_index[[1L]] <= top_index + 1L) return(x)
  mid_index <- mid_index[[1L]]
  header_lines <- lines[seq.int(top_index + 1L, mid_index - 1L)]
  repeated_header <- c(
    "\\endfirsthead",
    "\\toprule",
    header_lines,
    "\\midrule",
    "\\endhead"
  )
  lines <- append(lines, repeated_header, after = mid_index)
  paste(lines, collapse = "\n")
}
.mira_report_table <-
function(x, caption = NULL, cfg) {
  if (.mira_report_is_empty(x)) {
    cat("*No rows are available for display.*\n\n")
    return(invisible(NULL))
  }

  table_data <- tryCatch(
    .mira_report_prepare_table(x, digits = cfg$digits),
    error = function(e) NULL
  )
  if (is.null(table_data)) {
    .mira_report_text_block(x)
    return(invisible(NULL))
  }

  original_n <- nrow(table_data)
  if (is.finite(cfg$max_table_rows) && original_n > cfg$max_table_rows) {
    table_data <- table_data[seq_len(cfg$max_table_rows), , drop = FALSE]
  }

  is_latex <- knitr::is_latex_output()
  is_html <- knitr::is_html_output()
  output_format <- if (is_latex) "latex" else if (is_html) "html" else "pipe"
  max_columns <- max(2L, as.integer(cfg$max_table_columns))
  if (is_latex) max_columns <- min(max_columns, 8L)
  if (!is_latex && !is_html) max_columns <- min(max_columns, 7L)
  width_budget <- if (is_latex) 76 else if (is_html) Inf else 68
  blocks <- .mira_report_partition_columns(
    table_data, max_columns = max_columns, width_budget = width_budget
  )
  old_kable_options <- options(knitr.kable.NA = "--")
  on.exit(options(old_kable_options), add = TRUE)

  for (i in seq_along(blocks)) {
    block_caption <- caption
    if (length(blocks) > 1L && !is.null(caption)) {
      block_caption <- sprintf("%s (block %d of %d)", caption, i, length(blocks))
    }
    block <- table_data[, blocks[[i]], drop = FALSE]
    common <- list(
      x = block,
      format = output_format,
      digits = cfg$digits,
      row.names = FALSE,
      col.names = gsub("[._]+", " ", names(block)),
      caption = block_caption,
      label = NA,
      escape = TRUE
    )
    if (is_latex) {
      common$booktabs <- TRUE
      common$longtable <- TRUE
      # LaTeX alignment specifications are accepted directly by knitr. This
      # avoids rewriting the \begin{longtable} line after rendering and keeps
      # every brace under knitr's own table construction.
      common$align <- .mira_report_latex_alignment(block)
    } else if (is_html) {
      common$table.attr <- "class=\"mira-table\""
    }
    render_error <- NULL
    rendered <- tryCatch(
      do.call(knitr::kable, common),
      error = function(e) {
        render_error <<- conditionMessage(e)
        NULL
      }
    )
    if (is.null(rendered)) {
      .mira_report_callout(
        paste0(
          "The table could not be formatted by knitr (`",
          render_error,
          "`). Its contents are shown as text so that the report can continue."
        ),
        type = "warning", title = "Table formatting fallback"
      )
      .mira_report_text_block(block)
    } else {
      rendered_text <- as.character(rendered)
      if (is_latex) {
        rendered_text <- .mira_report_repeat_longtable_header(rendered_text)

        # In a booktabs longtable, a first data cell beginning with "(" or
        # "[" can be parsed as an optional trimming argument of \toprule or
        # \midrule.  The resulting TeX error mentions \cmrsideswitch.  A
        # strategically placed \relax terminates that argument scan while
        # leaving the printed table unchanged.
        rendered_text <- gsub(
          "\\\\toprule(?=\\s*(?:\\(|\\[))",
          "\\\\toprule\\\\relax",
          rendered_text,
          perl = TRUE
        )
        rendered_text <- gsub(
          "\\\\midrule(?=\\s*(?:\\(|\\[))",
          "\\\\midrule\\\\relax",
          rendered_text,
          perl = TRUE
        )

        # Technical identifiers and object paths often contain underscores or
        # dollar signs. They are kept verbatim, but explicit break points let
        # LaTeX wrap them inside fixed-width columns instead of crossing the
        # page margin.
        rendered_text <- .mira_report_latex_breaks(rendered_text)
        rendered_text <- paste0(
          "\\begingroup\n",
          "\\footnotesize\n",
          "\\setlength{\\tabcolsep}{3pt}\n",
          "\\setlength{\\LTleft}{0pt plus 1fill}\n",
          "\\setlength{\\LTright}{0pt plus 1fill}\n",
          "\\renewcommand{\\arraystretch}{1.14}\n",
          rendered_text,
          "\n\\endgroup"
        )
      }
      if (is_html) cat("<div class=\"mira-table-wrap\">\n")
      cat(rendered_text, "\n", sep = "")
      if (is_html) cat("</div>\n")
      cat("\n")
    }
  }

  if (is.finite(cfg$max_table_rows) && original_n > cfg$max_table_rows) {
    cat(sprintf(
      "*Showing %d of %d rows because of the explicit `max_table_rows` limit.*\n\n",
      nrow(table_data), original_n
    ))
  }
  invisible(table_data)
}
.mira_report_prepare_namespace <-
function(x) {
  package_map <- list(
    ggplot2 = c("ggplot", "ggplot2::ggplot"),
    lme4 = c("merMod", "lmerMod", "glmerMod", "summary.merMod"),
    nlme = c("lme", "gls"),
    geepack = c("geeglm", "gee", "summary.geeglm"),
    emmeans = c("emmGrid", "emm_list", "summary_emm"),
    afex = c("afex_aov"),
    clubSandwich = c("coef_test_clubSandwich", "Wald_test_clubSandwich"),
    performance = c("performance_model", "r2_nakagawa")
  )
  cls <- class(x)
  for (package in names(package_map)) {
    if (any(cls %in% package_map[[package]])) {
      suppressWarnings(requireNamespace(package, quietly = TRUE))
    }
  }
  invisible(NULL)
}
.mira_report_text_block <-
function(x) {
  .mira_report_prepare_namespace(x)
  text <- tryCatch(
    {
      if (inherits(x, c("lm", "glm", "merMod", "lmerMod", "glmerMod"))) {
        capture.output(print(summary(x)))
      } else {
        capture.output(print(x))
      }
    },
    error = function(e) paste("The object could not be printed:", conditionMessage(e))
  )
  text <- gsub("```", "'''", enc2utf8(text), fixed = TRUE)
  cat("```text\n", paste(text, collapse = "\n"), "\n```\n\n", sep = "")
}
.mira_report_callout <-
function(text, type = "note", title = "Note") {
  allowed <- c("note", "tip", "warning", "important", "caution")
  if (!type %in% allowed) type <- "note"
  cat(sprintf("::: {.callout-%s appearance=\"simple\"}\n", type))
  cat("**", title, ".** ", text, "\n", sep = "")
  cat(":::\n\n")
}
.mira_report_is_plot <-
function(x) {
  inherits(x, c("ggplot", "ggplot2::ggplot", "recordedplot", "grob", "gTree"))
}
.mira_report_plot <-
function(x, name = "Figure") {
  .mira_report_prepare_namespace(x)
  printed <- tryCatch({
    if (inherits(x, "recordedplot")) {
      grDevices::replayPlot(x)
    } else if (inherits(x, c("grob", "gTree"))) {
      grid::grid.newpage()
      grid::grid.draw(x)
    } else {
      print(x)
    }
    TRUE
  }, error = function(e) {
    .mira_report_callout(conditionMessage(e), type = "warning",
                         title = paste("Figure could not be rendered:", name))
    FALSE
  })
  if (printed) {
    caption_map <- c(
      boxplot = "Distribution of observed values by time point, with individual observations and confidence intervals for the mean.",
      trajectories = "Individual longitudinal trajectories.",
      trajectory = "Individual longitudinal trajectories.",
      spaghetti = "Individual longitudinal trajectories.",
      mean = "Longitudinal mean profile.",
      mean_ci = "Longitudinal mean profile with confidence intervals.",
      change = "Distribution of within-participant changes.",
      change_from_baseline = "Change from baseline during follow-up.",
      change_ci = "Mean longitudinal change with confidence intervals.",
      missingness = "Percentage of missing observations by time point.",
      correlation = "Correlation structure across time points.",
      arm = "Longitudinal pattern stratified by arm.",
      arm_mean_ci = "Longitudinal mean profile by arm with confidence intervals.",
      arm_boxplot = "Distribution of observed values by arm and time point.",
      arm_change = "Longitudinal change stratified by arm.",
      arm_change_ci = "Mean change by arm with confidence intervals.",
      arm_difference_ci = "Estimated arm differences with confidence intervals.",
      arm_missingness = "Percentage of missing observations by arm and time point.",
      response = "Frequency of individual response directions."
    )
    matching_keys <- names(caption_map)[vapply(names(caption_map), function(k) {
      grepl(tolower(k), tolower(name), fixed = TRUE)
    }, logical(1L))]
    if (length(matching_keys) > 1L) {
      matching_keys <- matching_keys[order(nchar(matching_keys), decreasing = TRUE)]
    }
    caption <- if (length(matching_keys) > 0L) {
      caption_map[[matching_keys[[1L]]]]
    } else {
      paste0("Plot `", name, "` produced by MIRA.")
    }
    cat("**Figure.** ", caption, "\n\n", sep = "")
  }
  invisible(printed)
}
.mira_report_is_plain_list <-
function(x) {
  is.list(x) && !is.data.frame(x) && !.mira_report_is_plot(x) &&
    !.mira_report_is_model(x)
}
.mira_report_is_model <-
function(x) {
  inherits(x, c("lm", "glm", "merMod", "lme", "gls", "geeglm",
                "summary.geeglm", "emmGrid", "emm_list", "afex_aov", "htest"))
}
.mira_report_extract_outcomes <-
function(result) {
  if (inherits(result, "mira_detect")) return(list())
  if (is.list(result$outcomes) && length(result$outcomes)) {
    out <- result$outcomes
  } else {
    out <- list(result)
  }
  nms <- names(out)
  if (is.null(nms)) nms <- rep("", length(out))
  for (i in seq_along(out)) {
    if (is.na(nms[i]) || !nzchar(nms[i])) {
      candidate <- out[[i]]$outcome
      nms[i] <- if (is.character(candidate) && length(candidate) == 1L &&
                     !is.na(candidate) && nzchar(candidate)) candidate else paste0("outcome_", i)
    }
  }
  names(out) <- make.unique(nms)
  out
}
.mira_report_section <-
function(key) {
  known <- c(settings = "methods", time_vars = "methods",
             time_labels = "methods", overview = "overview", descriptives = "descriptives",
             missing = "missingness", change = "change", arm_analysis = "arms",
             correlations = "correlations", variability = "variability",
             model = "models", advanced_tests = "models", advanced_models = "models",
             robustness = "robustness", effect_sizes = "robustness",
             multiplicity = "robustness", sensitivity = "robustness",
             outliers = "outliers", trajectories = "trajectories",
             plots = "figures", diagnostics = "diagnostics")
  if (key %in% names(known)) unname(known[[key]]) else key
}
.mira_report_selected <-
function(cfg, section) {
  "all" %in% cfg$sections || section %in% cfg$sections
}
.mira_report_node_has_result <-
function(node) {
  if (is.null(node)) return(FALSE)
  if (node$kind %in% c("table", "plot", "model")) return(TRUE)
  if (node$kind == "group") {
    return(any(vapply(node$children, .mira_report_node_has_result, logical(1L))))
  }
  FALSE
}
.mira_report_prune <-
function(x, key, path, cfg, depth = 0L,
                               parent_performed = FALSE) {
  if (depth > cfg$max_depth || is.null(x) || is.function(x) ||
      is.environment(x) || inherits(x, "mira_info_error")) return(NULL)
  # These fields describe the machinery or source data, not analytical output.
  if (key %in% c("call", "outcomes", "long_data",
                 "performed", "enabled", "package", "status", "object_path",
                 "plot_error", "disabled", "failed_outcomes")) return(NULL)
  if (identical(key, "data") &&
      grepl("(^|\\$)(model|advanced_tests|advanced_models|robustness)(\\$|$)",
            path)) return(NULL)
  if (grepl("(^|\\$)overview(\\$|$)", path) &&
      identical(key, "arm_analysis")) return(NULL)
  if (.mira_report_is_plot(x)) {
    return(list(kind = "plot", key = key, path = path, value = x))
  }
  if (is.data.frame(x) || is.matrix(x) || is.table(x)) {
    tab <- x
    if (key %in% c("pearson", "spearman") && all(is.na(tab))) return(NULL)
    # Keep unavailable rows alongside fitted models, so the comparison shows
    # which structures were attempted without presenting an all-empty table.
    if (identical(key, "model_comparison") && is.data.frame(tab) &&
        "object_path" %in% names(tab)) {
      fitted <- !is.na(tab$object_path) & nzchar(tab$object_path)
      if (!any(fitted)) return(NULL)
    }
    dimensions <- if (is.table(tab)) dim(as.data.frame(tab)) else dim(tab)
    if (length(dimensions) < 2L || dimensions[2L] == 0L) return(NULL)
    if (dimensions[1L] == 0L && !parent_performed &&
        !grepl("(^|\\$)outliers(\\$|$)", path)) return(NULL)
    return(list(kind = "table", key = key, path = path, value = tab,
                original_class = paste(class(x), collapse = "/")))
  }
  if (.mira_report_is_model(x)) {
    return(list(kind = "model", key = key, path = path, value = x))
  }
  if (is.list(x)) {
    if (!length(x)) return(NULL)
    if (grepl("(^|\\$)multiplicity(\\$|$)", path) &&
        "table" %in% names(x) && is.null(x$table)) return(NULL)
    if (identical(key, "arm_analysis") && !isTRUE(x$enabled)) return(NULL)
    if (identical(key, "model") && "fitted_model" %in% names(x) &&
        is.null(x$fitted_model) && is.null(x$error) &&
        !length(x$warnings) && !length(x$covariates_skipped)) return(NULL)
    performed <- isTRUE(x$performed)
    failed <- identical(x$performed, FALSE)
    active <- if (failed) FALSE else performed || parent_performed
    nms <- names(x)
    if (is.null(nms)) nms <- rep("", length(x))
    children <- list()
    for (i in seq_along(x)) {
      child_key <- if (is.na(nms[i]) || !nzchar(nms[i])) paste0("item_", i) else nms[i]
      if (child_key %in% c("model", "object") && "model" %in% nms &&
          "object" %in% nms && !is.null(x$model) &&
          identical(x$model, x$object) && identical(child_key, "object")) next
      child_value <- x[[i]]
      if (identical(child_key, "summary") &&
          inherits(child_value, "summary.merMod") &&
          (inherits(x$fitted_model, "merMod") || inherits(x$model, "merMod"))) next
      if (identical(key, "nlme") && identical(child_key, "comparison") &&
          is.data.frame(child_value) && "model" %in% names(child_value)) {
        fitted <- nms[vapply(x, function(item) {
          is.list(item) && isTRUE(item$performed)
        }, logical(1L))]
        child_value <- child_value[child_value$model %in% fitted, , drop = FALSE]
      }
      child <- .mira_report_prune(child_value, child_key, paste0(path, "$", child_key),
                                  cfg, depth + 1L, active)
      aliases <- c("summary", "tidy", "fixed_effects", "coefficient_tests")
      if (!is.null(child) && child$key %in% aliases &&
          any(vapply(children, function(previous) {
            previous$key %in% aliases && previous$kind == child$kind &&
              previous$kind == "table" && identical(previous$value, child$value)
          }, logical(1L)))) next
      if (!is.null(child)) children[[child_key]] <- child
    }
    if (!length(children)) return(NULL)
    if (key %in% c("outliers", "plots", "multiplicity") &&
        !any(vapply(children, .mira_report_node_has_result, logical(1L)))) return(NULL)
    return(list(kind = "group", key = key, path = path, children = children))
  }
  if (inherits(x, "formula") || is.call(x) || is.expression(x) || is.name(x)) {
    return(list(kind = "formula", key = key, path = path, value = x))
  }
  if (is.atomic(x) && length(x) > 0L && !all(is.na(x))) {
    return(list(kind = "value", key = key, path = path, value = x))
  }
  NULL
}
.mira_report_order <-
function(keys) {
  scientific <- c("settings", "time_vars", "time_labels", "overview",
                  "descriptives", "missing", "change", "arm_analysis",
                  "correlations", "variability", "trajectories", "model",
                  "advanced_tests", "advanced_models", "robustness", "effect_sizes",
                  "multiplicity", "sensitivity", "outliers", "plots", "diagnostics")
  c(scientific[scientific %in% keys], keys[!keys %in% scientific])
}
.mira_report_outcome_tree <-
function(outcome, cfg) {
  if (inherits(outcome, "mira_info_error") || !is.list(outcome)) return(list())
  skip <- c("call", "version", "config", "data_overview",
            "detected_variables", "outcome", "outcome_display",
            "outcomes", "long_data", "plot_error")
  keys <- .mira_report_order(setdiff(names(outcome), skip))
  nodes <- list()
  for (key in keys) {
    if (!.mira_report_selected(cfg, .mira_report_section(key))) next
    component <- outcome[[key]]
    if (identical(key, "outliers") && is.list(component)) {
      desc <- outcome$descriptives
      if (is.list(component$by_time) && is.data.frame(desc) &&
          all(c("time", "n") %in% names(desc))) {
        for (name in names(component$by_time)) {
          tab <- component$by_time[[name]]
          count <- desc$n[match(name, desc$time)]
          if (is.data.frame(tab) && nrow(tab) == 0L &&
              length(count) == 1L && !is.na(count) && count < 4L) {
            component$by_time[[name]] <- NULL
          }
        }
      }
      changes <- outcome$change
      if (is.list(component$change) && is.data.frame(changes) &&
          all(c("from", "to", "n") %in% names(changes))) {
        comparison <- paste(changes$from, changes$to, sep = "_to_")
        for (name in names(component$change)) {
          tab <- component$change[[name]]
          count <- changes$n[match(name, comparison)]
          if (is.data.frame(tab) && nrow(tab) == 0L &&
              length(count) == 1L && !is.na(count) && count < 4L) {
            component$change[[name]] <- NULL
          }
        }
      }
    }
    node <- .mira_report_prune(component, key, paste0("outcome$", key), cfg)
    if (!is.null(node)) nodes[[key]] <- node
  }
  if (.mira_report_selected(cfg, "diagnostics")) {
    diagnostic <- nodes$diagnostics
    if (is.null(diagnostic)) diagnostic <- list(kind = "group", key = "diagnostics",
                                               path = "outcome$diagnostics", children = list())
    model <- outcome$model
    if (is.list(model) && is.null(model$fitted_model) && is.null(nodes$model)) {
      for (key in c("error", "warnings")) {
        node <- .mira_report_prune(model[[key]], paste0("model_", key),
                                   paste0("outcome$model$", key), cfg)
        if (!is.null(node)) diagnostic$children[[paste0("model_", key)]] <- node
      }
    }
    plot_error <- .mira_report_prune(outcome$plot_error, "figure_error",
                                     "outcome$plot_error", cfg)
    if (!is.null(plot_error)) diagnostic$children$figure_error <- plot_error
    if (length(diagnostic$children)) nodes$diagnostics <- diagnostic
  }
  nodes
}
.mira_report_context_tree <-
function(result, cfg) {
  keys <- c(config = "methods", data_overview = "data_quality",
            detected_variables = "data_quality", diagnostics = "diagnostics")
  nodes <- list()
  for (key in names(keys)) {
    if (!.mira_report_selected(cfg, keys[[key]])) next
    # Outcome diagnostics are already shown inside each outcome for single output.
    if (key == "diagnostics" && !inherits(result, "mira_info_multi") &&
        !inherits(result, "mira_detect")) next
    if (is.null(result[[key]])) next
    node <- .mira_report_prune(result[[key]], key, paste0("result$", key), cfg)
    if (!is.null(node)) nodes[[key]] <- node
  }
  if (inherits(result, "mira_info_multi") &&
      .mira_report_selected(cfg, "diagnostics")) {
    errors <- list()
    for (name in names(result$outcomes)) {
      item <- result$outcomes[[name]]
      if (!inherits(item, "mira_info_error")) next
      node <- .mira_report_prune(item$error, name,
                                 paste0("result$outcomes$", name, "$error"), cfg)
      if (!is.null(node)) errors[[name]] <- node
    }
    if (length(errors)) {
      diagnostic <- nodes$diagnostics
      if (is.null(diagnostic)) diagnostic <- list(kind = "group", key = "diagnostics",
                                                 path = "result$diagnostics", children = list())
      diagnostic$children$outcome_errors <- list(kind = "group", key = "outcome_errors",
                                                  path = "result$outcomes", children = errors)
      nodes$diagnostics <- diagnostic
    }
  }
  nodes
}
.mira_report_render_node <-
function(node, level, cfg) {
  label <- .mira_report_human_name(node$key)
  if (node$kind == "value") {
    cat("**", label, ":** ", .mira_report_inline(node$value), "\n\n", sep = "")
    return(invisible(NULL))
  }
  .mira_report_heading(label, level)
  if (node$kind == "group") {
    for (child in node$children) .mira_report_render_node(child, min(6L, level + 1L), cfg)
  } else if (node$kind == "table") {
    rows <- if (is.table(node$value)) nrow(as.data.frame(node$value)) else nrow(node$value)
    if (rows == 0L) {
      if (grepl("(^|\\$)outliers(\\$|$)", node$path)) {
        cat("**Outliers detected:** 0.\n\n")
      } else {
        cat("**Observed results:** 0 rows.\n\n")
      }
    } else {
      .mira_report_table(node$value, caption = label, cfg = cfg)
    }
    for (attribute in c("note", "interpretation")) {
      value <- attr(node$value, attribute, exact = TRUE)
      if (!is.null(value) && length(value)) {
        cat("*", .mira_report_human_name(attribute), ":* ",
            .mira_report_inline(value), "\n\n", sep = "")
      }
    }
  } else if (node$kind == "plot") {
    .mira_report_plot(node$value, node$key)
  } else if (node$kind %in% c("model", "formula")) {
    .mira_report_text_block(node$value)
  }
  invisible(NULL)
}
.mira_report_render_abstract <-
function(outcomes, cfg) {
  if (!.mira_report_selected(cfg, "abstract")) return(invisible(NULL))
  available <- names(outcomes)[vapply(outcomes, function(x) {
    length(.mira_report_outcome_tree(x, cfg)) > 0L
  }, logical(1L))]
  if (!length(available)) return(invisible(NULL))
  .mira_report_heading("Abstract", 1L)
  cat("This report presents the observed analytical results for ",
      paste(available, collapse = ", "), ".\n\n", sep = "")
  invisible(NULL)
}
.mira_report_render_conclusions <-
function(outcomes, cfg) {
  if (!.mira_report_selected(cfg, "conclusions")) return(invisible(NULL))
  rows <- list()
  for (name in names(outcomes)) {
    desc <- outcomes[[name]]$descriptives
    if (!is.data.frame(desc) || nrow(desc) < 2L || !"mean" %in% names(desc)) next
    first <- suppressWarnings(as.numeric(desc$mean[1L]))
    last <- suppressWarnings(as.numeric(desc$mean[nrow(desc)]))
    if (!is.finite(first) || !is.finite(last)) next
    rows[[length(rows) + 1L]] <- data.frame(outcome = name,
      mean_first_visit = first, mean_last_visit = last,
      descriptive_difference = last - first, stringsAsFactors = FALSE)
  }
  if (!length(rows)) return(invisible(NULL))
  .mira_report_heading("Descriptive summary", 1L)
  .mira_report_table(do.call(rbind, rows),
                     caption = "Differences between observed means", cfg = cfg)
  invisible(NULL)
}
.mira_report_knit <-
function(payload) {
  result <- payload$result
  cfg <- payload$report
  if (knitr::is_html_output()) .mira_report_html_style()
  outcomes <- .mira_report_extract_outcomes(result)
  context <- .mira_report_context_tree(result, cfg)
  .mira_report_render_abstract(outcomes, cfg)
  if (length(context)) {
    .mira_report_heading("Analysis context", 1L)
    for (node in context) .mira_report_render_node(node, 2L, cfg)
  }
  for (i in seq_along(outcomes)) {
    nodes <- .mira_report_outcome_tree(outcomes[[i]], cfg)
    if (!length(nodes)) next
    name <- names(outcomes)[i]
    display <- outcomes[[i]]$outcome_display
    if (is.null(display) || !length(display) || is.na(display[1L])) display <- name
    .mira_report_heading(paste0("Outcome: ", display[1L]), 1L)
    for (node in nodes) .mira_report_render_node(node, 2L, cfg)
  }
  .mira_report_render_conclusions(outcomes, cfg)
  extra <- payload$extra_objects
  if (.mira_report_selected(cfg, "appendix") && length(extra)) {
    extra_nodes <- lapply(names(extra), function(key) {
      .mira_report_prune(extra[[key]], key, paste0("extra_objects$", key), cfg)
    })
    extra_nodes <- Filter(Negate(is.null), extra_nodes)
    if (length(extra_nodes)) {
      .mira_report_heading("Supplementary objects", 1L)
      for (node in extra_nodes) .mira_report_render_node(node, 2L, cfg)
    }
  }
  invisible(NULL)
}
.mira_report_html_style <-
function() {
  cat("```{=html}\n<style>\nbody{line-height:1.58;color:#14213d} h1,h2,h3{color:#14213d} .mira-table-wrap{overflow-x:auto} .mira-table{width:100%;font-size:.88rem}\n</style>\n```\n\n")
}
