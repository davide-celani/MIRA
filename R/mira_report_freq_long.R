# ============================================================
# MIRA_REPORT_FREQ
# Dynamic, publication-oriented Quarto reporting for mira_info_long()
#
# Typical use:
#   source("mira_info_long.R")
#   source("mira_report_freq_long.R")
#
#   fit <- mira_info_long(data = study_data, outcomes = c("score_a", "score_b"),
#                    analyses = "all", verbose = FALSE)
#   mira_report_freq_long(fit, format = c("html", "pdf"))
#
# Or run mira_info_long() and build the report in one call:
#   mira_report_freq_long(
#     data = study_data,
#     outcomes = "score_a",
#     arm = "treatment",
#     covariates = c("age", "sex"),
#     analyses = "all",
#     format = "html"
#   )
#
# Required for rendering: Quarto CLI, and the R packages quarto and knitr.
# PDF additionally requires a working LaTeX installation.
# ============================================================

.mira_report_or_long <- function(x, y) {
  if (is.null(x) || length(x) == 0L) y else x
}

.mira_report_is_empty_long <- function(x) {
  if (is.null(x)) return(TRUE)
  if (is.table(x)) return(length(x) == 0L)
  if (is.data.frame(x) || is.matrix(x)) return(nrow(x) == 0L)
  if (is.atomic(x) || is.list(x)) return(length(x) == 0L)
  FALSE
}

.mira_report_heading_long <- function(text, level = 1L) {
  text <- gsub("[\r\n]+", " ", as.character(text)[1L])
  level <- max(1L, min(6L, as.integer(level)[1L]))
  cat("\n", strrep("#", level), " ", text, "\n\n", sep = "")
}

.mira_report_human_name_long <- function(x) {
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

.mira_report_format_p_long <- function(x, digits = 3L) {
  x <- suppressWarnings(as.numeric(x))
  vapply(x, function(value) {
    if (is.na(value)) return(NA_character_)
    if (!is.finite(value)) return(as.character(value))
    threshold <- 10^(-digits)
    if (value >= 0 && value < threshold) {
      paste0("<", formatC(threshold, format = "f", digits = digits))
    } else {
      formatC(value, format = "f", digits = digits)
    }
  }, character(1L))
}

.mira_report_format_number_long <- function(x, digits = 3L) {
  if (length(x) == 0L || is.na(x[1L])) return("not estimable")
  value <- suppressWarnings(as.numeric(x[1L]))
  if (!is.finite(value)) return(as.character(value))
  formatC(value, format = "fg", digits = digits, flag = "#")
}

.mira_report_inline_long <- function(x) {
  if (is.null(x) || length(x) == 0L) return("not available")
  if (inherits(x, c("Date", "POSIXct", "POSIXlt"))) {
    return(paste(as.character(x), collapse = ", "))
  }
  if (is.atomic(x)) {
    values <- as.character(x)
    values[is.na(values)] <- "NA"
    return(paste(values, collapse = ", "))
  }
  paste(utils::capture.output(utils::str(x, give.attr = FALSE, vec.len = 8L)),
        collapse = " ")
}

.mira_report_prepare_table_long <- function(x, digits = 3L, format_p = TRUE) {
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
        value = vapply(x, .mira_report_inline_long, character(1L)),
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
          .mira_report_inline_long(flattened[i, , drop = TRUE])
        }, character(1L))
      } else {
        x[[j]] <- rep(.mira_report_inline_long(column), row_n)
      }
    } else if (inherits(column, c("POSIXct", "POSIXlt", "Date"))) {
      x[[j]] <- as.character(column)
    } else if (is.factor(column)) {
      x[[j]] <- as.character(column)
    } else if (is.list(column)) {
      row_n <- nrow(x)
      if (length(column) == row_n) {
        x[[j]] <- vapply(column, .mira_report_inline_long, character(1L))
      } else {
        x[[j]] <- rep(.mira_report_inline_long(column), row_n)
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
          .mira_report_inline_long(column[[i]])
        }, character(1L))
      } else {
        x[[j]] <- rep(.mira_report_inline_long(column), row_n)
      }
    }
  }

  if (format_p) {
    p_columns <- grepl(
      "((^|[._])p([._]|$)|(^|[._])p[._]?value([._]|$)|^pr\\()",
      names(x), ignore.case = TRUE, perl = TRUE
    )
    for (j in which(p_columns)) {
      if (is.numeric(x[[j]]) || is.integer(x[[j]])) {
        x[[j]] <- .mira_report_format_p_long(x[[j]], digits = digits)
      }
    }
  }
  x
}

.mira_report_numeric_column_long <- function(x, name = "") {
  if (is.numeric(x) || is.integer(x)) return(TRUE)
  if (grepl(
    paste0(
      "(^|[._])(n|count|percent|pct|mean|median|sd|se|ci|estimate|",
      "statistic|df|p|pvalue)([._]|$)|(^|[._])p[._]?value([._]|$)|^pr\\("
    ),
    name, ignore.case = TRUE, perl = TRUE
  )) return(TRUE)
  values <- trimws(as.character(x))
  values <- values[!is.na(values) & nzchar(values) & values != "--"]
  if (length(values) == 0L) return(FALSE)
  values <- sub("^[<>]=?", "", values)
  mean(!is.na(suppressWarnings(as.numeric(values)))) >= 0.9
}

.mira_report_column_widths_long <- function(x) {
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
    numeric_column <- .mira_report_numeric_column_long(x[[j]], column_name)
    cap <- if (numeric_column) 14 else 34
    max(4, min(cap, max(header, typical, na.rm = TRUE)))
  }, numeric(1L))
}

.mira_report_partition_columns_long <- function(x, max_columns, width_budget = Inf) {
  column_n <- ncol(x)
  if (column_n == 0L) return(list(integer()))
  widths <- .mira_report_column_widths_long(x)
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

.mira_report_latex_alignment_long <- function(x) {
  widths <- .mira_report_column_widths_long(x)
  if (length(widths) == 0L) return(NULL)
  # Reserve the exact inter-column padding (2 * tabcolsep per column) while
  # using nearly all of the available text width. The constants correspond to
  # A4 with 25 mm margins and the 3 pt tabcolsep set around every table.
  usable_fraction <- max(0.82, 0.965 - 0.0135 * length(widths))
  fractions <- usable_fraction * widths / sum(widths)
  vapply(seq_along(x), function(j) {
    direction <- if (.mira_report_numeric_column_long(x[[j]], names(x)[[j]])) {
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

.mira_report_latex_breaks_long <- function(x) {
  for (token in c("\\_", "\\$")) {
    pieces <- strsplit(x, token, fixed = TRUE)[[1L]]
    if (length(pieces) > 1L) {
      x <- paste(pieces, collapse = paste0(token, "\\allowbreak{}"))
    }
  }
  x
}

.mira_report_repeat_longtable_header_long <- function(x) {
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

.mira_report_table_long <- function(x, caption = NULL, cfg) {
  if (.mira_report_is_empty_long(x)) {
    cat("*No rows are available for display.*\n\n")
    return(invisible(NULL))
  }

  table_data <- tryCatch(
    .mira_report_prepare_table_long(x, digits = cfg$digits),
    error = function(e) NULL
  )
  if (is.null(table_data)) {
    .mira_report_text_block_long(x)
    return(invisible(NULL))
  }

  original_n <- nrow(table_data)
  row_limit <- if (is.finite(cfg$max_table_rows) &&
                   cfg$max_table_rows < .Machine$integer.max) {
    max(1L, as.integer(cfg$max_table_rows))
  } else {
    Inf
  }
  if (is.finite(row_limit) && original_n > row_limit) {
    table_data <- table_data[seq_len(row_limit), , drop = FALSE]
  }

  is_latex <- knitr::is_latex_output()
  is_html <- knitr::is_html_output()
  output_format <- if (is_latex) "latex" else if (is_html) "html" else "pipe"
  max_columns <- max(2L, as.integer(cfg$max_table_columns))
  if (is_latex) max_columns <- min(max_columns, 8L)
  if (!is_latex && !is_html) max_columns <- min(max_columns, 7L)
  width_budget <- if (is_latex) 76 else if (is_html) Inf else 68
  blocks <- .mira_report_partition_columns_long(
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
      common$align <- .mira_report_latex_alignment_long(block)
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
      .mira_report_callout_long(
        paste0(
          "The table could not be formatted by knitr (`",
          render_error,
          "`). Its contents are shown as text so that the report can continue."
        ),
        type = "warning", title = "Table formatting fallback"
      )
      .mira_report_text_block_long(block)
    } else {
      rendered_text <- as.character(rendered)
      if (is_latex) {
        rendered_text <- .mira_report_repeat_longtable_header_long(rendered_text)

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
        rendered_text <- .mira_report_latex_breaks_long(rendered_text)
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

  if (is.finite(row_limit) && original_n > row_limit) {
    cat(sprintf(
      "*Showing %d of %d rows because of the explicit `max_table_rows` limit.*\n\n",
      nrow(table_data), original_n
    ))
  }
  invisible(table_data)
}

.mira_report_prepare_namespace_long <- function(x) {
  package_map <- list(
    ggplot2 = c("ggplot", "ggplot2::ggplot"),
    lme4 = c("merMod", "lmerMod", "glmerMod", "summary.merMod"),
    lmerTest = c("lmerModLmerTest", "summary.lmerModLmerTest"),
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

.mira_report_text_block_long <- function(x) {
  .mira_report_prepare_namespace_long(x)
  text <- tryCatch(
    {
      if (inherits(x, c("lm", "glm", "merMod", "lmerMod", "glmerMod"))) {
        utils::capture.output(print(summary(x)))
      } else {
        utils::capture.output(print(x))
      }
    },
    error = function(e) paste("The object could not be printed:", conditionMessage(e))
  )
  text <- gsub("```", "'''", enc2utf8(text), fixed = TRUE)
  cat("```text\n", paste(text, collapse = "\n"), "\n```\n\n", sep = "")
}

.mira_report_callout_long <- function(text, type = "note", title = "Note") {
  allowed <- c("note", "tip", "warning", "important", "caution")
  if (!type %in% allowed) type <- "note"
  cat(sprintf("::: {.callout-%s appearance=\"simple\"}\n", type))
  cat("**", title, ".** ", text, "\n", sep = "")
  cat(":::\n\n")
}

.mira_report_is_plot_long <- function(x) {
  inherits(x, c("ggplot", "ggplot2::ggplot", "recordedplot", "grob", "gTree"))
}

.mira_report_plot_long <- function(x, name = "Figure") {
  .mira_report_prepare_namespace_long(x)
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
    .mira_report_callout_long(conditionMessage(e), type = "warning",
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

.mira_report_is_plain_list_long <- function(x) {
  is.list(x) && !is.data.frame(x) && !.mira_report_is_plot_long(x) &&
    !.mira_report_is_model_long(x)
}

.mira_report_is_model_long <- function(x) {
  inherits(x, c("lm", "glm", "merMod", "lme", "gls", "geeglm",
                "summary.geeglm", "emmGrid", "emm_list", "afex_aov", "htest"))
}

.mira_report_extract_outcomes_long <- function(result) {
  if (inherits(result, "mira_detect_long")) return(list())
  if (inherits(result, "mira_info_multi_long") &&
      (!is.list(result$outcomes) || !length(result$outcomes))) {
    return(list())
  }
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

.mira_report_section_long <- function(key) {
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

.mira_report_selected_long <- function(cfg, section) {
  "all" %in% cfg$sections || section %in% cfg$sections
}

.mira_report_node_has_result_long <- function(node) {
  if (is.null(node)) return(FALSE)
  if (node$kind %in% c("table", "plot", "model")) return(TRUE)
  if (node$kind == "group") {
    return(any(vapply(node$children, .mira_report_node_has_result_long, logical(1L))))
  }
  FALSE
}

.mira_report_prune_long <- function(x, key, path, cfg, depth = 0L,
                               parent_performed = FALSE) {
  if (depth > cfg$max_depth || is.null(x) || is.function(x) ||
      is.environment(x) || inherits(x, "mira_info_error_long")) return(NULL)
  # These fields describe the machinery or source data, not analytical output.
  if (key %in% c("call", "outcomes", "long_data", "reason_skipped",
                 "performed", "enabled", "package", "status", "object_path",
                 "plot_error", "disabled", "failed_outcomes")) return(NULL)
  if (identical(key, "data") &&
      grepl("(^|\\$)(model|advanced_tests|advanced_models|robustness)(\\$|$)",
            path)) return(NULL)
  if (grepl("(^|\\$)overview(\\$|$)", path) &&
      identical(key, "arm_analysis")) return(NULL)
  if (.mira_report_is_plot_long(x)) {
    return(list(kind = "plot", key = key, path = path, value = x))
  }
  if (is.data.frame(x) || is.matrix(x) || is.table(x)) {
    tab <- x
    if (key %in% c("pearson", "spearman") && all(is.na(tab))) return(NULL)
    # Show only models with a fitted object in the comparison table.
    if (identical(key, "model_comparison") && is.data.frame(tab) &&
        "object_path" %in% names(tab)) {
      tab <- tab[
        !is.na(tab$object_path) & nzchar(trimws(as.character(tab$object_path))),
        , drop = FALSE
      ]
    }
    dimensions <- if (is.table(tab)) dim(as.data.frame(tab)) else dim(tab)
    if (length(dimensions) < 2L || dimensions[2L] == 0L) return(NULL)
    if (dimensions[1L] == 0L && !parent_performed &&
        !grepl("(^|\\$)outliers(\\$|$)", path)) return(NULL)
    return(list(kind = "table", key = key, path = path, value = tab,
                original_class = paste(class(x), collapse = "/")))
  }
  if (.mira_report_is_model_long(x)) {
    return(list(kind = "model", key = key, path = path, value = x))
  }
  if (is.list(x)) {
    if (!length(x)) return(NULL)
    if ("performed" %in% names(x) && identical(x$performed, FALSE)) return(NULL)
    if (grepl("(^|\\$)multiplicity(\\$|$)", path) &&
        "table" %in% names(x) && is.null(x$table)) return(NULL)
    if (identical(key, "nlme") && "performed" %in% names(x) &&
        !isTRUE(x$performed)) return(NULL)
    if (identical(key, "arm_analysis") && !isTRUE(x$enabled)) return(NULL)
    if (identical(key, "model") && "fitted_model" %in% names(x) &&
        is.null(x$fitted_model)) return(NULL)
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
          is.list(item) && isTRUE(item$performed) &&
            (!is.null(item$model) || !is.null(item$object))
        }, logical(1L))]
        child_value <- child_value[child_value$model %in% fitted, , drop = FALSE]
      }
      child <- .mira_report_prune_long(child_value, child_key, paste0(path, "$", child_key),
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
    if ((failed || key %in% c("outliers", "plots", "multiplicity")) &&
        !any(vapply(children, .mira_report_node_has_result_long, logical(1L)))) return(NULL)
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

.mira_report_order_long <- function(keys) {
  scientific <- c("settings", "time_vars", "time_labels", "overview",
                  "descriptives", "missing", "change", "arm_analysis",
                  "correlations", "variability", "trajectories", "model",
                  "advanced_tests", "advanced_models", "robustness", "effect_sizes",
                  "multiplicity", "sensitivity", "outliers", "plots", "diagnostics")
  c(scientific[scientific %in% keys], keys[!keys %in% scientific])
}

.mira_report_outcome_tree_long <- function(outcome, cfg) {
  if (inherits(outcome, "mira_info_error_long") || !is.list(outcome)) return(list())
  skip <- c("call", "version", "config", "data_overview",
            "detected_variables", "outcome", "outcome_display",
            "outcomes", "long_data", "plot_error")
  keys <- .mira_report_order_long(setdiff(names(outcome), skip))
  nodes <- list()
  for (key in keys) {
    if (!.mira_report_selected_long(cfg, .mira_report_section_long(key))) next
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
    node <- .mira_report_prune_long(component, key, paste0("outcome$", key), cfg)
    if (!is.null(node)) nodes[[key]] <- node
  }
  if (.mira_report_selected_long(cfg, "diagnostics")) {
    diagnostic <- nodes$diagnostics
    if (is.null(diagnostic)) diagnostic <- list(kind = "group", key = "diagnostics",
                                               path = "outcome$diagnostics", children = list())
    model <- outcome$model
    if (is.list(model) && is.null(model$fitted_model)) {
      for (key in c("error", "warnings")) {
        node <- .mira_report_prune_long(model[[key]], paste0("model_", key),
                                   paste0("outcome$model$", key), cfg)
        if (!is.null(node)) diagnostic$children[[paste0("model_", key)]] <- node
      }
    }
    plot_error <- .mira_report_prune_long(outcome$plot_error, "figure_error",
                                     "outcome$plot_error", cfg)
    if (!is.null(plot_error)) diagnostic$children$figure_error <- plot_error
    if (length(diagnostic$children)) nodes$diagnostics <- diagnostic
  }
  nodes
}

.mira_report_context_tree_long <- function(result, cfg) {
  keys <- c(config = "methods", data_overview = "data_quality",
            detected_variables = "data_quality", diagnostics = "diagnostics")
  nodes <- list()
  for (key in names(keys)) {
    if (!.mira_report_selected_long(cfg, keys[[key]])) next
    # Outcome diagnostics are already shown inside each outcome for single output.
    if (key == "diagnostics" && !inherits(result, "mira_info_multi_long") &&
        !inherits(result, "mira_detect_long")) next
    if (is.null(result[[key]])) next
    node <- .mira_report_prune_long(result[[key]], key, paste0("result$", key), cfg)
    if (!is.null(node)) nodes[[key]] <- node
  }
  if (inherits(result, "mira_info_multi_long") &&
      .mira_report_selected_long(cfg, "diagnostics")) {
    errors <- list()
    for (name in names(result$outcomes)) {
      item <- result$outcomes[[name]]
      if (!inherits(item, "mira_info_error_long")) next
      node <- .mira_report_prune_long(item$error, name,
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

.mira_report_render_node_long <- function(node, level, cfg) {
  label <- .mira_report_human_name_long(node$key)
  if (node$kind == "value") {
    cat("**", label, ":** ", .mira_report_inline_long(node$value), "\n\n", sep = "")
    return(invisible(NULL))
  }
  .mira_report_heading_long(label, level)
  if (node$kind == "group") {
    for (child in node$children) .mira_report_render_node_long(child, min(6L, level + 1L), cfg)
  } else if (node$kind == "table") {
    rows <- if (is.table(node$value)) nrow(as.data.frame(node$value)) else nrow(node$value)
    if (rows == 0L) {
      if (grepl("(^|\\$)outliers(\\$|$)", node$path)) {
        cat("**Outliers detected:** 0.\n\n")
      } else {
        cat("**Observed results:** 0 rows.\n\n")
      }
    } else {
      .mira_report_table_long(node$value, caption = label, cfg = cfg)
    }
    for (attribute in c("note", "interpretation")) {
      value <- attr(node$value, attribute, exact = TRUE)
      if (!is.null(value) && length(value)) {
        cat("*", .mira_report_human_name_long(attribute), ":* ",
            .mira_report_inline_long(value), "\n\n", sep = "")
      }
    }
  } else if (node$kind == "plot") {
    .mira_report_plot_long(node$value, node$key)
  } else if (node$kind %in% c("model", "formula")) {
    .mira_report_text_block_long(node$value)
  }
  invisible(NULL)
}

.mira_report_render_abstract_long <- function(outcomes, cfg) {
  if (!.mira_report_selected_long(cfg, "abstract")) return(invisible(NULL))
  availability_cfg <- cfg
  availability_cfg$sections <- "all"
  available <- names(outcomes)[vapply(outcomes, function(x) {
    length(.mira_report_outcome_tree_long(x, availability_cfg)) > 0L
  }, logical(1L))]
  if (!length(available)) return(invisible(NULL))
  .mira_report_heading_long("Abstract", 1L)
  cat("This report presents the observed analytical results for ",
      paste(available, collapse = ", "), ".\n\n", sep = "")
  invisible(NULL)
}

.mira_report_render_conclusions_long <- function(outcomes, cfg) {
  if (!.mira_report_selected_long(cfg, "conclusions")) return(invisible(NULL))
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
  .mira_report_heading_long("Descriptive summary", 1L)
  .mira_report_table_long(do.call(rbind, rows),
                     caption = "Differences between observed means", cfg = cfg)
  invisible(NULL)
}

.mira_report_knit_long <- function(payload) {
  result <- payload$result
  cfg <- payload$report
  if (knitr::is_html_output()) .mira_report_html_style_long()
  outcomes <- .mira_report_extract_outcomes_long(result)
  context <- .mira_report_context_tree_long(result, cfg)
  .mira_report_render_abstract_long(outcomes, cfg)
  if (length(context)) {
    .mira_report_heading_long("Analysis context", 1L)
    for (node in context) .mira_report_render_node_long(node, 2L, cfg)
  }
  for (i in seq_along(outcomes)) {
    nodes <- .mira_report_outcome_tree_long(outcomes[[i]], cfg)
    if (!length(nodes)) next
    name <- names(outcomes)[i]
    display <- outcomes[[i]]$outcome_display
    if (is.null(display) || !length(display) || is.na(display[1L])) display <- name
    .mira_report_heading_long(paste0("Outcome: ", display[1L]), 1L)
    for (node in nodes) .mira_report_render_node_long(node, 2L, cfg)
  }
  .mira_report_render_conclusions_long(outcomes, cfg)
  extra <- payload$extra_objects
  if (.mira_report_selected_long(cfg, "appendix") && length(extra)) {
    extra_nodes <- lapply(names(extra), function(key) {
      .mira_report_prune_long(extra[[key]], key, paste0("extra_objects$", key), cfg)
    })
    extra_nodes <- Filter(Negate(is.null), extra_nodes)
    if (length(extra_nodes)) {
      .mira_report_heading_long("Supplementary objects", 1L)
      for (node in extra_nodes) .mira_report_render_node_long(node, 2L, cfg)
    }
  }
  invisible(NULL)
}

.mira_report_collect_tables_long <- function(result, cfg) {
  # Display filtering never changes the analytical archive on disk.
  cfg$sections <- "all"
  entries <- list()
  add_long <- function(node, outcome = "") {
    if (node$kind == "table") {
      entries[[length(entries) + 1L]] <<- list(outcome = outcome,
        path = node$path, value = node$value,
        original_class = node$original_class)
    } else if (node$kind == "group") {
      for (child in node$children) add_long(child, outcome)
    }
  }
  for (node in .mira_report_context_tree_long(result, cfg)) add_long(node)
  outcomes <- .mira_report_extract_outcomes_long(result)
  for (name in names(outcomes)) {
    for (node in .mira_report_outcome_tree_long(outcomes[[name]], cfg)) add_long(node, name)
  }
  entries
}

.mira_report_file_slug_long <- function(x) {
  x <- iconv(enc2utf8(as.character(x)[1L]), to = "ASCII//TRANSLIT", sub = "")
  x <- tolower(gsub("[^A-Za-z0-9]+", "_", x))
  x <- gsub("(^_+|_+$)", "", x)
  if (is.na(x) || !nzchar(x)) x <- "item"
  x <- substr(x, 1L, 90L)
  if (grepl("^(con|prn|aux|nul|com[1-9]|lpt[1-9])$", x)) x <- paste0("item_", x)
  x
}

.mira_report_csv_data_long <- function(x) {
  # A data frame can contain matrix-valued columns. write.csv() eventually
  # calls as.matrix.data.frame() on those columns, which can fail when their
  # expanded column names do not match the matrix width. Expand each matrix
  # column explicitly, retaining its values and keeping numeric p-values raw.
  if (is.table(x) || is.matrix(x)) {
    out <- .mira_report_prepare_table_long(x, format_p = FALSE)
    if (is.matrix(x) && ".row" %in% names(out)) {
      names(out)[names(out) == ".row"] <- "row"
    }
  } else {
    input <- as.data.frame(x, stringsAsFactors = FALSE, check.names = FALSE)
    pieces <- list()
    for (i in seq_along(input)) {
      key <- names(input)[i]
      if (is.na(key) || !nzchar(key)) key <- paste0("column_", i)
      column <- input[[i]]
      dimensions <- dim(column)
      if (length(dimensions) == 2L && dimensions[1L] == nrow(input) &&
          dimensions[2L] > 0L) {
        labels <- colnames(column)
        if (is.null(labels)) labels <- paste0("column_", seq_len(dimensions[2L]))
        for (j in seq_len(dimensions[2L])) {
          name <- paste0(key, "_", labels[j])
          value <- column[, j, drop = TRUE]
          if (inherits(value, "AsIs")) class(value) <- NULL
          if (!is.atomic(value) || !is.null(dim(value))) {
            value <- vapply(seq_len(nrow(input)), function(row) {
              .mira_report_inline_long(column[row, j, drop = TRUE])
            }, character(1L))
          }
          pieces[[length(pieces) + 1L]] <- value
          names(pieces)[length(pieces)] <- name
        }
      } else if (is.list(column) || !is.null(dimensions)) {
        normalized <- .mira_report_prepare_table_long(input[i], format_p = FALSE)
        pieces[[length(pieces) + 1L]] <- normalized[[1L]]
        names(pieces)[length(pieces)] <- key
      } else {
        pieces[[length(pieces) + 1L]] <- column
        names(pieces)[length(pieces)] <- key
      }
    }
    names(pieces) <- make.unique(names(pieces), sep = "_")
    out <- as.data.frame(pieces, stringsAsFactors = FALSE, check.names = FALSE)
  }
  names(out) <- make.unique(names(out), sep = "_")
  valid <- vapply(out, function(column) {
    is.atomic(column) && is.null(dim(column)) && length(column) == nrow(out)
  }, logical(1L))
  if (!all(valid)) stop("A result table contains columns that cannot be written to CSV.",
                        call. = FALSE)
  out
}

.mira_report_export_plan_long <- function(result, cfg, output_dir) {
  tables <- .mira_report_collect_tables_long(result, cfg)
  used <- character(0)
  for (i in seq_along(tables)) {
    path_keys <- strsplit(sub("^([^$]+)\\$", "", tables[[i]]$path),
                          "$", fixed = TRUE)[[1L]]
    base <- paste(vapply(path_keys, .mira_report_file_slug_long, character(1L)), collapse = "_")
    if (!nzchar(base)) base <- "result"
    base <- sub("_+$", "", substr(base, 1L, 90L))
    dir_part <- if (nzchar(tables[[i]]$outcome) &&
                    length(.mira_report_extract_outcomes_long(result)) > 1L) {
      substr(.mira_report_file_slug_long(tables[[i]]$outcome), 1L, 60L)
    } else ""
    relative <- if (nzchar(dir_part)) {
      file.path("results", dir_part, paste0(base, ".csv"))
    } else {
      file.path("results", paste0(base, ".csv"))
    }
    candidate <- relative
    serial <- 2L
    while (tolower(candidate) %in% used) {
      candidate <- file.path("results", dir_part, paste0(base, "_", serial, ".csv"))
      serial <- serial + 1L
    }
    used <- c(used, tolower(candidate))
    tables[[i]]$file <- gsub("\\\\", "/", candidate)
    tables[[i]]$absolute <- file.path(output_dir, candidate)
    tables[[i]]$csv <- .mira_report_csv_data_long(tables[[i]]$value)
  }
  tables
}

.mira_report_manifest_long <- function(plan) {
  if (!length(plan)) return(data.frame(outcome = character(), path = character(),
    file = character(), original_class = character(), rows = integer(),
    columns = integer(), stringsAsFactors = FALSE))
  do.call(rbind, lapply(plan, function(entry) data.frame(
    outcome = entry$outcome, path = entry$path, file = entry$file,
    original_class = entry$original_class, rows = nrow(entry$csv),
    columns = ncol(entry$csv), stringsAsFactors = FALSE)))
}

.mira_report_html_style_long <- function() {
  cat("```{=html}\n<style>\nbody{line-height:1.58;color:#14213d} h1,h2,h3{color:#14213d} .mira-table-wrap{overflow-x:auto} .mira-table{width:100%;font-size:.88rem}\n</style>\n```\n\n")
}

.mira_report_runtime_names_long <- function() {
  c(".mira_report_or_long", ".mira_report_is_empty_long", ".mira_report_heading_long",
    ".mira_report_human_name_long", ".mira_report_format_p_long",
    ".mira_report_format_number_long", ".mira_report_inline_long",
    ".mira_report_prepare_table_long", ".mira_report_numeric_column_long",
    ".mira_report_column_widths_long", ".mira_report_partition_columns_long",
    ".mira_report_latex_alignment_long", ".mira_report_latex_breaks_long",
    ".mira_report_repeat_longtable_header_long", ".mira_report_table_long",
    ".mira_report_prepare_namespace_long", ".mira_report_text_block_long",
    ".mira_report_callout_long", ".mira_report_is_plot_long", ".mira_report_plot_long",
    ".mira_report_is_plain_list_long", ".mira_report_is_model_long",
    ".mira_report_extract_outcomes_long", ".mira_report_section_long",
    ".mira_report_selected_long", ".mira_report_node_has_result_long",
    ".mira_report_prune_long", ".mira_report_order_long",
    ".mira_report_outcome_tree_long", ".mira_report_context_tree_long",
    ".mira_report_render_node_long", ".mira_report_render_abstract_long",
    ".mira_report_render_conclusions_long", ".mira_report_knit_long",
    ".mira_report_html_style_long")
}
.mira_report_slug_long <- function(x) {
  x <- iconv(enc2utf8(as.character(x)[1L]), to = "ASCII//TRANSLIT", sub = "")
  x <- tolower(gsub("[^A-Za-z0-9]+", "_", x))
  x <- gsub("(^_+|_+$)", "", x)
  if (is.na(x) || !nzchar(x)) x <- "mira_report"
  x <- substr(x, 1L, 90L)
  if (grepl("^(con|prn|aux|nul|com[1-9]|lpt[1-9])$", x)) {
    x <- paste0("mira_", x)
  }
  x
}

.mira_report_path_is_within_long <- function(path, root) {
  if (!is.character(path) || length(path) != 1L || is.na(path) ||
      !file.exists(path) || !dir.exists(root)) return(FALSE)
  normalized_path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  normalized_root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  if (.Platform$OS.type == "windows") {
    normalized_path <- tolower(normalized_path)
    normalized_root <- tolower(normalized_root)
  }
  identical(normalized_path, normalized_root) ||
    startsWith(normalized_path, paste0(sub("/+$", "", normalized_root), "/"))
}

.mira_report_safe_stale_files_long <- function(paths, output_dir) {
  paths <- paths[file.exists(paths)]
  if (!length(paths)) return(character())
  results_dir <- file.path(output_dir, "results")
  if (!.mira_report_path_is_within_long(results_dir, output_dir)) return(character())
  paths[vapply(paths, .mira_report_path_is_within_long, logical(1L),
               root = results_dir)]
}

.mira_report_yaml_quote_long <- function(x) {
  x <- gsub("[\r\n]+", " ", as.character(x)[1L])
  paste0("'", gsub("'", "''", x, fixed = TRUE), "'")
}

.mira_report_r_quote_long <- function(x) {
  encodeString(as.character(x)[1L], quote = "\"")
}

.mira_report_format_lines_long <- function(formats) {
  lines <- "format:"
  if ("html" %in% formats) {
    lines <- c(
      lines,
      "  html:",
      "    theme: cosmo",
      "    toc: true",
      "    toc-depth: 3",
      "    toc-location: left",
      "    number-sections: true",
      "    embed-resources: true",
      "    smooth-scroll: true",
      "    fig-cap-location: bottom",
      "    tbl-cap-location: top"
    )
  }
  if ("pdf" %in% formats) {
    lines <- c(
      lines,
      "  pdf:",
      "    documentclass: article",
      "    papersize: a4",
      "    geometry: margin=25mm",
      "    toc: true",
      "    toc-depth: 3",
      "    number-sections: true",
      "    number-depth: 3",
      "    colorlinks: true",
      "    fig-cap-location: bottom",
      "    tbl-cap-location: top"
    )
  }
  if ("docx" %in% formats) {
    lines <- c(
      lines,
      "  docx:",
      "    toc: true",
      "    toc-depth: 3",
      "    number-sections: true",
      "    fig-cap-location: bottom",
      "    tbl-cap-location: top"
    )
  }
  lines
}

.mira_report_qmd_lines_long <- function(title, subtitle, author, date, formats,
                                   runtime_file, payload_file) {
  yaml <- c("---", paste0("title: ", .mira_report_yaml_quote_long(title)))
  if (!is.null(subtitle) && nzchar(subtitle)) {
    yaml <- c(yaml, paste0("subtitle: ", .mira_report_yaml_quote_long(subtitle)))
  }
  if (!is.null(author) && length(author) > 0L) {
    if (length(author) == 1L) {
      yaml <- c(yaml, paste0("author: ", .mira_report_yaml_quote_long(author)))
    } else {
      yaml <- c(yaml, "author:", paste0("  - ", vapply(
        author, .mira_report_yaml_quote_long, character(1L)
      )))
    }
  }
  yaml <- c(
    yaml,
    paste0("date: ", .mira_report_yaml_quote_long(as.character(date))),
    "lang: en",
    .mira_report_format_lines_long(formats),
    "execute:",
    "  echo: false",
    "  warning: false",
    "  message: false",
    "  error: false",
    "  freeze: false",
    "---",
    "",
    "```{r}",
    "#| label: setup",
    "#| include: false",
    paste0("source(", .mira_report_r_quote_long(runtime_file),
           ", local = knitr::knit_global())"),
    paste0("payload <- readRDS(", .mira_report_r_quote_long(payload_file), ")"),
    "mira_plot_device <- if (knitr::is_latex_output() && capabilities('cairo')) 'cairo_pdf' else 'png'",
    "knitr::opts_chunk$set(fig.width = 7.0, fig.height = 4.6, dpi = 320,",
    "                      dev = mira_plot_device, fig.align = 'center',",
    "                      out.width = '100%')",
    "```",
    "",
    "```{r}",
    "#| label: mira-report",
    "#| results: asis",
    ".mira_report_knit_long(payload)",
    "```",
    ""
  )
  yaml
}

.mira_report_validate_scalar_flag_long <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop(sprintf("%s must be TRUE or FALSE.", name), call. = FALSE)
  }
  invisible(TRUE)
}

#' Build a reproducible Quarto report from MIRA frequentist results
#'
#' Creates a publication-oriented report from an existing \code{mira_info_long},
#' \code{mira_info_multi_long}, or \code{mira_detect_long} object, or first runs
#' [mira_info_long()] on supplied wide-format data. In every mode the function writes
#' a reproducibility bundle containing Quarto source, a serialized payload with
#' the original MIRA object, the report runtime, a machine-readable table
#' archive, and a manifest. It can then render any combination of HTML, PDF,
#' and DOCX output.
#'
#' @param x \code{NULL}; an object inheriting from \code{mira_info_long},
#'   \code{mira_info_multi_long}, or \code{mira_detect_long}; or a data frame. A data frame
#'   supplied through \code{x} is treated exactly as if it had been supplied
#'   through \code{data}. Do not supply both \code{x} and \code{data}. A
#'   standalone \code{mira_info_error_long} object is not a valid input, although
#'   failed outcomes nested inside a valid multi-outcome result are supported.
#' @param data An optional data frame to analyse before building the report.
#'   When supplied, the function evaluates
#'   \code{mira_info_long(data = data, ...)} and reports the returned object. See
#'   [mira_info_long()] for the required wide longitudinal structure and the complete
#'   analysis interface.
#' @param ... Arguments passed to [mira_info_long()] only when \code{data} is used.
#'   Typical arguments include \code{id}, \code{outcomes}, \code{time_vars},
#'   \code{time_labels}, \code{arm}, \code{covariates}, \code{analyses}, and
#'   \code{p_adjust_method}. If \code{verbose} is not supplied, it is set to
#'   \code{FALSE}. Supplying \code{...} with a precomputed \code{x}, or
#'   repeating \code{data} inside \code{...}, is an error.
#' @param output_dir A non-empty, length-one path to the bundle directory.
#'   Missing directories are created recursively. The default is
#'   \code{file.path(getwd(), "mira_analyses", "mira_report_flong")}. The
#'   normalized absolute path is available in the returned
#'   \code{output_dir} component.
#' @param output_file A non-empty character scalar used as a filename stem; do
#'   not include an output-format extension. For portability, the value is
#'   transliterated to ASCII, converted to lower case, reduced to letters,
#'   numbers, and underscores, and limited to 90 characters. Empty normalized
#'   stems fall back to \code{"mira_report"}; Windows-reserved stems receive a
#'   \code{"mira_"} prefix. For example, \code{"Primary Report.html"} produces
#'   the stem \code{"primary_report_html"}.
#' @param format A non-empty character vector containing one or more of
#'   \code{"html"}, \code{"pdf"}, \code{"docx"}, or \code{"all"}. Matching is
#'   case-insensitive and duplicates are removed while preserving order.
#'   \code{"all"} expands to HTML, PDF, and DOCX, in that order. With
#'   \code{render = FALSE}, this argument still controls the formats declared in
#'   the Quarto source and the paths listed in \code{expected_files}.
#' @param title A non-empty character scalar used as the Quarto document title.
#' @param subtitle \code{NULL} or a non-missing character scalar. \code{NULL}
#'   and \code{""} omit the subtitle from the Quarto metadata.
#' @param author \code{NULL} or a character vector of non-missing, non-blank
#'   author names. One name is written as scalar Quarto metadata; multiple names
#'   are written as a YAML list. \code{character(0)} omits the field.
#' @param date A non-missing object of length one that can be represented by
#'   \code{as.character()}, normally a \code{Date}, date-time, or character
#'   value. The resulting text is written to the Quarto metadata.
#' @param sections A non-empty character vector selecting visible report
#'   content. Matching is case-insensitive and duplicate values are removed.
#'   Use \code{"all"} or any combination of \code{"abstract"},
#'   \code{"methods"}, \code{"data_quality"}, \code{"overview"},
#'   \code{"descriptives"}, \code{"frequencies"}, \code{"missingness"},
#'   \code{"change"}, \code{"arms"}, \code{"correlations"},
#'   \code{"variability"}, \code{"models"}, \code{"robustness"},
#'   \code{"outliers"}, \code{"trajectories"}, \code{"figures"},
#'   \code{"diagnostics"}, \code{"conclusions"}, and \code{"appendix"}.
#'   Selection affects only the visible report and does not reorder it: the
#'   payload remains complete, and all eligible MIRA result tables are still
#'   exported to CSV.
#' @param include_complete_output A non-missing logical scalar retained for API
#'   compatibility. The current renderer already traverses all eligible
#'   scientific results without adding a duplicate complete-output appendix, so
#'   this value does not change visible content. It is recorded in
#'   \code{report_config}; the full MIRA object is always preserved in the
#'   payload, whether this argument is \code{TRUE} or \code{FALSE}.
#' @param extra_objects A list of supplementary objects, such as fitted models,
#'   to preserve in the payload and optionally print under the appendix.
#'   Prefer a named list, for example
#'   \code{list(adjusted_model = fitted_model)}. Empty or missing names become
#'   \code{extra_<position>} and duplicate names are made unique. These objects
#'   are visible only when \code{"appendix"} or \code{"all"} is selected and
#'   are not included in the CSV archive.
#' @param max_table_rows A positive numeric scalar or \code{Inf}, controlling
#'   the maximum number of rows displayed in each report table. The default,
#'   \code{Inf}, displays every row. A finite value is converted to an effective
#'   whole-row limit during rendering, with at least one row shown. Truncation
#'   is explicitly noted below the displayed table and never truncates the CSV
#'   export.
#' @param max_table_columns A finite numeric scalar greater than or equal to
#'   two, converted to an integer for rendering. It is an upper bound per
#'   displayed table block, not a request to discard columns. Wide tables are
#'   split into blocks and the first one or two columns, typically identifiers,
#'   are repeated. Width heuristics can produce smaller blocks; PDF blocks are
#'   capped at eight columns and non-HTML, non-LaTeX blocks at seven.
#' @param digits A finite numeric scalar from 1 through 10, converted to an
#'   integer. It controls numeric display precision in tables. Recognized
#'   numeric p-value columns below \code{10^-digits} are displayed as less than
#'   that threshold. CSV values are exported without this display formatting.
#' @param max_depth A finite numeric scalar greater than or equal to one,
#'   converted to an integer. This is a safety limit for recursive traversal of
#'   nested result and supplementary objects. Deeper content is omitted from
#'   the visible report. For MIRA results, the limit also applies to CSV
#'   traversal; supplementary objects are never exported to CSV. All original
#'   content remains in the payload.
#' @param render A non-missing logical scalar. If \code{TRUE}, render every
#'   requested format with [quarto::quarto_render()]. If \code{FALSE}, skip
#'   Quarto but still write the QMD source, runtime script, payload RDS,
#'   eligible-table CSV archive, and manifest. This mode is useful for
#'   inspection, testing, deferred rendering, and environments without the
#'   Quarto CLI.
#' @param open A non-missing logical scalar. If \code{TRUE} and at least one
#'   requested format renders successfully, open the first successful file in
#'   normalized request order using [utils::browseURL()]. It has no effect when
#'   \code{render = FALSE} or every render fails. Browser-opening errors are
#'   ignored. The default is \code{interactive()}.
#' @param overwrite A non-missing logical scalar. If \code{FALSE}, the call
#'   stops before writing bundle files when a protected source, requested
#'   rendered output, manifest, or planned CSV file already exists. If
#'   \code{TRUE}, those artifacts may be replaced. Obsolete CSV files listed in
#'   a previous valid manifest are removed only when their paths are safely
#'   contained below this bundle's \code{results} directory; unrelated files
#'   are not cleaned.
#' @param quiet A non-missing logical scalar passed to
#'   [quarto::quarto_render()] for each requested format. It controls Quarto
#'   console output, not validation warnings or the descriptor printed by this
#'   function.
#'
#' @details
#' \code{mira_report_freq_long()} has two input workflows. With a precomputed MIRA
#' object, it reports that object without refitting or recomputing its analyses.
#' With \code{data}, it first calls [mira_info_long()] and then reports the result.
#' The supplied or newly created MIRA object is serialized unchanged in the
#' payload. A \code{mira_detect_long} input produces analysis-context content only;
#' it has no outcome analyses to render.
#'
#' The function validates and normalizes the requested metadata, sections, and
#' formats; creates the output directory; plans all paths and CSV exports;
#' protects existing artifacts according to \code{overwrite}; writes the
#' reproducibility bundle; and, if requested, invokes Quarto once per format.
#' It prints a compact artifact summary before returning invisibly.
#'
#' The visible renderer is deliberately selective. It presents meaningful
#' analytical results while omitting machinery and source-data fields, disabled
#' or unperformed modules, empty non-result structures, all-missing correlation
#' matrices, duplicate equivalent model tables, and models without a fitted
#' object. Model and plotting failures can instead appear under diagnostics.
#' This presentation policy does not modify the original MIRA object and does
#' not remove anything from the serialized payload.
#'
#' Tables are adapted to the target medium. HTML uses horizontally scrollable
#' blocks, PDF uses repeated-header long tables with width-aware columns, and
#' DOCX uses portable pipe-table output. If \code{knitr} cannot format an
#' individual table, its contents fall back to a text block. If a stored plot
#' cannot be drawn, the report receives a warning callout and continues.
#'
#' @section Report sections:
#' Sections are emitted only when the corresponding eligible content exists.
#' Selecting \code{"all"} enables every route below.
#'
#' \describe{
#'   \item{\code{abstract}}{A generated statement naming outcomes with
#'     reportable analytical content.}
#'   \item{\code{methods}}{Resolved configuration and outcome settings,
#'     including source time variables and display labels.}
#'   \item{\code{data_quality}}{Dataset overview and detected-variable
#'     information.}
#'   \item{\code{overview}}{Outcome-level sample, completeness, and design
#'     summaries.}
#'   \item{\code{descriptives}}{Per-timepoint descriptive statistics.}
#'   \item{\code{frequencies}}{A result component named
#'     \code{frequencies}, when present; standard MIRA results need not contain
#'     such a component.}
#'   \item{\code{missingness}}{Availability summaries by timepoint and
#'     participant.}
#'   \item{\code{change}}{Within-participant change summaries and tests.}
#'   \item{\code{arms}}{Enabled treatment-arm descriptives and comparisons.}
#'   \item{\code{correlations}}{Estimable Pearson and Spearman matrices and
#'     pairwise sample sizes.}
#'   \item{\code{variability}}{Within- and between-participant variability and
#'     ICC results.}
#'   \item{\code{models}}{The primary model, advanced tests, advanced models,
#'     fitted objects, summaries, and model comparisons that have results.}
#'   \item{\code{robustness}}{Robust inference, effect sizes, multiplicity
#'     families, and sensitivity results.}
#'   \item{\code{outliers}}{Eligible IQR-based diagnostic outlier tables.}
#'   \item{\code{trajectories}}{Participant-level trajectory summaries.}
#'   \item{\code{figures}}{Stored MIRA plots that can be drawn.}
#'   \item{\code{diagnostics}}{Configuration diagnostics, relevant model and
#'     plotting errors or warnings, and failed outcomes from multi-outcome
#'     analyses.}
#'   \item{\code{conclusions}}{A generated descriptive table comparing the
#'     first and last observed means, when there are at least two descriptive
#'     rows and the means in the first and last rows are finite. It is a
#'     descriptive summary, not an inferential or causal conclusion.}
#'   \item{\code{appendix}}{Eligible objects supplied through
#'     \code{extra_objects}.}
#' }
#'
#' @section Reproducibility bundle:
#' For normalized stem \code{<stem>}, the bundle has the following layout:
#'
#'     <output_dir>/
#'       <stem>.qmd
#'       <stem>-runtime.R
#'       <stem>-payload.rds
#'       results_manifest.csv
#'       results/
#'         *.csv
#'         <outcome>/*.csv
#'       <stem>.html
#'       <stem>.pdf
#'       <stem>.docx
#'
#' The rendered paths at the bottom of the tree are targets, not a guarantee of
#' current success. Use \code{rendered} and \code{files} as the authoritative
#' status. In particular, with \code{overwrite = TRUE}, a pre-existing rendered
#' file can remain on disk after a skipped or failed render. The QMD uses English
#' document language, hides code and ordinary messages, installs a table of
#' contents, and numbers sections. HTML uses the Cosmo theme and embeds
#' resources; PDF uses an A4 article layout with 25-mm margins; DOCX uses
#' Quarto's standard document output.
#'
#' The runtime script contains the report helpers needed by the QMD. The payload
#' is an RDS list with exactly three top-level components:
#' \code{result}, \code{extra_objects}, and \code{report}. The last component
#' records normalized section and table settings plus a UTC creation timestamp.
#' Keep the QMD, runtime, and payload together when moving or rerendering the
#' source. For example, a deferred render can use
#' \code{quarto::quarto_render(report$source, execute_dir =
#' report$output_dir)}.
#' Deferred rendering still requires \pkg{knitr}, the Quarto CLI, and any R
#' packages needed to print or draw stored model and plot classes. The runtime
#' script supplies the MIRA-specific report helpers, not those dependencies.
#'
#' Because \code{results_manifest.csv} and the \code{results} directory have
#' bundle-level names, use a separate \code{output_dir} for each independently
#' managed report even when the \code{output_file} stems differ.
#'
#' @section CSV archive and manifest:
#' Every eligible tabular result found across all MIRA report sections is
#' exported independently of \code{sections}, \code{max_table_rows}, and
#' \code{max_table_columns}. Consequently, presentation choices cannot silently
#' truncate the analytical archive. P-values remain raw numeric values in CSV
#' files. Matrix-valued columns are expanded where possible, list-valued cells
#' are converted to plain-text representations, and file slugs plus numeric
#' suffixes prevent case-insensitive name collisions. Multi-outcome results are
#' partitioned into outcome-specific subdirectories.
#'
#' \code{results_manifest.csv} contains one row per exported table and exactly
#' these columns:
#'
#' \describe{
#'   \item{\code{outcome}}{Outcome key, or an empty value for shared context.}
#'   \item{\code{path}}{Traversal path of the table inside the MIRA object.}
#'   \item{\code{file}}{Portable path to the CSV, relative to
#'     \code{output_dir}.}
#'   \item{\code{original_class}}{Slash-separated class vector of the source
#'     table or matrix.}
#'   \item{\code{rows}, \code{columns}}{Dimensions of the exported CSV data.}
#' }
#'
#' Supplementary \code{extra_objects} are stored in the payload and can be
#' printed in the appendix, but they are not traversed for CSV export.
#'
#' @section Rendering requirements and failures:
#' Rendering requires the R packages \pkg{knitr} and \pkg{quarto} and a working
#' Quarto CLI. PDF additionally requires a compatible TeX installation. HTML
#' is generally the least dependency-intensive rendered target. DOCX rendering
#' does not require Microsoft Word, but opening or editing the result normally
#' requires compatible document software.
#'
#' Missing rendering packages cause an error after source artifacts have been
#' written. Once Quarto rendering begins, each requested format is attempted
#' independently. A per-format failure is captured in \code{render_errors},
#' contributes to one aggregate warning, and does not discard successful
#' formats or any source files. A Quarto call that returns without creating its
#' expected output is also treated as a failure. Because the bundle spans
#' multiple files, writing is not transactional: an unexpected filesystem
#' failure can leave a partial bundle for inspection or removal.
#'
#' @section Data confidentiality:
#' Treat \code{output_dir} as analytical data, not merely as a presentation
#' folder. The payload can contain longitudinal source data, participant
#' identifiers, fitted objects, diagnostics, and skipped-module reasons. CSV
#' files can also contain participant-level identifiers or diagnostic rows.
#' Review the complete bundle before sharing it and apply the same access,
#' retention, and de-identification controls used for the source analysis.
#'
#' @return Invisibly returns a named list with classes
#'   \code{c("mira_report_freq_long", "list")}. The object is also printed
#'   automatically. Its components are:
#'
#' \describe{
#'   \item{\code{call}}{The matched call.}
#'   \item{\code{mira_result}}{The supplied MIRA object or the object created by
#'     the internal [mira_info_long()] call.}
#'   \item{\code{formats}}{Normalized requested formats.}
#'   \item{\code{rendered}}{A named logical vector indicating whether each
#'     requested format was successfully created. Values are \code{FALSE} when
#'     rendering was skipped.}
#'   \item{\code{render_errors}}{A named character vector. A failed format has
#'     its captured message; successful or deliberately unrendered formats have
#'     \code{NA}.}
#'   \item{\code{files}}{Named absolute paths to successfully rendered files
#'     only. It has length zero when no format was rendered successfully.}
#'   \item{\code{expected_files}}{Named absolute target paths for every
#'     requested format, whether or not the files exist.}
#'   \item{\code{source}}{Absolute path to the generated QMD file.}
#'   \item{\code{payload}}{Absolute path to the compressed payload RDS file.}
#'   \item{\code{runtime}}{Absolute path to the generated R runtime script.}
#'   \item{\code{results_manifest}}{Absolute path to
#'     \code{results_manifest.csv}.}
#'   \item{\code{result_files}}{Absolute paths to all CSV files planned and
#'     written for the current bundle.}
#'   \item{\code{output_dir}}{Normalized absolute bundle-directory path.}
#'   \item{\code{report_config}}{Normalized section and display settings and
#'     the UTC bundle-creation timestamp.}
#' }
#'
#' The \code{print.mira_report_freq_long()} method reports the output directory,
#' source filename, successful formats, and formats not rendered successfully,
#' and returns its argument invisibly.
#'
#' @seealso [mira_info_long()] for the analysis that supplies report content,
#'   [mira_detect_long()] for configuration-only input, [quarto::quarto_render()] for
#'   manual rendering, [base::readRDS()] for inspecting the payload, and
#'   [utils::read.csv()] for inspecting the manifest and table archive.
#'
#' @examples
#' # ------------------------------------------------------------------
#' # Example 1: create and inspect a complete source bundle without Quarto
#' # ------------------------------------------------------------------
#' set.seed(202602)
#' n_subjects <- 24L
#' baseline <- rnorm(n_subjects, mean = 50, sd = 7)
#' report_data <- data.frame(
#'   subject_id = sprintf("S%03d", seq_len(n_subjects)),
#'   score_t0 = baseline,
#'   score_t1 = baseline - 1.5 + rnorm(n_subjects, sd = 2),
#'   score_t2 = baseline - 3.0 + rnorm(n_subjects, sd = 2)
#' )
#'
#' analysis <- mira_info_long(
#'   report_data,
#'   analyses = "none",
#'   verbose = FALSE
#' )
#'
#' report_dir <- tempfile("mira-report-")
#' report <- mira_report_freq_long(
#'   analysis,
#'   output_dir = report_dir,
#'   output_file = "Example Report",
#'   format = "html",
#'   sections = c(
#'     "abstract", "overview", "descriptives", "change", "conclusions"
#'   ),
#'   render = FALSE,
#'   open = FALSE
#' )
#'
#' # The stem is normalized, and all reproducibility sources exist even though
#' # the HTML file was not rendered.
#' basename(report$source)
#' file.exists(c(
#'   report$source, report$payload, report$runtime, report$results_manifest
#' ))
#' report$rendered
#' report$expected_files
#'
#' # Section selection controls the document, not the machine-readable archive.
#' manifest <- utils::read.csv(
#'   report$results_manifest,
#'   stringsAsFactors = FALSE
#' )
#' manifest[, c("outcome", "path", "file", "rows", "columns")]
#'
#' # The payload preserves the original result rather than a display-truncated
#' # copy, so exact objects remain available for audit or deferred rendering.
#' payload <- readRDS(report$payload)
#' identical(payload$result$descriptives, analysis$descriptives)
#' unlink(report_dir, recursive = TRUE)
#'
#' \dontrun{
#' # ------------------------------------------------------------------
#' # Example 2: analyse raw data and build the report in one call
#' # ------------------------------------------------------------------
#' direct_report <- mira_report_freq_long(
#'   data = report_data,
#'   id = "subject_id",
#'   outcomes = "score",
#'   analyses = "none",
#'   output_dir = tempfile("mira-direct-"),
#'   output_file = "direct_analysis",
#'   format = "html",
#'   render = TRUE,
#'   open = TRUE
#' )
#'
#' # ------------------------------------------------------------------
#' # Example 3: request every output format with publication metadata
#' # ------------------------------------------------------------------
#' all_formats <- mira_report_freq_long(
#'   analysis,
#'   output_dir = tempfile("mira-all-formats-"),
#'   output_file = "primary_longitudinal_analysis",
#'   format = "all",
#'   title = "Primary Longitudinal Analysis",
#'   subtitle = "Frequentist analysis set",
#'   author = c("A. Analyst", "B. Statistician"),
#'   date = as.Date("2026-09-15"),
#'   render = TRUE,
#'   open = FALSE
#' )
#' all_formats$rendered
#' all_formats$render_errors
#' all_formats$files
#'
#' # ------------------------------------------------------------------
#' # Example 4: make a focused report while retaining every table in CSV
#' # ------------------------------------------------------------------
#' focused <- mira_report_freq_long(
#'   analysis,
#'   output_dir = tempfile("mira-focused-"),
#'   format = "html",
#'   sections = c(
#'     "methods", "data_quality", "descriptives", "missingness",
#'     "models", "diagnostics"
#'   ),
#'   max_table_rows = 50,
#'   max_table_columns = 8,
#'   digits = 4,
#'   render = FALSE,
#'   open = FALSE
#' )
#' focused_manifest <- utils::read.csv(focused$results_manifest)
#' file.exists(file.path(focused$output_dir, focused_manifest$file))
#'
#' # ------------------------------------------------------------------
#' # Example 5: report multiple outcomes
#' # ------------------------------------------------------------------
#' multi_data <- transform(
#'   report_data,
#'   fatigue_t0 = 30 + rnorm(n_subjects, sd = 5),
#'   fatigue_t1 = 28 + rnorm(n_subjects, sd = 5),
#'   fatigue_t2 = 25 + rnorm(n_subjects, sd = 5)
#' )
#' multi_analysis <- mira_info_long(
#'   multi_data,
#'   outcomes = c("score", "fatigue"),
#'   analyses = "none",
#'   verbose = FALSE
#' )
#' multi_report <- mira_report_freq_long(
#'   multi_analysis,
#'   output_dir = tempfile("mira-multi-"),
#'   format = "html",
#'   render = FALSE,
#'   open = FALSE
#' )
#' multi_manifest <- utils::read.csv(multi_report$results_manifest)
#' split(multi_manifest$file[nzchar(multi_manifest$outcome)],
#'       multi_manifest$outcome[nzchar(multi_manifest$outcome)])
#'
#' # ------------------------------------------------------------------
#' # Example 6: preserve a fitted object for a future appendix
#' # ------------------------------------------------------------------
#' auxiliary_model <- stats::lm(score_t2 ~ score_t0, data = report_data)
#' with_appendix <- mira_report_freq_long(
#'   analysis,
#'   output_dir = tempfile("mira-appendix-"),
#'   sections = c("descriptives", "change", "appendix"),
#'   extra_objects = list(baseline_adjusted_model = auxiliary_model),
#'   format = "html",
#'   render = FALSE,
#'   open = FALSE
#' )
#' names(readRDS(with_appendix$payload)$extra_objects)
#'
#' # ------------------------------------------------------------------
#' # Example 7: document detection decisions without running analyses
#' # ------------------------------------------------------------------
#' detected <- mira_detect_long(report_data, verbose = FALSE)
#' detection_report <- mira_report_freq_long(
#'   detected,
#'   output_dir = tempfile("mira-detection-"),
#'   sections = c("methods", "data_quality", "diagnostics"),
#'   format = "html",
#'   render = FALSE,
#'   open = FALSE
#' )
#'
#' # ------------------------------------------------------------------
#' # Example 8: deliberately replace a previously generated bundle
#' # ------------------------------------------------------------------
#' stable_dir <- tempfile("mira-stable-report-")
#' first_report <- mira_report_freq_long(
#'   analysis,
#'   output_dir = stable_dir,
#'   render = FALSE,
#'   open = FALSE
#' )
#' updated_report <- mira_report_freq_long(
#'   analysis,
#'   output_dir = stable_dir,
#'   sections = c("overview", "descriptives", "change"),
#'   render = FALSE,
#'   open = FALSE,
#'   overwrite = TRUE
#' )
#'
#' # A source-only bundle can also be rendered later while keeping execution
#' # relative to its companion runtime and payload files.
#' quarto::quarto_render(
#'   input = updated_report$source,
#'   output_format = "html",
#'   execute_dir = updated_report$output_dir
#' )
#' }
#' @export
mira_report_freq_long <- function(
    x = NULL,
    data = NULL,
    ...,
    output_dir = NULL,
    output_file = "mira_report",
    format = c("html", "pdf"),
    title = "MIRA: Longitudinal Analysis",
    subtitle = "Dynamic and reproducible statistical report",
    author = NULL,
    date = Sys.Date(),
    sections = "all",
    include_complete_output = TRUE,
    extra_objects = list(),
    max_table_rows = Inf,
    max_table_columns = 12L,
    digits = 3L,
    max_depth = 50L,
    render = TRUE,
    open = interactive(),
    overwrite = FALSE,
    quiet = TRUE) {

  for (flag in c("include_complete_output", "render", "open", "overwrite", "quiet")) {
    .mira_report_validate_scalar_flag_long(get(flag, inherits = FALSE), flag)
  }

  for (metadata_name in c("output_file", "title")) {
    value <- get(metadata_name, inherits = FALSE)
    if (!is.character(value) || length(value) != 1L || is.na(value) ||
        !nzchar(trimws(value))) {
      stop(sprintf("%s must be a non-empty string.", metadata_name),
           call. = FALSE)
    }
  }
  if (!is.null(subtitle) &&
      (!is.character(subtitle) || length(subtitle) != 1L || is.na(subtitle))) {
    stop("subtitle must be NULL or a single string.", call. = FALSE)
  }
  if (!is.null(author) &&
      (!is.character(author) || anyNA(author) || any(!nzchar(trimws(author))))) {
    stop("author must be NULL or a vector of non-empty strings.", call. = FALSE)
  }
  if (length(date) != 1L || is.na(date)) {
    stop("date must have length one and must not be missing.", call. = FALSE)
  }

  allowed_sections <- c(
    "all", "abstract", "methods", "data_quality", "overview",
    "descriptives", "frequencies", "missingness", "change", "arms",
    "correlations", "variability", "models", "robustness", "outliers",
    "trajectories", "figures", "diagnostics", "conclusions", "appendix"
  )
  if (!is.character(sections) || length(sections) == 0L || anyNA(sections)) {
    stop("sections must be 'all' or a vector of section names.", call. = FALSE)
  }
  sections <- unique(tolower(sections))
  invalid_sections <- setdiff(sections, allowed_sections)
  if (length(invalid_sections) > 0L) {
    stop(sprintf(
      "Unrecognized sections: %s. Allowed values: %s.",
      paste(invalid_sections, collapse = ", "),
      paste(allowed_sections, collapse = ", ")
    ), call. = FALSE)
  }

  allowed_formats <- c("html", "pdf", "docx")
  if (!is.character(format) || length(format) == 0L || anyNA(format)) {
    stop("format must contain html, pdf, docx, or all.", call. = FALSE)
  }
  formats <- unique(tolower(format))
  invalid_formats <- setdiff(formats, c(allowed_formats, "all"))
  if (length(invalid_formats) > 0L) {
    stop(sprintf("Unrecognized formats: %s.",
                 paste(invalid_formats, collapse = ", ")), call. = FALSE)
  }
  if ("all" %in% formats) formats <- allowed_formats

  if (!is.numeric(max_table_rows) || length(max_table_rows) != 1L ||
      is.na(max_table_rows) || max_table_rows <= 0) {
    stop("max_table_rows must be a positive number or Inf.", call. = FALSE)
  }
  if (!is.numeric(max_table_columns) || length(max_table_columns) != 1L ||
      is.na(max_table_columns) || !is.finite(max_table_columns) ||
      max_table_columns < 2 || max_table_columns > .Machine$integer.max) {
    stop("max_table_columns must be a finite integer >= 2.", call. = FALSE)
  }
  if (!is.numeric(digits) || length(digits) != 1L || is.na(digits) ||
      !is.finite(digits) || digits < 1 || digits > 10) {
    stop("digits must be an integer between 1 and 10.", call. = FALSE)
  }
  if (!is.numeric(max_depth) || length(max_depth) != 1L || is.na(max_depth) ||
      !is.finite(max_depth) || max_depth < 1 || max_depth > .Machine$integer.max) {
    stop("max_depth must be a finite integer >= 1.", call. = FALSE)
  }
  if (!is.list(extra_objects)) {
    stop("extra_objects must be a list, preferably named.", call. = FALSE)
  }
  if (length(extra_objects) > 0L) {
    extra_names <- names(extra_objects)
    if (is.null(extra_names)) extra_names <- rep("", length(extra_objects))
    missing_names <- is.na(extra_names) | !nzchar(extra_names)
    extra_names[missing_names] <- paste0("extra_", which(missing_names))
    names(extra_objects) <- make.unique(extra_names)
  }

  dots <- list(...)
  if (is.data.frame(x)) {
    if (!is.null(data)) {
      stop("Supply data only once, through x or data.", call. = FALSE)
    }
    data <- x
    x <- NULL
  }

  if (!is.null(data)) {
    if (!is.null(x)) {
      stop("Supply a MIRA object in x or a data frame in data, but not both.",
           call. = FALSE)
    }
    if (!is.data.frame(data)) stop("data must be a data frame.", call. = FALSE)
    mira_fun <- get0("mira_info_long", mode = "function", inherits = TRUE)
    if (is.null(mira_fun)) {
      stop("mira_info_long() is not available; load the file that defines it first.",
           call. = FALSE)
    }
    if ("data" %in% names(dots)) {
      stop("Do not repeat data in ...; use the data argument.", call. = FALSE)
    }
    if (!"verbose" %in% names(dots)) dots$verbose <- FALSE
    x <- do.call(mira_fun, c(list(data = data), dots))
  } else if (length(dots) > 0L) {
    stop("Arguments in ... are allowed only when data is supplied.", call. = FALSE)
  }

  if (is.null(x)) {
    stop("Supply x (the output of mira_info_long) or data.", call. = FALSE)
  }
  accepted <- inherits(x, c("mira_info_long", "mira_info_multi_long", "mira_detect_long"))
  if (!accepted) {
    stop(
      "x must inherit from mira_info_long, mira_info_multi_long, or mira_detect_long.",
      call. = FALSE
    )
  }

  output_file <- .mira_report_slug_long(output_file)
  if (is.null(output_dir)) {
    output_dir <- file.path(getwd(), "mira_analyses", "mira_report_flong")
  }
  if (!is.character(output_dir) || length(output_dir) != 1L ||
      is.na(output_dir) || !nzchar(trimws(output_dir))) {
    stop("output_dir must be a non-empty path.", call. = FALSE)
  }
  if (!dir.exists(output_dir)) {
    created <- dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    if (!created && !dir.exists(output_dir)) {
      stop(sprintf("Cannot create directory: %s", output_dir), call. = FALSE)
    }
  }
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)

  stem <- output_file
  qmd_name <- paste0(stem, ".qmd")
  payload_name <- paste0(stem, "-payload.rds")
  runtime_name <- paste0(stem, "-runtime.R")
  qmd_path <- file.path(output_dir, qmd_name)
  payload_path <- file.path(output_dir, payload_name)
  runtime_path <- file.path(output_dir, runtime_name)
  extensions <- c(html = ".html", pdf = ".pdf", docx = ".docx")
  output_paths <- stats::setNames(
    file.path(output_dir, paste0(stem, unname(extensions[formats]))), formats
  )

  report_config <- list(
    sections = sections,
    include_complete_output = include_complete_output,
    max_table_rows = max_table_rows,
    max_table_columns = as.integer(max_table_columns),
    digits = as.integer(digits),
    max_depth = as.integer(max_depth),
    created_at = base::format(Sys.time(), tz = "UTC", usetz = TRUE)
  )
  payload <- list(
    result = x,
    extra_objects = extra_objects,
    report = report_config
  )

  export_plan <- .mira_report_export_plan_long(x, report_config, output_dir)
  manifest <- .mira_report_manifest_long(export_plan)
  manifest_path <- file.path(output_dir, "results_manifest.csv")
  planned_csv <- vapply(export_plan, function(entry) entry$absolute, character(1L))
  protected_paths <- c(qmd_path, payload_path, runtime_path,
                       unname(output_paths), manifest_path, planned_csv)
  existing <- protected_paths[file.exists(protected_paths)]
  if (length(existing) > 0L && !overwrite) {
    stop(sprintf(
      "Output files already exist. Use overwrite=TRUE or a different directory: %s",
      paste(basename(existing), collapse = ", ")
    ), call. = FALSE)
  }

  previous_csv <- character()
  if (overwrite && file.exists(manifest_path)) {
    previous <- tryCatch(utils::read.csv(manifest_path, stringsAsFactors = FALSE),
                         error = function(e) NULL)
    if (!is.null(previous) && "file" %in% names(previous)) {
      safe <- grepl("^results/([A-Za-z0-9_-]+/)*[A-Za-z0-9_-]+\\.csv$",
                    previous$file)
      previous_csv <- file.path(output_dir, previous$file[safe])
    }
  }

  saveRDS(payload, payload_path, compress = "xz")
  runtime_names <- .mira_report_runtime_names_long()
  runtime_env <- environment(mira_report_freq_long)
  missing_runtime <- runtime_names[!vapply(runtime_names, exists, logical(1L),
                                           envir = runtime_env, inherits = TRUE)]
  if (length(missing_runtime) > 0L) {
    stop(sprintf("Missing internal helpers: %s.",
                 paste(missing_runtime, collapse = ", ")), call. = FALSE)
  }
  dump(runtime_names, file = runtime_path, envir = runtime_env)

  qmd_lines <- .mira_report_qmd_lines_long(
    title = title, subtitle = subtitle, author = author, date = date,
    formats = formats, runtime_file = runtime_name, payload_file = payload_name
  )
  writeLines(qmd_lines, qmd_path, useBytes = TRUE)

  for (entry in export_plan) {
    if (!dir.exists(dirname(entry$absolute))) {
      dir.create(dirname(entry$absolute), recursive = TRUE, showWarnings = FALSE)
    }
    tryCatch(
      utils::write.csv(entry$csv, entry$absolute, row.names = FALSE,
                       fileEncoding = "UTF-8", na = ""),
      error = function(e) stop(sprintf(
        "Could not export table at %s to %s: %s",
        entry$path, entry$file, conditionMessage(e)
      ), call. = FALSE)
    )
  }
  utils::write.csv(manifest, manifest_path, row.names = FALSE,
                   fileEncoding = "UTF-8", na = "")
  stale <- .mira_report_safe_stale_files_long(
    setdiff(previous_csv, planned_csv), output_dir
  )
  if (length(stale)) unlink(stale)

  rendered <- stats::setNames(rep(FALSE, length(formats)), formats)
  render_errors <- stats::setNames(rep(NA_character_, length(formats)), formats)
  if (render) {
    if (!requireNamespace("knitr", quietly = TRUE)) {
      stop(
        paste0("The R package 'knitr' is required. Source files were created in: ",
               output_dir),
        call. = FALSE
      )
    }
    if (!requireNamespace("quarto", quietly = TRUE)) {
      stop(
        paste0("The R package 'quarto' is required. Source files were created in: ",
               output_dir),
        call. = FALSE
      )
    }

    for (fmt in formats) {
      target_name <- basename(output_paths[[fmt]])
      error_message <- tryCatch({
        quarto::quarto_render(
          input = qmd_path,
          output_format = fmt,
          output_file = target_name,
          execute_dir = output_dir,
          quiet = quiet,
          as_job = FALSE
        )
        NULL
      }, error = function(e) conditionMessage(e))
      if (is.null(error_message) && !file.exists(output_paths[[fmt]])) {
        error_message <- sprintf(
          "Quarto returned without creating the expected %s file.", fmt
        )
      }
      rendered[[fmt]] <- is.null(error_message) && file.exists(output_paths[[fmt]])
      if (!is.null(error_message)) render_errors[[fmt]] <- error_message
    }

    failed <- names(rendered)[!rendered]
    if (length(failed) > 0L) {
      details <- paste(
        sprintf("%s: %s", failed, render_errors[failed]),
        collapse = " | "
      )
      warning(sprintf(
        "Rendering did not complete for %s. Source files were retained. %s",
        paste(failed, collapse = ", "), details
      ), call. = FALSE)
    }
  }

  result <- structure(
    list(
      call = match.call(),
      mira_result = x,
      formats = formats,
      rendered = rendered,
      render_errors = render_errors,
      files = output_paths[rendered],
      expected_files = output_paths,
      source = qmd_path,
      payload = payload_path,
      runtime = runtime_path,
      results_manifest = manifest_path,
      result_files = planned_csv,
      output_dir = output_dir,
      report_config = report_config
    ),
    class = c("mira_report_freq_long", "list")
  )

  if (open && any(rendered)) {
    first_file <- unname(output_paths[which(rendered)[1L]])
    try(utils::browseURL(first_file), silent = TRUE)
  }
  print(result)
  invisible(result)
}

#' Print a MIRA report descriptor
#'
#' Prints the bundle directory, generated Quarto source, successfully rendered
#' formats, and requested formats that were not rendered successfully. The
#' method does not inspect or rerender report contents.
#'
#' @param x An object returned by [mira_report_freq_long()].
#' @param ... Reserved for compatibility with the \code{print} generic;
#'   currently ignored.
#'
#' @return \code{x}, invisibly.
#' @seealso [mira_report_freq_long()]
#' @method print mira_report_freq_long
#' @export
print.mira_report_freq_long <- function(x, ...) {
  cat("MIRA Quarto report\n")
  cat(sprintf("Directory: %s\n", x$output_dir))
  cat(sprintf("Quarto source: %s\n", basename(x$source)))
  if (any(x$rendered)) {
    cat("Rendered:\n")
    for (name in names(x$rendered)[x$rendered]) {
      cat(sprintf("  - %s: %s\n", name, basename(x$expected_files[[name]])))
    }
  }
  if (any(!x$rendered)) {
    label <- if (all(!x$rendered)) "Not rendered" else "Not rendered successfully"
    cat(sprintf("%s: %s\n", label, paste(names(x$rendered)[!x$rendered], collapse = ", ")))
  }
  invisible(x)
}
