# ============================================================
# MIRA_REPORT_FREQ
# Dynamic, publication-oriented Quarto reporting for mira_info()
#
# Typical use:
#   source("mira_info(2).R")
#   source("mira_report_freq.R")
#
#   fit <- mira_info(data = dati, outcomes = c("score_a", "score_b"),
#                    analyses = "all", verbose = FALSE)
#   mira_report_freq(fit, format = c("html", "pdf"))
#
# Or run mira_info() and build the report in one call:
#   mira_report_freq(
#     data = dati,
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

.mira_report_or <- function(x, y) {
  if (is.null(x) || length(x) == 0L) y else x
}

.mira_report_is_empty <- function(x) {
  if (is.null(x)) return(TRUE)
  if (is.data.frame(x) || is.matrix(x) || is.table(x)) return(nrow(x) == 0L)
  if (is.atomic(x) || is.list(x)) return(length(x) == 0L)
  FALSE
}

.mira_report_heading <- function(text, level = 1L) {
  text <- gsub("[\r\n]+", " ", as.character(text)[1L])
  level <- max(1L, min(6L, as.integer(level)[1L]))
  cat("\n", strrep("#", level), " ", text, "\n\n", sep = "")
}

.mira_report_human_name <- function(x) {
  labels <- c(
    call = "Chiamata analitica",
    version = "Versione MIRA",
    settings = "Impostazioni",
    overview = "Quadro generale",
    config = "Configurazione risolta",
    data_overview = "Struttura del dataset",
    detected_variables = "Variabili rilevate",
    diagnostics = "Diagnostica e adattamenti",
    descriptives = "Statistiche descrittive",
    missing = "Completezza e dati non disponibili",
    by_time = "Disponibilità per visita",
    by_patient = "Disponibilità per soggetto",
    change = "Cambiamenti entro soggetto",
    arm_analysis = "Analisi per braccio",
    baseline_balance = "Bilanciamento al basale",
    missingness = "Missingness per braccio",
    time_omnibus = "Confronti omnibus tra bracci per visita",
    time_pairwise = "Confronti a coppie tra bracci per visita",
    change_descriptives = "Cambiamento per braccio",
    change_omnibus = "Confronti omnibus del cambiamento",
    change_pairwise = "Confronti a coppie del cambiamento",
    correlations = "Correlazioni tra visite",
    pearson = "Correlazione di Pearson",
    spearman = "Correlazione di Spearman",
    pairwise_n = "Numerosità disponibili a coppie",
    variability = "Variabilità e ICC",
    trajectories = "Traiettorie individuali",
    model = "Modello misto principale",
    fitted_model = "Oggetto del modello stimato",
    summary = "Sintesi del modello",
    anova = "Tabella ANOVA del modello",
    global_time_test = "Test globale del tempo",
    global_arm_test = "Test globale del braccio",
    arm_time_interaction_test = "Test globale interazione braccio per tempo",
    fixed_parameters = "Parametri a effetti fissi",
    advanced_tests = "Test longitudinali avanzati",
    advanced_models = "Modelli longitudinali avanzati",
    robustness = "Inferenza robusta",
    effect_sizes = "Dimensioni dell'effetto",
    multiplicity = "Famiglie di molteplicità",
    sensitivity = "Analisi di sensibilità",
    outliers = "Valori anomali diagnostici",
    plots = "Figure",
    plot_error = "Diagnostica delle figure",
    long_data = "Dataset longitudinale analitico",
    outcomes = "Outcome",
    warnings = "Avvertenze",
    adaptation = "Adattamenti automatici",
    failed_outcomes = "Outcome non completati",
    emmeans = "Medie marginali stimate e contrasti",
    rm_anova = "ANOVA per misure ripetute",
    friedman = "Test di Friedman",
    random_slope = "Modello misto con pendenza casuale",
    nlme = "Modelli NLME",
    gee = "Modelli GEE",
    model_comparison = "Confronto descrittivo dei modelli",
    club_sandwich = "Inferenza cluster-robust CR2",
    compound_symmetry = "Correlazione compound symmetry",
    ar1 = "Correlazione AR(1)",
    independence = "Correlazione di lavoro indipendente",
    exchangeable = "Correlazione di lavoro scambiabile",
    data = "Dati analitici del modulo",
    formula = "Formula del modello",
    error = "Errore del modulo",
    note = "Nota",
    boxplot = "Distribuzione per visita",
    mean_ci = "Profilo medio con intervalli di confidenza",
    change_from_baseline = "Cambiamento rispetto al basale",
    change_ci = "Cambiamento medio con intervalli di confidenza",
    correlation_heatmap = "Matrice di correlazione",
    arm_mean_ci = "Profilo medio per braccio con intervalli di confidenza",
    arm_boxplot = "Distribuzione per braccio e visita",
    arm_change = "Cambiamento per braccio",
    arm_change_ci = "Cambiamento medio per braccio con intervalli di confidenza",
    arm_difference_ci = "Differenze tra bracci con intervalli di confidenza",
    arm_missingness = "Dati non disponibili per braccio",
    response = "Risposta clinica individuale"
  )
  key <- as.character(x)[1L]
  if (key %in% names(labels)) return(unname(labels[[key]]))
  text <- gsub("[._]+", " ", key)
  text <- trimws(text)
  if (!nzchar(text)) return("Elemento senza nome")
  paste0(toupper(substr(text, 1L, 1L)), substr(text, 2L, nchar(text)))
}

.mira_report_format_p <- function(x, digits = 3L) {
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

.mira_report_format_number <- function(x, digits = 3L) {
  if (length(x) == 0L || is.na(x[1L])) return("non stimabile")
  value <- suppressWarnings(as.numeric(x[1L]))
  if (!is.finite(value)) return(as.character(value))
  formatC(value, format = "fg", digits = digits, flag = "#")
}

.mira_report_inline <- function(x) {
  if (is.null(x) || length(x) == 0L) return("non disponibile")
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

.mira_report_prepare_table <- function(x, digits = 3L) {
  if (is.table(x)) x <- as.data.frame(x, stringsAsFactors = FALSE)

  if (is.matrix(x)) {
    row_labels <- rownames(x)
    x <- as.data.frame(x, stringsAsFactors = FALSE, check.names = FALSE)
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
        parametro = nms,
        valore = vapply(x, .mira_report_inline, character(1L)),
        stringsAsFactors = FALSE
      )
    } else if (is.atomic(x)) {
      nms <- names(x)
      if (is.null(nms)) nms <- seq_along(x)
      x <- data.frame(
        elemento = as.character(nms),
        valore = as.character(x),
        stringsAsFactors = FALSE
      )
    } else {
      stop("L'oggetto non è convertibile in tabella.", call. = FALSE)
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

  p_columns <- grepl(
    "(^p$|(^|[._])p([._]|$)|p.value|p_value|pvalue|^pr\\()",
    names(x), ignore.case = TRUE, perl = TRUE
  )
  for (j in which(p_columns)) {
    if (is.numeric(x[[j]]) || is.integer(x[[j]])) {
      x[[j]] <- .mira_report_format_p(x[[j]], digits = digits)
    }
  }
  x
}

.mira_report_numeric_column <- function(x, name = "") {
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

.mira_report_column_widths <- function(x) {
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

.mira_report_partition_columns <- function(x, max_columns, width_budget = Inf) {
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

.mira_report_latex_alignment <- function(x) {
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

.mira_report_latex_breaks <- function(x) {
  for (token in c("\\_", "\\$")) {
    pieces <- strsplit(x, token, fixed = TRUE)[[1L]]
    if (length(pieces) > 1L) {
      x <- paste(pieces, collapse = paste0(token, "\\allowbreak{}"))
    }
  }
  x
}

.mira_report_repeat_longtable_header <- function(x) {
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

.mira_report_table <- function(x, caption = NULL, cfg) {
  if (.mira_report_is_empty(x)) {
    cat("*Nessuna riga disponibile o stimabile.*\n\n")
    return(invisible(NULL))
  }

  table_data <- tryCatch(
    .mira_report_prepare_table(x, digits = cfg$digits),
    error = function(e) NULL
  )
  if (is.null(table_data)) {
    cat("*Oggetto non rappresentabile come tabella; viene riportato in forma testuale.*\n\n")
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
      block_caption <- sprintf("%s (blocco %d di %d)", caption, i, length(blocks))
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
          "La tabella non è stata convertita dal formatter di knitr (`",
          render_error,
          "`). Il contenuto viene mantenuto in forma testuale per non interrompere il report."
        ),
        type = "warning", title = "Fallback tabellare"
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
      "*Sono visualizzate %d di %d righe per il limite esplicito `max_table_rows`.*\n\n",
      nrow(table_data), original_n
    ))
  }
  invisible(table_data)
}

.mira_report_prepare_namespace <- function(x) {
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

.mira_report_text_block <- function(x) {
  .mira_report_prepare_namespace(x)
  text <- tryCatch(
    {
      if (inherits(x, c("lm", "glm", "merMod", "lmerMod", "glmerMod"))) {
        capture.output(print(summary(x)))
      } else {
        capture.output(print(x))
      }
    },
    error = function(e) paste("Impossibile stampare l'oggetto:", conditionMessage(e))
  )
  text <- gsub("```", "'''", enc2utf8(text), fixed = TRUE)
  cat("```text\n", paste(text, collapse = "\n"), "\n```\n\n", sep = "")
}

.mira_report_callout <- function(text, type = "note", title = "Nota") {
  allowed <- c("note", "tip", "warning", "important", "caution")
  if (!type %in% allowed) type <- "note"
  cat(sprintf("::: {.callout-%s appearance=\"simple\"}\n", type))
  cat("**", title, ".** ", text, "\n", sep = "")
  cat(":::\n\n")
}

.mira_report_is_plot <- function(x) {
  inherits(x, c("ggplot", "ggplot2::ggplot", "recordedplot", "grob", "gTree"))
}

.mira_report_plot <- function(x, name = "Figura") {
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
                         title = paste("Figura non renderizzata:", name))
    FALSE
  })
  if (printed) {
    caption_map <- c(
      boxplot = "Distribuzione dei valori osservati per visita, con osservazioni individuali e intervalli di confidenza della media.",
      trajectories = "Traiettorie individuali longitudinali.",
      trajectory = "Traiettorie individuali longitudinali.",
      spaghetti = "Traiettorie individuali longitudinali.",
      mean = "Profilo medio longitudinale.",
      mean_ci = "Profilo medio longitudinale con intervalli di confidenza.",
      change = "Distribuzione dei cambiamenti entro soggetto.",
      change_from_baseline = "Cambiamento rispetto al basale lungo il follow-up.",
      change_ci = "Cambiamento medio longitudinale con intervalli di confidenza.",
      missingness = "Percentuale di osservazioni non disponibili per visita.",
      correlation = "Struttura di correlazione tra visite.",
      arm = "Andamento longitudinale stratificato per braccio.",
      arm_mean_ci = "Profilo medio longitudinale per braccio con intervalli di confidenza.",
      arm_boxplot = "Distribuzione dei valori osservati per braccio e visita.",
      arm_change = "Cambiamento longitudinale stratificato per braccio.",
      arm_change_ci = "Cambiamento medio per braccio con intervalli di confidenza.",
      arm_difference_ci = "Differenze stimate tra bracci con intervalli di confidenza.",
      arm_missingness = "Percentuale di osservazioni non disponibili per braccio e visita.",
      response = "Frequenza delle direzioni di risposta individuale."
    )
    matching_keys <- names(caption_map)[vapply(names(caption_map), function(k) {
      grepl(k, name, ignore.case = TRUE, fixed = TRUE)
    }, logical(1L))]
    if (length(matching_keys) > 1L) {
      matching_keys <- matching_keys[order(nchar(matching_keys), decreasing = TRUE)]
    }
    caption <- if (length(matching_keys) > 0L) {
      caption_map[[matching_keys[[1L]]]]
    } else {
      paste0("Output grafico `", name, "` prodotto da MIRA.")
    }
    cat("**Figura.** ", caption, "\n\n", sep = "")
  }
  invisible(printed)
}

.mira_report_is_plain_list <- function(x) {
  if (!is.list(x)) return(FALSE)
  cls <- class(x)
  is.null(cls) || identical(cls, "list") ||
    inherits(x, c("mira_info", "mira_info_multi", "mira_info_error",
                  "mira_detect", "summary.mira_info", "summary.mira_info_multi"))
}

.mira_report_render_any <- function(x, name, level, cfg, depth = 0L) {
  if (depth > cfg$max_depth) {
    .mira_report_callout(
      sprintf("Profondità massima di sicurezza raggiunta nel percorso `%s`.", name),
      type = "warning", title = "Output annidato"
    )
    return(invisible(NULL))
  }

  if (.mira_report_is_empty(x)) {
    cat("*Output non disponibile, non richiesto o non stimabile.*\n\n")
    return(invisible(NULL))
  }

  if (.mira_report_is_plot(x)) {
    .mira_report_plot(x, name)
    return(invisible(NULL))
  }

  if (is.data.frame(x) || is.matrix(x) || is.table(x)) {
    .mira_report_table(x, caption = .mira_report_human_name(name), cfg = cfg)
    note <- attr(x, "note", exact = TRUE)
    interpretation <- attr(x, "interpretation", exact = TRUE)
    if (!is.null(note)) cat("*Nota:* ", .mira_report_inline(note), "\n\n", sep = "")
    if (!is.null(interpretation)) {
      cat("*Interpretazione:* ", .mira_report_inline(interpretation), "\n\n", sep = "")
    }
    return(invisible(NULL))
  }

  if (inherits(x, "formula") || is.call(x) || is.expression(x) || is.name(x)) {
    .mira_report_text_block(x)
    return(invisible(NULL))
  }

  if (is.atomic(x)) {
    if (length(x) == 1L) {
      cat("**Valore:** `", .mira_report_inline(x), "`\n\n", sep = "")
    } else {
      .mira_report_table(x, caption = .mira_report_human_name(name), cfg = cfg)
    }
    return(invisible(NULL))
  }

  if (.mira_report_is_plain_list(x)) {
    if (inherits(x, "mira_info_error")) {
      .mira_report_callout(
        .mira_report_or(x$error, "Errore analitico non specificato."),
        type = "warning", title = "Analisi non completata"
      )
      return(invisible(NULL))
    }

    nms <- names(x)
    if (is.null(nms)) nms <- paste0("item_", seq_along(x))
    for (i in seq_along(x)) {
      item_name <- nms[[i]]
      if (level <= 6L) {
        .mira_report_heading(.mira_report_human_name(item_name), level)
      } else {
        cat("**", .mira_report_human_name(item_name), "**\n\n", sep = "")
      }
      .mira_report_render_any(
        x[[i]], name = item_name, level = min(6L, level + 1L),
        cfg = cfg, depth = depth + 1L
      )
    }
    return(invisible(NULL))
  }

  .mira_report_text_block(x)
  invisible(NULL)
}

.mira_report_extract_outcomes <- function(result) {
  if (inherits(result, "mira_detect")) return(list())
  if (!is.null(result$outcomes) && is.list(result$outcomes) &&
      length(result$outcomes) > 0L) {
    outcomes <- result$outcomes
  } else {
    outcomes <- list(result)
  }
  nms <- names(outcomes)
  if (is.null(nms)) nms <- rep("", length(outcomes))
  for (i in seq_along(outcomes)) {
    if (!nzchar(nms[[i]])) {
      nms[[i]] <- .mira_report_or(outcomes[[i]]$outcome, paste0("outcome_", i))
    }
  }
  names(outcomes) <- make.unique(nms)
  outcomes
}

.mira_report_first_success <- function(outcomes) {
  if (length(outcomes) == 0L) return(NULL)
  ok <- !vapply(outcomes, inherits, logical(1L), what = "mira_info_error")
  if (!any(ok)) return(NULL)
  outcomes[[which(ok)[1L]]]
}

.mira_report_scalar_table <- function(x) {
  if (is.null(x) || !is.list(x)) return(data.frame())
  nms <- names(x)
  if (is.null(nms)) nms <- paste0("item_", seq_along(x))
  keep <- vapply(x, function(value) {
    is.null(value) || is.atomic(value) || inherits(value, c("Date", "POSIXct", "POSIXlt"))
  }, logical(1L))
  data.frame(
    parametro = vapply(nms[keep], .mira_report_human_name, character(1L)),
    valore = vapply(x[keep], .mira_report_inline, character(1L)),
    stringsAsFactors = FALSE
  )
}

.mira_report_frequency <- function(x, label = "categoria") {
  x <- as.character(x)
  x[is.na(x) | !nzchar(x)] <- "Non disponibile"
  counts <- sort(table(x, useNA = "no"), decreasing = TRUE)
  total <- sum(counts)
  out <- data.frame(
    categoria = names(counts),
    n = as.integer(counts),
    percentuale = if (total > 0L) as.numeric(counts) / total * 100 else numeric(length(counts)),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  stats::setNames(out, c(label, "n", "percentuale"))
}

.mira_report_object_index <- function(x, root = "result", max_depth = 50L) {
  rows <- list()
  counter <- 0L

  walk <- function(object, path, depth) {
    counter <<- counter + 1L
    dimensions <- tryCatch({
      if (!is.null(dim(object))) paste(dim(object), collapse = " x ")
      else as.character(length(object))
    }, error = function(e) "non disponibile")
    status <- NA_character_
    object_status <- if (is.list(object)) object[["status", exact = TRUE]] else NULL
    if (!is.null(object_status) && length(object_status) == 1L) {
      status <- as.character(object_status)
    } else if (.mira_report_is_empty(object)) {
      status <- "vuoto/non disponibile"
    }
    size <- tryCatch(
      format(utils::object.size(object), units = "auto", standard = "SI"),
      error = function(e) NA_character_
    )
    rows[[counter]] <<- data.frame(
      percorso = path,
      classe = paste(class(object), collapse = "/"),
      dimensione = dimensions,
      memoria = size,
      stato = status,
      stringsAsFactors = FALSE
    )

    if (depth >= max_depth || !.mira_report_is_plain_list(object) ||
        .mira_report_is_empty(object)) return(invisible(NULL))
    nms <- names(object)
    if (is.null(nms)) nms <- paste0("[[", seq_along(object), "]]" )
    for (i in seq_along(object)) {
      separator <- if (grepl("^\\[\\[", nms[[i]])) "" else "$"
      walk(object[[i]], paste0(path, separator, nms[[i]]), depth + 1L)
    }
    invisible(NULL)
  }

  walk(x, root, 0L)
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

.mira_report_has_section <- function(cfg, section) {
  "all" %in% cfg$sections || section %in% cfg$sections
}

.mira_report_longitudinal_map <- function(config) {
  if (is.null(config$time_vars)) return(data.frame())
  vars <- config$time_vars
  labs <- config$time_labels
  if (!is.list(vars)) {
    outcome_name <- .mira_report_or(config$outcomes, "outcome")
    vars <- stats::setNames(list(vars), as.character(outcome_name)[1L])
  }
  if (!is.list(labs)) labs <- list(labs)

  rows <- list()
  counter <- 0L
  outcome_names <- names(vars)
  if (is.null(outcome_names)) outcome_names <- paste0("outcome_", seq_along(vars))
  for (i in seq_along(vars)) {
    outcome <- outcome_names[[i]]
    variables <- as.character(vars[[i]])
    labels <- labs[[outcome]]
    if (is.null(labels) && length(labs) >= i) labels <- labs[[i]]
    if (is.null(labels)) labels <- variables
    labels <- as.character(labels)
    if (length(labels) != length(variables)) labels <- rep_len(labels, length(variables))
    for (j in seq_along(variables)) {
      counter <- counter + 1L
      rows[[counter]] <- data.frame(
        outcome = outcome,
        ordine = j,
        variabile = variables[[j]],
        etichetta_visita = labels[[j]],
        stringsAsFactors = FALSE
      )
    }
  }
  if (length(rows) == 0L) return(data.frame())
  do.call(rbind, rows)
}

.mira_report_abstract_rows <- function(outcomes) {
  rows <- list()
  counter <- 0L
  for (name in names(outcomes)) {
    outcome <- outcomes[[name]]
    counter <- counter + 1L
    if (inherits(outcome, "mira_info_error")) {
      rows[[counter]] <- data.frame(
        outcome = name, stato = "non completato", soggetti = NA_integer_,
        visite = NA_integer_, profili_completi_pct = NA_real_,
        media_basale = NA_real_, media_finale = NA_real_,
        cambiamento_basale_finale = NA_real_, p_aggiustato = NA_real_,
        stringsAsFactors = FALSE
      )
      next
    }

    desc <- outcome$descriptives
    change <- outcome$change
    baseline_mean <- final_mean <- change_mean <- p_adjusted <- NA_real_
    if (is.data.frame(desc) && nrow(desc) > 0L && "mean" %in% names(desc)) {
      baseline_mean <- desc$mean[[1L]]
      final_mean <- desc$mean[[nrow(desc)]]
    }
    if (is.data.frame(change) && nrow(change) > 0L) {
      index <- seq_len(nrow(change))
      if (all(c("from", "to") %in% names(change)) &&
          length(outcome$time_vars) >= 2L) {
        target <- which(
          change$from == outcome$time_vars[[1L]] &
            change$to == outcome$time_vars[[length(outcome$time_vars)]]
        )
        if (length(target) > 0L) index <- target[[1L]]
      }
      index <- index[[1L]]
      if ("mean_change" %in% names(change)) change_mean <- change$mean_change[[index]]
      p_candidates <- c("paired_t_p_adj", "wilcoxon_p_adj", "paired_t_p", "wilcoxon_p")
      p_name <- p_candidates[p_candidates %in% names(change)]
      if (length(p_name) > 0L) p_adjusted <- change[[p_name[[1L]]]][[index]]
    }

    overview <- outcome$overview
    rows[[counter]] <- data.frame(
      outcome = .mira_report_or(outcome$outcome_display, name),
      stato = "completato",
      soggetti = .mira_report_or(overview$n_patients, NA_integer_),
      visite = .mira_report_or(overview$n_timepoints, length(outcome$time_vars)),
      profili_completi_pct = .mira_report_or(overview$complete_profiles_pct, NA_real_),
      media_basale = baseline_mean,
      media_finale = final_mean,
      cambiamento_basale_finale = change_mean,
      p_aggiustato = p_adjusted,
      stringsAsFactors = FALSE
    )
  }
  if (length(rows) == 0L) return(data.frame())
  do.call(rbind, rows)
}

.mira_report_key_narrative <- function(outcome, name, cfg) {
  if (inherits(outcome, "mira_info_error")) {
    .mira_report_callout(
      .mira_report_or(outcome$error, "Analisi non completata."),
      type = "warning", title = paste("Outcome", name)
    )
    return(invisible(NULL))
  }

  desc <- outcome$descriptives
  if (!is.data.frame(desc) || nrow(desc) == 0L) return(invisible(NULL))
  display <- .mira_report_or(outcome$outcome_display, name)
  first_label <- if ("label" %in% names(desc)) desc$label[[1L]] else "prima visita"
  last_label <- if ("label" %in% names(desc)) desc$label[[nrow(desc)]] else "ultima visita"
  first_mean <- if ("mean" %in% names(desc)) desc$mean[[1L]] else NA_real_
  last_mean <- if ("mean" %in% names(desc)) desc$mean[[nrow(desc)]] else NA_real_
  cat(sprintf(
    "Per **%s**, la media osservata e **%s** alla visita %s e **%s** alla visita finale %s. ",
    display,
    .mira_report_format_number(first_mean, cfg$digits), first_label,
    .mira_report_format_number(last_mean, cfg$digits), last_label
  ))

  change <- outcome$change
  if (is.data.frame(change) && nrow(change) > 0L &&
      "mean_change" %in% names(change)) {
    index <- nrow(change)
    if (all(c("from", "to") %in% names(change)) && length(outcome$time_vars) >= 2L) {
      target <- which(
        change$from == outcome$time_vars[[1L]] &
          change$to == outcome$time_vars[[length(outcome$time_vars)]]
      )
      if (length(target) > 0L) index <- target[[1L]]
    }
    p_candidates <- c("paired_t_p_adj", "wilcoxon_p_adj", "paired_t_p", "wilcoxon_p")
    p_name <- p_candidates[p_candidates %in% names(change)]
    p_text <- if (length(p_name) > 0L) {
      p_value <- change[[p_name[[1L]]]][[index]]
      p_label <- if (is.na(p_value)) "non stimabile" else {
        .mira_report_format_p(p_value, cfg$digits)
      }
      paste0("; `", p_name[[1L]], "` = ", p_label)
    } else ""
    cat(sprintf(
      "Il cambiamento medio appaiato basale-finale e **%s**%s. ",
      .mira_report_format_number(change$mean_change[[index]], cfg$digits), p_text
    ))
  }
  cat("Le stime sono descrittive salvo dove il test e indicato esplicitamente; l'interpretazione clinica dipende dalla scala e dal protocollo.\n\n")
  invisible(NULL)
}

.mira_report_render_abstract <- function(result, outcomes, cfg) {
  .mira_report_heading("Sintesi", 1L)
  if (inherits(result, "mira_detect")) {
    cat("L'esecuzione ricevuta e di tipo `mira_detect`: il documento descrive la configurazione individuata, senza risultati inferenziali.\n\n")
    return(invisible(NULL))
  }

  completed <- sum(!vapply(outcomes, inherits, logical(1L), what = "mira_info_error"))
  failed <- length(outcomes) - completed
  cat(sprintf(
    paste0("Il report sintetizza in modo dinamico **%d outcome**, dei quali **%d** completati",
           "%s. Le sezioni e le tabelle sono generate esclusivamente dagli output presenti ",
           "nell'oggetto `mira_info`; i moduli non richiesti o non stimabili sono dichiarati come tali.\n\n"),
    length(outcomes), completed,
    if (failed > 0L) paste0(" e **", failed, "** non completati") else ""
  ))
  .mira_report_table(
    .mira_report_abstract_rows(outcomes),
    caption = "Sintesi descrittiva degli outcome", cfg = cfg
  )
  for (name in names(outcomes)) {
    .mira_report_key_narrative(outcomes[[name]], name, cfg)
  }
  .mira_report_callout(
    paste0(
      "Le frasi automatiche sono un ausilio alla lettura, non sostituiscono il piano statistico, ",
      "la verifica delle assunzioni, la rilevanza clinica o la revisione di un analista. ",
      "Un valore p non misura la dimensione ne l'importanza dell'effetto."
    ),
    type = "important", title = "Interpretazione responsabile"
  )
  invisible(NULL)
}

.mira_report_render_methods <- function(result, outcomes, cfg) {
  .mira_report_heading("Metodi e configurazione analitica", 1L)
  config <- result$config
  first <- .mira_report_first_success(outcomes)
  settings <- if (!is.null(first)) first$settings else NULL

  cat(
    paste0(
      "MIRA analizza dati longitudinali in formato wide, costruisce internamente la forma long ",
      "e conserva nell'output la provenienza delle scelte automatiche e manuali. Questo documento ",
      "non ricalcola i risultati quando riceve un oggetto `mira_info`: li presenta e li indicizza ",
      "in modo riproducibile. Se viene fornito `data`, `mira_report_freq()` esegue prima `mira_info()` ",
      "con gli argomenti ricevuti in `...`.\n\n"
    )
  )

  method_rows <- data.frame(
    elemento = c(
      "Versione MIRA", "ID", "Outcome", "Braccio", "Braccio di riferimento",
      "Covariate", "Alpha", "Livello di confidenza", "Correzione primaria",
      "Direzione del miglioramento", "Soglia di stabilita",
      "Gestione Inf/-Inf"
    ),
    valore = c(
      .mira_report_or(result$version, "non disponibile"),
      .mira_report_or(config$id, "non disponibile"),
      .mira_report_inline(config$outcomes),
      .mira_report_or(config$arm, "nessuno"),
      .mira_report_or(config$reference_arm, "nessuno"),
      .mira_report_inline(config$covariates),
      .mira_report_inline(settings$alpha),
      .mira_report_inline(settings$confidence_level),
      .mira_report_or(settings$p_adjust_method, "non disponibile"),
      .mira_report_inline(config$improvement_direction),
      .mira_report_inline(config$stable_threshold),
      .mira_report_or(settings$non_finite_handling, "documentata nell'output")
    ),
    stringsAsFactors = FALSE
  )
  .mira_report_table(method_rows, caption = "Parametri principali dell'analisi", cfg = cfg)
  .mira_report_callout(
    paste0(
      "Con `include_complete_output = TRUE` il documento può contenere identificativi ",
      "individuali e dati analitici a livello di soggetto. Applicare le regole di ",
      "riservatezza, minimizzazione e condivisione previste dal proprio contesto."
    ),
    type = "caution", title = "Riservatezza"
  )

  longitudinal_map <- .mira_report_longitudinal_map(config)
  if (nrow(longitudinal_map) > 0L) {
    .mira_report_heading("Mappa outcome-visite", 2L)
    .mira_report_table(
      longitudinal_map,
      caption = "Variabili longitudinali effettivamente utilizzate", cfg = cfg
    )
  }

  if (!is.null(config$analyses)) {
    .mira_report_heading("Moduli richiesti e attivati", 2L)
    .mira_report_table(
      .mira_report_scalar_table(config$analyses),
      caption = "Configurazione dei moduli analitici", cfg = cfg
    )
  }

  if (!is.null(result$call)) {
    .mira_report_heading("Chiamata riproducibile", 2L)
    .mira_report_text_block(result$call)
  }

  .mira_report_heading("Criteri di lettura", 2L)
  glossary <- data.frame(
    indicatore = c(
      "CI", "p_raw", "p_adj", "Cohen's dz", "Hedges g", "ICC",
      "Kendall's W", "CR2/HTZ", "AIC/BIC", "QIC/CIC"
    ),
    interpretazione = c(
      "Intervallo di confidenza al livello impostato nell'analisi.",
      "Valore p non corretto; va letto nella famiglia di test pertinente.",
      "Valore p corretto con il metodo configurato, Holm per impostazione predefinita.",
      "Cambiamento medio standardizzato entro soggetto.",
      "Differenza standardizzata tra gruppi con correzione per piccoli campioni.",
      "Quota della variabilità attribuibile alle differenze tra soggetti nel modello indicato.",
      "Dimensione dell'effetto basata sui ranghi per il test di Friedman.",
      "Inferenza cluster-robust con correzione per piccoli campioni.",
      "Criteri descrittivi per modelli likelihood-based; valori inferiori favoriscono il compromesso adattamento-complessita.",
      "Criteri analoghi per modelli GEE; non sono direttamente intercambiabili con AIC/BIC."
    ),
    stringsAsFactors = FALSE
  )
  .mira_report_table(glossary, caption = "Glossario statistico essenziale", cfg = cfg)
  invisible(NULL)
}

.mira_report_render_data_quality <- function(result, cfg) {
  .mira_report_heading("Qualità, struttura e selezione dei dati", 1L)

  overview <- result$data_overview
  if (!is.null(overview)) {
    scalar <- .mira_report_scalar_table(overview)
    if (nrow(scalar) > 0L) {
      .mira_report_table(scalar, caption = "Dimensioni del dataset", cfg = cfg)
    }
    if (is.data.frame(overview$column_profile)) {
      .mira_report_heading("Profilo delle colonne", 2L)
      .mira_report_table(
        overview$column_profile,
        caption = "Tipo, cardinalita e missingness delle variabili", cfg = cfg
      )
    }
  }

  if (!is.null(result$config$auto_detected)) {
    .mira_report_heading("Scelte automatiche", 2L)
    .mira_report_table(
      .mira_report_scalar_table(result$config$auto_detected),
      caption = "Elementi determinati automaticamente", cfg = cfg
    )
  }
  if (!is.null(result$config$specified_manually)) {
    .mira_report_heading("Scelte specificate dall'utente", 2L)
    .mira_report_table(
      .mira_report_scalar_table(result$config$specified_manually),
      caption = "Elementi specificati manualmente", cfg = cfg
    )
  }

  detected <- result$detected_variables
  if (!is.null(detected)) {
    .mira_report_heading("Dettaglio del rilevamento", 2L)
    .mira_report_render_any(detected, "detected_variables", 3L, cfg)
  }

  warnings <- result$diagnostics$warnings
  if (!is.null(warnings) && length(warnings) > 0L) {
    warning_text <- paste0("- ", as.character(warnings), collapse = "\n")
    .mira_report_callout(warning_text, type = "warning", title = "Avvertenze MIRA")
  } else {
    .mira_report_callout(
      "L'oggetto non contiene avvertenze globali registrate.",
      type = "tip", title = "Controlli automatici"
    )
  }

  if (!is.null(result$diagnostics$adaptation)) {
    .mira_report_heading("Adattamenti per disponibilità dei dati", 2L)
    .mira_report_render_any(result$diagnostics$adaptation, "adaptation", 3L, cfg)
  }
  invisible(NULL)
}

.mira_report_render_frequencies <- function(outcome, cfg) {
  rendered <- FALSE
  if (is.data.frame(outcome$trajectories) && nrow(outcome$trajectories) > 0L) {
    if ("direction" %in% names(outcome$trajectories)) {
      .mira_report_heading("Direzione osservata", 3L)
      .mira_report_table(
        .mira_report_frequency(outcome$trajectories$direction, "direzione"),
        caption = "Frequenze del cambiamento basale-finale", cfg = cfg
      )
      rendered <- TRUE
    }
    if ("clinical_direction" %in% names(outcome$trajectories) &&
        any(!is.na(outcome$trajectories$clinical_direction))) {
      .mira_report_heading("Direzione clinica", 3L)
      .mira_report_table(
        .mira_report_frequency(outcome$trajectories$clinical_direction,
                               "classificazione_clinica"),
        caption = "Frequenze secondo la direzione clinica configurata", cfg = cfg
      )
      rendered <- TRUE
    }
  }

  missing_patient <- outcome$missing$by_patient
  if (is.data.frame(missing_patient) && nrow(missing_patient) > 0L &&
      "unavailable_n" %in% names(missing_patient)) {
    .mira_report_heading("Numero di visite non disponibili per soggetto", 3L)
    .mira_report_table(
      .mira_report_frequency(missing_patient$unavailable_n, "visite_non_disponibili"),
      caption = "Distribuzione della completezza individuale", cfg = cfg
    )
    rendered <- TRUE
  }

  counts <- outcome$arm_analysis$counts
  if (!is.null(counts) && length(counts) > 0L) {
    count_table <- data.frame(
      braccio = names(counts), n = as.integer(counts),
      percentuale = as.numeric(counts) / sum(counts) * 100,
      stringsAsFactors = FALSE
    )
    .mira_report_heading("Composizione dei bracci", 3L)
    .mira_report_table(count_table, caption = "Frequenze per braccio", cfg = cfg)
    rendered <- TRUE
  }
  if (!rendered) cat("*Nessuna tabella di frequenza derivabile dagli output disponibili.*\n\n")
  invisible(rendered)
}

.mira_report_render_figures <- function(outcome, cfg) {
  plots <- outcome$plots
  if (is.null(plots) || length(plots) == 0L) {
    reason <- .mira_report_or(
      outcome$plot_error,
      "Le figure non sono state richieste oppure non erano stimabili per questo outcome."
    )
    .mira_report_callout(reason, type = "note", title = "Figure non disponibili")
    return(invisible(NULL))
  }
  plot_names <- names(plots)
  if (is.null(plot_names)) plot_names <- paste0("figura_", seq_along(plots))
  for (i in seq_along(plots)) {
    .mira_report_heading(.mira_report_human_name(plot_names[[i]]), 3L)
    .mira_report_plot(plots[[i]], plot_names[[i]])
  }
  if (!is.null(outcome$plot_error)) {
    .mira_report_callout(outcome$plot_error, type = "warning",
                         title = "Diagnostica delle figure")
  }
  invisible(NULL)
}

.mira_report_render_models <- function(outcome, cfg) {
  model <- outcome$model
  if (is.null(model)) {
    cat("*Il modulo dei modelli non è presente nell'output.*\n\n")
  } else {
    .mira_report_heading("Modello principale", 3L)
    model_order <- c(
      "fixed_parameters", "anova", "global_time_test", "global_arm_test",
      "arm_time_interaction_test", "singular", "converged", "warnings", "error",
      "covariates_requested", "covariates_used", "covariates_skipped", "summary",
      "fitted_model"
    )
    remaining <- setdiff(names(model), model_order)
    ordered_names <- c(model_order[model_order %in% names(model)], remaining)
    for (name in ordered_names) {
      .mira_report_heading(.mira_report_human_name(name), 4L)
      .mira_report_render_any(model[[name]], name, 5L, cfg)
    }
  }

  if (!is.null(outcome$advanced_tests)) {
    .mira_report_heading("Test avanzati", 3L)
    .mira_report_render_any(outcome$advanced_tests, "advanced_tests", 4L, cfg)
  }
  if (!is.null(outcome$advanced_models)) {
    .mira_report_heading("Modelli avanzati", 3L)
    .mira_report_render_any(outcome$advanced_models, "advanced_models", 4L, cfg)
  }
  invisible(NULL)
}

.mira_report_render_outcome <- function(outcome, outcome_name, index, cfg) {
  cat("\n\\newpage\n\n")
  display <- if (inherits(outcome, "mira_info_error")) outcome_name else {
    .mira_report_or(outcome$outcome_display, outcome_name)
  }
  .mira_report_heading(sprintf("Outcome %d: %s", index, display), 1L)

  if (inherits(outcome, "mira_info_error")) {
    .mira_report_callout(
      .mira_report_or(outcome$error, "Analisi non completata."),
      type = "warning", title = "Outcome non analizzato"
    )
    return(invisible(NULL))
  }

  .mira_report_key_narrative(outcome, outcome_name, cfg)

  if (.mira_report_has_section(cfg, "overview")) {
    .mira_report_heading("Popolazione analitica", 2L)
    .mira_report_table(
      .mira_report_scalar_table(outcome$overview),
      caption = "Quadro generale dell'outcome", cfg = cfg
    )
  }

  if (.mira_report_has_section(cfg, "descriptives")) {
    .mira_report_heading("Statistiche descrittive", 2L)
    cat(
      "Le statistiche sono calcolate sui valori finiti disponibili a ogni visita. Gli intervalli di confidenza riguardano la media osservata.\n\n"
    )
    .mira_report_table(
      outcome$descriptives,
      caption = paste("Statistiche descrittive complete per", display), cfg = cfg
    )
  }

  if (.mira_report_has_section(cfg, "frequencies")) {
    .mira_report_heading("Frequenze e classificazioni", 2L)
    .mira_report_render_frequencies(outcome, cfg)
  }

  if (.mira_report_has_section(cfg, "missingness")) {
    .mira_report_heading("Completezza e dati non disponibili", 2L)
    cat(
      "MIRA distingue i valori originariamente mancanti dai valori non finiti; entrambi risultano non disponibili per l'analisi.\n\n"
    )
    .mira_report_heading("Per visita", 3L)
    .mira_report_table(
      outcome$missing$by_time,
      caption = "Missingness e disponibilità per visita", cfg = cfg
    )
    .mira_report_heading("Per soggetto", 3L)
    .mira_report_table(
      outcome$missing$by_patient,
      caption = "Completezza del profilo longitudinale per soggetto", cfg = cfg
    )
  }

  if (.mira_report_has_section(cfg, "change")) {
    .mira_report_heading("Cambiamenti entro soggetto", 2L)
    cat(
      paste0(
        "Ogni riga usa soltanto i soggetti con entrambe le visite disponibili. ",
        "Il segno è calcolato come valore successivo meno valore precedente; `cohens_dz` ",
        "standardizza il cambiamento mediante la deviazione standard delle differenze.\n\n"
      )
    )
    .mira_report_table(
      outcome$change,
      caption = "Tutti i confronti longitudinali a coppie", cfg = cfg
    )
  }

  if (.mira_report_has_section(cfg, "arms")) {
    .mira_report_heading("Analisi per braccio", 2L)
    arm <- outcome$arm_analysis
    if (is.null(arm) || !isTRUE(arm$enabled)) {
      cat("*Analisi per braccio non attiva o non stimabile.*\n\n")
      if (!is.null(arm)) {
        .mira_report_table(
          .mira_report_scalar_table(arm),
          caption = "Stato del modulo per braccio", cfg = cfg
        )
      }
    } else {
      cat(
        paste0(
          "Le differenze sono presentate con la codifica documentata nelle colonne ",
          "(`arm_b - arm_a`). I confronti al basale descrivono l'equilibrio osservato; ",
          "non devono essere usati per ridefinire post hoc il modello primario.\n\n"
        )
      )
      arm_order <- c(
        "counts", "descriptives", "baseline_balance", "missingness",
        "time_omnibus", "time_pairwise", "change_descriptives",
        "change_omnibus", "change_pairwise"
      )
      for (name in arm_order[arm_order %in% names(arm)]) {
        .mira_report_heading(.mira_report_human_name(name), 3L)
        .mira_report_render_any(arm[[name]], name, 4L, cfg)
      }
    }
  }

  if (.mira_report_has_section(cfg, "correlations")) {
    .mira_report_heading("Correlazioni tra visite", 2L)
    cat(
      paste0(
        "Le matrici usano osservazioni complete a coppie. La matrice delle numerosità ",
        "deve essere letta insieme ai coefficienti perché il campione può variare tra celle.\n\n"
      )
    )
    .mira_report_render_any(outcome$correlations, "correlations", 3L, cfg)
  }

  if (.mira_report_has_section(cfg, "variability")) {
    .mira_report_heading("Variabilità e affidabilità", 2L)
    .mira_report_table(
      outcome$variability,
      caption = "Componenti di variabilità e coefficienti ICC", cfg = cfg
    )
  }

  if (.mira_report_has_section(cfg, "models")) {
    .mira_report_heading("Inferenza longitudinale e modelli", 2L)
    .mira_report_callout(
      paste0(
        "Convergenza, singolarità, struttura di correlazione, numerosità effettiva e ",
        "assunzioni devono essere valutate prima di interpretare i coefficienti. I criteri ",
        "di confronto fra famiglie di modelli non sono sempre direttamente comparabili."
      ),
      type = "note", title = "Controllo del modello"
    )
    .mira_report_render_models(outcome, cfg)
  }

  if (.mira_report_has_section(cfg, "robustness")) {
    .mira_report_heading("Robustezza, effect size e molteplicità", 2L)
    for (name in c("robustness", "effect_sizes", "multiplicity", "sensitivity")) {
      .mira_report_heading(.mira_report_human_name(name), 3L)
      .mira_report_render_any(outcome[[name]], name, 4L, cfg)
    }
  }

  if (.mira_report_has_section(cfg, "outliers")) {
    .mira_report_heading("Valori anomali diagnostici", 2L)
    cat(
      "I flag IQR sono indicatori diagnostici. Coerentemente con `mira_info`, non comportano esclusioni automatiche.\n\n"
    )
    .mira_report_render_any(outcome$outliers, "outliers", 3L, cfg)
  }

  if (.mira_report_has_section(cfg, "trajectories")) {
    .mira_report_heading("Traiettorie individuali", 2L)
    .mira_report_table(
      outcome$trajectories,
      caption = "Risultati individuali basale-finale", cfg = cfg
    )
  }

  if (.mira_report_has_section(cfg, "figures")) {
    .mira_report_heading("Figure", 2L)
    .mira_report_render_figures(outcome, cfg)
  }

  if (.mira_report_has_section(cfg, "diagnostics")) {
    .mira_report_heading("Diagnostica specifica dell'outcome", 2L)
    .mira_report_render_any(outcome$diagnostics, "diagnostics", 3L, cfg)
  }
  invisible(NULL)
}

.mira_report_render_conclusions <- function(outcomes, cfg) {
  cat("\n\\clearpage\n\n")
  .mira_report_heading("Conclusioni e limiti interpretativi", 1L)
  completed <- outcomes[!vapply(outcomes, inherits, logical(1L),
                                what = "mira_info_error")]
  if (length(completed) == 0L) {
    .mira_report_callout(
      "Nessun outcome è stato completato; non è possibile produrre una sintesi conclusiva.",
      type = "warning", title = "Conclusioni non disponibili"
    )
    return(invisible(NULL))
  }

  conclusion_rows <- lapply(names(completed), function(name) {
    outcome <- completed[[name]]
    desc <- outcome$descriptives
    baseline <- final <- delta <- NA_real_
    direction <- "non stimabile"
    if (is.data.frame(desc) && nrow(desc) > 0L && "mean" %in% names(desc)) {
      baseline <- desc$mean[[1L]]
      final <- desc$mean[[nrow(desc)]]
      delta <- final - baseline
      threshold <- .mira_report_or(outcome$settings$stable_threshold, 0)
      direction <- if (!is.finite(delta)) "non stimabile" else if (delta > threshold) {
        "aumento"
      } else if (delta < -threshold) {
        "diminuzione"
      } else {
        "stabile entro soglia"
      }
    }
    data.frame(
      outcome = .mira_report_or(outcome$outcome_display, name),
      media_basale = baseline,
      media_finale = final,
      differenza_descrittiva = delta,
      direzione = direction,
      direzione_clinica_configurata = .mira_report_or(
        outcome$settings$improvement_direction, "unknown"
      ),
      stringsAsFactors = FALSE
    )
  })
  conclusion_table <- do.call(rbind, conclusion_rows)
  rownames(conclusion_table) <- NULL
  .mira_report_table(
    conclusion_table,
    caption = "Sintesi descrittiva finale degli outcome", cfg = cfg
  )
  cat(
    paste0(
      "La direzione riportata in tabella descrive la variazione delle medie osservate e non, ",
      "da sola, evidenza di efficacia. Le conclusioni sostantive devono integrare intervalli di ",
      "confidenza, dimensioni dell'effetto, analisi di sensibilità, pattern di missingness, ",
      "qualità dell'adattamento e rilevanza clinica definita a priori.\n\n"
    )
  )
  limitations <- c(
    "Il report riflette le analisi presenti nell'oggetto: un modulo assente non equivale a un risultato nullo.",
    "Le analisi automatiche non sostituiscono un estimand e un piano statistico prespecificati.",
    "Missingness informativa, misure non comparabili tra visite e violazioni delle assunzioni possono modificare l'interpretazione.",
    "I confronti molteplici vanno letti nelle famiglie documentate; non selezionare post hoc il solo metodo più favorevole.",
    "Associazioni e differenze osservate non autorizzano automaticamente conclusioni causali."
  )
  cat(paste0("- ", limitations, collapse = "\n"), "\n\n", sep = "")
  invisible(NULL)
}

.mira_report_render_appendix <- function(result, extra_objects, cfg) {
  cat("\n\\clearpage\n\n")
  .mira_report_heading("Appendici di audit e riproducibilità", 1L)

  .mira_report_heading("Indice completo dell'oggetto", 2L)
  cat(
    paste0(
      "L'indice seguente elenca ogni componente attraversabile dell'oggetto ricevuto, ",
      "con classe, dimensione e stato. Gli oggetti di modello complessi sono trattati come ",
      "unità atomiche e stampati tramite il rispettivo metodo.\n\n"
    )
  )
  index <- .mira_report_object_index(result, root = "result", max_depth = cfg$max_depth)
  .mira_report_table(index, caption = "Inventario completo degli output MIRA", cfg = cfg)

  if (length(extra_objects) > 0L) {
    .mira_report_heading("Indice degli oggetti supplementari", 2L)
    extra_index <- do.call(rbind, lapply(names(extra_objects), function(name) {
      .mira_report_object_index(
        extra_objects[[name]], root = paste0("extra_objects$", name),
        max_depth = cfg$max_depth
      )
    }))
    rownames(extra_index) <- NULL
    .mira_report_table(extra_index, caption = "Inventario degli oggetti supplementari", cfg = cfg)
  }

  if (isTRUE(cfg$include_complete_output)) {
    .mira_report_heading("Output integrale", 2L)
    .mira_report_callout(
      paste0(
        "Questa appendice è deliberatamente estesa: mantiene tabelle, vettori, diagnostica, ",
        "oggetti di modello e componenti non ancora conosciuti da versioni future di MIRA. ",
        "La duplicazione di alcuni elementi già discussi nel testo principale serve all'audit."
      ),
      type = "note", title = "Completezza dell'appendice"
    )
    .mira_report_render_any(result, "result", 3L, cfg)

    if (length(extra_objects) > 0L) {
      .mira_report_heading("Oggetti supplementari integrali", 2L)
      for (name in names(extra_objects)) {
        .mira_report_heading(.mira_report_human_name(name), 3L)
        .mira_report_render_any(extra_objects[[name]], name, 4L, cfg)
      }
    }
  }

  .mira_report_heading("Ambiente di calcolo", 2L)
  .mira_report_text_block(utils::sessionInfo())
  invisible(NULL)
}

.mira_report_html_style <- function() {
  cat(
    paste0(
      "```{=html}\n",
      "<style>\n",
      ":root { --mira-ink:#14213d; --mira-accent:#1f6f8b; --mira-soft:#eef5f7; }\n",
      "body { color:var(--mira-ink); line-height:1.58; }\n",
      "h1,h2,h3 { color:var(--mira-ink); letter-spacing:-0.015em; }\n",
      "h1 { border-bottom:2px solid var(--mira-accent); padding-bottom:.32rem; }\n",
      ".mira-table { width:100%; font-size:.88rem; margin:1rem 0 1.6rem 0; }\n",
      ".mira-table-wrap { width:100%; overflow-x:auto; margin:1rem 0 1.6rem 0; }\n",
      ".mira-table-wrap .mira-table { margin:0; }\n",
      ".mira-table thead { background:var(--mira-ink); color:white; }\n",
      ".mira-table tbody tr:nth-child(even) { background:#f7f9fb; }\n",
      ".mira-table td,.mira-table th { padding:.42rem .55rem; vertical-align:top; }\n",
      ".callout { border-radius:4px; }\n",
      "code { color:#7b2d26; }\n",
      "@media print { .sidebar, #TOC { display:none!important; } body { font-size:10.5pt; } }\n",
      "</style>\n",
      "```\n\n"
    )
  )
}

.mira_report_knit <- function(payload) {
  result <- payload$result
  cfg <- payload$report
  extra_objects <- .mira_report_or(payload$extra_objects, list())

  old_options <- options(
    width = max(120L, getOption("width", 80L)),
    max.print = max(1000000L, getOption("max.print", 99999L)),
    scipen = 4
  )
  on.exit(options(old_options), add = TRUE)

  if (knitr::is_html_output()) .mira_report_html_style()
  outcomes <- .mira_report_extract_outcomes(result)

  if (.mira_report_has_section(cfg, "abstract")) {
    .mira_report_render_abstract(result, outcomes, cfg)
  }
  if (.mira_report_has_section(cfg, "methods")) {
    .mira_report_render_methods(result, outcomes, cfg)
  }
  if (.mira_report_has_section(cfg, "data_quality")) {
    .mira_report_render_data_quality(result, cfg)
  }

  outcome_sections <- c(
    "overview", "descriptives", "frequencies", "missingness", "change",
    "arms", "correlations", "variability", "models", "robustness",
    "outliers", "trajectories", "figures", "diagnostics"
  )
  render_outcomes <- "all" %in% cfg$sections ||
    any(outcome_sections %in% cfg$sections)
  if (render_outcomes && length(outcomes) > 0L) {
    for (i in seq_along(outcomes)) {
      .mira_report_render_outcome(
        outcomes[[i]], names(outcomes)[[i]], index = i, cfg = cfg
      )
    }
  }

  if (.mira_report_has_section(cfg, "conclusions") && length(outcomes) > 0L) {
    .mira_report_render_conclusions(outcomes, cfg)
  }

  if (.mira_report_has_section(cfg, "appendix")) {
    .mira_report_render_appendix(result, extra_objects, cfg)
  }
  invisible(NULL)
}

.mira_report_runtime_names <- function() {
  c(
    ".mira_report_or", ".mira_report_is_empty", ".mira_report_heading",
    ".mira_report_human_name", ".mira_report_format_p",
    ".mira_report_format_number", ".mira_report_inline",
    ".mira_report_prepare_table", ".mira_report_numeric_column",
    ".mira_report_column_widths", ".mira_report_partition_columns",
    ".mira_report_latex_alignment", ".mira_report_latex_breaks",
    ".mira_report_repeat_longtable_header", ".mira_report_table",
    ".mira_report_prepare_namespace", ".mira_report_text_block",
    ".mira_report_callout",
    ".mira_report_is_plot", ".mira_report_plot",
    ".mira_report_is_plain_list", ".mira_report_render_any",
    ".mira_report_extract_outcomes", ".mira_report_first_success",
    ".mira_report_scalar_table", ".mira_report_frequency",
    ".mira_report_object_index", ".mira_report_has_section",
    ".mira_report_longitudinal_map", ".mira_report_abstract_rows",
    ".mira_report_key_narrative", ".mira_report_render_abstract",
    ".mira_report_render_methods", ".mira_report_render_data_quality",
    ".mira_report_render_frequencies", ".mira_report_render_figures",
    ".mira_report_render_models", ".mira_report_render_outcome",
    ".mira_report_render_conclusions",
    ".mira_report_render_appendix", ".mira_report_html_style",
    ".mira_report_knit"
  )
}

.mira_report_slug <- function(x) {
  x <- iconv(as.character(x)[1L], to = "ASCII//TRANSLIT", sub = "")
  x <- tolower(gsub("[^A-Za-z0-9]+", "-", x))
  x <- gsub("(^-+|-+$)", "", x)
  if (is.na(x) || !nzchar(x)) x <- "mira-report"
  x
}

.mira_report_yaml_quote <- function(x) {
  x <- gsub("[\r\n]+", " ", as.character(x)[1L])
  paste0("'", gsub("'", "''", x, fixed = TRUE), "'")
}

.mira_report_r_quote <- function(x) {
  x <- gsub("\\\\", "\\\\\\\\", as.character(x)[1L])
  x <- gsub("\"", "\\\\\"", x, fixed = TRUE)
  paste0("\"", x, "\"")
}

.mira_report_format_lines <- function(formats) {
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

.mira_report_qmd_lines <- function(title, subtitle, author, date, formats,
                                   runtime_file, payload_file) {
  yaml <- c("---", paste0("title: ", .mira_report_yaml_quote(title)))
  if (!is.null(subtitle) && nzchar(subtitle)) {
    yaml <- c(yaml, paste0("subtitle: ", .mira_report_yaml_quote(subtitle)))
  }
  if (!is.null(author) && length(author) > 0L) {
    if (length(author) == 1L) {
      yaml <- c(yaml, paste0("author: ", .mira_report_yaml_quote(author)))
    } else {
      yaml <- c(yaml, "author:", paste0("  - ", vapply(
        author, .mira_report_yaml_quote, character(1L)
      )))
    }
  }
  yaml <- c(
    yaml,
    paste0("date: ", .mira_report_yaml_quote(as.character(date))),
    "lang: it",
    .mira_report_format_lines(formats),
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
    paste0("source(", .mira_report_r_quote(runtime_file),
           ", local = knitr::knit_global())"),
    paste0("payload <- readRDS(", .mira_report_r_quote(payload_file), ")"),
    "mira_plot_device <- if (knitr::is_latex_output() && capabilities('cairo')) 'cairo_pdf' else 'png'",
    "knitr::opts_chunk$set(fig.width = 7.0, fig.height = 4.6, dpi = 320,",
    "                      dev = mira_plot_device, fig.align = 'center',",
    "                      out.width = '100%')",
    "```",
    "",
    "```{r}",
    "#| label: mira-report",
    "#| results: asis",
    ".mira_report_knit(payload)",
    "```",
    ""
  )
  yaml
}

.mira_report_validate_scalar_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop(sprintf("%s deve essere TRUE o FALSE.", name), call. = FALSE)
  }
  invisible(TRUE)
}

#' Create a dynamic, exhaustive Quarto report from mira_info()
#'
#' @param x An existing `mira_info`, `mira_info_multi`, or `mira_detect` object.
#'   A data.frame is also accepted and is treated as `data`.
#' @param data Optional data.frame. If supplied, `mira_info(data, ...)` is run first.
#' @param ... Arguments forwarded unchanged to `mira_info()` when `data` is used.
#' @param output_dir Directory in which the Quarto source, RDS payload, runtime
#'   helper, and rendered reports are written. The default creates a timestamped
#'   directory under the current working directory.
#' @param output_file Base filename without extension.
#' @param format One or more of `"html"`, `"pdf"`, `"docx"`, or `"all"`.
#' @param title,subtitle,author,date Report metadata.
#' @param sections `"all"` or any subset of the documented report sections.
#' @param include_complete_output If TRUE, append a recursive rendering of the
#'   complete result object. This is TRUE by default to avoid silent omissions.
#' @param extra_objects Optional named list of additional fitted objects (for
#'   example `list(mira_fit = fit)`) to inventory and include in the appendix.
#' @param max_table_rows Maximum rows printed per table. Defaults to Inf: there is
#'   no arbitrary row or page limit. A finite value must be explicitly requested.
#' @param max_table_columns Maximum columns per displayed block. Wide tables are
#'   split into blocks; columns are never dropped.
#' @param digits Number of display digits.
#' @param max_depth Safety depth for recursive object traversal.
#' @param render If TRUE, render with Quarto; if FALSE, only create reproducible
#'   `.qmd`, `.R`, and `.rds` source files.
#' @param open Open the first rendered report in an interactive session.
#' @param overwrite Allow existing source or output files to be replaced.
#' @param quiet Passed to `quarto::quarto_render()`.
#'
#' @return Invisibly, an object of class `mira_report_freq` containing the MIRA
#'   result and all generated paths.
mira_report_freq <- function(
    x = NULL,
    data = NULL,
    ...,
    output_dir = NULL,
    output_file = "mira_report",
    format = c("html", "pdf"),
    title = "MIRA: analisi longitudinale",
    subtitle = "Report statistico dinamico e riproducibile",
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
    .mira_report_validate_scalar_flag(get(flag, inherits = FALSE), flag)
  }

  for (metadata_name in c("output_file", "title")) {
    value <- get(metadata_name, inherits = FALSE)
    if (!is.character(value) || length(value) != 1L || is.na(value) ||
        !nzchar(trimws(value))) {
      stop(sprintf("%s deve essere una stringa non vuota.", metadata_name),
           call. = FALSE)
    }
  }
  if (!is.null(subtitle) &&
      (!is.character(subtitle) || length(subtitle) != 1L || is.na(subtitle))) {
    stop("subtitle deve essere NULL o una sola stringa.", call. = FALSE)
  }
  if (!is.null(author) &&
      (!is.character(author) || anyNA(author) || any(!nzchar(trimws(author))))) {
    stop("author deve essere NULL o un vettore di stringhe non vuote.", call. = FALSE)
  }
  if (length(date) != 1L || is.na(date)) {
    stop("date deve avere lunghezza uno e non essere mancante.", call. = FALSE)
  }

  allowed_sections <- c(
    "all", "abstract", "methods", "data_quality", "overview",
    "descriptives", "frequencies", "missingness", "change", "arms",
    "correlations", "variability", "models", "robustness", "outliers",
    "trajectories", "figures", "diagnostics", "conclusions", "appendix"
  )
  if (!is.character(sections) || length(sections) == 0L || anyNA(sections)) {
    stop("sections deve essere 'all' o un vettore di nomi di sezione.", call. = FALSE)
  }
  sections <- unique(tolower(sections))
  invalid_sections <- setdiff(sections, allowed_sections)
  if (length(invalid_sections) > 0L) {
    stop(sprintf(
      "Sezioni non riconosciute: %s. Valori ammessi: %s.",
      paste(invalid_sections, collapse = ", "),
      paste(allowed_sections, collapse = ", ")
    ), call. = FALSE)
  }

  allowed_formats <- c("html", "pdf", "docx")
  if (!is.character(format) || length(format) == 0L || anyNA(format)) {
    stop("format deve contenere html, pdf, docx oppure all.", call. = FALSE)
  }
  formats <- unique(tolower(format))
  if ("all" %in% formats) formats <- allowed_formats
  invalid_formats <- setdiff(formats, allowed_formats)
  if (length(invalid_formats) > 0L) {
    stop(sprintf("Formati non riconosciuti: %s.",
                 paste(invalid_formats, collapse = ", ")), call. = FALSE)
  }

  if (!is.numeric(max_table_rows) || length(max_table_rows) != 1L ||
      is.na(max_table_rows) || max_table_rows <= 0) {
    stop("max_table_rows deve essere un numero positivo oppure Inf.", call. = FALSE)
  }
  if (!is.numeric(max_table_columns) || length(max_table_columns) != 1L ||
      is.na(max_table_columns) || !is.finite(max_table_columns) ||
      max_table_columns < 2) {
    stop("max_table_columns deve essere un intero finito >= 2.", call. = FALSE)
  }
  if (!is.numeric(digits) || length(digits) != 1L || is.na(digits) ||
      !is.finite(digits) || digits < 1 || digits > 10) {
    stop("digits deve essere un intero tra 1 e 10.", call. = FALSE)
  }
  if (!is.numeric(max_depth) || length(max_depth) != 1L || is.na(max_depth) ||
      !is.finite(max_depth) || max_depth < 1) {
    stop("max_depth deve essere un intero finito >= 1.", call. = FALSE)
  }
  if (!is.list(extra_objects)) {
    stop("extra_objects deve essere una lista, preferibilmente nominata.", call. = FALSE)
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
      stop("Fornire i dati una sola volta: in x oppure in data.", call. = FALSE)
    }
    data <- x
    x <- NULL
  }

  if (!is.null(data)) {
    if (!is.null(x)) {
      stop("Fornire un oggetto MIRA in x oppure un data.frame in data, non entrambi.",
           call. = FALSE)
    }
    if (!is.data.frame(data)) stop("data deve essere un data.frame.", call. = FALSE)
    mira_fun <- get0("mira_info", mode = "function", inherits = TRUE)
    if (is.null(mira_fun)) {
      stop("mira_info() non è disponibile: caricare prima il file che la definisce.",
           call. = FALSE)
    }
    if ("data" %in% names(dots)) {
      stop("Non ripetere data dentro ...; usare l'argomento data.", call. = FALSE)
    }
    if (!"verbose" %in% names(dots)) dots$verbose <- FALSE
    x <- do.call(mira_fun, c(list(data = data), dots))
  } else if (length(dots) > 0L) {
    stop("Gli argomenti in ... sono ammessi solo quando viene fornito data.", call. = FALSE)
  }

  if (is.null(x)) {
    stop("Fornire x (output di mira_info) oppure data.", call. = FALSE)
  }
  accepted <- inherits(x, c("mira_info", "mira_info_multi", "mira_detect"))
  if (!accepted) {
    stop(
      "x deve ereditare da mira_info, mira_info_multi o mira_detect.",
      call. = FALSE
    )
  }

  output_file <- .mira_report_slug(output_file)
  if (is.null(output_dir)) {
    timestamp <- base::format(Sys.time(), "%Y%m%d-%H%M%S")
    output_dir <- file.path(getwd(), paste0(output_file, "-", timestamp))
  }
  if (!is.character(output_dir) || length(output_dir) != 1L ||
      is.na(output_dir) || !nzchar(output_dir)) {
    stop("output_dir deve essere un percorso non vuoto.", call. = FALSE)
  }
  if (!dir.exists(output_dir)) {
    created <- dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    if (!created && !dir.exists(output_dir)) {
      stop(sprintf("Impossibile creare la directory: %s", output_dir), call. = FALSE)
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

  protected_paths <- c(qmd_path, payload_path, runtime_path, unname(output_paths))
  existing <- protected_paths[file.exists(protected_paths)]
  if (length(existing) > 0L && !overwrite) {
    stop(sprintf(
      "Esistono già file di destinazione. Usare overwrite=TRUE o un'altra directory: %s",
      paste(basename(existing), collapse = ", ")
    ), call. = FALSE)
  }

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

  saveRDS(payload, payload_path, compress = "xz")
  runtime_names <- .mira_report_runtime_names()
  runtime_env <- environment(mira_report_freq)
  missing_runtime <- runtime_names[!vapply(runtime_names, exists, logical(1L),
                                           envir = runtime_env, inherits = TRUE)]
  if (length(missing_runtime) > 0L) {
    stop(sprintf("Helper interni mancanti: %s.",
                 paste(missing_runtime, collapse = ", ")), call. = FALSE)
  }
  dump(runtime_names, file = runtime_path, envir = runtime_env)

  qmd_lines <- .mira_report_qmd_lines(
    title = title, subtitle = subtitle, author = author, date = date,
    formats = formats, runtime_file = runtime_name, payload_file = payload_name
  )
  writeLines(qmd_lines, qmd_path, useBytes = TRUE)

  rendered <- stats::setNames(rep(FALSE, length(formats)), formats)
  render_errors <- stats::setNames(rep(NA_character_, length(formats)), formats)
  if (render) {
    if (!requireNamespace("knitr", quietly = TRUE)) {
      stop(
        paste0("Il pacchetto R 'knitr' e necessario. I file sorgente sono stati creati in: ",
               output_dir),
        call. = FALSE
      )
    }
    if (!requireNamespace("quarto", quietly = TRUE)) {
      stop(
        paste0("Il pacchetto R 'quarto' e necessario. I file sorgente sono stati creati in: ",
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
        "Rendering non completato per %s. Sorgenti conservati. %s",
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
      output_dir = output_dir,
      report_config = report_config
    ),
    class = c("mira_report_freq", "list")
  )

  if (open && any(rendered)) {
    first_file <- unname(output_paths[which(rendered)[1L]])
    try(utils::browseURL(first_file), silent = TRUE)
  }
  print(result)
  invisible(result)
}

#' @export
print.mira_report_freq <- function(x, ...) {
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
