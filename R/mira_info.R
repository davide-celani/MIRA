

.mira_key <- function(x) {
  tolower(gsub("[^[:alnum:]]+", "_", trimws(as.character(x))))
}

.mira_compact_unique <- function(x) {
  unique(x[!is.na(x) & nzchar(x)])
}

.mira_validate_name <- function(x, data, argument, allow_null = FALSE) {
  if (allow_null && is.null(x)) return(NULL)
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(trimws(x))) {
    stop(sprintf("%s must name exactly one variable.", argument), call. = FALSE)
  }
  if (!x %in% names(data)) {
    stop(sprintf("The variable '%s' specified by %s does not exist in the data set.", x, argument),
         call. = FALSE)
  }
  x
}

.mira_column_profile <- function(data) {
  n <- nrow(data)
  data.frame(
    variable = names(data),
    class = vapply(data, function(x) paste(class(x), collapse = "/"), character(1L)),
    numeric = vapply(data, is.numeric, logical(1L)),
    unique_n = vapply(data, function(x) length(unique(x[!is.na(x)])), integer(1L)),
    missing_n = vapply(data, function(x) sum(is.na(x)), integer(1L)),
    missing_pct = vapply(data, function(x) {
      if (n == 0L) NA_real_ else sum(is.na(x)) / n * 100
    }, numeric(1L)),
    stringsAsFactors = FALSE
  )
}

.mira_time_order <- function(token) {
  key <- .mira_key(token)
  if (grepl("^(baseline|base|bl)$", key)) return(0)
  number <- regmatches(key, regexpr("[0-9]+(?:[.][0-9]+)?", key, perl = TRUE))
  if (length(number) == 0L || !nzchar(number)) return(NA_real_)
  suppressWarnings(as.numeric(number))
}

.mira_parse_longitudinal_name <- function(variable, variable_pattern = "auto") {
  if (is.function(variable_pattern)) {
    parsed <- variable_pattern(variable)
    if (is.null(parsed) || length(parsed) == 0L) return(NULL)
    if (is.data.frame(parsed)) parsed <- as.list(parsed[1L, , drop = FALSE])
    if (!is.list(parsed) || is.null(parsed$outcome)) {
      stop(
        "The variable_pattern function must return NULL or a list containing at least 'outcome'.",
        call. = FALSE
      )
    }
    label <- if (is.null(parsed$time_label)) variable else as.character(parsed$time_label)[1L]
    order <- if (is.null(parsed$time_order)) .mira_time_order(label) else {
      suppressWarnings(as.numeric(parsed$time_order)[1L])
    }
    return(list(
      outcome = as.character(parsed$outcome)[1L],
      time_label = label,
      time_order = order,
      pattern = "custom_function"
    ))
  }

  if (!is.character(variable_pattern) || length(variable_pattern) != 1L ||
      is.na(variable_pattern) || !nzchar(variable_pattern)) {
    stop("variable_pattern must be 'auto', a regular expression, or a function.", call. = FALSE)
  }

  if (!identical(variable_pattern, "auto")) {
    match <- regexec(variable_pattern, variable)
    pieces <- regmatches(variable, match)[[1L]]
    if (length(pieces) < 3L) return(NULL)
    return(list(
      outcome = pieces[[2L]],
      time_label = pieces[[3L]],
      time_order = .mira_time_order(pieces[[3L]]),
      pattern = "custom_regex"
    ))
  }

  patterns <- list(
    t = "^(.+)[._-]([tT][0-9]+)$",
    time = "^(.+)[._-]((time|visit|vis|v)[._-]?[0-9]+)$",
    baseline = "^(.+)[._-](baseline|base|BL|bl)$",
    month = "^(.+)[._-]((month|months|mo|m|M)[._-]?[0-9]+)$",
    week = "^(.+)[._-]((week|weeks|wk|w|W)[._-]?[0-9]+)$",
    day = "^(.+)[._-]((day|days|d|D)[._-]?[0-9]+)$",
    followup = "^(.+)[._-]((follow[_-]?up|followup|fu|FU)[._-]?[0-9]+)$",
    numeric = "^(.+)[._-]([0-9]+)$"
  )

  for (pattern_name in names(patterns)) {
    match <- regexec(patterns[[pattern_name]], variable, ignore.case = TRUE)
    pieces <- regmatches(variable, match)[[1L]]
    if (length(pieces) >= 3L && nzchar(pieces[[2L]])) {
      return(list(
        outcome = pieces[[2L]],
        time_label = pieces[[3L]],
        time_order = .mira_time_order(pieces[[3L]]),
        pattern = pattern_name
      ))
    }
  }
  NULL
}

.mira_make_longitudinal_group <- function(data,
                                          variables,
                                          outcome = NULL,
                                          variable_pattern = "auto",
                                          automatic = FALSE) {
  variables <- as.character(variables)
  if (length(variables) < 2L || anyNA(variables) || any(!nzchar(variables))) {
    stop("Each time_vars group must contain at least two non-empty names.", call. = FALSE)
  }
  if (anyDuplicated(variables)) {
    stop("time_vars contains duplicate names.", call. = FALSE)
  }
  missing <- setdiff(variables, names(data))
  if (length(missing) > 0L) {
    stop(sprintf("Longitudinal variables not found: %s.", paste(missing, collapse = ", ")),
         call. = FALSE)
  }
  non_numeric <- variables[!vapply(data[variables], is.numeric, logical(1L))]
  if (length(non_numeric) > 0L) {
    stop(sprintf("Longitudinal variables must be numeric: %s.",
                 paste(non_numeric, collapse = ", ")), call. = FALSE)
  }

  parsed <- lapply(variables, .mira_parse_longitudinal_name,
                   variable_pattern = variable_pattern)
  parsed_outcomes <- vapply(parsed, function(x) {
    if (is.null(x)) NA_character_ else x$outcome
  }, character(1L))
  parsed_labels <- vapply(seq_along(parsed), function(i) {
    if (is.null(parsed[[i]])) variables[[i]] else parsed[[i]]$time_label
  }, character(1L))
  parsed_order <- vapply(parsed, function(x) {
    if (is.null(x)) NA_real_ else x$time_order
  }, numeric(1L))

  if (is.null(outcome) || !nzchar(outcome)) {
    candidates <- .mira_compact_unique(parsed_outcomes)
    if (length(candidates) == 1L) {
      outcome <- candidates[[1L]]
    } else {
      common <- variables[[1L]]
      if (length(variables) > 1L) {
        chars <- strsplit(variables, "", fixed = TRUE)
        min_len <- min(lengths(chars))
        common_len <- 0L
        if (min_len > 0L) {
          for (k in seq_len(min_len)) {
            if (length(unique(vapply(chars, `[[`, character(1L), k))) == 1L) {
              common_len <- k
            } else break
          }
        }
        common <- if (common_len > 0L) substr(common, 1L, common_len) else ""
      }
      common <- sub("[._-]+$", "", common)
      outcome <- if (nzchar(common)) common else "outcome"
    }
  }

  # Explicit variables remain authoritative. When every suffix has a clear,
  # unique order, reproduce v3's chronological sorting; otherwise preserve the
  # user's order and expose the ambiguity in diagnostics.
  clear_order <- all(is.finite(parsed_order)) && !anyDuplicated(parsed_order)
  ord <- if (clear_order) order(parsed_order, variables) else seq_along(variables)

  list(
    outcome = as.character(outcome)[1L],
    variables = variables[ord],
    input_variables = variables,
    sort_index = ord,
    time_labels = parsed_labels[ord],
    time_order = parsed_order[ord],
    automatic = automatic,
    ambiguous_order = !clear_order,
    patterns = vapply(parsed[ord], function(x) {
      if (is.null(x)) NA_character_ else x$pattern
    }, character(1L))
  )
}

.mira_detect_longitudinal_variables <- function(data, variable_pattern = "auto") {
  parsed_rows <- list()
  counter <- 1L

  for (position in seq_along(data)) {
    parsed <- .mira_parse_longitudinal_name(names(data)[[position]], variable_pattern)
    if (is.null(parsed)) next
    parsed_rows[[counter]] <- data.frame(
      outcome = parsed$outcome,
      outcome_key = .mira_key(parsed$outcome),
      variable = names(data)[[position]],
      time_label = parsed$time_label,
      time_order = parsed$time_order,
      position = position,
      pattern = parsed$pattern,
      numeric = is.numeric(data[[position]]),
      stringsAsFactors = FALSE
    )
    counter <- counter + 1L
  }

  empty_table <- data.frame(
    outcome = character(0), outcome_key = character(0), variable = character(0),
    time_label = character(0), time_order = numeric(0), position = integer(0),
    pattern = character(0), numeric = logical(0), stringsAsFactors = FALSE
  )
  parsed_table <- if (length(parsed_rows) == 0L) empty_table else do.call(rbind, parsed_rows)
  numeric_table <- parsed_table[parsed_table$numeric, , drop = FALSE]
  non_numeric <- parsed_table$variable[!parsed_table$numeric]

  groups <- list()
  ambiguous_groups <- list()
  if (nrow(numeric_table) > 0L) {
    for (key in unique(numeric_table$outcome_key)) {
      rows <- numeric_table[numeric_table$outcome_key == key, , drop = FALSE]
      if (nrow(rows) < 2L) next
      duplicate_label <- anyDuplicated(tolower(rows$time_label)) > 0L
      finite_order <- is.finite(rows$time_order)
      duplicate_order <- anyDuplicated(rows$time_order[finite_order]) > 0L
      ambiguous <- duplicate_label || duplicate_order
      fallback_order <- rows$time_order
      if (any(!is.finite(fallback_order))) {
        max_known <- if (any(is.finite(fallback_order))) max(fallback_order[is.finite(fallback_order)]) else 0
        fallback_order[!is.finite(fallback_order)] <-
          max_known + seq_len(sum(!is.finite(fallback_order)))
      }
      ord <- order(fallback_order, rows$position)
      group <- list(
        outcome = rows$outcome[[1L]],
        variables = rows$variable[ord],
        input_variables = rows$variable,
        sort_index = ord,
        time_labels = rows$time_label[ord],
        time_order = rows$time_order[ord],
        automatic = TRUE,
        ambiguous_order = ambiguous,
        patterns = rows$pattern[ord]
      )
      if (ambiguous) {
        ambiguous_groups[[group$outcome]] <- group
      } else {
        groups[[group$outcome]] <- group
      }
    }
  }

  list(
    groups = groups,
    ambiguous_groups = ambiguous_groups,
    parsed = parsed_table,
    non_numeric_matches = non_numeric
  )
}

.mira_detect_id <- function(data, id = NULL, exclude = character(0)) {
  if (!is.null(id)) {
    selected <- .mira_validate_name(id, data, "id")
    return(list(
      selected = selected, automatic = FALSE, source = "manual",
      candidates = data.frame(), alternatives = character(0), warnings = character(0)
    ))
  }

  keys <- .mira_key(names(data))
  priority <- c(
    patient_id = 120, subject_id = 120, participant_id = 120,
    usubjid = 120, subjid = 115, patid = 115, record_id = 110,
    study_id = 105, patient = 100, subject = 100, participant = 100,
    id = 95
  )
  semantic <- unname(priority[keys])
  semantic[is.na(semantic) & grepl("_id$", keys)] <- 70
  semantic[is.na(semantic)] <- 0

  observed_n <- vapply(data, function(x) sum(!is.na(x)), integer(1L))
  unique_n <- vapply(data, function(x) length(unique(x[!is.na(x)])), integer(1L))
  missing_n <- vapply(data, function(x) sum(is.na(x)), integer(1L))
  duplicated_n <- observed_n - unique_n
  structural <- observed_n > 0L & duplicated_n == 0L & missing_n == 0L
  score <- semantic + ifelse(structural, 20, 0)
  keep <- names(data) %in% setdiff(names(data), exclude) & (semantic > 0 | structural)

  candidates <- data.frame(
    variable = names(data)[keep],
    semantic_score = semantic[keep],
    score = score[keep],
    observed_n = observed_n[keep],
    unique_n = unique_n[keep],
    missing_n = missing_n[keep],
    duplicated_n = duplicated_n[keep],
    structurally_valid = structural[keep],
    stringsAsFactors = FALSE
  )
  if (nrow(candidates) > 0L) {
    candidates <- candidates[order(-candidates$score, candidates$variable), , drop = FALSE]
    rownames(candidates) <- NULL
  }

  eligible <- candidates[candidates$structurally_valid, , drop = FALSE]
  warnings <- character(0)
  selected <- NULL
  source <- "generated"
  alternatives <- character(0)
  if (nrow(eligible) > 0L) {
    best_score <- max(eligible$score)
    best <- eligible$variable[eligible$score == best_score]
    # Structural uniqueness is useful evidence but, without a semantic ID name,
    # it is not sufficient to promote an arbitrary measurement to subject ID.
    if (length(best) == 1L && best_score >= 90) {
      selected <- best[[1L]]
      source <- "auto"
      alternatives <- setdiff(eligible$variable, selected)
      if (length(alternatives) > 0L) {
        warnings <- c(warnings, sprintf(
          "ID detected as '%s'; plausible alternatives: %s. Specify id= to override.",
          selected, paste(alternatives, collapse = ", ")
        ))
      }
    } else {
      alternatives <- eligible$variable
      warnings <- c(warnings, sprintf(
        "Ambiguous ID (%s): an internal row identifier will be used. Specify id= to choose.",
        paste(best, collapse = ", ")
      ))
    }
  } else {
    warnings <- c(warnings,
                  "No unique ID was detected: an internal row identifier will be used.")
  }

  list(
    selected = selected,
    automatic = TRUE,
    source = source,
    candidates = candidates,
    alternatives = alternatives,
    warnings = warnings
  )
}

.mira_detect_arm <- function(data, arm = NULL, exclude = character(0)) {
  if (!is.null(arm)) {
    selected <- .mira_validate_name(arm, data, "arm")
    return(list(
      selected = selected, automatic = FALSE, source = "manual",
      candidates = data.frame(), alternatives = character(0), warnings = character(0)
    ))
  }

  keys <- .mira_key(names(data))
  priority <- c(
    arm = 120, treatment_arm = 115, treatment_group = 110,
    randomized_group = 108, randomised_group = 108, treatment = 100,
    intervention = 95, group = 85, cohort = 75
  )
  semantic <- unname(priority[keys])
  semantic[is.na(semantic)] <- 0
  unique_n <- vapply(data, function(x) length(unique(x[!is.na(x)])), integer(1L))
  categorical <- vapply(data, function(x) is.factor(x) || is.character(x) || is.logical(x),
                        logical(1L))
  plausible_structure <- unique_n >= 2L & unique_n <= max(10L, floor(nrow(data) / 4L)) &
    (categorical | unique_n <= 10L)
  keep <- !names(data) %in% exclude & semantic > 0 & plausible_structure
  rejected <- names(data)[!names(data) %in% exclude & semantic > 0 & !plausible_structure]
  candidates <- data.frame(
    variable = names(data)[keep], semantic_score = semantic[keep],
    levels_n = unique_n[keep], stringsAsFactors = FALSE
  )
  if (nrow(candidates) > 0L) {
    candidates <- candidates[order(-candidates$semantic_score, candidates$variable), , drop = FALSE]
    rownames(candidates) <- NULL
  }

  warnings <- character(0)
  selected <- NULL
  alternatives <- character(0)
  if (nrow(candidates) > 0L) {
    best_score <- max(candidates$semantic_score)
    best <- candidates$variable[candidates$semantic_score == best_score]
    if (length(best) == 1L) {
      selected <- best[[1L]]
      alternatives <- setdiff(candidates$variable, selected)
      if (length(alternatives) > 0L) {
        warnings <- c(warnings, sprintf(
          "Treatment-arm variable detected as '%s'; alternatives: %s. Specify arm= to override.",
          selected, paste(alternatives, collapse = ", ")
        ))
      }
    } else {
      alternatives <- candidates$variable
      warnings <- c(warnings, sprintf(
        "Ambiguous treatment-arm variable (%s): between-group analyses remain disabled until arm= is specified.",
        paste(best, collapse = ", ")
      ))
    }
  }
  if (length(rejected) > 0L) {
    warnings <- c(warnings, sprintf(
      paste0("Treatment-arm candidates were excluded because they lacked at least two usable groups ",
             "or had implausible cardinality: %s."),
      paste(rejected, collapse = ", ")
    ))
  }

  list(
    selected = selected, automatic = TRUE,
    source = if (is.null(selected)) "none" else "auto",
    candidates = candidates, alternatives = alternatives, warnings = warnings
  )
}

.mira_classify_covariates <- function(data, variables) {
  numeric <- character(0)
  categorical <- character(0)
  for (variable in variables) {
    x <- data[[variable]]
    unique_n <- length(unique(x[!is.na(x)]))
    if (unique_n < 2L) next
    if (is.factor(x) || is.character(x) || is.logical(x) ||
        (is.numeric(x) && unique_n <= 5L)) {
      categorical <- c(categorical, variable)
    } else if (is.numeric(x)) {
      numeric <- c(numeric, variable)
    }
  }
  list(numeric = numeric, categorical = categorical)
}

.mira_detect_covariates <- function(data, covariates = NULL, exclude = character(0)) {
  reserved <- c(
    "patient", "outcome", "time", "time_index", "time_label", "arm", "value",
    "patient_factor", "time_factor", "arm_factor"
  )
  reserved_present <- intersect(names(data), reserved)
  if (!is.null(covariates) && !(length(covariates) == 1L && identical(covariates, "auto"))) {
    if (!is.character(covariates) || anyNA(covariates) || any(!nzchar(covariates))) {
      stop("covariates must be NULL, 'auto', or a vector of names.", call. = FALSE)
    }
    if (anyDuplicated(covariates)) stop("covariates contains duplicate names.", call. = FALSE)
    missing <- setdiff(covariates, names(data))
    if (length(missing) > 0L) {
      stop(sprintf("Covariates not found: %s.", paste(missing, collapse = ", ")),
           call. = FALSE)
    }
    overlap <- intersect(covariates, exclude)
    if (length(overlap) > 0L) {
      stop(sprintf("These covariates are already used as the ID, treatment arm, or outcome: %s.",
                   paste(overlap, collapse = ", ")), call. = FALSE)
    }
    reserved_overlap <- intersect(covariates, reserved_present)
    if (length(reserved_overlap) > 0L) {
      stop(sprintf(
        "Covariate names reserved for internal long-format and model data: %s.",
        paste(reserved_overlap, collapse = ", ")
      ), call. = FALSE)
    }
    classes <- .mira_classify_covariates(data, covariates)
    return(list(
      selected = covariates, numeric = classes$numeric,
      categorical = classes$categorical, detected = covariates,
      detected_numeric = classes$numeric,
      detected_categorical = classes$categorical,
      automatic = FALSE, warnings = character(0)
    ))
  }

  available <- setdiff(names(data), unique(c(exclude, reserved_present)))
  # Character columns with many distinct values are more likely identifiers or
  # free text than adjustment variables and are therefore not auto-selected.
  candidate <- available[vapply(data[available], function(x) {
    unique_n <- length(unique(x[!is.na(x)]))
    if (unique_n < 2L) return(FALSE)
    if (is.numeric(x)) return(TRUE)
    (is.factor(x) || is.character(x) || is.logical(x)) &&
      unique_n <= max(10L, floor(nrow(data) / 4L))
  }, logical(1L))]
  candidate <- candidate[!grepl("(^id$|_id$|^id_)", .mira_key(candidate))]
  classes <- .mira_classify_covariates(data, candidate)

  df_cost <- vapply(candidate, function(variable) {
    x <- data[[variable]]
    if (variable %in% classes$categorical) {
      max(1L, length(unique(x[!is.na(x)])) - 1L)
    } else 1L
  }, integer(1L))
  max_auto_df <- floor(nrow(data) / 10L)
  warnings <- character(0)
  selected <- candidate
  if (sum(df_cost) > max_auto_df) {
    selected <- character(0)
    if (length(candidate) > 0L) {
      warnings <- sprintf(
        paste0("Detected %d candidate covariates (%s), but their complexity exceeds ",
               "the conservative data-driven limit (%d df): they will not be included automatically ",
               "in the model. Specify covariates= to choose."),
        length(candidate), paste(candidate, collapse = ", "), max_auto_df
      )
    }
  }

  selected_classes <- .mira_classify_covariates(data, selected)
  list(
    selected = selected,
    numeric = selected_classes$numeric,
    categorical = selected_classes$categorical,
    detected = candidate,
    detected_numeric = classes$numeric,
    detected_categorical = classes$categorical,
    automatic = TRUE,
    warnings = warnings
  )
}

# ============================================================
# ADVANCED LONGITUDINAL INFERENCE HELPERS
# ============================================================

.mira_skipped_module <- function(package = NA_character_, reason) {
  list(
    performed = FALSE,
    package = package,
    object = NULL,
    summary = NULL,
    tidy = NULL,
    warnings = character(0),
    error = NULL,
    reason_skipped = reason
  )
}

.mira_eval <- function(expr) {
  captured_warnings <- character(0)
  captured_error <- NULL
  value <- tryCatch(
    withCallingHandlers(
      force(expr),
      warning = function(w) {
        captured_warnings <<- unique(c(captured_warnings, conditionMessage(w)))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      captured_error <<- conditionMessage(e)
      NULL
    }
  )
  list(value = value, warnings = captured_warnings, error = captured_error)
}

.mira_quote_formula_name <- function(x) {
  paste0("`", gsub("`", "", x, fixed = TRUE), "`")
}

.mira_formula_text <- function(x) {
  if (is.null(x)) return(NA_character_)
  paste(deparse(x, width.cutoff = 500L), collapse = " ")
}

.mira_make_fixed_formula <- function(arm_enabled, covariates = character(0)) {
  core <- if (isTRUE(arm_enabled)) "time_factor * arm_factor" else "time_factor"
  covariate_terms <- if (length(covariates) > 0L) {
    vapply(covariates, .mira_quote_formula_name, character(1L))
  } else character(0)
  stats::as.formula(paste(
    "value ~",
    paste(c(core, covariate_terms), collapse = " + ")
  ))
}

.mira_make_lmer_formula <- function(core, covariates, random_term) {
  covariate_terms <- if (length(covariates) > 0L) {
    vapply(covariates, .mira_quote_formula_name, character(1L))
  } else character(0)
  fixed_terms <- c(core, covariate_terms)
  fixed_terms <- fixed_terms[nzchar(fixed_terms)]
  if (length(fixed_terms) == 0L) fixed_terms <- "1"
  stats::as.formula(paste(
    "value ~",
    paste(c(fixed_terms, random_term), collapse = " + ")
  ))
}

.mira_as_numeric <- function(x) {
  if (is.null(x)) return(numeric(0))
  if (is.numeric(x)) return(as.numeric(x))
  x <- trimws(as.character(x))
  x <- sub("^[<=>]+", "", x)
  x <- gsub("[^0-9eE+.-]", "", x)
  suppressWarnings(as.numeric(x))
}

.mira_column_name <- function(x, candidates = character(0), contains = character(0)) {
  if (is.null(x) || length(names(x)) == 0L) return(NA_character_)
  clean <- tolower(gsub("[^[:alnum:]]+", "", names(x)))
  candidate_clean <- tolower(gsub("[^[:alnum:]]+", "", candidates))
  exact <- match(candidate_clean, clean, nomatch = 0L)
  exact <- exact[exact > 0L]
  if (length(exact) > 0L) return(names(x)[exact[[1L]]])
  for (pattern in contains) {
    idx <- grep(pattern, clean, perl = TRUE)
    if (length(idx) > 0L) return(names(x)[idx[[1L]]])
  }
  NA_character_
}

.mira_adjust_family <- function(table, p_col = "p_raw", primary = "holm") {
  if (is.null(table)) return(NULL)
  table <- as.data.frame(table, stringsAsFactors = FALSE)
  if (!p_col %in% names(table)) {
    detected <- .mira_column_name(
      table,
      candidates = c("p.value", "p_value", "p", "Pr(>|t|)", "Pr(>|z|)",
                     "Pr(>|W|)", "Pr(>F)")
    )
    if (!is.na(detected)) {
      table$p_raw <- .mira_as_numeric(table[[detected]])
      p_col <- "p_raw"
    } else {
      table$p_raw <- rep(NA_real_, nrow(table))
      p_col <- "p_raw"
    }
  }
  p <- .mira_as_numeric(table[[p_col]])
  adjust <- function(method) {
    out <- rep(NA_real_, length(p))
    ok <- is.finite(p)
    if (any(ok)) out[ok] <- stats::p.adjust(p[ok], method = method)
    out
  }
  table$p_raw <- p
  table$p_primary <- adjust(primary)
  table$primary_method <- rep(primary, nrow(table))
  table$p_bonferroni <- adjust("bonferroni")
  table$p_holm <- adjust("holm")
  table$p_bh <- adjust("BH")
  table$p_by <- adjust("BY")
  table
}

.mira_emm_summary <- function(object, primary, adjust = "none", confidence_level = 0.95) {
  if (is.null(object)) return(NULL)
  tab <- as.data.frame(
    summary(object, infer = c(TRUE, TRUE), adjust = adjust, level = confidence_level),
    stringsAsFactors = FALSE
  )
  if (identical(adjust, "none")) .mira_adjust_family(tab, primary = primary) else tab
}

.mira_prepare_advanced_data <- function(long_data,
                                        time_levels,
                                        arm_enabled,
                                        arm_levels,
                                        covariates,
                                        categorical_covariates) {
  out <- long_data[
    !is.na(long_data$patient) & is.finite(long_data$value),
    ,
    drop = FALSE
  ]
  out$patient_factor <- factor(out$patient)
  out$time_factor <- factor(as.character(out$time_label), levels = time_levels)
  out$time_ordinal <- as.numeric(out$time_index) - 1
  if (isTRUE(arm_enabled)) {
    out$arm_factor <- factor(as.character(out$arm), levels = arm_levels)
  }
  for (variable in covariates) {
    value <- out[[variable]]
    if (is.numeric(value)) value[!is.finite(value)] <- NA_real_
    if (variable %in% categorical_covariates && !is.factor(value)) value <- factor(value)
    if (is.character(value)) {
      value[!is.na(value) & !nzchar(trimws(value))] <- NA_character_
      value <- factor(value)
    }
    if (is.logical(value)) value <- factor(value)
    if (is.factor(value)) value <- droplevels(value)
    out[[variable]] <- value
  }
  required <- c(
    "patient_factor", "time_factor", "time_ordinal", "value",
    if (isTRUE(arm_enabled)) "arm_factor" else character(0),
    covariates
  )
  out <- out[stats::complete.cases(out[required]), , drop = FALSE]
  out$patient_factor <- droplevels(out$patient_factor)
  out$time_factor <- droplevels(out$time_factor)
  if (isTRUE(arm_enabled)) out$arm_factor <- droplevels(out$arm_factor)
  out <- out[order(out$patient_factor, out$time_index), , drop = FALSE]
  rownames(out) <- NULL
  out
}

.mira_run_emmeans <- function(model, arm_enabled, primary_method, confidence_level) {
  if (is.null(model)) {
    return(.mira_skipped_module(
      "emmeans",
      "The compatible mixed model was not available."
    ))
  }
  if (!requireNamespace("emmeans", quietly = TRUE)) {
    return(.mira_skipped_module(
      "emmeans",
      "Optional package 'emmeans' is not installed."
    ))
  }

  warnings <- character(0)
  errors <- list()
  safe <- function(name, expr) {
    evaluated <- .mira_eval(expr)
    warnings <<- unique(c(warnings, evaluated$warnings))
    if (!is.null(evaluated$error)) errors[[name]] <<- evaluated$error
    evaluated$value
  }

  time_emm <- safe("time", emmeans::emmeans(model, ~ time_factor))
  if (is.null(time_emm)) {
    return(list(
      performed = FALSE, package = "emmeans", object = NULL, summary = NULL,
      tidy = NULL, warnings = warnings,
      error = if (length(errors) > 0L) unname(errors[[1L]]) else "Unknown emmeans error.",
      reason_skipped = "Time marginal means could not be estimated.", errors = errors
    ))
  }

  time_tidy <- safe(
    "time_summary",
    as.data.frame(
      summary(time_emm, infer = c(TRUE, FALSE), level = confidence_level),
      stringsAsFactors = FALSE
    )
  )
  time_pairwise <- safe(
    "time_pairwise",
    emmeans::contrast(time_emm, method = "pairwise", adjust = "none")
  )
  time_pairwise_tidy <- safe(
    "time_pairwise_summary",
    .mira_emm_summary(time_pairwise, primary_method, confidence_level = confidence_level)
  )
  time_pairwise_tukey <- safe(
    "time_pairwise_tukey",
    .mira_emm_summary(
      time_pairwise, primary_method, adjust = "tukey",
      confidence_level = confidence_level
    )
  )

  time_grid <- safe(
    "time_grid",
    as.data.frame(time_emm, stringsAsFactors = FALSE)
  )
  observed_time_levels <- if (!is.null(time_grid) && "time_factor" %in% names(time_grid)) {
    unique(as.character(time_grid$time_factor))
  } else character(0)
  n_time <- length(observed_time_levels)

  baseline_methods <- if (n_time >= 2L) {
    lapply(seq.int(2L, n_time), function(k) {
      coefficient <- numeric(n_time)
      coefficient[[1L]] <- -1
      coefficient[[k]] <- 1
      coefficient
    })
  } else list()
  if (length(baseline_methods) > 0L) {
    names(baseline_methods) <- paste0(
      observed_time_levels[-1L], " - ", observed_time_levels[[1L]]
    )
  }
  baseline_followup <- safe(
    "baseline_followup",
    if (length(baseline_methods) > 0L) {
      emmeans::contrast(time_emm, method = baseline_methods, adjust = "none")
    } else NULL
  )
  baseline_followup_tidy <- safe(
    "baseline_followup_summary",
    .mira_emm_summary(
      baseline_followup, primary_method, confidence_level = confidence_level
    )
  )
  baseline_followup_dunnett <- safe(
    "baseline_followup_dunnett",
    .mira_emm_summary(
      baseline_followup, primary_method, adjust = "dunnettx",
      confidence_level = confidence_level
    )
  )

  baseline_final_method <- NULL
  if (n_time >= 2L) {
    coefficient <- numeric(n_time)
    coefficient[[1L]] <- -1
    coefficient[[n_time]] <- 1
    baseline_final_method <- list("final - baseline" = coefficient)
  }
  baseline_final <- safe(
    "baseline_final",
    if (!is.null(baseline_final_method)) {
      emmeans::contrast(time_emm, method = baseline_final_method, adjust = "none")
    } else NULL
  )
  baseline_final_tidy <- safe(
    "baseline_final_summary",
    .mira_emm_summary(
      baseline_final, primary_method, confidence_level = confidence_level
    )
  )

  consecutive_methods <- if (n_time >= 2L) {
    lapply(seq.int(2L, n_time), function(k) {
      coefficient <- numeric(n_time)
      coefficient[[k - 1L]] <- -1
      coefficient[[k]] <- 1
      coefficient
    })
  } else list()
  if (length(consecutive_methods) > 0L) {
    names(consecutive_methods) <- paste0(
      observed_time_levels[-1L], " - ", observed_time_levels[-n_time]
    )
  }
  consecutive <- safe(
    "consecutive",
    if (length(consecutive_methods) > 0L) {
      emmeans::contrast(time_emm, method = consecutive_methods, adjust = "none")
    } else NULL
  )
  consecutive_tidy <- safe(
    "consecutive_summary",
    .mira_emm_summary(
      consecutive, primary_method, confidence_level = confidence_level
    )
  )

  trends <- safe(
    "ordinal_polynomial_trends",
    if (n_time >= 3L) {
      emmeans::contrast(time_emm, method = "poly", adjust = "none")
    } else NULL
  )
  if (!is.null(trends)) {
    trend_grid <- safe(
      "ordinal_polynomial_trend_grid",
      as.data.frame(trends, stringsAsFactors = FALSE)
    )
    if (!is.null(trend_grid) && "contrast" %in% names(trend_grid)) {
      allowed <- c("linear", "quadratic", if (n_time >= 4L) "cubic")
      keep <- tolower(as.character(trend_grid$contrast)) %in% allowed
      subsetted_trends <- safe("ordinal_polynomial_trend_subset", trends[keep])
      if (!is.null(subsetted_trends)) trends <- subsetted_trends
    }
  }
  trends_tidy <- safe(
    "ordinal_polynomial_trends_summary",
    .mira_emm_summary(trends, primary_method, confidence_level = confidence_level)
  )
  if (!is.null(trends_tidy) && "contrast" %in% names(trends_tidy)) {
    allowed <- c("linear", "quadratic", if (n_time >= 4L) "cubic")
    trends_tidy <- trends_tidy[tolower(as.character(trends_tidy$contrast)) %in% allowed, , drop = FALSE]
  }

  arm_emm <- arm_tidy <- arm_pairwise <- arm_pairwise_tidy <- NULL
  arm_pairwise_tukey <- arm_time_emm <- arm_time_tidy <- NULL
  simple_arm <- simple_arm_tidy <- simple_time <- simple_time_tidy <- NULL
  interaction <- interaction_tidy <- NULL

  if (isTRUE(arm_enabled)) {
    arm_emm <- safe("arm", emmeans::emmeans(model, ~ arm_factor))
    arm_tidy <- safe(
      "arm_summary",
      if (!is.null(arm_emm)) {
        as.data.frame(
          summary(arm_emm, infer = c(TRUE, FALSE), level = confidence_level),
          stringsAsFactors = FALSE
        )
      } else NULL
    )
    arm_pairwise <- safe(
      "arm_pairwise",
      if (!is.null(arm_emm)) {
        emmeans::contrast(arm_emm, method = "pairwise", adjust = "none")
      } else NULL
    )
    arm_pairwise_tidy <- safe(
      "arm_pairwise_summary",
      .mira_emm_summary(
        arm_pairwise, primary_method, confidence_level = confidence_level
      )
    )
    arm_pairwise_tukey <- safe(
      "arm_pairwise_tukey",
      .mira_emm_summary(
        arm_pairwise, primary_method, adjust = "tukey",
        confidence_level = confidence_level
      )
    )

    arm_time_emm <- safe(
      "arm_time",
      emmeans::emmeans(model, ~ arm_factor * time_factor)
    )
    arm_time_tidy <- safe(
      "arm_time_summary",
      if (!is.null(arm_time_emm)) {
        as.data.frame(
          summary(arm_time_emm, infer = c(TRUE, FALSE), level = confidence_level),
          stringsAsFactors = FALSE
        )
      } else NULL
    )

    arm_by_time <- safe(
      "arm_by_time",
      emmeans::emmeans(model, ~ arm_factor | time_factor)
    )
    simple_arm <- safe(
      "simple_arm",
      if (!is.null(arm_by_time)) {
        emmeans::contrast(arm_by_time, method = "pairwise", adjust = "none")
      } else NULL
    )
    simple_arm_tidy <- safe(
      "simple_arm_summary",
      .mira_emm_summary(simple_arm, primary_method, confidence_level = confidence_level)
    )

    time_by_arm <- safe(
      "time_by_arm",
      emmeans::emmeans(model, ~ time_factor | arm_factor)
    )
    simple_time <- safe(
      "simple_time",
      if (!is.null(time_by_arm)) {
        emmeans::contrast(time_by_arm, method = "pairwise", adjust = "none")
      } else NULL
    )
    simple_time_tidy <- safe(
      "simple_time_summary",
      .mira_emm_summary(simple_time, primary_method, confidence_level = confidence_level)
    )

    interaction <- safe(
      "interaction_contrasts",
      if (!is.null(arm_time_emm)) {
        emmeans::contrast(
          arm_time_emm,
          interaction = c(arm_factor = "pairwise", time_factor = "pairwise"),
          by = NULL,
          adjust = "none"
        )
      } else NULL
    )
    interaction_tidy <- safe(
      "interaction_contrasts_summary",
      .mira_emm_summary(interaction, primary_method, confidence_level = confidence_level)
    )
  }

  list(
    performed = TRUE,
    package = "emmeans",
    object = time_emm,
    summary = time_tidy,
    tidy = time_tidy,
    time = list(
      object = time_emm,
      tidy = time_tidy,
      pairwise = time_pairwise,
      pairwise_tidy = time_pairwise_tidy,
      pairwise_tukey = time_pairwise_tukey
    ),
    arm = list(
      object = arm_emm,
      tidy = arm_tidy,
      pairwise = arm_pairwise,
      pairwise_tidy = arm_pairwise_tidy,
      pairwise_tukey = arm_pairwise_tukey
    ),
    arm_time = list(
      object = arm_time_emm,
      tidy = arm_time_tidy,
      simple_arm = simple_arm,
      simple_arm_tidy = simple_arm_tidy,
      simple_time = simple_time,
      simple_time_tidy = simple_time_tidy,
      interaction = interaction,
      interaction_tidy = interaction_tidy
    ),
    baseline_followup = list(
      object = baseline_followup,
      tidy = baseline_followup_tidy,
      dunnett = baseline_followup_dunnett,
      dunnett_method = "emmeans dunnettx approximation"
    ),
    baseline_final = list(object = baseline_final, tidy = baseline_final_tidy),
    consecutive = list(object = consecutive, tidy = consecutive_tidy),
    trends = list(
      object = trends,
      tidy = trends_tidy,
      scale = "ordered timepoint levels",
      note = paste0(
        "Polynomial contrasts use the order of timepoints; they do not assume ",
        "equally spaced chronological time."
      )
    ),
    warnings = warnings,
    error = NULL,
    errors = errors,
    reason_skipped = NULL
  )
}

.mira_standardize_afex <- function(table) {
  if (is.null(table)) return(NULL)
  out <- as.data.frame(table, stringsAsFactors = FALSE)
  effect_col <- .mira_column_name(out, candidates = c("Effect", "term"))
  effect <- if (!is.na(effect_col)) as.character(out[[effect_col]]) else rownames(out)
  if (is.null(effect) || length(effect) != nrow(out)) effect <- rep(NA_character_, nrow(out))
  stat_col <- .mira_column_name(out, candidates = c("F", "F.value", "F value"))
  p_col <- .mira_column_name(out, candidates = c("Pr(>F)", "p", "p.value"))
  num_df_col <- .mira_column_name(out, candidates = c("num Df", "num_df", "Df"))
  den_df_col <- .mira_column_name(out, candidates = c("den Df", "den_df"))
  pes_col <- .mira_column_name(out, candidates = c("pes", "partial eta squared"), contains = "pes")
  ges_col <- .mira_column_name(out, candidates = c("ges", "generalized eta squared"), contains = "ges")
  value <- function(column) {
    if (is.na(column)) rep(NA_real_, nrow(out)) else .mira_as_numeric(out[[column]])
  }
  data.frame(
    effect = effect,
    statistic = value(stat_col),
    df_num = value(num_df_col),
    df_den = value(den_df_col),
    p_value = value(p_col),
    partial_eta_squared = value(pes_col),
    generalized_eta_squared = value(ges_col),
    stringsAsFactors = FALSE
  )
}

.mira_run_rm_anova <- function(data,
                               arm_enabled,
                               covariates,
                               categorical_covariates) {
  if (!requireNamespace("afex", quietly = TRUE)) {
    return(.mira_skipped_module(
      "afex",
      "Optional package 'afex' is not installed."
    ))
  }

  between <- c(
    if (isTRUE(arm_enabled)) "arm_factor" else character(0),
    intersect(covariates, categorical_covariates)
  )
  numeric_covariates <- setdiff(covariates, categorical_covariates)
  numeric_covariates <- numeric_covariates[vapply(
    numeric_covariates,
    function(variable) is.numeric(data[[variable]]),
    logical(1L)
  )]
  observed <- unique(c(intersect(covariates, categorical_covariates), numeric_covariates))

  fitted <- .mira_eval(
    afex::aov_ez(
      id = "patient_factor",
      dv = "value",
      data = data,
      between = if (length(between) > 0L) between else NULL,
      within = "time_factor",
      covariate = if (length(numeric_covariates) > 0L) numeric_covariates else NULL,
      observed = if (length(observed) > 0L) observed else NULL,
      type = 3,
      factorize = FALSE,
      include_aov = TRUE,
      anova_table = list(correction = "none", es = "pes")
    )
  )
  if (is.null(fitted$value)) {
    return(list(
      performed = FALSE, package = "afex", object = NULL, summary = NULL,
      tidy = NULL, warnings = fitted$warnings, error = fitted$error,
      reason_skipped = "Repeated-measures ANOVA could not be estimated."
    ))
  }

  object <- fitted$value
  raw_pes_eval <- .mira_eval(stats::anova(object, correction = "none", es = "pes"))
  raw_ges_eval <- .mira_eval(stats::anova(object, correction = "none", es = "ges"))
  gg_eval <- .mira_eval(stats::anova(object, correction = "GG", es = "pes"))
  hf_eval <- .mira_eval(stats::anova(object, correction = "HF", es = "pes"))
  summary_eval <- .mira_eval(summary(object))

  raw <- .mira_standardize_afex(raw_pes_eval$value)
  ges <- .mira_standardize_afex(raw_ges_eval$value)
  gg <- .mira_standardize_afex(gg_eval$value)
  hf <- .mira_standardize_afex(hf_eval$value)
  if (!is.null(raw) && !is.null(ges)) {
    raw$generalized_eta_squared <- ges$generalized_eta_squared[
      match(raw$effect, ges$effect)
    ]
  }
  tidy <- raw
  if (!is.null(tidy) && !is.null(gg)) {
    tidy$df_num_gg <- gg$df_num[match(tidy$effect, gg$effect)]
    tidy$df_den_gg <- gg$df_den[match(tidy$effect, gg$effect)]
    tidy$p_gg <- gg$p_value[match(tidy$effect, gg$effect)]
  }
  if (!is.null(tidy) && !is.null(hf)) {
    tidy$df_num_hf <- hf$df_num[match(tidy$effect, hf$effect)]
    tidy$df_den_hf <- hf$df_den[match(tidy$effect, hf$effect)]
    tidy$p_hf <- hf$p_value[match(tidy$effect, hf$effect)]
  }

  full_summary <- summary_eval$value
  warnings <- unique(c(
    fitted$warnings, raw_pes_eval$warnings, raw_ges_eval$warnings,
    gg_eval$warnings, hf_eval$warnings, summary_eval$warnings
  ))
  errors <- Filter(Negate(is.null), list(
    raw_pes = raw_pes_eval$error,
    raw_ges = raw_ges_eval$error,
    green_house_geisser = gg_eval$error,
    huynh_feldt = hf_eval$error,
    summary = summary_eval$error
  ))

  rm_subjects <- tryCatch(
    as.integer(stats::nobs(object$lm)),
    error = function(e) {
      observed_by_subject <- tapply(
        as.character(data$time_factor),
        data$patient_factor,
        function(x) length(unique(x))
      )
      as.integer(sum(observed_by_subject == nlevels(data$time_factor)))
    }
  )

  list(
    performed = TRUE,
    package = "afex",
    object = object,
    summary = full_summary,
    tidy = tidy,
    tables = list(
      uncorrected_partial_eta = raw_pes_eval$value,
      uncorrected_generalized_eta = raw_ges_eval$value,
      greenhouse_geisser = gg_eval$value,
      huynh_feldt = hf_eval$value
    ),
    sphericity = list(
      mauchly = tryCatch(full_summary[["sphericity.tests"]], error = function(e) NULL),
      corrections = tryCatch(full_summary[["pval.adjustments"]], error = function(e) NULL)
    ),
    n = rm_subjects,
    n_subjects = rm_subjects,
    covariates = covariates,
    warnings = warnings,
    error = NULL,
    errors = errors,
    reason_skipped = NULL
  )
}

.mira_rank_biserial_paired <- function(delta) {
  delta <- delta[is.finite(delta) & delta != 0]
  if (length(delta) == 0L) return(NA_real_)
  ranks <- rank(abs(delta), ties.method = "average")
  positive <- sum(ranks[delta > 0])
  negative <- sum(ranks[delta < 0])
  denominator <- positive + negative
  if (!is.finite(denominator) || denominator == 0) NA_real_ else
    (positive - negative) / denominator
}

.mira_run_friedman <- function(analysis_data, time_vars, time_labels, primary_method) {
  complete <- analysis_data[
    stats::complete.cases(analysis_data[time_vars]),
    time_vars,
    drop = FALSE
  ]
  global_eval <- .mira_eval(stats::friedman.test(as.matrix(complete)))
  test <- global_eval$value
  n_complete <- nrow(complete)
  k <- length(time_vars)
  kendall_w <- if (!is.null(test) && n_complete > 0L && k > 1L) {
    as.numeric(test$statistic) / (n_complete * (k - 1L))
  } else NA_real_

  tests <- list()
  rows <- list()
  warnings <- global_eval$warnings
  errors <- list(global = global_eval$error)
  counter <- 1L
  for (i in seq_len(length(time_vars) - 1L)) {
    for (j in seq.int(i + 1L, length(time_vars))) {
      from <- time_vars[[i]]
      to <- time_vars[[j]]
      keep <- stats::complete.cases(analysis_data[[from]], analysis_data[[to]])
      x <- analysis_data[[from]][keep]
      y <- analysis_data[[to]][keep]
      evaluated <- .mira_eval(
        stats::wilcox.test(y, x, paired = TRUE, exact = FALSE)
      )
      key <- paste(from, to, sep = "_to_")
      tests[[key]] <- evaluated$value
      warnings <- unique(c(warnings, evaluated$warnings))
      errors[[key]] <- evaluated$error
      rows[[counter]] <- data.frame(
        from = from,
        to = to,
        from_label = unname(time_labels[from]),
        to_label = unname(time_labels[to]),
        n = sum(keep),
        statistic = if (!is.null(evaluated$value)) {
          as.numeric(evaluated$value$statistic)
        } else NA_real_,
        p_raw = if (!is.null(evaluated$value)) evaluated$value$p.value else NA_real_,
        rank_biserial = .mira_rank_biserial_paired(y - x),
        stringsAsFactors = FALSE
      )
      counter <- counter + 1L
    }
  }
  posthoc <- if (length(rows) > 0L) do.call(rbind, rows) else data.frame()
  posthoc <- .mira_adjust_family(posthoc, p_col = "p_raw", primary = primary_method)

  tidy <- data.frame(
    statistic = if (!is.null(test)) as.numeric(test$statistic) else NA_real_,
    df = if (!is.null(test)) as.numeric(test$parameter) else NA_real_,
    p_raw = if (!is.null(test)) test$p.value else NA_real_,
    kendalls_w = kendall_w,
    n = n_complete,
    n_timepoints = k,
    stringsAsFactors = FALSE
  )

  list(
    performed = !is.null(test),
    package = "stats",
    object = test,
    test = test,
    summary = tidy,
    tidy = tidy,
    posthoc = list(tests = tests, tidy = posthoc),
    warnings = warnings,
    error = global_eval$error,
    errors = Filter(Negate(is.null), errors),
    reason_skipped = if (is.null(test)) {
      "Friedman test could not be computed on complete repeated profiles."
    } else NULL
  )
}

.mira_effect_indices <- function(coefficient_names, question) {
  has_time <- grepl("time_factor", coefficient_names, fixed = TRUE)
  has_arm <- grepl("arm_factor", coefficient_names, fixed = TRUE)
  if (identical(question, "TIME_X_ARM")) {
    which(has_time & has_arm)
  } else if (identical(question, "TIME")) {
    which(has_time)
  } else if (identical(question, "ARM")) {
    which(has_arm)
  } else integer(0)
}

.mira_wald_chisq <- function(beta, covariance, indices) {
  if (length(indices) == 0L) return(NULL)
  beta <- as.numeric(beta)
  covariance <- as.matrix(covariance)
  constraint <- diag(length(beta))[indices, , drop = FALSE]
  estimate <- as.numeric(constraint %*% beta)
  variance <- constraint %*% covariance %*% t(constraint)
  decomposition <- tryCatch(
    eigen((variance + t(variance)) / 2, symmetric = TRUE),
    error = function(e) NULL
  )
  if (is.null(decomposition)) return(NULL)
  tolerance <- max(abs(decomposition$values), 1) * sqrt(.Machine$double.eps)
  keep <- decomposition$values > tolerance
  rank <- sum(keep)
  if (rank == 0L) return(NULL)
  rotated <- crossprod(decomposition$vectors[, keep, drop = FALSE], estimate)
  statistic <- sum(as.numeric(rotated)^2 / decomposition$values[keep])
  list(
    statistic = statistic,
    df = rank,
    p_value = stats::pchisq(statistic, df = rank, lower.tail = FALSE),
    constraints = constraint
  )
}

.mira_gee_coefficients <- function(model, confidence_level) {
  model_summary <- summary(model)
  table <- as.data.frame(model_summary$coefficients, stringsAsFactors = FALSE)
  table$term <- rownames(table)
  rownames(table) <- NULL
  table <- table[c("term", setdiff(names(table), "term"))]
  estimate_col <- .mira_column_name(table, candidates = "Estimate")
  se_col <- .mira_column_name(
    table,
    candidates = c("Std.err", "Std. Error", "San.se"),
    contains = c("stderr", "sanse")
  )
  estimate <- if (!is.na(estimate_col)) .mira_as_numeric(table[[estimate_col]]) else
    rep(NA_real_, nrow(table))
  robust_se <- if (!is.na(se_col)) .mira_as_numeric(table[[se_col]]) else
    rep(NA_real_, nrow(table))
  critical <- stats::qnorm(1 - (1 - confidence_level) / 2)
  table$robust_se <- robust_se
  table$conf_low <- estimate - critical * robust_se
  table$conf_high <- estimate + critical * robust_se
  table
}

.mira_run_gee <- function(data, fixed_formula, confidence_level) {
  if (!requireNamespace("geepack", quietly = TRUE)) {
    skipped <- .mira_skipped_module(
      "geepack",
      "Optional package 'geepack' is not installed."
    )
    return(list(
      performed = FALSE,
      package = "geepack",
      independence = skipped,
      exchangeable = skipped,
      ar1 = skipped,
      warnings = character(0),
      error = NULL,
      reason_skipped = skipped$reason_skipped
    ))
  }

  structures <- c(independence = "independence", exchangeable = "exchangeable", ar1 = "ar1")
  output <- list()
  all_warnings <- character(0)

  for (name in names(structures)) {
    correlation <- structures[[name]]
    fitted <- .mira_eval(
      geepack::geeglm(
        formula = fixed_formula,
        family = stats::gaussian(),
        data = data,
        id = data$patient_factor,
        waves = data$time_index,
        corstr = correlation,
        std.err = "san.se"
      )
    )
    all_warnings <- unique(c(all_warnings, fitted$warnings))
    if (is.null(fitted$value)) {
      output[[name]] <- list(
        performed = FALSE, package = "geepack", model = NULL, object = NULL,
        summary = NULL, coefficients = NULL, effect_tests = NULL, qic = NULL,
        correlation_structure = correlation, formula = fixed_formula,
        converged = FALSE, warnings = fitted$warnings, error = fitted$error,
        reason_skipped = "GEE estimation failed for this working correlation."
      )
      next
    }

    model <- fitted$value
    summary_eval <- .mira_eval(summary(model))
    coefficient_eval <- .mira_eval(.mira_gee_coefficients(model, confidence_level))
    qic_eval <- .mira_eval(geepack::QIC(model))
    anova_eval <- .mira_eval(stats::anova(model))
    beta_eval <- .mira_eval(stats::coef(model))
    vcov_eval <- .mira_eval(stats::vcov(model))
    effect_rows <- if (!is.null(beta_eval$value) && !is.null(vcov_eval$value)) {
      lapply(c("TIME", "ARM", "TIME_X_ARM"), function(question) {
        test <- .mira_wald_chisq(
          beta_eval$value,
          vcov_eval$value,
          .mira_effect_indices(names(beta_eval$value), question)
        )
        if (is.null(test)) return(NULL)
        data.frame(
          question = question,
          statistic = test$statistic,
          df = test$df,
          p_raw = test$p_value,
          stringsAsFactors = FALSE
        )
      })
    } else list()
    effect_rows <- Filter(Negate(is.null), effect_rows)
    effect_tests <- if (length(effect_rows) > 0L) do.call(rbind, effect_rows) else
      data.frame(question = character(0), statistic = numeric(0), df = numeric(0),
                 p_raw = numeric(0), stringsAsFactors = FALSE)
    warnings <- unique(c(
      fitted$warnings, summary_eval$warnings, coefficient_eval$warnings,
      qic_eval$warnings, anova_eval$warnings, beta_eval$warnings, vcov_eval$warnings
    ))
    all_warnings <- unique(c(all_warnings, warnings))
    output[[name]] <- list(
      performed = TRUE,
      package = "geepack",
      model = model,
      object = model,
      summary = summary_eval$value,
      coefficients = coefficient_eval$value,
      robust_standard_errors = if (!is.null(coefficient_eval$value)) {
        coefficient_eval$value[c("term", "robust_se", "conf_low", "conf_high")]
      } else NULL,
      wald = anova_eval$value,
      effect_tests = effect_tests,
      qic = qic_eval$value,
      correlation_structure = correlation,
      waves = "time_index (ordered timepoint index)",
      formula = fixed_formula,
      converged = tryCatch(isTRUE(model$geese$error == 0), error = function(e) NA),
      warnings = warnings,
      error = NULL,
      errors = Filter(Negate(is.null), list(
        summary = summary_eval$error,
        coefficients = coefficient_eval$error,
        qic = qic_eval$error,
        anova = anova_eval$error,
        coefficient_vector = beta_eval$error,
        covariance = vcov_eval$error
      )),
      reason_skipped = NULL
    )
  }

  list(
    performed = any(vapply(output, function(x) isTRUE(x$performed), logical(1L))),
    package = "geepack",
    independence = output$independence,
    exchangeable = output$exchangeable,
    ar1 = output$ar1,
    warnings = all_warnings,
    error = NULL,
    reason_skipped = NULL
  )
}

.mira_lmer_diagnostics <- function(model) {
  if (is.null(model)) {
    return(list(converged = FALSE, singular = NA, variance_components = NULL))
  }
  optinfo <- tryCatch(model@optinfo, error = function(e) NULL)
  messages <- if (!is.null(optinfo)) optinfo$conv$lme4$messages else NULL
  optimizer_ok <- if (!is.null(optinfo) && !is.null(optinfo$conv$opt)) {
    isTRUE(optinfo$conv$opt == 0)
  } else TRUE
  list(
    converged = optimizer_ok && is.null(messages),
    singular = tryCatch(lme4::isSingular(model, tol = 1e-4), error = function(e) NA),
    variance_components = tryCatch(
      as.data.frame(lme4::VarCorr(model)),
      error = function(e) NULL
    )
  )
}

.mira_run_random_slope <- function(data,
                                   arm_enabled,
                                   covariates,
                                   mixed_intercept_model) {
  if (!requireNamespace("lme4", quietly = TRUE)) {
    return(.mira_skipped_module(
      "lme4",
      "Optional package 'lme4' is not installed."
    ))
  }

  random_term <- "(1 + time_ordinal | patient_factor)"
  full_core <- if (isTRUE(arm_enabled)) "time_factor * arm_factor" else "time_factor"
  formula <- .mira_make_lmer_formula(full_core, covariates, random_term)
  fitted <- .mira_eval(
    if (requireNamespace("lmerTest", quietly = TRUE)) {
      lmerTest::lmer(
        formula, data = data, REML = TRUE, na.action = stats::na.omit
      )
    } else {
      lme4::lmer(
        formula, data = data, REML = TRUE, na.action = stats::na.omit
      )
    }
  )
  if (is.null(fitted$value)) {
    return(list(
      performed = FALSE, package = "lme4", model = NULL, object = NULL,
      summary = NULL, tidy = NULL, formula = formula, warnings = fitted$warnings,
      error = fitted$error,
      reason_skipped = "The random-slope mixed model could not be estimated."
    ))
  }

  model <- fitted$value
  diagnostics <- .mira_lmer_diagnostics(model)
  summary_eval <- .mira_eval(summary(model))
  anova_eval <- .mira_eval(stats::anova(model))
  fixed_eval <- .mira_eval(as.data.frame(coef(summary(model)), stringsAsFactors = FALSE))
  if (!is.null(fixed_eval$value)) {
    fixed_eval$value$term <- rownames(fixed_eval$value)
    rownames(fixed_eval$value) <- NULL
    fixed_eval$value <- fixed_eval$value[c("term", setdiff(names(fixed_eval$value), "term"))]
  }
  vc <- diagnostics$variance_components
  random_correlation <- NA_real_
  if (!is.null(vc)) {
    correlation_rows <- vc$grp == "patient_factor" & !is.na(vc$var2)
    if (any(correlation_rows)) random_correlation <- vc$sdcor[which(correlation_rows)[[1L]]]
  }

  warnings <- unique(c(
    fitted$warnings, summary_eval$warnings, anova_eval$warnings, fixed_eval$warnings
  ))
  errors <- Filter(Negate(is.null), list(
    summary = summary_eval$error,
    anova = anova_eval$error,
    fixed_effects = fixed_eval$error
  ))

  comparison_eval <- .mira_eval(
    if (!is.null(mixed_intercept_model)) {
      stats::anova(mixed_intercept_model, model, refit = FALSE)
    } else NULL
  )
  warnings <- unique(c(warnings, comparison_eval$warnings))
  if (!is.null(comparison_eval$error)) errors$random_effects_comparison <- comparison_eval$error

  ml_eval <- .mira_eval(
    lme4::lmer(formula, data = data, REML = FALSE, na.action = stats::na.omit)
  )
  warnings <- unique(c(warnings, ml_eval$warnings))
  if (!is.null(ml_eval$error)) errors$full_ml <- ml_eval$error

  fit_reduced <- function(name, core) {
    if (is.null(ml_eval$value)) return(NULL)
    reduced_formula <- .mira_make_lmer_formula(core, covariates, random_term)
    reduced_eval <- .mira_eval(
      lme4::lmer(
        reduced_formula, data = data, REML = FALSE, na.action = stats::na.omit
      )
    )
    warnings <<- unique(c(warnings, reduced_eval$warnings))
    if (!is.null(reduced_eval$error)) {
      errors[[paste0(name, "_model")]] <<- reduced_eval$error
      return(NULL)
    }
    test_eval <- .mira_eval(stats::anova(reduced_eval$value, ml_eval$value))
    warnings <<- unique(c(warnings, test_eval$warnings))
    if (!is.null(test_eval$error)) errors[[name]] <<- test_eval$error
    test_eval$value
  }

  global_time <- global_arm <- interaction <- NULL
  if (isTRUE(arm_enabled)) {
    interaction <- fit_reduced("interaction", "time_factor + arm_factor")
    global_time <- fit_reduced("global_time", "arm_factor")
    global_arm <- fit_reduced("global_arm", "time_factor")
  } else {
    global_time <- fit_reduced("global_time", "1")
  }

  list(
    performed = TRUE,
    package = if (inherits(model, "lmerModLmerTest")) "lmerTest/lme4" else "lme4",
    model = model,
    object = model,
    summary = summary_eval$value,
    tidy = fixed_eval$value,
    fixed_effects = fixed_eval$value,
    fixed_effect_tests = anova_eval$value,
    variance_components = vc,
    random_intercept_slope_correlation = random_correlation,
    formula = formula,
    time_scale = "ordered timepoint index",
    note = paste0(
      "The random slope uses the ordered timepoint index and does not imply ",
      "equally spaced chronological time."
    ),
    AIC = tryCatch(stats::AIC(model), error = function(e) NA_real_),
    BIC = tryCatch(stats::BIC(model), error = function(e) NA_real_),
    logLik = tryCatch(as.numeric(stats::logLik(model)), error = function(e) NA_real_),
    converged = diagnostics$converged,
    singular = diagnostics$singular,
    random_effects_lrt = comparison_eval$value,
    global_time_test = global_time,
    global_arm_test = global_arm,
    arm_time_interaction_test = interaction,
    warnings = warnings,
    error = NULL,
    errors = errors,
    reason_skipped = NULL
  )
}

.mira_run_nlme <- function(data, fixed_formula) {
  if (!requireNamespace("nlme", quietly = TRUE)) {
    skipped <- .mira_skipped_module(
      "nlme",
      "Optional package 'nlme' is not installed."
    )
    return(list(
      performed = FALSE, package = "nlme", compound_symmetry = skipped,
      ar1 = skipped, comparison = NULL, warnings = character(0), error = NULL,
      reason_skipped = skipped$reason_skipped
    ))
  }

  structures <- list(
    compound_symmetry = nlme::corCompSymm(form = ~1 | patient_factor),
    ar1 = nlme::corAR1(form = ~time_index | patient_factor)
  )
  output <- list()
  all_warnings <- character(0)
  for (name in names(structures)) {
    fitted <- .mira_eval(
      nlme::lme(
        fixed = fixed_formula,
        random = ~1 | patient_factor,
        correlation = structures[[name]],
        data = data,
        method = "REML",
        na.action = stats::na.omit,
        control = nlme::lmeControl(returnObject = TRUE)
      )
    )
    all_warnings <- unique(c(all_warnings, fitted$warnings))
    if (is.null(fitted$value)) {
      output[[name]] <- list(
        performed = FALSE, package = "nlme", model = NULL, object = NULL,
        summary = NULL, coefficients = NULL, formula = fixed_formula,
        correlation_structure = if (name == "compound_symmetry") {
          "compound symmetry"
        } else "AR(1)",
        converged = FALSE, warnings = fitted$warnings, error = fitted$error,
        reason_skipped = "The nlme correlation model could not be estimated."
      )
      next
    }
    model <- fitted$value
    summary_eval <- .mira_eval(summary(model))
    coefficients <- if (!is.null(summary_eval$value)) summary_eval$value$tTable else NULL
    correlation_parameter <- tryCatch(
      as.numeric(coef(model$modelStruct$corStruct, unconstrained = FALSE)),
      error = function(e) NA_real_
    )
    warnings <- unique(c(fitted$warnings, summary_eval$warnings))
    all_warnings <- unique(c(all_warnings, warnings))
    output[[name]] <- list(
      performed = TRUE,
      package = "nlme",
      model = model,
      object = model,
      summary = summary_eval$value,
      coefficients = coefficients,
      formula = fixed_formula,
      correlation_structure = if (name == "compound_symmetry") {
        "compound symmetry"
      } else "AR(1)",
      time_scale = if (name == "ar1") "ordered timepoint index" else NA_character_,
      correlation_parameter = correlation_parameter,
      AIC = tryCatch(stats::AIC(model), error = function(e) NA_real_),
      BIC = tryCatch(stats::BIC(model), error = function(e) NA_real_),
      logLik = tryCatch(as.numeric(stats::logLik(model)), error = function(e) NA_real_),
      converged = !any(grepl(
        "convergence|iteration limit|did not converge",
        warnings,
        ignore.case = TRUE
      )),
      warnings = warnings,
      error = NULL,
      errors = Filter(Negate(is.null), list(summary = summary_eval$error)),
      reason_skipped = NULL
    )
  }

  comparison_rows <- lapply(names(output), function(name) {
    item <- output[[name]]
    data.frame(
      model = name,
      correlation_structure = item$correlation_structure,
      AIC = if (isTRUE(item$performed)) item$AIC else NA_real_,
      BIC = if (isTRUE(item$performed)) item$BIC else NA_real_,
      logLik = if (isTRUE(item$performed)) item$logLik else NA_real_,
      converged = if (isTRUE(item$performed)) item$converged else FALSE,
      stringsAsFactors = FALSE
    )
  })

  list(
    performed = any(vapply(output, function(x) isTRUE(x$performed), logical(1L))),
    package = "nlme",
    compound_symmetry = output$compound_symmetry,
    ar1 = output$ar1,
    comparison = do.call(rbind, comparison_rows),
    warnings = all_warnings,
    error = NULL,
    reason_skipped = NULL
  )
}

.mira_run_robust_inference <- function(model, confidence_level) {
  if (is.null(model)) {
    return(.mira_skipped_module(
      "clubSandwich",
      "The random-intercept mixed model was not available."
    ))
  }
  if (!requireNamespace("clubSandwich", quietly = TRUE)) {
    return(.mira_skipped_module(
      "clubSandwich",
      "Optional package 'clubSandwich' is not installed."
    ))
  }
  if (!requireNamespace("lme4", quietly = TRUE)) {
    return(.mira_skipped_module(
      "clubSandwich",
      "Package 'lme4' is required for robust mixed-model inference."
    ))
  }

  covariance_eval <- .mira_eval(clubSandwich::vcovCR(model, type = "CR2"))
  if (is.null(covariance_eval$value)) {
    return(list(
      performed = FALSE, package = "clubSandwich", object = NULL,
      summary = NULL, tidy = NULL, warnings = covariance_eval$warnings,
      error = covariance_eval$error,
      reason_skipped = "CR2 covariance estimation failed."
    ))
  }
  covariance <- covariance_eval$value
  coefficient_eval <- .mira_eval(
    clubSandwich::coef_test(model, vcov = covariance, test = "Satterthwaite")
  )
  confidence_eval <- .mira_eval(
    clubSandwich::conf_int(
      model,
      vcov = covariance,
      level = confidence_level,
      test = "Satterthwaite"
    )
  )

  coefficients <- lme4::fixef(model)
  effect_tests <- list()
  warnings <- unique(c(
    covariance_eval$warnings, coefficient_eval$warnings, confidence_eval$warnings
  ))
  errors <- Filter(Negate(is.null), list(
    coefficients = coefficient_eval$error,
    confidence_intervals = confidence_eval$error
  ))
  for (question in c("TIME", "ARM", "TIME_X_ARM")) {
    indices <- .mira_effect_indices(names(coefficients), question)
    if (length(indices) == 0L) next
    constraint <- diag(length(coefficients))[indices, , drop = FALSE]
    evaluated <- .mira_eval(
      clubSandwich::Wald_test(
        model,
        constraints = constraint,
        vcov = covariance,
        test = "HTZ",
        tidy = TRUE
      )
    )
    warnings <- unique(c(warnings, evaluated$warnings))
    if (!is.null(evaluated$error)) errors[[question]] <- evaluated$error
    effect_tests[[question]] <- evaluated$value
  }

  list(
    performed = TRUE,
    package = "clubSandwich",
    object = covariance,
    covariance = covariance,
    covariance_estimator = "CR2",
    summary = coefficient_eval$value,
    tidy = coefficient_eval$value,
    coefficient_tests = coefficient_eval$value,
    confidence_intervals = confidence_eval$value,
    effect_tests = effect_tests,
    test = "Satterthwaite coefficient tests; HTZ joint Wald tests",
    warnings = warnings,
    error = NULL,
    errors = errors,
    reason_skipped = NULL
  )
}

.mira_run_r2 <- function(model) {
  if (is.null(model)) {
    return(.mira_skipped_module(
      "performance",
      "The mixed model was not available."
    ))
  }
  if (!requireNamespace("performance", quietly = TRUE)) {
    return(.mira_skipped_module(
      "performance",
      "Optional package 'performance' is not installed."
    ))
  }
  evaluated <- .mira_eval(performance::r2_nakagawa(model))
  if (is.null(evaluated$value)) {
    return(list(
      performed = FALSE, package = "performance", object = NULL,
      summary = NULL, tidy = NULL, warnings = evaluated$warnings,
      error = evaluated$error, reason_skipped = "Nakagawa R-squared could not be computed."
    ))
  }
  values <- unlist(evaluated$value, use.names = TRUE)
  marginal_idx <- grep("marginal", names(values), ignore.case = TRUE)
  conditional_idx <- grep("conditional", names(values), ignore.case = TRUE)
  tidy <- data.frame(
    marginal_r2 = if (length(marginal_idx) > 0L) {
      as.numeric(values[marginal_idx[[1L]]])
    } else NA_real_,
    conditional_r2 = if (length(conditional_idx) > 0L) {
      as.numeric(values[conditional_idx[[1L]]])
    } else NA_real_,
    stringsAsFactors = FALSE
  )
  list(
    performed = TRUE, package = "performance", object = evaluated$value,
    summary = evaluated$value, tidy = tidy, warnings = evaluated$warnings,
    error = NULL, reason_skipped = NULL
  )
}

.mira_qic_value <- function(qic, name) {
  if (is.null(qic)) return(NA_real_)
  if (is.matrix(qic) || is.data.frame(qic)) {
    column <- match(tolower(name), tolower(colnames(qic)), nomatch = 0L)
    if (column > 0L) return(.mira_as_numeric(qic[1L, column])[[1L]])
  }
  values <- unlist(qic, use.names = TRUE)
  index <- match(tolower(name), tolower(names(values)), nomatch = 0L)
  if (index > 0L) .mira_as_numeric(values[[index]])[[1L]] else NA_real_
}

.mira_build_model_comparison <- function(data,
                                         fixed_formula,
                                         mixed_model,
                                         mixed_converged,
                                         mixed_singular,
                                         random_slope,
                                         nlme_models,
                                         gee_models) {
  default_n <- nrow(data)
  default_subjects <- length(unique(data$patient_factor))
  safe_nobs <- function(model) {
    if (is.null(model)) return(default_n)
    value <- tryCatch(as.integer(stats::nobs(model)), error = function(e) NA_integer_)
    if (length(value) == 0L || is.na(value[[1L]])) default_n else value[[1L]]
  }
  likelihood <- function(model, fun) {
    if (is.null(model)) return(NA_real_)
    tryCatch(as.numeric(fun(model)), error = function(e) NA_real_)
  }
  row <- function(model_name,
                  model,
                  formula,
                  converged,
                  singular,
                  correlation_structure,
                  object_path,
                  likelihood_based = TRUE,
                  qic = NA_real_,
                  cic = NA_real_) {
    data.frame(
      model = model_name,
      formula = .mira_formula_text(formula),
      n = safe_nobs(model),
      n_subjects = default_subjects,
      AIC = if (likelihood_based) likelihood(model, stats::AIC) else NA_real_,
      BIC = if (likelihood_based) likelihood(model, stats::BIC) else NA_real_,
      logLik = if (likelihood_based) likelihood(model, stats::logLik) else NA_real_,
      QIC = qic,
      CIC = cic,
      converged = if (length(converged) == 0L) NA else as.logical(converged)[[1L]],
      singular = if (length(singular) == 0L) NA else as.logical(singular)[[1L]],
      correlation_structure = correlation_structure,
      object_path = if (is.null(model)) NA_character_ else object_path,
      stringsAsFactors = FALSE
    )
  }

  ri_formula <- if (!is.null(mixed_model)) {
    stats::formula(mixed_model)
  } else {
    stats::as.formula(paste(
      .mira_formula_text(fixed_formula),
      "+ (1 | patient_factor)"
    ))
  }
  random_slope_formula <- if (!is.null(random_slope$formula)) {
    random_slope$formula
  } else {
    stats::as.formula(paste(
      .mira_formula_text(fixed_formula),
      "+ (1 + time_ordinal | patient_factor)"
    ))
  }
  rows <- list(
    row(
      "mixed_random_intercept", mixed_model, ri_formula, mixed_converged,
      mixed_singular, "random intercept", "$model$fitted_model"
    ),
    row(
      "mixed_random_slope", random_slope$model, random_slope_formula,
      random_slope$converged, random_slope$singular,
      "random intercept + ordinal-time slope",
      "$advanced_models$random_slope$model"
    ),
    row(
      "nlme_compound_symmetry", nlme_models$compound_symmetry$model,
      fixed_formula, nlme_models$compound_symmetry$converged, NA,
      "compound symmetry", "$advanced_models$nlme$compound_symmetry$model"
    ),
    row(
      "nlme_ar1", nlme_models$ar1$model, fixed_formula,
      nlme_models$ar1$converged, NA, "AR(1)",
      "$advanced_models$nlme$ar1$model"
    )
  )

  for (name in c("independence", "exchangeable", "ar1")) {
    item <- gee_models[[name]]
    rows[[length(rows) + 1L]] <- row(
      paste0("gee_", name), item$model, fixed_formula, item$converged, NA,
      name, paste0("$advanced_models$gee$", name, "$model"),
      likelihood_based = FALSE,
      qic = .mira_qic_value(item$qic, "QIC"),
      cic = .mira_qic_value(item$qic, "CIC")
    )
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  attr(out, "note") <- paste0(
    "Likelihood criteria are descriptive and do not identify a uniquely true model. ",
    "GEE rows use QIC/CIC and intentionally leave AIC, BIC, and logLik unavailable."
  )
  out
}

.mira_sensitivity_row <- function(question,
                                  method,
                                  statistic = NA_real_,
                                  df = NA_character_,
                                  p_raw = NA_real_,
                                  effect_size = NA_real_,
                                  effect_size_type = NA_character_,
                                  n = NA_integer_,
                                  object_path,
                                  note = NA_character_) {
  data.frame(
    question = question,
    method = method,
    statistic = as.numeric(statistic)[[1L]],
    df = as.character(df)[[1L]],
    p_raw = as.numeric(p_raw)[[1L]],
    effect_size = as.numeric(effect_size)[[1L]],
    effect_size_type = as.character(effect_size_type)[[1L]],
    n = as.integer(n)[[1L]],
    object_path = object_path,
    note = as.character(note)[[1L]],
    stringsAsFactors = FALSE
  )
}

.mira_lrt_sensitivity <- function(test, question, method, object_path, n, note = NA_character_) {
  if (is.null(test)) return(NULL)
  table <- tryCatch(as.data.frame(test), error = function(e) NULL)
  if (is.null(table) || nrow(table) == 0L) return(NULL)
  p_col <- .mira_column_name(
    table,
    candidates = c("Pr(>Chisq)", "p.value", "p"),
    contains = c("prchisq", "pvalue")
  )
  statistic_col <- .mira_column_name(
    table,
    candidates = c("Chisq", "L.Ratio", "F", "Deviance"),
    contains = c("chisq", "lratio")
  )
  df_col <- .mira_column_name(
    table,
    candidates = c("Chi Df", "Df", "numDF"),
    contains = c("chidf", "numdf")
  )
  p <- if (!is.na(p_col)) .mira_as_numeric(table[[p_col]]) else rep(NA_real_, nrow(table))
  candidate <- which(is.finite(p))
  index <- if (length(candidate) > 0L) tail(candidate, 1L) else nrow(table)
  .mira_sensitivity_row(
    question = question,
    method = method,
    statistic = if (!is.na(statistic_col)) .mira_as_numeric(table[[statistic_col]])[[index]] else NA_real_,
    df = if (!is.na(df_col)) .mira_as_numeric(table[[df_col]])[[index]] else NA_real_,
    p_raw = p[[index]],
    n = n,
    object_path = object_path,
    note = note
  )
}

.mira_match_effect <- function(table, question) {
  if (is.null(table) || nrow(table) == 0L || !"effect" %in% names(table)) return(integer(0))
  effect <- tolower(gsub("[^[:alnum:]]+", "", table$effect))
  has_time <- grepl("timefactor|time", effect)
  has_arm <- grepl("armfactor|arm", effect)
  if (question == "TIME_X_ARM") which(has_time & has_arm) else if (question == "TIME") {
    which(has_time & !has_arm)
  } else if (question == "ARM") which(has_arm & !has_time) else integer(0)
}

.mira_build_sensitivity <- function(mixed_tests,
                                    rm_anova,
                                    friedman,
                                    gee,
                                    random_slope,
                                    robust,
                                    n_default) {
  rows <- list()
  add <- function(value) {
    if (!is.null(value)) rows[[length(rows) + 1L]] <<- value
  }

  add(.mira_lrt_sensitivity(
    mixed_tests$global_time, "TIME", "Mixed model: random-intercept LRT",
    "$model$global_time_test", n_default
  ))
  add(.mira_lrt_sensitivity(
    mixed_tests$global_arm, "ARM", "Mixed model: random-intercept LRT",
    "$model$global_arm_test", n_default
  ))
  add(.mira_lrt_sensitivity(
    mixed_tests$interaction, "TIME_X_ARM", "Mixed model: random-intercept LRT",
    "$model$arm_time_interaction_test", n_default
  ))

  if (isTRUE(rm_anova$performed) && !is.null(rm_anova$tidy)) {
    for (question in c("TIME", "ARM", "TIME_X_ARM")) {
      index <- .mira_match_effect(rm_anova$tidy, question)
      if (length(index) == 0L) next
      r <- rm_anova$tidy[index[[1L]], , drop = FALSE]
      df_raw <- paste(r$df_num, r$df_den, sep = "/")
      add(.mira_sensitivity_row(
        question, "RM-ANOVA (uncorrected)", r$statistic, df_raw, r$p_value,
        r$partial_eta_squared, "partial eta squared", rm_anova$n_subjects,
        "$advanced_tests$rm_anova$object",
        "Complete-profile repeated-measures analysis."
      ))
      if (question != "ARM" && "p_gg" %in% names(r)) {
        add(.mira_sensitivity_row(
          question, "RM-ANOVA (Greenhouse-Geisser)", r$statistic,
          paste(r$df_num_gg, r$df_den_gg, sep = "/"), r$p_gg,
          r$partial_eta_squared, "partial eta squared", rm_anova$n_subjects,
          "$advanced_tests$rm_anova$tables$greenhouse_geisser",
          "Degrees of freedom and p-value corrected for sphericity."
        ))
      }
      if (question != "ARM" && "p_hf" %in% names(r)) {
        add(.mira_sensitivity_row(
          question, "RM-ANOVA (Huynh-Feldt)", r$statistic,
          paste(r$df_num_hf, r$df_den_hf, sep = "/"), r$p_hf,
          r$partial_eta_squared, "partial eta squared", rm_anova$n_subjects,
          "$advanced_tests$rm_anova$tables$huynh_feldt",
          "Degrees of freedom and p-value corrected for sphericity."
        ))
      }
    }
  }

  if (isTRUE(friedman$performed)) {
    f <- friedman$tidy[1L, ]
    add(.mira_sensitivity_row(
      "TIME", "Friedman test", f$statistic, f$df, f$p_raw,
      f$kendalls_w, "Kendall's W", f$n,
      "$advanced_tests$friedman$test",
      "Rank-based complete-profile global test."
    ))
  }

  for (structure in c("independence", "exchangeable", "ar1")) {
    item <- gee[[structure]]
    if (is.null(item) || !isTRUE(item$performed) || is.null(item$effect_tests)) next
    for (i in seq_len(nrow(item$effect_tests))) {
      r <- item$effect_tests[i, ]
      add(.mira_sensitivity_row(
        r$question, paste0("GEE (", structure, ")"), r$statistic, r$df, r$p_raw,
        NA_real_, NA_character_, tryCatch(stats::nobs(item$model), error = function(e) n_default),
        paste0("$advanced_models$gee$", structure, "$model"),
        "Robust sandwich Wald test under the stated working correlation."
      ))
    }
  }

  add(.mira_lrt_sensitivity(
    random_slope$global_time_test, "TIME", "Mixed model: random-slope LRT",
    "$advanced_models$random_slope$global_time_test", n_default,
    "Random slope uses the ordinal timepoint index."
  ))
  add(.mira_lrt_sensitivity(
    random_slope$global_arm_test, "ARM", "Mixed model: random-slope LRT",
    "$advanced_models$random_slope$global_arm_test", n_default,
    "Random slope uses the ordinal timepoint index."
  ))
  add(.mira_lrt_sensitivity(
    random_slope$arm_time_interaction_test, "TIME_X_ARM",
    "Mixed model: random-slope LRT",
    "$advanced_models$random_slope$arm_time_interaction_test", n_default,
    "Random slope uses the ordinal timepoint index."
  ))

  if (isTRUE(robust$performed) && length(robust$effect_tests) > 0L) {
    for (question in names(robust$effect_tests)) {
      table <- robust$effect_tests[[question]]
      if (is.null(table)) next
      table <- tryCatch(
        as.data.frame(table, stringsAsFactors = FALSE),
        error = function(e) NULL
      )
      if (is.null(table)) next
      if (nrow(table) == 0L) next
      f_col <- .mira_column_name(table, candidates = c("Fstat", "F"))
      dfn_col <- .mira_column_name(table, candidates = c("df_num", "df num"))
      dfd_col <- .mira_column_name(table, candidates = c("df_denom", "df denom"))
      p_col <- .mira_column_name(table, candidates = c("p_val", "p.value", "p"))
      add(.mira_sensitivity_row(
        question, "Mixed model: CR2/HTZ robust Wald test",
        if (!is.na(f_col)) .mira_as_numeric(table[[f_col]])[[1L]] else NA_real_,
        paste(
          if (!is.na(dfn_col)) .mira_as_numeric(table[[dfn_col]])[[1L]] else NA_real_,
          if (!is.na(dfd_col)) .mira_as_numeric(table[[dfd_col]])[[1L]] else NA_real_,
          sep = "/"
        ),
        if (!is.na(p_col)) .mira_as_numeric(table[[p_col]])[[1L]] else NA_real_,
        NA_real_, NA_character_, n_default,
        paste0("$robustness$club_sandwich$effect_tests$", question),
        "CR2 covariance with small-sample HTZ correction."
      ))
    }
  }

  if (length(rows) == 0L) {
    return(data.frame(
      question = character(0), method = character(0), statistic = numeric(0),
      df = character(0), p_raw = numeric(0), effect_size = numeric(0),
      effect_size_type = character(0), n = integer(0), object_path = character(0),
      note = character(0), stringsAsFactors = FALSE
    ))
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  attr(out, "interpretation") <- paste0(
    "Use this table to compare sensitivity across methods. It is not a majority vote ",
    "and does not designate an effect as true based on the number of small p-values."
  )
  out
}

.mira_multiplicity_family <- function(table, primary_method, object_path) {
  list(
    primary_method = primary_method,
    available_adjustments = c("raw", "bonferroni", "holm", "BH", "BY"),
    table = table,
    object_path = if (is.null(table)) NA_character_ else object_path
  )
}

.mira_run_advanced_inference <- function(analysis_data,
                                         long_data,
                                         time_vars,
                                         time_labels,
                                         arm_enabled,
                                         arm_levels,
                                         covariates,
                                         categorical_covariates,
                                         primary_method,
                                         confidence_level,
                                         model_enabled,
                                         mixed_model,
                                         mixed_converged,
                                         mixed_singular,
                                         global_time_test,
                                         global_arm_test,
                                         interaction_test,
                                         change,
                                         icc_model) {
  time_levels <- unname(time_labels[time_vars])
  advanced_data <- .mira_prepare_advanced_data(
    long_data = long_data,
    time_levels = time_levels,
    arm_enabled = arm_enabled,
    arm_levels = arm_levels,
    covariates = covariates,
    categorical_covariates = categorical_covariates
  )
  fixed_formula <- .mira_make_fixed_formula(arm_enabled, covariates)
  friedman <- .mira_run_friedman(
    analysis_data, time_vars, time_labels, primary_method
  )

  if (isTRUE(model_enabled)) {
    emmeans <- .mira_run_emmeans(
      mixed_model, arm_enabled, primary_method, confidence_level
    )
    rm_anova <- .mira_run_rm_anova(
      advanced_data, arm_enabled, covariates, categorical_covariates
    )
    gee <- .mira_run_gee(advanced_data, fixed_formula, confidence_level)
    random_slope <- .mira_run_random_slope(
      advanced_data, arm_enabled, covariates, mixed_model
    )
    nlme_models <- .mira_run_nlme(advanced_data, fixed_formula)
    robust <- .mira_run_robust_inference(mixed_model, confidence_level)
  } else {
    reason <- "Advanced model-based analyses were disabled by model = FALSE."
    emmeans <- .mira_skipped_module("emmeans", reason)
    rm_anova <- .mira_skipped_module("afex", reason)
    random_slope <- .mira_skipped_module("lme4", reason)
    robust <- .mira_skipped_module("clubSandwich", reason)
    gee_item <- .mira_skipped_module("geepack", reason)
    gee <- list(
      performed = FALSE, package = "geepack", independence = gee_item,
      exchangeable = gee_item, ar1 = gee_item, warnings = character(0),
      error = NULL, reason_skipped = reason
    )
    nlme_item <- .mira_skipped_module("nlme", reason)
    nlme_models <- list(
      performed = FALSE, package = "nlme", compound_symmetry = nlme_item,
      ar1 = nlme_item, comparison = NULL, warnings = character(0),
      error = NULL, reason_skipped = reason
    )
  }

  random_intercept_r2 <- .mira_run_r2(mixed_model)
  random_slope_r2 <- .mira_run_r2(random_slope$model)
  model_comparison <- .mira_build_model_comparison(
    advanced_data, fixed_formula, mixed_model, mixed_converged, mixed_singular,
    random_slope, nlme_models, gee
  )
  sensitivity <- .mira_build_sensitivity(
    mixed_tests = list(
      global_time = global_time_test,
      global_arm = global_arm_test,
      interaction = interaction_test
    ),
    rm_anova = rm_anova,
    friedman = friedman,
    gee = gee,
    random_slope = random_slope,
    robust = robust,
    n_default = nrow(advanced_data)
  )

  paired_effects <- if (!is.null(change) && nrow(change) > 0L) {
    change[c("from", "to", "from_label", "to_label", "n", "cohens_dz")]
  } else data.frame()
  friedman_effects <- if (!is.null(friedman$posthoc$tidy) &&
                          nrow(friedman$posthoc$tidy) > 0L) {
    friedman$posthoc$tidy[c(
      "from", "to", "from_label", "to_label", "n", "rank_biserial"
    )]
  } else data.frame()
  rm_effects <- if (isTRUE(rm_anova$performed) && !is.null(rm_anova$tidy)) {
    rm_anova$tidy[c(
      "effect", "partial_eta_squared", "generalized_eta_squared"
    )]
  } else data.frame()

  emmeans_ok <- isTRUE(emmeans$performed)
  multiplicity <- list(
    primary_method = primary_method,
    note = "P-values are adjusted within scientifically distinct contrast families.",
    time_pairwise = .mira_multiplicity_family(
      if (emmeans_ok) emmeans$time$pairwise_tidy else NULL,
      primary_method,
      "$advanced_tests$emmeans$time$pairwise"
    ),
    baseline_vs_followup = .mira_multiplicity_family(
      if (emmeans_ok) emmeans$baseline_followup$tidy else NULL,
      primary_method,
      "$advanced_tests$emmeans$baseline_followup$object"
    ),
    consecutive_time = .mira_multiplicity_family(
      if (emmeans_ok) emmeans$consecutive$tidy else NULL,
      primary_method,
      "$advanced_tests$emmeans$consecutive$object"
    ),
    ordinal_trends = .mira_multiplicity_family(
      if (emmeans_ok) emmeans$trends$tidy else NULL,
      primary_method,
      "$advanced_tests$emmeans$trends$object"
    ),
    arm_pairwise = .mira_multiplicity_family(
      if (emmeans_ok) emmeans$arm$pairwise_tidy else NULL,
      primary_method,
      "$advanced_tests$emmeans$arm$pairwise"
    ),
    simple_effects = list(
      arm_within_time = .mira_multiplicity_family(
        if (emmeans_ok) emmeans$arm_time$simple_arm_tidy else NULL,
        primary_method,
        "$advanced_tests$emmeans$arm_time$simple_arm"
      ),
      time_within_arm = .mira_multiplicity_family(
        if (emmeans_ok) emmeans$arm_time$simple_time_tidy else NULL,
        primary_method,
        "$advanced_tests$emmeans$arm_time$simple_time"
      )
    ),
    interaction_contrasts = .mira_multiplicity_family(
      if (emmeans_ok) emmeans$arm_time$interaction_tidy else NULL,
      primary_method,
      "$advanced_tests$emmeans$arm_time$interaction"
    ),
    friedman_posthoc = .mira_multiplicity_family(
      friedman$posthoc$tidy,
      primary_method,
      "$advanced_tests$friedman$posthoc$tests"
    )
  )

  list(
    advanced_tests = list(
      emmeans = emmeans,
      rm_anova = rm_anova,
      friedman = friedman
    ),
    advanced_models = list(
      random_slope = random_slope,
      nlme = nlme_models,
      gee = gee,
      model_comparison = model_comparison
    ),
    robustness = list(club_sandwich = robust),
    effect_sizes = list(
      paired_cohens_dz = paired_effects,
      paired_rank_biserial = friedman_effects,
      rm_anova = rm_effects,
      friedman = list(
        kendalls_w = friedman$tidy$kendalls_w,
        object_path = if (!is.null(friedman$test)) {
          "$advanced_tests$friedman$test"
        } else NA_character_
      ),
      mixed_models = list(
        random_intercept = random_intercept_r2,
        random_slope = random_slope_r2,
        ICC = icc_model
      )
    ),
    multiplicity = multiplicity,
    sensitivity = sensitivity,
    data = advanced_data,
    formula = fixed_formula
  )
}

.mira_failed_advanced <- function(error, warnings = character(0)) {
  reason <- "Advanced longitudinal inference was unavailable because its protected wrapper failed."
  skipped <- .mira_skipped_module(NA_character_, reason)
  skipped$error <- error
  skipped$warnings <- warnings
  empty_sensitivity <- data.frame(
    question = character(0), method = character(0), statistic = numeric(0),
    df = character(0), p_raw = numeric(0), effect_size = numeric(0),
    effect_size_type = character(0), n = integer(0), object_path = character(0),
    note = character(0), stringsAsFactors = FALSE
  )
  list(
    advanced_tests = list(emmeans = skipped, rm_anova = skipped, friedman = skipped),
    advanced_models = list(
      random_slope = skipped,
      nlme = list(compound_symmetry = skipped, ar1 = skipped, comparison = NULL),
      gee = list(independence = skipped, exchangeable = skipped, ar1 = skipped),
      model_comparison = data.frame()
    ),
    robustness = list(club_sandwich = skipped),
    effect_sizes = list(),
    multiplicity = list(),
    sensitivity = empty_sensitivity,
    data = NULL,
    formula = NULL,
    error = error,
    warnings = warnings
  )
}

.mira_analyse_single_outcome <- function(data,
                                         id = "patient",
                                         time_vars = NULL,
                                         time_labels = NULL,
                                         arm = NULL,
                                         reference_arm = NULL,
                                         arm_tests = TRUE,
                                         alpha = 0.05,
                                         plots = TRUE,
                                         model = TRUE,
                                         outliers = TRUE,
                                         correlations = TRUE,
                                         verbose = TRUE,
                                         p_adjust_method = "holm",
                                         improvement_direction = c("unknown", "higher", "lower"),
                                         stable_threshold = 0,
                                         strict_id = TRUE,
                                         .outcome_name = NULL,
                                         .covariates = character(0),
                                         .categorical_covariates = character(0)) {

  # ----------------------------------------------------------
  # 0. INPUT VALIDATION
  # ----------------------------------------------------------

  if (!is.data.frame(data)) {
    stop("data must be a data.frame.", call. = FALSE)
  }

  if (nrow(data) == 0L) {
    stop("The data set contains no observations.", call. = FALSE)
  }

  if (!is.character(id) || length(id) != 1L || is.na(id) || !nzchar(id)) {
    stop("id must name exactly one variable.", call. = FALSE)
  }

  if (!id %in% names(data)) {
    stop(sprintf("The ID variable '%s' does not exist in the data set.", id), call. = FALSE)
  }

  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) ||
      alpha <= 0 || alpha >= 1) {
    stop("alpha must be a number strictly between 0 and 1.", call. = FALSE)
  }

  scalar_flag <- function(x, name) {
    if (!is.logical(x) || length(x) != 1L || is.na(x)) {
      stop(sprintf("%s must be TRUE or FALSE.", name), call. = FALSE)
    }
    invisible(TRUE)
  }

  scalar_flag(plots, "plots")
  scalar_flag(model, "model")
  scalar_flag(outliers, "outliers")
  scalar_flag(correlations, "correlations")
  scalar_flag(verbose, "verbose")
  scalar_flag(strict_id, "strict_id")
  scalar_flag(arm_tests, "arm_tests")

  improvement_direction <- match.arg(improvement_direction)

  if (!is.numeric(stable_threshold) || length(stable_threshold) != 1L ||
      !is.finite(stable_threshold) || stable_threshold < 0) {
    stop("stable_threshold must be a finite number >= 0.", call. = FALSE)
  }

  if (!is.character(p_adjust_method) || length(p_adjust_method) != 1L ||
      is.na(p_adjust_method) || !p_adjust_method %in% p.adjust.methods) {
    stop(
      sprintf(
        "p_adjust_method must be one of: %s.",
        paste(p.adjust.methods, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  # Wide-format repeated-measures data should have one row per subject.
  id_values <- data[[id]]
  missing_id_n <- sum(is.na(id_values))
  duplicated_id_n <- sum(duplicated(id_values[!is.na(id_values)]))

  if (strict_id && missing_id_n > 0L) {
    stop(
      sprintf("The ID variable contains %d missing values.", missing_id_n),
      call. = FALSE
    )
  }

  if (strict_id && duplicated_id_n > 0L) {
    stop(
      sprintf(
        paste0(
          "The ID variable contains %d duplicate values. ",
          "mira_info() assumes one row per subject in wide format. ",
          "Correct the duplicates, or use strict_id = FALSE only if duplicate IDs are intentional."
        ),
        duplicated_id_n
      ),
      call. = FALSE
    )
  }

  if (!strict_id && missing_id_n > 0L) {
    warning(sprintf("There are %d missing IDs.", missing_id_n), call. = FALSE)
  }

  if (!strict_id && duplicated_id_n > 0L) {
    warning(
      sprintf(
        "%d duplicated IDs will be treated as observations from the same subject.",
        duplicated_id_n
      ),
      call. = FALSE
    )
  }

  # ----------------------------------------------------------
  # 0B. TREATMENT ARM DETECTION / VALIDATION
  # ----------------------------------------------------------

  # Detection is performed by the public configuration layer. A NULL value here
  # is therefore authoritative and keeps between-group analyses disabled.

  arm_available <- !is.null(arm)
  arm_levels <- character(0)
  arm_counts <- NULL
  missing_arm_n <- 0L

  if (arm_available) {

    if (!is.character(arm) || length(arm) != 1L || is.na(arm) || !nzchar(arm)) {
      stop("arm must be NULL or the name of exactly one variable.", call. = FALSE)
    }

    if (!arm %in% names(data)) {
      stop(sprintf("The treatment-arm variable '%s' does not exist in the data set.", arm), call. = FALSE)
    }

    arm_raw <- data[[arm]]
    missing_arm_n <- sum(is.na(arm_raw) | !nzchar(trimws(as.character(arm_raw))))

    arm_character <- as.character(arm_raw)
    arm_character[is.na(arm_raw) | !nzchar(trimws(arm_character))] <- NA_character_

    observed_arms <- unique(arm_character[!is.na(arm_character)])

    if (length(observed_arms) < 2L) {
      warning(
        sprintf(
          "Variable '%s' contains fewer than two observed groups; treatment-arm analyses will be disabled.",
          arm
        ),
        call. = FALSE
      )
      arm_tests <- FALSE
    }

    if (is.null(reference_arm)) {
      if (length(observed_arms) > 0L) {
        reference_arm <- observed_arms[[1L]]
      }
    } else {
      reference_arm <- as.character(reference_arm)

      if (length(reference_arm) != 1L || is.na(reference_arm) || !nzchar(reference_arm)) {
        stop("reference_arm must identify exactly one non-empty group.", call. = FALSE)
      }

      if (!reference_arm %in% observed_arms) {
        stop(
          sprintf(
            "reference_arm='%s' is not present in the variable '%s'. Observed groups: %s.",
            reference_arm,
            arm,
            paste(observed_arms, collapse = ", ")
          ),
          call. = FALSE
        )
      }
    }

    if (length(observed_arms) > 0L) {
      arm_levels <- c(reference_arm, setdiff(observed_arms, reference_arm))
      arm_counts <- table(
        factor(
          arm_character,
          levels = arm_levels
        ),
        useNA = "no"
      )
    }

    if (missing_arm_n > 0L) {
      warning(
        sprintf(
          "%d subjects have a missing treatment arm and will be excluded only from between-group analyses.",
          missing_arm_n
        ),
        call. = FALSE
      )
    }
  }

  # ----------------------------------------------------------
  # 1. DETECT OUTCOME + TIMEPOINTS
  # ----------------------------------------------------------

  supplied_time_vars <- !is.null(time_vars)

  if (is.null(time_vars)) {
    detected <- .mira_detect_longitudinal_variables(data, "auto")$groups
    if (length(detected) != 1L) {
      stop(
        paste0(
          "The internal analyzer requires exactly one outcome. Use mira_info() ",
          "for automatic multi-outcome configuration."
        ),
        call. = FALSE
      )
    }
    detected_group <- detected[[1L]]
    time_vars <- detected_group$variables
    outcome_name <- detected_group$outcome
    detected_time_labels <- detected_group$time_labels
  } else {
    if (!is.character(time_vars) || length(time_vars) < 2L || anyNA(time_vars)) {
      stop("time_vars must contain at least two variable names.", call. = FALSE)
    }

    if (anyDuplicated(time_vars)) {
      stop("time_vars contains duplicate names.", call. = FALSE)
    }

    missing_vars <- setdiff(time_vars, names(data))
    if (length(missing_vars) > 0L) {
      stop(
        sprintf(
          "The following variables do not exist in the data set: %s",
          paste(missing_vars, collapse = ", ")
        ),
        call. = FALSE
      )
    }

    parsed <- lapply(time_vars, .mira_parse_longitudinal_name, variable_pattern = "auto")
    parsed_names <- vapply(parsed, function(x) {
      if (is.null(x)) NA_character_ else x$outcome
    }, character(1L))
    parsed_labels <- vapply(seq_along(parsed), function(i) {
      if (is.null(parsed[[i]])) time_vars[[i]] else parsed[[i]]$time_label
    }, character(1L))
    parsed_order <- vapply(parsed, function(x) {
      if (is.null(x)) NA_real_ else x$time_order
    }, numeric(1L))

    if (!is.null(.outcome_name)) {
      if (!is.character(.outcome_name) || length(.outcome_name) != 1L ||
          is.na(.outcome_name) || !nzchar(.outcome_name)) {
        stop(".outcome_name must be exactly one non-empty label.", call. = FALSE)
      }
      outcome_name <- .outcome_name
    } else {
      inferred <- .mira_compact_unique(parsed_names)
      outcome_name <- if (length(inferred) == 1L) inferred[[1L]] else "outcome"
    }

    clear_order <- all(is.finite(parsed_order)) && !anyDuplicated(parsed_order)
    ord <- if (clear_order) order(parsed_order, time_vars) else seq_along(time_vars)
    if (!is.null(time_labels)) {
      if (length(time_labels) != length(time_vars)) {
        stop("time_labels must have the same length as time_vars.", call. = FALSE)
      }
      time_labels <- time_labels[ord]
    }
    time_vars <- time_vars[ord]
    detected_time_labels <- parsed_labels[ord]
  }

  outcome_display <- toupper(outcome_name)

  if (!is.character(.covariates) || anyNA(.covariates) || any(!nzchar(.covariates))) {
    stop(".covariates must be a valid vector of names.", call. = FALSE)
  }
  missing_covariates <- setdiff(.covariates, names(data))
  if (length(missing_covariates) > 0L) {
    stop(sprintf("Covariates not found: %s.", paste(missing_covariates, collapse = ", ")),
         call. = FALSE)
  }
  if (!is.character(.categorical_covariates) || anyNA(.categorical_covariates) ||
      any(!nzchar(.categorical_covariates)) ||
      length(setdiff(.categorical_covariates, .covariates)) > 0L) {
    stop(".categorical_covariates must be a subset of .covariates.",
         call. = FALSE)
  }

  # ----------------------------------------------------------
  # 2. TIME VARIABLE VALIDATION
  # ----------------------------------------------------------

  non_numeric <- time_vars[
    !vapply(data[time_vars], is.numeric, logical(1L))
  ]

  if (length(non_numeric) > 0L) {
    stop(
      sprintf(
        "Longitudinal variables must be numeric: %s",
        paste(non_numeric, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  # ----------------------------------------------------------
  # 3. TIME LABELS
  # ----------------------------------------------------------

  if (is.null(time_labels)) {
    time_labels <- detected_time_labels
  } else {
    if (length(time_labels) != length(time_vars)) {
      stop("time_labels must have the same length as time_vars.", call. = FALSE)
    }
    time_labels <- as.character(time_labels)
  }

  if (anyNA(time_labels) || any(!nzchar(trimws(time_labels)))) {
    stop("time_labels cannot contain NA or empty labels.", call. = FALSE)
  }

  if (anyDuplicated(time_labels)) {
    stop("time_labels must contain unique labels.", call. = FALSE)
  }

  names(time_labels) <- time_vars

  # ----------------------------------------------------------
  # 4. CLEAN ANALYSIS COPY + HELPERS
  # ----------------------------------------------------------

  raw_data <- data
  analysis_data <- data

  non_finite_by_time <- vapply(
    raw_data[time_vars],
    function(x) sum(!is.na(x) & !is.finite(x)),
    integer(1L)
  )

  total_non_finite <- sum(non_finite_by_time)

  if (total_non_finite > 0L) {
    warning(
      sprintf(
        "%d Inf/-Inf values in longitudinal variables will be treated as NA in analyses.",
        total_non_finite
      ),
      call. = FALSE
    )
  }

  analysis_data[time_vars] <- lapply(
    analysis_data[time_vars],
    function(x) {
      x[!is.finite(x)] <- NA_real_
      x
    }
  )

  safe_mean <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0L) NA_real_ else mean(x)
  }

  safe_sd <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) < 2L) NA_real_ else stats::sd(x)
  }

  safe_var <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) < 2L) NA_real_ else stats::var(x)
  }

  safe_median <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0L) NA_real_ else stats::median(x)
  }

  safe_quantile <- function(x, p) {
    x <- x[is.finite(x)]
    if (length(x) == 0L) {
      NA_real_
    } else {
      as.numeric(stats::quantile(x, probs = p, names = FALSE, type = 7))
    }
  }

  mean_ci <- function(x) {
    x <- x[is.finite(x)]
    n <- length(x)

    if (n == 0L) {
      return(c(mean = NA_real_, se = NA_real_, lower = NA_real_, upper = NA_real_))
    }

    m <- mean(x)
    if (n < 2L) {
      return(c(mean = m, se = NA_real_, lower = NA_real_, upper = NA_real_))
    }

    s <- stats::sd(x)
    se <- s / sqrt(n)
    crit <- stats::qt(1 - alpha / 2, df = n - 1L)
    c(mean = m, se = se, lower = m - crit * se, upper = m + crit * se)
  }

  safe_t_p <- function(x, y) {
    if (length(x) < 2L || length(y) < 2L) return(NA_real_)
    z <- tryCatch(
      stats::t.test(y, x, paired = TRUE, conf.level = 1 - alpha),
      error = function(e) NULL
    )
    if (is.null(z)) NA_real_ else unname(z$p.value)
  }

  safe_wilcox_p <- function(x, y) {
    if (length(x) < 1L || length(y) < 1L) return(NA_real_)
    z <- tryCatch(
      suppressWarnings(stats::wilcox.test(y, x, paired = TRUE, exact = FALSE)),
      error = function(e) NULL
    )
    if (is.null(z)) NA_real_ else unname(z$p.value)
  }

  adjust_p <- function(p) {
    out <- rep(NA_real_, length(p))
    ok <- is.finite(p)
    if (any(ok)) out[ok] <- stats::p.adjust(p[ok], method = p_adjust_method)
    out
  }

  safe_welch_test <- function(x, y) {
    x <- x[is.finite(x)]
    y <- y[is.finite(y)]

    if (length(x) < 2L || length(y) < 2L) {
      return(list(
        p.value = NA_real_,
        estimate = NA_real_,
        conf.int = c(NA_real_, NA_real_)
      ))
    }

    z <- tryCatch(
      stats::t.test(y, x, paired = FALSE, var.equal = FALSE, conf.level = 1 - alpha),
      error = function(e) NULL
    )

    if (is.null(z)) {
      return(list(
        p.value = NA_real_,
        estimate = NA_real_,
        conf.int = c(NA_real_, NA_real_)
      ))
    }

    list(
      p.value = unname(z$p.value),
      estimate = mean(y) - mean(x),
      conf.int = unname(z$conf.int)
    )
  }

  safe_unpaired_wilcox_p <- function(x, y) {
    x <- x[is.finite(x)]
    y <- y[is.finite(y)]

    if (length(x) < 1L || length(y) < 1L) return(NA_real_)

    z <- tryCatch(
      suppressWarnings(
        stats::wilcox.test(y, x, paired = FALSE, exact = FALSE)
      ),
      error = function(e) NULL
    )

    if (is.null(z)) NA_real_ else unname(z$p.value)
  }

  hedges_g <- function(x, y) {
    x <- x[is.finite(x)]
    y <- y[is.finite(y)]

    nx <- length(x)
    ny <- length(y)

    if (nx < 2L || ny < 2L) return(NA_real_)

    vx <- stats::var(x)
    vy <- stats::var(y)
    df <- nx + ny - 2L

    if (!is.finite(vx) || !is.finite(vy) || df <= 0L) return(NA_real_)

    pooled_var <- ((nx - 1L) * vx + (ny - 1L) * vy) / df

    if (!is.finite(pooled_var) || pooled_var <= 0) return(NA_real_)

    d <- (mean(y) - mean(x)) / sqrt(pooled_var)

    # Small-sample correction.
    J <- if (df > 1L) 1 - 3 / (4 * df - 1) else 1

    J * d
  }

  safe_welch_anova_p <- function(value, group) {
    keep <- is.finite(value) & !is.na(group)
    value <- value[keep]
    group <- droplevels(factor(group[keep]))

    if (length(value) < 3L || nlevels(group) < 2L) return(NA_real_)

    z <- tryCatch(
      stats::oneway.test(value ~ group, var.equal = FALSE),
      error = function(e) NULL
    )

    if (is.null(z)) NA_real_ else unname(z$p.value)
  }

  safe_kruskal_p <- function(value, group) {
    keep <- is.finite(value) & !is.na(group)
    value <- value[keep]
    group <- droplevels(factor(group[keep]))

    if (length(value) < 2L || nlevels(group) < 2L) return(NA_real_)

    z <- tryCatch(
      stats::kruskal.test(value ~ group),
      error = function(e) NULL
    )

    if (is.null(z)) NA_real_ else unname(z$p.value)
  }

  improvement_counts <- function(delta) {
    increased <- delta > stable_threshold
    decreased <- delta < -stable_threshold
    stable <- !(increased | decreased)

    if (improvement_direction == "higher") {
      improved <- increased
      worsened <- decreased
    } else if (improvement_direction == "lower") {
      improved <- decreased
      worsened <- increased
    } else {
      improved <- rep(NA, length(delta))
      worsened <- rep(NA, length(delta))
    }

    list(
      increased_n = sum(increased),
      decreased_n = sum(decreased),
      stable_n = sum(stable),
      improved_n = if (improvement_direction == "unknown") NA_integer_ else sum(improved),
      worsened_n = if (improvement_direction == "unknown") NA_integer_ else sum(worsened)
    )
  }

  pct <- function(n, denom) {
    if (length(n) == 0L || is.na(n) || denom <= 0L) NA_real_ else n / denom * 100
  }

  confidence_level <- 1 - alpha
  confidence_percent <- 100 * confidence_level

  # ----------------------------------------------------------
  # 5. OVERVIEW
  # ----------------------------------------------------------

  n_rows <- nrow(analysis_data)
  n_patients <- length(unique(id_values[!is.na(id_values)]))
  complete_profiles <- sum(stats::complete.cases(analysis_data[time_vars]))

  # ----------------------------------------------------------
  # 6. DESCRIPTIVE STATISTICS
  # ----------------------------------------------------------

  descriptive_list <- lapply(time_vars, function(v) {
    raw_x <- raw_data[[v]]
    x <- analysis_data[[v]]
    finite_x <- x[is.finite(x)]
    n <- length(finite_x)
    na_n <- sum(is.na(raw_x))
    nonfinite_n <- sum(!is.na(raw_x) & !is.finite(raw_x))
    unavailable_n <- sum(is.na(x))
    ci <- mean_ci(x)
    q1 <- safe_quantile(x, 0.25)
    q3 <- safe_quantile(x, 0.75)
    m <- safe_mean(x)
    s <- safe_sd(x)

    data.frame(
      time = v,
      label = unname(time_labels[v]),
      n = n,
      missing = na_n,
      missing_pct = na_n / n_rows * 100,
      non_finite = nonfinite_n,
      non_finite_pct = nonfinite_n / n_rows * 100,
      unavailable = unavailable_n,
      unavailable_pct = unavailable_n / n_rows * 100,
      mean = m,
      sd = s,
      variance = safe_var(x),
      se = unname(ci["se"]),
      ci_lower = unname(ci["lower"]),
      ci_upper = unname(ci["upper"]),
      median = safe_median(x),
      q1 = q1,
      q3 = q3,
      iqr = if (is.na(q1) || is.na(q3)) NA_real_ else q3 - q1,
      min = if (n > 0L) min(finite_x) else NA_real_,
      max = if (n > 0L) max(finite_x) else NA_real_,
      cv_percent = if (!is.na(m) && !is.na(s) && m != 0) s / abs(m) * 100 else NA_real_,
      stringsAsFactors = FALSE
    )
  })

  descriptives <- do.call(rbind, descriptive_list)
  rownames(descriptives) <- NULL

  # Compatibility / convenience aliases.
  descriptives[[paste0(outcome_name, "_mean")]] <- descriptives$mean
  descriptives[[paste0(outcome_name, "_sd")]] <- descriptives$sd
  descriptives[[paste0(outcome_name, "_median")]] <- descriptives$median
  descriptives[[paste0(outcome_name, "_ci_lower")]] <- descriptives$ci_lower
  descriptives[[paste0(outcome_name, "_ci_upper")]] <- descriptives$ci_upper

  # ----------------------------------------------------------
  # 7. MISSINGNESS / DATA AVAILABILITY
  # ----------------------------------------------------------

  missing_by_patient <- data.frame(
    patient = analysis_data[[id]],
    missing_n = rowSums(is.na(raw_data[time_vars])),
    non_finite_n = rowSums(vapply(
      raw_data[time_vars],
      function(x) !is.na(x) & !is.finite(x),
      logical(n_rows)
    )),
    unavailable_n = rowSums(is.na(analysis_data[time_vars])),
    stringsAsFactors = FALSE
  )

  missing_by_patient$unavailable_pct <-
    missing_by_patient$unavailable_n / length(time_vars) * 100

  missing_summary <- descriptives[c(
    "time", "label", "missing", "missing_pct",
    "non_finite", "non_finite_pct", "unavailable", "unavailable_pct"
  )]

  names(missing_summary)[names(missing_summary) == "missing"] <- "missing_n"
  names(missing_summary)[names(missing_summary) == "non_finite"] <- "non_finite_n"
  names(missing_summary)[names(missing_summary) == "unavailable"] <- "unavailable_n"

  # ----------------------------------------------------------
  # 8. PAIRWISE CHANGE ANALYSIS
  # ----------------------------------------------------------

  change_list <- list()
  counter <- 1L

  for (i in seq_len(length(time_vars) - 1L)) {
    for (j in seq.int(i + 1L, length(time_vars))) {
      v1 <- time_vars[[i]]
      v2 <- time_vars[[j]]
      x <- analysis_data[[v1]]
      y <- analysis_data[[v2]]
      keep <- stats::complete.cases(x, y)
      x_complete <- x[keep]
      y_complete <- y[keep]
      delta <- y_complete - x_complete
      n <- length(delta)

      if (n == 0L) next

      delta_mean <- mean(delta)
      delta_sd <- if (n >= 2L) stats::sd(delta) else NA_real_
      delta_se <- if (n >= 2L) delta_sd / sqrt(n) else NA_real_

      if (n >= 2L && is.finite(delta_se)) {
        tcrit <- stats::qt(1 - alpha / 2, df = n - 1L)
        ci_low <- delta_mean - tcrit * delta_se
        ci_high <- delta_mean + tcrit * delta_se
      } else {
        ci_low <- NA_real_
        ci_high <- NA_real_
      }

      effect_size <- if (n >= 2L && is.finite(delta_sd) && delta_sd > 0) {
        delta_mean / delta_sd
      } else {
        NA_real_
      }

      counts <- improvement_counts(delta)

      change_list[[counter]] <- data.frame(
        from = v1,
        to = v2,
        from_label = unname(time_labels[v1]),
        to_label = unname(time_labels[v2]),
        n = n,
        mean_from = mean(x_complete),
        mean_to = mean(y_complete),
        mean_change = delta_mean,
        sd_change = delta_sd,
        se_change = delta_se,
        ci_lower = ci_low,
        ci_upper = ci_high,
        cohens_dz = effect_size,
        increased_n = counts$increased_n,
        increased_pct = pct(counts$increased_n, n),
        decreased_n = counts$decreased_n,
        decreased_pct = pct(counts$decreased_n, n),
        stable_n = counts$stable_n,
        stable_pct = pct(counts$stable_n, n),
        improved_n = counts$improved_n,
        improved_pct = pct(counts$improved_n, n),
        worsened_n = counts$worsened_n,
        worsened_pct = pct(counts$worsened_n, n),
        paired_t_p = safe_t_p(x_complete, y_complete),
        wilcoxon_p = safe_wilcox_p(x_complete, y_complete),
        stringsAsFactors = FALSE
      )

      counter <- counter + 1L
    }
  }

  if (length(change_list) > 0L) {
    change <- do.call(rbind, change_list)
    rownames(change) <- NULL
    change$paired_t_p_adj <- adjust_p(change$paired_t_p)
    change$wilcoxon_p_adj <- adjust_p(change$wilcoxon_p)
  } else {
    change <- data.frame(
      from = character(0), to = character(0),
      from_label = character(0), to_label = character(0),
      n = integer(0), mean_from = numeric(0), mean_to = numeric(0),
      mean_change = numeric(0), sd_change = numeric(0), se_change = numeric(0),
      ci_lower = numeric(0), ci_upper = numeric(0), cohens_dz = numeric(0),
      increased_n = integer(0), increased_pct = numeric(0),
      decreased_n = integer(0), decreased_pct = numeric(0),
      stable_n = integer(0), stable_pct = numeric(0),
      improved_n = integer(0), improved_pct = numeric(0),
      worsened_n = integer(0), worsened_pct = numeric(0),
      paired_t_p = numeric(0), paired_t_p_adj = numeric(0),
      wilcoxon_p = numeric(0), wilcoxon_p_adj = numeric(0),
      stringsAsFactors = FALSE
    )
  }

  if (nrow(change) > 0L) {
    change[[paste0(outcome_name, "_change")]] <- change$mean_change
    change[[paste0(outcome_name, "_mean_from")]] <- change$mean_from
    change[[paste0(outcome_name, "_mean_to")]] <- change$mean_to
  }

  # ----------------------------------------------------------
  # 9. CORRELATIONS + PAIRWISE N
  # ----------------------------------------------------------

  correlation_pearson <- NULL
  correlation_spearman <- NULL
  correlation_n <- NULL

  if (correlations) {
    cor_data <- analysis_data[time_vars]

    correlation_pearson <- tryCatch(
      stats::cor(cor_data, use = "pairwise.complete.obs", method = "pearson"),
      error = function(e) NULL
    )

    correlation_spearman <- tryCatch(
      stats::cor(cor_data, use = "pairwise.complete.obs", method = "spearman"),
      error = function(e) NULL
    )

    correlation_n <- matrix(
      NA_integer_,
      nrow = length(time_vars),
      ncol = length(time_vars),
      dimnames = list(time_vars, time_vars)
    )

    for (i in seq_along(time_vars)) {
      for (j in seq_along(time_vars)) {
        correlation_n[i, j] <- sum(stats::complete.cases(
          cor_data[[i]], cor_data[[j]]
        ))
      }
    }
  }

  # ----------------------------------------------------------
  # 10. LONG DATA
  # ----------------------------------------------------------

  long_list <- lapply(seq_along(time_vars), function(k) {
    v <- time_vars[[k]]
    row_block <- data.frame(
      patient = analysis_data[[id]],
      outcome = outcome_name,
      time = v,
      time_index = k,
      time_label = unname(time_labels[v]),
      arm = if (arm_available) as.character(analysis_data[[arm]]) else NA_character_,
      value = analysis_data[[v]],
      stringsAsFactors = FALSE
    )
    if (length(.covariates) > 0L) {
      row_block <- cbind(row_block, analysis_data[.covariates])
    }
    row_block
  })

  long_data <- do.call(rbind, long_list)
  rownames(long_data) <- NULL
  long_data$time_label <- factor(
    long_data$time_label,
    levels = unname(time_labels[time_vars])
  )

  if (arm_available) {
    long_data$arm <- factor(
      long_data$arm,
      levels = arm_levels
    )
  }

  # ----------------------------------------------------------
  # 10B. ARM-SPECIFIC EXPLORATORY ANALYSES
  # ----------------------------------------------------------

  arm_descriptives <- NULL
  arm_time_omnibus <- NULL
  arm_time_pairwise <- NULL
  arm_change_descriptives <- NULL
  arm_change_omnibus <- NULL
  arm_change_pairwise <- NULL
  arm_missingness <- NULL
  baseline_balance <- NULL

  if (arm_available && arm_tests && length(arm_levels) >= 2L) {

    # --------------------------------------------------------
    # Descriptives by arm and time
    # --------------------------------------------------------

    arm_desc_list <- list()
    a_counter <- 1L

    for (g in arm_levels) {
      for (k in seq_along(time_vars)) {

        v <- time_vars[[k]]
        idx <- !is.na(analysis_data[[arm]]) &
          as.character(analysis_data[[arm]]) == g

        x_all <- analysis_data[[v]][idx]
        x <- x_all[is.finite(x_all)]
        group_n <- sum(idx)
        n <- length(x)
        unavailable_n <- sum(!is.finite(x_all) | is.na(x_all))

        ci <- mean_ci(x)

        arm_desc_list[[a_counter]] <- data.frame(
          arm = g,
          time = v,
          time_index = k,
          time_label = unname(time_labels[v]),
          group_n = group_n,
          n = n,
          unavailable_n = unavailable_n,
          unavailable_pct = if (group_n > 0L) unavailable_n / group_n * 100 else NA_real_,
          mean = safe_mean(x),
          sd = safe_sd(x),
          se = unname(ci["se"]),
          ci_lower = unname(ci["lower"]),
          ci_upper = unname(ci["upper"]),
          median = safe_median(x),
          q1 = safe_quantile(x, 0.25),
          q3 = safe_quantile(x, 0.75),
          min = if (n > 0L) min(x) else NA_real_,
          max = if (n > 0L) max(x) else NA_real_,
          stringsAsFactors = FALSE
        )

        a_counter <- a_counter + 1L
      }
    }

    arm_descriptives <- do.call(rbind, arm_desc_list)
    rownames(arm_descriptives) <- NULL


    # --------------------------------------------------------
    # Missingness / availability by arm
    # --------------------------------------------------------

    arm_missing_list <- lapply(seq_along(time_vars), function(k) {

      v <- time_vars[[k]]
      group <- factor(
        as.character(analysis_data[[arm]]),
        levels = arm_levels
      )
      unavailable <- !is.finite(analysis_data[[v]]) | is.na(analysis_data[[v]])

      keep <- !is.na(group)

      tab <- table(
        group[keep],
        unavailable[keep]
      )

      test_p <- NA_real_
      test_used <- NA_character_

      if (nrow(tab) >= 2L && ncol(tab) >= 2L) {

        expected <- tryCatch(
          suppressWarnings(stats::chisq.test(tab)$expected),
          error = function(e) NULL
        )

        if (!is.null(expected) && any(expected < 5)) {
          ft <- tryCatch(stats::fisher.test(tab), error = function(e) NULL)
          if (!is.null(ft)) {
            test_p <- unname(ft$p.value)
            test_used <- "Fisher"
          }
        } else {
          ct <- tryCatch(
            suppressWarnings(stats::chisq.test(tab, correct = FALSE)),
            error = function(e) NULL
          )
          if (!is.null(ct)) {
            test_p <- unname(ct$p.value)
            test_used <- "Chi-square"
          }
        }
      }

      data.frame(
        time = v,
        time_index = k,
        time_label = unname(time_labels[v]),
        test = test_used,
        p_value = test_p,
        stringsAsFactors = FALSE
      )
    })

    arm_missingness <- do.call(rbind, arm_missing_list)
    arm_missingness$p_value_adj <- adjust_p(arm_missingness$p_value)


    # --------------------------------------------------------
    # Omnibus arm test at each timepoint
    # --------------------------------------------------------

    omnibus_list <- lapply(seq_along(time_vars), function(k) {
      v <- time_vars[[k]]
      value <- analysis_data[[v]]
      group <- factor(
        as.character(analysis_data[[arm]]),
        levels = arm_levels
      )

      data.frame(
        time = v,
        time_index = k,
        time_label = unname(time_labels[v]),
        n = sum(is.finite(value) & !is.na(group)),
        n_arms = length(unique(as.character(group[is.finite(value) & !is.na(group)]))),
        welch_anova_p = safe_welch_anova_p(value, group),
        kruskal_p = safe_kruskal_p(value, group),
        stringsAsFactors = FALSE
      )
    })

    arm_time_omnibus <- do.call(rbind, omnibus_list)
    arm_time_omnibus$welch_anova_p_adj <- adjust_p(arm_time_omnibus$welch_anova_p)
    arm_time_omnibus$kruskal_p_adj <- adjust_p(arm_time_omnibus$kruskal_p)


    # --------------------------------------------------------
    # Pairwise arm comparisons at each timepoint
    # Difference is arm_b - arm_a.
    # --------------------------------------------------------

    pair_list <- list()
    p_counter <- 1L
    arm_pairs <- utils::combn(arm_levels, 2L, simplify = FALSE)

    for (k in seq_along(time_vars)) {

      v <- time_vars[[k]]

      for (pair in arm_pairs) {

        a <- pair[[1L]]
        b <- pair[[2L]]

        xa <- analysis_data[[v]][as.character(analysis_data[[arm]]) == a]
        xb <- analysis_data[[v]][as.character(analysis_data[[arm]]) == b]

        xa <- xa[is.finite(xa)]
        xb <- xb[is.finite(xb)]

        wt <- safe_welch_test(xa, xb)

        pair_list[[p_counter]] <- data.frame(
          time = v,
          time_index = k,
          time_label = unname(time_labels[v]),
          arm_a = a,
          arm_b = b,
          n_a = length(xa),
          n_b = length(xb),
          mean_a = safe_mean(xa),
          mean_b = safe_mean(xb),
          mean_difference_b_minus_a = wt$estimate,
          ci_lower = wt$conf.int[[1L]],
          ci_upper = wt$conf.int[[2L]],
          hedges_g = hedges_g(xa, xb),
          welch_t_p = wt$p.value,
          wilcoxon_p = safe_unpaired_wilcox_p(xa, xb),
          stringsAsFactors = FALSE
        )

        p_counter <- p_counter + 1L
      }
    }

    arm_time_pairwise <- do.call(rbind, pair_list)
    arm_time_pairwise$welch_t_p_adj <- adjust_p(arm_time_pairwise$welch_t_p)
    arm_time_pairwise$wilcoxon_p_adj <- adjust_p(arm_time_pairwise$wilcoxon_p)

    baseline_balance <- arm_time_pairwise[
      arm_time_pairwise$time_index == 1L,
      ,
      drop = FALSE
    ]


    # --------------------------------------------------------
    # Baseline-to-follow-up change by arm
    # --------------------------------------------------------

    change_desc_list <- list()
    change_omnibus_list <- list()
    change_pair_list <- list()
    cd_counter <- 1L
    cp_counter <- 1L

    baseline_v <- time_vars[[1L]]

    for (k in seq.int(2L, length(time_vars))) {

      follow_v <- time_vars[[k]]
      baseline_values <- analysis_data[[baseline_v]]
      follow_values <- analysis_data[[follow_v]]
      delta_all <- follow_values - baseline_values

      # Arm-specific change descriptives
      for (g in arm_levels) {

        idx <- as.character(analysis_data[[arm]]) == g
        delta <- delta_all[idx]
        delta <- delta[is.finite(delta)]
        n <- length(delta)
        ci <- mean_ci(delta)
        counts <- improvement_counts(delta)

        change_desc_list[[cd_counter]] <- data.frame(
          arm = g,
          from = baseline_v,
          to = follow_v,
          from_label = unname(time_labels[baseline_v]),
          to_label = unname(time_labels[follow_v]),
          time_index = k,
          n = n,
          mean_change = safe_mean(delta),
          sd_change = safe_sd(delta),
          se_change = unname(ci["se"]),
          ci_lower = unname(ci["lower"]),
          ci_upper = unname(ci["upper"]),
          median_change = safe_median(delta),
          improved_n = counts$improved_n,
          improved_pct = pct(counts$improved_n, n),
          worsened_n = counts$worsened_n,
          worsened_pct = pct(counts$worsened_n, n),
          stable_n = counts$stable_n,
          stable_pct = pct(counts$stable_n, n),
          stringsAsFactors = FALSE
        )

        cd_counter <- cd_counter + 1L
      }

      group <- factor(
        as.character(analysis_data[[arm]]),
        levels = arm_levels
      )

      change_omnibus_list[[k - 1L]] <- data.frame(
        from = baseline_v,
        to = follow_v,
        to_label = unname(time_labels[follow_v]),
        time_index = k,
        n = sum(is.finite(delta_all) & !is.na(group)),
        welch_anova_p = safe_welch_anova_p(delta_all, group),
        kruskal_p = safe_kruskal_p(delta_all, group),
        stringsAsFactors = FALSE
      )

      for (pair in arm_pairs) {

        a <- pair[[1L]]
        b <- pair[[2L]]

        da <- delta_all[as.character(analysis_data[[arm]]) == a]
        db <- delta_all[as.character(analysis_data[[arm]]) == b]

        da <- da[is.finite(da)]
        db <- db[is.finite(db)]

        wt <- safe_welch_test(da, db)

        change_pair_list[[cp_counter]] <- data.frame(
          from = baseline_v,
          to = follow_v,
          to_label = unname(time_labels[follow_v]),
          time_index = k,
          arm_a = a,
          arm_b = b,
          n_a = length(da),
          n_b = length(db),
          mean_change_a = safe_mean(da),
          mean_change_b = safe_mean(db),
          difference_in_change_b_minus_a = wt$estimate,
          ci_lower = wt$conf.int[[1L]],
          ci_upper = wt$conf.int[[2L]],
          hedges_g = hedges_g(da, db),
          welch_t_p = wt$p.value,
          wilcoxon_p = safe_unpaired_wilcox_p(da, db),
          stringsAsFactors = FALSE
        )

        cp_counter <- cp_counter + 1L
      }
    }

    arm_change_descriptives <- do.call(rbind, change_desc_list)
    rownames(arm_change_descriptives) <- NULL

    arm_change_omnibus <- do.call(rbind, change_omnibus_list)
    rownames(arm_change_omnibus) <- NULL
    arm_change_omnibus$welch_anova_p_adj <- adjust_p(arm_change_omnibus$welch_anova_p)
    arm_change_omnibus$kruskal_p_adj <- adjust_p(arm_change_omnibus$kruskal_p)

    arm_change_pairwise <- do.call(rbind, change_pair_list)
    rownames(arm_change_pairwise) <- NULL
    arm_change_pairwise$welch_t_p_adj <- adjust_p(arm_change_pairwise$welch_t_p)
    arm_change_pairwise$wilcoxon_p_adj <- adjust_p(arm_change_pairwise$wilcoxon_p)
  }

  # ----------------------------------------------------------
  # 11. WITHIN / BETWEEN SUBJECT VARIABILITY
  # ----------------------------------------------------------

  individual_means <- apply(
    analysis_data[time_vars],
    1L,
    safe_mean
  )

  between_sd <- safe_sd(individual_means)
  grand_mean <- safe_mean(long_data$value)

  patient_means <- tapply(long_data$value, long_data$patient, safe_mean)
  patient_key <- as.character(long_data$patient)
  long_data$patient_mean <- unname(patient_means[patient_key])
  residuals_within <- long_data$value - long_data$patient_mean
  within_sd <- safe_sd(residuals_within)

  # Balanced-data one-way ICC(1,1) using only complete profiles.
  icc_anova <- NA_real_
  ms_between <- NA_real_
  ms_within <- NA_real_
  df_between <- NA_real_
  df_within <- NA_real_

  complete_profile_idx <- stats::complete.cases(analysis_data[time_vars]) &
    !is.na(analysis_data[[id]])

  complete_wide <- analysis_data[complete_profile_idx, c(id, time_vars), drop = FALSE]

  if (nrow(complete_wide) >= 2L && length(time_vars) >= 2L) {
    k_actual <- length(time_vars)
    y_matrix <- as.matrix(complete_wide[time_vars])
    subject_means <- rowMeans(y_matrix)
    grand <- mean(y_matrix)
    ss_between <- k_actual * sum((subject_means - grand)^2)
    ss_within <- sum((y_matrix - subject_means)^2)
    df_between <- nrow(y_matrix) - 1L
    df_within <- nrow(y_matrix) * (k_actual - 1L)

    if (df_between > 0L && df_within > 0L) {
      ms_between <- ss_between / df_between
      ms_within <- ss_within / df_within
      denom <- ms_between + (k_actual - 1L) * ms_within
      if (is.finite(denom) && denom != 0) {
        icc_anova <- (ms_between - ms_within) / denom
      }
    }
  }

  # ----------------------------------------------------------
  # 12. MIXED MODEL + MODEL-BASED ICC + GLOBAL TIME TEST
  # ----------------------------------------------------------

  mixed_model <- NULL
  model_summary <- NULL
  model_anova <- NULL
  global_time_test <- NULL
  global_arm_test <- NULL
  arm_time_interaction_test <- NULL
  model_error <- NULL
  model_warnings <- character(0)
  model_singular <- NA
  model_converged <- NA
  icc_model <- NA_real_
  model_covariates_requested <- .covariates
  model_covariates_used <- character(0)
  model_covariates_skipped <- character(0)
  model_fixed_parameters <- NA_integer_
  model_data <- NULL

  if (model) {
    if (!requireNamespace("lme4", quietly = TRUE)) {
      model_error <- "Package 'lme4' is not installed."
    } else {
      model_data <- long_data[
        !is.na(long_data$patient) & is.finite(long_data$value),
        ,
        drop = FALSE
      ]

      model_data$patient_factor <- factor(model_data$patient)
      model_data$time_factor <- factor(
        as.character(model_data$time_label),
        levels = unname(time_labels[time_vars])
      )

      model_uses_arm <- arm_available && arm_tests && length(arm_levels) >= 2L
      if (model_uses_arm) {
        model_data$arm_factor <- factor(
          as.character(model_data$arm),
          levels = arm_levels
        )
      }

      for (variable in .covariates) {
        x <- model_data[[variable]]
        if (is.numeric(x)) x[!is.finite(x)] <- NA_real_
        if (variable %in% .categorical_covariates && !is.factor(x)) x <- factor(x)
        if (is.character(x)) {
          x[!is.na(x) & !nzchar(trimws(x))] <- NA_character_
          x <- factor(x)
        }
        if (is.factor(x)) x <- droplevels(x)
        model_data[[variable]] <- x
        observed <- x[!is.na(x)]
        if (length(unique(observed)) >= 2L) {
          model_covariates_used <- c(model_covariates_used, variable)
        } else {
          model_covariates_skipped <- c(model_covariates_skipped, variable)
        }
      }

      required_model_columns <- c(
        "patient_factor", "time_factor", "value",
        if (model_uses_arm) "arm_factor" else character(0),
        model_covariates_used
      )
      model_data <- model_data[
        stats::complete.cases(model_data[required_model_columns]),
        ,
        drop = FALSE
      ]
      model_data$patient_factor <- droplevels(model_data$patient_factor)
      model_data$time_factor <- droplevels(model_data$time_factor)
      if (model_uses_arm) model_data$arm_factor <- droplevels(model_data$arm_factor)
      post_filter_skipped <- model_covariates_used[vapply(
        model_covariates_used,
        function(variable) length(unique(model_data[[variable]][
          !is.na(model_data[[variable]])
        ])) < 2L,
        logical(1L)
      )]
      if (length(post_filter_skipped) > 0L) {
        model_covariates_used <- setdiff(model_covariates_used, post_filter_skipped)
        model_covariates_skipped <- unique(c(model_covariates_skipped, post_filter_skipped))
      }
      covariate_df <- sum(vapply(model_covariates_used, function(variable) {
        x <- model_data[[variable]]
        if (is.factor(x) || is.character(x) || is.logical(x)) {
          max(1L, length(unique(x[!is.na(x)])) - 1L)
        } else 1L
      }, integer(1L)))
      core_df <- if (model_uses_arm) {
        nlevels(model_data$time_factor) * nlevels(model_data$arm_factor)
      } else {
        nlevels(model_data$time_factor)
      }
      model_fixed_parameters <- as.integer(core_df + covariate_df)

      if (nrow(model_data) < 3L || nlevels(model_data$patient_factor) < 2L ||
          nlevels(model_data$time_factor) < 2L ||
          nrow(model_data) <= model_fixed_parameters + 1L) {
        model_error <- paste0(
          "Insufficient data to estimate a mixed-effects model after applying ",
          "completeness requirements for the outcome, treatment arm, and covariates."
        )
      } else {
        fit_with_warnings <- function(expr) {
          withCallingHandlers(
            expr,
            warning = function(w) {
              model_warnings <<- unique(c(model_warnings, conditionMessage(w)))
              invokeRestart("muffleWarning")
            }
          )
        }

        quote_formula_name <- function(x) {
          paste0("`", gsub("`", "", x, fixed = TRUE), "`")
        }
        covariate_terms <- vapply(
          model_covariates_used,
          quote_formula_name,
          character(1L)
        )
        make_model_formula <- function(core) {
          stats::as.formula(paste(
            "value ~",
            paste(c(core, covariate_terms, "(1 | patient_factor)"), collapse = " + ")
          ))
        }

        model_formula <- if (model_uses_arm && nlevels(model_data$arm_factor) >= 2L) {
          make_model_formula("time_factor * arm_factor")
        } else {
          make_model_formula("time_factor")
        }

        mixed_model <- tryCatch(
          fit_with_warnings(
            if (requireNamespace("lmerTest", quietly = TRUE)) {
              lmerTest::lmer(
                model_formula,
                data = model_data,
                REML = TRUE,
                na.action = stats::na.omit
              )
            } else {
              lme4::lmer(
                model_formula,
                data = model_data,
                REML = TRUE,
                na.action = stats::na.omit
              )
            }
          ),
          error = function(e) {
            model_error <<- conditionMessage(e)
            NULL
          }
        )

        if (!is.null(mixed_model)) {
          model_summary <- summary(mixed_model)
          model_anova <- tryCatch(stats::anova(mixed_model), error = function(e) NULL)
          model_singular <- tryCatch(
            lme4::isSingular(mixed_model, tol = 1e-4),
            error = function(e) NA
          )

          optinfo <- tryCatch(mixed_model@optinfo, error = function(e) NULL)
          conv_messages <- if (!is.null(optinfo)) optinfo$conv$lme4$messages else NULL
          opt_ok <- if (!is.null(optinfo) && !is.null(optinfo$conv$opt)) {
            isTRUE(optinfo$conv$opt == 0)
          } else {
            TRUE
          }
          model_converged <- opt_ok && is.null(conv_messages)

          vc <- tryCatch(as.data.frame(lme4::VarCorr(mixed_model)), error = function(e) NULL)
          if (!is.null(vc)) {
            subject_var <- vc$vcov[vc$grp == "patient_factor" & is.na(vc$var2)]
            residual_var <- vc$vcov[vc$grp == "Residual"]
            if (length(subject_var) >= 1L && length(residual_var) >= 1L) {
              denom <- subject_var[[1L]] + residual_var[[1L]]
              if (is.finite(denom) && denom > 0) {
                icc_model <- subject_var[[1L]] / denom
              }
            }
          }

          # Robust ML likelihood-ratio tests.
          if (model_uses_arm &&
              "arm_factor" %in% names(model_data) &&
              nlevels(model_data$arm_factor) >= 2L) {

            full_ml <- tryCatch(
              fit_with_warnings(
                lme4::lmer(
                  make_model_formula("time_factor * arm_factor"),
                  data = model_data,
                  REML = FALSE,
                  na.action = stats::na.omit
                )
              ),
              error = function(e) NULL
            )

            additive_ml <- tryCatch(
              fit_with_warnings(
                lme4::lmer(
                  make_model_formula("time_factor + arm_factor"),
                  data = model_data,
                  REML = FALSE,
                  na.action = stats::na.omit
                )
              ),
              error = function(e) NULL
            )

            no_time_ml <- tryCatch(
              fit_with_warnings(
                lme4::lmer(
                  make_model_formula("arm_factor"),
                  data = model_data,
                  REML = FALSE,
                  na.action = stats::na.omit
                )
              ),
              error = function(e) NULL
            )

            no_arm_ml <- tryCatch(
              fit_with_warnings(
                lme4::lmer(
                  make_model_formula("time_factor"),
                  data = model_data,
                  REML = FALSE,
                  na.action = stats::na.omit
                )
              ),
              error = function(e) NULL
            )

            if (!is.null(full_ml) && !is.null(no_time_ml)) {
              global_time_test <- tryCatch(
                stats::anova(no_time_ml, full_ml),
                error = function(e) NULL
              )
            }

            if (!is.null(full_ml) && !is.null(no_arm_ml)) {
              global_arm_test <- tryCatch(
                stats::anova(no_arm_ml, full_ml),
                error = function(e) NULL
              )
            }

            if (!is.null(full_ml) && !is.null(additive_ml)) {
              arm_time_interaction_test <- tryCatch(
                stats::anova(additive_ml, full_ml),
                error = function(e) NULL
              )
            }

          } else {

            global_time_test <- tryCatch({
              full_ml <- fit_with_warnings(
                lme4::lmer(
                  make_model_formula("time_factor"),
                  data = model_data,
                  REML = FALSE,
                  na.action = stats::na.omit
                )
              )

              null_ml <- fit_with_warnings(
                lme4::lmer(
                  make_model_formula("1"),
                  data = model_data,
                  REML = FALSE,
                  na.action = stats::na.omit
                )
              )

              stats::anova(null_ml, full_ml)
            }, error = function(e) NULL)
          }
        }
      }
    }
  }

  variability <- data.frame(
    outcome = outcome_name,
    grand_mean = grand_mean,
    between_subject_sd = between_sd,
    within_subject_sd = within_sd,
    ICC_anova_complete_profiles = icc_anova,
    ICC_model = icc_model,
    ms_between = ms_between,
    ms_within = ms_within,
    df_between = df_between,
    df_within = df_within,
    stringsAsFactors = FALSE
  )

  # Backward-friendly alias: prefer model ICC if available.
  variability$ICC <- if (is.finite(icc_model)) icc_model else icc_anova
  variability[[paste0(outcome_name, "_grand_mean")]] <- variability$grand_mean
  variability[[paste0(outcome_name, "_between_subject_sd")]] <- variability$between_subject_sd
  variability[[paste0(outcome_name, "_within_subject_sd")]] <- variability$within_subject_sd

  # ----------------------------------------------------------
  # 12B. ADVANCED LONGITUDINAL INFERENCE
  # ----------------------------------------------------------

  advanced_covariates <- if (isTRUE(model) &&
                             requireNamespace("lme4", quietly = TRUE)) {
    model_covariates_used
  } else {
    .covariates
  }
  advanced_eval <- .mira_eval(
    .mira_run_advanced_inference(
      analysis_data = analysis_data,
      long_data = if (!is.null(model_data)) model_data else long_data,
      time_vars = time_vars,
      time_labels = time_labels,
      arm_enabled = arm_available && arm_tests && length(arm_levels) >= 2L,
      arm_levels = arm_levels,
      covariates = advanced_covariates,
      categorical_covariates = intersect(
        .categorical_covariates,
        advanced_covariates
      ),
      primary_method = p_adjust_method,
      confidence_level = confidence_level,
      model_enabled = model,
      mixed_model = mixed_model,
      mixed_converged = model_converged,
      mixed_singular = model_singular,
      global_time_test = global_time_test,
      global_arm_test = global_arm_test,
      interaction_test = arm_time_interaction_test,
      change = change,
      icc_model = icc_model
    )
  )
  advanced_results <- if (is.null(advanced_eval$value)) {
    .mira_failed_advanced(advanced_eval$error, advanced_eval$warnings)
  } else {
    advanced_eval$value
  }

  # ----------------------------------------------------------
  # 13. OUTLIERS (IQR FLAGS; NOT AUTOMATIC EXCLUSIONS)
  # ----------------------------------------------------------

  outlier_results <- list()
  change_outliers <- list()

  if (outliers) {
    for (v in time_vars) {
      x <- analysis_data[[v]]
      finite_x <- x[is.finite(x)]

      if (length(finite_x) < 4L) {
        outlier_results[[v]] <- data.frame(
          row = integer(0), patient = analysis_data[[id]][integer(0)],
          value = numeric(0), lower_bound = numeric(0), upper_bound = numeric(0),
          stringsAsFactors = FALSE
        )
        next
      }

      q1 <- as.numeric(stats::quantile(finite_x, 0.25, names = FALSE))
      q3 <- as.numeric(stats::quantile(finite_x, 0.75, names = FALSE))
      iqr_value <- q3 - q1
      lower <- q1 - 1.5 * iqr_value
      upper <- q3 + 1.5 * iqr_value
      idx <- which(is.finite(x) & (x < lower | x > upper))

      outlier_results[[v]] <- data.frame(
        row = idx,
        patient = analysis_data[[id]][idx],
        value = x[idx],
        lower_bound = rep(lower, length(idx)),
        upper_bound = rep(upper, length(idx)),
        stringsAsFactors = FALSE
      )
    }

    if (nrow(change) > 0L) {
      for (r in seq_len(nrow(change))) {
        v1 <- change$from[[r]]
        v2 <- change$to[[r]]
        keep <- stats::complete.cases(analysis_data[[v1]], analysis_data[[v2]])
        delta <- analysis_data[[v2]][keep] - analysis_data[[v1]][keep]
        patients <- analysis_data[[id]][keep]
        nm <- paste(v1, v2, sep = "_to_")

        if (length(delta) < 4L) {
          change_outliers[[nm]] <- data.frame(
            patient = patients[integer(0)], change = numeric(0),
            lower_bound = numeric(0), upper_bound = numeric(0),
            stringsAsFactors = FALSE
          )
          next
        }

        q1 <- as.numeric(stats::quantile(delta, 0.25, names = FALSE))
        q3 <- as.numeric(stats::quantile(delta, 0.75, names = FALSE))
        iqr_value <- q3 - q1
        lower <- q1 - 1.5 * iqr_value
        upper <- q3 + 1.5 * iqr_value
        idx <- which(delta < lower | delta > upper)

        change_outliers[[nm]] <- data.frame(
          patient = patients[idx],
          change = delta[idx],
          lower_bound = rep(lower, length(idx)),
          upper_bound = rep(upper, length(idx)),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  # ----------------------------------------------------------
  # 14. PATIENT TRAJECTORIES
  # ----------------------------------------------------------

  baseline <- analysis_data[[time_vars[[1L]]]]
  final <- analysis_data[[time_vars[[length(time_vars)]]]]
  absolute_change <- final - baseline

  trajectory_summary <- data.frame(
    patient = analysis_data[[id]],
    outcome = outcome_name,
    arm = if (arm_available) as.character(analysis_data[[arm]]) else NA_character_,
    baseline = baseline,
    final = final,
    absolute_change = absolute_change,
    relative_change_percent = ifelse(
      is.finite(baseline) & baseline != 0 & is.finite(final),
      (final - baseline) / abs(baseline) * 100,
      NA_real_
    ),
    stringsAsFactors = FALSE
  )

  trajectory_summary$direction <- ifelse(
    is.na(trajectory_summary$absolute_change),
    NA_character_,
    ifelse(
      trajectory_summary$absolute_change > stable_threshold,
      "Increase",
      ifelse(
        trajectory_summary$absolute_change < -stable_threshold,
        "Decrease",
        "Stable"
      )
    )
  )

  trajectory_summary$clinical_direction <- if (improvement_direction == "unknown") {
    rep(NA_character_, nrow(trajectory_summary))
  } else if (improvement_direction == "higher") {
    ifelse(
      is.na(trajectory_summary$direction), NA_character_,
      ifelse(trajectory_summary$direction == "Increase", "Improved",
             ifelse(trajectory_summary$direction == "Decrease", "Worsened", "Stable"))
    )
  } else {
    ifelse(
      is.na(trajectory_summary$direction), NA_character_,
      ifelse(trajectory_summary$direction == "Decrease", "Improved",
             ifelse(trajectory_summary$direction == "Increase", "Worsened", "Stable"))
    )
  }

  trajectory_summary[[paste0(outcome_name, "_baseline")]] <- trajectory_summary$baseline
  trajectory_summary[[paste0(outcome_name, "_final")]] <- trajectory_summary$final
  trajectory_summary[[paste0(outcome_name, "_change")]] <- trajectory_summary$absolute_change
  trajectory_summary[[paste0(outcome_name, "_relative_change_percent")]] <-
    trajectory_summary$relative_change_percent

  # ----------------------------------------------------------
  # 15. PLOTS (ggplot2 ONLY; PUBLICATION-READY)
  # ----------------------------------------------------------

  plots_list <- list()
  plot_error <- NULL

  if (plots) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      plot_error <- "Package 'ggplot2' is not installed."
    } else {
      ci_text <- paste0(formatC(confidence_percent, format = "fg", digits = 4), "% CI")

      # Shared paper-oriented theme. Plot titles/subtitles are deliberately omitted
      # so figures can be pasted directly into manuscripts and captioned externally.
      paper_theme <-
        ggplot2::theme_minimal(base_size = 11) +
        ggplot2::theme(
          plot.title = ggplot2::element_blank(),
          plot.subtitle = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          legend.position = "bottom",
          legend.title = ggplot2::element_text(size = 10),
          axis.title = ggplot2::element_text(size = 10),
          axis.text = ggplot2::element_text(size = 9)
        )

      # Boxplot summary uses the already validated descriptives table.
      boxplot_summary <- descriptives
      boxplot_summary$time_label_n <- paste0(
        boxplot_summary$label,
        "\n(N = ",
        boxplot_summary$n,
        ")"
      )

      long_plot <- long_data[is.finite(long_data$value), , drop = FALSE]
      long_plot$time_label <- factor(
        long_plot$time_label,
        levels = unname(time_labels[time_vars])
      )

      # 1) Distribution at each timepoint with raw observations + mean CI.
      plots_list$boxplot <-
        ggplot2::ggplot(long_plot, ggplot2::aes(x = time_label, y = value)) +
        ggplot2::geom_boxplot(
          width = 0.55, fill = "grey92", color = "grey30",
          linewidth = 0.7, outlier.shape = NA, na.rm = TRUE
        ) +
        ggplot2::geom_jitter(
          width = 0.10, height = 0, alpha = 0.22,
          size = 1.5, color = "grey35", na.rm = TRUE
        ) +
        ggplot2::geom_errorbar(
          data = boxplot_summary,
          ggplot2::aes(x = label, ymin = ci_lower, ymax = ci_upper),
          inherit.aes = FALSE, width = 0.10, linewidth = 0.9,
          color = "#2C7BB6", na.rm = TRUE
        ) +
        ggplot2::geom_point(
          data = boxplot_summary,
          ggplot2::aes(x = label, y = mean),
          inherit.aes = FALSE, shape = 21, size = 3.8,
          stroke = 1.1, fill = "white", color = "#2C7BB6", na.rm = TRUE
        ) +
        ggplot2::scale_x_discrete(labels = boxplot_summary$time_label_n) +
        ggplot2::labs(x = NULL, y = outcome_display) +
        paper_theme

      # Build patient-to-patient consecutive segments without dplyr.
      segment_list <- list()
      s_counter <- 1L
      split_long <- split(long_plot, long_plot$patient, drop = TRUE)

      for (patient_data in split_long) {
        patient_data <- patient_data[order(patient_data$time_index), , drop = FALSE]
        if (nrow(patient_data) < 2L) next

        for (k in seq_len(nrow(patient_data) - 1L)) {
          dlt <- patient_data$value[[k + 1L]] - patient_data$value[[k]]
          direction <- if (dlt > stable_threshold) {
            "Increase"
          } else if (dlt < -stable_threshold) {
            "Decrease"
          } else {
            "Stable"
          }

          segment_list[[s_counter]] <- data.frame(
            patient = patient_data$patient[[k]],
            x = patient_data$time_index[[k]],
            xend = patient_data$time_index[[k + 1L]],
            y = patient_data$value[[k]],
            yend = patient_data$value[[k + 1L]],
            trajectory_direction = direction,
            stringsAsFactors = FALSE
          )
          s_counter <- s_counter + 1L
        }
      }

      if (length(segment_list) > 0L) {
        trajectory_segments <- do.call(rbind, segment_list)
      } else {
        trajectory_segments <- data.frame(
          patient = character(0), x = numeric(0), xend = numeric(0),
          y = numeric(0), yend = numeric(0), trajectory_direction = character(0),
          stringsAsFactors = FALSE
        )
      }

      trajectory_summary_plot <- data.frame(
        time_index = seq_along(time_vars),
        time_label = unname(time_labels[time_vars]),
        n = descriptives$n,
        mean = descriptives$mean,
        ci_lower = descriptives$ci_lower,
        ci_upper = descriptives$ci_upper,
        stringsAsFactors = FALSE
      )
      trajectory_summary_plot$time_label_n <- paste0(
        trajectory_summary_plot$time_label,
        "\n(N = ", trajectory_summary_plot$n, ")"
      )

      # 2) Individual longitudinal trajectories + population mean and CI.
      plots_list$spaghetti <-
        ggplot2::ggplot() +
        ggplot2::geom_ribbon(
          data = trajectory_summary_plot,
          ggplot2::aes(x = time_index, ymin = ci_lower, ymax = ci_upper),
          inherit.aes = FALSE, alpha = 0.15, na.rm = TRUE
        ) +
        ggplot2::geom_segment(
          data = trajectory_segments,
          ggplot2::aes(
            x = x, xend = xend, y = y, yend = yend,
            group = patient, color = trajectory_direction
          ),
          alpha = 0.22, linewidth = 0.6, na.rm = TRUE
        ) +
        ggplot2::geom_point(
          data = long_plot,
          ggplot2::aes(x = time_index, y = value),
          color = "grey40", alpha = 0.20, size = 1.4, na.rm = TRUE
        ) +
        ggplot2::geom_line(
          data = trajectory_summary_plot,
          ggplot2::aes(x = time_index, y = mean),
          linewidth = 1.3, color = "black", na.rm = TRUE
        ) +
        ggplot2::geom_point(
          data = trajectory_summary_plot,
          ggplot2::aes(x = time_index, y = mean),
          shape = 21, size = 3.8, stroke = 1.1,
          fill = "white", color = "black", na.rm = TRUE
        ) +
        ggplot2::scale_color_manual(
          values = c(Increase = "#2C7BB6", Decrease = "#D7191C", Stable = "grey60"),
          breaks = c("Increase", "Decrease", "Stable"),
          name = "Consecutive change"
        ) +
        ggplot2::scale_x_continuous(
          breaks = trajectory_summary_plot$time_index,
          labels = trajectory_summary_plot$time_label_n
        ) +
        ggplot2::labs(x = NULL, y = outcome_display) +
        paper_theme

      # 3) Mean outcome and CI over time.
      mean_ci_plot_data <- descriptives
      mean_ci_plot_data$label <- factor(
        mean_ci_plot_data$label,
        levels = unname(time_labels[time_vars])
      )

      plots_list$mean_ci <-
        ggplot2::ggplot(
          mean_ci_plot_data,
          ggplot2::aes(x = label, y = mean, group = 1)
        ) +
        ggplot2::geom_line(color = "#2C7FB8", linewidth = 0.8, na.rm = TRUE) +
        ggplot2::geom_point(size = 3, color = "#2C7FB8", na.rm = TRUE) +
        ggplot2::geom_errorbar(
          ggplot2::aes(ymin = ci_lower, ymax = ci_upper),
          width = 0.1, na.rm = TRUE
        ) +
        ggplot2::labs(x = "Time", y = paste("Mean", outcome_display)) +
        paper_theme

      trajectory_plot_data <- trajectory_summary[
        is.finite(trajectory_summary$absolute_change),
        ,
        drop = FALSE
      ]

      # 4) Baseline-to-final individual change distribution.
      plots_list$change <-
        ggplot2::ggplot(
          trajectory_plot_data,
          ggplot2::aes(x = "All patients", y = absolute_change)
        ) +
        ggplot2::geom_boxplot(alpha = 0.7, na.rm = TRUE) +
        ggplot2::geom_jitter(width = 0.08, alpha = 0.30, na.rm = TRUE) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
        ggplot2::labs(x = NULL, y = paste("Final - Baseline", outcome_display)) +
        paper_theme

      # 5) NEW: distribution of change from baseline at every follow-up.
      baseline_change_list <- lapply(seq.int(2L, length(time_vars)), function(k) {
        delta <- analysis_data[[time_vars[[k]]]] - analysis_data[[time_vars[[1L]]]]
        keep <- is.finite(delta)
        data.frame(
          patient = analysis_data[[id]][keep],
          time_index = k,
          time_label = unname(time_labels[time_vars[[k]]]),
          change_from_baseline = delta[keep],
          stringsAsFactors = FALSE
        )
      })
      baseline_change_long <- do.call(rbind, baseline_change_list)
      baseline_change_long$time_label <- factor(
        baseline_change_long$time_label,
        levels = unname(time_labels[time_vars[-1L]])
      )

      plots_list$change_from_baseline <-
        ggplot2::ggplot(
          baseline_change_long,
          ggplot2::aes(x = time_label, y = change_from_baseline)
        ) +
        ggplot2::geom_boxplot(
          width = 0.55, outlier.shape = NA, fill = "grey92", color = "grey30",
          na.rm = TRUE
        ) +
        ggplot2::geom_jitter(
          width = 0.10, alpha = 0.25, size = 1.4, color = "grey35", na.rm = TRUE
        ) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
        ggplot2::labs(x = "Time", y = paste("Change from baseline", outcome_display)) +
        paper_theme

      # 6) NEW: forest-style plot of mean change from baseline and CI.
      baseline_change_ci <- change[
        change$from == time_vars[[1L]],
        c("to", "to_label", "n", "mean_change", "ci_lower", "ci_upper"),
        drop = FALSE
      ]

      if (nrow(baseline_change_ci) > 0L) {
        baseline_change_ci$to_label <- factor(
          baseline_change_ci$to_label,
          levels = rev(unname(time_labels[time_vars[-1L]]))
        )

        plots_list$change_ci <-
          ggplot2::ggplot(
            baseline_change_ci,
            ggplot2::aes(x = to_label, y = mean_change)
          ) +
          ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
          ggplot2::geom_errorbar(
            ggplot2::aes(ymin = ci_lower, ymax = ci_upper),
            width = 0.12, linewidth = 0.7, na.rm = TRUE
          ) +
          ggplot2::geom_point(size = 2.8, na.rm = TRUE) +
          ggplot2::coord_flip() +
          ggplot2::labs(
            x = NULL,
            y = paste("Mean change from baseline (", ci_text, ")", sep = "")
          ) +
          paper_theme
      }

      # 7) NEW: data availability / missingness by timepoint.
      missing_plot_data <- missing_summary
      missing_plot_data$label <- factor(
        missing_plot_data$label,
        levels = unname(time_labels[time_vars])
      )

      plots_list$missingness <-
        ggplot2::ggplot(
          missing_plot_data,
          ggplot2::aes(x = label, y = unavailable_pct)
        ) +
        ggplot2::geom_col(width = 0.62, fill = "grey55") +
        ggplot2::geom_text(
          ggplot2::aes(label = sprintf("%.1f%%", unavailable_pct)),
          vjust = -0.35, size = 3, na.rm = TRUE
        ) +
        ggplot2::scale_y_continuous(
          limits = c(0, max(c(5, missing_plot_data$unavailable_pct * 1.15), na.rm = TRUE)),
          expand = ggplot2::expansion(mult = c(0, 0.03))
        ) +
        ggplot2::labs(x = "Time", y = "Unavailable observations (%)") +
        paper_theme

      # 8) NEW: Pearson correlation heatmap across repeated timepoints.
      if (correlations && !is.null(correlation_pearson)) {
        cor_mat <- correlation_pearson
        cor_idx <- expand.grid(
          row = seq_len(nrow(cor_mat)),
          col = seq_len(ncol(cor_mat)),
          KEEP.OUT.ATTRS = FALSE,
          stringsAsFactors = FALSE
        )
        cor_plot_data <- data.frame(
          x = unname(time_labels[time_vars])[cor_idx$col],
          y = unname(time_labels[time_vars])[cor_idx$row],
          r = cor_mat[cbind(cor_idx$row, cor_idx$col)],
          stringsAsFactors = FALSE
        )
        cor_plot_data$x <- factor(cor_plot_data$x, levels = unname(time_labels[time_vars]))
        cor_plot_data$y <- factor(
          cor_plot_data$y,
          levels = rev(unname(time_labels[time_vars]))
        )

        plots_list$correlation_heatmap <-
          ggplot2::ggplot(
            cor_plot_data,
            ggplot2::aes(x = x, y = y, fill = r)
          ) +
          ggplot2::geom_tile(color = "white", linewidth = 0.5) +
          ggplot2::geom_text(
            ggplot2::aes(label = ifelse(is.finite(r), sprintf("%.2f", r), "")),
            size = 3
          ) +
          ggplot2::scale_fill_gradient2(
            low = "#2166AC", mid = "white", high = "#B2182B",
            midpoint = 0, limits = c(-1, 1), na.value = "grey95",
            name = "Pearson r"
          ) +
          ggplot2::coord_equal() +
          ggplot2::labs(x = NULL, y = NULL) +
          paper_theme +
          ggplot2::theme(panel.grid = ggplot2::element_blank())
      }

      # 9) NEW: responder/direction percentages for baseline-to-final change.
      response_variable <- if (improvement_direction == "unknown") "direction" else "clinical_direction"
      response_values <- trajectory_summary[[response_variable]]
      response_values <- response_values[!is.na(response_values)]

      if (length(response_values) > 0L) {
        response_tab <- table(response_values)
        response_plot_data <- data.frame(
          category = names(response_tab),
          n = as.integer(response_tab),
          percent = as.integer(response_tab) / sum(response_tab) * 100,
          stringsAsFactors = FALSE
        )

        preferred_order <- if (improvement_direction == "unknown") {
          c("Increase", "Stable", "Decrease")
        } else {
          c("Improved", "Stable", "Worsened")
        }
        response_plot_data$category <- factor(
          response_plot_data$category,
          levels = preferred_order[preferred_order %in% response_plot_data$category]
        )

        plots_list$response <-
          ggplot2::ggplot(
            response_plot_data,
            ggplot2::aes(x = category, y = percent)
          ) +
          ggplot2::geom_col(width = 0.62, fill = "grey55") +
          ggplot2::geom_text(
            ggplot2::aes(label = sprintf("%.1f%%\n(n=%d)", percent, n)),
            vjust = -0.25, size = 3
          ) +
          ggplot2::scale_y_continuous(
            limits = c(0, max(c(10, response_plot_data$percent * 1.18), na.rm = TRUE)),
            expand = ggplot2::expansion(mult = c(0, 0.03))
          ) +
          ggplot2::labs(x = NULL, y = "Patients (%)") +
          paper_theme
      }

      if (arm_available && arm_tests && !is.null(arm_descriptives)) {

        arm_plot_data <- arm_descriptives
        arm_plot_data$time_label <- factor(
          arm_plot_data$time_label,
          levels = unname(time_labels[time_vars])
        )
        arm_plot_data$arm <- factor(
          arm_plot_data$arm,
          levels = arm_levels
        )

        # 10) Mean and CI by treatment arm over time.
        plots_list$arm_mean_ci <-
          ggplot2::ggplot(
            arm_plot_data,
            ggplot2::aes(
              x = time_label,
              y = mean,
              group = arm,
              linetype = arm,
              shape = arm
            )
          ) +
          ggplot2::geom_line(linewidth = 0.8, na.rm = TRUE) +
          ggplot2::geom_point(size = 2.8, na.rm = TRUE) +
          ggplot2::geom_errorbar(
            ggplot2::aes(ymin = ci_lower, ymax = ci_upper),
            width = 0.08,
            na.rm = TRUE
          ) +
          ggplot2::labs(
            x = "Time",
            y = paste("Mean", outcome_display),
            linetype = "Arm",
            shape = "Arm"
          ) +
          paper_theme

        arm_long_plot <- long_plot[
          !is.na(long_plot$arm),
          ,
          drop = FALSE
        ]

        # 11) Distribution by arm and time.
        arm_long_plot$arm <- factor(
          arm_long_plot$arm,
          levels = arm_levels
        )

        arm_colours <- stats::setNames(
          grDevices::hcl.colors(
            n = length(arm_levels),
            palette = "Dark 3"
          ),
          arm_levels
        )

        plots_list$arm_boxplot <-
          ggplot2::ggplot(
            arm_long_plot,
            ggplot2::aes(
              x = time_label,
              y = value,
              group = interaction(time_label, arm)
            )
          ) +
          ggplot2::geom_boxplot(
            ggplot2::aes(
              fill = arm,
              colour = arm
            ),
            position = ggplot2::position_dodge(width = 0.72),
            width = 0.58,
            linewidth = 0.65,
            alpha = 0.45,
            outlier.shape = NA,
            na.rm = TRUE
          ) +
          ggplot2::geom_point(
            ggplot2::aes(colour = arm),
            position = ggplot2::position_jitterdodge(
              jitter.width = 0.08,
              jitter.height = 0,
              dodge.width = 0.72,
              seed = 123
            ),
            shape = 16,
            alpha = 0.42,
            size = 1.45,
            na.rm = TRUE
          ) +
          ggplot2::scale_fill_manual(
            values = arm_colours,
            name = "Arm",
            drop = FALSE
          ) +
          ggplot2::scale_colour_manual(
            values = arm_colours,
            name = "Arm",
            drop = FALSE
          ) +
          ggplot2::scale_y_continuous(
            expand = ggplot2::expansion(mult = c(0.04, 0.08))
          ) +
          ggplot2::labs(
            x = "Time",
            y = outcome_display
          ) +
          ggplot2::guides(
            colour = "none",
            fill = ggplot2::guide_legend(
              nrow = max(
                1L,
                ceiling(length(arm_levels) / 5L)
              ),
              byrow = TRUE,
              override.aes = list(
                alpha = 0.75,
                linewidth = 0.7
              )
            )
          ) +
          paper_theme +
          ggplot2::theme(
            panel.grid.major.x = ggplot2::element_blank(),
            panel.grid.major.y = ggplot2::element_line(
              colour = "grey88",
              linewidth = 0.35
            ),
            axis.line = ggplot2::element_line(
              colour = "grey25",
              linewidth = 0.45
            ),
            axis.ticks = ggplot2::element_line(
              colour = "grey25",
              linewidth = 0.4
            ),
            legend.key = ggplot2::element_blank()
          )

        arm_change_plot_data <- trajectory_summary[
          is.finite(trajectory_summary$absolute_change) &
            !is.na(trajectory_summary$arm),
          ,
          drop = FALSE
        ]

        arm_change_plot_data$arm <- factor(
          arm_change_plot_data$arm,
          levels = arm_levels
        )

        # 12) Baseline-to-final change distribution by arm.
        plots_list$arm_change <-
          ggplot2::ggplot(
            arm_change_plot_data,
            ggplot2::aes(x = arm, y = absolute_change)
          ) +
          ggplot2::geom_boxplot(outlier.shape = NA, na.rm = TRUE) +
          ggplot2::geom_jitter(
            width = 0.08,
            alpha = 0.30,
            na.rm = TRUE
          ) +
          ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
          ggplot2::labs(
            x = "Treatment arm",
            y = paste("Final - Baseline", outcome_display)
          ) +
          paper_theme

        # 13) NEW: mean change from baseline and CI over time by arm.
        if (!is.null(arm_change_descriptives) && nrow(arm_change_descriptives) > 0L) {
          arm_change_ci_data <- arm_change_descriptives
          arm_change_ci_data$to_label <- factor(
            arm_change_ci_data$to_label,
            levels = unname(time_labels[time_vars[-1L]])
          )
          arm_change_ci_data$arm <- factor(
            arm_change_ci_data$arm,
            levels = arm_levels
          )

          plots_list$arm_change_ci <-
            ggplot2::ggplot(
              arm_change_ci_data,
              ggplot2::aes(
                x = to_label,
                y = mean_change,
                group = arm,
                linetype = arm,
                shape = arm
              )
            ) +
            ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
            ggplot2::geom_line(linewidth = 0.8, na.rm = TRUE) +
            ggplot2::geom_point(size = 2.8, na.rm = TRUE) +
            ggplot2::geom_errorbar(
              ggplot2::aes(ymin = ci_lower, ymax = ci_upper),
              width = 0.08, na.rm = TRUE
            ) +
            ggplot2::labs(
              x = "Time",
              y = paste("Mean change from baseline (", ci_text, ")", sep = ""),
              linetype = "Arm",
              shape = "Arm"
            ) +
            paper_theme
        }

        # 14) NEW: between-arm mean differences with CI at each timepoint.
        if (!is.null(arm_time_pairwise) && nrow(arm_time_pairwise) > 0L) {
          arm_difference_data <- arm_time_pairwise
          arm_difference_data$comparison <- paste(
            arm_difference_data$arm_b,
            "-",
            arm_difference_data$arm_a
          )
          arm_difference_data$time_label <- factor(
            arm_difference_data$time_label,
            levels = rev(unname(time_labels[time_vars]))
          )

          plots_list$arm_difference_ci <-
            ggplot2::ggplot(
              arm_difference_data,
              ggplot2::aes(
                x = time_label,
                y = mean_difference_b_minus_a
              )
            ) +
            ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
            ggplot2::geom_errorbar(
              ggplot2::aes(ymin = ci_lower, ymax = ci_upper),
              width = 0.12, na.rm = TRUE
            ) +
            ggplot2::geom_point(size = 2.7, na.rm = TRUE) +
            ggplot2::coord_flip() +
            ggplot2::facet_wrap(~comparison, scales = "free_y") +
            ggplot2::labs(
              x = NULL,
              y = paste("Mean difference (", ci_text, ")", sep = "")
            ) +
            paper_theme
        }

        # 15) NEW: unavailable observations (%) by arm and time.
        arm_missing_plot_list <- list()
        am_counter <- 1L
        for (g in arm_levels) {
          for (k in seq_along(time_vars)) {
            v <- time_vars[[k]]
            idx <- !is.na(analysis_data[[arm]]) & as.character(analysis_data[[arm]]) == g
            n_group <- sum(idx)
            unavailable_n <- sum(is.na(analysis_data[[v]][idx]))
            arm_missing_plot_list[[am_counter]] <- data.frame(
              arm = g,
              time_label = unname(time_labels[v]),
              unavailable_pct = if (n_group > 0L) unavailable_n / n_group * 100 else NA_real_,
              stringsAsFactors = FALSE
            )
            am_counter <- am_counter + 1L
          }
        }
        arm_missing_plot_data <- do.call(rbind, arm_missing_plot_list)
        arm_missing_plot_data$arm <- factor(arm_missing_plot_data$arm, levels = arm_levels)
        arm_missing_plot_data$time_label <- factor(
          arm_missing_plot_data$time_label,
          levels = unname(time_labels[time_vars])
        )

        plots_list$arm_missingness <-
          ggplot2::ggplot(
            arm_missing_plot_data,
            ggplot2::aes(
              x = time_label,
              y = unavailable_pct,
              group = arm,
              linetype = arm,
              shape = arm
            )
          ) +
          ggplot2::geom_line(linewidth = 0.8, na.rm = TRUE) +
          ggplot2::geom_point(size = 2.7, na.rm = TRUE) +
          ggplot2::labs(
            x = "Time",
            y = "Unavailable observations (%)",
            linetype = "Arm",
            shape = "Arm"
          ) +
          paper_theme
      }
    }
  }

  # ----------------------------------------------------------
  # 16. FINAL RESULT
  # ----------------------------------------------------------

  overview <- list(
    outcome = outcome_name,
    outcome_display = outcome_display,
    id = id,
    n_rows = n_rows,
    n_patients = n_patients,
    n_timepoints = length(time_vars),
    timepoints = time_vars,
    time_labels = time_labels,
    supplied_time_vars = supplied_time_vars,
    missing_ids = missing_id_n,
    duplicated_ids = duplicated_id_n,
    complete_profiles = complete_profiles,
    complete_profiles_pct = complete_profiles / n_rows * 100,
    non_finite_values = total_non_finite,
    covariates = .covariates,
    arm_variable = if (arm_available) arm else NULL,
    reference_arm = if (arm_available) reference_arm else NULL,
    arm_levels = arm_levels,
    arm_counts = arm_counts,
    missing_arm = missing_arm_n,
    arm_analysis = arm_available && arm_tests && length(arm_levels) >= 2L
  )

  result <- list(
    call = match.call(),
    version = "4.0.0",
    settings = list(
      alpha = alpha,
      confidence_level = confidence_level,
      confidence_percent = confidence_percent,
      p_adjust_method = p_adjust_method,
      improvement_direction = improvement_direction,
      stable_threshold = stable_threshold,
      covariates = .covariates,
      arm_variable = if (arm_available) arm else NULL,
      reference_arm = if (arm_available) reference_arm else NULL,
      arm_tests = arm_tests,
      strict_id = strict_id,
      non_finite_handling = "Inf/-Inf treated as NA for analysis"
    ),
    overview = overview,
    outcome = outcome_name,
    outcome_display = outcome_display,
    time_vars = time_vars,
    time_labels = time_labels,
    descriptives = descriptives,
    missing = list(by_time = missing_summary, by_patient = missing_by_patient),
    change = change,
    arm_analysis = list(
      enabled = arm_available && arm_tests && length(arm_levels) >= 2L,
      arm_variable = if (arm_available) arm else NULL,
      reference_arm = if (arm_available) reference_arm else NULL,
      levels = arm_levels,
      counts = arm_counts,
      descriptives = arm_descriptives,
      baseline_balance = baseline_balance,
      missingness = arm_missingness,
      time_omnibus = arm_time_omnibus,
      time_pairwise = arm_time_pairwise,
      change_descriptives = arm_change_descriptives,
      change_omnibus = arm_change_omnibus,
      change_pairwise = arm_change_pairwise
    ),
    correlations = list(
      pearson = correlation_pearson,
      spearman = correlation_spearman,
      pairwise_n = correlation_n
    ),
    variability = variability,
    trajectories = trajectory_summary,
    model = list(
      fitted_model = mixed_model,
      summary = model_summary,
      anova = model_anova,
      global_time_test = global_time_test,
      global_arm_test = global_arm_test,
      arm_time_interaction_test = arm_time_interaction_test,
      singular = model_singular,
      converged = model_converged,
      warnings = model_warnings,
      error = model_error,
      covariates_requested = model_covariates_requested,
      covariates_used = model_covariates_used,
      covariates_skipped = model_covariates_skipped,
      fixed_parameters = model_fixed_parameters
    ),
    advanced_tests = advanced_results$advanced_tests,
    advanced_models = advanced_results$advanced_models,
    robustness = advanced_results$robustness,
    effect_sizes = advanced_results$effect_sizes,
    multiplicity = advanced_results$multiplicity,
    sensitivity = advanced_results$sensitivity,
    outliers = list(
      note = "IQR flags are diagnostic flags and are not automatically excluded from analyses.",
      by_time = if (outliers) outlier_results else NULL,
      change = if (outliers) change_outliers else NULL
    ),
    plots = plots_list,
    plot_error = plot_error,
    long_data = long_data
  )

  class(result) <- c("mira_info", "list")

  if (verbose) print(result)
  invisible(result)
}


# ============================================================
# DATA-DRIVEN CONFIGURATION LAYER
# ============================================================

.mira_name_groups <- function(groups) {
  if (length(groups) == 0L) return(groups)
  labels <- vapply(groups, function(x) x$outcome, character(1L))
  labels <- make.unique(labels, sep = "_")
  for (i in seq_along(groups)) groups[[i]]$outcome <- labels[[i]]
  names(groups) <- labels
  groups
}

.mira_match_group <- function(request, groups) {
  if (length(groups) == 0L) return(NA_integer_)
  hit <- which(.mira_key(names(groups)) == .mira_key(request))
  if (length(hit) == 1L) hit else NA_integer_
}

.mira_groups_from_time_vars <- function(data,
                                        time_vars,
                                        outcomes = NULL,
                                        variable_pattern = "auto") {
  groups <- list()

  if (is.list(time_vars) && !is.data.frame(time_vars)) {
    if (length(time_vars) == 0L) stop("time_vars cannot be an empty list.", call. = FALSE)
    supplied_names <- names(time_vars)
    has_names <- !is.null(supplied_names) && all(nzchar(supplied_names))
    outcome_hints <- outcomes
    if (!is.null(outcome_hints) && length(outcome_hints) == 1L &&
        identical(outcome_hints, "auto")) outcome_hints <- NULL

    if (!has_names && !is.null(outcome_hints) && length(outcome_hints) != length(time_vars)) {
      stop(
        "An unnamed time_vars list requires one outcome per element or explicit element names.",
        call. = FALSE
      )
    }

    for (i in seq_along(time_vars)) {
      hint <- if (has_names) supplied_names[[i]] else {
        if (!is.null(outcome_hints)) outcome_hints[[i]] else NULL
      }
      groups[[i]] <- .mira_make_longitudinal_group(
        data = data,
        variables = time_vars[[i]],
        outcome = hint,
        variable_pattern = variable_pattern,
        automatic = FALSE
      )
    }
    return(.mira_name_groups(groups))
  }

  if (!is.character(time_vars) || length(time_vars) < 2L || anyNA(time_vars)) {
    stop("time_vars must be a vector containing at least two names or a list by outcome.",
         call. = FALSE)
  }
  missing <- setdiff(time_vars, names(data))
  if (length(missing) > 0L) {
    stop(sprintf("Longitudinal variables not found: %s.", paste(missing, collapse = ", ")),
         call. = FALSE)
  }

  parsed <- lapply(time_vars, .mira_parse_longitudinal_name,
                   variable_pattern = variable_pattern)
  parsed_outcomes <- vapply(parsed, function(x) {
    if (is.null(x)) NA_character_ else x$outcome
  }, character(1L))
  parsed_keys <- .mira_compact_unique(.mira_key(parsed_outcomes))

  if (all(!is.na(parsed_outcomes)) && length(parsed_keys) > 1L) {
    for (key in parsed_keys) {
      idx <- which(.mira_key(parsed_outcomes) == key)
      groups[[length(groups) + 1L]] <- .mira_make_longitudinal_group(
        data = data,
        variables = time_vars[idx],
        outcome = parsed_outcomes[idx][[1L]],
        variable_pattern = variable_pattern,
        automatic = FALSE
      )
    }
  } else {
    hint <- NULL
    if (!is.null(outcomes) && length(outcomes) == 1L && !identical(outcomes, "auto") &&
        !outcomes %in% names(data)) hint <- outcomes
    groups[[1L]] <- .mira_make_longitudinal_group(
      data = data,
      variables = time_vars,
      outcome = hint,
      variable_pattern = variable_pattern,
      automatic = FALSE
    )
  }
  .mira_name_groups(groups)
}

.mira_apply_time_labels <- function(groups, time_labels) {
  if (is.null(time_labels)) return(groups)

  if (is.list(time_labels) && !is.data.frame(time_labels)) {
    label_names <- names(time_labels)
    if (length(groups) > 1L && (is.null(label_names) || any(!nzchar(label_names)))) {
      stop("For multiple outcomes, time_labels must be a named list.", call. = FALSE)
    }
    for (i in seq_along(groups)) {
      labels <- if (length(groups) == 1L && (is.null(label_names) || !nzchar(label_names[[1L]]))) {
        time_labels[[1L]]
      } else {
        idx <- which(.mira_key(label_names) == .mira_key(names(groups)[[i]]))
        if (length(idx) != 1L) {
          stop(sprintf("time_labels does not contain exactly one element for '%s'.",
                       names(groups)[[i]]), call. = FALSE)
        }
        time_labels[[idx]]
      }
      if (length(labels) != length(groups[[i]]$variables)) {
        stop(sprintf("time_labels for '%s' must contain %d elements.",
                     names(groups)[[i]], length(groups[[i]]$variables)), call. = FALSE)
      }
      if (!is.null(names(labels)) && all(groups[[i]]$variables %in% names(labels))) {
        labels <- labels[groups[[i]]$variables]
      } else if (!isTRUE(groups[[i]]$automatic) && !is.null(groups[[i]]$sort_index)) {
        labels <- labels[groups[[i]]$sort_index]
      }
      groups[[i]]$time_labels <- as.character(labels)
    }
  } else {
    if (!is.atomic(time_labels)) stop("time_labels must be a vector or a list.", call. = FALSE)
    labels <- as.character(time_labels)
    names(labels) <- names(time_labels)
    if (!is.null(names(time_labels)) &&
        all(unique(unlist(lapply(groups, `[[`, "variables"))) %in% names(time_labels))) {
      for (i in seq_along(groups)) {
        groups[[i]]$time_labels <- labels[groups[[i]]$variables]
      }
    } else {
      if (length(groups) != 1L) {
        stop(
          "For multiple outcomes, provide time_labels as a named list or a vector named by column.",
          call. = FALSE
        )
      }
      if (length(labels) != length(groups[[1L]]$variables)) {
        stop("time_labels must have the same length as time_vars.", call. = FALSE)
      }
      if (!isTRUE(groups[[1L]]$automatic) && !is.null(groups[[1L]]$sort_index)) {
        labels <- labels[groups[[1L]]$sort_index]
      }
      groups[[1L]]$time_labels <- labels
    }
  }

  for (i in seq_along(groups)) {
    labels <- groups[[i]]$time_labels
    if (anyNA(labels) || any(!nzchar(trimws(labels))) || anyDuplicated(labels)) {
      stop(sprintf("time_labels for '%s' must contain only non-empty, unique labels.",
                   names(groups)[[i]]), call. = FALSE)
    }
  }
  groups
}

.mira_resolve_outcome_groups <- function(data,
                                         detected,
                                         outcomes = NULL,
                                         time_vars = NULL,
                                         time_labels = NULL,
                                         variable_pattern = "auto") {
  warnings <- character(0)
  outcomes_auto <- is.null(outcomes) ||
    (length(outcomes) == 1L && is.character(outcomes) && identical(outcomes, "auto"))
  outcomes_input_auto <- outcomes_auto

  if (!outcomes_auto && (!is.character(outcomes) || anyNA(outcomes) ||
                         any(!nzchar(outcomes)) || anyDuplicated(.mira_key(outcomes)))) {
    stop("outcomes must be NULL, 'auto', or a vector of unique names.", call. = FALSE)
  }

  # Direct column specification is accepted through outcomes for convenience.
  outcomes_are_columns <- !outcomes_auto && length(outcomes) >= 2L &&
    all(outcomes %in% names(data))
  if (is.null(time_vars) && outcomes_are_columns) {
    time_vars <- outcomes
    outcomes <- NULL
    outcomes_auto <- TRUE
  }

  if (!is.null(time_vars)) {
    groups <- .mira_groups_from_time_vars(
      data, time_vars, outcomes = outcomes, variable_pattern = variable_pattern
    )
  } else {
    groups <- .mira_name_groups(detected$groups)
    if (length(detected$ambiguous_groups) > 0L) {
      warnings <- c(warnings, sprintf(
        paste0("Longitudinal groups with ambiguous timepoints were excluded from automatic detection: %s. ",
               "Specify time_vars and, if needed, time_labels."),
        paste(names(detected$ambiguous_groups), collapse = ", ")
      ))
    }
  }

  if (length(groups) == 0L) {
    stop(
      paste0(
        "No outcome with at least two numeric longitudinal columns was detected. ",
        "Specify time_vars or variable_pattern."
      ),
      call. = FALSE
    )
  }

  if (!outcomes_auto) {
    # A single explicitly supplied manual group can be relabelled by outcomes.
    if (!is.null(time_vars) && length(groups) == 1L && length(outcomes) == 1L) {
      groups[[1L]]$outcome <- outcomes[[1L]]
      names(groups) <- outcomes[[1L]]
    } else {
      selected <- list()
      unavailable <- character(0)
      for (requested in outcomes) {
        idx <- .mira_match_group(requested, groups)
        if (is.na(idx)) {
          ambiguous_idx <- .mira_match_group(requested, detected$ambiguous_groups)
          if (!is.na(ambiguous_idx)) {
            stop(sprintf(
              "Outcome '%s' has ambiguous timepoints; specify time_vars/time_labels manually.",
              requested
            ), call. = FALSE)
          }
          unavailable <- c(unavailable, requested)
        } else {
          selected[[requested]] <- groups[[idx]]
          selected[[requested]]$outcome <- requested
        }
      }
      if (length(unavailable) > 0L) {
        stop(sprintf(
          "Outcomes not detected: %s. Available outcomes: %s.",
          paste(unavailable, collapse = ", "), paste(names(groups), collapse = ", ")
        ), call. = FALSE)
      }
      groups <- selected
    }
  }

  groups <- .mira_apply_time_labels(groups, time_labels)
  ambiguous_manual <- names(groups)[vapply(groups, function(x) isTRUE(x$ambiguous_order), logical(1L))]
  if (!is.null(time_vars) && is.null(time_labels) && length(ambiguous_manual) > 0L) {
    warnings <- c(warnings, sprintf(
      paste0("Timepoint order could not be inferred for %s: the time_vars order was retained. ",
             "Specify time_labels to make the mapping explicit."),
      paste(ambiguous_manual, collapse = ", ")
    ))
  }

  list(
    groups = groups,
    warnings = warnings,
    outcomes_auto = outcomes_input_auto,
    time_vars_manual = !is.null(time_vars)
  )
}

.mira_resolve_direction <- function(value, outcome) {
  if (is.null(value)) value <- "auto"
  if (!is.character(value) || anyNA(value)) {
    stop("improvement_direction must be 'auto', 'higher', 'lower', 'unknown', or a named vector.",
         call. = FALSE)
  }
  selected <- value
  if (!is.null(names(value)) && any(nzchar(names(value)))) {
    if (any(!nzchar(names(value)))) {
      stop("When multiple values are supplied, improvement_direction must be named by outcome.", call. = FALSE)
    }
    idx <- which(.mira_key(names(value)) == .mira_key(outcome))
    if (length(idx) != 1L) {
      stop(sprintf("improvement_direction does not contain exactly one value for '%s'.", outcome),
           call. = FALSE)
    }
    selected <- value[[idx]]
  } else if (length(value) > 1L) {
    stop("When multiple values are supplied, improvement_direction must be named by outcome.", call. = FALSE)
  } else {
    selected <- value[[1L]]
  }
  if (!selected %in% c("auto", "higher", "lower", "unknown")) {
    stop("Allowed values for improvement_direction: auto, higher, lower, unknown.",
         call. = FALSE)
  }
  if (selected != "auto") {
    return(list(value = selected, automatic = FALSE, reason = "user"))
  }

  key <- .mira_key(outcome)
  higher <- c("bcva", "visual_acuity", "va")
  lower <- c("cmt", "crt", "iop")
  if (key %in% higher) {
    list(value = "higher", automatic = TRUE, reason = "internal_outcome_map")
  } else if (key %in% lower) {
    list(value = "lower", automatic = TRUE, reason = "internal_outcome_map")
  } else {
    list(value = "unknown", automatic = TRUE, reason = "no_reliable_mapping")
  }
}

.mira_resolve_threshold <- function(value, outcome) {
  if (is.null(value) || (is.character(value) && length(value) == 1L &&
                         identical(value, "auto"))) {
    return(list(value = 0, automatic = TRUE, reason = "safe_zero_fallback"))
  }
  if (!is.numeric(value) || anyNA(value) || any(!is.finite(value))) {
    stop("stable_threshold must be NULL, 'auto', or a numeric value >= 0.", call. = FALSE)
  }
  selected <- value
  if (!is.null(names(value)) && any(nzchar(names(value)))) {
    if (any(!nzchar(names(value)))) {
      stop("When multiple values are supplied, stable_threshold must be named by outcome.", call. = FALSE)
    }
    idx <- which(.mira_key(names(value)) == .mira_key(outcome))
    if (length(idx) != 1L) {
      stop(sprintf("stable_threshold does not contain exactly one value for '%s'.", outcome),
           call. = FALSE)
    }
    selected <- value[[idx]]
  } else if (length(value) > 1L) {
    stop("When multiple values are supplied, stable_threshold must be named by outcome.", call. = FALSE)
  } else selected <- value[[1L]]
  if (selected < 0) stop("stable_threshold must be >= 0.", call. = FALSE)
  list(value = unname(selected), automatic = FALSE, reason = "user")
}

.mira_choose_reference_arm <- function(data, arm, reference_arm = NULL) {
  if (is.null(arm)) {
    if (!is.null(reference_arm)) {
      stop("reference_arm requires an arm variable.", call. = FALSE)
    }
    return(list(value = NULL, automatic = is.null(reference_arm), warnings = character(0)))
  }

  raw <- as.character(data[[arm]])
  raw[is.na(raw) | !nzchar(trimws(raw))] <- NA_character_
  observed <- unique(raw[!is.na(raw)])
  if (length(observed) == 0L) {
    return(list(value = NULL, automatic = is.null(reference_arm),
                warnings = sprintf("Treatment-arm variable '%s' contains no observed groups.", arm)))
  }

  if (!is.null(reference_arm)) {
    if (!is.character(reference_arm) || length(reference_arm) != 1L ||
        is.na(reference_arm) || !nzchar(reference_arm)) {
      stop("reference_arm must identify exactly one non-empty group.", call. = FALSE)
    }
    if (!reference_arm %in% observed) {
      stop(sprintf("reference_arm='%s' is not present in '%s'. Observed groups: %s.",
                   reference_arm, arm, paste(observed, collapse = ", ")), call. = FALSE)
    }
    return(list(value = reference_arm, automatic = FALSE, warnings = character(0)))
  }

  control_pattern <- "^(control|ctrl|placebo|standard|usual[ _-]?care|comparator|sham|untreated)$"
  control_like <- observed[grepl(control_pattern, observed, ignore.case = TRUE, perl = TRUE)]
  warnings <- character(0)
  if (length(control_like) == 1L) {
    chosen <- control_like[[1L]]
  } else if (length(control_like) > 1L) {
    control_counts <- table(factor(raw, levels = control_like), useNA = "no")
    chosen <- control_like[which.max(as.numeric(control_counts))]
    warnings <- sprintf(
      "Multiple plausible reference arms (%s): the largest group '%s' was selected.",
      paste(control_like, collapse = ", "), chosen
    )
  } else {
    counts <- table(factor(raw, levels = observed), useNA = "no")
    chosen <- observed[which.max(as.numeric(counts))]
  }
  list(value = chosen, automatic = TRUE, warnings = warnings)
}

.mira_build_analysis_config <- function(data,
                                        id = NULL,
                                        outcomes = NULL,
                                        time_vars = NULL,
                                        time_labels = NULL,
                                        arm = NULL,
                                        reference_arm = NULL,
                                        covariates = NULL,
                                        variable_pattern = "auto",
                                        improvement_direction = "auto",
                                        stable_threshold = "auto",
                                        strict_id = TRUE) {
  if (!is.data.frame(data)) stop("data must be a data.frame.", call. = FALSE)
  if (nrow(data) == 0L) stop("The data set contains no observations.", call. = FALSE)
  if (anyDuplicated(names(data))) {
    stop("The data set contains duplicate column names; rename them before analysis.",
         call. = FALSE)
  }
  if (!is.logical(strict_id) || length(strict_id) != 1L || is.na(strict_id)) {
    stop("strict_id must be TRUE or FALSE.", call. = FALSE)
  }

  detected_longitudinal <- .mira_detect_longitudinal_variables(data, variable_pattern)
  resolved <- .mira_resolve_outcome_groups(
    data = data,
    detected = detected_longitudinal,
    outcomes = outcomes,
    time_vars = time_vars,
    time_labels = time_labels,
    variable_pattern = variable_pattern
  )
  groups <- resolved$groups

  parsed_long_vars <- unique(detected_longitudinal$parsed$variable[
    detected_longitudinal$parsed$numeric
  ])
  selected_long_vars <- unique(unlist(lapply(groups, `[[`, "variables"), use.names = FALSE))
  all_long_vars <- unique(c(parsed_long_vars, selected_long_vars))

  id_detection <- .mira_detect_id(data, id = id, exclude = all_long_vars)
  analysis_data <- data
  id_generated <- FALSE
  selected_id <- id_detection$selected
  if (is.null(selected_id)) {
    selected_id <- ".mira_subject_id"
    while (selected_id %in% names(analysis_data)) selected_id <- paste0(selected_id, "_")
    analysis_data[[selected_id]] <- seq_len(nrow(analysis_data))
    id_generated <- TRUE
  }
  id_validation_warnings <- character(0)
  if (!id_generated) {
    id_values <- analysis_data[[selected_id]]
    missing_id_n <- sum(is.na(id_values))
    duplicated_id_n <- sum(duplicated(id_values[!is.na(id_values)]))
    if (strict_id && missing_id_n > 0L) {
      stop(sprintf("The ID variable '%s' contains %d missing values.",
                   selected_id, missing_id_n), call. = FALSE)
    }
    if (strict_id && duplicated_id_n > 0L) {
      stop(sprintf(
        paste0("The ID variable '%s' contains %d duplicate values; wide format requires ",
               "one row per subject. Use strict_id=FALSE only if duplicate IDs are intentional."),
        selected_id, duplicated_id_n
      ), call. = FALSE)
    }
    if (!strict_id && missing_id_n > 0L) {
      id_validation_warnings <- c(id_validation_warnings, sprintf(
        "ID variable '%s' contains %d missing values.", selected_id, missing_id_n
      ))
    }
    if (!strict_id && duplicated_id_n > 0L) {
      id_validation_warnings <- c(id_validation_warnings, sprintf(
        "ID variable '%s' contains %d duplicates.", selected_id, duplicated_id_n
      ))
    }
  }

  arm_detection <- .mira_detect_arm(
    data,
    arm = arm,
    exclude = unique(c(all_long_vars, selected_id))
  )
  selected_arm <- arm_detection$selected
  reference <- .mira_choose_reference_arm(data, selected_arm, reference_arm)

  covariate_detection <- .mira_detect_covariates(
    data,
    covariates = covariates,
    exclude = unique(c(
      all_long_vars, selected_id, selected_arm,
      id_detection$candidates$variable,
      arm_detection$candidates$variable
    ))
  )

  directions <- lapply(names(groups), function(outcome) {
    .mira_resolve_direction(improvement_direction, outcome)
  })
  names(directions) <- names(groups)
  thresholds <- lapply(names(groups), function(outcome) {
    .mira_resolve_threshold(stable_threshold, outcome)
  })
  names(thresholds) <- names(groups)

  warnings <- c(
    resolved$warnings,
    id_detection$warnings,
    id_validation_warnings,
    arm_detection$warnings,
    reference$warnings,
    covariate_detection$warnings
  )
  if (length(detected_longitudinal$non_numeric_matches) > 0L) {
    warnings <- c(warnings, sprintf(
      "Columns with longitudinal names but non-numeric types were excluded: %s.",
      paste(detected_longitudinal$non_numeric_matches, collapse = ", ")
    ))
  }
  warnings <- unique(warnings[!is.na(warnings) & nzchar(warnings)])

  time_vars_config <- lapply(groups, `[[`, "variables")
  time_labels_config <- lapply(groups, function(group) {
    stats::setNames(as.character(group$time_labels), group$variables)
  })

  config <- list(
    id = selected_id,
    id_generated = id_generated,
    outcomes = names(groups),
    time_vars = time_vars_config,
    time_labels = time_labels_config,
    arm = selected_arm,
    reference_arm = reference$value,
    covariates = covariate_detection$selected,
    covariates_detected = covariate_detection$detected,
    covariates_numeric = covariate_detection$numeric,
    covariates_categorical = covariate_detection$categorical,
    improvement_direction = vapply(directions, `[[`, character(1L), "value"),
    stable_threshold = vapply(thresholds, `[[`, numeric(1L), "value"),
    variable_pattern = if (is.function(variable_pattern)) "<custom function>" else variable_pattern,
    auto_detected = list(
      id = is.null(id),
      outcomes = resolved$outcomes_auto,
      time_vars = !resolved$time_vars_manual,
      time_labels = is.null(time_labels),
      arm = is.null(arm),
      reference_arm = is.null(reference_arm),
      covariates = is.null(covariates) || identical(covariates, "auto"),
      improvement_direction = vapply(directions, `[[`, logical(1L), "automatic"),
      stable_threshold = vapply(thresholds, `[[`, logical(1L), "automatic")
    ),
    specified_manually = list(
      id = !is.null(id), outcomes = !resolved$outcomes_auto,
      time_vars = resolved$time_vars_manual, time_labels = !is.null(time_labels),
      arm = !is.null(arm), reference_arm = !is.null(reference_arm),
      covariates = !is.null(covariates) && !identical(covariates, "auto")
    ),
    sources = list(
      id = if (id_generated) "generated_row_id" else id_detection$source,
      arm = arm_detection$source,
      improvement_direction = vapply(directions, `[[`, character(1L), "reason"),
      stable_threshold = vapply(thresholds, `[[`, character(1L), "reason")
    ),
    alternatives = list(
      id = id_detection$alternatives,
      arm = arm_detection$alternatives,
      ambiguous_longitudinal = names(detected_longitudinal$ambiguous_groups)
    )
  )

  detected_variables <- list(
    id_candidates = id_detection$candidates,
    longitudinal = lapply(detected_longitudinal$groups, `[[`, "variables"),
    longitudinal_map = detected_longitudinal$parsed,
    ambiguous_longitudinal = lapply(detected_longitudinal$ambiguous_groups, `[[`, "variables"),
    outcomes = names(detected_longitudinal$groups),
    arm_candidates = arm_detection$candidates,
    covariates_numeric = covariate_detection$detected_numeric,
    covariates_categorical = covariate_detection$detected_categorical
  )

  data_overview <- list(
    n_rows = nrow(data),
    n_columns = ncol(data),
    column_profile = .mira_column_profile(data)
  )

  list(
    data = analysis_data,
    groups = groups,
    config = config,
    data_overview = data_overview,
    detected_variables = detected_variables,
    warnings = warnings
  )
}






#' Inspect longitudinal configuration without running analyses
#'
#' Resolves subject identifiers, longitudinal outcome groups, ordered
#' timepoints, treatment arm and reference level, adjustment covariates,
#' improvement directions, and stability thresholds using the same validation
#' and detection rules as [mira_info()]. No outcome-specific statistical
#' analysis is performed.
#'
#' `mira_detect()` provides the configuration-only workflow directly. The same
#' result can be requested from `mira_info()` with `inspect_only = TRUE`, which
#' additionally resolves the optional-analysis switches accepted by that
#' function.
#'
#' @inheritParams mira_info
#'
#' @return Invisibly returns an object with classes `mira_detect` and `list`.
#'   It contains `call`, the output-schema `version`, resolved `config`, a
#'   `data_overview`, `detected_variables`, and configuration `diagnostics`.
#'   If `verbose = TRUE`, the corresponding detection report is printed first.
#'
#' @seealso [mira_info()]
#'
#' @export
mira_detect <- function(data,
                        id = NULL,
                        outcomes = NULL,
                        time_vars = NULL,
                        time_labels = NULL,
                        arm = NULL,
                        reference_arm = NULL,
                        covariates = NULL,
                        variable_pattern = "auto",
                        improvement_direction = "auto",
                        stable_threshold = "auto",
                        strict_id = TRUE,
                        verbose = TRUE) {
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("verbose must be TRUE or FALSE.", call. = FALSE)
  }
  built <- .mira_build_analysis_config(
    data = data, id = id, outcomes = outcomes, time_vars = time_vars,
    time_labels = time_labels, arm = arm, reference_arm = reference_arm,
    covariates = covariates, variable_pattern = variable_pattern,
    improvement_direction = improvement_direction,
    stable_threshold = stable_threshold, strict_id = strict_id
  )
  for (message in built$warnings) warning(message, call. = FALSE)
  result <- list(
    call = match.call(),
    version = "4.0.0",
    config = built$config,
    data_overview = built$data_overview,
    detected_variables = built$detected_variables,
    diagnostics = list(warnings = built$warnings)
  )
  class(result) <- c("mira_detect", "list")
  if (verbose) print(result)
  invisible(result)
}

#' Analyse wide-format longitudinal and repeated-measures outcomes
#'
#' Performs analyses of one or more numeric longitudinal outcomes stored in wide format,
#' with one measurement column per timepoint. Subject identifiers, outcome
#' groups, time order, treatment arm, and adjustment covariates can be detected
#' from the data or specified explicitly. The function computes descriptive,
#' change, missingness, variability, and trajectory summaries and can add
#' treatment-group comparisons, correlations, diagnostic flags, plots, and a
#' suite of longitudinal models. Results are returned as structured objects for
#' programmatic use; printing is only a presentation layer.
#'
#' @param data A non-empty `data.frame` in wide longitudinal format. Repeated
#'   measurements must be numeric columns, normally with one row per subject.
#'   Column names must be unique. Other columns may contain a subject identifier,
#'   treatment arm, or time-invariant covariates. `Inf` and `-Inf` in selected
#'   longitudinal columns are retained in the supplied data but treated as
#'   unavailable (`NA`) for analysis.
#' @param id `NULL` (the default) or a length-one character string naming the
#'   subject-ID column. An explicit name must exist and takes precedence over
#'   detection. With `NULL`, the function searches semantically ID-like names and
#'   requires the selected column to be complete and unique. If there is no
#'   single safe candidate, it creates an internal row ID named
#'   `.mira_subject_id` (with additional underscores if needed) and records this
#'   in `config$id_generated`.
#' @param time_vars `NULL` (the default), a character vector of at least two
#'   longitudinal column names, or a list of such vectors, one per outcome. All
#'   selected columns must exist, be numeric, and be distinct within each outcome
#'   group. A named list uses its element names as outcome names. For an unnamed
#'   list, `outcomes` may supply one name per element; otherwise names are
#'   inferred from the columns where possible. A character vector whose names
#'   all parse into more than one outcome family is split into those families.
#'   With `NULL`, all unambiguous automatically detected outcome groups are used.
#'   Explicit groups are ordered by parsed time values when every value is
#'   finite and unique; otherwise their supplied order is retained.
#' @param time_labels `NULL` (the default), an atomic vector, or a list of
#'   display labels. Labels must be non-missing, non-empty, unique within an
#'   outcome, and have the same length as the corresponding `time_vars` group.
#'   For multiple outcomes, use a list named by outcome or a single vector named
#'   by every selected longitudinal column. For one outcome, an unnamed vector
#'   or an unnamed one-element list is also accepted. A vector named by source
#'   column is matched by name. Otherwise, when explicit `time_vars` are
#'   reordered from parsed suffixes, the labels are reordered by the same
#'   permutation. `NULL` uses parsed suffixes, or the source column name when
#'   parsing did not supply a label.
#' @param arm `NULL` (the default) or a length-one character string naming a
#'   treatment or grouping column. An explicit name must exist and overrides
#'   detection. With `NULL`, only an unambiguous, plausibly categorical column
#'   with an arm-like name is selected; otherwise between-group analyses remain
#'   disabled. The arm is still carried into long-format and trajectory output
#'   when selected, even if `arm_tests` is false.
#' @param reference_arm `NULL` (the default) or one non-empty character value
#'   identifying an observed level of `arm`. It is an error to supply this
#'   without a resolved arm or to name an unobserved level. With `NULL`, a single
#'   control-like level (`control`, `ctrl`, `placebo`, `standard`, `comparator`,
#'   `sham`, `untreated`, or `usual care` with an optional space, underscore, or
#'   hyphen, all matched case-insensitively) is preferred. If several such levels
#'   exist, the largest is used; if none exists, the largest observed arm is
#'   used. Ties follow the order in which levels occur in `data`. The chosen
#'   level is first in model and comparison level order.
#' @param arm_tests A non-missing logical scalar, defaulting to `TRUE`. If
#'   `TRUE`, and a resolved arm has at least two usable groups, compute
#'   time-by-arm descriptives, missingness and balance diagnostics, between-arm
#'   tests, change comparisons, and arm terms in requested models. It can be
#'   overridden by `analyses`.
#' @param alpha A finite numeric scalar strictly between 0 and 1; the default is
#'   `0.05`. The confidence level is `1 - alpha` for mean, change, Welch,
#'   marginal-mean, and supported model-based confidence intervals. It is not
#'   used to remove results or to define clinical improvement.
#' @param plots A non-missing logical scalar, defaulting to `TRUE`, requesting
#'   stored `ggplot2` plots. Plotting is automatically disabled for an outcome
#'   with no finite values and can be overridden by `analyses`. If `ggplot2` is
#'   unavailable, analysis continues, `plots` is empty, and `plot_error` records
#'   the reason.
#' @param model A non-missing logical scalar, defaulting to `TRUE`, requesting
#'   the random-intercept mixed model and the optional model-based longitudinal
#'   modules. It is automatically disabled for an outcome with fewer than three
#'   finite values in total or no timepoint containing at least two distinct
#'   finite values, and it can be overridden by `analyses`. The base-R Friedman
#'   module is attempted independently of this switch.
#' @param outliers A non-missing logical scalar, defaulting to `TRUE`, requesting
#'   Tukey 1.5-IQR flags for values at each timepoint and for every available
#'   pairwise change. At least four finite values are required for a set of flags.
#'   Flagged observations are diagnostic only and are never removed
#'   automatically. This argument can be overridden by `analyses`.
#' @param correlations A non-missing logical scalar, defaulting to `TRUE`,
#'   requesting Pearson and Spearman correlations between timepoint columns,
#'   using pairwise-complete observations, together with the corresponding
#'   pairwise sample-size matrix. It is automatically disabled unless at least
#'   two timepoints each contain at least two finite observations and at least
#'   two distinct finite values, and can be overridden by `analyses`.
#' @param verbose A non-missing logical scalar, defaulting to `TRUE`. If `TRUE`,
#'   the appropriate `print()` method is called before the result is returned
#'   invisibly. Assigning the call, for example `fit <- mira_info(data)`, still
#'   captures the complete object. If `FALSE`, no automatic report is printed;
#'   warnings are unaffected.
#' @param p_adjust_method A length-one character string in
#'   `stats::p.adjust.methods`, defaulting to `"holm"`, passed to
#'   [stats::p.adjust()] as the primary adjustment. Adjustments are made within
#'   the implemented comparison families, not across every test produced by the
#'   function. Unadjusted advanced contrast summaries are augmented with raw,
#'   primary, Bonferroni, Holm, Benjamini-Hochberg, and Benjamini-Yekutieli
#'   values where estimable.
#' @param improvement_direction `NULL`, `"auto"` (the default), `"higher"`,
#'   `"lower"`, or `"unknown"`; for multiple outcomes it may be a character
#'   vector named by outcome. Each named outcome must match exactly once after
#'   case-and-punctuation normalization; an unnamed vector of length greater
#'   than one is invalid. `"higher"` and `"lower"` determine whether positive
#'   or negative change is labelled improved; `"unknown"` reports only increase,
#'   decrease, and stability.
#'   `"auto"` (also selected by `NULL`) maps `bcva`, `visual_acuity`, and `va` to
#'   higher-is-better and `cmt`, `crt`, and `iop` to lower-is-better after name
#'   normalization; all other outcome names fall back to `"unknown"`. This
#'   setting changes labels and counts, not estimates or hypothesis tests, and
#'   does not validate clinical importance.
#' @param stable_threshold `NULL`, `"auto"` (the default), a finite non-negative
#'   numeric scalar, or a finite non-negative numeric vector named by outcome.
#'   Each named outcome must match exactly once after case-and-punctuation
#'   normalization; an unnamed numeric vector of length greater than one is
#'   invalid. `NULL` and `"auto"` use zero. A change greater than the threshold
#'   is an increase, a change less than its negative is a decrease, and all
#'   intervening values are stable. The threshold affects response
#'   classifications and trajectory labels, not statistical tests; it is not
#'   estimated or clinically validated by the function.
#' @param strict_id A non-missing logical scalar, defaulting to `TRUE`. With
#'   `TRUE`, missing or duplicated values in a selected, non-generated ID column
#'   are fatal because wide input is assumed to contain one row per subject. With
#'   `FALSE`, these conditions generate warnings and duplicated IDs are treated
#'   as records from the same subject by subject-indexed models. Use `FALSE` only
#'   when that interpretation is intentional.
#' @param outcomes `NULL` (the default), `"auto"`, or a character vector of
#'   outcome names, unique after case-and-punctuation normalization. `NULL` and
#'   `"auto"` select all unambiguous detected groups. Explicit names select
#'   detected groups by normalized name, label an explicitly supplied single
#'   `time_vars` group, or provide labels for an unnamed `time_vars` list. As a
#'   compatibility convenience, when `time_vars` is `NULL` and `outcomes`
#'   contains at least two values that are all column names, those values are
#'   treated as `time_vars` rather than as outcome-family names. The
#'   implementation also accepts `character(0)` as an explicit selection and
#'   returns an empty `mira_info_multi` container.
#' @param covariates `NULL` (the default), `"auto"`, or a character vector of
#'   column names; `character(0)` explicitly selects none. Explicit columns must
#'   exist, be unique, and not overlap selected or detected longitudinal
#'   variables, the chosen ID or arm, retained ID/arm candidates, or reserved
#'   internal columns (`patient`, `outcome`, `time`, `time_index`, `time_label`,
#'   `arm`, `value`, `patient_factor`, `time_factor`, and `arm_factor`). With
#'   `NULL` or `"auto"`, eligible numeric and low-cardinality categorical columns
#'   are detected after excluding ID-, arm-, and longitudinal-like variables.
#'   Every candidate needs at least two distinct non-missing values; factor,
#'   character, and logical candidates may have at most
#'   `max(10, floor(nrow(data) / 4))` such values,
#'   whereas numeric candidates have no corresponding upper limit. Numeric
#'   columns with at most five distinct observed values are classified as
#'   categorical. The fixed-effect cost is one degree of freedom for a numeric
#'   covariate and `number of observed levels - 1` for a categorical covariate.
#'   Automatic covariates are used only when the combined cost does not exceed
#'   `floor(nrow(data) / 10)`; otherwise none are selected and a warning records
#'   the candidates.
#' @param analyses `NULL` (the default) or a character vector selecting optional
#'   modules from `"plots"`, `"model"`, `"outliers"`, `"correlations"`, and
#'   `"arm_tests"`. With `NULL`, the five individual logical arguments retain
#'   their supplied values. Matching is case-insensitive and duplicate names are
#'   removed. When non-NULL, this vector completely replaces the five
#'   corresponding logical switches: omitted modules are set to false. `"all"`
#'   selects every optional module; `"none"` or `character(0)` selects none. The
#'   implementation expands `"all"` before checking `"none"`; consequently
#'   `c("all", "none")` selects all, whereas any `"none"` without `"all"` clears
#'   the selection. Remaining unknown names are errors. The logical arguments
#'   are still validated before this override. Core summaries and the Friedman
#'   branch are not disabled by `analyses = "none"`.
#' @param variable_pattern Exactly `"auto"` (the default), a non-empty
#'   regular-expression string, or a parser function called once per column
#'   name. With `"auto"`, the supported suffix rules are described under
#'   **Input data and longitudinal naming**. A custom regular-expression match
#'   is recognized only when it yields at least two capture groups: the first is
#'   the outcome name and the second is the time label; time order is derived
#'   from that label. A custom function must return `NULL` for a non-match or a
#'   list (a one-row data frame is also accepted) containing `outcome`, with
#'   optional `time_label` and `time_order`. A missing `time_label` defaults to
#'   the column name, and a missing `time_order` is derived from the label.
#'   Custom regular expressions use the default case-sensitive, non-Perl behavior
#'   of `regexec()`. The first value of each parser field is used, with
#'   `time_order` coerced to numeric.
#' @param inspect_only A non-missing logical scalar, defaulting to `FALSE`. If
#'   `TRUE`, perform input validation and configuration resolution, emit
#'   configuration warnings, and return a `mira_detect` object without running
#'   any outcome analysis. The resolved `analyses` settings are included in
#'   `config`. `verbose = TRUE` prints the detection report.
#'
#' @details
#' `mira_info()` first validates global arguments and resolves the analysis
#' configuration. It groups repeated-measure columns by parsed outcome name,
#' orders their timepoints, resolves or generates an ID, selects an arm and
#' reference level when possible, and selects covariates. Explicit `id`,
#' `time_vars`, `time_labels`, `arm`, `reference_arm`, `outcomes`, and
#' `covariates` values take precedence over inference. The resolved decisions,
#' alternatives, and automatic/manual provenance are retained in `config` and
#' `detected_variables`.
#'
#' Consequently, `mira_info(data)` selects every unambiguous detected outcome,
#' chooses or generates an ID, attempts arm and covariate detection, uses the
#' outcome-name map for improvement direction with a zero stability threshold,
#' requests all five optional modules, and prints the resulting report. Every
#' requested module remains conditional on usable data and its dependency.
#'
#' Unless `inspect_only = TRUE`, each outcome is analysed separately. The wide
#' measurements are converted to an outcome-specific long table while the wide
#' copy is used for timepoint summaries and pairwise changes. Optional analyses
#' are adapted per outcome: a requested model, correlation matrix, plot suite,
#' or arm analysis can be disabled when the outcome lacks the minimum usable
#' data. These adaptations are recorded under `diagnostics$adaptation`.
#'
#' For a single outcome, an analysis error stops the call. With multiple
#' outcomes, each outcome is protected separately: successful outcomes remain
#' available, failures become `mira_info_error` entries, and one warning names
#' the failed outcomes. Missing optional packages and most advanced-model
#' failures instead produce structured skipped modules with `performed = FALSE`,
#' an error or `reason_skipped`, and any captured warnings.
#'
#' @section Input data and longitudinal naming:
#' Automatic parsing requires a separator (`.`, `_`, or `-`) between an outcome
#' stem and a recognized suffix. Matching is case-insensitive for the automatic
#' rules. The implemented suffixes are:
#'
#' - `tN`;
#' - `timeN`, `visitN`, `visN`, or `vN`, with an optional separator before `N`;
#' - `baseline`, `base`, or `BL`;
#' - `monthN`, `monthsN`, `moN`, or `mN`;
#' - `weekN`, `weeksN`, `wkN`, or `wN`;
#' - `dayN`, `daysN`, or `dN`;
#' - `followupN`, `follow-upN`, or `fuN`; and
#' - a bare numeric suffix `N`.
#'
#' Here `N` is a sequence of digits. The time/visit, month, week, day, and
#' follow-up forms also allow an optional `.`, `_`, or `-` immediately before
#' `N`; `tN` requires the letter and digits to be adjacent. Examples include
#' `score_t0`, `score_BL`, `score_M6`, `score_visit_2`, and `score_12`. Baseline
#' labels receive order zero; otherwise the first number in the suffix determines
#' order. At least two numeric columns are required per automatically detected
#' group. A group is excluded from automatic analysis when time labels repeat
#' (case-insensitively) or finite time orders repeat. It remains visible under
#' `detected_variables$ambiguous_longitudinal` and can be supplied explicitly.
#' Columns that match a longitudinal name rule but are not numeric are excluded
#' and reported in a warning and in the parsed longitudinal map.
#'
#' @section Automatic detection and explicit overrides:
#' ID detection ranks complete, unique columns with names such as `patient_id`,
#' `subject_id`, `participant_id`, `usubjid`, `subjid`, `patid`, `record_id`,
#' `study_id`, `patient`, `subject`, `participant`, `id`, or another name ending
#' in `_id`. Structural uniqueness alone is not enough to promote an arbitrary
#' column. Ambiguous or unsafe detection therefore creates a row ID rather than
#' guessing.
#'
#' Arm detection considers arm-like names including `arm`, `treatment_arm`,
#' `treatment_group`, randomized/randomised group, `treatment`, `intervention`,
#' `group`, and `cohort`. A candidate must have from two through
#' `max(10, floor(nrow(data) / 4))` distinct non-missing values and a unique
#' highest semantic rank. Character, factor, and logical columns qualify
#' throughout that range; other column types qualify only with at most ten
#' distinct non-missing values. Ambiguity leaves the arm unset. Detection
#' warnings report excluded or alternative candidates without changing explicit
#' choices.
#'
#' @section Core and optional analyses:
#' The following outcome-specific results are computed independently of the
#' optional switches:
#'
#' - per-timepoint counts, missing and non-finite counts, mean, standard
#'   deviation, variance, standard error, Student-t confidence interval, median,
#'   quartiles, range, and coefficient of variation (`100 * SD / abs(mean)`,
#'   unavailable when the mean is zero);
#' - all ordered timepoint-pair changes (`to - from`) using subjects observed at
#'   both timepoints, with a mean-change t interval, Cohen's `dz`, two-sided
#'   paired t test, asymptotic paired Wilcoxon signed-rank test (`exact = FALSE`),
#'   and direction/response counts;
#' - availability summaries; baseline-to-final trajectories with relative change
#'   `100 * (final - baseline) / abs(baseline)` when baseline is finite and
#'   nonzero; pooled within-subject residual standard deviation; between-subject
#'   standard deviation of subject means; and a complete-profile one-way
#'   ICC(1,1);
#' - a Friedman test on complete repeated profiles, Kendall's W, and all
#'   asymptotic pairwise Wilcoxon signed-rank post-hoc tests with paired
#'   rank-biserial correlations computed after removing zero changes. This
#'   branch is stored under `advanced_tests$friedman` and is attempted even when
#'   `model = FALSE`.
#'
#' The complete-profile ANOVA ICC is
#' `(MS_between - MS_within) / (MS_between + (k - 1) * MS_within)`, where `k` is
#' the number of selected timepoints. The model ICC, when available, is the
#' random-intercept variance divided by the sum of random-intercept and residual
#' variances; `variability$ICC` prefers this estimate and otherwise uses the
#' ANOVA estimate.
#'
#' Requested correlations, IQR flags, plots, arm analyses, and model-based
#' methods are added conditionally. The primary mixed model, when estimable and
#' `lme4` is installed, is a Gaussian random-intercept model fitted by REML with
#' categorical time, additive selected covariates, and, when enabled, the full
#' time-by-arm interaction. `lmerTest::lmer()` is used when `lmerTest` is
#' installed; otherwise `lme4::lmer()` is used. Separate maximum-likelihood
#' nested-model comparisons test time and, where applicable, arm and the
#' time-by-arm interaction. The output also reports convergence messages,
#' singularity, variance components through the fitted object, and a
#' random-intercept/residual ICC.
#'
#' When `model` is enabled, the function conditionally attempts:
#'
#' - estimated marginal means and pairwise, baseline-versus-follow-up,
#'   baseline-versus-final, consecutive-time, arm, simple, and interaction
#'   contrasts with `emmeans`, plus linear, quadratic, and, when available,
#'   cubic ordered polynomial contrasts; Tukey summaries are retained for
#'   ordinary pairwise contrasts and a `dunnettx` summary for baseline-versus-
#'   follow-up contrasts, labelled in the result as an approximation. Polynomial
#'   trends use timepoint order and do not assert equal chronological spacing;
#' - type-III repeated-measures ANOVA with `afex`, including uncorrected,
#'   Greenhouse-Geisser, and Huynh-Feldt tables, Mauchly information when
#'   available, and partial/generalized eta-squared output;
#' - Gaussian GEE fits with `geepack` under independence, exchangeable, and
#'   AR(1) working correlations, sandwich standard errors, and joint Wald tests;
#' - a random-intercept-and-ordinal-time-slope model with `lme4`, including
#'   likelihood-ratio sensitivity tests and random-effect diagnostics; because
#'   the slope uses the ordered index, it does not represent actual unequal
#'   chronological spacing;
#' - REML random-intercept `nlme` fits with compound-symmetry and AR(1) residual
#'   correlations, with AR(1) indexed by ordered timepoint number; and
#' - `clubSandwich` CR2 covariance estimates with Satterthwaite coefficient
#'   tests and HTZ joint Wald tests.
#'
#' For GEE and CR2 joint tests, `TIME` selects every fixed-effect coefficient
#' whose name contains the time factor (including interaction coefficients),
#' `ARM` analogously selects every arm coefficient, and `TIME_X_ARM` selects only
#' coefficients containing both. Thus these are joint coefficient tests, not
#' comparisons of a single reference-cell coefficient.
#'
#' The mixed models specify a conditional mean with subject-specific random
#' effects, whereas GEE describes a population-averaged marginal mean under the
#' stated working correlation. These are different inferential targets even when
#' their numerical estimates are similar.
#'
#' A model-comparison table reports likelihood criteria for likelihood-based
#' models and QIC/CIC for GEE without treating those criteria as a single common
#' scale. `performance`, when installed, supplies Nakagawa marginal and
#' conditional R-squared values for available mixed models. The sensitivity
#' table aligns raw tests of time, arm, and interaction across methods; it is a
#' comparison aid, not a voting rule for scientific conclusions, and the rows
#' can use different analysis samples or null hypotheses.
#'
#' @section Treatment-arm analyses:
#' Arm analyses require a resolved arm, `arm_tests = TRUE`, at least two observed
#' groups, and usable outcome values. Empty strings and missing arm values are
#' excluded from between-group work but do not remove subjects from analyses
#' that do not use arm. The first ordered timepoint is treated as baseline.
#'
#' The arm module provides time-by-arm descriptives; differential-unavailability
#' tests using Pearson chi-squared tests or Fisher's exact tests when any expected
#' count is below five; Welch and Kruskal-Wallis omnibus tests at each timepoint;
#' and two-sided pairwise Welch t and asymptotic unpaired Wilcoxon tests
#' (`exact = FALSE`). Pairwise tables report `arm_b - arm_a`, a Welch interval,
#' and Hedges' g. `baseline_balance` is the first-timepoint subset of these
#' exploratory pairwise comparisons. For every follow-up, analogous descriptives
#' and between-arm tests are computed for
#' change from baseline. These comparisons are associational summaries and do
#' not by themselves establish a causal treatment effect.
#'
#' @section Covariate adjustment:
#' Selected covariates enter model formulas as additive fixed effects. Core
#' descriptive and pairwise-change tables remain unadjusted; adjustment applies
#' only to the model-based results. In the repeated-measures ANOVA, categorical
#' covariates are between-subject factors and numeric covariates are supplied as
#' covariates. The same resolved covariates are used, when possible, by the
#' primary mixed model, marginal means, random-slope model, GEE, and `nlme` fits.
#' Character and logical values are converted to factors; non-finite numeric
#' covariate values and blank character values become missing, and covariates
#' with fewer than two usable
#' values are omitted from the primary mixed model and recorded in
#' `model$covariates_skipped`. Adjustment is statistical and does not imply that
#' the selected variables suffice for causal control.
#'
#' @section Missing data:
#' The function performs no imputation. Descriptives use all finite values at
#' each timepoint. Output distinguishes raw `NA` counts (`missing`), `Inf` or
#' `-Inf` counts (`non_finite`), and their combined analysis-unavailable count.
#' Change estimates and paired tests use a separate pairwise-complete subset for
#' each comparison; correlations also use pairwise-complete observations.
#' Baseline-to-final trajectories require both endpoints, and arm change
#' comparisons require baseline and the relevant
#' follow-up. The complete-profile ANOVA ICC and Friedman test require every
#' selected timepoint. The `afex` analysis receives the prepared repeated-
#' measures data and applies its own repeated-measures completeness requirements.
#' Mixed, GEE, and residual-correlation models use long rows complete for
#' outcome, ID, required arm, and retained covariates, so their analysis samples
#' can differ from the descriptive and pairwise samples. Module-specific sample
#' sizes and prepared data are retained where available.
#'
#' @section Multiplicity:
#' The primary `p_adjust_method` is applied separately to paired t-test and
#' paired Wilcoxon families across timepoint pairs. Within the arm module it is
#' applied separately to missingness tests across time, Welch and Kruskal-Wallis
#' omnibus families, and Welch and Wilcoxon pairwise families for observed values
#' and changes. Advanced contrast and Friedman post-hoc tables define additional
#' scientifically distinct families. Their tables preserve raw p-values and
#' store primary, Bonferroni, Holm, BH, and BY adjustments. Global
#' likelihood-ratio, ANOVA, GEE, and sensitivity-table p-values are not pooled
#' into a single adjustment family.
#'
#' @section Optional dependencies and failure handling:
#' Base summaries and rank tests use R's `stats` package. Optional capabilities
#' require `ggplot2` (plots), `lme4` (mixed models), `lmerTest` (used in place of
#' `lme4::lmer()` when present), `emmeans` (marginal means and contrasts), `afex`
#' (repeated-measures ANOVA), `geepack` (GEE), `nlme` (residual correlation
#' models), `clubSandwich` (CR2 inference), and `performance` (Nakagawa
#' R-squared). Except for plotting, each advanced module records `performed`,
#' package, warnings, error, and/or `reason_skipped` information. Absence of one
#' package does not prevent independent modules from being attempted. Absence of
#' `ggplot2` is recorded in `plot_error`. Captured model-fitting warnings are
#' diagnostics, not evidence that a fit is scientifically adequate.
#'
#' @section Warnings and errors:
#' Invalid argument types or values, missing named columns, duplicated data-frame
#' column names, duplicate or non-numeric time variables, invalid labels,
#' unresolved requested outcomes, and strict-ID violations stop configuration.
#' Automatic ambiguity, generated IDs, excluded detection candidates,
#' non-finite outcome values, missing arm values, and conservative covariate
#' suppression generate warnings. Fewer than two observed arm groups disable
#' arm analysis with a warning. Most unavailable or failed optional modules are
#' represented in the result rather than stopping it. An unexpected analysis
#' error is fatal for a single outcome but is isolated as a `mira_info_error`
#' entry during a multi-outcome call.
#'
#' @return Invisibly returns an object whose class and structure depend on the
#'   requested mode.
#'
#' With `inspect_only = TRUE`, the result has classes `mira_detect` and `list`
#' and contains:
#'
#' \describe{
#'   \item{call, version}{The matched call and output-schema version, currently
#'     `"4.0.0"`.}
#'   \item{config}{Resolved ID, outcome groups, ordered time variables and
#'     labels, arm/reference, covariates, improvement directions, thresholds,
#'     optional-analysis switches, provenance, sources, and alternatives.}
#'   \item{data_overview}{Row and column counts plus a column profile containing
#'     class, numeric status, distinct-value count, and missingness.}
#'   \item{detected_variables}{ID and arm candidates, the parsed longitudinal
#'     map, accepted and ambiguous groups, and detected covariate classes.}
#'   \item{diagnostics}{Configuration and detection warnings.}
#' }
#'
#' A successful single-outcome analysis has classes `mira_info` and `list`.
#' Its important top-level components are:
#'
#' \describe{
#'   \item{call, version}{The matched public call and output-schema version,
#'     currently `"4.0.0"`.}
#'   \item{settings, overview}{Outcome-specific settings and sample,
#'     completeness, timepoint, covariate, and arm summaries.}
#'   \item{outcome, outcome_display, time_vars, time_labels}{Outcome identifiers
#'     (including an upper-case display form) and the ordered source-to-display
#'     timepoint map.}
#'   \item{descriptives}{One row per timepoint with availability, location,
#'     dispersion, and mean-confidence-interval statistics.}
#'   \item{missing}{Lists `by_time` and `by_patient` availability tables; the
#'     latter has one row per wide input row and carries the resolved ID value.}
#'   \item{change}{All estimable ordered timepoint-pair changes, confidence
#'     intervals, Cohen's `dz`, direction counts, raw tests, and adjusted tests.}
#'   \item{arm_analysis}{An `enabled` flag, arm metadata, and conditional
#'     descriptives, baseline balance, missingness, timepoint comparisons, and
#'     baseline-to-follow-up change comparisons. Disabled results are `NULL`.}
#'   \item{correlations}{Pearson and Spearman matrices and `pairwise_n`; values
#'     are `NULL` when correlations were disabled or adapted away.}
#'   \item{variability}{Within- and between-subject standard deviations,
#'     complete-profile ANOVA ICC, model ICC, the preferred available `ICC`, and
#'     ANOVA mean-square components.}
#'   \item{trajectories}{Baseline, final, absolute and relative change, numeric
#'     direction, clinical direction, and arm when available, with one row per
#'     wide input row.}
#'   \item{model}{The primary fitted model and summary, fixed-effect ANOVA,
#'     global likelihood-ratio tests, convergence/singularity diagnostics,
#'     warnings/error, and requested, used, or skipped covariates.}
#'   \item{advanced_tests}{Structured `emmeans`, repeated-measures ANOVA, and
#'     Friedman modules.}
#'   \item{advanced_models}{The random-slope, `nlme`, and GEE modules plus a
#'     model-comparison table and paths to available fit objects.}
#'   \item{robustness}{The `club_sandwich` CR2 module.}
#'   \item{effect_sizes}{Paired Cohen's `dz`, paired rank-biserial correlations,
#'     ANOVA eta-squared measures, Kendall's W, mixed-model R-squared modules,
#'     and model ICC where available.}
#'   \item{multiplicity}{Contrast-family tables, available adjustments, and
#'     source-object paths.}
#'   \item{sensitivity}{A method-by-scientific-question table for time, arm, and
#'     interaction tests, with raw p-values and source paths.}
#'   \item{outliers}{A note and optional IQR-flag tables `by_time` and `change`.}
#'   \item{plots, plot_error}{Named `ggplot` objects and an optional dependency
#'     message. Possible names are `boxplot`, `spaghetti`, `mean_ci`, `change`,
#'     `change_from_baseline`, `change_ci`, `missingness`,
#'     `correlation_heatmap`, `response`, `arm_mean_ci`, `arm_boxplot`,
#'     `arm_change`, `arm_change_ci`, `arm_difference_ci`, and
#'     `arm_missingness`, subject to enabled modules and available data.}
#'   \item{long_data}{Outcome-specific long-format data with subject, source
#'     time variable, ordered index and label, arm, value, and selected
#'     covariates.}
#'   \item{config, data_overview, detected_variables, diagnostics}{Resolved
#'     configuration, detection evidence, and per-outcome adaptations.}
#'   \item{outcomes}{A one-element outcome-indexed list, permitting the same
#'     access convention as multi-outcome results.}
#' }
#'
#' A multi-outcome result has classes `mira_info_multi`, `mira_info`, and `list`.
#' Its top level contains `call`, `version`, `config`, `data_overview`,
#' `detected_variables`, `diagnostics`, and the outcome-indexed `outcomes` list.
#' Each successful element of `outcomes` is an outcome-specific `mira_info`
#' object with the components above. A failed element has classes
#' `mira_info_error` and `list` and contains `outcome` and `error`. Use
#' `names(result$outcomes)` and `result$outcomes[[name]]` for programmatic access.
#'
#' All modes are returned invisibly. `print()` provides detection, compact
#' multi-outcome, error, or full single-outcome reports as appropriate;
#' `summary()` has methods for single- and multi-outcome results. `plot()` prints
#' and invisibly returns a stored plot; for a multi-outcome object, specify
#' `outcome` when more than one successful outcome is available. Requesting a
#' plot that was not generated is an error.
#'
#' @seealso [mira_detect()] for configuration inspection without analysis,
#'   [base::summary()], [graphics::plot()], [stats::p.adjust()],
#'   [stats::friedman.test()], and [stats::wilcox.test()].
#'
#' @examples
#' # ------------------------------------------------------------------
#' # Example 1: Minimal automatic workflow
#' # ------------------------------------------------------------------
#' # The input is wide: one row per subject and one numeric column per visit.
#' # The names subject_id and score_t0, score_t1, score_t2 are deliberately
#' # conventional, so ID, outcome family, and time order can all be inferred.
#' set.seed(202601)
#' n_wide <- 48L
#' score_baseline <- rnorm(n_wide, mean = 50, sd = 7)
#' wide <- data.frame(
#'   subject_id = sprintf("S%03d", seq_len(n_wide)),
#'   score_t0 = score_baseline,
#'   score_t1 = score_baseline - 2 + rnorm(n_wide, sd = 1.8),
#'   score_t2 = score_baseline - 5 + rnorm(n_wide, sd = 2.0)
#' )
#'
#' # analyses = "none" keeps this first workflow fast. It disables optional
#' # modules, but not descriptives, changes, missingness, variability,
#' # trajectories, or the base-R Friedman branch. Holm adjustment is applied
#' # separately within each implemented comparison family.
#' auto_fit <- mira_info(
#'   wide,
#'   analyses = "none",
#'   p_adjust_method = "holm",
#'   verbose = FALSE
#' )
#'
#' # Start with sample size and completeness, then verify the ordered mapping
#' # from source columns to display labels.
#' auto_fit$overview[c(
#'   "n_patients", "n_timepoints", "complete_profiles",
#'   "complete_profiles_pct"
#' )]
#' data.frame(
#'   order = seq_along(auto_fit$time_vars),
#'   variable = auto_fit$time_vars,
#'   label = unname(auto_fit$time_labels)
#' )
#' utils::head(auto_fit$long_data[, c(
#'   "patient", "outcome", "time", "time_index", "time_label", "value"
#' )])
#'
#' # Per-timepoint availability and descriptive statistics.
#' auto_fit$descriptives[, c(
#'   "label", "n", "missing", "mean", "sd", "ci_lower", "ci_upper"
#' )]
#'
#' # Every ordered to-minus-from comparison is in $change. The adjusted columns
#' # belong to the paired-t and paired-Wilcoxon families, respectively.
#' auto_fit$change[, c(
#'   "from_label", "to_label", "n", "mean_change", "cohens_dz",
#'   "paired_t_p", "paired_t_p_adj", "wilcoxon_p_adj"
#' )]
#'
#' # Select the clinically common baseline-to-final row without relying on row
#' # position. This remains reliable if more visits are added later.
#' baseline_final <- auto_fit$change[
#'   auto_fit$change$from == auto_fit$time_vars[[1L]] &
#'     auto_fit$change$to == auto_fit$time_vars[[length(auto_fit$time_vars)]],
#'   ,
#'   drop = FALSE
#' ]
#' baseline_final[, c(
#'   "n", "mean_from", "mean_to", "mean_change", "ci_lower", "ci_upper"
#' )]
#'
#' # Variability and the preferred available ICC are kept in a one-row table.
#' auto_fit$variability[, c(
#'   "between_subject_sd", "within_subject_sd",
#'   "ICC_anova_complete_profiles", "ICC_model", "ICC"
#' )]
#'
#' # $trajectories contains one baseline-to-final record per input row.
#' utils::head(auto_fit$trajectories[, c(
#'   "patient", "baseline", "final", "absolute_change", "direction"
#' )])
#' table(auto_fit$trajectories$direction, useNA = "ifany")
#'
#' # Friedman is attempted with base R even though model fitting was disabled.
#' auto_fit$advanced_tests$friedman$tidy[, c(
#'   "statistic", "df", "p_raw", "kendalls_w", "n", "n_timepoints"
#' )]
#' auto_fit$advanced_tests$friedman$posthoc$tidy[, c(
#'   "from_label", "to_label", "n", "p_raw", "p_holm", "rank_biserial"
#' )]
#'
#' # ------------------------------------------------------------------
#' # Example 2: Audit what mira_info() detected before fitting anything
#' # ------------------------------------------------------------------
#' # This reusable trial has a genuine time-by-arm difference plus two
#' # time-invariant covariates. Eighty subjects and four visits also provide a
#' # stable basis for the model examples below.
#' set.seed(202602)
#' n_trial <- 80L
#' trial_arm <- rep(c("Control", "Treatment"), each = n_trial / 2L)
#' trial_site <- factor(rep(c("North", "South", "East", "West"), length.out = n_trial))
#' trial_age <- rep(50:69, length.out = n_trial) +
#'   sample(-1:1, n_trial, replace = TRUE)
#' treatment_indicator <- as.numeric(trial_arm == "Treatment")
#' subject_intercept <- rnorm(n_trial, sd = 5.5)
#' subject_slope <- rnorm(n_trial, sd = 0.8)
#' site_shift <- c(North = 0, South = 1.0, East = -0.8, West = 0.5)
#' trial_baseline <- 55 + 0.12 * (trial_age - 60) +
#'   unname(site_shift[as.character(trial_site)]) + subject_intercept +
#'   rnorm(n_trial, sd = 1.4)
#' trial <- data.frame(
#'   subject_id = sprintf("T%03d", seq_len(n_trial)),
#'   arm = factor(trial_arm, levels = c("Control", "Treatment")),
#'   age = trial_age,
#'   site = trial_site,
#'   outcome_t0 = trial_baseline,
#'   outcome_t1 = trial_baseline + subject_slope - 1 -
#'     3 * treatment_indicator + rnorm(n_trial, sd = 1.5),
#'   outcome_t2 = trial_baseline + 2 * subject_slope - 2 -
#'     6 * treatment_indicator + rnorm(n_trial, sd = 1.5),
#'   outcome_t3 = trial_baseline + 3 * subject_slope - 3 -
#'     9 * treatment_indicator + rnorm(n_trial, sd = 1.5)
#' )
#'
#' # inspect_only validates and resolves configuration but performs no outcome
#' # analysis. The requested optional switches are still recorded in $config.
#' detected <- mira_info(
#'   trial,
#'   analyses = c("model", "arm_tests"),
#'   inspect_only = TRUE,
#'   verbose = FALSE
#' )
#'
#' # These are the principal resolved choices. Control is selected as the
#' # reference because its label is recognized as control-like.
#' detected$config[c(
#'   "id", "id_generated", "outcomes", "arm", "reference_arm", "covariates"
#' )]
#' detected$config$time_vars[["outcome"]]
#' detected$config$time_labels[["outcome"]]
#' detected$config$analyses
#'
#' # These fields distinguish inferred choices from explicit overrides and
#' # record why an ID, arm, direction, or threshold was chosen.
#' detected$config$auto_detected
#' detected$config$specified_manually
#' detected$config$sources
#' detected$config$alternatives
#'
#' # The parsed map is the compact audit trail for outcome grouping and time
#' # ordering. Candidate tables and covariate classes explain the other choices.
#' detected$detected_variables$longitudinal_map[, c(
#'   "outcome", "variable", "time_label", "time_order", "pattern", "numeric"
#' )]
#' detected$detected_variables$id_candidates
#' detected$detected_variables$arm_candidates
#' detected$detected_variables$covariates_numeric
#' detected$detected_variables$covariates_categorical
#'
#' # Empty vectors are informative: they mean that no unresolved alternatives or
#' # diagnostic warnings remained. Non-empty entries explain what needs review.
#' detected$config$alternatives
#' detected$diagnostics$warnings
#'
#' # ------------------------------------------------------------------
#' # Example 3: Take control when column names do not encode time
#' # ------------------------------------------------------------------
#' set.seed(202603)
#' n_explicit <- 42L
#' symptom_baseline <- rnorm(n_explicit, mean = 30, sd = 5)
#' explicit_data <- data.frame(
#'   case_id = sprintf("E%03d", seq_len(n_explicit)),
#'   pre = symptom_baseline,
#'   week4 = symptom_baseline - 2 + rnorm(n_explicit, sd = 1.2),
#'   final = symptom_baseline - 4 + rnorm(n_explicit, sd = 1.4)
#' )
#'
#' # None of pre, week4, and final contains an outcome stem separated from a
#' # recognized suffix. Supply the ID, ordered columns, display labels, and
#' # outcome name explicitly. The supplied order defines baseline and final;
#' # stable_threshold = 1 treats changes within +/-1 symptom point as stable.
#' explicit_fit <- mira_info(
#'   explicit_data,
#'   id = "case_id",
#'   time_vars = c("pre", "week4", "final"),
#'   time_labels = c("Baseline", "Week 4", "Week 12"),
#'   outcomes = "symptom",
#'   covariates = character(0),
#'   analyses = "none",
#'   improvement_direction = "lower",
#'   stable_threshold = 1,
#'   strict_id = TRUE,
#'   verbose = FALSE
#' )
#'
#' data.frame(
#'   order = seq_along(explicit_fit$time_vars),
#'   variable = explicit_fit$time_vars,
#'   label = unname(explicit_fit$time_labels)
#' )
#' explicit_fit$config$outcomes
#' explicit_fit$config$specified_manually[c(
#'   "id", "outcomes", "time_vars", "time_labels"
#' )]
#' explicit_fit$change[
#'   explicit_fit$change$from == "pre" & explicit_fit$change$to == "final",
#'   c(
#'     "from_label", "to_label", "n", "mean_change",
#'     "improved_n", "stable_n", "worsened_n", "paired_t_p_adj"
#'   ),
#'   drop = FALSE
#' ]
#'
#' # Wide format means one row per subject. strict_id = TRUE is the protective
#' # default; these quick checks should both be zero before analysis.
#' c(
#'   missing_ids = sum(is.na(explicit_data$case_id)),
#'   duplicated_ids = sum(duplicated(explicit_data$case_id))
#' )
#'
#' # ------------------------------------------------------------------
#' # Example 4: Parse a non-standard naming convention with a regex
#' # ------------------------------------------------------------------
#' set.seed(202604)
#' n_regex <- 45L
#' pain_baseline <- rnorm(n_regex, mean = 7, sd = 1.3)
#' regex_data <- data.frame(
#'   record_id = sprintf("R%03d", seq_len(n_regex)),
#'   pain__OBS_12 = pain_baseline - 3.0 + rnorm(n_regex, sd = 0.7),
#'   pain__OBS_0 = pain_baseline,
#'   pain__OBS_4 = pain_baseline - 1.4 + rnorm(n_regex, sd = 0.6)
#' )
#'
#' # The first capture group is the outcome and the second is the time label.
#' # Numeric order is derived from that second group. The input columns are
#' # intentionally out of order so the resolved order is visible.
#' obs_pattern <- "^(.+)__OBS_([0-9]+)$"
#' regex_detected <- mira_info(
#'   regex_data,
#'   variable_pattern = obs_pattern,
#'   analyses = "none",
#'   inspect_only = TRUE,
#'   verbose = FALSE
#' )
#' regex_detected$detected_variables$longitudinal_map[, c(
#'   "outcome", "variable", "time_label", "time_order", "pattern"
#' )]
#' regex_detected$config$time_vars[["pain"]]
#' regex_detected$config$time_labels[["pain"]]
#'
#' regex_fit <- mira_info(
#'   regex_data,
#'   variable_pattern = obs_pattern,
#'   outcomes = "pain",
#'   analyses = "none",
#'   verbose = FALSE
#' )
#' data.frame(
#'   order = seq_along(regex_fit$time_vars),
#'   variable = regex_fit$time_vars,
#'   label = unname(regex_fit$time_labels)
#' )
#' regex_fit$change[
#'   regex_fit$change$from == regex_fit$time_vars[[1L]] &
#'     regex_fit$change$to == regex_fit$time_vars[[length(regex_fit$time_vars)]],
#'   c("from_label", "to_label", "n", "mean_change", "paired_t_p_adj"),
#'   drop = FALSE
#' ]
#'
#' # ------------------------------------------------------------------
#' # Example 5: Treatment-arm comparisons without fitting a model
#' # ------------------------------------------------------------------
#' # Explicit arguments make the scientific comparison unambiguous. A change
#' # within +/-2 outcome points is classified as stable for this illustration.
#' arm_fit <- mira_info(
#'   trial,
#'   id = "subject_id",
#'   outcomes = "outcome",
#'   arm = "arm",
#'   reference_arm = "Control",
#'   covariates = character(0),
#'   analyses = "arm_tests",
#'   p_adjust_method = "holm",
#'   improvement_direction = "lower",
#'   stable_threshold = 2,
#'   verbose = FALSE
#' )
#' arm_fit$arm_analysis[c("enabled", "arm_variable", "reference_arm", "levels")]
#'
#' # Values are compared between arms at every visit. The column name states the
#' # direction exactly: mean_difference_b_minus_a is arm_b minus arm_a.
#' arm_fit$arm_analysis$time_pairwise[, c(
#'   "time_label", "arm_a", "arm_b", "n_a", "n_b",
#'   "mean_difference_b_minus_a", "ci_lower", "ci_upper",
#'   "welch_t_p", "welch_t_p_adj"
#' )]
#'
#' # Baseline balance is the first-timepoint subset; omnibus p-values form their
#' # own across-time families and therefore have their own adjusted columns.
#' arm_fit$arm_analysis$baseline_balance[, c(
#'   "arm_a", "arm_b", "mean_difference_b_minus_a", "welch_t_p_adj"
#' )]
#' arm_fit$arm_analysis$time_omnibus[, c(
#'   "time_label", "n", "welch_anova_p", "welch_anova_p_adj",
#'   "kruskal_p_adj"
#' )]
#'
#' # Change comparisons are baseline-to-follow-up differences. Here a negative
#' # difference_in_change_b_minus_a means a larger reduction in Treatment.
#' arm_fit$arm_analysis$change_pairwise[, c(
#'   "to_label", "arm_a", "arm_b", "n_a", "n_b",
#'   "mean_change_a", "mean_change_b", "difference_in_change_b_minus_a",
#'   "ci_lower", "ci_upper", "welch_t_p_adj"
#' )]
#' arm_fit$arm_analysis$change_descriptives[, c(
#'   "arm", "to_label", "n", "mean_change",
#'   "improved_n", "stable_n", "worsened_n"
#' )]
#'
#' \donttest{
#' # ------------------------------------------------------------------
#' # Example 6: Covariate-adjusted longitudinal model
#' # ------------------------------------------------------------------
#' # lme4 is the required engine for the primary mixed model. Other advanced
#' # packages are optional and record a structured skip when unavailable.
#' if (requireNamespace("lme4", quietly = TRUE)) {
#'   adjusted_fit <- mira_info(
#'     trial,
#'     id = "subject_id",
#'     outcomes = "outcome",
#'     arm = "arm",
#'     reference_arm = "Control",
#'     covariates = c("age", "site"),
#'     analyses = c("model", "arm_tests"),
#'     p_adjust_method = "holm",
#'     improvement_direction = "lower",
#'     stable_threshold = 2,
#'     verbose = FALSE
#'   )
#'
#'   # Requested, retained, and skipped covariates are distinguished explicitly.
#'   # Covariate adjustment changes model-based estimates; it does not make the
#'   # treatment comparison causal.
#'   adjusted_fit$model[c(
#'     "covariates_requested", "covariates_used", "covariates_skipped",
#'     "converged", "singular", "warnings", "error"
#'   )]
#'   if (!is.null(adjusted_fit$model$fitted_model)) {
#'     stats::formula(adjusted_fit$model$fitted_model)
#'   }
#'
#'   # These likelihood-ratio tables answer different global questions. They are
#'   # not combined into the contrast-level multiplicity families.
#'   adjusted_fit$model$global_time_test
#'   adjusted_fit$model$global_arm_test
#'   adjusted_fit$model$arm_time_interaction_test
#'   adjusted_fit$variability[, c("ICC_model", "ICC")]
#'
#'   # ----------------------------------------------------------------
#'   # Example 7: Estimated marginal means and targeted contrasts
#'   # ----------------------------------------------------------------
#'   if (requireNamespace("emmeans", quietly = TRUE) &&
#'       isTRUE(adjusted_fit$advanced_tests$emmeans$performed)) {
#'     emm <- adjusted_fit$advanced_tests$emmeans
#'
#'     # This named list is a compact path index. It shows where each estimand or
#'     # contrast family lives without printing every table.
#'     emm_tables <- list(
#'       time_means = emm$time$tidy,
#'       all_time_pairs = emm$time$pairwise_tidy,
#'       baseline_vs_followup = emm$baseline_followup$tidy,
#'       baseline_vs_final = emm$baseline_final$tidy,
#'       consecutive_time = emm$consecutive$tidy,
#'       arm_pairs = emm$arm$pairwise_tidy,
#'       arm_within_time = emm$arm_time$simple_arm_tidy,
#'       time_within_arm = emm$arm_time$simple_time_tidy,
#'       interaction_contrasts = emm$arm_time$interaction_tidy,
#'       ordinal_trends = emm$trends$tidy
#'     )
#'     vapply(
#'       emm_tables,
#'       function(table) if (is.null(table)) 0L else nrow(table),
#'       integer(1L)
#'     )
#'
#'     # Estimated marginal means are adjusted for age and site because those
#'     # covariates were included in the fitted model.
#'     emm$time$tidy
#'
#'     compact_contrasts <- function(table) {
#'       if (is.null(table)) return(NULL)
#'       keep <- intersect(
#'         c("contrast", "estimate", "SE", "df", "p_raw", "p_holm"),
#'         names(table)
#'       )
#'       table[, keep, drop = FALSE]
#'     }
#'
#'     # Pairwise enumerates every pair of visits. Baseline-versus-final is one
#'     # prespecified endpoint contrast, stored separately for direct retrieval.
#'     compact_contrasts(emm$time$pairwise_tidy)
#'     compact_contrasts(emm$baseline_final$tidy)
#'
#'     # Simple arm effects compare arms within each visit. The path index above
#'     # also exposes time-within-arm and interaction contrasts.
#'     compact_contrasts(emm$arm_time$simple_arm_tidy)
#'
#'     # With four ordered visits, linear, quadratic, and cubic trends can be
#'     # returned. These use visit order, not assumed equal chronological spacing.
#'     compact_contrasts(emm$trends$tidy)
#'   }
#'
#'   # ----------------------------------------------------------------
#'   # Example 8: Advanced modules and sensitivity analysis
#'   # ----------------------------------------------------------------
#'   # Check performed before attempting to read module-specific tables. This
#'   # remains safe when some optional packages are not installed or a fit fails.
#'   advanced_status <- c(
#'     rm_anova = isTRUE(adjusted_fit$advanced_tests$rm_anova$performed),
#'     friedman = isTRUE(adjusted_fit$advanced_tests$friedman$performed),
#'     gee_independence = isTRUE(
#'       adjusted_fit$advanced_models$gee$independence$performed
#'     ),
#'     gee_exchangeable = isTRUE(
#'       adjusted_fit$advanced_models$gee$exchangeable$performed
#'     ),
#'     gee_ar1 = isTRUE(adjusted_fit$advanced_models$gee$ar1$performed),
#'     random_slope = isTRUE(
#'       adjusted_fit$advanced_models$random_slope$performed
#'     ),
#'     nlme_compound_symmetry = isTRUE(
#'       adjusted_fit$advanced_models$nlme$compound_symmetry$performed
#'     ),
#'     nlme_ar1 = isTRUE(adjusted_fit$advanced_models$nlme$ar1$performed),
#'     robust_cr2 = isTRUE(
#'       adjusted_fit$robustness$club_sandwich$performed
#'     )
#'   )
#'   advanced_status
#'
#'   if (advanced_status[["rm_anova"]]) {
#'     adjusted_fit$advanced_tests$rm_anova$tidy
#'   }
#'   adjusted_fit$advanced_tests$friedman$tidy
#'
#'   if (advanced_status[["gee_exchangeable"]]) {
#'     adjusted_fit$advanced_models$gee$exchangeable$effect_tests
#'   }
#'   if (advanced_status[["random_slope"]]) {
#'     adjusted_fit$advanced_models$random_slope[c(
#'       "formula", "converged", "singular",
#'       "random_intercept_slope_correlation"
#'     )]
#'   }
#'   if (advanced_status[["nlme_compound_symmetry"]] ||
#'       advanced_status[["nlme_ar1"]]) {
#'     adjusted_fit$advanced_models$nlme$comparison
#'   }
#'   if (advanced_status[["robust_cr2"]]) {
#'     utils::head(adjusted_fit$robustness$club_sandwich$coefficient_tests)
#'   }
#'
#'   # Likelihood criteria and GEE QIC/CIC are retained on different scales.
#'   model_comparison <- adjusted_fit$advanced_models$model_comparison
#'   if (is.data.frame(model_comparison) && nrow(model_comparison) > 0L) {
#'     model_comparison[, c(
#'       "model", "AIC", "BIC", "QIC", "CIC", "converged", "singular",
#'       "object_path"
#'     )]
#'   }
#'
#'   # Effect-size and multiplicity objects point back to their source families.
#'   paired_dz <- adjusted_fit$effect_sizes$paired_cohens_dz
#'   if (!is.null(paired_dz)) {
#'     paired_dz[, c("from_label", "to_label", "n", "cohens_dz")]
#'   }
#'   friedman_family <- adjusted_fit$multiplicity$friedman_posthoc$table
#'   if (!is.null(friedman_family)) {
#'     friedman_family[, c(
#'       "from_label", "to_label", "p_raw", "p_holm", "p_bh"
#'     )]
#'   }
#'   if (isTRUE(adjusted_fit$advanced_tests$emmeans$performed)) {
#'     time_family <- adjusted_fit$multiplicity$time_pairwise$table
#'     if (!is.null(time_family)) {
#'       time_family[, c("contrast", "p_raw", "p_holm", "p_bh")]
#'     }
#'   }
#'
#'   # Sensitivity rows align methods by scientific question and give the exact
#'   # source path. They may use different samples or hypotheses, so this table
#'   # must not be interpreted as a vote among p-values.
#'   adjusted_fit$sensitivity[, c(
#'     "question", "method", "statistic", "df", "p_raw", "effect_size",
#'     "effect_size_type", "n", "object_path"
#'   )]
#' }
#' }
#'
#' # ------------------------------------------------------------------
#' # Example 9: Multiple outcomes with different improvement directions
#' # ------------------------------------------------------------------
#' set.seed(202609)
#' n_multi <- 54L
#' response_group <- rep(c("Improved", "Stable", "Worsened"), each = 18L)
#' symptom_baseline_multi <- rnorm(n_multi, mean = 40, sd = 5)
#' function_baseline_multi <- rnorm(n_multi, mean = 55, sd = 7)
#' symptom_delta <- unname(
#'   c(Improved = -6, Stable = 0, Worsened = 4)[response_group]
#' ) + rnorm(n_multi, sd = 0.15)
#' function_delta <- unname(
#'   c(Improved = 8, Stable = 0, Worsened = -5)[response_group]
#' ) + rnorm(n_multi, sd = 0.20)
#' multi_data <- data.frame(
#'   subject_id = sprintf("M%03d", seq_len(n_multi)),
#'   symptom_t0 = symptom_baseline_multi,
#'   symptom_t1 = symptom_baseline_multi + 0.5 * symptom_delta +
#'     rnorm(n_multi, sd = 0.15),
#'   symptom_t2 = symptom_baseline_multi + symptom_delta,
#'   function_t0 = function_baseline_multi,
#'   function_t1 = function_baseline_multi + 0.5 * function_delta +
#'     rnorm(n_multi, sd = 0.20),
#'   function_t2 = function_baseline_multi + function_delta
#' )
#'
#' # Symptom is better when lower; function is better when higher. The named
#' # thresholds mean that changes within +/-1.5 symptom points or +/-2 function
#' # points are treated as stable. These are user-supplied practical thresholds,
#' # not values estimated or clinically validated by mira_info().
#' directions <- stats::setNames(
#'   c("lower", "higher"), c("symptom", "function")
#' )
#' practical_thresholds <- stats::setNames(
#'   c(1.5, 2), c("symptom", "function")
#' )
#' multi_fit <- mira_info(
#'   multi_data,
#'   outcomes = c("symptom", "function"),
#'   improvement_direction = directions,
#'   stable_threshold = practical_thresholds,
#'   analyses = "none",
#'   verbose = FALSE
#' )
#'
#' # A multi-outcome result is an outcome-indexed container. Each element is a
#' # complete single-outcome mira_info object with the same internal paths.
#' names(multi_fit$outcomes)
#' multi_fit$config$improvement_direction
#' multi_fit$config$stable_threshold
#' multi_fit$outcomes[["symptom"]]$settings[c(
#'   "improvement_direction", "stable_threshold"
#' )]
#' summary(multi_fit)
#'
#' # Extract the same result for every outcome instead of exploring nested lists.
#' baseline_to_final <- function(result) {
#'   result$change[
#'     result$change$from == result$time_vars[[1L]] &
#'       result$change$to == result$time_vars[[length(result$time_vars)]],
#'     c(
#'       "from_label", "to_label", "n", "mean_change",
#'       "improved_n", "stable_n", "worsened_n"
#'     ),
#'     drop = FALSE
#'   ]
#' }
#' lapply(
#'   multi_fit$outcomes,
#'   function(result) result$descriptives[, c("label", "n", "mean", "sd")]
#' )
#' lapply(multi_fit$outcomes, baseline_to_final)
#'
#' # Compare response labels under a zero threshold and the practical thresholds.
#' # With a zero threshold, even negligible simulated deviations count as change.
#' multi_zero <- mira_info(
#'   multi_data,
#'   outcomes = c("symptom", "function"),
#'   improvement_direction = directions,
#'   stable_threshold = stats::setNames(c(0, 0), c("symptom", "function")),
#'   analyses = "none",
#'   verbose = FALSE
#' )
#' response_counts <- function(result) {
#'   lapply(
#'     result$outcomes,
#'     function(outcome_result) {
#'       table(outcome_result$trajectories$clinical_direction, useNA = "ifany")
#'     }
#'   )
#' }
#' list(
#'   zero_threshold = response_counts(multi_zero),
#'   practical_threshold = response_counts(multi_fit)
#' )
#'
#' \donttest{
#' # The multi-outcome plot method needs outcome = when several outcomes succeed.
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   multi_plot_fit <- mira_info(
#'     multi_data,
#'     outcomes = c("symptom", "function"),
#'     improvement_direction = directions,
#'     stable_threshold = practical_thresholds,
#'     analyses = "plots",
#'     verbose = FALSE
#'   )
#'   plot(multi_plot_fit, outcome = "symptom", which = "mean_ci")
#' }
#' }
#'
#' # ------------------------------------------------------------------
#' # Example 10: Missing follow-up data and effective sample sizes
#' # ------------------------------------------------------------------
#' missing_data <- wide
#' missing_data$score_t1[c(4, 9, 15, 22)] <- NA_real_
#' missing_data$score_t2[c(2, 4, 11, 19, 28, 37)] <- NA_real_
#' missing_fit <- mira_info(
#'   missing_data,
#'   analyses = "correlations",
#'   verbose = FALSE
#' )
#'
#' # mira_info() does not impute. Descriptives use all finite observations at
#' # each visit, while complete_profiles requires all selected visits.
#' missing_fit$overview[c(
#'   "n_patients", "complete_profiles", "complete_profiles_pct"
#' )]
#' missing_fit$descriptives[, c(
#'   "label", "n", "missing", "non_finite", "unavailable"
#' )]
#' missing_fit$missing$by_time
#'
#' # Identify affected subjects without printing all rows.
#' affected <- missing_fit$missing$by_patient$unavailable_n > 0L
#' missing_fit$missing$by_patient[affected, , drop = FALSE]
#'
#' # Each paired comparison has its own pairwise-complete N. Friedman instead
#' # uses complete profiles across all visits, reported in its tidy row.
#' missing_fit$change[, c(
#'   "from_label", "to_label", "n", "mean_change", "paired_t_p_adj"
#' )]
#' missing_fit$advanced_tests$friedman$tidy[, c("n", "n_timepoints", "p_raw")]
#'
#' # Per-outcome adaptation records show usable counts and any requested module
#' # disabled because of insufficient data.
#' missing_fit$diagnostics$adaptation[["score"]]
#'
#' # ------------------------------------------------------------------
#' # Example 11: Correlations, outlier diagnostics, and public plot methods
#' # ------------------------------------------------------------------
#' diagnostic_data <- wide
#' diagnostic_data$score_t2[[1L]] <- diagnostic_data$score_t2[[1L]] + 35
#' diagnostic_fit <- mira_info(
#'   diagnostic_data,
#'   analyses = c("correlations", "outliers", "plots"),
#'   verbose = FALSE
#' )
#'
#' diagnostic_fit$correlations$pearson
#' diagnostic_fit$correlations$pairwise_n
#'
#' # Outliers are Tukey 1.5-IQR flags only; they are never removed automatically.
#' vapply(diagnostic_fit$outliers$by_time, nrow, integer(1L))
#' diagnostic_fit$outliers$by_time[["score_t2"]]
#' change_flag_counts <- vapply(
#'   diagnostic_fit$outliers$change, nrow, integer(1L)
#' )
#' change_flag_counts
#' most_flagged_change <- names(change_flag_counts)[which.max(change_flag_counts)]
#' diagnostic_fit$outliers$change[[most_flagged_change]]
#'
#' # Use plot() rather than reaching into $plots when displaying a stored plot.
#' # If ggplot2 is absent, the analysis still succeeds and $plot_error explains
#' # why no plot was stored.
#' diagnostic_fit$plot_error
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   names(diagnostic_fit$plots)
#'   plot(diagnostic_fit, which = "mean_ci")
#'   plot(diagnostic_fit, which = "spaghetti")
#'   plot(diagnostic_fit, which = "missingness")
#' }
#'
#' # ------------------------------------------------------------------
#' # Example 12: Select exactly which optional modules to request
#' # ------------------------------------------------------------------
#' # inspect_only makes this comparison cheap: no statistical module is run.
#' cfg_none <- mira_info(
#'   wide, analyses = "none", inspect_only = TRUE, verbose = FALSE
#' )
#' cfg_arm <- mira_info(
#'   trial, analyses = "arm_tests", inspect_only = TRUE, verbose = FALSE
#' )
#' cfg_model_cor <- mira_info(
#'   wide,
#'   analyses = c("model", "correlations"),
#'   inspect_only = TRUE,
#'   verbose = FALSE
#' )
#' cfg_all <- mira_info(
#'   trial, analyses = "all", inspect_only = TRUE, verbose = FALSE
#' )
#'
#' module_switches <- function(configuration) {
#'   unlist(
#'     configuration$config$analyses[c(
#'       "plots", "model", "outliers", "correlations", "arm_tests"
#'     )],
#'     use.names = TRUE
#'   )
#' }
#' rbind(
#'   none = module_switches(cfg_none),
#'   arm_tests = module_switches(cfg_arm),
#'   model_correlations = module_switches(cfg_model_cor),
#'   all = module_switches(cfg_all)
#' )
#'
#' # Core summaries remain enabled under analyses = "none". On a real analysis,
#' # the base-R Friedman branch is also attempted independently of model =.
#' cfg_none$config$analyses$core
#' @export
mira_info <- function(data,
                      id = NULL,
                      time_vars = NULL,
                      time_labels = NULL,
                      arm = NULL,
                      reference_arm = NULL,
                      arm_tests = TRUE,
                      alpha = 0.05,
                      plots = TRUE,
                      model = TRUE,
                      outliers = TRUE,
                      correlations = TRUE,
                      verbose = TRUE,
                      p_adjust_method = "holm",
                      improvement_direction = "auto",
                      stable_threshold = "auto",
                      strict_id = TRUE,
                      outcomes = NULL,
                      covariates = NULL,
                      analyses = NULL,
                      variable_pattern = "auto",
                      inspect_only = FALSE) {
  scalar_flag <- function(x, name) {
    if (!is.logical(x) || length(x) != 1L || is.na(x)) {
      stop(sprintf("%s must be TRUE or FALSE.", name), call. = FALSE)
    }
  }
  flags <- list(
    arm_tests = arm_tests, plots = plots, model = model, outliers = outliers,
    correlations = correlations, verbose = verbose, strict_id = strict_id,
    inspect_only = inspect_only
  )
  for (flag_name in names(flags)) {
    scalar_flag(flags[[flag_name]], flag_name)
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) ||
      alpha <= 0 || alpha >= 1) {
    stop("alpha must be a number strictly between 0 and 1.", call. = FALSE)
  }
  if (!is.character(p_adjust_method) || length(p_adjust_method) != 1L ||
      is.na(p_adjust_method) || !p_adjust_method %in% p.adjust.methods) {
    stop(sprintf("p_adjust_method must be one of: %s.",
                 paste(p.adjust.methods, collapse = ", ")), call. = FALSE)
  }

  optional_analyses <- c("plots", "model", "outliers", "correlations", "arm_tests")
  if (!is.null(analyses)) {
    if (!is.character(analyses) || anyNA(analyses)) {
      stop("analyses must be NULL or a vector of names.", call. = FALSE)
    }
    analyses <- unique(tolower(analyses))
    if ("all" %in% analyses) analyses <- optional_analyses
    if ("none" %in% analyses) analyses <- character(0)
    invalid <- setdiff(analyses, optional_analyses)
    if (length(invalid) > 0L) {
      stop(sprintf("Unrecognized optional analyses: %s. Allowed values: %s.",
                   paste(invalid, collapse = ", "),
                   paste(optional_analyses, collapse = ", ")), call. = FALSE)
    }
    plots <- "plots" %in% analyses
    model <- "model" %in% analyses
    outliers <- "outliers" %in% analyses
    correlations <- "correlations" %in% analyses
    arm_tests <- "arm_tests" %in% analyses
  }

  built <- .mira_build_analysis_config(
    data = data, id = id, outcomes = outcomes, time_vars = time_vars,
    time_labels = time_labels, arm = arm, reference_arm = reference_arm,
    covariates = covariates, variable_pattern = variable_pattern,
    improvement_direction = improvement_direction,
    stable_threshold = stable_threshold, strict_id = strict_id
  )
  built$config$analyses <- list(
    core = c("descriptives", "missingness", "change", "variability", "trajectories"),
    plots = plots, model = model, outliers = outliers,
    correlations = correlations, arm_tests = arm_tests
  )
  for (message in built$warnings) warning(message, call. = FALSE)

  if (inspect_only) {
    result <- list(
      call = match.call(), version = "4.0.0", config = built$config,
      data_overview = built$data_overview,
      detected_variables = built$detected_variables,
      diagnostics = list(warnings = built$warnings)
    )
    class(result) <- c("mira_detect", "list")
    if (verbose) print(result)
    return(invisible(result))
  }

  outcome_results <- list()
  adaptations <- list()
  multiple <- length(built$groups) > 1L

  for (outcome_name in names(built$groups)) {
    group <- built$groups[[outcome_name]]
    availability <- vapply(group$variables, function(variable) {
      x <- built$data[[variable]]
      sum(!is.na(x) & is.finite(x))
    }, integer(1L))
    all_values <- unlist(lapply(group$variables, function(variable) {
      x <- built$data[[variable]]
      x[!is.na(x) & is.finite(x)]
    }), use.names = FALSE)
    variable_timepoints <- sum(vapply(group$variables, function(variable) {
      x <- built$data[[variable]]
      x <- x[!is.na(x) & is.finite(x)]
      length(unique(x)) >= 2L
    }, logical(1L)))

    outcome_model <- model && length(all_values) >= 3L && variable_timepoints >= 1L
    outcome_correlations <- correlations && sum(availability >= 2L) >= 2L &&
      variable_timepoints >= 2L
    outcome_plots <- plots && length(all_values) > 0L
    outcome_arm_tests <- arm_tests && length(all_values) > 0L &&
      !is.null(built$config$arm) &&
      length(unique(built$data[[built$config$arm]][
        !is.na(built$data[[built$config$arm]])
      ])) >= 2L

    disabled <- character(0)
    if (model && !outcome_model) disabled <- c(disabled, "model: insufficient data/variability")
    if (correlations && !outcome_correlations) {
      disabled <- c(disabled, "correlations: fewer than two estimable timepoints")
    }
    if (plots && !outcome_plots) disabled <- c(disabled, "plots: no finite outcome values")
    if (arm_tests && !outcome_arm_tests) disabled <- c(disabled, "arm_tests: no usable arm")
    adaptations[[outcome_name]] <- list(
      available_by_time = stats::setNames(availability, group$variables),
      variable_timepoints = variable_timepoints,
      disabled = disabled
    )

    analyse <- function() {
      .mira_analyse_single_outcome(
        data = built$data,
        id = built$config$id,
        time_vars = group$variables,
        time_labels = group$time_labels,
        arm = built$config$arm,
        reference_arm = built$config$reference_arm,
        arm_tests = outcome_arm_tests,
        alpha = alpha,
        plots = outcome_plots,
        model = outcome_model,
        outliers = outliers,
        correlations = outcome_correlations,
        verbose = FALSE,
        p_adjust_method = p_adjust_method,
        improvement_direction = built$config$improvement_direction[[outcome_name]],
        stable_threshold = built$config$stable_threshold[[outcome_name]],
        strict_id = strict_id,
        .outcome_name = outcome_name,
        .covariates = built$config$covariates,
        .categorical_covariates = built$config$covariates_categorical
      )
    }

    outcome_result <- if (multiple) {
      tryCatch(
        analyse(),
        error = function(e) structure(
          list(outcome = outcome_name, error = conditionMessage(e)),
          class = c("mira_info_error", "list")
        )
      )
    } else analyse()

    if (!inherits(outcome_result, "mira_info_error")) {
      outcome_config <- built$config
      outcome_config$outcomes <- outcome_name
      outcome_config$time_vars <- group$variables
      outcome_config$time_labels <- stats::setNames(group$time_labels, group$variables)
      outcome_config$improvement_direction <-
        built$config$improvement_direction[[outcome_name]]
      outcome_config$stable_threshold <- built$config$stable_threshold[[outcome_name]]
      outcome_result$config <- outcome_config
      outcome_result$data_overview <- built$data_overview
      outcome_result$detected_variables <- built$detected_variables
      outcome_result$diagnostics <- list(
        warnings = built$warnings,
        adaptation = adaptations[[outcome_name]]
      )
      outcome_result$call <- match.call()
      outcome_result$version <- "4.0.0"
    }
    outcome_results[[outcome_name]] <- outcome_result
  }

  failed <- names(outcome_results)[vapply(outcome_results, inherits, logical(1L),
                                          what = "mira_info_error")]
  if (length(failed) > 0L) {
    warning(sprintf("Analyses were not completed for: %s. See the $outcomes component.",
                    paste(failed, collapse = ", ")), call. = FALSE)
  }

  if (length(outcome_results) == 1L && length(failed) == 0L) {
    result <- outcome_results[[1L]]
    result$config <- built$config
    result$data_overview <- built$data_overview
    result$detected_variables <- built$detected_variables
    result$diagnostics <- list(warnings = built$warnings, adaptation = adaptations)
    result$outcomes <- outcome_results
    result$call <- match.call()
    class(result) <- c("mira_info", "list")
  } else {
    result <- list(
      call = match.call(),
      version = "4.0.0",
      config = built$config,
      data_overview = built$data_overview,
      detected_variables = built$detected_variables,
      outcomes = outcome_results,
      diagnostics = list(warnings = built$warnings, adaptation = adaptations,
                         failed_outcomes = failed)
    )
    class(result) <- c("mira_info_multi", "mira_info", "list")
  }

  if (verbose) print(result)
  invisible(result)
}

.mira_output_component_descriptions <- c(
  call = "Matched function call that created the object.",
  version = "MIRA output-schema version recorded in this object.",
  settings = "Outcome-specific analysis settings actually applied.",
  overview = "Outcome-specific sample, timepoint, completeness, and arm overview.",
  outcome = "Machine-readable outcome identifier.",
  outcome_display = "Human-readable outcome label used in output and plots.",
  time_vars = "Ordered source-column names for the repeated measurements.",
  time_labels = "Named display labels for the ordered timepoint columns.",
  descriptives = "Per-timepoint descriptive statistics and confidence intervals.",
  missing = "Availability summaries by timepoint and by subject.",
  change = "All pairwise longitudinal changes, effect sizes, and adjusted tests.",
  arm_analysis = "Arm-specific descriptives, balance, missingness, and comparisons.",
  correlations = "Pearson and Spearman matrices with pairwise sample sizes.",
  variability = "Within- and between-subject variability with ICC estimates.",
  trajectories = "Subject-level baseline-to-final changes and direction labels.",
  model = "Mixed-effects model, diagnostics, ANOVA, and global tests.",
  advanced_tests = "RM-ANOVA, Friedman, emmeans, and longitudinal contrasts.",
  advanced_models = "Random-slope, nlme, and GEE fits with model comparison.",
  robustness = "Cluster-robust CR2 inference for the mixed-effects model.",
  effect_sizes = "Method-specific longitudinal effect-size results.",
  multiplicity = "Multiplicity-adjustment families and source-object paths.",
  sensitivity = "Effect-aligned inference comparison across methods.",
  outliers = "IQR diagnostic flags by timepoint and change comparison.",
  plots = "Available ggplot objects indexed by plot name.",
  plot_error = "Plot-generation error message, or NULL when no error occurred.",
  long_data = "Long-format analysis data used by downstream methods.",
  config = "Resolved configuration and automatic/manual selection provenance.",
  data_overview = "Dataset-wide dimensions and column-level profile.",
  detected_variables = "Detection candidates and longitudinal-variable map.",
  diagnostics = "Detection warnings, adaptations, and failed outcomes when applicable.",
  outcomes = "Outcome-indexed results, one entry per requested outcome."
)

.mira_print_output_guide <- function(x, title = "AVAILABLE OUTPUT ELEMENTS") {
  cat("\n", title, "\n", sep = "")
  cat(strrep("-", 84L), "\n", sep = "")

  components <- names(x)
  if (is.null(components) || length(components) == 0L) {
    cat("No named output elements are available.\n")
    return(invisible(NULL))
  }

  for (component in components) {
    description <- unname(.mira_output_component_descriptions[component])
    if (length(description) == 0L || is.na(description) || !nzchar(description)) {
      description <- "Additional output component."
    }
    cat(sprintf("  $%-20s %s\n", component, description))
  }

  invisible(NULL)
}

#' @export
print.mira_detect <- function(x, ...) {
  cat("MIRA INFO\n")
  cat(strrep("=", 84L), "\n", sep = "")
  cat("\nDETECTION OVERVIEW\n")
  cat(strrep("-", 84L), "\n", sep = "")
  cat(sprintf("Rows: %d | Columns: %d\n",
              x$data_overview$n_rows, x$data_overview$n_columns))
  cat(sprintf("ID: %s%s\n", x$config$id,
              if (isTRUE(x$config$id_generated)) " (generated safely)" else ""))
  cat(sprintf("Outcomes: %s\n", paste(x$config$outcomes, collapse = ", ")))
  for (outcome in x$config$outcomes) {
    variables <- x$config$time_vars[[outcome]]
    labels <- unname(x$config$time_labels[[outcome]])
    cat(sprintf("  %s: %s\n", outcome,
                paste(sprintf("%s [%s]", variables, labels), collapse = ", ")))
  }
  cat(sprintf("Arm: %s\n", if (is.null(x$config$arm)) "none" else x$config$arm))
  if (!is.null(x$config$arm)) {
    cat(sprintf("Reference arm: %s\n",
                if (is.null(x$config$reference_arm)) "none" else x$config$reference_arm))
  }
  cat(sprintf("Selected covariates: %s\n",
              if (length(x$config$covariates) == 0L) "none" else
                paste(x$config$covariates, collapse = ", ")))
  cat(sprintf("Detected numeric covariates: %s\n",
              if (length(x$detected_variables$covariates_numeric) == 0L) "none" else
                paste(x$detected_variables$covariates_numeric, collapse = ", ")))
  cat(sprintf("Detected categorical covariates: %s\n",
              if (length(x$detected_variables$covariates_categorical) == 0L) "none" else
                paste(x$detected_variables$covariates_categorical, collapse = ", ")))
  if (length(x$diagnostics$warnings) > 0L) {
    cat("Diagnostics:\n")
    cat(paste0("  - ", x$diagnostics$warnings, collapse = "\n"), "\n")
  }
  .mira_print_output_guide(x)
  invisible(x)
}

#' @export
print.mira_info_multi <- function(x, ...) {
  cat("MIRA INFO\n")
  cat(strrep("=", 84L), "\n", sep = "")
  cat("\nMULTI-OUTCOME OVERVIEW\n")
  cat(strrep("-", 84L), "\n", sep = "")
  cat(sprintf("Rows: %d | ID: %s | Outcomes: %d\n",
              x$data_overview$n_rows, x$config$id, length(x$outcomes)))
  cat(sprintf("Arm: %s | Covariates: %s\n",
              if (is.null(x$config$arm)) "none" else x$config$arm,
              if (length(x$config$covariates) == 0L) "none" else
                paste(x$config$covariates, collapse = ", ")))
  rows <- lapply(names(x$outcomes), function(outcome) {
    result <- x$outcomes[[outcome]]
    if (inherits(result, "mira_info_error")) {
      return(data.frame(outcome = outcome, timepoints = NA_integer_, subjects = NA_integer_,
                        complete_profiles = NA_integer_, status = result$error,
                        stringsAsFactors = FALSE))
    }
    data.frame(
      outcome = outcome,
      timepoints = result$overview$n_timepoints,
      subjects = result$overview$n_patients,
      complete_profiles = result$overview$complete_profiles,
      status = "ok",
      stringsAsFactors = FALSE
    )
  })
  print(do.call(rbind, rows), row.names = FALSE)
  .mira_print_output_guide(x)
  cat("\nAccess an outcome-specific result through $outcomes[[\"OUTCOME_NAME\"]]\n")
  invisible(x)
}

#' @export
print.mira_info_error <- function(x, ...) {
  cat("MIRA INFO\n")
  cat(strrep("=", 84L), "\n", sep = "")
  cat("\nANALYSIS ERROR\n")
  cat(strrep("-", 84L), "\n", sep = "")
  cat(sprintf("Outcome %s: analysis not completed - %s\n", x$outcome, x$error))
  invisible(x)
}


# ============================================================
# PRINT METHOD
# ============================================================

#' @export
print.mira_info <- function(x,
                            digits = 3,
                            max_rows = 20,
                            correlations = TRUE,
                            model = TRUE,
                            outliers = TRUE,
                            trajectories = TRUE,
                            missingness = TRUE,
                            plots = TRUE,
                            ...) {

  if (!is.numeric(digits) || length(digits) != 1L || is.na(digits) || digits < 0) {
    stop("digits must be an integer >= 0.", call. = FALSE)
  }
  digits <- as.integer(digits)

  if (!is.numeric(max_rows) || length(max_rows) != 1L || is.na(max_rows) || max_rows <= 0) {
    stop("max_rows must be > 0 or Inf.", call. = FALSE)
  }

  line <- function(char = "-", n = 84L) cat(strrep(char, n), "\n", sep = "")

  section <- function(title) {
    cat("\n", title, "\n", sep = "")
    line()
  }

  fmt_num <- function(z, d = digits) {
    ifelse(is.na(z), "NA", formatC(z, format = "f", digits = d))
  }

  fmt_pct <- function(z, d = 1L) {
    ifelse(is.na(z), "NA", paste0(formatC(z, format = "f", digits = d), "%"))
  }

  fmt_p <- function(p) {
    ifelse(
      is.na(p),
      "NA",
      ifelse(p < 0.001, "<0.001", formatC(p, format = "f", digits = 3L))
    )
  }

  has_rows <- function(z) {
    rows <- tryCatch(nrow(z), error = function(e) NULL)
    !is.null(rows) && length(rows) == 1L && !is.na(rows) && rows > 0L
  }

  fmt_text <- function(z, fallback = "none") {
    if (is.null(z) || length(z) == 0L) return(fallback)
    z <- as.character(z)
    z <- z[!is.na(z) & nzchar(z)]
    if (length(z) == 0L) fallback else paste(z, collapse = ", ")
  }

  extract_test_p <- function(test) {
    if (is.null(test)) return(NA_real_)
    table <- tryCatch(as.data.frame(test), error = function(e) NULL)
    if (is.null(table) || nrow(table) == 0L) return(NA_real_)
    p_col <- .mira_column_name(
      table,
      candidates = c("Pr(>Chisq)", "p.value", "p_value", "p"),
      contains = c("prchisq", "pvalue")
    )
    if (is.na(p_col)) return(NA_real_)
    values <- .mira_as_numeric(table[[p_col]])
    available <- which(is.finite(values))
    if (length(available) == 0L) NA_real_ else values[[tail(available, 1L)]]
  }

  effect_label <- function(d) {
    if (length(d) == 0L || is.na(d)) return(NA_character_)
    a <- abs(d)
    if (a < 0.20) "negligible" else if (a < 0.50) "small" else if (a < 0.80) "moderate" else "large"
  }

  limit_table <- function(z, label = "rows") {
    if (!is.data.frame(z) || nrow(z) == 0L) return(z)
    if (is.finite(max_rows) && nrow(z) > max_rows) {
      cat(sprintf("Showing the first %d of %d %s. Set max_rows = Inf in print() to display all rows.\n",
                  as.integer(max_rows), nrow(z), label))
      return(utils::head(z, as.integer(max_rows)))
    }
    z
  }

  ov <- x$overview
  d <- x$descriptives
  conf_pct <- x$settings$confidence_percent
  conf_label <- paste0(formatC(conf_pct, format = "fg", digits = 4L), "% CI")

  baseline_final <- NULL
  if (!is.null(x$change) && nrow(x$change) > 0L && length(ov$timepoints) >= 2L) {
    baseline_final <- x$change[
      x$change$from == ov$timepoints[[1L]] &
        x$change$to == ov$timepoints[[length(ov$timepoints)]],
      ,
      drop = FALSE
    ]
  }

  cat("MIRA INFO\n")
  line("=")

  # ------------------------------------------------------------------
  # KEY RESULTS
  # ------------------------------------------------------------------
  section("KEY RESULTS")
  cat(sprintf("Outcome: %s\n", ov$outcome_display))
  cat(sprintf(
    "Subjects: %d | Rows: %d | Timepoints: %d | Complete profiles: %d/%d (%s)\n",
    ov$n_patients, ov$n_rows, ov$n_timepoints,
    ov$complete_profiles, ov$n_rows, fmt_pct(ov$complete_profiles_pct)
  ))

  if (isTRUE(ov$arm_analysis)) {
    cat(sprintf(
      "Treatment arms: %d | Reference: %s | Arm variable: %s\n",
      length(ov$arm_levels), ov$reference_arm, ov$arm_variable
    ))
  }

  if (!is.null(d) && nrow(d) > 0L) {
    baseline_row <- d[1L, , drop = FALSE]
    cat(sprintf(
      "Baseline (%s): N=%d | mean=%s | SD=%s\n",
      baseline_row$label, baseline_row$n,
      fmt_num(baseline_row$mean), fmt_num(baseline_row$sd)
    ))
    if (nrow(d) > 1L) {
      final_row <- d[nrow(d), , drop = FALSE]
      cat(sprintf(
        "Final (%s): N=%d | mean=%s | SD=%s\n",
        final_row$label, final_row$n,
        fmt_num(final_row$mean), fmt_num(final_row$sd)
      ))
    }
  }

  if (!is.null(x$missing$by_time) && nrow(x$missing$by_time) > 0L) {
    total_cells <- ov$n_rows * ov$n_timepoints
    unavailable_cells <- sum(x$missing$by_time$unavailable_n, na.rm = TRUE)
    available_cells <- total_cells - unavailable_cells
    cat(sprintf(
      "Available longitudinal cells: %d/%d (%s)\n",
      available_cells, total_cells,
      fmt_pct(if (total_cells > 0L) available_cells / total_cells * 100 else NA_real_)
    ))
  }

  if (!is.null(baseline_final) && nrow(baseline_final) == 1L) {
    bf_key <- baseline_final[1L, , drop = FALSE]
    cat(sprintf(
      "Baseline-to-final change: mean=%s | %s [%s, %s] | Cohen dz=%s (%s)\n",
      fmt_num(bf_key$mean_change), conf_label,
      fmt_num(bf_key$ci_lower), fmt_num(bf_key$ci_upper),
      fmt_num(bf_key$cohens_dz), effect_label(bf_key$cohens_dz)
    ))
    cat(sprintf(
      "Adjusted paired tests: t=%s | Wilcoxon=%s\n",
      fmt_p(bf_key$paired_t_p_adj), fmt_p(bf_key$wilcoxon_p_adj)
    ))
    if (x$settings$improvement_direction != "unknown") {
      cat(sprintf(
        "Clinical direction: improved=%s | stable=%s | worsened=%s\n",
        fmt_pct(bf_key$improved_pct), fmt_pct(bf_key$stable_pct),
        fmt_pct(bf_key$worsened_pct)
      ))
    }
  }

  if (model) {
    if (!is.null(x$model$fitted_model)) {
      cat(sprintf(
        "Mixed model: available | converged=%s | singular=%s\n",
        as.character(x$model$converged), as.character(x$model$singular)
      ))
    } else if (!is.null(x$model$error)) {
      cat("Mixed model: unavailable | ", x$model$error, "\n", sep = "")
    } else {
      cat("Mixed model: not available\n")
    }

    global_tests <- character(0)
    if (!is.null(x$model$global_time_test)) {
      global_tests <- c(
        global_tests,
        paste0("time p=", fmt_p(extract_test_p(x$model$global_time_test)))
      )
    }
    if (!is.null(x$model$global_arm_test)) {
      global_tests <- c(
        global_tests,
        paste0("arm p=", fmt_p(extract_test_p(x$model$global_arm_test)))
      )
    }
    if (!is.null(x$model$arm_time_interaction_test)) {
      global_tests <- c(
        global_tests,
        paste0(
          "time-by-arm p=",
          fmt_p(extract_test_p(x$model$arm_time_interaction_test))
        )
      )
    }
    if (length(global_tests) > 0L) {
      cat("Global mixed-model tests: ", paste(global_tests, collapse = " | "), "\n", sep = "")
    }
  }

  if (outliers && !is.null(x$outliers$by_time)) {
    time_flags <- sum(vapply(x$outliers$by_time, nrow, integer(1L)))
    change_flags <- if (!is.null(x$outliers$change)) {
      sum(vapply(x$outliers$change, nrow, integer(1L)))
    } else 0L
    cat(sprintf(
      "IQR diagnostic flags: timepoint=%d | change=%d\n",
      time_flags, change_flags
    ))
  }

  # ------------------------------------------------------------------
  # ANALYSIS SETTINGS
  # ------------------------------------------------------------------
  section("ANALYSIS SETTINGS")
  cat(sprintf("Outcome: %s | ID variable: %s\n", ov$outcome_display, ov$id))
  cat(sprintf("Confidence level: %s | alpha: %s | p-value adjustment: %s\n",
              conf_label, fmt_num(x$settings$alpha), x$settings$p_adjust_method))
  cat(sprintf("Improvement direction: %s | Stable threshold: %s\n",
              x$settings$improvement_direction, fmt_num(x$settings$stable_threshold)))
  cat(sprintf("Strict ID checks: %s | Non-finite handling: %s\n",
              as.character(x$settings$strict_id), x$settings$non_finite_handling))
  if (!is.null(x$config)) {
    cat(sprintf("Selected covariates: %s\n",
                if (length(x$config$covariates) == 0L) "none" else
                  paste(x$config$covariates, collapse = ", ")))
  }

  if (isTRUE(ov$arm_analysis)) {
    cat(sprintf(
      "Arm variable: %s | Reference arm: %s | Groups: %d\n",
      ov$arm_variable,
      ov$reference_arm,
      length(ov$arm_levels)
    ))
  }

  if (!is.null(x$config$auto_detected)) {
    automatic <- x$config$auto_detected
    mode <- function(value) if (isTRUE(value)) "automatic" else "manual"
    id_mode <- if (isTRUE(x$config$id_generated)) "generated" else mode(automatic$id)
    arm_mode <- if (is.null(x$config$arm)) "none" else mode(automatic$arm)
    cat(sprintf(
      "Selection mode: ID=%s | outcome=%s | time variables=%s | arm=%s | covariates=%s\n",
      id_mode, mode(automatic$outcomes), mode(automatic$time_vars),
      arm_mode, mode(automatic$covariates)
    ))
    cat("Full configuration provenance: $config\n")
    cat("Detection candidates and map:  $detected_variables\n")
  }

  # ------------------------------------------------------------------
  # DATASET OVERVIEW
  # ------------------------------------------------------------------
  section("DATASET OVERVIEW")
  cat(sprintf("Subjects: %d | Rows: %d | Timepoints: %d\n",
              ov$n_patients, ov$n_rows, ov$n_timepoints))
  cat(sprintf("Complete profiles: %d/%d (%s)\n",
              ov$complete_profiles, ov$n_rows, fmt_pct(ov$complete_profiles_pct)))
  cat(sprintf("Duplicated IDs: %d | Missing IDs: %d | Non-finite values: %d\n",
              ov$duplicated_ids, ov$missing_ids, ov$non_finite_values))

  if (!is.null(ov$arm_variable)) {
    cat(sprintf("Missing arm values: %d\n", ov$missing_arm))
    if (!is.null(ov$arm_counts)) {
      cat("Subjects by arm:\n")
      print(ov$arm_counts)
    }
  }

  tp <- data.frame(
    Variable = ov$timepoints,
    Label = unname(ov$time_labels[ov$timepoints]),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  cat("\nTimepoint map:\n")
  print(tp, row.names = FALSE)

  # ------------------------------------------------------------------
  # DESCRIPTIVES
  # ------------------------------------------------------------------
  section("DESCRIPTIVE STATISTICS")
  desc_print <- data.frame(
    Time = d$label,
    N = d$n,
    Missing = paste0(d$missing, " (", round(d$missing_pct, 1), "%)"),
    Non_finite = paste0(d$non_finite, " (", round(d$non_finite_pct, 1), "%)"),
    Unavailable = paste0(d$unavailable, " (", round(d$unavailable_pct, 1), "%)"),
    Mean = round(d$mean, digits),
    SD = round(d$sd, digits),
    SE = round(d$se, digits),
    Median = round(d$median, digits),
    Q1 = round(d$q1, digits),
    Q3 = round(d$q3, digits),
    IQR = round(d$iqr, digits),
    Min = round(d$min, digits),
    Max = round(d$max, digits),
    CV_pct = round(d$cv_percent, 1L),
    CI_low = round(d$ci_lower, digits),
    CI_high = round(d$ci_upper, digits),
    check.names = FALSE
  )
  print(desc_print, row.names = FALSE)
  cat(sprintf("Confidence intervals above: %s\n", conf_label))

  # ------------------------------------------------------------------
  # TREATMENT ARM EXPLORATION
  # ------------------------------------------------------------------
  if (isTRUE(ov$arm_analysis) && !is.null(x$arm_analysis)) {

    section("TREATMENT ARM EXPLORATION")

    cat(sprintf(
      "Reference arm: %s | Groups: %s\n",
      x$arm_analysis$reference_arm,
      paste(x$arm_analysis$levels, collapse = ", ")
    ))

    if (!is.null(x$arm_analysis$descriptives) &&
        nrow(x$arm_analysis$descriptives) > 0L) {

      ad <- x$arm_analysis$descriptives

      ad_print <- data.frame(
        Arm = ad$arm,
        Time = ad$time_label,
        Group_N = ad$group_n,
        Available_N = ad$n,
        Unavailable_pct = round(ad$unavailable_pct, 1L),
        Mean = round(ad$mean, digits),
        SD = round(ad$sd, digits),
        Median = round(ad$median, digits),
        CI_low = round(ad$ci_lower, digits),
        CI_high = round(ad$ci_upper, digits),
        check.names = FALSE
      )

      cat("\nDescriptives by arm and time:\n")
      print(limit_table(ad_print, "time-by-arm rows"), row.names = FALSE)
    }

    if (!is.null(x$arm_analysis$missingness) &&
        nrow(x$arm_analysis$missingness) > 0L) {

      am <- x$arm_analysis$missingness

      am_print <- data.frame(
        Time = am$time_label,
        Test = am$test,
        p = vapply(am$p_value, fmt_p, character(1L)),
        p_adj = vapply(am$p_value_adj, fmt_p, character(1L)),
        check.names = FALSE
      )

      cat("\nDifferential missingness across arms:\n")
      print(am_print, row.names = FALSE)
    }

    if (!is.null(x$arm_analysis$baseline_balance) &&
        nrow(x$arm_analysis$baseline_balance) > 0L) {

      bb <- x$arm_analysis$baseline_balance

      bb_print <- data.frame(
        Arm_A = bb$arm_a,
        Arm_B = bb$arm_b,
        N_A = bb$n_a,
        N_B = bb$n_b,
        Mean_A = round(bb$mean_a, digits),
        Mean_B = round(bb$mean_b, digits),
        Difference_B_minus_A = round(bb$mean_difference_b_minus_a, digits),
        CI_low = round(bb$ci_lower, digits),
        CI_high = round(bb$ci_upper, digits),
        Hedges_g = round(bb$hedges_g, digits),
        Welch_p = vapply(bb$welch_t_p, fmt_p, character(1L)),
        Wilcoxon_p = vapply(bb$wilcoxon_p, fmt_p, character(1L)),
        check.names = FALSE
      )

      cat("\nBaseline balance (exploratory; arm B - arm A):\n")
      print(limit_table(bb_print, "baseline comparisons"), row.names = FALSE)
    }

    if (!is.null(x$arm_analysis$time_omnibus) &&
        nrow(x$arm_analysis$time_omnibus) > 0L) {

      ao <- x$arm_analysis$time_omnibus

      ao_print <- data.frame(
        Time = ao$time_label,
        N = ao$n,
        Welch_ANOVA_p = vapply(ao$welch_anova_p, fmt_p, character(1L)),
        Welch_ANOVA_p_adj = vapply(ao$welch_anova_p_adj, fmt_p, character(1L)),
        Kruskal_p = vapply(ao$kruskal_p, fmt_p, character(1L)),
        Kruskal_p_adj = vapply(ao$kruskal_p_adj, fmt_p, character(1L)),
        check.names = FALSE
      )

      cat("\nOmnibus arm comparison at each timepoint:\n")
      print(ao_print, row.names = FALSE)
    }

    if (!is.null(x$arm_analysis$time_pairwise) &&
        nrow(x$arm_analysis$time_pairwise) > 0L) {

      ap <- x$arm_analysis$time_pairwise

      ap_print <- data.frame(
        Time = ap$time_label,
        Arm_A = ap$arm_a,
        Arm_B = ap$arm_b,
        N_A = ap$n_a,
        N_B = ap$n_b,
        Difference_B_minus_A = round(ap$mean_difference_b_minus_a, digits),
        CI_low = round(ap$ci_lower, digits),
        CI_high = round(ap$ci_upper, digits),
        Hedges_g = round(ap$hedges_g, digits),
        Welch_p = vapply(ap$welch_t_p, fmt_p, character(1L)),
        Welch_p_adj = vapply(ap$welch_t_p_adj, fmt_p, character(1L)),
        Wilcoxon_p = vapply(ap$wilcoxon_p, fmt_p, character(1L)),
        Wilcoxon_p_adj = vapply(ap$wilcoxon_p_adj, fmt_p, character(1L)),
        check.names = FALSE
      )

      cat("\nPairwise arm comparisons by timepoint:\n")
      print(limit_table(ap_print, "arm pairwise comparisons"), row.names = FALSE)
    }

    if (!is.null(x$arm_analysis$change_descriptives) &&
        nrow(x$arm_analysis$change_descriptives) > 0L) {

      acd <- x$arm_analysis$change_descriptives

      acd_print <- data.frame(
        Arm = acd$arm,
        Follow_up = acd$to_label,
        N = acd$n,
        Mean_change = round(acd$mean_change, digits),
        SD_change = round(acd$sd_change, digits),
        SE_change = round(acd$se_change, digits),
        Median_change = round(acd$median_change, digits),
        CI_low = round(acd$ci_lower, digits),
        CI_high = round(acd$ci_upper, digits),
        Stable_pct = round(acd$stable_pct, 1L),
        check.names = FALSE
      )

      if (x$settings$improvement_direction != "unknown") {
        acd_print$Improved_pct <- round(acd$improved_pct, 1L)
        acd_print$Worsened_pct <- round(acd$worsened_pct, 1L)
        acd_print <- acd_print[c(
          "Arm", "Follow_up", "N", "Mean_change", "SD_change", "SE_change",
          "Median_change", "CI_low", "CI_high", "Improved_pct", "Stable_pct",
          "Worsened_pct"
        )]
      }

      cat("\nChange from baseline by arm and follow-up:\n")
      print(limit_table(acd_print, "arm change descriptive rows"), row.names = FALSE)
    }

    if (!is.null(x$arm_analysis$change_omnibus) &&
        nrow(x$arm_analysis$change_omnibus) > 0L) {

      co <- x$arm_analysis$change_omnibus

      co_print <- data.frame(
        Follow_up = co$to_label,
        N = co$n,
        Welch_ANOVA_p = vapply(co$welch_anova_p, fmt_p, character(1L)),
        Welch_ANOVA_p_adj = vapply(co$welch_anova_p_adj, fmt_p, character(1L)),
        Kruskal_p = vapply(co$kruskal_p, fmt_p, character(1L)),
        Kruskal_p_adj = vapply(co$kruskal_p_adj, fmt_p, character(1L)),
        check.names = FALSE
      )

      cat("\nOmnibus between-arm comparison of baseline-to-follow-up change:\n")
      print(co_print, row.names = FALSE)
    }

    if (!is.null(x$arm_analysis$change_pairwise) &&
        nrow(x$arm_analysis$change_pairwise) > 0L) {

      cp <- x$arm_analysis$change_pairwise

      cp_print <- data.frame(
        Follow_up = cp$to_label,
        Arm_A = cp$arm_a,
        Arm_B = cp$arm_b,
        N_A = cp$n_a,
        N_B = cp$n_b,
        Mean_change_A = round(cp$mean_change_a, digits),
        Mean_change_B = round(cp$mean_change_b, digits),
        Difference_in_change_B_minus_A =
          round(cp$difference_in_change_b_minus_a, digits),
        CI_low = round(cp$ci_lower, digits),
        CI_high = round(cp$ci_upper, digits),
        Hedges_g = round(cp$hedges_g, digits),
        Welch_p = vapply(cp$welch_t_p, fmt_p, character(1L)),
        Welch_p_adj = vapply(cp$welch_t_p_adj, fmt_p, character(1L)),
        Wilcoxon_p = vapply(cp$wilcoxon_p, fmt_p, character(1L)),
        Wilcoxon_p_adj = vapply(cp$wilcoxon_p_adj, fmt_p, character(1L)),
        check.names = FALSE
      )

      cat("\nBetween-arm comparison of baseline-to-follow-up change:\n")
      print(limit_table(cp_print, "change comparisons"), row.names = FALSE)
    }

    cat(
      "\nInterpretation note: these are exploratory frequentist summaries; ",
      "the primary longitudinal treatment inference can remain the Bayesian MIRA model.\n",
      sep = ""
    )
  }

  # ------------------------------------------------------------------
  # MISSINGNESS
  # ------------------------------------------------------------------
  if (missingness && !is.null(x$missing)) {
    section("MISSINGNESS / DATA AVAILABILITY")

    if (!is.null(x$missing$by_time) && nrow(x$missing$by_time) > 0L) {
      mbt <- x$missing$by_time
      mbt_print <- data.frame(
        Time = mbt$label,
        Missing_n = mbt$missing_n,
        Missing_pct = round(mbt$missing_pct, 1L),
        Non_finite_n = mbt$non_finite_n,
        Non_finite_pct = round(mbt$non_finite_pct, 1L),
        Unavailable_n = mbt$unavailable_n,
        Unavailable_pct = round(mbt$unavailable_pct, 1L),
        check.names = FALSE
      )
      cat("By timepoint:\n")
      print(mbt_print, row.names = FALSE)
    }

    if (!is.null(x$missing$by_patient) && nrow(x$missing$by_patient) > 0L) {
      mbp <- x$missing$by_patient
      affected <- mbp[mbp$unavailable_n > 0L | mbp$missing_n > 0L | mbp$non_finite_n > 0L, , drop = FALSE]
      cat(sprintf("\nSubjects with >=1 unavailable longitudinal value: %d/%d (%s)\n",
                  sum(mbp$unavailable_n > 0L), nrow(mbp),
                  fmt_pct(mean(mbp$unavailable_n > 0L) * 100)))
      cat(sprintf("Total unavailable cells: %d of %d (%s)\n",
                  sum(mbp$unavailable_n), nrow(mbp) * ov$n_timepoints,
                  fmt_pct(sum(mbp$unavailable_n) / (nrow(mbp) * ov$n_timepoints) * 100)))

      if (nrow(affected) == 0L) {
        cat("No subjects have missing or non-finite longitudinal values.\n")
      } else {
        cat("\nAffected subjects:\n")
        affected <- limit_table(affected, "affected subjects")
        print(affected, row.names = FALSE)
      }
    }
  }

  # ------------------------------------------------------------------
  # LONGITUDINAL CHANGE
  # ------------------------------------------------------------------
  section("LONGITUDINAL CHANGE")

  if (is.null(x$change) || nrow(x$change) == 0L) {
    cat("No longitudinal comparisons available.\n")
  } else {
    ch <- x$change
    baseline_final <- ch[
      ch$from == ov$timepoints[[1L]] &
        ch$to == ov$timepoints[[length(ov$timepoints)]],
      ,
      drop = FALSE
    ]

    if (nrow(baseline_final) == 1L) {
      bf <- baseline_final[1L, ]
      cat(sprintf("Baseline -> Final: %s -> %s (N=%d)\n",
                  bf$from_label, bf$to_label, bf$n))
      cat(sprintf("Mean: %s -> %s | Mean change: %s | SD change: %s | SE change: %s\n",
                  fmt_num(bf$mean_from), fmt_num(bf$mean_to), fmt_num(bf$mean_change),
                  fmt_num(bf$sd_change), fmt_num(bf$se_change)))
      cat(sprintf("%s: [%s, %s] | Cohen dz: %s (%s)\n",
                  conf_label, fmt_num(bf$ci_lower), fmt_num(bf$ci_upper),
                  fmt_num(bf$cohens_dz), effect_label(bf$cohens_dz)))
      cat(sprintf("Increase / decrease / stable: %s / %s / %s\n",
                  fmt_pct(bf$increased_pct), fmt_pct(bf$decreased_pct), fmt_pct(bf$stable_pct)))

      if (x$settings$improvement_direction != "unknown") {
        cat(sprintf("Improved / worsened / stable: %s / %s / %s\n",
                    fmt_pct(bf$improved_pct), fmt_pct(bf$worsened_pct), fmt_pct(bf$stable_pct)))
      } else {
        cat("Clinical improvement/worsening not classified because improvement_direction='unknown'.\n")
      }

      cat(sprintf("Paired t p: %s (adjusted %s) | Wilcoxon p: %s (adjusted %s)\n",
                  fmt_p(bf$paired_t_p), fmt_p(bf$paired_t_p_adj),
                  fmt_p(bf$wilcoxon_p), fmt_p(bf$wilcoxon_p_adj)))
    }

    pairwise_print <- data.frame(
      From = ch$from_label,
      To = ch$to_label,
      N = ch$n,
      Mean_from = round(ch$mean_from, digits),
      Mean_to = round(ch$mean_to, digits),
      Mean_change = round(ch$mean_change, digits),
      SD_change = round(ch$sd_change, digits),
      SE_change = round(ch$se_change, digits),
      CI_low = round(ch$ci_lower, digits),
      CI_high = round(ch$ci_upper, digits),
      Cohen_dz = round(ch$cohens_dz, digits),
      Increase_pct = round(ch$increased_pct, 1L),
      Decrease_pct = round(ch$decreased_pct, 1L),
      Stable_pct = round(ch$stable_pct, 1L),
      t_p = vapply(ch$paired_t_p, fmt_p, character(1L)),
      t_p_adj = vapply(ch$paired_t_p_adj, fmt_p, character(1L)),
      Wilcoxon_p = vapply(ch$wilcoxon_p, fmt_p, character(1L)),
      Wilcoxon_p_adj = vapply(ch$wilcoxon_p_adj, fmt_p, character(1L)),
      check.names = FALSE
    )

    if (x$settings$improvement_direction != "unknown") {
      pairwise_print$Improved_pct <- round(ch$improved_pct, 1L)
      pairwise_print$Worsened_pct <- round(ch$worsened_pct, 1L)
    }

    cat("\nALL PAIRWISE COMPARISONS\n")
    pairwise_print <- limit_table(pairwise_print, "pairwise comparisons")
    print(pairwise_print, row.names = FALSE)
  }

  # ------------------------------------------------------------------
  # CORRELATIONS
  # ------------------------------------------------------------------
  if (correlations) {
    section("CORRELATIONS")

    if (!is.null(x$correlations$pearson)) {
      cat("Pearson correlation matrix:\n")
      print(round(x$correlations$pearson, digits))
    } else {
      cat("Pearson correlations not available.\n")
    }

    if (!is.null(x$correlations$spearman)) {
      cat("\nSpearman correlation matrix:\n")
      print(round(x$correlations$spearman, digits))
    } else {
      cat("\nSpearman correlations not available.\n")
    }

    if (!is.null(x$correlations$pairwise_n)) {
      cat("\nPairwise N matrix:\n")
      print(x$correlations$pairwise_n)
    }
  }

  # ------------------------------------------------------------------
  # VARIABILITY / ICC
  # ------------------------------------------------------------------
  section("WITHIN- / BETWEEN-SUBJECT VARIABILITY")
  v <- x$variability[1L, ]
  ratio <- if (is.finite(v$between_subject_sd) && v$between_subject_sd != 0) {
    v$within_subject_sd / v$between_subject_sd
  } else {
    NA_real_
  }
  cat(sprintf("Grand mean: %s\n", fmt_num(v$grand_mean)))
  cat(sprintf("Between-subject SD: %s\n", fmt_num(v$between_subject_sd)))
  cat(sprintf("Within-subject SD:  %s\n", fmt_num(v$within_subject_sd)))
  cat(sprintf("Within / Between ratio: %s\n", fmt_num(ratio)))
  cat(sprintf("ICC (preferred available estimate): %s\n", fmt_num(v$ICC)))
  cat(sprintf("Model-based ICC: %s\n", fmt_num(v$ICC_model)))
  cat(sprintf("Complete-profile ANOVA ICC: %s\n", fmt_num(v$ICC_anova_complete_profiles)))
  if (all(c("ms_between", "ms_within", "df_between", "df_within") %in% names(v))) {
    cat(sprintf("ANOVA components - MS between: %s | MS within: %s | df: %s / %s\n",
                fmt_num(v$ms_between), fmt_num(v$ms_within),
                fmt_num(v$df_between, 0L), fmt_num(v$df_within, 0L)))
  }

  # ------------------------------------------------------------------
  # MIXED MODEL
  # ------------------------------------------------------------------
  if (model) {
    section("MIXED-EFFECTS MODEL")

    if (!is.null(x$model$error)) {
      cat("Model not estimated: ", x$model$error, "\n", sep = "")
    } else if (is.null(x$model$fitted_model)) {
      cat("Mixed model not available.\n")
    } else {
      fit <- x$model$fitted_model
      model_formula <- tryCatch(
        paste(deparse(stats::formula(fit)), collapse = " "),
        error = function(e) NA_character_
      )
      if (!is.na(model_formula) && nzchar(model_formula)) {
        cat("Formula: ", model_formula, "\n", sep = "")
      }
      cat(sprintf("Covariates requested: %s\n",
                  fmt_text(x$model$covariates_requested)))
      cat(sprintf("Covariates retained:  %s\n",
                  fmt_text(x$model$covariates_used)))
      if (length(x$model$covariates_skipped) > 0L) {
        cat(sprintf("Covariates skipped (insufficient variation): %s\n",
                    paste(x$model$covariates_skipped, collapse = ", ")))
      }
      if (length(x$model$fixed_parameters) == 1L &&
          is.finite(x$model$fixed_parameters)) {
        cat(sprintf("Fixed-effect parameters: %d\n", x$model$fixed_parameters))
      }
      cat(sprintf("Converged: %s | Singular: %s\n",
                  as.character(x$model$converged), as.character(x$model$singular)))

      nobs_fit <- tryCatch(stats::nobs(fit), error = function(e) NA_integer_)
      ll <- tryCatch(as.numeric(stats::logLik(fit)), error = function(e) NA_real_)
      aic <- tryCatch(stats::AIC(fit), error = function(e) NA_real_)
      bic <- tryCatch(stats::BIC(fit), error = function(e) NA_real_)
      cat(sprintf("Observations used: %s | logLik: %s | AIC: %s | BIC: %s\n",
                  ifelse(is.na(nobs_fit), "NA", as.character(nobs_fit)),
                  fmt_num(ll), fmt_num(aic), fmt_num(bic)))

      if (length(x$model$warnings) > 0L) {
        cat("\nModel warnings:\n")
        cat(paste0("  - ", x$model$warnings, collapse = "\n"), "\n")
      }

      coef_table <- tryCatch(as.data.frame(coef(summary(fit))), error = function(e) NULL)
      if (!is.null(coef_table)) {
        cat("\nFixed effects:\n")
        coef_table$Term <- rownames(coef_table)
        rownames(coef_table) <- NULL
        coef_table <- coef_table[c("Term", setdiff(names(coef_table), "Term"))]
        numeric_cols <- vapply(coef_table, is.numeric, logical(1L))
        coef_table[numeric_cols] <- lapply(coef_table[numeric_cols], round, digits = digits)
        print(coef_table, row.names = FALSE)
      }

      if (!is.null(x$model$anova)) {
        cat("\nREML model ANOVA / fixed-effect test:\n")
        print(x$model$anova)
      }

      cat("\nFull fitted model object: $model$fitted_model\n")
      cat("Full model summary:       $model$summary\n")
    }
  }

  # ------------------------------------------------------------------
  # GLOBAL ARM / TIME TESTS
  # ------------------------------------------------------------------

  if (model && isTRUE(ov$arm_analysis)) {

    section("GLOBAL TIME-BY-ARM TESTS")

    if (!is.null(x$model$arm_time_interaction_test)) {
      cat("Likelihood-ratio test of the time-by-arm interaction (interaction model vs additive model):\n")
      print(x$model$arm_time_interaction_test)
    } else {
      cat("Time-by-arm interaction test not available.\n")
    }

    if (!is.null(x$model$global_arm_test)) {
      cat("\nGlobal contribution of arm (model without arm vs full time-by-arm model):\n")
      print(x$model$global_arm_test)
    }

    if (!is.null(x$model$global_time_test)) {
      cat("\nGlobal contribution of time (model without time vs full time-by-arm model):\n")
      print(x$model$global_time_test)
    }
  }

  # ------------------------------------------------------------------
  # GLOBAL TIME TEST
  # ------------------------------------------------------------------
  if (model && !isTRUE(ov$arm_analysis)) {
    section("GLOBAL TEST OF TIME")
    if (!is.null(x$model$global_time_test)) {
      cat("Likelihood-ratio test: ML full model with time vs. random-intercept model without time.\n")
      print(x$model$global_time_test)
    } else if (!is.null(x$model$error)) {
      cat("Global test unavailable because the mixed model was not estimated.\n")
    } else {
      cat("Global time test not available.\n")
    }
  }

  # ------------------------------------------------------------------
  # ADVANCED LONGITUDINAL INFERENCE
  # ------------------------------------------------------------------
  if (model && (!is.null(x$advanced_tests) || !is.null(x$advanced_models))) {
    section("ADVANCED LONGITUDINAL INFERENCE")

    rm <- x$advanced_tests$rm_anova
    cat("RM-ANOVA:\n")
    if (!is.null(rm) && isTRUE(rm$performed) && !is.null(rm$tidy)) {
      for (question in c("TIME", "ARM", "TIME_X_ARM")) {
        idx <- .mira_match_effect(rm$tidy, question)
        if (length(idx) == 0L) next
        row <- rm$tidy[idx[[1L]], , drop = FALSE]
        gg <- if ("p_gg" %in% names(row)) fmt_p(row$p_gg) else "NA"
        hf <- if ("p_hf" %in% names(row)) fmt_p(row$p_hf) else "NA"
        cat(sprintf("  %s: F=%s | p=%s | GG p=%s | HF p=%s\n",
                    question, fmt_num(row$statistic), fmt_p(row$p_value), gg, hf))
      }
      cat("  Full object: $advanced_tests$rm_anova$object\n")
    } else {
      reason <- if (!is.null(rm$reason_skipped)) rm$reason_skipped else "not available"
      cat("  Not available: ", reason, "\n", sep = "")
    }

    friedman <- x$advanced_tests$friedman
    cat("\nFriedman:\n")
    if (!is.null(friedman) && isTRUE(friedman$performed)) {
      row <- friedman$tidy[1L, ]
      cat(sprintf("  statistic=%s | df=%s | p=%s | Kendall W=%s\n",
                  fmt_num(row$statistic), fmt_num(row$df, 0L),
                  fmt_p(row$p_raw), fmt_num(row$kendalls_w)))
      cat("  Full test: $advanced_tests$friedman$test\n")
    } else {
      reason <- if (!is.null(friedman$reason_skipped)) {
        friedman$reason_skipped
      } else "not available"
      cat("  Not available: ", reason, "\n", sep = "")
    }

    cat("\nGEE robust Wald inference:\n")
    gee <- x$advanced_models$gee
    for (structure in c("independence", "exchangeable", "ar1")) {
      item <- gee[[structure]]
      if (!is.null(item) && isTRUE(item$performed)) {
        time_row <- item$effect_tests[item$effect_tests$question == "TIME", , drop = FALSE]
        arm_row <- item$effect_tests[item$effect_tests$question == "ARM", , drop = FALSE]
        interaction_row <- item$effect_tests[
          item$effect_tests$question == "TIME_X_ARM", , drop = FALSE
        ]
        gee_tests <- c(
          if (nrow(time_row) > 0L) paste0("TIME p=", fmt_p(time_row$p_raw[[1L]])),
          if (nrow(arm_row) > 0L) paste0("ARM p=", fmt_p(arm_row$p_raw[[1L]])),
          if (nrow(interaction_row) > 0L) {
            paste0("TIME-BY-ARM p=", fmt_p(interaction_row$p_raw[[1L]]))
          }
        )
        if (length(gee_tests) == 0L) gee_tests <- "no estimable global effects"
        cat(sprintf("  %-12s %s\n", structure, paste(gee_tests, collapse = " | ")))
      } else {
        cat(sprintf("  %-12s not available\n", structure))
      }
    }
    cat("  Models: $advanced_models$gee\n")

    cat("\nMixed models:\n")
    random_slope <- x$advanced_models$random_slope
    cat(sprintf("  random intercept: %s\n",
                if (!is.null(x$model$fitted_model)) "available" else "not available"))
    if (!is.null(random_slope) && isTRUE(random_slope$performed)) {
      cat(sprintf("  random slope: available | converged=%s | singular=%s\n",
                  as.character(random_slope$converged), as.character(random_slope$singular)))
      cat("  Full object: $advanced_models$random_slope$model\n")
    } else {
      cat("  random slope: not available\n")
    }

    cat("\nCorrelation structures:\n")
    nlme_models <- x$advanced_models$nlme
    for (name in c("compound_symmetry", "ar1")) {
      item <- nlme_models[[name]]
      label <- if (name == "compound_symmetry") "CS" else "AR1"
      if (!is.null(item) && isTRUE(item$performed)) {
        cat(sprintf("  %s: available | AIC=%s | BIC=%s\n",
                    label, fmt_num(item$AIC), fmt_num(item$BIC)))
      } else {
        cat(sprintf("  %s: not available\n", label))
      }
    }
    cat("  Models: $advanced_models$nlme\n")

    comparison <- x$advanced_models$model_comparison
    if (has_rows(comparison) && "object_path" %in% names(comparison)) {
      available_models <- comparison[
        !is.na(comparison$object_path) & nzchar(comparison$object_path),
        ,
        drop = FALSE
      ]
      if (nrow(available_models) > 0L) {
        comparison_print <- data.frame(
          Model = available_models$model,
          N = available_models$n,
          Subjects = available_models$n_subjects,
          AIC = round(available_models$AIC, digits),
          BIC = round(available_models$BIC, digits),
          QIC = round(available_models$QIC, digits),
          CIC = round(available_models$CIC, digits),
          Converged = available_models$converged,
          Singular = available_models$singular,
          check.names = FALSE
        )
        optional_columns <- c("AIC", "BIC", "QIC", "CIC", "Converged", "Singular")
        keep_optional <- optional_columns[vapply(
          comparison_print[optional_columns],
          function(column) any(!is.na(column)),
          logical(1L)
        )]
        comparison_print <- comparison_print[c("Model", "N", "Subjects", keep_optional)]
        cat("\nAvailable-model comparison:\n")
        print(comparison_print, row.names = FALSE)
        comparison_note <- attr(comparison, "note")
        if (!is.null(comparison_note) && nzchar(comparison_note)) {
          cat("Note: ", comparison_note, "\n", sep = "")
        }
        cat("Full comparison and model paths: $advanced_models$model_comparison\n")
      }
    }

    cat("\nRobust mixed-model inference:\n")
    robust <- x$robustness$club_sandwich
    if (!is.null(robust) && isTRUE(robust$performed)) {
      cat("  CR2 covariance with Satterthwaite/HTZ tests available.\n")
      cat("  Full results: $robustness$club_sandwich\n")
    } else {
      cat("  Not available.\n")
    }

    emmeans <- x$advanced_tests$emmeans
    cat("\nKey marginal contrasts:\n")
    if (!is.null(emmeans) && isTRUE(emmeans$performed)) {
      emmeans_baseline_final <- emmeans$baseline_final$tidy
      baseline_followup <- emmeans$baseline_followup$tidy
      if (!is.null(emmeans_baseline_final) && nrow(emmeans_baseline_final) > 0L) {
        cat(sprintf("  baseline vs final: estimate=%s | primary-adjusted p=%s\n",
                    fmt_num(emmeans_baseline_final$estimate[[1L]]),
                    fmt_p(emmeans_baseline_final$p_primary[[1L]])))
      }
      cat(sprintf("  baseline vs follow-ups: %d contrast(s)\n",
                  if (is.null(baseline_followup)) 0L else nrow(baseline_followup)))
      cat("  Full objects: $advanced_tests$emmeans\n")
    } else {
      cat("  Not available.\n")
    }

    sensitivity <- x$sensitivity
    cat("\nSensitivity:\n")
    if (is.data.frame(sensitivity) && nrow(sensitivity) > 0L) {
      available_p <- sensitivity$p_raw[is.finite(sensitivity$p_raw)]
      range_text <- if (length(available_p) > 0L) {
        paste0(fmt_p(min(available_p)), " to ", fmt_p(max(available_p)))
      } else "NA"
      cat(sprintf("  %d method/effect rows; raw p-value range: %s\n",
                  nrow(sensitivity), range_text))
      sensitivity_print <- data.frame(
        Effect = sensitivity$question,
        Method = sensitivity$method,
        Statistic = round(sensitivity$statistic, digits),
        df = sensitivity$df,
        p = vapply(sensitivity$p_raw, fmt_p, character(1L)),
        Effect_size = round(sensitivity$effect_size, digits),
        Effect_size_type = sensitivity$effect_size_type,
        N = sensitivity$n,
        check.names = FALSE
      )
      print(limit_table(sensitivity_print, "sensitivity rows"), row.names = FALSE)
      cat("  Full comparison, source paths, and notes: $sensitivity\n")
    } else {
      cat("  No comparable inference rows available.\n")
    }
  }

  # ------------------------------------------------------------------
  # ADDITIONAL EFFECT SIZES
  # ------------------------------------------------------------------
  rm_effects <- if (model) x$effect_sizes$rm_anova else NULL
  rank_effects <- x$effect_sizes$paired_rank_biserial
  r2_rows <- list()

  add_r2_row <- function(label, item, object_path) {
    if (!model || is.null(item) || !isTRUE(item$performed) || !has_rows(item$tidy)) {
      return(NULL)
    }
    data.frame(
      Model = label,
      Marginal_R2 = round(item$tidy$marginal_r2[[1L]], digits),
      Conditional_R2 = round(item$tidy$conditional_r2[[1L]], digits),
      Object = object_path,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }

  if (!is.null(x$effect_sizes$mixed_models)) {
    r2_rows <- Filter(Negate(is.null), list(
      add_r2_row(
        "Random intercept",
        x$effect_sizes$mixed_models$random_intercept,
        "$effect_sizes$mixed_models$random_intercept"
      ),
      add_r2_row(
        "Random slope",
        x$effect_sizes$mixed_models$random_slope,
        "$effect_sizes$mixed_models$random_slope"
      )
    ))
  }

  if (has_rows(rm_effects) || has_rows(rank_effects) || length(r2_rows) > 0L) {
    section("ADDITIONAL EFFECT SIZES")

    if (has_rows(rm_effects)) {
      rm_effects_print <- data.frame(
        Effect = rm_effects$effect,
        Partial_eta_squared = round(rm_effects$partial_eta_squared, digits),
        Generalized_eta_squared = round(rm_effects$generalized_eta_squared, digits),
        check.names = FALSE
      )
      cat("RM-ANOVA effect sizes:\n")
      print(rm_effects_print, row.names = FALSE)
      cat("Full results: $effect_sizes$rm_anova\n")
    }

    if (has_rows(rank_effects)) {
      rank_effects_print <- data.frame(
        From = rank_effects$from_label,
        To = rank_effects$to_label,
        N = rank_effects$n,
        Rank_biserial = round(rank_effects$rank_biserial, digits),
        check.names = FALSE
      )
      cat("\nPaired rank-biserial effect sizes:\n")
      print(limit_table(rank_effects_print, "rank-biserial comparisons"), row.names = FALSE)
      cat("Full results: $effect_sizes$paired_rank_biserial\n")
    }

    if (length(r2_rows) > 0L) {
      cat("\nNakagawa R-squared:\n")
      print(do.call(rbind, r2_rows), row.names = FALSE)
    }

    cat("Paired Cohen dz values are reported under LONGITUDINAL CHANGE and stored in ",
        "$effect_sizes$paired_cohens_dz\n", sep = "")
  }

  # ------------------------------------------------------------------
  # MULTIPLICITY CONTROL
  # ------------------------------------------------------------------
  multiplicity <- x$multiplicity
  family_entries <- list()
  if (!is.null(multiplicity)) {
    if (model) {
      family_entries <- c(family_entries, list(
        list("Time pairwise contrasts", multiplicity$time_pairwise,
             "$multiplicity$time_pairwise$table"),
        list("Baseline vs follow-up", multiplicity$baseline_vs_followup,
             "$multiplicity$baseline_vs_followup$table"),
        list("Consecutive timepoints", multiplicity$consecutive_time,
             "$multiplicity$consecutive_time$table"),
        list("Ordinal trends", multiplicity$ordinal_trends,
             "$multiplicity$ordinal_trends$table"),
        list("Arm pairwise contrasts", multiplicity$arm_pairwise,
             "$multiplicity$arm_pairwise$table"),
        list("Arm within time", multiplicity$simple_effects$arm_within_time,
             "$multiplicity$simple_effects$arm_within_time$table"),
        list("Time within arm", multiplicity$simple_effects$time_within_arm,
             "$multiplicity$simple_effects$time_within_arm$table"),
        list("Interaction contrasts", multiplicity$interaction_contrasts,
             "$multiplicity$interaction_contrasts$table")
      ))
    }
    family_entries <- c(family_entries, list(
      list("Friedman post hoc", multiplicity$friedman_posthoc,
           "$multiplicity$friedman_posthoc$table")
    ))
  }

  family_rows <- Filter(Negate(is.null), lapply(family_entries, function(entry) {
    family <- entry[[2L]]
    if (is.null(family) || !has_rows(family$table)) return(NULL)
    data.frame(
      Family = entry[[1L]],
      Comparisons = nrow(family$table),
      Primary_adjustment = fmt_text(family$primary_method, "not recorded"),
      Table = entry[[3L]],
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }))

  if (length(family_rows) > 0L) {
    section("MULTIPLICITY CONTROL")
    family_summary <- do.call(rbind, family_rows)
    rownames(family_summary) <- NULL
    print(family_summary, row.names = FALSE)
    available_adjustments <- unique(unlist(lapply(family_entries, function(entry) {
      family <- entry[[2L]]
      if (is.null(family) || !has_rows(family$table)) character(0) else
        family$available_adjustments
    }), use.names = FALSE))
    if (length(available_adjustments) > 0L) {
      cat("Stored adjustments: ", paste(available_adjustments, collapse = ", "), "\n", sep = "")
    }
    if (!is.null(multiplicity$note) && nzchar(multiplicity$note)) {
      cat("Note: ", multiplicity$note, "\n", sep = "")
    }
    cat("Each family also records its source object in $object_path; full details: $multiplicity\n")
  }

  # ------------------------------------------------------------------
  # OUTLIERS
  # ------------------------------------------------------------------
  if (outliers) {
    section("OUTLIER FLAGS")

    if (is.null(x$outliers$by_time)) {
      cat("Outlier detection was disabled.\n")
    } else {
      counts <- vapply(x$outliers$by_time, nrow, integer(1L))
      cat(sprintf("Timepoint IQR flags: %d total\n", sum(counts)))
      if (length(counts) > 0L) print(counts)

      flagged_time <- do.call(rbind, lapply(names(x$outliers$by_time), function(nm) {
        z <- x$outliers$by_time[[nm]]
        if (nrow(z) == 0L) return(NULL)
        z$time <- nm
        z[, c("time", setdiff(names(z), "time")), drop = FALSE]
      }))

      if (!is.null(flagged_time) && nrow(flagged_time) > 0L) {
        cat("\nDetailed timepoint flags:\n")
        flagged_time <- limit_table(flagged_time, "timepoint outlier flags")
        print(flagged_time, row.names = FALSE)
      }

      if (!is.null(x$outliers$change)) {
        change_counts <- vapply(x$outliers$change, nrow, integer(1L))
        cat(sprintf("\nChange IQR flags: %d total\n", sum(change_counts)))
        if (length(change_counts) > 0L) print(change_counts)

        flagged_change <- do.call(rbind, lapply(names(x$outliers$change), function(nm) {
          z <- x$outliers$change[[nm]]
          if (nrow(z) == 0L) return(NULL)
          z$comparison <- nm
          z[, c("comparison", setdiff(names(z), "comparison")), drop = FALSE]
        }))

        if (!is.null(flagged_change) && nrow(flagged_change) > 0L) {
          cat("\nDetailed change flags:\n")
          flagged_change <- limit_table(flagged_change, "change outlier flags")
          print(flagged_change, row.names = FALSE)
        }
      }

      cat("\nNote: IQR flags are diagnostic flags and are not automatic exclusion criteria.\n")
    }
  }

  # ------------------------------------------------------------------
  # INDIVIDUAL TRAJECTORIES
  # ------------------------------------------------------------------
  if (trajectories && !is.null(x$trajectories)) {
    section("INDIVIDUAL TRAJECTORIES")
    tr <- x$trajectories
    finite_change <- tr$absolute_change[is.finite(tr$absolute_change)]

    cat(sprintf("Subjects with baseline-final data: %d/%d (%s)\n",
                length(finite_change), nrow(tr), fmt_pct(length(finite_change) / nrow(tr) * 100)))

    if (length(finite_change) > 0L) {
      cat(sprintf("Mean individual change: %s | SD: %s | Median: %s\n",
                  fmt_num(mean(finite_change)), fmt_num(stats::sd(finite_change)),
                  fmt_num(stats::median(finite_change))))
      cat(sprintf("Range: [%s, %s]\n", fmt_num(min(finite_change)), fmt_num(max(finite_change))))
    }

    direction_tab <- table(tr$direction, useNA = "no")
    if (length(direction_tab) > 0L) {
      cat("\nNumeric direction counts:\n")
      print(direction_tab)
    }

    if (x$settings$improvement_direction != "unknown") {
      clinical_tab <- table(tr$clinical_direction, useNA = "no")
      if (length(clinical_tab) > 0L) {
        cat("\nClinical direction counts:\n")
        print(clinical_tab)
      }
    }

    trajectory_columns <- c(
      "patient",
      if (isTRUE(ov$arm_analysis)) "arm",
      "baseline", "final", "absolute_change", "relative_change_percent",
      "direction", "clinical_direction"
    )
    tr_print <- tr[trajectory_columns]
    numeric_cols <- vapply(tr_print, is.numeric, logical(1L))
    tr_print[numeric_cols] <- lapply(tr_print[numeric_cols], round, digits = digits)
    cat("\nPatient-level trajectories:\n")
    tr_print <- limit_table(tr_print, "patient trajectories")
    print(tr_print, row.names = FALSE)
  }

  # ------------------------------------------------------------------
  # DIAGNOSTICS AND AUTOMATIC ADAPTATIONS
  # ------------------------------------------------------------------
  diagnostic_warnings <- x$diagnostics$warnings
  diagnostic_warnings <- as.character(diagnostic_warnings)
  diagnostic_warnings <- unique(
    diagnostic_warnings[!is.na(diagnostic_warnings) & nzchar(diagnostic_warnings)]
  )

  adaptation <- x$diagnostics$adaptation
  is_adaptation_item <- function(item) {
    is.list(item) &&
      all(c("available_by_time", "variable_timepoints", "disabled") %in% names(item))
  }
  adaptation_items <- if (is_adaptation_item(adaptation)) {
    stats::setNames(
      list(adaptation),
      fmt_text(x$outcome_display, fmt_text(x$outcome, "Outcome"))
    )
  } else if (is.list(adaptation)) {
    adaptation
  } else {
    list()
  }
  if (length(adaptation_items) > 0L && is.null(names(adaptation_items))) {
    names(adaptation_items) <- paste0("Outcome ", seq_along(adaptation_items))
  }

  adaptation_lines <- unlist(lapply(seq_along(adaptation_items), function(i) {
    item <- adaptation_items[[i]]
    if (!is_adaptation_item(item)) return(character(0))
    disabled <- as.character(item$disabled)
    disabled <- disabled[!is.na(disabled) & nzchar(disabled)]
    if (length(disabled) == 0L) return(character(0))
    label <- names(adaptation_items)[[i]]
    sprintf("%s: %s", label, paste(disabled, collapse = "; "))
  }), use.names = FALSE)

  if (length(diagnostic_warnings) > 0L || length(adaptation_lines) > 0L) {
    section("DIAGNOSTICS AND AUTOMATIC ADAPTATIONS")
    if (length(diagnostic_warnings) > 0L) {
      cat("Configuration and detection warnings:\n")
      cat(paste0("  - ", diagnostic_warnings, collapse = "\n"), "\n")
    }
    if (length(adaptation_lines) > 0L) {
      if (length(diagnostic_warnings) > 0L) cat("\n")
      cat("Analyses automatically disabled for insufficient usable data:\n")
      cat(paste0("  - ", adaptation_lines, collapse = "\n"), "\n")
    }
    cat("Full diagnostic record: $diagnostics\n")
  }

  # ------------------------------------------------------------------
  # PLOTS
  # ------------------------------------------------------------------
  if (plots) {
    section("PLOTS AVAILABLE")
    plot_descriptions <- c(
      boxplot = "Outcome distributions by timepoint with observations and mean CIs.",
      spaghetti = "Individual trajectories with the population mean and confidence band.",
      mean_ci = "Mean outcome and confidence interval across timepoints.",
      change = "Distribution of individual baseline-to-final changes.",
      change_from_baseline = "Change-from-baseline distributions at each follow-up.",
      change_ci = "Mean changes from baseline with confidence intervals.",
      missingness = "Unavailable-observation percentage by timepoint.",
      correlation_heatmap = "Pearson correlation heatmap across repeated timepoints.",
      response = "Baseline-to-final direction or clinical-response percentages.",
      arm_mean_ci = "Mean outcome and confidence intervals over time by arm.",
      arm_boxplot = "Outcome distributions by treatment arm and timepoint.",
      arm_change = "Baseline-to-final change distributions by treatment arm.",
      arm_change_ci = "Mean changes from baseline with confidence intervals by arm.",
      arm_difference_ci = "Between-arm mean differences with confidence intervals.",
      arm_missingness = "Unavailable-observation percentage by arm and timepoint."
    )
    if (length(x$plots) == 0L) {
      cat("No plot objects available.\n")
    } else {
      available_plots <- names(x$plots)[vapply(x$plots, function(z) !is.null(z), logical(1L))]
      if (length(available_plots) == 0L) {
        cat("No plot objects available.\n")
      } else {
        for (nm in available_plots) {
          description <- unname(plot_descriptions[nm])
          if (length(description) == 0L || is.na(description)) {
            description <- "Additional available plot."
          }
          cat(sprintf("  $plots$%-24s %s\n", nm, description))
        }
      }
    }
    if (!is.null(x$plot_error)) cat("Plot note: ", x$plot_error, "\n", sep = "")
  }

  # ------------------------------------------------------------------
  # COMPLETE OUTPUT GUIDE
  # ------------------------------------------------------------------
  .mira_print_output_guide(x, title = "COMPLETE OUTPUT GUIDE")
  cat("\nUse names() on the saved object to verify its current top-level components.\n")
  cat("Set max_rows = Inf in print() to display every patient, comparison, and flag row.\n")
  line("=")
  invisible(x)
}


# ============================================================
# SUMMARY METHOD
# ============================================================

#' @export
summary.mira_info <- function(object, ...) {
  ov <- object$overview
  ch <- object$change

  baseline_final <- if (nrow(ch) > 0L) {
    ch[
      ch$from == ov$timepoints[[1L]] &
        ch$to == ov$timepoints[[length(ov$timepoints)]],
      ,
      drop = FALSE
    ]
  } else {
    ch
  }

  advanced_summary <- NULL
  if (!is.null(object$advanced_tests) || !is.null(object$advanced_models)) {
    rm <- object$advanced_tests$rm_anova
    friedman <- object$advanced_tests$friedman
    emmeans <- object$advanced_tests$emmeans
    random_slope <- object$advanced_models$random_slope
    robust <- object$robustness$club_sandwich
    gee_available <- vapply(
      c("independence", "exchangeable", "ar1"),
      function(name) isTRUE(object$advanced_models$gee[[name]]$performed),
      logical(1L)
    )
    advanced_summary <- list(
      rm_anova = if (!is.null(rm) && isTRUE(rm$performed)) rm$tidy else NULL,
      friedman = if (!is.null(friedman) && isTRUE(friedman$performed)) {
        friedman$tidy
      } else NULL,
      gee_available = gee_available,
      random_slope = list(
        performed = !is.null(random_slope) && isTRUE(random_slope$performed),
        converged = if (!is.null(random_slope)) random_slope$converged else NA,
        singular = if (!is.null(random_slope)) random_slope$singular else NA
      ),
      robust_inference = !is.null(robust) && isTRUE(robust$performed),
      baseline_final_emmean = if (!is.null(emmeans) && isTRUE(emmeans$performed)) {
        emmeans$baseline_final$tidy
      } else NULL,
      sensitivity = object$sensitivity
    )
  }

  out <- list(
    outcome = object$outcome,
    n_patients = ov$n_patients,
    n_timepoints = ov$n_timepoints,
    complete_profiles = ov$complete_profiles,
    complete_profiles_pct = ov$complete_profiles_pct,
    descriptives = object$descriptives,
    baseline_final = baseline_final,
    variability = object$variability,
    arm_analysis = object$arm_analysis,
    global_time_test = object$model$global_time_test,
    global_arm_test = object$model$global_arm_test,
    arm_time_interaction_test = object$model$arm_time_interaction_test,
    model_singular = object$model$singular,
    model_converged = object$model$converged,
    advanced = advanced_summary
  )

  class(out) <- c("summary.mira_info", "list")
  out
}

#' @method print summary.mira_info
#' @export
print.summary.mira_info <- function(x, digits = 3, ...) {
  cat("MIRA INFO\n")
  cat(strrep("=", 84L), "\n", sep = "")
  cat("\nSUMMARY\n")
  cat(strrep("-", 84L), "\n", sep = "")
  cat(sprintf("Outcome: %s\n", x$outcome))
  cat(sprintf("Subjects: %d | Timepoints: %d | Complete profiles: %.1f%%\n",
              x$n_patients, x$n_timepoints, x$complete_profiles_pct))

  if (nrow(x$baseline_final) == 1L) {
    bf <- x$baseline_final[1L, ]
    cat(sprintf("Baseline-final mean change: %.*f | adjusted paired-t p: %s\n",
                digits, bf$mean_change,
                ifelse(is.na(bf$paired_t_p_adj), "NA", format.pval(bf$paired_t_p_adj, digits = 3L))))
  }

  cat(sprintf("ICC: %s\n",
              ifelse(is.na(x$variability$ICC[[1L]]), "NA",
                     formatC(x$variability$ICC[[1L]], format = "f", digits = digits))))

  if (!is.null(x$arm_time_interaction_test)) {
    cat("Time-by-arm interaction test available in $arm_time_interaction_test\n")
  }

  if (!is.null(x$global_arm_test)) {
    cat("Global arm test available in $global_arm_test\n")
  }

  if (!is.null(x$global_time_test)) {
    cat("Global time test available in $global_time_test\n")
  }

  if (!is.null(x$advanced)) {
    cat("Advanced longitudinal inference:\n")
    cat(sprintf("  RM-ANOVA: %s | Friedman: %s | GEE structures: %d/3\n",
                if (!is.null(x$advanced$rm_anova)) "available" else "not available",
                if (!is.null(x$advanced$friedman)) "available" else "not available",
                sum(x$advanced$gee_available)))
    cat(sprintf("  Random slope: %s | Robust CR2: %s\n",
                if (isTRUE(x$advanced$random_slope$performed)) "available" else "not available",
                if (isTRUE(x$advanced$robust_inference)) "available" else "not available"))
    sensitivity_n <- if (is.data.frame(x$advanced$sensitivity)) {
      nrow(x$advanced$sensitivity)
    } else 0L
    cat(sprintf("  Sensitivity rows: %d\n", sensitivity_n))
    cat("  Full sensitivity table: $advanced$sensitivity\n")
  }

  invisible(x)
}


# ============================================================
# PLOT METHOD
# ============================================================

#' @export
plot.mira_info <- function(
    x,
    which = c(
      "boxplot",
      "spaghetti",
      "mean_ci",
      "change",
      "change_from_baseline",
      "change_ci",
      "missingness",
      "correlation_heatmap",
      "response",
      "arm_mean_ci",
      "arm_boxplot",
      "arm_change",
      "arm_change_ci",
      "arm_difference_ci",
      "arm_missingness"
    ),
    ...
) {
  which <- match.arg(which)

  if (length(x$plots) == 0L || is.null(x$plots[[which]])) {
    stop(
      sprintf(
        "Plot '%s' is not available. Run mira_info(..., plots = TRUE) with ggplot2 installed.",
        which
      ),
      call. = FALSE
    )
  }

  print(x$plots[[which]])
  invisible(x$plots[[which]])
}


# ============================================================
# MULTI-OUTCOME SUMMARY AND PLOT METHODS
# ============================================================

#' @export
summary.mira_info_multi <- function(object, ...) {
  summaries <- lapply(object$outcomes, function(result) {
    if (inherits(result, "mira_info_error")) return(result)
    summary.mira_info(result, ...)
  })
  out <- list(
    version = object$version,
    config = object$config,
    outcomes = summaries,
    failed_outcomes = object$diagnostics$failed_outcomes
  )
  class(out) <- c("summary.mira_info_multi", "list")
  out
}

#' @method print summary.mira_info_multi
#' @export
print.summary.mira_info_multi <- function(x, digits = 3, ...) {
  cat(sprintf("mira_info multi-outcome summary - %d outcomes\n", length(x$outcomes)))
  for (outcome in names(x$outcomes)) {
    cat("\n")
    value <- x$outcomes[[outcome]]
    if (inherits(value, "mira_info_error")) {
      print(value)
    } else {
      print.summary.mira_info(value, digits = digits, ...)
    }
  }
  invisible(x)
}

#' @export
plot.mira_info_multi <- function(x, outcome = NULL, which = "boxplot", ...) {
  available <- names(x$outcomes)[!vapply(x$outcomes, inherits, logical(1L),
                                         what = "mira_info_error")]
  if (length(available) == 0L) {
    stop("No outcome has available plots.", call. = FALSE)
  }
  if (is.null(outcome)) {
    if (length(available) > 1L) {
      stop(sprintf("Specify outcome=. Available values: %s.",
                   paste(available, collapse = ", ")), call. = FALSE)
    }
    outcome <- available[[1L]]
  }
  idx <- which(.mira_key(available) == .mira_key(outcome))
  if (length(idx) != 1L) {
    stop(sprintf("Outcome '%s' is not available. Available values: %s.",
                 outcome, paste(available, collapse = ", ")), call. = FALSE)
  }
  plot.mira_info(x$outcomes[[available[[idx]]]], which = which, ...)
}
