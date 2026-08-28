# ============================================================
# MIRA_INFO v2.1 COMPLETE REPORT
# Robust longitudinal analysis for wide-format repeated measures
#
# Expected longitudinal column names:
#   BCVA_t0, BCVA_t1, BCVA_t2, ...
#   CMT_t0,  CMT_t1,  CMT_t2,  ...
#
# Main improvements vs. the original version:
#   - consistent handling of Inf/-Inf as unavailable observations
#   - safe ID validation for wide longitudinal data
#   - correct reordering of manual time_labels together with time_vars
#   - generic increase/decrease/stable classification
#   - optional direction of clinical improvement (higher/lower/unknown)
#   - multiplicity-adjusted pairwise p-values (Holm by default)
#   - pairwise sample-size matrix for correlations
#   - model diagnostics (warnings, convergence, singularity)
#   - likelihood-ratio global time test when lme4 is available
#   - model-based ICC in addition to balanced-data ANOVA ICC
#   - plotting requires ggplot2 only (no hidden dplyr dependency)
#   - dynamic confidence-level labels (not hard-coded to 95%)
#   - compact print(), summary(), and plot() S3 methods
# ============================================================

mira_info <- function(data,
                      id = "patient",
                      time_vars = NULL,
                      time_labels = NULL,
                      alpha = 0.05,
                      plots = TRUE,
                      model = TRUE,
                      outliers = TRUE,
                      correlations = TRUE,
                      verbose = TRUE,
                      p_adjust_method = "holm",
                      improvement_direction = c("unknown", "higher", "lower"),
                      stable_threshold = 0,
                      strict_id = TRUE) {

  # ----------------------------------------------------------
  # 0. INPUT VALIDATION
  # ----------------------------------------------------------

  if (!is.data.frame(data)) {
    stop("data deve essere un data.frame.", call. = FALSE)
  }

  if (nrow(data) == 0L) {
    stop("Il dataset non contiene osservazioni.", call. = FALSE)
  }

  if (!is.character(id) || length(id) != 1L || is.na(id) || !nzchar(id)) {
    stop("id deve essere il nome di una sola variabile.", call. = FALSE)
  }

  if (!id %in% names(data)) {
    stop(sprintf("La variabile ID '%s' non esiste nel dataset.", id), call. = FALSE)
  }

  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) ||
      alpha <= 0 || alpha >= 1) {
    stop("alpha deve essere un numero compreso tra 0 e 1.", call. = FALSE)
  }

  scalar_flag <- function(x, name) {
    if (!is.logical(x) || length(x) != 1L || is.na(x)) {
      stop(sprintf("%s deve essere TRUE o FALSE.", name), call. = FALSE)
    }
    invisible(TRUE)
  }

  scalar_flag(plots, "plots")
  scalar_flag(model, "model")
  scalar_flag(outliers, "outliers")
  scalar_flag(correlations, "correlations")
  scalar_flag(verbose, "verbose")
  scalar_flag(strict_id, "strict_id")

  improvement_direction <- match.arg(improvement_direction)

  if (!is.numeric(stable_threshold) || length(stable_threshold) != 1L ||
      !is.finite(stable_threshold) || stable_threshold < 0) {
    stop("stable_threshold deve essere un numero finito >= 0.", call. = FALSE)
  }

  if (!is.character(p_adjust_method) || length(p_adjust_method) != 1L ||
      is.na(p_adjust_method) || !p_adjust_method %in% p.adjust.methods) {
    stop(
      sprintf(
        "p_adjust_method deve essere uno tra: %s.",
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
      sprintf("La variabile ID contiene %d valori mancanti.", missing_id_n),
      call. = FALSE
    )
  }

  if (strict_id && duplicated_id_n > 0L) {
    stop(
      sprintf(
        paste0(
          "La variabile ID contiene %d duplicati. ",
          "mira_info() assume una riga per soggetto nel formato wide. ",
          "Correggi i duplicati oppure usa strict_id = FALSE consapevolmente."
        ),
        duplicated_id_n
      ),
      call. = FALSE
    )
  }

  if (!strict_id && missing_id_n > 0L) {
    warning(sprintf("Sono presenti %d ID mancanti.", missing_id_n), call. = FALSE)
  }

  if (!strict_id && duplicated_id_n > 0L) {
    warning(
      sprintf(
        "%d ID duplicati: saranno trattati come osservazioni dello stesso soggetto.",
        duplicated_id_n
      ),
      call. = FALSE
    )
  }

  # ----------------------------------------------------------
  # 1. DETECT OUTCOME + TIMEPOINTS
  # ----------------------------------------------------------

  supplied_time_vars <- !is.null(time_vars)

  if (is.null(time_vars)) {
    matched <- grep("^(.+)_t([0-9]+)$", names(data), value = TRUE)

    if (length(matched) < 2L) {
      stop(
        paste0(
          "Non sono state trovate almeno due variabili longitudinali nel formato ",
          "OUTCOME_tTIME (es. BCVA_t0, BCVA_t1)."
        ),
        call. = FALSE
      )
    }

    outcome_names <- sub("_t[0-9]+$", "", matched)
    unique_outcomes <- unique(outcome_names)

    if (length(unique_outcomes) != 1L) {
      stop(
        paste0(
          "Sono stati rilevati più outcome: ",
          paste(unique_outcomes, collapse = ", "),
          ". Specifica time_vars per analizzarne uno alla volta."
        ),
        call. = FALSE
      )
    }

    time_numbers <- as.numeric(sub("^.+_t([0-9]+)$", "\\1", matched))
    ord <- order(time_numbers, matched)
    time_vars <- matched[ord]
    outcome_name <- unique_outcomes[[1L]]

  } else {
    if (!is.character(time_vars) || length(time_vars) < 2L || anyNA(time_vars)) {
      stop("time_vars deve contenere almeno due nomi di variabili.", call. = FALSE)
    }

    if (anyDuplicated(time_vars)) {
      stop("time_vars contiene nomi duplicati.", call. = FALSE)
    }

    missing_vars <- setdiff(time_vars, names(data))
    if (length(missing_vars) > 0L) {
      stop(
        sprintf(
          "Le seguenti variabili non esistono nel dataset: %s",
          paste(missing_vars, collapse = ", ")
        ),
        call. = FALSE
      )
    }

    invalid <- time_vars[!grepl("^(.+)_t([0-9]+)$", time_vars)]
    if (length(invalid) > 0L) {
      stop(
        sprintf(
          "Le seguenti variabili non rispettano OUTCOME_tTIME: %s",
          paste(invalid, collapse = ", ")
        ),
        call. = FALSE
      )
    }

    outcome_names <- sub("_t[0-9]+$", "", time_vars)
    unique_outcomes <- unique(outcome_names)

    if (length(unique_outcomes) != 1L) {
      stop(
        sprintf(
          "time_vars contiene più outcome: %s",
          paste(unique_outcomes, collapse = ", ")
        ),
        call. = FALSE
      )
    }

    outcome_name <- unique_outcomes[[1L]]
    time_numbers <- as.numeric(sub("^.+_t([0-9]+)$", "\\1", time_vars))
    ord <- order(time_numbers, time_vars)

    # IMPORTANT: reorder labels together with manually supplied variables.
    if (!is.null(time_labels)) {
      if (length(time_labels) != length(time_vars)) {
        stop("time_labels deve avere la stessa lunghezza di time_vars.", call. = FALSE)
      }
      time_labels <- time_labels[ord]
    }

    time_vars <- time_vars[ord]
  }

  if (anyDuplicated(as.numeric(sub("^.+_t([0-9]+)$", "\\1", time_vars)))) {
    stop("Sono presenti timepoint numerici duplicati.", call. = FALSE)
  }

  outcome_display <- toupper(outcome_name)
  detected_time_labels <- sub("^.+_(t[0-9]+)$", "\\1", time_vars)

  # ----------------------------------------------------------
  # 2. TIME VARIABLE VALIDATION
  # ----------------------------------------------------------

  non_numeric <- time_vars[
    !vapply(data[time_vars], is.numeric, logical(1L))
  ]

  if (length(non_numeric) > 0L) {
    stop(
      sprintf(
        "Le variabili longitudinali devono essere numeriche: %s",
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
      stop("time_labels deve avere la stessa lunghezza di time_vars.", call. = FALSE)
    }
    time_labels <- as.character(time_labels)
  }

  if (anyNA(time_labels) || any(!nzchar(trimws(time_labels)))) {
    stop("time_labels non può contenere NA o etichette vuote.", call. = FALSE)
  }

  if (anyDuplicated(time_labels)) {
    stop("time_labels deve contenere etichette univoche.", call. = FALSE)
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
        "%d valori Inf/-Inf nelle variabili longitudinali saranno trattati come NA nelle analisi.",
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
    data.frame(
      patient = analysis_data[[id]],
      outcome = outcome_name,
      time = v,
      time_index = k,
      time_label = unname(time_labels[v]),
      value = analysis_data[[v]],
      stringsAsFactors = FALSE
    )
  })

  long_data <- do.call(rbind, long_list)
  rownames(long_data) <- NULL
  long_data$time_label <- factor(
    long_data$time_label,
    levels = unname(time_labels[time_vars])
  )

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
  model_error <- NULL
  model_warnings <- character(0)
  model_singular <- NA
  model_converged <- NA
  icc_model <- NA_real_

  if (model) {
    if (!requireNamespace("lme4", quietly = TRUE)) {
      model_error <- "Il pacchetto 'lme4' non è installato."
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

      if (nrow(model_data) < 3L || nlevels(model_data$patient_factor) < 2L ||
          nlevels(model_data$time_factor) < 2L) {
        model_error <- "Dati insufficienti per stimare un mixed-effects model."
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

        model_formula <- value ~ time_factor + (1 | patient_factor)

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

          # Robust global time test: ML likelihood-ratio full vs. no-time model.
          global_time_test <- tryCatch({
            full_ml <- fit_with_warnings(
              lme4::lmer(
                value ~ time_factor + (1 | patient_factor),
                data = model_data,
                REML = FALSE,
                na.action = stats::na.omit
              )
            )
            null_ml <- fit_with_warnings(
              lme4::lmer(
                value ~ 1 + (1 | patient_factor),
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
  # 15. PLOTS (ggplot2 ONLY)
  # ----------------------------------------------------------

  plots_list <- list()
  plot_error <- NULL

  if (plots) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      plot_error <- "Il pacchetto 'ggplot2' non è installato."
    } else {
      ci_text <- paste0(formatC(confidence_percent, format = "fg", digits = 4), "% CI")

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
        ggplot2::labs(
          x = NULL,
          y = outcome_display,
          title = paste("Distribution of", outcome_display, "Across Timepoints"),
          subtitle = paste("Individual observations, IQR, mean and", ci_text)
        ) +
        ggplot2::theme_minimal(base_size = 12)

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
        ggplot2::labs(
          x = NULL,
          y = outcome_display,
          title = paste("Individual", outcome_display, "Longitudinal Trajectories"),
          subtitle = paste("Subject trajectories, population mean and", ci_text)
        ) +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(legend.position = "bottom")

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
        ggplot2::labs(
          x = "Time", y = paste("Mean", outcome_display),
          title = paste("Mean", outcome_display, "and", ci_text)
        ) +
        ggplot2::theme_minimal()

      trajectory_plot_data <- trajectory_summary[
        is.finite(trajectory_summary$absolute_change),
        ,
        drop = FALSE
      ]

      plots_list$change <-
        ggplot2::ggplot(
          trajectory_plot_data,
          ggplot2::aes(x = "All patients", y = absolute_change)
        ) +
        ggplot2::geom_boxplot(fill = "#59A14F", alpha = 0.7, na.rm = TRUE) +
        ggplot2::geom_jitter(width = 0.08, alpha = 0.30, na.rm = TRUE) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
        ggplot2::labs(
          x = NULL,
          y = paste("Final - Baseline", outcome_display),
          title = paste("Individual", outcome_display, "Change")
        ) +
        ggplot2::theme_minimal()
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
    non_finite_values = total_non_finite
  )

  result <- list(
    call = match.call(),
    version = "2.1.0",
    settings = list(
      alpha = alpha,
      confidence_level = confidence_level,
      confidence_percent = confidence_percent,
      p_adjust_method = p_adjust_method,
      improvement_direction = improvement_direction,
      stable_threshold = stable_threshold,
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
      singular = model_singular,
      converged = model_converged,
      warnings = model_warnings,
      error = model_error
    ),
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
# PRINT METHOD
# ============================================================

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
    stop("digits deve essere un intero >= 0.", call. = FALSE)
  }
  digits <- as.integer(digits)

  if (!is.numeric(max_rows) || length(max_rows) != 1L || is.na(max_rows) || max_rows <= 0) {
    stop("max_rows deve essere > 0 oppure Inf.", call. = FALSE)
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

  effect_label <- function(d) {
    if (length(d) == 0L || is.na(d)) return(NA_character_)
    a <- abs(d)
    if (a < 0.20) "negligible" else if (a < 0.50) "small" else if (a < 0.80) "moderate" else "large"
  }

  limit_table <- function(z, label = "rows") {
    if (!is.data.frame(z) || nrow(z) == 0L) return(z)
    if (is.finite(max_rows) && nrow(z) > max_rows) {
      cat(sprintf("Showing first %d of %d %s. Use print(result, max_rows = Inf) for all.\n",
                  as.integer(max_rows), nrow(z), label))
      return(utils::head(z, as.integer(max_rows)))
    }
    z
  }

  ov <- x$overview
  d <- x$descriptives
  conf_pct <- x$settings$confidence_percent
  conf_label <- paste0(formatC(conf_pct, format = "fg", digits = 4L), "% CI")

  cat("\n")
  line("=")
  cat(sprintf("MIRA INFO v%s — COMPLETE REPORT — %s\n", x$version, ov$outcome_display))
  line("=")

  # ------------------------------------------------------------------
  # ANALYSIS SETTINGS
  # ------------------------------------------------------------------
  section("ANALYSIS SETTINGS")
  cat(sprintf("Outcome: %s | ID variable: %s\n", ov$outcome_display, ov$id))
  cat(sprintf("Confidence level: %s | alpha: %s | p-adjust: %s\n",
              conf_label, fmt_num(x$settings$alpha), x$settings$p_adjust_method))
  cat(sprintf("Improvement direction: %s | Stable threshold: %s\n",
              x$settings$improvement_direction, fmt_num(x$settings$stable_threshold)))
  cat(sprintf("Strict ID checks: %s | Non-finite handling: %s\n",
              as.character(x$settings$strict_id), x$settings$non_finite_handling))

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
    cat(sprintf("ANOVA components — MS between: %s | MS within: %s | df: %s / %s\n",
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

      cat("\nFull fitted model object: result$model$fitted_model\n")
      cat("Full model summary:       result$model$summary\n")
    }
  }

  # ------------------------------------------------------------------
  # GLOBAL TIME TEST
  # ------------------------------------------------------------------
  section("GLOBAL TEST OF TIME")
  if (!is.null(x$model$global_time_test)) {
    cat("Likelihood-ratio test: ML full model with time vs. random-intercept model without time.\n")
    print(x$model$global_time_test)
  } else if (!is.null(x$model$error)) {
    cat("Global test unavailable because the mixed model was not estimated.\n")
  } else {
    cat("Global time test not available.\n")
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

    tr_print <- tr[c(
      "patient", "baseline", "final", "absolute_change",
      "relative_change_percent", "direction", "clinical_direction"
    )]
    numeric_cols <- vapply(tr_print, is.numeric, logical(1L))
    tr_print[numeric_cols] <- lapply(tr_print[numeric_cols], round, digits = digits)
    cat("\nPatient-level trajectories:\n")
    tr_print <- limit_table(tr_print, "patient trajectories")
    print(tr_print, row.names = FALSE)
  }

  # ------------------------------------------------------------------
  # PLOTS
  # ------------------------------------------------------------------
  if (plots) {
    section("PLOTS AVAILABLE")
    if (length(x$plots) == 0L) {
      cat("No plot objects available.\n")
    } else {
      available_plots <- names(x$plots)[vapply(x$plots, function(z) !is.null(z), logical(1L))]
      if (length(available_plots) == 0L) {
        cat("No plot objects available.\n")
      } else {
        for (nm in available_plots) {
          cat(sprintf("  %-12s -> plot(result, which = \"%s\")  /  result$plots$%s\n", nm, nm, nm))
        }
      }
    }
    if (!is.null(x$plot_error)) cat("Plot note: ", x$plot_error, "\n", sep = "")
  }

  # ------------------------------------------------------------------
  # COMPLETE OUTPUT GUIDE
  # ------------------------------------------------------------------
  section("COMPLETE OUTPUT GUIDE")
  cat(sprintf("  $outcome        Detected outcome (%s)\n", ov$outcome_display))
  cat("  $overview       Dataset structure, IDs, completeness and timepoints\n")
  cat("  $settings       Confidence level, p-adjustment and clinical direction settings\n")
  cat("  $descriptives   Detailed statistics for every timepoint\n")
  cat("  $missing        Missing/non-finite/unavailable data by timepoint and subject\n")
  cat("  $change         All longitudinal pairwise comparisons and adjusted tests\n")
  cat("  $correlations   Pearson, Spearman and pairwise sample-size matrices\n")
  cat("  $variability    Within/between-subject variability and ICC estimates\n")
  cat("  $trajectories   Subject-level baseline-to-final changes and direction\n")
  cat("  $model          Mixed-effects model, diagnostics, ANOVA and global time test\n")
  cat("  $outliers       Timepoint and change IQR diagnostic flags\n")
  cat("  $plots          ggplot objects generated by mira_info()\n")
  cat("  $long_data      Long-format analysis dataset\n")
  cat("\nUse names(result) for all top-level components.\n")
  cat("Use print(result, max_rows = Inf) to print every patient/comparison/flag.\n")
  line("=")
  invisible(x)
}


# ============================================================
# SUMMARY METHOD
# ============================================================

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

  out <- list(
    outcome = object$outcome,
    n_patients = ov$n_patients,
    n_timepoints = ov$n_timepoints,
    complete_profiles = ov$complete_profiles,
    complete_profiles_pct = ov$complete_profiles_pct,
    descriptives = object$descriptives,
    baseline_final = baseline_final,
    variability = object$variability,
    global_time_test = object$model$global_time_test,
    model_singular = object$model$singular,
    model_converged = object$model$converged
  )

  class(out) <- c("summary.mira_info", "list")
  out
}


print.summary.mira_info <- function(x, digits = 3, ...) {
  cat(sprintf("mira_info summary — %s\n", toupper(x$outcome)))
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

  if (!is.null(x$global_time_test)) {
    cat("Global time test available in $global_time_test.\n")
  }

  invisible(x)
}


# ============================================================
# PLOT METHOD
# ============================================================

plot.mira_info <- function(x,
                           which = c("boxplot", "spaghetti", "mean_ci", "change"),
                           ...) {
  which <- match.arg(which)

  if (length(x$plots) == 0L || is.null(x$plots[[which]])) {
    stop(
      sprintf(
        "Il grafico '%s' non è disponibile. Esegui mira_info(..., plots = TRUE) con ggplot2 installato.",
        which
      ),
      call. = FALSE
    )
  }

  print(x$plots[[which]])
  invisible(x$plots[[which]])
}
