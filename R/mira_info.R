# ============================================================
# MIRA_INFO v3.0 COMPLETE REPORT
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
#   - automatic treatment-arm detection when a column named `arm` exists
#   - arm x time descriptives and baseline balance diagnostics
#   - Welch/Kruskal omnibus arm tests at each timepoint
#   - multiplicity-adjusted pairwise arm comparisons with Hedges g
#   - between-arm comparisons of baseline-to-follow-up change
#   - differential missingness tests across treatment arms
#   - mixed-effects arm x time interaction and global likelihood-ratio tests
#   - arm-specific exploratory plots
# ============================================================

mira_info <- function(data,
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
  scalar_flag(arm_tests, "arm_tests")

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
  # 0B. TREATMENT ARM DETECTION / VALIDATION
  # ----------------------------------------------------------

  # Backward compatible behaviour:
  # - if arm is NULL and a column named "arm" exists, use it automatically;
  # - if no arm column is available, all original single-cohort analyses remain available.
  if (is.null(arm) && "arm" %in% names(data)) {
    arm <- "arm"
  }

  arm_available <- !is.null(arm)
  arm_levels <- character(0)
  arm_counts <- NULL
  missing_arm_n <- 0L

  if (arm_available) {

    if (!is.character(arm) || length(arm) != 1L || is.na(arm) || !nzchar(arm)) {
      stop("arm deve essere NULL oppure il nome di una sola variabile.", call. = FALSE)
    }

    if (!arm %in% names(data)) {
      stop(sprintf("La variabile di trattamento '%s' non esiste nel dataset.", arm), call. = FALSE)
    }

    arm_raw <- data[[arm]]
    missing_arm_n <- sum(is.na(arm_raw) | !nzchar(trimws(as.character(arm_raw))))

    arm_character <- as.character(arm_raw)
    arm_character[is.na(arm_raw) | !nzchar(trimws(arm_character))] <- NA_character_

    observed_arms <- unique(arm_character[!is.na(arm_character)])

    if (length(observed_arms) < 2L) {
      warning(
        sprintf(
          "La variabile '%s' contiene meno di due gruppi osservati; le analisi tra arm saranno disabilitate.",
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
        stop("reference_arm deve identificare un singolo gruppo non vuoto.", call. = FALSE)
      }

      if (!reference_arm %in% observed_arms) {
        stop(
          sprintf(
            "reference_arm='%s' non è presente nella variabile '%s'. Gruppi osservati: %s.",
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
          "%d soggetti hanno arm mancante e saranno esclusi solo dalle analisi tra gruppi.",
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
    data.frame(
      patient = analysis_data[[id]],
      outcome = outcome_name,
      time = v,
      time_index = k,
      time_label = unname(time_labels[v]),
      arm = if (arm_available) as.character(analysis_data[[arm]]) else NA_character_,
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
    # Descriptives by arm x time
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

      if (arm_available) {
        model_data$arm_factor <- factor(
          as.character(model_data$arm),
          levels = arm_levels
        )
        model_data <- model_data[!is.na(model_data$arm_factor), , drop = FALSE]
        model_data$arm_factor <- droplevels(model_data$arm_factor)
      }

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

        model_formula <- if (arm_available && arm_tests &&
                             "arm_factor" %in% names(model_data) &&
                             nlevels(model_data$arm_factor) >= 2L) {
          value ~ time_factor * arm_factor + (1 | patient_factor)
        } else {
          value ~ time_factor + (1 | patient_factor)
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
          if (arm_available && arm_tests &&
              "arm_factor" %in% names(model_data) &&
              nlevels(model_data$arm_factor) >= 2L) {

            full_ml <- tryCatch(
              fit_with_warnings(
                lme4::lmer(
                  value ~ time_factor * arm_factor + (1 | patient_factor),
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
                  value ~ time_factor + arm_factor + (1 | patient_factor),
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
                  value ~ arm_factor + (1 | patient_factor),
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
                  value ~ time_factor + (1 | patient_factor),
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
      plot_error <- "Il pacchetto 'ggplot2' non è installato."
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
    arm_variable = if (arm_available) arm else NULL,
    reference_arm = if (arm_available) reference_arm else NULL,
    arm_levels = arm_levels,
    arm_counts = arm_counts,
    missing_arm = missing_arm_n,
    arm_analysis = arm_available && arm_tests && length(arm_levels) >= 2L
  )

  result <- list(
    call = match.call(),
    version = "3.0.0",
    settings = list(
      alpha = alpha,
      confidence_level = confidence_level,
      confidence_percent = confidence_percent,
      p_adjust_method = p_adjust_method,
      improvement_direction = improvement_direction,
      stable_threshold = stable_threshold,
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

  if (isTRUE(ov$arm_analysis)) {
    cat(sprintf(
      "Arm variable: %s | Reference arm: %s | Groups: %d\n",
      ov$arm_variable,
      ov$reference_arm,
      length(ov$arm_levels)
    ))
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
      print(limit_table(ad_print, "arm x time rows"), row.names = FALSE)
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
        Mean_A = round(bb$mean_a, digits),
        Mean_B = round(bb$mean_b, digits),
        Difference_B_minus_A = round(bb$mean_difference_b_minus_a, digits),
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
        Difference_B_minus_A = round(ap$mean_difference_b_minus_a, digits),
        CI_low = round(ap$ci_lower, digits),
        CI_high = round(ap$ci_upper, digits),
        Hedges_g = round(ap$hedges_g, digits),
        Welch_p_adj = vapply(ap$welch_t_p_adj, fmt_p, character(1L)),
        Wilcoxon_p_adj = vapply(ap$wilcoxon_p_adj, fmt_p, character(1L)),
        check.names = FALSE
      )

      cat("\nPairwise arm comparisons by timepoint:\n")
      print(limit_table(ap_print, "arm pairwise comparisons"), row.names = FALSE)
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
        Mean_change_A = round(cp$mean_change_a, digits),
        Mean_change_B = round(cp$mean_change_b, digits),
        Difference_in_change_B_minus_A =
          round(cp$difference_in_change_b_minus_a, digits),
        Hedges_g = round(cp$hedges_g, digits),
        Welch_p_adj = vapply(cp$welch_t_p_adj, fmt_p, character(1L)),
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
  # GLOBAL ARM / TIME TESTS
  # ------------------------------------------------------------------

  if (isTRUE(ov$arm_analysis)) {

    section("GLOBAL ARM × TIME TESTS")

    if (!is.null(x$model$arm_time_interaction_test)) {
      cat("Likelihood-ratio test of the arm × time interaction (interaction model vs additive model):\n")
      print(x$model$arm_time_interaction_test)
    } else {
      cat("Arm × time interaction test not available.\n")
    }

    if (!is.null(x$model$global_arm_test)) {
      cat("\nGlobal contribution of arm (model without arm vs full arm × time model):\n")
      print(x$model$global_arm_test)
    }

    if (!is.null(x$model$global_time_test)) {
      cat("\nGlobal contribution of time (model without time vs full arm × time model):\n")
      print(x$model$global_time_test)
    }
  }

  # ------------------------------------------------------------------
  # GLOBAL TIME TEST
  # ------------------------------------------------------------------
  if (!isTRUE(ov$arm_analysis)) {
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
    arm_analysis = object$arm_analysis,
    global_time_test = object$model$global_time_test,
    global_arm_test = object$model$global_arm_test,
    arm_time_interaction_test = object$model$arm_time_interaction_test,
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

  if (!is.null(x$arm_time_interaction_test)) {
    cat("Arm × time interaction test available in $arm_time_interaction_test.\n")
  }

  if (!is.null(x$global_arm_test)) {
    cat("Global arm test available in $global_arm_test.\n")
  }

  if (!is.null(x$global_time_test)) {
    cat("Global time test available in $global_time_test.\n")
  }

  invisible(x)
}


# ============================================================
# PLOT METHOD
# ============================================================

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
        "Il grafico '%s' non è disponibile. Esegui mira_info(..., plots = TRUE) con ggplot2 installato.",
        which
      ),
      call. = FALSE
    )
  }

  print(x$plots[[which]])
  invisible(x$plots[[which]])
}
