# ============================================================
# MIRA_INFO
# Comprehensive Longitudinal Data Analysis
#
# Automatic outcome detection
#
# Expected column structure:
#
#   BCVA_t0
#   BCVA_t1
#   BCVA_t2
#
# or:
#
#   CMT_t0
#   CMT_t1
#   CMT_t2
#
# or:
#
#   IOP_t0
#   IOP_t1
#   IOP_t2
#
# The function automatically detects:
#
#   outcome = BCVA / CMT / IOP / ...
#   timepoints = t0 / t1 / t2 / ...
#
# ============================================================


# ============================================================
# MAIN FUNCTION
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
                      verbose = TRUE) {


  # ----------------------------------------------------------
  # 0. CHECK INPUT
  # ----------------------------------------------------------

  if (!is.data.frame(data)) {

    stop(
      "data deve essere un data.frame."
    )
  }


  if (nrow(data) == 0) {

    stop(
      "Il dataset non contiene osservazioni."
    )
  }


  if (!is.character(id) ||
      length(id) != 1) {

    stop(
      "id deve essere il nome di una sola variabile."
    )
  }


  if (!id %in% names(data)) {

    stop(
      sprintf(
        "La variabile ID '%s' non esiste nel dataset.",
        id
      )
    )
  }


  if (!is.numeric(alpha) ||
      length(alpha) != 1 ||
      is.na(alpha) ||
      alpha <= 0 ||
      alpha >= 1) {

    stop(
      "alpha deve essere un numero compreso tra 0 e 1."
    )
  }


  # ----------------------------------------------------------
  # 1. IDENTIFICAZIONE AUTOMATICA OUTCOME + TIMEPOINT
  # ----------------------------------------------------------

  #
  # La struttura attesa è:
  #
  # OUTCOME_t0
  # OUTCOME_t1
  # OUTCOME_t2
  #
  # ecc.
  #

  if (is.null(time_vars)) {


    # --------------------------------------------------------
    # Cerca tutte le colonne che rispettano:
    #
    # qualcosa_tnumero
    #
    # Esempi:
    # BCVA_t0
    # CMT_t1
    # IOP_t2
    # --------------------------------------------------------

    matched <- grep(
      "^(.+)_t([0-9]+)$",
      names(data),
      value = TRUE
    )


    if (length(matched) == 0) {

      stop(
        paste0(
          "Non sono state trovate variabili longitudinali.\n\n",
          "Il dataset deve contenere colonne con struttura:\n\n",
          "  BCVA_t0\n",
          "  BCVA_t1\n",
          "  BCVA_t2\n\n",
          "oppure:\n\n",
          "  CMT_t0\n",
          "  CMT_t1\n",
          "  CMT_t2\n\n",
          "oppure qualsiasi altro OUTCOME_tTIME."
        )
      )
    }


    # --------------------------------------------------------
    # Estrae il nome dell'outcome
    # --------------------------------------------------------

    outcome_names <- sub(
      "_t[0-9]+$",
      "",
      matched
    )


    unique_outcomes <- unique(
      outcome_names
    )


    # --------------------------------------------------------
    # Deve esserci un solo outcome
    # --------------------------------------------------------

    if (length(unique_outcomes) > 1) {

      stop(
        paste0(
          "Sono stati rilevati più outcome nello stesso dataset:\n\n",
          paste(
            unique_outcomes,
            collapse = ", "
          ),
          "\n\n",
          "mira_info() analizza un outcome alla volta.\n",
          "Usa un dataset contenente un solo outcome longitudinale ",
          "oppure seleziona esplicitamente le colonne tramite ",
          "time_vars."
        )
      )
    }


    outcome_name <- unique_outcomes[1]


    # --------------------------------------------------------
    # Estrae il numero del timepoint
    # --------------------------------------------------------

    time_numbers <- as.numeric(
      sub(
        "^.+_t([0-9]+)$",
        "\\1",
        matched
      )
    )


    # --------------------------------------------------------
    # Ordina cronologicamente
    #
    # t0, t1, t2, ..., t10
    #
    # --------------------------------------------------------

    time_vars <- matched[
      order(time_numbers)
    ]


  } else {


    # --------------------------------------------------------
    # TIME_VARS FORNITO MANUALMENTE
    # --------------------------------------------------------

    if (!is.character(time_vars) ||
        length(time_vars) < 2) {

      stop(
        "time_vars deve contenere almeno due variabili longitudinali."
      )
    }


    # --------------------------------------------------------
    # Controlla che le variabili esistano
    # --------------------------------------------------------

    missing_vars <- setdiff(
      time_vars,
      names(data)
    )


    if (length(missing_vars) > 0) {

      stop(
        paste0(
          "Le seguenti variabili non esistono nel dataset: ",
          paste(
            missing_vars,
            collapse = ", "
          )
        )
      )
    }


    # --------------------------------------------------------
    # Controlla formato OUTCOME_tTIME
    # --------------------------------------------------------

    invalid_time_vars <- time_vars[
      !grepl(
        "^(.+)_t([0-9]+)$",
        time_vars
      )
    ]


    if (length(invalid_time_vars) > 0) {

      stop(
        paste0(
          "Le seguenti variabili non rispettano il formato ",
          "'OUTCOME_tTIME': ",
          paste(
            invalid_time_vars,
            collapse = ", "
          )
        )
      )
    }


    # --------------------------------------------------------
    # Estrae outcome
    # --------------------------------------------------------

    outcome_names <- sub(
      "_t[0-9]+$",
      "",
      time_vars
    )


    unique_outcomes <- unique(
      outcome_names
    )


    if (length(unique_outcomes) > 1) {

      stop(
        paste0(
          "time_vars contiene più outcome:\n\n",
          paste(
            unique_outcomes,
            collapse = ", "
          ),
          "\n\n",
          "mira_info() analizza un solo outcome alla volta."
        )
      )
    }


    outcome_name <- unique_outcomes[1]


    # --------------------------------------------------------
    # Ordina timepoint
    # --------------------------------------------------------

    time_numbers <- as.numeric(
      sub(
        "^.+_t([0-9]+)$",
        "\\1",
        time_vars
      )
    )


    time_vars <- time_vars[
      order(time_numbers)
    ]
  }


  # ----------------------------------------------------------
  # Nome visualizzato
  #
  # BCVA -> BCVA
  # cmt  -> CMT
  # iop  -> IOP
  # ----------------------------------------------------------

  outcome_display <- toupper(
    outcome_name
  )


  # ----------------------------------------------------------
  # Timepoint names
  #
  # BCVA_t0 -> t0
  # BCVA_t1 -> t1
  # BCVA_t2 -> t2
  # ----------------------------------------------------------

  detected_time_labels <- sub(
    "^.+_(t[0-9]+)$",
    "\\1",
    time_vars
  )


  # ----------------------------------------------------------
  # 2. CONTROLLO TIME VARIABLES
  # ----------------------------------------------------------

  missing_vars <- setdiff(
    time_vars,
    names(data)
  )


  if (length(missing_vars) > 0) {

    stop(
      paste0(
        "Le seguenti variabili non esistono nel dataset: ",
        paste(
          missing_vars,
          collapse = ", "
        )
      )
    )
  }


  non_numeric <- time_vars[
    !vapply(
      data[time_vars],
      is.numeric,
      logical(1)
    )
  ]


  if (length(non_numeric) > 0) {

    stop(
      paste0(
        "Le variabili longitudinali devono essere numeriche: ",
        paste(
          non_numeric,
          collapse = ", "
        )
      )
    )
  }


  # ----------------------------------------------------------
  # 3. TIME LABELS
  # ----------------------------------------------------------

  if (is.null(time_labels)) {

    time_labels <- detected_time_labels

  } else {

    if (length(time_labels) != length(time_vars)) {

      stop(
        "time_labels deve avere la stessa lunghezza di time_vars."
      )
    }


    if (!is.character(time_labels)) {

      time_labels <- as.character(
        time_labels
      )
    }
  }


  names(time_labels) <- time_vars


  # ----------------------------------------------------------
  # 4. HELPER FUNCTIONS
  # ----------------------------------------------------------

  mean_ci <- function(x,
                      alpha = 0.05) {

    x <- x[
      is.finite(x)
    ]

    n <- length(x)


    if (n == 0) {

      return(
        c(
          mean = NA_real_,
          se = NA_real_,
          lower = NA_real_,
          upper = NA_real_
        )
      )
    }


    m <- mean(x)


    if (n < 2) {

      return(
        c(
          mean = m,
          se = NA_real_,
          lower = NA_real_,
          upper = NA_real_
        )
      )
    }


    s <- sd(x)

    se <- s / sqrt(n)


    crit <- qt(
      1 - alpha / 2,
      df = n - 1
    )


    c(
      mean = m,
      se = se,
      lower = m - crit * se,
      upper = m + crit * se
    )
  }


  safe_mean <- function(x) {

    x <- x[
      is.finite(x)
    ]


    if (length(x) == 0) {

      return(
        NA_real_
      )
    }


    mean(x)
  }


  safe_sd <- function(x) {

    x <- x[
      is.finite(x)
    ]


    if (length(x) < 2) {

      return(
        NA_real_
      )
    }


    sd(x)
  }


  safe_median <- function(x) {

    x <- x[
      is.finite(x)
    ]


    if (length(x) == 0) {

      return(
        NA_real_
      )
    }


    median(x)
  }


  safe_quantile <- function(x,
                            p) {

    x <- x[
      is.finite(x)
    ]


    if (length(x) == 0) {

      return(
        NA_real_
      )
    }


    as.numeric(
      quantile(
        x,
        probs = p,
        names = FALSE
      )
    )
  }


  safe_var <- function(x) {

    x <- x[
      is.finite(x)
    ]


    if (length(x) < 2) {

      return(
        NA_real_
      )
    }


    var(x)
  }


  # ----------------------------------------------------------
  # 5. BASIC OVERVIEW
  # ----------------------------------------------------------

  n_rows <- nrow(data)


  patient_values <- data[[id]]


  n_patients <- length(
    unique(
      patient_values[
        !is.na(patient_values)
      ]
    )
  )


  duplicated_ids <- sum(
    duplicated(
      patient_values[
        !is.na(patient_values)
      ]
    )
  )


  # ----------------------------------------------------------
  # Complete profile
  # ----------------------------------------------------------

  complete_cases <- sum(
    complete.cases(
      data[time_vars]
    )
  )


  # ----------------------------------------------------------
  # 6. DESCRIPTIVE STATISTICS
  # ----------------------------------------------------------

  descriptive_list <- lapply(

    time_vars,

    function(v) {


      x <- data[[v]]


      finite_x <- x[
        is.finite(x)
      ]


      n <- length(
        finite_x
      )


      missing_n <- sum(
        is.na(x)
      )


      non_finite_n <- sum(
        !is.na(x) &
          !is.finite(x)
      )


      ci <- mean_ci(
        x,
        alpha = alpha
      )


      m <- safe_mean(x)

      s <- safe_sd(x)

      variance <- safe_var(x)


      q1 <- safe_quantile(
        x,
        0.25
      )


      q3 <- safe_quantile(
        x,
        0.75
      )


      cv <- if (
        !is.na(m) &&
        !is.na(s) &&
        m != 0
      ) {

        s / abs(m) * 100

      } else {

        NA_real_
      }


      data.frame(

        time = v,

        label = unname(
          time_labels[v]
        ),

        n = n,

        missing = missing_n,

        missing_pct =
          missing_n /
          n_rows *
          100,

        non_finite =
          non_finite_n,

        mean = m,

        sd = s,

        variance = variance,

        se =
          unname(
            ci["se"]
          ),

        ci_lower =
          unname(
            ci["lower"]
          ),

        ci_upper =
          unname(
            ci["upper"]
          ),

        median =
          safe_median(x),

        q1 = q1,

        q3 = q3,

        iqr =
          if (
            !is.na(q1) &&
            !is.na(q3)
          ) {

            q3 - q1

          } else {

            NA_real_
          },

        min =
          if (n > 0) {

            min(finite_x)

          } else {

            NA_real_
          },

        max =
          if (n > 0) {

            max(finite_x)

          } else {

            NA_real_
          },

        cv_percent = cv,

        stringsAsFactors = FALSE
      )
    }
  )


  descriptives <- do.call(
    rbind,
    descriptive_list
  )


  rownames(
    descriptives
  ) <- NULL


  # ----------------------------------------------------------
  # Outcome-specific descriptive column names
  #
  # BCVA_mean
  # BCVA_sd
  # BCVA_median
  # ...
  # ----------------------------------------------------------

  descriptives[[paste0(
    outcome_name,
    "_mean"
  )]] <-
    descriptives$mean


  descriptives[[paste0(
    outcome_name,
    "_sd"
  )]] <-
    descriptives$sd


  descriptives[[paste0(
    outcome_name,
    "_median"
  )]] <-
    descriptives$median


  descriptives[[paste0(
    outcome_name,
    "_ci_lower"
  )]] <-
    descriptives$ci_lower


  descriptives[[paste0(
    outcome_name,
    "_ci_upper"
  )]] <-
    descriptives$ci_upper


  # ----------------------------------------------------------
  # 7. MISSINGNESS
  # ----------------------------------------------------------

  missing_by_patient <- data.frame(

    patient =
      data[[id]],

    missing_n =
      rowSums(
        is.na(
          data[time_vars]
        )
      ),

    missing_pct =
      rowMeans(
        is.na(
          data[time_vars]
        )
      ) * 100,

    stringsAsFactors = FALSE
  )


  missing_summary <- data.frame(

    time = time_vars,

    label =
      unname(
        time_labels[time_vars]
      ),

    missing_n =
      vapply(
        data[time_vars],
        function(x) {
          sum(is.na(x))
        },
        numeric(1)
      ),

    missing_pct =
      vapply(
        data[time_vars],
        function(x) {
          mean(is.na(x)) * 100
        },
        numeric(1)
      ),

    stringsAsFactors = FALSE
  )


  # ----------------------------------------------------------
  # 8. CHANGE ANALYSIS
  # ----------------------------------------------------------

  change_list <- list()

  counter <- 1


  for (i in seq_len(
    length(time_vars) - 1
  )) {


    for (j in (i + 1):length(time_vars)) {


      v1 <- time_vars[i]

      v2 <- time_vars[j]


      x <- data[[v1]]

      y <- data[[v2]]


      keep <- complete.cases(
        x,
        y
      )


      x_complete <- x[
        keep
      ]


      y_complete <- y[
        keep
      ]


      delta <- y_complete -
        x_complete


      finite <- is.finite(
        delta
      )


      delta <-
        delta[finite]


      x_complete <-
        x_complete[finite]


      y_complete <-
        y_complete[finite]


      n <- length(delta)


      if (n == 0) {

        next
      }


      delta_mean <-
        mean(delta)


      delta_sd <-

        if (n >= 2) {

          sd(delta)

        } else {

          NA_real_
        }


      delta_se <-

        if (n >= 2) {

          delta_sd /
            sqrt(n)

        } else {

          NA_real_
        }


      if (n >= 2) {

        tcrit <- qt(
          1 - alpha / 2,
          df = n - 1
        )


        ci_low <-
          delta_mean -
          tcrit *
          delta_se


        ci_high <-
          delta_mean +
          tcrit *
          delta_se

      } else {

        ci_low <- NA_real_

        ci_high <- NA_real_
      }


      # ------------------------------------------------------
      # Cohen's dz
      # ------------------------------------------------------

      effect_size <-

        if (
          n >= 2 &&
          !is.na(delta_sd) &&
          delta_sd > 0
        ) {

          delta_mean /
            delta_sd

        } else {

          NA_real_
        }


      # ------------------------------------------------------
      # Paired t-test
      # ------------------------------------------------------

      ttest <- tryCatch(

        t.test(
          y_complete,
          x_complete,
          paired = TRUE,
          conf.level = 1 - alpha
        ),

        error = function(e) {
          NULL
        }
      )


      p_value <-

        if (!is.null(ttest)) {

          ttest$p.value

        } else {

          NA_real_
        }


      # ------------------------------------------------------
      # Wilcoxon
      # ------------------------------------------------------

      wilcox <- tryCatch(

        suppressWarnings(

          wilcox.test(
            y_complete,
            x_complete,
            paired = TRUE,
            exact = FALSE
          )
        ),

        error = function(e) {
          NULL
        }
      )


      wilcox_p <-

        if (!is.null(wilcox)) {

          wilcox$p.value

        } else {

          NA_real_
        }


      # ------------------------------------------------------
      # Direction
      # ------------------------------------------------------

      improved <- sum(
        delta > 0
      )


      worsened <- sum(
        delta < 0
      )


      unchanged <- sum(
        delta == 0
      )


      # ------------------------------------------------------
      # Store result
      # ------------------------------------------------------

      change_list[[counter]] <-

        data.frame(

          from = v1,

          to = v2,

          from_label =
            unname(
              time_labels[v1]
            ),

          to_label =
            unname(
              time_labels[v2]
            ),

          n = n,

          mean_from =
            mean(x_complete),

          mean_to =
            mean(y_complete),

          mean_change =
            delta_mean,

          sd_change =
            delta_sd,

          se_change =
            delta_se,

          ci_lower =
            ci_low,

          ci_upper =
            ci_high,

          cohens_dz =
            effect_size,

          improved_n =
            improved,

          improved_pct =
            improved /
            n *
            100,

          worsened_n =
            worsened,

          worsened_pct =
            worsened /
            n *
            100,

          unchanged_n =
            unchanged,

          unchanged_pct =
            unchanged /
            n *
            100,

          paired_t_p =
            p_value,

          wilcoxon_p =
            wilcox_p,

          stringsAsFactors = FALSE
        )


      counter <- counter + 1
    }
  }


  # ----------------------------------------------------------
  # Combine change results
  # ----------------------------------------------------------

  if (length(change_list) > 0) {

    change <- do.call(
      rbind,
      change_list
    )


    rownames(change) <- NULL

  } else {

    change <- data.frame(

      from = character(0),

      to = character(0),

      from_label = character(0),

      to_label = character(0),

      n = numeric(0),

      mean_from = numeric(0),

      mean_to = numeric(0),

      mean_change = numeric(0),

      sd_change = numeric(0),

      se_change = numeric(0),

      ci_lower = numeric(0),

      ci_upper = numeric(0),

      cohens_dz = numeric(0),

      improved_n = numeric(0),

      improved_pct = numeric(0),

      worsened_n = numeric(0),

      worsened_pct = numeric(0),

      unchanged_n = numeric(0),

      unchanged_pct = numeric(0),

      paired_t_p = numeric(0),

      wilcoxon_p = numeric(0),

      stringsAsFactors = FALSE
    )
  }


  # ----------------------------------------------------------
  # Outcome-specific change columns
  # ----------------------------------------------------------

  if (nrow(change) > 0) {

    change[[paste0(
      outcome_name,
      "_change"
    )]] <-
      change$mean_change


    change[[paste0(
      outcome_name,
      "_mean_from"
    )]] <-
      change$mean_from


    change[[paste0(
      outcome_name,
      "_mean_to"
    )]] <-
      change$mean_to
  }


  # ----------------------------------------------------------
  # 9. CORRELATIONS
  # ----------------------------------------------------------

  correlation_pearson <- NULL

  correlation_spearman <- NULL


  if (correlations) {


    correlation_pearson <- tryCatch(

      cor(
        data[time_vars],
        use = "pairwise.complete.obs",
        method = "pearson"
      ),

      error = function(e) {
        NULL
      }
    )


    correlation_spearman <- tryCatch(

      cor(
        data[time_vars],
        use = "pairwise.complete.obs",
        method = "spearman"
      ),

      error = function(e) {
        NULL
      }
    )
  }


  # ----------------------------------------------------------
  # 10. LONG DATA
  # ----------------------------------------------------------

  long_list <- lapply(

    time_vars,

    function(v) {

      data.frame(

        patient =
          data[[id]],

        outcome =
          outcome_name,

        time = v,

        time_label =
          unname(
            time_labels[v]
          ),

        value =
          data[[v]],

        stringsAsFactors = FALSE
      )
    }
  )


  long_data <- do.call(
    rbind,
    long_list
  )


  rownames(long_data) <- NULL


  long_data$time_label <-

    factor(

      long_data$time_label,

      levels =
        unname(
          time_labels[time_vars]
        )
    )


  # ----------------------------------------------------------
  # 11. WITHIN / BETWEEN SUBJECT VARIABILITY
  # ----------------------------------------------------------

  individual_means <- apply(

    data[time_vars],

    1,

    function(x) {

      x <- x[
        is.finite(x)
      ]


      if (length(x) == 0) {

        return(
          NA_real_
        )
      }


      mean(x)
    }
  )


  between_sd <- safe_sd(
    individual_means
  )


  grand_mean <- safe_mean(
    long_data$value
  )


  # ----------------------------------------------------------
  # Patient means
  # ----------------------------------------------------------

  patient_means <- tapply(

    long_data$value,

    long_data$patient,

    safe_mean
  )


  patient_key <- as.character(
    long_data$patient
  )


  long_data$patient_mean <-

    unname(
      patient_means[
        patient_key
      ]
    )


  residuals_within <-

    long_data$value -
    long_data$patient_mean


  within_sd <- safe_sd(
    residuals_within
  )


  # ----------------------------------------------------------
  # 12. ICC
  # ----------------------------------------------------------

  icc <- NA_real_

  ms_between <- NA_real_

  ms_within <- NA_real_

  df_between <- NA_real_

  df_within <- NA_real_


  complete_long <- long_data[

    is.finite(
      long_data$value
    ) &

      !is.na(
        long_data$patient_mean
      ) &

      !is.na(
        long_data$patient
      ),

  ]


  if (
    nrow(complete_long) > 1
  ) {


    patient_counts <- table(
      complete_long$patient
    )


    if (

      length(
        unique(
          patient_counts
        )
      ) == 1 &&

      length(
        patient_counts
      ) > 1

    ) {


      k_actual <-

        as.numeric(
          unique(
            patient_counts
          )
        )


      means2 <- tapply(

        complete_long$value,

        complete_long$patient,

        mean
      )


      grand_mean2 <-
        mean(
          complete_long$value
        )


      ss_between <-

        k_actual *

        sum(

          (
            means2 -
              grand_mean2
          )^2
        )


      ss_within <- sum(

        (
          complete_long$value -

            means2[
              as.character(
                complete_long$patient
              )
            ]

        )^2
      )


      df_between <-
        length(means2) - 1


      df_within <-
        nrow(complete_long) -
        length(means2)


      if (
        df_between > 0 &&
        df_within > 0
      ) {


        ms_between <-
          ss_between /
          df_between


        ms_within <-
          ss_within /
          df_within


        denominator <-

          ms_between +

          (k_actual - 1) *
          ms_within


        if (
          is.finite(denominator) &&
          denominator != 0
        ) {

          icc <-

            (
              ms_between -
                ms_within
            ) /

            denominator
        }
      }
    }
  }


  # ----------------------------------------------------------
  # Variability output
  # ----------------------------------------------------------

  variability <- data.frame(

    outcome =
      outcome_name,

    grand_mean =
      grand_mean,

    between_subject_sd =
      between_sd,

    within_subject_sd =
      within_sd,

    ICC =
      icc,

    ms_between =
      ms_between,

    ms_within =
      ms_within,

    df_between =
      df_between,

    df_within =
      df_within,

    stringsAsFactors = FALSE
  )


  variability[[paste0(
    outcome_name,
    "_grand_mean"
  )]] <-
    variability$grand_mean


  variability[[paste0(
    outcome_name,
    "_between_subject_sd"
  )]] <-
    variability$between_subject_sd


  variability[[paste0(
    outcome_name,
    "_within_subject_sd"
  )]] <-
    variability$within_subject_sd


  # ----------------------------------------------------------
  # 13. MIXED MODEL
  # ----------------------------------------------------------

  mixed_model <- NULL

  model_summary <- NULL

  model_anova <- NULL

  model_error <- NULL


  if (model) {


    if (
      requireNamespace(
        "lme4",
        quietly = TRUE
      )
    ) {


      model_data <- long_data


      model_data$time_factor <-

        factor(

          model_data$time,

          levels =
            time_vars,

          labels =
            unname(
              time_labels[time_vars]
            )
        )


      model_data$patient_factor <-

        factor(
          model_data$patient
        )


      model_formula <-

        value ~
        time_factor +
        (1 | patient_factor)


      # ------------------------------------------------------
      # lmerTest
      # ------------------------------------------------------

      if (
        requireNamespace(
          "lmerTest",
          quietly = TRUE
        )
      ) {


        mixed_model <- tryCatch(

          lmerTest::lmer(

            model_formula,

            data = model_data,

            REML = TRUE,

            na.action = na.omit
          ),

          error = function(e) {

            model_error <<-
              conditionMessage(e)

            NULL
          }
        )


      } else {


        mixed_model <- tryCatch(

          lme4::lmer(

            model_formula,

            data = model_data,

            REML = TRUE,

            na.action = na.omit
          ),

          error = function(e) {

            model_error <<-
              conditionMessage(e)

            NULL
          }
        )
      }


      # ------------------------------------------------------
      # Model summary
      # ------------------------------------------------------

      if (!is.null(mixed_model)) {


        model_summary <-

          summary(
            mixed_model
          )


        model_anova <- tryCatch(

          anova(
            mixed_model
          ),

          error = function(e) {
            NULL
          }
        )
      }


    } else {


      model_error <-

        paste(
          "Il pacchetto 'lme4' non è installato."
        )
    }
  }


  # ----------------------------------------------------------
  # 14. OUTLIERS
  # ----------------------------------------------------------

  outlier_results <- list()

  change_outliers <- list()


  if (outliers) {


    # --------------------------------------------------------
    # OUTLIERS PER TIMEPOINT
    # --------------------------------------------------------

    for (v in time_vars) {


      x <- data[[v]]


      finite_x <- x[
        is.finite(x)
      ]


      if (length(finite_x) < 4) {


        outlier_results[[v]] <-

          data.frame(

            row = integer(0),

            patient =
              data[[id]][
                integer(0)
              ],

            value = numeric(0),

            lower_bound =
              numeric(0),

            upper_bound =
              numeric(0),

            stringsAsFactors =
              FALSE
          )


        next
      }


      q1 <- quantile(

        finite_x,

        0.25,

        names = FALSE
      )


      q3 <- quantile(

        finite_x,

        0.75,

        names = FALSE
      )


      iqr_value <-
        q3 - q1


      lower <-
        q1 -
        1.5 *
        iqr_value


      upper <-
        q3 +
        1.5 *
        iqr_value


      idx <- which(

        !is.na(x) &

          is.finite(x) &

          (
            x < lower |
              x > upper
          )
      )


      if (length(idx) == 0) {


        outlier_results[[v]] <-

          data.frame(

            row = integer(0),

            patient =
              data[[id]][
                integer(0)
              ],

            value = numeric(0),

            lower_bound =
              numeric(0),

            upper_bound =
              numeric(0),

            stringsAsFactors =
              FALSE
          )


      } else {


        outlier_results[[v]] <-

          data.frame(

            row = idx,

            patient =
              data[[id]][idx],

            value =
              x[idx],

            lower_bound =
              rep(
                lower,
                length(idx)
              ),

            upper_bound =
              rep(
                upper,
                length(idx)
              ),

            stringsAsFactors =
              FALSE
          )
      }
    }


    # --------------------------------------------------------
    # OUTLIERS DEL CAMBIAMENTO
    # --------------------------------------------------------

    if (
      nrow(change) > 0
    ) {


      for (i in seq_len(
        nrow(change)
      )) {


        v1 <- change$from[i]

        v2 <- change$to[i]


        x <- data[[v1]]

        y <- data[[v2]]


        keep <- complete.cases(
          x,
          y
        )


        delta <- y[keep] -
          x[keep]


        patients <-
          data[[id]][keep]


        finite <- is.finite(
          delta
        )


        delta <-
          delta[finite]


        patients <-
          patients[finite]


        name <- paste(

          v1,

          v2,

          sep = "_to_"
        )


        if (length(delta) < 4) {


          change_outliers[[name]] <-

            data.frame(

              patient =
                patients[
                  integer(0)
                ],

              change =
                numeric(0),

              lower_bound =
                numeric(0),

              upper_bound =
                numeric(0),

              stringsAsFactors =
                FALSE
            )


          next
        }


        q1 <- quantile(

          delta,

          0.25,

          names = FALSE
        )


        q3 <- quantile(

          delta,

          0.75,

          names = FALSE
        )


        iqr_value <-
          q3 - q1


        lower <-
          q1 -
          1.5 *
          iqr_value


        upper <-
          q3 +
          1.5 *
          iqr_value


        idx <- which(

          is.finite(delta) &

            (
              delta < lower |
                delta > upper
            )
        )


        if (length(idx) == 0) {


          change_outliers[[name]] <-

            data.frame(

              patient =
                patients[
                  integer(0)
                ],

              change =
                numeric(0),

              lower_bound =
                numeric(0),

              upper_bound =
                numeric(0),

              stringsAsFactors =
                FALSE
            )


        } else {


          change_outliers[[name]] <-

            data.frame(

              patient =
                patients[idx],

              change =
                delta[idx],

              lower_bound =
                rep(
                  lower,
                  length(idx)
                ),

              upper_bound =
                rep(
                  upper,
                  length(idx)
                ),

              stringsAsFactors =
                FALSE
            )
        }
      }
    }
  }


  # ----------------------------------------------------------
  # 15. PATIENT TRAJECTORIES
  # ----------------------------------------------------------

  baseline <- data[[
    time_vars[1]
  ]]


  final <- data[[
    time_vars[
      length(time_vars)
    ]
  ]]


  trajectory_summary <- data.frame(

    patient =
      data[[id]],

    outcome =
      outcome_name,

    baseline =
      baseline,

    final =
      final,

    stringsAsFactors =
      FALSE
  )


  trajectory_summary$absolute_change <-

    trajectory_summary$final -

    trajectory_summary$baseline


  trajectory_summary$relative_change_percent <-

    ifelse(

      !is.na(
        trajectory_summary$baseline
      ) &

        trajectory_summary$baseline != 0,

      (

        trajectory_summary$final -

          trajectory_summary$baseline

      ) /

        abs(
          trajectory_summary$baseline
        ) *

        100,

      NA_real_
    )


  # ----------------------------------------------------------
  # Outcome-specific trajectory columns
  # ----------------------------------------------------------

  trajectory_summary[[paste0(
    outcome_name,
    "_baseline"
  )]] <-
    trajectory_summary$baseline


  trajectory_summary[[paste0(
    outcome_name,
    "_final"
  )]] <-
    trajectory_summary$final


  trajectory_summary[[paste0(
    outcome_name,
    "_change"
  )]] <-
    trajectory_summary$absolute_change


  trajectory_summary[[paste0(
    outcome_name,
    "_relative_change_percent"
  )]] <-
    trajectory_summary$relative_change_percent


  # ----------------------------------------------------------
  # 16. PLOTS
  # ----------------------------------------------------------

  plots_list <- list()


  if (plots) {


    if (
      requireNamespace(
        "ggplot2",
        quietly = TRUE
      )
    ) {


      # ======================================================
      # BOXPLOT
      # ======================================================

      boxplot_summary <-

        long_data |>

        dplyr::filter(

          !is.na(value),

          !is.na(time_label)

        ) |>

        dplyr::group_by(

          time,

          time_label

        ) |>

        dplyr::summarise(

          n =
            dplyr::n_distinct(
              patient
            ),

          mean =
            mean(
              value,
              na.rm = TRUE
            ),

          se =
            stats::sd(
              value,
              na.rm = TRUE
            ) /

            sqrt(
              sum(
                !is.na(value)
              )
            ),

          ci_lower =

            mean -

            stats::qt(

              1 - alpha / 2,

              df =
                sum(
                  !is.na(value)
                ) - 1

            ) *

            se,

          ci_upper =

            mean +

            stats::qt(

              1 - alpha / 2,

              df =
                sum(
                  !is.na(value)
                ) - 1

            ) *

            se,

          median =
            stats::median(
              value,
              na.rm = TRUE
            ),

          .groups = "drop"
        )


      # ------------------------------------------------------
      # Label con N
      # ------------------------------------------------------

      boxplot_summary$time_label_n <-

        paste0(

          boxplot_summary$time_label,

          "\n(N = ",

          boxplot_summary$n,

          ")"
        )


      # ------------------------------------------------------
      # Ordine timepoint
      # ------------------------------------------------------

      long_data_boxplot <-

        long_data |>

        dplyr::mutate(

          time_label = factor(

            time_label,

            levels =
              boxplot_summary$time_label
          )
        )


      boxplot_summary$time_label <-

        factor(

          boxplot_summary$time_label,

          levels =
            boxplot_summary$time_label
        )


      # ------------------------------------------------------
      # BOXplot
      # ------------------------------------------------------

      plots_list$boxplot <-

        ggplot2::ggplot(

          long_data_boxplot,

          ggplot2::aes(

            x = time_label,

            y = value
          )
        ) +


        ggplot2::geom_boxplot(

          width = 0.55,

          fill = "grey92",

          color = "grey30",

          linewidth = 0.7,

          outlier.shape = NA,

          na.rm = TRUE
        ) +


        ggplot2::geom_jitter(

          width = 0.10,

          height = 0,

          alpha = 0.22,

          size = 1.5,

          color = "grey35",

          na.rm = TRUE
        ) +


        ggplot2::geom_errorbar(

          data =
            boxplot_summary,

          ggplot2::aes(

            x = time_label,

            ymin = ci_lower,

            ymax = ci_upper
          ),

          inherit.aes = FALSE,

          width = 0.10,

          linewidth = 1.0,

          color = "#2C7BB6",

          na.rm = TRUE
        ) +


        ggplot2::geom_point(

          data =
            boxplot_summary,

          ggplot2::aes(

            x = time_label,

            y = mean
          ),

          inherit.aes = FALSE,

          shape = 21,

          size = 3.8,

          stroke = 1.1,

          fill = "white",

          color = "#2C7BB6",

          na.rm = TRUE
        ) +


        ggplot2::geom_point(

          data =
            boxplot_summary,

          ggplot2::aes(

            x = time_label,

            y = median
          ),

          inherit.aes = FALSE,

          shape = 95,

          size = 7,

          color = "grey20",

          na.rm = TRUE
        ) +


        ggplot2::scale_x_discrete(

          labels =
            boxplot_summary$time_label_n
        ) +


        ggplot2::labs(

          x = NULL,

          y = outcome_display,

          title = paste(

            "Distribution of",

            outcome_display,

            "Across Timepoints"
          ),

          subtitle = paste(

            "Individual observations, interquartile range,",

            "population mean and",

            "95% confidence interval"
          )
        ) +


        ggplot2::theme_minimal(

          base_size = 12
        ) +


        ggplot2::theme(

          plot.title =

            ggplot2::element_text(

              face = "bold",

              size = 15
            ),

          plot.subtitle =

            ggplot2::element_text(

              size = 10,

              color = "grey35"
            ),

          axis.title.y =

            ggplot2::element_text(

              face = "bold"
            ),

          axis.text.x =

            ggplot2::element_text(

              angle = 0,

              hjust = 0.5,

              lineheight = 0.9
            ),

          axis.text.y =

            ggplot2::element_text(

              color = "grey25"
            ),

          panel.grid.minor =

            ggplot2::element_blank(),

          panel.grid.major.x =

            ggplot2::element_blank(),

          panel.grid.major.y =

            ggplot2::element_line(

              color = "grey90",

              linewidth = 0.4
            ),

          legend.position = "none",

          plot.margin =

            ggplot2::margin(

              12,

              18,

              12,

              12
            )
        )


      # ======================================================
      # SPAGHETTI PLOT
      # ======================================================

      stable_threshold <- 0


      long_data_spaghetti <- long_data


      long_data_spaghetti$time_index <-

        match(

          long_data_spaghetti$time,

          time_vars
        )


      trajectory_summary_plot <- data.frame(

        time = time_vars,

        time_index =
          seq_along(time_vars),

        time_label =
          unname(
            time_labels[time_vars]
          ),

        mean =
          descriptives$mean,

        se =
          descriptives$se,

        ci_lower =
          descriptives$ci_lower,

        ci_upper =
          descriptives$ci_upper,

        stringsAsFactors = FALSE
      )


      n_per_timepoint <-

        long_data_spaghetti |>

        dplyr::filter(

          !is.na(value),

          !is.na(time_index)
        ) |>

        dplyr::group_by(

          time_index
        ) |>

        dplyr::summarise(

          n =
            dplyr::n_distinct(
              patient
            ),

          .groups = "drop"
        )


      trajectory_summary_plot <-

        trajectory_summary_plot |>

        dplyr::left_join(

          n_per_timepoint,

          by = "time_index"
        )


      trajectory_summary_plot$time_label_n <-

        paste0(

          trajectory_summary_plot$time_label,

          "\n(N = ",

          trajectory_summary_plot$n,

          ")"
        )


      # ------------------------------------------------------
      # Direction
      # ------------------------------------------------------

      trajectory_segments <-

        long_data_spaghetti |>

        dplyr::filter(

          !is.na(value),

          !is.na(time_index)
        ) |>

        dplyr::arrange(

          patient,

          time_index
        ) |>

        dplyr::group_by(

          patient
        ) |>

        dplyr::mutate(

          next_time_index =
            dplyr::lead(
              time_index
            ),

          next_value =
            dplyr::lead(
              value
            ),

          change =
            next_value -
            value,

          trajectory_direction =

            dplyr::case_when(

              is.na(next_value) ~
                NA_character_,

              change >
                stable_threshold ~
                "Increase",

              change <
                -stable_threshold ~
                "Decrease",

              TRUE ~
                "Stable"
            )
        ) |>

        dplyr::ungroup()


      # ------------------------------------------------------
      # Spaghetti
      # ------------------------------------------------------

      plots_list$spaghetti <-

        ggplot2::ggplot() +


        ggplot2::geom_ribbon(

          data =
            trajectory_summary_plot,

          ggplot2::aes(

            x = time_index,

            ymin = ci_lower,

            ymax = ci_upper,

            group = 1
          ),

          inherit.aes = FALSE,

          alpha = 0.20,

          na.rm = TRUE
        ) +


        ggplot2::geom_segment(

          data =
            trajectory_segments,

          ggplot2::aes(

            x = time_index,

            xend =
              next_time_index,

            y = value,

            yend =
              next_value,

            group = patient,

            color =
              trajectory_direction
          ),

          alpha = 0.20,

          linewidth = 0.6,

          na.rm = TRUE
        ) +


        ggplot2::geom_point(

          data =
            long_data_spaghetti,

          ggplot2::aes(

            x = time_index,

            y = value
          ),

          color = "grey40",

          alpha = 0.20,

          size = 1.5,

          na.rm = TRUE
        ) +


        ggplot2::geom_line(

          data =
            trajectory_summary_plot,

          ggplot2::aes(

            x = time_index,

            y = mean,

            group = 1
          ),

          inherit.aes = FALSE,

          linewidth = 1.5,

          color = "black",

          na.rm = TRUE
        ) +


        ggplot2::geom_point(

          data =
            trajectory_summary_plot,

          ggplot2::aes(

            x = time_index,

            y = mean
          ),

          inherit.aes = FALSE,

          shape = 21,

          size = 4,

          stroke = 1.2,

          fill = "white",

          color = "black",

          na.rm = TRUE
        ) +


        ggplot2::scale_color_manual(

          values = c(

            "Increase" =
              "#2C7BB6",

            "Decrease" =
              "#D7191C",

            "Stable" =
              "grey60"
          ),

          breaks = c(

            "Increase",

            "Decrease",

            "Stable"
          ),

          labels = c(

            "↑ Increase",

            "↓ Decrease",

            "→ Stable"
          ),

          name =
            "Change between consecutive timepoints"
        ) +


        ggplot2::scale_x_continuous(

          breaks =
            trajectory_summary_plot$time_index,

          labels =
            trajectory_summary_plot$time_label_n,

          expand =
            ggplot2::expansion(

              mult =
                c(
                  0.03,
                  0.03
                )
            )
        ) +


        ggplot2::labs(

          x = NULL,

          y = outcome_display,

          title = paste(

            "Individual",

            outcome_display,

            "Longitudinal Trajectories"
          ),

          subtitle = paste(

            "Subject-level changes between consecutive timepoints,",

            "with population mean and 95% confidence interval"
          )
        ) +


        ggplot2::theme_minimal(

          base_size = 12
        ) +


        ggplot2::theme(

          plot.title =

            ggplot2::element_text(

              face = "bold",

              size = 15
            ),

          plot.subtitle =

            ggplot2::element_text(

              size = 10
            ),

          axis.title =

            ggplot2::element_text(

              face = "bold"
            ),

          axis.text.x =

            ggplot2::element_text(

              angle = 0,

              hjust = 0.5,

              lineheight = 0.9
            ),

          panel.grid.minor =

            ggplot2::element_blank(),

          panel.grid.major.x =

            ggplot2::element_blank(),

          legend.position =
            "bottom",

          legend.title =

            ggplot2::element_text(

              face = "bold"
            ),

          legend.text =

            ggplot2::element_text(

              size = 10
            ),

          plot.margin =

            ggplot2::margin(

              12,

              18,

              12,

              12
            )
        )


      # ======================================================
      # MEAN + CI
      # ======================================================

      plots_list$mean_ci <-

        ggplot2::ggplot(

          descriptives,

          ggplot2::aes(

            x = label,

            y = mean,

            group = 1
          )
        ) +


        ggplot2::geom_line(

          color = "#2C7FB8",

          linewidth = 0.8,

          na.rm = TRUE
        ) +


        ggplot2::geom_point(

          size = 3,

          color = "#2C7FB8",

          na.rm = TRUE
        ) +


        ggplot2::geom_errorbar(

          ggplot2::aes(

            ymin = ci_lower,

            ymax = ci_upper
          ),

          width = 0.1,

          na.rm = TRUE
        ) +


        ggplot2::labs(

          x = "Time",

          y = paste(
            "Mean",
            outcome_display
          ),

          title = paste(

            "Mean",

            outcome_display,

            "and 95% CI"
          )
        ) +


        ggplot2::theme_minimal()


      # ======================================================
      # INDIVIDUAL CHANGE
      # ======================================================

      plots_list$change <-

        ggplot2::ggplot(

          trajectory_summary,

          ggplot2::aes(

            x = "All patients",

            y = absolute_change
          )
        ) +


        ggplot2::geom_boxplot(

          fill = "#59A14F",

          alpha = 0.7,

          na.rm = TRUE
        ) +


        ggplot2::geom_jitter(

          width = 0.08,

          alpha = 0.30,

          na.rm = TRUE
        ) +


        ggplot2::geom_hline(

          yintercept = 0,

          linetype = "dashed"
        ) +


        ggplot2::labs(

          x = NULL,

          y = paste(

            "Final - Baseline",

            outcome_display
          ),

          title = paste(

            "Individual",

            outcome_display,

            "Change"
          )
        ) +


        ggplot2::theme_minimal()
    }
  }


  # ----------------------------------------------------------
  # 17. GLOBAL TEST OF TIME
  # ----------------------------------------------------------

  global_time_test <-
    model_anova


  # ----------------------------------------------------------
  # 18. OVERVIEW
  # ----------------------------------------------------------

  overview <- list(

    outcome =
      outcome_name,

    outcome_display =
      outcome_display,

    n_rows =
      n_rows,

    n_patients =
      n_patients,

    n_timepoints =
      length(time_vars),

    timepoints =
      time_vars,

    time_labels =
      time_labels,

    duplicated_ids =
      duplicated_ids,

    complete_profiles =
      complete_cases,

    complete_profiles_pct =

      complete_cases /
      n_rows *
      100
  )


  # ----------------------------------------------------------
  # 19. FINAL RESULT
  # ----------------------------------------------------------

  result <- list(

    call =
      match.call(),


    overview =
      overview,


    outcome =
      outcome_name,


    outcome_display =
      outcome_display,


    time_vars =
      time_vars,


    time_labels =
      time_labels,


    descriptives =
      descriptives,


    missing =
      list(

        by_time =
          missing_summary,

        by_patient =
          missing_by_patient
      ),


    change =
      change,


    correlations =
      list(

        pearson =
          correlation_pearson,

        spearman =
          correlation_spearman
      ),


    variability =
      variability,


    trajectories =
      trajectory_summary,


    model =
      list(

        fitted_model =
          mixed_model,

        summary =
          model_summary,

        anova =
          model_anova,

        global_time_test =
          global_time_test,

        error =
          model_error
      ),


    outliers =
      list(

        by_time =

          if (outliers)

            outlier_results

        else

          NULL,


        change =

          if (outliers)

            change_outliers

        else

          NULL
      ),


    plots =
      plots_list,


    long_data =
      long_data
  )


  class(result) <-
    "mira_info"


  # ----------------------------------------------------------
  # 20. PRINT
  # ----------------------------------------------------------

  if (verbose) {

    print(
      result
    )
  }


  invisible(
    result
  )
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
                            ...) {


  # ----------------------------------------------------------
  # Helper
  # ----------------------------------------------------------

  line <- function(
    char = "-",
    n = 72
  ) {

    cat(
      paste0(
        strrep(
          char,
          n
        ),
        "\n"
      )
    )
  }


  fmt_num <- function(
    z,
    digits = 3
  ) {

    ifelse(

      is.na(z),

      "NA",

      formatC(

        z,

        format = "f",

        digits = digits
      )
    )
  }


  fmt_pct <- function(
    z,
    digits = 1
  ) {

    ifelse(

      is.na(z),

      "NA",

      paste0(

        formatC(

          z,

          format = "f",

          digits = digits
        ),

        "%"
      )
    )
  }


  fmt_p <- function(p) {

    ifelse(

      is.na(p),

      "NA",

      ifelse(

        p < 0.001,

        "<0.001",

        formatC(

          p,

          format = "f",

          digits = 3
        )
      )
    )
  }


  effect_label <- function(d) {

    if (is.na(d)) {

      return(
        NA_character_
      )
    }


    a <- abs(d)


    if (a < 0.20) {

      "negligible"

    } else if (a < 0.50) {

      "small"

    } else if (a < 0.80) {

      "moderate"

    } else {

      "large"
    }
  }


  # ----------------------------------------------------------
  # Extract automatic outcome
  # ----------------------------------------------------------

  outcome_display <-

    if (!is.null(x$overview$outcome_display)) {

      x$overview$outcome_display

    } else if (!is.null(x$outcome_display)) {

      x$outcome_display

    } else {

      "OUTCOME"
    }


  # ----------------------------------------------------------
  # HEADER
  # ----------------------------------------------------------

  cat("\n")

  line("=")


  cat(

    sprintf(

      "                         MIRA INFO — %s\n",

      outcome_display
    )
  )


  cat(

    "        Comprehensive Longitudinal Dataset Snapshot\n"
  )


  line("=")


  # ==========================================================
  # DATASET OVERVIEW
  # ==========================================================

  cat(

    sprintf(

      "\nDATASET OVERVIEW — %s\n",

      outcome_display
    )
  )


  line()


  ov <- x$overview

  d <- x$descriptives


  total_cells <-

    ov$n_rows *
    ov$n_timepoints


  observed_cells <-

    sum(
      d$n
    )


  missing_cells <-

    sum(
      d$missing
    )


  non_finite_cells <-

    sum(
      d$non_finite
    )


  complete_profiles <-

    ov$complete_profiles


  incomplete_profiles <-

    ov$n_rows -
    complete_profiles


  # ----------------------------------------------------------
  # Outcome
  # ----------------------------------------------------------

  cat(

    sprintf(

      "Outcome:                     %s\n",

      outcome_display
    )
  )


  cat(

    sprintf(

      "ID variable:                 %s\n",

      as.character(
        x$call$id
      )
    )
  )


  cat(

    sprintf(

      "Rows:                        %s\n",

      ov$n_rows
    )
  )


  cat(

    sprintf(

      "Unique subjects:             %s\n",

      ov$n_patients
    )
  )


  cat(

    sprintf(

      "Timepoints:                  %s\n",

      ov$n_timepoints
    )
  )


  cat(

    sprintf(

      "Time variables:              %s\n",

      paste(

        ov$timepoints,

        collapse = ", "
      )
    )
  )


  cat(

    sprintf(

      "Time labels:                 %s\n",

      paste(

        unname(
          ov$time_labels
        ),

        collapse = ", "
      )
    )
  )


  cat(

    sprintf(

      "Duplicated IDs:              %s\n",

      ov$duplicated_ids
    )
  )


  # ==========================================================
  # DATA QUALITY
  # ==========================================================

  cat(

    sprintf(

      "\nDATA QUALITY — %s\n",

      outcome_display
    )
  )


  line()


  cat(

    sprintf(

      "Total longitudinal cells:    %s\n",

      total_cells
    )
  )


  cat(

    sprintf(

      "Finite observations:         %s (%s)\n",

      observed_cells,

      fmt_pct(

        observed_cells /
          total_cells *
          100
      )
    )
  )


  cat(

    sprintf(

      "Missing values:              %s (%s)\n",

      missing_cells,

      fmt_pct(

        missing_cells /
          total_cells *
          100
      )
    )
  )


  cat(

    sprintf(

      "Non-finite values:           %s (%s)\n",

      non_finite_cells,

      fmt_pct(

        non_finite_cells /
          total_cells *
          100
      )
    )
  )


  cat(

    sprintf(

      "Complete profiles:           %s (%s)\n",

      complete_profiles,

      fmt_pct(
        ov$complete_profiles_pct
      )
    )
  )


  cat(

    sprintf(

      "Incomplete profiles:         %s (%s)\n",

      incomplete_profiles,

      fmt_pct(

        incomplete_profiles /
          ov$n_rows *
          100
      )
    )
  )


  # ==========================================================
  # DESCRIPTIVE STATISTICS
  # ==========================================================

  cat(

    sprintf(

      "\nDESCRIPTIVE STATISTICS FOR %s BY TIMEPOINT\n",

      outcome_display
    )
  )


  line()


  desc_print <- data.frame(

    Time =
      d$label,

    N =
      d$n,

    Missing =

      paste0(

        d$missing,

        " (",

        formatC(

          d$missing_pct,

          format = "f",

          digits = 1
        ),

        "%)"
      ),

    Mean =
      round(
        d$mean,
        digits
      ),

    SD =
      round(
        d$sd,
        digits
      ),

    Median =
      round(
        d$median,
        digits
      ),

    IQR =
      round(
        d$iqr,
        digits
      ),

    Min =
      round(
        d$min,
        digits
      ),

    Max =
      round(
        d$max,
        digits
      ),

    CV_pct =
      round(
        d$cv_percent,
        1
      ),

    CI_lower =
      round(
        d$ci_lower,
        digits
      ),

    CI_upper =
      round(
        d$ci_upper,
        digits
      ),

    check.names = FALSE
  )


  print(

    desc_print,

    row.names = FALSE
  )


  # ==========================================================
  # MISSING DATA
  # ==========================================================

  cat(

    sprintf(

      "\nMISSING DATA FOR %s BY TIMEPOINT\n",

      outcome_display
    )
  )


  line()


  miss <- x$missing$by_time


  miss_print <- data.frame(

    Time =
      miss$label,

    Missing_n =
      miss$missing_n,

    Missing_pct =
      round(
        miss$missing_pct,
        1
      )
  )


  print(

    miss_print,

    row.names = FALSE
  )


  # ==========================================================
  # LONGITUDINAL CHANGE
  # ==========================================================

  cat(

    sprintf(

      "\nLONGITUDINAL CHANGE — %s\n",

      outcome_display
    )
  )


  line()


  if (nrow(x$change) > 0) {


    ch <- x$change


    # --------------------------------------------------------
    # Baseline -> Final
    # --------------------------------------------------------

    baseline_final <-

      ch[

        ch$from ==
          ov$timepoints[1] &

          ch$to ==
          ov$timepoints[
            length(
              ov$timepoints
            )
          ],

        ,

        drop = FALSE
      ]


    if (nrow(baseline_final) == 1) {


      bf <-
        baseline_final[1, ]


      cat(

        sprintf(

          "Baseline -> Final:          %s -> %s\n",

          bf$from_label,

          bf$to_label
        )
      )


      cat(

        sprintf(

          "Paired observations:        %s\n",

          bf$n
        )
      )


      cat(

        sprintf(

          "Mean %s change:             %s\n",

          outcome_display,

          fmt_num(

            bf$mean_change,

            digits
          )
        )
      )


      cat(

        sprintf(

          "%s%% CI:                     [%s, %s]\n",

          100 * (1 - 0.05),

          fmt_num(

            bf$ci_lower,

            digits
          ),

          fmt_num(

            bf$ci_upper,

            digits
          )
        )
      )


      cat(

        sprintf(

          "Cohen's dz:                 %s (%s)\n",

          fmt_num(

            bf$cohens_dz,

            digits
          ),

          effect_label(
            bf$cohens_dz
          )
        )
      )


      cat(

        sprintf(

          "Improved / worsened / same: %s / %s / %s\n",

          fmt_pct(
            bf$improved_pct
          ),

          fmt_pct(
            bf$worsened_pct
          ),

          fmt_pct(
            bf$unchanged_pct
          )
        )
      )


      cat(

        sprintf(

          "Paired t-test p:            %s\n",

          fmt_p(
            bf$paired_t_p
          )
        )
      )


      cat(

        sprintf(

          "Wilcoxon p:                 %s\n",

          fmt_p(
            bf$wilcoxon_p
          )
        )
      )
    }


    # --------------------------------------------------------
    # Pairwise comparisons
    # --------------------------------------------------------

    cat(

      sprintf(

        "\nPAIRWISE %s TIMEPOINT COMPARISONS\n",

        outcome_display
      )
    )


    change_print <- data.frame(

      From =
        ch$from_label,

      To =
        ch$to_label,

      N =
        ch$n,

      Mean_change =
        round(
          ch$mean_change,
          digits
        ),

      CI_lower =
        round(
          ch$ci_lower,
          digits
        ),

      CI_upper =
        round(
          ch$ci_upper,
          digits
        ),

      Cohen_dz =
        round(
          ch$cohens_dz,
          digits
        ),

      Improved_pct =
        round(
          ch$improved_pct,
          1
        ),

      Worsened_pct =
        round(
          ch$worsened_pct,
          1
        ),

      t_p =

        vapply(

          ch$paired_t_p,

          fmt_p,

          character(1)
        ),

      Wilcoxon_p =

        vapply(

          ch$wilcoxon_p,

          fmt_p,

          character(1)
        ),

      check.names = FALSE
    )


    if (
      nrow(change_print) >
      max_rows
    ) {


      cat(

        sprintf(

          "\nShowing first %s of %s comparisons.\n",

          max_rows,

          nrow(change_print)
        )
      )


      change_print <-

        head(

          change_print,

          max_rows
        )
    }


    print(

      change_print,

      row.names = FALSE
    )


  } else {


    cat(

      "No longitudinal comparisons available.\n"
    )
  }


  # ==========================================================
  # CORRELATIONS
  # ==========================================================

  if (

    correlations &&

    !is.null(
      x$correlations$pearson
    )

  ) {


    cat(

      sprintf(

        "\nCORRELATIONS — %s\n",

        outcome_display
      )
    )


    line()


    cor_mat <-

      x$correlations$pearson


    if (ncol(cor_mat) >= 2) {


      upper_values <-

        cor_mat[
          upper.tri(cor_mat)
        ]


      finite_cor <-

        upper_values[
          is.finite(
            upper_values
          )
        ]


      if (length(finite_cor) > 0) {


        cat(

          sprintf(

            "Pearson correlation range: %s to %s\n",

            fmt_num(
              min(
                finite_cor
              ),
              digits
            ),

            fmt_num(
              max(
                finite_cor
              ),
              digits
            )
          )
        )


        cat(

          sprintf(

            "Median correlation:        %s\n",

            fmt_num(

              median(
                finite_cor
              ),

              digits
            )
          )
        )
      }
    }


    cat(

      sprintf(

        "\nPearson correlation matrix for %s:\n",

        outcome_display
      )
    )


    print(

      round(
        cor_mat,
        digits
      )
    )
  }


  # ==========================================================
  # VARIABILITY
  # ==========================================================

  cat(

    sprintf(

      "\nVARIABILITY AND SUBJECT DEPENDENCE — %s\n",

      outcome_display
    )
  )


  line()


  v <- x$variability


  ratio <-

    if (

      !is.na(
        v$between_subject_sd
      ) &&

      v$between_subject_sd != 0

    ) {


      v$within_subject_sd /
        v$between_subject_sd


    } else {


      NA_real_
    }


  cat(

    sprintf(

      "Grand mean %s:              %s\n",

      outcome_display,

      fmt_num(

        v$grand_mean,

        digits
      )
    )
  )


  cat(

    sprintf(

      "Between-subject SD:          %s\n",

      fmt_num(

        v$between_subject_sd,

        digits
      )
    )
  )


  cat(

    sprintf(

      "Within-subject SD:           %s\n",

      fmt_num(

        v$within_subject_sd,

        digits
      )
    )
  )


  cat(

    sprintf(

      "Within / Between ratio:      %s\n",

      fmt_num(

        ratio,

        digits
      )
    )
  )


  cat(

    sprintf(

      "ICC:                         %s\n",

      fmt_num(

        v$ICC,

        digits
      )
    )
  )


  # ==========================================================
  # MIXED MODEL
  # ==========================================================

  if (model) {


    cat(

      sprintf(

        "\nMIXED-EFFECTS MODEL — %s\n",

        outcome_display
      )
    )


    line()


    if (
      !is.null(
        x$model$fitted_model
      )
    ) {


      cat(

        sprintf(

          "Model: %s ~ time + (1 | subject)\n",

          outcome_display
        )
      )


      if (
        !is.null(
          x$model$anova
        )
      ) {


        cat(

          sprintf(

            "\nGlobal test of %s over time:\n",

            outcome_display
          )
        )


        print(
          x$model$anova
        )
      }


      cat(

        "\nFull model output available in:\n"
      )


      cat(

        "  result$model$summary\n"
      )


    } else if (
      !is.null(
        x$model$error
      )
    ) {


      cat(

        "Model not estimated: ",

        x$model$error,

        "\n",

        sep = ""
      )
    }
  }


  # ==========================================================
  # OUTLIERS
  # ==========================================================

  if (

    outliers &&

    !is.null(
      x$outliers$by_time
    )

  ) {


    cat(

      sprintf(

        "\nOUTLIER SUMMARY — %s\n",

        outcome_display
      )
    )


    line()


    outlier_counts <-

      vapply(

        x$outliers$by_time,

        nrow,

        integer(1)
      )


    total_outliers <-

      sum(
        outlier_counts
      )


    cat(

      sprintf(

        "Total %s timepoint outliers: %s\n",

        outcome_display,

        total_outliers
      )
    )


    for (
      nm in names(
        outlier_counts
      )
    ) {


      cat(

        sprintf(

          "  %-25s %s\n",

          nm,

          outlier_counts[nm]
        )
      )
    }


    if (
      !is.null(
        x$outliers$change
      )
    ) {


      change_outlier_counts <-

        vapply(

          x$outliers$change,

          nrow,

          integer(1)
        )


      cat(

        sprintf(

          "Total %s change outliers:    %s\n",

          outcome_display,

          sum(
            change_outlier_counts
          )
        )
      )
    }
  }


  # ==========================================================
  # TRAJECTORIES
  # ==========================================================

  if (
    !is.null(
      x$trajectories
    )
  ) {


    tr <-
      x$trajectories


    finite_change <-

      tr$absolute_change[

        is.finite(
          tr$absolute_change
        )
      ]


    if (
      length(finite_change) > 0
    ) {


      cat(

        sprintf(

          "\nINDIVIDUAL %s TRAJECTORIES\n",

          outcome_display
        )
      )


      line()


      cat(

        sprintf(

          "Subjects with baseline-final data: %s\n",

          length(
            finite_change
          )
        )
      )


      cat(

        sprintf(

          "Mean individual %s change:           %s\n",

          outcome_display,

          fmt_num(

            mean(
              finite_change
            ),

            digits
          )
        )
      )


      cat(

        sprintf(

          "Median individual %s change:         %s\n",

          outcome_display,

          fmt_num(

            median(
              finite_change
            ),

            digits
          )
        )
      )


      cat(

        sprintf(

          "Range of individual %s change:       [%s, %s]\n",

          outcome_display,

          fmt_num(

            min(
              finite_change
            ),

            digits
          ),

          fmt_num(

            max(
              finite_change
            ),

            digits
          )
        )
      )
    }
  }


  # ==========================================================
  # OUTPUT GUIDE
  # ==========================================================

  cat("\n")

  line("=")


  cat(

    sprintf(

      "COMPLETE OUTPUT — %s\n",

      outcome_display
    )
  )


  line()


  cat(

    "Use names(result) to inspect all available components.\n\n"
  )


  cat(

    "Key components:\n"
  )


  cat(

    sprintf(

      "  outcome        Automatically detected outcome (%s)\n",

      outcome_display
    )
  )


  cat(

    "  overview       Dataset structure and completeness\n"
  )


  cat(

    sprintf(

      "  descriptives   Detailed %s statistics by timepoint\n",

      outcome_display
    )
  )


  cat(

    "  missing        Missingness by timepoint and subject\n"
  )


  cat(

    sprintf(

      "  change         All %s longitudinal pairwise comparisons\n",

      outcome_display
    )
  )


  cat(

    sprintf(

      "  correlations   %s Pearson and Spearman correlation matrices\n",

      outcome_display
    )
  )


  cat(

    sprintf(

      "  variability    %s within/between-subject variability and ICC\n",

      outcome_display
    )
  )


  cat(

    sprintf(

      "  trajectories    %s subject-level baseline-to-final changes\n",

      outcome_display
    )
  )


  cat(

    sprintf(

      "  model          %s mixed-effects model and global time test\n",

      outcome_display
    )
  )


  cat(

    sprintf(

      "  outliers       Detected %s value and change outliers\n",

      outcome_display
    )
  )


  cat(

    "  plots          Generated visualization objects\n"
  )


  cat(

    sprintf(

      "  long_data      Long-format %s dataset used internally\n",

      outcome_display
    )
  )


  line("=")


  cat("\n")


  invisible(x)
}
