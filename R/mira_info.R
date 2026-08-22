# ============================================================
# MIRA_INFO
# Comprehensive longitudinal data analysis
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
    stop("data deve essere un data.frame.")
  }

  if (nrow(data) == 0) {
    stop("Il dataset non contiene osservazioni.")
  }

  if (!is.character(id) || length(id) != 1) {
    stop("id deve essere il nome di una sola variabile.")
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

    stop("alpha deve essere un numero compreso tra 0 e 1.")
  }


  # ----------------------------------------------------------
  # 1. IDENTIFICAZIONE TIMEPOINT
  # ----------------------------------------------------------

  if (is.null(time_vars)) {

    # Prima cerca t0, t1, t2, ...
    time_vars <- grep(
      "^t[0-9]+$",
      names(data),
      value = TRUE
    )

    # Se non trova t0/t1/t2...
    # cerca variabili numeriche
    if (length(time_vars) == 0) {

      numeric_vars <- names(data)[
        vapply(
          data,
          is.numeric,
          logical(1)
        )
      ]

      time_vars <- setdiff(
        numeric_vars,
        id
      )
    }
  }

  if (!is.character(time_vars) ||
      length(time_vars) < 2) {

    stop(
      paste0(
        "Sono necessarie almeno due variabili longitudinali.\n",
        "Variabili trovate: ",
        paste(time_vars, collapse = ", ")
      )
    )
  }


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
        paste(missing_vars, collapse = ", ")
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
        paste(non_numeric, collapse = ", ")
      )
    )
  }


  # ----------------------------------------------------------
  # 3. TIME LABELS
  # ----------------------------------------------------------

  if (is.null(time_labels)) {

    time_labels <- time_vars

  } else {

    if (length(time_labels) != length(time_vars)) {

      stop(
        "time_labels deve avere la stessa lunghezza di time_vars."
      )
    }

    if (!is.character(time_labels)) {
      time_labels <- as.character(time_labels)
    }
  }

  names(time_labels) <- time_vars


  # ----------------------------------------------------------
  # 4. HELPER FUNCTIONS
  # ----------------------------------------------------------

  mean_ci <- function(x, alpha = 0.05) {

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

    x <- x[is.finite(x)]

    if (length(x) == 0) {
      return(NA_real_)
    }

    mean(x)
  }


  safe_sd <- function(x) {

    x <- x[is.finite(x)]

    if (length(x) < 2) {
      return(NA_real_)
    }

    sd(x)
  }


  safe_median <- function(x) {

    x <- x[is.finite(x)]

    if (length(x) == 0) {
      return(NA_real_)
    }

    median(x)
  }


  safe_quantile <- function(x, p) {

    x <- x[is.finite(x)]

    if (length(x) == 0) {
      return(NA_real_)
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

    x <- x[is.finite(x)]

    if (length(x) < 2) {
      return(NA_real_)
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


  # Complete profile:
  # il paziente deve avere tutte le misure disponibili

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

      n <- length(finite_x)

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
          missing_n / n_rows * 100,

        non_finite = non_finite_n,

        mean = m,

        sd = s,

        variance = variance,

        se = unname(
          ci["se"]
        ),

        ci_lower = unname(
          ci["lower"]
        ),

        ci_upper = unname(
          ci["upper"]
        ),

        median = safe_median(x),

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

  rownames(descriptives) <- NULL


  # ----------------------------------------------------------
  # 7. MISSINGNESS
  # ----------------------------------------------------------

  missing_by_patient <- data.frame(

    patient = data[[id]],

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

    label = unname(
      time_labels[time_vars]
    ),

    missing_n = vapply(
      data[time_vars],
      function(x) {
        sum(is.na(x))
      },
      numeric(1)
    ),

    missing_pct = vapply(
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

      keep <- complete.cases(x, y)

      x_complete <- x[keep]

      y_complete <- y[keep]

      delta <- y_complete - x_complete

      # Elimina eventuali Inf
      finite <- is.finite(delta)

      delta <- delta[finite]

      x_complete <- x_complete[finite]

      y_complete <- y_complete[finite]

      n <- length(delta)

      if (n == 0) {
        next
      }

      delta_mean <- mean(delta)

      delta_sd <-
        if (n >= 2) {
          sd(delta)
        } else {
          NA_real_
        }

      delta_se <-
        if (n >= 2) {
          delta_sd / sqrt(n)
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
          tcrit * delta_se

        ci_high <-
          delta_mean +
          tcrit * delta_se

      } else {

        ci_low <- NA_real_

        ci_high <- NA_real_
      }


      # Cohen's dz

      effect_size <-
        if (
          n >= 2 &&
          !is.na(delta_sd) &&
          delta_sd > 0
        ) {

          delta_mean / delta_sd

        } else {

          NA_real_
        }


      # Paired t-test

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


      # Wilcoxon

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


      improved <- sum(
        delta > 0
      )

      worsened <- sum(
        delta < 0
      )

      unchanged <- sum(
        delta == 0
      )


      change_list[[counter]] <-
        data.frame(

          from = v1,

          to = v2,

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
            improved / n * 100,

          worsened_n =
            worsened,

          worsened_pct =
            worsened / n * 100,

          unchanged_n =
            unchanged,

          unchanged_pct =
            unchanged / n * 100,

          paired_t_p =
            p_value,

          wilcoxon_p =
            wilcox_p,

          stringsAsFactors = FALSE
        )

      counter <- counter + 1
    }
  }


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

        patient = data[[id]],

        time = v,

        value = data[[v]],

        stringsAsFactors = FALSE
      )
    }
  )

  long_data <- do.call(
    rbind,
    long_list
  )

  rownames(long_data) <- NULL


  long_data$time_label <- factor(
    long_data$time,
    levels = time_vars,
    labels = time_labels
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
        return(NA_real_)
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


  # Patient means

  patient_means <- tapply(
    long_data$value,
    long_data$patient,
    safe_mean
  )


  # Match patient means back to long data

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
    is.finite(long_data$value) &
      !is.na(long_data$patient_mean) &
      !is.na(long_data$patient),
  ]


  if (
    nrow(complete_long) > 1
  ) {

    patient_counts <- table(
      complete_long$patient
    )

    # ICC classico ANOVA:
    # richiede dati bilanciati

    if (
      length(unique(patient_counts)) == 1 &&
      length(patient_counts) > 1
    ) {

      k_actual <-
        as.numeric(
          unique(patient_counts)
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


  variability <- data.frame(

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
          levels = time_vars,
          labels = time_labels
        )

      model_data$patient_factor <-
        factor(
          model_data$patient
        )


      model_formula <-
        value ~ time_factor +
        (1 | patient_factor)


      # Se lmerTest è disponibile,
      # utilizziamo direttamente lmerTest

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

            lower_bound = numeric(0),

            upper_bound = numeric(0),

            stringsAsFactors = FALSE
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
        1.5 * iqr_value

      upper <-
        q3 +
        1.5 * iqr_value


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

            lower_bound = numeric(0),

            upper_bound = numeric(0),

            stringsAsFactors = FALSE
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

            stringsAsFactors = FALSE
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
          1.5 * iqr_value

        upper <-
          q3 +
          1.5 * iqr_value


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

  baseline <- data[[time_vars[1]]]

  final <- data[[time_vars[length(time_vars)]]]


  trajectory_summary <- data.frame(

    patient =
      data[[id]],

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

      # ------------------------------------------------------
      # BOXPLOT
      # ------------------------------------------------------

      plots_list$boxplot <-

        ggplot2::ggplot(

          long_data,

          ggplot2::aes(
            x = time_label,
            y = value
          )
        ) +

        ggplot2::geom_boxplot(
          fill = "#4C78A8",
          alpha = 0.7,
          na.rm = TRUE
        ) +

        ggplot2::geom_jitter(
          width = 0.08,
          alpha = 0.25,
          na.rm = TRUE
        ) +

        ggplot2::labs(

          x = "Time",

          y = "Value",

          title =
            "Distribution by timepoint"
        ) +

        ggplot2::theme_minimal()


      # ------------------------------------------------------
      # SPAGHETTI
      # ------------------------------------------------------

      plots_list$spaghetti <-

        ggplot2::ggplot(

          long_data,

          ggplot2::aes(

            x = time_label,

            y = value,

            group = patient
          )
        ) +

        ggplot2::geom_line(

          alpha = 0.15,

          na.rm = TRUE
        ) +

        ggplot2::geom_point(

          alpha = 0.20,

          na.rm = TRUE
        ) +

        ggplot2::stat_summary(

          ggplot2::aes(
            group = 1
          ),

          fun = mean,

          geom = "line",

          color = "red",

          linewidth = 1.5,

          na.rm = TRUE
        ) +

        ggplot2::stat_summary(

          ggplot2::aes(
            group = 1
          ),

          fun = mean,

          geom = "point",

          color = "red",

          size = 3,

          na.rm = TRUE
        ) +

        ggplot2::labs(

          x = "Time",

          y = "Value",

          title =
            "Individual trajectories"
        ) +

        ggplot2::theme_minimal()


      # ------------------------------------------------------
      # MEAN + CI
      # ------------------------------------------------------

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

          y = "Mean",

          title =
            "Mean and 95% CI"
        ) +

        ggplot2::theme_minimal()


      # ------------------------------------------------------
      # INDIVIDUAL CHANGE
      # ------------------------------------------------------

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

          y = "Final - Baseline",

          title =
            "Individual change"
        ) +

        ggplot2::theme_minimal()
    }
  }


  # ----------------------------------------------------------
  # 17. GLOBAL TEST OF TIME
  # ----------------------------------------------------------

  global_time_test <- model_anova


  # ----------------------------------------------------------
  # 18. OVERVIEW
  # ----------------------------------------------------------

  overview <- list(

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
    print(result)
  }


  invisible(result)
}



# ============================================================
# PRINT METHOD
# ============================================================

print.mira_info <- function(x, ...) {

  cat("\n")

  cat(
    "============================================================\n"
  )

  cat(
    "                    MIRA INFO REPORT\n"
  )

  cat(
    "============================================================\n"
  )


  # ----------------------------------------------------------
  # DATASET
  # ----------------------------------------------------------

  cat("\nDATASET\n")

  cat(
    "------------------------------------------------------------\n"
  )

  cat(
    "Patients:       ",
    x$overview$n_patients,
    "\n"
  )

  cat(
    "Rows:           ",
    x$overview$n_rows,
    "\n"
  )

  cat(
    "Timepoints:     ",
    x$overview$n_timepoints,
    "\n"
  )

  cat(
    "Time variables: ",
    paste(
      x$overview$timepoints,
      collapse = ", "
    ),
    "\n"
  )

  cat(
    "Complete:       ",
    sprintf(
      "%.1f%%",
      x$overview$complete_profiles_pct
    ),
    "\n"
  )

  cat(
    "Duplicated IDs:  ",
    x$overview$duplicated_ids,
    "\n"
  )


  # ----------------------------------------------------------
  # DESCRIPTIVES
  # ----------------------------------------------------------

  cat("\n\nDESCRIPTIVE STATISTICS\n")

  cat(
    "------------------------------------------------------------\n"
  )

  d <- x$descriptives


  print(

    data.frame(

      time =
        d$time,

      n =
        d$n,

      missing =
        d$missing,

      mean =
        round(d$mean, 3),

      sd =
        round(d$sd, 3),

      median =
        round(d$median, 3),

      q1 =
        round(d$q1, 3),

      q3 =
        round(d$q3, 3),

      min =
        round(d$min, 3),

      max =
        round(d$max, 3),

      ci_lower =
        round(d$ci_lower, 3),

      ci_upper =
        round(d$ci_upper, 3),

      stringsAsFactors =
        FALSE
    ),

    row.names = FALSE
  )


  # ----------------------------------------------------------
  # CHANGE
  # ----------------------------------------------------------

  cat("\n\nLONGITUDINAL CHANGES\n")

  cat(
    "------------------------------------------------------------\n"
  )


  if (
    nrow(x$change) > 0
  ) {

    ch <- x$change


    print(

      data.frame(

        from =
          ch$from,

        to =
          ch$to,

        n =
          ch$n,

        mean_change =
          round(
            ch$mean_change,
            4
          ),

        ci_lower =
          round(
            ch$ci_lower,
            4
          ),

        ci_upper =
          round(
            ch$ci_upper,
            4
          ),

        cohens_dz =
          round(
            ch$cohens_dz,
            4
          ),

        improved_pct =
          round(
            ch$improved_pct,
            1
          ),

        worsened_pct =
          round(
            ch$worsened_pct,
            1
          ),

        paired_t_p =
          signif(
            ch$paired_t_p,
            4
          ),

        wilcoxon_p =
          signif(
            ch$wilcoxon_p,
            4
          ),

        stringsAsFactors =
          FALSE
      ),

      row.names = FALSE
    )

  } else {

    cat(
      "Nessun confronto longitudinale disponibile.\n"
    )
  }


  # ----------------------------------------------------------
  # CORRELATIONS
  # ----------------------------------------------------------

  if (
    !is.null(
      x$correlations$pearson
    )
  ) {

    cat("\n\nPEARSON CORRELATIONS\n")

    cat(
      "------------------------------------------------------------\n"
    )

    print(
      round(
        x$correlations$pearson,
        3
      )
    )
  }


  # ----------------------------------------------------------
  # VARIABILITY
  # ----------------------------------------------------------

  cat("\n\nVARIABILITY\n")

  cat(
    "------------------------------------------------------------\n"
  )

  print(

    data.frame(

      grand_mean =
        round(
          x$variability$grand_mean,
          4
        ),

      between_subject_sd =
        round(
          x$variability$between_subject_sd,
          4
        ),

      within_subject_sd =
        round(
          x$variability$within_subject_sd,
          4
        ),

      ICC =
        round(
          x$variability$ICC,
          4
        ),

      stringsAsFactors =
        FALSE
    ),

    row.names = FALSE
  )


  # ----------------------------------------------------------
  # MIXED MODEL
  # ----------------------------------------------------------

  if (
    !is.null(
      x$model$fitted_model
    )
  ) {

    cat("\n\nMIXED MODEL\n")

    cat(
      "------------------------------------------------------------\n"
    )

    cat(
      "Model: value ~ time + (1 | patient)\n"
    )

    print(
      summary(
        x$model$fitted_model
      )
    )


    if (
      !is.null(
        x$model$anova
      )
    ) {

      cat("\n\nGLOBAL TEST OF TIME\n")

      cat(
        "------------------------------------------------------------\n"
      )

      print(
        x$model$anova
      )
    }

  } else if (
    !is.null(
      x$model$error
    )
  ) {

    cat("\n\nMIXED MODEL\n")

    cat(
      "------------------------------------------------------------\n"
    )

    cat(
      "Model non stimato: ",
      x$model$error,
      "\n"
    )
  }


  # ----------------------------------------------------------
  # OUTLIERS
  # ----------------------------------------------------------

  if (
    !is.null(
      x$outliers$by_time
    )
  ) {

    cat("\n\nOUTLIERS\n")

    cat(
      "------------------------------------------------------------\n"
    )


    for (
      nm in names(
        x$outliers$by_time
      )
    ) {

      tmp <-
        x$outliers$by_time[[nm]]

      cat(

        nm,

        ": ",

        nrow(tmp),

        " outlier(s)\n",

        sep = ""
      )
    }
  }


  # ----------------------------------------------------------
  # END
  # ----------------------------------------------------------

  cat("\n")

  cat(
    "============================================================\n"
  )

  cat(
    "Use names(result) to explore the complete output.\n"
  )

  cat(
    "============================================================\n\n"
  )


  invisible(x)
}
