#' Research-oriented summary of a MIRA Bayesian fit
#'
#' Produces a comprehensive research-oriented summary of a fitted
#' MIRA longitudinal Student-t mixed-effects model.
#'
#' The function is fully dynamic with respect to the number of
#' longitudinal measurement occasions K and the number of subjects S.
#'
#' @param fit A fitted CmdStanMCMC object.
#' @param stan_data Optional Stan data list returned by
#'   [mira_prepare_data()]. If supplied, it is used to recover
#'   N, S, K, time_value, y and MCID.
#' @param meaningful_change Optional MCID. If NULL, it is taken
#'   from `stan_data`.
#' @param y Optional observed outcome vector. If NULL, it is taken
#'   from `stan_data`.
#' @param credible_level Width of posterior credible intervals.
#' @param responder_thresholds Posterior probability thresholds
#'   used for responder classification.
#'
#' @return A list containing population summaries, longitudinal
#'   changes, clinical quantities, individual summaries, responder
#'   analyses, heterogeneity, posterior predictive checks,
#'   log-likelihood information, diagnostics, model information
#'   and posterior draws.
#'
#' @export
mira_summary <- function(
    fit,
    stan_data = NULL,
    meaningful_change = NULL,
    y = NULL,
    credible_level = 0.90,
    responder_thresholds = c(
      0.50,
      0.80,
      0.95
    )
) {

  # ============================================================
  # VALIDATION
  # ============================================================

  if (!inherits(fit, "CmdStanMCMC")) {

    stop(
      "`fit` must be a fitted CmdStanMCMC object.",
      call. = FALSE
    )
  }


  if (
    length(credible_level) != 1 ||
    !is.numeric(credible_level) ||
    !is.finite(credible_level) ||
    credible_level <= 0 ||
    credible_level >= 1
  ) {

    stop(
      "`credible_level` must be between 0 and 1.",
      call. = FALSE
    )
  }


  if (
    length(responder_thresholds) == 0 ||
    any(!is.finite(responder_thresholds)) ||
    any(
      responder_thresholds <= 0 |
      responder_thresholds >= 1
    )
  ) {

    stop(
      "`responder_thresholds` must contain values between 0 and 1.",
      call. = FALSE
    )
  }


  # ============================================================
  # RECOVER DATA
  # ============================================================

  if (!is.null(stan_data)) {

    if (!is.list(stan_data)) {

      stop(
        "`stan_data` must be a list.",
        call. = FALSE
      )
    }


    if (is.null(y) && "y" %in% names(stan_data)) {

      y <- stan_data$y
    }


    if (
      is.null(meaningful_change) &&
      "meaningful_change" %in% names(stan_data)
    ) {

      meaningful_change <-
        stan_data$meaningful_change
    }
  }


  if (is.null(y)) {

    stop(
      "Observed outcome `y` must be supplied either directly ",
      "or through `stan_data`.",
      call. = FALSE
    )
  }


  if (
    !is.numeric(y) ||
    length(y) == 0 ||
    any(!is.finite(y))
  ) {

    stop(
      "`y` must be a non-empty numeric vector containing finite values.",
      call. = FALSE
    )
  }


  if (is.null(meaningful_change)) {

    stop(
      "`meaningful_change` must be supplied either directly ",
      "or through `stan_data`.",
      call. = FALSE
    )
  }


  if (
    length(meaningful_change) != 1 ||
    !is.numeric(meaningful_change) ||
    !is.finite(meaningful_change)
  ) {

    stop(
      "`meaningful_change` must be a single finite numeric value.",
      call. = FALSE
    )
  }


  # ============================================================
  # EXTRACT ALL POSTERIOR DRAWS
  #
  # This is intentional:
  #
  # the current model contains many dynamic quantities whose
  # names depend on K and S.
  #
  # ============================================================

  draws_raw <- tryCatch(

    fit$draws(
      format = "draws_matrix"
    ),

    error = function(e) {

      stop(
        paste0(
          "Could not extract posterior draws.\n\n",
          "Original error: ",
          conditionMessage(e)
        ),
        call. = FALSE
      )
    }
  )


  draws <- as.data.frame(
    draws_raw
  )


  # ============================================================
  # HELPERS
  # ============================================================

  get_cols <- function(
    prefix
  ) {

    grep(
      paste0("^", prefix, "\\["),
      names(draws),
      value = TRUE
    )
  }


  posterior_stats <- function(
    x
  ) {

    x <- as.numeric(x)

    alpha <-
      (1 - credible_level) / 2

    q <- stats::quantile(
      x,
      probs = c(
        alpha,
        1 - alpha
      ),
      names = FALSE
    )

    c(
      mean = mean(x),
      median = stats::median(x),
      sd = stats::sd(x),
      mad = stats::mad(x),
      lower = q[1],
      upper = q[2],
      CrI_width = q[2] - q[1]
    )
  }


  summarize_vector <- function(
    x,
    parameter
  ) {

    s <- posterior_stats(x)

    data.frame(

      parameter = parameter,

      mean = unname(s["mean"]),

      median = unname(s["median"]),

      sd = unname(s["sd"]),

      mad = unname(s["mad"]),

      lower = unname(s["lower"]),

      upper = unname(s["upper"]),

      CrI_width = unname(s["CrI_width"]),

      stringsAsFactors = FALSE
    )
  }


  summarize_columns <- function(
    columns,
    prefix
  ) {

    if (length(columns) == 0) {

      return(
        tibble::tibble()
      )
    }


    result <- lapply(

      seq_along(columns),

      function(i) {

        s <-
          posterior_stats(
            draws[[columns[i]]]
          )

        data.frame(

          index = i,

          mean = unname(s["mean"]),

          median = unname(s["median"]),

          sd = unname(s["sd"]),

          mad = unname(s["mad"]),

          lower = unname(s["lower"]),

          upper = unname(s["upper"]),

          CrI_width =
            unname(s["CrI_width"]),

          stringsAsFactors = FALSE
        )
      }
    )


    result <-
      dplyr::bind_rows(result)


    result
  }


  # ============================================================
  # MODEL DIMENSIONS
  # ============================================================

  mu_time_cols <-
    get_cols("mu_time")


  K <-
    length(mu_time_cols)


  if (K < 2) {

    stop(
      "Could not determine at least two measurement occasions ",
      "from posterior variable `mu_time`.",
      call. = FALSE
    )
  }


  individual_change_cols <-
    get_cols(
      "individual_change_from_baseline"
    )


  if (length(individual_change_cols) == 0) {

    stop(
      "The fitted model does not contain ",
      "`individual_change_from_baseline`.",
      call. = FALSE
    )
  }


  # ------------------------------------------------------------
  # Infer S from individual slope
  # ------------------------------------------------------------

  individual_slope_cols <-
    get_cols(
      "individual_slope"
    )


  S <-
    length(individual_slope_cols)


  if (S < 1) {

    stop(
      "Could not determine the number of subjects.",
      call. = FALSE
    )
  }


  # ============================================================
  # TIME VALUES
  # ============================================================

  if (
    !is.null(stan_data) &&
    "time_value" %in% names(stan_data)
  ) {

    time_value <-
      as.numeric(
        stan_data$time_value
      )

  } else {

    time_value <-
      seq_len(K) - 1
  }


  if (length(time_value) != K) {

    stop(
      "`time_value` must contain exactly K values.",
      call. = FALSE
    )
  }


  # ============================================================
  # POPULATION PARAMETERS
  # ============================================================

  population_scalar_parameters <- c(

    "sigma",

    "nu",

    "sigma_intercept",

    "sigma_slope",

    "rho_subject"
  )


  population_scalar_parameters <-
    population_scalar_parameters[
      population_scalar_parameters %in%
        names(draws)
    ]


  population_summary <-
    dplyr::bind_rows(

      lapply(

        population_scalar_parameters,

        function(parameter) {

          summarize_vector(
            draws[[parameter]],
            parameter
          )
        }
      )
    )


  # ============================================================
  # POPULATION TIME MEANS
  # ============================================================

  mu_time_summary <-
    summarize_columns(
      mu_time_cols,
      "mu_time"
    )


  mu_time_summary$time <-
    seq_len(K)


  mu_time_summary$time_value <-
    time_value


  mu_time_summary <-
    mu_time_summary[
      c(
        "index",
        "time",
        "time_value",
        "mean",
        "median",
        "sd",
        "mad",
        "lower",
        "upper",
        "CrI_width"
      )
    ]


  # ============================================================
  # POPULATION CHANGES FROM BASELINE
  # ============================================================

  change_cols <-
    get_cols(
      "change_from_baseline"
    )


  population_change_summary <-
    summarize_columns(
      change_cols,
      "change_from_baseline"
    )


  population_change_summary$time <-
    seq_len(K)


  population_change_summary$time_value <-
    time_value


  population_change_summary$P_positive <-
    vapply(

      change_cols,

      function(column) {

        mean(
          draws[[column]] > 0
        )
      },

      numeric(1)
    )


  population_change_summary$P_negative <-
    vapply(

      change_cols,

      function(column) {

        mean(
          draws[[column]] < 0
        )
      },

      numeric(1)
    )


  population_change_summary$P_MCID <-
    vapply(

      change_cols,

      function(column) {

        mean(
          draws[[column]] >=
            meaningful_change
        )
      },

      numeric(1)
    )


  population_change_summary$distance_from_MCID_mean <-
    vapply(

      change_cols,

      function(column) {

        mean(
          draws[[column]] -
            meaningful_change
        )
      },

      numeric(1)
    )


  population_change_summary <-
    population_change_summary[
      c(
        "index",
        "time",
        "time_value",
        "mean",
        "median",
        "sd",
        "mad",
        "lower",
        "upper",
        "CrI_width",
        "P_positive",
        "P_negative",
        "P_MCID",
        "distance_from_MCID_mean"
      )
    ]


  # ============================================================
  # CONSECUTIVE POPULATION CHANGES
  # ============================================================

  consecutive_cols <-
    get_cols(
      "consecutive_change"
    )


  consecutive_summary <-
    summarize_columns(
      consecutive_cols,
      "consecutive_change"
    )


  if (length(consecutive_cols) > 0) {

    consecutive_summary$from_time <-
      seq_len(K - 1)

    consecutive_summary$to_time <-
      seq(2, K)

    consecutive_summary$from_time_value <-
      time_value[seq_len(K - 1)]

    consecutive_summary$to_time_value <-
      time_value[seq(2, K)]

    consecutive_summary$P_positive <-
      vapply(

        consecutive_cols,

        function(column) {

          mean(
            draws[[column]] > 0
          )
        },

        numeric(1)
      )


    consecutive_summary$P_negative <-
      vapply(

        consecutive_cols,

        function(column) {

          mean(
            draws[[column]] < 0
          )
        },

        numeric(1)
      )
  }


  # ============================================================
  # POPULATION SLOPES
  # ============================================================

  slope_cols <-
    get_cols(
      "population_slope_from_baseline"
    )


  population_slope_summary <-
    summarize_columns(
      slope_cols,
      "population_slope_from_baseline"
    )


  if (length(slope_cols) > 0) {

    population_slope_summary$time <-
      seq_len(K)

    population_slope_summary$time_value <-
      time_value
  }


  # ============================================================
  # STANDARDIZED CHANGE
  # ============================================================

  standardized_cols <-
    get_cols(
      "standardized_change_from_baseline"
    )


  standardized_change_summary <-
    summarize_columns(
      standardized_cols,
      "standardized_change_from_baseline"
    )


  if (length(standardized_cols) > 0) {

    standardized_change_summary$time <-
      seq_len(K)

    standardized_change_summary$time_value <-
      time_value
  }


  # ============================================================
  # POPULATION CLINICAL RESPONSE
  # ============================================================

  population_any_cols <-
    get_cols(
      "population_any_improvement"
    )


  population_MCID_cols <-
    get_cols(
      "population_meaningful_improvement"
    )


  population_distance_cols <-
    get_cols(
      "population_change_minus_MCID"
    )


  population_clinical <- tibble::tibble(

    time =
      seq_len(K),

    time_value =
      time_value,

    P_any_improvement =
      vapply(

        population_any_cols,

        function(column) {

          mean(
            draws[[column]]
          )
        },

        numeric(1)
      ),

    P_meaningful_improvement =
      vapply(

        population_MCID_cols,

        function(column) {

          mean(
            draws[[column]]
          )
        },

        numeric(1)
      ),

    mean_change_minus_MCID =
      vapply(

        population_distance_cols,

        function(column) {

          mean(
            draws[[column]]
          )
        },

        numeric(1)
      ),

    P_change_minus_MCID_ge_0 =
      vapply(

        population_distance_cols,

        function(column) {

          mean(
            draws[[column]] >= 0
          )
        },

        numeric(1)
      )
  )


  # ============================================================
  # INDIVIDUAL SLOPES
  # ============================================================

  individual_slope_summary <-
    summarize_columns(
      individual_slope_cols,
      "individual_slope"
    )


  individual_slope_summary$subject <-
    seq_len(S)


  individual_slope_summary <-
    individual_slope_summary[
      c(
        "subject",
        "mean",
        "median",
        "sd",
        "mad",
        "lower",
        "upper",
        "CrI_width"
      )
    ]


  # ============================================================
  # INDIVIDUAL CHANGE FROM BASELINE
  # ============================================================

  individual_change_summary <-
    summarize_columns(
      individual_change_cols,
      "individual_change_from_baseline"
    )


  # Parse [time,subject] from Stan column names

  individual_indices <-
    regmatches(

      individual_change_cols,

      regexec(
        "\\[([0-9]+),([0-9]+)\\]",
        individual_change_cols
      )
    )


  individual_change_summary$time <-
    vapply(
      individual_indices,
      function(x) as.integer(x[2]),
      integer(1)
    )


  individual_change_summary$subject <-
    vapply(
      individual_indices,
      function(x) as.integer(x[3]),
      integer(1)
    )


  individual_change_summary$time_value <-
    time_value[
      individual_change_summary$time
    ]


  # ============================================================
  # INDIVIDUAL CLINICAL PROBABILITIES
  # ============================================================

  individual_any_cols <-
    get_cols(
      "individual_any_improvement"
    )


  individual_MCID_cols <-
    get_cols(
      "individual_meaningful_improvement"
    )


  individual_distance_cols <-
    get_cols(
      "individual_change_minus_MCID"
    )


  individual_negative_cols <-
    individual_change_cols


  individual_clinical_summary <-
    data.frame(

      subject = integer(0),

      time = integer(0),

      time_value = numeric(0),

      P_improvement = numeric(0),

      P_MCID = numeric(0),

      P_negative_change = numeric(0),

      mean_change_minus_MCID = numeric(0),

      P_change_minus_MCID_ge_0 = numeric(0),

      stringsAsFactors = FALSE
    )


  if (length(individual_any_cols) > 0) {

    indices <-
      regmatches(

        individual_any_cols,

        regexec(
          "\\[([0-9]+),([0-9]+)\\]",
          individual_any_cols
        )
      )


    individual_clinical_summary <-
      data.frame(

        subject =
          vapply(
            indices,
            function(x) as.integer(x[3]),
            integer(1)
          ),

        time =
          vapply(
            indices,
            function(x) as.integer(x[2]),
            integer(1)
          ),

        P_improvement =
          vapply(
            individual_any_cols,
            function(column) {

              mean(
                draws[[column]]
              )
            },
            numeric(1)
          ),

        P_MCID =
          vapply(
            individual_MCID_cols,
            function(column) {

              mean(
                draws[[column]]
              )
            },
            numeric(1)
          ),

        P_negative_change =
          vapply(
            individual_negative_cols,
            function(column) {

              mean(
                draws[[column]] < 0
              )
            },
            numeric(1)
          ),

        mean_change_minus_MCID =
          vapply(
            individual_distance_cols,
            function(column) {

              mean(
                draws[[column]]
              )
            },
            numeric(1)
          ),

        P_change_minus_MCID_ge_0 =
          vapply(
            individual_distance_cols,
            function(column) {

              mean(
                draws[[column]] >= 0
              )
            },
            numeric(1)
          ),

        stringsAsFactors = FALSE
      )


    individual_clinical_summary$time_value <-
      time_value[
        individual_clinical_summary$time
      ]


    individual_clinical_summary <-
      individual_clinical_summary[
        order(
          individual_clinical_summary$subject,
          individual_clinical_summary$time
        ),
      ]
  }


  # ============================================================
  # INDIVIDUAL RESPONDER CLASSIFICATION
  # ============================================================

  if (nrow(individual_clinical_summary) > 0) {

    individual_clinical_summary$responder_class <-

      cut(

        individual_clinical_summary$P_MCID,

        breaks = c(
          -Inf,
          0.20,
          0.50,
          0.80,
          0.95,
          Inf
        ),

        labels = c(
          "very_unlikely",
          "uncertain",
          "probable",
          "high_probability",
          "very_high_probability"
        ),

        right = FALSE
      )


    individual_clinical_summary$response_80 <-
      individual_clinical_summary$P_MCID >= 0.80


    individual_clinical_summary$response_95 <-
      individual_clinical_summary$P_MCID >= 0.95
  }


  # ============================================================
  # RESPONDER PROPORTION
  # ============================================================

  responder_proportion_cols <-
    get_cols(
      "population_responder_proportion"
    )


  responder_proportion_summary <-
    summarize_columns(
      responder_proportion_cols,
      "population_responder_proportion"
    )


  if (length(responder_proportion_cols) > 0) {

    responder_proportion_summary$time <-
      seq_len(K)

    responder_proportion_summary$time_value <-
      time_value
  }


  # ============================================================
  # RESPONDER SUMMARY BY TIME
  # ============================================================

  responder_summary <- tibble::tibble(

    time =
      seq_len(K),

    time_value =
      time_value,

    expected_responder_proportion =
      NA_real_,

    proportion_P_MCID_ge_50 =
      NA_real_,

    proportion_P_MCID_ge_80 =
      NA_real_,

    proportion_P_MCID_ge_95 =
      NA_real_
  )


  if (nrow(individual_clinical_summary) > 0) {

    for (k in seq_len(K)) {

      x <-
        individual_clinical_summary$P_MCID[
          individual_clinical_summary$time == k
        ]


      responder_summary$expected_responder_proportion[k] <-
        mean(x)


      responder_summary$proportion_P_MCID_ge_50[k] <-
        mean(x >= 0.50)


      responder_summary$proportion_P_MCID_ge_80[k] <-
        mean(x >= 0.80)


      responder_summary$proportion_P_MCID_ge_95[k] <-
        mean(x >= 0.95)
    }
  }


  # ============================================================
  # HETEROGENEITY
  # ============================================================

  heterogeneity_parameters <- c(

    "sigma_intercept",

    "sigma_slope",

    "rho_subject"
  )


  heterogeneity_parameters <-
    heterogeneity_parameters[
      heterogeneity_parameters %in%
        names(draws)
    ]


  heterogeneity <-
    dplyr::bind_rows(

      lapply(

        heterogeneity_parameters,

        function(parameter) {

          summarize_vector(
            draws[[parameter]],
            parameter
          )
        }
      )
    )


  if ("rho_subject" %in% names(draws)) {

    heterogeneity$P_positive <-
      NA_real_

    heterogeneity$P_negative <-
      NA_real_


    i <-
      heterogeneity$parameter ==
      "rho_subject"


    heterogeneity$P_positive[i] <-
      mean(
        draws$rho_subject > 0
      )


    heterogeneity$P_negative[i] <-
      mean(
        draws$rho_subject < 0
      )
  }


  # ============================================================
  # CLINICAL SUMMARY
  # ============================================================

  clinical_summary <-
    tibble::tibble(

      measure = c(

        "MCID",

        "Final time point",

        "Final population change",

        "P(final population change > 0)",

        "P(final population change >= MCID)",

        "P(final population change < MCID)",

        "Mean final distance from MCID",

        "P(final distance from MCID >= 0)"
      ),

      estimate = c(

        meaningful_change,

        K,

        mean(draws[[change_cols[K]]]),

        mean(draws[[change_cols[K]]] > 0),

        mean(draws[[change_cols[K]]] >= meaningful_change),

        mean(draws[[change_cols[K]]] < meaningful_change),

        mean(draws[[change_cols[K]]] -  meaningful_change),

        mean(draws[[change_cols[K]]] >= meaningful_change
        )
      )
    )


  # ============================================================
  # STANDARDIZED FINAL CHANGE
  # ============================================================

  final_standardized_change <- draws[[standardized_cols[K]]]


  clinical_summary <-
    tibble::add_row(

      clinical_summary,

      measure =
        "Standardized final population change",

      estimate =
        mean(
          final_standardized_change
        )
    )


  # ============================================================
  # PPC
  # ============================================================

  y_rep_cols <-
    get_cols(
      "y_rep"
    )


  ppc <- NULL


  if (length(y_rep_cols) > 0) {

    y_rep_matrix <-
      as.matrix(
        draws[
          ,
          y_rep_cols,
          drop = FALSE
        ]
      )


    predictive_mean <-
      rowMeans(
        y_rep_matrix
      )


    predictive_sd <-
      apply(
        y_rep_matrix,
        1,
        stats::sd
      )


    predictive_median <-
      apply(
        y_rep_matrix,
        1,
        stats::median
      )


    predictive_q05 <-
      apply(
        y_rep_matrix,
        1,
        stats::quantile,
        probs = 0.05
      )


    predictive_q95 <-
      apply(
        y_rep_matrix,
        1,
        stats::quantile,
        probs = 0.95
      )


    observed <- c(

      mean =
        mean(y),

      sd =
        stats::sd(y),

      median =
        stats::median(y),

      q05 =
        as.numeric(
          stats::quantile(
            y,
            0.05
          )
        ),

      q95 =
        as.numeric(
          stats::quantile(
            y,
            0.95
          )
        )
    )


    ppc <- list(

      observed = observed,

      posterior_predictive = list(

        mean =
          posterior_stats(
            predictive_mean
          ),

        sd =
          posterior_stats(
            predictive_sd
          ),

        median =
          posterior_stats(
            predictive_median
          ),

        q05 =
          posterior_stats(
            predictive_q05
          ),

        q95 =
          posterior_stats(
            predictive_q95
          )
      ),

      bayesian_p_values = c(

        mean =
          mean(
            predictive_mean >=
              mean(y)
          ),

        sd =
          mean(
            predictive_sd >=
              stats::sd(y)
          ),

        median =
          mean(
            predictive_median >=
              stats::median(y)
          )
      ),

      predictive_draws =
        y_rep_matrix
    )
  }


  # ============================================================
  # LOG LIKELIHOOD / LOO INFORMATION
  # ============================================================

  log_lik_cols <-
    get_cols(
      "log_lik"
    )


  loo_information <- NULL


  if (length(log_lik_cols) > 0) {

    log_lik_matrix <-
      as.matrix(
        draws[
          ,
          log_lik_cols,
          drop = FALSE
        ]
      )


    loo_information <- list(

      n_observations =
        length(log_lik_cols),

      mean_total_log_lik =
        mean(
          rowSums(
            log_lik_matrix
          )
        ),

      pointwise_mean_log_lik =
        colMeans(
          log_lik_matrix
        ),

      draws =
        log_lik_matrix
    )
  }


  # ============================================================
  # MCMC DIAGNOSTICS
  # ============================================================

  diagnostics <- tryCatch({

    s <-
      fit$summary()


    d <-
      data.frame(

        parameter =
          s$variable,

        Rhat =
          s$rhat,

        ESS_bulk =
          s$ess_bulk,

        ESS_tail =
          s$ess_tail,

        stringsAsFactors = FALSE
      )


    d$Rhat_ok <-
      d$Rhat < 1.01


    d$ESS_bulk_ok <-
      d$ESS_bulk >= 400


    d$ESS_tail_ok <-
      d$ESS_tail >= 400


    list(

      parameters = d,

      max_Rhat =
        max(
          d$Rhat,
          na.rm = TRUE
        ),

      min_ESS_bulk =
        min(
          d$ESS_bulk,
          na.rm = TRUE
        ),

      min_ESS_tail =
        min(
          d$ESS_tail,
          na.rm = TRUE
        )
    )

  }, error = function(e) {

    list(
      error =
        conditionMessage(e)
    )
  })


  # ============================================================
  # GLOBAL QUALITY FLAGS
  # ============================================================

  quality_flags <- list(

    Rhat_ok =
      if (
        !is.null(diagnostics$max_Rhat)
      ) {

        diagnostics$max_Rhat < 1.01

      } else {

        NA
      },

    ESS_bulk_ok =
      if (
        !is.null(diagnostics$min_ESS_bulk)
      ) {

        diagnostics$min_ESS_bulk >= 400

      } else {

        NA
      },

    ESS_tail_ok =
      if (
        !is.null(diagnostics$min_ESS_tail)
      ) {

        diagnostics$min_ESS_tail >= 400

      } else {

        NA
      }
  )


  # ============================================================
  # MODEL INFORMATION
  # ============================================================

  model_information <- list(

    n_observations =
      length(y),

    n_subjects =
      S,

    n_time_points =
      K,

    time_value =
      time_value,

    MCID =
      meaningful_change,

    credible_level =
      credible_level,

    responder_thresholds =
      responder_thresholds,

    observed_mean =
      mean(y),

    observed_sd =
      stats::sd(y),

    observed_median =
      stats::median(y),

    observed_q05 =
      as.numeric(
        stats::quantile(
          y,
          0.05
        )
      ),

    observed_q95 =
      as.numeric(
        stats::quantile(
          y,
          0.95
        )
      )
  )


  # ============================================================
  # ALL MODEL VARIABLES
  #
  # Useful for research/debugging:
  # explicitly expose which variables were actually generated
  # by Stan.
  # ============================================================

  variable_inventory <- list(

    population =
      c(
        mu_time_cols,
        change_cols,
        consecutive_cols,
        slope_cols,
        standardized_cols
      ),

    population_clinical =
      c(
        population_any_cols,
        population_MCID_cols,
        population_distance_cols,
        responder_proportion_cols
      ),

    individual =
      c(
        individual_slope_cols,
        individual_change_cols
      ),

    individual_clinical =
      c(
        individual_any_cols,
        individual_MCID_cols,
        individual_distance_cols
      ),

    observation_level =
      c(
        log_lik_cols,
        y_rep_cols
      )
  )


  # ============================================================
  # RETURN
  # ============================================================

  list(

    # ----------------------------------------------------------
    # Population-level results
    # ----------------------------------------------------------

    population =
      population_summary,

    population_time_means =
      mu_time_summary,

    population_change =
      population_change_summary,

    population_consecutive_change =
      consecutive_summary,

    population_slope =
      population_slope_summary,

    population_standardized_change =
      standardized_change_summary,

    population_clinical =
      population_clinical,

    # ----------------------------------------------------------
    # Individual results
    # ----------------------------------------------------------

    individual_slope =
      individual_slope_summary,

    individual_change =
      individual_change_summary,

    individual_clinical =
      individual_clinical_summary,

    # ----------------------------------------------------------
    # Responders
    # ----------------------------------------------------------

    responder_proportion =
      responder_proportion_summary,

    responders =
      responder_summary,

    # ----------------------------------------------------------
    # Clinical interpretation
    # ----------------------------------------------------------

    clinical =
      clinical_summary,

    # ----------------------------------------------------------
    # Heterogeneity
    # ----------------------------------------------------------

    heterogeneity =
      heterogeneity,

    # ----------------------------------------------------------
    # Posterior predictive checks
    # ----------------------------------------------------------

    ppc =
      ppc,

    # ----------------------------------------------------------
    # LOO / likelihood
    # ----------------------------------------------------------

    loo =
      loo_information,

    # ----------------------------------------------------------
    # Diagnostics
    # ----------------------------------------------------------

    diagnostics =
      diagnostics,

    quality_flags =
      quality_flags,

    # ----------------------------------------------------------
    # Model metadata
    # ----------------------------------------------------------

    model_information =
      model_information,

    variable_inventory =
      variable_inventory,

    # ----------------------------------------------------------
    # Raw posterior draws
    # ----------------------------------------------------------

    draws =
      draws
  )
}
