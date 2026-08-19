#' Research-oriented summary of a MIRA Bayesian fit
#'
#' @param fit A fitted CmdStanMCMC object.
#' @param meaningful_change Minimum clinically important difference.
#' @param y Observed outcome vector.
#' @param credible_level Width of the posterior credible interval.
#' @param responder_thresholds Posterior probability thresholds.
#'
#' @return A list containing population, clinical, individual,
#' responders, heterogeneity, PPC, LOO information and diagnostics.
#'
#' @export
mira_summary <- function(
    fit,
    meaningful_change,
    y,
    credible_level = 0.90,
    responder_thresholds = c(0.50, 0.80, 0.95)
) {

  # ============================================================
  # VALIDATION
  # ============================================================

  if (!inherits(fit, "CmdStanMCMC")) {
    stop(
      "'fit' must be a fitted CmdStanMCMC object.",
      call. = FALSE
    )
  }


  if (
    length(meaningful_change) != 1 ||
    !is.numeric(meaningful_change) ||
    !is.finite(meaningful_change)
  ) {
    stop(
      "'meaningful_change' must be a single finite numeric value.",
      call. = FALSE
    )
  }


  if (
    !is.numeric(y) ||
    length(y) == 0 ||
    any(!is.finite(y))
  ) {
    stop(
      "'y' must be a non-empty numeric vector containing finite values.",
      call. = FALSE
    )
  }


  if (
    length(credible_level) != 1 ||
    credible_level <= 0 ||
    credible_level >= 1
  ) {
    stop(
      "'credible_level' must be between 0 and 1.",
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
      "'responder_thresholds' must contain values between 0 and 1.",
      call. = FALSE
    )
  }


  # ============================================================
  # HELPER
  # ============================================================

  posterior_stats <- function(x) {

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


  get_cols <- function(
    x,
    prefix
  ) {

    grep(
      paste0("^", prefix, "\\["),
      names(x),
      value = TRUE
    )
  }


  # ============================================================
  # VARIABLES REQUIRED BY THE CURRENT STAN MODEL
  # ============================================================

  draw_variables <- c(

    # population
    "mu_time",
    "change_T1",
    "change_T2",
    "change_T2_vs_T1",
    "population_slope",
    "sigma_intercept",
    "sigma_slope",
    "rho_subject",
    "sigma",
    "nu",

    # clinical population
    "population_meaningful_improvement_T2",
    "population_any_improvement_T2",
    "population_change_minus_MCID_T2",

    # individual
    "individual_change_T1",
    "individual_change_T2",
    "individual_change_T2_vs_T1",
    "individual_any_improvement_T2",
    "individual_meaningful_improvement_T2",
    "individual_change_minus_MCID_T2",

    # predictive / likelihood
    "log_lik",
    "y_rep"
  )


  # ============================================================
  # EXTRACT DRAWS
  # ============================================================

  draws_raw <- tryCatch(

    fit$draws(
      variables = draw_variables,
      format = "draws_matrix"
    ),

    error = function(e) {

      stop(
        paste0(
          "Could not extract posterior draws. ",
          "This usually means that the fitted Stan model ",
          "does not contain one or more required generated ",
          "quantities.\n\nOriginal error: ",
          conditionMessage(e)
        ),
        call. = FALSE
      )
    }
  )


  draws <-
    as.data.frame(draws_raw)


  # ============================================================
  # CHECK EXTRACTED VARIABLES
  # ============================================================

  required_scalar_variables <- c(

    "change_T1",
    "change_T2",
    "change_T2_vs_T1",
    "population_slope",
    "sigma_intercept",
    "sigma_slope",
    "rho_subject",
    "sigma",
    "nu",
    "population_change_minus_MCID_T2"
  )


  missing_scalars <-
    setdiff(
      required_scalar_variables,
      names(draws)
    )


  if (length(missing_scalars) > 0) {

    stop(
      paste0(
        "The fitted Stan model is missing: ",
        paste(
          missing_scalars,
          collapse = ", "
        ),
        ".\n\n",
        "Make sure `fit` was generated using the current ",
        "MIRA Stan model."
      ),
      call. = FALSE
    )
  }


  # ============================================================
  # POPULATION SUMMARY
  # ============================================================

  population_targets <- c(

    "change_T1",
    "change_T2",
    "change_T2_vs_T1",
    "population_slope",
    "sigma_intercept",
    "sigma_slope",
    "rho_subject",
    "sigma",
    "nu",
    "population_change_minus_MCID_T2"
  )


  population_summary <- do.call(

    rbind,

    lapply(

      population_targets,

      function(parameter) {

        s <-
          posterior_stats(
            draws[[parameter]]
          )

        data.frame(

          parameter = parameter,

          mean = s["mean"],

          median = s["median"],

          sd = s["sd"],

          mad = s["mad"],

          lower = s["lower"],

          upper = s["upper"],

          CrI_width = s["CrI_width"],

          stringsAsFactors = FALSE
        )
      }
    )
  )


  rownames(population_summary) <-
    NULL


  # ============================================================
  # POPULATION PROBABILITIES
  # ============================================================

  population_summary$P_positive <-
    NA_real_


  population_summary$P_MCID <-
    NA_real_


  population_summary$P_negative <-
    NA_real_


  population_summary$P_positive[
    population_summary$parameter ==
      "change_T1"
  ] <-
    mean(
      draws$change_T1 > 0
    )


  population_summary$P_positive[
    population_summary$parameter ==
      "change_T2"
  ] <-
    mean(
      draws$change_T2 > 0
    )


  population_summary$P_positive[
    population_summary$parameter ==
      "change_T2_vs_T1"
  ] <-
    mean(
      draws$change_T2_vs_T1 > 0
    )


  population_summary$P_negative[
    population_summary$parameter ==
      "change_T1"
  ] <-
    mean(
      draws$change_T1 < 0
    )


  population_summary$P_negative[
    population_summary$parameter ==
      "change_T2"
  ] <-
    mean(
      draws$change_T2 < 0
    )


  population_summary$P_negative[
    population_summary$parameter ==
      "change_T2_vs_T1"
  ] <-
    mean(
      draws$change_T2_vs_T1 < 0
    )


  population_summary$P_MCID[
    population_summary$parameter ==
      "change_T2"
  ] <-
    mean(
      draws$change_T2 >=
        meaningful_change
    )


  # ============================================================
  # CLINICAL SUMMARY
  # ============================================================

  probability_any_improvement <-
    mean(
      draws$change_T2 > 0
    )


  probability_MCID <-
    mean(
      draws$change_T2 >=
        meaningful_change
    )


  distance_MCID <-
    draws$population_change_minus_MCID_T2


  clinical_summary <-
    tibble::tibble(

      measure = c(

        "MCID",

        "P(population change T2 > 0)",

        "P(population change T2 >= MCID)",

        "P(population change T2 < MCID)",

        "Mean population distance from MCID",

        "P(population distance from MCID >= 0)"
      ),

      estimate = c(

        meaningful_change,

        probability_any_improvement,

        probability_MCID,

        1 - probability_MCID,

        mean(distance_MCID),

        mean(distance_MCID >= 0)
      )
    )


  # ============================================================
  # STANDARDIZED POPULATION CHANGE
  # ============================================================

  standardized_change_T2 <-
    draws$change_T2 /
    stats::sd(y)


  standardized_stats <-
    posterior_stats(
      standardized_change_T2
    )


  clinical_summary <-
    tibble::add_row(

      clinical_summary,

      measure =
        "Standardized population change T2",

      estimate =
        standardized_stats["mean"]
    )


  # ============================================================
  # INDIVIDUAL VARIABLES
  # ============================================================

  individual_T1_cols <-
    get_cols(
      draws,
      "individual_change_T1"
    )


  individual_T2_cols <-
    get_cols(
      draws,
      "individual_change_T2"
    )


  individual_T2_vs_T1_cols <-
    get_cols(
      draws,
      "individual_change_T2_vs_T1"
    )


  individual_any_cols <-
    get_cols(
      draws,
      "individual_any_improvement_T2"
    )


  individual_MCID_cols <-
    get_cols(
      draws,
      "individual_meaningful_improvement_T2"
    )


  individual_distance_cols <-
    get_cols(
      draws,
      "individual_change_minus_MCID_T2"
    )


  required_groups <- list(

    individual_change_T1 =
      individual_T1_cols,

    individual_change_T2 =
      individual_T2_cols,

    individual_change_T2_vs_T1 =
      individual_T2_vs_T1_cols,

    individual_any_improvement_T2 =
      individual_any_cols,

    individual_meaningful_improvement_T2 =
      individual_MCID_cols,

    individual_change_minus_MCID_T2 =
      individual_distance_cols
  )


  missing_groups <-
    names(required_groups)[
      vapply(
        required_groups,
        length,
        integer(1)
      ) == 0
    ]


  if (length(missing_groups) > 0) {

    stop(
      paste0(
        "Missing individual generated quantities: ",
        paste(
          missing_groups,
          collapse = ", "
        )
      ),
      call. = FALSE
    )
  }


  # ============================================================
  # NUMBER OF SUBJECTS
  # ============================================================

  n_subjects <-
    length(individual_T2_cols)


  # ============================================================
  # INDIVIDUAL SUMMARY
  # ============================================================

  individual_summary <-
    tibble::tibble(
      subject = seq_len(n_subjects)
    )


  add_stats <- function(
    data,
    columns,
    prefix
  ) {

    sm <- t(
      vapply(
        columns,
        function(column) {

          posterior_stats(
            draws[[column]]
          )
        },
        numeric(7)
      )
    )


    data[[paste0(prefix, "_mean")]] <-
      sm[, "mean"]


    data[[paste0(prefix, "_median")]] <-
      sm[, "median"]


    data[[paste0(prefix, "_sd")]] <-
      sm[, "sd"]


    data[[paste0(prefix, "_mad")]] <-
      sm[, "mad"]


    data[[paste0(prefix, "_lower")]] <-
      sm[, "lower"]


    data[[paste0(prefix, "_upper")]] <-
      sm[, "upper"]


    data[[paste0(prefix, "_CrI_width")]] <-
      sm[, "CrI_width"]


    data
  }


  individual_summary <-
    add_stats(
      individual_summary,
      individual_T1_cols,
      "change_T1"
    )


  individual_summary <-
    add_stats(
      individual_summary,
      individual_T2_cols,
      "change_T2"
    )


  individual_summary <-
    add_stats(
      individual_summary,
      individual_T2_vs_T1_cols,
      "change_T2_vs_T1"
    )


  individual_summary <-
    add_stats(
      individual_summary,
      individual_distance_cols,
      "distance_from_MCID"
    )


  # ============================================================
  # INDIVIDUAL CLINICAL PROBABILITIES
  # ============================================================

  individual_summary$P_improvement_T2 <-
    vapply(
      individual_any_cols,
      function(column) {

        mean(
          draws[[column]]
        )
      },
      numeric(1)
    )


  individual_summary$P_MCID_T2 <-
    vapply(
      individual_MCID_cols,
      function(column) {

        mean(
          draws[[column]]
        )
      },
      numeric(1)
    )


  individual_summary$P_negative_change_T2 <-
    vapply(
      individual_T2_cols,
      function(column) {

        mean(
          draws[[column]] < 0
        )
      },
      numeric(1)
    )


  # ============================================================
  # RESPONDER CLASSIFICATION
  # ============================================================

  individual_summary$responder_class <-
    cut(

      individual_summary$P_MCID_T2,

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


  individual_summary$high_probability_response_80 <-
    individual_summary$P_MCID_T2 >=
    0.80


  individual_summary$high_probability_response_95 <-
    individual_summary$P_MCID_T2 >=
    0.95


  # ============================================================
  # RESPONDER SUMMARY
  # ============================================================

  responders_summary <-
    tibble::tibble(

      posterior_probability_threshold =
        responder_thresholds,

      proportion_subjects =
        vapply(

          responder_thresholds,

          function(threshold) {

            mean(
              individual_summary$P_MCID_T2 >=
                threshold
            )
          },

          numeric(1)
        )
    )


  # Expected proportion of responders

  expected_responder_proportion <-
    mean(
      individual_summary$P_MCID_T2
    )


  # Add it as a separate quantity

  responders_summary <-
    tibble::add_row(

      responders_summary,

      posterior_probability_threshold =
        NA_real_,

      proportion_subjects =
        expected_responder_proportion
    )


  # ============================================================
  # HETEROGENEITY
  # ============================================================

  heterogeneity_parameters <- c(
    "sigma_intercept",
    "sigma_slope",
    "rho_subject"
  )


  heterogeneity <-
    do.call(

      rbind,

      lapply(

        heterogeneity_parameters,

        function(parameter) {

          s <-
            posterior_stats(
              draws[[parameter]]
            )


          data.frame(

            parameter = parameter,

            mean = s["mean"],

            median = s["median"],

            sd = s["sd"],

            mad = s["mad"],

            lower = s["lower"],

            upper = s["upper"],

            CrI_width =
              s["CrI_width"],

            stringsAsFactors = FALSE
          )
        }
      )
    )


  rownames(heterogeneity) <-
    NULL


  heterogeneity$P_positive <-
    NA_real_


  heterogeneity$P_negative <-
    NA_real_


  heterogeneity$P_positive[
    heterogeneity$parameter ==
      "rho_subject"
  ] <-
    mean(
      draws$rho_subject > 0
    )


  heterogeneity$P_negative[
    heterogeneity$parameter ==
      "rho_subject"
  ] <-
    mean(
      draws$rho_subject < 0
    )


  # ============================================================
  # PPC
  # ============================================================

  y_rep_cols <-
    get_cols(
      draws,
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
      )
    )
  }


  # ============================================================
  # LOG LIKELIHOOD
  # ============================================================

  log_lik_cols <-
    get_cols(
      draws,
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
      fit$summary(
        variables =
          population_targets
      )


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
      n_subjects,

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

    expected_responder_proportion =
      expected_responder_proportion
  )


  # ============================================================
  # RETURN
  # ============================================================

  list(

    population =
      population_summary,

    clinical =
      clinical_summary,

    individual =
      individual_summary,

    responders =
      responders_summary,

    heterogeneity =
      heterogeneity,

    ppc =
      ppc,

    loo =
      loo_information,

    diagnostics =
      diagnostics,

    quality_flags =
      quality_flags,

    model_information =
      model_information,

    draws =
      draws
  )
}
