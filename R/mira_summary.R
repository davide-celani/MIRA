#' Research-oriented summary of a MIRA Bayesian fit
#'
#' Summarizes the current multi-arm MIRA longitudinal Student-t model.
#' The function is dynamic in the number of treatment arms G, measurement
#' occasions K and subjects S, and is aligned with the variables generated
#' by the current Stan model.
#'
#' @param fit A fitted CmdStanMCMC object.
#' @param stan_data Optional data list used to fit the model. Supplying it is
#'   strongly recommended because it provides arm membership, time values,
#'   direction of improvement and clinical thresholds.
#' @param meaningful_change Optional fixed MCID fallback for backwards
#'   compatibility. In the current model the posterior draws of `mcid` are
#'   used whenever available.
#' @param y Optional observed outcome vector. If NULL, it is recovered from
#'   `stan_data`.
#' @param credible_level Width of posterior credible intervals.
#' @param responder_thresholds Posterior probability thresholds used for
#'   existing-subject responder classification.
#'
#' @return A list containing population trajectories, changes, new-subject
#'   responder estimands, individual summaries, treatment contrasts,
#'   time-specific gender and age-threshold contrasts (with group sizes when
#'   `stan_data` is supplied), heterogeneity,
#'   posterior predictive checks, log-likelihood information,
#'   MCMC diagnostics, model information and raw posterior draws.
#'
#' @export
mira_summary <- function(
    fit,
    stan_data = NULL,
    meaningful_change = NULL,
    y = NULL,
    credible_level = 0.90,
    responder_thresholds = c(0.50, 0.80, 0.95)
) {

  # ============================================================
  # VALIDATION
  # ============================================================

  if (!inherits(fit, "CmdStanMCMC")) {
    stop("`fit` must be a fitted CmdStanMCMC object.", call. = FALSE)
  }

  if (length(credible_level) != 1 ||
      !is.numeric(credible_level) ||
      !is.finite(credible_level) ||
      credible_level <= 0 ||
      credible_level >= 1) {
    stop("`credible_level` must be between 0 and 1.", call. = FALSE)
  }

  if (length(responder_thresholds) == 0 ||
      !is.numeric(responder_thresholds) ||
      any(!is.finite(responder_thresholds)) ||
      any(responder_thresholds <= 0 | responder_thresholds >= 1)) {
    stop("`responder_thresholds` must contain values between 0 and 1.", call. = FALSE)
  }

  if (!is.null(stan_data) && !is.list(stan_data)) {
    stop("`stan_data` must be NULL or a list.", call. = FALSE)
  }

  if (is.null(y) && !is.null(stan_data) && "y" %in% names(stan_data)) {
    y <- stan_data$y
  }

  if (!is.null(y) &&
      (!is.numeric(y) || length(y) == 0 || any(!is.finite(y)))) {
    stop("`y` must be a non-empty numeric vector of finite values.", call. = FALSE)
  }

  # ============================================================
  # POSTERIOR DRAWS
  # ============================================================

  draws_raw <- tryCatch(
    fit$draws(format = "draws_matrix"),
    error = function(e) {
      stop(
        "Could not extract posterior draws. Original error: ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )

  draws <- as.data.frame(draws_raw)

  if (nrow(draws) < 1) {
    stop("No posterior draws were found.", call. = FALSE)
  }

  # ============================================================
  # HELPERS
  # ============================================================

  alpha <- (1 - credible_level) / 2

  posterior_stats <- function(x) {
    x <- as.numeric(x)
    q <- stats::quantile(
      x,
      probs = c(alpha, 1 - alpha),
      names = FALSE,
      na.rm = TRUE
    )

    c(
      mean = mean(x, na.rm = TRUE),
      median = stats::median(x, na.rm = TRUE),
      sd = stats::sd(x, na.rm = TRUE),
      mad = stats::mad(x, na.rm = TRUE),
      lower = q[1],
      upper = q[2],
      CrI_width = q[2] - q[1]
    )
  }

  summarize_vector <- function(x, parameter) {
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

  get_cols <- function(prefix) {
    grep(
      paste0("^", prefix, "\\["),
      names(draws),
      value = TRUE
    )
  }

  parse_indices <- function(columns) {
    if (length(columns) == 0) {
      return(list())
    }

    lapply(columns, function(x) {
      inside <- sub("^.*\\[", "", x)
      inside <- sub("\\]$", "", inside)
      as.integer(strsplit(inside, ",", fixed = TRUE)[[1]])
    })
  }

  summarize_indexed <- function(prefix, index_names) {
    columns <- get_cols(prefix)

    if (length(columns) == 0) {
      out <- data.frame(stringsAsFactors = FALSE)
      for (nm in index_names) out[[nm]] <- integer(0)
      out$mean <- numeric(0)
      out$median <- numeric(0)
      out$sd <- numeric(0)
      out$mad <- numeric(0)
      out$lower <- numeric(0)
      out$upper <- numeric(0)
      out$CrI_width <- numeric(0)
      return(out)
    }

    idx <- parse_indices(columns)

    if (any(vapply(idx, length, integer(1)) != length(index_names))) {
      stop(
        "Unexpected index structure for posterior variable `",
        prefix,
        "`.",
        call. = FALSE
      )
    }

    stats_list <- lapply(columns, function(col) posterior_stats(draws[[col]]))
    stats_matrix <- do.call(rbind, stats_list)

    index_matrix <- do.call(rbind, idx)
    colnames(index_matrix) <- index_names

    out <- data.frame(
      index_matrix,
      mean = stats_matrix[, "mean"],
      median = stats_matrix[, "median"],
      sd = stats_matrix[, "sd"],
      mad = stats_matrix[, "mad"],
      lower = stats_matrix[, "lower"],
      upper = stats_matrix[, "upper"],
      CrI_width = stats_matrix[, "CrI_width"],
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

    out
  }

  probability_indexed <- function(prefix, index_names, probability_name = "probability") {
    out <- summarize_indexed(prefix, index_names)
    if (nrow(out) > 0) {
      out[[probability_name]] <- out$mean
    } else {
      out[[probability_name]] <- numeric(0)
    }
    out
  }

  first_existing <- function(x, fallback = NULL) {
    if (length(x) > 0) x[[1]] else fallback
  }

  # ============================================================
  # MODEL DIMENSIONS
  # ============================================================

  population_mean_cols <- get_cols("population_mean")
  population_mean_idx <- parse_indices(population_mean_cols)

  if (!is.null(stan_data) && all(c("G", "K") %in% names(stan_data))) {
    G <- as.integer(stan_data$G)
    K <- as.integer(stan_data$K)
  } else {
    if (length(population_mean_idx) == 0) {
      stop(
        "Could not determine G and K. Supply `stan_data` or ensure ",
        "`population_mean[g,k]` is present in the fitted model.",
        call. = FALSE
      )
    }
    G <- max(vapply(population_mean_idx, function(x) x[1], integer(1)))
    K <- max(vapply(population_mean_idx, function(x) x[2], integer(1)))
  }

  individual_change_cols <- get_cols("individual_change_from_baseline")
  individual_change_idx <- parse_indices(individual_change_cols)

  if (!is.null(stan_data) && "S" %in% names(stan_data)) {
    S <- as.integer(stan_data$S)
  } else if (length(individual_change_idx) > 0) {
    S <- max(vapply(individual_change_idx, function(x) x[2], integer(1)))
  } else {
    stop(
      "Could not determine the number of subjects. Supply `stan_data` or ",
      "ensure `individual_change_from_baseline[k,s]` is present.",
      call. = FALSE
    )
  }

  if (!is.null(stan_data) && "N" %in% names(stan_data)) {
    N <- as.integer(stan_data$N)
  } else {
    N <- length(get_cols("log_lik"))
    if (N < 1 && !is.null(y)) N <- length(y)
  }

  # ============================================================
  # METADATA
  # ============================================================

  time_value <- if (!is.null(stan_data) && "time_value" %in% names(stan_data)) {
    as.numeric(stan_data$time_value)
  } else {
    seq_len(K) - 1
  }

  if (length(time_value) != K) {
    stop("`time_value` must contain exactly K values.", call. = FALSE)
  }

  direction <- if (!is.null(stan_data) && "direction" %in% names(stan_data)) {
    as.integer(stan_data$direction)
  } else {
    NA_integer_
  }

  arm_labels <- if (!is.null(stan_data) &&
                    "arm_labels" %in% names(stan_data) &&
                    length(stan_data$arm_labels) == G) {
    as.character(stan_data$arm_labels)
  } else {
    paste0("arm_", seq_len(G))
  }

  subject_arm <- if (!is.null(stan_data) &&
                     "arm" %in% names(stan_data) &&
                     length(stan_data$arm) == S) {
    as.integer(stan_data$arm)
  } else {
    rep(NA_integer_, S)
  }

  subject_male <- if (!is.null(stan_data) &&
                      "male" %in% names(stan_data) &&
                      length(stan_data$male) == S) {
    as.integer(stan_data$male)
  } else {
    rep(NA_integer_, S)
  }

  subject_age_above_threshold <- if (!is.null(stan_data) &&
                                     "age_above_threshold" %in% names(stan_data) &&
                                     length(stan_data$age_above_threshold) == S) {
    as.integer(stan_data$age_above_threshold)
  } else {
    rep(NA_integer_, S)
  }

  gender_labels <- if (!is.null(stan_data) &&
                       "gender_labels" %in% names(stan_data) &&
                       length(stan_data$gender_labels) >= 2) {
    c(
      reference = as.character(stan_data$gender_labels[[1]]),
      comparison = as.character(stan_data$gender_labels[[2]])
    )
  } else {
    c(reference = "Female", comparison = "Male")
  }

  age_threshold <- if (!is.null(stan_data) &&
                       "age_threshold" %in% names(stan_data) &&
                       length(stan_data$age_threshold) >= 1) {
    as.numeric(stan_data$age_threshold)[1]
  } else {
    NA_real_
  }

  age_group_labels <- if (!is.null(stan_data) &&
                          "age_group_labels" %in% names(stan_data) &&
                          length(stan_data$age_group_labels) >= 2) {
    c(
      reference = as.character(stan_data$age_group_labels[[1]]),
      comparison = as.character(stan_data$age_group_labels[[2]])
    )
  } else if (is.finite(age_threshold)) {
    c(
      reference = paste0("<=", age_threshold),
      comparison = paste0(">", age_threshold)
    )
  } else {
    c(reference = "age_reference", comparison = "age_above_threshold")
  }

  between_arm_threshold <- if (!is.null(stan_data) &&
                               "meaningful_between_arm_difference" %in% names(stan_data)) {
    as.numeric(stan_data$meaningful_between_arm_difference)[1]
  } else {
    NA_real_
  }

  # Current model: MCID is a parameter. Fixed meaningful_change is fallback only.
  if ("mcid" %in% names(draws)) {
    mcid_draws <- as.numeric(draws$mcid)
  } else if (!is.null(meaningful_change)) {
    if (length(meaningful_change) != 1 ||
        !is.numeric(meaningful_change) ||
        !is.finite(meaningful_change)) {
      stop("`meaningful_change` must be one finite numeric value.", call. = FALSE)
    }
    mcid_draws <- rep(as.numeric(meaningful_change), nrow(draws))
  } else {
    stop(
      "The fitted model does not contain posterior variable `mcid`. ",
      "For legacy fits, supply `meaningful_change`.",
      call. = FALSE
    )
  }

  mcid_summary <- summarize_vector(mcid_draws, "mcid")

  # ============================================================
  # SCALAR POPULATION / VARIANCE PARAMETERS
  # ============================================================

  scalar_parameters <- c(
    "baseline_mean",
    "beta_time",
    "tau_common",
    "arm_baseline_sd",
    "gender_baseline_effect",
    "beta_gender_time",
    "tau_gender",
    "age_baseline_effect",
    "beta_age_time",
    "tau_age",
    "sigma",
    "nu",
    "sigma_intercept",
    "sigma_slope",
    "rho_subject",
    "residual_sd",
    "mcid"
  )

  scalar_parameters <- scalar_parameters[scalar_parameters %in% names(draws)]

  population_summary <- do.call(
    rbind,
    lapply(scalar_parameters, function(parameter) {
      summarize_vector(draws[[parameter]], parameter)
    })
  )

  if (is.null(population_summary)) {
    population_summary <- data.frame()
  }

  beta_treatment_summary <- summarize_indexed("beta_treatment", "contrast")
  tau_treatment_summary <- summarize_indexed("tau_treatment", "contrast")
  arm_baseline_offset_summary <- summarize_indexed("arm_baseline_offset", "contrast")

  for (obj_name in c("beta_treatment_summary", "tau_treatment_summary", "arm_baseline_offset_summary")) {
    obj <- get(obj_name)
    if (nrow(obj) > 0) {
      obj$treatment_arm <- obj$contrast + 1L
      obj$treatment_label <- arm_labels[obj$treatment_arm]
      obj$reference_arm <- 1L
      obj$reference_label <- arm_labels[1]
      assign(obj_name, obj)
    }
  }

  # ============================================================
  # POPULATION TRAJECTORIES BY ARM
  # ============================================================

  population_time_means <- summarize_indexed("population_mean", c("arm", "time"))

  if (nrow(population_time_means) > 0) {
    population_time_means$arm_label <- arm_labels[population_time_means$arm]
    population_time_means$time_value <- time_value[population_time_means$time]
    population_time_means <- population_time_means[order(
      population_time_means$arm,
      population_time_means$time
    ), ]
  }

  population_change <- summarize_indexed(
    "population_change_from_baseline",
    c("arm", "time")
  )

  directional_change_cols <- get_cols("directional_population_change")
  directional_change_idx <- parse_indices(directional_change_cols)

  if (nrow(population_change) > 0) {
    population_change$arm_label <- arm_labels[population_change$arm]
    population_change$time_value <- time_value[population_change$time]
    population_change$P_positive_raw <- NA_real_
    population_change$P_negative_raw <- NA_real_
    population_change$P_improvement <- NA_real_
    population_change$P_responder <- NA_real_
    population_change$mean_directional_change_minus_mcid <- NA_real_

    raw_change_cols <- get_cols("population_change_from_baseline")
    raw_change_idx <- parse_indices(raw_change_cols)

    for (i in seq_len(nrow(population_change))) {
      g <- population_change$arm[i]
      k <- population_change$time[i]

      raw_pos <- which(vapply(raw_change_idx, function(x) all(x == c(g, k)), logical(1)))
      dir_pos <- which(vapply(directional_change_idx, function(x) all(x == c(g, k)), logical(1)))

      if (length(raw_pos) == 1) {
        x <- draws[[raw_change_cols[raw_pos]]]
        population_change$P_positive_raw[i] <- mean(x > 0)
        population_change$P_negative_raw[i] <- mean(x < 0)
      }

      if (length(dir_pos) == 1) {
        d <- as.numeric(draws[[directional_change_cols[dir_pos]]])
        population_change$P_improvement[i] <- mean(d > 0)
        population_change$P_responder[i] <- mean(d >= mcid_draws)
        population_change$mean_directional_change_minus_mcid[i] <- mean(d - mcid_draws)
      }
    }

    population_change <- population_change[order(
      population_change$arm,
      population_change$time
    ), ]
  }

  population_standardized_change <- summarize_indexed(
    "standardized_population_change",
    c("arm", "time")
  )

  if (nrow(population_standardized_change) > 0) {
    population_standardized_change$arm_label <-
      arm_labels[population_standardized_change$arm]
    population_standardized_change$time_value <-
      time_value[population_standardized_change$time]
  }

  directional_population_change <- summarize_indexed(
    "directional_population_change",
    c("arm", "time")
  )

  if (nrow(directional_population_change) > 0) {
    directional_population_change$arm_label <-
      arm_labels[directional_population_change$arm]
    directional_population_change$time_value <-
      time_value[directional_population_change$time]
  }

  # ============================================================
  # GENDER AND AGE-THRESHOLD DIFFERENCES BY TIME
  # ============================================================

  add_time_metadata <- function(x) {
    if (nrow(x) > 0) {
      x$time_value <- time_value[x$time]
      x <- x[order(x$time), ]
    }
    x
  }

  add_time_probability <- function(x, prefix, output_name, predicate) {
    if (nrow(x) == 0) {
      x[[output_name]] <- numeric(0)
      return(x)
    }

    columns <- get_cols(prefix)
    idx <- parse_indices(columns)

    x[[output_name]] <- vapply(
      x$time,
      function(k) {
        pos <- which(vapply(
          idx,
          function(z) length(z) == 1 && z[1] == k,
          logical(1)
        ))

        if (length(pos) != 1) return(NA_real_)
        mean(predicate(as.numeric(draws[[columns[pos]]])), na.rm = TRUE)
      },
      numeric(1)
    )

    x
  }

  gender_level_difference <- summarize_indexed(
    "male_vs_female_difference",
    "time"
  )
  gender_change_difference <- summarize_indexed(
    "male_vs_female_change_difference",
    "time"
  )
  gender_directional_change_difference <- summarize_indexed(
    "directional_male_vs_female_change_difference",
    "time"
  )

  gender_level_difference <- add_time_metadata(gender_level_difference)
  gender_change_difference <- add_time_metadata(gender_change_difference)
  gender_directional_change_difference <- add_time_metadata(
    gender_directional_change_difference
  )

  if (nrow(gender_level_difference) > 0) {
    gender_level_difference$comparison_group <- gender_labels[["comparison"]]
    gender_level_difference$reference_group <- gender_labels[["reference"]]
    gender_level_difference <- add_time_probability(
      gender_level_difference,
      "male_vs_female_difference",
      "P_comparison_minus_reference_gt_0",
      function(z) z > 0
    )
  }

  if (nrow(gender_change_difference) > 0) {
    gender_change_difference$comparison_group <- gender_labels[["comparison"]]
    gender_change_difference$reference_group <- gender_labels[["reference"]]
    gender_change_difference <- add_time_probability(
      gender_change_difference,
      "male_vs_female_change_difference",
      "P_change_difference_gt_0",
      function(z) z > 0
    )
  }

  if (nrow(gender_directional_change_difference) > 0) {
    gender_directional_change_difference$comparison_group <-
      gender_labels[["comparison"]]
    gender_directional_change_difference$reference_group <-
      gender_labels[["reference"]]
    gender_directional_change_difference <- add_time_probability(
      gender_directional_change_difference,
      "directional_male_vs_female_change_difference",
      "P_comparison_has_better_change",
      function(z) z > 0
    )
  }

  age_level_difference <- summarize_indexed(
    "older_vs_younger_difference",
    "time"
  )
  age_change_difference <- summarize_indexed(
    "older_vs_younger_change_difference",
    "time"
  )
  age_directional_change_difference <- summarize_indexed(
    "directional_older_vs_younger_change_difference",
    "time"
  )

  age_level_difference <- add_time_metadata(age_level_difference)
  age_change_difference <- add_time_metadata(age_change_difference)
  age_directional_change_difference <- add_time_metadata(
    age_directional_change_difference
  )

  if (nrow(age_level_difference) > 0) {
    age_level_difference$comparison_group <- age_group_labels[["comparison"]]
    age_level_difference$reference_group <- age_group_labels[["reference"]]
    age_level_difference$age_threshold <- age_threshold
    age_level_difference <- add_time_probability(
      age_level_difference,
      "older_vs_younger_difference",
      "P_comparison_minus_reference_gt_0",
      function(z) z > 0
    )
  }

  if (nrow(age_change_difference) > 0) {
    age_change_difference$comparison_group <- age_group_labels[["comparison"]]
    age_change_difference$reference_group <- age_group_labels[["reference"]]
    age_change_difference$age_threshold <- age_threshold
    age_change_difference <- add_time_probability(
      age_change_difference,
      "older_vs_younger_change_difference",
      "P_change_difference_gt_0",
      function(z) z > 0
    )
  }

  if (nrow(age_directional_change_difference) > 0) {
    age_directional_change_difference$comparison_group <-
      age_group_labels[["comparison"]]
    age_directional_change_difference$reference_group <-
      age_group_labels[["reference"]]
    age_directional_change_difference$age_threshold <- age_threshold
    age_directional_change_difference <- add_time_probability(
      age_directional_change_difference,
      "directional_older_vs_younger_change_difference",
      "P_comparison_has_better_change",
      function(z) z > 0
    )
  }

  covariate_effects <- list(
    gender = list(
      level_difference = gender_level_difference,
      change_difference = gender_change_difference,
      directional_change_difference = gender_directional_change_difference
    ),
    age_threshold = list(
      threshold = age_threshold,
      level_difference = age_level_difference,
      change_difference = age_change_difference,
      directional_change_difference = age_directional_change_difference
    )
  )

  # ============================================================
  # NEW-SUBJECT RESPONDER ESTIMANDS
  # ============================================================

  latent_any_improvement <- summarize_indexed(
    "latent_new_subject_any_improvement_prob",
    c("arm", "time")
  )

  latent_responder <- summarize_indexed(
    "latent_new_subject_responder_prob",
    c("arm", "time")
  )

  new_subject_latent_responder <- probability_indexed(
    "new_subject_latent_responder_draw",
    c("arm", "time"),
    "posterior_predictive_probability"
  )

  new_subject_predictive_responder <- probability_indexed(
    "new_subject_predictive_responder_draw",
    c("arm", "time"),
    "posterior_predictive_probability"
  )

  new_subject_latent_change <- summarize_indexed(
    "new_subject_latent_change_draw",
    c("arm", "time")
  )

  new_subject_predictive_change <- summarize_indexed(
    "new_subject_predictive_change_draw",
    c("arm", "time")
  )

  add_arm_time_metadata <- function(x) {
    if (nrow(x) > 0) {
      x$arm_label <- arm_labels[x$arm]
      x$time_value <- time_value[x$time]
      x$gender_profile <- gender_labels[["reference"]]
      x$age_group_profile <- age_group_labels[["reference"]]
      x <- x[order(x$arm, x$time), ]
    }
    x
  }

  latent_any_improvement <- add_arm_time_metadata(latent_any_improvement)
  latent_responder <- add_arm_time_metadata(latent_responder)
  new_subject_latent_responder <- add_arm_time_metadata(new_subject_latent_responder)
  new_subject_predictive_responder <- add_arm_time_metadata(new_subject_predictive_responder)
  new_subject_latent_change <- add_arm_time_metadata(new_subject_latent_change)
  new_subject_predictive_change <- add_arm_time_metadata(new_subject_predictive_change)

  # ============================================================
  # INDIVIDUAL CHANGE AND CLINICAL RESPONSE
  # ============================================================

  individual_change <- summarize_indexed(
    "individual_change_from_baseline",
    c("time", "subject")
  )

  individual_directional_change <- summarize_indexed(
    "individual_directional_change",
    c("time", "subject")
  )

  individual_any <- probability_indexed(
    "individual_any_improvement_draw",
    c("time", "subject"),
    "P_improvement"
  )

  individual_responder <- probability_indexed(
    "individual_meaningful_responder_draw",
    c("time", "subject"),
    "P_MCID"
  )

  individual_distance <- summarize_indexed(
    "individual_change_minus_mcid",
    c("time", "subject")
  )

  add_subject_metadata <- function(x) {
    if (nrow(x) > 0) {
      x$time_value <- time_value[x$time]
      x$arm <- subject_arm[x$subject]
      x$arm_label <- ifelse(
        is.na(x$arm),
        NA_character_,
        arm_labels[x$arm]
      )
      x$male <- subject_male[x$subject]
      x$gender_label <- ifelse(
        is.na(x$male),
        NA_character_,
        ifelse(
          x$male == 1L,
          gender_labels[["comparison"]],
          gender_labels[["reference"]]
        )
      )
      x$age_above_threshold <- subject_age_above_threshold[x$subject]
      x$age_group_label <- ifelse(
        is.na(x$age_above_threshold),
        NA_character_,
        ifelse(
          x$age_above_threshold == 1L,
          age_group_labels[["comparison"]],
          age_group_labels[["reference"]]
        )
      )
      x <- x[order(x$subject, x$time), ]
    }
    x
  }

  individual_change <- add_subject_metadata(individual_change)
  individual_directional_change <- add_subject_metadata(individual_directional_change)
  individual_any <- add_subject_metadata(individual_any)
  individual_responder <- add_subject_metadata(individual_responder)
  individual_distance <- add_subject_metadata(individual_distance)

  individual_clinical <- data.frame()

  if (nrow(individual_responder) > 0) {
    individual_clinical <- individual_responder[
      , c(
        "time", "subject", "P_MCID", "time_value", "arm", "arm_label",
        "male", "gender_label", "age_above_threshold", "age_group_label"
      ),
      drop = FALSE
    ]

    if (nrow(individual_any) > 0) {
      key <- paste(individual_any$time, individual_any$subject, sep = ":")
      target <- paste(individual_clinical$time, individual_clinical$subject, sep = ":")
      m <- match(target, key)
      individual_clinical$P_improvement <- individual_any$P_improvement[m]
    } else {
      individual_clinical$P_improvement <- NA_real_
    }

    if (nrow(individual_distance) > 0) {
      key <- paste(individual_distance$time, individual_distance$subject, sep = ":")
      target <- paste(individual_clinical$time, individual_clinical$subject, sep = ":")
      m <- match(target, key)
      individual_clinical$mean_change_minus_MCID <- individual_distance$mean[m]
      individual_clinical$P_change_minus_MCID_ge_0 <- vapply(
        seq_len(nrow(individual_clinical)),
        function(i) {
          col_name <- paste0(
            "individual_change_minus_mcid[",
            individual_clinical$time[i],
            ",",
            individual_clinical$subject[i],
            "]"
          )
          if (col_name %in% names(draws)) mean(draws[[col_name]] >= 0) else NA_real_
        },
        numeric(1)
      )
    } else {
      individual_clinical$mean_change_minus_MCID <- NA_real_
      individual_clinical$P_change_minus_MCID_ge_0 <- NA_real_
    }

    individual_clinical$responder_class <- cut(
      individual_clinical$P_MCID,
      breaks = c(-Inf, 0.20, 0.50, 0.80, 0.95, Inf),
      labels = c(
        "very_unlikely",
        "uncertain",
        "probable",
        "high_probability",
        "very_high_probability"
      ),
      right = FALSE
    )

    for (threshold in responder_thresholds) {
      nm <- paste0("response_", sprintf("%02d", round(100 * threshold)))
      individual_clinical[[nm]] <- individual_clinical$P_MCID >= threshold
    }
  }

  # Existing-subject responder overview by arm and time.
  responder_summary <- data.frame()

  if (nrow(individual_clinical) > 0) {
    combinations <- unique(individual_clinical[, c("arm", "arm_label", "time", "time_value")])
    combinations <- combinations[order(combinations$arm, combinations$time), ]

    rows <- lapply(seq_len(nrow(combinations)), function(i) {
      ix <- individual_clinical$arm == combinations$arm[i] &
        individual_clinical$time == combinations$time[i]
      p <- individual_clinical$P_MCID[ix]

      out <- data.frame(
        arm = combinations$arm[i],
        arm_label = combinations$arm_label[i],
        time = combinations$time[i],
        time_value = combinations$time_value[i],
        expected_responder_proportion = mean(p),
        stringsAsFactors = FALSE
      )

      for (threshold in responder_thresholds) {
        nm <- paste0("proportion_P_MCID_ge_", sprintf("%02d", round(100 * threshold)))
        out[[nm]] <- mean(p >= threshold)
      }

      out
    })

    responder_summary <- do.call(rbind, rows)
  }

  # ============================================================
  # TREATMENT EFFECTS VS REFERENCE ARM
  # ============================================================

  treatment_effects <- summarize_indexed(
    "treatment_change_difference",
    c("contrast", "time")
  )

  directional_treatment <- summarize_indexed(
    "directional_treatment_benefit",
    c("contrast", "time")
  )

  positive_treatment <- probability_indexed(
    "treatment_benefit_positive_draw",
    c("contrast", "time"),
    "P_benefit_positive"
  )

  meaningful_treatment <- probability_indexed(
    "treatment_benefit_meaningful_draw",
    c("contrast", "time"),
    "P_benefit_meaningful"
  )

  responder_uplift <- summarize_indexed(
    "latent_responder_probability_difference",
    c("contrast", "time")
  )

  if (nrow(treatment_effects) > 0) {
    treatment_effects$treatment_arm <- treatment_effects$contrast + 1L
    treatment_effects$treatment_label <- arm_labels[treatment_effects$treatment_arm]
    treatment_effects$reference_arm <- 1L
    treatment_effects$reference_label <- arm_labels[1]
    treatment_effects$time_value <- time_value[treatment_effects$time]

    treatment_key <- paste(treatment_effects$contrast, treatment_effects$time, sep = ":")

    add_from_table <- function(base, table, value_column, output_name) {
      if (nrow(table) == 0) {
        base[[output_name]] <- NA_real_
        return(base)
      }
      key <- paste(table$contrast, table$time, sep = ":")
      base[[output_name]] <- table[[value_column]][match(treatment_key, key)]
      base
    }

    treatment_effects <- add_from_table(
      treatment_effects, directional_treatment, "mean", "mean_directional_benefit"
    )
    treatment_effects <- add_from_table(
      treatment_effects, positive_treatment, "P_benefit_positive", "P_benefit_positive"
    )
    treatment_effects <- add_from_table(
      treatment_effects, meaningful_treatment, "P_benefit_meaningful", "P_benefit_meaningful"
    )
    treatment_effects <- add_from_table(
      treatment_effects, responder_uplift, "mean", "mean_responder_probability_difference"
    )

    treatment_effects <- treatment_effects[order(
      treatment_effects$contrast,
      treatment_effects$time
    ), ]
  }

  treatment_responder_uplift <- responder_uplift
  if (nrow(treatment_responder_uplift) > 0) {
    treatment_responder_uplift$treatment_arm <- treatment_responder_uplift$contrast + 1L
    treatment_responder_uplift$treatment_label <-
      arm_labels[treatment_responder_uplift$treatment_arm]
    treatment_responder_uplift$reference_label <- arm_labels[1]
    treatment_responder_uplift$time_value <-
      time_value[treatment_responder_uplift$time]
  }

  # ============================================================
  # HETEROGENEITY
  # ============================================================

  heterogeneity_parameters <- c(
    "sigma_intercept",
    "sigma_slope",
    "rho_subject",
    "arm_baseline_sd",
    "tau_common",
    "tau_gender",
    "tau_age"
  )

  heterogeneity_parameters <- heterogeneity_parameters[
    heterogeneity_parameters %in% names(draws)
  ]

  heterogeneity <- do.call(
    rbind,
    lapply(heterogeneity_parameters, function(parameter) {
      summarize_vector(draws[[parameter]], parameter)
    })
  )

  if (is.null(heterogeneity)) heterogeneity <- data.frame()

  if (nrow(heterogeneity) > 0) {
    heterogeneity$P_positive <- NA_real_
    heterogeneity$P_negative <- NA_real_

    for (i in seq_len(nrow(heterogeneity))) {
      parameter <- heterogeneity$parameter[i]
      heterogeneity$P_positive[i] <- mean(draws[[parameter]] > 0)
      heterogeneity$P_negative[i] <- mean(draws[[parameter]] < 0)
    }
  }

  # ============================================================
  # CLINICAL FINAL-TIME SUMMARY
  # ============================================================

  clinical_summary <- data.frame()

  if (nrow(population_change) > 0) {
    final_rows <- population_change$time == K
    final_change <- population_change[final_rows, , drop = FALSE]

    clinical_summary <- data.frame(
      arm = final_change$arm,
      arm_label = final_change$arm_label,
      final_time = K,
      final_time_value = time_value[K],
      gender_profile = gender_labels[["reference"]],
      age_group_profile = age_group_labels[["reference"]],
      final_population_change_mean = final_change$mean,
      P_final_improvement = final_change$P_improvement,
      P_final_responder = final_change$P_responder,
      mean_final_directional_change_minus_mcid =
        final_change$mean_directional_change_minus_mcid,
      stringsAsFactors = FALSE
    )
  }

  # ============================================================
  # POSTERIOR PREDICTIVE CHECKS
  # ============================================================

  y_rep_cols <- get_cols("y_rep")
  ppc <- NULL

  if (length(y_rep_cols) > 0 && !is.null(y)) {
    y_rep_matrix <- as.matrix(draws[, y_rep_cols, drop = FALSE])

    predictive_mean <- rowMeans(y_rep_matrix)
    predictive_sd <- apply(y_rep_matrix, 1, stats::sd)
    predictive_median <- apply(y_rep_matrix, 1, stats::median)
    predictive_q05 <- apply(y_rep_matrix, 1, stats::quantile, probs = 0.05)
    predictive_q95 <- apply(y_rep_matrix, 1, stats::quantile, probs = 0.95)

    observed <- c(
      mean = mean(y),
      sd = stats::sd(y),
      median = stats::median(y),
      q05 = as.numeric(stats::quantile(y, 0.05)),
      q95 = as.numeric(stats::quantile(y, 0.95))
    )

    ppc <- list(
      observed = observed,
      posterior_predictive = list(
        mean = posterior_stats(predictive_mean),
        sd = posterior_stats(predictive_sd),
        median = posterior_stats(predictive_median),
        q05 = posterior_stats(predictive_q05),
        q95 = posterior_stats(predictive_q95)
      ),
      bayesian_p_values = c(
        mean = mean(predictive_mean >= mean(y)),
        sd = mean(predictive_sd >= stats::sd(y)),
        median = mean(predictive_median >= stats::median(y))
      ),
      predictive_draws = y_rep_matrix
    )
  }

  # ============================================================
  # LOG-LIKELIHOOD INFORMATION
  # ============================================================

  log_lik_cols <- get_cols("log_lik")
  loo_information <- NULL

  if (length(log_lik_cols) > 0) {
    log_lik_matrix <- as.matrix(draws[, log_lik_cols, drop = FALSE])

    loo_information <- list(
      n_observations = length(log_lik_cols),
      mean_total_log_lik = mean(rowSums(log_lik_matrix)),
      pointwise_mean_log_lik = colMeans(log_lik_matrix),
      draws = log_lik_matrix
    )
  }

  # ============================================================
  # MCMC DIAGNOSTICS
  # ============================================================

  diagnostics <- tryCatch({
    s <- fit$summary()

    d <- data.frame(
      parameter = s$variable,
      Rhat = s$rhat,
      ESS_bulk = s$ess_bulk,
      ESS_tail = s$ess_tail,
      stringsAsFactors = FALSE
    )

    d$Rhat_ok <- is.na(d$Rhat) | d$Rhat < 1.01
    d$ESS_bulk_ok <- is.na(d$ESS_bulk) | d$ESS_bulk >= 400
    d$ESS_tail_ok <- is.na(d$ESS_tail) | d$ESS_tail >= 400

    finite_rhat <- d$Rhat[is.finite(d$Rhat)]
    finite_bulk <- d$ESS_bulk[is.finite(d$ESS_bulk)]
    finite_tail <- d$ESS_tail[is.finite(d$ESS_tail)]

    list(
      parameters = d,
      max_Rhat = if (length(finite_rhat)) max(finite_rhat) else NA_real_,
      min_ESS_bulk = if (length(finite_bulk)) min(finite_bulk) else NA_real_,
      min_ESS_tail = if (length(finite_tail)) min(finite_tail) else NA_real_
    )
  }, error = function(e) {
    list(error = conditionMessage(e))
  })

  quality_flags <- list(
    Rhat_ok = if (!is.null(diagnostics$max_Rhat) && is.finite(diagnostics$max_Rhat)) {
      diagnostics$max_Rhat < 1.01
    } else {
      NA
    },
    ESS_bulk_ok = if (!is.null(diagnostics$min_ESS_bulk) && is.finite(diagnostics$min_ESS_bulk)) {
      diagnostics$min_ESS_bulk >= 400
    } else {
      NA
    },
    ESS_tail_ok = if (!is.null(diagnostics$min_ESS_tail) && is.finite(diagnostics$min_ESS_tail)) {
      diagnostics$min_ESS_tail >= 400
    } else {
      NA
    }
  )

  # ============================================================
  # MODEL INFORMATION
  # ============================================================

  gender_counts <- if (all(!is.na(subject_male))) {
    stats::setNames(
      c(sum(subject_male == 0L), sum(subject_male == 1L)),
      c(gender_labels[["reference"]], gender_labels[["comparison"]])
    )
  } else {
    stats::setNames(
      c(NA_integer_, NA_integer_),
      c(gender_labels[["reference"]], gender_labels[["comparison"]])
    )
  }

  age_group_counts <- if (all(!is.na(subject_age_above_threshold))) {
    stats::setNames(
      c(
        sum(subject_age_above_threshold == 0L),
        sum(subject_age_above_threshold == 1L)
      ),
      c(age_group_labels[["reference"]], age_group_labels[["comparison"]])
    )
  } else {
    stats::setNames(
      c(NA_integer_, NA_integer_),
      c(age_group_labels[["reference"]], age_group_labels[["comparison"]])
    )
  }

  model_information <- list(
    n_observations = N,
    n_subjects = S,
    n_time_points = K,
    n_arms = G,
    arm_labels = arm_labels,
    reference_arm = 1L,
    reference_label = arm_labels[1],
    gender_labels = gender_labels,
    gender_reference = gender_labels[["reference"]],
    gender_comparison = gender_labels[["comparison"]],
    gender_counts = gender_counts,
    age_threshold = age_threshold,
    age_group_labels = age_group_labels,
    age_reference = age_group_labels[["reference"]],
    age_comparison = age_group_labels[["comparison"]],
    age_group_counts = age_group_counts,
    population_reference_profile = paste0(
      gender_labels[["reference"]],
      ", age ",
      age_group_labels[["reference"]]
    ),
    time_value = time_value,
    direction = direction,
    direction_interpretation = if (is.na(direction)) {
      NA_character_
    } else if (direction == 1L) {
      "higher outcome = better"
    } else {
      "lower outcome = better"
    },
    mcid = mcid_summary,
    mcid_prior_mean = if (!is.null(stan_data) && "mcid_prior_mean" %in% names(stan_data)) {
      stan_data$mcid_prior_mean
    } else {
      NA_real_
    },
    mcid_prior_sd = if (!is.null(stan_data) && "mcid_prior_sd" %in% names(stan_data)) {
      stan_data$mcid_prior_sd
    } else {
      NA_real_
    },
    meaningful_between_arm_difference = between_arm_threshold,
    covariate_prior_structure = list(
      gender = c(
        "gender_baseline_prior_sd",
        "beta_gender_prior_sd",
        "tau_gender_prior_rate"
      ),
      age_threshold = c(
        "age_baseline_prior_sd",
        "beta_age_prior_sd",
        "tau_age_prior_rate"
      )
    ),
    credible_level = credible_level,
    responder_thresholds = responder_thresholds,
    observed_mean = if (!is.null(y)) mean(y) else NA_real_,
    observed_sd = if (!is.null(y)) stats::sd(y) else NA_real_,
    observed_median = if (!is.null(y)) stats::median(y) else NA_real_
  )

  # ============================================================
  # VARIABLE INVENTORY
  # ============================================================

  variable_inventory <- list(
    population = c(
      get_cols("population_mean"),
      get_cols("population_change_from_baseline"),
      get_cols("directional_population_change"),
      get_cols("standardized_population_change")
    ),
    covariates = c(
      get_cols("gender_effect"),
      get_cols("age_threshold_effect"),
      get_cols("male_vs_female_difference"),
      get_cols("male_vs_female_change_difference"),
      get_cols("directional_male_vs_female_change_difference"),
      get_cols("older_vs_younger_difference"),
      get_cols("older_vs_younger_change_difference"),
      get_cols("directional_older_vs_younger_change_difference")
    ),
    new_subject = c(
      get_cols("latent_new_subject_any_improvement_prob"),
      get_cols("latent_new_subject_responder_prob"),
      get_cols("new_subject_latent_change_draw"),
      get_cols("new_subject_latent_responder_draw"),
      get_cols("new_subject_predictive_change_draw"),
      get_cols("new_subject_predictive_responder_draw")
    ),
    individual = c(
      get_cols("individual_change_from_baseline"),
      get_cols("individual_directional_change"),
      get_cols("individual_any_improvement_draw"),
      get_cols("individual_meaningful_responder_draw"),
      get_cols("individual_change_minus_mcid")
    ),
    treatment = c(
      get_cols("treatment_change_difference"),
      get_cols("directional_treatment_benefit"),
      get_cols("treatment_benefit_positive_draw"),
      get_cols("treatment_benefit_meaningful_draw"),
      get_cols("latent_responder_probability_difference")
    ),
    observation_level = c(
      get_cols("log_lik"),
      get_cols("y_rep")
    )
  )

  # ============================================================
  # RETURN
  # ============================================================

  list(
    # Core population parameters
    population = population_summary,
    beta_treatment = beta_treatment_summary,
    tau_treatment = tau_treatment_summary,
    arm_baseline_offset = arm_baseline_offset_summary,

    # Population trajectories
    population_time_means = population_time_means,
    population_change = population_change,
    population_directional_change = directional_population_change,
    population_standardized_change = population_standardized_change,

    # Time-specific adjusted gender and age-threshold contrasts
    covariate_effects = covariate_effects,
    gender_effects = covariate_effects$gender,
    age_threshold_effects = covariate_effects$age_threshold,

    # Kept for backwards-compatible list structure; these estimands are
    # not explicitly generated by the new Stan model.
    population_consecutive_change = data.frame(),
    population_slope = data.frame(),

    # New-subject clinical estimands
    population_clinical = list(
      latent_any_improvement = latent_any_improvement,
      latent_responder = latent_responder,
      latent_responder_draw = new_subject_latent_responder,
      predictive_responder_draw = new_subject_predictive_responder,
      latent_change_draw = new_subject_latent_change,
      predictive_change_draw = new_subject_predictive_change
    ),

    # Existing subjects
    individual_change = individual_change,
    individual_directional_change = individual_directional_change,
    individual_clinical = individual_clinical,
    individual_slope = data.frame(),

    # Responders
    responder_proportion = latent_responder,
    responders = responder_summary,

    # Treatment effects vs reference arm
    treatment_effects = treatment_effects,
    treatment_responder_uplift = treatment_responder_uplift,

    # Clinical interpretation
    clinical = clinical_summary,
    mcid = mcid_summary,

    # Heterogeneity
    heterogeneity = heterogeneity,

    # Model checking
    ppc = ppc,
    loo = loo_information,
    diagnostics = diagnostics,
    quality_flags = quality_flags,

    # Metadata and debugging
    model_information = model_information,
    variable_inventory = variable_inventory,
    draws = draws
  )
}
