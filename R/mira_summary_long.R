#' Research-oriented summary of a MIRA Bayesian fit
#'
#' Summarizes the outcome-adaptive multi-arm MIRA longitudinal model. The
#' function is dynamic in likelihood/link, treatment arms G, measurement
#' occasions K and subjects S. Natural-scale outcomes and clinical changes are
#' reported for every family; log-scale and ratio estimands are additionally
#' exposed for log-normal CMT fits.
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
#' @param verbose Logical. If TRUE, print the compact Bayesian report when
#'   the summary object is created.
#'
#' @return A list containing population trajectories, changes, new-subject
#'   responder estimands, individual summaries, treatment contrasts,
#'   time-specific effects for every selected encoded covariate term (with
#'   readable design-matrix metadata when available), heterogeneity,
#'   posterior predictive checks, log-likelihood information,
#'   MCMC diagnostics, model information and raw posterior draws.
#'
#' @export
mira_summary_long <- function(
    fit,
    stan_data = NULL,
    meaningful_change = NULL,
    y = NULL,
    credible_level = 0.90,
    responder_thresholds = c(0.50, 0.80, 0.95),
    verbose = TRUE
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

  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("`verbose` must be TRUE or FALSE.", call. = FALSE)
  }

  fit_info <- attr(fit, "mira_fit_info", exact = TRUE)
  if (!is.null(fit_info) && !is.list(fit_info)) {
    fit_info <- NULL
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

  posterior_stats_long <- function(x) {
    x <- as.numeric(x)
    x <- x[is.finite(x)]

    if (length(x) == 0L) {
      return(c(
        mean = NA_real_,
        median = NA_real_,
        sd = NA_real_,
        mad = NA_real_,
        lower = NA_real_,
        upper = NA_real_,
        CrI_width = NA_real_
      ))
    }

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

  summarize_vector_long <- function(x, parameter) {
    s <- posterior_stats_long(x)
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

  get_cols_long <- function(prefix) {
    grep(
      paste0("^", prefix, "\\["),
      names(draws),
      value = TRUE
    )
  }

  parse_indices_long <- function(columns) {
    if (length(columns) == 0) {
      return(list())
    }

    lapply(columns, function(x) {
      inside <- sub("^.*\\[", "", x)
      inside <- sub("\\]$", "", inside)
      as.integer(strsplit(inside, ",", fixed = TRUE)[[1]])
    })
  }

  summarize_indexed_long <- function(prefix, index_names) {
    columns <- get_cols_long(prefix)

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

    idx <- parse_indices_long(columns)

    if (any(vapply(idx, length, integer(1)) != length(index_names))) {
      stop(
        "Unexpected index structure for posterior variable `",
        prefix,
        "`.",
        call. = FALSE
      )
    }

    stats_list <- lapply(columns, function(col) posterior_stats_long(draws[[col]]))
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

  probability_indexed_long <- function(prefix, index_names, probability_name = "probability") {
    out <- summarize_indexed_long(prefix, index_names)
    if (nrow(out) > 0) {
      out[[probability_name]] <- out$mean
    } else {
      out[[probability_name]] <- numeric(0)
    }
    out
  }

  first_existing_long <- function(x, fallback = NULL) {
    if (length(x) > 0) x[[1]] else fallback
  }

  # ============================================================
  # MODEL DIMENSIONS
  # ============================================================

  population_mean_cols <- get_cols_long("population_mean")
  population_mean_idx <- parse_indices_long(population_mean_cols)

  validate_dimension_long <- function(x, name, lower = 1L) {
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
        x < lower || x > .Machine$integer.max || x != as.integer(x)) {
      stop(
        "`stan_data$", name, "` must be one integer >= ", lower, ".",
        call. = FALSE
      )
    }
    as.integer(x)
  }

  if (!is.null(stan_data) && all(c("G", "K") %in% names(stan_data))) {
    G <- validate_dimension_long(stan_data$G, "G", 2L)
    K <- validate_dimension_long(stan_data$K, "K", 2L)
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

  individual_change_cols <- get_cols_long("individual_change_from_baseline")
  individual_change_idx <- parse_indices_long(individual_change_cols)

  if (!is.null(stan_data) && "S" %in% names(stan_data)) {
    S <- validate_dimension_long(stan_data$S, "S")
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
    N <- validate_dimension_long(stan_data$N, "N")
  } else {
    N <- length(get_cols_long("log_lik"))
    if (N < 1 && !is.null(y)) N <- length(y)
  }

  covariate_effect_cols <- get_cols_long("covariate_effect")
  covariate_effect_idx <- parse_indices_long(covariate_effect_cols)

  if (!is.null(stan_data) && "P" %in% names(stan_data)) {
    P <- validate_dimension_long(stan_data$P, "P", 0L)
  } else if (!is.null(fit_info) && !is.null(fit_info$P)) {
    P_value <- fit_info$P
    if (!is.numeric(P_value) || length(P_value) != 1L ||
        !is.finite(P_value) || P_value < 0L ||
        P_value != as.integer(P_value)) {
      stop("Attached fit metadata contain an invalid covariate dimension P.", call. = FALSE)
    }
    P <- as.integer(P_value)
  } else if (length(covariate_effect_idx) > 0L) {
    P <- max(vapply(covariate_effect_idx, function(x) x[[1L]], integer(1L)))
  } else {
    P <- 0L
  }

  # When both sources are available, never let caller-supplied `stan_data`
  # relabel posterior draws from a different same-sized design matrix.
  if (!is.null(stan_data) && !is.null(fit_info) && !is.null(fit_info$P)) {
    fit_P <- fit_info$P
    if (!is.numeric(fit_P) || length(fit_P) != 1L ||
        !is.finite(fit_P) || fit_P < 0L || fit_P != as.integer(fit_P)) {
      stop("Attached fit metadata contain an invalid covariate dimension P.", call. = FALSE)
    }
    if (as.integer(fit_P) != P) {
      stop(
        "`stan_data$P` does not match the covariate dimension stored in `fit`.",
        call. = FALSE
      )
    }

    stan_covariate_names <- if (!is.null(stan_data$covariate_names)) {
      as.character(stan_data$covariate_names)
    } else if (!is.null(stan_data$X) && !is.null(colnames(stan_data$X))) {
      as.character(colnames(stan_data$X))
    } else {
      NULL
    }
    fit_covariate_names <- if (!is.null(fit_info$covariate_names)) {
      as.character(fit_info$covariate_names)
    } else {
      NULL
    }
    if (!is.null(stan_covariate_names) && !is.null(fit_covariate_names) &&
        !identical(stan_covariate_names, fit_covariate_names)) {
      stop(
        "Covariate names/order in `stan_data` do not match the metadata ",
        "stored in `fit`.",
        call. = FALSE
      )
    }
  }

  # ============================================================
  # METADATA
  # ============================================================

  time_value <- if (!is.null(stan_data) && "time_value" %in% names(stan_data)) {
    if (!is.numeric(stan_data$time_value)) {
      stop("`stan_data$time_value` must be numeric.", call. = FALSE)
    }
    as.numeric(stan_data$time_value)
  } else {
    seq_len(K) - 1
  }

  if (length(time_value) != K || any(!is.finite(time_value)) ||
      anyDuplicated(time_value) || is.unsorted(time_value, strictly = TRUE)) {
    stop(
      "`time_value` must contain exactly K finite, strictly increasing values.",
      call. = FALSE
    )
  }

  direction <- if (!is.null(stan_data) && "direction" %in% names(stan_data)) {
    direction_value <- stan_data$direction
    if (!is.numeric(direction_value) || length(direction_value) != 1L ||
        !is.finite(direction_value) || !(direction_value %in% c(-1, 1))) {
      stop("`stan_data$direction` must be exactly +1 or -1.", call. = FALSE)
    }
    as.integer(direction_value)
  } else {
    NA_integer_
  }

  likelihood_id <- if (!is.null(stan_data) &&
                       "likelihood_id" %in% names(stan_data)) {
    as.integer(stan_data$likelihood_id)[1L]
  } else if (!is.null(fit_info) && !is.null(fit_info$likelihood_id)) {
    as.integer(fit_info$likelihood_id)[1L]
  } else {
    # Legacy fits were Student-t identity models.
    1L
  }

  if (!(likelihood_id %in% 1:3)) {
    stop("Could not determine a valid MIRA likelihood id.", call. = FALSE)
  }

  likelihood <- c("student_t", "gaussian", "lognormal")[[likelihood_id]]
  modeling_scale <- if (likelihood_id == 3L) "log" else "identity"
  outcome_name <- if (!is.null(stan_data) && !is.null(stan_data$outcome_name)) {
    as.character(stan_data$outcome_name)[1L]
  } else if (!is.null(fit_info) && !is.null(fit_info$outcome_name)) {
    as.character(fit_info$outcome_name)[1L]
  } else {
    NA_character_
  }
  outcome <- if (!is.null(stan_data) && !is.null(stan_data$outcome)) {
    as.character(stan_data$outcome)[1L]
  } else if (!is.null(fit_info) && !is.null(fit_info$outcome)) {
    as.character(fit_info$outcome)[1L]
  } else {
    outcome_name
  }

  fit_bounds <- if (!is.null(fit_info) &&
                    !is.null(fit_info$outcome_bounds) &&
                    length(fit_info$outcome_bounds) == 2L) {
    as.numeric(fit_info$outcome_bounds)
  } else {
    c(NA_real_, NA_real_)
  }
  has_lower_bound <- if (!is.null(stan_data) &&
                         !is.null(stan_data$has_lower_bound)) {
    identical(as.integer(stan_data$has_lower_bound)[1L], 1L)
  } else {
    is.finite(fit_bounds[[1L]])
  }
  has_upper_bound <- if (!is.null(stan_data) &&
                         !is.null(stan_data$has_upper_bound)) {
    identical(as.integer(stan_data$has_upper_bound)[1L], 1L)
  } else {
    is.finite(fit_bounds[[2L]])
  }
  lower_bound <- if (has_lower_bound) {
    if (!is.null(stan_data) && !is.null(stan_data$outcome_lower_bound)) {
      as.numeric(stan_data$outcome_lower_bound)[1L]
    } else {
      fit_bounds[[1L]]
    }
  } else {
    NA_real_
  }
  upper_bound <- if (has_upper_bound) {
    if (!is.null(stan_data) && !is.null(stan_data$outcome_upper_bound)) {
      as.numeric(stan_data$outcome_upper_bound)[1L]
    } else {
      fit_bounds[[2L]]
    }
  } else {
    NA_real_
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

  subject_labels <- if (!is.null(stan_data) &&
                        "subject_labels" %in% names(stan_data) &&
                        length(stan_data$subject_labels) == S) {
    as.character(stan_data$subject_labels)
  } else {
    as.character(seq_len(S))
  }

  metadata_value_long <- function(name, default = NULL) {
    if (!is.null(stan_data) && !is.null(stan_data[[name]])) {
      return(stan_data[[name]])
    }
    if (!is.null(fit_info) && !is.null(fit_info[[name]])) {
      return(fit_info[[name]])
    }
    default
  }

  covariate_names <- metadata_value_long("covariate_names")
  if (is.null(covariate_names) && !is.null(stan_data$X) &&
      !is.null(colnames(stan_data$X))) {
    covariate_names <- colnames(stan_data$X)
  }
  if (is.null(covariate_names)) {
    covariate_names <- paste0("covariate_", seq_len(P))
  }
  covariate_names <- as.character(covariate_names)
  if (length(covariate_names) != P || anyNA(covariate_names) ||
      any(!nzchar(covariate_names)) || anyDuplicated(covariate_names)) {
    stop(
      "Covariate metadata must provide exactly P unique, non-missing encoded names.",
      call. = FALSE
    )
  }

  covariate_field_long <- function(name, default) {
    value <- metadata_value_long(name)
    if (is.null(value) || length(value) != P) default else value
  }

  covariate_original_names <- as.character(covariate_field_long(
    "covariate_original_names", covariate_names
  ))
  covariate_labels <- as.character(covariate_field_long(
    "covariate_labels", covariate_names
  ))
  covariate_types <- as.character(covariate_field_long(
    "covariate_types", rep("unspecified", P)
  ))
  covariate_reference_levels <- as.character(covariate_field_long(
    "covariate_reference_levels", rep(NA_character_, P)
  ))
  covariate_centers <- as.numeric(covariate_field_long(
    "covariate_centers", rep(0, P)
  ))
  covariate_scales <- as.numeric(covariate_field_long(
    "covariate_scales", rep(1, P)
  ))

  covariate_map <- metadata_value_long("covariate_map")
  covariate_metadata_object <- metadata_value_long("covariate_metadata")
  if ((!is.data.frame(covariate_map) || nrow(covariate_map) != P) &&
      is.list(covariate_metadata_object) &&
      is.data.frame(covariate_metadata_object$columns)) {
    covariate_map <- covariate_metadata_object$columns
  }
  if (!is.data.frame(covariate_map) || nrow(covariate_map) != P) {
    covariate_map <- data.frame(
      index = seq_len(P),
      name = covariate_names,
      label = covariate_labels,
      original_name = covariate_original_names,
      original_label = covariate_original_names,
      type = covariate_types,
      encoding = rep("unspecified", P),
      level = rep(NA_character_, P),
      reference_level = covariate_reference_levels,
      center = covariate_centers,
      scale = covariate_scales,
      unit = rep(NA_character_, P),
      stringsAsFactors = FALSE
    )
  }

  required_map_defaults <- list(
    index = seq_len(P),
    name = covariate_names,
    label = covariate_labels,
    original_name = covariate_original_names,
    original_label = covariate_original_names,
    type = covariate_types,
    encoding = rep("unspecified", P),
    level = rep(NA_character_, P),
    reference_level = covariate_reference_levels,
    center = covariate_centers,
    scale = covariate_scales,
    unit = rep(NA_character_, P)
  )
  for (nm in names(required_map_defaults)) {
    if (!nm %in% names(covariate_map)) {
      covariate_map[[nm]] <- required_map_defaults[[nm]]
    }
  }
  covariate_map$index <- as.integer(covariate_map$index)
  covariate_map$name <- as.character(covariate_map$name)
  if (!identical(covariate_map$index, seq_len(P)) ||
      !identical(covariate_map$name, covariate_names)) {
    stop(
      "`covariate_map` indices and names must match the P columns of X in order.",
      call. = FALSE
    )
  }

  subject_X <- NULL
  if (!is.null(stan_data) && !is.null(stan_data$X)) {
    subject_X <- stan_data$X
    if (!is.matrix(subject_X) || !is.numeric(subject_X) ||
        !identical(as.integer(dim(subject_X)), c(S, P)) ||
        any(!is.finite(subject_X))) {
      stop("`stan_data$X` must be a finite numeric S by P matrix.", call. = FALSE)
    }
    colnames(subject_X) <- covariate_names
  } else if (P == 0L) {
    subject_X <- matrix(numeric(0), nrow = S, ncol = 0L)
  }

  individual_covariates <- data.frame(
    subject = seq_len(S),
    subject_id = subject_labels,
    stringsAsFactors = FALSE
  )
  if (!is.null(subject_X) && P > 0L) {
    individual_covariates <- cbind(
      individual_covariates,
      as.data.frame(subject_X, check.names = FALSE, stringsAsFactors = FALSE)
    )
  }

  covariate_distributions <- metadata_value_long("covariate_distributions")
  if ((!is.data.frame(covariate_distributions) ||
       nrow(covariate_distributions) != P) && !is.null(subject_X)) {
    covariate_distributions <- data.frame(
      index = seq_len(P),
      name = covariate_names,
      n = rep.int(S, P),
      n_unique = integer(P),
      mean = numeric(P),
      sd = numeric(P),
      min = numeric(P),
      q25 = numeric(P),
      median = numeric(P),
      q75 = numeric(P),
      max = numeric(P),
      n_zero = integer(P),
      n_one = integer(P),
      stringsAsFactors = FALSE
    )
    if (P > 0L) {
      for (j in seq_len(P)) {
        xj <- as.numeric(subject_X[, j])
        qj <- stats::quantile(xj, c(0.25, 0.50, 0.75), names = FALSE)
        covariate_distributions$n_unique[j] <- length(unique(xj))
        covariate_distributions$mean[j] <- mean(xj)
        covariate_distributions$sd[j] <- stats::sd(xj)
        covariate_distributions$min[j] <- min(xj)
        covariate_distributions$q25[j] <- qj[[1L]]
        covariate_distributions$median[j] <- qj[[2L]]
        covariate_distributions$q75[j] <- qj[[3L]]
        covariate_distributions$max[j] <- max(xj)
        covariate_distributions$n_zero[j] <- sum(xj == 0)
        covariate_distributions$n_one[j] <- sum(xj == 1)
      }
    }
  }
  if (!is.data.frame(covariate_distributions)) {
    covariate_distributions <- data.frame()
  }

  reference_profile_default <- "X = 0 on the encoded design-matrix scale"
  if (is.list(covariate_metadata_object) &&
      is.list(covariate_metadata_object$reference_profile)) {
    reference_description <-
      covariate_metadata_object$reference_profile$description
    if (is.character(reference_description) &&
        length(reference_description) == 1L &&
        !is.na(reference_description) && nzchar(reference_description)) {
      reference_profile_default <- reference_description
    }
  }
  population_reference_profile <- metadata_value_long(
    "population_reference_profile",
    reference_profile_default
  )

  between_arm_threshold <- if (!is.null(stan_data) &&
                               "meaningful_between_arm_difference" %in% names(stan_data)) {
    as.numeric(stan_data$meaningful_between_arm_difference)[1]
  } else {
    NA_real_
  }

  prior_names <- c(
    "baseline_prior_mean", "baseline_prior_sd",
    "beta_time_prior_mean", "beta_time_prior_sd",
    "tau_common_prior_rate",
    "beta_treatment_prior_sd", "tau_treatment_prior_rate",
    "arm_baseline_sd_prior_rate",
    "covariate_baseline_prior_sd", "beta_covariate_prior_sd",
    "tau_covariate_prior_rate",
    "sigma_intercept_prior_rate", "sigma_slope_prior_rate",
    "sigma_prior_rate", "nu_prior_shape", "nu_prior_rate"
  )

  prior_data <- NULL
  if (!is.null(fit_info) && is.list(fit_info$prior_data) &&
      all(prior_names %in% names(fit_info$prior_data))) {
    prior_data <- fit_info$prior_data[prior_names]
  } else if (!is.null(stan_data) && all(prior_names %in% names(stan_data))) {
    prior_data <- stan_data[prior_names]
  }

  # Current model: MCID is a parameter. Fixed meaningful_change is fallback only.
  if ("mcid" %in% names(draws)) {
    mcid_draws <- as.numeric(draws$mcid)
  } else if (!is.null(meaningful_change)) {
    if (length(meaningful_change) != 1 ||
        !is.numeric(meaningful_change) ||
        !is.finite(meaningful_change) || meaningful_change < 0) {
      stop("`meaningful_change` must be one finite non-negative numeric value.", call. = FALSE)
    }
    mcid_draws <- rep(as.numeric(meaningful_change), nrow(draws))
  } else {
    stop(
      "The fitted model does not contain posterior variable `mcid`. ",
      "For legacy fits, supply `meaningful_change`.",
      call. = FALSE
    )
  }

  mcid_summary <- summarize_vector_long(mcid_draws, "mcid")

  # ============================================================
  # SCALAR POPULATION / VARIANCE PARAMETERS
  # ============================================================

  scalar_parameters <- c(
    "baseline_mean",
    "beta_time",
    "tau_common",
    "arm_baseline_sd",
    "sigma",
    "nu",
    "nu_value",
    "sigma_intercept",
    "sigma_slope",
    "rho_subject",
    "residual_sd",
    "residual_cv",
    "mcid"
  )

  if (likelihood_id != 1L) {
    scalar_parameters <- setdiff(scalar_parameters, c("nu", "nu_value"))
  }
  if (likelihood_id != 3L) {
    scalar_parameters <- setdiff(scalar_parameters, "residual_cv")
  }

  scalar_parameters <- scalar_parameters[scalar_parameters %in% names(draws)]

  population_summary <- do.call(
    rbind,
    lapply(scalar_parameters, function(parameter) {
      summarize_vector_long(draws[[parameter]], parameter)
    })
  )

  if (is.null(population_summary)) {
    population_summary <- data.frame()
  }

  beta_treatment_summary <- summarize_indexed_long("beta_treatment", "contrast")
  tau_treatment_summary <- summarize_indexed_long("tau_treatment", "contrast")
  arm_baseline_offset_summary <- summarize_indexed_long("arm_baseline_offset", "contrast")

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

  population_time_means <- summarize_indexed_long("population_mean", c("arm", "time"))

  if (nrow(population_time_means) > 0) {
    population_time_means$arm_label <- arm_labels[population_time_means$arm]
    population_time_means$time_value <- time_value[population_time_means$time]
    population_time_means <- population_time_means[order(
      population_time_means$arm,
      population_time_means$time
    ), ]
  }

  add_population_metadata_long <- function(x) {
    if (nrow(x) > 0) {
      x$arm_label <- arm_labels[x$arm]
      x$time_value <- time_value[x$time]
      x <- x[order(x$arm, x$time), , drop = FALSE]
    }
    x
  }

  population_model_location <- add_population_metadata_long(
    summarize_indexed_long("population_model_location", c("arm", "time"))
  )
  population_median <- add_population_metadata_long(
    summarize_indexed_long("population_median", c("arm", "time"))
  )
  population_ratio <- add_population_metadata_long(
    summarize_indexed_long("population_ratio_from_baseline", c("arm", "time"))
  )
  population_percent_change <- add_population_metadata_long(
    summarize_indexed_long(
      "population_percent_change_from_baseline",
      c("arm", "time")
    )
  )

  population_change <- summarize_indexed_long(
    "population_change_from_baseline",
    c("arm", "time")
  )

  directional_change_cols <- get_cols_long("directional_population_change")
  directional_change_idx <- parse_indices_long(directional_change_cols)

  if (nrow(population_change) > 0) {
    population_change$arm_label <- arm_labels[population_change$arm]
    population_change$time_value <- time_value[population_change$time]
    population_change$P_positive_raw <- NA_real_
    population_change$P_negative_raw <- NA_real_
    population_change$P_improvement <- NA_real_
    population_change$P_responder <- NA_real_
    population_change$mean_directional_change_minus_mcid <- NA_real_

    raw_change_cols <- get_cols_long("population_change_from_baseline")
    raw_change_idx <- parse_indices_long(raw_change_cols)

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

  population_standardized_change <- summarize_indexed_long(
    "standardized_population_change",
    c("arm", "time")
  )

  if (nrow(population_standardized_change) > 0) {
    population_standardized_change$arm_label <-
      arm_labels[population_standardized_change$arm]
    population_standardized_change$time_value <-
      time_value[population_standardized_change$time]
  }

  directional_population_change <- summarize_indexed_long(
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
  # DYNAMIC COVARIATE PARAMETERS AND TIME-SPECIFIC EFFECTS
  # ============================================================

  covariate_effect_interpretation_long <- function(j) {
    row <- covariate_map[j, , drop = FALSE]
    level <- as.character(row$level[[1L]])
    reference <- as.character(row$reference_level[[1L]])
    encoding <- tolower(as.character(row$encoding[[1L]]))
    scale <- suppressWarnings(as.numeric(row$scale[[1L]]))
    center <- suppressWarnings(as.numeric(row$center[[1L]]))
    unit <- as.character(row$unit[[1L]])

    if (!is.na(level) && nzchar(level) &&
        !is.na(reference) && nzchar(reference)) {
      return(paste0(level, " versus ", reference, " (encoded 0 -> 1)"))
    }
    if (grepl("dummy|treatment|binary|logical", encoding) &&
        !is.na(reference) && nzchar(reference)) {
      return(paste0("encoded 0 -> 1 relative to ", reference))
    }
    if (is.finite(scale) && abs(scale - 1) > sqrt(.Machine$double.eps)) {
      center_text <- if (is.finite(center)) paste0("; centered at ", signif(center, 6L)) else ""
      return(paste0(
        "+1 encoded unit (= +", signif(scale, 6L),
        " original-scale units", center_text, ")"
      ))
    }
    if (!is.na(unit) && nzchar(unit)) {
      return(paste0("+1 original-scale unit (", unit, ")"))
    }
    "+1 original-scale unit"
  }

  covariate_interpretations <- if (P > 0L) {
    vapply(seq_len(P), covariate_effect_interpretation_long, character(1L))
  } else {
    character(0)
  }

  add_covariate_metadata_long <- function(x, include_time = FALSE) {
    if (!"covariate" %in% names(x)) {
      stop("Internal error: covariate-indexed summary lacks `covariate`.", call. = FALSE)
    }
    idx <- x$covariate
    x$encoded_term <- covariate_names[idx]
    x$original_covariate <- as.character(covariate_map$original_name[idx])
    x$label <- as.character(covariate_map$label[idx])
    x$type <- as.character(covariate_map$type[idx])
    x$encoding <- as.character(covariate_map$encoding[idx])
    x$level <- as.character(covariate_map$level[idx])
    x$reference_level <- as.character(covariate_map$reference_level[idx])
    x$center <- suppressWarnings(as.numeric(covariate_map$center[idx]))
    x$scale <- suppressWarnings(as.numeric(covariate_map$scale[idx]))
    x$effect_interpretation <- covariate_interpretations[idx]
    if (isTRUE(include_time)) {
      x$time_value <- time_value[x$time]
      x <- x[order(x$covariate, x$time), , drop = FALSE]
    } else {
      x <- x[order(x$covariate), , drop = FALSE]
    }
    rownames(x) <- NULL
    x
  }

  add_covariate_probability_long <- function(x, prefix, output_name, predicate) {
    if (nrow(x) == 0L) {
      x[[output_name]] <- numeric(0)
      return(x)
    }
    columns <- get_cols_long(prefix)
    idx <- parse_indices_long(columns)
    x[[output_name]] <- vapply(
      seq_len(nrow(x)),
      function(i) {
        wanted <- c(x$covariate[[i]], x$time[[i]])
        pos <- which(vapply(
          idx,
          function(z) length(z) == 2L && all(z == wanted),
          logical(1L)
        ))
        if (length(pos) != 1L) return(NA_real_)
        mean(predicate(as.numeric(draws[[columns[[pos]]]])), na.rm = TRUE)
      },
      numeric(1L)
    )
    x
  }

  covariate_baseline_summary <- add_covariate_metadata_long(
    summarize_indexed_long("covariate_baseline_effect", "covariate")
  )
  beta_covariate_summary <- add_covariate_metadata_long(
    summarize_indexed_long("beta_covariate_time", "covariate")
  )
  tau_covariate_summary <- add_covariate_metadata_long(
    summarize_indexed_long("tau_covariate", "covariate")
  )

  covariate_model_effect <- add_covariate_metadata_long(
    summarize_indexed_long("covariate_effect", c("covariate", "time")),
    include_time = TRUE
  )
  covariate_level_difference <- add_covariate_metadata_long(
    summarize_indexed_long(
      "covariate_population_mean_difference",
      c("covariate", "time")
    ),
    include_time = TRUE
  )
  covariate_change_difference <- add_covariate_metadata_long(
    summarize_indexed_long(
      "covariate_population_mean_change_difference",
      c("covariate", "time")
    ),
    include_time = TRUE
  )
  covariate_directional_change_difference <- add_covariate_metadata_long(
    summarize_indexed_long(
      "directional_covariate_population_mean_change_difference",
      c("covariate", "time")
    ),
    include_time = TRUE
  )

  covariate_level_difference <- add_covariate_probability_long(
    covariate_level_difference,
    "covariate_population_mean_difference",
    "P_level_difference_gt_0",
    function(z) z > 0
  )
  covariate_change_difference <- add_covariate_probability_long(
    covariate_change_difference,
    "covariate_population_mean_change_difference",
    "P_change_difference_gt_0",
    function(z) z > 0
  )
  covariate_directional_change_difference <- add_covariate_probability_long(
    covariate_directional_change_difference,
    "directional_covariate_population_mean_change_difference",
    "P_directional_change_gt_0",
    function(z) z > 0
  )

  add_effect_stats_long <- function(base, extra, prefix) {
    statistic_names <- c("mean", "median", "sd", "mad", "lower", "upper", "CrI_width")
    key <- paste(base$covariate, base$time, sep = ":")
    extra_key <- paste(extra$covariate, extra$time, sep = ":")
    matched <- match(key, extra_key)
    for (nm in statistic_names) {
      base[[paste0(prefix, "_", nm)]] <- extra[[nm]][matched]
    }
    base
  }

  covariate_effects <- covariate_level_difference
  if (nrow(covariate_effects) == 0L && nrow(covariate_model_effect) > 0L) {
    covariate_effects <- covariate_model_effect
    covariate_effects$estimand_scale <- rep(
      "model/link scale", nrow(covariate_effects)
    )
  } else {
    covariate_effects$estimand_scale <- rep(
      "natural-outcome mean difference", nrow(covariate_effects)
    )
  }
  if (nrow(covariate_effects) > 0L) {
    covariate_effects <- add_effect_stats_long(
      covariate_effects, covariate_model_effect, "model_scale"
    )
    covariate_effects <- add_effect_stats_long(
      covariate_effects, covariate_change_difference, "change_from_baseline"
    )
    covariate_effects <- add_effect_stats_long(
      covariate_effects,
      covariate_directional_change_difference,
      "directional_change_from_baseline"
    )
    change_key <- paste(
      covariate_change_difference$covariate,
      covariate_change_difference$time,
      sep = ":"
    )
    directional_key <- paste(
      covariate_directional_change_difference$covariate,
      covariate_directional_change_difference$time,
      sep = ":"
    )
    base_key <- paste(covariate_effects$covariate, covariate_effects$time, sep = ":")
    covariate_effects$P_change_difference_gt_0 <-
      covariate_change_difference$P_change_difference_gt_0[
        match(base_key, change_key)
      ]
    covariate_effects$P_directional_change_gt_0 <-
      covariate_directional_change_difference$P_directional_change_gt_0[
        match(base_key, directional_key)
      ]
  }


  # ============================================================
  # ALL PAIRWISE TEMPORAL CONTRASTS
  # ============================================================
  #
  # Derived directly from posterior draws already generated by Stan.
  # For K time points, return all unique forward comparisons:
  # t0 -> t1, t0 -> t2, ..., t(K-2) -> t(K-1).
  # Reverse comparisons are omitted because they contain the same
  # information with the opposite sign.
  # ============================================================

  time_pairs <- if (K >= 2L) {
    utils::combn(seq_len(K), 2L, simplify = FALSE)
  } else {
    list()
  }

  time_label_pairwise_long <- function(k) paste0("t", as.integer(k) - 1L)

  indexed_draw_column_long <- function(prefix, wanted_index) {
    columns <- get_cols_long(prefix)
    idx <- parse_indices_long(columns)

    pos <- which(vapply(
      idx,
      function(z) {
        length(z) == length(wanted_index) &&
          all(z == wanted_index)
      },
      logical(1)
    ))

    if (length(pos) != 1L) return(NA_character_)
    columns[[pos]]
  }

  posterior_pairwise_stats_long <- function(z) {
    s <- posterior_stats_long(z)
    list(
      mean = unname(s["mean"]),
      median = unname(s["median"]),
      sd = unname(s["sd"]),
      mad = unname(s["mad"]),
      lower = unname(s["lower"]),
      upper = unname(s["upper"]),
      CrI_width = unname(s["CrI_width"])
    )
  }

  direction_available <- length(direction) == 1L && !is.na(direction)

  # ------------------------------------------------------------
  # Generic population change: every arm x every unique time pair
  # ------------------------------------------------------------

  population_pairwise_rows <- list()
  pair_counter <- 1L

  for (g in seq_len(G)) {
    for (pair in time_pairs) {
      from <- pair[[1L]]
      to <- pair[[2L]]

      col_from <- indexed_draw_column_long("population_mean", c(g, from))
      col_to <- indexed_draw_column_long("population_mean", c(g, to))

      if (is.na(col_from) || is.na(col_to)) next

      z <- as.numeric(draws[[col_to]]) - as.numeric(draws[[col_from]])
      s <- posterior_pairwise_stats_long(z)

      directional_z <- if (direction_available) {
        direction * z
      } else {
        rep(NA_real_, length(z))
      }

      population_pairwise_rows[[pair_counter]] <- data.frame(
        arm = g,
        arm_label = arm_labels[g],
        from = from,
        to = to,
        from_label = time_label_pairwise_long(from),
        to_label = time_label_pairwise_long(to),
        from_time_value = time_value[from],
        to_time_value = time_value[to],
        elapsed = time_value[to] - time_value[from],
        mean = s$mean,
        median = s$median,
        sd = s$sd,
        mad = s$mad,
        lower = s$lower,
        upper = s$upper,
        CrI_width = s$CrI_width,
        mean_directional_change =
          if (direction_available) mean(directional_z) else NA_real_,
        P_positive_raw = mean(z > 0),
        P_negative_raw = mean(z < 0),
        P_improvement =
          if (direction_available) mean(directional_z > 0) else NA_real_,
        P_MCID =
          if (direction_available) mean(directional_z >= mcid_draws) else NA_real_,
        stringsAsFactors = FALSE
      )

      pair_counter <- pair_counter + 1L
    }
  }

  population_pairwise_change <- if (length(population_pairwise_rows) > 0L) {
    do.call(rbind, population_pairwise_rows)
  } else {
    data.frame()
  }

  if (nrow(population_pairwise_change) > 0L) {
    rownames(population_pairwise_change) <- NULL
    population_pairwise_change <- population_pairwise_change[
      order(
        population_pairwise_change$arm,
        population_pairwise_change$from,
        population_pairwise_change$to
      ),
      ,
      drop = FALSE
    ]
  }

  population_consecutive_change <- if (nrow(population_pairwise_change) > 0L) {
    population_pairwise_change[
      population_pairwise_change$to == population_pairwise_change$from + 1L,
      ,
      drop = FALSE
    ]
  } else {
    data.frame()
  }

  # ------------------------------------------------------------
  # Covariates: all pairwise differences in their time-varying
  # natural-scale effects, one table for every encoded term.
  # ------------------------------------------------------------

  covariate_pairwise_rows <- list()
  pair_counter <- 1L

  if (P > 0L) {
    for (p in seq_len(P)) {
      for (pair in time_pairs) {
        from <- pair[[1L]]
        to <- pair[[2L]]

        col_from <- indexed_draw_column_long(
          "covariate_population_mean_difference", c(p, from)
        )
        col_to <- indexed_draw_column_long(
          "covariate_population_mean_difference", c(p, to)
        )

        if (is.na(col_from) || is.na(col_to)) next

        z <- as.numeric(draws[[col_to]]) - as.numeric(draws[[col_from]])
        s <- posterior_pairwise_stats_long(z)
        directional_z <- if (direction_available) {
          direction * z
        } else {
          rep(NA_real_, length(z))
        }

        covariate_pairwise_rows[[pair_counter]] <- data.frame(
          covariate = p,
          encoded_term = covariate_names[p],
          original_covariate = as.character(covariate_map$original_name[p]),
          label = as.character(covariate_map$label[p]),
          type = as.character(covariate_map$type[p]),
          encoding = as.character(covariate_map$encoding[p]),
          level = as.character(covariate_map$level[p]),
          reference_level = as.character(covariate_map$reference_level[p]),
          effect_interpretation = covariate_interpretations[p],
          from = from,
          to = to,
          from_label = time_label_pairwise_long(from),
          to_label = time_label_pairwise_long(to),
          from_time_value = time_value[from],
          to_time_value = time_value[to],
          elapsed = time_value[to] - time_value[from],
          mean = s$mean,
          median = s$median,
          sd = s$sd,
          mad = s$mad,
          lower = s$lower,
          upper = s$upper,
          CrI_width = s$CrI_width,
          P_change_difference_gt_0 = mean(z > 0),
          mean_directional_change_difference =
            if (direction_available) mean(directional_z) else NA_real_,
          P_directional_change_gt_0 =
            if (direction_available) mean(directional_z > 0) else NA_real_,
          stringsAsFactors = FALSE
        )
        pair_counter <- pair_counter + 1L
      }
    }
  }

  covariate_pairwise_change <- if (length(covariate_pairwise_rows) > 0L) {
    do.call(rbind, covariate_pairwise_rows)
  } else {
    data.frame()
  }

  if (nrow(covariate_pairwise_change) > 0L) {
    rownames(covariate_pairwise_change) <- NULL
    covariate_pairwise_change <- covariate_pairwise_change[
      order(
        covariate_pairwise_change$covariate,
        covariate_pairwise_change$from,
        covariate_pairwise_change$to
      ),
      ,
      drop = FALSE
    ]
  }

  covariate_consecutive_change <- if (nrow(covariate_pairwise_change) > 0L) {
    covariate_pairwise_change[
      covariate_pairwise_change$to == covariate_pairwise_change$from + 1L,
      ,
      drop = FALSE
    ]
  } else {
    data.frame()
  }

  # ------------------------------------------------------------
  # Treatment: adjusted level difference at every time and
  # difference in temporal change for every unique time pair
  # ------------------------------------------------------------

  treatment_level_rows <- list()
  treatment_pairwise_rows <- list()
  level_counter <- 1L
  pair_counter <- 1L

  if (G >= 2L) {
    for (g in seq.int(2L, G)) {

      for (k in seq_len(K)) {
        col_ref <- indexed_draw_column_long("population_mean", c(1L, k))
        col_trt <- indexed_draw_column_long("population_mean", c(g, k))

        if (is.na(col_ref) || is.na(col_trt)) next

        z <- as.numeric(draws[[col_trt]]) - as.numeric(draws[[col_ref]])
        s <- posterior_pairwise_stats_long(z)
        directional_z <- if (direction_available) {
          direction * z
        } else {
          rep(NA_real_, length(z))
        }

        treatment_level_rows[[level_counter]] <- data.frame(
          contrast = g - 1L,
          treatment_arm = g,
          treatment_label = arm_labels[g],
          reference_arm = 1L,
          reference_label = arm_labels[1L],
          time = k,
          time_label = time_label_pairwise_long(k),
          time_value = time_value[k],
          mean = s$mean,
          median = s$median,
          sd = s$sd,
          mad = s$mad,
          lower = s$lower,
          upper = s$upper,
          CrI_width = s$CrI_width,
          mean_directional_difference =
            if (direction_available) mean(directional_z) else NA_real_,
          P_treatment_better =
            if (direction_available) mean(directional_z > 0) else NA_real_,
          stringsAsFactors = FALSE
        )

        level_counter <- level_counter + 1L
      }

      for (pair in time_pairs) {
        from <- pair[[1L]]
        to <- pair[[2L]]

        ref_from <- indexed_draw_column_long("population_mean", c(1L, from))
        ref_to <- indexed_draw_column_long("population_mean", c(1L, to))
        trt_from <- indexed_draw_column_long("population_mean", c(g, from))
        trt_to <- indexed_draw_column_long("population_mean", c(g, to))

        if (anyNA(c(ref_from, ref_to, trt_from, trt_to))) next

        ref_change <- as.numeric(draws[[ref_to]]) - as.numeric(draws[[ref_from]])
        trt_change <- as.numeric(draws[[trt_to]]) - as.numeric(draws[[trt_from]])
        z <- trt_change - ref_change
        s <- posterior_pairwise_stats_long(z)

        directional_z <- if (direction_available) {
          direction * z
        } else {
          rep(NA_real_, length(z))
        }

        treatment_pairwise_rows[[pair_counter]] <- data.frame(
          contrast = g - 1L,
          treatment_arm = g,
          treatment_label = arm_labels[g],
          reference_arm = 1L,
          reference_label = arm_labels[1L],
          from = from,
          to = to,
          from_label = time_label_pairwise_long(from),
          to_label = time_label_pairwise_long(to),
          from_time_value = time_value[from],
          to_time_value = time_value[to],
          elapsed = time_value[to] - time_value[from],
          mean = s$mean,
          median = s$median,
          sd = s$sd,
          mad = s$mad,
          lower = s$lower,
          upper = s$upper,
          CrI_width = s$CrI_width,
          mean_directional_benefit =
            if (direction_available) mean(directional_z) else NA_real_,
          P_benefit_positive =
            if (direction_available) mean(directional_z > 0) else NA_real_,
          P_benefit_meaningful =
            if (
              direction_available &&
              length(between_arm_threshold) == 1L &&
              is.finite(between_arm_threshold)
            ) {
              mean(directional_z >= between_arm_threshold)
            } else {
              NA_real_
            },
          stringsAsFactors = FALSE
        )

        pair_counter <- pair_counter + 1L
      }
    }
  }

  treatment_level_difference <- if (length(treatment_level_rows) > 0L) {
    do.call(rbind, treatment_level_rows)
  } else {
    data.frame()
  }

  treatment_pairwise_change <- if (length(treatment_pairwise_rows) > 0L) {
    do.call(rbind, treatment_pairwise_rows)
  } else {
    data.frame()
  }

  if (nrow(treatment_level_difference) > 0L) {
    rownames(treatment_level_difference) <- NULL
  }

  if (nrow(treatment_pairwise_change) > 0L) {
    rownames(treatment_pairwise_change) <- NULL
    treatment_pairwise_change <- treatment_pairwise_change[
      order(
        treatment_pairwise_change$contrast,
        treatment_pairwise_change$from,
        treatment_pairwise_change$to
      ),
      ,
      drop = FALSE
    ]
  }

  treatment_consecutive_change <- if (nrow(treatment_pairwise_change) > 0L) {
    treatment_pairwise_change[
      treatment_pairwise_change$to == treatment_pairwise_change$from + 1L,
      ,
      drop = FALSE
    ]
  } else {
    data.frame()
  }

  # ============================================================
  # NEW-SUBJECT RESPONDER ESTIMANDS
  # ============================================================

  latent_any_improvement <- probability_indexed_long(
    "new_subject_latent_any_improvement_draw",
    c("arm", "time"),
    "probability"
  )

  latent_responder <- probability_indexed_long(
    "new_subject_latent_responder_draw",
    c("arm", "time"),
    "probability"
  )

  new_subject_latent_responder <- probability_indexed_long(
    "new_subject_latent_responder_draw",
    c("arm", "time"),
    "posterior_predictive_probability"
  )

  new_subject_predictive_responder <- probability_indexed_long(
    "new_subject_predictive_responder_draw",
    c("arm", "time"),
    "posterior_predictive_probability"
  )

  new_subject_latent_change <- summarize_indexed_long(
    "new_subject_latent_change_draw",
    c("arm", "time")
  )

  new_subject_predictive_change <- summarize_indexed_long(
    "new_subject_predictive_change_draw",
    c("arm", "time")
  )

  add_arm_time_metadata_long <- function(x) {
    if (nrow(x) > 0) {
      x$arm_label <- arm_labels[x$arm]
      x$time_value <- time_value[x$time]
      x$covariate_profile <- population_reference_profile
      x <- x[order(x$arm, x$time), ]
    }
    x
  }

  latent_any_improvement <- add_arm_time_metadata_long(latent_any_improvement)
  latent_responder <- add_arm_time_metadata_long(latent_responder)
  new_subject_latent_responder <- add_arm_time_metadata_long(new_subject_latent_responder)
  new_subject_predictive_responder <- add_arm_time_metadata_long(new_subject_predictive_responder)
  new_subject_latent_change <- add_arm_time_metadata_long(new_subject_latent_change)
  new_subject_predictive_change <- add_arm_time_metadata_long(new_subject_predictive_change)

  # ============================================================
  # INDIVIDUAL CHANGE AND CLINICAL RESPONSE
  # ============================================================

  individual_change <- summarize_indexed_long(
    "individual_change_from_baseline",
    c("time", "subject")
  )

  individual_directional_change <- summarize_indexed_long(
    "individual_directional_change",
    c("time", "subject")
  )

  individual_any <- probability_indexed_long(
    "individual_any_improvement_draw",
    c("time", "subject"),
    "P_improvement"
  )

  individual_responder <- probability_indexed_long(
    "individual_meaningful_responder_draw",
    c("time", "subject"),
    "P_MCID"
  )

  individual_distance <- summarize_indexed_long(
    "individual_change_minus_mcid",
    c("time", "subject")
  )

  add_subject_metadata_long <- function(x) {
    if (nrow(x) > 0) {
      x$subject_id <- subject_labels[x$subject]
      x$time_value <- time_value[x$time]
      x$arm <- subject_arm[x$subject]
      x$arm_label <- ifelse(
        is.na(x$arm),
        NA_character_,
        arm_labels[x$arm]
      )
      x <- x[order(x$subject, x$time), ]
    }
    x
  }

  individual_change <- add_subject_metadata_long(individual_change)
  individual_directional_change <- add_subject_metadata_long(individual_directional_change)
  individual_any <- add_subject_metadata_long(individual_any)
  individual_responder <- add_subject_metadata_long(individual_responder)
  individual_distance <- add_subject_metadata_long(individual_distance)

  individual_clinical <- data.frame()

  if (nrow(individual_responder) > 0) {
    individual_clinical <- individual_responder[
      , c(
        "time", "subject", "subject_id", "P_MCID", "time_value", "arm", "arm_label"
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

  treatment_effects <- summarize_indexed_long(
    "treatment_change_difference",
    c("contrast", "time")
  )

  directional_treatment <- summarize_indexed_long(
    "directional_treatment_benefit",
    c("contrast", "time")
  )

  positive_treatment <- probability_indexed_long(
    "treatment_benefit_positive_draw",
    c("contrast", "time"),
    "P_benefit_positive"
  )

  meaningful_treatment <- probability_indexed_long(
    "treatment_benefit_meaningful_draw",
    c("contrast", "time"),
    "P_benefit_meaningful"
  )

  responder_uplift <- summarize_indexed_long(
    "latent_responder_probability_difference",
    c("contrast", "time")
  )

  treatment_ratio_of_ratios <- summarize_indexed_long(
    "treatment_ratio_of_ratios",
    c("contrast", "time")
  )

  if (nrow(treatment_ratio_of_ratios) > 0) {
    treatment_ratio_of_ratios$treatment_arm <-
      treatment_ratio_of_ratios$contrast + 1L
    treatment_ratio_of_ratios$treatment_label <-
      arm_labels[treatment_ratio_of_ratios$treatment_arm]
    treatment_ratio_of_ratios$reference_label <- arm_labels[1L]
    treatment_ratio_of_ratios$time_value <-
      time_value[treatment_ratio_of_ratios$time]
  }

  if (nrow(treatment_effects) > 0) {
    treatment_effects$treatment_arm <- treatment_effects$contrast + 1L
    treatment_effects$treatment_label <- arm_labels[treatment_effects$treatment_arm]
    treatment_effects$reference_arm <- 1L
    treatment_effects$reference_label <- arm_labels[1]
    treatment_effects$time_value <- time_value[treatment_effects$time]

    treatment_key <- paste(treatment_effects$contrast, treatment_effects$time, sep = ":")

    add_from_table_long <- function(base, table, value_column, output_name) {
      if (nrow(table) == 0) {
        base[[output_name]] <- NA_real_
        return(base)
      }
      key <- paste(table$contrast, table$time, sep = ":")
      base[[output_name]] <- table[[value_column]][match(treatment_key, key)]
      base
    }

    treatment_effects <- add_from_table_long(
      treatment_effects, directional_treatment, "mean", "mean_directional_benefit"
    )
    treatment_effects <- add_from_table_long(
      treatment_effects, positive_treatment, "P_benefit_positive", "P_benefit_positive"
    )
    treatment_effects <- add_from_table_long(
      treatment_effects, meaningful_treatment, "P_benefit_meaningful", "P_benefit_meaningful"
    )
    treatment_effects <- add_from_table_long(
      treatment_effects, responder_uplift, "mean", "mean_responder_probability_difference"
    )
    treatment_effects <- add_from_table_long(
      treatment_effects,
      treatment_ratio_of_ratios,
      "mean",
      "mean_ratio_of_ratios"
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
    "tau_common"
  )

  heterogeneity_parameters <- heterogeneity_parameters[
    heterogeneity_parameters %in% names(draws)
  ]

  heterogeneity <- do.call(
    rbind,
    lapply(heterogeneity_parameters, function(parameter) {
      summarize_vector_long(draws[[parameter]], parameter)
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
      covariate_profile = population_reference_profile,
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

  y_rep_cols <- get_cols_long("y_rep")
  ppc <- NULL

  if (length(y_rep_cols) > 0 && !is.null(y)) {
    y_rep_indices <- vapply(
      parse_indices_long(y_rep_cols),
      function(z) z[[1L]],
      integer(1L)
    )
    y_rep_cols <- y_rep_cols[order(y_rep_indices)]
    y_rep_matrix <- as.matrix(draws[, y_rep_cols, drop = FALSE])

    safe_skewness_long <- function(z) {
      z <- as.numeric(z)
      s <- stats::sd(z)
      n <- length(z)
      if (n < 3L || !is.finite(s) || s <= 0) return(NA_real_)
      n / ((n - 1) * (n - 2)) * sum(((z - mean(z)) / s)^3)
    }

    predictive_mean <- rowMeans(y_rep_matrix)
    predictive_sd <- apply(y_rep_matrix, 1, stats::sd)
    predictive_median <- apply(y_rep_matrix, 1, stats::median)
    predictive_q05 <- apply(y_rep_matrix, 1, stats::quantile, probs = 0.05)
    predictive_q95 <- apply(y_rep_matrix, 1, stats::quantile, probs = 0.95)
    predictive_min <- apply(y_rep_matrix, 1, min)
    predictive_max <- apply(y_rep_matrix, 1, max)
    predictive_skewness <- apply(y_rep_matrix, 1, safe_skewness_long)
    predictive_cv <- predictive_sd / pmax(abs(predictive_mean), 1e-12)

    observed <- c(
      mean = mean(y),
      sd = stats::sd(y),
      median = stats::median(y),
      q05 = as.numeric(stats::quantile(y, 0.05)),
      q95 = as.numeric(stats::quantile(y, 0.95)),
      min = min(y),
      max = max(y),
      skewness = safe_skewness_long(y),
      coefficient_of_variation = stats::sd(y) / max(abs(mean(y)), 1e-12)
    )

    posterior_predictive <- list(
      mean = posterior_stats_long(predictive_mean),
      sd = posterior_stats_long(predictive_sd),
      median = posterior_stats_long(predictive_median),
      q05 = posterior_stats_long(predictive_q05),
      q95 = posterior_stats_long(predictive_q95),
      min = posterior_stats_long(predictive_min),
      max = posterior_stats_long(predictive_max),
      skewness = posterior_stats_long(predictive_skewness),
      coefficient_of_variation = posterior_stats_long(predictive_cv)
    )

    lag1_observed <- NA_real_
    lag1_predictive <- NULL
    if (!is.null(stan_data) &&
        !is.null(stan_data$subject) &&
        !is.null(stan_data$time) &&
        length(stan_data$subject) == length(y) &&
        length(stan_data$time) == length(y)) {
      order_index <- order(stan_data$subject, stan_data$time)
      ordered_subject <- as.integer(stan_data$subject)[order_index]
      ordered_time <- as.integer(stan_data$time)[order_index]
      adjacent <- which(
        ordered_subject[-length(ordered_subject)] ==
          ordered_subject[-1L] &
          ordered_time[-length(ordered_time)] + 1L == ordered_time[-1L]
      )

      if (length(adjacent) >= 3L) {
        from_index <- order_index[adjacent]
        to_index <- order_index[adjacent + 1L]
        lag1_observed <- stats::cor(y[from_index], y[to_index])
        lag1_predictive <- posterior_stats_long(apply(
          y_rep_matrix,
          1,
          function(z) stats::cor(z[from_index], z[to_index])
        ))
      }
    }

    support_violation <- matrix(FALSE, nrow = nrow(y_rep_matrix), ncol = N)
    if (has_lower_bound) {
      support_violation <- support_violation | y_rep_matrix < lower_bound
    }
    if (has_upper_bound) {
      support_violation <- support_violation | y_rep_matrix > upper_bound
    }
    if (likelihood_id == 3L) {
      support_violation <- support_violation | y_rep_matrix <= 0
    }

    family_checks <- list(
      likelihood = likelihood,
      modeling_scale = modeling_scale,
      lower_boundary = if (has_lower_bound) lower_bound else NA_real_,
      upper_boundary = if (has_upper_bound) upper_bound else NA_real_,
      observed_lower_boundary_fraction = if (has_lower_bound) {
        mean(y == lower_bound)
      } else {
        NA_real_
      },
      observed_upper_boundary_fraction = if (has_upper_bound) {
        mean(y == upper_bound)
      } else {
        NA_real_
      },
      predictive_lower_boundary_fraction = if (has_lower_bound) {
        posterior_stats_long(rowMeans(y_rep_matrix == lower_bound))
      } else {
        NULL
      },
      predictive_upper_boundary_fraction = if (has_upper_bound) {
        posterior_stats_long(rowMeans(y_rep_matrix == upper_bound))
      } else {
        NULL
      },
      predictive_support_violation_fraction = posterior_stats_long(
        rowMeans(support_violation)
      ),
      observed_longitudinal_lag1_correlation = lag1_observed,
      predictive_longitudinal_lag1_correlation = lag1_predictive
    )

    if (likelihood_id == 3L) {
      log_y <- log(y)
      log_y_rep <- log(y_rep_matrix)
      family_checks$log_scale <- list(
        observed = c(
          mean = mean(log_y),
          sd = stats::sd(log_y),
          skewness = safe_skewness_long(log_y)
        ),
        posterior_predictive = list(
          mean = posterior_stats_long(rowMeans(log_y_rep)),
          sd = posterior_stats_long(apply(log_y_rep, 1, stats::sd)),
          skewness = posterior_stats_long(apply(log_y_rep, 1, safe_skewness_long))
        )
      )
    }

    by_time <- NULL
    if (!is.null(stan_data) &&
        !is.null(stan_data$time) &&
        length(stan_data$time) == length(y)) {
      by_time <- lapply(seq_len(K), function(k) {
        ix <- which(as.integer(stan_data$time) == k)
        observed_k <- y[ix]
        replicated_k <- y_rep_matrix[, ix, drop = FALSE]

        list(
          time = k,
          time_value = time_value[k],
          observed = c(
            mean = mean(observed_k),
            sd = stats::sd(observed_k),
            skewness = safe_skewness_long(observed_k)
          ),
          posterior_predictive = list(
            mean = posterior_stats_long(rowMeans(replicated_k)),
            sd = posterior_stats_long(apply(replicated_k, 1, stats::sd)),
            skewness = posterior_stats_long(
              apply(replicated_k, 1, safe_skewness_long)
            )
          )
        )
      })
    }

    ppc <- list(
      observed = observed,
      posterior_predictive = posterior_predictive,
      bayesian_p_values = c(
        mean = mean(predictive_mean >= mean(y)),
        sd = mean(predictive_sd >= stats::sd(y)),
        median = mean(predictive_median >= stats::median(y)),
        skewness = mean(predictive_skewness >= safe_skewness_long(y), na.rm = TRUE)
      ),
      family_checks = family_checks,
      by_time = by_time,
      predictive_draws = y_rep_matrix
    )
  }

  # ============================================================
  # LOG-LIKELIHOOD INFORMATION
  # ============================================================

  log_lik_cols <- get_cols_long("log_lik")
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

  model_information <- list(
    outcome = outcome,
    outcome_name = outcome_name,
    likelihood = likelihood,
    likelihood_id = likelihood_id,
    link = if (likelihood_id == 3L) "log" else "identity",
    modeling_scale = modeling_scale,
    clinical_estimand_scale = "natural outcome units",
    population_mean_definition = if (likelihood_id == 3L) {
      paste0(
        "marginal uncensored natural-scale mean integrating Gaussian subject ",
        "effects and log-normal observation variance, then mapped to any ",
        "declared observable endpoints"
      )
    } else if (has_lower_bound || has_upper_bound) {
      "bounded latent location (not the exact mean of the censored distribution)"
    } else {
      "natural-scale population location/mean"
    },
    standardized_change_definition = if (likelihood_id == 3L) {
      "change in log location divided by residual log-SD"
    } else {
      "natural-scale change divided by residual SD"
    },
    outcome_bounds = c(lower = lower_bound, upper = upper_bound),
    boundary_strategy = if (likelihood_id == 3L &&
                            !has_upper_bound &&
                            (!has_lower_bound || lower_bound <= 0)) {
      "strictly positive log-normal support"
    } else if (likelihood_id == 3L) {
      "strictly positive log-normal support with endpoint censoring"
    } else if (has_lower_bound || has_upper_bound) {
      "endpoint censoring"
    } else {
      "none"
    },
    outcome_diagnostics = if (!is.null(stan_data) &&
                              !is.null(stan_data$outcome_diagnostics)) {
      stan_data$outcome_diagnostics
    } else {
      NULL
    },
    n_observations = N,
    n_subjects = S,
    n_time_points = K,
    n_arms = G,
    n_covariates = P,
    P = P,
    subject_labels = subject_labels,
    arm_labels = arm_labels,
    reference_arm = 1L,
    reference_label = arm_labels[1],
    covariate_names = covariate_names,
    covariate_original_names = covariate_original_names,
    covariate_labels = covariate_labels,
    covariate_types = covariate_types,
    covariate_reference_levels = covariate_reference_levels,
    covariate_centers = covariate_centers,
    covariate_scales = covariate_scales,
    covariate_map = covariate_map,
    covariate_metadata = metadata_value_long("covariate_metadata"),
    covariates_requested = metadata_value_long("covariates_requested"),
    covariates_selected = metadata_value_long("covariates_selected"),
    covariate_distributions = covariate_distributions,
    population_reference_profile = population_reference_profile,
    covariate_effect_definition = paste0(
      "Each term is evaluated as a +1 design-matrix-unit change while all ",
      "other covariates are fixed at X = 0. Natural-scale mean differences ",
      "are computed at the reference treatment arm."
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
    prior_source = if (!is.null(fit_info) && !is.null(fit_info$prior_source)) {
      fit_info$prior_source
    } else {
      NA_character_
    },
    prior_profile = if (!is.null(fit_info) && !is.null(fit_info$prior_profile)) {
      fit_info$prior_profile
    } else {
      NA_character_
    },
    prior_data = prior_data,
    covariate_prior_structure = list(
      baseline = "covariate_baseline_prior_sd",
      linear_time = "beta_covariate_prior_sd",
      nonlinear_time = "tau_covariate_prior_rate",
      dimension = P,
      terms = covariate_names
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
      get_cols_long("population_model_location"),
      get_cols_long("population_median"),
      get_cols_long("population_mean"),
      get_cols_long("population_change_from_baseline"),
      get_cols_long("directional_population_change"),
      get_cols_long("standardized_population_change"),
      get_cols_long("population_ratio_from_baseline"),
      get_cols_long("population_percent_change_from_baseline")
    ),
    covariates = c(
      get_cols_long("covariate_baseline_effect"),
      get_cols_long("beta_covariate_time"),
      get_cols_long("tau_covariate"),
      get_cols_long("covariate_effect"),
      get_cols_long("covariate_population_mean_difference"),
      get_cols_long("covariate_population_mean_change_difference"),
      get_cols_long("directional_covariate_population_mean_change_difference")
    ),
    new_subject = c(
      get_cols_long("new_subject_latent_any_improvement_draw"),
      get_cols_long("new_subject_latent_change_draw"),
      get_cols_long("new_subject_latent_responder_draw"),
      get_cols_long("new_subject_predictive_any_improvement_draw"),
      get_cols_long("new_subject_predictive_change_draw"),
      get_cols_long("new_subject_predictive_responder_draw")
    ),
    individual = c(
      get_cols_long("individual_change_from_baseline"),
      get_cols_long("individual_directional_change"),
      get_cols_long("individual_any_improvement_draw"),
      get_cols_long("individual_meaningful_responder_draw"),
      get_cols_long("individual_change_minus_mcid")
    ),
    treatment = c(
      get_cols_long("treatment_change_difference"),
      get_cols_long("directional_treatment_benefit"),
      get_cols_long("treatment_benefit_positive_draw"),
      get_cols_long("treatment_benefit_meaningful_draw"),
      get_cols_long("latent_responder_probability_difference"),
      get_cols_long("treatment_ratio_of_ratios")
    ),
    observation_level = c(
      get_cols_long("log_lik"),
      get_cols_long("y_rep")
    )
  )

  tag_covariate_parameter_long <- function(x, component) {
    x$component <- rep(component, nrow(x))
    x
  }
  covariate_parameter_summary <- do.call(
    rbind,
    list(
      tag_covariate_parameter_long(covariate_baseline_summary, "baseline"),
      tag_covariate_parameter_long(beta_covariate_summary, "linear_time"),
      tag_covariate_parameter_long(tau_covariate_summary, "rw1_scale")
    )
  )
  rownames(covariate_parameter_summary) <- NULL

  add_legacy_group_columns_long <- function(x) {
    if (!is.data.frame(x)) return(x)
    x$comparison_group <- ifelse(
      is.na(x$level) | !nzchar(x$level), x$label, x$level
    )
    x$reference_group <- x$reference_level
    if ("P_level_difference_gt_0" %in% names(x)) {
      x$P_comparison_minus_reference_gt_0 <- x$P_level_difference_gt_0
    }
    if ("P_directional_change_gt_0" %in% names(x)) {
      x$P_comparison_has_better_change <- x$P_directional_change_gt_0
    }
    x
  }

  legacy_covariate_alias_long <- function(index, threshold = NULL) {
    if (length(index) != 1L) return(NULL)
    level <- add_legacy_group_columns_long(
      covariate_level_difference[
        covariate_level_difference$covariate == index,
        ,
        drop = FALSE
      ]
    )
    from_baseline <- add_legacy_group_columns_long(
      covariate_change_difference[
        covariate_change_difference$covariate == index,
        ,
        drop = FALSE
      ]
    )
    directional <- add_legacy_group_columns_long(
      covariate_directional_change_difference[
        covariate_directional_change_difference$covariate == index,
        ,
        drop = FALSE
      ]
    )
    pairwise <- add_legacy_group_columns_long(
      covariate_pairwise_change[
        covariate_pairwise_change$covariate == index,
        ,
        drop = FALSE
      ]
    )
    family <- list(
      level = level,
      change = pairwise,
      from_baseline = from_baseline,
      directional_from_baseline = directional,
      consecutive = if (nrow(pairwise) > 0L) {
        pairwise[pairwise$to == pairwise$from + 1L, , drop = FALSE]
      } else {
        data.frame()
      }
    )
    effects <- list(
      level_difference = level,
      change_difference = from_baseline,
      directional_change_difference = directional
    )
    if (!is.null(threshold)) {
      family$threshold <- threshold
      effects$threshold <- threshold
    }
    list(family = family, effects = effects)
  }

  normalized_original_names <- tolower(trimws(covariate_map$original_name))
  has_clean_binary_contrast <- !is.na(covariate_map$level) &
    nzchar(covariate_map$level) &
    !is.na(covariate_map$reference_level) &
    nzchar(covariate_map$reference_level)
  gender_index <- which(
    normalized_original_names %in% c("gender", "sex") & has_clean_binary_contrast
  )
  age_threshold_index <- which(
    (
      grepl("age.*threshold|age_above_threshold", normalized_original_names) |
        grepl("threshold", tolower(covariate_map$encoding))
    ) & has_clean_binary_contrast
  )
  gender_legacy <- legacy_covariate_alias_long(gender_index)
  age_threshold_legacy <- legacy_covariate_alias_long(age_threshold_index)

  # ============================================================
  # RETURN
  # ============================================================

  result <- list(
    # Core population parameters
    population = population_summary,
    beta_treatment = beta_treatment_summary,
    tau_treatment = tau_treatment_summary,
    arm_baseline_offset = arm_baseline_offset_summary,

    # Population trajectories
    population_model_location = population_model_location,
    population_median = population_median,
    population_time_means = population_time_means,
    population_change = population_change,
    population_directional_change = directional_population_change,
    population_standardized_change = population_standardized_change,
    population_ratio_from_baseline = population_ratio,
    population_percent_change_from_baseline = population_percent_change,

    # ----------------------------------------------------------
    # Family-oriented interface
    # ----------------------------------------------------------
    change = population_pairwise_change,
    change_from_baseline = population_change,
    change_consecutive = population_consecutive_change,

    covariates = list(
      P = P,
      names = covariate_names,
      map = covariate_map,
      metadata = metadata_value_long("covariate_metadata"),
      distributions = covariate_distributions,
      reference_profile = population_reference_profile,
      parameters = list(
        baseline = covariate_baseline_summary,
        linear_time = beta_covariate_summary,
        rw1_scale = tau_covariate_summary
      ),
      parameter_summary = covariate_parameter_summary,
      model_scale = covariate_model_effect,
      level = covariate_level_difference,
      from_baseline = covariate_change_difference,
      directional_from_baseline = covariate_directional_change_difference,
      pairwise_change = covariate_pairwise_change,
      consecutive_change = covariate_consecutive_change
    ),

    treatment = list(
      level = treatment_level_difference,
      change = treatment_pairwise_change,
      from_baseline = treatment_effects,
      consecutive = treatment_consecutive_change,
      responder_uplift = treatment_responder_uplift,
      ratio_of_ratios = treatment_ratio_of_ratios
    ),

    # Primary generic covariate table and metadata-derived legacy aliases.
    covariate_effects = covariate_effects,
    covariate_parameters = covariate_parameter_summary,
    covariate_pairwise_change = covariate_pairwise_change,
    gender = if (!is.null(gender_legacy)) gender_legacy$family else NULL,
    age = if (!is.null(age_threshold_legacy)) age_threshold_legacy$family else NULL,
    gender_effects = if (!is.null(gender_legacy)) gender_legacy$effects else NULL,
    age_threshold_effects = if (!is.null(age_threshold_legacy)) {
      age_threshold_legacy$effects
    } else {
      NULL
    },
    population_consecutive_change = population_consecutive_change,
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
    individual_covariates = individual_covariates,
    individual_slope = data.frame(),

    # Responders
    responder_proportion = latent_responder,
    responders = responder_summary,

    # Treatment effects vs reference arm
    treatment_effects = treatment_effects,
    treatment_responder_uplift = treatment_responder_uplift,
    treatment_ratio_of_ratios = treatment_ratio_of_ratios,

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

  class(result) <- c("mira_summary_long", "list")

  if (verbose) print(result)
  invisible(result)
}


# ============================================================
# PRINT METHOD
# ============================================================

#' Print a complete MIRA Bayesian summary
#'
#' Prints the clinically and statistically relevant results from
#' [mira_summary_long()] and finishes with a guide to the detailed `$` components.
#' The report is printed automatically by `mira_summary_long()` when `verbose = TRUE`.
#'
#' @param x A `mira_summary_long` object.
#' @param digits Number of decimal places used for continuous quantities.
#' @param max_rows Maximum number of rows printed for any one table.
#' @param trajectories Print adjusted population trajectories by arm.
#' @param treatment Print treatment effects versus the reference arm.
#' @param covariates Print time-specific effects for selected encoded
#'   covariate terms.
#' @param responders Print responder summaries.
#' @param heterogeneity Print random-effect / heterogeneity parameters.
#' @param ppc Print posterior predictive checks.
#' @param diagnostics Print MCMC diagnostics.
#' @param ... Unused.
#'
#' @return Invisibly returns `x`.
#' @export
print.mira_summary_long <- function(
    x,
    digits = 3,
    max_rows = Inf,
    trajectories = TRUE,
    treatment = TRUE,
    covariates = TRUE,
    responders = TRUE,
    heterogeneity = TRUE,
    ppc = TRUE,
    diagnostics = TRUE,
    ...
) {

  if (!is.numeric(digits) || length(digits) != 1L || is.na(digits) || digits < 0) {
    stop("`digits` must be one integer >= 0.", call. = FALSE)
  }
  digits <- as.integer(digits)

  if (!is.numeric(max_rows) || length(max_rows) != 1L || is.na(max_rows) || max_rows <= 0) {
    stop("`max_rows` must be > 0 or Inf.", call. = FALSE)
  }

  flags <- c(
    trajectories = trajectories,
    treatment = treatment,
    covariates = covariates,
    responders = responders,
    heterogeneity = heterogeneity,
    ppc = ppc,
    diagnostics = diagnostics
  )
  if (anyNA(flags) || !all(vapply(as.list(flags), is.logical, logical(1L)))) {
    stop("Section switches must be TRUE or FALSE.", call. = FALSE)
  }

  line_long <- function(char = "-", n = 92L) {
    cat(strrep(char, n), "\n", sep = "")
  }

  section_long <- function(title) {
    cat("\n", title, "\n", sep = "")
    line_long()
  }

  fmt_num_long <- function(z, d = digits) {
    if (length(z) == 0L || is.na(z) || !is.finite(z)) return("NA")
    formatC(z, format = "f", digits = d)
  }

  fmt_prob_long <- function(z, d = 3L) {
    if (length(z) == 0L || is.na(z) || !is.finite(z)) return("NA")
    formatC(z, format = "f", digits = d)
  }

  fmt_bool_long <- function(z) {
    if (length(z) == 0L || is.na(z)) "NA" else if (isTRUE(z)) "OK" else "CHECK"
  }

  fmt_interval_long <- function(mean, lower, upper, d = digits) {
    paste0(fmt_num_long(mean, d), " [", fmt_num_long(lower, d), ", ", fmt_num_long(upper, d), "]")
  }

  limit_table_long <- function(z, label = "rows") {
    if (!is.data.frame(z) || nrow(z) == 0L) return(z)
    if (is.finite(max_rows) && nrow(z) > max_rows) {
      cat(sprintf(
        "Showing first %d of %d %s. Use print(x, max_rows = Inf) for all.\n",
        as.integer(max_rows), nrow(z), label
      ))
      return(utils::head(z, as.integer(max_rows)))
    }
    z
  }

  safe_df_long <- function(z) {
    is.data.frame(z) && nrow(z) > 0L
  }

  info <- x$model_information
  if (is.null(info)) info <- list()

  cred_pct <- if (!is.null(info$credible_level) && is.finite(info$credible_level)) {
    100 * info$credible_level
  } else {
    NA_real_
  }
  cri_label <- if (is.finite(cred_pct)) {
    paste0(formatC(cred_pct, format = "fg", digits = 4L), "% CrI")
  } else {
    "CrI"
  }

  cat("\n")
  line_long("=")
  cat("MIRA BAYESIAN LONGITUDINAL MODEL - COMPLETE REPORT\n")
  line_long("=")

  # ------------------------------------------------------------
  # MODEL OVERVIEW
  # ------------------------------------------------------------
  section_long("MODEL OVERVIEW")

  if (!is.null(info$outcome) || !is.null(info$likelihood)) {
    cat(sprintf(
      "Outcome: %s | Likelihood: %s | Link/modeling scale: %s\n",
      ifelse(is.null(info$outcome), "NA", as.character(info$outcome)),
      ifelse(is.null(info$likelihood), "NA", as.character(info$likelihood)),
      ifelse(is.null(info$modeling_scale), "NA", as.character(info$modeling_scale))
    ))
  }

  if (!is.null(info$outcome_bounds) && any(is.finite(info$outcome_bounds))) {
    bound_text <- paste0(
      ifelse(is.finite(info$outcome_bounds[[1L]]), info$outcome_bounds[[1L]], "-Inf"),
      " to ",
      ifelse(is.finite(info$outcome_bounds[[2L]]), info$outcome_bounds[[2L]], "Inf")
    )
    cat("Observable support: ", bound_text,
        " | Boundary handling: ", info$boundary_strategy, "\n", sep = "")
  }

  cat(sprintf(
    "Subjects: %s | Observations: %s | Timepoints: %s | Arms: %s\n",
    ifelse(is.null(info$n_subjects), "NA", as.character(info$n_subjects)),
    ifelse(is.null(info$n_observations), "NA", as.character(info$n_observations)),
    ifelse(is.null(info$n_time_points), "NA", as.character(info$n_time_points)),
    ifelse(is.null(info$n_arms), "NA", as.character(info$n_arms))
  ))

  if (!is.null(info$arm_labels)) {
    cat(sprintf(
      "Arms: %s | Reference: %s\n",
      paste(info$arm_labels, collapse = ", "),
      ifelse(is.null(info$reference_label), "NA", info$reference_label)
    ))
  }

  if (!is.null(info$time_value)) {
    cat("Measurement times: ", paste(info$time_value, collapse = ", "), "\n", sep = "")
  }

  cat(sprintf(
    "Improvement direction: %s\n",
    ifelse(is.null(info$direction_interpretation), "NA", info$direction_interpretation)
  ))

  if (!is.null(info$population_reference_profile)) {
    cat("Population reference profile: ", info$population_reference_profile, "\n", sep = "")
  }

  if (!is.null(info$n_covariates) && info$n_covariates == 0L) {
    cat("Selected covariates: none (P = 0).\n")
  } else if (is.data.frame(info$covariate_map) && nrow(info$covariate_map) > 0L) {
    cat("Selected encoded covariate terms: ", nrow(info$covariate_map), "\n", sep = "")
    keep <- intersect(
      c(
        "index", "name", "label", "original_name", "type", "encoding",
        "level", "reference_level", "center", "scale", "unit"
      ),
      names(info$covariate_map)
    )
    print(info$covariate_map[, keep, drop = FALSE], row.names = FALSE)
  }

  if (safe_df_long(x$mcid)) {
    m <- x$mcid[1L, ]
    cat(sprintf(
      "MCID posterior: %s %s | Prior mean: %s | Prior SD: %s\n",
      fmt_interval_long(m$mean, m$lower, m$upper),
      cri_label,
      fmt_num_long(info$mcid_prior_mean),
      fmt_num_long(info$mcid_prior_sd)
    ))
  }

  if (!is.null(info$meaningful_between_arm_difference)) {
    cat("Clinically meaningful between-arm threshold: ",
        fmt_num_long(info$meaningful_between_arm_difference), "\n", sep = "")
  }

  # ------------------------------------------------------------
  # MCMC DIAGNOSTICS
  # ------------------------------------------------------------
  if (diagnostics) {
    section_long("MCMC DIAGNOSTICS")

    d <- x$diagnostics
    q <- x$quality_flags

    if (is.null(d)) {
      cat("Diagnostics not available.\n")
    } else if (!is.null(d$error)) {
      cat("Diagnostics error: ", d$error, "\n", sep = "")
    } else {
      cat(sprintf(
        "Max R-hat: %s | Min ESS bulk: %s | Min ESS tail: %s\n",
        fmt_num_long(d$max_Rhat, 4L), fmt_num_long(d$min_ESS_bulk, 0L), fmt_num_long(d$min_ESS_tail, 0L)
      ))
      if (!is.null(q)) {
        cat(sprintf(
          "R-hat < 1.01: %s | ESS bulk >= 400: %s | ESS tail >= 400: %s\n",
          fmt_bool_long(q$Rhat_ok), fmt_bool_long(q$ESS_bulk_ok), fmt_bool_long(q$ESS_tail_ok)
        ))
      }

      if (safe_df_long(d$parameters)) {
        bad <- d$parameters[
          !(d$parameters$Rhat_ok & d$parameters$ESS_bulk_ok & d$parameters$ESS_tail_ok),
          ,
          drop = FALSE
        ]
        if (nrow(bad) > 0L) {
          cat(sprintf("Parameters requiring attention: %d\n", nrow(bad)))
        } else {
          cat("All summarized parameters pass the stored R-hat / ESS thresholds.\n")
        }
      }
    }
  }

  # ------------------------------------------------------------
  # CORE PARAMETERS
  # ------------------------------------------------------------
  if (safe_df_long(x$population)) {
    section_long("CORE POSTERIOR PARAMETERS")

    core_names <- c(
      "baseline_mean", "beta_time", "sigma_intercept", "sigma_slope",
      "rho_subject", "sigma", "residual_sd", "residual_cv", "nu_value", "mcid"
    )
    core <- x$population[x$population$parameter %in% core_names, , drop = FALSE]
    if (nrow(core) > 0L) {
      core_print <- data.frame(
        Parameter = core$parameter,
        Mean = round(core$mean, digits),
        Median = round(core$median, digits),
        SD = round(core$sd, digits),
        CrI_low = round(core$lower, digits),
        CrI_high = round(core$upper, digits),
        check.names = FALSE
      )
      print(core_print, row.names = FALSE)
      cat("Intervals: ", cri_label, "\n", sep = "")
    }
  }

  # ------------------------------------------------------------
  # POPULATION TRAJECTORIES
  # ------------------------------------------------------------
  if (trajectories) {
    section_long("ADJUSTED POPULATION TRAJECTORIES")

    pt <- x$population_time_means
    if (!safe_df_long(pt)) {
      cat("Population trajectories not available.\n")
    } else {
      pt_print <- data.frame(
        Arm = pt$arm_label,
        Time = pt$time_value,
        Mean = round(pt$mean, digits),
        CrI_low = round(pt$lower, digits),
        CrI_high = round(pt$upper, digits),
        check.names = FALSE
      )
      print(limit_table_long(pt_print, "trajectory rows"), row.names = FALSE)
      if (!is.null(info$population_reference_profile)) {
        cat("Covariate profile: ", info$population_reference_profile, "\n", sep = "")
      }
    }
  }

  # ------------------------------------------------------------
  # ALL PAIRWISE POPULATION CHANGES
  # ------------------------------------------------------------
  section_long("ALL PAIRWISE POPULATION CHANGES")

  pc <- x$change
  if (!safe_df_long(pc)) {
    cat("Pairwise population changes not available.\n")
  } else {
    pc_print <- data.frame(
      Arm = pc$arm_label,
      From = pc$from_label,
      To = pc$to_label,
      From_time = pc$from_time_value,
      To_time = pc$to_time_value,
      Change = round(pc$mean, digits),
      CrI_low = round(pc$lower, digits),
      CrI_high = round(pc$upper, digits),
      Directional_change = round(pc$mean_directional_change, digits),
      P_improvement = round(pc$P_improvement, 3L),
      P_MCID = round(pc$P_MCID, 3L),
      check.names = FALSE
    )
    print(limit_table_long(pc_print, "pairwise population-change rows"), row.names = FALSE)
    cat("All unique forward time pairs are shown; reverse pairs are identical with opposite sign.\n")
  }

  # ------------------------------------------------------------
  # TREATMENT EFFECTS
  # ------------------------------------------------------------
  if (treatment) {
    section_long("TREATMENT EFFECTS VS REFERENCE ARM")

    te <- x$treatment_effects
    if (!safe_df_long(te)) {
      cat("Treatment contrasts not available.\n")
    } else {
      te_print <- data.frame(
        Treatment = te$treatment_label,
        Reference = te$reference_label,
        Time = te$time_value,
        Change_difference = round(te$mean, digits),
        CrI_low = round(te$lower, digits),
        CrI_high = round(te$upper, digits),
        Directional_benefit = round(te$mean_directional_benefit, digits),
        P_benefit = round(te$P_benefit_positive, 3L),
        P_meaningful = round(te$P_benefit_meaningful, 3L),
        Responder_uplift = round(te$mean_responder_probability_difference, 3L),
        check.names = FALSE
      )
      if (identical(info$likelihood, "lognormal") &&
          "mean_ratio_of_ratios" %in% names(te)) {
        te_print$Ratio_of_ratios <- round(te$mean_ratio_of_ratios, digits)
      }
      print(limit_table_long(te_print, "treatment-effect rows"), row.names = FALSE)
      cat("Change_difference is treatment minus reference in change from baseline.\n")
      cat("Directional_benefit > 0 favors treatment according to the declared improvement direction.\n")
      if (identical(info$likelihood, "lognormal")) {
        cat("Ratio_of_ratios is the treatment/reference ratio of temporal ratios.\n")
      }
    }

    tp <- if (!is.null(x$treatment)) x$treatment$change else NULL
    if (safe_df_long(tp)) {
      cat("\nAll pairwise treatment differences in temporal change:\n")
      tp_print <- data.frame(
        Treatment = tp$treatment_label,
        Reference = tp$reference_label,
        From = tp$from_label,
        To = tp$to_label,
        Change_difference = round(tp$mean, digits),
        CrI_low = round(tp$lower, digits),
        CrI_high = round(tp$upper, digits),
        Directional_benefit = round(tp$mean_directional_benefit, digits),
        P_benefit = round(tp$P_benefit_positive, 3L),
        P_meaningful = round(tp$P_benefit_meaningful, 3L),
        check.names = FALSE
      )
      print(limit_table_long(tp_print, "pairwise treatment-change rows"), row.names = FALSE)
    }
  }

  # ------------------------------------------------------------
  # SELECTED COVARIATE EFFECTS
  # ------------------------------------------------------------
  if (covariates) {
    section_long("SELECTED COVARIATE EFFECTS")

    ce <- x$covariate_effects
    if (!safe_df_long(ce)) {
      if (!is.null(info$P) && info$P == 0L) {
        cat("No additional covariates were selected (P = 0).\n")
      } else {
        cat("Covariate effects are not available.\n")
      }
    } else {
      value_or_na_long <- function(data, name) {
        if (name %in% names(data)) data[[name]] else rep(NA_real_, nrow(data))
      }
      ce_print <- data.frame(
        Encoded_term = ce$encoded_term,
        Covariate = ce$original_covariate,
        Interpretation = ce$effect_interpretation,
        Time = ce$time_value,
        Natural_difference = round(ce$mean, digits),
        CrI_low = round(ce$lower, digits),
        CrI_high = round(ce$upper, digits),
        Model_scale_effect = round(value_or_na_long(ce, "model_scale_mean"), digits),
        Change_from_baseline = round(
          value_or_na_long(ce, "change_from_baseline_mean"), digits
        ),
        Directional_change = round(
          value_or_na_long(ce, "directional_change_from_baseline_mean"), digits
        ),
        P_directional_change_gt_0 = round(
          value_or_na_long(ce, "P_directional_change_gt_0"), 3L
        ),
        check.names = FALSE
      )
      print(limit_table_long(ce_print, "covariate-effect rows"), row.names = FALSE)
      cat(
        "Natural differences compare X[p] = 1 with X[p] = 0 at the reference ",
        "arm while all other encoded covariates remain zero.\n",
        sep = ""
      )
    }

    cp <- x$covariate_parameters
    if (safe_df_long(cp)) {
      cat("\nCovariate trajectory parameters on the model/link scale:\n")
      cp_print <- data.frame(
        Encoded_term = cp$encoded_term,
        Component = cp$component,
        Mean = round(cp$mean, digits),
        SD = round(cp$sd, digits),
        CrI_low = round(cp$lower, digits),
        CrI_high = round(cp$upper, digits),
        check.names = FALSE
      )
      print(limit_table_long(cp_print, "covariate-parameter rows"), row.names = FALSE)
    }

    cc <- x$covariate_pairwise_change
    if (safe_df_long(cc)) {
      cat("\nAll pairwise changes in covariate effects:\n")
      cc_print <- data.frame(
        Encoded_term = cc$encoded_term,
        From = cc$from_label,
        To = cc$to_label,
        Change_difference = round(cc$mean, digits),
        CrI_low = round(cc$lower, digits),
        CrI_high = round(cc$upper, digits),
        Directional_difference = round(
          cc$mean_directional_change_difference, digits
        ),
        P_directional_gt_0 = round(cc$P_directional_change_gt_0, 3L),
        check.names = FALSE
      )
      print(limit_table_long(cc_print, "pairwise covariate-change rows"), row.names = FALSE)
    }

    cd <- if (!is.null(x$covariates)) x$covariates$distributions else NULL
    if (safe_df_long(cd)) {
      cat("\nEncoded covariate distributions across subjects:\n")
      keep <- intersect(
        c("name", "n", "n_unique", "mean", "sd", "min", "median", "max", "n_zero", "n_one"),
        names(cd)
      )
      cd_print <- cd[, keep, drop = FALSE]
      numeric_cols <- vapply(cd_print, is.numeric, logical(1L))
      cd_print[numeric_cols] <- lapply(cd_print[numeric_cols], round, digits = digits)
      print(limit_table_long(cd_print, "covariate-distribution rows"), row.names = FALSE)
    }
  }

  # ------------------------------------------------------------
  # CLINICAL / RESPONDER INFORMATION
  # ------------------------------------------------------------
  if (responders) {
    section_long("CLINICAL AND RESPONDER SUMMARY")

    cl <- x$clinical
    if (safe_df_long(cl)) {
      cl_print <- data.frame(
        Arm = cl$arm_label,
        Final_time = cl$final_time_value,
        Mean_change = round(cl$final_population_change_mean, digits),
        P_improvement = round(cl$P_final_improvement, 3L),
        P_responder = round(cl$P_final_responder, 3L),
        Directional_change_minus_MCID = round(cl$mean_final_directional_change_minus_mcid, digits),
        check.names = FALSE
      )
      cat("Final-time population summary:\n")
      print(cl_print, row.names = FALSE)
    }

    rp <- x$responder_proportion
    if (safe_df_long(rp)) {
      rp_print <- data.frame(
        Arm = rp$arm_label,
        Time = rp$time_value,
        Latent_new_subject_P_responder = round(rp$mean, 3L),
        CrI_low = round(rp$lower, 3L),
        CrI_high = round(rp$upper, 3L),
        check.names = FALSE
      )
      cat("\nNew-subject latent responder probability:\n")
      print(limit_table_long(rp_print, "responder rows"), row.names = FALSE)
    }

    er <- x$responders
    if (safe_df_long(er)) {
      keep <- c("arm_label", "time_value", "expected_responder_proportion")
      extra <- grep("^proportion_P_MCID_ge_", names(er), value = TRUE)
      er_print <- er[, c(keep, extra), drop = FALSE]
      names(er_print)[1:3] <- c("Arm", "Time", "Expected_responder_proportion")
      numeric_cols <- vapply(er_print, is.numeric, logical(1L))
      er_print[numeric_cols] <- lapply(er_print[numeric_cols], round, digits = 3L)
      cat("\nExisting-subject responder overview:\n")
      print(limit_table_long(er_print, "existing-subject responder rows"), row.names = FALSE)
    }
  }

  # ------------------------------------------------------------
  # HETEROGENEITY
  # ------------------------------------------------------------
  if (heterogeneity) {
    section_long("HETEROGENEITY / VARIANCE COMPONENTS")

    h <- x$heterogeneity
    if (!safe_df_long(h)) {
      cat("Heterogeneity summaries not available.\n")
    } else {
      h_print <- data.frame(
        Parameter = h$parameter,
        Mean = round(h$mean, digits),
        SD = round(h$sd, digits),
        CrI_low = round(h$lower, digits),
        CrI_high = round(h$upper, digits),
        check.names = FALSE
      )
      print(h_print, row.names = FALSE)
    }
  }

  # ------------------------------------------------------------
  # POSTERIOR PREDICTIVE CHECKS
  # ------------------------------------------------------------
  if (ppc) {
    section_long("POSTERIOR PREDICTIVE CHECKS")

    z <- x$ppc
    if (is.null(z)) {
      cat("Posterior predictive checks not available.\n")
    } else {
      stat_names <- intersect(
        c(
          "mean", "sd", "median", "q05", "q95", "min", "max",
          "skewness", "coefficient_of_variation"
        ),
        names(z$posterior_predictive)
      )

      if (length(stat_names) > 0L) {
        rows <- lapply(stat_names, function(nm) {
          pp <- z$posterior_predictive[[nm]]
          data.frame(
            Statistic = nm,
            Observed = as.numeric(z$observed[[nm]]),
            Predictive_mean = as.numeric(pp[["mean"]]),
            Predictive_CrI_low = as.numeric(pp[["lower"]]),
            Predictive_CrI_high = as.numeric(pp[["upper"]]),
            stringsAsFactors = FALSE
          )
        })
        pp_print <- do.call(rbind, rows)
        numeric_cols <- vapply(pp_print, is.numeric, logical(1L))
        pp_print[numeric_cols] <- lapply(pp_print[numeric_cols], round, digits = digits)
        print(pp_print, row.names = FALSE)
      }

      if (!is.null(z$bayesian_p_values)) {
        cat("Bayesian predictive p-values: ")
        cat(
          paste0(
            names(z$bayesian_p_values), "=",
            vapply(z$bayesian_p_values, fmt_prob_long, character(1L)),
            collapse = " | "
          ),
          "\n"
        )
      }

      if (!is.null(z$family_checks)) {
        fc <- z$family_checks
        if (!is.null(fc$predictive_support_violation_fraction)) {
          cat(
            "Predictive support-violation fraction: ",
            fmt_interval_long(
              fc$predictive_support_violation_fraction[["mean"]],
              fc$predictive_support_violation_fraction[["lower"]],
              fc$predictive_support_violation_fraction[["upper"]],
              4L
            ),
            "\n",
            sep = ""
          )
        }
        if (!is.null(fc$predictive_lower_boundary_fraction)) {
          cat(
            "Lower-bound mass (observed / predictive mean): ",
            fmt_num_long(fc$observed_lower_boundary_fraction, 4L), " / ",
            fmt_num_long(fc$predictive_lower_boundary_fraction[["mean"]], 4L),
            "\n",
            sep = ""
          )
        }
        if (!is.null(fc$predictive_upper_boundary_fraction)) {
          cat(
            "Upper-bound mass (observed / predictive mean): ",
            fmt_num_long(fc$observed_upper_boundary_fraction, 4L), " / ",
            fmt_num_long(fc$predictive_upper_boundary_fraction[["mean"]], 4L),
            "\n",
            sep = ""
          )
        }
        if (!is.null(fc$predictive_longitudinal_lag1_correlation)) {
          cat(
            "Longitudinal lag-1 correlation (observed / predictive mean): ",
            fmt_num_long(fc$observed_longitudinal_lag1_correlation, 3L), " / ",
            fmt_num_long(
              fc$predictive_longitudinal_lag1_correlation[["mean"]],
              3L
            ),
            "\n",
            sep = ""
          )
        }
      }
      cat("Values near 0.5 indicate good centering; values near 0 or 1 can indicate lack of fit.\n")
    }
  }

  # ------------------------------------------------------------
  # COMPLETE OUTPUT GUIDE
  # ------------------------------------------------------------
  section_long("COMPLETE OUTPUT GUIDE")

  cat("Use the following `$` components on the object returned by mira_summary_long().\n")
  cat("For example, if you used `res <- mira_summary_long(...)`, read `$treatment_effects` as `res$treatment_effects`.\n\n")

  guide <- c(
    "$change" = "ALL unique population time-to-time changes for every arm (t0->t1, t0->t2, ...)",
    "$change_from_baseline" = "Population changes restricted to baseline -> each follow-up",
    "$change_consecutive" = "Population changes restricted to consecutive visits",
    "$covariate_effects" = "Readable term-by-time natural/model-scale effects and changes from baseline",
    "$covariate_parameters" = "Baseline, linear-time and RW1-scale parameters for every encoded term",
    "$covariate_pairwise_change" = "ALL pairwise changes in each encoded covariate effect",
    "$covariates$map" = "Mapping from encoded terms to original variables, levels, references and scaling",
    "$covariates$distributions" = "Subject-level distributions of encoded design-matrix columns",
    "$covariates$model_scale" = "Time-specific +1-design-unit effects on the model/link scale",
    "$covariates$level" = "Natural-scale +1-design-unit mean differences at each time",
    "$covariates$from_baseline" = "Covariate differences in change from baseline",
    "$covariates$directional_from_baseline" = "Clinically oriented covariate differences in change",
    "$treatment$level" = "Treatment - reference adjusted level difference at every time",
    "$treatment$change" = "ALL pairwise treatment-vs-reference differences in temporal change",
    "$treatment$from_baseline" = "Treatment-vs-reference differences in change from baseline",
    "$treatment$consecutive" = "Treatment-vs-reference differences in consecutive-visit changes",
    "$model_information" = "Outcome family/link/support plus study setup, direction, MCID and group counts",
    "$population" = "Core posterior parameters with posterior mean, SD and credible intervals",
    "$beta_treatment" = "Treatment slope parameters versus the reference arm",
    "$tau_treatment" = "Treatment-specific RW1 trajectory scales",
    "$arm_baseline_offset" = "Posterior baseline imbalance between treatment arms",
    "$population_model_location" = "Adjusted trajectory on the active identity/log modeling scale",
    "$population_median" = "Natural-scale conditional median/location for each arm and time",
    "$population_time_means" = "Natural-scale adjusted mean/location trajectory for each arm and time",
    "$population_change" = "Absolute natural-scale adjusted change from baseline",
    "$population_directional_change" = "Change from baseline oriented so positive values mean clinical improvement",
    "$population_standardized_change" = "Standardized change on the family-appropriate modeling scale",
    "$population_ratio_from_baseline" = "CMT natural-scale ratio to baseline (one for identity models)",
    "$population_percent_change_from_baseline" = "CMT percent change from baseline (zero for identity models)",
    "$treatment_effects" = "Primary treatment-versus-reference change contrasts, probabilities of benefit and responder uplift",
    "$treatment_responder_uplift" = "Difference between arms in latent responder probability",
    "$treatment_ratio_of_ratios" = "CMT treatment/reference ratio of temporal ratios",
    "$clinical" = "Final-time population change, probability of improvement and probability of MCID response",
    "$population_clinical$latent_any_improvement" = "Latent probability of any improvement for a new subject",
    "$population_clinical$latent_responder" = "Latent probability of clinically meaningful response for a new subject",
    "$population_clinical$predictive_responder_draw" = "Posterior-predictive responder draws for a new subject",
    "$individual_change" = "Subject-specific posterior change from baseline",
    "$individual_directional_change" = "Subject-specific change oriented according to the clinical direction",
    "$individual_clinical" = "Subject-level posterior responder probabilities and responder classifications",
    "$individual_covariates" = "Subject lookup table containing the active encoded design-matrix values",
    "$responders" = "Existing-subject responder overview by arm and time",
    "$heterogeneity" = "Random-intercept, random-slope, correlation and trajectory-heterogeneity summaries",
    "$ppc" = "Family-specific PPCs, boundary/support checks, visit checks and predictive p-values",
    "$loo" = "Pointwise log-likelihood information for model comparison / LOO workflows",
    "$diagnostics" = "MCMC diagnostic summary",
    "$diagnostics$parameters" = "R-hat, bulk ESS and tail ESS for every monitored parameter",
    "$quality_flags" = "Quick diagnostic pass/fail flags",
    "$variable_inventory" = "Names of Stan generated quantities grouped by purpose",
    "$draws" = "Complete posterior draws; use only when a custom posterior calculation is needed"
  )

  if (!is.null(x$gender_effects)) {
    guide <- c(
      guide,
      "$gender_effects" = "Metadata-derived legacy alias for a single binary gender/sex term"
    )
  }
  if (!is.null(x$age_threshold_effects)) {
    guide <- c(
      guide,
      "$age_threshold_effects" = "Metadata-derived legacy alias for an explicit age-threshold term"
    )
  }

  width <- max(nchar(names(guide)))
  for (nm in names(guide)) {
    cat(sprintf("  %-*s  %s\n", width, nm, unname(guide[[nm]])))
  }

  cat("\nUseful commands:\n")
  cat("  names(object)                         list all top-level components\n")
  cat("  object$change                         inspect every population time-to-time change\n")
  cat("  object$treatment$change               inspect every pairwise treatment change contrast\n")
  cat("  object$covariate_effects               inspect all term-by-time covariate effects\n")
  cat("  object$covariates$map                 inspect encoding, references and scaling\n")
  cat("  object$covariate_pairwise_change       inspect pairwise covariate-effect changes\n")
  cat("  object$treatment_effects               primary baseline treatment contrasts\n")
  cat("  object$diagnostics$parameters         inspect detailed MCMC diagnostics\n")
  cat("  print(object, max_rows = Inf)         print every row of the report tables\n")
  cat("  print(object, ppc = FALSE)            skip PPC output when a shorter print is desired\n")

  line_long("=")

  invisible(x)
}
