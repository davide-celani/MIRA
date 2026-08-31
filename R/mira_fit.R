#' Fit MIRA longitudinal treatment model
#'
#' Fits the MIRA longitudinal Student-t mixed-effects model with
#' treatment-, gender-, and age-threshold-specific trajectories using CmdStan.
#'
#' @param stan_data Data prepared for the MIRA Stan model. The list may
#'   contain additional R-side metadata (for example `mean_y`, `sd_y`, or
#'   `arm_labels`, gender labels, or the age threshold); only variables required
#'   by Stan are passed to CmdStan. The current model requires subject-level
#'   binary indicators `male` (0 = Female reference, 1 = Male) and
#'   `age_above_threshold` (0 = age <= threshold, 1 = age > threshold).
#' @param prior A `mira_prior` object, or a named list containing the Stan
#'   prior fields required by the current model, including dedicated priors
#'   for gender- and age-threshold trajectories. If `NULL`, `mira_prior()`
#'   is called with profile = "default", unless all prior fields are already
#'   present in `stan_data`.
#' @param chains Number of MCMC chains.
#' @param parallel_chains Number of parallel chains.
#' @param iter_warmup Number of warmup iterations.
#' @param iter_sampling Number of sampling iterations.
#' @param seed Random seed.
#' @param refresh Number of iterations between progress messages.
#' @param stan_file Optional path to the Stan file. If `NULL`, MIRA first
#'   looks for `inst/stan/gaussian_longitudinal_gender_age.stan`, then the
#'   legacy MIRA Stan filenames for backwards compatibility.
#'
#' @return A CmdStanMCMC object.
#'
#' @export
mira_fit <- function(
    stan_data,
    prior = NULL,
    chains = 4,
    parallel_chains = chains,
    iter_warmup = 1000,
    iter_sampling = 3000,
    seed = 123,
    refresh = 100,
    stan_file = NULL
) {

  # ------------------------------------------------------------
  # Data object
  # ------------------------------------------------------------

  if (!is.list(stan_data)) {
    stop("`stan_data` must be a list.", call. = FALSE)
  }

  model_data_names <- c(
    "N", "S", "K", "G",
    "y", "subject", "time", "arm",
    "male", "age_above_threshold", "time_value",
    "direction",
    "mcid_prior_mean", "mcid_prior_sd",
    "meaningful_between_arm_difference"
  )

  missing_data <- setdiff(model_data_names, names(stan_data))

  if (length(missing_data) > 0) {
    stop(
      "Missing data required by the new MIRA Stan model: ",
      paste(missing_data, collapse = ", "),
      call. = FALSE
    )
  }

  # ------------------------------------------------------------
  # Dimensions and indices
  # ------------------------------------------------------------

  scalar_integer_names <- c("N", "S", "K", "G")

  for (nm in scalar_integer_names) {
    x <- stan_data[[nm]]
    if (length(x) != 1 || !is.numeric(x) || !is.finite(x) || x != as.integer(x)) {
      stop("`", nm, "` must be one finite integer.", call. = FALSE)
    }
  }

  if (stan_data$N < 1) stop("`N` must be >= 1.", call. = FALSE)
  if (stan_data$S < 1) stop("`S` must be >= 1.", call. = FALSE)
  if (stan_data$K < 2) stop("`K` must be >= 2.", call. = FALSE)
  if (stan_data$G < 2) stop("`G` must be >= 2 for the current treatment model.", call. = FALSE)

  if (!is.numeric(stan_data$y) ||
      length(stan_data$y) != stan_data$N ||
      any(!is.finite(stan_data$y))) {
    stop("`y` must contain exactly N finite numeric values.", call. = FALSE)
  }

  if (length(stan_data$subject) != stan_data$N) {
    stop("`length(subject)` must equal `N`.", call. = FALSE)
  }

  if (length(stan_data$time) != stan_data$N) {
    stop("`length(time)` must equal `N`.", call. = FALSE)
  }

  if (length(stan_data$arm) != stan_data$S) {
    stop("`length(arm)` must equal `S`.", call. = FALSE)
  }

  if (any(!is.finite(stan_data$subject)) ||
      any(stan_data$subject != as.integer(stan_data$subject)) ||
      any(stan_data$subject < 1) ||
      any(stan_data$subject > stan_data$S)) {
    stop("`subject` must contain integers between 1 and S.", call. = FALSE)
  }

  if (any(!is.finite(stan_data$time)) ||
      any(stan_data$time != as.integer(stan_data$time)) ||
      any(stan_data$time < 1) ||
      any(stan_data$time > stan_data$K)) {
    stop("`time` must contain integers between 1 and K.", call. = FALSE)
  }

  if (any(!is.finite(stan_data$arm)) ||
      any(stan_data$arm != as.integer(stan_data$arm)) ||
      any(stan_data$arm < 1) ||
      any(stan_data$arm > stan_data$G)) {
    stop("`arm` must contain integers between 1 and G.", call. = FALSE)
  }

  # ------------------------------------------------------------
  # Subject-level gender and age-group indicators
  # ------------------------------------------------------------

  validate_binary_subject_indicator <- function(x, name) {
    if (!is.numeric(x) || length(x) != stan_data$S) {
      stop(
        "`", name, "` must be a numeric/integer vector of length S.",
        call. = FALSE
      )
    }

    if (any(!is.finite(x)) ||
        any(x != as.integer(x)) ||
        any(!x %in% c(0, 1))) {
      stop(
        "`", name, "` must contain exactly S binary integer values (0/1).",
        call. = FALSE
      )
    }
  }

  validate_binary_subject_indicator(stan_data$male, "male")
  validate_binary_subject_indicator(
    stan_data$age_above_threshold,
    "age_above_threshold"
  )

  if (length(unique(stan_data$male)) < 2) {
    warning(
      "`male` contains only one observed category; the gender-by-time effect ",
      "will be weakly/non-identified by these data.",
      call. = FALSE
    )
  }

  if (length(unique(stan_data$age_above_threshold)) < 2) {
    warning(
      "`age_above_threshold` contains only one observed group; the age-group-by-time ",
      "effect will be weakly/non-identified by these data.",
      call. = FALSE
    )
  }

  if (!is.numeric(stan_data$time_value) ||
      length(stan_data$time_value) != stan_data$K ||
      any(!is.finite(stan_data$time_value)) ||
      is.unsorted(stan_data$time_value, strictly = TRUE)) {
    stop(
      "`time_value` must contain exactly K finite, strictly increasing values.",
      call. = FALSE
    )
  }

  # ------------------------------------------------------------
  # Clinical inputs
  # ------------------------------------------------------------

  if (length(stan_data$direction) != 1 ||
      !is.numeric(stan_data$direction) ||
      !is.finite(stan_data$direction) ||
      !(stan_data$direction %in% c(-1, 1))) {
    stop("`direction` must be exactly +1 or -1.", call. = FALSE)
  }

  positive_scalar <- function(x, name, allow_zero = FALSE) {
    ok <- length(x) == 1 && is.numeric(x) && is.finite(x)
    ok <- ok && if (allow_zero) x >= 0 else x > 0
    if (!ok) {
      comparator <- if (allow_zero) ">= 0" else "> 0"
      stop("`", name, "` must be one finite numeric value ", comparator, ".", call. = FALSE)
    }
  }

  positive_scalar(stan_data$mcid_prior_mean, "mcid_prior_mean", allow_zero = TRUE)
  positive_scalar(stan_data$mcid_prior_sd, "mcid_prior_sd")
  positive_scalar(
    stan_data$meaningful_between_arm_difference,
    "meaningful_between_arm_difference",
    allow_zero = TRUE
  )

  # ------------------------------------------------------------
  # Prior data required by Stan
  # ------------------------------------------------------------

  prior_names <- c(
    "baseline_prior_mean",
    "baseline_prior_sd",
    "beta_time_prior_mean",
    "beta_time_prior_sd",
    "tau_common_prior_rate",
    "beta_treatment_prior_sd",
    "tau_treatment_prior_rate",
    "arm_baseline_sd_prior_rate",
    "gender_baseline_prior_sd",
    "beta_gender_prior_sd",
    "tau_gender_prior_rate",
    "age_baseline_prior_sd",
    "beta_age_prior_sd",
    "tau_age_prior_rate",
    "sigma_intercept_prior_rate",
    "sigma_slope_prior_rate",
    "sigma_prior_rate",
    "nu_prior_shape",
    "nu_prior_rate"
  )

  # Allow fully assembled Stan data as an advanced use case.
  if (all(prior_names %in% names(stan_data))) {

    stan_prior_data <- stan_data[prior_names]

  } else {

    if (is.null(prior)) {
      prior <- mira_prior(
        stan_data,
        profile = "default"
      )
    }

    # A direct named list is useful for development/testing of the new model.
    if (is.list(prior) && all(prior_names %in% names(prior))) {
      stan_prior_data <- prior[prior_names]
    } else {
      mira_validate_prior(prior)
      stan_prior_data <- mira_prior_stan_data(prior)
    }
  }

  missing_prior <- setdiff(prior_names, names(stan_prior_data))

  if (length(missing_prior) > 0) {
    stop(
      "The prior specification is not compatible with the new Stan model. ",
      "Missing Stan prior fields: ",
      paste(missing_prior, collapse = ", "),
      ". Update `mira_prior()` / `mira_prior_stan_data()` so that gender and age priors are included, or pass a complete named prior list.",
      call. = FALSE
    )
  }

  stan_prior_data <- stan_prior_data[prior_names]

  finite_prior_names <- c("baseline_prior_mean", "beta_time_prior_mean")
  for (nm in finite_prior_names) {
    x <- stan_prior_data[[nm]]
    if (length(x) != 1 || !is.numeric(x) || !is.finite(x)) {
      stop("`", nm, "` must be one finite numeric value.", call. = FALSE)
    }
  }

  positive_prior_names <- setdiff(prior_names, finite_prior_names)
  for (nm in positive_prior_names) {
    positive_scalar(stan_prior_data[[nm]], nm)
  }

  # Pass only objects declared in the Stan data block. This lets stan_data
  # safely retain R-side metadata such as mean_y, sd_y and arm labels.
  sampling_data <- c(
    stan_data[model_data_names],
    stan_prior_data
  )

  sampling_data$N <- as.integer(sampling_data$N)
  sampling_data$S <- as.integer(sampling_data$S)
  sampling_data$K <- as.integer(sampling_data$K)
  sampling_data$G <- as.integer(sampling_data$G)
  sampling_data$subject <- as.integer(sampling_data$subject)
  sampling_data$time <- as.integer(sampling_data$time)
  sampling_data$arm <- as.integer(sampling_data$arm)
  sampling_data$male <- as.integer(sampling_data$male)
  sampling_data$age_above_threshold <- as.integer(
    sampling_data$age_above_threshold
  )
  sampling_data$direction <- as.integer(sampling_data$direction)

  # ------------------------------------------------------------
  # Stan model file
  # ------------------------------------------------------------

  if (is.null(stan_file)) {

    candidates <- c(
      "gaussian_longitudinal_gender_age.stan",
      "mira_longitudinal.stan",
      "gaussian_longitudinal.stan"
    )

    candidate_paths <- vapply(
      candidates,
      function(x) system.file("stan", x, package = "MIRA"),
      character(1)
    )

    existing <- candidate_paths[nzchar(candidate_paths)]

    # Also support sourcing the function from the package project before
    # installation, when system.file() cannot yet resolve the package path.
    if (length(existing) == 0) {
      local_candidates <- file.path("inst", "stan", candidates)
      existing <- local_candidates[file.exists(local_candidates)]
    }

    if (length(existing) == 0) {
      stop(
        "Could not find the MIRA Stan model in `inst/stan`. Expected ",
        "`gaussian_longitudinal_gender_age.stan` (preferred), ",
        "`mira_longitudinal.stan`, or `gaussian_longitudinal.stan`.",
        call. = FALSE
      )
    }

    stan_file <- existing[[1]]

  } else {

    if (length(stan_file) != 1 || !is.character(stan_file) || !nzchar(stan_file)) {
      stop("`stan_file` must be NULL or one non-empty character path.", call. = FALSE)
    }

    if (!file.exists(stan_file)) {
      packaged_file <- system.file("stan", stan_file, package = "MIRA")
      if (!nzchar(packaged_file)) {
        stop("Could not find Stan file: ", stan_file, call. = FALSE)
      }
      stan_file <- packaged_file
    }
  }

  message("Compiling MIRA Stan model: ", basename(stan_file))

  model <- cmdstanr::cmdstan_model(
    stan_file,
    quiet = TRUE
  )

  # ------------------------------------------------------------
  # Initial values
  # ------------------------------------------------------------

  y <- as.numeric(stan_data$y)
  baseline_y <- y[stan_data$time == 1]

  mean_y <- if (!is.null(stan_data$mean_y)) {
    as.numeric(stan_data$mean_y)[1]
  } else {
    mean(y)
  }

  sd_y <- if (!is.null(stan_data$sd_y)) {
    as.numeric(stan_data$sd_y)[1]
  } else {
    stats::sd(y)
  }

  if (!is.finite(mean_y)) mean_y <- mean(y)
  if (!is.finite(sd_y) || sd_y <= 0) sd_y <- max(abs(mean_y) * 0.1, 1)

  baseline_init <- if (length(baseline_y) > 0) mean(baseline_y) else mean_y
  elapsed <- max(stan_data$time_value) - min(stan_data$time_value)
  elapsed_safe <- max(elapsed, 1e-6)

  slope_scale <- max(sd_y / elapsed_safe, 0.01)
  rw_scale <- max(sd_y / sqrt(elapsed_safe) / 10, 0.01)

  init <- function() {
    list(
      baseline_mean = baseline_init,
      beta_time = 0,
      z_common_step = rep(0, stan_data$K - 1),
      tau_common = rw_scale,
      beta_treatment = rep(0, stan_data$G - 1),
      z_treatment_step = matrix(
        0,
        nrow = stan_data$G - 1,
        ncol = stan_data$K - 1
      ),
      tau_treatment = rep(rw_scale, stan_data$G - 1),
      z_arm_baseline = rep(0, stan_data$G - 1),
      arm_baseline_sd = max(sd_y / 10, 0.01),

      # Gender-by-time trajectory: Male - Female.
      gender_baseline_effect = 0,
      beta_gender_time = 0,
      z_gender_step = rep(0, stan_data$K - 1),
      tau_gender = rw_scale,

      # Age-group-by-time trajectory: age > threshold - age <= threshold.
      age_baseline_effect = 0,
      beta_age_time = 0,
      z_age_step = rep(0, stan_data$K - 1),
      tau_age = rw_scale,

      z_subject = matrix(
        0,
        nrow = 2,
        ncol = stan_data$S
      ),
      sigma_subject = c(
        max(sd_y / 2, 0.1),
        max(slope_scale / 2, 0.01)
      ),
      L_subject = diag(2),
      sigma = max(sd_y / 2, 0.1),
      nu = 10,
      mcid = max(stan_data$mcid_prior_mean, 1e-6)
    )
  }

  # ------------------------------------------------------------
  # Sampling
  # ------------------------------------------------------------

  message("Sampling posterior...")

  fit <- model$sample(
    data = sampling_data,
    chains = chains,
    parallel_chains = parallel_chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    seed = seed,
    init = init,
    refresh = refresh
  )

  return(fit)
}
