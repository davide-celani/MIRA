#' Create a MIRA prior specification
#'
#' Creates prior specifications for the MIRA longitudinal multi-arm
#' Student-t model with treatment-, gender-, and age-threshold-specific
#' longitudinal effects.
#'
#' Priors are scaled using the observed outcome SD and the total elapsed
#' follow-up time, so the specification remains applicable with any
#' number of measurement occasions K >= 2 and with irregular time spacing.
#'
#' @param stan_data Data prepared by [mira_prepare_data()].
#' @param profile Character string specifying a predefined prior profile.
#'   Available profiles are `"default"`, `"weak"`, `"regularized"`,
#'   `"informative"`, and `"custom"`.
#' @param baseline_mean Prior mean for the reference-arm baseline mean.
#' @param baseline_sd Prior SD for the reference-arm baseline mean.
#' @param beta_time_mean Prior mean for the reference-arm global time slope.
#' @param beta_time_sd Prior SD for the reference-arm global time slope.
#' @param tau_common_rate Exponential rate for continuous-time RW1
#'   deviations of the reference trajectory.
#' @param beta_treatment_sd Prior SD for treatment slopes relative to the
#'   reference arm.
#' @param tau_treatment_rate Exponential rate for treatment-specific RW1
#'   deviations.
#' @param arm_baseline_sd_rate Exponential rate for baseline imbalance SD
#'   between treatment arms.
#' @param gender_baseline_sd Prior SD for the Male - Female baseline
#'   difference. For predefined profiles, defaults to the same scale as
#'   `baseline_sd`; for `profile = "custom"`, NULL inherits `baseline_sd`.
#' @param beta_gender_sd Prior SD for the Male - Female global time-slope
#'   difference. For predefined profiles, defaults to the treatment-slope
#'   scale; for `profile = "custom"`, NULL inherits `beta_treatment_sd`.
#' @param tau_gender_rate Exponential rate for the gender-specific RW1 scale.
#'   For `profile = "custom"`, NULL inherits `tau_treatment_rate`.
#' @param age_baseline_sd Prior SD for the baseline difference between
#'   subjects above versus at/below the age threshold. For predefined profiles,
#'   defaults to the same scale as `baseline_sd`; for `profile = "custom"`,
#'   NULL inherits `baseline_sd`.
#' @param beta_age_sd Prior SD for the age-group global time-slope difference.
#'   For predefined profiles, defaults to the treatment-slope scale; for
#'   `profile = "custom"`, NULL inherits `beta_treatment_sd`.
#' @param tau_age_rate Exponential rate for the age-group-specific RW1 scale.
#'   For `profile = "custom"`, NULL inherits `tau_treatment_rate`.
#' @param sigma_intercept_rate Exponential rate for subject random-intercept SD.
#' @param sigma_slope_rate Exponential rate for subject random-slope SD.
#' @param sigma_rate Exponential rate for residual Student-t scale.
#' @param nu_shape Shape parameter for the Gamma prior on Student-t degrees
#'   of freedom.
#' @param nu_rate Rate parameter for the Gamma prior on Student-t degrees
#'   of freedom.
#'
#' @return An object of class `mira_prior`.
#'
#' @export
mira_prior <- function(
    stan_data,
    profile = c(
      "default",
      "weak",
      "regularized",
      "informative",
      "custom"
    ),
    baseline_mean = NULL,
    baseline_sd = NULL,
    beta_time_mean = NULL,
    beta_time_sd = NULL,
    tau_common_rate = NULL,
    beta_treatment_sd = NULL,
    tau_treatment_rate = NULL,
    arm_baseline_sd_rate = NULL,
    gender_baseline_sd = NULL,
    beta_gender_sd = NULL,
    tau_gender_rate = NULL,
    age_baseline_sd = NULL,
    beta_age_sd = NULL,
    tau_age_rate = NULL,
    sigma_intercept_rate = NULL,
    sigma_slope_rate = NULL,
    sigma_rate = NULL,
    nu_shape = NULL,
    nu_rate = NULL
) {

  # ============================================================
  # DATA VALIDATION
  # ============================================================

  if (!is.list(stan_data)) {
    stop("`stan_data` must be a list.", call. = FALSE)
  }

  required <- c(
    "y",
    "time",
    "time_value",
    "mean_y",
    "sd_y"
  )

  missing <- setdiff(required, names(stan_data))

  if (length(missing) > 0) {
    stop(
      "Missing data required to construct MIRA priors: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  profile <- match.arg(profile)

  y <- as.numeric(stan_data$y)
  time <- as.integer(stan_data$time)
  time_value <- as.numeric(stan_data$time_value)
  mean_y <- as.numeric(stan_data$mean_y)
  sd_y <- as.numeric(stan_data$sd_y)

  if (length(y) < 1 || any(!is.finite(y))) {
    stop("`stan_data$y` must contain finite numeric values.", call. = FALSE)
  }

  if (length(time) != length(y) ||
      any(!is.finite(time)) ||
      any(time < 1)) {
    stop(
      "`stan_data$time` must contain one valid time index per observation.",
      call. = FALSE
    )
  }

  if (length(time_value) < 2 ||
      any(!is.finite(time_value)) ||
      is.unsorted(time_value, strictly = TRUE)) {
    stop(
      "`stan_data$time_value` must contain at least two finite, strictly increasing values.",
      call. = FALSE
    )
  }

  if (length(mean_y) != 1 || !is.finite(mean_y)) {
    stop(
      "`stan_data$mean_y` must be a single finite numeric value.",
      call. = FALSE
    )
  }

  if (length(sd_y) != 1 || !is.finite(sd_y) || sd_y <= 0) {
    stop(
      "`stan_data$sd_y` must be a single positive finite numeric value.",
      call. = FALSE
    )
  }

  # ============================================================
  # DATA-ADAPTIVE REFERENCE SCALES
  # ============================================================

  baseline_y <- y[time == 1]

  if (length(baseline_y) < 1) {
    stop(
      "No baseline observations (`time == 1`) were found.",
      call. = FALSE
    )
  }

  baseline_center <- mean(baseline_y)

  time_span <- max(time_value) - min(time_value)

  if (!is.finite(time_span) || time_span <= 0) {
    stop(
      "`time_value` must span a positive amount of time.",
      call. = FALSE
    )
  }

  # beta_* and subject slope have units outcome / time.
  slope_scale <- sd_y / time_span

  # tau_* multiplies sqrt(dt), therefore has units outcome / sqrt(time).
  rw_scale <- sd_y / sqrt(time_span)

  # Numerical guards only; these do not materially alter ordinary data.
  outcome_scale <- max(sd_y, 1e-8)
  slope_scale <- max(slope_scale, 1e-8)
  rw_scale <- max(rw_scale, 1e-8)

  # ============================================================
  # PREDEFINED PROFILES
  # ============================================================

  if (profile == "default") {

    baseline_mean <- baseline_center
    baseline_sd <- 2 * outcome_scale

    beta_time_mean <- 0
    beta_time_sd <- 2 * slope_scale

    tau_common_rate <- 1 / rw_scale

    beta_treatment_sd <- 2 * slope_scale
    tau_treatment_rate <- 1 / rw_scale

    arm_baseline_sd_rate <- 1 / outcome_scale

    # Covariate effects are centered at zero. Dedicated scales default to
    # the analogous baseline/treatment scales for this profile.
    gender_baseline_sd <- baseline_sd
    beta_gender_sd <- beta_treatment_sd
    tau_gender_rate <- tau_treatment_rate

    age_baseline_sd <- baseline_sd
    beta_age_sd <- beta_treatment_sd
    tau_age_rate <- tau_treatment_rate

    sigma_intercept_rate <- 1 / outcome_scale
    sigma_slope_rate <- 1 / slope_scale

    sigma_rate <- 1 / outcome_scale

    nu_shape <- 2
    nu_rate <- 0.1
  }


  if (profile == "weak") {

    baseline_mean <- baseline_center
    baseline_sd <- 5 * outcome_scale

    beta_time_mean <- 0
    beta_time_sd <- 5 * slope_scale

    tau_common_rate <- 0.5 / rw_scale

    beta_treatment_sd <- 5 * slope_scale
    tau_treatment_rate <- 0.5 / rw_scale

    arm_baseline_sd_rate <- 0.5 / outcome_scale

    gender_baseline_sd <- baseline_sd
    beta_gender_sd <- beta_treatment_sd
    tau_gender_rate <- tau_treatment_rate

    age_baseline_sd <- baseline_sd
    beta_age_sd <- beta_treatment_sd
    tau_age_rate <- tau_treatment_rate

    sigma_intercept_rate <- 0.5 / outcome_scale
    sigma_slope_rate <- 0.5 / slope_scale

    sigma_rate <- 0.5 / outcome_scale

    nu_shape <- 2
    nu_rate <- 0.05
  }


  if (profile == "regularized") {

    baseline_mean <- baseline_center
    baseline_sd <- 1.5 * outcome_scale

    beta_time_mean <- 0
    beta_time_sd <- 1 * slope_scale

    tau_common_rate <- 2 / rw_scale

    beta_treatment_sd <- 1 * slope_scale
    tau_treatment_rate <- 2 / rw_scale

    arm_baseline_sd_rate <- 2 / outcome_scale

    gender_baseline_sd <- baseline_sd
    beta_gender_sd <- beta_treatment_sd
    tau_gender_rate <- tau_treatment_rate

    age_baseline_sd <- baseline_sd
    beta_age_sd <- beta_treatment_sd
    tau_age_rate <- tau_treatment_rate

    sigma_intercept_rate <- 2 / outcome_scale
    sigma_slope_rate <- 2 / slope_scale

    sigma_rate <- 2 / outcome_scale

    nu_shape <- 2
    nu_rate <- 0.1
  }


  if (profile == "informative") {

    baseline_mean <- baseline_center
    baseline_sd <- 0.75 * outcome_scale

    beta_time_mean <- 0
    beta_time_sd <- 0.5 * slope_scale

    tau_common_rate <- 5 / rw_scale

    beta_treatment_sd <- 0.5 * slope_scale
    tau_treatment_rate <- 5 / rw_scale

    arm_baseline_sd_rate <- 5 / outcome_scale

    gender_baseline_sd <- baseline_sd
    beta_gender_sd <- beta_treatment_sd
    tau_gender_rate <- tau_treatment_rate

    age_baseline_sd <- baseline_sd
    beta_age_sd <- beta_treatment_sd
    tau_age_rate <- tau_treatment_rate

    sigma_intercept_rate <- 5 / outcome_scale
    sigma_slope_rate <- 5 / slope_scale

    sigma_rate <- 5 / outcome_scale

    nu_shape <- 5
    nu_rate <- 0.5
  }


  # ============================================================
  # CUSTOM PROFILE
  # ============================================================

  if (profile == "custom") {

    if (is.null(baseline_mean)) {
      baseline_mean <- baseline_center
    }

    if (is.null(beta_time_mean)) {
      beta_time_mean <- 0
    }

    required_custom <- c(
      "baseline_sd",
      "beta_time_sd",
      "tau_common_rate",
      "beta_treatment_sd",
      "tau_treatment_rate",
      "arm_baseline_sd_rate",
      "sigma_intercept_rate",
      "sigma_slope_rate",
      "sigma_rate",
      "nu_shape",
      "nu_rate"
    )

    custom_values <- mget(
      required_custom,
      envir = environment(),
      inherits = FALSE
    )

    missing_custom <- names(custom_values)[
      vapply(custom_values, is.null, logical(1))
    ]

    if (length(missing_custom) > 0) {
      stop(
        "For `profile = 'custom'`, supply: ",
        paste(missing_custom, collapse = ", "),
        ".",
        call. = FALSE
      )
    }

    # Backwards-compatible defaults for the newly introduced covariate priors.
    # They remain separately represented in the prior object and Stan data, so
    # users can override them independently when desired.
    if (is.null(gender_baseline_sd)) gender_baseline_sd <- baseline_sd
    if (is.null(beta_gender_sd)) beta_gender_sd <- beta_treatment_sd
    if (is.null(tau_gender_rate)) tau_gender_rate <- tau_treatment_rate

    if (is.null(age_baseline_sd)) age_baseline_sd <- baseline_sd
    if (is.null(beta_age_sd)) beta_age_sd <- beta_treatment_sd
    if (is.null(tau_age_rate)) tau_age_rate <- tau_treatment_rate
  }

  # ============================================================
  # VALIDATION
  # ============================================================

  finite_parameters <- list(
    baseline_mean = baseline_mean,
    beta_time_mean = beta_time_mean
  )

  for (name in names(finite_parameters)) {
    value <- finite_parameters[[name]]

    if (length(value) != 1 ||
        !is.numeric(value) ||
        !is.finite(value)) {
      stop(
        "`", name, "` must be a single finite numeric value.",
        call. = FALSE
      )
    }
  }

  positive_parameters <- list(
    baseline_sd = baseline_sd,
    beta_time_sd = beta_time_sd,
    tau_common_rate = tau_common_rate,
    beta_treatment_sd = beta_treatment_sd,
    tau_treatment_rate = tau_treatment_rate,
    arm_baseline_sd_rate = arm_baseline_sd_rate,
    gender_baseline_sd = gender_baseline_sd,
    beta_gender_sd = beta_gender_sd,
    tau_gender_rate = tau_gender_rate,
    age_baseline_sd = age_baseline_sd,
    beta_age_sd = beta_age_sd,
    tau_age_rate = tau_age_rate,
    sigma_intercept_rate = sigma_intercept_rate,
    sigma_slope_rate = sigma_slope_rate,
    sigma_rate = sigma_rate,
    nu_shape = nu_shape,
    nu_rate = nu_rate
  )

  for (name in names(positive_parameters)) {
    value <- positive_parameters[[name]]

    if (length(value) != 1 ||
        !is.numeric(value) ||
        !is.finite(value) ||
        value <= 0) {
      stop(
        "`", name, "` must be a single positive finite numeric value.",
        call. = FALSE
      )
    }
  }

  # ============================================================
  # BUILD PRIOR OBJECT
  # ============================================================

  prior <- list(
    baseline_prior_mean = as.numeric(baseline_mean),
    baseline_prior_sd = as.numeric(baseline_sd),

    beta_time_prior_mean = as.numeric(beta_time_mean),
    beta_time_prior_sd = as.numeric(beta_time_sd),

    tau_common_prior_rate = as.numeric(tau_common_rate),

    beta_treatment_prior_sd = as.numeric(beta_treatment_sd),
    tau_treatment_prior_rate = as.numeric(tau_treatment_rate),

    arm_baseline_sd_prior_rate = as.numeric(arm_baseline_sd_rate),

    # Dedicated priors for Male - Female trajectory.
    gender_baseline_prior_sd = as.numeric(gender_baseline_sd),
    beta_gender_prior_sd = as.numeric(beta_gender_sd),
    tau_gender_prior_rate = as.numeric(tau_gender_rate),

    # Dedicated priors for age > threshold - age <= threshold trajectory.
    age_baseline_prior_sd = as.numeric(age_baseline_sd),
    beta_age_prior_sd = as.numeric(beta_age_sd),
    tau_age_prior_rate = as.numeric(tau_age_rate),

    sigma_intercept_prior_rate = as.numeric(sigma_intercept_rate),
    sigma_slope_prior_rate = as.numeric(sigma_slope_rate),

    sigma_prior_rate = as.numeric(sigma_rate),

    nu_prior_shape = as.numeric(nu_shape),
    nu_prior_rate = as.numeric(nu_rate),

    profile = profile,

    # R-side metadata: useful for inspection but not sent to Stan.
    reference_scales = list(
      baseline_center = baseline_center,
      outcome_scale = outcome_scale,
      slope_scale = slope_scale,
      rw_scale = rw_scale,
      time_span = time_span
    )
  )

  class(prior) <- "mira_prior"

  prior
}


#' Print a MIRA prior specification
#'
#' @param x A `mira_prior` object.
#' @param ... Additional arguments.
#'
#' @export
print.mira_prior <- function(
    x,
    ...
) {

  mira_validate_prior(x)

  cat("\n")
  cat("MIRA prior specification\n")
  cat("========================\n\n")

  cat("Profile: ", x$profile, "\n\n", sep = "")

  cat("Reference-arm baseline:\n")
  cat(
    "  baseline_mean ~ Normal(",
    x$baseline_prior_mean,
    ", ",
    x$baseline_prior_sd,
    ")\n\n",
    sep = ""
  )

  cat("Reference-arm global time trend:\n")
  cat(
    "  beta_time ~ Normal(",
    x$beta_time_prior_mean,
    ", ",
    x$beta_time_prior_sd,
    ")\n\n",
    sep = ""
  )

  cat("Reference trajectory RW1 scale:\n")
  cat(
    "  tau_common ~ Exponential(",
    x$tau_common_prior_rate,
    ")\n\n",
    sep = ""
  )

  cat("Treatment slope deviations:\n")
  cat(
    "  beta_treatment[g] ~ Normal(0, ",
    x$beta_treatment_prior_sd,
    ")\n\n",
    sep = ""
  )

  cat("Treatment trajectory RW1 scales:\n")
  cat(
    "  tau_treatment[g] ~ Exponential(",
    x$tau_treatment_prior_rate,
    ")\n\n",
    sep = ""
  )

  cat("Baseline arm imbalance SD:\n")
  cat(
    "  arm_baseline_sd ~ Exponential(",
    x$arm_baseline_sd_prior_rate,
    ")\n\n",
    sep = ""
  )

  cat("Gender (Male - Female) baseline effect:\n")
  cat(
    "  gender_baseline_effect ~ Normal(0, ",
    x$gender_baseline_prior_sd,
    ")\n\n",
    sep = ""
  )

  cat("Gender time-slope difference:\n")
  cat(
    "  beta_gender_time ~ Normal(0, ",
    x$beta_gender_prior_sd,
    ")\n\n",
    sep = ""
  )

  cat("Gender trajectory RW1 scale:\n")
  cat(
    "  tau_gender ~ Exponential(",
    x$tau_gender_prior_rate,
    ")\n\n",
    sep = ""
  )

  cat("Age-threshold baseline effect (> threshold - <= threshold):\n")
  cat(
    "  age_baseline_effect ~ Normal(0, ",
    x$age_baseline_prior_sd,
    ")\n\n",
    sep = ""
  )

  cat("Age-threshold time-slope difference:\n")
  cat(
    "  beta_age_time ~ Normal(0, ",
    x$beta_age_prior_sd,
    ")\n\n",
    sep = ""
  )

  cat("Age-threshold trajectory RW1 scale:\n")
  cat(
    "  tau_age ~ Exponential(",
    x$tau_age_prior_rate,
    ")\n\n",
    sep = ""
  )

  cat("Subject random intercept SD:\n")
  cat(
    "  sigma_intercept ~ Exponential(",
    x$sigma_intercept_prior_rate,
    ")\n\n",
    sep = ""
  )

  cat("Subject random slope SD:\n")
  cat(
    "  sigma_slope ~ Exponential(",
    x$sigma_slope_prior_rate,
    ")\n\n",
    sep = ""
  )

  cat("Residual Student-t scale:\n")
  cat(
    "  sigma ~ Exponential(",
    x$sigma_prior_rate,
    ")\n\n",
    sep = ""
  )

  cat("Student-t degrees of freedom:\n")
  cat(
    "  nu ~ Gamma(",
    x$nu_prior_shape,
    ", ",
    x$nu_prior_rate,
    ")\n",
    sep = ""
  )

  cat("\n")

  invisible(x)
}


#' Validate a MIRA prior specification
#'
#' @param prior A `mira_prior` object.
#'
#' @return Invisibly returns TRUE.
#'
#' @export
mira_validate_prior <- function(
    prior
) {

  if (!inherits(prior, "mira_prior")) {
    stop(
      "`prior` must be created with `mira_prior()`.",
      call. = FALSE
    )
  }

  required <- c(
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

  missing <- setdiff(required, names(prior))

  if (length(missing) > 0) {
    stop(
      "Invalid MIRA prior. Missing fields: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  finite_parameters <- c(
    "baseline_prior_mean",
    "beta_time_prior_mean"
  )

  positive_parameters <- setdiff(
    required,
    finite_parameters
  )

  for (name in finite_parameters) {
    value <- prior[[name]]

    if (length(value) != 1 ||
        !is.numeric(value) ||
        !is.finite(value)) {
      stop(
        "Invalid MIRA prior field `",
        name,
        "`: expected one finite numeric value.",
        call. = FALSE
      )
    }
  }

  for (name in positive_parameters) {
    value <- prior[[name]]

    if (length(value) != 1 ||
        !is.numeric(value) ||
        !is.finite(value) ||
        value <= 0) {
      stop(
        "Invalid MIRA prior field `",
        name,
        "`: expected one positive finite numeric value.",
        call. = FALSE
      )
    }
  }

  invisible(TRUE)
}


#' Convert MIRA priors to Stan data
#'
#' @param prior A `mira_prior` object.
#'
#' @return A named list suitable for CmdStan.
#'
#' @export
mira_prior_stan_data <- function(
    prior
) {

  mira_validate_prior(prior)

  list(
    baseline_prior_mean =
      prior$baseline_prior_mean,

    baseline_prior_sd =
      prior$baseline_prior_sd,

    beta_time_prior_mean =
      prior$beta_time_prior_mean,

    beta_time_prior_sd =
      prior$beta_time_prior_sd,

    tau_common_prior_rate =
      prior$tau_common_prior_rate,

    beta_treatment_prior_sd =
      prior$beta_treatment_prior_sd,

    tau_treatment_prior_rate =
      prior$tau_treatment_prior_rate,

    arm_baseline_sd_prior_rate =
      prior$arm_baseline_sd_prior_rate,

    gender_baseline_prior_sd =
      prior$gender_baseline_prior_sd,

    beta_gender_prior_sd =
      prior$beta_gender_prior_sd,

    tau_gender_prior_rate =
      prior$tau_gender_prior_rate,

    age_baseline_prior_sd =
      prior$age_baseline_prior_sd,

    beta_age_prior_sd =
      prior$beta_age_prior_sd,

    tau_age_prior_rate =
      prior$tau_age_prior_rate,

    sigma_intercept_prior_rate =
      prior$sigma_intercept_prior_rate,

    sigma_slope_prior_rate =
      prior$sigma_slope_prior_rate,

    sigma_prior_rate =
      prior$sigma_prior_rate,

    nu_prior_shape =
      prior$nu_prior_shape,

    nu_prior_rate =
      prior$nu_prior_rate
  )
}
