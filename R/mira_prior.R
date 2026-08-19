#' Create a MIRA prior specification
#'
#' Creates a prior specification for the MIRA longitudinal
#' Student-t model.
#'
#' @param stan_data Data prepared by [mira_prepare_data()].
#' @param profile Character string specifying a predefined prior
#'   profile. Available profiles are `"default"`, `"weak"`,
#'   `"regularized"`, `"informative"`, and `"custom"`.
#' @param mu_time_mean Prior mean for the population time-specific means.
#' @param mu_time_sd Prior standard deviation for the population
#'   time-specific means.
#' @param sigma_intercept_rate Rate parameter for the exponential
#'   prior on the random-intercept standard deviation.
#' @param sigma_slope_rate Rate parameter for the exponential
#'   prior on the random-slope standard deviation.
#' @param sigma_rate Rate parameter for the exponential prior on
#'   the residual standard deviation.
#' @param nu_shape Shape parameter for the Gamma prior on Student-t
#'   degrees of freedom.
#' @param nu_rate Rate parameter for the Gamma prior on Student-t
#'   degrees of freedom.
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
    mu_time_mean = NULL,
    mu_time_sd = NULL,
    sigma_intercept_rate = NULL,
    sigma_slope_rate = NULL,
    sigma_rate = NULL,
    nu_shape = NULL,
    nu_rate = NULL
) {

  # ------------------------------------------------------------
  # Check data
  # ------------------------------------------------------------

  if (!is.list(stan_data)) {
    stop("`stan_data` must be a list.", call. = FALSE)
  }

  required <- c(
    "mean_y",
    "sd_y"
  )

  missing <- setdiff(
    required,
    names(stan_data)
  )

  if (length(missing) > 0) {
    stop(
      "Missing Stan data: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }


  # ------------------------------------------------------------
  # Profile
  # ------------------------------------------------------------

  profile <- match.arg(profile)


  # ------------------------------------------------------------
  # Scale
  # ------------------------------------------------------------

  mean_y <- stan_data$mean_y
  sd_y <- stan_data$sd_y


  if (!is.finite(mean_y)) {
    stop(
      "`stan_data$mean_y` must be finite.",
      call. = FALSE
    )
  }

  if (!is.finite(sd_y) || sd_y <= 0) {
    stop(
      "`stan_data$sd_y` must be positive and finite.",
      call. = FALSE
    )
  }


  # ------------------------------------------------------------
  # Default profile
  #
  # This reproduces the priors from the original model.
  # ------------------------------------------------------------

  if (profile == "default") {

    mu_time_mean <- mean_y
    mu_time_sd <- 5 * sd_y

    sigma_intercept_rate <- 1 / sd_y
    sigma_slope_rate <- 1 / sd_y

    sigma_rate <- 1 / sd_y

    nu_shape <- 2
    nu_rate <- 0.1
  }


  # ------------------------------------------------------------
  # Weak profile
  # ------------------------------------------------------------

  if (profile == "weak") {

    mu_time_mean <- mean_y
    mu_time_sd <- 10 * sd_y

    sigma_intercept_rate <- 0.5 / sd_y
    sigma_slope_rate <- 0.5 / sd_y

    sigma_rate <- 0.5 / sd_y

    nu_shape <- 2
    nu_rate <- 0.05
  }


  # ------------------------------------------------------------
  # Regularized profile
  # ------------------------------------------------------------

  if (profile == "regularized") {

    mu_time_mean <- mean_y
    mu_time_sd <- 2 * sd_y

    sigma_intercept_rate <- 2 / sd_y
    sigma_slope_rate <- 2 / sd_y

    sigma_rate <- 2 / sd_y

    nu_shape <- 2
    nu_rate <- 0.1
  }


  # ------------------------------------------------------------
  # Informative profile
  # ------------------------------------------------------------

  if (profile == "informative") {

    mu_time_mean <- mean_y
    mu_time_sd <- 1 * sd_y

    sigma_intercept_rate <- 5 / sd_y
    sigma_slope_rate <- 5 / sd_y

    sigma_rate <- 5 / sd_y

    nu_shape <- 5
    nu_rate <- 0.5
  }


  # ------------------------------------------------------------
  # Custom profile
  # ------------------------------------------------------------

  if (profile == "custom") {

    if (is.null(mu_time_mean)) {
      mu_time_mean <- mean_y
    }

    if (is.null(mu_time_sd)) {
      stop(
        "`mu_time_sd` must be supplied for `profile = 'custom'`.",
        call. = FALSE
      )
    }

    if (is.null(sigma_intercept_rate)) {
      stop(
        "`sigma_intercept_rate` must be supplied for ",
        "`profile = 'custom'`.",
        call. = FALSE
      )
    }

    if (is.null(sigma_slope_rate)) {
      stop(
        "`sigma_slope_rate` must be supplied for ",
        "`profile = 'custom'`.",
        call. = FALSE
      )
    }

    if (is.null(sigma_rate)) {
      stop(
        "`sigma_rate` must be supplied for ",
        "`profile = 'custom'`.",
        call. = FALSE
      )
    }

    if (is.null(nu_shape)) {
      stop(
        "`nu_shape` must be supplied for ",
        "`profile = 'custom'`.",
        call. = FALSE
      )
    }

    if (is.null(nu_rate)) {
      stop(
        "`nu_rate` must be supplied for ",
        "`profile = 'custom'`.",
        call. = FALSE
      )
    }
  }


  # ------------------------------------------------------------
  # Validation
  # ------------------------------------------------------------

  values <- list(
    mu_time_mean = mu_time_mean,
    mu_time_sd = mu_time_sd,
    sigma_intercept_rate = sigma_intercept_rate,
    sigma_slope_rate = sigma_slope_rate,
    sigma_rate = sigma_rate,
    nu_shape = nu_shape,
    nu_rate = nu_rate
  )

  for (name in names(values)) {

    value <- values[[name]]

    if (
      length(value) != 1 ||
      !is.numeric(value) ||
      !is.finite(value)
    ) {

      stop(
        "`",
        name,
        "` must be a single finite numeric value.",
        call. = FALSE
      )
    }
  }


  if (mu_time_sd <= 0) {
    stop(
      "`mu_time_sd` must be > 0.",
      call. = FALSE
    )
  }

  if (sigma_intercept_rate <= 0) {
    stop(
      "`sigma_intercept_rate` must be > 0.",
      call. = FALSE
    )
  }

  if (sigma_slope_rate <= 0) {
    stop(
      "`sigma_slope_rate` must be > 0.",
      call. = FALSE
    )
  }

  if (sigma_rate <= 0) {
    stop(
      "`sigma_rate` must be > 0.",
      call. = FALSE
    )
  }

  if (nu_shape <= 0) {
    stop(
      "`nu_shape` must be > 0.",
      call. = FALSE
    )
  }

  if (nu_rate <= 0) {
    stop(
      "`nu_rate` must be > 0.",
      call. = FALSE
    )
  }


  # ------------------------------------------------------------
  # Build object
  # ------------------------------------------------------------

  prior <- list(

    mu_time_prior_mean =
      mu_time_mean,

    mu_time_prior_sd =
      mu_time_sd,

    sigma_intercept_prior_rate =
      sigma_intercept_rate,

    sigma_slope_prior_rate =
      sigma_slope_rate,

    sigma_prior_rate =
      sigma_rate,

    nu_prior_shape =
      nu_shape,

    nu_prior_rate =
      nu_rate,

    profile =
      profile
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

  cat("\n")
  cat("MIRA prior specification\n")
  cat("------------------------\n")

  cat(
    "Profile: ",
    x$profile,
    "\n\n",
    sep = ""
  )

  cat(
    "Population means:\n",
    "  mu_time ~ Normal(",
    x$mu_time_prior_mean,
    ", ",
    x$mu_time_prior_sd,
    ")\n\n",
    sep = ""
  )

  cat(
    "Random intercept:\n",
    "  sigma_intercept ~ Exponential(",
    x$sigma_intercept_prior_rate,
    ")\n\n",
    sep = ""
  )

  cat(
    "Random slope:\n",
    "  sigma_slope ~ Exponential(",
    x$sigma_slope_prior_rate,
    ")\n\n",
    sep = ""
  )

  cat(
    "Residual SD:\n",
    "  sigma ~ Exponential(",
    x$sigma_prior_rate,
    ")\n\n",
    sep = ""
  )

  cat(
    "Student-t degrees of freedom:\n",
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
mira_validate_prior <- function(prior) {

  if (!inherits(prior, "mira_prior")) {

    stop(
      "`prior` must be created with `mira_prior()`.",
      call. = FALSE
    )
  }

  required <- c(
    "mu_time_prior_mean",
    "mu_time_prior_sd",
    "sigma_intercept_prior_rate",
    "sigma_slope_prior_rate",
    "sigma_prior_rate",
    "nu_prior_shape",
    "nu_prior_rate"
  )

  missing <- setdiff(
    required,
    names(prior)
  )

  if (length(missing) > 0) {

    stop(
      "Invalid MIRA prior. Missing fields: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
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
mira_prior_stan_data <- function(prior) {

  mira_validate_prior(prior)

  list(

    mu_time_prior_mean =
      prior$mu_time_prior_mean,

    mu_time_prior_sd =
      prior$mu_time_prior_sd,

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
