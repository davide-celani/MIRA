#' Fit MIRA longitudinal model
#'
#' Fits the longitudinal Student-t mixed-effects model using CmdStan.
#'
#' @param stan_data Data prepared for the MIRA Stan model.
#' @param prior A `mira_prior` object. If omitted, the default
#'   MIRA prior specification is used.
#' @param chains Number of MCMC chains.
#' @param parallel_chains Number of parallel chains.
#' @param iter_warmup Number of warmup iterations.
#' @param iter_sampling Number of sampling iterations.
#' @param seed Random seed.
#' @param refresh Number of iterations between progress messages.
#'
#' @return A CmdStanMCMC object.
#'
#' @export
mira_fit <- function(
    stan_data,
    prior = mira_prior(
      stan_data,
      profile = "default"
    ),
    chains = 4,
    parallel_chains = chains,
    iter_warmup = 1000,
    iter_sampling = 3000,
    seed = 123,
    refresh = 100
) {

  # ------------------------------------------------------------
  # Check data
  # ------------------------------------------------------------

  if (!is.list(stan_data)) {
    stop(
      "`stan_data` must be a list.",
      call. = FALSE
    )
  }


  # ------------------------------------------------------------
  # Required data for the Stan model
  # ------------------------------------------------------------

  required <- c(
    "N",
    "S",
    "K",
    "y",
    "subject",
    "time",
    "time_value",
    "meaningful_change",
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
  # Basic dimension checks
  # ------------------------------------------------------------

  if (stan_data$N < 1) {
    stop(
      "`N` must be >= 1.",
      call. = FALSE
    )
  }


  if (stan_data$S < 1) {
    stop(
      "`S` must be >= 1.",
      call. = FALSE
    )
  }


  if (stan_data$K < 2) {
    stop(
      "`K` must be >= 2.",
      call. = FALSE
    )
  }


  if (length(stan_data$y) != stan_data$N) {
    stop(
      "`length(y)` must equal `N`.",
      call. = FALSE
    )
  }


  if (length(stan_data$subject) != stan_data$N) {
    stop(
      "`length(subject)` must equal `N`.",
      call. = FALSE
    )
  }


  if (length(stan_data$time) != stan_data$N) {
    stop(
      "`length(time)` must equal `N`.",
      call. = FALSE
    )
  }


  if (length(stan_data$time_value) != stan_data$K) {
    stop(
      "`length(time_value)` must equal `K`.",
      call. = FALSE
    )
  }


  if (any(stan_data$subject < 1) ||
      any(stan_data$subject > stan_data$S)) {

    stop(
      "`subject` must contain integers between 1 and S.",
      call. = FALSE
    )
  }


  if (any(stan_data$time < 1) ||
      any(stan_data$time > stan_data$K)) {

    stop(
      "`time` must contain integers between 1 and K.",
      call. = FALSE
    )
  }


  if (anyDuplicated(stan_data$time_value) > 0) {
    stop(
      "`time_value` must contain distinct measurement times.",
      call. = FALSE
    )
  }


  # ------------------------------------------------------------
  # Validate prior
  # ------------------------------------------------------------

  mira_validate_prior(prior)


  # ------------------------------------------------------------
  # Stan model
  # ------------------------------------------------------------

  stan_file <- system.file(
    "stan",
    "gaussian_longitudinal.stan",
    package = "MIRA"
  )


  if (stan_file == "") {
    stop(
      "Could not find ",
      "`gaussian_longitdinal.stan` ",
      "in the MIRA package.",
      call. = FALSE
    )
  }


  # ------------------------------------------------------------
  # Compile Stan model
  # ------------------------------------------------------------

  message(
    "Compiling MIRA Stan model: ",
    "gaussian longitudinal"
  )


  model <- cmdstanr::cmdstan_model(
    stan_file,
    quiet = TRUE
  )


  # ------------------------------------------------------------
  # Add prior specification to Stan data
  # ------------------------------------------------------------

  stan_prior_data <- mira_prior_stan_data(
    prior
  )


  # Check for duplicated names before combining
  duplicated_names <- intersect(
    names(stan_data),
    names(stan_prior_data)
  )


  if (length(duplicated_names) > 0) {
    stop(
      "Duplicated data names between `stan_data` ",
      "and prior specification: ",
      paste(duplicated_names, collapse = ", "),
      call. = FALSE
    )
  }


  stan_data <- c(
    stan_data,
    stan_prior_data
  )


  # ------------------------------------------------------------
  # Initial values
  # ------------------------------------------------------------

  init <- function() {

    list(

      # Population mean at every measurement occasion
      mu_time = rep(
        stan_data$mean_y,
        stan_data$K
      ),


      # Non-centered subject-level random effects
      z_subject = matrix(
        0,
        nrow = 2,
        ncol = stan_data$S
      ),


      # Random-effect standard deviations:
      # 1 = intercept
      # 2 = slope
      sigma_subject = c(

        max(
          stan_data$sd_y,
          0.1
        ),

        max(
          stan_data$sd_y / 10,
          0.01
        )
      ),


      # Initial correlation matrix = identity
      L_subject = diag(2),


      # Residual SD
      sigma = max(
        stan_data$sd_y,
        0.1
      ),


      # Student-t degrees of freedom
      nu = 10
    )
  }


  # ------------------------------------------------------------
  # Sampling
  # ------------------------------------------------------------

  message(
    "Sampling posterior..."
  )


  fit <- model$sample(

    data = stan_data,

    chains = chains,

    parallel_chains = parallel_chains,

    iter_warmup = iter_warmup,

    iter_sampling = iter_sampling,

    seed = seed,

    init = init,

    refresh = refresh
  )


  # ------------------------------------------------------------
  # Return
  # ------------------------------------------------------------

  return(fit)
}
