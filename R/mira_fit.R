#' Fit MIRA longitudinal model
#'
#' Fits the longitudinal Student-t model using CmdStan.
#'
#' @param stan_data Data prepared for Stan.
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


  required <- c(
    "N",
    "S",
    "y",
    "subject",
    "time",
    "time_value",
    "mean_y",
    "sd_y",
    "meaningful_change"
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
      "`gaussian_longitudinal.stan` ",
      "in the MIRA package.",
      call. = FALSE
    )
  }


  # ------------------------------------------------------------
  # Compile
  # ------------------------------------------------------------

  message(
    "Compiling MIRA Stan model: ",
    "gaussian_longitudinal"
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


  stan_data <- c(
    stan_data,
    stan_prior_data
  )


  # ------------------------------------------------------------
  # Initial values
  # ------------------------------------------------------------

  init <- function() {

    list(

      mu_time = rep(
        stan_data$mean_y,
        3
      ),

      z_subject = matrix(
        0,
        nrow = 2,
        ncol = stan_data$S
      ),

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

      L_subject = diag(2),

      sigma = max(
        stan_data$sd_y,
        0.1
      ),

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

    chains =
      chains,

    parallel_chains =
      parallel_chains,

    iter_warmup =
      iter_warmup,

    iter_sampling =
      iter_sampling,

    seed =
      seed,

    init =
      init,

    refresh =
      refresh
  )


  # ------------------------------------------------------------
  # Return
  # ------------------------------------------------------------

  return(fit)
}
