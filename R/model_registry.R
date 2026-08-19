#' Get the Stan file associated with a MIRA model
#'
#' @param model Character string identifying the MIRA model.
#'
#' @return Path to the Stan file.
#'
#' @keywords internal
mira_stan_file <- function(model = "gaussian_longitudinal") {

  registry <- list(

    gaussian_longitudinal =
      "gaussian_longitudinal.stan"

  )

  if (!is.character(model) || length(model) != 1) {
    stop(
      "`model` must be a single character string.",
      call. = FALSE
    )
  }

  if (!model %in% names(registry)) {
    stop(
      "Unknown MIRA model: '",
      model,
      "'. Available models: ",
      paste(names(registry), collapse = ", "),
      call. = FALSE
    )
  }

  stan_file <- system.file(
    "stan",
    registry[[model]],
    package = "MIRA"
  )

  # Development mode
  if (stan_file == "") {
    stan_file <- file.path(
      "inst",
      "stan",
      registry[[model]]
    )
  }

  if (!file.exists(stan_file)) {
    stop(
      "Stan model not found: ",
      stan_file,
      call. = FALSE
    )
  }

  stan_file
}
