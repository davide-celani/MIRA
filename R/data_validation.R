#' Validate data supplied to a MIRA Stan model
#'
#' @param stan_data Named list containing the data passed to Stan.
#' @param model Character string identifying the MIRA model.
#'
#' @return Invisibly returns TRUE.
#'
#' @keywords internal
mira_validate_stan_data <- function(
    stan_data,
    model = "gaussian_longitudinal"
) {

  # ---------------------------------------------------------
  # Basic check
  # ---------------------------------------------------------

  if (!is.list(stan_data)) {
    stop(
      "`stan_data` must be a list.",
      call. = FALSE
    )
  }


  # ---------------------------------------------------------
  # Required variables by model
  # ---------------------------------------------------------

  required <- switch(

    model,

    gaussian_longitudinal = c(
      "N",
      "S",
      "y",
      "subject",
      "time",
      "mean_y"
    ),

    student_t_longitudinal = c(
      "N",
      "S",
      "y",
      "subject",
      "time",
      "time_value",
      "mean_y",
      "sd_y"
    ),

    stop(
      "No Stan data specification available for model '",
      model,
      "'.",
      call. = FALSE
    )
  )


  # ---------------------------------------------------------
  # Check missing variables
  # ---------------------------------------------------------

  missing <- setdiff(
    required,
    names(stan_data)
  )

  if (length(missing) > 0) {

    stop(
      "Missing Stan data for model '",
      model,
      "': ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }


  # ---------------------------------------------------------
  # Extract core variables
  # ---------------------------------------------------------

  N <- stan_data$N
  S <- stan_data$S
  y <- stan_data$y
  subject <- stan_data$subject
  time <- stan_data$time


  # ---------------------------------------------------------
  # Validate N
  # ---------------------------------------------------------

  if (
    length(N) != 1 ||
    !is.numeric(N) ||
    is.na(N) ||
    N < 1 ||
    N != as.integer(N)
  ) {

    stop(
      "`N` must be a single positive integer.",
      call. = FALSE
    )
  }


  # ---------------------------------------------------------
  # Validate S
  # ---------------------------------------------------------

  if (
    length(S) != 1 ||
    !is.numeric(S) ||
    is.na(S) ||
    S < 1 ||
    S != as.integer(S)
  ) {

    stop(
      "`S` must be a single positive integer.",
      call. = FALSE
    )
  }


  # ---------------------------------------------------------
  # Validate outcome
  # ---------------------------------------------------------

  if (length(y) != N) {

    stop(
      "`length(y)` must equal `N`.",
      call. = FALSE
    )
  }

  if (!is.numeric(y)) {

    stop(
      "`y` must be numeric.",
      call. = FALSE
    )
  }

  if (anyNA(y)) {

    stop(
      "`y` contains missing values.",
      call. = FALSE
    )
  }


  # ---------------------------------------------------------
  # Validate subject
  # ---------------------------------------------------------

  if (length(subject) != N) {

    stop(
      "`length(subject)` must equal `N`.",
      call. = FALSE
    )
  }

  if (!is.numeric(subject) && !is.integer(subject)) {

    stop(
      "`subject` must be numeric or integer.",
      call. = FALSE
    )
  }

  if (anyNA(subject)) {

    stop(
      "`subject` contains missing values.",
      call. = FALSE
    )
  }

  if (any(subject < 1 | subject > S)) {

    stop(
      "`subject` contains values outside 1:S.",
      call. = FALSE
    )
  }


  # ---------------------------------------------------------
  # Validate time
  # ---------------------------------------------------------

  if (length(time) != N) {

    stop(
      "`length(time)` must equal `N`.",
      call. = FALSE
    )
  }

  if (anyNA(time)) {

    stop(
      "`time` contains missing values.",
      call. = FALSE
    )
  }

  if (any(time < 1 | time > 3)) {

    stop(
      "`time` must contain values between 1 and 3.",
      call. = FALSE
    )
  }


  # ---------------------------------------------------------
  # Model-specific validation
  # ---------------------------------------------------------

  if (model == "gaussian_longitudinal") {

    if (
      length(stan_data$mean_y) != 1 ||
      !is.numeric(stan_data$mean_y) ||
      is.na(stan_data$mean_y)
    ) {

      stop(
        "`mean_y` must be a single non-missing numeric value.",
        call. = FALSE
      )
    }
  }


  if (model == "student_t_longitudinal") {

    time_value <- stan_data$time_value
    sd_y <- stan_data$sd_y

    if (length(time_value) != 3) {

      stop(
        "`time_value` must contain exactly 3 values.",
        call. = FALSE
      )
    }

    if (!is.numeric(time_value)) {

      stop(
        "`time_value` must be numeric.",
        call. = FALSE
      )
    }

    if (anyNA(time_value)) {

      stop(
        "`time_value` contains missing values.",
        call. = FALSE
      )
    }

    if (
      length(sd_y) != 1 ||
      !is.numeric(sd_y) ||
      is.na(sd_y) ||
      sd_y <= 0
    ) {

      stop(
        "`sd_y` must be a single positive numeric value.",
        call. = FALSE
      )
    }
  }


  invisible(TRUE)
}
