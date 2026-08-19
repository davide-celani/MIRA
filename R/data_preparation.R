#' Prepare Longitudinal Data for the MIRA Stan Model
#'
#' Prepares and validates longitudinal data for fitting the MIRA
#' (Multilevel Individual Response Analysis) Bayesian model implemented
#' in Stan.
#'
#' The function expects the input data in wide format, with one row per
#' subject and three repeated measurements (`t0`, `t1`, and `t2`). It
#' validates the subject identifiers, outcome variables, time points,
#' missing values, and clinical threshold, and then converts the data
#' from wide to long format.
#'
#' The resulting long-format representation follows the ordering required
#' by the MIRA Stan model:
#'
#' \itemize{
#'   \item observations at `t0` for subjects 1, ..., S;
#'   \item observations at `t1` for subjects 1, ..., S;
#'   \item observations at `t2` for subjects 1, ..., S.
#' }
#'
#' Subject identifiers are internally represented using consecutive
#' integer indices from 1 to `S`, while the actual measurement times are
#' supplied separately through `time_value`.
#'
#' The function does not perform any imputation or transformation of the
#' outcome variables. Missing values are therefore not currently
#' supported and result in an error.
#'
#' @param data A data frame containing one row per subject. The data frame
#'   must contain the following columns:
#'   \describe{
#'     \item{patient}{A unique subject identifier. Missing and duplicated
#'       identifiers are not allowed.}
#'     \item{t0}{Numeric outcome measurement at the baseline time point.}
#'     \item{t1}{Numeric outcome measurement at the first follow-up time
#'       point.}
#'     \item{t2}{Numeric outcome measurement at the second follow-up time
#'       point.}
#'   }
#'
#' @param time_value A numeric vector of length three containing the
#'   actual time values corresponding to `t0`, `t1`, and `t2`,
#'   respectively. Values must be finite, distinct, and strictly
#'   increasing. The default is `c(0, 6, 12)`.
#'
#' @param meaningful_change A single finite numeric value defining the
#'   minimum clinically important difference (MCID). This threshold is
#'   passed to the Stan model and is used to define clinically meaningful
#'   individual improvement. The default is `5`.
#'
#' @return A named list containing data formatted for direct use with
#'   the MIRA Stan model. The returned object contains:
#'   \describe{
#'     \item{N}{Total number of observations, equal to `3 * S`.}
#'     \item{S}{Number of subjects.}
#'     \item{y}{Numeric vector containing the outcome measurements in
#'       long format.}
#'     \item{subject}{Integer vector identifying the subject associated
#'       with each observation.}
#'     \item{time}{Integer vector identifying the measurement occasion
#'       (`1`, `2`, or `3`) associated with each observation.}
#'     \item{time_value}{Numeric vector containing the actual values of
#'       the three measurement times.}
#'     \item{mean_y}{Mean of all observed outcome measurements.}
#'     \item{sd_y}{Standard deviation of all observed outcome
#'       measurements.}
#'     \item{meaningful_change}{The MCID supplied through the
#'       `meaningful_change` argument.}
#'   }
#'
#' @details
#' The function performs several validation steps before constructing the
#' Stan data object. These include:
#' \itemize{
#'   \item verification that `data` is a data frame;
#'   \item verification that all required columns are present;
#'   \item verification that at least one subject is available;
#'   \item verification that subject identifiers are non-missing and
#'     unique;
#'   \item verification that exactly three valid measurement times are
#'     supplied;
#'   \item verification that measurement times are distinct and strictly
#'     increasing;
#'   \item verification that the MCID is a single finite numeric value;
#'   \item verification that all outcome variables are numeric;
#'   \item verification that no missing values are present;
#'   \item verification that the outcome standard deviation is positive
#'     and finite.
#' }
#'
#' No statistical model is fitted by this function. Its purpose is
#' restricted to data validation, restructuring, and construction of
#' the data list required by the MIRA Stan model.
#'
#' @section Data ordering:
#' The returned observations are ordered by measurement occasion rather
#' than by subject. For `S` subjects, the ordering is:
#'
#' \preformatted{
#' t0: subject 1, subject 2, ..., subject S
#' t1: subject 1, subject 2, ..., subject S
#' t2: subject 1, subject 2, ..., subject S
#' }
#'
#' This ordering is consistent with the construction of the `subject`
#' and `time` index vectors and must be preserved when passing the
#' resulting list to Stan.
#'
#' @section Missing data:
#' Missing outcome measurements are not currently supported. If any
#' missing value is detected in `patient`, `t0`, `t1`, or `t2`, the
#' function stops with an informative error.
#'
#' @section Clinical interpretation:
#' The `meaningful_change` argument represents the minimum clinically
#' important difference used by the MIRA model to classify individual
#' changes as clinically meaningful. The threshold is not used to modify
#' or transform the observed outcome data.
#'
#' @section Validation:
#' The function stops with an error if any input fails the required
#' validation criteria. This is intended to prevent invalid data from
#' being passed to Stan and to ensure consistency between the R data
#' preparation step and the Stan data block.
#'
#' @examples
#' data <- data.frame(
#'   patient = 1:3,
#'   t0 = c(50, 55, 48),
#'   t1 = c(54, 58, 51),
#'   t2 = c(60, 63, 56)
#' )
#'
#' stan_data <- mira_prepare_data(
#'   data = data,
#'   time_value = c(0, 6, 12),
#'   meaningful_change = 5
#' )
#'
#' str(stan_data)
#'
#' @references
#' Stan Development Team. Stan Reference Manual.
#' \url{https://mc-stan.org/docs/reference-manual/}
#'
#' Gelman, A., Hill, J., & Vehtari, A. (2020).
#' Regression and Other Stories. Cambridge University Press.
#'
#' McElreath, R. (2020).
#' Statistical Rethinking: A Bayesian Course with Examples in R and
#' Stan (2nd ed.). CRC Press.
#'
#' @seealso
#' \code{\link{mira_prepare_data}}
#'
#' @export
mira_prepare_data <- function(
    data,
    time_value = c(0, 6, 12),
    meaningful_change = 5
) {

  # ============================================================
  # BASIC VALIDATION
  # ============================================================

  if (!is.data.frame(data)) {

    stop(
      "`data` must be a data frame.",
      call. = FALSE
    )
  }


  required_columns <- c(
    "patient",
    "t0",
    "t1",
    "t2"
  )


  missing_columns <-
    setdiff(
      required_columns,
      names(data)
    )


  if (length(missing_columns) > 0) {

    stop(
      paste0(
        "Missing required columns: ",
        paste(
          missing_columns,
          collapse = ", "
        )
      ),
      call. = FALSE
    )
  }


  # ============================================================
  # NUMBER OF SUBJECTS
  # ============================================================

  S <- nrow(data)


  if (S < 1) {

    stop(
      "`data` must contain at least one subject.",
      call. = FALSE
    )
  }


  # ============================================================
  # PATIENT IDENTIFIER
  # ============================================================

  if (anyNA(data$patient)) {

    stop(
      "`patient` contains missing values.",
      call. = FALSE
    )
  }


  if (anyDuplicated(data$patient)) {

    stop(
      paste0(
        "Each row of `data` must represent one unique patient. ",
        "Duplicated patient identifiers were found."
      ),
      call. = FALSE
    )
  }


  # ============================================================
  # TIME VALUES
  # ============================================================

  if (
    !is.numeric(time_value) ||
    length(time_value) != 3 ||
    any(!is.finite(time_value))
  ) {

    stop(
      "`time_value` must contain exactly three finite numeric values.",
      call. = FALSE
    )
  }


  if (anyDuplicated(time_value)) {

    stop(
      "`time_value` must contain three distinct time points.",
      call. = FALSE
    )
  }


  # Require chronological ordering

  if (is.unsorted(time_value, strictly = TRUE)) {

    stop(
      "`time_value` must be strictly increasing.",
      call. = FALSE
    )
  }


  # ============================================================
  # MCID
  # ============================================================

  if (
    !is.numeric(meaningful_change) ||
    length(meaningful_change) != 1 ||
    !is.finite(meaningful_change)
  ) {

    stop(
      "`meaningful_change` must be a single finite numeric value.",
      call. = FALSE
    )
  }


  # ============================================================
  # OUTCOME VARIABLES
  # ============================================================

  outcome_variables <- c(
    "t0",
    "t1",
    "t2"
  )


  non_numeric_outcomes <-
    outcome_variables[
      !vapply(
        data[outcome_variables],
        is.numeric,
        logical(1)
      )
    ]


  if (length(non_numeric_outcomes) > 0) {

    stop(
      paste0(
        "The following outcome columns must be numeric: ",
        paste(
          non_numeric_outcomes,
          collapse = ", "
        )
      ),
      call. = FALSE
    )
  }


  # ============================================================
  # MISSING VALUES
  # ============================================================

  if (
    anyNA(
      data[
        c(
          "patient",
          "t0",
          "t1",
          "t2"
        )
      ]
    )
  ) {

    stop(
      paste0(
        "Missing values are not currently supported in ",
        "the MIRA model."
      ),
      call. = FALSE
    )
  }


  # ============================================================
  # LONG-FORMAT OUTCOME
  #
  # Order:
  #
  # t0: subject 1 ... subject S
  # t1: subject 1 ... subject S
  # t2: subject 1 ... subject S
  #
  # This matches the construction of `subject` and `time`.
  # ============================================================

  y <- c(
    data$t0,
    data$t1,
    data$t2
  )


  # ============================================================
  # SUBJECT INDEX
  # ============================================================

  subject <- rep(
    seq_len(S),
    times = 3
  )


  # ============================================================
  # TIME INDEX
  # ============================================================

  time <- rep(
    seq_len(3),
    each = S
  )


  # ============================================================
  # NUMBER OF OBSERVATIONS
  # ============================================================

  N <- length(y)


  # Sanity check

  if (N != 3 * S) {

    stop(
      "Internal error: N must equal 3 * S.",
      call. = FALSE
    )
  }


  # ============================================================
  # OUTCOME MOMENTS
  # ============================================================

  mean_y <- mean(y)

  sd_y <- stats::sd(y)


  if (
    !is.finite(mean_y) ||
    !is.finite(sd_y) ||
    sd_y <= 0
  ) {

    stop(
      "The outcome standard deviation must be positive and finite.",
      call. = FALSE
    )
  }


  # ============================================================
  # RETURN
  # ============================================================

  list(

    # ----------------------------------------------------------
    # Dimensions
    # ----------------------------------------------------------

    N = N,

    S = S,


    # ----------------------------------------------------------
    # Outcome
    # ----------------------------------------------------------

    y = y,


    # ----------------------------------------------------------
    # Indices
    # ----------------------------------------------------------

    subject = subject,

    time = time,


    # ----------------------------------------------------------
    # Actual time values
    # ----------------------------------------------------------

    time_value =
      as.numeric(time_value),


    # ----------------------------------------------------------
    # Outcome moments
    # ----------------------------------------------------------

    mean_y =
      mean_y,

    sd_y =
      sd_y,


    # ----------------------------------------------------------
    # Clinical threshold
    # ----------------------------------------------------------

    meaningful_change =
      meaningful_change
  )
}
