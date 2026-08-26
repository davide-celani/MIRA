#' Prepare Longitudinal Data for the MIRA Stan Model
#'
#' Prepares and validates longitudinal data for the MIRA
#' (Multilevel Individual Response Analysis) Bayesian model implemented
#' in Stan.
#'
#' Longitudinal measurement columns are detected automatically using
#' the naming convention:
#'
#'   <outcome>_t0
#'   <outcome>_t1
#'   <outcome>_t2
#'   ...
#'   <outcome>_tK
#'
#' The number of measurement occasions is therefore determined
#' automatically from the input data.
#'
#' @param data A data frame containing one row per subject. It must
#'   contain a `patient` identifier and at least two longitudinal
#'   measurement columns following the `<outcome>_t0`, ..., `<outcome>_tK`
#'   naming convention.
#'
#' @param time_value Numeric vector containing the actual measurement
#'   times corresponding to t0, t1, ..., tK.
#'
#' @param meaningful_change Single finite numeric value defining the
#'   minimum clinically important difference (MCID).
#'
#' @return A named list containing data directly suitable for CmdStan.
#'
#' @export
mira_prepare_data <- function(
    data,
    time_value,
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


  # ============================================================
  # PATIENT IDENTIFIER
  # ============================================================

  if (!"patient" %in% names(data)) {

    stop(
      "The data frame must contain a `patient` column.",
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
  # PATIENT VALIDATION
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
  # AUTOMATIC DETECTION OF LONGITUDINAL COLUMNS
  #
  # Expected:
  #
  # IOP_t0
  # IOP_t1
  # IOP_t2
  # IOP_t3
  #
  # or:
  #
  # BCVA_t0
  # BCVA_t1
  #
  # ============================================================

  measurement_columns <- names(data)[
    grepl(
      "_t[0-9]+$",
      names(data)
    )
  ]


  # ============================================================
  # MINIMUM NUMBER OF TIME POINTS
  # ============================================================

  if (length(measurement_columns) < 2) {

    stop(
      paste0(
        "At least two longitudinal measurement columns are required. ",
        "Expected columns named `<outcome>_t0`, `<outcome>_t1`, ..., ",
        "`<outcome>_tK`."
      ),
      call. = FALSE
    )
  }


  # ============================================================
  # EXTRACT TIME INDICES
  # ============================================================

  time_indices <- as.integer(
    sub(
      "^.*_t([0-9]+)$",
      "\\1",
      measurement_columns
    )
  )


  if (anyNA(time_indices)) {

    stop(
      "Could not determine the time index from the longitudinal column names.",
      call. = FALSE
    )
  }


  # ============================================================
  # CHECK DUPLICATE TIME INDICES
  # ============================================================

  if (anyDuplicated(time_indices)) {

    duplicated_times <- unique(
      time_indices[
        duplicated(time_indices)
      ]
    )

    stop(
      paste0(
        "Multiple measurement columns correspond to the same time index: ",
        paste(
          duplicated_times,
          collapse = ", "
        ),
        "."
      ),
      call. = FALSE
    )
  }


  # ============================================================
  # ORDER COLUMNS CHRONOLOGICALLY
  # ============================================================

  order_time <- order(time_indices)

  measurement_columns <- measurement_columns[
    order_time
  ]

  time_indices <- time_indices[
    order_time
  ]


  # ============================================================
  # REQUIRE t0
  # ============================================================

  if (time_indices[1] != 0) {

    stop(
      paste0(
        "The first longitudinal measurement must be `t0`. ",
        "Detected first time index: t",
        time_indices[1],
        "."
      ),
      call. = FALSE
    )
  }


  # ============================================================
  # REQUIRE CONSECUTIVE TIME INDICES
  #
  # Valid:
  #
  # t0 t1
  # t0 t1 t2
  # t0 t1 t2 t3
  # t0 ... tK
  #
  # Invalid:
  #
  # t0 t2
  # t0 t1 t3
  #
  # ============================================================

  expected_time_indices <- seq(
    from = 0,
    to = length(time_indices) - 1
  )


  if (!identical(
    time_indices,
    expected_time_indices
  )) {

    stop(
      paste0(
        "Longitudinal measurement columns must have consecutive ",
        "time indices starting at t0. Detected: ",
        paste(
          paste0("t", time_indices),
          collapse = ", "
        ),
        ". Expected: ",
        paste(
          paste0("t", expected_time_indices),
          collapse = ", "
        ),
        "."
      ),
      call. = FALSE
    )
  }


  # ============================================================
  # NUMBER OF TIME POINTS
  # ============================================================

  K <- length(measurement_columns)


  # ============================================================
  # CHECK SAME OUTCOME
  #
  # Valid:
  #
  # IOP_t0 IOP_t1 IOP_t2 IOP_t3
  #
  # Invalid:
  #
  # IOP_t0 IOP_t1 BCVA_t2
  #
  # ============================================================

  outcome_names <- sub(
    "_t[0-9]+$",
    "",
    measurement_columns
  )


  if (length(unique(outcome_names)) != 1) {

    stop(
      paste0(
        "All longitudinal measurement columns must refer to ",
        "the same outcome. Detected prefixes: ",
        paste(
          unique(outcome_names),
          collapse = ", "
        ),
        "."
      ),
      call. = FALSE
    )
  }


  # ============================================================
  # TIME VALUES
  # ============================================================

  if (
    !is.numeric(time_value) ||
    length(time_value) != K ||
    any(!is.finite(time_value))
  ) {

    stop(
      paste0(
        "`time_value` must contain exactly ",
        K,
        " finite numeric values corresponding to ",
        "t0 ... t",
        K - 1,
        "."
      ),
      call. = FALSE
    )
  }


  # ============================================================
  # DISTINCT TIME VALUES
  # ============================================================

  if (anyDuplicated(time_value)) {

    stop(
      "`time_value` must contain distinct measurement times.",
      call. = FALSE
    )
  }


  # ============================================================
  # CHRONOLOGICAL TIME VALUES
  # ============================================================

  if (is.unsorted(
    time_value,
    strictly = TRUE
  )) {

    stop(
      "`time_value` must be strictly increasing.",
      call. = FALSE
    )
  }


  # ============================================================
  # MCID VALIDATION
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
  # OUTCOME TYPE VALIDATION
  # ============================================================

  non_numeric_outcomes <- measurement_columns[
    !vapply(
      data[measurement_columns],
      is.numeric,
      logical(1)
    )
  ]


  if (length(non_numeric_outcomes) > 0) {

    stop(
      paste0(
        "The following longitudinal measurement columns ",
        "must be numeric: ",
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
          measurement_columns
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
  # CONSTRUCT LONG-FORMAT OUTCOME
  #
  # For:
  #
  # IOP_t0 IOP_t1 IOP_t2 IOP_t3
  #
  # the resulting order is:
  #
  # t0: patient 1 ... patient S
  # t1: patient 1 ... patient S
  # t2: patient 1 ... patient S
  # t3: patient 1 ... patient S
  #
  # ============================================================

  y <- unlist(
    data[measurement_columns],
    use.names = FALSE
  )


  # ============================================================
  # SUBJECT INDEX
  # ============================================================

  subject <- rep(
    seq_len(S),
    times = K
  )


  # ============================================================
  # TIME INDEX
  # ============================================================

  time <- rep(
    seq_len(K),
    each = S
  )


  # ============================================================
  # NUMBER OF OBSERVATIONS
  # ============================================================

  N <- length(y)


  # ============================================================
  # INTERNAL CONSISTENCY
  # ============================================================

  if (N != K * S) {

    stop(
      "Internal error: N must equal K * S.",
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
  # RETURN DATA FOR STAN
  #
  # IMPORTANT:
  #
  # Only objects actually required by Stan are returned.
  #
  # `outcome_name` and `measurement_columns` are deliberately
  # NOT returned here.
  #
  # ============================================================

  list(

    # ----------------------------------------------------------
    # Dimensions
    # ----------------------------------------------------------

    N = as.integer(N),

    S = as.integer(S),

    K = as.integer(K),


    # ----------------------------------------------------------
    # Outcome
    # ----------------------------------------------------------

    y = as.numeric(y),


    # ----------------------------------------------------------
    # Indices
    # ----------------------------------------------------------

    subject = as.integer(subject),

    time = as.integer(time),


    # ----------------------------------------------------------
    # Actual measurement times
    # ----------------------------------------------------------

    time_value = as.numeric(
      time_value
    ),


    # ----------------------------------------------------------
    # Outcome moments
    # ----------------------------------------------------------

    mean_y = as.numeric(
      mean_y
    ),

    sd_y = as.numeric(
      sd_y
    ),


    # ----------------------------------------------------------
    # Clinical threshold
    # ----------------------------------------------------------

    meaningful_change = as.numeric(
      meaningful_change
    )
  )
}
