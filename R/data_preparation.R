#' Prepare longitudinal data for the MIRA treatment model
#'
#' Prepares and validates longitudinal data for the MIRA Bayesian
#' longitudinal Student-t mixed-effects model with treatment-specific
#' trajectories.
#'
#' Longitudinal measurement columns are detected automatically using
#' the naming convention `<outcome>_t0`, ..., `<outcome>_tK`.
#' The number of measurement occasions is therefore determined from
#' the input data and may be any K >= 2.
#'
#' @param data A data frame with one row per subject. It must contain a
#'   `patient` identifier, a treatment-arm column, and at least two
#'   longitudinal measurement columns named `<outcome>_t0`, ..., `<outcome>_tK`.
#' @param time_value Numeric vector of actual measurement times corresponding
#'   to t0, ..., tK. Values must be finite and strictly increasing.
#' @param meaningful_change Non-negative numeric value giving the prior mean
#'   of the minimum clinically important difference (MCID).
#' @param meaningful_change_sd Positive numeric value giving the prior SD of
#'   the uncertain MCID. This must be externally specified; the outcome data
#'   should not be used to identify MCID uncertainty.
#' @param direction Direction of clinical improvement. Use `1` or
#'   `"higher"` when higher outcome values are better; use `-1` or
#'   `"lower"` when lower outcome values are better.
#' @param meaningful_between_arm_difference Non-negative threshold defining
#'   a clinically meaningful between-arm difference in change. By default it
#'   is set equal to `meaningful_change`.
#' @param arm_column Name of the treatment-arm column in `data`.
#' @param reference_arm Value identifying the reference arm. If NULL, the
#'   first observed arm is used and a warning is emitted.
#' @param gender_column Name of the sex/gender column in `data`.
#' @param female_label Value identifying the Female reference category.
#' @param male_label Value identifying the Male category. The Stan model
#'   estimates the Male - Female difference at each measurement occasion.
#' @param age_column Name of the age column in `data`.
#' @param age_threshold Numeric age threshold used to define two groups.
#'   Subjects with age > `age_threshold` are coded as the older group;
#'   subjects with age <= `age_threshold` are the reference group.
#'
#' @return A named list containing the variables required by the MIRA Stan
#'   model plus R-side metadata (`mean_y`, `sd_y`, `arm_labels`,
#'   `outcome_name`, `measurement_columns`, gender coding, age threshold,
#'   and observed group counts used for interpretation of covariate effects).
#'
#' @export
mira_prepare_data <- function(
    data,
    time_value,
    meaningful_change = 5,
    meaningful_change_sd,
    direction,
    meaningful_between_arm_difference = meaningful_change,
    arm_column = "arm",
    reference_arm = NULL,
    gender_column = "gender",
    female_label = "Female",
    male_label = "Male",
    age_column = "age",
    age_threshold = 60
) {

  # ============================================================
  # BASIC VALIDATION
  # ============================================================

  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.", call. = FALSE)
  }

  if (!is.character(arm_column) || length(arm_column) != 1 || !nzchar(arm_column)) {
    stop("`arm_column` must be one non-empty character string.", call. = FALSE)
  }

  if (!is.character(gender_column) || length(gender_column) != 1 || !nzchar(gender_column)) {
    stop("`gender_column` must be one non-empty character string.", call. = FALSE)
  }

  if (!is.character(age_column) || length(age_column) != 1 || !nzchar(age_column)) {
    stop("`age_column` must be one non-empty character string.", call. = FALSE)
  }

  if (!"patient" %in% names(data)) {
    stop("The data frame must contain a `patient` column.", call. = FALSE)
  }

  if (!arm_column %in% names(data)) {
    stop(
      "The data frame must contain the treatment-arm column `",
      arm_column,
      "`.",
      call. = FALSE
    )
  }

  if (!gender_column %in% names(data)) {
    stop(
      "The data frame must contain the gender column `",
      gender_column,
      "`.",
      call. = FALSE
    )
  }

  if (!age_column %in% names(data)) {
    stop(
      "The data frame must contain the age column `",
      age_column,
      "`.",
      call. = FALSE
    )
  }

  S <- nrow(data)

  if (S < 1) {
    stop("`data` must contain at least one subject.", call. = FALSE)
  }

  # ============================================================
  # PATIENT IDENTIFIER
  # ============================================================

  if (anyNA(data$patient)) {
    stop("`patient` contains missing values.", call. = FALSE)
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
  # TREATMENT ARMS
  # ============================================================

  arm_raw <- data[[arm_column]]

  if (anyNA(arm_raw)) {
    stop("`", arm_column, "` contains missing values.", call. = FALSE)
  }

  arm_chr <- as.character(arm_raw)

  if (any(!nzchar(arm_chr))) {
    stop("`", arm_column, "` contains empty treatment labels.", call. = FALSE)
  }

  observed_arms <- unique(arm_chr)

  if (length(observed_arms) < 2) {
    stop(
      "The current MIRA treatment model requires at least two treatment arms.",
      call. = FALSE
    )
  }

  if (is.null(reference_arm)) {
    reference_arm_chr <- observed_arms[1]

    warning(
      "`reference_arm` was not supplied; using the first observed arm (`",
      reference_arm_chr,
      "`) as the reference group.",
      call. = FALSE
    )
  } else {
    if (length(reference_arm) != 1 || is.na(reference_arm)) {
      stop("`reference_arm` must identify exactly one non-missing arm.", call. = FALSE)
    }

    reference_arm_chr <- as.character(reference_arm)

    if (!reference_arm_chr %in% observed_arms) {
      stop(
        "`reference_arm` = `",
        reference_arm_chr,
        "` was not found in `",
        arm_column,
        "`.",
        call. = FALSE
      )
    }
  }

  arm_labels <- c(
    reference_arm_chr,
    observed_arms[observed_arms != reference_arm_chr]
  )

  arm <- match(arm_chr, arm_labels)
  G <- length(arm_labels)

  if (anyNA(arm) || any(arm < 1L) || any(arm > G)) {
    stop("Internal error while encoding treatment arms.", call. = FALSE)
  }

  # ============================================================
  # GENDER: FEMALE REFERENCE, MALE INDICATOR
  # ============================================================

  if (
    length(female_label) != 1 ||
    is.na(female_label) ||
    !nzchar(as.character(female_label))
  ) {
    stop("`female_label` must identify one non-missing category.", call. = FALSE)
  }

  if (
    length(male_label) != 1 ||
    is.na(male_label) ||
    !nzchar(as.character(male_label))
  ) {
    stop("`male_label` must identify one non-missing category.", call. = FALSE)
  }

  female_label_chr <- as.character(female_label)
  male_label_chr <- as.character(male_label)

  if (identical(female_label_chr, male_label_chr)) {
    stop("`female_label` and `male_label` must be different.", call. = FALSE)
  }

  gender_raw <- data[[gender_column]]

  if (anyNA(gender_raw)) {
    stop("`", gender_column, "` contains missing values.", call. = FALSE)
  }

  gender_chr <- as.character(gender_raw)
  allowed_gender <- c(female_label_chr, male_label_chr)

  if (any(!gender_chr %in% allowed_gender)) {
    bad_gender <- unique(gender_chr[!gender_chr %in% allowed_gender])

    stop(
      paste0(
        "`", gender_column, "` contains unsupported categories: ",
        paste(bad_gender, collapse = ", "),
        ". Expected only `", female_label_chr, "` and `", male_label_chr, "`."
      ),
      call. = FALSE
    )
  }

  if (!all(allowed_gender %in% unique(gender_chr))) {
    stop(
      paste0(
        "Both gender categories are required to estimate the Male - Female ",
        "difference. Expected `", female_label_chr, "` and `", male_label_chr, "`."
      ),
      call. = FALSE
    )
  }

  male <- as.integer(gender_chr == male_label_chr)

  # ============================================================
  # AGE THRESHOLD GROUP
  # ============================================================

  age_raw <- data[[age_column]]

  if (!is.numeric(age_raw)) {
    stop("`", age_column, "` must be numeric.", call. = FALSE)
  }

  if (anyNA(age_raw) || any(!is.finite(age_raw))) {
    stop("`", age_column, "` must contain only finite non-missing values.", call. = FALSE)
  }

  if (
    length(age_threshold) != 1 ||
    !is.numeric(age_threshold) ||
    !is.finite(age_threshold)
  ) {
    stop("`age_threshold` must be one finite numeric value.", call. = FALSE)
  }

  age_above_threshold <- as.integer(age_raw > age_threshold)

  if (length(unique(age_above_threshold)) < 2) {
    stop(
      paste0(
        "`age_threshold` = ", age_threshold,
        " does not create two observed age groups. Choose a threshold with ",
        "at least one subject on each side."
      ),
      call. = FALSE
    )
  }

  # ============================================================
  # LONGITUDINAL MEASUREMENT COLUMNS
  # ============================================================

  measurement_columns <- names(data)[
    grepl("_t[0-9]+$", names(data))
  ]

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

  time_indices <- as.integer(
    sub("^.*_t([0-9]+)$", "\\1", measurement_columns)
  )

  if (anyNA(time_indices)) {
    stop(
      "Could not determine the time index from the longitudinal column names.",
      call. = FALSE
    )
  }

  if (anyDuplicated(time_indices)) {
    duplicated_times <- unique(time_indices[duplicated(time_indices)])

    stop(
      paste0(
        "Multiple measurement columns correspond to the same time index: ",
        paste(duplicated_times, collapse = ", "),
        "."
      ),
      call. = FALSE
    )
  }

  order_time <- order(time_indices)
  measurement_columns <- measurement_columns[order_time]
  time_indices <- time_indices[order_time]

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

  expected_time_indices <- seq.int(0L, length(time_indices) - 1L)

  if (!identical(time_indices, expected_time_indices)) {
    stop(
      paste0(
        "Longitudinal measurement columns must have consecutive ",
        "time indices starting at t0. Detected: ",
        paste(paste0("t", time_indices), collapse = ", "),
        ". Expected: ",
        paste(paste0("t", expected_time_indices), collapse = ", "),
        "."
      ),
      call. = FALSE
    )
  }

  K <- length(measurement_columns)

  outcome_names <- sub("_t[0-9]+$", "", measurement_columns)

  if (length(unique(outcome_names)) != 1) {
    stop(
      paste0(
        "All longitudinal measurement columns must refer to the same outcome. ",
        "Detected prefixes: ",
        paste(unique(outcome_names), collapse = ", "),
        "."
      ),
      call. = FALSE
    )
  }

  outcome_name <- unique(outcome_names)

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
        " finite numeric values corresponding to t0 ... t",
        K - 1,
        "."
      ),
      call. = FALSE
    )
  }

  if (anyDuplicated(time_value)) {
    stop("`time_value` must contain distinct measurement times.", call. = FALSE)
  }

  if (is.unsorted(time_value, strictly = TRUE)) {
    stop("`time_value` must be strictly increasing.", call. = FALSE)
  }

  # ============================================================
  # DIRECTION OF IMPROVEMENT
  # ============================================================

  if (is.character(direction)) {
    if (length(direction) != 1 || is.na(direction)) {
      stop("`direction` must identify exactly one improvement direction.", call. = FALSE)
    }

    direction_key <- tolower(trimws(direction))

    if (direction_key %in% c("higher", "higher_better", "increase", "increasing")) {
      direction <- 1L
    } else if (direction_key %in% c("lower", "lower_better", "decrease", "decreasing")) {
      direction <- -1L
    } else {
      stop(
        "Character `direction` must be `higher` or `lower` ",
        "(or a supported synonym).",
        call. = FALSE
      )
    }
  }

  if (
    length(direction) != 1 ||
    !is.numeric(direction) ||
    !is.finite(direction) ||
    !(direction %in% c(-1, 1))
  ) {
    stop("`direction` must be exactly +1 or -1.", call. = FALSE)
  }

  direction <- as.integer(direction)

  # ============================================================
  # CLINICAL THRESHOLDS
  # ============================================================

  if (
    length(meaningful_change) != 1 ||
    !is.numeric(meaningful_change) ||
    !is.finite(meaningful_change) ||
    meaningful_change < 0
  ) {
    stop(
      "`meaningful_change` must be one finite non-negative numeric value.",
      call. = FALSE
    )
  }

  if (missing(meaningful_change_sd)) {
    stop(
      paste0(
        "`meaningful_change_sd` must be supplied for the new model. ",
        "It defines the externally informed prior uncertainty of the MCID."
      ),
      call. = FALSE
    )
  }

  if (
    length(meaningful_change_sd) != 1 ||
    !is.numeric(meaningful_change_sd) ||
    !is.finite(meaningful_change_sd) ||
    meaningful_change_sd <= 0
  ) {
    stop(
      "`meaningful_change_sd` must be one finite positive numeric value.",
      call. = FALSE
    )
  }

  if (
    length(meaningful_between_arm_difference) != 1 ||
    !is.numeric(meaningful_between_arm_difference) ||
    !is.finite(meaningful_between_arm_difference) ||
    meaningful_between_arm_difference < 0
  ) {
    stop(
      paste0(
        "`meaningful_between_arm_difference` must be one finite ",
        "non-negative numeric value."
      ),
      call. = FALSE
    )
  }

  # ============================================================
  # OUTCOME VALIDATION
  # ============================================================

  non_numeric_outcomes <- measurement_columns[
    !vapply(data[measurement_columns], is.numeric, logical(1))
  ]

  if (length(non_numeric_outcomes) > 0) {
    stop(
      paste0(
        "The following longitudinal measurement columns must be numeric: ",
        paste(non_numeric_outcomes, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  if (anyNA(data[c("patient", arm_column, gender_column, age_column, measurement_columns)])) {
    stop(
      "Missing values are not currently supported in the MIRA model.",
      call. = FALSE
    )
  }

  # ============================================================
  # CONSTRUCT LONG-FORMAT DATA
  # ============================================================

  # unlist(data[measurement_columns]) is column-major:
  # t0 for all subjects, then t1 for all subjects, ..., tK.
  y <- unlist(
    data[measurement_columns],
    use.names = FALSE
  )

  subject <- rep(
    seq_len(S),
    times = K
  )

  time <- rep(
    seq_len(K),
    each = S
  )

  N <- length(y)

  if (N != K * S) {
    stop("Internal error: N must equal K * S.", call. = FALSE)
  }

  # ============================================================
  # OUTCOME MOMENTS FOR R-SIDE PRIOR/INITIALIZATION HELPERS
  # ============================================================

  mean_y <- mean(y)
  sd_y <- stats::sd(y)

  if (!is.finite(mean_y) || !is.finite(sd_y) || sd_y <= 0) {
    stop(
      "The outcome standard deviation must be positive and finite.",
      call. = FALSE
    )
  }

  # ============================================================
  # RETURN
  # ============================================================

  list(
    # Stan model data
    N = as.integer(N),
    S = as.integer(S),
    K = as.integer(K),
    G = as.integer(G),
    y = as.numeric(y),
    subject = as.integer(subject),
    time = as.integer(time),
    arm = as.integer(arm),
    male = as.integer(male),
    age_above_threshold = as.integer(age_above_threshold),
    time_value = as.numeric(time_value),
    direction = as.integer(direction),
    mcid_prior_mean = as.numeric(meaningful_change),
    mcid_prior_sd = as.numeric(meaningful_change_sd),
    meaningful_between_arm_difference = as.numeric(
      meaningful_between_arm_difference
    ),

    # R-side metadata; mira_fit() does not pass these to Stan
    mean_y = as.numeric(mean_y),
    sd_y = as.numeric(sd_y),
    meaningful_change = as.numeric(meaningful_change),
    arm_labels = arm_labels,
    reference_arm = reference_arm_chr,
    gender_column = gender_column,
    gender_labels = c(reference = female_label_chr, comparison = male_label_chr),
    gender_counts = stats::setNames(
      c(sum(male == 0L), sum(male == 1L)),
      c(female_label_chr, male_label_chr)
    ),
    age_column = age_column,
    age_threshold = as.numeric(age_threshold),
    age_group_labels = c(
      reference = paste0("<=", age_threshold),
      comparison = paste0(">", age_threshold)
    ),
    age_group_counts = stats::setNames(
      c(sum(age_above_threshold == 0L), sum(age_above_threshold == 1L)),
      c(paste0("<=", age_threshold), paste0(">", age_threshold))
    ),
    outcome_name = outcome_name,
    measurement_columns = measurement_columns
  )
}
