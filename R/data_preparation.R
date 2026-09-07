# Internal helpers are deliberately kept in this file so that
# `mira_prepare_data()` also works when the script is sourced outside an R
# package.  The Stan model uses the numeric likelihood id, while users work
# with readable names.
.mira_prepare_outcome_key <- function(x) {
  value <- toupper(trimws(as.character(x)))

  if (value %in% c("BCVA", "VA", "ETDRS")) return("BCVA")
  if (value %in% c("CMT", "CST", "CSFT")) return("CMT")
  "generic"
}

.mira_prepare_likelihood <- function(likelihood, outcome) {
  if (length(likelihood) != 1L || is.na(likelihood) || !nzchar(likelihood)) {
    stop("`likelihood` must be one non-empty character value.", call. = FALSE)
  }

  value <- tolower(trimws(as.character(likelihood)))
  aliases <- c(
    student_t = "student_t",
    student = "student_t",
    robust = "student_t",
    gaussian = "gaussian",
    normal = "gaussian",
    lognormal = "lognormal",
    log_normal = "lognormal"
  )

  if (value == "auto") {
    if (identical(outcome, "CMT")) return("lognormal")
    if (identical(outcome, "BCVA")) return("student_t")

    warning(
      "The outcome is not a recognized BCVA/CMT variable; using the robust ",
      "identity-scale Student-t model. Set `likelihood` explicitly after ",
      "checking the outcome support and empirical distribution.",
      call. = FALSE
    )
    return("student_t")
  }

  if (!(value %in% names(aliases))) {
    stop(
      "Unsupported `likelihood`: ", likelihood,
      ". Use 'auto', 'student_t', 'gaussian', or 'lognormal'.",
      call. = FALSE
    )
  }

  unname(aliases[[value]])
}

.mira_sample_skewness <- function(x) {
  x <- as.numeric(x)
  n <- length(x)
  s <- stats::sd(x)
  if (n < 3L || !is.finite(s) || s <= 0) return(NA_real_)
  n / ((n - 1) * (n - 2)) * sum(((x - mean(x)) / s)^3)
}

#' Prepare longitudinal data for the MIRA treatment model
#'
#' Prepares and validates longitudinal data for the outcome-adaptive MIRA
#' Bayesian model. By default BCVA uses a censored robust Student-t model on
#' its letter-score scale, whereas CMT uses a log-normal model with a log link.
#' Gaussian identity and explicit family overrides are available for
#' sensitivity analyses.
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
#' @param outcome Outcome to select when the data contain one or more
#'   longitudinal series. `"auto"` is allowed only when there is one outcome
#'   prefix. Recognized aliases are BCVA/VA/ETDRS and CMT/CST/CSFT; any other
#'   value is treated as an exact column prefix and as a generic outcome.
#' @param likelihood Observation model. `"auto"` selects Student-t for BCVA,
#'   log-normal for CMT, and a warned Student-t fallback for an unknown
#'   continuous outcome. Explicit options are `"student_t"`, `"gaussian"`,
#'   and `"lognormal"`.
#' @param outcome_bounds Optional numeric vector `c(lower, upper)` defining
#'   observable bounds. Use `NA` for an absent endpoint. BCVA defaults to
#'   `c(0, 100)` and is modeled as a censored continuous score. CMT defaults
#'   to the positive support of the log-normal distribution.
#' @param meaningful_change Positive numeric value giving the prior mean of
#'   the minimum clinically important difference (MCID). Together with
#'   `meaningful_change_sd`, it defines a Gamma prior with exactly these
#'   moments. When omitted, BCVA defaults to 5 ETDRS letters. CMT and generic
#'   outcomes require an explicit value because no universal threshold exists.
#' @param meaningful_change_sd Positive numeric value giving the prior SD of
#'   the uncertain MCID. This must be externally specified; the outcome data
#'   should not be used to identify MCID uncertainty.
#' @param direction Direction of clinical improvement. `"auto"` selects
#'   higher-is-better for BCVA and lower-is-better for CMT. Otherwise use `1`
#'   or `"higher"`, or `-1` or `"lower"`.
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
#' @return A named list containing the variables required by Stan plus
#'   outcome-family metadata and empirical diagnostics. All clinical changes
#'   and MCID values remain on the natural outcome scale even when the latent
#'   trajectory is modeled on the log scale.
#'
#' @export
mira_prepare_data <- function(
    data,
    time_value,
    outcome = "auto",
    likelihood = "auto",
    outcome_bounds = NULL,
    meaningful_change = NULL,
    meaningful_change_sd,
    direction = "auto",
    meaningful_between_arm_difference = NULL,
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

  all_measurement_columns <- names(data)[
    grepl("_t[0-9]+$", names(data))
  ]

  if (length(all_measurement_columns) < 2) {
    stop(
      paste0(
        "At least two longitudinal measurement columns are required. ",
        "Expected columns named `<outcome>_t0`, `<outcome>_t1`, ..., ",
        "`<outcome>_tK`."
      ),
      call. = FALSE
    )
  }

  if (length(outcome) != 1L || is.na(outcome) || !nzchar(outcome)) {
    stop("`outcome` must be one non-empty character value.", call. = FALSE)
  }

  all_outcome_names <- sub("_t[0-9]+$", "", all_measurement_columns)
  available_outcomes <- unique(all_outcome_names)
  requested_outcome <- trimws(as.character(outcome))

  if (tolower(requested_outcome) == "auto") {
    if (length(available_outcomes) != 1L) {
      stop(
        paste0(
          "Multiple longitudinal outcomes were detected: ",
          paste(available_outcomes, collapse = ", "),
          ". Select one with `outcome`, for example `outcome = 'BCVA'` ",
          "or `outcome = 'CMT'`."
        ),
        call. = FALSE
      )
    }
    outcome_name <- available_outcomes[[1L]]
  } else {
    exact_match <- available_outcomes[
      tolower(available_outcomes) == tolower(requested_outcome)
    ]
    requested_key <- .mira_prepare_outcome_key(requested_outcome)

    if (length(exact_match) == 1L) {
      outcome_name <- exact_match[[1L]]
    } else if (requested_key %in% c("BCVA", "CMT")) {
      alias_match <- available_outcomes[
        vapply(
          available_outcomes,
          function(x) identical(.mira_prepare_outcome_key(x), requested_key),
          logical(1L)
        )
      ]

      if (length(alias_match) != 1L) {
        stop(
          "`outcome = '", requested_outcome,
          "'` matches ", length(alias_match),
          " available prefixes. Available outcomes: ",
          paste(available_outcomes, collapse = ", "),
          ". Use the exact prefix to disambiguate.",
          call. = FALSE
        )
      }
      outcome_name <- alias_match[[1L]]
    } else if (tolower(requested_outcome) == "generic" &&
               length(available_outcomes) == 1L) {
      outcome_name <- available_outcomes[[1L]]
    } else {
      stop(
        "Could not find longitudinal columns for `outcome = '",
        requested_outcome, "'`. Available outcomes: ",
        paste(available_outcomes, collapse = ", "), ".",
        call. = FALSE
      )
    }
  }

  measurement_columns <- all_measurement_columns[
    tolower(all_outcome_names) == tolower(outcome_name)
  ]

  if (length(measurement_columns) < 2L) {
    stop(
      "The selected outcome `", outcome_name,
      "` must have at least two measurement occasions.",
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

  outcome_class <- .mira_prepare_outcome_key(outcome_name)
  if (tolower(requested_outcome) == "generic") {
    outcome_class <- "generic"
  } else if (tolower(requested_outcome) != "auto" &&
             .mira_prepare_outcome_key(requested_outcome) %in% c("BCVA", "CMT")) {
    outcome_class <- .mira_prepare_outcome_key(requested_outcome)
  }

  likelihood <- .mira_prepare_likelihood(likelihood, outcome_class)
  if (outcome_class == "BCVA" && likelihood == "lognormal") {
    stop(
      "A log-normal model is incompatible with bounded ETDRS letter scores. ",
      "Use `student_t` (primary) or `gaussian` (sensitivity).",
      call. = FALSE
    )
  }
  likelihood_id <- c(student_t = 1L, gaussian = 2L, lognormal = 3L)[[likelihood]]
  modeling_scale <- if (likelihood == "lognormal") "log" else "identity"

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

    if (direction_key == "auto") {
      if (outcome_class == "BCVA") {
        direction <- 1L
      } else if (outcome_class == "CMT") {
        direction <- -1L
        warning(
          "`direction = 'auto'` assumes that lower CMT is anatomically ",
          "better. Verify this for the disease and time horizon (for example, ",
          "retinal thinning/atrophy can invalidate a monotone lower-is-better ",
          "interpretation); otherwise set `direction` explicitly.",
          call. = FALSE
        )
      } else {
        stop(
          "`direction = 'auto'` is available only for recognized BCVA or ",
          "CMT outcomes. Specify 'higher' or 'lower' for this outcome.",
          call. = FALSE
        )
      }
    } else if (direction_key %in% c("higher", "higher_better", "increase", "increasing")) {
      direction <- 1L
    } else if (direction_key %in% c("lower", "lower_better", "decrease", "decreasing")) {
      direction <- -1L
    } else {
      stop(
        "Character `direction` must be `auto`, `higher`, or `lower` ",
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

  if (is.null(meaningful_change)) {
    if (outcome_class == "BCVA") {
      meaningful_change <- 5
      warning(
        "`meaningful_change` was not supplied; using 5 ETDRS letters as a ",
        "candidate responder threshold. This is not treated as a universal ",
        "BCVA MCID: override it when the protocol defines another threshold.",
        call. = FALSE
      )
    } else {
      stop(
        "`meaningful_change` must be supplied for ", outcome_class,
        ". In particular, CMT has no universal MCID; provide an externally ",
        "justified threshold in micrometres.",
        call. = FALSE
      )
    }
  }

  if (
    length(meaningful_change) != 1 ||
    !is.numeric(meaningful_change) ||
    !is.finite(meaningful_change) ||
    meaningful_change <= 0
  ) {
    stop(
      "`meaningful_change` must be one finite positive numeric value.",
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

  if (is.null(meaningful_between_arm_difference)) {
    meaningful_between_arm_difference <- meaningful_change
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

  non_finite_outcomes <- measurement_columns[
    vapply(
      data[measurement_columns],
      function(x) any(!is.finite(x)),
      logical(1)
    )
  ]

  if (length(non_finite_outcomes) > 0) {
    stop(
      paste0(
        "The following longitudinal measurement columns contain missing or ",
        "non-finite values: ",
        paste(non_finite_outcomes, collapse = ", "),
        ". Missing values are not currently supported by the MIRA Stan model."
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

  # Observable support. Identity-scale models use censoring at supplied
  # endpoints; the log-normal family additionally enforces strict positivity.
  if (is.null(outcome_bounds)) {
    outcome_bounds <- if (outcome_class == "BCVA") {
      c(0, 100)
    } else if (outcome_class == "CMT") {
      c(0, NA_real_)
    } else {
      c(NA_real_, NA_real_)
    }
  }

  if (is.logical(outcome_bounds) && length(outcome_bounds) == 2L &&
      all(is.na(outcome_bounds))) {
    outcome_bounds <- as.numeric(outcome_bounds)
  }

  if (!is.numeric(outcome_bounds) || length(outcome_bounds) != 2L) {
    stop(
      "`outcome_bounds` must be NULL or numeric `c(lower, upper)`; use NA ",
      "for a missing endpoint.",
      call. = FALSE
    )
  }

  finite_or_missing <- is.na(outcome_bounds) | is.finite(outcome_bounds)
  if (!all(finite_or_missing)) {
    stop("Finite bounds or NA must be used in `outcome_bounds`.", call. = FALSE)
  }

  has_lower_bound <- as.integer(!is.na(outcome_bounds[[1L]]))
  has_upper_bound <- as.integer(!is.na(outcome_bounds[[2L]]))
  outcome_lower_bound <- if (has_lower_bound == 1L) {
    as.numeric(outcome_bounds[[1L]])
  } else {
    0
  }
  outcome_upper_bound <- if (has_upper_bound == 1L) {
    as.numeric(outcome_bounds[[2L]])
  } else {
    0
  }

  if (has_lower_bound == 1L && has_upper_bound == 1L &&
      outcome_lower_bound >= outcome_upper_bound) {
    stop("The lower outcome bound must be smaller than the upper bound.", call. = FALSE)
  }

  if (has_lower_bound == 1L && any(y < outcome_lower_bound)) {
    stop(
      "Observed `", outcome_name, "` values fall below the declared lower ",
      "bound of ", outcome_lower_bound, ".",
      call. = FALSE
    )
  }

  if (has_upper_bound == 1L && any(y > outcome_upper_bound)) {
    stop(
      "Observed `", outcome_name, "` values exceed the declared upper ",
      "bound of ", outcome_upper_bound, ".",
      call. = FALSE
    )
  }

  if (likelihood == "lognormal" && any(y <= 0)) {
    stop(
      "The log-normal model requires every observed outcome to be strictly ",
      "positive. Zeros require a hurdle/zero-inflated model, not an automatic ",
      "log transformation.",
      call. = FALSE
    )
  }

  if (outcome_class == "BCVA" && any(abs(y - round(y)) > sqrt(.Machine$double.eps))) {
    stop(
      "The automatic BCVA branch expects integer ETDRS letter scores. For ",
      "logMAR or another visual-acuity scale, use `outcome = 'generic'` (or ",
      "a non-BCVA column prefix) and explicitly set `likelihood`, `direction`, ",
      "`outcome_bounds`, and the clinical thresholds.",
      call. = FALSE
    )
  }

  # ============================================================
  # OUTCOME MOMENTS FOR R-SIDE PRIOR/INITIALIZATION HELPERS
  # ============================================================

  mean_y <- mean(y)
  sd_y <- stats::sd(y)

  y_model <- if (likelihood == "lognormal") log(y) else y
  mean_model_y <- mean(y_model)
  sd_model_y <- stats::sd(y_model)

  if (!is.finite(mean_y) || !is.finite(sd_y) || sd_y <= 0) {
    stop(
      "The outcome standard deviation must be positive and finite.",
      call. = FALSE
    )
  }

  if (!is.finite(mean_model_y) || !is.finite(sd_model_y) || sd_model_y <= 0) {
    stop(
      "The outcome standard deviation on the selected modeling scale must ",
      "be positive and finite.",
      call. = FALSE
    )
  }

  lower_hits <- if (has_lower_bound == 1L) {
    sum(y == outcome_lower_bound)
  } else {
    0L
  }
  upper_hits <- if (has_upper_bound == 1L) {
    sum(y == outcome_upper_bound)
  } else {
    0L
  }

  outcome_diagnostics <- list(
    n = N,
    min = min(y),
    q05 = unname(stats::quantile(y, 0.05)),
    median = stats::median(y),
    mean = mean_y,
    q95 = unname(stats::quantile(y, 0.95)),
    max = max(y),
    sd = sd_y,
    coefficient_of_variation = if (mean_y != 0) sd_y / abs(mean_y) else NA_real_,
    skewness = .mira_sample_skewness(y),
    log_skewness = if (all(y > 0)) .mira_sample_skewness(log(y)) else NA_real_,
    lower_boundary_count = as.integer(lower_hits),
    upper_boundary_count = as.integer(upper_hits),
    lower_boundary_fraction = lower_hits / N,
    upper_boundary_fraction = upper_hits / N,
    likelihood = likelihood,
    modeling_scale = modeling_scale
  )

  # ============================================================
  # RETURN
  # ============================================================

  list(
    # Stan model data
    N = as.integer(N),
    S = as.integer(S),
    K = as.integer(K),
    G = as.integer(G),
    likelihood_id = as.integer(likelihood_id),
    has_lower_bound = as.integer(has_lower_bound),
    outcome_lower_bound = as.numeric(outcome_lower_bound),
    has_upper_bound = as.integer(has_upper_bound),
    outcome_upper_bound = as.numeric(outcome_upper_bound),
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
    mean_model_y = as.numeric(mean_model_y),
    sd_model_y = as.numeric(sd_model_y),
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
    subject_labels = as.character(data$patient),
    outcome_name = outcome_name,
    outcome = outcome_class,
    likelihood = likelihood,
    modeling_scale = modeling_scale,
    outcome_bounds = c(
      lower = if (has_lower_bound == 1L) outcome_lower_bound else NA_real_,
      upper = if (has_upper_bound == 1L) outcome_upper_bound else NA_real_
    ),
    boundary_strategy = if (likelihood == "lognormal" &&
                            has_upper_bound == 0L &&
                            (has_lower_bound == 0L || outcome_lower_bound <= 0)) {
      "strictly positive log-normal support"
    } else if (likelihood == "lognormal") {
      "strictly positive log-normal support with endpoint censoring"
    } else if (has_lower_bound == 1L || has_upper_bound == 1L) {
      "endpoint censoring"
    } else {
      "none"
    },
    change_scale = "absolute natural-outcome units",
    outcome_diagnostics = outcome_diagnostics,
    available_outcomes = available_outcomes,
    measurement_columns = measurement_columns
  )
}
