
mira_create_ophthalmology_data <- function(
    n_per_arm = 50L,
    seed = 2026L,
    time_months = c(0, 3, 5, 12, 15)
) {

  if (length(n_per_arm) != 1L ||
      !is.numeric(n_per_arm) ||
      !is.finite(n_per_arm) ||
      n_per_arm != floor(n_per_arm) ||
      n_per_arm < 8L) {
    stop("`n_per_arm` must be one integer greater than or equal to 8.",
         call. = FALSE)
  }

  if (length(seed) != 1L ||
      !is.numeric(seed) ||
      !is.finite(seed) ||
      seed != floor(seed)) {
    stop("`seed` must be one finite integer.", call. = FALSE)
  }

  if (!is.numeric(time_months) ||
      length(time_months) < 2L ||
      any(!is.finite(time_months)) ||
      time_months[1L] != 0 ||
      is.unsorted(time_months, strictly = TRUE)) {
    stop(
      "`time_months` must contain at least two strictly increasing values, starting at 0.",
      call. = FALSE
    )
  }

  n_per_arm <- as.integer(n_per_arm)
  seed <- as.integer(seed)
  time_months <- as.numeric(time_months)

  clip <- function(x, lower, upper) {
    pmin(pmax(x, lower), upper)
  }

  arm_labels <- c("control", "ranibizumab", "aflibercept")
  treatment_labels <- c(
    control = "Bevacizumab",
    ranibizumab = "Ranibizumab",
    aflibercept = "Aflibercept"
  )

  # Approximate mean 12-month changes reported in Protocol T.
  bcva_gain_12m <- c(
    control = 9.7,
    ranibizumab = 11.2,
    aflibercept = 13.3
  )

  cmt_reduction_12m <- c(
    control = 101,
    ranibizumab = 147,
    aflibercept = 169
  )

  set.seed(seed)

  n <- length(arm_labels) * n_per_arm
  arm <- sample(rep(arm_labels, each = n_per_arm), size = n, replace = FALSE)

  age <- round(clip(stats::rnorm(n, mean = 61, sd = 10), 25, 85))
  gender <- ifelse(stats::runif(n) < 0.47, "Female", "Male")
  study_eye <- sample(c("OD", "OS"), n, replace = TRUE)

  diabetes_duration <- round(
    clip(stats::rnorm(n, mean = 17, sd = 10), 1, 45),
    1
  )
  hba1c <- round(clip(stats::rnorm(n, mean = 7.8, sd = 1.3), 5.2, 13), 1)

  # Baseline BCVA and CMT are moderately negatively associated: thicker
  # maculae tend to have worse visual acuity, without treating CMT as a
  # deterministic surrogate for vision.
  z_bcva <- stats::rnorm(n)
  z_cmt <- -0.35 * z_bcva + sqrt(1 - 0.35^2) * stats::rnorm(n)

  bcva_0 <- clip(
    64.8 - 0.08 * (age - 61) - 0.35 * (hba1c - 7.8) + 11.3 * z_bcva,
    24,
    78
  )

  cmt_0 <- clip(
    412 + 5 * (hba1c - 7.8) + 130 * z_cmt,
    250,
    750
  )

  # Correlated subject-level response heterogeneity. Better anatomic response
  # is only modestly associated with larger BCVA gain, consistent with the
  # limited individual-level correlation reported in Protocol T analyses.
  z_response_bcva <- stats::rnorm(n)
  z_response_cmt <-
    0.35 * z_response_bcva + sqrt(1 - 0.35^2) * stats::rnorm(n)

  bcva_gain_i <-
    unname(bcva_gain_12m[arm]) +
    0.18 * (64.8 - bcva_0) -
    0.04 * (age - 61) -
    0.20 * (hba1c - 7.8) +
    4.5 * z_response_bcva

  cmt_reduction_i <-
    unname(cmt_reduction_12m[arm]) +
    0.30 * (cmt_0 - 412) -
    3 * (hba1c - 7.8) +
    45 * z_response_cmt

  # Smooth early response, normalized to 1 at month 12. A small attenuation
  # after month 12 avoids imposing continued linear improvement.
  response_fraction <- function(month) {
    early <- (1 - exp(-month / 3.5)) / (1 - exp(-12 / 3.5))
    ifelse(month <= 12, early, exp(-0.005 * (month - 12)))
  }

  bcva <- matrix(NA_real_, nrow = n, ncol = length(time_months))
  cmt <- matrix(NA_real_, nrow = n, ncol = length(time_months))

  bcva[, 1L] <- bcva_0
  cmt[, 1L] <- cmt_0

  if (length(time_months) > 1L) {
    for (k in 2:length(time_months)) {
      fraction <- response_fraction(time_months[k])

      # Correlated visit-level deviations preserve realistic within-person
      # noise without making one outcome deterministically predict the other.
      visit_bcva <- stats::rnorm(n)
      visit_cmt <-
        0.25 * visit_bcva + sqrt(1 - 0.25^2) * stats::rnorm(n)

      bcva[, k] <-
        bcva_0 + fraction * bcva_gain_i + 2.2 * visit_bcva

      cmt[, k] <-
        cmt_0 - fraction * cmt_reduction_i - 18 * visit_cmt
    }
  }

  # ETDRS letters and OCT thickness are stored as integer-valued clinical
  # measurements. Bounds prevent impossible simulated values.
  bcva <- round(clip(bcva, 0, 100))
  cmt <- round(clip(cmt, 150, 800))

  colnames(bcva) <- paste0("BCVA_t", seq_along(time_months) - 1L)
  colnames(cmt) <- paste0("CMT_t", seq_along(time_months) - 1L)

  data <- data.frame(
    patient = sprintf("DME%04d", seq_len(n)),
    arm = factor(arm, levels = arm_labels),
    treatment = unname(treatment_labels[arm]),
    gender = factor(gender, levels = c("Female", "Male")),
    age = age,
    study_eye = study_eye,
    diabetes_duration_years = diabetes_duration,
    hba1c_percent = hba1c,
    synthetic = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  data <- cbind(data, as.data.frame(bcva), as.data.frame(cmt))

  attr(data, "time_months") <- time_months
  attr(data, "time_labels") <- ifelse(
    time_months == 0,
    "Baseline",
    paste("Month", format(time_months, trim = TRUE, scientific = FALSE))
  )
  attr(data, "reference_arm") <- "control"
  attr(data, "outcome_units") <- c(
    BCVA = "ETDRS letters",
    CMT = "micrometres"
  )
  attr(data, "simulation_truth") <- list(
    population = "Synthetic DME anti-VEGF cohort",
    treatment_mapping = treatment_labels,
    bcva_gain_12m = bcva_gain_12m,
    cmt_reduction_12m = cmt_reduction_12m,
    baseline_bcva_mean = 64.8,
    baseline_bcva_sd = 11.3,
    baseline_cmt_mean = 412,
    baseline_cmt_sd = 130
  )
  attr(data, "literature") <- c(
    protocol_t_one_year = "doi:10.1056/NEJMoa1414264",
    protocol_t_two_year = "doi:10.1016/j.ophtha.2016.02.022",
    bcva_cmt_association = "doi:10.1001/jamaophthalmol.2019.1963"
  )

  data
}


#' Select one ophthalmology outcome for the univariate MIRA model
#'
#' The current Stan model accepts one longitudinal outcome at a time. This
#' helper retains the metadata and exactly one set of `*_t*` columns.
#'
#' @param data A data frame created by [mira_create_ophthalmology_data()].
#' @param outcome Either `"BCVA"` or `"CMT"`.
#'
#' @return A data frame ready for `mira_info()` or `mira_prepare_data()`.
#'
#' @export
mira_select_ophthalmology_outcome <- function(
    data,
    outcome = c("BCVA", "CMT")
) {

  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.", call. = FALSE)
  }

  outcome <- match.arg(toupper(outcome), c("BCVA", "CMT"))
  outcome_columns <- grep(
    paste0("^", outcome, "_t[0-9]+$"),
    names(data),
    value = TRUE
  )

  if (length(outcome_columns) < 2L) {
    stop(
      "Could not find at least two columns named `",
      outcome,
      "_t0`, `",
      outcome,
      "_t1`, ... .",
      call. = FALSE
    )
  }

  all_outcome_columns <- grep("^(BCVA|CMT)_t[0-9]+$", names(data), value = TRUE)
  metadata_columns <- setdiff(names(data), all_outcome_columns)
  selected <- data[c(metadata_columns, outcome_columns)]

  custom_attributes <- setdiff(
    names(attributes(data)),
    c("names", "row.names", "class")
  )
  for (name in custom_attributes) {
    attr(selected, name) <- attr(data, name)
  }
  attr(selected, "outcome_name") <- outcome
  attr(selected, "outcome_unit") <- attr(data, "outcome_units")[[outcome]]

  selected
}



