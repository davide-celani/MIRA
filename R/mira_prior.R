.mira_prior_fields <- c(
  "baseline_prior_mean",
  "baseline_prior_sd",
  "beta_time_prior_mean",
  "beta_time_prior_sd",
  "tau_common_prior_rate",
  "beta_treatment_prior_sd",
  "tau_treatment_prior_rate",
  "arm_baseline_sd_prior_rate",
  "gender_baseline_prior_sd",
  "beta_gender_prior_sd",
  "tau_gender_prior_rate",
  "age_baseline_prior_sd",
  "beta_age_prior_sd",
  "tau_age_prior_rate",
  "sigma_intercept_prior_rate",
  "sigma_slope_prior_rate",
  "sigma_prior_rate",
  "nu_prior_shape",
  "nu_prior_rate"
)


.mira_prior_user_fields <- c(
  "baseline_mean",
  "baseline_sd",
  "beta_time_mean",
  "beta_time_sd",
  "tau_common_rate",
  "beta_treatment_sd",
  "tau_treatment_rate",
  "arm_baseline_sd_rate",
  "gender_baseline_sd",
  "beta_gender_sd",
  "tau_gender_rate",
  "age_baseline_sd",
  "beta_age_sd",
  "tau_age_rate",
  "sigma_intercept_rate",
  "sigma_slope_rate",
  "sigma_rate",
  "nu_shape",
  "nu_rate"
)


.mira_prior_field_map <- stats::setNames(
  .mira_prior_user_fields,
  .mira_prior_fields
)


.mira_normalize_outcome <- function(
    outcome,
    allow_auto = TRUE,
    unknown_as_generic = FALSE
) {
  if (length(outcome) != 1L || is.na(outcome) || !nzchar(outcome)) {
    stop("`outcome` must be one non-empty character value.", call. = FALSE)
  }

  value <- toupper(trimws(as.character(outcome)))

  if (allow_auto && value == "AUTO") return("auto")
  if (value %in% c("BCVA", "VA", "ETDRS")) return("BCVA")
  if (value %in% c("CMT", "CST", "CSFT")) return("CMT")
  if (value %in% c("GENERIC", "OTHER")) return("generic")

  if (unknown_as_generic) return("generic")

  stop(
    "Unsupported `outcome`: ", outcome,
    ". Use 'auto', 'BCVA', 'CMT', or 'generic'.",
    call. = FALSE
  )
}


.mira_normalize_informativeness <- function(informativeness) {
  if (length(informativeness) < 1L ||
      is.na(informativeness[1L]) ||
      !nzchar(informativeness[1L])) {
    stop(
      "`informativeness` must be one non-empty character value.",
      call. = FALSE
    )
  }

  value <- tolower(trimws(as.character(informativeness[1L])))

  if (value %in% c("standard", "normal", "default")) return("standard")
  if (value %in% c("weak", "less_informative", "less-informative")) {
    return("weak")
  }
  if (value %in% c("informative", "strong")) return("informative")
  if (value == "custom") return("custom")

  stop(
    "Unsupported `informativeness`: ", informativeness[1L],
    ". Use 'standard', 'weak', 'informative', or 'custom'.",
    call. = FALSE
  )
}


.mira_apply_informativeness <- function(values, informativeness) {
  if (informativeness == "standard" || informativeness == "custom") {
    return(values)
  }

  # Weak priors have wider Normal distributions and larger expected scales.
  # Informative priors have narrower Normal distributions and smaller scales.
  spread_factor <- switch(
    informativeness,
    weak = 2.5,
    informative = 0.5
  )

  normal_sd_fields <- c(
    "baseline_sd",
    "beta_time_sd",
    "beta_treatment_sd",
    "gender_baseline_sd",
    "beta_gender_sd",
    "age_baseline_sd",
    "beta_age_sd"
  )

  exponential_rate_fields <- c(
    "tau_common_rate",
    "tau_treatment_rate",
    "arm_baseline_sd_rate",
    "tau_gender_rate",
    "tau_age_rate",
    "sigma_intercept_rate",
    "sigma_slope_rate",
    "sigma_rate"
  )

  values[normal_sd_fields] <- lapply(
    values[normal_sd_fields],
    function(x) x * spread_factor
  )

  # For Exponential(rate), the mean is 1 / rate. Dividing the rate therefore
  # makes the corresponding scale prior wider.
  values[exponential_rate_fields] <- lapply(
    values[exponential_rate_fields],
    function(x) x / spread_factor
  )

  if (informativeness == "weak") {
    values$nu_shape <- 2
    values$nu_rate <- 0.1
  } else {
    values$nu_shape <- 5
    values$nu_rate <- 0.5
  }

  values
}


.mira_normalize_custom_prior <- function(custom_prior) {
  if (is.null(custom_prior)) return(list())

  if (!is.list(custom_prior) ||
      is.null(names(custom_prior)) ||
      anyNA(names(custom_prior)) ||
      any(!nzchar(names(custom_prior))) ||
      anyDuplicated(names(custom_prior))) {
    stop(
      "`custom_prior` must be a named list with unique, non-empty names.",
      call. = FALSE
    )
  }

  normalized_names <- names(custom_prior)
  stan_names <- normalized_names %in% names(.mira_prior_field_map)
  normalized_names[stan_names] <- unname(
    .mira_prior_field_map[normalized_names[stan_names]]
  )

  unknown_names <- setdiff(normalized_names, .mira_prior_user_fields)
  if (length(unknown_names) > 0L) {
    stop(
      "Unknown fields in `custom_prior`: ",
      paste(unknown_names, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  if (anyDuplicated(normalized_names)) {
    stop(
      "`custom_prior` specifies the same parameter more than once.",
      call. = FALSE
    )
  }

  names(custom_prior) <- normalized_names
  custom_prior
}


#' Create BCVA, CMT, or custom priors for the MIRA model
#'
#' `outcome` and `informativeness` are independent choices. For example,
#' BCVA can use weak, standard, informative, or custom priors, and the same is
#' true for CMT. The returned object contains exactly the 19 prior fields used
#' by the current Stan model, plus R-side metadata.
#'
#' The BCVA and CMT centres are calibrated to DRCR.net Protocol T summaries.
#' BCVA is expressed in ETDRS letters and CMT in micrometres. The CMT profile
#' uses central subfield thickness as a practical proxy; OCT device and
#' segmentation conventions can require different values.
#'
#' @param stan_data Data prepared by [mira_prepare_data()].
#' @param outcome Outcome scale: `"auto"`, `"BCVA"`, `"CMT"`, or
#'   `"generic"`. With `"auto"`, `stan_data$outcome_name` is used.
#' @param informativeness Prior strength:
#'   * `"weak"`: broad, less informative priors;
#'   * `"standard"`: clinically scaled default priors (`"normal"` is an
#'     accepted alias);
#'   * `"informative"`: narrower priors;
#'   * `"custom"`: standard priors with user-supplied replacements.
#' @param time_unit Unit used by `stan_data$time_value`. Clinical slope and
#'   continuous-time RW1 priors are defined per month and converted to this
#'   unit.
#' @param custom_prior Optional named list of replacements. It may use the
#'   short argument names, such as `baseline_mean`, or the exact Stan names,
#'   such as `baseline_prior_mean`. Unspecified fields inherit the selected
#'   outcome's standard prior.
#' @param profile Optional legacy profile. Supported for compatibility with
#'   existing code: `"auto"`, `"bcva"`, `"cmt"`, `"default"`, `"weak"`,
#'   `"regularized"`, `"informative"`, and `"custom"`.
#' @param baseline_mean Prior mean for the reference-arm baseline.
#' @param baseline_sd Prior SD for the reference-arm baseline.
#' @param beta_time_mean Prior mean for the reference-arm time slope.
#' @param beta_time_sd Prior SD for the reference-arm time slope.
#' @param tau_common_rate Exponential rate for the common RW1 scale.
#' @param beta_treatment_sd Prior SD for treatment slope differences.
#' @param tau_treatment_rate Exponential rate for treatment RW1 scales.
#' @param arm_baseline_sd_rate Exponential rate for baseline arm imbalance SD.
#' @param gender_baseline_sd Prior SD for the gender baseline difference.
#' @param beta_gender_sd Prior SD for the gender slope difference.
#' @param tau_gender_rate Exponential rate for the gender RW1 scale.
#' @param age_baseline_sd Prior SD for the age-group baseline difference.
#' @param beta_age_sd Prior SD for the age-group slope difference.
#' @param tau_age_rate Exponential rate for the age-group RW1 scale.
#' @param sigma_intercept_rate Exponential rate for random-intercept SD.
#' @param sigma_slope_rate Exponential rate for random-slope SD.
#' @param sigma_rate Exponential rate for the residual Student-t scale.
#' @param nu_shape Shape of the Gamma prior for Student-t degrees of freedom.
#' @param nu_rate Rate of the Gamma prior for Student-t degrees of freedom.
#'
#' @return An object of class `mira_prior`.
#'
#' @export
mira_prior <- function(
    stan_data,
    outcome = c("auto", "BCVA", "CMT", "generic"),
    informativeness = c("standard", "weak", "informative", "custom"),
    time_unit = c("months", "years", "weeks", "days"),
    custom_prior = NULL,
    profile = NULL,
    baseline_mean = NULL,
    baseline_sd = NULL,
    beta_time_mean = NULL,
    beta_time_sd = NULL,
    tau_common_rate = NULL,
    beta_treatment_sd = NULL,
    tau_treatment_rate = NULL,
    arm_baseline_sd_rate = NULL,
    gender_baseline_sd = NULL,
    beta_gender_sd = NULL,
    tau_gender_rate = NULL,
    age_baseline_sd = NULL,
    beta_age_sd = NULL,
    tau_age_rate = NULL,
    sigma_intercept_rate = NULL,
    sigma_slope_rate = NULL,
    sigma_rate = NULL,
    nu_shape = NULL,
    nu_rate = NULL
) {

  outcome_was_missing <- missing(outcome)
  informativeness_was_missing <- missing(informativeness)

  if (!is.list(stan_data)) {
    stop("`stan_data` must be a list.", call. = FALSE)
  }

  required_data <- c("y", "time", "time_value", "mean_y", "sd_y")
  missing_data <- setdiff(required_data, names(stan_data))

  if (length(missing_data) > 0L) {
    stop(
      "Missing data required to construct priors: ",
      paste(missing_data, collapse = ", "),
      call. = FALSE
    )
  }

  y <- stan_data$y
  time <- stan_data$time
  time_value <- stan_data$time_value
  mean_y <- stan_data$mean_y
  sd_y <- stan_data$sd_y

  if (!is.numeric(y) || length(y) < 1L || any(!is.finite(y))) {
    stop("`stan_data$y` must contain finite numeric values.", call. = FALSE)
  }

  if (!is.numeric(time) ||
      length(time) != length(y) ||
      any(!is.finite(time)) ||
      any(time != floor(time)) ||
      any(time < 1L)) {
    stop("`stan_data$time` must contain one positive integer per observation.",
         call. = FALSE)
  }

  if (!is.numeric(time_value) ||
      length(time_value) < 2L ||
      any(!is.finite(time_value)) ||
      is.unsorted(time_value, strictly = TRUE)) {
    stop("`stan_data$time_value` must be finite and strictly increasing.",
         call. = FALSE)
  }

  if (any(time > length(time_value))) {
    stop("A time index exceeds `length(stan_data$time_value)`.", call. = FALSE)
  }

  if (length(mean_y) != 1L || !is.numeric(mean_y) || !is.finite(mean_y)) {
    stop("`stan_data$mean_y` must be one finite numeric value.", call. = FALSE)
  }

  if (length(sd_y) != 1L ||
      !is.numeric(sd_y) ||
      !is.finite(sd_y) ||
      sd_y <= 0) {
    stop("`stan_data$sd_y` must be one positive finite numeric value.",
         call. = FALSE)
  }

  legacy_profile <- NULL
  if (!is.null(profile)) {
    if (length(profile) != 1L || is.na(profile) || !nzchar(profile)) {
      stop("`profile` must be NULL or one non-empty character value.",
           call. = FALSE)
    }

    legacy_profile <- tolower(trimws(as.character(profile)))
    allowed_profiles <- c(
      "auto", "bcva", "cmt", "default", "weak",
      "regularized", "informative", "custom"
    )

    if (!legacy_profile %in% allowed_profiles) {
      stop(
        "Unsupported legacy `profile`: ", profile, ".",
        call. = FALSE
      )
    }

    if (legacy_profile %in% c("bcva", "cmt")) {
      legacy_outcome <- toupper(legacy_profile)
      if (!outcome_was_missing &&
          .mira_normalize_outcome(outcome[1L]) != legacy_outcome) {
        stop("Legacy `profile` conflicts with `outcome`.", call. = FALSE)
      }
      outcome <- legacy_outcome
      if (informativeness_was_missing) informativeness <- "standard"
    }

    legacy_information <- switch(
      legacy_profile,
      default = "standard",
      weak = "weak",
      regularized = "informative",
      informative = "informative",
      custom = "custom",
      NULL
    )

    if (!is.null(legacy_information)) {
      if (!informativeness_was_missing &&
          .mira_normalize_informativeness(informativeness) !=
          legacy_information) {
        stop("Legacy `profile` conflicts with `informativeness`.",
             call. = FALSE)
      }
      informativeness <- legacy_information
    }
  }

  requested_outcome <- .mira_normalize_outcome(outcome[1L])
  informativeness <- .mira_normalize_informativeness(informativeness)
  time_unit <- match.arg(
    tolower(time_unit),
    c("months", "years", "weeks", "days")
  )

  data_outcome <- if (!is.null(stan_data$outcome_name)) {
    .mira_normalize_outcome(
      stan_data$outcome_name,
      allow_auto = FALSE,
      unknown_as_generic = TRUE
    )
  } else {
    "generic"
  }

  resolved_outcome <- if (requested_outcome == "auto") {
    data_outcome
  } else {
    requested_outcome
  }

  if (requested_outcome != "auto" &&
      data_outcome != "generic" &&
      requested_outcome != data_outcome) {
    stop(
      "`outcome = '", requested_outcome,
      "'` does not match `stan_data$outcome_name = '",
      data_outcome, "'`.",
      call. = FALSE
    )
  }

  y <- as.numeric(y)
  time <- as.integer(time)
  time_value <- as.numeric(time_value)
  baseline_y <- y[time == 1L]

  if (length(baseline_y) < 1L) {
    stop("No baseline observations (`time == 1`) were found.", call. = FALSE)
  }

  baseline_center <- mean(baseline_y)
  time_span <- max(time_value) - min(time_value)
  outcome_scale <- max(as.numeric(sd_y), 1e-8)
  slope_scale <- max(outcome_scale / time_span, 1e-8)
  rw_scale <- max(outcome_scale / sqrt(time_span), 1e-8)

  months_per_unit <- c(
    months = 1,
    years = 12,
    weeks = 12 / 52.1775,
    days = 12 / 365.25
  )[[time_unit]]

  slope_unit <- months_per_unit
  rw_unit <- sqrt(months_per_unit)

  if (resolved_outcome == "BCVA") {
    values <- list(
      baseline_mean = 65,
      baseline_sd = 15,
      beta_time_mean = 0.8 * slope_unit,
      beta_time_sd = 0.6 * slope_unit,
      tau_common_rate = 1 / (1.5 * rw_unit),
      beta_treatment_sd = 0.5 * slope_unit,
      tau_treatment_rate = 1 / (1.0 * rw_unit),
      arm_baseline_sd_rate = 1 / 4,
      gender_baseline_sd = 3,
      beta_gender_sd = 0.25 * slope_unit,
      tau_gender_rate = 1 / (1.0 * rw_unit),
      age_baseline_sd = 5,
      beta_age_sd = 0.25 * slope_unit,
      tau_age_rate = 1 / (1.0 * rw_unit),
      sigma_intercept_rate = 1 / 10,
      sigma_slope_rate = 1 / (0.6 * slope_unit),
      sigma_rate = 1 / 3,
      nu_shape = 3,
      nu_rate = 0.2
    )
  } else if (resolved_outcome == "CMT") {
    values <- list(
      baseline_mean = 412,
      baseline_sd = 130,
      beta_time_mean = -8 * slope_unit,
      beta_time_sd = 8 * slope_unit,
      tau_common_rate = 1 / (20 * rw_unit),
      beta_treatment_sd = 6 * slope_unit,
      tau_treatment_rate = 1 / (15 * rw_unit),
      arm_baseline_sd_rate = 1 / 30,
      gender_baseline_sd = 25,
      beta_gender_sd = 3 * slope_unit,
      tau_gender_rate = 1 / (15 * rw_unit),
      age_baseline_sd = 30,
      beta_age_sd = 3 * slope_unit,
      tau_age_rate = 1 / (15 * rw_unit),
      sigma_intercept_rate = 1 / 60,
      sigma_slope_rate = 1 / (6 * slope_unit),
      sigma_rate = 1 / 25,
      nu_shape = 3,
      nu_rate = 0.2
    )
  } else {
    values <- list(
      baseline_mean = baseline_center,
      baseline_sd = 2 * outcome_scale,
      beta_time_mean = 0,
      beta_time_sd = 2 * slope_scale,
      tau_common_rate = 1 / rw_scale,
      beta_treatment_sd = 2 * slope_scale,
      tau_treatment_rate = 1 / rw_scale,
      arm_baseline_sd_rate = 1 / outcome_scale,
      gender_baseline_sd = 2 * outcome_scale,
      beta_gender_sd = 2 * slope_scale,
      tau_gender_rate = 1 / rw_scale,
      age_baseline_sd = 2 * outcome_scale,
      beta_age_sd = 2 * slope_scale,
      tau_age_rate = 1 / rw_scale,
      sigma_intercept_rate = 1 / outcome_scale,
      sigma_slope_rate = 1 / slope_scale,
      sigma_rate = 1 / outcome_scale,
      nu_shape = 3,
      nu_rate = 0.2
    )
  }

  values <- .mira_apply_informativeness(values, informativeness)

  custom_values <- .mira_normalize_custom_prior(custom_prior)
  for (name in names(custom_values)) {
    values[[name]] <- custom_values[[name]]
  }

  direct_overrides <- list(
    baseline_mean = baseline_mean,
    baseline_sd = baseline_sd,
    beta_time_mean = beta_time_mean,
    beta_time_sd = beta_time_sd,
    tau_common_rate = tau_common_rate,
    beta_treatment_sd = beta_treatment_sd,
    tau_treatment_rate = tau_treatment_rate,
    arm_baseline_sd_rate = arm_baseline_sd_rate,
    gender_baseline_sd = gender_baseline_sd,
    beta_gender_sd = beta_gender_sd,
    tau_gender_rate = tau_gender_rate,
    age_baseline_sd = age_baseline_sd,
    beta_age_sd = beta_age_sd,
    tau_age_rate = tau_age_rate,
    sigma_intercept_rate = sigma_intercept_rate,
    sigma_slope_rate = sigma_slope_rate,
    sigma_rate = sigma_rate,
    nu_shape = nu_shape,
    nu_rate = nu_rate
  )

  direct_names <- names(direct_overrides)[
    !vapply(direct_overrides, is.null, logical(1))
  ]

  for (name in direct_names) {
    values[[name]] <- direct_overrides[[name]]
  }

  customized_fields <- unique(c(names(custom_values), direct_names))
  if (informativeness == "custom" && length(customized_fields) == 0L) {
    stop(
      "With `informativeness = 'custom'`, supply `custom_prior` or at least one explicit prior argument.",
      call. = FALSE
    )
  }

  finite_names <- c("baseline_mean", "beta_time_mean")
  positive_names <- setdiff(.mira_prior_user_fields, finite_names)

  for (name in finite_names) {
    value <- values[[name]]
    if (length(value) != 1L ||
        !is.numeric(value) ||
        !is.finite(value)) {
      stop("`", name, "` must be one finite numeric value.", call. = FALSE)
    }
  }

  for (name in positive_names) {
    value <- values[[name]]
    if (length(value) != 1L ||
        !is.numeric(value) ||
        !is.finite(value) ||
        value <= 0) {
      stop("`", name, "` must be one positive finite numeric value.",
           call. = FALSE)
    }
  }

  prior <- list(
    baseline_prior_mean = as.numeric(values$baseline_mean),
    baseline_prior_sd = as.numeric(values$baseline_sd),
    beta_time_prior_mean = as.numeric(values$beta_time_mean),
    beta_time_prior_sd = as.numeric(values$beta_time_sd),
    tau_common_prior_rate = as.numeric(values$tau_common_rate),
    beta_treatment_prior_sd = as.numeric(values$beta_treatment_sd),
    tau_treatment_prior_rate = as.numeric(values$tau_treatment_rate),
    arm_baseline_sd_prior_rate = as.numeric(values$arm_baseline_sd_rate),
    gender_baseline_prior_sd = as.numeric(values$gender_baseline_sd),
    beta_gender_prior_sd = as.numeric(values$beta_gender_sd),
    tau_gender_prior_rate = as.numeric(values$tau_gender_rate),
    age_baseline_prior_sd = as.numeric(values$age_baseline_sd),
    beta_age_prior_sd = as.numeric(values$beta_age_sd),
    tau_age_prior_rate = as.numeric(values$tau_age_rate),
    sigma_intercept_prior_rate = as.numeric(values$sigma_intercept_rate),
    sigma_slope_prior_rate = as.numeric(values$sigma_slope_rate),
    sigma_prior_rate = as.numeric(values$sigma_rate),
    nu_prior_shape = as.numeric(values$nu_shape),
    nu_prior_rate = as.numeric(values$nu_rate),
    profile = paste(tolower(resolved_outcome), informativeness, sep = "_"),
    outcome = resolved_outcome,
    informativeness = informativeness,
    time_unit = time_unit,
    customized_fields = customized_fields,
    legacy_profile = legacy_profile,
    reference_scales = list(
      baseline_center = baseline_center,
      outcome_scale = outcome_scale,
      slope_scale = slope_scale,
      rw_scale = rw_scale,
      time_span = time_span,
      months_per_time_unit = months_per_unit
    ),
    literature = if (resolved_outcome %in% c("BCVA", "CMT")) c(
      protocol_t_one_year = "doi:10.1056/NEJMoa1414264",
      protocol_t_two_year = "doi:10.1016/j.ophtha.2016.02.022",
      bcva_cmt_association = "doi:10.1001/jamaophthalmol.2019.1963"
    ) else NULL
  )

  class(prior) <- "mira_prior"
  mira_validate_prior(prior)
  prior
}


#' Create BCVA priors
#'
#' @param stan_data Data prepared from BCVA columns.
#' @param informativeness `"standard"`, `"weak"`, `"informative"`, or
#'   `"custom"`.
#' @param time_unit Unit used by `stan_data$time_value`.
#' @param ... Optional prior replacements passed to [mira_prior()].
#'
#' @return A `mira_prior` object.
#'
#' @export
mira_prior_bcva <- function(
    stan_data,
    informativeness = "standard",
    time_unit = "months",
    ...
) {
  mira_prior(
    stan_data = stan_data,
    outcome = "BCVA",
    informativeness = informativeness,
    time_unit = time_unit,
    ...
  )
}


#' Create CMT priors
#'
#' @param stan_data Data prepared from CMT columns.
#' @param informativeness `"standard"`, `"weak"`, `"informative"`, or
#'   `"custom"`.
#' @param time_unit Unit used by `stan_data$time_value`.
#' @param ... Optional prior replacements passed to [mira_prior()].
#'
#' @return A `mira_prior` object.
#'
#' @export
mira_prior_cmt <- function(
    stan_data,
    informativeness = "standard",
    time_unit = "months",
    ...
) {
  mira_prior(
    stan_data = stan_data,
    outcome = "CMT",
    informativeness = informativeness,
    time_unit = time_unit,
    ...
  )
}


#' Validate a MIRA prior specification
#'
#' @param prior A `mira_prior` object or complete named list.
#'
#' @return Invisibly returns `TRUE`.
#'
#' @export
mira_validate_prior <- function(prior) {
  if (!is.list(prior)) {
    stop("`prior` must be a list.", call. = FALSE)
  }

  if (is.null(names(prior)) ||
      anyNA(names(prior)) ||
      any(!nzchar(names(prior))) ||
      anyDuplicated(names(prior))) {
    stop("`prior` must have unique, non-empty field names.", call. = FALSE)
  }

  missing_fields <- setdiff(.mira_prior_fields, names(prior))
  if (length(missing_fields) > 0L) {
    stop(
      "Invalid MIRA prior. Missing fields: ",
      paste(missing_fields, collapse = ", "),
      call. = FALSE
    )
  }

  finite_fields <- c("baseline_prior_mean", "beta_time_prior_mean")
  positive_fields <- setdiff(.mira_prior_fields, finite_fields)

  for (name in finite_fields) {
    value <- prior[[name]]
    if (length(value) != 1L ||
        !is.numeric(value) ||
        !is.finite(value)) {
      stop("Invalid prior field `", name, "`: expected one finite number.",
           call. = FALSE)
    }
  }

  for (name in positive_fields) {
    value <- prior[[name]]
    if (length(value) != 1L ||
        !is.numeric(value) ||
        !is.finite(value) ||
        value <= 0) {
      stop(
        "Invalid prior field `", name,
        "`: expected one positive finite number.",
        call. = FALSE
      )
    }
  }

  invisible(TRUE)
}


#' Convert MIRA priors to Stan data
#'
#' @param prior A `mira_prior` object or complete named list.
#'
#' @return A list containing exactly the 19 prior fields declared in Stan.
#'
#' @export
mira_prior_stan_data <- function(prior) {
  mira_validate_prior(prior)
  stan_prior <- lapply(.mira_prior_fields, function(name) prior[[name]])
  names(stan_prior) <- .mira_prior_fields
  stan_prior
}


#' Print a MIRA prior specification
#'
#' @param x A `mira_prior` object.
#' @param ... Additional arguments (unused).
#'
#' @export
print.mira_prior <- function(x, ...) {
  mira_validate_prior(x)

  cat("\nMIRA prior specification\n")
  cat("========================\n")
  if (!is.null(x$outcome)) cat("Outcome: ", x$outcome, "\n", sep = "")
  if (!is.null(x$informativeness)) {
    cat("Informativeness: ", x$informativeness, "\n", sep = "")
  }
  if (!is.null(x$time_unit)) cat("Time unit: ", x$time_unit, "\n", sep = "")
  if (!is.null(x$customized_fields) && length(x$customized_fields) > 0L) {
    cat(
      "Customized: ", paste(x$customized_fields, collapse = ", "), "\n",
      sep = ""
    )
  }
  cat("\n")
  cat(
    "Baseline ~ Normal(", x$baseline_prior_mean, ", ",
    x$baseline_prior_sd, ")\n",
    sep = ""
  )
  cat(
    "Reference slope ~ Normal(", x$beta_time_prior_mean, ", ",
    x$beta_time_prior_sd, ")\n",
    sep = ""
  )
  cat(
    "Treatment slope SD: ", x$beta_treatment_prior_sd, "\n",
    "Residual scale ~ Exponential(", x$sigma_prior_rate, ")\n",
    "Student-t nu ~ Gamma(", x$nu_prior_shape, ", ",
    x$nu_prior_rate, "), truncated at 2\n",
    sep = ""
  )

  invisible(x)
}
