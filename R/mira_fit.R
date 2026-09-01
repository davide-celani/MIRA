#' Fit MIRA longitudinal treatment model
#'
#' Fits the MIRA longitudinal Student-t mixed-effects model with
#' treatment-, gender-, and age-threshold-specific trajectories using CmdStan.
#'
#' @param stan_data Data prepared for the MIRA Stan model. The list may
#'   contain additional R-side metadata (for example `mean_y`, `sd_y`, or
#'   `arm_labels`, gender labels, or the age threshold); only variables required
#'   by Stan are passed to CmdStan. The current model requires subject-level
#'   binary indicators `male` (0 = Female reference, 1 = Male) and
#'   `age_above_threshold` (0 = age <= threshold, 1 = age > threshold).
#' @param prior A `mira_prior` object, or a named list containing the Stan
#'   prior fields required by the current model, including dedicated priors
#'   for gender- and age-threshold trajectories. If `NULL`, `mira_prior()`
#'   is called with profile = "default", unless all prior fields are already
#'   present in `stan_data`. A non-NULL `prior` takes precedence over prior
#'   fields embedded in `stan_data`.
#' @param chains Number of MCMC chains.
#' @param parallel_chains Number of parallel chains.
#' @param iter_warmup Number of warmup iterations.
#' @param iter_sampling Number of sampling iterations.
#' @param seed Random seed.
#' @param refresh Number of iterations between progress messages.
#' @param verbose Logical. If TRUE (default), print a compact MIRA fit report
#'   automatically after successful sampling. The fitted CmdStanMCMC object is
#'   still returned invisibly and can be assigned normally.
#' @param stan_file Optional path to the Stan file. If `NULL`, MIRA first
#'   looks for `inst/stan/gaussian_longitudinal.stan`, then legacy MIRA Stan
#'   filenames for backwards compatibility.
#'
#' @return A CmdStanMCMC object.
#'
#' @export
mira_fit <- function(
    stan_data,
    prior = NULL,
    chains = 4,
    parallel_chains = chains,
    iter_warmup = 1000,
    iter_sampling = 3000,
    seed = 123,
    refresh = 100,
    stan_file = NULL,
    verbose = TRUE
) {

  # ------------------------------------------------------------
  # Data object
  # ------------------------------------------------------------

  if (!is.list(stan_data)) {
    stop("`stan_data` must be a list.", call. = FALSE)
  }

  positive_integer <- function(x, name, allow_zero = FALSE) {
    lower <- if (allow_zero) 0L else 1L
    ok <- is.numeric(x) && length(x) == 1L && is.finite(x) &&
      x >= lower && x <= .Machine$integer.max && x == as.integer(x)
    if (!ok) {
      stop(
        "`", name, "` must be one integer ",
        if (allow_zero) ">= 0." else ">= 1.",
        call. = FALSE
      )
    }
    invisible(TRUE)
  }

  positive_integer(chains, "chains")
  positive_integer(parallel_chains, "parallel_chains")
  positive_integer(iter_warmup, "iter_warmup", allow_zero = TRUE)
  positive_integer(iter_sampling, "iter_sampling")
  positive_integer(refresh, "refresh", allow_zero = TRUE)
  positive_integer(seed, "seed", allow_zero = TRUE)

  if (parallel_chains > chains) {
    stop("`parallel_chains` cannot be larger than `chains`.", call. = FALSE)
  }

  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("`verbose` must be TRUE or FALSE.", call. = FALSE)
  }

  model_data_names <- c(
    "N", "S", "K", "G",
    "y", "subject", "time", "arm",
    "male", "age_above_threshold", "time_value",
    "direction",
    "mcid_prior_mean", "mcid_prior_sd",
    "meaningful_between_arm_difference"
  )

  missing_data <- setdiff(model_data_names, names(stan_data))

  if (length(missing_data) > 0) {
    stop(
      "Missing data required by the new MIRA Stan model: ",
      paste(missing_data, collapse = ", "),
      call. = FALSE
    )
  }

  # ------------------------------------------------------------
  # Dimensions and indices
  # ------------------------------------------------------------

  scalar_integer_names <- c("N", "S", "K", "G")

  for (nm in scalar_integer_names) {
    x <- stan_data[[nm]]
    if (length(x) != 1 || !is.numeric(x) || !is.finite(x) ||
        x < 0 || x > .Machine$integer.max || x != floor(x)) {
      stop("`", nm, "` must be one finite integer.", call. = FALSE)
    }
  }

  if (stan_data$N < 1) stop("`N` must be >= 1.", call. = FALSE)
  if (stan_data$S < 1) stop("`S` must be >= 1.", call. = FALSE)
  if (stan_data$K < 2) stop("`K` must be >= 2.", call. = FALSE)
  if (stan_data$G < 2) stop("`G` must be >= 2 for the current treatment model.", call. = FALSE)

  if (!is.numeric(stan_data$y) ||
      length(stan_data$y) != stan_data$N ||
      any(!is.finite(stan_data$y))) {
    stop("`y` must contain exactly N finite numeric values.", call. = FALSE)
  }

  if (length(stan_data$subject) != stan_data$N) {
    stop("`length(subject)` must equal `N`.", call. = FALSE)
  }

  if (length(stan_data$time) != stan_data$N) {
    stop("`length(time)` must equal `N`.", call. = FALSE)
  }

  if (length(stan_data$arm) != stan_data$S) {
    stop("`length(arm)` must equal `S`.", call. = FALSE)
  }

  if (!is.numeric(stan_data$subject) ||
      any(!is.finite(stan_data$subject)) ||
      any(stan_data$subject != floor(stan_data$subject)) ||
      any(stan_data$subject < 1) ||
      any(stan_data$subject > stan_data$S)) {
    stop("`subject` must contain integers between 1 and S.", call. = FALSE)
  }

  if (!is.numeric(stan_data$time) ||
      any(!is.finite(stan_data$time)) ||
      any(stan_data$time != floor(stan_data$time)) ||
      any(stan_data$time < 1) ||
      any(stan_data$time > stan_data$K)) {
    stop("`time` must contain integers between 1 and K.", call. = FALSE)
  }

  if (!is.numeric(stan_data$arm) ||
      any(!is.finite(stan_data$arm)) ||
      any(stan_data$arm != floor(stan_data$arm)) ||
      any(stan_data$arm < 1) ||
      any(stan_data$arm > stan_data$G)) {
    stop("`arm` must contain integers between 1 and G.", call. = FALSE)
  }

  # ------------------------------------------------------------
  # Subject-level gender and age-group indicators
  # ------------------------------------------------------------

  validate_binary_subject_indicator <- function(x, name) {
    if (!is.numeric(x) || length(x) != stan_data$S) {
      stop(
        "`", name, "` must be a numeric/integer vector of length S.",
        call. = FALSE
      )
    }

    if (any(!is.finite(x)) ||
        any(x != floor(x)) ||
        any(!x %in% c(0, 1))) {
      stop(
        "`", name, "` must contain exactly S binary integer values (0/1).",
        call. = FALSE
      )
    }
  }

  validate_binary_subject_indicator(stan_data$male, "male")
  validate_binary_subject_indicator(
    stan_data$age_above_threshold,
    "age_above_threshold"
  )

  if (length(unique(stan_data$male)) < 2) {
    warning(
      "`male` contains only one observed category; the gender-by-time effect ",
      "will be weakly/non-identified by these data.",
      call. = FALSE
    )
  }

  if (length(unique(stan_data$age_above_threshold)) < 2) {
    warning(
      "`age_above_threshold` contains only one observed group; the age-group-by-time ",
      "effect will be weakly/non-identified by these data.",
      call. = FALSE
    )
  }

  if (!is.numeric(stan_data$time_value) ||
      length(stan_data$time_value) != stan_data$K ||
      any(!is.finite(stan_data$time_value)) ||
      anyDuplicated(stan_data$time_value) ||
      is.unsorted(stan_data$time_value, strictly = TRUE)) {
    stop(
      "`time_value` must contain exactly K finite, strictly increasing values.",
      call. = FALSE
    )
  }

  # ------------------------------------------------------------
  # Clinical inputs
  # ------------------------------------------------------------

  if (length(stan_data$direction) != 1 ||
      !is.numeric(stan_data$direction) ||
      !is.finite(stan_data$direction) ||
      !(stan_data$direction %in% c(-1, 1))) {
    stop("`direction` must be exactly +1 or -1.", call. = FALSE)
  }

  positive_scalar <- function(x, name, allow_zero = FALSE) {
    ok <- length(x) == 1 && is.numeric(x) && is.finite(x)
    ok <- ok && if (allow_zero) x >= 0 else x > 0
    if (!ok) {
      comparator <- if (allow_zero) ">= 0" else "> 0"
      stop("`", name, "` must be one finite numeric value ", comparator, ".", call. = FALSE)
    }
  }

  positive_scalar(stan_data$mcid_prior_mean, "mcid_prior_mean")
  positive_scalar(stan_data$mcid_prior_sd, "mcid_prior_sd")
  positive_scalar(
    stan_data$meaningful_between_arm_difference,
    "meaningful_between_arm_difference",
    allow_zero = TRUE
  )

  # ------------------------------------------------------------
  # Prior data required by Stan
  # ------------------------------------------------------------

  prior_names <- c(
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

  embedded_prior_names <- intersect(prior_names, names(stan_data))

  # An explicit `prior` argument takes precedence over any embedded fields.
  if (!is.null(prior)) {

    mira_validate_prior(prior)
    stan_prior_data <- mira_prior_stan_data(prior)
    prior_source <- "argument"

  } else if (all(prior_names %in% names(stan_data))) {

    stan_prior_data <- stan_data[prior_names]
    prior_source <- "embedded in stan_data"

  } else {

    if (length(embedded_prior_names) > 0L) {
      stop(
        "`stan_data` contains an incomplete embedded prior specification. ",
        "Missing fields: ",
        paste(setdiff(prior_names, embedded_prior_names), collapse = ", "),
        ". Supply a complete `prior` argument or remove the partial fields.",
        call. = FALSE
      )
    }

    prior <- mira_prior(
      stan_data,
      profile = "default"
    )
    stan_prior_data <- mira_prior_stan_data(prior)
    prior_source <- "automatic default"
  }

  missing_prior <- setdiff(prior_names, names(stan_prior_data))

  if (length(missing_prior) > 0) {
    stop(
      "The prior specification is not compatible with the new Stan model. ",
      "Missing Stan prior fields: ",
      paste(missing_prior, collapse = ", "),
      ". Update `mira_prior()` / `mira_prior_stan_data()` so that gender and age priors are included, or pass a complete named prior list.",
      call. = FALSE
    )
  }

  stan_prior_data <- stan_prior_data[prior_names]

  finite_prior_names <- c("baseline_prior_mean", "beta_time_prior_mean")
  for (nm in finite_prior_names) {
    x <- stan_prior_data[[nm]]
    if (length(x) != 1 || !is.numeric(x) || !is.finite(x)) {
      stop("`", nm, "` must be one finite numeric value.", call. = FALSE)
    }
  }

  positive_prior_names <- setdiff(prior_names, finite_prior_names)
  for (nm in positive_prior_names) {
    positive_scalar(stan_prior_data[[nm]], nm)
  }

  # Pass only objects declared in the Stan data block. This lets stan_data
  # safely retain R-side metadata such as mean_y, sd_y and arm labels.
  sampling_data <- c(
    stan_data[model_data_names],
    stan_prior_data
  )

  sampling_data$N <- as.integer(sampling_data$N)
  sampling_data$S <- as.integer(sampling_data$S)
  sampling_data$K <- as.integer(sampling_data$K)
  sampling_data$G <- as.integer(sampling_data$G)
  sampling_data$subject <- as.integer(sampling_data$subject)
  sampling_data$time <- as.integer(sampling_data$time)
  sampling_data$arm <- as.integer(sampling_data$arm)
  sampling_data$male <- as.integer(sampling_data$male)
  sampling_data$age_above_threshold <- as.integer(
    sampling_data$age_above_threshold
  )
  sampling_data$direction <- as.integer(sampling_data$direction)

  # ------------------------------------------------------------
  # Stan model file
  # ------------------------------------------------------------

  if (is.null(stan_file)) {

    candidates <- c(
      "gaussian_longitudinal.stan",
      "gaussian_longitudinal_gender_age.stan",
      "mira_longitudinal.stan"
    )

    candidate_paths <- vapply(
      candidates,
      function(x) system.file("stan", x, package = "MIRA"),
      character(1)
    )

    existing <- candidate_paths[nzchar(candidate_paths)]

    # Also support sourcing the function from the package project before
    # installation, when system.file() cannot yet resolve the package path.
    if (length(existing) == 0) {
      local_candidates <- file.path("inst", "stan", candidates)
      existing <- local_candidates[file.exists(local_candidates)]
    }

    if (length(existing) == 0) {
      stop(
        "Could not find the MIRA Stan model in `inst/stan`. Expected ",
        "`gaussian_longitudinal.stan` (preferred), ",
        "`gaussian_longitudinal_gender_age.stan`, or `mira_longitudinal.stan`.",
        call. = FALSE
      )
    }

    stan_file <- existing[[1]]

  } else {

    if (length(stan_file) != 1 || !is.character(stan_file) || !nzchar(stan_file)) {
      stop("`stan_file` must be NULL or one non-empty character path.", call. = FALSE)
    }

    if (!file.exists(stan_file)) {
      packaged_file <- system.file("stan", stan_file, package = "MIRA")
      if (!nzchar(packaged_file)) {
        stop("Could not find Stan file: ", stan_file, call. = FALSE)
      }
      stan_file <- packaged_file
    }
  }

  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    stop(
      "Package `cmdstanr` is required to compile and fit the MIRA model.",
      call. = FALSE
    )
  }

  message("Compiling MIRA Stan model: ", basename(stan_file))

  compile_started <- Sys.time()

  model <- cmdstanr::cmdstan_model(
    stan_file,
    quiet = TRUE
  )

  compile_elapsed_seconds <- as.numeric(
    difftime(Sys.time(), compile_started, units = "secs")
  )

  # ------------------------------------------------------------
  # Initial values
  # ------------------------------------------------------------

  y <- as.numeric(stan_data$y)
  baseline_y <- y[stan_data$time == 1]

  mean_y <- if (!is.null(stan_data$mean_y)) {
    as.numeric(stan_data$mean_y)[1]
  } else {
    mean(y)
  }

  sd_y <- if (!is.null(stan_data$sd_y)) {
    as.numeric(stan_data$sd_y)[1]
  } else {
    stats::sd(y)
  }

  if (!is.finite(mean_y)) mean_y <- mean(y)
  if (!is.finite(sd_y) || sd_y <= 0) sd_y <- max(abs(mean_y) * 0.1, 1)

  baseline_init <- if (length(baseline_y) > 0) mean(baseline_y) else mean_y
  elapsed <- max(stan_data$time_value) - min(stan_data$time_value)
  elapsed_safe <- max(elapsed, 1e-6)

  slope_scale <- max(sd_y / elapsed_safe, 0.01)
  rw_scale <- max(sd_y / sqrt(elapsed_safe) / 10, 0.01)

  init <- function() {
    list(
      baseline_mean = baseline_init,
      beta_time = 0,
      z_common_step = rep(0, stan_data$K - 1),
      tau_common = rw_scale,
      beta_treatment = rep(0, stan_data$G - 1),
      z_treatment_step = matrix(
        0,
        nrow = stan_data$G - 1,
        ncol = stan_data$K - 1
      ),
      tau_treatment = rep(rw_scale, stan_data$G - 1),
      z_arm_baseline = rep(0, stan_data$G - 1),
      arm_baseline_sd = max(sd_y / 10, 0.01),

      # Gender-by-time trajectory: Male - Female.
      gender_baseline_effect = 0,
      beta_gender_time = 0,
      z_gender_step = rep(0, stan_data$K - 1),
      tau_gender = rw_scale,

      # Age-group-by-time trajectory: age > threshold - age <= threshold.
      age_baseline_effect = 0,
      beta_age_time = 0,
      z_age_step = rep(0, stan_data$K - 1),
      tau_age = rw_scale,

      z_subject = matrix(
        0,
        nrow = 2,
        ncol = stan_data$S
      ),
      sigma_subject = c(
        max(sd_y / 2, 0.1),
        max(slope_scale / 2, 0.01)
      ),
      L_subject = diag(2),
      sigma = max(sd_y / 2, 0.1),
      nu = 10,
      mcid = max(stan_data$mcid_prior_mean, 1e-6)
    )
  }

  # ------------------------------------------------------------
  # Sampling
  # ------------------------------------------------------------

  message("Sampling posterior...")

  sampling_started <- Sys.time()

  fit <- model$sample(
    data = sampling_data,
    chains = chains,
    parallel_chains = parallel_chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    seed = seed,
    init = init,
    refresh = refresh
  )

  sampling_elapsed_seconds <- as.numeric(
    difftime(Sys.time(), sampling_started, units = "secs")
  )

  # ------------------------------------------------------------
  # Attach lightweight MIRA metadata to the CmdStanMCMC object.
  #
  # The object remains a CmdStanMCMC object, so existing code such as
  # mira_summary(fit, ...) and all CmdStanR $methods continue to work.
  # ------------------------------------------------------------

  arm_labels <- if (!is.null(stan_data$arm_labels) &&
                    length(stan_data$arm_labels) == stan_data$G) {
    as.character(stan_data$arm_labels)
  } else {
    paste0("arm_", seq_len(stan_data$G))
  }

  arm_counts <- table(
    factor(
      stan_data$arm,
      levels = seq_len(stan_data$G),
      labels = arm_labels
    )
  )

  gender_counts <- c(
    Female = sum(stan_data$male == 0L),
    Male = sum(stan_data$male == 1L)
  )

  age_threshold <- if (!is.null(stan_data$age_threshold)) {
    as.numeric(stan_data$age_threshold)[1]
  } else {
    NA_real_
  }

  age_counts <- c(
    at_or_below_threshold = sum(stan_data$age_above_threshold == 0L),
    above_threshold = sum(stan_data$age_above_threshold == 1L)
  )

  prior_profile <- if (!is.null(prior) && !is.null(prior$profile)) {
    as.character(prior$profile)[1]
  } else if (identical(prior_source, "embedded in stan_data")) {
    "embedded in stan_data"
  } else {
    "named/custom prior list"
  }

  mira_fit_info <- list(
    model_file = basename(stan_file),
    model_path = tryCatch(
      normalizePath(stan_file, winslash = "/", mustWork = FALSE),
      error = function(e) stan_file
    ),
    n_observations = as.integer(stan_data$N),
    n_subjects = as.integer(stan_data$S),
    n_time_points = as.integer(stan_data$K),
    n_arms = as.integer(stan_data$G),
    arm_labels = arm_labels,
    arm_counts = arm_counts,
    gender_counts = gender_counts,
    age_threshold = age_threshold,
    age_counts = age_counts,
    time_value = as.numeric(stan_data$time_value),
    direction = as.integer(stan_data$direction),
    direction_interpretation = if (stan_data$direction == 1L) {
      "higher outcome = better"
    } else {
      "lower outcome = better"
    },
    mcid_prior_mean = as.numeric(stan_data$mcid_prior_mean),
    mcid_prior_sd = as.numeric(stan_data$mcid_prior_sd),
    meaningful_between_arm_difference =
      as.numeric(stan_data$meaningful_between_arm_difference),
    prior_profile = prior_profile,
    prior_source = prior_source,
    prior_data = stan_prior_data,
    chains = as.integer(chains),
    parallel_chains = as.integer(parallel_chains),
    iter_warmup = as.integer(iter_warmup),
    iter_sampling = as.integer(iter_sampling),
    post_warmup_draws = as.integer(chains * iter_sampling),
    seed = as.integer(seed),
    compile_elapsed_seconds = compile_elapsed_seconds,
    sampling_elapsed_seconds = sampling_elapsed_seconds
  )

  attr(fit, "mira_fit_info") <- mira_fit_info
  class(fit) <- unique(c("mira_fit", class(fit)))

  if (isTRUE(verbose)) {
    print(fit)
  }

  invisible(fit)
}


#' Print a fitted MIRA model
#'
#' Compact technical report for a fitted MIRA CmdStanMCMC object. The report
#' focuses on model setup, sampling quality and core posterior parameters.
#' Detailed clinical estimands are intentionally left to [mira_summary()].
#'
#' @param x A fitted object returned by [mira_fit()].
#' @param digits Number of digits used for numeric output.
#' @param parameters Logical; print core posterior parameter summaries.
#' @param diagnostics Logical; print MCMC diagnostics.
#' @param guide Logical; print a guide to useful CmdStanR methods and the next
#'   MIRA analysis step.
#' @param ... Additional arguments (currently ignored).
#'
#' @return Invisibly returns `x`.
#'
#' @export
print.mira_fit <- function(
    x,
    digits = 3,
    parameters = TRUE,
    diagnostics = TRUE,
    guide = TRUE,
    ...
) {

  if (!inherits(x, "CmdStanMCMC")) {
    stop("`x` must inherit from CmdStanMCMC.", call. = FALSE)
  }

  if (!is.numeric(digits) || length(digits) != 1L ||
      is.na(digits) || digits < 0) {
    stop("`digits` must be one integer >= 0.", call. = FALSE)
  }

  digits <- as.integer(digits)

  info <- attr(x, "mira_fit_info")

  if (is.null(info)) {
    info <- list()
  }

  line <- function(char = "-", n = 92L) {
    cat(strrep(char, n), "\n", sep = "")
  }

  section <- function(title) {
    cat("\n", title, "\n", sep = "")
    line()
  }

  fmt_num <- function(z, d = digits) {
    if (length(z) == 0L || is.na(z) || !is.finite(z)) {
      return("NA")
    }
    formatC(z, format = "f", digits = d)
  }

  fmt_time <- function(seconds) {
    if (length(seconds) == 0L || is.na(seconds) || !is.finite(seconds)) {
      return("NA")
    }
    if (seconds < 60) {
      return(sprintf("%.1f sec", seconds))
    }
    if (seconds < 3600) {
      return(sprintf("%.1f min", seconds / 60))
    }
    sprintf("%.2f h", seconds / 3600)
  }

  fmt_flag <- function(ok) {
    if (length(ok) == 0L || is.na(ok)) return("UNKNOWN")
    if (isTRUE(ok)) "OK" else "CHECK"
  }

  safe_metadata <- tryCatch(x$metadata(), error = function(e) NULL)

  # ============================================================
  # HEADER
  # ============================================================

  cat("\n")
  line("=")
  cat("MIRA BAYESIAN LONGITUDINAL MODEL - FIT REPORT\n")
  line("=")

  # ============================================================
  # MODEL / DATA
  # ============================================================

  section("MODEL AND DATA")

  if (!is.null(info$model_file)) {
    cat(sprintf("Stan model: %s\n", info$model_file))
  }

  if (!is.null(info$n_subjects)) {
    cat(sprintf(
      "Subjects: %d | Observations: %d | Timepoints: %d | Arms: %d\n",
      info$n_subjects,
      info$n_observations,
      info$n_time_points,
      info$n_arms
    ))
  }

  if (!is.null(info$time_value)) {
    cat("Time values: ", paste(info$time_value, collapse = ", "), "\n", sep = "")
  }

  if (!is.null(info$direction_interpretation)) {
    cat("Clinical direction: ", info$direction_interpretation, "\n", sep = "")
  }

  if (!is.null(info$prior_profile)) {
    cat("Prior profile: ", info$prior_profile, "\n", sep = "")
  }

  if (!is.null(info$prior_source)) {
    cat("Prior source: ", info$prior_source, "\n", sep = "")
  }

  if (!is.null(info$mcid_prior_mean)) {
    cat(sprintf(
      "MCID Gamma prior: mean %s | SD %s | Between-arm meaningful threshold: %s\n",
      fmt_num(info$mcid_prior_mean),
      fmt_num(info$mcid_prior_sd),
      fmt_num(info$meaningful_between_arm_difference)
    ))
  }

  if (!is.null(info$arm_counts)) {
    cat("\nSubjects by treatment arm:\n")
    print(info$arm_counts)
  }

  if (!is.null(info$gender_counts)) {
    cat("\nSubjects by gender (Female reference):\n")
    print(info$gender_counts)
  }

  if (!is.null(info$age_counts)) {
    if (!is.null(info$age_threshold) && is.finite(info$age_threshold)) {
      names(info$age_counts) <- c(
        paste0("age <= ", info$age_threshold),
        paste0("age > ", info$age_threshold)
      )
    }
    cat("\nSubjects by age-threshold group:\n")
    print(info$age_counts)
  }

  # ============================================================
  # SAMPLING
  # ============================================================

  section("SAMPLING")

  if (!is.null(info$chains)) {
    cat(sprintf(
      "Chains: %d | Parallel chains: %d | Warmup: %d | Sampling: %d per chain\n",
      info$chains,
      info$parallel_chains,
      info$iter_warmup,
      info$iter_sampling
    ))
    cat(sprintf(
      "Post-warmup draws: %d | Seed: %d\n",
      info$post_warmup_draws,
      info$seed
    ))
    cat(sprintf(
      "Compilation: %s | Sampling: %s\n",
      fmt_time(info$compile_elapsed_seconds),
      fmt_time(info$sampling_elapsed_seconds)
    ))
  } else if (!is.null(safe_metadata)) {
    cat("CmdStanR metadata are available via fit$metadata().\n")
  }

  # ============================================================
  # DIAGNOSTICS
  # ============================================================

  if (isTRUE(diagnostics)) {
    section("MCMC DIAGNOSTICS")

    diagnostic_table <- tryCatch(
      x$diagnostic_summary(),
      error = function(e) NULL
    )

    divergences <- NA_real_
    treedepth_hits <- NA_real_
    min_ebfmi <- NA_real_

    if (!is.null(diagnostic_table)) {
      diagnostic_df <- as.data.frame(diagnostic_table)
      nm <- names(diagnostic_df)

      div_col <- grep("diverg", nm, ignore.case = TRUE, value = TRUE)[1]
      tree_col <- grep("treedepth", nm, ignore.case = TRUE, value = TRUE)[1]
      ebfmi_col <- grep("ebfmi|bfmi", nm, ignore.case = TRUE, value = TRUE)[1]

      if (!is.na(div_col) && nzchar(div_col)) {
        divergences <- sum(diagnostic_df[[div_col]], na.rm = TRUE)
      }

      if (!is.na(tree_col) && nzchar(tree_col)) {
        treedepth_hits <- sum(diagnostic_df[[tree_col]], na.rm = TRUE)
      }

      if (!is.na(ebfmi_col) && nzchar(ebfmi_col)) {
        z <- diagnostic_df[[ebfmi_col]]
        z <- z[is.finite(z)]
        if (length(z) > 0L) min_ebfmi <- min(z)
      }
    }

    model_param_names <- NULL
    if (!is.null(safe_metadata) && !is.null(safe_metadata$model_params)) {
      model_param_names <- safe_metadata$model_params
    }

    diagnostic_parameter_summary <- tryCatch(
      {
        if (!is.null(model_param_names) && length(model_param_names) > 0L) {
          x$summary(variables = model_param_names)
        } else {
          NULL
        }
      },
      error = function(e) NULL
    )

    max_rhat <- NA_real_
    min_ess_bulk <- NA_real_
    min_ess_tail <- NA_real_

    if (!is.null(diagnostic_parameter_summary)) {
      ds <- as.data.frame(diagnostic_parameter_summary)

      if ("rhat" %in% names(ds)) {
        z <- ds$rhat[is.finite(ds$rhat)]
        if (length(z) > 0L) max_rhat <- max(z)
      }

      if ("ess_bulk" %in% names(ds)) {
        z <- ds$ess_bulk[is.finite(ds$ess_bulk)]
        if (length(z) > 0L) min_ess_bulk <- min(z)
      }

      if ("ess_tail" %in% names(ds)) {
        z <- ds$ess_tail[is.finite(ds$ess_tail)]
        if (length(z) > 0L) min_ess_tail <- min(z)
      }
    }

    cat(sprintf(
      "Divergences: %s [%s] | Max-treedepth hits: %s [%s]\n",
      ifelse(is.finite(divergences), as.character(divergences), "NA"),
      fmt_flag(if (is.finite(divergences)) divergences == 0 else NA),
      ifelse(is.finite(treedepth_hits), as.character(treedepth_hits), "NA"),
      fmt_flag(if (is.finite(treedepth_hits)) treedepth_hits == 0 else NA)
    ))

    cat(sprintf(
      "Max R-hat: %s [%s] | Min bulk ESS: %s [%s] | Min tail ESS: %s [%s]\n",
      fmt_num(max_rhat),
      fmt_flag(if (is.finite(max_rhat)) max_rhat < 1.01 else NA),
      fmt_num(min_ess_bulk, 0L),
      fmt_flag(if (is.finite(min_ess_bulk)) min_ess_bulk >= 400 else NA),
      fmt_num(min_ess_tail, 0L),
      fmt_flag(if (is.finite(min_ess_tail)) min_ess_tail >= 400 else NA)
    ))

    cat(sprintf(
      "Minimum E-BFMI: %s [%s]\n",
      fmt_num(min_ebfmi),
      fmt_flag(if (is.finite(min_ebfmi)) min_ebfmi > 0.30 else NA)
    ))

    checks <- c(
      if (is.finite(divergences)) divergences == 0 else NA,
      if (is.finite(treedepth_hits)) treedepth_hits == 0 else NA,
      if (is.finite(max_rhat)) max_rhat < 1.01 else NA,
      if (is.finite(min_ess_bulk)) min_ess_bulk >= 400 else NA,
      if (is.finite(min_ess_tail)) min_ess_tail >= 400 else NA,
      if (is.finite(min_ebfmi)) min_ebfmi > 0.30 else NA
    )

    known_checks <- checks[!is.na(checks)]

    overall <- if (length(known_checks) == 0L) {
      "DIAGNOSTICS UNAVAILABLE"
    } else if (all(known_checks)) {
      "NO OBVIOUS MCMC PROBLEMS DETECTED"
    } else {
      "CHECK MCMC DIAGNOSTICS BEFORE INTERPRETATION"
    }

    cat("Overall: ", overall, "\n", sep = "")
  }

  # ============================================================
  # CORE PARAMETERS
  # ============================================================

  if (isTRUE(parameters)) {
    section("CORE POSTERIOR PARAMETERS")

    core_candidates <- c(
      "baseline_mean",
      "beta_time",
      "tau_common",
      "beta_treatment",
      "tau_treatment",
      "arm_baseline_sd",
      "gender_baseline_effect",
      "beta_gender_time",
      "tau_gender",
      "age_baseline_effect",
      "beta_age_time",
      "tau_age",
      "sigma_subject",
      "sigma",
      "nu",
      "mcid"
    )

    available_core <- core_candidates

    if (!is.null(safe_metadata) && !is.null(safe_metadata$model_params)) {
      available_core <- intersect(core_candidates, safe_metadata$model_params)
    }

    core_summary <- tryCatch(
      x$summary(variables = available_core),
      error = function(e) NULL
    )

    if (is.null(core_summary) || nrow(core_summary) == 0L) {
      cat("Core parameter summary unavailable. Use fit$summary() for the full CmdStanR summary.\n")
    } else {
      cs <- as.data.frame(core_summary)

      wanted <- intersect(
        c(
          "variable", "mean", "median", "sd",
          "q5", "q95", "rhat", "ess_bulk", "ess_tail"
        ),
        names(cs)
      )

      cs <- cs[, wanted, drop = FALSE]

      numeric_cols <- vapply(cs, is.numeric, logical(1L))
      cs[numeric_cols] <- lapply(
        cs[numeric_cols],
        function(z) round(z, digits)
      )

      names(cs)[names(cs) == "q5"] <- "CrI_5%"
      names(cs)[names(cs) == "q95"] <- "CrI_95%"

      print(cs, row.names = FALSE)
      cat("Credible interval columns above use CmdStanR's default 5%-95% posterior quantiles.\n")
    }
  }

  # ============================================================
  # GUIDE
  # ============================================================

  if (isTRUE(guide)) {
    section("FIT OUTPUT GUIDE")

    cat(
      "  fit$summary()              Full posterior summary of selected/all variables\n",
      "  fit$draws()                Posterior draws in posterior package formats\n",
      "  fit$diagnostic_summary()   Divergences, treedepth and E-BFMI by chain\n",
      "  fit$cmdstan_diagnose()     Full CmdStan diagnostic report\n",
      "  fit$metadata()             CmdStan model and sampling metadata\n",
      "  fit$time()                 CmdStan timing information\n",
      "  attr(fit, 'mira_fit_info') MIRA data/sampling metadata attached to this fit\n",
      "\n",
      "  mira_summary(fit, stan_data = stan_data)\n",
      "      -> treatment, longitudinal change, gender, age, responder and clinical estimands\n",
      sep = ""
    )
  }

  line("=")

  invisible(x)
}
