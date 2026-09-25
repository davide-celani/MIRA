.mira_covariate_initial_values_long <- function(P, K, rw_scale) {
  # A 0 x (K - 1) R matrix is serialized to JSON as `[]`, which loses the
  # second dimension and cannot initialize Stan's matrix[0, K - 1]. Empty
  # Stan parameters need no explicit values, so omit the whole covariate
  # initialization block when P = 0.
  if (P == 0L) {
    return(list())
  }

  list(
    covariate_baseline_effect = rep(0, P),
    beta_covariate_time = rep(0, P),
    z_covariate_step = matrix(
      0,
      nrow = P,
      ncol = K - 1L
    ),
    tau_covariate = rep(rw_scale, P)
  )
}


#' Fit MIRA longitudinal treatment model
#'
#' Fits the outcome-adaptive MIRA longitudinal mixed-effects model with
#' treatment- and user-selected covariate-specific trajectories using CmdStan.
#' The likelihood and link are selected by [mira_data_long()]: censored
#' Student-t/Gaussian identity models and a positive log-normal/log-link model
#' share a coherent natural-scale clinical output interface.
#'
#' @param stan_data Data prepared for the MIRA Stan model. The list may
#'   contain additional R-side metadata (for example `mean_y`, `sd_y`,
#'   `arm_labels`, and `covariate_metadata`); only variables required by Stan
#'   are passed to CmdStan. Subject-level covariates are supplied through `P`
#'   and the `S` by `P` design matrix `X`. `P = 0` is supported.
#' @param prior A `mira_prior_long` object, or a named list containing the Stan
#'   prior fields required by the selected model, including prior vectors for
#'   the active encoded covariate terms. If `NULL`, `mira_prior_long()`
#'   automatically selects the outcome and uses standard priors, unless all
#'   prior fields are already present in `stan_data`. A non-NULL `prior` takes
#'   precedence over prior fields embedded in `stan_data`.
#' @param chains Number of MCMC chains.
#' @param parallel_chains Number of parallel chains.
#' @param iter_warmup Number of warmup iterations.
#' @param iter_sampling Number of sampling iterations.
#' @param seed Random seed.
#' @param refresh Number of iterations between progress messages.
#' @param adapt_delta Target average acceptance probability used during NUTS
#'   adaptation. Values closer to one reduce integration error and can remove
#'   divergent transitions, at the cost of additional computation.
#' @param step_size Positive initial integrator step size. It is adapted during
#'   warmup. A conservative default is used to avoid invalid first-iteration
#'   proposals for positive scales and Cholesky factors.
#' @param max_treedepth Maximum NUTS tree depth.
#' @param metric Euclidean metric used by NUTS: `"diag_e"`, `"dense_e"`, or
#'   `"unit_e"`. The default diagonal metric is appropriate for this
#'   non-centred hierarchical model.
#' @param rhat_threshold Upper diagnostic threshold for rank-normalized R-hat.
#'   Values below this threshold pass the automatic fit check.
#' @param ess_threshold Minimum acceptable bulk and tail effective sample size.
#' @param ebfmi_threshold Lower diagnostic threshold for E-BFMI.
#' @param treedepth_tolerance Proportion of post-warmup transitions allowed to
#'   hit `max_treedepth` before the report raises an efficiency warning. A
#'   positive proportion not exceeding this value is labelled `MINOR`: it does
#'   not by itself make the posterior fit invalid.
#' @param show_messages Logical. If `TRUE`, retain CmdStan informational
#'   messages, including rejected warmup proposals.
#' @param verbose Logical. If TRUE (default), print a compact MIRA fit report
#'   automatically after successful sampling. The fitted CmdStanMCMC object is
#'   still returned invisibly and can be assigned normally.
#' @param stan_file Optional path to the Stan file. If `NULL`, MIRA first
#'   looks for `inst/stan/mira_longitudinal.stan`.
#'
#' @return A CmdStanMCMC object.
#'
#' @export
mira_fit_long <- function(
    stan_data,
    prior = NULL,
    chains = 4,
    parallel_chains = chains,
    iter_warmup = 2000,
    iter_sampling = 3000,
    seed = 123,
    refresh = 100,
    adapt_delta = 0.95,
    step_size = 0.1,
    max_treedepth = 12,
    metric = c("diag_e", "dense_e", "unit_e"),
    show_messages = TRUE,
    stan_file = NULL,
    verbose = TRUE,
    rhat_threshold = 1.01,
    ess_threshold = 400,
    ebfmi_threshold = 0.30,
    treedepth_tolerance = 0.01
) {

  # ------------------------------------------------------------
  # Data object
  # ------------------------------------------------------------

  if (!is.list(stan_data)) {
    stop("`stan_data` must be a list.", call. = FALSE)
  }

  positive_integer_long <- function(x, name, allow_zero = FALSE) {
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

  positive_integer_long(chains, "chains")
  positive_integer_long(parallel_chains, "parallel_chains")
  positive_integer_long(iter_warmup, "iter_warmup", allow_zero = TRUE)
  positive_integer_long(iter_sampling, "iter_sampling")
  positive_integer_long(refresh, "refresh", allow_zero = TRUE)
  positive_integer_long(seed, "seed", allow_zero = TRUE)
  positive_integer_long(max_treedepth, "max_treedepth")

  if (parallel_chains > chains) {
    stop("`parallel_chains` cannot be larger than `chains`.", call. = FALSE)
  }

  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("`verbose` must be TRUE or FALSE.", call. = FALSE)
  }

  if (!is.logical(show_messages) ||
      length(show_messages) != 1L ||
      is.na(show_messages)) {
    stop("`show_messages` must be TRUE or FALSE.", call. = FALSE)
  }

  if (!is.numeric(adapt_delta) ||
      length(adapt_delta) != 1L ||
      !is.finite(adapt_delta) ||
      adapt_delta <= 0 ||
      adapt_delta >= 1) {
    stop("`adapt_delta` must be one finite number strictly between 0 and 1.",
         call. = FALSE)
  }

  if (!is.numeric(step_size) ||
      length(step_size) != 1L ||
      !is.finite(step_size) ||
      step_size <= 0) {
    stop("`step_size` must be one positive finite number.", call. = FALSE)
  }

  if (!is.numeric(rhat_threshold) ||
      length(rhat_threshold) != 1L ||
      !is.finite(rhat_threshold) ||
      rhat_threshold <= 1) {
    stop("`rhat_threshold` must be one finite number greater than 1.",
         call. = FALSE)
  }

  if (!is.numeric(ess_threshold) ||
      length(ess_threshold) != 1L ||
      !is.finite(ess_threshold) ||
      ess_threshold <= 0) {
    stop("`ess_threshold` must be one positive finite number.", call. = FALSE)
  }

  if (!is.numeric(ebfmi_threshold) ||
      length(ebfmi_threshold) != 1L ||
      !is.finite(ebfmi_threshold) ||
      ebfmi_threshold <= 0) {
    stop("`ebfmi_threshold` must be one positive finite number.",
         call. = FALSE)
  }

  if (!is.numeric(treedepth_tolerance) ||
      length(treedepth_tolerance) != 1L ||
      !is.finite(treedepth_tolerance) ||
      treedepth_tolerance < 0 ||
      treedepth_tolerance > 1) {
    stop("`treedepth_tolerance` must be one finite number between 0 and 1.",
         call. = FALSE)
  }

  metric <- match.arg(metric)

  model_data_names <- c(
    "N", "S", "K", "G", "P",
    "likelihood_id",
    "has_lower_bound", "outcome_lower_bound",
    "has_upper_bound", "outcome_upper_bound",
    "y", "subject", "time", "arm", "X", "time_value",
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

  scalar_integer_names <- c("N", "S", "K", "G", "P")

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
  if (stan_data$P < 0) stop("`P` must be >= 0.", call. = FALSE)

  if (!is.numeric(stan_data$likelihood_id) ||
      length(stan_data$likelihood_id) != 1L ||
      !is.finite(stan_data$likelihood_id) ||
      stan_data$likelihood_id != as.integer(stan_data$likelihood_id) ||
      !(stan_data$likelihood_id %in% 1:3)) {
    stop(
      "`likelihood_id` must be 1 (Student-t), 2 (Gaussian), or 3 (log-normal).",
      call. = FALSE
    )
  }

  for (nm in c("has_lower_bound", "has_upper_bound")) {
    x <- stan_data[[nm]]
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
        x != as.integer(x) || !(x %in% c(0, 1))) {
      stop("`", nm, "` must be exactly 0 or 1.", call. = FALSE)
    }
  }

  for (nm in c("outcome_lower_bound", "outcome_upper_bound")) {
    x <- stan_data[[nm]]
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x)) {
      stop("`", nm, "` must be one finite numeric value.", call. = FALSE)
    }
  }

  if (stan_data$has_lower_bound == 1L &&
      stan_data$has_upper_bound == 1L &&
      stan_data$outcome_lower_bound >= stan_data$outcome_upper_bound) {
    stop("The lower outcome bound must be smaller than the upper bound.", call. = FALSE)
  }

  if (!is.numeric(stan_data$y) ||
      length(stan_data$y) != stan_data$N ||
      any(!is.finite(stan_data$y))) {
    stop("`y` must contain exactly N finite numeric values.", call. = FALSE)
  }

  if (stan_data$has_lower_bound == 1L &&
      any(stan_data$y < stan_data$outcome_lower_bound)) {
    stop("Observed `y` falls below `outcome_lower_bound`.", call. = FALSE)
  }

  if (stan_data$has_upper_bound == 1L &&
      any(stan_data$y > stan_data$outcome_upper_bound)) {
    stop("Observed `y` exceeds `outcome_upper_bound`.", call. = FALSE)
  }

  if (stan_data$likelihood_id == 3L && any(stan_data$y <= 0)) {
    stop("The log-normal model requires strictly positive `y`.", call. = FALSE)
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
  # Subject-level covariate design matrix
  # ------------------------------------------------------------

  if (!is.matrix(stan_data$X) || !is.numeric(stan_data$X)) {
    stop("`X` must be a numeric matrix with dimensions S by P.", call. = FALSE)
  }

  expected_x_dim <- c(as.integer(stan_data$S), as.integer(stan_data$P))
  if (!identical(as.integer(dim(stan_data$X)), expected_x_dim)) {
    stop(
      "`X` must have dimensions S by P (expected ",
      expected_x_dim[[1L]], " by ", expected_x_dim[[2L]], "; got ",
      paste(dim(stan_data$X), collapse = " by "), ").",
      call. = FALSE
    )
  }

  if (any(!is.finite(stan_data$X))) {
    stop("Every element of `X` must be finite.", call. = FALSE)
  }

  covariate_names <- if (!is.null(stan_data$covariate_names)) {
    as.character(stan_data$covariate_names)
  } else if (!is.null(colnames(stan_data$X))) {
    as.character(colnames(stan_data$X))
  } else {
    paste0("covariate_", seq_len(stan_data$P))
  }

  if (length(covariate_names) != stan_data$P ||
      anyNA(covariate_names) || any(!nzchar(covariate_names)) ||
      anyDuplicated(covariate_names)) {
    stop(
      "`covariate_names` must contain exactly P unique, non-missing names.",
      call. = FALSE
    )
  }

  if (!is.null(colnames(stan_data$X)) &&
      !identical(as.character(colnames(stan_data$X)), covariate_names)) {
    stop(
      "Column names of `X` must match `covariate_names` in the same order.",
      call. = FALSE
    )
  }

  metadata_vectors <- c(
    "covariate_names", "covariate_original_names", "covariate_labels",
    "covariate_types", "covariate_reference_levels", "covariate_centers",
    "covariate_scales"
  )
  for (nm in intersect(metadata_vectors, names(stan_data))) {
    if (length(stan_data[[nm]]) != stan_data$P) {
      stop(
        "`", nm, "` must have exactly P elements.",
        call. = FALSE
      )
    }
  }

  if (!is.null(stan_data$covariate_centers) &&
      any(!is.finite(stan_data$covariate_centers))) {
    stop("Every `covariate_centers` value must be finite.", call. = FALSE)
  }
  if (!is.null(stan_data$covariate_scales) &&
      (any(!is.finite(stan_data$covariate_scales)) ||
       any(stan_data$covariate_scales <= 0))) {
    stop("Every `covariate_scales` value must be positive and finite.", call. = FALSE)
  }

  validate_covariate_table_long <- function(metadata_terms, label) {
    if (!is.data.frame(metadata_terms) || nrow(metadata_terms) != stan_data$P) {
      stop("`", label, "` must be a data frame with exactly P rows.", call. = FALSE)
    }
    metadata_name_field <- intersect(
      c("name", "covariate_name", "encoded_name", "term"),
      names(metadata_terms)
    )
    if (length(metadata_name_field) > 0L) {
      metadata_names <- as.character(metadata_terms[[metadata_name_field[[1L]]]])
      if (!identical(metadata_names, covariate_names)) {
        stop(
          "Covariate names in `", label,
          "` must match `covariate_names` in order.",
          call. = FALSE
        )
      }
    }
    invisible(TRUE)
  }

  if (!is.null(stan_data$covariate_map)) {
    validate_covariate_table_long(stan_data$covariate_map, "covariate_map")
  }

  if (!is.null(stan_data$covariate_metadata)) {
    metadata_terms <- stan_data$covariate_metadata
    if (is.list(metadata_terms) && !is.data.frame(metadata_terms)) {
      if (!is.null(metadata_terms$columns)) {
        metadata_terms <- metadata_terms$columns
      } else if (!is.null(metadata_terms$terms)) {
        metadata_terms <- metadata_terms$terms
      }
    }
    if (is.data.frame(metadata_terms)) {
      validate_covariate_table_long(metadata_terms, "covariate_metadata")
    }
  }

  if (stan_data$P > 0L) {
    constant_covariates <- vapply(
      seq_len(stan_data$P),
      function(j) length(unique(stan_data$X[, j])) < 2L,
      logical(1L)
    )
    if (any(constant_covariates)) {
      warning(
        "The following active design-matrix columns are constant and their ",
        "covariate trajectories will be weakly/non-identified: ",
        paste(covariate_names[constant_covariates], collapse = ", "),
        ".",
        call. = FALSE
      )
    }
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

  positive_scalar_long <- function(x, name, allow_zero = FALSE) {
    ok <- length(x) == 1 && is.numeric(x) && is.finite(x)
    ok <- ok && if (allow_zero) x >= 0 else x > 0
    if (!ok) {
      comparator <- if (allow_zero) ">= 0" else "> 0"
      stop("`", name, "` must be one finite numeric value ", comparator, ".", call. = FALSE)
    }
  }

  positive_scalar_long(stan_data$mcid_prior_mean, "mcid_prior_mean")
  positive_scalar_long(stan_data$mcid_prior_sd, "mcid_prior_sd")
  positive_scalar_long(
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
    "covariate_baseline_prior_sd",
    "beta_covariate_prior_sd",
    "tau_covariate_prior_rate",
    "sigma_intercept_prior_rate",
    "sigma_slope_prior_rate",
    "sigma_prior_rate",
    "nu_prior_shape",
    "nu_prior_rate"
  )

  covariate_prior_names <- c(
    "covariate_baseline_prior_sd",
    "beta_covariate_prior_sd",
    "tau_covariate_prior_rate"
  )

  validate_prior_covariate_identity_long <- function(
      prior_values,
      metadata_names = NULL,
      source = "prior"
  ) {
    named_prior_fields <- lapply(
      covariate_prior_names,
      function(name) names(prior_values[[name]])
    )
    named_prior_fields <- Filter(Negate(is.null), named_prior_fields)
    if (length(named_prior_fields) > 1L &&
        !all(vapply(
          named_prior_fields[-1L],
          identical,
          logical(1L),
          named_prior_fields[[1L]]
        ))) {
      stop(
        "Names on the generic covariate prior vectors must match in order.",
        call. = FALSE
      )
    }

    prior_covariate_names <- if (!is.null(metadata_names)) {
      as.character(metadata_names)
    } else if (length(named_prior_fields) > 0L) {
      named_prior_fields[[1L]]
    } else {
      NULL
    }
    if (!is.null(prior_covariate_names) &&
        !identical(prior_covariate_names, covariate_names)) {
      stop(
        "The ", source, " covariate names/order do not match ",
        "`stan_data$covariate_names`. Rebuild the prior from these data ",
        "or provide values in the active encoded-term order.",
        call. = FALSE
      )
    }
    invisible(TRUE)
  }

  embedded_prior_names <- intersect(prior_names, names(stan_data))

  # An explicit `prior` argument takes precedence over any embedded fields.
  if (!is.null(prior)) {

    mira_validate_prior_long(prior)

    prior_P <- length(prior$covariate_baseline_prior_sd)
    if (prior_P != stan_data$P) {
      stop(
        "The supplied prior was built for P = ", prior_P,
        ", but `stan_data` has P = ", stan_data$P, ".",
        call. = FALSE
      )
    }

    # `mira_prior_stan_data_long()` intentionally removes R-side names before
    # serialization. Check identity/order first so that a prior built for a
    # different same-size design matrix cannot silently target the wrong term.
    validate_prior_covariate_identity_long(
      prior,
      metadata_names = prior$covariate_names,
      source = "supplied prior"
    )

    if (!is.null(prior$likelihood)) {
      expected_likelihood <- c("student_t", "gaussian", "lognormal")[[
        as.integer(stan_data$likelihood_id)
      ]]
      supplied_likelihood <- tolower(as.character(prior$likelihood)[1L])
      if (!identical(supplied_likelihood, expected_likelihood)) {
        stop(
          "The prior was built for `", supplied_likelihood,
          "`, but `stan_data` requests `", expected_likelihood, "`.",
          call. = FALSE
        )
      }
    }

    stan_prior_data <- mira_prior_stan_data_long(prior)
    prior_source <- "argument"

  } else if (all(prior_names %in% names(stan_data))) {

    stan_prior_data <- stan_data[prior_names]
    validate_prior_covariate_identity_long(
      stan_prior_data,
      source = "embedded prior"
    )
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

    prior <- mira_prior_long(
      stan_data = stan_data,
      outcome = "auto",
      informativeness = "standard"
    )
    stan_prior_data <- mira_prior_stan_data_long(prior)
    prior_source <- "automatic default"
  }

  missing_prior <- setdiff(prior_names, names(stan_prior_data))

  if (length(missing_prior) > 0) {
    stop(
      "The prior specification is not compatible with the MIRA Stan model. ",
      "Missing Stan prior fields: ",
      paste(missing_prior, collapse = ", "),
      ". Update `mira_prior_long()` / `mira_prior_stan_data_long()` so that generic ",
      "covariate priors are included, or pass a complete named prior list.",
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

  positive_prior_names <- setdiff(
    prior_names,
    c(finite_prior_names, covariate_prior_names)
  )
  for (nm in positive_prior_names) {
    positive_scalar_long(stan_prior_data[[nm]], nm)
  }

  for (nm in covariate_prior_names) {
    x <- stan_prior_data[[nm]]
    if (!is.numeric(x) || length(x) != stan_data$P ||
        any(!is.finite(x)) || any(x <= 0)) {
      stop(
        "`", nm, "` must be a positive finite numeric vector of length P ",
        "(use `numeric(0)` when P = 0).",
        call. = FALSE
      )
    }
    stan_prior_data[[nm]] <- as.numeric(x)
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
  sampling_data$P <- as.integer(sampling_data$P)
  sampling_data$likelihood_id <- as.integer(sampling_data$likelihood_id)
  sampling_data$has_lower_bound <- as.integer(sampling_data$has_lower_bound)
  sampling_data$has_upper_bound <- as.integer(sampling_data$has_upper_bound)
  sampling_data$subject <- as.integer(sampling_data$subject)
  sampling_data$time <- as.integer(sampling_data$time)
  sampling_data$arm <- as.integer(sampling_data$arm)
  storage.mode(sampling_data$X) <- "double"
  sampling_data$direction <- as.integer(sampling_data$direction)

  # ------------------------------------------------------------
  # Stan model file
  # ------------------------------------------------------------

  if (is.null(stan_file)) {

    candidates <- "mira_longitudinal.stan"

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
        "`mira_longitudinal.stan`.",
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
  model_y <- if (stan_data$likelihood_id == 3L) log(y) else y
  baseline_y <- model_y[stan_data$time == 1]

  mean_y <- if (!is.null(stan_data$mean_model_y)) {
    as.numeric(stan_data$mean_model_y)[1]
  } else {
    mean(model_y)
  }

  sd_y <- if (!is.null(stan_data$sd_model_y)) {
    as.numeric(stan_data$sd_model_y)[1]
  } else {
    stats::sd(model_y)
  }

  if (!is.finite(mean_y)) mean_y <- mean(model_y)
  if (!is.finite(sd_y) || sd_y <= 0) {
    sd_y <- if (stan_data$likelihood_id == 3L) 0.2 else max(abs(mean_y) * 0.1, 1)
  }

  baseline_init <- if (length(baseline_y) > 0) mean(baseline_y) else mean_y
  elapsed <- max(stan_data$time_value) - min(stan_data$time_value)
  elapsed_safe <- max(elapsed, 1e-6)

  scale_floor <- if (stan_data$likelihood_id == 3L) 0.001 else 0.01
  slope_scale <- max(sd_y / elapsed_safe, scale_floor)
  rw_scale <- max(sd_y / sqrt(elapsed_safe) / 10, scale_floor)

  init_long <- function() {
    initial_values <- c(
      list(
        baseline_mean = baseline_init,
        beta_time = as.numeric(stan_prior_data$beta_time_prior_mean),
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
        arm_baseline_sd = max(sd_y / 10, scale_floor)
      ),
      .mira_covariate_initial_values_long(
        P = stan_data$P,
        K = stan_data$K,
        rw_scale = rw_scale
      ),
      list(
        z_subject = matrix(
          0,
          nrow = 2,
          ncol = stan_data$S
        ),
        sigma_subject = c(
          max(sd_y / 2, scale_floor),
          max(slope_scale / 2, scale_floor)
        ),
        L_subject = diag(2),
        sigma = max(sd_y / 2, scale_floor),
        mcid = max(stan_data$mcid_prior_mean, 1e-6)
      )
    )

    # Name nu in every init list so CmdStanR does not report the inactive
    # vector[0] as missing for Gaussian and log-normal models.
    initial_values$nu <- if (stan_data$likelihood_id == 1L) {
      array(10, dim = 1L)
    } else {
      numeric(0)
    }
    initial_values
  }

  # ------------------------------------------------------------
  # Sampling
  # ------------------------------------------------------------

  message("Sampling posterior...")

  sampling_started <- Sys.time()

  sample_model_long <- function() {
    model$sample(
      data = sampling_data,
      chains = chains,
      parallel_chains = parallel_chains,
      iter_warmup = iter_warmup,
      iter_sampling = iter_sampling,
      seed = seed,
      init = init_long,
      refresh = refresh,
      adapt_delta = adapt_delta,
      step_size = step_size,
      max_treedepth = max_treedepth,
      metric = metric,
      show_messages = show_messages
    )
  }

  fit <- if (stan_data$P == 0L) {
    # CmdStanR reports intentionally omitted zero-dimensional parameters as
    # partial initial values. Silence only that notice for this sample call
    # and restore the user's option even if sampling fails.
    local({
      old_options <- options(cmdstanr_warn_inits = FALSE)
      on.exit(options(old_options), add = TRUE)
      sample_model_long()
    })
  } else {
    sample_model_long()
  }

  sampling_elapsed_seconds <- as.numeric(
    difftime(Sys.time(), sampling_started, units = "secs")
  )

  # ------------------------------------------------------------
  # Attach lightweight MIRA metadata to the CmdStanMCMC object.
  #
  # The object remains a CmdStanMCMC object, so existing code such as
  # mira_summary_long(fit, ...) and all CmdStanR $methods continue to work.
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

  covariate_field_long <- function(name, default) {
    value <- stan_data[[name]]
    if (is.null(value) || length(value) != stan_data$P) default else value
  }

  covariate_map <- stan_data$covariate_map
  if (!is.data.frame(covariate_map) || nrow(covariate_map) != stan_data$P) {
    covariate_map <- data.frame(
      index = seq_len(stan_data$P),
      name = covariate_names,
      label = as.character(covariate_field_long("covariate_labels", covariate_names)),
      original_name = as.character(covariate_field_long(
        "covariate_original_names", covariate_names
      )),
      type = as.character(covariate_field_long(
        "covariate_types", rep("unspecified", stan_data$P)
      )),
      reference_level = as.character(covariate_field_long(
        "covariate_reference_levels", rep(NA_character_, stan_data$P)
      )),
      center = as.numeric(covariate_field_long(
        "covariate_centers", rep(0, stan_data$P)
      )),
      scale = as.numeric(covariate_field_long(
        "covariate_scales", rep(1, stan_data$P)
      )),
      stringsAsFactors = FALSE
    )
  }

  if (!"index" %in% names(covariate_map)) {
    covariate_map$index <- seq_len(stan_data$P)
  }
  if (!"name" %in% names(covariate_map)) {
    covariate_map$name <- covariate_names
  }

  covariate_distributions <- data.frame(
    index = seq_len(stan_data$P),
    name = covariate_names,
    n = rep.int(stan_data$S, stan_data$P),
    n_unique = integer(stan_data$P),
    mean = numeric(stan_data$P),
    sd = numeric(stan_data$P),
    min = numeric(stan_data$P),
    q25 = numeric(stan_data$P),
    median = numeric(stan_data$P),
    q75 = numeric(stan_data$P),
    max = numeric(stan_data$P),
    n_zero = integer(stan_data$P),
    n_one = integer(stan_data$P),
    stringsAsFactors = FALSE
  )

  if (stan_data$P > 0L) {
    for (j in seq_len(stan_data$P)) {
      xj <- as.numeric(stan_data$X[, j])
      qj <- stats::quantile(xj, c(0.25, 0.50, 0.75), names = FALSE)
      covariate_distributions$n_unique[j] <- length(unique(xj))
      covariate_distributions$mean[j] <- mean(xj)
      covariate_distributions$sd[j] <- stats::sd(xj)
      covariate_distributions$min[j] <- min(xj)
      covariate_distributions$q25[j] <- qj[[1L]]
      covariate_distributions$median[j] <- qj[[2L]]
      covariate_distributions$q75[j] <- qj[[3L]]
      covariate_distributions$max[j] <- max(xj)
      covariate_distributions$n_zero[j] <- sum(xj == 0)
      covariate_distributions$n_one[j] <- sum(xj == 1)
    }
  }

  prior_profile <- if (!is.null(prior) && !is.null(prior$profile)) {
    as.character(prior$profile)[1]
  } else if (identical(prior_source, "embedded in stan_data")) {
    "embedded in stan_data"
  } else {
    "named/custom prior list"
  }

  population_reference_profile <- "X = 0 on the encoded design-matrix scale"
  if (is.list(stan_data$covariate_metadata) &&
      is.list(stan_data$covariate_metadata$reference_profile)) {
    reference_description <-
      stan_data$covariate_metadata$reference_profile$description
    if (is.character(reference_description) &&
        length(reference_description) == 1L &&
        !is.na(reference_description) && nzchar(reference_description)) {
      population_reference_profile <- reference_description
    }
  }

  mira_fit_info <- list(
    model_file = basename(stan_file),
    model_path = tryCatch(
      normalizePath(stan_file, winslash = "/", mustWork = FALSE),
      error = function(e) stan_file
    ),
    outcome = if (!is.null(stan_data$outcome)) {
      as.character(stan_data$outcome)[1L]
    } else if (!is.null(stan_data$outcome_name)) {
      as.character(stan_data$outcome_name)[1L]
    } else {
      "generic"
    },
    outcome_name = if (!is.null(stan_data$outcome_name)) {
      as.character(stan_data$outcome_name)[1L]
    } else {
      NA_character_
    },
    likelihood = c("student_t", "gaussian", "lognormal")[[
      as.integer(stan_data$likelihood_id)
    ]],
    likelihood_id = as.integer(stan_data$likelihood_id),
    modeling_scale = if (stan_data$likelihood_id == 3L) "log" else "identity",
    clinical_estimand_scale = "absolute natural-outcome units",
    outcome_bounds = c(
      lower = if (stan_data$has_lower_bound == 1L) {
        stan_data$outcome_lower_bound
      } else {
        NA_real_
      },
      upper = if (stan_data$has_upper_bound == 1L) {
        stan_data$outcome_upper_bound
      } else {
        NA_real_
      }
    ),
    boundary_strategy = if (stan_data$likelihood_id == 3L &&
                            stan_data$has_upper_bound == 0L &&
                            (stan_data$has_lower_bound == 0L ||
                             stan_data$outcome_lower_bound <= 0)) {
      "strictly positive log-normal support"
    } else if (stan_data$likelihood_id == 3L) {
      "strictly positive log-normal support with endpoint censoring"
    } else if (stan_data$has_lower_bound == 1L ||
               stan_data$has_upper_bound == 1L) {
      "endpoint censoring"
    } else {
      "none"
    },
    outcome_diagnostics = stan_data$outcome_diagnostics,
    n_observations = as.integer(stan_data$N),
    n_subjects = as.integer(stan_data$S),
    n_time_points = as.integer(stan_data$K),
    n_arms = as.integer(stan_data$G),
    n_covariates = as.integer(stan_data$P),
    P = as.integer(stan_data$P),
    arm_labels = arm_labels,
    arm_counts = arm_counts,
    covariate_names = covariate_names,
    covariate_original_names = as.character(covariate_field_long(
      "covariate_original_names", covariate_names
    )),
    covariate_labels = as.character(covariate_field_long(
      "covariate_labels", covariate_names
    )),
    covariate_types = as.character(covariate_field_long(
      "covariate_types", rep("unspecified", stan_data$P)
    )),
    covariate_reference_levels = covariate_field_long(
      "covariate_reference_levels", rep(NA_character_, stan_data$P)
    ),
    covariate_centers = as.numeric(covariate_field_long(
      "covariate_centers", rep(0, stan_data$P)
    )),
    covariate_scales = as.numeric(covariate_field_long(
      "covariate_scales", rep(1, stan_data$P)
    )),
    covariate_map = covariate_map,
    covariate_metadata = stan_data$covariate_metadata,
    covariates_requested = stan_data$covariates_requested,
    covariates_selected = stan_data$covariates_selected,
    covariate_distributions = covariate_distributions,
    population_reference_profile = population_reference_profile,
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
    adapt_delta = as.numeric(adapt_delta),
    initial_step_size = as.numeric(step_size),
    max_treedepth = as.integer(max_treedepth),
    metric = metric,
    rhat_threshold = as.numeric(rhat_threshold),
    ess_threshold = as.numeric(ess_threshold),
    ebfmi_threshold = as.numeric(ebfmi_threshold),
    treedepth_tolerance = as.numeric(treedepth_tolerance),
    show_messages = show_messages,
    compile_elapsed_seconds = compile_elapsed_seconds,
    sampling_elapsed_seconds = sampling_elapsed_seconds
  )

  attr(fit, "mira_fit_info") <- mira_fit_info
  class(fit) <- unique(c("mira_fit_long", class(fit)))

  if (isTRUE(verbose)) {
    print(fit)
  }

  invisible(fit)
}


#' Print a fitted MIRA model
#'
#' Compact technical report for a fitted MIRA CmdStanMCMC object. The report
#' focuses on model setup, sampling quality and core posterior parameters.
#' Detailed clinical estimands are intentionally left to [mira_summary_long()].
#'
#' @param x A fitted object returned by [mira_fit_long()].
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
print.mira_fit_long <- function(
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

  line_long <- function(char = "-", n = 92L) {
    cat(strrep(char, n), "\n", sep = "")
  }

  section_long <- function(title) {
    cat("\n", title, "\n", sep = "")
    line_long()
  }

  fmt_num_long <- function(z, d = digits) {
    if (length(z) == 0L || is.na(z) || !is.finite(z)) {
      return("NA")
    }
    formatC(z, format = "f", digits = d)
  }

  fmt_time_long <- function(seconds) {
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

  fmt_flag_long <- function(ok) {
    if (length(ok) == 0L || is.na(ok)) return("UNKNOWN")
    if (isTRUE(ok)) "OK" else "CHECK"
  }

  diagnostic_limit_long <- function(value, fallback) {
    if (length(value) == 1L && is.numeric(value) && is.finite(value)) {
      return(as.numeric(value))
    }
    fallback
  }

  safe_metadata <- tryCatch(x$metadata(), error = function(e) NULL)

  # ============================================================
  # HEADER
  # ============================================================

  cat("\n")
  line_long("=")
  cat("MIRA BAYESIAN LONGITUDINAL MODEL - FIT REPORT\n")
  line_long("=")

  # ============================================================
  # MODEL / DATA
  # ============================================================

  section_long("MODEL AND DATA")

  if (!is.null(info$model_file)) {
    cat(sprintf("Stan model: %s\n", info$model_file))
  }

  if (!is.null(info$outcome) || !is.null(info$likelihood)) {
    cat(
      "Outcome: ", if (!is.null(info$outcome)) info$outcome else "unknown",
      " | Likelihood: ",
      if (!is.null(info$likelihood)) info$likelihood else "unknown",
      " | Modeling scale: ",
      if (!is.null(info$modeling_scale)) info$modeling_scale else "unknown",
      "\n",
      sep = ""
    )
  }

  if (!is.null(info$outcome_bounds) && any(is.finite(info$outcome_bounds))) {
    bound_text <- paste0(
      if (is.finite(info$outcome_bounds[[1L]])) info$outcome_bounds[[1L]] else "-Inf",
      ", ",
      if (is.finite(info$outcome_bounds[[2L]])) info$outcome_bounds[[2L]] else "Inf"
    )
    cat("Observable support: [", bound_text, "] via ",
        info$boundary_strategy, "\n", sep = "")
  }

  if (!is.null(info$n_subjects)) {
    cat(sprintf(
      "Subjects: %d | Observations: %d | Timepoints: %d | Arms: %d | Covariate terms: %d\n",
      info$n_subjects,
      info$n_observations,
      info$n_time_points,
      info$n_arms,
      if (!is.null(info$n_covariates)) info$n_covariates else 0L
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
      fmt_num_long(info$mcid_prior_mean),
      fmt_num_long(info$mcid_prior_sd),
      fmt_num_long(info$meaningful_between_arm_difference)
    ))
  }

  if (!is.null(info$arm_counts)) {
    cat("\nSubjects by treatment arm:\n")
    print(info$arm_counts)
  }

  if (!is.null(info$n_covariates) && info$n_covariates == 0L) {
    cat("Active covariates: none (time-only adjustment profile; P = 0).\n")
  } else if (is.data.frame(info$covariate_map) && nrow(info$covariate_map) > 0L) {
    cat("\nActive encoded covariate terms:\n")
    keep <- intersect(
      c(
        "index", "name", "label", "original_name", "type", "encoding",
        "level", "reference_level", "center", "scale", "unit",
        "n_reference", "n_comparison"
      ),
      names(info$covariate_map)
    )
    print(info$covariate_map[, keep, drop = FALSE], row.names = FALSE)
  }

  # ============================================================
  # SAMPLING
  # ============================================================

  section_long("SAMPLING")

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
    if (!is.null(info$adapt_delta)) {
      cat(sprintf(
        paste0(
          "NUTS controls: adapt_delta=%s | initial step_size=%s | ",
          "max_treedepth=%d | metric=%s\n"
        ),
        fmt_num_long(info$adapt_delta, 2L),
        fmt_num_long(info$initial_step_size, 3L),
        info$max_treedepth,
        info$metric
      ))
    }
    cat(sprintf(
      "Compilation: %s | Sampling: %s\n",
      fmt_time_long(info$compile_elapsed_seconds),
      fmt_time_long(info$sampling_elapsed_seconds)
    ))
  } else if (!is.null(safe_metadata)) {
    cat("CmdStanR metadata are available via fit$metadata().\n")
  }

  # ============================================================
  # DIAGNOSTICS
  # ============================================================

  if (isTRUE(diagnostics)) {
    section_long("MCMC DIAGNOSTICS")

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

    rhat_limit <- diagnostic_limit_long(info$rhat_threshold, 1.01)
    ess_limit <- diagnostic_limit_long(info$ess_threshold, 400)
    ebfmi_limit <- diagnostic_limit_long(info$ebfmi_threshold, 0.30)
    treedepth_limit <- diagnostic_limit_long(info$treedepth_tolerance, 0.01)

    post_warmup_draws <- diagnostic_limit_long(info$post_warmup_draws, NA_real_)
    treedepth_rate <- if (is.finite(treedepth_hits) &&
                          is.finite(post_warmup_draws) &&
                          post_warmup_draws > 0) {
      treedepth_hits / post_warmup_draws
    } else {
      NA_real_
    }

    treedepth_status <- if (!is.finite(treedepth_hits)) {
      "UNKNOWN"
    } else if (treedepth_hits == 0) {
      "OK"
    } else if (is.finite(treedepth_rate) &&
               treedepth_rate <= treedepth_limit) {
      "MINOR"
    } else {
      "CHECK"
    }

    treedepth_text <- if (is.finite(treedepth_hits)) {
      if (is.finite(treedepth_rate)) {
        sprintf("%d (%.2f%%)", as.integer(treedepth_hits), 100 * treedepth_rate)
      } else {
        as.character(treedepth_hits)
      }
    } else {
      "NA"
    }

    cat(sprintf(
      "Divergences: %s [%s] | Max-treedepth hits: %s [%s]\n",
      ifelse(is.finite(divergences), as.character(divergences), "NA"),
      fmt_flag_long(if (is.finite(divergences)) divergences == 0 else NA),
      treedepth_text,
      treedepth_status
    ))

    cat(sprintf(
      "Max R-hat: %s [%s] | Min bulk ESS: %s [%s] | Min tail ESS: %s [%s]\n",
      fmt_num_long(max_rhat),
      fmt_flag_long(if (is.finite(max_rhat)) max_rhat < rhat_limit else NA),
      fmt_num_long(min_ess_bulk, 0L),
      fmt_flag_long(if (is.finite(min_ess_bulk)) min_ess_bulk >= ess_limit else NA),
      fmt_num_long(min_ess_tail, 0L),
      fmt_flag_long(if (is.finite(min_ess_tail)) min_ess_tail >= ess_limit else NA)
    ))

    cat(sprintf(
      "Minimum E-BFMI: %s [%s]\n",
      fmt_num_long(min_ebfmi),
      fmt_flag_long(if (is.finite(min_ebfmi)) min_ebfmi > ebfmi_limit else NA)
    ))

    cat(sprintf(
      paste0(
        "Thresholds: R-hat < %.3f | ESS >= %.0f | E-BFMI > %.2f | ",
        "minor treedepth <= %.2f%%\n"
      ),
      rhat_limit,
      ess_limit,
      ebfmi_limit,
      100 * treedepth_limit
    ))

    critical_checks <- c(
      if (is.finite(divergences)) divergences == 0 else NA,
      if (is.finite(max_rhat)) max_rhat < rhat_limit else NA,
      if (is.finite(min_ess_bulk)) min_ess_bulk >= ess_limit else NA,
      if (is.finite(min_ess_tail)) min_ess_tail >= ess_limit else NA,
      if (is.finite(min_ebfmi)) min_ebfmi > ebfmi_limit else NA
    )

    known_critical_checks <- critical_checks[!is.na(critical_checks)]

    overall <- if (length(known_critical_checks) == 0L) {
      "DIAGNOSTICS UNAVAILABLE"
    } else if (!all(known_critical_checks)) {
      "CHECK MCMC DIAGNOSTICS BEFORE INTERPRETATION"
    } else if (identical(treedepth_status, "CHECK")) {
      "CRITICAL DIAGNOSTICS OK; CHECK TREEDEPTH EFFICIENCY"
    } else if (identical(treedepth_status, "MINOR")) {
      "NO CRITICAL MCMC PROBLEMS; MINOR TREEDEPTH WARNING"
    } else {
      "NO OBVIOUS MCMC PROBLEMS DETECTED"
    }

    cat("Overall: ", overall, "\n", sep = "")

    if (is.finite(divergences) && divergences > 0) {
      cat(
        "Action: divergences affect posterior reliability; consider a higher ",
        "adapt_delta and inspect model geometry.\n",
        sep = ""
      )
    }

    if (identical(treedepth_status, "MINOR")) {
      cat(
        "Note: the small treedepth proportion is an efficiency warning only; ",
        "the critical diagnostics determine interpretability.\n",
        sep = ""
      )
    } else if (identical(treedepth_status, "CHECK")) {
      cat(
        "Action: inspect posterior geometry before only increasing ",
        "max_treedepth; frequent hits can make sampling very slow.\n",
        sep = ""
      )
    }
  }

  # ============================================================
  # CORE PARAMETERS
  # ============================================================

  if (isTRUE(parameters)) {
    section_long("CORE POSTERIOR PARAMETERS")

    core_candidates <- c(
      "baseline_mean",
      "beta_time",
      "tau_common",
      "beta_treatment",
      "tau_treatment",
      "arm_baseline_sd",
      "covariate_baseline_effect",
      "beta_covariate_time",
      "tau_covariate",
      "sigma_subject",
      "sigma",
      "nu",
      "nu_value",
      "residual_sd",
      "residual_cv",
      "mcid"
    )

    if (!is.null(info$likelihood) && info$likelihood != "student_t") {
      core_candidates <- setdiff(core_candidates, c("nu", "nu_value"))
    }
    if (!is.null(info$likelihood) && info$likelihood != "lognormal") {
      core_candidates <- setdiff(core_candidates, "residual_cv")
    }
    if (!is.null(info$n_covariates) && info$n_covariates == 0L) {
      core_candidates <- setdiff(
        core_candidates,
        c(
          "covariate_baseline_effect",
          "beta_covariate_time",
          "tau_covariate"
        )
      )
    }

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
    section_long("FIT OUTPUT GUIDE")

    cat(
      "  fit$summary()              Full posterior summary of selected/all variables\n",
      "  fit$draws()                Posterior draws in posterior package formats\n",
      "  fit$diagnostic_summary()   Divergences, treedepth and E-BFMI by chain\n",
      "  fit$cmdstan_diagnose()     Full CmdStan diagnostic report\n",
      "  fit$metadata()             CmdStan model and sampling metadata\n",
      "  fit$time()                 CmdStan timing information\n",
      "  attr(fit, 'mira_fit_info') MIRA data/sampling metadata attached to this fit\n",
      "\n",
      "  mira_summary_long(fit, stan_data = stan_data)\n",
      "      -> treatment, longitudinal change, selected-covariate, responder and clinical estimands\n",
      sep = ""
    )
  }

  line_long("=")

  invisible(x)
}
