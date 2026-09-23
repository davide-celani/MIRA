dynamic_covariate_fixture <- function() {
  out <- data.frame(
    patient = sprintf("P%02d", seq_len(8L)),
    arm = rep(c("Control", "Treatment"), 4L),
    age = c(45, 50, 55, 60, 65, 70, 75, 80),
    gender = factor(
      rep(c("Female", "Male"), 4L),
      levels = c("Female", "Male")
    ),
    smoking = factor(
      c("Never", "Former", "Current", "Never",
        "Former", "Current", "Never", "Former"),
      levels = c("Never", "Former", "Current")
    ),
    baseline_score = c(10, 13, 11, 16, 15, 19, 18, 22),
    integer_score = as.integer(c(2, 4, 3, 5, 7, 6, 9, 8)),
    eligible_flag = c(FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, TRUE),
    region = c("North", "South", "West", "North",
               "South", "West", "North", "South"),
    unselected_marker = c(101, 103, 102, 108, 105, 107, 104, 106),
    BCVA_t0 = c(60, 61, 62, 63, 64, 65, 66, 67),
    BCVA_t1 = c(61, 62, 64, 64, 66, 66, 68, 69),
    BCVA_t2 = c(62, 63, 65, 66, 67, 68, 69, 71),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  # Columns with internal/identifier semantics must never enter auto mode.
  out[["_row_tag"]] <- seq_len(nrow(out))
  out$subject_id <- paste0("external-", seq_len(nrow(out)))
  out
}


prepare_dynamic_covariates <- function(covariates = "auto", ...) {
  mira_prepare_data(
    data = dynamic_covariate_fixture(),
    time_value = c(0, 1, 3),
    outcome = "BCVA",
    likelihood = "gaussian",
    meaningful_change = 5,
    meaningful_change_sd = 1,
    direction = "higher",
    meaningful_between_arm_difference = 5,
    reference_arm = "Control",
    covariates = covariates,
    ...
  )
}


dynamic_covariate_stan_file <- function() {
  installed <- system.file(
    "stan", "gaussian_longitudinal.stan", package = "MIRA"
  )
  if (nzchar(installed)) return(installed)

  testthat::test_path("..", "..", "inst", "stan", "gaussian_longitudinal.stan")
}


mock_dynamic_covariate_fit <- function(stan_data, n_draws = 20L) {
  variable_names <- "mcid"
  if (stan_data$P > 0L) {
    indexed_names <- unlist(lapply(seq_len(stan_data$P), function(p) {
      unlist(lapply(seq_len(stan_data$K), function(k) {
        paste0(
          c(
            "covariate_effect",
            "covariate_population_mean_difference",
            "covariate_population_mean_change_difference",
            "directional_covariate_population_mean_change_difference"
          ),
          "[", p, ",", k, "]"
        )
      }), use.names = FALSE)
    }), use.names = FALSE)
    variable_names <- c(variable_names, indexed_names)
  }

  draws <- vapply(seq_along(variable_names), function(j) {
    seq(-0.2, 0.2, length.out = n_draws) + j / 100
  }, numeric(n_draws))
  colnames(draws) <- variable_names
  draws[, "mcid"] <- seq(4.5, 5.5, length.out = n_draws)

  fit <- list(
    draws = function(format = "draws_matrix", ...) draws,
    summary = function(...) {
      data.frame(
        variable = colnames(draws),
        rhat = rep(1, ncol(draws)),
        ess_bulk = rep(1000, ncol(draws)),
        ess_tail = rep(1000, ncol(draws)),
        stringsAsFactors = FALSE
      )
    }
  )
  class(fit) <- c("CmdStanMCMC", "list")
  fit
}


dynamic_covariate_sampling_data <- function(stan_data, prior) {
  model_fields <- c(
    "N", "S", "K", "G", "P", "likelihood_id",
    "has_lower_bound", "outcome_lower_bound",
    "has_upper_bound", "outcome_upper_bound",
    "y", "subject", "time", "arm", "X", "time_value", "direction",
    "mcid_prior_mean", "mcid_prior_sd",
    "meaningful_between_arm_difference"
  )
  c(stan_data[model_fields], mira_prior_stan_data(prior))
}


test_that("A: NULL and character(0) produce a complete P = 0 pipeline", {
  none <- prepare_dynamic_covariates(character(0))
  null <- prepare_dynamic_covariates(NULL)

  for (stan_data in list(none, null)) {
    expect_identical(stan_data$P, 0L)
    expect_true(is.matrix(stan_data$X))
    expect_identical(dim(stan_data$X), c(stan_data$S, 0L))
    expect_length(stan_data$covariate_names, 0L)
    expect_length(stan_data$covariates_selected, 0L)
    expect_equal(nrow(stan_data$covariate_map), 0L)
    expect_length(
      stan_data$covariate_metadata$reference_profile$design_values,
      0L
    )

    prior <- mira_prior(stan_data)
    prior_data <- mira_prior_stan_data(prior)
    expect_length(prior_data$covariate_baseline_prior_sd, 0L)
    expect_length(prior_data$beta_covariate_prior_sd, 0L)
    expect_length(prior_data$tau_covariate_prior_rate, 0L)
  }

  expect_false(none$covariate_metadata$selection$request_was_null)
  expect_true(null$covariate_metadata$selection$request_was_null)
})


test_that("A: the generic summary has a typed empty P = 0 result", {
  stan_data <- prepare_dynamic_covariates(character(0))
  fit <- mock_dynamic_covariate_fit(stan_data)
  result <- mira_summary(fit, stan_data = stan_data, verbose = FALSE)

  expect_s3_class(result, "mira_summary")
  expect_true(is.data.frame(result$covariate_effects))
  expect_equal(nrow(result$covariate_effects), 0L)
  expect_identical(result$covariates$P, 0L)
  expect_length(result$covariates$names, 0L)
  expect_equal(nrow(result$covariates$map), 0L)
  expect_null(result$gender_effects)
  expect_null(result$age_threshold_effects)
})


test_that("B and I: an explicit one-variable subset controls P and X", {
  stan_data <- prepare_dynamic_covariates("age")

  expect_identical(stan_data$P, 1L)
  expect_identical(colnames(stan_data$X), "age")
  expect_identical(stan_data$covariates_selected, "age")
  expect_identical(stan_data$covariate_original_names, "age")
  expect_false(any(c("gender", "smoking", "baseline_score") %in%
                     stan_data$covariate_original_names))
  expect_false(any(grepl(
    "gender|smoking|baseline_score|unselected_marker",
    colnames(stan_data$X),
    ignore.case = TRUE
  )))
})


test_that("C: several selected variables determine encoded P", {
  stan_data <- prepare_dynamic_covariates(c("age", "gender", "smoking"))

  expect_identical(stan_data$P, 4L)
  expect_identical(
    colnames(stan_data$X),
    c("age", "genderMale", "smokingFormer", "smokingCurrent")
  )
  expect_identical(
    stan_data$covariates_selected,
    c("age", "gender", "smoking")
  )
  expect_identical(dim(stan_data$X), c(stan_data$S, stan_data$P))
  expect_true(all(is.finite(stan_data$X)))
})


test_that("D: auto selects all and only eligible subject-level variables", {
  stan_data <- prepare_dynamic_covariates("auto")
  expected <- c(
    "age", "gender", "smoking", "baseline_score", "integer_score",
    "eligible_flag", "region", "unselected_marker"
  )

  expect_identical(stan_data$covariates_selected, expected)
  expect_setequal(unique(stan_data$covariate_original_names), expected)
  expect_false(any(c(
    "patient", "arm", "BCVA_t0", "BCVA_t1", "BCVA_t2",
    "_row_tag", "subject_id"
  ) %in% stan_data$covariate_original_names))
  expect_identical(stan_data$P, as.integer(ncol(stan_data$X)))
  expect_gt(stan_data$P, length(expected))
})


test_that("E: multilevel factors use treatment coding and retain metadata", {
  stan_data <- prepare_dynamic_covariates("smoking")
  map <- stan_data$covariate_map

  expect_identical(stan_data$P, 2L)
  expect_identical(colnames(stan_data$X), c("smokingFormer", "smokingCurrent"))
  expect_identical(map$original_name, rep("smoking", 2L))
  expect_identical(map$reference_level, rep("Never", 2L))
  expect_identical(map$level, c("Former", "Current"))
  expect_identical(map$encoding, rep("treatment", 2L))
  expect_true(all(stan_data$X %in% c(0, 1)))

  never <- dynamic_covariate_fixture()$smoking == "Never"
  expect_true(all(rowSums(stan_data$X[never, , drop = FALSE]) == 0))
})


test_that("E and I: readable selected names survive the generic summary", {
  stan_data <- prepare_dynamic_covariates("smoking")
  fit <- mock_dynamic_covariate_fit(stan_data)
  result <- mira_summary(fit, stan_data = stan_data, verbose = FALSE)

  expect_s3_class(result, "mira_summary")
  expect_true(is.data.frame(result$covariate_effects))
  expect_setequal(
    unique(result$covariate_effects$encoded_term),
    c("smokingFormer", "smokingCurrent")
  )
  expect_identical(
    unique(result$covariate_effects$original_covariate),
    "smoking"
  )
  expect_match(
    result$covariates$reference_profile,
    "reference category"
  )
  expect_false(any(grepl(
    "age|gender|unselected_marker",
    result$covariate_effects$encoded_term,
    ignore.case = TRUE
  )))
})


test_that("summary rejects same-P metadata from a different design matrix", {
  age_data <- prepare_dynamic_covariates("age")
  gender_data <- prepare_dynamic_covariates("gender")
  fit <- mock_dynamic_covariate_fit(age_data)
  attr(fit, "mira_fit_info") <- list(
    P = age_data$P,
    covariate_names = age_data$covariate_names
  )

  expect_identical(age_data$P, gender_data$P)
  expect_error(
    mira_summary(fit, stan_data = gender_data, verbose = FALSE),
    "Covariate names/order.*do not match"
  )
})


test_that("F: continuous age is centered/scaled, never dichotomized", {
  raw <- dynamic_covariate_fixture()
  stan_data <- prepare_dynamic_covariates("age")
  map <- stan_data$covariate_map

  expect_identical(map$original_name, "age")
  expect_identical(map$type, "numeric")
  expect_identical(map$encoding, "center_scale")
  expect_equal(map$center, mean(raw$age))
  expect_equal(map$scale, stats::sd(raw$age))
  expect_equal(
    as.numeric(stan_data$X[, "age"]),
    as.numeric(scale(raw$age))
  )
  expect_false(any(grepl("threshold|above", c(colnames(stan_data$X), map$name))))
  expect_match(map$unit, "SD")
})


test_that("G: invalid and reserved selections fail clearly", {
  expect_error(
    prepare_dynamic_covariates("variabile_inesistente"),
    "not found.*variabile_inesistente|variabile_inesistente.*not found"
  )
  expect_error(
    prepare_dynamic_covariates("patient"),
    "Reserved.*patient|patient.*cannot be covariate"
  )
  expect_error(
    prepare_dynamic_covariates(c("auto", "age")),
    "auto.*cannot be combined"
  )
})


test_that("generic prior vectors follow active encoded columns and names", {
  zero <- prepare_dynamic_covariates(character(0))
  one <- prepare_dynamic_covariates("age")
  many <- prepare_dynamic_covariates(c("age", "gender", "smoking"))

  for (stan_data in list(zero, one, many)) {
    prior <- mira_prior(stan_data)
    prior_data <- mira_prior_stan_data(prior)
    for (field in c(
      "covariate_baseline_prior_sd",
      "beta_covariate_prior_sd",
      "tau_covariate_prior_rate"
    )) {
      expect_length(prior_data[[field]], stan_data$P)
      expect_true(all(is.finite(prior_data[[field]])))
      expect_true(all(prior_data[[field]] > 0))
      # R-side priors retain readable matching names; the Stan-only payload
      # intentionally strips attributes before JSON serialization.
      expect_identical(names(prior[[field]]), stan_data$covariate_names)
      expect_null(names(prior_data[[field]]))
    }
  }

  named_override <- stats::setNames(c(2, 3, 4, 5), many$covariate_names)
  direct <- mira_prior(many, covariate_baseline_sd = named_override)
  expect_equal(direct$covariate_baseline_prior_sd, named_override)

  nested <- mira_prior(
    many,
    custom_prior = list(
      covariates = list(
        age = list(baseline_sd = 2.5),
        smoking = list(tau_rate = 7)
      )
    )
  )
  expect_equal(unname(nested$covariate_baseline_prior_sd["age"]), 2.5)
  expect_equal(
    unname(nested$tau_covariate_prior_rate[
      c("smokingFormer", "smokingCurrent")
    ]),
    c(7, 7)
  )
})


test_that("a prior cannot be reused for different encoded terms with the same P", {
  age_data <- prepare_dynamic_covariates("age")
  gender_data <- prepare_dynamic_covariates("gender")
  age_prior <- mira_prior(age_data)

  expect_identical(age_data$P, gender_data$P)
  expect_false(identical(
    age_data$covariate_names,
    gender_data$covariate_names
  ))
  expect_error(
    mira_fit(
      gender_data,
      prior = age_prior,
      chains = 1,
      parallel_chains = 1,
      iter_warmup = 0,
      iter_sampling = 1,
      refresh = 0,
      stan_file = "not_reached.stan",
      verbose = FALSE
    ),
    "prior covariate names/order do not match"
  )

  embedded <- age_data
  embedded_prior <- mira_prior_stan_data(age_prior)
  embedded_prior$covariate_baseline_prior_sd <- stats::setNames(
    embedded_prior$covariate_baseline_prior_sd,
    "different_term"
  )
  embedded_prior$beta_covariate_prior_sd <- stats::setNames(
    embedded_prior$beta_covariate_prior_sd,
    "different_term"
  )
  embedded_prior$tau_covariate_prior_rate <- stats::setNames(
    embedded_prior$tau_covariate_prior_rate,
    "different_term"
  )
  embedded[names(embedded_prior)] <- embedded_prior

  expect_error(
    mira_fit(
      embedded,
      chains = 1,
      parallel_chains = 1,
      iter_warmup = 0,
      iter_sampling = 1,
      refresh = 0,
      stan_file = "not_reached.stan",
      verbose = FALSE
    ),
    "embedded prior covariate names/order do not match"
  )
})


test_that("H and I: Stan is generic, legacy-free, and statically P = 0 safe", {
  stan_file <- dynamic_covariate_stan_file()
  expect_true(file.exists(stan_file))
  stan_lines <- readLines(stan_file, warn = FALSE)
  stan_code <- paste(stan_lines, collapse = "\n")

  expect_match(stan_code, "int<lower=0> P;", fixed = TRUE)
  expect_match(stan_code, "matrix[S, P] X;", fixed = TRUE)
  expect_match(stan_code, "vector[P] covariate_baseline_effect;", fixed = TRUE)
  expect_match(stan_code, "vector[P] beta_covariate_time;", fixed = TRUE)
  expect_match(stan_code, "matrix[P, K - 1] z_covariate_step;", fixed = TRUE)
  expect_match(stan_code, "vector<lower=0>[P] tau_covariate;", fixed = TRUE)
  expect_match(stan_code, "matrix[P, K] covariate_effect;", fixed = TRUE)

  legacy_identifiers <- c(
    "male", "age_above_threshold", "gender_effect",
    "age_threshold_effect", "gender_baseline_effect",
    "age_baseline_effect", "beta_gender_time", "beta_age_time",
    "z_gender_step", "z_age_step", "tau_gender", "tau_age",
    "male_vs_female_difference", "older_vs_younger_difference"
  )
  for (identifier in legacy_identifiers) {
    expect_false(
      grepl(paste0("\\b", identifier, "\\b"), stan_code, perl = TRUE),
      info = paste("legacy Stan identifier remains:", identifier)
    )
  }

  p_loops <- grep("for \\(p in 1:P\\)", stan_lines)
  expect_gt(length(p_loops), 0L)
  for (line_number in p_loops) {
    preceding <- stan_lines[seq.int(max(1L, line_number - 2L), line_number - 1L)]
    expect_true(
      any(grepl("if \\(P > 0\\)", preceding)),
      info = paste("unguarded 1:P loop at Stan line", line_number)
    )
  }
})


test_that("H: opt-in CmdStan smoke test accepts P = 0 data", {
  run_smoke <- tolower(Sys.getenv("MIRA_RUN_CMDSTAN_TESTS", unset = "false")) %in%
    c("1", "true", "yes")
  skip_if_not(run_smoke, "Set MIRA_RUN_CMDSTAN_TESTS=true for CmdStan smoke tests")
  skip_if_not_installed("cmdstanr")

  cmdstan_path <- tryCatch(
    cmdstanr::cmdstan_path(),
    error = function(e) NA_character_
  )
  skip_if(
    length(cmdstan_path) != 1L || is.na(cmdstan_path) ||
      !nzchar(cmdstan_path) || !dir.exists(cmdstan_path),
    "CmdStan is not configured or its installation directory is inaccessible"
  )

  cmdstan_version <- tryCatch(
    cmdstanr::cmdstan_version(error_on_NA = FALSE),
    error = function(e) NA_character_
  )
  skip_if(
    length(cmdstan_version) < 1L || anyNA(cmdstan_version),
    "CmdStan is not configured"
  )

  stan_data <- prepare_dynamic_covariates(character(0))
  prior <- mira_prior(stan_data)
  sampling_data <- dynamic_covariate_sampling_data(stan_data, prior)
  expect_identical(sampling_data$P, 0L)
  expect_identical(dim(sampling_data$X), c(sampling_data$S, 0L))

  temporary_stan <- file.path(tempdir(), "mira_dynamic_covariates_smoke.stan")
  expect_true(file.copy(
    dynamic_covariate_stan_file(), temporary_stan, overwrite = TRUE
  ))
  model <- cmdstanr::cmdstan_model(temporary_stan, quiet = TRUE)
  expect_s3_class(model, "CmdStanModel")

  fit <- model$sample(
    data = sampling_data,
    seed = 20260923,
    chains = 1,
    parallel_chains = 1,
    iter_warmup = 20,
    iter_sampling = 1,
    refresh = 0
  )
  expect_s3_class(fit, "CmdStanMCMC")
})
