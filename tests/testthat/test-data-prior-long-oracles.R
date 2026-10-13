data_prior_oracle_fixture <- function() {
  data.frame(
    patient = c("P3", "P1", "P2"),
    arm = c("Treatment", "Control", "Treatment"),
    score_t2 = c(9, 18, 36),
    age = c(30, 40, 50),
    score_t0 = c(1, 2, 4),
    eligible = c(FALSE, TRUE, FALSE),
    score_t1 = c(3, 6, 12),
    stringsAsFactors = FALSE
  )
}

data_prior_oracle_prepare <- function(data = data_prior_oracle_fixture(), ...) {
  defaults <- list(data = data, time_value = c(0, 2, 9), outcome = "generic",
    likelihood = "gaussian", direction = "higher", meaningful_change = 6,
    meaningful_change_sd = 2, reference_arm = "Control", covariates = character())
  additions <- list(...)
  defaults[names(additions)] <- additions
  do.call(mira_data_long, defaults)
}

test_that("Stan observations retain subject, visit, arm and irregular-time association", {
  x <- data_prior_oracle_prepare()
  expect_identical(c(x$N, x$S, x$K, x$G, x$P), c(9L, 3L, 3L, 2L, 0L))
  expect_identical(x$measurement_columns, c("score_t0", "score_t1", "score_t2"))
  expect_identical(x$y, c(1, 2, 4, 3, 6, 12, 9, 18, 36))
  expect_identical(x$subject, c(1L, 2L, 3L, 1L, 2L, 3L, 1L, 2L, 3L))
  expect_identical(x$time, c(1L, 1L, 1L, 2L, 2L, 2L, 3L, 3L, 3L))
  expect_identical(x$time_value, c(0, 2, 9))
  expect_identical(x$subject_labels, c("P3", "P1", "P2"))
  expect_identical(x$arm_labels, c("Control", "Treatment"))
  expect_identical(x$arm, c(2L, 1L, 2L))
  expect_identical(x$reference_arm, "Control")
  expect_identical(x$mcid_prior_mean, 6)
  expect_identical(x$mcid_prior_sd, 2)
  expect_identical(x$meaningful_between_arm_difference, 6)
  # The nine hand-listed observations sum to 91 and their squares to 1911.
  # Sample variance = (1911 - 91^2 / 9) / 8 = 4459 / 36.
  expect_equal(x$mean_y, 91 / 9, tolerance = 1e-12)
  expect_equal(x$sd_y, sqrt(4459) / 6, tolerance = 1e-12)
  expect_identical(x$mean_model_y, x$mean_y)
  expect_identical(x$sd_model_y, x$sd_y)
})

test_that("permuting rows and columns preserves keyed Stan observations with fixed reference", {
  original <- data_prior_oracle_prepare()
  raw <- data_prior_oracle_fixture()
  permuted <- data_prior_oracle_prepare(raw[c(3L, 1L, 2L), rev(seq_along(raw))])
  # With explicit reference and outcome, permutations only relabel internal
  # subject indices: aligning by patient and visit must recover identical data.
  keyed_observations <- function(x) {
    out <- data.frame(patient = x$subject_labels[x$subject], visit = x$time,
      value = x$y, arm = x$arm_labels[x$arm[x$subject]])
    out <- out[order(out$patient, out$visit), , drop = FALSE]
    rownames(out) <- NULL
    out
  }
  expect_identical(keyed_observations(permuted), keyed_observations(original))
  expect_identical(permuted$time_value, original$time_value)
  expect_identical(permuted$measurement_columns, original$measurement_columns)
  expect_identical(permuted$mcid_prior_mean, original$mcid_prior_mean)
})

test_that("lognormal preparation separates natural and log moments and retains endpoints", {
  x <- data_prior_oracle_prepare(likelihood = "lognormal", outcome_bounds = c(1, 36))
  expect_identical(x$likelihood_id, 3L)
  expect_identical(x$y, c(1, 2, 4, 3, 6, 12, 9, 18, 36))
  expect_equal(x$mean_y, 91 / 9, tolerance = 1e-12)
  expect_equal(x$sd_y, sqrt(4459) / 6, tolerance = 1e-12)
  # log(y) = i*log(2) + j*log(3), the complete i,j = 0,1,2 grid.
  # The sample log variance is (3/4)*(log(2)^2 + log(3)^2).
  expect_equal(x$mean_model_y, log(6), tolerance = 1e-12)
  expect_equal(x$sd_model_y, sqrt(3 / 4 * (log(2)^2 + log(3)^2)), tolerance = 1e-12)
  expect_identical(c(x$has_lower_bound, x$has_upper_bound), c(1L, 1L))
  expect_identical(c(x$outcome_lower_bound, x$outcome_upper_bound), c(1, 36))
  expect_identical(x$outcome_diagnostics$lower_boundary_count, 1L)
  expect_identical(x$outcome_diagnostics$upper_boundary_count, 1L)
  expect_equal(x$outcome_diagnostics$lower_boundary_fraction, 1 / 9, tolerance = 1e-12)
  expect_equal(x$outcome_diagnostics$upper_boundary_fraction, 1 / 9, tolerance = 1e-12)
  expect_identical(c(x$mcid_prior_mean, x$mcid_prior_sd), c(6, 2))
  expect_identical(x$change_scale, "absolute natural-outcome units")
})

test_that("Bayesian preparation rejects unsupported missing and infinite observations", {
  for (bad in c(NA_real_, NaN, Inf, -Inf)) {
    raw <- data_prior_oracle_fixture()
    raw$score_t1[2L] <- bad
    expect_error(data_prior_oracle_prepare(raw),
      "measurement columns contain missing or non-finite values: score_t1",
      info = paste("Invalid observation:", format(bad)))
  }
  raw <- data_prior_oracle_fixture()
  raw$patient[2L] <- raw$patient[1L]
  expect_error(data_prior_oracle_prepare(raw), "Duplicated patient identifiers")
  raw$patient[2L] <- NA_character_
  expect_error(data_prior_oracle_prepare(raw), "`patient` contains missing values")
})

data_prior_oracle_prior <- function() {
  x <- data_prior_oracle_prepare(covariates = c("age", "eligible"))
  mira_prior_long(x, baseline_mean = -3, baseline_sd = 11, beta_time_mean = -2,
    beta_time_sd = 12, tau_common_rate = 13, beta_treatment_sd = 14,
    tau_treatment_rate = 15, arm_baseline_sd_rate = 16,
    covariate_baseline_sd = c(2, 3), beta_covariate_sd = c(4, 5),
    tau_covariate_rate = c(6, 7), sigma_intercept_rate = 17,
    sigma_slope_rate = 18, sigma_rate = 19, nu_shape = 20, nu_rate = 21)
}

test_that("prior conversion independently maps every declared Stan field without metadata", {
  prior <- data_prior_oracle_prior()
  expected <- list(baseline_prior_mean = -3, baseline_prior_sd = 11,
    beta_time_prior_mean = -2, beta_time_prior_sd = 12, tau_common_prior_rate = 13,
    beta_treatment_prior_sd = 14, tau_treatment_prior_rate = 15,
    arm_baseline_sd_prior_rate = 16, covariate_baseline_prior_sd = c(2, 3),
    beta_covariate_prior_sd = c(4, 5), tau_covariate_prior_rate = c(6, 7),
    sigma_intercept_prior_rate = 17, sigma_slope_prior_rate = 18,
    sigma_prior_rate = 19, nu_prior_shape = 20, nu_prior_rate = 21)
  expect_identical(mira_prior_stan_data_long(prior), expected)
  expect_identical(names(prior$covariate_baseline_prior_sd), prior$covariate_names)
  expect_identical(unname(prior$covariate_baseline_prior_sd), c(2, 3))
  expect_identical(mira_validate_prior_long(prior), TRUE)
})

test_that("prior validation rejects invalid values and coherent-dimension metadata mismatches", {
  prior <- data_prior_oracle_prior()
  for (bad in list(NA_real_, NaN, Inf, -Inf, TRUE, numeric(), c(1, 2))) {
    broken <- prior
    broken$baseline_prior_mean <- bad
    expect_error(mira_validate_prior_long(broken), "Invalid prior field `baseline_prior_mean`")
  }
  for (bad in list(0, -1, NA_real_, NaN, Inf, -Inf, TRUE)) {
    broken <- prior
    broken$sigma_prior_rate <- bad
    expect_error(mira_validate_prior_long(broken), "Invalid prior field `sigma_prior_rate`")
  }
  for (bad in list(c(0, 2), c(-1, 2), c(NA_real_, 2), c(NaN, 2), c(Inf, 2), c(-Inf, 2))) {
    broken <- prior
    broken$beta_covariate_prior_sd <- bad
    expect_error(mira_validate_prior_long(broken), "Invalid prior field `beta_covariate_prior_sd`")
  }
  broken <- prior
  broken$covariate_baseline_prior_sd <- 1
  expect_error(mira_validate_prior_long(broken), "three generic covariate prior vectors must have the same length")
  broken <- prior
  broken$P <- 3L
  expect_error(mira_validate_prior_long(broken), "Prior metadata `P` does not match")
  broken <- prior
  broken$covariate_names <- rev(broken$covariate_names)
  expect_error(mira_validate_prior_long(broken), "Names of prior field .* do not match `covariate_names`")
  broken <- prior
  broken$nu_prior_shape <- NULL
  expect_error(mira_prior_stan_data_long(broken), "Missing fields: nu_prior_shape")
})
