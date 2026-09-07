functions {
  // Point masses at declared endpoints represent observable floor/ceiling
  // censoring. With no endpoint this reduces to the ordinary density.
  real censored_student_t_lpdf(
      real y,
      real nu,
      real mu,
      real sigma,
      int has_lower_bound,
      real lower_bound,
      int has_upper_bound,
      real upper_bound) {
    if (has_lower_bound == 1 && y <= lower_bound)
      return student_t_lcdf(lower_bound | nu, mu, sigma);

    if (has_upper_bound == 1 && y >= upper_bound)
      return student_t_lccdf(upper_bound | nu, mu, sigma);

    return student_t_lpdf(y | nu, mu, sigma);
  }

  real censored_normal_lpdf(
      real y,
      real mu,
      real sigma,
      int has_lower_bound,
      real lower_bound,
      int has_upper_bound,
      real upper_bound) {
    if (has_lower_bound == 1 && y <= lower_bound)
      return normal_lcdf(lower_bound | mu, sigma);

    if (has_upper_bound == 1 && y >= upper_bound)
      return normal_lccdf(upper_bound | mu, sigma);

    return normal_lpdf(y | mu, sigma);
  }

  real censored_lognormal_lpdf(
      real y,
      real mu_log,
      real sigma_log,
      int has_lower_bound,
      real lower_bound,
      int has_upper_bound,
      real upper_bound) {
    if (has_lower_bound == 1 && y <= lower_bound)
      return lognormal_lcdf(lower_bound | mu_log, sigma_log);

    if (has_upper_bound == 1 && y >= upper_bound)
      return lognormal_lccdf(upper_bound | mu_log, sigma_log);

    return lognormal_lpdf(y | mu_log, sigma_log);
  }

  real apply_observable_bounds(
      real value,
      int has_lower_bound,
      real lower_bound,
      int has_upper_bound,
      real upper_bound) {
    real bounded_value = value;

    if (has_lower_bound == 1)
      bounded_value = fmax(bounded_value, lower_bound);

    if (has_upper_bound == 1)
      bounded_value = fmin(bounded_value, upper_bound);

    return bounded_value;
  }

  // For log-normal outcomes this is the conditional median; for identity
  // models it is the latent location. Observable bounds are then applied.
  real natural_latent_location(
      real eta,
      int likelihood_id,
      int has_lower_bound,
      real lower_bound,
      int has_upper_bound,
      real upper_bound) {
    real value = likelihood_id == 3 ? exp(eta) : eta;

    return apply_observable_bounds(
      value,
      has_lower_bound,
      lower_bound,
      has_upper_bound,
      upper_bound
    );
  }

  // For the log-normal family this integrates over observation noise and
  // Gaussian subject effects. For censored identity models it is the bounded
  // latent location (the exact censored Student-t mean has no simple form).
  real natural_population_mean(
      real eta,
      real sigma,
      real random_effect_variance,
      int likelihood_id,
      int has_lower_bound,
      real lower_bound,
      int has_upper_bound,
      real upper_bound) {
    real value;

    if (likelihood_id == 3)
      value = exp(eta + 0.5 * (square(sigma) + random_effect_variance));
    else
      value = eta;

    return apply_observable_bounds(
      value,
      has_lower_bound,
      lower_bound,
      has_upper_bound,
      upper_bound
    );
  }

  real censored_outcome_rng(
      int likelihood_id,
      real eta,
      real sigma,
      real nu,
      int has_lower_bound,
      real lower_bound,
      int has_upper_bound,
      real upper_bound) {
    real value;

    if (likelihood_id == 1)
      value = student_t_rng(nu, eta, sigma);
    else if (likelihood_id == 2)
      value = normal_rng(eta, sigma);
    else
      value = lognormal_rng(eta, sigma);

    return apply_observable_bounds(
      value,
      has_lower_bound,
      lower_bound,
      has_upper_bound,
      upper_bound
    );
  }
}

data {
  // ============================================================
  // DIMENSIONS AND OUTCOME FAMILY
  // ============================================================
  int<lower=1> N;
  int<lower=1> S;
  int<lower=2> K;
  int<lower=2> G;

  // 1 = censored Student-t with identity link
  // 2 = censored Gaussian with identity link
  // 3 = censored log-normal with log link
  int<lower=1, upper=3> likelihood_id;

  int<lower=0, upper=1> has_lower_bound;
  real outcome_lower_bound;
  int<lower=0, upper=1> has_upper_bound;
  real outcome_upper_bound;

  vector[N] y;
  array[N] int<lower=1, upper=S> subject;
  array[N] int<lower=1, upper=K> time;
  array[S] int<lower=1, upper=G> arm;
  array[S] int<lower=0, upper=1> male;
  array[S] int<lower=0, upper=1> age_above_threshold;
  vector[K] time_value;
  int<lower=-1, upper=1> direction;

  // Clinical thresholds are always in natural outcome units.
  real<lower=0> mcid_prior_mean;
  real<lower=0> mcid_prior_sd;
  real<lower=0> meaningful_between_arm_difference;

  // Externally configurable priors on the active link scale.
  real baseline_prior_mean;
  real<lower=0> baseline_prior_sd;
  real beta_time_prior_mean;
  real<lower=0> beta_time_prior_sd;
  real<lower=0> tau_common_prior_rate;
  real<lower=0> beta_treatment_prior_sd;
  real<lower=0> tau_treatment_prior_rate;
  real<lower=0> arm_baseline_sd_prior_rate;
  real<lower=0> gender_baseline_prior_sd;
  real<lower=0> beta_gender_prior_sd;
  real<lower=0> tau_gender_prior_rate;
  real<lower=0> age_baseline_prior_sd;
  real<lower=0> beta_age_prior_sd;
  real<lower=0> tau_age_prior_rate;
  real<lower=0> sigma_intercept_prior_rate;
  real<lower=0> sigma_slope_prior_rate;
  real<lower=0> sigma_prior_rate;
  real<lower=0> nu_prior_shape;
  real<lower=0> nu_prior_rate;
}

transformed data {
  vector[K - 1] dt;
  real mcid_prior_shape;
  real mcid_prior_rate;
  int<lower=0, upper=1> number_student_t_parameters;

  number_student_t_parameters = likelihood_id == 1;

  if (direction == 0)
    reject("direction must be +1 or -1");
  if (mcid_prior_mean <= 0)
    reject("mcid_prior_mean must be strictly positive");
  if (mcid_prior_sd <= 0)
    reject("mcid_prior_sd must be strictly positive");
  if (has_lower_bound == 1 && has_upper_bound == 1 &&
      outcome_lower_bound >= outcome_upper_bound)
    reject("outcome_lower_bound must be smaller than outcome_upper_bound");

  for (i in 1:N) {
    if (has_lower_bound == 1 && y[i] < outcome_lower_bound)
      reject("An observed outcome is below outcome_lower_bound");
    if (has_upper_bound == 1 && y[i] > outcome_upper_bound)
      reject("An observed outcome is above outcome_upper_bound");
    if (likelihood_id == 3 && y[i] <= 0)
      reject("The log-normal likelihood requires strictly positive outcomes");
  }

  mcid_prior_shape = square(mcid_prior_mean / mcid_prior_sd);
  mcid_prior_rate = mcid_prior_mean / square(mcid_prior_sd);

  for (k in 2:K) {
    if (time_value[k] <= time_value[k - 1])
      reject("time_value must be strictly increasing");
    dt[k - 1] = time_value[k] - time_value[k - 1];
  }
}

parameters {
  // Reference-arm trajectory on the active link scale.
  real baseline_mean;
  real beta_time;
  vector[K - 1] z_common_step;
  real<lower=0> tau_common;

  // Treatment change is anchored to zero at baseline.
  vector[G - 1] beta_treatment;
  matrix[G - 1, K - 1] z_treatment_step;
  vector<lower=0>[G - 1] tau_treatment;
  vector[G - 1] z_arm_baseline;
  real<lower=0> arm_baseline_sd;

  // Male - Female trajectory.
  real gender_baseline_effect;
  real beta_gender_time;
  vector[K - 1] z_gender_step;
  real<lower=0> tau_gender;

  // Age above - age at/below threshold trajectory.
  real age_baseline_effect;
  real beta_age_time;
  vector[K - 1] z_age_step;
  real<lower=0> tau_age;

  // Correlated subject random intercept and slope on the link scale.
  matrix[2, S] z_subject;
  vector<lower=0>[2] sigma_subject;
  cholesky_factor_corr[2] L_subject;

  // Natural-scale residual scale for identity models; log-SD for CMT.
  real<lower=0> sigma;

  // Dimension zero outside the Student-t branch, so Gaussian/log-normal
  // models do not sample an unidentified degrees-of-freedom parameter.
  vector<lower=2>[number_student_t_parameters] nu;

  // Externally informed MCID on the natural outcome scale.
  real<lower=0> mcid;
}

transformed parameters {
  vector[K] mu_reference;
  matrix[G - 1, K] treatment_change;
  vector[G - 1] arm_baseline_offset;
  vector[K] gender_effect;
  vector[K] age_threshold_effect;
  matrix[2, S] b_subject;

  mu_reference[1] = baseline_mean;
  for (k in 2:K) {
    mu_reference[k] =
      mu_reference[k - 1]
      + beta_time * dt[k - 1]
      + tau_common * sqrt(dt[k - 1]) * z_common_step[k - 1];
  }

  for (g in 1:(G - 1)) {
    treatment_change[g, 1] = 0;
    for (k in 2:K) {
      treatment_change[g, k] =
        treatment_change[g, k - 1]
        + beta_treatment[g] * dt[k - 1]
        + tau_treatment[g] * sqrt(dt[k - 1])
          * z_treatment_step[g, k - 1];
    }
  }

  arm_baseline_offset = arm_baseline_sd * z_arm_baseline;

  gender_effect[1] = gender_baseline_effect;
  for (k in 2:K) {
    gender_effect[k] =
      gender_effect[k - 1]
      + beta_gender_time * dt[k - 1]
      + tau_gender * sqrt(dt[k - 1]) * z_gender_step[k - 1];
  }

  age_threshold_effect[1] = age_baseline_effect;
  for (k in 2:K) {
    age_threshold_effect[k] =
      age_threshold_effect[k - 1]
      + beta_age_time * dt[k - 1]
      + tau_age * sqrt(dt[k - 1]) * z_age_step[k - 1];
  }

  b_subject = diag_pre_multiply(sigma_subject, L_subject) * z_subject;
}

model {
  baseline_mean ~ normal(baseline_prior_mean, baseline_prior_sd);
  beta_time ~ normal(beta_time_prior_mean, beta_time_prior_sd);
  z_common_step ~ std_normal();
  tau_common ~ exponential(tau_common_prior_rate);

  beta_treatment ~ normal(0, beta_treatment_prior_sd);
  to_vector(z_treatment_step) ~ std_normal();
  tau_treatment ~ exponential(tau_treatment_prior_rate);
  z_arm_baseline ~ std_normal();
  arm_baseline_sd ~ exponential(arm_baseline_sd_prior_rate);

  gender_baseline_effect ~ normal(0, gender_baseline_prior_sd);
  beta_gender_time ~ normal(0, beta_gender_prior_sd);
  z_gender_step ~ std_normal();
  tau_gender ~ exponential(tau_gender_prior_rate);

  age_baseline_effect ~ normal(0, age_baseline_prior_sd);
  beta_age_time ~ normal(0, beta_age_prior_sd);
  z_age_step ~ std_normal();
  tau_age ~ exponential(tau_age_prior_rate);

  to_vector(z_subject) ~ std_normal();
  sigma_subject[1] ~ exponential(sigma_intercept_prior_rate);
  sigma_subject[2] ~ exponential(sigma_slope_prior_rate);
  L_subject ~ lkj_corr_cholesky(2);

  sigma ~ exponential(sigma_prior_rate);
  if (number_student_t_parameters == 1)
    nu[1] ~ gamma(nu_prior_shape, nu_prior_rate);

  mcid ~ gamma(mcid_prior_shape, mcid_prior_rate);

  // Family-specific observation model.
  for (i in 1:N) {
    int s = subject[i];
    int k = time[i];
    int g = arm[s];
    real eta = mu_reference[k];

    if (g > 1)
      eta += arm_baseline_offset[g - 1] + treatment_change[g - 1, k];

    eta +=
      male[s] * gender_effect[k]
      + age_above_threshold[s] * age_threshold_effect[k]
      + b_subject[1, s]
      + b_subject[2, s] * (time_value[k] - time_value[1]);

    if (likelihood_id == 1) {
      target += censored_student_t_lpdf(
        y[i] |
        nu[1], eta, sigma,
        has_lower_bound, outcome_lower_bound,
        has_upper_bound, outcome_upper_bound
      );
    } else if (likelihood_id == 2) {
      target += censored_normal_lpdf(
        y[i] |
        eta, sigma,
        has_lower_bound, outcome_lower_bound,
        has_upper_bound, outcome_upper_bound
      );
    } else {
      target += censored_lognormal_lpdf(
        y[i] |
        eta, sigma,
        has_lower_bound, outcome_lower_bound,
        has_upper_bound, outcome_upper_bound
      );
    }
  }
}

generated quantities {
  // Family and variance quantities.
  real nu_value = 0;
  real sigma_intercept = sigma_subject[1];
  real sigma_slope = sigma_subject[2];
  matrix[2, 2] subject_correlation =
    multiply_lower_tri_self_transpose(L_subject);
  matrix[2, 2] subject_cholesky =
    diag_pre_multiply(sigma_subject, L_subject);
  matrix[2, 2] subject_covariance =
    multiply_lower_tri_self_transpose(subject_cholesky);
  real rho_subject = subject_correlation[1, 2];
  real residual_sd;
  real residual_cv = 0;

  // Population estimands.
  matrix[G, K] population_model_location;
  matrix[G, K] population_median;
  matrix[G, K] population_mean;
  matrix[G, K] population_change_from_baseline;
  matrix[G, K] directional_population_change;
  matrix[G, K] standardized_population_change;
  matrix[G, K] population_ratio_from_baseline;
  matrix[G, K] population_percent_change_from_baseline;

  vector[K] male_vs_female_difference;
  vector[K] male_vs_female_change_difference;
  vector[K] older_vs_younger_difference;
  vector[K] older_vs_younger_change_difference;
  vector[K] directional_male_vs_female_change_difference;
  vector[K] directional_older_vs_younger_change_difference;

  // Posterior means of these binary draws estimate new-subject response
  // probabilities while retaining the complete generative model.
  matrix[G, K] new_subject_latent_any_improvement_draw;
  matrix[G, K] new_subject_latent_change_draw;
  matrix[G, K] new_subject_latent_responder_draw;
  matrix[G, K] new_subject_predictive_any_improvement_draw;
  matrix[G, K] new_subject_predictive_change_draw;
  matrix[G, K] new_subject_predictive_responder_draw;

  matrix[K, S] individual_change_from_baseline;
  matrix[K, S] individual_directional_change;
  matrix[K, S] individual_any_improvement_draw;
  matrix[K, S] individual_meaningful_responder_draw;
  matrix[K, S] individual_change_minus_mcid;

  // Baseline-standardized treatment effects on natural-scale change.
  matrix[G - 1, K] treatment_change_difference;
  matrix[G - 1, K] directional_treatment_benefit;
  matrix[G - 1, K] treatment_benefit_positive_draw;
  matrix[G - 1, K] treatment_benefit_meaningful_draw;
  matrix[G - 1, K] latent_responder_probability_difference;
  matrix[G - 1, K] treatment_ratio_of_ratios;

  vector[N] log_lik;
  vector[N] y_rep;

  if (number_student_t_parameters == 1) {
    nu_value = nu[1];
    residual_sd = sigma * sqrt(nu[1] / (nu[1] - 2));
  } else {
    // For log-normal fits this is the residual SD on log(Y).
    residual_sd = sigma;
  }

  if (likelihood_id == 3)
    residual_cv = sqrt(exp(square(sigma)) - 1);

  // Natural-scale population trajectories. For log-normal outcomes,
  // population_mean integrates over subject effects and observation noise;
  // population_median is the zero-random-effect conditional median.
  for (g in 1:G) {
    real baseline_offset = 0;

    if (g > 1)
      baseline_offset = arm_baseline_offset[g - 1];

    for (k in 1:K) {
      real elapsed = time_value[k] - time_value[1];
      real treatment_effect = 0;
      real eta;
      real random_effect_variance =
        subject_covariance[1, 1]
        + 2 * elapsed * subject_covariance[1, 2]
        + square(elapsed) * subject_covariance[2, 2];

      if (g > 1)
        treatment_effect = treatment_change[g - 1, k];

      eta = mu_reference[k] + baseline_offset + treatment_effect;
      population_model_location[g, k] = eta;
      population_median[g, k] = natural_latent_location(
        eta,
        likelihood_id,
        has_lower_bound, outcome_lower_bound,
        has_upper_bound, outcome_upper_bound
      );
      population_mean[g, k] = natural_population_mean(
        eta,
        sigma,
        fmax(random_effect_variance, 0),
        likelihood_id,
        has_lower_bound, outcome_lower_bound,
        has_upper_bound, outcome_upper_bound
      );
    }

    for (k in 1:K) {
      population_change_from_baseline[g, k] =
        population_mean[g, k] - population_mean[g, 1];
      directional_population_change[g, k] =
        direction * population_change_from_baseline[g, k];

      if (likelihood_id == 3) {
        standardized_population_change[g, k] =
          (population_model_location[g, k]
           - population_model_location[g, 1]) / sigma;
        population_ratio_from_baseline[g, k] =
          population_mean[g, k] / population_mean[g, 1];
        population_percent_change_from_baseline[g, k] =
          100 * (population_ratio_from_baseline[g, k] - 1);
      } else {
        standardized_population_change[g, k] =
          population_change_from_baseline[g, k] / residual_sd;
        population_ratio_from_baseline[g, k] = 1;
        population_percent_change_from_baseline[g, k] = 0;
      }
    }
  }

  // Natural-scale covariate contrasts at the reference arm.
  for (k in 1:K) {
    real elapsed = time_value[k] - time_value[1];
    real random_effect_variance =
      subject_covariance[1, 1]
      + 2 * elapsed * subject_covariance[1, 2]
      + square(elapsed) * subject_covariance[2, 2];
    real reference_mean = natural_population_mean(
      mu_reference[k], sigma, fmax(random_effect_variance, 0), likelihood_id,
      has_lower_bound, outcome_lower_bound,
      has_upper_bound, outcome_upper_bound
    );
    real male_mean = natural_population_mean(
      mu_reference[k] + gender_effect[k],
      sigma, fmax(random_effect_variance, 0), likelihood_id,
      has_lower_bound, outcome_lower_bound,
      has_upper_bound, outcome_upper_bound
    );
    real older_mean = natural_population_mean(
      mu_reference[k] + age_threshold_effect[k],
      sigma, fmax(random_effect_variance, 0), likelihood_id,
      has_lower_bound, outcome_lower_bound,
      has_upper_bound, outcome_upper_bound
    );

    male_vs_female_difference[k] = male_mean - reference_mean;
    older_vs_younger_difference[k] = older_mean - reference_mean;
  }

  for (k in 1:K) {
    male_vs_female_change_difference[k] =
      male_vs_female_difference[k] - male_vs_female_difference[1];
    older_vs_younger_change_difference[k] =
      older_vs_younger_difference[k] - older_vs_younger_difference[1];
    directional_male_vs_female_change_difference[k] =
      direction * male_vs_female_change_difference[k];
    directional_older_vs_younger_change_difference[k] =
      direction * older_vs_younger_change_difference[k];
  }

  // New-subject latent and posterior-predictive responder draws. A full
  // correlated random intercept/slope draw is required for absolute changes
  // under the CMT log link, where the intercept does not cancel.
  for (g in 1:G) {
    vector[2] new_b = multi_normal_cholesky_rng(
      rep_vector(0, 2), subject_cholesky
    );
    real baseline_eta = population_model_location[g, 1] + new_b[1];
    real baseline_latent = natural_latent_location(
      baseline_eta,
      likelihood_id,
      has_lower_bound, outcome_lower_bound,
      has_upper_bound, outcome_upper_bound
    );
    real baseline_observed;
    real nu_for_rng = 10;

    if (number_student_t_parameters == 1)
      nu_for_rng = nu[1];

    baseline_observed = censored_outcome_rng(
      likelihood_id, baseline_eta, sigma, nu_for_rng,
      has_lower_bound, outcome_lower_bound,
      has_upper_bound, outcome_upper_bound
    );

    new_subject_latent_any_improvement_draw[g, 1] = 0;
    new_subject_latent_change_draw[g, 1] = 0;
    new_subject_latent_responder_draw[g, 1] = 0;
    new_subject_predictive_any_improvement_draw[g, 1] = 0;
    new_subject_predictive_change_draw[g, 1] = 0;
    new_subject_predictive_responder_draw[g, 1] = 0;

    for (k in 2:K) {
      real elapsed = time_value[k] - time_value[1];
      real followup_eta =
        population_model_location[g, k]
        + new_b[1]
        + new_b[2] * elapsed;
      real followup_latent = natural_latent_location(
        followup_eta,
        likelihood_id,
        has_lower_bound, outcome_lower_bound,
        has_upper_bound, outcome_upper_bound
      );
      real followup_observed = censored_outcome_rng(
        likelihood_id, followup_eta, sigma, nu_for_rng,
        has_lower_bound, outcome_lower_bound,
        has_upper_bound, outcome_upper_bound
      );
      real latent_change = followup_latent - baseline_latent;
      real predictive_change = followup_observed - baseline_observed;

      new_subject_latent_any_improvement_draw[g, k] =
        direction * latent_change > 0;
      new_subject_latent_change_draw[g, k] = latent_change;
      new_subject_latent_responder_draw[g, k] =
        direction * latent_change >= mcid;
      new_subject_predictive_any_improvement_draw[g, k] =
        direction * predictive_change > 0;
      new_subject_predictive_change_draw[g, k] = predictive_change;
      new_subject_predictive_responder_draw[g, k] =
        direction * predictive_change >= mcid;
    }
  }

  // Existing-subject latent changes on the natural outcome scale.
  for (s in 1:S) {
    int g = arm[s];
    real baseline_eta =
      population_model_location[g, 1]
      + male[s] * gender_effect[1]
      + age_above_threshold[s] * age_threshold_effect[1]
      + b_subject[1, s];
    real baseline_latent = natural_latent_location(
      baseline_eta,
      likelihood_id,
      has_lower_bound, outcome_lower_bound,
      has_upper_bound, outcome_upper_bound
    );

    for (k in 1:K) {
      real elapsed = time_value[k] - time_value[1];
      real eta =
        population_model_location[g, k]
        + male[s] * gender_effect[k]
        + age_above_threshold[s] * age_threshold_effect[k]
        + b_subject[1, s]
        + b_subject[2, s] * elapsed;
      real latent_outcome = natural_latent_location(
        eta,
        likelihood_id,
        has_lower_bound, outcome_lower_bound,
        has_upper_bound, outcome_upper_bound
      );
      real latent_change = latent_outcome - baseline_latent;
      real directional_change = direction * latent_change;

      individual_change_from_baseline[k, s] = latent_change;
      individual_directional_change[k, s] = directional_change;
      individual_any_improvement_draw[k, s] = directional_change > 0;
      individual_meaningful_responder_draw[k, s] =
        directional_change >= mcid;
      individual_change_minus_mcid[k, s] = directional_change - mcid;
    }
  }

  // Treatment effects at a common baseline (arm imbalance excluded).
  // Absolute differences are primary; ratio-of-ratios is also returned for
  // the log-link model.
  for (g in 1:(G - 1)) {
    for (k in 1:K) {
      real elapsed = time_value[k] - time_value[1];
      real random_effect_variance_k =
        subject_covariance[1, 1]
        + 2 * elapsed * subject_covariance[1, 2]
        + square(elapsed) * subject_covariance[2, 2];
      real random_effect_variance_0 = subject_covariance[1, 1];
      real reference_0 = natural_population_mean(
        mu_reference[1], sigma, random_effect_variance_0, likelihood_id,
        has_lower_bound, outcome_lower_bound,
        has_upper_bound, outcome_upper_bound
      );
      real reference_k = natural_population_mean(
        mu_reference[k], sigma, fmax(random_effect_variance_k, 0), likelihood_id,
        has_lower_bound, outcome_lower_bound,
        has_upper_bound, outcome_upper_bound
      );
      real treatment_0 = reference_0;
      real treatment_k = natural_population_mean(
        mu_reference[k] + treatment_change[g, k],
        sigma, fmax(random_effect_variance_k, 0), likelihood_id,
        has_lower_bound, outcome_lower_bound,
        has_upper_bound, outcome_upper_bound
      );
      real difference_in_change =
        (treatment_k - treatment_0) - (reference_k - reference_0);
      real directional_difference = direction * difference_in_change;

      treatment_change_difference[g, k] = difference_in_change;
      directional_treatment_benefit[g, k] = directional_difference;
      treatment_benefit_positive_draw[g, k] = directional_difference > 0;
      treatment_benefit_meaningful_draw[g, k] =
        directional_difference >= meaningful_between_arm_difference;
      latent_responder_probability_difference[g, k] =
        new_subject_latent_responder_draw[g + 1, k]
        - new_subject_latent_responder_draw[1, k];

      if (likelihood_id == 3)
        treatment_ratio_of_ratios[g, k] = exp(treatment_change[g, k]);
      else
        treatment_ratio_of_ratios[g, k] = 1;
    }
  }

  // Pointwise log likelihood and exact family-specific replication.
  for (i in 1:N) {
    int s = subject[i];
    int k = time[i];
    int g = arm[s];
    real eta = mu_reference[k];
    real nu_for_rng = 10;

    if (g > 1)
      eta += arm_baseline_offset[g - 1] + treatment_change[g - 1, k];

    eta +=
      male[s] * gender_effect[k]
      + age_above_threshold[s] * age_threshold_effect[k]
      + b_subject[1, s]
      + b_subject[2, s] * (time_value[k] - time_value[1]);

    if (likelihood_id == 1) {
      log_lik[i] = censored_student_t_lpdf(
        y[i] |
        nu[1], eta, sigma,
        has_lower_bound, outcome_lower_bound,
        has_upper_bound, outcome_upper_bound
      );
      nu_for_rng = nu[1];
    } else if (likelihood_id == 2) {
      log_lik[i] = censored_normal_lpdf(
        y[i] |
        eta, sigma,
        has_lower_bound, outcome_lower_bound,
        has_upper_bound, outcome_upper_bound
      );
    } else {
      log_lik[i] = censored_lognormal_lpdf(
        y[i] |
        eta, sigma,
        has_lower_bound, outcome_lower_bound,
        has_upper_bound, outcome_upper_bound
      );
    }

    y_rep[i] = censored_outcome_rng(
      likelihood_id, eta, sigma, nu_for_rng,
      has_lower_bound, outcome_lower_bound,
      has_upper_bound, outcome_upper_bound
    );
  }
}
