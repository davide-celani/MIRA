data {
  // ============================================================
  // DIMENSIONS
  // ============================================================
  int<lower=1> N;                  // observations
  int<lower=1> S;                  // subjects
  int<lower=2> K;                  // measurement occasions
  int<lower=2> G;                  // treatment groups; group 1 = reference

  // ============================================================
  // OBSERVATIONS AND INDICES
  // ============================================================
  vector[N] y;
  array[N] int<lower=1, upper=S> subject;
  array[N] int<lower=1, upper=K> time;
  array[S] int<lower=1, upper=G> arm;

  // Strictly increasing actual times, e.g. 0, 3, 5, 12
  vector[K] time_value;

  // +1: higher outcome = better
  // -1: lower outcome = better
  int<lower=-1, upper=1> direction;

  // ============================================================
  // CLINICAL THRESHOLDS
  // ============================================================
  // MCID is treated as externally informed but uncertain.
  // Its posterior equals this prior because no likelihood term
  // identifies it directly; uncertainty is propagated into all
  // responder quantities.
  real<lower=0> mcid_prior_mean;
  real<lower=0> mcid_prior_sd;

  // Threshold for a clinically relevant BETWEEN-ARM difference.
  // This is intentionally distinct from the individual-level MCID.
  real<lower=0> meaningful_between_arm_difference;

  // ============================================================
  // PRIORS: POPULATION TRAJECTORY
  // ============================================================
  real baseline_prior_mean;
  real<lower=0> baseline_prior_sd;

  real beta_time_prior_mean;
  real<lower=0> beta_time_prior_sd;

  real<lower=0> tau_common_prior_rate;

  // Treatment trajectory priors (deviation from reference arm)
  real<lower=0> beta_treatment_prior_sd;
  real<lower=0> tau_treatment_prior_rate;

  // Baseline arm imbalance (reference arm fixed to zero)
  real<lower=0> arm_baseline_sd_prior_rate;

  // ============================================================
  // PRIORS: SUBJECT RANDOM EFFECTS
  // ============================================================
  real<lower=0> sigma_intercept_prior_rate;
  real<lower=0> sigma_slope_prior_rate;

  // ============================================================
  // PRIORS: RESIDUAL MODEL
  // ============================================================
  real<lower=0> sigma_prior_rate;
  real<lower=0> nu_prior_shape;
  real<lower=0> nu_prior_rate;
}

transformed data {
  vector[K - 1] dt;

  if (direction == 0)
    reject("direction must be +1 or -1");

  for (k in 2:K) {
    if (time_value[k] <= time_value[k - 1])
      reject("time_value must be strictly increasing");

    dt[k - 1] = time_value[k] - time_value[k - 1];
  }
}

parameters {
  // ============================================================
  // REFERENCE-ARM POPULATION TRAJECTORY
  // ============================================================
  real baseline_mean;
  real beta_time;

  // Continuous-time RW1 deviations around the global linear trend.
  // Increment variance scales with elapsed time.
  vector[K - 1] z_common_step;
  real<lower=0> tau_common;

  // ============================================================
  // TREATMENT-BY-TIME TRAJECTORIES
  // ============================================================
  // Each non-reference arm has an overall treatment slope plus
  // time-local deviations. Treatment change is anchored to zero
  // at baseline, so treatment effects are effects on CHANGE.
  vector[G - 1] beta_treatment;
  matrix[G - 1, K - 1] z_treatment_step;
  vector<lower=0>[G - 1] tau_treatment;

  // Allows chance baseline imbalance without contaminating
  // the treatment effect on change.
  vector[G - 1] z_arm_baseline;
  real<lower=0> arm_baseline_sd;

  // ============================================================
  // SUBJECT RANDOM INTERCEPT + RANDOM SLOPE
  // ============================================================
  matrix[2, S] z_subject;
  vector<lower=0>[2] sigma_subject;
  cholesky_factor_corr[2] L_subject;

  // ============================================================
  // ROBUST RESIDUAL MODEL
  // ============================================================
  real<lower=0> sigma;
  real<lower=2> nu;

  // ============================================================
  // EXTERNALLY INFORMED UNCERTAIN MCID
  // ============================================================
  real<lower=0> mcid;
}

transformed parameters {
  vector[K] mu_reference;
  matrix[G - 1, K] treatment_change;
  vector[G - 1] arm_baseline_offset;
  matrix[2, S] b_subject;

  // ------------------------------------------------------------
  // Reference-arm trajectory:
  // global linear component + continuous-time random-walk detail.
  // ------------------------------------------------------------
  mu_reference[1] = baseline_mean;

  for (k in 2:K) {
    mu_reference[k] =
      mu_reference[k - 1]
      + beta_time * dt[k - 1]
      + tau_common * sqrt(dt[k - 1]) * z_common_step[k - 1];
  }

  // ------------------------------------------------------------
  // Treatment effects on CHANGE, anchored at baseline = 0.
  // ------------------------------------------------------------
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

  // ------------------------------------------------------------
  // Correlated subject random effects, non-centered.
  // Time is centered at baseline, so b_subject[1, s] is the
  // subject deviation at baseline and b_subject[2, s] controls
  // subject-specific longitudinal departure.
  // ------------------------------------------------------------
  b_subject =
    diag_pre_multiply(sigma_subject, L_subject) * z_subject;
}

model {
  vector[N] mu;

  // ============================================================
  // PRIORS: REFERENCE TRAJECTORY
  // ============================================================
  baseline_mean ~ normal(baseline_prior_mean, baseline_prior_sd);
  beta_time ~ normal(beta_time_prior_mean, beta_time_prior_sd);

  z_common_step ~ std_normal();
  tau_common ~ exponential(tau_common_prior_rate);

  // ============================================================
  // PRIORS: TREATMENT TRAJECTORIES
  // ============================================================
  beta_treatment ~ normal(0, beta_treatment_prior_sd);
  to_vector(z_treatment_step) ~ std_normal();
  tau_treatment ~ exponential(tau_treatment_prior_rate);

  z_arm_baseline ~ std_normal();
  arm_baseline_sd ~ exponential(arm_baseline_sd_prior_rate);

  // ============================================================
  // PRIORS: SUBJECT RANDOM EFFECTS
  // ============================================================
  to_vector(z_subject) ~ std_normal();
  sigma_subject[1] ~ exponential(sigma_intercept_prior_rate);
  sigma_subject[2] ~ exponential(sigma_slope_prior_rate);
  L_subject ~ lkj_corr_cholesky(2);

  // ============================================================
  // PRIORS: ROBUST RESIDUAL MODEL
  // ============================================================
  sigma ~ exponential(sigma_prior_rate);
  nu ~ gamma(nu_prior_shape, nu_prior_rate);

  // External MCID uncertainty; outcome data do not identify MCID.
  mcid ~ normal(mcid_prior_mean, mcid_prior_sd);

  // ============================================================
  // OBSERVATION MODEL
  // ============================================================
  for (i in 1:N) {
    int s = subject[i];
    int k = time[i];
    int g = arm[s];
    real population_mu = mu_reference[k];

    if (g > 1) {
      population_mu +=
        arm_baseline_offset[g - 1]
        + treatment_change[g - 1, k];
    }

    mu[i] =
      population_mu
      + b_subject[1, s]
      + b_subject[2, s] * (time_value[k] - time_value[1]);
  }

  y ~ student_t(nu, mu, sigma);
}

generated quantities {
  // ============================================================
  // BASIC MODEL QUANTITIES
  // ============================================================
  real sigma_intercept = sigma_subject[1];
  real sigma_slope = sigma_subject[2];
  real rho_subject =
    multiply_lower_tri_self_transpose(L_subject)[1, 2];

  // True SD of Student-t residuals (nu > 2), not merely its scale.
  real residual_sd = sigma * sqrt(nu / (nu - 2));

  // ============================================================
  // POPULATION TRAJECTORIES AND CHANGE
  // ============================================================
  matrix[G, K] population_mean;
  matrix[G, K] population_change_from_baseline;
  matrix[G, K] directional_population_change;
  matrix[G, K] standardized_population_change;

  // ============================================================
  // RESPONDER ESTIMANDS FOR A NEW SUBJECT
  // ============================================================
  // Analytic probability for a NEW latent subject, integrating
  // over the random slope but not residual observation noise.
  matrix[G, K] latent_new_subject_any_improvement_prob;
  matrix[G, K] latent_new_subject_responder_prob;

  // One posterior predictive draw per MCMC draw. Posterior means of
  // the *_responder_draw objects estimate responder probabilities.
  matrix[G, K] new_subject_latent_change_draw;
  matrix[G, K] new_subject_latent_responder_draw;
  matrix[G, K] new_subject_predictive_change_draw;
  matrix[G, K] new_subject_predictive_responder_draw;

  // ============================================================
  // SUBJECT-SPECIFIC LATENT CHANGE
  // ============================================================
  matrix[K, S] individual_change_from_baseline;
  matrix[K, S] individual_directional_change;
  matrix[K, S] individual_any_improvement_draw;
  matrix[K, S] individual_meaningful_responder_draw;
  matrix[K, S] individual_change_minus_mcid;

  // ============================================================
  // TREATMENT EFFECTS ON CHANGE (vs reference arm)
  // ============================================================
  matrix[G - 1, K] treatment_change_difference;
  matrix[G - 1, K] directional_treatment_benefit;
  matrix[G - 1, K] treatment_benefit_positive_draw;
  matrix[G - 1, K] treatment_benefit_meaningful_draw;
  matrix[G - 1, K] latent_responder_probability_difference;

  // ============================================================
  // OBSERVATION-LEVEL DIAGNOSTICS
  // ============================================================
  vector[N] log_lik;
  vector[N] y_rep;

  // ------------------------------------------------------------
  // Population means and changes by arm.
  // ------------------------------------------------------------
  for (g in 1:G) {
    real baseline_offset = 0;

    if (g > 1)
      baseline_offset = arm_baseline_offset[g - 1];

    for (k in 1:K) {
      real trt = 0;

      if (g > 1)
        trt = treatment_change[g - 1, k];

      population_mean[g, k] =
        mu_reference[k] + baseline_offset + trt;
    }

    for (k in 1:K) {
      population_change_from_baseline[g, k] =
        population_mean[g, k] - population_mean[g, 1];

      directional_population_change[g, k] =
        direction * population_change_from_baseline[g, k];

      standardized_population_change[g, k] =
        population_change_from_baseline[g, k] / residual_sd;
    }
  }

  // ------------------------------------------------------------
  // New-subject latent responder probabilities.
  // Random intercept cancels in change; random slope remains.
  // ------------------------------------------------------------
  for (g in 1:G) {
    latent_new_subject_any_improvement_prob[g, 1] = 0;
    latent_new_subject_responder_prob[g, 1] = 0;

    for (k in 2:K) {
      real elapsed = time_value[k] - time_value[1];
      real mean_directional_change =
        directional_population_change[g, k];
      real sd_latent_change = abs(elapsed) * sigma_slope + 1e-12;

      latent_new_subject_any_improvement_prob[g, k] =
        exp(normal_lccdf(0 | mean_directional_change, sd_latent_change));

      latent_new_subject_responder_prob[g, k] =
        exp(normal_lccdf(mcid | mean_directional_change, sd_latent_change));
    }
  }

  // ------------------------------------------------------------
  // Posterior predictive responder draws for one new subject per arm.
  // The same random slope and baseline residual are reused across
  // times to preserve within-subject longitudinal dependence.
  // ------------------------------------------------------------
  for (g in 1:G) {
    real new_random_slope = normal_rng(0, sigma_slope);
    real baseline_error = student_t_rng(nu, 0, sigma);

    new_subject_latent_change_draw[g, 1] = 0;
    new_subject_latent_responder_draw[g, 1] = 0;
    new_subject_predictive_change_draw[g, 1] = 0;
    new_subject_predictive_responder_draw[g, 1] = 0;

    for (k in 2:K) {
      real elapsed = time_value[k] - time_value[1];
      real latent_change =
        population_change_from_baseline[g, k]
        + new_random_slope * elapsed;
      real followup_error = student_t_rng(nu, 0, sigma);
      real predictive_change =
        latent_change + followup_error - baseline_error;

      new_subject_latent_change_draw[g, k] = latent_change;
      new_subject_latent_responder_draw[g, k] =
        (direction * latent_change >= mcid) ? 1 : 0;

      new_subject_predictive_change_draw[g, k] = predictive_change;
      new_subject_predictive_responder_draw[g, k] =
        (direction * predictive_change >= mcid) ? 1 : 0;
    }
  }

  // ------------------------------------------------------------
  // Existing-subject latent responder draws.
  // These are conditional on each subject's posterior random slope.
  // ------------------------------------------------------------
  for (s in 1:S) {
    int g = arm[s];

    for (k in 1:K) {
      real elapsed = time_value[k] - time_value[1];
      real latent_change =
        population_change_from_baseline[g, k]
        + b_subject[2, s] * elapsed;
      real dchange = direction * latent_change;

      individual_change_from_baseline[k, s] = latent_change;
      individual_directional_change[k, s] = dchange;
      individual_any_improvement_draw[k, s] =
        (dchange > 0) ? 1 : 0;
      individual_meaningful_responder_draw[k, s] =
        (dchange >= mcid) ? 1 : 0;
      individual_change_minus_mcid[k, s] = dchange - mcid;
    }
  }

  // ------------------------------------------------------------
  // Treatment effects vs reference arm.
  // Separate continuous effect, positive-effect event,
  // clinically meaningful between-arm event, and responder uplift.
  // ------------------------------------------------------------
  for (g in 1:(G - 1)) {
    for (k in 1:K) {
      real diff =
        population_change_from_baseline[g + 1, k]
        - population_change_from_baseline[1, k];
      real directional_diff = direction * diff;

      treatment_change_difference[g, k] = diff;
      directional_treatment_benefit[g, k] = directional_diff;
      treatment_benefit_positive_draw[g, k] =
        (directional_diff > 0) ? 1 : 0;
      treatment_benefit_meaningful_draw[g, k] =
        (directional_diff >= meaningful_between_arm_difference)
        ? 1 : 0;

      latent_responder_probability_difference[g, k] =
        latent_new_subject_responder_prob[g + 1, k]
        - latent_new_subject_responder_prob[1, k];
    }
  }

  // ------------------------------------------------------------
  // Log-likelihood and posterior predictive observations.
  // ------------------------------------------------------------
  for (i in 1:N) {
    int s = subject[i];
    int k = time[i];
    int g = arm[s];
    real population_mu = mu_reference[k];
    real mu_i;

    if (g > 1) {
      population_mu +=
        arm_baseline_offset[g - 1]
        + treatment_change[g - 1, k];
    }

    mu_i =
      population_mu
      + b_subject[1, s]
      + b_subject[2, s] * (time_value[k] - time_value[1]);

    log_lik[i] = student_t_lpdf(y[i] | nu, mu_i, sigma);
    y_rep[i] = student_t_rng(nu, mu_i, sigma);
  }
}
