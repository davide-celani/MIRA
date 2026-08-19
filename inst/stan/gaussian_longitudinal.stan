data {

  // ============================================================
  // DATA
  // ============================================================

  int<lower=1> N;
  int<lower=1> S;

  vector[N] y;

  array[N] int<lower=1, upper=S> subject;
  array[N] int<lower=1, upper=3> time;

  vector[3] time_value;

  // ------------------------------------------------------------
  // Clinical threshold
  // ------------------------------------------------------------

  real meaningful_change;


  // ============================================================
  // PRIOR HYPERPARAMETERS
  // ============================================================

  real mu_time_prior_mean;
  real<lower=0> mu_time_prior_sd;

  real<lower=0> sigma_intercept_prior_rate;

  real<lower=0> sigma_slope_prior_rate;

  real<lower=0> sigma_prior_rate;

  real<lower=0> nu_prior_shape;
  real<lower=0> nu_prior_rate;
}


parameters {

  // ============================================================
  // POPULATION-LEVEL TRAJECTORY
  // ============================================================

  vector[3] mu_time;


  // ============================================================
  // SUBJECT-LEVEL RANDOM EFFECTS
  // ============================================================

  matrix[2, S] z_subject;

  vector<lower=0>[2] sigma_subject;

  cholesky_factor_corr[2] L_subject;


  // ============================================================
  // RESIDUAL DISTRIBUTION
  // ============================================================

  real<lower=0> sigma;

  real<lower=2> nu;
}


transformed parameters {

  // ============================================================
  // SUBJECT RANDOM EFFECTS
  // ============================================================

  matrix[2, S] b_subject;

  b_subject =
    diag_pre_multiply(
      sigma_subject,
      L_subject
    ) * z_subject;
}


model {

  vector[N] mu;


  // ============================================================
  // POPULATION TRAJECTORY
  // ============================================================

  mu_time ~ normal(
    mu_time_prior_mean,
    mu_time_prior_sd
  );


  // ============================================================
  // HIERARCHICAL SUBJECT EFFECTS
  // ============================================================

  to_vector(z_subject) ~ std_normal();


  sigma_subject[1] ~ exponential(
    sigma_intercept_prior_rate
  );


  sigma_subject[2] ~ exponential(
    sigma_slope_prior_rate
  );


  L_subject ~ lkj_corr_cholesky(2);


  // ============================================================
  // RESIDUAL VARIATION
  // ============================================================

  sigma ~ exponential(
    sigma_prior_rate
  );


  // ============================================================
  // STUDENT-T DEGREES OF FREEDOM
  // ============================================================

  nu ~ gamma(
    nu_prior_shape,
    nu_prior_rate
  );


  // ============================================================
  // OBSERVATION MODEL
  // ============================================================

  for (i in 1:N) {

    mu[i] =
      mu_time[time[i]]
      +
      b_subject[1, subject[i]]
      +
      b_subject[2, subject[i]]
        * time_value[time[i]];
  }


  y ~ student_t(
    nu,
    mu,
    sigma
  );
}


generated quantities {

  // ============================================================
  // POPULATION CONTRASTS
  // ============================================================

  real change_T1;
  real change_T2;
  real change_T2_vs_T1;


  // ============================================================
  // POPULATION SLOPE
  // ============================================================

  real population_slope;


  // ============================================================
  // VARIANCE COMPONENTS
  // ============================================================

  real sigma_intercept;
  real sigma_slope;
  real rho_subject;


  // ============================================================
  // STANDARDIZED EFFECT
  // ============================================================

  real standardized_change_T2;


  // ============================================================
  // POPULATION CLINICAL QUANTITIES
  // ============================================================

  real population_meaningful_improvement_T2;
  real population_any_improvement_T2;

  real population_change_minus_MCID_T2;


  // ============================================================
  // INDIVIDUAL EFFECTS
  // ============================================================

  vector[S] individual_slope;

  vector[S] individual_change_T1;
  vector[S] individual_change_T2;
  vector[S] individual_change_T2_vs_T1;


  // ============================================================
  // INDIVIDUAL CLINICAL RESPONSE
  // ============================================================

  vector[S] individual_any_improvement_T2;

  vector[S] individual_meaningful_improvement_T2;

  vector[S] individual_change_minus_MCID_T2;


  // ============================================================
  // POPULATION RESPONDER PROPORTION
  // ============================================================

  real population_responder_proportion_T2;


  // ============================================================
  // LOG LIKELIHOOD
  // ============================================================

  vector[N] log_lik;


  // ============================================================
  // POSTERIOR PREDICTIVE OBSERVATIONS
  // ============================================================

  vector[N] y_rep;


  // ============================================================
  // POPULATION CONTRASTS
  // ============================================================

  change_T1 =
    mu_time[2] - mu_time[1];


  change_T2 =
    mu_time[3] - mu_time[1];


  change_T2_vs_T1 =
    mu_time[3] - mu_time[2];


  // ============================================================
  // POPULATION SLOPE
  // ============================================================

  population_slope =
    (mu_time[3] - mu_time[1])
    /
    (time_value[3] - time_value[1]);


  // ============================================================
  // VARIANCE COMPONENTS
  // ============================================================

  sigma_intercept =
    sigma_subject[1];


  sigma_slope =
    sigma_subject[2];


  rho_subject =
    multiply_lower_tri_self_transpose(
      L_subject
    )[1, 2];


  // ============================================================
  // STANDARDIZED EFFECT
  // ============================================================

  standardized_change_T2 =
    change_T2 / sigma;


  // ============================================================
  // POPULATION CLINICAL QUANTITIES
  // ============================================================

  population_any_improvement_T2 =
    (change_T2 > 0)
      ? 1 : 0;


  population_meaningful_improvement_T2 =
    (change_T2 >= meaningful_change)
      ? 1 : 0;


  population_change_minus_MCID_T2 =
    change_T2 - meaningful_change;


  // ============================================================
  // SUBJECT-SPECIFIC QUANTITIES
  // ============================================================

  for (s in 1:S) {

    real mu_T0;
    real mu_T1;
    real mu_T2;


    // ----------------------------------------------------------
    // Individual slope
    // ----------------------------------------------------------

    individual_slope[s] =
      population_slope
      +
      b_subject[2, s];


    // ----------------------------------------------------------
    // Subject-specific expected values
    // ----------------------------------------------------------

    mu_T0 =
      mu_time[1]
      +
      b_subject[1, s]
      +
      b_subject[2, s]
        * time_value[1];


    mu_T1 =
      mu_time[2]
      +
      b_subject[1, s]
      +
      b_subject[2, s]
        * time_value[2];


    mu_T2 =
      mu_time[3]
      +
      b_subject[1, s]
      +
      b_subject[2, s]
        * time_value[3];


    // ----------------------------------------------------------
    // Individual changes
    // ----------------------------------------------------------

    individual_change_T1[s] =
      mu_T1 - mu_T0;


    individual_change_T2[s] =
      mu_T2 - mu_T0;


    individual_change_T2_vs_T1[s] =
      mu_T2 - mu_T1;


    // ----------------------------------------------------------
    // Clinical response
    // ----------------------------------------------------------

    individual_any_improvement_T2[s] =
      (individual_change_T2[s] > 0)
        ? 1 : 0;


    individual_meaningful_improvement_T2[s] =
      (individual_change_T2[s] >= meaningful_change)
        ? 1 : 0;


    individual_change_minus_MCID_T2[s] =
      individual_change_T2[s]
      - meaningful_change;
  }


  // ============================================================
  // POPULATION RESPONDER PROPORTION
  // ============================================================

  population_responder_proportion_T2 =
    mean(
      individual_meaningful_improvement_T2
    );


  // ============================================================
  // OBSERVATION-LEVEL QUANTITIES
  // ============================================================

  for (i in 1:N) {

    real mu_i;


    mu_i =
      mu_time[time[i]]
      +
      b_subject[1, subject[i]]
      +
      b_subject[2, subject[i]]
        * time_value[time[i]];


    log_lik[i] =
      student_t_lpdf(
        y[i] |
        nu,
        mu_i,
        sigma
      );


    y_rep[i] =
      student_t_rng(
        nu,
        mu_i,
        sigma
      );
  }
}
