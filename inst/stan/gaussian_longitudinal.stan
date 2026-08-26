data {

  // ============================================================
  // DATA DIMENSIONS
  // ============================================================

  int<lower=1> N;
  int<lower=1> S;
  int<lower=2> K;


  // ============================================================
  // OBSERVATIONS
  // ============================================================

  vector[N] y;


  // ============================================================
  // SUBJECT INDEX
  // ============================================================

  array[N] int<lower=1, upper=S> subject;


  // ============================================================
  // TIME INDEX
  // ============================================================

  array[N] int<lower=1, upper=K> time;


  // ============================================================
  // ACTUAL TIME VALUES
  //
  // Example:
  // c(0, 3, 5, 12)
  // ============================================================

  vector[K] time_value;


  // ============================================================
  // CLINICAL THRESHOLD
  // ============================================================

  real meaningful_change;


  // ============================================================
  // PRIOR HYPERPARAMETERS
  // ============================================================

  real mu_time_prior_mean;

  real<lower=0> mu_time_prior_sd;


  real<lower=0>
    sigma_intercept_prior_rate;


  real<lower=0>
    sigma_slope_prior_rate;


  real<lower=0>
    sigma_prior_rate;


  real<lower=0>
    nu_prior_shape;


  real<lower=0>
    nu_prior_rate;
}


parameters {

  // ============================================================
  // POPULATION-LEVEL TIME MEANS
  //
  // One mean for every measurement occasion.
  //
  // K = 2  -> mu_time[1:2]
  // K = 4  -> mu_time[1:4]
  // K = 5  -> mu_time[1:5]
  // ============================================================

  vector[K] mu_time;


  // ============================================================
  // SUBJECT-LEVEL RANDOM EFFECTS
  //
  // Row 1 = random intercept
  // Row 2 = random slope
  // ============================================================

  matrix[2, S] z_subject;


  vector<lower=0>[2] sigma_subject;


  cholesky_factor_corr[2] L_subject;


  // ============================================================
  // RESIDUAL STANDARD DEVIATION
  // ============================================================

  real<lower=0> sigma;


  // ============================================================
  // STUDENT-T DEGREES OF FREEDOM
  // ============================================================

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
  // POPULATION TIME MEANS
  // ============================================================

  mu_time ~ normal(
    mu_time_prior_mean,
    mu_time_prior_sd
  );


  // ============================================================
  // NON-CENTERED RANDOM EFFECTS
  // ============================================================

  to_vector(z_subject) ~ std_normal();


  // ============================================================
  // RANDOM INTERCEPT SD
  // ============================================================

  sigma_subject[1] ~ exponential(
    sigma_intercept_prior_rate
  );


  // ============================================================
  // RANDOM SLOPE SD
  // ============================================================

  sigma_subject[2] ~ exponential(
    sigma_slope_prior_rate
  );


  // ============================================================
  // RANDOM EFFECT CORRELATION
  // ============================================================

  L_subject ~ lkj_corr_cholesky(2);


  // ============================================================
  // RESIDUAL SD
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
  //
  // For each observation:
  //
  // population mean at that time
  // + subject random intercept
  // + subject random slope * actual time
  //
  // This works for ANY K.
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


  // ============================================================
  // LIKELIHOOD
  // ============================================================

  y ~ student_t(
    nu,
    mu,
    sigma
  );
}


generated quantities {

  // ============================================================
  // POPULATION CHANGE FROM BASELINE
  //
  // change_from_baseline[k]
  //
  // k = 1 -> 0 by definition
  // k = 2 -> time 2 - baseline
  // k = 3 -> time 3 - baseline
  // ...
  // ============================================================

  vector[K] change_from_baseline;


  // ============================================================
  // POPULATION CHANGE BETWEEN CONSECUTIVE TIMES
  //
  // consecutive_change[1] =
  //   time 2 - time 1
  //
  // consecutive_change[2] =
  //   time 3 - time 2
  //
  // etc.
  // ============================================================

  vector[K - 1] consecutive_change;


  // ============================================================
  // POPULATION SLOPE BETWEEN BASELINE AND EACH TIME
  // ============================================================

  vector[K] population_slope_from_baseline;


  // ============================================================
  // STANDARDIZED CHANGE FROM BASELINE
  // ============================================================

  vector[K] standardized_change_from_baseline;


  // ============================================================
  // POPULATION CLINICAL IMPROVEMENT
  // ============================================================

  vector[K] population_any_improvement;

  vector[K] population_meaningful_improvement;


  // ============================================================
  // POPULATION CHANGE MINUS MCID
  // ============================================================

  vector[K] population_change_minus_MCID;


  // ============================================================
  // SUBJECT-SPECIFIC SLOPES
  // ============================================================

  vector[S] individual_slope;


  // ============================================================
  // SUBJECT-SPECIFIC CHANGES FROM BASELINE
  //
  // Matrix:
  //
  // row = time
  // col = subject
  // ============================================================

  matrix[K, S] individual_change_from_baseline;


  // ============================================================
  // SUBJECT-SPECIFIC CONSECUTIVE CHANGES
  // ============================================================

  matrix[K - 1, S] individual_consecutive_change;


  // ============================================================
  // INDIVIDUAL CLINICAL RESPONSE
  // ============================================================

  matrix[K, S] individual_any_improvement;

  matrix[K, S] individual_meaningful_improvement;


  // ============================================================
  // INDIVIDUAL CHANGE MINUS MCID
  // ============================================================

  matrix[K, S] individual_change_minus_MCID;


  // ============================================================
  // POPULATION RESPONDER PROPORTION
  // ============================================================

  vector[K] population_responder_proportion;


  // ============================================================
  // RANDOM EFFECT PARAMETERS
  // ============================================================

  real sigma_intercept;

  real sigma_slope;

  real rho_subject;


  // ============================================================
  // LOG LIKELIHOOD
  // ============================================================

  vector[N] log_lik;


  // ============================================================
  // POSTERIOR PREDICTIVE OBSERVATIONS
  // ============================================================

  vector[N] y_rep;


  // ============================================================
  // POPULATION CHANGES
  // ============================================================

  for (k in 1:K) {

    change_from_baseline[k] =
      mu_time[k] - mu_time[1];
  }


  // ============================================================
  // CONSECUTIVE CHANGES
  // ============================================================

  for (k in 1:(K - 1)) {

    consecutive_change[k] =
      mu_time[k + 1]
      -
      mu_time[k];
  }


  // ============================================================
  // POPULATION SLOPES
  //
  // Slope from baseline to each time point.
  //
  // At baseline the slope is set to 0.
  // ============================================================

  population_slope_from_baseline[1] = 0;


  for (k in 2:K) {

    population_slope_from_baseline[k] =
      (
        mu_time[k]
        -
        mu_time[1]
      )
      /
      (
        time_value[k]
        -
        time_value[1]
      );
  }


  // ============================================================
  // STANDARDIZED POPULATION CHANGE
  // ============================================================

  for (k in 1:K) {

    standardized_change_from_baseline[k] =
      change_from_baseline[k]
      /
      sigma;
  }


  // ============================================================
  // POPULATION CLINICAL QUANTITIES
  // ============================================================

  for (k in 1:K) {

    population_any_improvement[k] =
      (change_from_baseline[k] > 0)
      ? 1 : 0;


    population_meaningful_improvement[k] =
      (
        change_from_baseline[k]
        >= meaningful_change
      )
      ? 1 : 0;


    population_change_minus_MCID[k] =
      change_from_baseline[k]
      -
      meaningful_change;
  }


  // ============================================================
  // RANDOM EFFECT PARAMETERS
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
  // INDIVIDUAL SLOPES
  // ============================================================

  for (s in 1:S) {

    individual_slope[s] =
      population_slope_from_baseline[K]
      +
      b_subject[2, s];
  }


  // ============================================================
  // INDIVIDUAL CHANGES
  // ============================================================

  for (s in 1:S) {

    // ----------------------------------------------------------
    // Changes from baseline
    // ----------------------------------------------------------

    for (k in 1:K) {

      real mu_k;
      real mu_baseline;


      mu_k =
        mu_time[k]
        +
        b_subject[1, s]
        +
        b_subject[2, s]
          * time_value[k];


      mu_baseline =
        mu_time[1]
        +
        b_subject[1, s]
        +
        b_subject[2, s]
          * time_value[1];


      individual_change_from_baseline[k, s] =
        mu_k
        -
        mu_baseline;


      // --------------------------------------------------------
      // Clinical response
      // --------------------------------------------------------

      individual_any_improvement[k, s] =
        (
          individual_change_from_baseline[k, s]
          > 0
        )
        ? 1 : 0;


      individual_meaningful_improvement[k, s] =
        (
          individual_change_from_baseline[k, s]
          >= meaningful_change
        )
        ? 1 : 0;


      individual_change_minus_MCID[k, s] =
        individual_change_from_baseline[k, s]
        -
        meaningful_change;
    }


    // ----------------------------------------------------------
    // Consecutive changes
    // ----------------------------------------------------------

    for (k in 1:(K - 1)) {

      individual_consecutive_change[k, s] =
        individual_change_from_baseline[k + 1, s]
        -
        individual_change_from_baseline[k, s];
    }
  }


  // ============================================================
  // POPULATION RESPONDER PROPORTION
  // ============================================================

  for (k in 1:K) {

    population_responder_proportion[k] =
      mean(
        individual_meaningful_improvement[k]
      );
  }


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


    // ----------------------------------------------------------
    // Log likelihood
    // ----------------------------------------------------------

    log_lik[i] =
      student_t_lpdf(
        y[i] |
        nu,
        mu_i,
        sigma
      );


    // ----------------------------------------------------------
    // Posterior predictive observation
    // ----------------------------------------------------------

    y_rep[i] =
      student_t_rng(
        nu,
        mu_i,
        sigma
      );
  }
}
