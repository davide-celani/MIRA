generated quantities {

  // ============================================================
  // 1. POPULATION EFFECTS
  // ============================================================

  real change_T1;
  real change_T2;
  real change_T2_vs_T1;

  real population_slope;

  real standardized_change_T1;
  real standardized_change_T2;


  // ============================================================
  // 2. POPULATION CLINICAL RESPONSE
  // ============================================================

  int population_any_improvement_T1;
  int population_any_improvement_T2;

  int population_meaningful_improvement_T1;
  int population_meaningful_improvement_T2;

  int population_meaningful_deterioration_T1;
  int population_meaningful_deterioration_T2;

  real population_distance_from_MCID_T1;
  real population_distance_from_MCID_T2;


  // ============================================================
  // 3. CLINICAL STATES
  // ============================================================

  int population_clinical_benefit_T2;
  int population_clinical_no_change_T2;
  int population_clinical_harm_T2;


  // ============================================================
  // 4. VARIANCE COMPONENTS
  // ============================================================

  real sigma_intercept;
  real sigma_slope;
  real rho_subject;


  // ============================================================
  // 5. INDIVIDUAL TRAJECTORIES
  // ============================================================

  vector[S] individual_slope;

  vector[S] individual_change_T1;
  vector[S] individual_change_T2;
  vector[S] individual_change_T2_vs_T1;


  // ============================================================
  // 6. INDIVIDUAL CLINICAL RESPONSE
  // ============================================================

  vector[S] individual_any_improvement_T1;
  vector[S] individual_any_improvement_T2;

  vector[S] individual_meaningful_improvement_T1;
  vector[S] individual_meaningful_improvement_T2;

  vector[S] individual_meaningful_deterioration_T1;
  vector[S] individual_meaningful_deterioration_T2;


  // ============================================================
  // 7. LONGITUDINAL RESPONSE PATTERNS
  // ============================================================

  vector[S] individual_early_responder;
  vector[S] individual_late_responder;
  vector[S] individual_sustained_responder;
  vector[S] individual_loss_of_response;


  // ============================================================
  // 8. DISTANCE FROM CLINICAL THRESHOLD
  // ============================================================

  vector[S] individual_distance_from_MCID_T1;
  vector[S] individual_distance_from_MCID_T2;


  // ============================================================
  // 9. POPULATION RESPONDER QUANTITIES
  // ============================================================

  real responder_proportion_T1;
  real responder_proportion_T2;

  real deterioration_proportion_T1;
  real deterioration_proportion_T2;

  real early_responder_proportion;
  real late_responder_proportion;
  real sustained_responder_proportion;
  real loss_of_response_proportion;


  // ============================================================
  // 10. FUTURE-PATIENT QUANTITIES
  // ============================================================

  vector[2] z_new;
  vector[2] b_new;

  real new_patient_change_T2;
  int new_patient_responder_T2;
  int new_patient_deterioration_T2;


  // ============================================================
  // 11. MODEL DIAGNOSTICS
  // ============================================================

  vector[N] log_lik;
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
  // STANDARDIZED EFFECTS
  // ============================================================

  standardized_change_T1 =
    change_T1 / sigma;

  standardized_change_T2 =
    change_T2 / sigma;


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
  // POPULATION CLINICAL RESPONSE
  // ============================================================

  population_any_improvement_T1 =
    (change_T1 > 0);

  population_any_improvement_T2 =
    (change_T2 > 0);


  population_meaningful_improvement_T1 =
    (change_T1 >= meaningful_change);

  population_meaningful_improvement_T2 =
    (change_T2 >= meaningful_change);


  population_meaningful_deterioration_T1 =
    (change_T1 <= -meaningful_change);

  population_meaningful_deterioration_T2 =
    (change_T2 <= -meaningful_change);


  population_distance_from_MCID_T1 =
    change_T1 - meaningful_change;

  population_distance_from_MCID_T2 =
    change_T2 - meaningful_change;


  // ============================================================
  // POPULATION CLINICAL STATES
  // ============================================================

  population_clinical_benefit_T2 =
    (change_T2 >= meaningful_change);

  population_clinical_harm_T2 =
    (change_T2 <= -meaningful_change);

  population_clinical_no_change_T2 =
    (fabs(change_T2) < meaningful_change);


  // ============================================================
  // SUBJECT-SPECIFIC QUANTITIES
  // ============================================================

  for (s in 1:S) {

    real mu_T0;
    real mu_T1;
    real mu_T2;


    individual_slope[s] =
      population_slope
      +
      b_subject[2, s];


    mu_T0 =
      mu_time[1]
      + b_subject[1, s]
      + b_subject[2, s] * time_value[1];

    mu_T1 =
      mu_time[2]
      + b_subject[1, s]
      + b_subject[2, s] * time_value[2];

    mu_T2 =
      mu_time[3]
      + b_subject[1, s]
      + b_subject[2, s] * time_value[3];


    // ----------------------------------------------------------
    // CHANGES
    // ----------------------------------------------------------

    individual_change_T1[s] =
      mu_T1 - mu_T0;

    individual_change_T2[s] =
      mu_T2 - mu_T0;

    individual_change_T2_vs_T1[s] =
      mu_T2 - mu_T1;


    // ----------------------------------------------------------
    // IMPROVEMENT
    // ----------------------------------------------------------

    individual_any_improvement_T1[s] =
      (individual_change_T1[s] > 0);

    individual_any_improvement_T2[s] =
      (individual_change_T2[s] > 0);


    individual_meaningful_improvement_T1[s] =
      (individual_change_T1[s] >= meaningful_change);

    individual_meaningful_improvement_T2[s] =
      (individual_change_T2[s] >= meaningful_change);


    // ----------------------------------------------------------
    // DETERIORATION
    // ----------------------------------------------------------

    individual_meaningful_deterioration_T1[s] =
      (individual_change_T1[s] <= -meaningful_change);

    individual_meaningful_deterioration_T2[s] =
      (individual_change_T2[s] <= -meaningful_change);


    // ----------------------------------------------------------
    // DISTANCE FROM MCID
    // ----------------------------------------------------------

    individual_distance_from_MCID_T1[s] =
      individual_change_T1[s] - meaningful_change;

    individual_distance_from_MCID_T2[s] =
      individual_change_T2[s] - meaningful_change;


    // ----------------------------------------------------------
    // LONGITUDINAL RESPONSE PATTERNS
    // ----------------------------------------------------------

    individual_early_responder[s] =
      (individual_change_T1[s] >= meaningful_change);


    individual_late_responder[s] =
      (individual_change_T1[s] < meaningful_change)
      &&
      (individual_change_T2[s] >= meaningful_change);


    individual_sustained_responder[s] =
      (individual_change_T1[s] >= meaningful_change)
      &&
      (individual_change_T2[s] >= meaningful_change);


    individual_loss_of_response[s] =
      (individual_change_T1[s] >= meaningful_change)
      &&
      (individual_change_T2[s] < meaningful_change);
  }


  // ============================================================
  // OBSERVED POPULATION RESPONSE PROPORTIONS
  // ============================================================

  responder_proportion_T1 =
    mean(individual_meaningful_improvement_T1);

  responder_proportion_T2 =
    mean(individual_meaningful_improvement_T2);


  deterioration_proportion_T1 =
    mean(individual_meaningful_deterioration_T1);

  deterioration_proportion_T2 =
    mean(individual_meaningful_deterioration_T2);


  early_responder_proportion =
    mean(individual_early_responder);

  late_responder_proportion =
    mean(individual_late_responder);

  sustained_responder_proportion =
    mean(individual_sustained_responder);

  loss_of_response_proportion =
    mean(individual_loss_of_response);


  // ============================================================
  // NEW PATIENT FROM THE POPULATION
  // ============================================================

  z_new[1] = normal_rng(0, 1);
  z_new[2] = normal_rng(0, 1);

  b_new =
    diag_pre_multiply(
      sigma_subject,
      L_subject
    ) * z_new;


  new_patient_change_T2 =
    change_T2
    +
    b_new[2]
      * (time_value[3] - time_value[1]);


  new_patient_responder_T2 =
    (new_patient_change_T2 >= meaningful_change);


  new_patient_deterioration_T2 =
    (new_patient_change_T2 <= -meaningful_change);


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
