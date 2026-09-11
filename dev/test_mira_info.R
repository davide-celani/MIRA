

devtools::load_all()



# ------------------------------------------------------------------
# Example 1: Minimal automatic workflow
# ------------------------------------------------------------------
# The input is wide: one row per subject and one numeric column per visit.
# The names subject_id and score_t0, score_t1, score_t2 are deliberately
# conventional, so ID, outcome family, and time order can all be inferred.
set.seed(202601)
n_wide <- 48L
score_baseline <- rnorm(n_wide, mean = 50, sd = 7)
wide <- data.frame(
  subject_id = sprintf("S%03d", seq_len(n_wide)),
  score_t0 = score_baseline,
  score_t1 = score_baseline - 2 + rnorm(n_wide, sd = 1.8),
  score_t2 = score_baseline - 5 + rnorm(n_wide, sd = 2.0)
)

# analyses = "none" keeps this first workflow fast. It disables optional
# modules, but not descriptives, changes, missingness, variability,
# trajectories, or the base-R Friedman branch. Holm adjustment is applied
# separately within each implemented comparison family.
auto_fit <- mira_info(
  wide,
  analyses = "none",
  p_adjust_method = "holm",
  verbose = FALSE
)

# Start with sample size and completeness, then verify the ordered mapping
# from source columns to display labels.
auto_fit$overview[c(
  "n_patients", "n_timepoints", "complete_profiles",
  "complete_profiles_pct"
)]
data.frame(
  order = seq_along(auto_fit$time_vars),
  variable = auto_fit$time_vars,
  label = unname(auto_fit$time_labels)
)
utils::head(auto_fit$long_data[, c(
  "patient", "outcome", "time", "time_index", "time_label", "value"
)])

# Per-timepoint availability and descriptive statistics.
auto_fit$descriptives[, c(
  "label", "n", "missing", "mean", "sd", "ci_lower", "ci_upper"
)]

# Every ordered to-minus-from comparison is in $change. The adjusted columns
# belong to the paired-t and paired-Wilcoxon families, respectively.
auto_fit$change[, c(
  "from_label", "to_label", "n", "mean_change", "cohens_dz",
  "paired_t_p", "paired_t_p_adj", "wilcoxon_p_adj"
)]

# Select the clinically common baseline-to-final row without relying on row
# position. This remains reliable if more visits are added later.
baseline_final <- auto_fit$change[
  auto_fit$change$from == auto_fit$time_vars[[1L]] &
    auto_fit$change$to == auto_fit$time_vars[[length(auto_fit$time_vars)]],
  ,
  drop = FALSE
]
baseline_final[, c(
  "n", "mean_from", "mean_to", "mean_change", "ci_lower", "ci_upper"
)]

# Variability and the preferred available ICC are kept in a one-row table.
auto_fit$variability[, c(
  "between_subject_sd", "within_subject_sd",
  "ICC_anova_complete_profiles", "ICC_model", "ICC"
)]

# $trajectories contains one baseline-to-final record per input row.
utils::head(auto_fit$trajectories[, c(
  "patient", "baseline", "final", "absolute_change", "direction"
)])
table(auto_fit$trajectories$direction, useNA = "ifany")

# Friedman is attempted with base R even though model fitting was disabled.
auto_fit$advanced_tests$friedman$tidy[, c(
  "statistic", "df", "p_raw", "kendalls_w", "n", "n_timepoints"
)]
auto_fit$advanced_tests$friedman$posthoc$tidy[, c(
  "from_label", "to_label", "n", "p_raw", "p_holm", "rank_biserial"
)]

# ------------------------------------------------------------------
# Example 2: Audit what mira_info() detected before fitting anything
# ------------------------------------------------------------------
# This reusable trial has a genuine time-by-arm difference plus two
# time-invariant covariates. Eighty subjects and four visits also provide a
# stable basis for the model examples below.
set.seed(202602)
n_trial <- 80L
trial_arm <- rep(c("Control", "Treatment"), each = n_trial / 2L)
trial_site <- factor(rep(c("North", "South", "East", "West"), length.out = n_trial))
trial_age <- rep(50:69, length.out = n_trial) +
  sample(-1:1, n_trial, replace = TRUE)
treatment_indicator <- as.numeric(trial_arm == "Treatment")
subject_intercept <- rnorm(n_trial, sd = 5.5)
subject_slope <- rnorm(n_trial, sd = 0.8)
site_shift <- c(North = 0, South = 1.0, East = -0.8, West = 0.5)
trial_baseline <- 55 + 0.12 * (trial_age - 60) +
  unname(site_shift[as.character(trial_site)]) + subject_intercept +
  rnorm(n_trial, sd = 1.4)
trial <- data.frame(
  subject_id = sprintf("T%03d", seq_len(n_trial)),
  arm = factor(trial_arm, levels = c("Control", "Treatment")),
  age = trial_age,
  site = trial_site,
  outcome_t0 = trial_baseline,
  outcome_t1 = trial_baseline + subject_slope - 1 -
    3 * treatment_indicator + rnorm(n_trial, sd = 1.5),
  outcome_t2 = trial_baseline + 2 * subject_slope - 2 -
    6 * treatment_indicator + rnorm(n_trial, sd = 1.5),
  outcome_t3 = trial_baseline + 3 * subject_slope - 3 -
    9 * treatment_indicator + rnorm(n_trial, sd = 1.5)
)

# inspect_only validates and resolves configuration but performs no outcome
# analysis. The requested optional switches are still recorded in $config.
detected <- mira_info(
  trial,
  analyses = c("model", "arm_tests"),
  inspect_only = TRUE,
  verbose = FALSE
)

# These are the principal resolved choices. Control is selected as the
# reference because its label is recognized as control-like.
detected$config[c(
  "id", "id_generated", "outcomes", "arm", "reference_arm", "covariates"
)]
detected$config$time_vars[["outcome"]]
detected$config$time_labels[["outcome"]]
detected$config$analyses

# These fields distinguish inferred choices from explicit overrides and
# record why an ID, arm, direction, or threshold was chosen.
detected$config$auto_detected
detected$config$specified_manually
detected$config$sources
detected$config$alternatives

# The parsed map is the compact audit trail for outcome grouping and time
# ordering. Candidate tables and covariate classes explain the other choices.
detected$detected_variables$longitudinal_map[, c(
  "outcome", "variable", "time_label", "time_order", "pattern", "numeric"
)]
detected$detected_variables$id_candidates
detected$detected_variables$arm_candidates
detected$detected_variables$covariates_numeric
detected$detected_variables$covariates_categorical

# Empty vectors are informative: they mean that no unresolved alternatives or
# diagnostic warnings remained. Non-empty entries explain what needs review.
detected$config$alternatives
detected$diagnostics$warnings

# ------------------------------------------------------------------
# Example 3: Take control when column names do not encode time
# ------------------------------------------------------------------
set.seed(202603)
n_explicit <- 42L
symptom_baseline <- rnorm(n_explicit, mean = 30, sd = 5)
explicit_data <- data.frame(
  case_id = sprintf("E%03d", seq_len(n_explicit)),
  pre = symptom_baseline,
  week4 = symptom_baseline - 2 + rnorm(n_explicit, sd = 1.2),
  final = symptom_baseline - 4 + rnorm(n_explicit, sd = 1.4)
)

# None of pre, week4, and final contains an outcome stem separated from a
# recognized suffix. Supply the ID, ordered columns, display labels, and
# outcome name explicitly. The supplied order defines baseline and final;
# stable_threshold = 1 treats changes within +/-1 symptom point as stable.
explicit_fit <- mira_info(
  explicit_data,
  id = "case_id",
  time_vars = c("pre", "week4", "final"),
  time_labels = c("Baseline", "Week 4", "Week 12"),
  outcomes = "symptom",
  covariates = character(0),
  analyses = "none",
  improvement_direction = "lower",
  stable_threshold = 1,
  strict_id = TRUE,
  verbose = FALSE
)

data.frame(
  order = seq_along(explicit_fit$time_vars),
  variable = explicit_fit$time_vars,
  label = unname(explicit_fit$time_labels)
)
explicit_fit$config$outcomes
explicit_fit$config$specified_manually[c(
  "id", "outcomes", "time_vars", "time_labels"
)]
explicit_fit$change[
  explicit_fit$change$from == "pre" & explicit_fit$change$to == "final",
  c(
    "from_label", "to_label", "n", "mean_change",
    "improved_n", "stable_n", "worsened_n", "paired_t_p_adj"
  ),
  drop = FALSE
]

# Wide format means one row per subject. strict_id = TRUE is the protective
# default; these quick checks should both be zero before analysis.
c(
  missing_ids = sum(is.na(explicit_data$case_id)),
  duplicated_ids = sum(duplicated(explicit_data$case_id))
)

# ------------------------------------------------------------------
# Example 4: Parse a non-standard naming convention with a regex
# ------------------------------------------------------------------
set.seed(202604)
n_regex <- 45L
pain_baseline <- rnorm(n_regex, mean = 7, sd = 1.3)
regex_data <- data.frame(
  record_id = sprintf("R%03d", seq_len(n_regex)),
  pain__OBS_12 = pain_baseline - 3.0 + rnorm(n_regex, sd = 0.7),
  pain__OBS_0 = pain_baseline,
  pain__OBS_4 = pain_baseline - 1.4 + rnorm(n_regex, sd = 0.6)
)

# The first capture group is the outcome and the second is the time label.
# Numeric order is derived from that second group. The input columns are
# intentionally out of order so the resolved order is visible.
obs_pattern <- "^(.+)__OBS_([0-9]+)$"
regex_detected <- mira_info(
  regex_data,
  variable_pattern = obs_pattern,
  analyses = "none",
  inspect_only = TRUE,
  verbose = FALSE
)
regex_detected$detected_variables$longitudinal_map[, c(
  "outcome", "variable", "time_label", "time_order", "pattern"
)]
regex_detected$config$time_vars[["pain"]]
regex_detected$config$time_labels[["pain"]]

regex_fit <- mira_info(
  regex_data,
  variable_pattern = obs_pattern,
  outcomes = "pain",
  analyses = "none",
  verbose = FALSE
)
data.frame(
  order = seq_along(regex_fit$time_vars),
  variable = regex_fit$time_vars,
  label = unname(regex_fit$time_labels)
)
regex_fit$change[
  regex_fit$change$from == regex_fit$time_vars[[1L]] &
    regex_fit$change$to == regex_fit$time_vars[[length(regex_fit$time_vars)]],
  c("from_label", "to_label", "n", "mean_change", "paired_t_p_adj"),
  drop = FALSE
]

# ------------------------------------------------------------------
# Example 5: Treatment-arm comparisons without fitting a model
# ------------------------------------------------------------------
# Explicit arguments make the scientific comparison unambiguous. A change
# within +/-2 outcome points is classified as stable for this illustration.
arm_fit <- mira_info(
  trial,
  id = "subject_id",
  outcomes = "outcome",
  arm = "arm",
  reference_arm = "Control",
  covariates = character(0),
  analyses = "arm_tests",
  p_adjust_method = "holm",
  improvement_direction = "lower",
  stable_threshold = 2,
  verbose = FALSE
)
arm_fit$arm_analysis[c("enabled", "arm_variable", "reference_arm", "levels")]

# Values are compared between arms at every visit. The column name states the
# direction exactly: mean_difference_b_minus_a is arm_b minus arm_a.
arm_fit$arm_analysis$time_pairwise[, c(
  "time_label", "arm_a", "arm_b", "n_a", "n_b",
  "mean_difference_b_minus_a", "ci_lower", "ci_upper",
  "welch_t_p", "welch_t_p_adj"
)]

# Baseline balance is the first-timepoint subset; omnibus p-values form their
# own across-time families and therefore have their own adjusted columns.
arm_fit$arm_analysis$baseline_balance[, c(
  "arm_a", "arm_b", "mean_difference_b_minus_a", "welch_t_p_adj"
)]
arm_fit$arm_analysis$time_omnibus[, c(
  "time_label", "n", "welch_anova_p", "welch_anova_p_adj",
  "kruskal_p_adj"
)]

# Change comparisons are baseline-to-follow-up differences. Here a negative
# difference_in_change_b_minus_a means a larger reduction in Treatment.
arm_fit$arm_analysis$change_pairwise[, c(
  "to_label", "arm_a", "arm_b", "n_a", "n_b",
  "mean_change_a", "mean_change_b", "difference_in_change_b_minus_a",
  "ci_lower", "ci_upper", "welch_t_p_adj"
)]
arm_fit$arm_analysis$change_descriptives[, c(
  "arm", "to_label", "n", "mean_change",
  "improved_n", "stable_n", "worsened_n"
)]


# ------------------------------------------------------------------
# Example 6: Covariate-adjusted longitudinal model
# ------------------------------------------------------------------
# lme4 is the required engine for the primary mixed model. Other advanced
# packages are optional and record a structured skip when unavailable.
if (requireNamespace("lme4", quietly = TRUE)) {
  adjusted_fit <- mira_info(
    trial,
    id = "subject_id",
    outcomes = "outcome",
    arm = "arm",
    reference_arm = "Control",
    covariates = c("age", "site"),
    analyses = c("model", "arm_tests"),
    p_adjust_method = "holm",
    improvement_direction = "lower",
    stable_threshold = 2,
    verbose = FALSE
  )

  # Requested, retained, and skipped covariates are distinguished explicitly.
  # Covariate adjustment changes model-based estimates; it does not make the
  # treatment comparison causal.
  adjusted_fit$model[c(
    "covariates_requested", "covariates_used", "covariates_skipped",
    "converged", "singular", "warnings", "error"
  )]
  if (!is.null(adjusted_fit$model$fitted_model)) {
    stats::formula(adjusted_fit$model$fitted_model)
  }

  # These likelihood-ratio tables answer different global questions. They are
  # not combined into the contrast-level multiplicity families.
  adjusted_fit$model$global_time_test
  adjusted_fit$model$global_arm_test
  adjusted_fit$model$arm_time_interaction_test
  adjusted_fit$variability[, c("ICC_model", "ICC")]

  # ----------------------------------------------------------------
  # Example 7: Estimated marginal means and targeted contrasts
  # ----------------------------------------------------------------
  if (requireNamespace("emmeans", quietly = TRUE) &&
      isTRUE(adjusted_fit$advanced_tests$emmeans$performed)) {
    emm <- adjusted_fit$advanced_tests$emmeans

    # This named list is a compact path index. It shows where each estimand or
    # contrast family lives without printing every table.
    emm_tables <- list(
      time_means = emm$time$tidy,
      all_time_pairs = emm$time$pairwise_tidy,
      baseline_vs_followup = emm$baseline_followup$tidy,
      baseline_vs_final = emm$baseline_final$tidy,
      consecutive_time = emm$consecutive$tidy,
      arm_pairs = emm$arm$pairwise_tidy,
      arm_within_time = emm$arm_time$simple_arm_tidy,
      time_within_arm = emm$arm_time$simple_time_tidy,
      interaction_contrasts = emm$arm_time$interaction_tidy,
      ordinal_trends = emm$trends$tidy
    )
    vapply(
      emm_tables,
      function(table) if (is.null(table)) 0L else nrow(table),
      integer(1L)
    )

    # Estimated marginal means are adjusted for age and site because those
    # covariates were included in the fitted model.
    emm$time$tidy

    compact_contrasts <- function(table) {
      if (is.null(table)) return(NULL)
      keep <- intersect(
        c("contrast", "estimate", "SE", "df", "p_raw", "p_holm"),
        names(table)
      )
      table[, keep, drop = FALSE]
    }

    # Pairwise enumerates every pair of visits. Baseline-versus-final is one
    # prespecified endpoint contrast, stored separately for direct retrieval.
    compact_contrasts(emm$time$pairwise_tidy)
    compact_contrasts(emm$baseline_final$tidy)

    # Simple arm effects compare arms within each visit. The path index above
    # also exposes time-within-arm and interaction contrasts.
    compact_contrasts(emm$arm_time$simple_arm_tidy)

    # With four ordered visits, linear, quadratic, and cubic trends can be
    # returned. These use visit order, not assumed equal chronological spacing.
    compact_contrasts(emm$trends$tidy)
  }

  # ----------------------------------------------------------------
  # Example 8: Advanced modules and sensitivity analysis
  # ----------------------------------------------------------------
  # Check performed before attempting to read module-specific tables. This
  # remains safe when some optional packages are not installed or a fit fails.
  advanced_status <- c(
    rm_anova = isTRUE(adjusted_fit$advanced_tests$rm_anova$performed),
    friedman = isTRUE(adjusted_fit$advanced_tests$friedman$performed),
    gee_independence = isTRUE(
      adjusted_fit$advanced_models$gee$independence$performed
    ),
    gee_exchangeable = isTRUE(
      adjusted_fit$advanced_models$gee$exchangeable$performed
    ),
    gee_ar1 = isTRUE(adjusted_fit$advanced_models$gee$ar1$performed),
    random_slope = isTRUE(
      adjusted_fit$advanced_models$random_slope$performed
    ),
    nlme_compound_symmetry = isTRUE(
      adjusted_fit$advanced_models$nlme$compound_symmetry$performed
    ),
    nlme_ar1 = isTRUE(adjusted_fit$advanced_models$nlme$ar1$performed),
    robust_cr2 = isTRUE(
      adjusted_fit$robustness$club_sandwich$performed
    )
  )
  advanced_status

  if (advanced_status[["rm_anova"]]) {
    adjusted_fit$advanced_tests$rm_anova$tidy
  }
  adjusted_fit$advanced_tests$friedman$tidy

  if (advanced_status[["gee_exchangeable"]]) {
    adjusted_fit$advanced_models$gee$exchangeable$effect_tests
  }
  if (advanced_status[["random_slope"]]) {
    adjusted_fit$advanced_models$random_slope[c(
      "formula", "converged", "singular",
      "random_intercept_slope_correlation"
    )]
  }
  if (advanced_status[["nlme_compound_symmetry"]] ||
      advanced_status[["nlme_ar1"]]) {
    adjusted_fit$advanced_models$nlme$comparison
  }
  if (advanced_status[["robust_cr2"]]) {
    utils::head(adjusted_fit$robustness$club_sandwich$coefficient_tests)
  }

  # Likelihood criteria and GEE QIC/CIC are retained on different scales.
  model_comparison <- adjusted_fit$advanced_models$model_comparison
  if (is.data.frame(model_comparison) && nrow(model_comparison) > 0L) {
    model_comparison[, c(
      "model", "AIC", "BIC", "QIC", "CIC", "converged", "singular",
      "object_path"
    )]
  }

  # Effect-size and multiplicity objects point back to their source families.
  paired_dz <- adjusted_fit$effect_sizes$paired_cohens_dz
  if (!is.null(paired_dz)) {
    paired_dz[, c("from_label", "to_label", "n", "cohens_dz")]
  }
  friedman_family <- adjusted_fit$multiplicity$friedman_posthoc$table
  if (!is.null(friedman_family)) {
    friedman_family[, c(
      "from_label", "to_label", "p_raw", "p_holm", "p_bh"
    )]
  }
  if (isTRUE(adjusted_fit$advanced_tests$emmeans$performed)) {
    time_family <- adjusted_fit$multiplicity$time_pairwise$table
    if (!is.null(time_family)) {
      time_family[, c("contrast", "p_raw", "p_holm", "p_bh")]
    }
  }

  # Sensitivity rows align methods by scientific question and give the exact
  # source path. They may use different samples or hypotheses, so this table
  # must not be interpreted as a vote among p-values.
  adjusted_fit$sensitivity[, c(
    "question", "method", "statistic", "df", "p_raw", "effect_size",
    "effect_size_type", "n", "object_path"
  )]
}


# ------------------------------------------------------------------
# Example 9: Multiple outcomes with different improvement directions
# ------------------------------------------------------------------
set.seed(202609)
n_multi <- 54L
response_group <- rep(c("Improved", "Stable", "Worsened"), each = 18L)
symptom_baseline_multi <- rnorm(n_multi, mean = 40, sd = 5)
function_baseline_multi <- rnorm(n_multi, mean = 55, sd = 7)
symptom_delta <- unname(
  c(Improved = -6, Stable = 0, Worsened = 4)[response_group]
) + rnorm(n_multi, sd = 0.15)
function_delta <- unname(
  c(Improved = 8, Stable = 0, Worsened = -5)[response_group]
) + rnorm(n_multi, sd = 0.20)
multi_data <- data.frame(
  subject_id = sprintf("M%03d", seq_len(n_multi)),
  symptom_t0 = symptom_baseline_multi,
  symptom_t1 = symptom_baseline_multi + 0.5 * symptom_delta +
    rnorm(n_multi, sd = 0.15),
  symptom_t2 = symptom_baseline_multi + symptom_delta,
  function_t0 = function_baseline_multi,
  function_t1 = function_baseline_multi + 0.5 * function_delta +
    rnorm(n_multi, sd = 0.20),
  function_t2 = function_baseline_multi + function_delta
)

# Symptom is better when lower; function is better when higher. The named
# thresholds mean that changes within +/-1.5 symptom points or +/-2 function
# points are treated as stable. These are user-supplied practical thresholds,
# not values estimated or clinically validated by mira_info().
directions <- stats::setNames(
  c("lower", "higher"), c("symptom", "function")
)
practical_thresholds <- stats::setNames(
  c(1.5, 2), c("symptom", "function")
)
multi_fit <- mira_info(
  multi_data,
  outcomes = c("symptom", "function"),
  improvement_direction = directions,
  stable_threshold = practical_thresholds,
  analyses = "none",
  verbose = FALSE
)

# A multi-outcome result is an outcome-indexed container. Each element is a
# complete single-outcome mira_info object with the same internal paths.
names(multi_fit$outcomes)
multi_fit$config$improvement_direction
multi_fit$config$stable_threshold
multi_fit$outcomes[["symptom"]]$settings[c(
  "improvement_direction", "stable_threshold"
)]
summary(multi_fit)

# Extract the same result for every outcome instead of exploring nested lists.
baseline_to_final <- function(result) {
  result$change[
    result$change$from == result$time_vars[[1L]] &
      result$change$to == result$time_vars[[length(result$time_vars)]],
    c(
      "from_label", "to_label", "n", "mean_change",
      "improved_n", "stable_n", "worsened_n"
    ),
    drop = FALSE
  ]
}
lapply(
  multi_fit$outcomes,
  function(result) result$descriptives[, c("label", "n", "mean", "sd")]
)
lapply(multi_fit$outcomes, baseline_to_final)

# Compare response labels under a zero threshold and the practical thresholds.
# With a zero threshold, even negligible simulated deviations count as change.
multi_zero <- mira_info(
  multi_data,
  outcomes = c("symptom", "function"),
  improvement_direction = directions,
  stable_threshold = stats::setNames(c(0, 0), c("symptom", "function")),
  analyses = "none",
  verbose = FALSE
)
response_counts <- function(result) {
  lapply(
    result$outcomes,
    function(outcome_result) {
      table(outcome_result$trajectories$clinical_direction, useNA = "ifany")
    }
  )
}
list(
  zero_threshold = response_counts(multi_zero),
  practical_threshold = response_counts(multi_fit)
)


# The multi-outcome plot method needs outcome = when several outcomes succeed.
if (requireNamespace("ggplot2", quietly = TRUE)) {
  multi_plot_fit <- mira_info(
    multi_data,
    outcomes = c("symptom", "function"),
    improvement_direction = directions,
    stable_threshold = practical_thresholds,
    analyses = "plots",
    verbose = FALSE
  )
  plot(multi_plot_fit, outcome = "symptom", which = "mean_ci")
}


# ------------------------------------------------------------------
# Example 10: Missing follow-up data and effective sample sizes
# ------------------------------------------------------------------
missing_data <- wide
missing_data$score_t1[c(4, 9, 15, 22)] <- NA_real_
missing_data$score_t2[c(2, 4, 11, 19, 28, 37)] <- NA_real_
missing_fit <- mira_info(
  missing_data,
  analyses = "correlations",
  verbose = FALSE
)

# mira_info() does not impute. Descriptives use all finite observations at
# each visit, while complete_profiles requires all selected visits.
missing_fit$overview[c(
  "n_patients", "complete_profiles", "complete_profiles_pct"
)]
missing_fit$descriptives[, c(
  "label", "n", "missing", "non_finite", "unavailable"
)]
missing_fit$missing$by_time

# Identify affected subjects without printing all rows.
affected <- missing_fit$missing$by_patient$unavailable_n > 0L
missing_fit$missing$by_patient[affected, , drop = FALSE]

# Each paired comparison has its own pairwise-complete N. Friedman instead
# uses complete profiles across all visits, reported in its tidy row.
missing_fit$change[, c(
  "from_label", "to_label", "n", "mean_change", "paired_t_p_adj"
)]
missing_fit$advanced_tests$friedman$tidy[, c("n", "n_timepoints", "p_raw")]

# Per-outcome adaptation records show usable counts and any requested module
# disabled because of insufficient data.
missing_fit$diagnostics$adaptation[["score"]]

# ------------------------------------------------------------------
# Example 11: Correlations, outlier diagnostics, and public plot methods
# ------------------------------------------------------------------
diagnostic_data <- wide
diagnostic_data$score_t2[[1L]] <- diagnostic_data$score_t2[[1L]] + 35
diagnostic_fit <- mira_info(
  diagnostic_data,
  analyses = c("correlations", "outliers", "plots"),
  verbose = FALSE
)

diagnostic_fit$correlations$pearson
diagnostic_fit$correlations$pairwise_n

# Outliers are Tukey 1.5-IQR flags only; they are never removed automatically.
vapply(diagnostic_fit$outliers$by_time, nrow, integer(1L))
diagnostic_fit$outliers$by_time[["score_t2"]]
change_flag_counts <- vapply(
  diagnostic_fit$outliers$change, nrow, integer(1L)
)
change_flag_counts
most_flagged_change <- names(change_flag_counts)[which.max(change_flag_counts)]
diagnostic_fit$outliers$change[[most_flagged_change]]

# Use plot() rather than reaching into $plots when displaying a stored plot.
# If ggplot2 is absent, the analysis still succeeds and $plot_error explains
# why no plot was stored.
diagnostic_fit$plot_error
if (requireNamespace("ggplot2", quietly = TRUE)) {
  names(diagnostic_fit$plots)
  plot(diagnostic_fit, which = "mean_ci")
  plot(diagnostic_fit, which = "spaghetti")
  plot(diagnostic_fit, which = "missingness")
}

# ------------------------------------------------------------------
# Example 12: Select exactly which optional modules to request
# ------------------------------------------------------------------
# inspect_only makes this comparison cheap: no statistical module is run.
cfg_none <- mira_info(
  wide, analyses = "none", inspect_only = TRUE, verbose = FALSE
)
cfg_arm <- mira_info(
  trial, analyses = "arm_tests", inspect_only = TRUE, verbose = FALSE
)
cfg_model_cor <- mira_info(
  wide,
  analyses = c("model", "correlations"),
  inspect_only = TRUE,
  verbose = FALSE
)
cfg_all <- mira_info(
  trial, analyses = "all", inspect_only = TRUE, verbose = FALSE
)

module_switches <- function(configuration) {
  unlist(
    configuration$config$analyses[c(
      "plots", "model", "outliers", "correlations", "arm_tests"
    )],
    use.names = TRUE
  )
}
rbind(
  none = module_switches(cfg_none),
  arm_tests = module_switches(cfg_arm),
  model_correlations = module_switches(cfg_model_cor),
  all = module_switches(cfg_all)
)

# Core summaries remain enabled under analyses = "none". On a real analysis,
# the base-R Friedman branch is also attempted independently of model =.
cfg_none$config$analyses$core





