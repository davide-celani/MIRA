

devtools::load_all()



# Ready-to-use objects. Source this file and start from `bcva_data` or
# `cmt_data`; `data` contains both outcomes for descriptive work.
data <- mira_create_ophthalmology_data()


outcome <- "BCVA"

# Seleziona automaticamente le colonne dell'outcome scelto
analysis_data <- mira_select_ophthalmology_outcome(
  data = data,
  outcome = outcome
)

# Impostazioni specifiche dell'outcome
direction <- "higher" # "higher" else "lower"
meaningful_change <- 5 # 5 else 50
meaningful_change_sd <- 1.5 # 1.5 else 15







result <- mira_info(
  data = data,
  time_vars = c("BCVA_t0", "BCVA_t1", "BCVA_t2")
)

report <- mira_report_freq(
  result,
  author = "Davide Celani",
  output_dir = "report_mira",
  output_file = "analisi_longitudinale",
  format = "pdf",
  quiet = FALSE,
  overwrite = TRUE
)












# ============================================================
# 1. ANALISI DI BASE
# ============================================================

# Solo timepoint
result <- mira_info(
  data = data,
  time_vars = c("BCVA_t0", "BCVA_t1", "BCVA_t2")
)

# Timepoint + età
result <- mira_info(
  data = data,
  time_vars = c("BCVA_t0", "BCVA_t1", "BCVA_t2"),
  covariates = "age"
)

# Timepoint + genere
result <- mira_info(
  data = data,
  time_vars = c("BCVA_t0", "BCVA_t1", "BCVA_t2"),
  covariates = "gender"
)

# Timepoint + età + genere
result <- mira_info(
  data = data,
  time_vars = c("BCVA_t0", "BCVA_t1", "BCVA_t2"),
  covariates = c("age", "gender")
)

# ID specificato manualmente
result <- mira_info(
  data = data,
  id = "patient",
  time_vars = c("BCVA_t0", "BCVA_t1", "BCVA_t2")
)

# ID + covariate
result <- mira_info(
  data = data,
  id = "patient",
  time_vars = c("BCVA_t0", "BCVA_t1", "BCVA_t2"),
  covariates = c("age", "gender")
)


# ============================================================
# 2. TRATTAMENTO E CONFRONTO TRA ARM
# ============================================================

# Timepoint + trattamento
result <- mira_info(
  data = data,
  time_vars = c("BCVA_t0", "BCVA_t1", "BCVA_t2"),
  arm = "treatment"
)

# Analisi completa con trattamento e covariate
result <- mira_info(
  data = data,
  id = "patient",
  time_vars = c(
    "BCVA_t0", "BCVA_t1", "BCVA_t2",
    "BCVA_t3", "BCVA_t4"
  ),
  arm = "treatment",
  reference_arm = "Aflibercept",
  covariates = c("age", "gender", "study_eye"),
  p_adjust_method = "holm",
  verbose = TRUE
)

# Stampa il report completo dopo l'elaborazione
print(result)

# Versione compatta
summary(result)


# ============================================================
# 3. OUTCOME AUTOMATICI E MULTI-OUTCOME
# ============================================================

# Un outcome con rilevamento automatico dei timepoint
result <- mira_info(
  data = data,
  outcomes = "BCVA"
)

# Più outcome rilevati automaticamente
result <- mira_info(
  data = data,
  outcomes = c("BCVA", "CMT", "IOP")
)

# Più outcome con colonne definite manualmente
result <- mira_info(
  data = data,
  time_vars = list(
    BCVA = c("BCVA_t0", "BCVA_t1", "BCVA_t2"),
    CMT  = c("CMT_t0",  "CMT_t1",  "CMT_t2"),
    IOP  = c("IOP_t0",  "IOP_t1",  "IOP_t2")
  )
)

# Multi-outcome con trattamento e covariate
result <- mira_info(
  data = data,
  id = "patient",
  outcomes = c("BCVA", "CMT"),
  arm = "treatment",
  reference_arm = "Aflibercept",
  covariates = c("age", "gender")
)


# ============================================================
# 4. TIME LABELS E NOMI NON STANDARD
# ============================================================

# Etichette temporali personalizzate
result <- mira_info(
  data = data,
  time_vars = c("BCVA_t0", "BCVA_t1", "BCVA_t2"),
  time_labels = c("Baseline", "Month 1", "Month 3")
)

# Nomi longitudinali non standard specificati manualmente
result <- mira_info(
  data = data,
  outcomes = "BCVA",
  time_vars = c(
    "BCVA_baseline",
    "BCVA_month1",
    "BCVA_month3"
  ),
  time_labels = c("Baseline", "Month 1", "Month 3")
)

# Regex personalizzata:
# gruppo 1 = nome outcome, gruppo 2 = timepoint
result <- mira_info(
  data = data,
  variable_pattern = "^(.+)_visit_([0-9]+)$"
)


# ============================================================
# 5. DIREZIONE CLINICA
# ============================================================

# BCVA: valori maggiori rappresentano miglioramento
result <- mira_info(
  data = data,
  time_vars = c("BCVA_t0", "BCVA_t1", "BCVA_t2"),
  improvement_direction = "higher",
  stable_threshold = 5
)

# CMT: valori minori rappresentano miglioramento
result <- mira_info(
  data = data,
  time_vars = c("CMT_t0", "CMT_t1", "CMT_t2"),
  improvement_direction = "lower",
  stable_threshold = 20
)

# Più outcome con direzioni e soglie differenti
result <- mira_info(
  data = data,
  outcomes = c("BCVA", "CMT"),
  improvement_direction = c(
    BCVA = "higher",
    CMT  = "lower"
  ),
  stable_threshold = c(
    BCVA = 5,
    CMT  = 20
  )
)


# ============================================================
# 6. SELEZIONE DELLE ANALISI
# ============================================================

# Modelli, correlazioni e outlier
result <- mira_info(
  data = data,
  time_vars = c("BCVA_t0", "BCVA_t1", "BCVA_t2"),
  analyses = c("model", "correlations", "outliers")
)

# Modelli longitudinali con trattamento
# arm_tests deve essere incluso per TIME × ARM
result <- mira_info(
  data = data,
  time_vars = c("BCVA_t0", "BCVA_t1", "BCVA_t2"),
  arm = "treatment",
  analyses = c("model", "arm_tests")
)

# Disabilita le analisi opzionali.
# Descrittive, cambiamenti e Friedman restano disponibili.
result <- mira_info(
  data = data,
  time_vars = c("BCVA_t0", "BCVA_t1", "BCVA_t2"),
  analyses = "none"
)

# Ispezione senza eseguire le analisi
config <- mira_info(
  data = data,
  inspect_only = TRUE
)

# Configurazione completamente automatica
result <- mira_info(data)





result$descriptives
result$change
result$missing
result$correlations
result$variability
result$trajectories
result$model
result$arm_analysis
result$plots
result$long_data

# Modello mixed originale
summary(result$model$fitted_model)
result$model$global_time_test
result$model$global_arm_test
result$model$arm_time_interaction_test

# Plot
result$plots$boxplot
result$plots$mean_ci
result$plots$arm_boxplot
result$plots$arm_change_ci

plot(result, which = "boxplot")
plot(result, which = "arm_boxplot")





# Stato del modulo
result$advanced_tests$emmeans$performed
result$advanced_tests$emmeans$error
result$advanced_tests$emmeans$reason_skipped

# Estimated marginal means del tempo
result$advanced_tests$emmeans$time$tidy
summary(result$advanced_tests$emmeans$time$object)

# Tutti i confronti tra timepoint
result$advanced_tests$emmeans$time$pairwise_tidy
summary(result$advanced_tests$emmeans$time$pairwise)

# Confronti Tukey
result$advanced_tests$emmeans$time$pairwise_tukey

# Baseline contro ogni follow-up
result$advanced_tests$emmeans$baseline_followup$tidy
summary(result$advanced_tests$emmeans$baseline_followup$object)

# Correzione Dunnett
result$advanced_tests$emmeans$baseline_followup$dunnett

# Baseline contro ultimo follow-up
result$advanced_tests$emmeans$baseline_final$tidy
summary(result$advanced_tests$emmeans$baseline_final$object)

# Confronti consecutivi
result$advanced_tests$emmeans$consecutive$tidy

# Trend lineare, quadratico e cubico
result$advanced_tests$emmeans$trends$tidy
summary(result$advanced_tests$emmeans$trends$object)





# Marginal means per arm
result$advanced_tests$emmeans$arm$tidy
summary(result$advanced_tests$emmeans$arm$object)

# Confronti globali tra arm
result$advanced_tests$emmeans$arm$pairwise_tidy

# Marginal means ARM × TIME
result$advanced_tests$emmeans$arm_time$tidy
summary(result$advanced_tests$emmeans$arm_time$object)

# Effetto dell'arm in ogni timepoint
result$advanced_tests$emmeans$arm_time$simple_arm_tidy
summary(result$advanced_tests$emmeans$arm_time$simple_arm)

# Effetto del tempo separatamente in ogni arm
result$advanced_tests$emmeans$arm_time$simple_time_tidy
summary(result$advanced_tests$emmeans$arm_time$simple_time)

# Contrasti dell'interazione TIME × ARM
result$advanced_tests$emmeans$arm_time$interaction_tidy
summary(result$advanced_tests$emmeans$arm_time$interaction)




rm_anova <- result$advanced_tests$rm_anova

rm_anova$performed
rm_anova$error
rm_anova$warnings

# Tabella compatta:
# F, df, p non corretto, GG, HF, eta squared
rm_anova$tidy

# Oggetto originale
rm_anova$object
summary(rm_anova$object)

# Tabelle specifiche
rm_anova$tables$uncorrected_partial_eta
rm_anova$tables$uncorrected_generalized_eta
rm_anova$tables$greenhouse_geisser
rm_anova$tables$huynh_feldt

# Mauchly e informazioni sulla sfericità
rm_anova$sphericity$mauchly
rm_anova$sphericity$corrections




friedman <- result$advanced_tests$friedman

# Test globale originale
friedman$test

# Statistica, p-value e Kendall's W
friedman$tidy

# Tutti i Wilcoxon paired post-hoc
friedman$posthoc$tidy

# Oggetti htest originali
friedman$posthoc$tests

# Un confronto specifico
friedman$posthoc$tests$BCVA_t0_to_BCVA_t1



gee <- result$advanced_models$gee

# Stato dei tre modelli
gee$independence$performed
gee$exchangeable$performed
gee$ar1$performed

# Modello con correlazione exchangeable
summary(gee$exchangeable$model)

# Coefficienti, SE robusti e intervalli di confidenza
gee$exchangeable$coefficients
gee$exchangeable$robust_standard_errors

# Test Wald globali di TIME, ARM e TIME × ARM
gee$exchangeable$effect_tests

# QIC e CIC
gee$exchangeable$qic

# Confronto tra strutture
gee$independence$qic
gee$exchangeable$qic
gee$ar1$qic





random_slope <- result$advanced_models$random_slope

random_slope$performed
random_slope$error
random_slope$warnings

# Oggetto originale
random_slope$model
summary(random_slope$model)

# Risultati estratti
random_slope$fixed_effects
random_slope$fixed_effect_tests
random_slope$variance_components
random_slope$random_intercept_slope_correlation

# AIC, BIC e log-likelihood
random_slope$AIC
random_slope$BIC
random_slope$logLik

# Convergenza
random_slope$converged
random_slope$singular

# Random intercept contro random slope
random_slope$random_effects_lrt

# Test globali
random_slope$global_time_test
random_slope$global_arm_test
random_slope$arm_time_interaction_test




nlme_models <- result$advanced_models$nlme

# Compound symmetry
summary(nlme_models$compound_symmetry$model)
nlme_models$compound_symmetry$coefficients
nlme_models$compound_symmetry$correlation_parameter

# AR(1)
summary(nlme_models$ar1$model)
nlme_models$ar1$coefficients
nlme_models$ar1$correlation_parameter

# Tabella comparativa CS vs AR(1)
nlme_models$comparison






robust <- result$robustness$club_sandwich

robust$performed
robust$error
robust$warnings

# Matrice di covarianza CR2
robust$covariance

# Test dei singoli coefficienti
robust$coefficient_tests

# Intervalli di confidenza robusti
robust$confidence_intervals

# Test globali robusti
robust$effect_tests$TIME
robust$effect_tests$ARM
robust$effect_tests$TIME_X_ARM


# Cohen's dz per confronti paired
result$effect_sizes$paired_cohens_dz

# Rank-biserial correlation dei Wilcoxon paired
result$effect_sizes$paired_rank_biserial

# Partial e generalized eta squared
result$effect_sizes$rm_anova

# Kendall's W
result$effect_sizes$friedman$kendalls_w

# R² marginale e condizionale del random-intercept model
result$effect_sizes$mixed_models$random_intercept$tidy

# R² marginale e condizionale del random-slope model
result$effect_sizes$mixed_models$random_slope$tidy

# ICC del mixed model
result$effect_sizes$mixed_models$ICC




# Il metodo scelto dall'utente resta quello primario
result$multiplicity$primary_method

# Confronti tra timepoint
time_p <- result$multiplicity$time_pairwise$table

time_p[, c(
  "contrast",
  "estimate",
  "p_raw",
  "p_primary",
  "p_bonferroni",
  "p_holm",
  "p_bh",
  "p_by"
)]

# Baseline contro follow-up
result$multiplicity$baseline_vs_followup$table

# Confronti consecutivi
result$multiplicity$consecutive_time$table

# Trend temporali
result$multiplicity$ordinal_trends$table

# Confronti tra arm
result$multiplicity$arm_pairwise$table

# Simple effects
result$multiplicity$simple_effects$arm_within_time$table
result$multiplicity$simple_effects$time_within_arm$table

# Contrasti di interazione
result$multiplicity$interaction_contrasts$table

# Wilcoxon post-hoc
result$multiplicity$friedman_posthoc$table





# Confronto compatto di tutti i modelli longitudinali
result$advanced_models$model_comparison

# Tutte le inferenze confrontabili
result$sensitivity

# Solo effetto globale del tempo
subset(result$sensitivity, question == "TIME")

# Solo effetto globale dell'arm
subset(result$sensitivity, question == "ARM")

# Solo interazione TIME × ARM
subset(result$sensitivity, question == "TIME_X_ARM")




result <- mira_info(
  data = data,
  outcomes = c("BCVA", "CMT"),
  arm = "treatment",
  covariates = c("age", "gender"),
  verbose = FALSE
)

# Estrazione di un outcome
bcva <- result$outcomes[["BCVA"]]
cmt  <- result$outcomes[["CMT"]]

# Analisi avanzate per BCVA
bcva$advanced_tests$rm_anova$tidy
bcva$advanced_models$gee$exchangeable$effect_tests
bcva$sensitivity

# Analisi avanzate per CMT
cmt$advanced_tests$friedman$tidy
cmt$advanced_models$model_comparison
cmt$sensitivity

# Sensitivity table per tutti gli outcome
sensitivity_by_outcome <- lapply(
  result$outcomes,
  function(x) {
    if (inherits(x, "mira_info_error")) {
      return(data.frame(error = x$error))
    }
    x$sensitivity
  }
)




result$advanced_tests$rm_anova[
  c("performed", "package", "error", "reason_skipped", "warnings")
]

result$advanced_models$gee$exchangeable[
  c("performed", "package", "error", "reason_skipped", "warnings")
]

result$robustness$club_sandwich[
  c("performed", "package", "error", "reason_skipped", "warnings")
]
