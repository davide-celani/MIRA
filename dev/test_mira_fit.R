# ============================================================
# MIRA - Development test script
# ============================================================

# Questo script serve per testare manualmente
# il workflow principale di MIRA durante lo sviluppo.
#
# Workflow:
#
# simulated data
#      ↓
# mira_prepare_data()
#      ↓
# Stan data
#      ↓
# mira_fit()
#      ↓
# posterior
#      ↓
# mira_summary()
#
# NON è un testthat test.
# È uno script di sviluppo.
# ============================================================

devtools::document()




# ------------------------------------------------------------
# 0. Carica la versione locale di MIRA
# ------------------------------------------------------------

warnings()
getwd()


devtools::load_all()

list_stan_models("MIRA")


# Usa "BCVA" oppure "CMT"


# Ready-to-use objects. Source this file and start from `bcva_data` or
# `cmt_data`; `data` contains both outcomes for descriptive work.
data <- mira_create_ophthalmology_data()
bcva_data <- mira_select_ophthalmology_outcome(data, "BCVA")
cmt_data <- mira_select_ophthalmology_outcome(data, "CMT")


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


# ============================================================
# EXPLORATORY ANALYSIS
# ============================================================

# Solo timepoint
resul <- mira_info(
  data,
  time_vars = c("BCVA_t0", "BCVA_t1", "BCVA_t2")
)

# Timepoint + age
result <- mira_info(
  data,
  time_vars = c("BCVA_t0", "BCVA_t1", "BCVA_t2"),
  covariates = "age"
)

# Timepoint + gender
result <- mira_info(
  data,
  time_vars = c("BCVA_t0", "BCVA_t1", "BCVA_t2"),
  covariates = "gender"
)

# Timepoint + age + gender
result <- mira_info(
  data,
  time_vars = c("BCVA_t0", "BCVA_t1", "BCVA_t2"),
  covariates = c("age", "gender")
)

# ID + timepoint
result <- mira_info(
  data,
  id = "patient",
  time_vars = c("BCVA_t0", "BCVA_t1", "BCVA_t2")
)

# ID + timepoint + age + gender
result <- mira_info(
  data,
  id = "patient",
  time_vars = c("BCVA_t0", "BCVA_t1", "BCVA_t2"),
  covariates = c("age", "gender")
)

# Timepoint + trattamento
result <- mira_info(
  data,
  time_vars = c("BCVA_t0", "BCVA_t1", "BCVA_t2"),
  arm = "treatment"
)

# ID + timepoint + trattamento + covariate
result <- mira_info(
  data,
  id = "patient",
  time_vars = c("BCVA_t0", "BCVA_t1", "BCVA_t2", "BCVA_t3","BCVA_t4"),
  arm = "treatment",
  reference_arm = "Aflibercept",
  covariates = c("age", "gender","study_eye")
)

# Un solo outcome selezionato, con rilevamento automatico dei timepoint
result <- mira_info(
  data,
  outcomes = "BCVA"
)

# Più outcome rilevati automaticamente
result <- mira_info(
  data,
  outcomes = c("BCVA", "CMT", "IOP")
)

# Più outcome con colonne specificate manualmente
result <- mira_info(
  data,
  time_vars = list(
    BCVA = c("BCVA_t0", "BCVA_t1", "BCVA_t2"),
    CMT  = c("CMT_t0", "CMT_t1", "CMT_t2")
  )
)

# Timepoint con etichette personalizzate
result <- mira_info(
  data,
  time_vars = c("BCVA_t0", "BCVA_t1", "BCVA_t2"),
  time_labels = c("Baseline", "Month 1", "Month 3")
)

# Nomi longitudinali non standard
result <- mira_info(
  data,
  outcomes = "BCVA",
  time_vars = c("BCVA_baseline", "BCVA_month1", "BCVA_month3"),
  time_labels = c("Baseline", "Month 1", "Month 3")
)

# Direzione del miglioramento e soglia di stabilità
result <- mira_info(
  data,
  time_vars = c("BCVA_t0", "BCVA_t1", "BCVA_t2"),
  improvement_direction = "higher",
  stable_threshold = 5
)

# Più outcome con direzioni diverse
result <- mira_info(
  data,
  outcomes = c("BCVA", "CMT"),
  improvement_direction = c(
    BCVA = "higher",
    CMT = "lower"
  ),
  stable_threshold = c(
    BCVA = 5,
    CMT = 20
  )
)

# Solo analisi selezionate
result <- mira_info(
  data,
  time_vars = c("BCVA_t0", "BCVA_t1", "BCVA_t2"),
  analyses = c("model", "correlations", "outliers")
)

# Tutte le analisi opzionali disabilitate
result <- mira_info(
  data,
  time_vars = c("BCVA_t0", "BCVA_t1", "BCVA_t2"),
  analyses = "none"
)

# Solo ispezione della configurazione
config <- mira_info(
  data,
  inspect_only = TRUE
)

# Chiamata completamente automatica
result <- mira_info(data)

res$descriptives
res$change
res$variability
res$model
res$plots
res$plots$arm_boxplot


# ============================================================
# PREPARE DATA
# ============================================================

stan_data <- mira_prepare_data(
  data = analysis_data,
  time_value = c(0, 3, 5, 12, 15),
  meaningful_change = meaningful_change,
  meaningful_change_sd = meaningful_change_sd,
  direction = direction,
  reference_arm = "control",
  age_threshold = 60
)


# ============================================================
# PRIORS
# ============================================================

prior <- mira_prior(stan_data, outcome = "BCVA")

print(prior)


# ============================================================
# TEST FIT
# ============================================================

fit <- mira_fit(
  stan_data = stan_data,
  prior = prior,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 2000,
  iter_sampling = 3000,
  seed = 123,
  refresh = 100,
  verbose = TRUE
)

fit$diagnostic_summary()
fit$cmdstan_diagnose()


# ============================================================
# SUMMARY
# ============================================================

mira_res <- mira_summary(
  fit = fit,
  stan_data = stan_data,
  verbose = TRUE
)

mira_res$population_time_means
mira_res$change
mira_res$diagnostics

# ------------------------------------------------------------
# 6. Verifica il fit
# ------------------------------------------------------------


print(class(fit))






# ------------------------------------------------------------
# 8. Controllo semplice dei parametri
# ------------------------------------------------------------

cat("\n")
cat("============================================\n")
cat("EXPECTED VALUES\n")
cat("============================================\n\n")

cat("True T0:", true_t0, "\n")
cat("True T1:", true_t1, "\n")
cat("True T2:", true_t2, "\n")

cat("\nTrue change T1 - T0:",
    true_t1 - true_t0,
    "\n")

cat("True change T2 - T0:",
    true_t2 - true_t0,
    "\n")

cat("True change T2 - T1:",
    true_t2 - true_t1,
    "\n")


# ------------------------------------------------------------
# 9. Fine
# ------------------------------------------------------------

cat("\n")
cat("============================================\n")
cat("MIRA DEVELOPMENT TEST COMPLETED\n")
cat("============================================\n\n")


print("first change done on git")
