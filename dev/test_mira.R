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

res <- mira_info(
  data = analysis_data,
  id = "patient",
  time_vars = paste0(outcome, "_t", 0:4),
  time_labels = c(
    "Baseline",
    "Month 3",
    "Month 5",
    "Month 12",
    "Month 15"
  ),
  arm = "arm",
  reference_arm = "control",
  improvement_direction = direction,
  plots = TRUE,
  model = TRUE,
  outliers = TRUE,
  correlations = TRUE,
  verbose = TRUE,
  p_adjust_method = "holm",
  stable_threshold = 0,
  strict_id = TRUE,
  arm_tests = TRUE
)

res$descriptives
res$change
res$variability
res$model
res$plots


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

prior <- mira_prior(
  stan_data = stan_data,
  outcome = outcome,
  informativeness = "standard"
)

print(prior)


# ============================================================
# TEST FIT
# ============================================================

fit <- mira_fit(
  stan_data = stan_data,
  prior = prior,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1500,
  iter_sampling = 2000,
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
mira_res$treatment_effects
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
