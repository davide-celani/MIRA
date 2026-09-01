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


# ============================================================
# MIRA INFO
# ============================================================

res <- mira_info(
  data = data,
  id = "patient",
  time_vars = c(
    "IOP_t0",
    "IOP_t1",
    "IOP_t2",
    "IOP_t3",
    "IOP_t4"
  ),
  time_labels = c(
    "Baseline",
    "Month 3",
    "Month 5",
    "Month 12",
    "Month 15"
  ),
  arm = "arm",
  reference_arm = "control",
  improvement_direction = "lower",
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

# Alcuni risultati esplorativi.
res$descriptives
res$change
res$variability
res$model

# ============================================================
# PREPARE DATA
# ============================================================

stan_data <- mira_prepare_data(
  data = data,
  time_value = c(0, 3, 5, 12, 15),
  meaningful_change = 2,
  meaningful_change_sd = 0.5,
  direction = "lower",
  reference_arm = "control",
  age_threshold = 60
)

# ============================================================
# PRIORS
# ============================================================

prior <- mira_prior(
  stan_data = stan_data,
  profile = "default"
)

print(prior)

# ============================================================
# TEST FIT
# ============================================================
# Poche iterazioni: serve solo a verificare che il modello compili e campioni.

fit <- mira_fit(
  stan_data = stan_data,
  prior = prior,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  seed = 123,
  refresh = 50,
  verbose = TRUE
)

# Diagnostica tecnica del test.
fit$diagnostic_summary()

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
