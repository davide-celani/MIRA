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
# ------------------------------------------------------------
# 1. Impostazioni
# ------------------------------------------------------------

set.seed(123)

n_patients <- 150


# ------------------------------------------------------------
# 2. Crea dataset simulato
# ------------------------------------------------------------




true_t0 <- 10
true_t1 <- 11
true_t2 <- 12

true_sigma_subject <- 1
true_sigma <- 0.5


# Effetto casuale del paziente

subject_effect <- rnorm(
  n_patients,
  mean = 0,
  sd = true_sigma_subject
)


# Genera le tre misure per ogni paziente

data <- data.frame(
  patient = seq_len(n_patients)
)

data$IOP_t0 <- rnorm(
  n_patients,
  mean = true_t0 + subject_effect,
  sd = true_sigma
)

data$IOP_t1 <- rnorm(
  n_patients,
  mean = true_t1 + subject_effect,
  sd = true_sigma
)

data$IOP_t2 <- rnorm(
  n_patients,
  mean = true_t2 + subject_effect,
  sd = true_sigma
)


data$IOP_t3 <- rnorm(
  n_patients,
  mean = true_t2 + subject_effect,
  sd = true_sigma
)

data$IOP_t4 <- rnorm(
  n_patients,
  mean = true_t1 + subject_effect,
  sd = true_sigma
)



res <- mira_info(
  data = data,
  id = "patient",
  time_vars = c("IOP_t0", "IOP_t1", "IOP_t2", "IOP_t3", "IOP_t4"),
  time_labels = c(
    "Baseline",
    "Time 1",
    "Time 2",
    "Time 3",
    "Time 4"
  ),
  plots = TRUE,
  model = TRUE,
  outliers = TRUE,
  correlations = TRUE,
  verbose = TRUE
)





res$descriptives
res$change
res$variability
res$trajectories
res$long_data
res$call
res$overview
res$missing
res$correlations
res$plots
res$plots$spaghetti
res$outliers
res$model
# ============================================================
# PREPARE DATA
# ============================================================

stan_data <- mira_prepare_data(
  data = data,
  time_value = c(0, 3, 5, 12, 15),
  meaningful_change = 1
)


# ============================================================
# PRIOR
# ============================================================


prior <- mira_prior(
  stan_data,
  profile = "default"
)

fit <- mira_fit(
  stan_data = stan_data,
  prior = prior,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 300,
  iter_sampling = 1000,
  seed = 123
)

# ============================================================
# STAN MODEL
# ============================================================

stan_file <- system.file(
  "stan",
  "gaussian_longitudinal.stan",
  package = "MIRA"
)




# ============================================================
# PRINT RESULTS
# ============================================================

summary <- mira_summary(
  fit = fit,
  stan_data = stan_data
)


summary$clinical
summary$responders



# ------------------------------------------------------------
# 6. Verifica il fit
# ------------------------------------------------------------


print(class(fit))


# ------------------------------------------------------------
# 7. Posterior summary
# ------------------------------------------------------------

summary_mira <- mira_summary(
  fit = fit,
  meaningful_change = stan_data$meaningful_change,
  y = stan_data$y
)



summary_mira$population
summary_mira$clinical
summary_mira$individual
summary_mira$responders
summary_mira$heterogeneity
summary_mira$draws
summary_mira$ppc
summary_mira$loo
summary_mira$diagnostics
summary_mira$quality_flags
summary_mira$model_information

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
