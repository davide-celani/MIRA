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


# ------------------------------------------------------------
# 1. Impostazioni
# ------------------------------------------------------------

set.seed(123)

n_patients <- 20


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

data$t0 <- rnorm(
  n_patients,
  mean = true_t0 + subject_effect,
  sd = true_sigma
)

data$t1 <- rnorm(
  n_patients,
  mean = true_t1 + subject_effect,
  sd = true_sigma
)

data$t2 <- rnorm(
  n_patients,
  mean = true_t2 + subject_effect,
  sd = true_sigma
)


# ------------------------------------------------------------
# Prepare data for Stan
# ------------------------------------------------------------

stan_data <- mira_prepare_data(
  data = data,
  time_value = c(0, 6, 12),
  meaningful_change = 1.5
)

# ------------------------------------------------------------
# Fit model
# ------------------------------------------------------------


prior <- mira_prior(
  stan_data,
  profile = "custom",
  mu_time_sd = 5 * stan_data$sd_y,
  sigma_intercept_rate = 1 / stan_data$sd_y,
  sigma_slope_rate = 1 / stan_data$sd_y,
  sigma_rate = 1 / stan_data$sd_y,
  nu_shape = 2,
  nu_rate = 0.1
)

prior <- mira_prior(
  stan_data,
  profile = "regularized"
)

print(prior)

fit <- mira_fit(
  stan_data,
  prior=prior
)


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
