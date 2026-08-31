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
# Parametri
# ------------------------------------------------------------

true_t0 <- 10
true_t1 <- 11
true_t2 <- 12

true_sigma_subject <- 1
true_sigma <- 0.5


# ------------------------------------------------------------
# Effetto casuale del paziente
# ------------------------------------------------------------

subject_effect <- rnorm(
  n_patients,
  mean = 0,
  sd = true_sigma_subject
)


# ------------------------------------------------------------
# Dataset
# ------------------------------------------------------------

data <- data.frame(

  patient = seq_len(n_patients),

  arm = rep(
    c("control", "treatment"),
    length.out = n_patients
  ),

  gender = sample(
    c("Female", "Male"),
    size = n_patients,
    replace = TRUE,
    prob = c(0.5, 0.5)
  ),

  age = round(
    rnorm(
      n_patients,
      mean = 60,
      sd = 10
    )
  )
)


# Limitiamo l'età a un range plausibile
data$age <- pmin(
  pmax(data$age, 30),
  85
)


# ------------------------------------------------------------
# Misure IOP
# ------------------------------------------------------------

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
  reference_arm = "control",
  improvement_direction = "lower"
)


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




res$descriptives
res$change
res$variability
res$trajectories
res$long_data
res$call
res$overview
res$missing
res$correlations
res$plots$boxplot
res$plots$spaghetti
res$outliers
res$model
# ============================================================
# PREPARE DATA
# ============================================================

stan_data <- mira_prepare_data(
  data = data,
  time_value = c(0, 1, 2, 3, 4),
  meaningful_change = 2,
  meaningful_change_sd = 0.5,
  direction = "lower",
  reference_arm = "control",
  age_threshold = 60
)

prior <- mira_prior(
  stan_data,
  profile = "default"
)

fit <- mira_fit(
  stan_data = stan_data,
  prior = prior
)


fit <- mira_fit(
  stan_data = stan_data,
  prior = prior,
  chains = 2,
  parallel_chains = 2,
  iter_warmup = 500,
  iter_sampling = 500,
  seed = 123,
  refresh = 100
)



summary <- mira_summary(
  fit,
  stan_data = stan_data
)

summary$population_time_means
summary$change_from_baseline
summary$gender$change


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
