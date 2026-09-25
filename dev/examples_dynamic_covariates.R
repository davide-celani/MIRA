# ============================================================
# MIRA - Esempi diretti: dal modello semplice a quello completo
# ============================================================
#
# Eseguire questo script dalla root del package MIRA:
#
#   source("dev/examples_dynamic_covariates.R")
#
# Ogni scenario usa direttamente il flusso pubblico del package:
#
#   mira_data_long()
#   mira_prior_long()
#   mira_fit_long()
#   diagnostica CmdStan
#   mira_summary_long()
#
# Non ci sono funzioni wrapper. Lo script esegue otto fit completi; per
# provare un solo scenario, eseguire interattivamente il relativo blocco.
# ============================================================

if (!file.exists("DESCRIPTION") || !file.exists("dev/data_creation.R")) {
  stop(
    "Eseguire `dev/examples_dynamic_covariates.R` dalla root del package MIRA.",
    call. = FALSE
  )
}

if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("Il package `devtools` è necessario per caricare MIRA localmente.",
       call. = FALSE)
}

devtools::load_all(quiet = TRUE)
source("dev/data_creation.R", local = FALSE)

data <- mira_create_ophthalmology_data(
  n_per_arm = 10,
  seed = 2026,
  time_months = c(0, 3, 5, 12, 15)
)

bcva_visits <- paste0("BCVA_t", 0:4)
cmt_visits <- paste0("CMT_t", 0:4)


# ============================================================
# BCVA 1 - Solo tempo e trattamento, nessuna covariata
# ============================================================

data_bcva_01 <- data[c(
  "patient",
  "arm",
  bcva_visits
)]

stan_data_bcva_01 <- mira_data_long(
  data = data_bcva_01,
  time_value = c(0, 3, 5, 12, 15),
  outcome = "BCVA",
  meaningful_change = 5,
  meaningful_change_sd = 1.5,
  direction = "higher",
  reference_arm = "control",
  covariates = character(0)
)

prior_bcva_01 <- mira_prior_long(
  stan_data_bcva_01,
  outcome = "BCVA"
)

print(prior_bcva_01)

fit_bcva_01 <- mira_fit_long(
  stan_data = stan_data_bcva_01,
  prior = prior_bcva_01,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  seed = 123,
  refresh = 100,
  verbose = TRUE
)

print(fit_bcva_01$diagnostic_summary())
fit_bcva_01$cmdstan_diagnose()

mira_res_bcva_01 <- mira_summary_long(
  fit = fit_bcva_01,
  stan_data = stan_data_bcva_01,
  verbose = TRUE
)


# ============================================================
# BCVA 2 - Tempo, trattamento ed età continua
# ============================================================
# `age_threshold` non è usato: nell'API corrente è deprecato e ignorato.
# Qui l'età entra direttamente come covariata continua.

data_bcva_02 <- data[c(
  "patient",
  "arm",
  "age",
  bcva_visits
)]

stan_data_bcva_02 <- mira_data_long(
  data = data_bcva_02,
  time_value = c(0, 3, 5, 12, 15),
  outcome = "BCVA",
  meaningful_change = 5,
  meaningful_change_sd = 1.5,
  direction = "higher",
  reference_arm = "control",
  covariates = "age"
)

prior_bcva_02 <- mira_prior_long(
  stan_data_bcva_02,
  outcome = "BCVA"
)

print(prior_bcva_02)

fit_bcva_02 <- mira_fit_long(
  stan_data = stan_data_bcva_02,
  prior = prior_bcva_02,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  seed = 124,
  refresh = 100,
  verbose = TRUE
)

print(fit_bcva_02$diagnostic_summary())
fit_bcva_02$cmdstan_diagnose()

mira_res_bcva_02 <- mira_summary_long(
  fit = fit_bcva_02,
  stan_data = stan_data_bcva_02,
  verbose = TRUE
)


# ============================================================
# BCVA 3 - Età e genere
# ============================================================

data_bcva_03 <- data[c(
  "patient",
  "arm",
  "age",
  "gender",
  bcva_visits
)]

stan_data_bcva_03 <- mira_data_long(
  data = data_bcva_03,
  time_value = c(0, 3, 5, 12, 15),
  outcome = "BCVA",
  meaningful_change = 5,
  meaningful_change_sd = 1.5,
  direction = "higher",
  reference_arm = "control",
  covariates = c("age", "gender"),
  covariate_reference_levels = c(
    gender = "Female"
  )
)

prior_bcva_03 <- mira_prior_long(
  stan_data_bcva_03,
  outcome = "BCVA"
)

print(prior_bcva_03)

fit_bcva_03 <- mira_fit_long(
  stan_data = stan_data_bcva_03,
  prior = prior_bcva_03,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  seed = 125,
  refresh = 100,
  verbose = TRUE
)

print(fit_bcva_03$diagnostic_summary())
fit_bcva_03$cmdstan_diagnose()

mira_res_bcva_03 <- mira_summary_long(
  fit = fit_bcva_03,
  stan_data = stan_data_bcva_03,
  verbose = TRUE
)


# ============================================================
# BCVA 4 - Modello clinico completo, covariate automatiche
# ============================================================

data_bcva_04 <- data[c(
  "patient",
  "arm",
  "age",
  "gender",
  "study_eye",
  "diabetes_duration_years",
  "hba1c_percent",
  bcva_visits
)]

stan_data_bcva_04 <- mira_data_long(
  data = data_bcva_04,
  time_value = c(0, 3, 5, 12, 15),
  outcome = "BCVA",
  meaningful_change = 5,
  meaningful_change_sd = 1.5,
  direction = "higher",
  reference_arm = "control",
  covariates = "auto",
  covariate_reference_levels = c(
    gender = "Female",
    study_eye = "OD"
  )
)

prior_bcva_04 <- mira_prior_long(
  stan_data_bcva_04,
  outcome = "BCVA"
)

print(prior_bcva_04)

fit_bcva_04 <- mira_fit_long(
  stan_data = stan_data_bcva_04,
  prior = prior_bcva_04,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  seed = 126,
  refresh = 100,
  verbose = TRUE
)

print(fit_bcva_04$diagnostic_summary())
fit_bcva_04$cmdstan_diagnose()

mira_res_bcva_04 <- mira_summary_long(
  fit = fit_bcva_04,
  stan_data = stan_data_bcva_04,
  verbose = TRUE
)


# ============================================================
# CMT 1 - Solo tempo e trattamento, nessuna covariata
# ============================================================

data_cmt_01 <- data[c(
  "patient",
  "arm",
  cmt_visits
)]

stan_data_cmt_01 <- mira_data_long(
  data = data_cmt_01,
  time_value = c(0, 3, 5, 12, 15),
  outcome = "CMT",
  meaningful_change = 50,
  meaningful_change_sd = 15,
  direction = "lower",
  reference_arm = "control",
  covariates = character(0)
)

prior_cmt_01 <- mira_prior_long(
  stan_data_cmt_01,
  outcome = "CMT"
)

print(prior_cmt_01)

fit_cmt_01 <- mira_fit_long(
  stan_data = stan_data_cmt_01,
  prior = prior_cmt_01,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  seed = 223,
  refresh = 100,
  verbose = TRUE
)

print(fit_cmt_01$diagnostic_summary())
fit_cmt_01$cmdstan_diagnose()

mira_res_cmt_01 <- mira_summary_long(
  fit = fit_cmt_01,
  stan_data = stan_data_cmt_01,
  verbose = TRUE
)


# ============================================================
# CMT 2 - Tempo, trattamento e HbA1c
# ============================================================

data_cmt_02 <- data[c(
  "patient",
  "arm",
  "hba1c_percent",
  cmt_visits
)]

stan_data_cmt_02 <- mira_data_long(
  data = data_cmt_02,
  time_value = c(0, 3, 5, 12, 15),
  outcome = "CMT",
  meaningful_change = 50,
  meaningful_change_sd = 15,
  direction = "lower",
  reference_arm = "control",
  covariates = "hba1c_percent"
)

prior_cmt_02 <- mira_prior_long(
  stan_data_cmt_02,
  outcome = "CMT"
)

print(prior_cmt_02)

fit_cmt_02 <- mira_fit_long(
  stan_data = stan_data_cmt_02,
  prior = prior_cmt_02,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  seed = 224,
  refresh = 100,
  verbose = TRUE
)

print(fit_cmt_02$diagnostic_summary())
fit_cmt_02$cmdstan_diagnose()

mira_res_cmt_02 <- mira_summary_long(
  fit = fit_cmt_02,
  stan_data = stan_data_cmt_02,
  verbose = TRUE
)


# ============================================================
# CMT 3 - Durata del diabete e occhio dello studio
# ============================================================

data_cmt_03 <- data[c(
  "patient",
  "arm",
  "diabetes_duration_years",
  "study_eye",
  cmt_visits
)]

stan_data_cmt_03 <- mira_data_long(
  data = data_cmt_03,
  time_value = c(0, 3, 5, 12, 15),
  outcome = "CMT",
  meaningful_change = 50,
  meaningful_change_sd = 15,
  direction = "lower",
  reference_arm = "control",
  covariates = c("diabetes_duration_years", "study_eye"),
  covariate_reference_levels = c(
    study_eye = "OD"
  )
)

prior_cmt_03 <- mira_prior_long(
  stan_data_cmt_03,
  outcome = "CMT"
)

print(prior_cmt_03)

fit_cmt_03 <- mira_fit_long(
  stan_data = stan_data_cmt_03,
  prior = prior_cmt_03,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  seed = 225,
  refresh = 100,
  verbose = TRUE
)

print(fit_cmt_03$diagnostic_summary())
fit_cmt_03$cmdstan_diagnose()

mira_res_cmt_03 <- mira_summary_long(
  fit = fit_cmt_03,
  stan_data = stan_data_cmt_03,
  verbose = TRUE
)


# ============================================================
# CMT 4 - Modello clinico completo, covariate automatiche
# ============================================================

data_cmt_04 <- data[c(
  "patient",
  "arm",
  "age",
  "gender",
  "study_eye",
  "diabetes_duration_years",
  "hba1c_percent",
  cmt_visits
)]

stan_data_cmt_04 <- mira_data_long(
  data = data_cmt_04,
  time_value = c(0, 3, 5, 12, 15),
  outcome = "CMT",
  meaningful_change = 50,
  meaningful_change_sd = 15,
  direction = "lower",
  reference_arm = "control",
  covariates = "auto",
  covariate_reference_levels = c(
    gender = "Female",
    study_eye = "OD"
  )
)

prior_cmt_04 <- mira_prior_long(
  stan_data_cmt_04,
  outcome = "CMT"
)

print(prior_cmt_04)

fit_cmt_04 <- mira_fit_long(
  stan_data = stan_data_cmt_04,
  prior = prior_cmt_04,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  seed = 226,
  refresh = 100,
  verbose = TRUE
)

print(fit_cmt_04$diagnostic_summary())
fit_cmt_04$cmdstan_diagnose()

mira_res_cmt_04 <- mira_summary_long(
  fit = fit_cmt_04,
  stan_data = stan_data_cmt_04,
  verbose = TRUE
)


# ============================================================
# Oggetti finali disponibili
# ============================================================
#
# BCVA:
#   stan_data_bcva_01, prior_bcva_01, fit_bcva_01, mira_res_bcva_01
#   stan_data_bcva_02, prior_bcva_02, fit_bcva_02, mira_res_bcva_02
#   stan_data_bcva_03, prior_bcva_03, fit_bcva_03, mira_res_bcva_03
#   stan_data_bcva_04, prior_bcva_04, fit_bcva_04, mira_res_bcva_04
#
# CMT:
#   stan_data_cmt_01, prior_cmt_01, fit_cmt_01, mira_res_cmt_01
#   stan_data_cmt_02, prior_cmt_02, fit_cmt_02, mira_res_cmt_02
#   stan_data_cmt_03, prior_cmt_03, fit_cmt_03, mira_res_cmt_03
#   stan_data_cmt_04, prior_cmt_04, fit_cmt_04, mira_res_cmt_04
