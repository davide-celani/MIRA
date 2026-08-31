# ============================================================
# DATASET DETERMINISTICO PER PARAMETER RECOVERY
# ============================================================

# ------------------------------------------------------------
# Struttura dello studio
# ------------------------------------------------------------

n_rep <- 20

# 2 arm x 2 gender x 2 gruppi età x 20 = 160 pazienti
n_patients <- 2 * 2 * 2 * n_rep

time_value <- 0:4
K <- length(time_value)

age_threshold <- 60


# ============================================================
# PARAMETRI VERI
# ============================================================

# ------------------------------------------------------------
# Traiettoria del gruppo di riferimento:
#
# control
# Female
# age <= 60
# ------------------------------------------------------------

true_reference <- c(
  15.0,   # t0
  14.7,   # t1
  14.2,   # t2
  13.9,   # t3
  13.7    # t4
)


# ------------------------------------------------------------
# Effetto treatment - control a ogni tempo
#
# Valori negativi = IOP più bassa nel treatment
# ------------------------------------------------------------

true_treatment_effect <- c(
  0.0,   # t0
  -0.6,   # t1
  -1.2,   # t2
  -1.7,   # t3
  -2.0    # t4
)


# ------------------------------------------------------------
# Effetto Male - Female
# ------------------------------------------------------------

true_male_effect <- c(
  0.8,    # t0
  0.7,    # t1
  0.5,    # t2
  0.4,    # t3
  0.3     # t4
)


# ------------------------------------------------------------
# Effetto age > 60 rispetto ad age <= 60
# ------------------------------------------------------------

true_older_effect <- c(
  1.1,    # t0
  1.0,    # t1
  0.9,    # t2
  0.8,    # t3
  0.7     # t4
)


# ------------------------------------------------------------
# Eterogeneità tra soggetti
# ------------------------------------------------------------

true_sigma_subject_intercept <- 1.0
true_sigma_subject_slope     <- 0.10


# ------------------------------------------------------------
# Residui Student-t
# ------------------------------------------------------------

true_sigma <- 0.35
true_nu    <- 8


# ============================================================
# DISEGNO COMPLETAMENTE BILANCIATO
# ============================================================

design <- expand.grid(

  rep = seq_len(n_rep),

  arm = c(
    "control",
    "treatment"
  ),

  gender = c(
    "Female",
    "Male"
  ),

  age_group = c(
    "younger",
    "older"
  ),

  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)


design$patient <- seq_len(nrow(design))


# ============================================================
# ETÀ DETERMINISTICA
# ============================================================

# Younger: 41, 42, ..., 60
#
# Older:
# 61, 62, ..., 80

design$age <- ifelse(
  design$age_group == "younger",
  40 + design$rep,
  60 + design$rep
)


# Controllo

table(
  design$age > age_threshold
)

table(
  design$arm,
  design$gender,
  design$age_group
)


# ============================================================
# RANDOM EFFECT "DETERMINISTICI"
# ============================================================
#
# Non utilizziamo rnorm().
#
# Prendiamo invece quantili prefissati della distribuzione
# normale. In questo modo ogni esecuzione produce ESATTAMENTE
# gli stessi soggetti.
# ============================================================

z_subject <- qnorm(
  (seq_len(n_rep) - 0.5) / n_rep
)

# Forziamo media = 0 e SD = 1 esattamente nel campione

z_subject <- as.numeric(
  scale(
    z_subject,
    center = TRUE,
    scale = TRUE
  )
)


# ------------------------------------------------------------
# Random intercept
# ------------------------------------------------------------

subject_intercept_rep <-
  true_sigma_subject_intercept *
  z_subject


# ------------------------------------------------------------
# Random slope
#
# Permutiamo deterministicamente i quantili per evitare che
# intercept e slope siano fortemente correlati.
# ------------------------------------------------------------

slope_index <- (
  ((seq_len(n_rep) - 1L) * 13L + 16L) %%
    n_rep
) + 1L


subject_slope_rep <-
  true_sigma_subject_slope *
  z_subject[slope_index]


# Assegniamo gli stessi pattern in ogni cella sperimentale.
# Questo è importante per non confondere treatment/gender/age
# con l'eterogeneità soggettiva.

design$subject_intercept <-
  subject_intercept_rep[design$rep]

design$subject_slope <-
  subject_slope_rep[design$rep]


# ============================================================
# RESIDUI DETERMINISTICI STUDENT-t
# ============================================================
#
# Anche qui NON usiamo rt().
#
# Usiamo quantili deterministici della Student-t.
# ============================================================

residual_a <- c(
  7L,
  13L,
  17L,
  9L,
  3L
)

residual_b <- c(
  3L,
  16L,
  8L,
  7L,
  11L
)


# ============================================================
# GENERAZIONE DELLE IOP
# ============================================================

for (k in seq_len(K)) {

  # ----------------------------------------------------------
  # Effetto treatment
  # ----------------------------------------------------------

  treatment_component <- ifelse(
    design$arm == "treatment",
    true_treatment_effect[k],
    0
  )


  # ----------------------------------------------------------
  # Effetto gender
  # ----------------------------------------------------------

  gender_component <- ifelse(
    design$gender == "Male",
    true_male_effect[k],
    0
  )


  # ----------------------------------------------------------
  # Effetto gruppo età
  # ----------------------------------------------------------

  age_component <- ifelse(
    design$age > age_threshold,
    true_older_effect[k],
    0
  )


  # ----------------------------------------------------------
  # Random intercept + random slope
  # ----------------------------------------------------------

  subject_component <-
    design$subject_intercept +
    design$subject_slope * time_value[k]


  # ----------------------------------------------------------
  # Residuo deterministico
  # ----------------------------------------------------------

  residual_index <- (
    (
      (design$rep - 1L) * residual_a[k] +
        residual_b[k]
    ) %%
      n_rep
  ) + 1L


  residual_probability <- (
    residual_index - 0.5
  ) / n_rep


  residual <- true_sigma * qt(
    residual_probability,
    df = true_nu
  )


  # ----------------------------------------------------------
  # Media vera
  # ----------------------------------------------------------

  mu <-
    true_reference[k] +
    treatment_component +
    gender_component +
    age_component +
    subject_component


  # ----------------------------------------------------------
  # Outcome osservato
  # ----------------------------------------------------------

  data_column <- paste0(
    "IOP_t",
    k - 1
  )


  design[[data_column]] <-
    mu +
    residual
}


# ============================================================
# DATASET FINALE DA PASSARE A MIRA
# ============================================================

data <- design[
  ,
  c(
    "patient",
    "arm",
    "gender",
    "age",
    paste0("IOP_t", 0:4)
  )
]


# ============================================================
# CONTROLLO
# ============================================================

head(data)

dim(data)

table(
  data$arm,
  data$gender,
  data$age > age_threshold
)
