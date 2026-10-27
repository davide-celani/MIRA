rank_longitudinal_call <- function(...) {
  # This also lets the standalone development runner source the two files
  # without changing an installed namespace. Normal tests resolve MIRA's env.
  get(".mira_run_rank_longitudinal_long", envir = environment(mira_info_long),
    inherits = FALSE)(...)
}

rank_longitudinal_small <- function() {
  data.frame(patient = rep(c("P3", "P1", "P2"), 2L),
    time = rep(c("score_t0", "score_t1"), each = 3L),
    time_index = rep(1:2, each = 3L),
    time_label = rep(c("Zbase", "Afinal"), each = 3L),
    value = c(1, 3, 3, 2, 4, NA_real_))
}

rank_longitudinal_arm_data <- function() {
  values <- rbind(c(1, 3, 4), c(4, 5, 7), c(6, NA, 8),
    c(2, 4, 8), c(5, 7, 9), c(8, 9, 11), c(10, 12, 14))
  data.frame(patient = rep(c("Z3", "C1", "M2", "A4", "B2", "X1", "F8"), each = 3L),
    time = rep(c("score_t0", "score_t1", "score_t2"), 7L),
    time_index = rep(1:3, 7L),
    time_label = rep(c("Zbase", "Amid", "Mfinal"), 7L),
    arm = rep(c(rep("control", 3L), rep("active", 4L)), each = 3L),
    value = as.numeric(t(values)),
    unused_covariate = c(NA, rep(11, 19L), NA))
}

test_that("rank relative effects, uncertainty and global tests have a manual missing-data oracle", {
  data <- rank_longitudinal_small()
  original <- data
  result <- rank_longitudinal_call(data, c("Zbase", "Afinal"), FALSE, character(), .90)
  expect_true(result$performed)
  expect_null(result$error)
  expect_null(result$reason_skipped)
  expect_s3_class(result$object, "nparld_fit")
  expect_identical(result$n_subjects, 3L)
  expect_identical(result$n_observations, 5L)
  effects <- result$effects
  expect_identical(nrow(effects), 2L)
  expect_identical(effects$time, c("score_t0", "score_t1"))
  expect_identical(effects$time_label, c("Zbase", "Afinal"))
  expect_identical(effects$n, c(3L, 2L))
  expect_identical(effects$missing_n, c(0L, 1L))
  expect_identical(effects$n_subjects, c(3L, 3L))
  # Of six cross-visit comparisons, two place an observation at the first
  # visit above the second. Self-distribution comparisons contribute 1/2.
  # With G=(F_base+F_final)/2, p=(5/12,7/12), not the allocation-weighted
  # values (13/30,3/5). The subject influence values for p_base are
  # (1/12,-1/6,1/12), after available-cell weighting 3/3 and 3/2.
  # Thus Var(p_base)=sum(influence^2)/(3*(3-1))=1/144; p_final=1-p_base.
  p <- c(5, 7) / 12
  expect_equal(effects$relative_effect, p, tolerance = 1e-12)
  expect_equal(effects$SE, rep(1 / 12, 2L), tolerance = 1e-12)
  expect_equal(unname(result$covariance), matrix(c(1, -1, -1, 1), 2L) / 144,
    tolerance = 1e-12)
  z <- stats::qnorm(.95)
  logit_se <- (1 / 12) / (p * (1 - p))
  expect_equal(effects$ci_lower, stats::plogis(stats::qlogis(p) - z * logit_se),
    tolerance = 1e-12)
  expect_equal(effects$ci_upper, stats::plogis(stats::qlogis(p) + z * logit_se),
    tolerance = 1e-12)
  # Native nparLD's field is sqrt(N) times estimator SE. The wrapper must not
  # mistake that asymptotic scaling for the SE displayed in MIRA/report.
  expect_equal(result$object$effects$std.error, rep(sqrt(3) / 12, 2L), tolerance = 1e-12)
  tests <- result$tests
  expect_identical(paste(tests$question, tests$method), c("TIME ATS", "TIME WTS"))
  expect_true(all(tests$estimable))
  expect_equal(tests$statistic, c(1, 1), tolerance = 1e-12)
  expect_equal(tests$df_num, c(1, 1), tolerance = 1e-12)
  expect_identical(tests$df_den, c(1e10, NA_real_))
  expect_equal(tests$p_raw, c(stats::pf(1, 1, 1e10, lower.tail = FALSE),
    stats::pchisq(1, 1, lower.tail = FALSE)), tolerance = 1e-12)
  expect_identical(data, original)
})

test_that("rank factorial inference preserves subject-visit keys and ignores nuisance missingness", {
  data <- rank_longitudinal_arm_data()
  visits <- c("Zbase", "Amid", "Mfinal")
  arms <- c("control", "active")
  result <- rank_longitudinal_call(data, visits, TRUE, arms, .90)
  expect_true(result$performed)
  expect_null(result$error)
  expect_identical(result$n_subjects, 7L)
  expect_identical(result$n_observations, 20L)
  expect_true(result$unadjusted)
  expect_identical(nrow(result$data), 21L)
  key <- function(x) paste(x$patient, x$time, sep = "::")
  expect_identical(anyDuplicated(key(result$data)), 0L)
  expect_identical(key(result$data), key(data))
  expect_identical(result$data$value, data$value)
  expect_false("unused_covariate" %in% names(result$data))
  effects <- result$effects
  expected <- data.frame(arm = rep(arms, each = 3L),
    time_label = rep(visits, 2L), n = c(3L, 2L, 3L, 4L, 4L, 4L))
  expect_identical(nrow(effects), 6L)
  expect_identical(anyDuplicated(effects[c("arm", "time_label")]), 0L)
  expect_identical(effects$arm, expected$arm)
  expect_identical(effects$time_label, expected$time_label)
  expect_identical(effects$n, expected$n)
  # Direct engine oracle verifies the wrapper's factorial extraction/term
  # association. It complements, not replaces, the independent arithmetic above.
  reference <- nparLD::nparLD(value ~ arm_factor * time_factor, data = result$data,
    subject = "patient_factor", effect = "unweighted", hypothesis = "H0p",
    cell.weights = "subjects", alpha = .10, CI.method = "logit",
    covariance = TRUE, contrast = NULL, perm.test = FALSE)
  expect_identical(result$object$effects, reference$effects)
  tests <- result$tests
  expect_identical(nrow(tests), 6L)
  expect_identical(anyDuplicated(tests[c("question", "method")]), 0L)
  expect_identical(tests$question, rep(c("TIME", "ARM", "TIME_X_ARM"), each = 2L))
  expect_identical(tests$term,
    rep(c("time_factor", "arm_factor", "arm_factor:time_factor"), each = 2L))
  for (i in seq_len(nrow(tests))) {
    table <- reference[[tests$method[[i]]]]
    term <- tests$term[[i]]
    expect_identical(sum(rownames(table) == term), 1L)
    expect_equal(c(tests$statistic[[i]], tests$df_num[[i]], tests$p_raw[[i]]),
      as.numeric(table[term, c("Statistic", "df", "p-value")]), tolerance = 1e-12)
  }
})

test_that("rank inference is invariant to record order and bijective subject renaming", {
  data <- rank_longitudinal_arm_data()
  visits <- c("Zbase", "Amid", "Mfinal")
  result <- rank_longitudinal_call(data, visits, TRUE, c("control", "active"), .90)
  permuted <- data[rev(seq_len(nrow(data))), , drop = FALSE]
  ids <- unique(permuted$patient)
  permuted$patient <- paste0("Renamed", match(permuted$patient, ids))
  changed <- rank_longitudinal_call(permuted, visits, TRUE, c("control", "active"), .90)
  expect_true(changed$performed)
  expect_equal(changed$effects, result$effects, tolerance = 1e-12)
  expect_equal(changed$tests, result$tests, tolerance = 1e-12)
  expect_equal(changed$covariance, result$covariance, tolerance = 1e-12)
})

test_that("rank sampling excludes wholly unobserved profiles without imputing partial profiles", {
  data <- rank_longitudinal_small()
  result <- rank_longitudinal_call(data, c("Zbase", "Afinal"), FALSE, character(), .90)
  absent <- data[1:2, ]
  absent$patient <- "unobserved"
  absent$time <- c("score_t0", "score_t1")
  absent$time_label <- c("Zbase", "Afinal")
  absent$time_index <- 1:2
  absent$value <- c(Inf, NaN)
  changed <- rank_longitudinal_call(rbind(data, absent), c("Zbase", "Afinal"),
    FALSE, character(), .90)
  expect_identical(changed$excluded_subjects, "unobserved")
  expect_identical(changed$n_subjects, 3L)
  expect_identical(changed$n_observations, 5L)
  expect_equal(changed$effects, result$effects, tolerance = 0)
  expect_equal(changed$tests, result$tests, tolerance = 0)
  expect_identical(changed$data, result$data)
  expect_true(is.na(changed$data$value[changed$data$patient == "P2" &
    changed$data$time == "score_t1"]))
})

test_that("rank adaptation exposes unsupported and insufficient designs without fake inference", {
  data <- rank_longitudinal_small()
  run <- function(x) rank_longitudinal_call(x, c("Zbase", "Afinal"), FALSE, character(), .90)
  duplicate <- run(rbind(data, data[1L, ]))
  expect_false(duplicate$performed)
  expect_match(duplicate$reason_skipped, "Duplicate subject-by-visit")
  expect_null(duplicate$object)
  expect_null(duplicate$error)
  incomplete <- run(data[-6L, ])
  expect_false(incomplete$performed)
  expect_match(incomplete$reason_skipped, "full subject-by-visit schedule")
  sparse <- data
  sparse$value[4L] <- NA_real_
  insufficient <- run(sparse)
  expect_false(insufficient$performed)
  expect_match(insufficient$reason_skipped, "two finite observations")
  tied <- data
  tied$value[] <- 7
  unavailable <- run(tied)
  expect_false(unavailable$performed)
  expect_match(unavailable$reason_skipped, "all finite responses are tied")
  expect_identical(nrow(unavailable$effects), 0L)
  expect_identical(nrow(unavailable$tests), 0L)
  arm <- rank_longitudinal_arm_data()
  arm$arm[2L] <- "active"
  invalid_arm <- rank_longitudinal_call(arm, c("Zbase", "Amid", "Mfinal"),
    TRUE, c("control", "active"), .90)
  expect_false(invalid_arm$performed)
  expect_match(invalid_arm$reason_skipped, "constant treatment arm")
})

test_that("rank dependency and engine failures remain structured", {
  data <- rank_longitudinal_small()
  unavailable <- testthat::with_mocked_bindings(
    rank_longitudinal_call(data, c("Zbase", "Afinal"), FALSE, character(), .90),
    requireNamespace = function(...) FALSE, .package = "base")
  expect_false(unavailable$performed)
  expect_match(unavailable$reason_skipped, "not installed")
  old <- testthat::with_mocked_bindings(
    rank_longitudinal_call(data, c("Zbase", "Afinal"), FALSE, character(), .90),
    packageVersion = function(...) package_version("2.2"), .package = "utils")
  expect_false(old$performed)
  expect_match(old$reason_skipped, "different API")
  failed <- testthat::with_mocked_bindings(
    rank_longitudinal_call(data, c("Zbase", "Afinal"), FALSE, character(), .90),
    nparLD = function(...) stop("controlled rank engine failure"), .package = "nparLD")
  expect_false(failed$performed)
  expect_identical(failed$error, "controlled rank engine failure")
  expect_null(failed$object)
  expect_identical(failed$data$value, c(1, 2, 3, 4, 3, NA_real_))
  expect_true(all(c("error", "reason_skipped", "object", "effects", "tests") %in% names(failed)))
})

test_that("rank inference does not turn nonfinite engine tests into a performed analysis", {
  data <- rank_longitudinal_small()
  original <- nparLD::nparLD
  result <- testthat::with_mocked_bindings(
    rank_longitudinal_call(data, c("Zbase", "Afinal"), FALSE, character(), .90),
    nparLD = function(...) {
      object <- original(...)
      object$ATS[, "p-value"] <- NaN
      object$WTS[, "p-value"] <- NaN
      object
    }, .package = "nparLD")
  expect_false(result$performed)
  expect_match(result$reason_skipped, "no finite estimable global test")
  expect_false(any(result$tests$estimable))
  expect_true(all(is.na(as.matrix(result$tests[c("statistic", "df_num", "df_den", "p_raw")]))))
  expect_true(all(grepl("engine returned", result$tests$note, fixed = TRUE)))
  expect_true(all(is.nan(result$object$ATS[, "p-value"])))
  expect_s3_class(result$object, "nparld_fit")
  expect_equal(result$effects$relative_effect, c(5, 7) / 12, tolerance = 1e-12)
})

test_that("static complete rank profiles retain valid arm inference without inventing time inference", {
  data <- rank_longitudinal_arm_data()
  # The groups overlap and differ; every subject's profile is exactly static.
  # Thus ARM remains testable while TIME and TIME_X_ARM have no within-subject
  # variation. This structural property does not depend on roundoff in nparLD.
  values <- c(1, 4, 6, 3, 5, 8, 11)
  data$value <- rep(values, each = 3L)
  result <- rank_longitudinal_call(data, c("Zbase", "Amid", "Mfinal"),
    TRUE, c("control", "active"), .90)
  expect_true(result$performed)
  expect_null(result$reason_skipped)
  arm <- result$tests[result$tests$question == "ARM", , drop = FALSE]
  time <- result$tests[result$tests$question != "ARM", , drop = FALSE]
  expect_identical(nrow(arm), 2L)
  expect_true(all(arm$estimable))
  expect_true(all(is.finite(arm$p_raw)))
  expect_identical(nrow(time), 4L)
  expect_false(any(time$estimable))
  expect_true(all(is.na(as.matrix(time[c("statistic", "df_num", "df_den", "p_raw")]))))
  expect_true(all(grepl("no within-subject variation", time$note, fixed = TRUE)))
  # Retain native engine evidence, including any finite p=1 due to numerical
  # regularization. MIRA's own rows must not advertise that as time inference.
  expect_s3_class(result$object, "nparld_fit")
  expect_true(all(c("time_factor", "arm_factor:time_factor") %in% rownames(result$object$ATS)))
  pooled <- rank_longitudinal_call(data, c("Zbase", "Amid", "Mfinal"),
    FALSE, character(), .90)
  expect_false(pooled$performed)
  expect_match(pooled$reason_skipped, "no within-subject variation")
  expect_false(any(pooled$tests$estimable))
})

test_that("deterministic rank fitting preserves the caller RNG state", {
  withr::local_preserve_seed()
  set.seed(20261026)
  seed <- .Random.seed
  kind <- RNGkind()
  result <- rank_longitudinal_call(rank_longitudinal_small(), c("Zbase", "Afinal"),
    FALSE, character(), .90)
  expect_true(result$performed)
  expect_identical(.Random.seed, seed)
  expect_identical(RNGkind(), kind)
  rm(".Random.seed", envir = .GlobalEnv)
  absent <- rank_longitudinal_call(rank_longitudinal_small(), c("Zbase", "Afinal"),
    FALSE, character(), .90)
  expect_true(absent$performed)
  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
  expect_identical(RNGkind(), kind)
})
