# Direct tests for the configuration and data-preparation helpers used by mira_info().

test_that("longitudinal names recognize documented time suffixes and separators", {
  cases <- list(
    c("score_t0", "score", "t0", "0", "t"),
    c("score.time1", "score", "time1", "1", "time"),
    c("score-visit2", "score", "visit2", "2", "time"),
    c("score_baseline", "score", "baseline", "0", "baseline"),
    c("score_month3", "score", "month3", "3", "month"),
    c("score.week4", "score", "week4", "4", "week"),
    c("score-day7", "score", "day7", "7", "day"),
    c("score_fu1", "score", "fu1", "1", "followup"),
    c("score_2", "score", "2", "2", "numeric"),
    c("SCORE_T5", "SCORE", "T5", "5", "t")
  )
  for (case in cases) {
    found <- MIRA:::.mira_parse_longitudinal_name(case[[1L]])
    expect_equal(found$outcome, case[[2L]], info = case[[1L]])
    expect_equal(found$time_label, case[[3L]], info = case[[1L]])
    expect_equal(found$time_order, as.numeric(case[[4L]]), info = case[[1L]])
    expect_equal(found$pattern, case[[5L]], info = case[[1L]])
  }
  expect_null(MIRA:::.mira_parse_longitudinal_name("age"))
  expect_null(MIRA:::.mira_parse_longitudinal_name("score_time"))
  expect_equal(MIRA:::.mira_time_order("BL"), 0)
  expect_true(is.na(MIRA:::.mira_time_order("followup")))
})

test_that("custom longitudinal regex and parser functions retain supplied metadata", {
  regex <- "^(.+)__at__(visit[0-9]+)$"
  hit <- MIRA:::.mira_parse_longitudinal_name("pain__at__visit12", regex)
  expect_equal(hit[c("outcome", "time_label", "time_order", "pattern")],
               list(outcome = "pain", time_label = "visit12",
                    time_order = 12, pattern = "custom_regex"))
  expect_null(MIRA:::.mira_parse_longitudinal_name("pain_t1", regex))
  expect_null(MIRA:::.mira_parse_longitudinal_name("pain_12", "^(.+)_([0-9]+)_x$"))

  parser <- function(x) {
    if (x == "ignore") return(NULL)
    list(outcome = "pain", time_label = "week6", time_order = 99)
  }
  hit <- MIRA:::.mira_parse_longitudinal_name("measurement", parser)
  expect_equal(hit$time_order, 99)
  expect_equal(hit$pattern, "custom_function")
  expect_null(MIRA:::.mira_parse_longitudinal_name("ignore", parser))

  dataframe_parser <- function(x) data.frame(outcome = "pain", time_label = "month3")
  hit <- MIRA:::.mira_parse_longitudinal_name("measurement", dataframe_parser)
  expect_equal(hit$time_order, 3)
  expect_equal(hit$pattern, "custom_function")
  expect_equal(MIRA:::.mira_parse_longitudinal_name("measurement",
    function(x) list(outcome = "pain"))$time_label, "measurement")

  expect_error(MIRA:::.mira_parse_longitudinal_name("x", function(x) "bad"),
               "list containing at least 'outcome'")
  expect_error(MIRA:::.mira_parse_longitudinal_name("x", function(x) list(time_label = "t1")),
               "list containing at least 'outcome'")
  expect_error(MIRA:::.mira_parse_longitudinal_name("x", 1), "variable_pattern")
  expect_error(MIRA:::.mira_parse_longitudinal_name("x", NA_character_), "variable_pattern")
})

test_that("explicit longitudinal groups sort clear time order and preserve ambiguous order", {
  dat <- data.frame(score_t0 = 1:6, score_t1 = 2:7, score_t2 = 3:8)
  ordered <- MIRA:::.mira_make_longitudinal_group(
    dat, c("score_t2", "score_t0", "score_t1"))
  expect_equal(ordered$variables, c("score_t0", "score_t1", "score_t2"))
  expect_equal(ordered$input_variables, c("score_t2", "score_t0", "score_t1"))
  expect_equal(ordered$sort_index, c(2L, 3L, 1L))
  expect_equal(ordered$outcome, "score")
  expect_false(ordered$ambiguous_order)

  ambiguous <- MIRA:::.mira_make_longitudinal_group(
    dat, c("score_t2", "score_t0"), variable_pattern = function(x) NULL)
  expect_equal(ambiguous$variables, c("score_t2", "score_t0"))
  expect_true(ambiguous$ambiguous_order)
  expect_equal(ambiguous$outcome, "score_t")

  expect_error(MIRA:::.mira_make_longitudinal_group(dat, "score_t0"), "at least two")
  expect_error(MIRA:::.mira_make_longitudinal_group(dat, c("score_t0", "score_t0")),
               "duplicate")
  expect_error(MIRA:::.mira_make_longitudinal_group(dat, c("score_t0", "missing")),
               "not found")
  dat$score_t3 <- letters[1:6]
  expect_error(MIRA:::.mira_make_longitudinal_group(dat, c("score_t0", "score_t3")),
               "must be numeric")
})

test_that("automatic longitudinal detection separates valid, ambiguous, and nonnumeric matches", {
  dat <- data.frame(
    pain_t2 = 1:8, pain_t0 = 2:9, pain_t1 = 3:10,
    fatigue_t1 = 1:8, fatigue_time1 = 2:9,
    vision_t0 = 1:8, vision_t1 = letters[1:8],
    solo_t0 = 1:8, unrelated = 1:8
  )
  result <- MIRA:::.mira_detect_longitudinal_variables(dat)
  expect_equal(result$groups$pain$variables, c("pain_t0", "pain_t1", "pain_t2"))
  expect_true("fatigue" %in% names(result$ambiguous_groups))
  expect_false("fatigue" %in% names(result$groups))
  expect_false("solo" %in% names(result$groups))
  expect_false("vision" %in% names(result$groups))
  expect_true("vision_t1" %in% result$non_numeric_matches)
  expect_false("unrelated" %in% result$parsed$variable)

  duplicate_labels <- data.frame(pain_t1 = 1:8, pain_T1 = 2:9)
  found <- MIRA:::.mira_detect_longitudinal_variables(duplicate_labels)
  expect_true("pain" %in% names(found$ambiguous_groups))
})

test_that("outcome resolution accepts explicit columns, grouped lists, and named labels", {
  dat <- data.frame(pain_t0 = 1:8, pain_t1 = 2:9,
                    vision_t0 = 3:10, vision_t1 = 4:11)
  detected <- MIRA:::.mira_detect_longitudinal_variables(dat)
  selected <- MIRA:::.mira_resolve_outcome_groups(
    dat, detected, outcomes = "vision")
  expect_named(selected$groups, "vision")
  expect_equal(selected$groups$vision$variables, c("vision_t0", "vision_t1"))

  by_columns <- MIRA:::.mira_resolve_outcome_groups(
    dat, detected, outcomes = c("pain_t1", "pain_t0"))
  expect_named(by_columns$groups, "pain")
  expect_equal(by_columns$groups$pain$variables, c("pain_t0", "pain_t1"))

  manual <- MIRA:::.mira_resolve_outcome_groups(
    dat, detected,
    time_vars = list(pain = c("pain_t1", "pain_t0"),
                     vision = c("vision_t0", "vision_t1")),
    time_labels = list(pain = c(pain_t0 = "Baseline", pain_t1 = "Follow-up"),
                       vision = c(vision_t0 = "Start", vision_t1 = "End")))
  expect_equal(unname(manual$groups$pain$time_labels), c("Baseline", "Follow-up"))
  expect_equal(unname(manual$groups$vision$time_labels), c("Start", "End"))
  expect_true(manual$time_vars_manual)

  expect_error(MIRA:::.mira_resolve_outcome_groups(dat, detected, outcomes = "absent"),
               "not detected")
  expect_error(MIRA:::.mira_resolve_outcome_groups(dat, detected,
    time_vars = list(c("pain_t0", "pain_t1"), c("vision_t0", "vision_t1")),
    time_labels = c("Start", "End")), "multiple outcomes")
  expect_error(MIRA:::.mira_resolve_outcome_groups(dat, detected,
    time_vars = c("pain_t0", "pain_t1"), time_labels = c("Same", "Same")),
    "unique labels")
})

test_that("ID detection requires a unique semantic candidate and records ambiguity", {
  dat <- data.frame(patient_id = sprintf("P%02d", 1:12),
                    measurement_t0 = 1:12, measurement_t1 = 2:13)
  chosen <- MIRA:::.mira_detect_id(dat, exclude = c("measurement_t0", "measurement_t1"))
  expect_equal(chosen$selected, "patient_id")
  expect_equal(chosen$source, "auto")
  expect_true(chosen$candidates$structurally_valid[chosen$candidates$variable == "patient_id"])

  manual <- MIRA:::.mira_detect_id(dat, id = "patient_id")
  expect_false(manual$automatic)
  expect_equal(manual$source, "manual")
  expect_error(MIRA:::.mira_detect_id(dat, id = "absent"), "does not exist")

  ambiguous <- dat
  ambiguous$subject_id <- sprintf("S%02d", 1:12)
  tie <- MIRA:::.mira_detect_id(ambiguous,
    exclude = c("measurement_t0", "measurement_t1"))
  expect_null(tie$selected)
  expect_equal(tie$source, "generated")
  expect_match(tie$warnings, "Ambiguous ID")
  expect_setequal(tie$alternatives, c("patient_id", "subject_id"))

  no_semantic <- data.frame(x = 1:12, y_t0 = 2:13, y_t1 = 3:14)
  generated <- MIRA:::.mira_detect_id(no_semantic, exclude = c("y_t0", "y_t1"))
  expect_null(generated$selected)
  expect_equal(generated$source, "generated")

  invalid <- dat
  invalid$patient_id[[2L]] <- invalid$patient_id[[1L]]
  invalid$patient_id[[3L]] <- NA_character_
  detected <- MIRA:::.mira_detect_id(invalid,
    exclude = c("measurement_t0", "measurement_t1"))
  expect_null(detected$selected)
  expect_false(detected$candidates$structurally_valid[
    detected$candidates$variable == "patient_id"])
})

test_that("arm detection checks names and cardinality while manual selection takes precedence", {
  dat <- data.frame(arm = rep(c("control", "active"), each = 6),
                    treatment_arm = rep(c("A", "B"), 6),
                    score_t0 = 1:12, score_t1 = 2:13)
  found <- MIRA:::.mira_detect_arm(dat, exclude = c("score_t0", "score_t1"))
  expect_equal(found$selected, "arm")
  expect_true("treatment_arm" %in% found$alternatives)
  expect_match(found$warnings, "alternatives")

  dat$arm <- sprintf("G%02d", 1:12)
  rejected <- MIRA:::.mira_detect_arm(dat, exclude = c("score_t0", "score_t1"))
  expect_equal(rejected$selected, "treatment_arm")
  expect_match(rejected$warnings, "implausible cardinality")

  dat$arm <- "only"
  manual <- MIRA:::.mira_detect_arm(dat, arm = "arm")
  expect_equal(manual$selected, "arm")
  expect_false(manual$automatic)
  expect_error(MIRA:::.mira_detect_arm(dat, arm = "missing"), "does not exist")
})

test_that("reference arms prefer controls and otherwise the largest observed group", {
  dat <- data.frame(arm = c("drug", "control", "drug", "control", "drug"))
  expect_equal(MIRA:::.mira_choose_reference_arm(dat, "arm")$value, "control")
  expect_equal(MIRA:::.mira_choose_reference_arm(dat, "arm", "drug")$value, "drug")
  expect_false(MIRA:::.mira_choose_reference_arm(dat, "arm", "drug")$automatic)
  expect_error(MIRA:::.mira_choose_reference_arm(dat, "arm", "absent"), "not present")
  expect_error(MIRA:::.mira_choose_reference_arm(dat, NULL, "control"), "requires an arm")

  no_control <- data.frame(arm = c("B", "A", "B", "A", "B"))
  expect_equal(MIRA:::.mira_choose_reference_arm(no_control, "arm")$value, "B")
  multiple_controls <- data.frame(arm = c("control", "placebo", "placebo", "drug"))
  ref <- MIRA:::.mira_choose_reference_arm(multiple_controls, "arm")
  expect_equal(ref$value, "placebo")
  expect_match(ref$warnings, "Multiple plausible")
})

test_that("covariate classification distinguishes continuous and categorical inputs", {
  dat <- data.frame(age = seq(20, 58, by = 2),
                    score = rep(1:4, each = 5),
                    sex = rep(c("F", "M"), 10),
                    smoker = rep(c(TRUE, FALSE), 10),
                    constant = 1)
  classified <- MIRA:::.mira_classify_covariates(dat, names(dat))
  expect_equal(classified$numeric, "age")
  expect_setequal(classified$categorical, c("score", "sex", "smoker"))
  expect_false("constant" %in% unlist(classified, use.names = FALSE))
})

test_that("covariate detection excludes identifiers and enforces manual validation", {
  dat <- data.frame(age = seq_len(40) + 20,
                    sex = rep(c("F", "M"), 20),
                    subject_id = sprintf("P%02d", 1:40),
                    notes = sprintf("text%02d", 1:40),
                    pain_t0 = 1:40, pain_t1 = 2:41)
  found <- MIRA:::.mira_detect_covariates(
    dat, exclude = c("pain_t0", "pain_t1"))
  expect_setequal(found$selected, c("age", "sex"))
  expect_false("subject_id" %in% found$detected)
  expect_false("notes" %in% found$detected)

  manual <- MIRA:::.mira_detect_covariates(dat, c("age", "sex"),
    exclude = c("pain_t0", "pain_t1"))
  expect_false(manual$automatic)
  expect_equal(manual$numeric, "age")
  expect_equal(manual$categorical, "sex")
  expect_error(MIRA:::.mira_detect_covariates(dat, c("age", "age")), "duplicate")
  expect_error(MIRA:::.mira_detect_covariates(dat, "absent"), "not found")
  expect_error(MIRA:::.mira_detect_covariates(dat, "pain_t0",
    exclude = "pain_t0"), "already used")
  dat$patient <- seq_len(nrow(dat))
  expect_error(MIRA:::.mira_detect_covariates(dat, "patient"), "reserved")
  expect_error(MIRA:::.mira_detect_covariates(dat, 1), "covariates must")
})

test_that("automatic covariates are withheld when adjustment degrees of freedom exceed limit", {
  dat <- data.frame(age = 1:20, weight = 21:40, height = 41:60,
                    sex = rep(c("F", "M"), 10))
  found <- MIRA:::.mira_detect_covariates(dat)
  expect_length(found$selected, 0L)
  expect_setequal(found$detected, names(dat))
  expect_match(found$warnings, "complexity exceeds")
  manual <- MIRA:::.mira_detect_covariates(dat, c("age", "weight", "height"))
  expect_equal(manual$selected, c("age", "weight", "height"))
})

test_that("direction and threshold resolution validate outcome-specific settings", {
  expect_equal(MIRA:::.mira_resolve_direction("auto", "BCVA")$value, "higher")
  expect_equal(MIRA:::.mira_resolve_direction("auto", "cmt")$value, "lower")
  expect_equal(MIRA:::.mira_resolve_direction("auto", "pain")$value, "unknown")
  expect_equal(MIRA:::.mira_resolve_direction("lower", "pain")$reason, "user")
  expect_equal(MIRA:::.mira_resolve_direction(
    c(Pain = "lower", BCVA = "higher"), "pain")$value, "lower")
  expect_error(MIRA:::.mira_resolve_direction("invalid", "pain"), "Allowed values")
  expect_error(MIRA:::.mira_resolve_direction(c("lower", "higher"), "pain"),
               "must be named")
  expect_error(MIRA:::.mira_resolve_direction(c(other = "lower"), "pain"),
               "does not contain exactly one")

  expect_equal(MIRA:::.mira_resolve_threshold("auto", "pain")$value, 0)
  expect_equal(MIRA:::.mira_resolve_threshold(NULL, "pain")$reason,
               "safe_zero_fallback")
  expect_equal(MIRA:::.mira_resolve_threshold(c(Pain = 1.5, BCVA = 2), "pain")$value,
               1.5)
  expect_error(MIRA:::.mira_resolve_threshold(-1, "pain"), "must be >= 0")
  expect_error(MIRA:::.mira_resolve_threshold(Inf, "pain"), "numeric value")
  expect_error(MIRA:::.mira_resolve_threshold(c(1, 2), "pain"), "must be named")
})

test_that("configuration builder records generated IDs, provenance, and column profiles", {
  dat <- data.frame(score_t0 = 1:12, score_t1 = 2:13,
                    .mira_subject_id = rep("occupied", 12),
                    check.names = FALSE)
  built <- MIRA:::.mira_build_analysis_config(dat)
  expect_true(built$config$id_generated)
  expect_equal(built$config$id, ".mira_subject_id_")
  expect_equal(built$data[[built$config$id]], seq_len(nrow(dat)))
  expect_equal(built$config$sources$id, "generated_row_id")
  expect_true(built$config$auto_detected$id)
  expect_equal(built$data_overview$n_rows, 12L)
  expect_equal(built$data_overview$n_columns, 3L)
  expect_setequal(built$data_overview$column_profile$variable, names(dat))
  expect_named(built$groups, "score")

  dat$patient_id <- sprintf("P%02d", 1:12)
  manual <- MIRA:::.mira_build_analysis_config(dat, id = "patient_id",
    improvement_direction = "higher", stable_threshold = 1)
  expect_false(manual$config$id_generated)
  expect_false(manual$config$auto_detected$id)
  expect_equal(manual$config$improvement_direction[["score"]], "higher")
  expect_equal(manual$config$stable_threshold[["score"]], 1)
})

test_that("configuration builder validates wide-format IDs and structural data", {
  dat <- data.frame(patient_id = c("p1", "p1", "p2", "p3"),
                    score_t0 = 1:4, score_t1 = 2:5)
  expect_error(MIRA:::.mira_build_analysis_config(dat, id = "patient_id"),
               "duplicate values")
  permissive <- MIRA:::.mira_build_analysis_config(dat, id = "patient_id",
                                                     strict_id = FALSE)
  expect_match(permissive$warnings, "duplicates")
  dat$patient_id[[4L]] <- NA_character_
  expect_error(MIRA:::.mira_build_analysis_config(dat, id = "patient_id"),
               "missing values")
  expect_error(MIRA:::.mira_build_analysis_config(data.frame()),
               "no observations")
  duplicate_names <- dat
  names(duplicate_names)[[3L]] <- "score_t0"
  expect_error(MIRA:::.mira_build_analysis_config(duplicate_names),
               "duplicate column names")
  expect_error(MIRA:::.mira_build_analysis_config(dat, strict_id = NA),
               "strict_id")
})

test_that("blank categorical covariates are excluded from advanced model data", {
  long <- data.frame(
    patient = c("p1", "p2"),
    time_label = c("Start", "Start"),
    time_index = c(1L, 1L),
    value = c(1, 2),
    sex = c("F", " ")
  )
  prepared <- MIRA:::.mira_prepare_advanced_data(
    long, "Start", FALSE, character(0), "sex", "sex")
  expect_equal(nrow(prepared), 1L)

  long$sex <- factor(long$sex)
  prepared_factor <- MIRA:::.mira_prepare_advanced_data(
    long, "Start", FALSE, character(0), "sex", "sex")
  expect_equal(nrow(prepared_factor), 1L)
})
