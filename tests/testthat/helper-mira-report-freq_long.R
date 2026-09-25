report_test_cfg <- function(sections = "all", ...) {
  defaults <- list(
    sections = sections,
    include_complete_output = TRUE,
    max_table_rows = Inf,
    max_table_columns = 12L,
    digits = 3L,
    max_depth = 50L
  )
  utils::modifyList(defaults, list(...))
}

report_test_mira <- function(outcome = "score") {
  pearson <- matrix(
    c(1, 0.5, 0.5, 1), nrow = 2L,
    dimnames = list(c("t0", "t1"), c("t0", "t1"))
  )
  result <- list(
    version = "4.0.0",
    config = list(
      id = "patient",
      outcomes = outcome,
      analyses = list(
        core = c("descriptives", "missingness", "change", "variability",
                 "trajectories"),
        plots = FALSE,
        model = FALSE,
        outliers = FALSE,
        correlations = TRUE,
        arm_tests = FALSE
      )
    ),
    data_overview = data.frame(
      records = 6L, variables = 3L, stringsAsFactors = FALSE
    ),
    detected_variables = list(id = "patient", outcomes = outcome),
    outcome = outcome,
    outcome_display = paste("Outcome", outcome),
    settings = list(alpha = 0.05, p_adjust_method = "holm"),
    overview = data.frame(
      metric = c("participants", "timepoints"), value = c(6, 2),
      stringsAsFactors = FALSE
    ),
    descriptives = data.frame(
      time = c("t0", "t1"), n = c(6L, 6L), mean = c(10.25, 11.75),
      sd = c(1.5, 1.25), p.value = c(0.0004, 0.05),
      check.names = FALSE, stringsAsFactors = FALSE
    ),
    missing = list(
      by_time = data.frame(
        time = c("t0", "t1"), n_missing = c(0L, 1L),
        percent_missing = c(0, 100 / 6), stringsAsFactors = FALSE
      )
    ),
    change = data.frame(
      from = "t0", to = "t1", n = 5L, mean_change = 1.5,
      p_value = 0.012345, check.names = FALSE, stringsAsFactors = FALSE
    ),
    correlations = list(
      pearson = pearson,
      spearman = pearson - matrix(c(0, 0.1, 0.1, 0), 2L),
      pairwise_n = matrix(c(6L, 5L, 5L, 5L), 2L,
                          dimnames = dimnames(pearson))
    ),
    variability = list(
      summary = data.frame(component = c("within", "between"),
                           variance = c(1.25, 2.5), stringsAsFactors = FALSE)
    ),
    trajectories = data.frame(
      patient = paste0("P", 1:3), baseline = c(9, 10, 11),
      final = c(10, 12, 12), change = c(1, 2, 1),
      stringsAsFactors = FALSE
    )
  )
  structure(result, class = c("mira_info_long", "list"))
}

report_test_multi <- function(include_error = FALSE) {
  score <- report_test_mira("score")
  fatigue <- report_test_mira("fatigue")
  fatigue$descriptives$mean <- c(20.5, 18.25)
  outcomes <- list(score = score, fatigue = fatigue)
  if (include_error) {
    outcomes$failed <- structure(
      list(outcome = "failed", error = "controlled outcome failure"),
      class = c("mira_info_error_long", "list")
    )
  }
  structure(
    list(
      version = "4.0.0",
      config = list(outcomes = names(outcomes)),
      data_overview = data.frame(records = 6L, variables = 5L),
      detected_variables = list(outcomes = names(outcomes)),
      outcomes = outcomes,
      diagnostics = list(
        warnings = character(),
        failed_outcomes = if (include_error) "failed" else character()
      )
    ),
    class = c("mira_info_multi_long", "mira_info_long", "list")
  )
}

report_test_detect <- function() {
  structure(
    list(
      version = "4.0.0",
      config = list(
        id = "patient",
        outcomes = "score",
        time_vars = c("score_t0", "score_t1"),
        analyses = list(core = TRUE, model = FALSE)
      ),
      data_overview = data.frame(records = 4L, variables = 3L),
      detected_variables = list(id = "patient", outcomes = "score"),
      diagnostics = list(warnings = character())
    ),
    class = c("mira_detect_long", "list")
  )
}

report_test_run <- function(..., .capture = TRUE) {
  if (!.capture) return(mira_report_freq_long(...))
  printed <- utils::capture.output(value <- mira_report_freq_long(...))
  attr(value, "test_printed_output") <- printed
  value
}

report_test_call <- function(args = list(), output_dir = tempfile("mira-report-test-")) {
  call_args <- list(
    x = report_test_mira(),
    output_dir = output_dir,
    format = "html",
    render = FALSE,
    open = FALSE
  )
  call_args[names(args)] <- args
  printed <- utils::capture.output(value <- do.call(mira_report_freq_long, call_args))
  attr(value, "test_printed_output") <- printed
  value
}

report_test_atomic_table <- function(x) {
  is.data.frame(x) && all(vapply(x, function(column) {
    is.atomic(column) && is.null(dim(column)) && length(column) == nrow(x)
  }, logical(1L)))
}

report_test_relpath <- function(path) {
  gsub("\\\\", "/", path)
}
