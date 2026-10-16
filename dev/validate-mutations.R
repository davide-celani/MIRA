# Targeted sensitivity checks on isolated copies of the current working tree.
# Rscript --vanilla dev/validate-mutations.R <existing-output-dir> [name1,name2,...]
# Rscript --vanilla dev/validate-mutations.R --list
# No dependencies are installed; .git is never copied or changed.
args <- commandArgs(trailingOnly = TRUE)
repo <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (!file.exists(file.path(repo, "DESCRIPTION"))) stop("Run from the MIRA repository root.")
Sys.setenv(MIRA_RUN_CMDSTAN_TESTS = "false", MIRA_RUN_QUARTO_TESTS = "false",
           MIRA_RUN_QUARTO_PDF_TESTS = "false",
           QUARTO_TEXLIVE_AUTO_INSTALL = "false", QUARTO_TEXLIVE_AUTO_UPDATE = "false")

# The marker is silent and external to the copied repository. It records the
# altered expression, not merely entry to its containing function.
marker <- 'base::cat("reached\\n", file = base::Sys.getenv("MIRA_MUTATION_MARKER"), append = TRUE)'
marked <- function(expression) paste0("{ ", marker, "; ", expression, " }")
mutant <- function(name, file, old, new, filter, relevant_test) {
  list(name = name, file = file, old = old, new = new,
       filter = filter, relevant_test = relevant_test)
}
info <- "R/mira_info_long.R"
report <- "R/mira_report_freq_long.R"
cases <- list(
  mutant("change-sign", info, "delta <- y_complete - x_complete",
    paste0("delta <- ", marked("x_complete - y_complete")),
    "statistical-oracles|sampling-inference|generated-properties",
    "pairwise change estimates|bounded finite-pair cases"),
  mutant("wrong-subject-pairing", info, "y_complete <- y[keep]",
    paste0("y_complete <- ", marked("rev(y[keep])")),
    "statistical-oracles|sampling-inference|generated-properties",
    "pairwise change estimates|bounded finite-pair cases"),
  mutant("wrong-paired-N", info,
    "delta <- y_complete - x_complete\n      n <- length(delta)",
    paste0("delta <- y_complete - x_complete\n      n <- ", marked("length(delta) + 1L")),
    "statistical-oracles|sampling-inference|generated-properties",
    "pairwise change estimates|bounded finite-pair cases"),
  mutant("missing-BH-family", info,
    "change$paired_t_p_adj <- adjust_p_long(change$paired_t_p)",
    paste0("change$paired_t_p_adj <- ", marked("change$paired_t_p")),
    "sampling-inference|generated-properties", "bounded finite-pair cases"),
  mutant("ignore-custom-chronology", info,
    "ord <- if (clear_order) order(parsed_order, time_vars) else seq_along(time_vars)",
    paste0("ord <- if (clear_order) ", marked("seq_along(time_vars)"),
           " else seq_along(time_vars)"),
    "configuration-contracts|properties", "custom chronology determines baseline and final"),
  mutant("retain-infinite-observations", info,
    "function(x) {\n      x[!is.finite(x)] <- NA_real_",
    paste0("function(x) {\n      ", marker, "\n      x[is.nan(x)] <- NA_real_"),
    "boundaries-pathological|sampling-inference|generated-properties",
    "mixed missing profiles|bounded finite-pair cases|singleton missing and nonfinite"),
  mutant("truncate-scientific-CSV", report,
    "utils::write.csv(entry$csv, entry$absolute",
    paste0("utils::write.csv(", marked("utils::head(entry$csv, 1L)"), ", entry$absolute"),
    "bundle-contracts|files", "display row limit never truncates|display settings and selected sections"),
  mutant("absolute-Wald-rank", info,
    "tolerance <- max(abs(decomposition$values)) * sqrt(.Machine$double.eps)",
    paste0("tolerance <- ", marked("max(abs(decomposition$values), 1) * sqrt(.Machine$double.eps)")),
    "scientific-regressions", "Wald rank and inference|rescaling an observed outcome"),
  mutant("GEE-covariate-hypothesis", info,
    ".mira_effect_indices_long(names(beta_eval$value), question, model = model)",
    marked(".mira_effect_indices_long(names(beta_eval$value), question)"),
    "scientific-regressions", "GEE scientific hypotheses are unchanged"),
  mutant("Stan-observation-pairing", "R/mira_data_long.R", "y = as.numeric(y)",
    paste0("y = ", marked("rev(as.numeric(y))")),
    "data-prior-long-oracles", "Stan observations retain subject, visit, arm"),
  mutant("false-render-success", report,
    "rendered[[fmt]] <- is.null(error_message) && file.exists(output_paths[[fmt]])",
    paste0("rendered[[fmt]] <- ", marked("file.exists(output_paths[[fmt]])")),
    "validation-regressions|rendering-contracts", "old rendered files cannot make a no-op renderer"),
  mutant("archive-section-loss", report, '  cfg$sections <- "all"',
    paste0('  cfg$sections <- ', marked('"descriptives"')),
    "bundle-contracts|files", "display settings and selected sections preserve every archived"),
  mutant("arm-sign", info, "estimate = mean(y) - mean(x),",
    paste0("estimate = ", marked("mean(x) - mean(y)"), ","),
    "audit-2026-10-16", "^asymmetric arm samples have independent Welch, Hedges and finite-family Holm oracles$"),
  mutant("parser-ragged", "tests/testthat/helper-mira-report-freq_long.R",
    "if (any(widths != length(headers))) {", paste0("if (", marked("FALSE"), ") {"),
    "html-parser", "^published HTML rejects malformed rows before recycling scientific cells$")
)
names(cases) <- vapply(cases, `[[`, character(1L), "name")
if (identical(args, "--list")) {
  cat(paste(names(cases), collapse = "\n"), "\n")
  quit(status = 0L)
}
if (length(args) < 1L || length(args) > 2L || !dir.exists(args[[1L]])) {
  stop("Supply an existing external output directory and optional comma-separated mutant names.")
}
# Default is the focused longitudinal campaign. Historical cases, including
# Stan, remain selectable but are not part of this scoped default profile.
selection <- if (length(args) == 2L) strsplit(args[[2L]], ",", fixed = TRUE)[[1L]] else {
  c("arm-sign", "missing-BH-family", "truncate-scientific-CSV",
    "false-render-success", "parser-ragged", "ignore-custom-chronology")
}
if (!length(selection) || anyDuplicated(selection) || any(!selection %in% names(cases))) {
  stop("Mutant names must be unique and listed by --list.")
}
cases <- cases[selection]
output <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
if (startsWith(paste0(output, "/"), paste0(repo, "/"))) {
  stop("Use an output directory OUTSIDE the repository to avoid recursive copies.")
}
runner <- file.path(output, "mutation-runner.R")
baseline <- file.path(output, "unmutated-working-tree")
summary_file <- file.path(output, "mutation-summary.csv")
planned_files <- unlist(lapply(c("unmutated", selection), function(name) {
  file.path(output, paste0(name, c(".log", "-tests.csv", "-failures.csv", "-reached.txt")))
}), use.names = FALSE)
if (any(file.exists(c(runner, baseline, summary_file, planned_files, file.path(output, selection))))) {
  stop("Use a fresh output directory; previous evidence will not be overwritten.")
}
writeLines(c(
  'args <- commandArgs(trailingOnly = TRUE)',
  'library(testthat)',
  'started <- proc.time()',
  'r <- test_local(args[[1]], filter = args[[2]], reporter = "silent", stop_on_failure = FALSE)',
  'tab <- as.data.frame(r)',
  'write.csv(tab[, !vapply(tab, is.list, logical(1)), drop = FALSE], args[[3]], row.names = FALSE)',
  'failures <- do.call(rbind, lapply(r, function(block) {',
  '  do.call(rbind, lapply(block$results, function(result) {',
  '    if (inherits(result, c("expectation_failure", "expectation_error")))',
  '      data.frame(test = block$test, class = class(result)[1L], message = result$message)',
  '  }))',
  '}))',
  'if (is.null(failures)) failures <- data.frame(test = character(), class = character(), message = character())',
  'write.csv(failures, sub("-tests.csv$", "-failures.csv", args[[3]]), row.names = FALSE)',
  'cat("blocks", nrow(tab), "failed", sum(tab$failed), "errors", sum(tab$error), "warnings", sum(tab$warning), "skips", sum(tab$skipped), "\\n")',
  'print(proc.time() - started)'
), runner)
copy_tree <- function(destination) {
  dir.create(destination)
  paths <- list.files(repo, all.files = TRUE, full.names = TRUE, no.. = TRUE)
  paths <- paths[basename(paths) != ".git"]
  if (!all(file.copy(paths, destination, recursive = TRUE, copy.mode = TRUE))) stop("Working tree copy failed.")
  if (dir.exists(file.path(destination, ".git"))) stop("Git directory copied unexpectedly.")
}
run_tests <- function(tree, filter, name) {
  csv <- file.path(output, paste0(name, "-tests.csv"))
  log <- file.path(output, paste0(name, ".log"))
  executable <- if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript"
  status <- system2(file.path(R.home("bin"), executable),
    c("--vanilla", shQuote(runner), shQuote(tree), shQuote(filter), shQuote(csv)),
    stdout = log, stderr = log)
  if (status != 0L || !file.exists(csv)) return(list(blocked = TRUE, log = log))
  list(blocked = FALSE, log = log, tab = utils::read.csv(csv, stringsAsFactors = FALSE),
    failures = utils::read.csv(sub("-tests.csv$", "-failures.csv", csv), stringsAsFactors = FALSE))
}
record <- function(case, classification, interpretation, run = NULL, reached = FALSE, relevant = 0L) {
  data.frame(mutant = case$name, file = case$file, filter = case$filter,
    relevant_test = case$relevant_test, reached = reached,
    assertion_failures = if (is.null(run$tab)) NA_integer_ else sum(run$tab$failed),
    relevant_failures = relevant,
    test_errors = if (is.null(run$tab)) NA_integer_ else sum(run$tab$error),
    warnings = if (is.null(run$tab)) NA_integer_ else sum(run$tab$warning),
    skips = if (is.null(run$tab)) NA_integer_ else sum(run$tab$skipped),
    detected = classification == "detected", classification = classification,
    interpretation = interpretation)
}

# A passing source suite elsewhere does not replace this copy-level baseline.
copy_tree(baseline)
union_filter <- paste(unique(vapply(cases, `[[`, character(1L), "filter")), collapse = "|")
Sys.unsetenv("MIRA_MUTATION_MARKER")
base_run <- run_tests(baseline, union_filter, "unmutated")
baseline_ok <- !base_run$blocked && nrow(base_run$tab) > 0L &&
  all(base_run$tab$failed == 0L & !base_run$tab$error & base_run$tab$warning == 0L)
if (baseline_ok) {
  baseline_ok <- all(vapply(cases, function(case) {
    any(grepl(case$relevant_test, base_run$tab$test, perl = TRUE))
  }, logical(1L)))
}
if (!baseline_ok) {
  results <- lapply(cases, record, classification = "blocked",
    interpretation = "Unmutated copied baseline failed or lacks a designated relevant test; inspect unmutated.log and CSV.")
  utils::write.csv(do.call(rbind, results), summary_file, row.names = FALSE)
  stop("Unmutated copied baseline did not pass; no mutants were run.")
}
results <- vector("list", length(cases))
for (i in seq_along(cases)) {
  case <- cases[[i]]
  destination <- file.path(output, case$name)
  copy_tree(destination)
  target <- file.path(destination, case$file)
  source <- paste(readLines(target, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  hits <- gregexpr(case$old, source, fixed = TRUE)[[1L]]
  if (length(hits) != 1L || hits[[1L]] < 0L) {
    results[[i]] <- record(case, "invalid", "Mutation target is absent or not unique.")
    next
  }
  writeLines(sub(case$old, case$new, source, fixed = TRUE), target, useBytes = TRUE)
  parse_error <- tryCatch({ parse(file = target); NULL }, error = conditionMessage)
  if (!is.null(parse_error)) {
    results[[i]] <- record(case, "invalid", paste("Mutated file does not parse:", parse_error))
    next
  }
  marker_file <- file.path(output, paste0(case$name, "-reached.txt"))
  Sys.setenv(MIRA_MUTATION_MARKER = marker_file)
  run <- run_tests(destination, case$filter, case$name)
  reached <- file.exists(marker_file) && file.info(marker_file)$size > 0L
  if (run$blocked) {
    results[[i]] <- record(case, "blocked", paste("Load/process infrastructure failed; inspect", run$log), run, reached)
  } else {
    relevant <- sum(run$failures$class == "expectation_failure" &
      grepl(case$relevant_test, run$failures$test, perl = TRUE))
    if (reached && relevant > 0L) {
      results[[i]] <- record(case, "detected", "Altered expression reached and a predeclared relevant assertion failed; inspect failure CSV.", run, reached, relevant)
    } else if (any(run$tab$error | run$tab$warning > 0L | run$tab$failed > 0L)) {
      results[[i]] <- record(case, "blocked", "Only errors/warnings or unrelated assertions failed; no relevant detection established.", run, reached, relevant)
    } else {
      results[[i]] <- record(case, if (reached) "survived" else "non-reached",
        if (reached) "Altered expression reached without relevant assertion failure; investigate the gap."
        else "Altered expression was not reached; no sensitivity conclusion.", run, reached, relevant)
    }
  }
  print(results[[i]], row.names = FALSE)
}
utils::write.csv(do.call(rbind, results), summary_file, row.names = FALSE)
if (any(vapply(results, function(x) x$classification %in% c("invalid", "blocked"), logical(1L)))) {
  quit(status = 1L)
}
