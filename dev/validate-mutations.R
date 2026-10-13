# Targeted sensitivity checks in isolated copies of the ENTIRE working tree,
# including local/untracked files and excluding .git. Never edits real sources.
# Run from repo root: Rscript --vanilla dev/validate-mutations.R <existing-output-dir>
# The companion runner is embedded below; only existing dependencies are used.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L || !dir.exists(args[[1L]])) stop("Supply one existing output directory.")
repo <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (!file.exists(file.path(repo, "DESCRIPTION"))) stop("Run from the MIRA repository root.")
output <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
Sys.setenv(MIRA_RUN_CMDSTAN_TESTS = "false", MIRA_RUN_QUARTO_TESTS = "false",
           MIRA_RUN_QUARTO_PDF_TESTS = "false",
           QUARTO_TEXLIVE_AUTO_INSTALL = "false", QUARTO_TEXLIVE_AUTO_UPDATE = "false")
if (startsWith(paste0(output, "/"), paste0(repo, "/"))) {
  stop("Use an output directory OUTSIDE the repository to avoid recursive copies.")
}
runner <- file.path(output, "mutation-runner.R")
writeLines(c(
  'args <- commandArgs(trailingOnly = TRUE)',
  'library(testthat)',
  'started <- proc.time()',
  'r <- test_local(args[[1]], filter = args[[2]], reporter = "silent", stop_on_failure = FALSE)',
  'tab <- as.data.frame(r)',
  'write.csv(tab[, !vapply(tab, is.list, logical(1)), drop=FALSE], args[[3]], row.names=FALSE)',
  'failures <- do.call(rbind, lapply(r, function(block) {',
  '  do.call(rbind, lapply(block$results, function(result) {',
  '    if (inherits(result, c("expectation_failure", "expectation_error")))',
  '      data.frame(test=block$test, class=class(result)[1L], message=result$message)',
  '  }))',
  '}))',
  'if (is.null(failures)) failures <- data.frame(test=character(), class=character(), message=character())',
  'write.csv(failures, sub("-tests.csv$", "-failures.csv", args[[3]]), row.names=FALSE)',
  'cat("blocks", nrow(tab), "failed", sum(tab$failed), "errors", sum(tab$error), "warnings", sum(tab$warning), "skips", sum(tab$skipped), "\\n")',
  'print(proc.time() - started)'
), runner)
mutant <- function(name, file, old, new, filter) list(name=name, file=file, old=old, new=new, filter=filter)
info <- "R/mira_info_long.R"
report <- "R/mira_report_freq_long.R"
cases <- list(
  mutant("change-sign", info, "delta <- y_complete - x_complete", "delta <- x_complete - y_complete", "statistical-oracles|sampling-inference|generated-properties"),
  mutant("wrong-subject-pairing", info, "y_complete <- y[keep]", "y_complete <- rev(y[keep])", "statistical-oracles|sampling-inference|generated-properties"),
  mutant("wrong-paired-N", info, "delta <- y_complete - x_complete\n      n <- length(delta)", "delta <- y_complete - x_complete\n      n <- length(delta) + 1L", "statistical-oracles|sampling-inference|generated-properties"),
  mutant("missing-BH-family", info, "change$paired_t_p_adj <- adjust_p_long(change$paired_t_p)", "change$paired_t_p_adj <- change$paired_t_p", "sampling-inference|generated-properties"),
  mutant("ignore-custom-chronology", info, "ord <- if (clear_order) order(parsed_order, time_vars) else seq_along(time_vars)", "ord <- seq_along(time_vars)", "configuration-contracts|properties"),
  mutant("retain-infinite-observations", info, "function(x) {\n      x[!is.finite(x)] <- NA_real_", "function(x) {\n      x[is.nan(x)] <- NA_real_", "boundaries-pathological|sampling-inference|generated-properties"),
  mutant("truncate-scientific-CSV", report, "utils::write.csv(entry$csv, entry$absolute", "utils::write.csv(utils::head(entry$csv, 1L), entry$absolute", "bundle-contracts|integration|files"),
  mutant("absolute-Wald-rank", info, "tolerance <- max(abs(decomposition$values)) * sqrt(.Machine$double.eps)", "tolerance <- max(abs(decomposition$values), 1) * sqrt(.Machine$double.eps)", "scientific-regressions"),
  mutant("GEE-covariate-hypothesis", info, ".mira_effect_indices_long(names(beta_eval$value), question, model = model)", ".mira_effect_indices_long(names(beta_eval$value), question)", "scientific-regressions"),
  mutant("Stan-observation-pairing", "R/mira_data_long.R", "y = as.numeric(y)", "y = rev(as.numeric(y))", "data-prior-long-oracles"),
  mutant("false-render-success", report, "rendered[[fmt]] <- is.null(error_message) && file.exists(output_paths[[fmt]])", "rendered[[fmt]] <- file.exists(output_paths[[fmt]])", "validation-regressions|rendering-contracts"),
  mutant("archive-section-loss", report, "  cfg$sections <- \"all\"", "  cfg$sections <- \"descriptives\"", "bundle-contracts|files|integration")
)
copy_tree <- function(destination) {
  dir.create(destination)
  paths <- list.files(repo, all.files=TRUE, full.names=TRUE, no..=TRUE)
  paths <- paths[basename(paths) != ".git"]
  if (!all(file.copy(paths, destination, recursive=TRUE, copy.mode=TRUE))) stop("Working tree copy failed.")
  if (dir.exists(file.path(destination, ".git"))) stop("Git directory copied unexpectedly.")
}
results <- vector("list", length(cases))
for (i in seq_along(cases)) {
  case <- cases[[i]]
  destination <- file.path(output, case$name)
  if (file.exists(destination)) stop("Output directory already contains mutant: ", case$name)
  copy_tree(destination)
  target <- file.path(destination, case$file)
  source <- paste(readLines(target, warn=FALSE, encoding="UTF-8"), collapse="\n")
  hits <- gregexpr(case$old, source, fixed=TRUE)[[1L]]
  if (length(hits) != 1L || hits[[1L]] < 0L) stop("Mutation target not unique: ", case$name)
  source <- sub(case$old, case$new, source, fixed=TRUE)
  writeLines(source, target, useBytes=TRUE)
  csv <- file.path(output, paste0(case$name, "-tests.csv"))
  log <- file.path(output, paste0(case$name, ".log"))
  executable <- if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript"
  status <- system2(file.path(R.home("bin"), executable),
    c("--vanilla", shQuote(runner), shQuote(destination), shQuote(case$filter), shQuote(csv)),
    stdout=log, stderr=log)
  if (status != 0L || !file.exists(csv)) stop("Infrastructure error for ", case$name, "; inspect ", log)
  tab <- utils::read.csv(csv, stringsAsFactors=FALSE)
  failures <- utils::read.csv(sub("-tests.csv$", "-failures.csv", csv), stringsAsFactors=FALSE)
  detected <- any(failures$class == "expectation_failure")
  results[[i]] <- data.frame(mutant=case$name, file=case$file, filter=case$filter,
    assertion_failures=sum(tab$failed), test_errors=sum(tab$error), warnings=sum(tab$warning),
    skips=sum(tab$skipped), detected=detected,
    interpretation=if (detected) "Relevant assertions failed; inspect failure evidence" else "SURVIVED: investigate equivalence/path/test gap")
  print(results[[i]], row.names=FALSE)
}
utils::write.csv(do.call(rbind, results), file.path(output, "mutation-summary.csv"), row.names=FALSE)
