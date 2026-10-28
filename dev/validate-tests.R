main <- function() {
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L || length(args) > 3L) {
  stop("Usage: validate-tests.R <repository> <existing-directory/output-prefix> [source|installed]")
}
repo <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
prefix <- args[[2]]
mode <- if (length(args) >= 3L) args[[3]] else "source"
if (!mode %in% c("source", "installed")) stop("Mode must be source or installed.")
if (!dir.exists(dirname(prefix))) stop("The output-prefix directory must already exist.")
# Keep development and installed-namespace validation on the package edition,
# including when this runner is launched without a testthat package context.
testthat::local_edition(3)
cat("Requested testthat edition: 3\n")
Sys.setenv(QUARTO_TEXLIVE_AUTO_INSTALL = "false",
           QUARTO_TEXLIVE_AUTO_UPDATE = "false",
           MIRA_RUN_CMDSTAN_TESTS = "false",
           MIRA_RUN_QUARTO_TESTS = "false",
           MIRA_RUN_QUARTO_PDF_TESTS = "false")
cat("Mode:", mode, "\nRepository:", repo, "\n")
cat("Started:", format(Sys.time(), tz = "Europe/Rome"), "\n")
print(sessionInfo())
ip <- installed.packages()
write.csv(ip[, c("Package", "Version", "LibPath")],
          paste0(prefix, "-packages.csv"), row.names = FALSE)
started <- proc.time()
if (mode == "source") {
  pkgload::load_all(repo, quiet = TRUE)
} else {
  library(MIRA)
  cat("Installed package path:", getNamespaceInfo(asNamespace("MIRA"), "path"), "\n")
  cat("Public exports:", paste(sort(getNamespaceExports("MIRA")), collapse = ", "), "\n")
}
set.seed(20261012)
snapshot <- function() list(seed = .Random.seed, rng = RNGkind(), wd = getwd(),
  devices = grDevices::dev.list(), options = options(),
  environment = Sys.getenv(c("TZ", "LANG", "LC_ALL", "CMDSTAN", "QUARTO_PATH"), unset = NA))
before <- snapshot()
reporter <- testthat::SummaryReporter$new(show_praise = FALSE)
if (mode == "source") {
  result <- testthat::test_local(repo, reporter = reporter, stop_on_failure = FALSE)
} else {
  result <- testthat::test_dir(file.path(repo, "tests", "testthat"),
    package = "MIRA", load_package = "installed", reporter = reporter,
    stop_on_failure = FALSE)
}
after <- snapshot()
changed_options <- names(before$options)[!vapply(names(before$options), function(name) {
  identical(before$options[[name]], after$options[[name]])
}, logical(1))]
state <- list(seed_restored = identical(before$seed, after$seed),
  rng_restored = identical(before$rng, after$rng),
  working_directory_restored = identical(before$wd, after$wd),
  devices_restored = identical(before$devices, after$devices),
  environment_restored = identical(before$environment, after$environment),
  changed_options = changed_options, new_options = setdiff(names(after$options), names(before$options)))
writeLines(capture.output(str(state)), paste0(prefix, "-state.txt"))
print(state)
elapsed <- proc.time() - started
saveRDS(result, paste0(prefix, "-results.rds"))
tab <- as.data.frame(result)
write.csv(tab[, !vapply(tab, is.list, logical(1)), drop = FALSE],
          paste0(prefix, "-results.csv"), row.names = FALSE)
cat("\nTotals:\n")
print(colSums(tab[, intersect(c("nb", "failed", "skipped", "error", "warning", "passed"), names(tab)), drop = FALSE]))
cat("Test blocks:", nrow(tab), "\n")
print(elapsed)
writeLines(capture.output(print(elapsed)), paste0(prefix, "-elapsed.txt"))
print(sessionInfo())
# First loading Imports may add documented package defaults; changing an
# existing option or any other recorded state is an unexpected leak.
state_failed <- !all(unlist(state[1:5])) || length(state$changed_options) > 0L
as.integer(any(tab$failed > 0 | tab$error | tab$warning > 0) || state_failed)
}
# Run deferred edition cleanup before quitting, while its namespace is alive.
quit(status = main())
