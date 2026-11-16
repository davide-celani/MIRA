# Infrastructure for explicitly enabled, extended validation profiles.
# Evidence is always written outside the package source directory.
mira_validation_campaign_enabled <- function(variable) {
  tolower(Sys.getenv(variable, unset = "false")) %in% c("true", "1", "yes")
}

mira_validation_campaign_root <- function(require_sources = FALSE) {
  candidate <- Sys.getenv("MIRA_VALIDATION_REPO_DIR", unset = "")
  if (!nzchar(candidate)) candidate <- testthat::test_path("..", "..")
  root <- normalizePath(candidate, winslash = "/", mustWork = TRUE)
  if (require_sources && !all(file.exists(file.path(root, c(
    "DESCRIPTION", "R/mira_info_long.R", "R/mira_report_freq_long.R",
    "tests/testthat/helper-mira-report-freq_long.R"
  ))))) {
    stop("Mutation validation requires the MIRA source repository; set MIRA_VALIDATION_REPO_DIR.",
         call. = FALSE)
  }
  root
}

mira_validation_campaign_output <- function(profile, repo) {
  base <- Sys.getenv("MIRA_VALIDATION_OUTPUT_DIR", unset = "")
  if (!nzchar(base)) base <- tempdir()
  if (!dir.exists(base)) {
    stop("MIRA_VALIDATION_OUTPUT_DIR must name an existing external directory.", call. = FALSE)
  }
  base <- normalizePath(base, winslash = "/", mustWork = TRUE)
  key <- function(path) if (.Platform$OS.type == "windows") tolower(path) else path
  if (startsWith(paste0(key(base), "/"), paste0(key(repo), "/"))) {
    stop("Validation output must be outside the source repository.", call. = FALSE)
  }
  output <- tempfile(paste0(profile, "-"), tmpdir = base)
  if (!dir.create(output)) stop("Cannot create validation output directory.", call. = FALSE)
  normalizePath(output, winslash = "/", mustWork = TRUE)
}
