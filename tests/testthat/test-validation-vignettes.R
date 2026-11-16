# The complete profile needs a source checkout as well as an installed MIRA.
# Set MIRA_VALIDATION_REPO_DIR when tests run against an installed/check copy.
validation_vignettes_repo <- function() {
  valid <- function(path) {
    description <- file.path(path, "DESCRIPTION")
    if (!file.exists(description) ||
        !file.exists(file.path(path, "R", "mira_info_long.R")) ||
        !file.exists(file.path(path, "vignettes", "MIRA.Rmd"))) return(FALSE)
    metadata <- tryCatch(read.dcf(description), error = function(e) NULL)
    !is.null(metadata) && "Package" %in% colnames(metadata) &&
      identical(unname(metadata[1L, "Package"]), "MIRA")
  }
  explicit <- Sys.getenv("MIRA_VALIDATION_REPO_DIR", "")
  if (nzchar(explicit)) {
    if (!valid(explicit)) stop("MIRA_VALIDATION_REPO_DIR must identify a MIRA source checkout.",
      call. = FALSE)
    return(normalizePath(explicit, winslash = "/", mustWork = TRUE))
  }
  candidates <- unique(c(testthat::test_path("..", ".."), getwd()))
  available <- candidates[vapply(candidates, valid, logical(1L))]
  if (!length(available)) return(NULL)
  normalizePath(available[[1L]], winslash = "/", mustWork = TRUE)
}

test_that("the frequentist vignette index lists every source article", {
  repo <- validation_vignettes_repo()
  if (is.null(repo)) skip("MIRA vignette sources are unavailable in this installed/check copy.")
  index <- readLines(file.path(repo, "vignettes", "MIRA.Rmd"),
    warn = FALSE, encoding = "UTF-8")
  links <- unlist(regmatches(index,
    gregexpr("frequentist-long-[A-Za-z0-9_-]+[.]html", index)))
  expected <- c("MIRA", unique(sub("[.]html$", "", links)))
  available <- tools::file_path_sans_ext(list.files(file.path(repo, "vignettes"),
    pattern = "^(MIRA|frequentist-long-.+)[.]Rmd$"))
  expect_setequal(expected, available)
})

test_that("the installed package executes the complete frequentist vignette collection", {
  if (!identical(Sys.getenv("MIRA_RUN_VIGNETTE_TESTS"), "true")) {
    skip("Set MIRA_RUN_VIGNETTE_TESTS=true for isolated rendering and both code modes.")
  }
  repo <- validation_vignettes_repo()
  if (is.null(repo)) stop("MIRA_RUN_VIGNETTE_TESTS requires a source checkout; set MIRA_VALIDATION_REPO_DIR.",
    call. = FALSE)
  runner <- file.path(repo, "tests", "testthat", "fixtures", "run-vignette-validation.R")
  if (!file.exists(runner)) stop("The vignette validation fixture is missing from the source checkout.",
    call. = FALSE)
  base <- Sys.getenv("MIRA_VALIDATION_OUTPUT_DIR", "")
  if (!nzchar(base)) base <- tempdir()
  # The runner rejects paths inside the checkout before creating directories.
  # Keep its output fresh and its failure log outside that guarded location.
  output <- tempfile(pattern = "vignettes-", tmpdir = path.expand(base))
  log <- withr::local_tempfile(pattern = "mira-vignettes-", fileext = ".log")
  executable <- file.path(R.home("bin"),
    if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
  # A fresh R process must see the same installed package libraries as the
  # calling suite, including an explicitly selected validation library.
  withr::local_envvar(R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep))
  withr::local_dir(repo)
  status <- system2(executable,
    c("--vanilla", shQuote(runner), shQuote(output)), stdout = log, stderr = log)
  diagnostics <- paste(readLines(log, warn = FALSE), collapse = "\n")
  expect_identical(as.integer(status), 0L, info = diagnostics)
  status_file <- file.path(output, "validation", "vignette-status.csv")
  expect_true(file.exists(status_file), info = diagnostics)
  if (file.exists(status_file)) {
    completed <- utils::read.csv(status_file, stringsAsFactors = FALSE)
    expect_true(all(completed$success), info = diagnostics)
    expect_setequal(unique(completed$mode), c("render", "extracted", "visible"))
  }
})
