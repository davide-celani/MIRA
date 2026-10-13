# Installed-package real Quarto profile. All software must already exist.
# Rscript --vanilla dev/validate-rendering.R <repository> <existing-output-directory>
# Set MIRA_RUN_QUARTO_PDF_TESTS=true explicitly to include PDF with existing TeX.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || !dir.exists(args[[2L]])) {
  stop("Supply repository and one existing output directory.")
}
repo <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
output <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
Sys.setenv(MIRA_RUN_CMDSTAN_TESTS = "false", MIRA_RUN_QUARTO_TESTS = "true",
           QUARTO_TEXLIVE_AUTO_INSTALL = "false", QUARTO_TEXLIVE_AUTO_UPDATE = "false")
library(MIRA)
cat("Installed package:", getNamespaceInfo(asNamespace("MIRA"), "path"), "\n")
cat("Quarto:", as.character(quarto::quarto_version()), "\n")
print(Sys.which(c("pdflatex", "xelatex", "lualatex")))
print(sessionInfo())
result <- testthat::test_dir(file.path(repo, "tests", "testthat"),
  filter = "mira-report-freq_long-(integration|rendering-contracts)",
  package = "MIRA", load_package = "installed",
  reporter = testthat::SummaryReporter$new(show_praise = FALSE),
  stop_on_failure = FALSE)
saveRDS(result, file.path(output, "quarto-results.rds"))
tab <- as.data.frame(result)
write.csv(tab[, !vapply(tab, is.list, logical(1)), drop = FALSE],
          file.path(output, "quarto-results.csv"), row.names = FALSE)
print(colSums(tab[, c("nb", "failed", "skipped", "error", "warning", "passed")]))
if (any(tab$failed > 0 | tab$error | tab$warning > 0)) quit(status = 1L)
