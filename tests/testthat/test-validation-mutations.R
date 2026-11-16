# Enable with MIRA_RUN_MUTATION_TESTS=true. MIRA_MUTATION_CASES accepts a comma-
# separated selection of the 14 cases listed by fixtures/run-mutation-validation.R
# --list. The default is the six-case focused longitudinal campaign.
test_that("extended mutation campaign classifies reached scientific defects", {
  skip_if_not(mira_validation_campaign_enabled("MIRA_RUN_MUTATION_TESTS"),
              "Set MIRA_RUN_MUTATION_TESTS=true for isolated mutation testing")
  repo <- mira_validation_campaign_root(require_sources = TRUE)
  worker <- file.path(repo, "tests", "testthat", "fixtures", "run-mutation-validation.R")
  expect_true(file.exists(worker), info = "Mutation worker must be available in the source tree")
  output <- mira_validation_campaign_output("mutations", repo)
  selection <- Sys.getenv("MIRA_MUTATION_CASES", unset = "")
  if (!nzchar(selection)) {
    selection <- paste(c("arm-sign", "missing-BH-family", "truncate-scientific-CSV",
      "false-render-success", "parser-ragged", "ignore-custom-chronology"), collapse = ",")
  }
  executable <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
  log <- file.path(output, "campaign.log")
  status <- system2(executable,
    c("--vanilla", shQuote(worker), shQuote(output), shQuote(selection), shQuote(repo)),
    stdout = log, stderr = log)
  message("Mutation evidence: ", output)
  expect_identical(status, 0L, info = paste("Inspect", log))
  summary_file <- file.path(output, "mutation-summary.csv")
  expect_true(file.exists(summary_file), info = paste("Inspect", log))
  if (!file.exists(summary_file)) stop("Mutation campaign did not produce a summary.")
  results <- utils::read.csv(summary_file, stringsAsFactors = FALSE)
  expect_identical(results$mutant, strsplit(selection, ",", fixed = TRUE)[[1L]])
  expect_true(all(results$classification %in% c("detected", "survived", "non-reached")),
              info = paste("Invalid or blocked campaigns require investigation; inspect", summary_file))
  # Survival/non-reach is a reported limitation, never evidence of detection.
  limitations <- results$classification != "detected"
  if (any(limitations)) {
    message("Mutation sensitivity limits: ", paste(
      paste0(results$mutant[limitations], "=", results$classification[limitations]),
      collapse = ", "))
  }
  print(results, row.names = FALSE)
})
