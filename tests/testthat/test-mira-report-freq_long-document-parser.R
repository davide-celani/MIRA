report_docx_parser_fixture <- function(headers = c("time", "n", "mean"),
                                       rows = list(c("Baseline", "6", "3.500"))) {
  row <- function(cells) paste0('<w:tr>', paste0('<w:tc><w:p><w:r><w:t>',
    cells, '</w:t></w:r></w:p></w:tc>', collapse = ""), '</w:tr>')
  paste0('<w:tbl><w:tblPr/><w:tblGrid/>', row(headers),
    paste0(vapply(rows, row, character(1L)), collapse = ""), '</w:tbl>')
}

test_that("DOCX verifier selects a unique scientific header signature and decodes cells once", {
  target <- report_docx_parser_fixture(rows = list(
    c("A &amp;lt; &#956;", "0", "-0.0010"), c("B", "NA", "")))
  # Text split across Word runs must preserve the original cell text.
  target <- sub('<w:t>B</w:t>', '<w:t>B</w:t></w:r><w:r><w:br/><w:t>C</w:t>',
    target, fixed = TRUE)
  decoy <- report_docx_parser_fixture(c("time", "other"), list(c("X", "999")))
  observed <- report_test_docx_table(paste0(decoy, target, decoy), c("time", "n", "mean"))
  expect_identical(names(observed), c("time", "n", "mean"))
  expect_identical(observed$time, c("A &lt; \u03bc", "B C"))
  expect_identical(observed$n, c("0", "NA"))
  expect_identical(observed$mean, c("-0.0010", ""))
})

test_that("DOCX verifier identifies split scientific blocks by their repeated visit keys", {
  counts <- report_docx_parser_fixture(c("time", "label", "n"),
    list(c("score_t0", "Baseline", "6"), c("score_t1", "Final", "5")))
  means <- report_docx_parser_fixture(c("time", "label", "mean"),
    list(c("score_t0", "Baseline", "3.500"), c("score_t1", "Final", "4.500")))
  xml <- paste0(counts, means)
  n_block <- report_test_docx_table(xml, c("time", "n"))
  mean_block <- report_test_docx_table(xml, c("time", "mean"))
  expect_identical(n_block$time, c("score_t0", "score_t1"))
  expect_identical(mean_block$time, n_block$time)
  expect_identical(n_block$n, c("6", "5"))
  expect_identical(mean_block$mean, c("3.500", "4.500"))
  expect_error(report_test_docx_table(xml, c("time", "n", "mean")),
    "expected exactly one table with headers [time, n, mean]; found 0", fixed = TRUE)
})

test_that("DOCX verifier diagnoses absent, duplicate and malformed scientific tables", {
  valid <- report_docx_parser_fixture(rows = list(c("0", "0", "0"), c("B", "6", "3.5")))
  expect_error(report_test_docx_table(valid, c("time", "estimate")),
    "expected exactly one table with headers [time, estimate]; found 0", fixed = TRUE)
  expect_error(report_test_docx_table(paste0(valid, valid), c("time", "mean")),
    "expected exactly one table with headers [time, mean]; found 2", fixed = TRUE)
  ragged <- report_docx_parser_fixture(rows = list("0", c("B", "6", "3.5")))
  expect_error(report_test_docx_table(ragged, c("time", "mean")),
    "data row width differs from the header", fixed = TRUE)
  expect_error(report_test_docx_table('<w:tbl></w:tbl>', "time"),
    "no header row", fixed = TRUE)
  missing_cell <- sub('</w:tc>', '', valid, fixed = TRUE)
  expect_error(report_test_docx_table(missing_cell, "time"),
    "malformed or missing cells", fixed = TRUE)
  missing_row <- sub('</w:tr>', '', valid, fixed = TRUE)
  expect_error(report_test_docx_table(missing_row, "time"),
    "malformed or missing rows", fixed = TRUE)
  expect_error(report_test_docx_table(sub('</w:tbl>', '', valid, fixed = TRUE), "time"),
    "malformed or missing table boundaries", fixed = TRUE)
  missing_text <- sub('</w:t>', '', valid, fixed = TRUE)
  expect_error(report_test_docx_table(missing_text, "time"), "malformed text nodes", fixed = TRUE)
  duplicate_header <- report_docx_parser_fixture(c("time", "time", "mean"))
  expect_error(report_test_docx_table(duplicate_header, "time"),
    "header names must be nonempty and unique", fixed = TRUE)
})
