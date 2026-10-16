test_that("contrast verifier rejects omission, replacement, reversal and empty tables", {
  visits <- c("score_t0", "score_t1", "score_t2")
  full <- data.frame(from = visits[c(1L, 1L, 2L)], to = visits[c(2L, 3L, 3L)])
  expect_no_error(report_test_contrast_keys(full, visits))
  for (index in list(c(1L, 3L), c(1L, 1L, 3L), integer(0))) {
    expect_error(report_test_contrast_keys(full[index, , drop = FALSE], visits),
      "every expected from/to key exactly once", fixed = TRUE)
  }
  wrong <- full
  wrong[1L, ] <- full[1L, 2:1]
  expect_error(report_test_contrast_keys(wrong, visits), "exactly once", fixed = TRUE)
  empty_html <- paste0('<table><caption>Within-participant changes</caption>',
    '<tr><th>from</th><th>to</th></tr></table>')
  empty <- report_test_html_table(empty_html, "Within-participant changes")
  expect_identical(nrow(empty), 0L)
  expect_error(report_test_contrast_keys(empty, visits), "exactly once", fixed = TRUE)
})

test_that("DOCX inference verifier rejects keyed cell corruption across separate blocks", {
  raw <- report_test_inference_data()[c("score_t0", "score_t1", "score_t2")]
  # Independently construct Word XML, not via MIRA's table preparation. Each
  # statistic has its own block, deliberately separating n and mean.
  xml_table <- function(table) {
    cell <- function(x) paste0("<w:tc><w:p><w:r><w:t>",
      gsub("<", "&lt;", x, fixed = TRUE), "</w:t></w:r></w:p></w:tc>")
    rows <- rbind(names(table), as.matrix(table))
    paste0("<w:tbl>", paste(apply(rows, 1L, function(row) {
      paste0("<w:tr>", paste(cell(row), collapse = ""), "</w:tr>")
    }), collapse = ""), "</w:tbl>")
  }
  expected <- vapply(raw, report_test_mean_oracle, numeric(6L))
  expect_identical(unname(expected["n", ]), c(5, 6, 5))
  tables <- lapply(rownames(expected), function(column) {
    tab <- data.frame(time = names(raw), stringsAsFactors = FALSE)
    tab[[column]] <- sprintf("%.3f", expected[column, ])
    tab
  })
  pairs <- utils::combn(names(raw), 2L)
  delta <- lapply(1:3, function(i) raw[[pairs[2L, i]]] - raw[[pairs[1L, i]]])
  changes <- vapply(delta, report_test_mean_oracle, numeric(6L))
  expect_identical(unname(changes["n", ]), c(5, 4, 5))
  columns <- c(n = "n", mean = "mean change", sd = "sd change", se = "se change",
    "ci lower" = "ci lower", "ci upper" = "ci upper")
  for (column in names(columns)) {
    tab <- data.frame(from = pairs[1L, ], to = pairs[2L, ], stringsAsFactors = FALSE)
    tab[[columns[[column]]]] <- sprintf("%.3f", changes[column, ])
    tables[[length(tables) + 1L]] <- tab
  }
  p <- 2 * stats::pt(-abs(changes["mean", ] / changes["se", ]), changes["n", ] - 1L)
  tab <- data.frame(from = pairs[1L, ], to = pairs[2L, ], stringsAsFactors = FALSE)
  tab[["paired t p"]] <- ifelse(p < .001, "<0.001", sprintf("%.3f", p))
  tables[[length(tables) + 1L]] <- tab
  encode <- function(tables) paste(vapply(tables, xml_table, character(1L)), collapse = "")
  expect_no_error(report_test_docx_inference(encode(tables), raw))
  # Every protected numerical field is changed independently by more than its
  # display rounding bound, while all headers, keys and other cells stay valid.
  for (i in seq_along(tables)) {
    bad <- tables
    field <- ncol(bad[[i]])
    bad[[i]][1L, field] <- "99.000"
    expect_error(report_test_docx_inference(encode(bad), raw),
      "Incorrect scientific cells", fixed = TRUE)
  }
  for (index in list(1:2, c(1L, 1L, 3L), integer(0))) {
    bad <- tables
    bad[[7L]] <- bad[[7L]][index, , drop = FALSE]
    expect_error(report_test_docx_inference(encode(bad), raw), "exactly once", fixed = TRUE)
  }
  # One complete-case analysis would use only four profiles for every visit.
  bad <- tables
  bad[[1L]]$n <- rep("4", 3L)
  expect_error(report_test_docx_inference(encode(bad), raw), "scientific cells: n", fixed = TRUE)
  expect_error(report_test_docx_inference(encode(c(tables, tables[1L])), raw),
    "expected exactly one table", fixed = TRUE)
})

test_that("PDF verifier requires readable pages and keyed scientific rows", {
  good <- paste("MIRA: Longitudinal Analysis Descriptive statistics",
    "time n mean sd p value", "t0 6 10.250 1.500 <0.001", "t1 6 11.750 1.250 0.050")
  expect_no_error(report_test_pdf_science(good))
  expect_no_error(report_test_pdf_science(paste("Descriptive statistics time n mean sd p value",
    "t0 6 10.25 1.50 <0.001", "t1 6 11.75 1.25 0.050")))
  expect_error(report_test_pdf_science(character()), "no readable pages", fixed = TRUE)
  expect_error(report_test_pdf_science(""), "no readable pages", fixed = TRUE)
  for (bad in c("%PDF-1.7", "MIRA: Longitudinal Analysis",
                sub("11.750", "99.750", good, fixed = TRUE),
                sub("t0", "not0", good, fixed = TRUE),
                sub("0.050", "0.050999", good, fixed = TRUE),
                sub("t1 6 11.750 1.250 0.050", "", good, fixed = TRUE),
                sub("t1 6 11.750 1.250 0.050", "t0 6 11.750 1.250 0.050", good, fixed = TRUE))) {
    expect_error(report_test_pdf_science(bad), "keyed scientific table cells", fixed = TRUE)
  }
})
