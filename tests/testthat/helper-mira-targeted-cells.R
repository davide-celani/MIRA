# Independent checks shared by the real document integrations and their
# adversarial, non-rendering verifier tests. No expected cell comes from MIRA.
report_test_contrast_keys <- function(change, visits) {
  expected <- apply(utils::combn(visits, 2L), 2L, paste, collapse = "::")
  if (!all(c("from", "to") %in% names(change))) {
    stop("Contrast table requires from/to keys.", call. = FALSE)
  }
  keys <- paste(change$from, change$to, sep = "::")
  counts <- tabulate(match(keys, expected), nbins = length(expected))
  if (nrow(change) != length(expected) || anyDuplicated(keys) ||
      !setequal(keys, expected) || any(counts != 1L)) {
    stop("Contrast table must contain every expected from/to key exactly once.",
         call. = FALSE)
  }
  invisible(keys)
}

report_test_inference_data <- function() {
  # First two visits reuse the independent descriptive oracle fixture from
  # test-mira-info_long-statistical-oracles.R. A missing third visit makes
  # descriptive n=(5,6,5) differ from pairwise n=(5,4,5).
  data.frame(patient = 1:6, score_t0 = c(1, 2, 4, 8, NA, 16),
    score_t1 = c(3, 5, 7, 9, 11, 13), score_t2 = c(4, 7, NA, 11, 15, 18))
}

report_test_mean_oracle <- function(x, alpha = 0.10) {
  # Same base-R mean/SD/Student interval oracle as statistical-oracles.R;
  # apply to finite observations for descriptives and paired differences.
  x <- x[is.finite(x)]
  n <- length(x)
  se <- stats::sd(x) / sqrt(n)
  half_width <- stats::qt(1 - alpha / 2, n - 1L) * se
  c(n = n, mean = mean(x), sd = stats::sd(x), se = se,
    "ci lower" = mean(x) - half_width, "ci upper" = mean(x) + half_width)
}

report_test_numeric_cells <- function(observed, expected, label, rounding = 0.0005000001) {
  # Three printed decimals: half a last-place unit plus floating arithmetic
  # slack, not a statistical/engine tolerance. Counts use rounding=0.
  numbers <- as.numeric(observed)
  if (length(numbers) != length(expected) || any(!is.finite(numbers)) ||
      any(abs(numbers - expected) > rounding)) {
    stop(paste0("Incorrect scientific cells: ", label, "."), call. = FALSE)
  }
  invisible(TRUE)
}

report_test_docx_inference <- function(xml, raw, alpha = 0.10) {
  visits <- names(raw)
  expected <- vapply(raw, report_test_mean_oracle, numeric(6L), alpha = alpha)
  for (column in rownames(expected)) {
    # Each signature identifies one block, irrespective of column partition.
    block <- report_test_docx_table(xml, c("time", column))
    if (!identical(block$time, visits)) stop("Incorrect descriptive visit keys.", call. = FALSE)
    report_test_numeric_cells(block[[column]], expected[column, ], column,
      rounding = if (column == "n") 0 else 0.0005000001)
  }
  pairs <- utils::combn(visits, 2L)
  differences <- lapply(seq_len(ncol(pairs)), function(i) {
    raw[[pairs[2L, i]]] - raw[[pairs[1L, i]]]
  })
  changes <- vapply(differences, report_test_mean_oracle, numeric(6L), alpha = alpha)
  columns <- c(n = "n", mean = "mean change", sd = "sd change",
    se = "se change", "ci lower" = "ci lower", "ci upper" = "ci upper")
  expected_keys <- apply(pairs, 2L, paste, collapse = "::")
  for (column in names(columns)) {
    block <- report_test_docx_table(xml, c("from", "to", columns[[column]]))
    keys <- report_test_contrast_keys(block, visits)
    report_test_numeric_cells(block[[columns[[column]]]],
      changes[column, match(keys, expected_keys)], columns[[column]],
      rounding = if (column == "n") 0 else 0.0005000001)
  }
  block <- report_test_docx_table(xml, c("from", "to", "paired t p"))
  keys <- report_test_contrast_keys(block, visits)
  p <- 2 * stats::pt(-abs(changes["mean", ] / changes["se", ]),
                     df = changes["n", ] - 1L)
  labels <- ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
  if (!identical(block[["paired t p"]], unname(labels[match(keys, expected_keys)]))) {
    stop("Incorrect scientific cells: paired t p.", call. = FALSE)
  }
  invisible(TRUE)
}

report_test_pdf_pages <- function(path) {
  if (requireNamespace("pdftools", quietly = TRUE)) {
    info <- pdftools::pdf_info(path)
    pages <- pdftools::pdf_text(path)
    if (info$pages < 1L || length(pages) != info$pages) stop("PDF page count mismatch.")
    return(pages)
  }
  commands <- Sys.which(c("pdfinfo", "pdftotext"))
  if (any(!nzchar(commands))) {
    stop(paste("PDF scientific inspection requires existing pdftools or both",
      "pdfinfo and pdftotext; tests never install them."), call. = FALSE)
  }
  output <- tempfile("mira-pdf-text-", fileext = ".txt")
  errors <- tempfile("mira-pdf-inspect-", fileext = ".log")
  info_file <- tempfile("mira-pdf-info-", fileext = ".txt")
  on.exit(unlink(c(output, errors, info_file)), add = TRUE)
  status <- system2(commands[["pdfinfo"]], shQuote(path), stdout = info_file, stderr = errors)
  if (status != 0L) stop("PDF inspection failed (pdfinfo).")
  info <- readLines(info_file, warn = FALSE)
  count <- grep("^Pages:[[:space:]]+[0-9]+", info, value = TRUE)
  if (length(count) != 1L) stop("PDF has no verified page count.")
  n <- as.integer(sub("^Pages:[[:space:]]+", "", count))
  status <- system2(commands[["pdftotext"]], c("-layout", shQuote(path), shQuote(output)),
                    stdout = errors, stderr = errors)
  if (status != 0L || !file.exists(output)) stop("PDF inspection failed (pdftotext).")
  text <- paste(readLines(output, warn = FALSE), collapse = "\n")
  pages <- strsplit(text, "\f", fixed = TRUE)[[1L]]
  if (length(pages) > n && !nzchar(trimws(tail(pages, 1L)))) pages <- head(pages, -1L)
  if (n < 1L || length(pages) != n) stop("PDF page count mismatch.")
  pages
}

report_test_pdf_science <- function(pages) {
  if (!length(pages) || anyNA(pages) || !any(nzchar(trimws(pages)))) {
    stop("PDF contains no readable pages.", call. = FALSE)
  }
  text <- paste0(" ", trimws(gsub("[[:space:]]+", " ", paste(pages, collapse = " "))), " ")
  # Keys and cells must occur together. knitr rounds to at most three places,
  # but does not pad every numeric column to exactly three decimal places.
  if (!grepl(" Descriptive statistics ", text, fixed = TRUE) ||
      !grepl(" time n mean sd p value ", text, fixed = TRUE)) {
    stop("PDF is missing keyed scientific table cells.", call. = FALSE)
  }
  tokens <- strsplit(trimws(text), " ", fixed = TRUE)[[1L]]
  visits <- c("t0", "t1")
  expected <- rbind(c(10.25, 1.5), c(11.75, 1.25))
  p <- c("<0.001", "0.050")
  for (i in seq_along(visits)) {
    hits <- which(tokens == visits[[i]])
    valid <- vapply(hits, function(index) {
      if (index + 4L > length(tokens)) return(FALSE)
      cells <- tokens[index + 1:4]
      if (cells[[1L]] != "6" || cells[[4L]] != p[[i]] ||
          !all(grepl("^[+-]?[0-9]+(?:\\.[0-9]+)?$", cells[2:3], perl = TRUE))) return(FALSE)
      all(abs(as.numeric(cells[2:3]) - expected[i, ]) <= 0.0005000001)
    }, logical(1L))
    if (sum(valid) != 1L) stop("PDF is missing keyed scientific table cells.", call. = FALSE)
  }
  invisible(TRUE)
}
