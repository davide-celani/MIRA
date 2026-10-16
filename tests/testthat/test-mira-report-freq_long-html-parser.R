report_parser_fixture <- function(caption) {
  paste0('<table class="mira-table">', caption,
    '<thead><tr><th>time</th><th>n</th><th>mean</th><th>p value</th></tr></thead>',
    '<tbody><tr><td>t0</td><td>6</td><td>3.500</td><td>&lt;0.001</td></tr>',
    '<tr><td>t1</td><td>5</td><td>4.500</td><td>0.123</td></tr></tbody></table>')
}

test_that("published HTML captions tolerate whitespace, attributes, and inline markup", {
  captions <- c(
    '<caption>\nDescriptive statistics\n</caption>',
    '<caption class="table-caption">  Descriptive statistics </caption>',
    '<caption id="tbl-description"><strong>Descriptive</strong>\n statistics</caption>'
  )
  for (caption in captions) {
    observed <- report_test_html_table(report_parser_fixture(caption), "Descriptive statistics")
    expect_identical(names(observed), c("time", "n", "mean", "p value"))
    expect_identical(observed$time, c("t0", "t1"))
    expect_identical(observed$n, c("6", "5"))
    expect_identical(observed$mean, c("3.500", "4.500"))
    expect_identical(observed[["p value"]], c("<0.001", "0.123"))
  }
})

test_that("caption and cell entities decode exactly once", {
  html <- report_parser_fixture(
    '<caption>A &amp; B &amp;lt; &#945;&#160;&#x03B2; &quot;quoted&quot;</caption>')
  html <- sub('<td>t0</td>', '<td>A &lt; B &amp; C &amp;quot;</td>', html, fixed = TRUE)
  observed <- report_test_html_table(html, 'A & B &lt; \u03b1 \u03b2 "quoted"')
  expect_identical(observed$time, c('A < B & C &quot;', 't1'))
  expect_identical(observed[["p value"]], c("<0.001", "0.123"))
})

test_that("missing, near-matching, and duplicate captions fail clearly", {
  table <- report_parser_fixture('<caption>Descriptive statistics</caption>')
  expect_error(report_test_html_table(table, "Within-participant changes"),
    "Expected exactly one HTML table with caption 'Within-participant changes'; found 0.",
    fixed = TRUE)
  block <- report_parser_fixture('<caption>Descriptive statistics (block 1 of 2)</caption>')
  expect_error(report_test_html_table(block, "Descriptive statistics"),
    "Expected exactly one HTML table with caption 'Descriptive statistics'; found 0.",
    fixed = TRUE)
  expect_error(report_test_html_table(paste0(table, table), "Descriptive statistics"),
    "Expected exactly one HTML table with caption 'Descriptive statistics'; found 2.",
    fixed = TRUE)
})

# The verifier must reject absent cells before rbind can recycle a short row.
# The literal zero makes that failure scientifically plausible (0, 0, 0).
test_that("published HTML rejects malformed rows before recycling scientific cells", {
  valid <- paste0('<table><caption>Results</caption>',
    '<thead><tr><th>id</th><th>n</th><th>mean</th></tr></thead>',
    '<tbody><tr><td>0</td><td>0</td><td>0</td></tr>',
    '<tr><td>B</td><td>6</td><td>3.5</td></tr></tbody></table>')
  ragged <- sub('<td>0</td><td>0</td><td>0</td>', '<td>0</td>', valid, fixed = TRUE)
  expect_error(report_test_html_table(ragged, "Results"),
    "data row 1 has 1 cells; expected 3", fixed = TRUE)
  unclosed <- sub('<td>3.5</td>', '<td>3.5', valid, fixed = TRUE)
  expect_error(report_test_html_table(unclosed, "Results"),
    "malformed or mismatched cells in row 3", fixed = TRUE)
  mismatched <- sub('<td>3.5</td>', '<td>3.5</th>', valid, fixed = TRUE)
  expect_error(report_test_html_table(mismatched, "Results"),
    "malformed or mismatched cells in row 3", fixed = TRUE)
  missing_row <- sub('</tr></tbody>', '</tbody>', valid, fixed = TRUE)
  expect_error(report_test_html_table(missing_row, "Results"),
    "malformed or missing rows", fixed = TRUE)
  missing_table <- sub('</table>', '', valid, fixed = TRUE)
  expect_error(report_test_html_table(missing_table, "Results"),
    "malformed or missing table boundaries", fixed = TRUE)
})

test_that("published HTML requires an unambiguous header before indexing", {
  valid <- report_parser_fixture('<caption>Descriptive statistics</caption>')
  no_header <- gsub('th>', 'td>', valid, fixed = TRUE)
  expect_error(report_test_html_table(no_header, "Descriptive statistics"),
    "first row must contain only header cells", fixed = TRUE)
  duplicate_header <- sub('<th>n</th>', '<th>time</th>', valid, fixed = TRUE)
  expect_error(report_test_html_table(duplicate_header, "Descriptive statistics"),
    "header names must be nonempty and unique", fixed = TRUE)
  empty_header <- sub('<th>n</th>', '<th> </th>', valid, fixed = TRUE)
  expect_error(report_test_html_table(empty_header, "Descriptive statistics"),
    "header names must be nonempty and unique", fixed = TRUE)
  empty_table <- '<table><caption>Descriptive statistics</caption></table>'
  expect_error(report_test_html_table(empty_table, "Descriptive statistics"),
    "no header row", fixed = TRUE)
  repeated_header <- sub('<td>t0</td>', '<th>t0</th>', valid, fixed = TRUE)
  expect_error(report_test_html_table(repeated_header, "Descriptive statistics"),
    "data row 1 must contain only data cells", fixed = TRUE)
  spanning <- sub('<td>t0</td>', '<td colspan="2">t0</td>', valid, fixed = TRUE)
  expect_error(report_test_html_table(spanning, "Descriptive statistics"),
    "rowspan and colspan are unsupported", fixed = TRUE)
})

test_that("published HTML preserves empty, NA, zero, signs, precision and table identity", {
  target <- paste0('<table><caption>Results</caption>',
    '<tr><th>id</th><th>empty</th><th>missing</th><th>value</th></tr>',
    '<tr><td><strong>A</strong> &#956;</td><td></td><td>NA</td><td>0</td></tr>',
    '<tr><td>B<br>C</td><td> </td><td>NaN</td><td>-0.0010</td></tr></table>')
  decoy <- report_parser_fixture('<caption>Results similar</caption>')
  observed <- report_test_html_table(paste0(decoy, target, decoy), "Results")
  expect_identical(names(observed), c("id", "empty", "missing", "value"))
  expect_identical(observed$id, c("A \u03bc", "B C"))
  expect_identical(observed$empty, c("", ""))
  expect_identical(observed$missing, c("NA", "NaN"))
  expect_identical(observed$value, c("0", "-0.0010"))
})
