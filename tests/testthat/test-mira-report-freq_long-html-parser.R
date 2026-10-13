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
