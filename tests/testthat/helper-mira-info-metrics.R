# Statistical oracle tests use canonical base names so the same assertions apply
# to explicitly named and inferred outcomes. Public naming/ordering contracts
# are tested separately without this view; no result values are changed here.
mira_test_info_base <- function(...) {
  result <- mira_info_long(...)
  if (inherits(result, "mira_info_multi_long")) {
    result$outcomes <- lapply(result$outcomes, MIRA:::.mira_metric_base_view_long)
    result
  } else MIRA:::.mira_metric_base_view_long(result)
}
