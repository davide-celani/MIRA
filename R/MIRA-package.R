#' MIRA: Bayesian Multilevel Inference for Longitudinal Data
#'
#' MIRA provides data preparation, configurable priors, Stan fitting, and
#' posterior summaries for continuous longitudinal outcomes. Frequentist
#' exploration and reproducible report bundles are also available.
#'
#' @details
#' Start with [mira_data_long()], configure priors with [mira_prior_long()],
#' fit with [mira_fit_long()], and inspect [mira_summary_long()]. Bayesian
#' fitting requires CmdStan and a C++ toolchain. For data exploration use
#' [mira_info_long()] and [mira_report_freq_long()].
#'
#' @seealso <https://github.com/davide-celani/MIRA>
#' @importFrom stats coef p.adjust.methods
#' @importFrom utils tail
#' @keywords internal
"_PACKAGE"

# Names evaluated in ggplot2 data masks rather than the package namespace.
utils::globalVariables(c(
  "category", "change_from_baseline", "ci_lower", "ci_upper", "label",
  "mean_change", "mean_difference_b_minus_a", "patient", "percent",
  "time_index", "time_label", "to_label", "trajectory_direction",
  "unavailable_pct", "value", "xend", "yend"
))
