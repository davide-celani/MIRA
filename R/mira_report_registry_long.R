#' List outcomes supported by MIRA
#'
#' Returns the canonical outcomes supported by the MIRA longitudinal model,
#' together with a concise description of each one.
#'
#' @return A data frame with one row per supported outcome and two columns:
#'   \describe{
#'   \item{outcome}{Canonical outcome name.}
#'   \item{description}{Plain-language description of the outcome.}
#'   }
#'
#' @examples
#' mira_outcome_registry_long()
#'
#' @export
mira_outcome_registry_long <- function() {
  data.frame(
    outcome = c("BCVA", "CMT", "generic"),
    description = c(
      "Best-corrected visual acuity, measured as an ETDRS letter score.",
      "Central macular thickness, typically measured in micrometres.",
      "A user-defined continuous longitudinal outcome."
    ),
    stringsAsFactors = FALSE
  )
}
