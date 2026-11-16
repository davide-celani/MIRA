# Repository metadata is intentionally absent from built/installed packages.
# Source CI requires it explicitly so deleting CITATION.cff cannot hide a check.
mira_repository_validation_root <- function() {
  requested <- Sys.getenv("MIRA_VALIDATION_REPO_DIR", "")
  root <- if (nzchar(requested)) requested else testthat::test_path("..", "..")
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  if (!file.exists(file.path(root, "CITATION.cff"))) {
    required <- tolower(Sys.getenv("MIRA_REQUIRE_REPOSITORY_METADATA", "false")) %in%
      c("1", "true", "yes")
    if (required || nzchar(requested)) {
      stop("Repository metadata checks require CITATION.cff in the source root: ",
           root, call. = FALSE)
    }
    testthat::skip(
      "Repository-only metadata: CITATION.cff is excluded from built/installed packages."
    )
  }
  root
}

mira_repository_metadata <- function() {
  root <- mira_repository_validation_root()
  if (!requireNamespace("yaml", quietly = TRUE)) {
    required <- tolower(Sys.getenv("MIRA_REQUIRE_REPOSITORY_METADATA", "false")) %in%
      c("1", "true", "yes")
    if (required || nzchar(Sys.getenv("MIRA_VALIDATION_REPO_DIR", ""))) {
      stop("Required repository metadata checks need the Suggests dependency 'yaml'.",
           call. = FALSE)
    }
    testthat::skip("Repository metadata checks need the optional yaml dependency.")
  }
  description <- as.list(read.dcf(file.path(root, "DESCRIPTION"))[1L, ])
  people <- eval(parse(text = description$`Authors@R`),
    envir = list2env(list(person = utils::person), parent = baseenv()))
  if (!inherits(people, "person")) {
    stop("Authors@R must evaluate to R person objects.", call. = FALSE)
  }
  citation <- utils::readCitationFile(file.path(root, "inst", "CITATION"),
                                    meta = description)
  if (length(citation) != 1L) {
    stop("Review validation when adding a preferred citation.", call. = FALSE)
  }
  list(root = root,
    cff = yaml::read_yaml(file.path(root, "CITATION.cff"), eval.expr = FALSE),
    description = description,
    authors = people[vapply(people, function(p) "aut" %in% p$role, logical(1))],
    citation = citation)
}

mira_metadata_valid_orcid <- function(id) {
  if (!is.character(id) || length(id) != 1L || is.na(id) ||
      !grepl("^[0-9]{4}-[0-9]{4}-[0-9]{4}-[0-9]{3}[0-9X]$", id)) {
    return(FALSE)
  }
  digits <- strsplit(gsub("-", "", id, fixed = TRUE), "", fixed = TRUE)[[1L]]
  total <- 0L
  for (digit in digits[seq_len(15L)]) total <- (total + as.integer(digit)) * 2L
  check <- (12L - total %% 11L) %% 11L
  identical(digits[16L], if (check == 10L) "X" else as.character(check))
}

mira_metadata_identity <- function(given, family, entity = NULL) {
  if (!is.null(entity)) return(paste0("entity\t", entity))
  if (length(given) < 1L || length(family) < 1L) {
    stop("Every personal author needs given and family names.", call. = FALSE)
  }
  paste(paste(given, collapse = " "), paste(family, collapse = " "), sep = "\t")
}

mira_metadata_orcid <- function(value, url = FALSE) {
  if (is.null(value)) return(NA_character_)
  if (!is.character(value) || length(value) != 1L || is.na(value)) {
    stop("An ORCID must be a single string.", call. = FALSE)
  }
  if (url) {
    if (!startsWith(value, "https://orcid.org/")) {
      stop("CFF authors need a full https://orcid.org/ URL.", call. = FALSE)
    }
    value <- sub("^https://orcid.org/", "", value)
  }
  if (!mira_metadata_valid_orcid(value)) {
    stop("Invalid ORCID (format/checksum): ", value, call. = FALSE)
  }
  value
}

mira_metadata_person_authors <- function(people) {
  data.frame(
    identity = vapply(people, function(p) {
      if (is.null(p$given)) mira_metadata_identity(NULL, NULL, p$family)
      else mira_metadata_identity(p$given, p$family)
    }, character(1)),
    orcid = vapply(people, function(p) mira_metadata_orcid(p$comment[["ORCID"]]),
                   character(1)),
    stringsAsFactors = FALSE)
}
