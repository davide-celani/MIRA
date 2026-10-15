# Run from the repository root: Rscript --vanilla dev/validate-metadata.R
# Uses yaml, already in Suggests. Schema validation remains in citation.yaml.
# Reads metadata only; does not install MIRA, contact Zenodo, or modify Git.
if (!requireNamespace("yaml", quietly = TRUE)) {
  stop("Install the existing Suggests dependency 'yaml' to validate metadata.")
}
if (file.exists(".zenodo.json")) {
  stop(".zenodo.json overrides CITATION.cff in Zenodo. Review the metadata ",
       "source and extend this check before introducing that file.")
}

assert <- function(ok, message) {
  if (!isTRUE(ok)) stop(message, call. = FALSE)
}
expected_orcid <- "0009-0006-3321-3928"
expected_author <- "Davide\tCelani"

valid_orcid <- function(id) {
  if (length(id) != 1L || is.na(id) ||
      !grepl("^[0-9]{4}-[0-9]{4}-[0-9]{4}-[0-9]{3}[0-9X]$", id)) {
    return(FALSE)
  }
  digits <- strsplit(gsub("-", "", id, fixed = TRUE), "", fixed = TRUE)[[1L]]
  total <- 0L
  for (digit in digits[seq_len(15L)]) total <- (total + as.integer(digit)) * 2L
  check <- (12L - total %% 11L) %% 11L
  identical(digits[16L], if (check == 10L) "X" else as.character(check))
}

identity <- function(given, family, entity = NULL) {
  if (!is.null(entity)) return(paste0("entity\t", entity))
  assert(length(given) >= 1L && length(family) >= 1L,
         "Every personal author needs given and family names.")
  paste(paste(given, collapse = " "), paste(family, collapse = " "), sep = "\t")
}
orcid_value <- function(value, url = FALSE) {
  if (is.null(value)) return(NA_character_)
  assert(is.character(value) && length(value) == 1L && !is.na(value),
         "An ORCID must be a single string.")
  if (url) {
    assert(startsWith(value, "https://orcid.org/"),
           "CFF authors need a full https://orcid.org/ URL.")
    value <- sub("^https://orcid.org/", "", value)
  }
  assert(valid_orcid(value), paste("Invalid ORCID (format/checksum):", value))
  value
}
from_person <- function(people) {
  data.frame(
    identity = vapply(people, function(p) {
      if (is.null(p$given)) identity(NULL, NULL, p$family)
      else identity(p$given, p$family)
    }, character(1)),
    orcid = vapply(people, function(p) orcid_value(p$comment[["ORCID"]]),
                   character(1)),
    stringsAsFactors = FALSE
  )
}
check_authors <- function(authors, source) {
  assert(nrow(authors) > 0L && !anyDuplicated(authors$identity),
         paste(source, "has no authors or has duplicate author identities."))
  assert(!anyDuplicated(authors$orcid[!is.na(authors$orcid)]),
         paste(source, "assigns the same ORCID to multiple authors."))
  target <- which(authors$identity == expected_author)
  assert(length(target) == 1L && identical(authors$orcid[target], expected_orcid),
         paste(source, "must associate the expected ORCID with Davide Celani."))
  authors <- authors[order(authors$identity), , drop = FALSE]
  rownames(authors) <- NULL
  authors
}

assert(valid_orcid(expected_orcid), "The configured ORCID fails its checksum.")
cff <- yaml::read_yaml("CITATION.cff", eval.expr = FALSE)
description <- as.list(read.dcf("DESCRIPTION")[1L, ])
people <- eval(parse(text = description$`Authors@R`),
               envir = list2env(list(person = utils::person), parent = baseenv()))
assert(inherits(people, "person"), "Authors@R must evaluate to R person objects.")
package_authors <- people[vapply(people, function(p) "aut" %in% p$role, logical(1))]
citation <- utils::readCitationFile("inst/CITATION", meta = description)
assert(length(citation) == 1L, "Review validation when adding a preferred citation.")

cff_authors <- check_authors(data.frame(
  identity = vapply(cff$authors, function(p) {
    identity(p$`given-names`, p$`family-names`, p$name)
  }, character(1)),
  orcid = vapply(cff$authors, function(p) orcid_value(p$orcid, url = TRUE),
                 character(1)),
  stringsAsFactors = FALSE
), "CITATION.cff")
r_authors <- check_authors(from_person(package_authors), "DESCRIPTION Authors@R")
bib_authors <- check_authors(from_person(citation$author), "inst/CITATION")
assert(identical(cff_authors, r_authors) && identical(cff_authors, bib_authors),
       "Author identities/ORCIDs differ between CFF, DESCRIPTION and R citation.")
assert(identical(cff$title, citation$title), "CFF and R citation titles differ.")
assert(identical(as.character(cff$version), description$Version),
       "CFF and DESCRIPTION versions differ.")
assert(identical(citation$note, sprintf("R package version %s", description$Version)),
       "R citation version note differs from DESCRIPTION.")
assert(identical(cff$doi, citation$doi), "CFF and R citation DOIs differ.")
assert(identical(cff$`repository-code`, description$URL) &&
         identical(cff$`repository-code`, citation$url),
       "Repository URLs differ between metadata files.")
cat("PASS: ORCID checksum, Davide Celani association, all authors/ORCIDs,\n",
    "      CFF/R citation title, version, DOI and repository URL agree.\n",
    "      CITATION.cff is the configured Zenodo source (.zenodo.json absent).\n",
    "      Run cffconvert --validate separately for the full CFF schema.\n",
    sep = "")
