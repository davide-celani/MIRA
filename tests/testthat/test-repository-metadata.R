test_that("CITATION.cff remains the repository's Zenodo metadata source", {
  root <- mira_repository_validation_root()
  expect_false(file.exists(file.path(root, ".zenodo.json")),
    info = ".zenodo.json overrides CITATION.cff; review and extend validation first.")
})

test_that("repository author identities and valid ORCIDs agree across citation formats", {
  metadata <- mira_repository_metadata()
  expected_orcid <- "0009-0006-3321-3928"
  expected_author <- "Davide\tCelani"
  expect_true(mira_metadata_valid_orcid(expected_orcid))
  authors <- list(
    CFF = data.frame(
      identity = vapply(metadata$cff$authors, function(p) {
        mira_metadata_identity(p$`given-names`, p$`family-names`, p$name)
      }, character(1)),
      orcid = vapply(metadata$cff$authors, function(p) {
        mira_metadata_orcid(p$orcid, url = TRUE)
      }, character(1)),
      stringsAsFactors = FALSE),
    DESCRIPTION = mira_metadata_person_authors(metadata$authors),
    CITATION = mira_metadata_person_authors(metadata$citation$author))

  for (source in names(authors)) {
    current <- authors[[source]]
    expect_true(nrow(current) > 0L, info = source)
    expect_identical(anyDuplicated(current$identity), 0L, info = source)
    expect_identical(anyDuplicated(current$orcid[!is.na(current$orcid)]), 0L,
                     info = source)
    target <- which(current$identity == expected_author)
    expect_identical(length(target), 1L, info = source)
    expect_identical(current$orcid[target], expected_orcid, info = source)
    current <- current[order(current$identity), , drop = FALSE]
    rownames(current) <- NULL
    authors[[source]] <- current
  }
  expect_identical(authors$CFF, authors$DESCRIPTION)
  expect_identical(authors$CFF, authors$CITATION)
})

test_that("repository citation title, version, DOI and source URL agree", {
  metadata <- mira_repository_metadata()
  expect_identical(metadata$cff$title, metadata$citation$title)
  expect_identical(as.character(metadata$cff$version), metadata$description$Version)
  expect_identical(metadata$citation$note,
    sprintf("R package version %s", metadata$description$Version))
  expect_identical(metadata$cff$doi, metadata$citation$doi)
  expect_identical(metadata$cff$`repository-code`, metadata$description$URL)
  expect_identical(metadata$cff$`repository-code`, metadata$citation$url)
})
