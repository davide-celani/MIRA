# Contributing to MIRA

Contributions to MIRA are welcome: reproducible bug reports, improvements to
documentation, tests, and proposals for statistical or computational changes.
Please follow the [Code of Conduct](CODE_OF_CONDUCT.md).

## Questions and proposals

Use the repository's GitHub Issues for public questions, bugs, and feature
requests. Search existing issues first. Describe the research question or
workflow that motivates a proposed feature, and discuss substantial model or
API changes before implementing them. The maintainer reviews contributions;
opening an issue does not guarantee an implementation date.

For confidential project matters, contact Davide Celani at
<davide.celani.00@gmail.com>.

## Reproducible bug reports

Include a small, runnable example, the result you expected, the actual result,
and any error messages. Include `sessionInfo()` and the MIRA version
(`packageVersion("MIRA")`). For Bayesian fitting issues, include the CmdStanR
and CmdStan versions, random seed, sampling settings, and available diagnostic
output. For rendering issues, include the Quarto version and requested format.

Use simulated or suitably anonymized data. Do not upload identifiable patient
data, credentials, or private study files to public issues.

## Development workflow

1. Fork the repository and create a branch for a focused change.
2. Install the package dependencies and the development packages `devtools`,
   `roxygen2`, and `testthat`. Follow the README for CmdStanR, CmdStan, and
   reporting tools needed by your change.
3. Work from the package directory. Load the source with
   `devtools::load_all()` and make the change.
4. Update roxygen comments when the public interface changes, then run
   `devtools::document()`. Generated `man/` files and `NAMESPACE` should stay
   consistent with the source.
5. Run relevant existing tests with `devtools::test()` and check the package
   with `devtools::check()` before opening a pull request. Report failures,
   warnings, and skipped checks in the pull request.

Tests are in `tests/testthat/`. The CmdStan smoke test is opt-in and needs a
working CmdStan installation and C++ toolchain:

```r
Sys.setenv(MIRA_RUN_CMDSTAN_TESTS = "true")
devtools::test(filter = "dynamic-covariates_long")
```

External integrations are disabled in the ordinary suite. A skipped test does
not establish that the corresponding feature works. Explicit rendering profiles
require Quarto; PDF rendering also requires LaTeX. Tests never install tools or
dependencies.

## Validation profiles

Functional and scientific checks are in `tests/testthat/`; the outer runner
below also audits the full suite's session state. The longer campaigns use
the same framework and are enabled explicitly:

| Profile | Environment flag set to `true` | `testthat::test_local()` filter |
| --- | --- | --- |
| CmdStan sampling | `MIRA_RUN_CMDSTAN_TESTS` | `dynamic-covariates_long` |
| Quarto HTML and DOCX | `MIRA_RUN_QUARTO_TESTS` | `mira-report-freq_long-(integration\|rendering-contracts)` |
| Quarto PDF | `MIRA_RUN_QUARTO_PDF_TESTS` | `mira-report-freq_long-rendering-contracts` |
| Extended Monte Carlo study | `MIRA_RUN_SCIENCE_TESTS` | `validation-paired-change-monte-carlo` |
| Targeted mutation campaign | `MIRA_RUN_MUTATION_TESTS` | `validation-mutations` |
| Rendered, extracted and visible vignette code | `MIRA_RUN_VIGNETTE_TESTS` | `validation-vignettes` |

For example, from the repository root:

```r
Sys.setenv(MIRA_RUN_SCIENCE_TESTS = "true")
testthat::test_local(filter = "validation-paired-change-monte-carlo",
                     stop_on_failure = TRUE, stop_on_warning = TRUE)
```

Set `MIRA_VALIDATION_OUTPUT_DIR` to an existing output directory outside the repository
to retain campaign evidence; otherwise campaigns use temporary directories.
Mutation and vignette campaigns need the source repository. Set
`MIRA_VALIDATION_REPO_DIR` explicitly when it cannot be inferred from the test
directory. `MIRA_MUTATION_CASES` selects comma-separated mutant names; its
default preserves the six-case campaign. Mutation summaries distinguish
detected, survived, non-reached, invalid and blocked cases: survival or
non-reach is a coverage limitation, not evidence of detection.

Vignette checks run each article's rendering, extracted R code and reader-visible
R code in separate clean R processes. `MIRA_RUN_ADVANCED_VIGNETTES=true` enables
advanced tutorial branches; `MIRA_RENDER_VIGNETTE_FORMATS=html,docx,pdf` requests
external report formats. `MIRA_VALIDATE_INSTALLED_VIGNETTES=true` additionally
checks installed vignette files and requires a package built with its vignettes.
HTML link/anchor checks use the optional `xml2` dependency.

Repository metadata consistency runs in the source suite. CI sets
`MIRA_REQUIRE_REPOSITORY_METADATA=true` so missing source metadata fails there;
tests of built packages skip this repository-only check because `CITATION.cff`
is deliberately excluded from the package archive.

For an ordinary full-suite run with RNG, working-directory, graphics-device,
environment and existing-option restoration checks, create an output directory
outside the repository and run:

```text
Rscript --vanilla tests/run-validation.R <repository> <output-directory>/source source
Rscript --vanilla tests/run-validation.R <repository> <output-directory>/installed installed
```

The installed mode uses the MIRA installation selected by the R library paths.
This outer runner disables all optional profiles, records CSV/RDS results and
session state, and exits nonzero on failures, warnings or unexpected state
changes. It is excluded from the built package to avoid automatic recursive
execution by `R CMD check`. The extended CI validation job runs the long
campaigns on pushes to `main`/`master` and on manual dispatch. The regular
checks and Quarto HTML/DOCX integration also run on pull requests.

Pushes to `main`/`master` run nine automatic checks: five R configurations,
real Quarto HTML/DOCX integration, CmdStan sampling, citation validation and
the extended campaigns. Pull requests retain the eight regular checks.
Before changing `.github/workflows/`, run `actionlint` from the repository root
to validate GitHub Actions expressions and context availability as well as YAML.

Workflow definitions belong in `.github/workflows/` and automated tests in
`tests/`. Development examples in `dev/` remain available in Git, while
`.Rbuildignore` excludes them from the distributed R package.

Before a release, keep `DESCRIPTION`, `CITATION.cff`, `inst/CITATION` and `NEWS.md`
consistent, run the applicable profiles and inspect the CI results. Verify the
native GitHub-Zenodo integration and the archived version/DOI after publishing.

## Pull requests

Explain the problem, the resulting behavior, and how you verified it. Include
a NEWS entry for a user-visible change. Add meaningful regression tests for
bug fixes and behavioral changes, and use synthetic data in examples.

Changes to the Stan model, priors, likelihoods, or posterior estimands should
state the statistical rationale and any implications for existing analyses.
Document the relevant model assumptions and diagnostics. Passing software
tests alone does not validate a statistical method.

If AI tools assisted the contribution, describe the tool, scope, and human
verification in the pull request. Contributors remain responsible for the
accuracy and licensing of submitted material.

By contributing, you agree that your contribution can be distributed under
MIRA's MIT license. You must have the right to submit the material.
