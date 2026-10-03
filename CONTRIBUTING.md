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

Some tests skip when external tools or optional analysis engines are
unavailable. A skipped test does not establish that the corresponding feature
works. Quarto rendering tests require Quarto; PDF rendering also requires
LaTeX.

## Pull requests

Explain the problem, the resulting behavior, and how you verified it. Include
a NEWS entry for a user-visible change. Add meaningful regression tests for
bug fixes and behavioral changes, and use synthetic data in examples.

Changes to the Stan model, priors, likelihoods, or posterior estimands should
state the statistical rationale and any implications for existing analyses.
Document the relevant model assumptions and diagnostics. Passing software
tests alone does not validate a statistical method.

If AI tools assisted the contribution, describe the tool, scope, and human
verification in the pull request; see [AI usage](dev/AI_USAGE.md) for the record
of this release-preparation session. Contributors remain responsible for the
accuracy and licensing of submitted material.

By contributing, you agree that your contribution can be distributed under
MIRA's MIT license. You must have the right to submit the material.
