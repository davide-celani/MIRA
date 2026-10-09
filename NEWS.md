# MIRA (development version)

- Preserve the chronology resolved by custom longitudinal-variable parsers in
  analysis results, including baseline/final labels and change directions.
- Support one-subject longitudinal data in the non-finite availability audit.
- Preserve terminal escaped underscores and dollar signs in the report's
  internal LaTeX wrapping helper.
- Expand deterministic statistical-oracle, configuration, report-bundle,
  publication-formatting, and rendering-failure regression tests.

# MIRA 0.0.1 (2026-10-03)

Initial release, archived on Zenodo with version DOI
[10.5281/zenodo.23120184](https://doi.org/10.5281/zenodo.23120184).
The concept DOI for all MIRA versions is
[10.5281/zenodo.23120183](https://doi.org/10.5281/zenodo.23120183).

## Initial capabilities

- Preparation and validation of wide longitudinal data with `mira_data_long()`,
  treatment-arm metadata, measurement times, outcome direction, and
  user-selected subject-level covariates.
- Bayesian hierarchical longitudinal models fitted through CmdStanR and the
  bundled Stan source. Student-t, Gaussian, and log-normal likelihood options
  support different continuous-outcome scales.
- Configurable prior objects and validation, including BCVA and CMT prior
  helpers and conversion to Stan data.
- Posterior summaries of population trajectories, individual changes,
  treatment contrasts, covariate effects, responder probabilities, clinically
  meaningful change, heterogeneity, posterior predictive quantities, and
  sampling diagnostics through `mira_summary_long()`.
- Longitudinal variable detection and frequentist exploratory analysis through
  `mira_detect_long()` and `mira_info_long()`.
- Quarto report bundles for frequentist analysis through
  `mira_report_freq_long()`.
- An outcome registry for BCVA, CMT, and generic continuous outcomes, and an
  existing `testthat` suite with opt-in CmdStan smoke testing.

## Release documentation

- Added software citation metadata for GitHub/Zenodo and R, installation and
  synthetic workflow examples, and release/JOSS preparation guides.
- Added package, citation, and CmdStan smoke-test workflows for GitHub Actions.
- Added contribution guidance, a project Code of Conduct, and GitHub forms
  for bug reports and feature proposals.
- Recorded AI assistance during release preparation in `dev/AI_USAGE.md`,
  with human review still pending.

## Installation and verification fixes

- Normalize accented report filenames consistently across operating systems.
- Allow machine-precision rounding differences in the correlation boundary
  test, matching the tolerance already used for Spearman correlations.
- Compile Stan executables in a writable session cache, allowing package
  installations and custom source directories to remain read-only.
- Added regression coverage for compilation-directory handling and expanded
  the opt-in CmdStan smoke test to Student-t, Gaussian, and log-normal models,
  each with zero or one encoded covariate.
- Declared missing base-R namespaces, imported previously unresolved helpers,
  removed unused mandatory dependencies, and declared ggplot2 data-mask names
  for package checks.
- Exclude repository governance, CI, citation metadata, and developer scripts
  from the built R package while retaining them in the source repository.

This inventory describes the initial release. R CMD check passed on Windows,
Linux, and macOS; citation validation and CmdStan execution smoke tests passed.
These checks do not establish completed statistical validation or JOSS
acceptance.
