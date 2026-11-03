# MIRA (development version)

- Audit and refine the complete help for `mira_info_long()` and
  `mira_report_freq_long()` against their implementation. Clarify model-switch
  resolution, analysis samples, native engine boundaries, rank uncertainty,
  MMRM contrast families and likelihood criteria, report diagnostics and
  archive depth. Extend the existing examples with direct KR/rank/MMRM result
  access, incomplete-data and unavailable-visit cases, and executable
  source-only report workflows; preserve the analytical implementation and API.

- Rebuild all frequentist vignette examples around ophthalmic studies with
  BCVA ETDRS letters, defined OCT central subfield thickness and IOP in mmHg.
  Declare one prespecified study eye per participant, scale-dependent direction,
  illustrative response bands, structure-function distinctions and synthetic
  strategy comparisons. Recompute scientific results and archive checks while
  preserving the analysis implementations and API.

- Complete the frequentist handbook's learning bridges: two visits before
  three, direct reading of means and paired changes, a minimal source-only
  archive, and a separate first-model article before adjustment and alternatives.
  Expose all required reader code; add executed reference, independent-missingness,
  missing-group, one-pair and transfer checkpoints. Validate rendered, extracted
  and reader-visible paths in separate R processes with a flat linked handbook.

- Refine the frequentist handbook with independent sample, interval, ICC,
  contrast-orientation and report-archive checks. Align the complete workflow's
  adjusted primary contrast with its stated difference-in-changes estimand,
  document the ANOVA adapter's distinct mean structure, and handle unavailable
  endpoint pairs explicitly in tutorial extraction code.

- Add a ten-article frequentist longitudinal handbook with deterministic
  datasets, explicit detection and inference assumptions, multi-outcome
  workflows, diagnostic checks and auditable source-only report examples.
  Register the standard knitr vignette builder and ordered pkgdown articles.

- Add joint Kenward-Roger F tests, unadjusted unweighted relative marginal
  effects with asymptotic ATS/WTS, and an unstructured Gaussian MMRM with
  REML/Satterthwaite inference, preserving all existing longitudinal methods.
  Declare `pbkrtest`, `nparLD >= 2.3.3`, and `mmrm >= 0.3.17` dependencies.
- Integrate the new methods into metric naming, print/summary methods,
  sensitivity comparisons, MMRM contrast multiplicity, and report CSV/RDS
  exports. Preserve native objects and expose requested non-estimability
  reasons without presenting disabled modules as fitted results.
- Preserve incomplete subject-by-visit schedules for the unadjusted rank
  module independently of model covariate filtering; retain MMRM's selected
  visit order and report native likelihood-criterion conventions explicitly.
- Add independent subject-level KR/MMRM and rank-effect oracles, adapter and
  sample-isolation checks, and real report-runtime/HTML table assertions.
- Restore the caller's numeric display precision after native object printing
  in reports, including GEE printers that change the global `digits` option.

- Repeat complete visit and treatment-arm identifiers in split report tables,
  preserving the identity of each between-arm comparison.
- Mark wholly non-estimable repeated-measures ANOVA as unavailable while
  retaining its raw fit and diagnostics, including its reason in reports.
- Validate explicit reference arms even when all group values are missing,
  and identify empty outcome settings with argument-specific errors.
- Strengthen marginal-contrast oracles and PDF duplicate-row detection;
  restore RNG isolation in test fixtures and repair prefixed-metric checks
  in the scientific and bundle validation scripts.
- Preserve literal backticks and backslashes in longitudinal adjustment
  covariate names, preventing distinct selected columns from collapsing into
  a different model specification.
- Reject empty rendered documents before replacing previous reports; retain
  explicit failure diagnostics and clean up temporary outputs for each format.
- Prevent response-derived and ordinal internal columns from being selected as
  adjustment covariates; document the complete reserved-name contract.
- Select GEE, CR2, and ANOVA sensitivity hypotheses by model terms rather than
  substrings in covariate names. Preserve the existing joint hypotheses.
- Use a relative covariance-rank threshold for measurement-unit invariant
  Wald inference, including very small finite outcome scales.
- Preserve Markdown pipe-table line breaks so DOCX renders scientific results
  as tables with intact cells; correct the unstratified missingness heading.
- Verify newly rendered report files before replacing previous outputs, escape
  literal scientific labels, and expose failed advanced-module diagnostics.
- Repair QTE argument-validation false positives and RNG cleanup; add independent
  data/prior, Friedman, paired-inference and transformation oracles, bounded
  generated cases, offline validation profiles and targeted mutation studies.
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
