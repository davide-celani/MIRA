# MIRA frequentist longitudinal handbook — implementation report

Date: 29 October 2026. Status: **implemented with specific validation limitations**.

Work was performed directly in `C:\Users\david\OneDrive\Desktop\MIRA\MIRA`,
starting from clean Git revision `0742f80f89d0f540d91821decfdfab38d979b6d8`.
The two supplied R-file references identify the authoritative repository files
themselves; there was no separate attached implementation to reconcile.
The requested tutorial collection is complete. External Quarto rendering and
complete pkgdown generation remain unverified for the environmental reasons below.

## A. Repository modifications

The ten files in the architecture and tutorial tables below were created.
They supply the handbook index and all nine requested scientific tutorials.
The remaining created or modified repository files are:

| File | Change and purpose |
|:--|:--|
| `DESCRIPTION` | Add `rmarkdown` to Suggests and `VignetteBuilder: knitr, rmarkdown`; retain the existing analysis/Stan dependency strategy. |
| `.Rbuildignore` | Exclude `_pkgdown.yml` and generated `docs/` from source bundles. |
| `.gitignore` | Ignore generated vignette HTML/R/asset directories, `inst/doc/`, and `docs/`. |
| `_pkgdown.yml` | Create minimal Bootstrap 5 configuration and order the articles by scientific progression. |
| `README.md` | Add links to all tutorials and installation/discovery instructions. |
| `NEWS.md` | Record the new handbook and documentation integration. |
| `dev/validate_frequentist_vignettes.R` | Create a reproducible installed-package render runner with an output-directory guard, status CSV and session record. |
| `dev/AI_USAGE.md` | Record this session's AI assistance and pending human review. |
| `dev/FREQUENTIST_VIGNETTES.md` | Preserve this implementation report in the repository. |

This is 19 repository files in total. `R/`, `man/`, `NAMESPACE`, `tests/` and
`inst/` retain their original tracked contents. No statistical behavior,
function signature, Bayesian implementation, Bayesian documentation or
Bayesian test was changed. Existing reference pages already contain detailed
argument contracts; the new tutorials complement them.

## B. Final architecture

```text
vignettes/
  MIRA.Rmd
  frequentist-long-getting-started.Rmd
  frequentist-long-data-and-detection.Rmd
  frequentist-long-descriptive-analysis.Rmd
  frequentist-long-treatment-groups.Rmd
  frequentist-long-missing-data.Rmd
  frequentist-long-multiple-outcomes.Rmd
  frequentist-long-advanced-models.Rmd
  frequentist-long-reporting.Rmd
  frequentist-long-complete-workflow.Rmd
```

All sources are directly under the standard directory and declare
`knitr::rmarkdown`, UTF-8 and `rmarkdown::html_vignette`. pkgdown has an index
group and a nine-article **Frequentist Longitudinal Analysis** group. Future
`bayesian-long-*` articles can be added as a separate group without moving
these files. There are no empty Bayesian articles or navigation entries.

## C. Tutorial-by-tutorial coverage

All ten articles rendered successfully against installed MIRA. The last complete
standalone render run took about 69 seconds in total. Times vary with concurrent
work. Each dataset is deterministic and generated locally with an explicit seed.

| Source | Audience, dataset and scientific content | Principal settings and outputs | Last render |
|:--|:--|:--|--:|
| `MIRA.Rmd` | All users; purpose, installation, dependency distinctions and learning sequence | Installed package version, vignette discovery and links | 1.66 s |
| `frequentist-long-getting-started.Rmd` | Beginners; 12-person walking-speed study, three visits | `mira_info_long()`, core analysis, classes, means/CIs, paired endpoint changes, trajectories and RDS preservation | 1.64 s |
| `frequentist-long-data-and-detection.Rmd` | Data preparers; 30-person score and 36-person pain/function examples | IDs, automatic/manual outcomes, time mappings, regex/function parsers, inspection, metric maps, ambiguous/missing/duplicate IDs and invalid types | 2.53 s |
| `frequentist-long-descriptive-analysis.Rmd` | Applied researchers; 64-person observational SBP cohort, four visits | Optional switches, pointwise CIs, paired changes, individual trajectories, correlations/pair counts, ICC estimators, outliers and plots | 3.60 s |
| `frequentist-long-treatment-groups.Rmd` | Comparative-study analysts; 96-person two-arm pain and 108-person three-arm function studies | Arm/reference selection, balance, visit differences, differences in change, Welch/rank tests, Hedges' g, Holm verification and plots | 2.58 s |
| `frequentist-long-missing-data.Rmd` | Analysts of incomplete studies; 82-person respiratory cohort and separate sparse instrument | Missing/non-finite availability, unequal groups, irregular scheduled visits, pairwise/complete-profile counts, adaptation/status and missingness assumptions | 2.51 s |
| `frequentist-long-multiple-outcomes.Rmd` | Analysts managing several endpoints; 72-person recovery cohort | Different schedules/scales, explicit/automatic outcomes, per-outcome directions/thresholds, safe extraction/error handling, threshold invariance and multi-outcome plotting | 3.53 s |
| `frequentist-long-advanced-models.Rmd` | Statistically trained readers; 96-person four-visit adjusted function study | One model bundle, numeric/categorical covariates, actual formulas, 13 module statuses, contrasts, diagnostics, sensitivity and model-comparison limitations | 18.16 s |
| `frequentist-long-reporting.Rmd` | Reproducible-report authors; 48-person rehabilitation study | Both input workflows, all report controls, metadata, sections, complete payload, CSV/manifest assertions, multi-outcome/configuration-only reports, overwrite and optional rendering | 2.92 s |
| `frequentist-long-complete-workflow.Rmd` | Experienced researchers; 120-person randomized fatigue study with observed-burden-dependent missingness | Inspection, unadjusted and adjusted inference, diagnostics, complete-profile comparison, plots, computed conclusions and verified report artifacts | 30.23 s |

Analysis tutorials use `mira_info_long()`; reporting and the complete workflow
also execute `mira_report_freq_long()`. Conclusions in the flagship article
are calculated from the result tables through inline R, including estimates,
pointwise intervals, participant counts and adjusted p-values.

## D. Public API coverage

The source-call audit in `api-coverage.csv` confirms code examples for all
**22 analysis parameters** and **21 reporting parameters**, including forwarded
analysis arguments through `...`. Positional data/result arguments are included.

Analysis: `data`, `id`, `time_vars`, `time_labels`, `arm`, `reference_arm`,
`arm_tests`, `alpha`, `plots`, `model`, `outliers`, `correlations`, `verbose`,
`p_adjust_method`, `improvement_direction`, `stable_threshold`, `strict_id`,
`outcomes`, `covariates`, `analyses`, `variable_pattern`, `inspect_only`.

Reporting: `x`, `data`, `...`, `output_dir`, `output_file`, `format`, `title`,
`subtitle`, `author`, `date`, `sections`, `include_complete_output`,
`extra_objects`, `max_table_rows`, `max_table_columns`, `digits`, `max_depth`,
`render`, `open`, `overwrite`, `quiet`.

Coverage means an argument is illustrated, not that every possible value and
validation branch is exhaustively tested. External `render = TRUE` examples
are opt-in; the public three-format rendering call was attempted but failed
before document generation. Browser opening is disabled in executable builds.

Key implementation distinctions are explicit: `analyses` overrides individual
switches; `"none"` retains core inference; declared outcome names prefix metrics;
the metric-name attribute resolves fields; module selection is global across
outcomes; time labels do not encode physical elapsed intervals; section/display
filters preserve eligible CSV exports, whereas `max_depth` can restrict traversal;
`include_complete_output` currently records a compatibility option without
altering visible content; the payload remains complete.

## E. Statistical execution and interpretation

Successfully demonstrated core methods include visit-specific summaries and
Student-t mean intervals, paired t and asymptotic signed-rank comparisons,
complete-profile Friedman tests, Kendall's W, paired Cohen's dz, between-arm
Welch and rank tests, Hedges' g, Pearson/Spearman correlations, pair counts,
ANOVA/model ICC and descriptive response classifications.

The advanced article's parsed HTML status table contains **13 successful
modules**: emmeans, repeated-measures ANOVA, Friedman, Kenward-Roger,
rank longitudinal inference, random slopes, three GEE correlations, two NLME
correlations, MMRM, and CR2. `module-status.csv` preserves their actual status
and warnings. The primary random-intercept model also fitted successfully.
Examples cover adjusted means and longitudinal contrasts, ML interaction
comparisons, GG/HF ANOVA corrections, REML/Satterthwaite MMRM with unstructured
covariance, unweighted relative marginal effects with ATS/WTS, cluster-robust
inference, model fit criteria, effect sizes and sensitivity views.

The articles distinguish visit differences, differences in change and joint
interaction hypotheses; pointwise confidence intervals and adjusted p-values;
conditional and population-averaged models; response classification and
clinical validation; and available-pair, complete-profile and model samples.
The missingness tutorial does not claim that its joint masking process is MAR.
The flagship generator deliberately uses fully observed screening burden for
MAR dropout, while explaining why an unadjusted MCAR rank analysis is not an
equally justified primary analysis for that generator.

## F. Validation evidence and commands

Environment: Windows 10 x64 build 19045, R 4.4.3, MIRA 0.0.1. Full dependency
versions and session information are in `evidence/`. Initial installation with
R 4.5.3 failed because that library lacked eight existing Imports (`afex`,
`clubSandwich`, `cmdstanr`, `geepack`, `lmerTest`, `mmrm`, `nparLD`, `quarto`).
The existing R 4.4 library supplied them; no dependency was newly installed.

| Validation | Observed outcome |
|:--|:--|
| Clean source installation and namespace loading | Passed |
| All ten standalone vignette renders, including assertions and figures | Passed |
| Extra advanced sensitivity branch | Passed with `MIRA_RUN_ADVANCED_VIGNETTES=true` |
| Public parameter/cross-reference audit | All 43 parameters illustrated; every local HTML article link resolves |
| Source-only reports | QMD/runtime/RDS/manifest/CSV existence, payload identity, archive independence, overwrite handling and multi-outcome exports verified |
| Existing complete source test suite | **11,912 pass, 0 failure/error/warning, 6 skips** |
| Full `R CMD check --as-cran --no-manual` | **0 errors, 0 warnings, 2 notes**; package tests independently confirm the same 11,912 passes and six skips |
| Final snapshot `R CMD check --as-cran --no-manual --no-tests` | **0 errors, 0 warnings, 1 note** (unable to verify current time); all final vignettes rebuilt successfully |
| `R CMD build` and installed discovery | Passed; all ten entries and HTML sources installed; generated index also present |
| Bundle identity | All 19 R/source-vignette files match the repository; all ten installed vignette sources match |
| Quarto HTML/PDF/DOCX | Attempted; all blocked by Windows process-spawn/pipe errors |
| pkgdown | Configuration parses and discovers ten ordered articles; full initialization blocked by CDN SSL failure |

The full check's two notes were inability to verify current time and failure
to find Pandoc on PATH for checking README/NEWS. Pandoc was available to
rmarkdown through `RSTUDIO_PANDOC`. The final check adds its directory to PATH.
The final snapshot check skips re-running tests already executed twice against
unchanged R/test sources, while rebuilding all final vignettes and checking
metadata, reference documentation, examples and installation. It includes the
last display/prose edits and the current manual's two-package builder declaration.
Neither check builds the PDF reference manual (`--no-manual`). Remote CRAN
incoming checks were disabled explicitly; network time verification remains on.

The six ordinary-suite skips are one opt-in CmdStan sampling profile and five
opt-in Quarto integration profiles. No analysis engine was skipped for absence
of an installed dependency. The direct rendering attempt is additional to these
ordinary-suite skips and is recorded as a failure, not a successful render.

Relevant actual commands (PowerShell; quoted executable paths are literal):

```powershell
$env:LC_ALL='C'
$env:LANG='C'
$env:R_LIBS='C:\Users\david\Documents\Codex\2026-10-29\files-mentioned-by-the-user-mira\work\library'
$env:RSTUDIO_PANDOC='C:\Program Files\RStudio\resources\app\bin\quarto\bin\tools'
$env:_R_CHECK_FORCE_SUGGESTS_='true'
$env:_R_CHECK_CRAN_INCOMING_REMOTE_='false'
$env:MIRA_RUN_CMDSTAN_TESTS='false'

# From the chat workspace:
& 'C:\Program Files\R\R-4.4.3\bin\R.exe' CMD INSTALL --library=work/library --no-multiarch 'C:\Users\david\OneDrive\Desktop\MIRA\MIRA'
& 'C:\Program Files\R\R-4.4.3\bin\Rscript.exe' --vanilla work/run-tests.R 'C:\Users\david\OneDrive\Desktop\MIRA\MIRA' 'C:\Users\david\Documents\Codex\2026-10-29\files-mentioned-by-the-user-mira\work\validation'

# From the repository:
& 'C:\Program Files\R\R-4.4.3\bin\Rscript.exe' --vanilla dev/validate_frequentist_vignettes.R 'C:\Users\david\Documents\Codex\2026-10-29\files-mentioned-by-the-user-mira\work\validation\final-vignettes'

# From work/ (full check) and work/final-build/ (final build/check):
& 'C:\Program Files\R\R-4.4.3\bin\R.exe' CMD build 'C:\Users\david\OneDrive\Desktop\MIRA\MIRA'
& 'C:\Program Files\R\R-4.4.3\bin\R.exe' CMD check --as-cran --no-manual MIRA_0.0.1.tar.gz
$env:PATH='C:\Program Files\RStudio\resources\app\bin\quarto\bin\tools;'+$env:PATH
& 'C:\Program Files\R\R-4.4.3\bin\R.exe' CMD check --as-cran --no-manual --no-tests MIRA_0.0.1.tar.gz

# Final bundle installation, from the chat workspace:
& 'C:\Program Files\R\R-4.4.3\bin\R.exe' CMD INSTALL --library=work/library --no-multiarch work/final-build/MIRA_0.0.1.tar.gz
```

The additional model run set `MIRA_RUN_ADVANCED_VIGNETTES=true` and rendered
the two advanced/workflow drafts using the same installed public package.
The corresponding log is retained. `audit-documentation.R` parses every
executable source call and the generated HTML, checking links and the 13-row
module table. Reports were generated only under temporary/workspace paths.

The public external-render probe used:

```r
set.seed(1)
d <- data.frame(id = 1:12, score_t0 = rnorm(12, 40, 5),
                score_t1 = rnorm(12, 38, 5))
a <- mira_info_long(d, analyses = "none", verbose = FALSE)
b <- mira_report_freq_long(a, output_dir = "<workspace>/work/render-probe",
  format = c("html", "pdf", "docx"), render = TRUE,
  open = FALSE, overwrite = TRUE, quiet = FALSE)
```

Here `<workspace>` was the absolute chat workspace given in the environment.
`QUARTO_PATH` pointed to RStudio's `quarto.exe`. A separate direct CLI attempt
with `quarto.exe render ... --to html --execute-dir ...` also failed spawning R.
The actual captured errors are preserved in `evidence/`.

For pkgdown, the executed public call was `pkgdown::build_articles()` with the
repository as `pkg`, `lazy = FALSE`, `quiet = FALSE`, `preview = FALSE` and
`override = list(destination = "<workspace>/work/pkgdown-articles")`.
`R_USER_CACHE_DIR` was redirected to `<workspace>/work/cache` after the default
cache path proved unwritable. Downloading standard Headroom JavaScript then
failed with an SSL connect error. Certificate verification was retained.

## G. Outstanding issues and validation limits

1. **External report rendering:** Quarto's R interface raises Windows system
   error 5 while creating a processx write pipe. Direct CLI execution also
   reports an invalid handle spawning R. No completed HTML, PDF or DOCX report
   is claimed. Source-only bundles and ordinary vignette HTML are validated.
   TinyTeX was found, but the PDF path failed before its use.
2. **Website generation:** pkgdown's config and article discovery pass; the
   complete website build needs a working TLS/CDN connection or valid cached
   standard assets. Site layout has not been certified.
3. **Existing sparse plotting defect:** an outcome with some finite values
   but no available baseline-to-follow-up pair can error while constructing
   empty change-plot data. The missing-data article documents it and requests
   only model/correlation diagnostics for the sparse example. Reproduction:

   ```r
   sparse <- data.frame(participant_id = sprintf("S%02d", 1:6),
     instrument_t0 = c(12, rep(NA_real_, 5)),
     instrument_t1 = rep(NA_real_, 6), instrument_t2 = rep(NA_real_, 6))
   mira_info_long(sparse, id = "participant_id", outcomes = "instrument",
     covariates = character(0), analyses = c("model", "correlations", "plots"),
     verbose = FALSE)
   ```

   Observed error: `arguments imply differing number of rows: 0, 1`.
   Removing `"plots"` permits the documented diagnostic container. The engine
   was left unchanged; no undocumented statistical adjustment was introduced.
4. **Platform scope:** validation used one Windows/R 4.4.3 environment and
   existing dependencies, not fresh provisioning, R 4.2/R 4.5, Linux/macOS,
   remote CI or the entire version matrix. CmdStan sampling was not requested.
5. **Presentation scope:** HTML tables, structure, unintended math and links
   were inspected programmatically and editorially. Browser screenshots were
   unavailable because the browser policy blocks local `file:` URLs; no
   pixel-level browser certification is asserted. Pandoc emits an external
   `--highlight-style` deprecation warning during vignette rendering.
6. **Scientific scope:** synthetic examples do not validate clinical thresholds,
   identify causal effects or guarantee missing-data ignorability. Independent
   statistical/editorial review and real-study validation remain necessary.

## H. CRAN and JOSS assessment

The documentation follows standard flat vignette discovery, registers its
engines/dependencies, uses portable temporary paths, avoids runtime downloads,
keeps external report rendering optional, and builds/rebuilds using R's own
tools. R Core's current development manual also identifies both knitr and rmarkdown for
this builder configuration. See [Writing R Extensions](https://www.stat.math.ethz.ch/CRAN/doc/manuals/r-devel/R-exts.html#The-DESCRIPTION-file).
CRAN's policy emphasizes portable code, limited checking resources and graceful
network behavior. The final package still needs other-platform checks and
resolution/explanation of environmental notes before submission; these local
results do not establish CRAN acceptance. See the [CRAN Repository Policy](https://mirror.metanet.ch/cran/web/packages/policies.html).

The handbook supplies a scientific audience, reproducible use cases, API
explanations, objective execution checks, assumptions and limitations. Those
address documentation/functionality aspects of the [JOSS review criteria](https://joss.readthedocs.io/en/latest/review_criteria.html)
and [review checklist](https://joss.readthedocs.io/en/latest/review_checklist.html).
License/citation/contribution/support files already exist. This task did not
establish research adoption, independent-user experience or manuscript readiness.

Current [JOSS submission guidance](https://joss.readthedocs.io/en/latest/submitting.html)
requires sustained public development exceeding six months and demonstrated
research use. The non-shallow local history begins on 19 August 2026, less
than six months before this work; public availability and any earlier history
would need separate evidence. A paper and verified research-impact account
remain separate deliverables.

JOSS also requires disclosure of AI tools/versions where known, their scope,
and the authors' confirmation of human review and validation. This session's
record describes OpenAI Codex assistance without inventing an exposed model
identifier. Human verification is pending and cannot be supplied by this agent.
See the [JOSS AI policy](https://joss.readthedocs.io/en/latest/policies.html#ai-usage-policy).

## I. Final implementation status

**Implemented with specific validation limitations.** The complete frequentist
longitudinal collection is written, executed, built and discoverable; its
scientific engine and Bayesian component are preserved. The remaining blockers
concern external rendered reports, full pkgdown generation, environmental check
notes, the disclosed existing sparse plotting defect, and validation beyond
this local platform. This is a documentation implementation outcome, not formal
CRAN/JOSS compliance or a substitute for human scientific review.

The companion chat output directory,
`C:\Users\david\Documents\Codex\2026-10-29\files-mentioned-by-the-user-mira\outputs`,
contains the linked HTML handbook and Rmd sources, complete parameter/module
CSVs, selected evidence logs, and the source package built from this working
tree. Relative evidence filenames in this report refer to that directory.
Repository edits remain local changes.
