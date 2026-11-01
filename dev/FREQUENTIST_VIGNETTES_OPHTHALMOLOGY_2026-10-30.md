# Ophthalmic frequentist vignette rewrite — 30 October 2026

All eleven frequentist R Markdown articles have been rewritten around synthetic
ophthalmic studies and executed again. New generators, measurement scales,
covariates, clinical questions and interpretations replace the earlier generic
examples. Statistical implementation, public APIs, tests, dependencies and
Bayesian/Stan material are unchanged. Earlier local changes and reports were
preserved. This is a completed documentation revision with recorded validation;
it is not a guarantee of perfection or an attestation of human clinical review.

## Read and reproduce

- Open `handbook/MIRA.html` for the linked HTML handbook. Its eleven Rmd sources
  are beside the HTML files and match the repository and installed package.
- `MIRA_0.0.1.tar.gz` is the newly built source archive, installed and checked
  during this phase. `bundle-source-identity.csv` records its source identity.
- `clinical-case-map.csv` describes every main study and the named edge and
  sensitivity cases. `article-summary.csv`, `api-coverage.csv` and
  `scientific-coverage.csv` connect the teaching material to fresh evidence.
- `CLINICAL_SOURCE_RECORD.md` links the primary and official measurement sources.
- `review/final-clinical-reader-review.md`, `html-audit/` and `figures/` contain
  independent AI reader review, generated-DOM checks and decoded plot images.
- `reproducibility/` contains the scripts used. `minimal-report-example/` is an
  additional synthetic BCVA **source-only** report bundle, independently checked.

The repository copy of this report is
`dev/FREQUENTIST_VIGNETTES_OPHTHALMOLOGY_2026-10-30.md`. The earlier implementation,
revision and closure reports retain their historical evidence. They do not
describe the new clinical datasets or supply this phase's execution results.

## Clinical redesign

| Article | New question and dataset |
|---|---|
| MIRA index | Ophthalmic reading routes; measurement definitions, direction, study-eye selection and bilateral-design limits. |
| Getting started | Twelve participants, one prespecified eligible right eye each: rounded BCVA ETDRS letters at baseline and week 4, an absent follow-up copy, then an added week 12. |
| Data and detection | Separate 30-person and 36-person retinal BCVA/import cases, defined OCT CMT fields and a separate four-person logMAR direction example. Malformed imports remain explicit diagnostic cases. |
| Descriptive analysis | Sixty-four people with glaucoma or an elevated-pressure history, one study eye each; IOP in mmHg over baseline/months 2, 5 and 9. Visit counts are 64/62/61/59, with 56 complete profiles. |
| Treatment groups | New 96-person two-strategy and 108-person three-strategy glaucoma IOP cohorts. Raw visit differences, differences in changes and reference reversal remain distinct. |
| Missing data | Eighty-two retinal-edema participants with OCT central-subfield thickness in micrometres; complete, independently masked, complex masked and missing-group copies, plus a sparse six-person edge case. |
| Multiple outcomes | Seventy-two retinal participants with BCVA letters, defined OCT CMT and IOP safety measurements, retaining different visit schedules and zero/one-pair endpoint cases. |
| First model | Eighty participants with diabetic macular edema; rounded study-eye BCVA at baseline/months 3 and 6, model-vs-paired contrast checks and conditional residual diagnostics. |
| Advanced models | Ninety-six diabetic-macular-edema participants, four BCVA visits, age, lens status and separately observed pre-enrolment screening CMT. Missing-age and covariate-free sensitivity copies are explicit. |
| Reporting | Forty-eight diabetic-macular-edema participants, three BCVA visits, plus a same-eye BCVA/CMT multi-report. Results, metadata, display controls, archive depth and overwrite contracts are checked. |
| Complete workflow | A new 120-person diabetic-macular-edema study with strategies A/B, BCVA at weeks 0/4/8/12, observed-covariate-dependent missingness, adjusted primary contrast and reversed week-8 transfer. |

Each analysis uses one study eye per person, selected before follow-up under
the stated hypothetical protocol. Person/eye mappings are explicit in the
model studies. Separate eye identifiers do not establish independence when
both eyes of a real person contribute. The material does not recommend
discarding fellow-eye information for all ophthalmic designs or bypassing
strict identifier checks.

BCVA is rounded ETDRS letters, with higher values favorable; the separate
logMAR import explicitly uses lower values. No universal letter/logMAR
conversion is asserted. In these examples `cmt` means mean OCT retinal
thickness in the central **1-mm-diameter subfield**, locally synonymous with
CST, measured in micrometres under a consistent fictional instrument and
segmentation protocol. A thickness decrease is anatomical change and does
not establish visual benefit. Lower IOP is a favorable teaching direction
only in the elevated-pressure scenarios. Retinal IOP safety uses `unknown`;
the name-based auto-lower demonstration is explicitly rejected for that task.

All values, effects, variances, schedules and missingness probabilities are
invented. Rounding and fixed-seed range checks do not constitute clinical
eligibility limits; no clipping is used to hide implausible draws. Gaussian
analyses of letter scores are stated as working models. Five letters,
25 micrometres and one/two mmHg are pedagogical classification bands where
configured, not validated MCIDs, efficacy thresholds or treatment decisions.

## Fresh execution and scientific checks

The phase uses R 4.4.3 and an isolated installation under
`work/ophthalmology/library`. The environment probe found all Imports available.
R 4.5.3 was also probed and lacks eight Imports plus the optional `pdftools`; it was not used to claim a
successful full validation. No dependencies or toolchains were installed.

The final runner completed **33/33 processes**: HTML rendering, extracted R
and all visible evaluated R chunks for each of the eleven exact sources.
Every mode ran in a new R process with the recorded library and source MD5.
The advanced opt-in branch was enabled; external report formats were disabled
for this ordinary matrix. The process-time sum was **394.18 seconds**.
`handbook/validation/vignette-status.csv`, its per-mode contexts/logs and
`evidence/final-freeze-sha256.csv` identify those exact inputs.

The unchanged test suite ran anew: **387 cases, 11,912 passed expectations,
zero failed/error/warning expectations, six explicit skips**, in **408.08
seconds**. The skips cover opt-in CmdStan/Quarto integrations; no skipped test
is reported as passing. The CSV and session record are in `evidence/`.

Raw-data checks retain finite visit, pair, complete-profile and fitted-key
denominators; independent reconstructions of means, variances, quantiles and
confidence intervals;
within-outcome multiplicity families; independent Welch algebra; sign and
reference reversals; metric-name-safe extraction; and full-object/archive
identity despite bounded report display. Intentional error fixtures are
labelled separately from valid analyses. `evidence/defensive-fixtures.csv`
records error, absent-table and nonestimable-interval fixtures; the validator
also asserts duplicate endpoint rejection and multi-outcome dispatch.

The advanced article requests thirteen real modules: **nine complete and four
partially usable**. The random-slope fits exist but their native comparison
rejects different data-call symbols. Three GEE fits exist but their native
ANOVA components fail to find `geese.fit` with namespace-only loading. These
states are printed and interpreted; they were not patched, hidden or called
complete. Sparse plotting without consecutive observed pairs remains an
explicitly reproduced native limitation.

The new independent native probe verifies afex's actual 96-eye sample, exact
source covariates and 96-by-4 response matrix, its expanded mean formula
including `arm_factor:lens_status`, the `~ time_factor` within design and
rank-three sum-contrast transformation. Native Friedman counts agree with
the raw complete profiles. This does not equate afex's hypotheses with the
additive mixed model. The NLME adapter fits `nlme::lme` with a random intercept
and residual covariance structure; a separate standalone GLS fit is not
claimed. The probe script, compact proof and log are retained in the new phase.

The complete workflow's model uses **440 measurements from 120 study eyes**.
The final raw paired endpoint uses **98 eyes** (48 Strategy A, 50 Strategy B);
22 participants contribute model information without that paired endpoint.
The complete-profile sensitivity has 84 people/336 measurements. These samples
and targets are not interchangeable.

- Raw week-12 difference in BCVA changes, **B minus A**: **+6.42 letters**, with
  pointwise Welch 95% interval **[4.2001, 8.6399]**.
- Adjusted primary difference in baseline-to-week-12 changes, **B minus A**:
  **+6.1447 letters**, pointwise 95% interval **[4.0897, 8.1997]**;
  Holm-adjusted p = **5.1288e-08** across six arm-by-visit-pair contrasts.
- Executed week-8 transfer, **A minus B**: **−6.3837 letters**, pointwise 95%
  interval **[−8.3830, −4.3844]**, Holm-adjusted p = **6.6387e-09**.

The model-based primary contrast is checked against all four native marginal
means and the explicit change algebra; the reversed transfer recomputes its
own endpoint and direction. A positive B-minus-A difference in BCVA gains
favors B under the declared higher-is-better scale. These generated effects provide teaching examples,
not evidence about a real intervention. All quoted values come from the new
final workflow ledgers, not earlier rehabilitation examples.

## HTML, package and archive verification

Independent AI review read **all 71 generated tables and 1,665 cells** and
viewed **all nine decoded PNGs**. Nine visible `div.figcaption` captions were
verified against their image alt text. An initial overly narrow caption
selector was corrected in the audit; no document change was necessary.
The generalized HTML audit reports **zero errors, zero warnings and ten
informational NULL contexts**, each corresponding to an explicitly inspected
diagnostic/status field. Plot axes, units, uncertainty, signs and denominators
were checked. Exact source/HTML/PNG hashes accompany the review.

Source build, archive installation and `R CMD check --as-cran --no-tests`
completed. The PDF and HTML reference manuals, examples, installed vignette
discovery and vignette rebuilding passed. Tests were run separately as above,
then reconciled with the archive's unchanged R/test sources. Build time was
155.38 seconds; check time 302.58 seconds. Result: **zero ERRORs, zero WARNINGs,
one NOTE: unable to verify current time**. Remote incoming checks were disabled;
the console also records unavailable repository-index access. This is local
validation, not a CRAN acceptance decision.

The archive has **98 checked source files**: 97 byte-identical and DESCRIPTION
semantically identical apart from the standard build-added fields. The
archive is **864,151 bytes**; compiled documentation is **1,423,707 bytes**.
Archive SHA256: `920b3c47ef9aa8d03312dba67746c3fd13735f020caa8fb8906091a1420d0d65`.
All eleven installed Rmds match repository and delivered sources, with twelve
installed HTML files including R's generated index. The generated R-help index
routes are verified by static mapping, not by treating them as standalone
`file://` links or starting a help server. Delivery links and exact hashes are
recorded under `evidence/`; all **143** source/flat/installed link checks pass.

## Explicit limits and scope

The ordinary vignette HTML, extracted code and source-only report bundles
pass. Fresh opt-in attempts to render the reporting vignette's HTML/DOCX/PDF
reports all fail at Windows `processx` write-pipe creation, system error 5
(access denied). A separate direct public `render = TRUE` HTML/DOCX attempt
records the same failure. Source files are retained and each failed format
is labelled. No successful external report rendering is claimed. The PDF
attempt disables LaTeX auto-installation; no external installer was invoked.

The environment's security policy blocked local-file browser inspection in
the preceding closure phase, with no permitted browser/server workaround.
It was not retried here. This phase inspected actual generated HTML/DOM,
all table cells, decoded PNGs and visible caption nodes. Browser viewport,
responsive layout and interactive usability have not been certified.

Only the eleven vignette sources, the focused README/NEWS/pkgdown article
group description, the appended AI-assistance record and this new technical
report changed in the package during this phase. The developer runner,
statistical R, Rd, tests, `inst`, dependencies, NAMESPACE, Bayesian/Stan
material, earlier reports and pre-existing local changes were preserved.
`evidence/scope-comparison.csv` compares every file with the phase-start
snapshot, including ignored files; no file was removed. Final source,
execution, bundle, installed and delivered identities are reconciled by the
delivery audit and SHA256 manifest.

No commit, staging, push, PR, release, publication or pkgdown site build was
performed. The redesign and cross-reviews were AI-assisted. Human
ophthalmologic, statistical and editorial review, and responsibility for any
accepted scientific or publication claims, remain with the maintainers.
