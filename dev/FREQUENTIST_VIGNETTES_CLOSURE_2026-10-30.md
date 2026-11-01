# Frequentist longitudinal handbook closure — 30 October 2026

Status: **Implemented with explicit validation limitations** (the requested
Italian state: **Implementato con limitazioni di validazione esplicite**).

The two-function documentation assignment is implemented as an executed,
auditable candidate for human scientific and editorial review. The final
collection contains eleven articles, including one substantive new first-model
bridge. This status is scoped to the local environment and the taught paths.
It is not JOSS approval, CRAN eligibility, human peer review or validation of
an arbitrary real-world analysis.

The working repository is `C:/Users/david/OneDrive/Desktop/MIRA/MIRA`.
Deliverables are rooted at
`C:/Users/david/Documents/Codex/2026-10-30/ese/outputs`.
Artifact paths below refer to that delivery root, including when this same
report is retained in repository `dev/`. Historical reports were read and
preserved; none has been rewritten to claim this phase's results.

## Starting state, authorization and preservation

The starting HEAD was `0742f80f89d0f540d91821decfdfab38d979b6d8`.
The working tree already contained modified DESCRIPTION, NEWS, README,
AI_USAGE and ignore files, plus untracked vignette/config/validation material.
It was not clean. Before edits, Git status, binary diff, untracked names,
HEAD and SHA256/lengths of all non-Git files were recorded in `evidence/`.
The prompt, historical report, API matrix, HTML collection, logs and
reproduction scripts from 29 October were available and consulted.

New work is identified against that observed state, not against an assumed
clean checkout. The preservation proof is `evidence/scope-comparison.csv` and
its JSON summary. Statistical R code, public signatures, NAMESPACE, man,
tests, Stan/Bayesian material, DESCRIPTION and ignore rules retain their
initial bytes. Only the frequentist articles, navigation and necessary
developer documentation/runner were changed or added. Existing local work
was not reset, restored or staged. No commit, push, PR, tag, release,
submission, deployment, external message or pkgdown build was performed.

There are **15 modified and two added files** relative to the captured initial
tree, with no removed or unexpected files. All 87 protected implementation,
reference, test, metadata and ignore-rule inputs are unchanged. The complete
per-file before/after list is `evidence/phase-changes.csv`.

| Changed support or source group | Purpose |
|:--|:--|
| Ten existing vignette Rmds | Preserve useful studies; repair progression, visibility, samples, extraction and archive explanations. Article-by-article details are in article-summary.csv. |
| New vignettes/frequentist-long-first-model.Rmd | Substantive bridge between paired change and adjusted model reference. |
| README.md | Add the first-model source link and reading-route entry. |
| NEWS.md | Record this documentation closure without rewriting historical entries. |
| _pkgdown.yml | Register the new article statically; no site build. |
| dev/validate_frequentist_vignettes.R | Isolated processes, architecture-derived collection, output guard, exact hashes, code modes, flat HTML and ledgers. |
| dev/AI_USAGE.md | Add a dated scope/validation-limit disclosure; earlier records preserved. |
| New dev/FREQUENTIST_VIGNETTES_CLOSURE_2026-10-30.md | Retain this new technical report separately from historical reports. |

## Diagnosis and final learning architecture

The existing collection contained sound explanations and valuable deterministic
studies. The main barriers were progression: three visits before the first
result was understood, required helpers hidden in setup, all adjusted model
machinery introduced together, and a first report call containing many controls.
The revision preserves the useful studies and scientific audits, while moving
conventional configuration/ID decisions ahead of optional parsing references.

| Route | Prerequisites and reading order | Intended stopping point |
|:--|:--|:--|
| Essential | Basic R data frames; Getting started; conventional Data and detection sections when importing data | A correct wide table, secure ID/order, mean versus paired change, one figure, RDS and a minimal source-only bundle. |
| Applied | Essential concepts; Descriptive analysis, Treatment groups, Missing data, Multiple outcomes as needed, Reporting | Named denominators/subtractions/families, missing/partial outputs, outcome settings and manifest-led archives. |
| Advanced | Regression plus group/missingness concepts; First model, Advanced models, Complete workflow | Model dependence and rows before adjustment, coherent contrasts, assumptions/alternatives, partial status, transfer and correct archived primary result. |

Getting started now changes two visits to one missing value and then adds a
third visit without silently replacing participants or earlier measurements.
Its first-result bridge defines means, pairs, subtraction, interval and p-value
before a method survey. Its first archive reuses the understood analysis with
a short source-only call. Reporting introduces metadata and controls later.
The new first-model article uses an explicitly new, adequate 80-person study
to explain random intercepts, baseline as response and one contrast. The real
public model bundle is requested; no nonexistent single-model selector is
invented. Advanced models then introduces a new adjusted two-arm question.

All eleven sources remain directly under `vignettes/`. The existing names
are stable; `frequentist-long-first-model.Rmd` is the only new article.
The index derives the expected validation collection. Static pkgdown links
include it; future Bayesian articles and R's generated installed index remain
possible without creating empty Bayesian material.

`article-summary.csv` provides audience, question, data, methods, visible
outputs, final execution times and limitations for every article.
`case-matrix.csv` records eighteen motivated capability decisions, novelty,
complexity, location, verification and limits. The additional cases address
reference reversal, independent missingness, one missing group, one missing
predictor, zero/one endpoint pairs and a capstone adaptation task. They are
not repeated fits added to reach a target count.

## Public contracts and scientific coverage

`api-coverage.csv` has all **43 current public formals**: 22 analysis and 21
reporting parameters. Each row records default, accepted types/values,
precedences, explanation, exact source/chunk examples, explicit branches,
ordinary/optional execution, numerical or contract oracle and unverified limits.
All defaults and all eleven final AST input hashes were cross-checked.
The AST inventories 61 public calls and 477 explicit argument records; those
counts are inventories, not evidence that every branch was reached.
`report-api-contract.md` gives a compact interpreted reporting reference.

Key distinctions retained include: non-NULL analyses overrides five switches;
none keeps core inference/Friedman; explicit declarations affect owned metric
names while native engine fields follow native contracts; multi-outcome class
must be checked first; module selection/covariates are global; direction and
threshold can be per outcome. Permissive IDs do not deduplicate participants.
`analyses=c("all","none")` follows the implementation's all-first expansion;
the unusual combination is documented in the matrix, not promoted as normal use.

`scientific-coverage.csv` records estimands/hypotheses, actual structures,
samples/units, subtraction, uncertainty, families, assumptions and locations.
The manual distinguishes observed means, paired changes, arm differences,
differences in changes, adjusted visit contrasts and joint hypotheses. It
does not infer model n from wide nrow or compare within-arm significance.
Pointwise intervals remain pointwise when p-values are adjusted. Families
are separate for paired tests, arm values/changes, omnibus tests, missingness,
post hoc and model contrasts; no study-wide or cross-outcome correction is
implied. Holm/Bonferroni and BH/BY retain their distinct error criteria.

Independent checks use raw finite values, paired participant identities,
participant–visit keys, raw means/differences, elementary t/Welch intervals,
ANOVA mean squares for ICC, actual finite-p adjustment and response thresholds.
The within-person SD uses all finite person-centered residuals with the
ordinary SD denominator, not an assumed pooled residual degrees of freedom.
Preferred ICC uses finite model ICC and otherwise the complete-profile ANOVA
estimate. Correlation-one with an exact 10-mmHg shift is an arithmetic example
of why association does not establish agreement.

Time labels do not supply physical or irregular individual time. Baseline is
a response in the repeated-response models. Automatic covariates satisfy
structural rules rather than a confounding strategy. The afex adapter expands
between factors and time-by-factor/covariate terms, unlike the additive
covariate mean structures of the other wrappers. KR/CR2/GEE/MMRM TIME and ARM
joint tests include the relevant interaction coefficients. Rank RTE/ATS/WTS
do not estimate Gaussian mean effects, and the implemented rank wrapper's
MCAR condition is explicit. Alternatives are not votes for a preferred p-value.

## Coherent primary workflow and transfer

The capstone declares
`(Program Week12 - Program Baseline) - (Care Week12 - Care Baseline)`
in fatigue points. Negative values favor the program. Its unadjusted
observed-pair estimate is -5.0070 points, CI [-7.1317, -2.8824], from 49 care
and 53 program pairs. The adjusted primary estimate is **-5.0579**,
pointwise 95% CI **[-7.0078, -3.1080]**, Holm p **3.4459e-06** within six
arm-by-visit-pair contrasts. Native double-negative labels are reconciled
through the four model-based means. This is the result extracted, interpreted
and preserved; a final-visit arm difference is not substituted for it.

The model uses 443 measurements from 120 people; 102 people have both primary
endpoints and 88 have complete four-visit profiles. Exact keys, not merely
counts, are audited. The distinct complete-profile result is labeled as such.
The executed transfer checkpoint asks the reader to use Week 8 and reverse
the subtraction: care change minus program change is **4.2289**, pointwise
CI **[2.3134, 6.1445]**, Holm p **9.4390e-05** in the existing six-contrast
family. It verifies the manifest-linked native CSV without changing the
primary endpoint or claiming new independent treatment evidence.

Four-cell contrast checks share the fitted model/emmeans engine and do not
independently validate its covariance or CI implementation. The raw unadjusted
Welch audit has a separate algebraic basis. Synthetic generators establish
their own mechanism by construction; they do not establish missingness or
causal identification in real data. The respiratory joint mask is not labeled
MAR. Its independent-mask variant is MCAR by construction; the capstone's
MAR masking uses an always-observed screening predictor.

## Modules, samples and known partial results

The advanced complete-data configuration has a primary model plus thirteen
containers with performed=TRUE. **Nine are completed and four are partially
usable**. The ledger retains requested context, package/version, warnings,
nested error text and interpretation. A top-level NULL error is not treated
as a clean bill of health.

| State in the final advanced example | Components |
|:--|:--|
| Completed | emmeans, RM-ANOVA, Friedman, KR, rank longitudinal, NLME CS/AR1, MMRM, CR2 |
| Fit available, comparison unavailable | Random slope: native LRT refuses distinct data-call symbols despite equal observations. |
| Fit/joint tests available, native ANOVA unavailable | Independence, exchangeable and AR1 GEE; namespace-only geepack 1.3.13 cannot resolve geese.fit. |

Final ledgers are under `handbook/validation/frequentist-long-advanced-models/render/`.
Model/sample/formula rows preserve 384 measurements from 96 people; RM-ANOVA
and Friedman counts have participant units. The focused missing-age variant
preserves all responses and descriptive baseline n=96 but has 95 model
participants/380 rows, verified with exact keys. It does not impute age or
establish ignorability. A separate native sample/formula probe records actual
afex structure and retained participant identities.

Nonrequested, insufficient-information, design-inapplicable, dependency/version,
calculation-error, partial and indeterminate states are explained separately.
The multi-outcome extractor distinguishes absent table, no endpoint pairs,
one pair, finite available contrast and duplicate keys. Explicitly constructed
defensive fixtures are tested separately in
`reproducibility/validate-defensive-extraction.R`; they are never represented
as ordinary MIRA analytic outputs or generated through a bug.

Native minimal reproductions reconfirm sparse plotting, random-slope data-call
comparison and GEE ANOVA limitations. Linking geepack in a diagnostic process
recovers a sequential ANOVA, not a replacement for MIRA's joint hypotheses or
a retroactive correction. `SOFTWARE_FOLLOWUP.md` gives separate concrete fixes
and regression proposals. No engine was patched to make a tutorial pass.

## Reporting and archive integrity

The path progresses from stored-result source-only reporting to direct data,
metadata, sections, controls, archive depth, overwrite and optional rendering.
Sources/runtime/payload/manifest and eligible CSVs are verified. Presentation
limits and digits do not truncate the numeric archive; max_depth can constrain
traversal/export while preserving the payload. include_complete_output is a
recorded compatibility flag without a current visible-content difference.
Extras are retained/appendix-visible when selected and excluded from CSV.
One independent directory per bundle is necessary even with different stems.
Writes are not fully transactional, and stale rendered documents may survive
a failed attempt. expected_files is a plan, not success evidence.

Before serialization, input identity is checked. After round trip, selected
tables, name maps, settings, native parameters and covariances are checked;
whole fitted-model identity and recursive pointer/environment equality are
not used as the integrity criterion. The manifest locates actual CSVs;
CSV does not preserve all native classes/attributes. Payloads may contain
data and identifiers. Automatic descriptive conclusions are not inference.

`minimal-source-bundle/` is a small retained synthetic example with a QMD,
runtime, payload and manifest-linked tables. Its tables/settings are separately
verified from the deterministic raw study. It contains no rendered HTML.
An initial standalone verification used the wrong payload key, was corrected
to result, and both the corrected creator and read-only delivered-bundle
validation pass. A cleanup attempt was automatically refused by sandbox
approval policy; no deletion occurred and the existing bundle was preserved
and validated without bypassing the restriction.

Real public render=TRUE dispatch was attempted for HTML/DOCX. The final
reporting vignette's separate deferred branch attempted HTML/PDF/DOCX after
disabling automatic LaTeX package installation in a QMD copy. **All formats
failed before rendering at Windows processx write-pipe creation, system error
5 (access denied).** Exact final reporting source MD5 is in the deferred
context. No external finished document is claimed. The source-only descriptor
does not update after external rendering. open remains FALSE throughout.
Ordinary rmarkdown vignette HTML and R reference manuals are different paths
and succeed. No CLI/toolchain installation or security bypass was attempted.

## Independent AI reading and generated-output inspection

Three agents had nonoverlapping authoring responsibilities; other agents
then read their reader-visible material before using code as an explanation.
The 12-checkpoint walkthrough covers preparation, identity/order, execution,
denominators, means/subtraction, CI/p-values, group tests, model contrasts,
partial results, mapped extraction, manifest archive and applicability limits.
Concrete corrections include visible libraries/helpers, the first-model bridge,
guarded dependency-version lookup, explicit reference-grid averaging, same-
configuration missing-group example and finite endpoint-interval state.
The additional missing-covariate case closed an actual coverage gap.

Final DOM review read **69 tables, 295 body rows and 1,639 cells including
headers**. All **63 promised literal scientific/audit captions** are present.
Header/row/value sequences, keys, p-value formatting, missing values and nested
errors were reviewed, not inferred from kable counts. A final three-arm table
initially exposed source row names 7/8/9; row.names=FALSE removed them. The
affected article's three modes were rerun and its source/hash/status replaced,
then the package was rebuilt and checked again. Other sources did not change.
`targeted-reconciliation.txt` documents the exact replacement.

The structural audit reports zero errors/warnings and ten informational NULL
output-context items, reviewed as explained absence/status output. All nine
extracted PNG figures were viewed in a contact sheet; first-model diagnostics,
workflow arm means and two-arm pain means were also viewed at native size.
Axes, units, pointwise intervals, ordinal spacing and monochrome group
distinctions were inspected. No image blankness or clipping was found at those
image sizes. The final table reviewer independently reconciled captions and
complete cells after the correction.

The in-app browser refused file URLs under its security policy. **Zero browser
pages/viewports and zero interactive link navigations were observed.** No
localhost, alternate browser or indirect workaround was attempted. Consequently
page wrapping, scrolling and responsive table layout remain visually unverified.
Actual flat and installed article target/anchor checks are filesystem/DOM checks.
R's generated installed index uses dynamic-help /library routes; those were
statically mapped to installed targets using the local R navigation contract,
without starting a help server or claiming a standalone file-index navigation.
Repository source navigation is supplied by README's Rmd links. HTML article
links in Rmd sources target generated/installed articles rather than raw source
checkout files. This is not a pkgdown website review.

These are AI cross-reviews and execution audits, not human usability studies,
success percentages, external-user trials or peer review. Review records are
retained in `evidence/reviews/` and defects in `defect-ledger.csv`.

## Environment, final validation and source reconciliation

Windows 10 x64 build 19045; R 4.4.3; MIRA 0.0.1. R 4.5.3 was inspected but
lacks eight current Imports in its own dependency library. R 4.4.3 supplies
all dependencies, including geepack 1.3.13, mmrm 0.3.17 and nparLD 2.3.3.
Full versions/libraries are in `evidence/environment-4.4.csv` and logs;
no dependencies were installed. Locale was C. Pandoc/Quarto are existing
RStudio tools; the existing TinyTeX path was used for reference-manual checks.

Initial source installation used `work/library`; final bundle installation
used the separate `work/bundle-library`. Each render, extracted-code run and
reader-visible run uses its own Rscript --vanilla process and staged input
outside the repository. The runner rejects inside-repository output before
mkdir; its negative check left no directory. It hashes source/code/tests,
records paths/versions/optional flags/times, captures compact ledgers and
checks actual relative HTML targets/anchors. No source/install tree is used
as an analysis output directory.

| Final check | Result and evidence |
|:--|:--|
| Eleven article HTML renders | 11/11, 140.39 seconds summed, advanced optional branch true. |
| Extracted R paths | 11/11 in distinct clean processes, 116.36 seconds. |
| Reader-visible code paths | 11/11 in distinct clean processes, 118.76 seconds; no hidden analytical setup needed. |
| Existing full suite | 11,912 passed, zero failure/error/test warning, six explicit opt-in skips; 423.51 seconds, 387 test cases. |
| Final R CMD build | Success, 114.03 seconds, outside source tree. |
| Source bundle installation/discovery | Success; eleven expected vignettes and twelve HTML files including R-generated index. |
| Source/bundle input identity | 98 files: 97 byte-identical, DESCRIPTION semantically equivalent across original DCF fields after standard normalization. |
| Installed Rmd identity | All eleven byte-identical to final repository sources. |
| Final R CMD check | --as-cran --no-tests, 265.01 seconds; zero errors, zero warnings, one current-time-verification NOTE. |
| PDF/HTML reference manuals and vignette rebuild | All passed; --no-manual was not used. |
| Real external Quarto | Failed per format at processx launch; kept separate from ordinary HTML/source-only success. |

The suite has 11,918 total expectations including skips; six are skipped
expectations, not passes. One is opt-in CmdStan sampling and five are real
Quarto integration profiles. The existing asserts were unchanged. Tests run
separately with installed MIRA against source tests; the final bundle comparison
establishes identical R and test inputs, so --no-tests in the final check does
not hide a different code revision. Source/Rd/examples and vignette rebuilds
were checked normally, including --run-donttest examples through --as-cran.
Remote CRAN incoming checks were explicitly disabled; force-suggests remained
true. Repository-index access messages and Pandoc's deprecated-highlight-style
console message are distinguished from test warnings and the check's formal
zero-warning result. Timings are observed local runs with some concurrent
validation, not performance benchmarks. The earlier candidate build/check is
an intermediate work record and is not the final source-package proof.

The delivered archive is 829,207 bytes;
its SHA256 is `750c8f0311e74a2ac0dfe030dfc469a6a4e199bc70e8e51750b983ce38852405`. Exact final sources,
per-mode hashes, bundle inputs, installed Rmds, table-review hashes and artifact
hashes form the reconciliation chain. The final display correction is fully
included. DESCRIPTION normalization is explicitly not called byte identity.
`evidence/delivery-link-checks.csv` checks repository README source links,
flat/installed article targets/anchors and statically mapped R-help index
routes. `_pkgdown.yml` is parsed
statically; no site initialization or construction was run.

Executed drivers, commands, environment flags and minimal reproductions are
under `reproducibility/`; tests/build/check/render logs and compact ledgers
are in `evidence/` and `handbook/validation/`. The final source archive includes
standard installed vignette sources and HTML. No fitted-model dump archive or
redundant release ZIP is delivered.

## JOSS/CRAN scope and remaining human responsibilities

`SOURCE_RECORD.md` records official URLs, access date and mirror fallbacks.
Current [JOSS criteria](https://joss.readthedocs.io/en/latest/review_criteria.html)
and [checklist](https://joss.readthedocs.io/en/latest/review_checklist.html)
support assessing installation/dependencies,
worked examples, documented functionality and reproducible verification.
This work supplies evidence for those dimensions. JOSS's broader requirements
also concern public development history, research use/impact, human design
contribution, licensing/community practice and the manuscript. This task does
not attest their adequacy or make an editorial recommendation.

The local DESCRIPTION still strongly imports CmdStanR and points to Stan
R-universe. The [current CRAN policy](https://stat.ethz.ch/CRAN/web/packages/policies.html)
expects strong dependencies from CRAN or
Bioconductor; Additional_repositories is not an automatic exception. A further
identity concern was confirmed: [Bioconductor already distributes an unrelated
MIRA package](https://bioconductor.org/packages/release/bioc/html/MIRA.html).
Its naming policy implication requires maintainer review.
Neither dependency strategy nor package name was changed here. Package size
and local check success do not settle eligibility, resource suitability on a
shared farm, current-R compatibility or cross-platform behavior.

Residual categories are explicit: three software/engine limitations above;
Quarto and browser validation environment limits; dependency/name submission
issues; study-specific statistical assumptions; and pending human review.
None invalidates the demonstrated basic summaries or the capstone's taught
primary contrast. No known uncorrected fact/API/sign/sample error remains in
those taught conclusions after these scoped checks; this is not a guarantee
of universal correctness or comprehension.

OpenAI Codex with a GPT-6-family agent and three parallel AI agents assisted
with source/test reading, synthetic cases, editing, code execution, primary
source research, audits and reporting. The exact deployment identifier is not
exposed and is not invented. PowerShell/R, Python structural/image tools,
hash/Git inspection and web retrieval were used. No human verification is
claimed. AI review records and a new dated AI_USAGE entry preserve scope
without reattributing historical work.

Under the [JOSS AI usage policy](https://joss.readthedocs.io/en/latest/policies.html#ai-usage-policy),
human authors must review/edit and validate accepted AI-assisted material,
confirm study-specific estimands/assumptions and consequential calculations,
assess independent usability and browser layout, resolve package identity/
dependency eligibility, and prepare accurate manuscript/AI disclosures.
Those are later human responsibilities, not tasks an agent can certify as
completed. The documentation candidate is ready for that review with the
explicit local validation limits stated above.
