# MIRA frequentist longitudinal handbook: substantive revision report

Revision date: 2026-10-29 (Europe/Rome). Package: MIRA 0.0.1.

This report distinguishes historical evidence, intermediate investigations, and
execution of the final revised sources. The final build/check and scope proof
is recorded separately in section E.

The assignment was a scientific, didactic and technical revision of the
existing frequentist longitudinal documentation. It was not an engine rewrite.
The work compares the tutorials with the public analysis/reporting functions,
improves concrete deficiencies, and adds independent checks where the
scientific meaning of a number matters. Execution success does not establish
statistical validity, formal CRAN compliance, or JOSS acceptance.

## A. Changes and preservation

The repository at `C:/Users/david/OneDrive/Desktop/MIRA/MIRA` already contained
modified `.Rbuildignore`, `.gitignore`, `DESCRIPTION`, `NEWS.md`, `README.md`,
and `dev/AI_USAGE.md`, together with an untracked vignette collection,
`_pkgdown.yml`, the historical implementation report and the validation runner.
Git status and diffs were examined before editing. The primary agent saved a pre-revision
snapshot and hashes outside the repository. These pre-existing changes were
treated as user work, not as a clean baseline to replace.

| File revised in this phase | Reason and substantive change |
|:--|:--|
| `vignettes/MIRA.Rmd` | Distinguish standard knitr/rmarkdown/Pandoc vignette building from Quarto report rendering, and describe multiple-outcome workflows without implying joint modeling. |
| `vignettes/frequentist-long-getting-started.Rmd` | Make participant versus measurement counts explicit; independently check visit counts, paired endpoint count and Student-t change interval; qualify overview as contextual rather than a universal analysis sample. |
| `vignettes/frequentist-long-data-and-detection.Rmd` | Actually analyze custom-regex data; verify unordered manual input with named labels; demonstrate absent-ID fallback; require the intended ambiguity error; verify automatic covariates and qualify constant/all-missing classification. |
| `vignettes/frequentist-long-descriptive-analysis.Rmd` | Correct the meaning/denominator of person-mean residual SD; independently reconstruct availability, paired change/interval and complete-profile ANOVA ICC; improve figures, units and p-value display. |
| `vignettes/frequentist-long-treatment-groups.Rmd` | Independently verify named arm contrasts, group samples and Welch intervals; retain pair identities in three-arm examples; define actual multiplicity families and dependence conditions; improve small-p-value display and group figures. |
| `vignettes/frequentist-long-missing-data.Rmd` | Verify paired samples by participant identity and arithmetic, add sample ledger/interval check, distinguish execution/estimability statuses, and separate the reproduced sparse plotting defect from correct nonplotting use. |
| `vignettes/frequentist-long-multiple-outcomes.Rmd` | Verify outcome-specific endpoint samples, changes and response labels; make extraction safe when a successful outcome has no final pairs; execute a sparse final-visit example with plotting disabled. |
| `vignettes/frequentist-long-advanced-models.Rmd` | Audit exact model subject-by-visit samples and formulas; focus contrast/status tables; explain afex's distinct mean-model expansion; center its numeric covariates after response generation; improve inference and diagnostic presentation. |
| `vignettes/frequentist-long-reporting.Rmd` | Independently check source stem, CSV values, incompatible-input errors and shallow archive behavior; retain the complete reporting contract and source-only routine builds. |
| `vignettes/frequentist-long-complete-workflow.Rmd` | Align the prespecified primary estimand with the extracted adjusted difference in changes, reconstruct its direction from four marginal means, audit endpoint/model participant sets and preserve the discussed analysis and contrasts in the report payload. |
| `dev/validate_frequentist_vignettes.R` | Guard output paths before directory creation; allow future Bayesian articles and an automatic installed index; record source hashes, selected installed library, status, sample and contrast evidence. |
| `dev/FREQUENTIST_VIGNETTES_REVISION_2026-10-29.md`, `dev/AI_USAGE.md`, `NEWS.md` | Record this phase separately, with actual final checks, limitations, AI assistance and pending human responsibilities. |

`dev/FREQUENTIST_VIGNETTES.md` is the historical implementation report and
remains unchanged by this revision. Its earlier timings, 11,912 test
expectations, six skips, build/install results and check notes are historical
claims to compare with current evidence; none is inherited as validation of
the new version. Its earlier pkgdown attempt is also historical. This phase
does not build or publish the website.

Material already correct was retained: the connected ten-article route;
deterministic local generators; the installed public API; explicit/inferred
outcome-name provenance; metric-name-safe extraction; later-minus-earlier
changes; pointwise intervals versus adjusted p-values; baseline as a repeated
response; ordinal versus physical time; and non-interchangeable missingness
assumptions and inferential targets.

No statistical engine source, public signature, Bayesian function, Stan model,
prior, posterior calculation, Bayesian example/documentation, reference Rd,
namespace or existing test was intentionally changed. MIRA was reinstalled into an isolated library; no third-party dependency was
installed and no dependency strategy was altered. No reset/restore, commit,
push, PR, tag, release, deployment or publication was performed.

The final phase-only scope/hash comparison and whitespace result are recorded
with the final package proof in section E. The original historical report and
pre-existing user changes remain distinct from this phase's revisions.

## B. Article-by-article teaching and execution coverage

The runtime field separates focused intermediate checks from the final
integrated run. Focused times are evidence for the reviewed copies at that
moment, not final package-build timings. The final integrated status and hash
manifest take precedence after all edits have been applied.

### MIRA: a guide to longitudinal analysis

| Dimension | Revised coverage |
|:--|:--|
| Audience | Researchers, clinical scientists and statistical programmers choosing an entry point. |
| Question | What scientific/data problems does the frequentist interface address, and which article should the reader use next? |
| Dataset | No analytical dataset; installed package version is queried. |
| Methods | Learning sequence and installation/rendering distinctions. |
| API | `mira_info_long()`, `mira_report_freq_long()` introduced; precise help pages and vignette discovery identified. |
| Output | Package version and links only to existing articles. |
| Interpretation | Repeated measurements remain correlated; execution does not justify scientific assumptions. CmdStanR, CmdStan, Quarto CLI, Pandoc and TeX have distinct roles. |
| Runtime | Focused revised render passed, 1.96 s. Final revised-source render: PASS, 1.66 s. |
| Residual limits | No fresh provisioning or network installation demonstration; Bayesian tutorials remain a future separate collection. |

### Getting started

| Dimension | Revised coverage |
|:--|:--|
| Audience | New users learning one complete first analysis. |
| Question | How do observed walking speeds and within-person changes vary over follow-up? |
| Dataset | Twelve synthetic participants, baseline/week 4/week 12, speed in m/s. |
| Methods | Core visit summaries, paired change and Student-t interval; no complex model on this small teaching sample. |
| API | Automatic data/ID/outcome detection; `analyses="none"`, `verbose=FALSE`; named components and RDS round trip. |
| Output | Class, overview/map, visit N/mean/SD/CI, named baseline-final change, trajectories and long rows. |
| Interpretation | Thirty-six measurements represent twelve independent participant profiles. Overview is not every estimator's sample. Paired mean/interval/count are checked independently; intervals are pointwise and assume independent participant differences. `none` retains core inference and Friedman. |
| Runtime | Focused revised render passed, 1.79 s. Final revised-source render: PASS, 1.53 s. |
| Residual limits | Small synthetic example; no clinical effect, calibrated individual treatment effect or model-selection claim. |

### Data preparation and detection

| Dimension | Revised coverage |
|:--|:--|
| Audience | Analysts importing wide data and resolving its scientific mapping. |
| Question | Which columns become participants, outcome visits, groups and covariates? |
| Dataset | Thirty symptom profiles; 36 rehabilitation profiles with pain/function and different naming conventions. |
| Methods | Automatic/manual detection, alternate suffixes, custom regex and function parser, validation and provenance. |
| API | `id`, `time_vars`, `time_labels`, `outcomes`, `variable_pattern`, `inspect_only`, `strict_id`, `covariates`, `analyses`, `verbose`; arm/reference discussed. |
| Output | Parsed maps, resolved config, names/types/provenance, actual parser analyses and controlled ID/type/time ambiguity errors. |
| Interpretation | Numeric suffixes establish order rather than a calendar unit. Named labels survive reordering. Generated row IDs do not establish real participant identity. Metric attributes distinguish declared from inferred names. Automatic covariates satisfy structural rules, not scientific adjustment criteria. |
| Runtime | Focused revised render passed, 2.44 s, including regex/function analyses and intended errors. Final revised-source render: PASS, 2.00 s. |
| Residual limits | Data dictionary, chronology, identity and adjustment justification remain investigator responsibilities; not every parser/value combination is tested. |

### Descriptive analysis

| Dimension | Revised coverage |
|:--|:--|
| Audience | Applied researchers who understand the wide-data introduction. |
| Question | How do visit distributions, paired changes, association, variation and reliability differ? |
| Dataset | Sixty-four synthetic observational SBP participants, four ordered visits with incomplete profiles; mmHg. |
| Methods | Means/intervals, paired change, trajectories, Pearson/Spearman correlation, person-mean residual SD, complete-profile ANOVA ICC, IQR flags, core Friedman branch. |
| API | Explicit ID/visits/labels/outcome; logical optional switches; `alpha`, direction/threshold, no covariate adjustment, mapping-based extraction and plot copies. |
| Output | Visit samples 64/62/61/59; baseline-final pairs 59; complete profiles 56; focused tables and two scientific figures. |
| Interpretation | Correlation does not measure agreement or exclude systematic change. Residual SD uses finite residuals about each participant's available mean with denominator N-1, not an ANOVA pooled-within estimator. ICC is independently reconstructed from a one-way ANOVA, with complete-profile assumptions/limits. |
| Runtime | Focused revised render passed, approximately 5.42 s; reconstructed ICC 0.8195281. Final revised-source render: PASS, 4.00 s. |
| Residual limits | No fitted-model ICC demonstration in this article; summaries do not identify a missingness mechanism. Threshold is illustrative. |

### Treatment groups

| Dimension | Revised coverage |
|:--|:--|
| Audience | Researchers interpreting comparative repeated-measure studies. |
| Question | Do observed visit values or baseline-follow-up changes differ between named groups? |
| Dataset | Ninety-six randomized synthetic pain participants, two groups of 48 and three visits; separate 108-person function study with three groups of 36. |
| Methods | Group descriptives, Welch/rank visit and change tests, omnibus and pairwise comparisons, Hedges' g, response counts and plots. |
| API | Explicit arm/reference, direct optional flags and `analyses` precedence, outcome/direction/threshold, `alpha`, `p_adjust_method`, no covariates. |
| Output | Named B-minus-A contrasts and independently checked raw means, sample sizes and Welch limits; all non-reference pairs retained for three arms. |
| Interpretation | Different within-group significance does not test a between-group difference. Baseline p-values do not certify randomization or select covariates. Holm's full family includes nine visit-by-pair tests, even when only three final pairs are displayed. Intervals remain pointwise. |
| Runtime | Focused revised render passed, approximately 3.39 s. Final revised-source render: PASS, 3.02 s. |
| Residual limits | Complete teaching data do not demonstrate differential dropout; adjusted/global model inference is deferred. No validated clinical threshold is asserted. |

### Missing data

| Dimension | Revised coverage |
|:--|:--|
| Audience | Researchers interpreting incomplete and analytically unavailable responses. |
| Question | Which participants contribute to each summary, pair, complete-profile procedure or model? |
| Dataset | Eighty-two synthetic respiratory participants, unequal groups, visits at weeks 0/2/9/18, missing baselines, gaps/dropout and one infinite value; separate sparse instrument. |
| Methods | Raw versus analytical availability, pairwise change/interval, correlation pair counts, complete profiles, group imbalance and module-status investigation. |
| API | Explicit wide mapping, arm/reference, `analyses`, `alpha`, Holm, config adaptation and diagnostics. |
| Output | Independently checked changes for every pair, exact endpoint participant identities, endpoint Student-t interval and enrollment/pair/complete-profile ledger. |
| Interpretation | Available-case and complete-case analyses differ. A baseline-conditioned dropout generator is not automatically jointly MAR once that baseline predictor is masked. Pair deletion does not remove selection bias. Status distinguishes unrequested, inapplicable, disabled, missing dependency, failed computation and partially usable output. |
| Runtime | Focused revised render passed, 6.38 s. Sparse plotting defect reproduced separately. Final revised-source render: PASS, 4.17 s; adaptation and skipped-model states remain explicit. |
| Residual limits | No imputation or MNAR identification; ordinary GEE is not a general MAR remedy. Sparse plotting defect remains an engine issue outside this scope. |

### Multiple outcomes

| Dimension | Revised coverage |
|:--|:--|
| Audience | Researchers managing several endpoints with different scales, schedules and favorable directions. |
| Question | How can outcome-specific results be configured, safely extracted and interpreted? |
| Dataset | Seventy-two synthetic participants: symptom (three visits), function (three), positive biomarker (two), outcome-specific endpoint availability; separate data copy with no final function visit. |
| Methods | Parallel outcome analyses, endpoint extraction, response classifications, threshold sensitivity and multi-outcome plotting. |
| API | Automatic/explicit outcome families, named visit/label lists, named directions/thresholds, global analysis settings, single/multi/error result classes. |
| Output | Endpoint counts/mean changes checked from original data, response labels checked independently, unavailable endpoints retained with n=0/NA/status. One pair is explicitly identified as insufficient for an interval. |
| Interpretation | The container is not a joint multivariate model. Scale/favorable direction/threshold require scientific justification. Classifications do not imply individual causal effects or account automatically for measurement error. Across-outcome multiplicity is not jointly corrected. |
| Runtime | Focused revised render passed, 4.21 s, including the unavailable-final-visit case. Final revised-source render: PASS, 3.69 s. |
| Residual limits | The defensive `mira_info_error_long` branch is explained but not deliberately induced by an ordinary failing-outcome tutorial; that branch must not be claimed as newly executed. |

### Advanced models

| Dimension | Revised coverage |
|:--|:--|
| Audience | Readers trained in regression and repeated-measures inference. |
| Question | Which longitudinal estimand, mean/covariance structure, sample and hypothesis does each implemented procedure represent? |
| Dataset | Ninety-six independent synthetic participants, two arms, four complete equally spaced visits; function score, age, sex and separate screening severity; random intercept/slope variation. Numeric covariates are centered after generating unchanged responses. |
| Methods | Primary/random-slope mixed models, emmeans, RM-ANOVA/GG/HF, Friedman, KR, CR2, three GEE fits, two NLME fits, MMRM, rank ATS/WTS, effect sizes/R-squared and sensitivity. |
| API | Model/arm analysis bundle, explicit covariates/reference/labels, pointwise confidence level, primary p-adjustment, statuses/diagnostics, opt-in no-adjustment refit. |
| Output | Thirteen module statuses; exact subject-by-visit sample ledger; formula ledger; named focused contrasts and coefficient-set hypotheses; residual panels and engine-specific criteria. |
| Interpretation | Methods are not interchangeable. TIME/ARM joint tests include appropriate interaction coefficients. Afex expands categorical between-subject factors factorially and time interactions with them/numeric covariates, unlike additive-covariate mixed/GEE/NLME/MMRM mean structures. Centering changes covariate reference, not generated outcomes. Ordinal slope/AR(1) is not physical calendar time. |
| Runtime | Final revised-source render passed in 24.22 s with the opt-in refit enabled. All 13 module containers were available: nine complete at the inspected level and four partially usable (random slope and three GEE native subresults). Nested errors are disclosed in D/F; no container is falsely treated as wholly successful. |
| Residual limits | Complete-data agreement of sample keys does not equate hypotheses; random-slope and correlation adequacy still require study-specific assessment. No new missing-data or arbitrary-formula engine interface is implied. |

### Reporting

| Dimension | Revised coverage |
|:--|:--|
| Audience | Analysts preparing auditable bundles and researchers tracing a reported number. |
| Question | What is preserved, exported, displayed, rendered and opened, and how can each be verified? |
| Dataset | Forty-eight synthetic rehabilitation symptom profiles; a second function outcome and auxiliary regression illustrate multi-outcome/supplementary objects. |
| Methods | Stored-result versus direct-data reporting; QMD/runtime/RDS/CSV/manifest inspection; source-only normal build, shallow traversal, incompatible inputs, regeneration and opt-in external rendering. |
| API | All 21 reporting parameters are illustrated; analysis arguments forwarded through `...` are distinguished from forbidden arguments with an existing result. |
| Output | Stem check, source existence/status, manifest-linked CSV values, exact simple-result payload, display/archive independence and depth behavior. |
| Interpretation | `expected_files` is a target list, not render success. Sources are not finished documents; old rendered files can survive a failed attempt. Display limits differ from `max_depth`. Compatibility-only `include_complete_output` does not alter visible content. `extra_objects` are preserved but excluded from CSV. |
| Runtime | Independent focused source-only render passed, 6.63 s. Final revised-source render: PASS, 3.46 s. External CLI rendering: current explicit HTML/DOCX/PDF attempts all failed at processx pipe creation with Windows Access denied (system error 5); sources were retained. |
| Residual limits | Serialized result may contain participant data/IDs; bundle governance is investigator responsibility. Native model identity across serialization requires scientific-component checks, not pointer/environment identity. |

### Complete research workflow

| Dimension | Revised coverage |
|:--|:--|
| Audience | Experienced analysts connecting one scientific question to a preserved analytical result. |
| Question | What is the adjusted Program-minus-Usual-care difference in mean baseline-to-week-12 fatigue changes? |
| Dataset | One 120-person synthetic randomized dataset used throughout, four visits, complete screening burden, burden-dependent follow-up loss and independent baseline masking. |
| Methods | Inspection, core/group summaries, primary repeated-response mixed model, marginal interaction contrast, diagnostics, model perspectives and an additional complete-profile descriptive comparison. |
| API | Explicit ID/visits/labels/outcome/arm/reference/covariates, all optional modules, direction/threshold/alpha/Holm, stored-result source-only report. |
| Output | Named unadjusted endpoint contrast; adjusted arm-by-time-pair interaction reconstructed from four marginal means; endpoint/model participant sets; exact archived scientific tables/settings/native fixed effects/covariance and both reported contrasts. |
| Interpretation | The earlier adjusted simple-arm endpoint contrast answered a different estimand. It is now replaced by the prespecified difference in changes. The subtraction is Program change minus Usual-care change; negative favors Program. Same-scale unadjusted and adjusted change contrasts can still use different participants and assumptions. Fully observed screening burden makes this generator's MAR statement defensible; MCAR rank inference is not thereby an equally justified primary analysis. |
| Runtime | Initial independent run found an overly strict whole-object RDS identity assertion after all DID checks passed. The primary agent replaced that assertion with meaningful component/native-parameter preservation checks; the revised extracted examples then ran without error. Final integrated render: PASS, 32.19 s; adjusted Program-minus-care change difference -5.0579 points, pointwise 95% CI [-7.0078, -3.1080], within-family Holm p=3.4459e-06. The model uses 443 observed long rows from 120 participants; the unadjusted endpoint comparison uses 49 Usual-care and 53 Program participants. |
| Residual limits | No power calculation, clinical-threshold validation, imputation, MNAR sensitivity or automatic causal identification. The generated report freezes the discussion's analysis; it does not validate its scientific assumptions. |

## C. Complete public API coverage

The companion matrix `api-coverage-ast.csv` has one row for each of 43 public
parameters (22 analysis, 21 reporting). `api-ast-calls.csv` retains call text,
article/chunk/line, branches and controlled-error context;
`api-ast-call-arguments.csv` retains explicit expressions and reporting
forwarded-argument names. `api-ast-input-md5.csv` identifies the audited source
copies. The inventory script parses executable chunks and resolves positional
and named arguments with the installed public formals; it does not evaluate
examples or equate textual presence with verification.

At the recorded inventory snapshot, all 43 parameters appeared explicitly in
at least one routine example. There were 50 public calls: 48 routine and two
opt-in calls, with 393 explicit parameter occurrences. The final-source AST refresh and annotated coverage matrix are supplied with the companion API evidence; their hashes identify the exact source snapshot. Display-only edits do not add new public calls or change argument coverage.

The two opt-in routes are the advanced no-adjustment sensitivity refit
(`MIRA_RUN_ADVANCED_VIGNETTES=true`) and actual external rendering
(`MIRA_RENDER_VIGNETTE_FORMATS`). Conditional inference/plot access also depends
on actual module availability; the final ledgers, not source inspection,
establish which outputs were obtained. Calls inside controlled error handlers
demonstrate expected rejection, not a successfully fitted analysis.

| Coverage dimension | Evidence and limit |
|:--|:--|
| Explained in prose | Manual article review; an inline-code parameter token alone is not adequate evidence of explanation. |
| Present in executable example | AST call/argument inventory, including positional and forwarded inputs. |
| Actually executed | All ten final Rmd executions, installed-path assertion and final source-hash record; optional and controlled-error routes are distinguished. |
| Only optional / not executed | Branch metadata and final environment flags; successful source-only reporting does not exercise Quarto rendering or document opening. |
| Independently validated | Raw-data counts/identities/means/intervals, ANOVA ICC, p-adjustment arithmetic, four-cell DID, CSV/payload/native-parameter checks. These verify specific behaviors, not all values of an argument. |

The complete parameter inventory is:

| Function | Public parameters |
|:--|:--|
| `mira_info_long()` | `data`, `id`, `time_vars`, `time_labels`, `arm`, `reference_arm`, `arm_tests`, `alpha`, `plots`, `model`, `outliers`, `correlations`, `verbose`, `p_adjust_method`, `improvement_direction`, `stable_threshold`, `strict_id`, `outcomes`, `covariates`, `analyses`, `variable_pattern`, `inspect_only` |
| `mira_report_freq_long()` | `x`, `data`, `...`, `output_dir`, `output_file`, `format`, `title`, `subtitle`, `author`, `date`, `sections`, `include_complete_output`, `extra_objects`, `max_table_rows`, `max_table_columns`, `digits`, `max_depth`, `render`, `open`, `overwrite`, `quiet` |

Important behaviors checked against current implementation include:

- `analyses` replaces the five optional logical switches; `none` retains
  core inference. There is no independent public switch selecting only the
  primary mixed fit from its implemented model bundle.
- Explicit names declare metric prefixes, while inferred outcomes retain base
  metric names. Structural fields/native engine fields follow their own
  contracts. Table attributes resolve MIRA-owned metric names.
- Per-outcome visit/label/direction/threshold settings are distinguished from
  global model/alpha/adjustment/module selection. The multi-result container
  performs separate analyses.
- `alpha` affects supported confidence intervals, not clinical importance or
  whether results appear. Adjusted p-values do not turn intervals simultaneous.
- `strict_id=FALSE` changes admissibility and interpretation; it is not a
  duplicate-removal procedure. Row-ID generation does not repair repeated real
  participants stored in multiple wide rows.
- Reporting's default requests HTML and PDF; `render` defaults TRUE and `open`
  to `interactive()`. Routine examples explicitly use FALSE for both.
- Section/display settings leave eligible CSV tables intact; traversal depth
  can omit tables. `include_complete_output` is retained for compatibility
  without the claimed visible-content effect of its name.

No claim is made that every admissible argument value, partial-match case,
combination, error branch, external rendering path or native-engine option is
exhaustively tested. Imported engines remain installation dependencies even
when a particular module checks availability at runtime.

## D. Statistical coverage, actual samples and hypotheses

| Procedure | Quantity/hypothesis and sample | New verification / interpretive boundary |
|:--|:--|:--|
| Visit descriptives and mean intervals | Available finite participants separately by visit | Counts reconstructed from raw matrices; between-person spread differs from mean uncertainty. |
| Paired change and t interval | Later-minus-earlier differences among participants observed at both selected visits | Raw differences, participant identities, pair counts and elementary t limits reconstructed independently. |
| Paired t/Wilcoxon and Friedman | Mean/rank questions are distinct; Friedman requires complete profiles, while post hoc pairs can use different available pairs | Core inference remains when optional modules are disabled; absence of an estimate is not evidence of the null. |
| Correlations and variability | Pair-specific association; sample SD about available participant means; complete-profile one-way ANOVA ICC | Pair counts and residual SD checked independently; ICC reconstructed from a base-R one-way ANOVA rather than copied from implementation. Correlation is not agreement. |
| Arm visit/change comparisons | Observed B-minus-A means versus B-minus-A mean paired changes; group and pair availability differ | Named final visits/groups, raw samples, mean differences and Welch limits checked. Within-arm significance is not a between-arm test. |
| Primary mixed model | Categorical visit * arm, additive selected covariates, participant random intercept; REML estimation and ML nested interaction comparison | Native formula/frame and exact source subject-by-visit identities. Baseline is a response; this does not automatically impose baseline equality or endpoint ANCOVA. |
| Random-slope model | Same fixed mean with participant slope on ordered visit index | Exact model-frame keys and formula; ordinal slope is not automatically a per-month/week coefficient. |
| emmeans | Marginal visit/arm means, simple effects and contrasts of contrasts | Workflow primary difference in changes selected by both labels and rebuilt from four cell means; applicable primary p-adjustment independently checked. |
| Repeated-measures ANOVA | afex type-III repeated-measures effects, GG/HF correction; complete profiles under engine rules | Afex's categorical factorial and time-by-covariate expansion differs from the other additive-covariate models. Centering makes zero the study-mean reference. Sphericity corrections do not repair selection bias. |
| Kenward-Roger | Approximate model-based joint Gaussian F hypotheses | Explicit coefficient sets: TIME/ARM may include interaction coefficients; not a single lower-order coefficient and not a sandwich correction. |
| CR2 | Participant-cluster robust covariance; Satterthwaite coefficient and HTZ joint tests | Primary-model sample and coefficient constraints; finite-cluster approximations still require adequate independent clusters/mean specification. |
| GEE | Gaussian marginal mean, independence/exchangeable/AR(1) working correlation | Exact three fit samples/clusters/waves; ordinary GEE is not a general MAR or MNAR method. QIC/CIC are not likelihood AIC/BIC. |
| NLME | Additive fixed mean, participant random intercept, CS or ordinal AR(1) residual structure | Actual data/formula sample audit; fixed formula alone does not describe random/covariance assumptions. |
| MMRM | Gaussian marginal mean with common-across-arm unstructured response covariance, REML/Satterthwaite | Exact source sample; no subject random effect or imputation. Available-response inference needs ignorable MAR and adequate mean/covariance model. Native information criteria are not automatically comparable across engines. |
| Rank longitudinal inference | Unweighted relative marginal effects and ATS/WTS for factorial cells | Exact finite-response keys; MCAR requirement for implemented incomplete-data inference. These are relative-effect hypotheses, not Gaussian mean coefficients or covariate-adjusted contrasts. |
| Effect sizes and R-squared | Paired dz, Hedges' g, Kendall W, eta-squared, model ICC/Nakagawa R-squared | Different denominators/estimands remain explicit; availability and diagnostics are retained. |

The advanced teaching dataset has 96 participants and 384 finite response
rows. Its agreement of fitted/data sample keys is deliberately checked, rather
than inferred from `nrow(study)`. KR, CR2 and emmeans inherit the primary model
frame. The complete-profile tests count participants rather than long rows.
The incomplete workflow separately compares actual endpoint/model participant
sets and uses observed predictors in the model eligibility rule.

**Final statistical evidence.** The final run records `module_status.csv`,
`sample_ledger.csv` and `formula_ledger.csv` under the advanced article's
validation directory; all nine model/data sample entries have 384 observations
and 96 participants, checked by exact subject-by-visit keys. The 13 module
containers comprise nine without recorded errors/warnings and four partial
containers: random slope and GEE independence/exchangeable/AR(1). Their
unavailable supplementary subresults are documented below; the fitted models
and independently constructed joint tests remain available. No top-level
fitting warning should be substituted for the named nested error record.

The complete workflow's final `adjusted_endpoint.csv` records the primary
adjusted difference in changes as -5.05789379058301 fatigue points, pointwise
95% interval [-7.00781296696341, -3.10797461420262], raw p=5.74315966691906e-07
and Holm p=3.44589580015144e-06 across six arm-by-visit-pair contrasts.
The raw group comparison is -5.00702734555245 points, interval
[-7.13168761565061, -2.88236707545428], Holm p=2.86400938033758e-05,
using 49 Usual-care and 53 Program endpoint-complete participants. The model
instead uses 443 finite rows and all 120 participants; 18 model participants
lack a complete endpoint pair. These are computed records, not manually
inserted values in the tutorial's conclusion.

**Four partially usable containers.** In the random-slope module,
`errors[["random_effects_comparison"]]` records
`all models must be fit to the same data object`; `random_effects_lrt` is
unavailable even though the fit is converged and nonsingular. The primary
and slope model calls use different data-expression names; a native comparison
guard rejects them despite identical observed keys. For each GEE working
correlation, `errors[["anova"]]` records
`object 'geese.fit' of mode 'function' was not found`. The namespace-only
native sequential ANOVA lookup fails; MIRA's separate three-row joint
`effect_tests` tables remain available. This is neither an absent engine nor
failure of the GEE fit itself. `performed=TRUE` does not certify all subresults.

Two independent minimal native reproductions and diagnostic workarounds were
executed in this phase (`partial_module_minimal_repro.R` and
`probe_partial_modules.R`; 5.0 s and 5.6 s). Refitting the original mixed-model
formulas on the same explicitly named data object allowed a separate native
comparison; attaching geepack in the temporary R process allowed its native
ANOVA on the same stored fits. Those workarounds do not modify the MIRA result,
resolve the integration defects, or make sequential ANOVA equivalent to
MIRA's joint TIME/ARM hypotheses. The variance-component boundary caveat
remains. `partial_module_findings.md` records mechanisms, sample proof and
interpretation limits; no engine patch was made.


Multiplicity is described by implemented family, not by an undifferentiated
claim of correction: core paired visit families, arm-value visit/pair families,
change families, separate omnibus/missingness families and advanced contrast
families remain distinct. The workflow's four-visit/two-arm interaction family
contains six arm-by-visit-pair contrasts; its selected primary contrast does
not inherit a simultaneous interval from the Holm p-value. Across outcomes,
methods and sensitivity analyses, a study-wide strategy remains the researcher's
responsibility. Holm/Bonferroni target familywise error, BH false discovery rate
under independence/appropriate positive dependence, and BY under arbitrary
dependence given valid individual p-values.

## E. Validation: environment, commands and evidence

### Historical, intermediate and final evidence

The unchanged historical implementation report is not evidence that the
revised sources pass. Focused checks in this revision informed the edits and
are labeled in B. Final revised-source execution is recorded by
`work/final-delivery-vignettes/vignette-status.csv`, its source-MD5 inventory
and session information. All ten final articles passed in **79.94 seconds**
in total: MIRA 1.66; getting started 1.53; data/detection 2.00; descriptive
4.00; groups 3.02; missing data 4.17; multiple outcomes 3.69; advanced 24.22;
reporting 3.46; complete workflow 32.19 seconds. These are local elapsed
times, not a performance guarantee for CRAN or another computer.

R 4.4.3 and 4.5.3 were inspected before selecting the environment. R 4.4.3
supplied all required engines; R 4.5.3's inspected library lacked afex,
clubSandwich, CmdStanR, geepack, lmerTest, mmrm, nparLD, quarto and pdftools.
No third-party dependency was installed to make the examples pass.

| Environment component | Recorded version |
|:--|:--|
| OS / R / architecture | Windows 10 x64 build 19045; R 4.4.3 (2025-02-28 UCRT), x86_64-w64-mingw32 |
| Isolated MIRA | 0.0.1 at `C:/Users/david/Documents/Codex/2026-10-29/il-x20/work/library/MIRA` |
| Primary engines | lme4 2.0-1; lmerTest 3.2-1; emmeans 2.0.3; afex 1.5-1; pbkrtest 0.5.4 |
| Other inference engines | clubSandwich 0.7.0; geepack 1.3.13; nlme 3.1-168; mmrm 0.3.17; nparLD 2.3.3; performance 0.15.3 |
| Data/graphics/rendering packages | ggplot2 3.5.2; knitr 1.50; rmarkdown 2.29; quarto R package 1.5.1; TMB 1.9.21 |
| Existing external tools | Quarto CLI 1.9.38; Pandoc 3.8.3; pdfTeX 1.40.29, TeX Live 2026 |
| Locale / time zone | C for validation processes; system code page 65001; Europe/Rome |

The environment inventories (`environment-4.4.3/`, `environment-4.5.3/`) and
final `session-info.txt` retain full package/library details. Quarto was
located under RStudio's existing resources; no tool was installed or security
check disabled. Standard R Markdown vignettes use knitr/rmarkdown/Pandoc and
do not require successful MIRA Quarto-report rendering.

### Completed execution and interpretation

| Validation | Actual result | Evidence |
|:--|:--|:--|
| Isolated source installation and path assertion | Passed; installed public namespace used | `logs/revised-install.log`; runner's installed-path assertion and `validation-context.txt` |
| Final ten frequentist Rmd articles | 10/10 PASS, 79.94 s total | `final-delivery-vignettes/vignette-status.csv`; per-file source hashes and session |
| Optional advanced no-adjustment refit | Executed in final advanced article | `MIRA_RUN_ADVANCED_VIGNETTES=true` recorded in context; article PASS |
| Model bundle | 13 available/performed containers: nine complete and four partial at the inspected level | Final module/sample/formula ledgers; details D/F |
| Source-only report examples | Passed manifest paths, CSV values, archive depth, error/overwrite contracts and scientific payload preservation | Reporting/workflow assertions; no external renderer used by routine examples |
| Independent reporting contracts | Compatibility flag TRUE/FALSE produced identical knitted content; positional data and default HTML/PDF normalization passed | `report-contract-probe.R`, `logs/report-contract-probe.log`, `report-contracts/compatibility-status.txt` |
| Existing test suite | 11,912 passing expectations; zero failures, errors or warnings; six explicit skips across 387 test cases; 393.38 s | `run-tests.R`, `logs/source-tests-installed.log`, `test-summary.csv`, `test-results.rds` |
| Explicit actual Quarto rendering | HTML, DOCX and PDF each attempted; all three failed before successful CLI startup with processx Access denied, system error 5; sources retained | `report-contracts/external-render-status.csv` and probe log |
| Source-tree output guard | Rejected forbidden output before directory creation | Independent guard test; no `vignettes/guard-must-not-be-created` created |
| Rendered output visibility | Corrected conditional chunks to explicitly print or separate significant outputs; static/HTML caption audit and real visual inspection performed | `audit_html.py`, structural JSON/CSV, figure contact sheet and browser observations; final confirmation with package proof below |

The first test-driver attempt used an unsuitable `test_dir` loading context
and was interrupted after missing installed-namespace symbols. That harness
setup was corrected to `load_package="installed"`; the results above are from
the complete corrected run, not from the interrupted attempt. No software
test failure is inferred from the abandoned harness.

The first bundle-identity assertion treated R CMD build's normal DESCRIPTION
canonicalization as a mismatch. The verifier was corrected to compare source
DESCRIPTION fields semantically after whitespace normalization, allow only
normal build-generated metadata, and retain byte identity checks for every
other audited file. This was a verification-assertion correction, not evidence
of a changed statistical source or a final bundle failure.

The six test skips are explicit gated integrations: one CmdStan smoke test
requires `MIRA_RUN_CMDSTAN_TESTS=true`; three HTML/DOCX rendering tests require
`MIRA_RUN_QUARTO_TESTS=true`; two PDF rendering tests require
`MIRA_RUN_QUARTO_PDF_TESTS=true`. Their absence does not validate those
integrations. A separate public-API report probe attempted all three real
formats and recorded their current processx failure.

### Exact command/script record

Paths below identify the actual local run. Supporting scripts retain their
full commands, assertions and process-local settings; they are evidence, not
an alternative documentation infrastructure inside the package. Working
directory is `C:/Users/david/Documents/Codex/2026-10-29/il-x20` unless stated
otherwise. Common environment: `LC_ALL=C`, `LANG=C`, `R_LIBS=C:/Users/david/Documents/Codex/2026-10-29/il-x20/work/library`,
and `RSTUDIO_PANDOC=C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools`.

```powershell
# From the chat workspace: isolated package installation.
& 'C:/Program Files/R/R-4.4.3/bin/R.exe' CMD INSTALL --library=work/library --no-multiarch 'C:/Users/david/OneDrive/Desktop/MIRA/MIRA'

# From the repository, with MIRA_VALIDATION_LIBRARY set to that library
# and MIRA_RUN_ADVANCED_VIGNETTES=true:
& 'C:/Program Files/R/R-4.4.3/bin/Rscript.exe' --vanilla dev/validate_frequentist_vignettes.R 'C:/Users/david/Documents/Codex/2026-10-29/il-x20/work/final-delivery-vignettes'

# Corrected complete source-test driver, against the installed namespace:
& 'C:/Program Files/R/R-4.4.3/bin/Rscript.exe' --vanilla work/run-tests.R 'C:/Users/david/OneDrive/Desktop/MIRA/MIRA' work/test-results.rds

# Public reporting-contract and separate rendering probe:
& 'C:/Program Files/R/R-4.4.3/bin/Rscript.exe' --vanilla work/report-contract-probe.R

# Native partial-module investigation, separate from final tutorial values:
& 'C:/Program Files/R/R-4.4.3/bin/Rscript.exe' --vanilla work/probe_partial_modules.R
& 'C:/Program Files/R/R-4.4.3/bin/Rscript.exe' --vanilla work/partial_module_minimal_repro.R

# Package build, bundle installation, verification and CRAN-style check:
& ./work/build-check.ps1
```

`build-check.ps1` executes R CMD build from `work/final-build`, installs the
resulting `MIRA_0.0.1.tar.gz` into the isolated library, invokes
`verify-bundle.R` with repository/archive/library/output arguments, then
runs **`R CMD check --as-cran --no-tests`**. The PDF reference manual remains
enabled; `--no-manual` is not used. `_R_CHECK_FORCE_SUGGESTS_=true` is set;
`_R_CHECK_CRAN_INCOMING_REMOTE_=false` disables remote incoming checks.
`MIRA_RUN_CMDSTAN_TESTS=false` and optional extra advanced refits are disabled
for this package-check profile; the latter was separately executed above.
The check's `--no-tests` omission is explicit: it must not be described as a
full unmodified R CMD check. Code/test identity with the separately executed
suite is established by the final bundle/source proof.

The structural audit command uses the existing bundled Python and reads final
HTML without executing any statistical analyses:

```powershell
& 'C:/Users/david/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe' work/audit_html.py --source-dir 'C:/Users/david/OneDrive/Desktop/MIRA/MIRA/vignettes' --html-root work/final-delivery-vignettes --output-dir work/html-audit-delivery --figures-dir work/final-delivery-figures
```

The audit checks ten expected articles, flat-handbook relative links, literal
table captions, headers/row consistency, p-values displayed as numeric zero,
unintended math, image decoding/nonblank content and printed NULL contexts.
Conditional captions are reconciled with the module ledger where possible.
The immediately preceding structural pass found ten articles, seven valid
images, 46 table captions, zero structural errors/warnings and nine NULL
findings attributable to legitimate diagnostic fields. The final-source
rerun and hashes are included with the final package proof below.

All seven figures were actually reviewed through a contact sheet, with
focused full-image inspection of scientific figures. The root reviewer also
observed all ten article entry views and selected advanced-model/workflow
tables in the browser at a 1280-pixel viewport. The complete pages had no
horizontal document overflow in that view; a missingness table is slightly
wider than the nominal text column but remains within the viewport. The
handbook's index-to-first-tutorial link was followed successfully. This is
selected visual coverage, not a claim of every section/viewport/browser
combination. Structural checks cover the complete collection; no pkgdown
site was built. Local preview used loopback HTTP, without publishing a site
or bypassing a browser security warning.

### Final package build/check, identity, size and scope proof

| Final check | Observed result |
|:--|:--|
| `R CMD build` | Passed in 94.23 seconds; `MIRA_0.0.1.tar.gz`, 745,661 bytes (0.711 MiB). |
| Source-bundle installation | Passed in the isolated library; loading and installed-path assertions passed. |
| Installed vignette discovery | All ten expected names, HTML files and Rmd sources found. Eleven HTML files include R's generated index; extra future articles are permitted. |
| Source/bundle comparison | 97 audited package-input files validated: 96 byte-identical; original DESCRIPTION fields semantically identical after R's standard DCF whitespace normalization. Only `NeedsCompilation`, `Packaged`, `Author`, and `Maintainer` were added by the build. Every installed vignette Rmd matches the final source bytes. |
| `R CMD check --as-cran --no-tests` | Exit 0 in 232.01 seconds: **0 errors, 0 warnings, 1 NOTE**, `unable to verify current time`. Rebuilding vignette outputs passed (86 s); PDF and HTML reference manuals both passed. |
| Scope of the check | Package tests were explicitly skipped by this invocation and were run separately: 11,912 passing expectations, six gated skips, 393.38 s. The audited R/test files in the bundle are identical to the tested sources. `--no-manual` was **not** used. Remote incoming checks were disabled; repository-index access warnings in the console log limit remote availability verification. |
| Compiled documentation | `inst/doc` totals 1,059,686 bytes (1.011 MiB), including compiled articles and their sources. These measurements support a local resource assessment, not an acceptance decision. |
| Final HTML audit | Ten byte-final articles, seven decoded nonblank figures, 46 table captions; zero structural errors or warnings. Nine informational NULL findings are legitimate absence-of-error/skip/alternative fields, not demonstrations of missing analytical results. The six previously suppressed advanced tables are now visible. |
| API/source execution reconciliation | 43 public parameters; 50 static calls (48 routine and two optional); 393 explicit argument occurrences. All ten input Rmd hashes match the final execution record. The manual API matrix distinguishes actual execution, optional attempts and independent checks. |
| Final diff and whitespace | Fourteen task-authored source/documentation changes: the ten Rmd files, NEWS, AI usage record, validation runner and this new phase report. `git diff --check` and a whitespace audit of all phase-only files pass. |

The phase-only patch and SHA256 comparison preserve the existing working
tree rather than compare only against HEAD. `R/`, `man/`, `NAMESPACE`,
`tests/`, `inst/`, DESCRIPTION, README, both ignore files, `_pkgdown.yml`, and
the historical `dev/FREQUENTIST_VIGNETTES.md` retain their pre-revision source
hashes. Ignored RStudio state under `.Rproj.user/` changed during observation;
the agent did not edit or restore it, and it is excluded from the source patch.
The raw comparison marks it separately from task-authored files. No statistical
or Bayesian source was changed. The new developer report is excluded from the
package by the pre-existing `^dev$` build rule, so adding this audit record
does not change the validated tarball.

Evidence is delivered in `evidence/`: `bundle-source-identity.csv`,
`scope-comparison.csv`, `phase-only.patch`, `vignette-status.csv`,
`test-summary.csv`, `module_status.csv`, sample/formula/contrast ledgers,
`check.log`, environment records and the final HTML audit. The source tarball,
HTML/Rmd handbook and reproducibility scripts accompany this report. Script
paths describe this environment; adapt the workspace paths before rerunning
them on another machine.

## F. Open issues and limits

| Category | Current issue and impact | Resolution / remaining action |
|:--|:--|:--|
| Documentation defect corrected | Complete workflow previously substituted an adjusted endpoint arm contrast for its declared change estimand. | Corrected named interaction contrast and four-marginal-mean direction check; final integrated execution passed. |
| Documentation defect corrected | Multi-outcome extractor assumed a final-change row whenever the outcome analysis succeeded. | Return explicit unavailable endpoint row when n=0; execute final-unavailable example without plotting; reject unexpected duplicate selected rows. |
| Documentation omission corrected | RM-ANOVA was insufficiently distinguished from additive-covariate model alternatives. | Explain afex factorial/time-by-covariate mean structure and center numeric predictors; exact final results recorded separately. |
| Validation-assertion defect corrected | Whole-object identity after RDS can fail for native model environments/pointers even when scientific content is preserved. | Assert unchanged analysis object before serialization and exact relevant tables/settings/native parameters after serialization. Avoid traversing entire native environment graphs. |
| Software defect confirmed | With finite values but no baseline-follow-up pairs, requesting plots can raise `arguments imply differing number of rows: 0, 1`. | Reproduced through installed public API, minimal example below. Nonplotting availability analysis remains usable; engine repair is separate and was not attempted. |
| Software integration defect confirmed | Random-slope fit is available but supplementary `random_effects_lrt` fails the native same-data-object guard. | Identical subject-by-visit keys and a standalone native reproduction verified. Same-named-data native refits provide a diagnostic comparison, retaining boundary-hypothesis caveats; separate engine repair needed. |
| Native integration limitation confirmed | GEE fits and standardized joint tests are available, but native sequential ANOVA fails to find `geese.fit` in the namespace-only context for three working correlations. | Standalone reproduction and temporary-process `library(geepack)` workaround verified; recovered sequential tests must not replace MIRA's different joint hypotheses. No engine/global attachment change made. |
| Presentation defect corrected | Several significant `kable()` or model expressions were nested in a single conditional block and not visibly rendered. | Separate visible-output chunks/explicit printing now expose the returned scientific tables; source-only execution alone was insufficient verification. Final caption/HTML audit accompanies delivery. |
| Environment / rendering | Historical Windows Quarto/processx and TLS failures are not automatically current results. | HTML/DOCX/PDF were explicitly attempted in this phase; all failed at processx write-pipe creation with Access denied (system error 5), before rendering. Source-only examples passed; no bypass or tool installation.. |
| Installation / dependencies | Some engines described as optional modules are Imports in the package metadata. | Documentation now distinguishes installation requirements from runtime module gating. No dependency strategy changed. |
| Statistical limits | Synthetic thresholds are not clinical validation; model execution does not establish assumptions, MAR, causal identification, measurement-error calibration or cross-method equivalence. | Human scientific review and study-specific protocol/sensitivity work remain necessary. |
| Platform/version limits | The recorded installed Windows/R4.4.3 engines do not establish behavior on current R release/devel, Linux/macOS or every supported minimum/current engine version. | Independent cross-platform/current-version checks remain separate. |
| Unexecuted branches | External rendering/opening, defensive failed-outcome branch, or absent-dependency branches may remain unexecuted. | CmdStan smoke and five gated rendering tests skipped by their explicit flags; three separate real-render attempts failed environmentally. Document opening (`open=TRUE`), induced defensive failed-outcome extraction, and absent-dependency branches were not newly executed. Disabled/sparse modules are explained separately from errors and partial results.. |

The sparse plotting reproduction used in this phase is deliberately separated
from ordinary tutorial execution:

```r
sparse <- data.frame(
  participant_id = sprintf("S%02d", 1:6),
  instrument_t0 = c(12, rep(NA_real_, 5)),
  instrument_t1 = rep(NA_real_, 6),
  instrument_t2 = rep(NA_real_, 6)
)
safe <- mira_info_long(sparse, id = "participant_id",
  covariates = character(0), analyses = "none", verbose = FALSE)
# The same public input with analyses = "plots" reproduces the defect.
```

`work/plot-defect-repro.log` records the new reproduction. This is an engine
limitation rather than correct plotting behavior, a missingness remedy or a
reason to hide tutorial errors. No engine change was made to make an example
pass. The final unavailable-outcome extractor uses a valid nonplotting result.

## G. CRAN/JOSS assessment and primary sources

Official policy and engine material was consulted online on 2026-10-29;
`work/policy_engine_sources.md` retains the full source/claim record. Current
policy and engine documentation inform interpretation, while local function
versions and executed objects establish what MIRA actually did.

R's [Writing R Extensions](https://cran.r-project.org/doc/manuals/r-release/R-exts.html)
distinguishes dependencies, source vignettes, installed documentation and
package build/check. The [CRAN Repository Policy](https://stat.ethz.ch/CRAN/web/packages/policies.html)
requires portability and restrained check-time resources and places conditions
on strong dependencies. No fixed ten-minute package-check limit is asserted.
The policy's general documentation/source-size guidance should be assessed
using the final measured sizes and times in E, not historical figures.

The review found current CRAN entries meeting the declared minima for
[`nparLD >= 2.3.3`](https://stat.ethz.ch/CRAN/web/packages/nparLD/index.html) and
[`mmrm >= 0.3.17`](https://stat.ethz.ch/CRAN/web/packages/mmrm/index.html).
Current mmrm documentation is not proof that local mmrm 0.3.17 behaves
identically to the current release; no update was installed. Official
[CmdStanR installation guidance](https://mc-stan.org/cmdstanr/) and the
[publisher repository](https://stan-dev.r-universe.dev/cmdstanr) demonstrate
distribution through Stan R-universe. Because MIRA declares CmdStanR in
Imports, the authors must resolve its repository eligibility against CRAN's
strong-dependency policy before submission. `Additional_repositories` is not
itself a waiver or acceptance decision. This review did not establish an
exhaustive negative inventory of every possible CRAN/Bioconductor package
listing; the distribution/dependency concern is a documented inference requiring
maintainer confirmation, not a silent change to Bayesian dependencies.

The handbook improves reproducible installation/use documentation, scientific
questions, examples, API contracts, output interpretation and independent
checks. Those support some [JOSS review criteria](https://joss.readthedocs.io/en/latest/review_criteria.html)
and [checklist](https://joss.readthedocs.io/en/latest/review_checklist.html)
dimensions. Current [submission requirements and AI policy](https://joss.readthedocs.io/en/latest/submitting.html)
also address sustained public development, research use, contribution,
licensing/open practices, manuscript and human responsibility. This task did
not verify those broader requirements. No editor/reviewer message was drafted
or sent as part of submission correspondence.

Method-specific primary references supporting consequential distinctions
include [afex formula construction](https://search.r-project.org/CRAN/refmans/afex/html/aov_car.html),
[emmeans interaction contrasts](https://rvlenth.github.io/emmeans/articles/interactions.html),
[Kenward-Roger restriction matrices](https://hojsgaard.github.io/pbkrtest/reference/kr-modcomp.html),
[CR2 joint inference](https://jepusto.github.io/clubSandwich/reference/Wald_test.html),
[MMRM 0.3.17 methodology](https://openpharma.github.io/mmrm/v0.3.17/articles/methodological_introduction.html),
[nparLD's authors' vignette](https://stat.ethz.ch/CRAN/web/packages/nparLD/vignettes/nparLD.html)
and [`p.adjust`](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/p.adjust.html).
The full source record covers the remaining engines. No invented DOI or
attribution is used.

### AI assistance and outstanding human work

OpenAI Codex, using a GPT-6-family model, assisted with source/test inspection,
scientific/editorial revision, synthetic examples, extraction, validation
scripts, execution, policy/engine research and reporting. The exact runtime
deployment identifier was not exposed and is not invented. Parallel AI agents
reviewed distinct article groups and conducted second readings across groups.
This was AI cross-review, not human peer review. Tools included local
PowerShell/R/package tools, source/hash/diff inspection, web retrieval of primary
documents and image inspection of generated figures. No site deployment,
publication or communication with external reviewers occurred.

Human authors remain responsible for reviewing/editing accepted prose,
confirming statistical interpretations and design choices, independently
checking consequential computations, assessing study-specific applicability,
and preparing truthful manuscript/AI disclosures. This report does not attest
that those human responsibilities have been completed. Public development
history, real research adoption, comparative need, licensing/contributor
permissions and manuscript readiness require separate evidence.

## H. Concluding implementation status

**Implemented with explicit validation limitations** (the requested status
“implementato con limitazioni di validazione esplicite”).

The ten existing articles were substantively revised, all final-source
articles executed successfully, the applicable optional advanced branch ran,
and the complete existing test suite passed with six explicitly gated skips.
The final package build/install/discovery/check and source-scope proof are
recorded in E. Successful execution is not presented as proof of statistical
validity or formal CRAN/JOSS acceptance.

The three confirmed software/integration limitations remain disclosed:
sparse-data plotting, the random-effects comparison, and native GEE ANOVA.
The latter two leave four model containers partially usable despite
`performed=TRUE`. Real Quarto rendering was attempted but remained blocked
by Windows processx pipe creation. Current-R, cross-platform, dependency-
eligibility and human scientific/editorial review remain separate work.
No statistical engine or Bayesian component was changed to make a tutorial
or check pass.

The new technical report is maintained separately from historical evidence in
`dev/FREQUENTIST_VIGNETTES_REVISION_2026-10-29.md`. User-facing deliverables
include the revised Rmd/HTML handbook, API matrix and current execution,
module, sample, formula, contrast, export/render and test evidence. Temporary
investigation artifacts remain outside the repository; deliverables are
collected in this chat's `outputs/` directory. No commit, push, PR, tag,
release, deployment, publication or pkgdown build was performed.
