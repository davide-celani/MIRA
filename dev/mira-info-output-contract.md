# Longitudinal output contract

This note records the scope, ownership and evidence behind the output of
`mira_info_long()`, `mira_detect_long()` and `mira_report_freq_long()`.
The executable naming and ordering specification is
`.mira_metric_families_long()` / `.mira_metric_spec_long()` in
`R/mira_info_long.R`. Statistical calculations precede output normalization.

## Name provenance

`config$outcome_names` records each outcome's key, exact declared name (or
`NULL`), source and display label. `outcomes` declaring a logical outcome and
nonempty names of `time_vars` list elements declare names. Column selection
through either argument, automatic detection and `variable_pattern` infer
names. A partially named list is handled per element. Selectors such as
`NULL` and `"auto"` are not declarations. Wrapper forwarding does not promote
an inferred name to a declaration. A display label is never the metric prefix.

Every registered metric uses its canonical base name when inferred, and
`paste0(declared, "_", base)` when declared. Case, underscores and words inside
the declared name are preserved. For example, the same rule gives
`Mobility_mean`, `daily_change_mean` and `mean` for the inferred case.
Matching-key collisions between declarations and collisions between final
metric fields and existing fields raise explicit errors before overwriting.
Names never depend on observed values.

`attr(table, "mira_metric_names")` maps base names to actual fields; it holds
references, not duplicate values. `attr(result, "mira_metric_paths")` bounds
traversal to registered components. Print methods use a transient base-name
view, and the report resolves fields through the map. Stored objects, CSVs and
RDS payloads retain canonical names.

## Component and consumer map

| Components | Producer / ownership | Treatment and consumers |
|---|---|---|
| `call`, `config`, `settings`, identities, time variables/labels | public resolution and configuration helpers | Context; original source names and declared outcome order retained. Detection, inspect, printing and report settings consume these. |
| `data_overview`, `detected_variables`, `diagnostics`, `overview` | detection/configuration and single-outcome analyzer | Dataset profiles, selection diagnostics, sample sizes, completeness and adaptation metadata stay unprefixed. No traversal into user data. |
| `missing`, descriptive availability fields, arm counts | single-outcome analyzer | Sample/availability metadata (`n`, `n_a`, `n_b`, `n_subjects`, missing/nonfinite/unavailable counts and percentages) stay generic. Counts and their percentages remain paired. |
| `descriptives`, `change`, arm descriptives/comparisons/missingness tests | single-outcome analyzer | Owned estimates, dispersion, intervals, effect sizes, response counts/percentages, test statistics, degrees of freedom and p-values are outcome metrics. Print, report conclusions, tables and exports use their field maps. |
| `variability`, `trajectories`, `outliers` | single-outcome analyzer | ICC/variance diagnostics, trajectory change magnitudes and diagnostic bounds are metrics. Subject IDs, rows, source observed values and direction/clinical-direction states remain unprefixed and aligned. |
| `correlations` | `stats::cor` plus pairwise sample counts | Method containers and matrices retain their source-variable dimensions; these are not metric column labels. |
| `model` | lme4/lmerTest native fits, summaries, ANOVA and tests | Native contracts remain unchanged; model status flags, warnings and design counts are metadata. |
| `advanced_tests` | emmeans, RM-ANOVA and Friedman helpers | MIRA standardized/adjusted tables receive metric naming. Native emmeans means and Tukey/Dunnett summaries, `emmGrid`, afex objects/tables, sphericity and `htest` objects keep their native classes, fields and attributes. |
| `advanced_models` | random-slope, nlme, GEE and comparison helpers | Owned extracted coefficient tables, GEE effect tests, model-fit scalars and comparison tables receive metric naming. Native model/summary/ANOVA objects, nlme coefficient matrices, variance-component tables and QIC vectors retain their contracts. |
| `robustness` | clubSandwich helper | Native covariance matrices, coefficient/interval/Wald tables retain their external contracts. Hypothesis diagnostics, omitted coefficients and failure explanations remain intact. |
| `effect_sizes`, `sensitivity` | advanced-inference helpers | Owned effect sizes, R-squared, preferred model ICC and sensitivity statistics receive metric names. Method, question, effect-size type, source path and sample size are metadata. |
| `multiplicity` | advanced-inference helper | Adjustment metadata and canonical table `object_path`; no second copy of the result table. Raw/primary/individual adjustment results remain at their owning contrast table. |
| `plots`, `long_data` | single-outcome analyzer | Native ggplot contracts and source response/covariate columns are preserved. The derived `patient_mean` support metric follows naming. Plot rows and IDs are not reordered. |
| multiple outcomes / failures | public wrapper | One ordered result per outcome; failed outcomes retain their error class and messages. A single result has no duplicate `$outcomes` copy. |
| report source, runtime, manifest, CSV, RDS | report builder and dumped runtime helpers | Canonical names/order flow through export; rendering resolves names from maps. User `extra_objects` and environment/package/dependency versions are preserved. |

## Proven duplicate inventory

Removal depends on assignments in producer code, not numeric equality.

| Path | Removed fields | Canonical field / evidence |
|---|---|---|
| `descriptives` | outcome aliases of mean, SD, median, CI endpoints | Direct copies of `mean`, `sd`, `median`, `ci_lower`, `ci_upper`; these bases now follow the naming rule. |
| `change` | outcome `change`, `mean_from`, `mean_to` aliases | Direct copies of `mean_change`, `mean_from`, `mean_to`; the full base name is preserved. |
| `variability` | outcome grand-mean/between-SD/within-SD aliases | Direct assignment from the corresponding generic fields. |
| `trajectories` | outcome baseline/final/change/relative-change aliases | Direct copies of `baseline`, `final`, `absolute_change`, `relative_change_percent`. |
| `advanced_tests$emmeans` | top-level `object`, `summary`, `tidy` | Same `time_emm` / `time_tidy` variables stored in `time$object` / `time$tidy`. |
| `advanced_tests$friedman` | `object`, `summary` | Same `test` / `tidy` variables retained as `test` / `tidy`. |
| `advanced_tests$rm_anova` | `n` | Both count fields assigned `rm_subjects`; retain `n_subjects`. |
| random-slope, nlme, GEE modules | `object` | Same fitted model assigned twice; retain `model`. |
| random-slope module | `tidy` | Both fields assigned `fixed_eval$value`; retain `fixed_effects`. |
| GEE modules | `robust_standard_errors` | Exact subset of `coefficients`, with no extra method/sample/context. |
| GEE owned coefficient tables | source SE column | `robust_se` is its numeric copy; native model summary retains its original contract. |
| `robustness$club_sandwich` | `object`, `summary`, `tidy` | Same covariance and coefficient evaluation objects; retain `covariance`, `coefficient_tests`. |
| R-squared modules | `summary` | Same `evaluated$value` as `object`; retain native `object` and standardized `tidy`. |
| owned multiplicity-adjusted tables | source raw p-value column | Numeric source copied to `p_raw`; raw and adjusted methods remain distinct. |
| multiplicity families | `table` | Same canonical contrast table passed directly; retain reference `object_path`. |
| single-outcome top level | `outcomes` | The analyzer result was repeated inside its own one-element outcome list. |

Consumers affected are the print/summary methods, report conclusions and table
formatting, report export/runtime, examples and schema/statistical tests.
There are no compatibility result aliases. The report no longer discards
arbitrary equal-value tables.

## Deliberately retained information

The preferred `ICC` selects a model estimate with ANOVA fallback, so neither
`ICC_model` nor `ICC_anova_complete_profiles` is its unconditional equivalent.
Raw, selected-primary and individually adjusted p-values describe different
inferential roles even when equal. SD, variance, SE, counts/percentages and
numeric/clinical direction remain distinct. Effect-size summaries and baseline
balance retain context useful for interpretation; native model summaries remain
available beside standardized extracts. IDs and metadata repeated across
tables make those tables interpretable on their own. No constant/NA/value-based
deduplication is performed. Selected versus detected covariates retain different
configuration roles even when manual selection makes them equal. Identical
skipped-module placeholders remain under distinct method names, so unavailability
does not erase which analyses were requested.

## Ordering and structural boundaries

Top-level order is context/configuration, data overview and availability,
descriptives, changes and trajectories, arm comparisons/correlations/variability,
models, effect sizes and multiplicity, robustness/sensitivity/outliers and
diagnostics, then plots/support data. Multi-outcome containers follow context
and diagnostics and preserve the user's outcome order. Paths are unchanged
except proven duplicate removals.

The report keeps arm comparisons immediately after overall changes, so readers
can compare population and arm results before individual trajectories. This
small presentation-order difference leaves the stored component order intact.

Column order is specified by family using canonical bases: identifiers and
time/comparison fields first, sample availability, estimates/dispersion and
SE/intervals (lower before upper), effect sizes, paired response counts and
percentages, and test results with raw/adjusted p-values adjacent. Prefixes do
not affect ordering. Rows, factor levels, reference arms, comparison orientation
and subject associations are retained. No absent module or missing result is
invented; NULL, empty tables and failed modules retain their meaning.

Function/output-schema version fields and historical revision labels are
removed from detection, analysis, summaries, report-generated content and
directly related documentation. Package/dependency/environment versions and
native or arbitrary user metadata are outside this removal.

## Verification strategy

Testthat covers provenance and arbitrary names, exact family schemas, semantic
deduplication, equal/NA-valued distinct statistics, empty/failed/disabled cases,
native boundaries, printing/summary, exports and isolated report runtime.
Statistical oracle tests resolve field maps before applying their original
numerical assertions. Before/after comparisons use saved working-source
results with fixed data/seeds and explicit removal/rename/permutation maps.
Bundle generation with `render = FALSE, open = FALSE` is separate from final
Quarto document rendering; isolated runtime/knitr tests exercise presentation
helpers without claiming Quarto/PDF rendering.
CSV round-trip tests allow `1e-12` for decimal serialization of floating-point
values; RDS and before/after scientific comparisons use exact values.
