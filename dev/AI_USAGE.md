# AI assistance record

## Ophthalmic frequentist examples: 2026-10-30

OpenAI Codex (GPT-6 family; exact deployment identifier not exposed) and three
parallel AI agents redesigned all eleven frequentist articles around new
synthetic BCVA, defined OCT central-subfield thickness and IOP studies. Work
covered study-eye definitions, measurement scales and directions, generators,
clinical covariates, named sensitivity/edge cases, interpretations, figures,
report bundles and primary/official ophthalmic measurement sources.

All 33 final-source execution modes passed. Fresh unchanged-suite results,
package build/install/check, native model/sample proof, 71-table/1,665-cell
and nine-figure/caption reviews, source/archive/installed identities and
clinical/API/scientific matrices are recorded separately in
FREQUENTIST_VIGNETTES_OPHTHALMOLOGY_2026-10-30.md. Historical reports and
pre-existing local edits retain their original evidence. Statistical R,
public signatures, tests, Rd, dependencies and Bayesian/Stan material did not
change. External Quarto and browser-inspection limits are explicitly recorded.

All data and effect sizes are invented for teaching; illustrative bands are
not validated MCIDs or clinical decision rules. The clinical and statistical
cross-reviews were AI-assisted, not human clinical review. Maintainers remain
responsible for scientific/editorial review and publication disclosures. No
commit, staging, push, PR, publication, release or pkgdown build was performed.

## Frequentist handbook closure: 2026-10-30

OpenAI Codex (GPT-6 family; exact deployment identifier not exposed) and three
parallel AI agents completed a documentation-only closure phase. Assistance
covered reader-first restructuring, a new first-model article, visible executable
helpers, focused missingness/reference/sample/transfer cases, API and scientific
matrices, isolated per-article/code validation, generated-table/image review,
package build/install/check and primary-source policy/engine research.

The final collection has eleven articles. Native sparse plotting, random-slope
comparison and GEE ANOVA limitations remain disclosed rather than patched.
External Quarto formats failed at Windows processx pipe creation; browser file
URL inspection was blocked. These limits are separate from passing ordinary
vignettes, source-only archives, unchanged tests and package-manual checks.
The new report is `FREQUENTIST_VIGNETTES_CLOSURE_2026-10-30.md`; historical
reports retain their original evidence. Statistical/Bayesian code, public
signatures, tests and dependency strategy were preserved. No commit, staging,
push, PR, publication, submission, deployment or pkgdown build was performed.

Cross-review was by AI agents, not human researchers. Human scientific,
editorial and usability review, responsibility for accepted outputs and a
truthful manuscript AI disclosure remain pending. No JOSS/CRAN acceptance or
human-review attestation is made.

## Substantive frequentist handbook revision: 2026-10-29, second phase

OpenAI Codex (GPT-6 family; exact deployment identifier not exposed) and three
parallel AI agents reviewed the existing ten-article handbook against the
analysis and reporting implementations, tests, and current primary sources.
They made targeted scientific and editorial revisions, added independent
sample/interval/ICC/contrast/export checks, audited rendered output, and
strengthened the developer validation runner. AI review is not human review.

This phase corrects the primary adjusted contrast in the complete workflow,
distinguishes afex's expanded mean structure, and documents safe extraction
when an otherwise valid outcome has no observed endpoint pairs. Statistical
source code, public signatures, tests, and Bayesian material were preserved.
The earlier implementation report below remains historical evidence. The
new phase's commands, final-source hashes, execution results and limitations
are recorded separately in `FREQUENTIST_VIGNETTES_REVISION_2026-10-29.md`.

The authors must still perform and attest their own scientific/editorial
review, take responsibility for design and interpretation, and reconcile the
paper's AI disclosure with the current JOSS policy. This record supplies no
human-review attestation or formal CRAN/JOSS acceptance judgment. No commit,
push, PR, release, deployment or pkgdown build was performed in this phase.

## Frequentist longitudinal handbook: 2026-10-29

OpenAI Codex assisted with inspection of the public longitudinal analysis and
reporting implementation and the design, writing, execution and review of a
ten-article R Markdown handbook. Assistance covered deterministic synthetic
studies, scientific explanations, metric-name-safe extraction, module-status
and report-artifact checks, standard vignette/pkgdown integration, and a
developer validation script. Parallel agents reviewed distinct tutorial
groups and reporting claims. The exact runtime/deployment model identifier
was not exposed and is not invented here.

The implementation report records the actual package build, installation,
vignette discovery, examples, tests, extended model execution and checks,
including environmental rendering and website limitations. Statistical
source files, function signatures, Bayesian examples and tests were preserved.
An existing sparse-data plotting defect was reproduced and documented rather
than changing statistical behavior to make a tutorial pass.

Human scientific and editorial review remains pending. This record does not
assert that the authors have completed the human verification required by
JOSS, nor that package checks establish scientific validity or eligibility.
The submitting authors must reconcile the eventual paper disclosure with
their own review, modifications and validation of the accepted material.
See `FREQUENTIST_VIGNETTES.md` for scope, commands, results and limitations.

## CmdStan CI timeout correction: 2026-10-26

OpenAI Codex inspected the supplied setup-r-dependencies log and public GitHub
run 40307665205 for commit 5de4bb3. GitHub explicitly reported the job's
30-minute limit was exceeded; setup-r consumed about 20 minutes and the
dependency step was cancelled before CmdStan installation or tests ran.
Assistance increased the job budget to 60 minutes, restricted dependency
types to the same check profile used elsewhere, and replaced the unnecessary
devtools wrapper with testthat::test_local using the same filter and failure
setting. Local validation is distinct from an observed successful rerun of
the modified workflow: YAML/R syntax and configuration checks passed, and the
same local test filter passed 161 expectations with zero failures/warnings
and one explicit CmdStan skip (the local sampling flag was false). The CI
sampling flag remains true. No commit, staging, push or remote rerun was performed;
human review remains pending.

## Longitudinal test audit: 2026-10-26

OpenAI Codex assisted with the requested local audit of `mira_info_long()` and
`mira_report_freq_long()`, using parallel statistical, report and suite reviews.
Assistance covered complete test-suite inventory, independent contrast oracles,
minimal fixes for RM-ANOVA availability, report comparison keys and empty-input
validation, PDF verifier sensitivity, RNG isolation, and repair of the scientific
and bundle validation scripts. Source and isolated installed-package tests,
build/check, explicit rendering profiles and controlled corruptions are recorded
in `TEST_AUDIT_LONG.md`, including environmental failures and remaining limits.
No commit, staging, push, release, PR or external communication was performed.
No new dependency was installed. Human review remains pending; this record
does not assert human verification, universal statistical validity or a JOSS
editorial outcome. Earlier dated reports retain their own provenance.

## Bibliographic identity session: 2026-10-14

OpenAI Codex inspected local citation metadata and existing release workflows,
consulted official CFF, R and Zenodo documentation, and added Davide Celani's
user-provided ORCID to CFF and R author metadata. Assistance included a local
consistency/checksum validator, its integration in the existing R check job,
and documentation of persistent identity versus release-specific metadata.
No scientific algorithms were changed. No staging, commit, push, tag, release,
remote configuration change or Zenodo deposit was performed. Human review
remains pending; private Zenodo integration settings and actual ingestion on
the next release are not certified by the local validation.

## Validation session: 2026-10-12

OpenAI Codex assisted with an audit of the local implementations, complete
test suite, dependency declarations, documentation and CI for the requested
validation of `mira_info_long()` and `mira_report_freq_long()`. The exact
runtime model identifier and deployment version were not exposed; this record
does not invent them.

Assistance included identifying invalid-argument tests that passed on an
unrelated missing-symbol error, correcting those QTE tests and their RNG
cleanup, replacing broad warning suppression with explicit checks, adding
hand-derived data-mapping and prior-conversion oracles, and preparing CI
profiles for offline tests and explicit rendering integrations. Codex also
consulted the official JOSS submission and review guidance. The validation
report records the other test and implementation changes, commands, outcomes
and limitations of this session; a configured workflow is not a completed CI
run, and mocked rendering is not real rendering.

Human review remains pending. The maintainer must review and validate all
AI-assisted outputs, their scientific interpretation and the reported
evidence before accepting the working-tree changes or making a submission.
This record does not assert human validation, public development history,
external adoption, or an editorial verdict. No submission or communication
with editors or reviewers was performed.

## Release-preparation session: 2026-10-03

OpenAI Codex assisted with review of the existing package and preparation of
release documentation and metadata in this session. The exact runtime model
identifier was not provided; this record does not assign a model name or
version.

The recorded scope includes drafting contribution guidance, the project Code
of Conduct, the initial changelog, GitHub issue forms, and this disclosure.
It also includes assistance with review and preparation of release metadata
and citation information, rewriting the README, drafting release and JOSS
preparation guides, and creating GitHub workflows. Code assistance comprised
a writable temporary compilation cache in `mira_fit_long()`, its regression
tests, expansion of the existing opt-in smoke test to all three likelihoods,
and package namespace/import cleanup. Follow-up assistance corrected the
R CMD check workflow input quoting, normalized accented report filenames for
macOS, and allowed machine-precision rounding in a correlation test.
After publication, Codex verified the public Zenodo record and both DOIs,
then updated citation metadata, DOI links, and release documentation in a
follow-up source commit.
Existing tests, examples, documentation
generation, and R package checks were executed and their results reviewed by
Codex. The maintainer should update this scope to match the final changes
accepted for publication.

This record describes assistance in this session only. It does not establish
the development history of the pre-existing R or Stan implementation, the
origin of earlier documentation or tests, or AI usage in a manuscript.

Human review is pending. The maintainer must verify descriptions against the
implementation, confirm authorship and citation metadata, and review the
reported validation results before publishing. No scientific validation or
successful package check should be inferred from AI assistance alone.

For future contributions, record the tool, date, scope, and human verification
when AI assistance materially contributes to a change. Do not list an AI
system as an author, and do not upload private participant data to an AI tool
without appropriate authorization.

## Longitudinal methodology expansion: 2026-10-26

OpenAI Codex assisted with primary-source method research and implementation
of joint Kenward-Roger tests, unweighted relative marginal effects with
asymptotic ATS/WTS, and an unstructured REML/Satterthwaite MMRM. Assistance
included dependency declarations, adaptive sample/estimability guards,
print/summary and report integration, independent mathematical test oracles,
mutation checks, documentation, and local source/installed package validation.
The existing longitudinal methods were retained. The evidence and limitations
are recorded in `dev/TEST_AUDIT_LONG.md`; its earlier audit counts describe
the preceding implementation, not the expanded suite.

Synthetic deterministic fixtures were used for statistical checks. No private
participant dataset was supplied for this expansion. Human scientific review
and acceptance of the local changes remain pending. No commit, staging,
push, publication, or communication with third parties was performed.
