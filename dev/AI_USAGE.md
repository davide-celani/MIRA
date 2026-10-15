# AI assistance record

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
