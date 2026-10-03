# Preparing MIRA for a future JOSS submission

The first Zenodo archive can be created independently of JOSS. This file tracks
the additional work for a software paper; it is not a submission or a claim of
eligibility. Check the [current JOSS requirements](https://joss.readthedocs.io/en/latest/submitting.html)
again when preparing the submission.

## Public development and research use

- Maintain more than six months of public, active development, with meaningful
  iteration over time. The earliest local commit is dated 2026-08-19; a commit
  date alone does not prove when the repository became public. Record the
  verified public date before planning submission.
- Preserve releases, changelog entries, issues, CI results, and contribution
  guidance. A sole-author project can demonstrate open maintenance without
  manufacturing contributions or activity.
- Document actual research use and the exact versions involved. Distinguish
  the earlier model used by a manuscript from a new package release.
- Obtain independent installation and usability feedback from a colleague.

## Scientific and technical work

- Execute compilation and sampling checks for all supported likelihoods.
- Validate parameter and contrast recovery in simulated datasets, including
  multiple time points, multiple arms, covariates, and bounded outcomes.
- Check responder estimands and uncertainty against independently derived
  calculations; document the external justification of clinical thresholds.
- Retain clear documentation of the bounded latent-location estimand and
  its distinction from an observed censored marginal mean.
- Validate Bayesian summaries beyond mocked fits and short smoke chains.
- Add a complete worked tutorial with diagnostics and prior sensitivity.
- Verify supported systems through GitHub Actions and test clean installation.
- Review internal helpers that still assume exactly three visits, notably
  `R/data_validation_long.R`, before extending or exposing those helpers.
- Add eye/patient nesting or crossover-specific structure only through an
  explicit design, implementation, and validation effort when needed.

## Software paper and disclosure

Prepare `paper/paper.md` and `paper/paper.bib` once scope, research evidence,
affiliation, and references are confirmed. The paper should explain the need,
alternatives, design decisions, and demonstrated research impact. Use
`dev/AI_USAGE.md` to assemble an accurate AI disclosure after human validation.
Follow the [JOSS review criteria](https://joss.readthedocs.io/en/latest/review_criteria.html)
and their current paper instructions.

Keep the software contribution separate from clinical findings. Reviewers will
evaluate the implemented software and its evidence, not prospective promises.
