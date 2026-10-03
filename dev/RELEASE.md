# Maintaining citable MIRA releases

Version `0.0.1` was published on 2026-10-03 as tag `v0.0.1`, targeting commit
`40ad4af0ed13bbd0fedd2efd4fd583fdaeaf0fac`, and is archived on
[Zenodo](https://zenodo.org/records/23120184).
Its version DOI is `10.5281/zenodo.23120184`; the concept DOI for all MIRA
versions is `10.5281/zenodo.23120183`.
The software title is **MIRA: Bayesian Multilevel Inference for Longitudinal
Data**, authored by Davide Celani and licensed under MIT.

New commits on `main` continue development without changing the published
tag or the archived snapshot. Adding DOI links and citation metadata in a
follow-up commit does not require another release.

## Prepare the next software version

Update `DESCRIPTION`, `CITATION.cff`, and `NEWS.md` together when preparing a
new version. Remove the previous release date and version DOI from the CFF,
and remove the previous version DOI from the R citation
in `inst/CITATION`; use the repository URL while the new archive is pending.
Do not pair a new version number with an older archive's DOI. Keep the concept
DOI for the general README link or badge. After publication, add the new
verified version DOI in a follow-up commit and preserve the published tag.

The steps below describe the release workflow; replace `0.0.1` and `v0.0.1`
with the version being prepared. The first release's identifiers above remain
its historical record.

## Validate the release candidate

1. Confirm the version in `DESCRIPTION` and `CITATION.cff` matches the planned
   tag `v0.0.1`. Keep the changelog entry unreleased until the publication date
   is fixed, then record that date and set `date-released` in the CFF.
2. Review the README against the public API, run its synthetic example, and
   verify installation from a clean R library.
3. Run `devtools::document()`, `devtools::test()`, and
   `devtools::check(args = "--no-manual")`. Resolve errors and warnings and
   explain any remaining notes. Check the GitHub workflows on each configured
   platform after pushing the preparation commit.
4. Configure CmdStan and run the opt-in smoke tests. Validate meaningful
   scientific behavior separately through simulation and reference calculations;
   smoke tests are only execution checks.
5. Run `cffconvert --validate` from the repository root. Confirm
   `citation("MIRA")` for an installed package returns the same author/title.
6. Review every changed file with `git diff` and exclude local reports,
   compiled models, private study data, credentials, and session files.
7. Review the recorded AI assistance in `dev/AI_USAGE.md` before publishing.

An unresolved failed check means the candidate still needs work. Do not mark
the release as validated solely because metadata has been added.

## Enable Zenodo before publishing the GitHub release

1. [Create a Zenodo account](https://help.zenodo.org/docs/get-started/create-an-account/).
   Use GitHub to connect the account if convenient. Add an ORCID only when it
   belongs to the author and its value has been verified.
2. In Zenodo's GitHub settings, authorize GitHub access and
   [enable the MIRA repository](https://help.zenodo.org/docs/github/enable-repository/).
   Confirm it is enabled before publishing the first release.
3. Use `CITATION.cff` as the metadata source. The present release does not need
   `.zenodo.json`. If both exist, Zenodo uses `.zenodo.json` and ignores the CFF.
   See [Zenodo metadata guidance](https://help.zenodo.org/docs/github/describe-software/citation-file/).

## Publish and verify the archive

1. Commit and push the reviewed release candidate and wait for CI results.
2. Publish a GitHub release using tag `v0.0.1`, targeted at the exact verified
   commit. Use a title such as `MIRA 0.0.1` and a concise description based on
   `NEWS.md`, including implemented capabilities and current limitations.
3. Zenodo archives the new release. Follow the
   [GitHub upload guide](https://help.zenodo.org/docs/github/archive-software/github-upload/)
   and inspect any reported metadata errors.
4. Open the Zenodo record and verify the exact title, author, MIT license,
   version, source URL, and archived files. Record both the version DOI and
   the concept DOI.

For a paper's reproducible analysis, cite the DOI of the version actually used.
The concept DOI identifies MIRA across versions. See
[Zenodo DOI versioning](https://support.zenodo.org/help/en-gb/1-upload-deposit/97-what-is-doi-versioning).

## Add the real DOI after publication

Add the verified version DOI to `CITATION.cff` under `doi` and to the
`bibentry()` in `inst/CITATION`, then update the README citation and DOI badge.
These additions go into a subsequent source commit; the original archive and
tag remain the release's historical snapshot. Zenodo already provides a DOI
citation for that archived snapshot. Do not move the published tag merely to
insert its newly assigned DOI.

When the R citation is eventually changed to a `Software` BibTeX entry or a
preferred software citation, check that the target journal accepts it. Until
then the existing `Manual` entry is compatible with R's citation mechanism.

## Earlier analyses

If an article used an earlier MIRA model, identify and preserve that exact
source separately before claiming full reproducibility. A new release DOI
documents the new snapshot; it does not establish that earlier analyses used
the identical code. The package citation and the analysis archive may therefore
be separate references.
