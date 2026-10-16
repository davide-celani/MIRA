# Reproducing the longitudinal validation

This guide accompanies `VALIDATION_REPORT_2026-10-16.md` and preserves the
historical procedures in `VALIDATION_REPORT_2026-10-12.md`. Recorded results
belong to the report's stated working tree, not automatically to current files.
Run from the repository
root with a supported R installation and **all DESCRIPTION Imports** already
installed. `testthat` and `withr` are needed for tests. These scripts do not
install packages, Quarto, TeX or CmdStan and do not need network access.
Choose output directories **outside the repository**, especially for mutation
copies. Do not point mutation output at an ancestor of the repository.

## Offline suite: source and installed package

The following PowerShell example uses R 4.4.3 and an existing dependency library.
Replace paths for another machine. `LC_ALL=C` and `LANG=C` avoid unsupported
Windows locale startup warnings; they are set only in this shell process.

```powershell
$rscript = 'C:\Program Files\R\R-4.4.3\bin\Rscript.exe'
$r = 'C:\Program Files\R\R-4.4.3\bin\R.exe'
$repo = (Get-Location).Path
$audit = 'C:\path\outside-repository\mira-validation'
New-Item -ItemType Directory -Path $audit -Force | Out-Null
$env:LC_ALL = 'C'
$env:LANG = 'C'
& $rscript --vanilla dev/validate-tests.R $repo "$audit/source" source
```

`validate-tests.R` explicitly disables all three external-engine integration
flags. It saves per-block CSV/RDS results, installed-package versions, timings
and a session-state audit; it exits nonzero on failures, errors, unexpected
warnings or changes to the recorded RNG, directory, graphics devices,
environment or existing options. Inspect added options separately: first
loading Imports can initialize their own package defaults. The summary
reporter's random praise is disabled so it does not consume the audited RNG.

Build and install a source archive, then run tests through the installed namespace:

```powershell
Push-Location $audit
& $r CMD build --no-build-vignettes --no-manual $repo
New-Item -ItemType Directory -Path "$audit/library" -Force | Out-Null
& $r CMD INSTALL "--library=$audit/library" MIRA_0.0.1.tar.gz
$env:R_LIBS = "$audit/library"
& $rscript --vanilla "$repo/dev/validate-tests.R" $repo "$audit/installed" installed
$env:MIRA_RUN_CMDSTAN_TESTS = 'false'
$env:MIRA_RUN_QUARTO_TESTS = 'false'
$env:MIRA_RUN_QUARTO_PDF_TESTS = 'false'
& $r CMD check --no-manual --no-vignettes MIRA_0.0.1.tar.gz
Pop-Location
```

Use the archive name produced by the current DESCRIPTION version. Check each
process exit code before proceeding. Installing into a temporary library does
not remove the need for Imports in the existing R library paths. This profile
does not simulate missing mandatory dependencies. It does not build the PDF
manual or perform the network-dependent CRAN incoming checks.
Standard R CMD check can still try configured repository indexes during
dependency checks. Unavailable-index messages are console warnings recorded
separately from test warnings and the formal final check status; they do not
trigger package installation. The ordinary test runner requires no network.

## Independent scientific checks

The ordinary suite includes hand-derived paired summaries, tie-corrected
Friedman statistics, independent cluster-sandwich arithmetic for independence
GEE, explicit Stan observation mapping and prior-conversion formulas. Model
engine comparisons validate wrapper behavior against external engine outputs;
they do not independently validate the external engine. Metamorphic checks freeze configuration and
compare aligned quantities under ID/row/column permutations, affine outcome
changes, sign/direction reversal and an entirely missing participant. Twelve
small generated cases use seed 20261012, restore RNG and report configuration.

The extended Monte Carlo study is separate from rapid tests:

```powershell
New-Item -ItemType Directory -Path "$audit/science" -Force | Out-Null
& $rscript --vanilla dev/validate-science.R "$audit/science"
```

It uses the installed package: 24 independent participants, normal baseline
and independent normal change, targets 0 and 2, 500 replications each, seed
2026101201. Every estimate and t interval is first checked against independent
arithmetic at absolute tolerance 1e-10. Bias, interval coverage and rejection
are reported with Monte Carlo standard errors. A deviation exceeding four
specified MCSE is an investigation signal, not a guarantee or a fast-suite
pass threshold. No result establishes adequacy under informative missingness
or validates Bayesian inference, clinical thresholds or all model engines.

Deterministic arithmetic uses approximately 1e-12 tolerances; fitted-model
invariance and sandwich comparisons use 1e-10 where observed numerical
precision permits. Existing engine-specific tolerances are retained. CSVs
preserve raw values; display rounding is tested separately.

## Report bundle, relocation and real rendering

```powershell
New-Item -ItemType Directory -Path "$audit/bundle" -Force | Out-Null
& $rscript --vanilla dev/validate-bundle.R "$audit/bundle"
```

Use a fresh output directory. This generates an actual analysis, checks a
hand-derived final contrast, full numeric CSVs and relative manifest, moves
the bundle, and knits the QMD with a fresh environment whose parent is
`baseenv()`. It checks analytical text and actual captured figures. This is
real **knitr execution**, and is not a Quarto HTML/DOCX/PDF render.

For the public Quarto integration, install Quarto separately beforehand:

```powershell
New-Item -ItemType Directory -Path "$audit/quarto" -Force | Out-Null
# Optional when the CLI is outside PATH: point at the existing executable.
$env:QUARTO_PATH = 'C:\path\to\quarto.exe'
$env:MIRA_RUN_QUARTO_PDF_TESTS = 'false'
& $rscript --vanilla dev/validate-rendering.R $repo "$audit/quarto"
```

HTML and DOCX are explicitly enabled by this script. To include PDF, put an
existing TeX engine on PATH and set `MIRA_RUN_QUARTO_PDF_TESTS=true` before
running into a separate output directory. Automatic TeX installation/update
is disabled. Enabled integrations must fail, rather than skip, when required
tools are unavailable. Tests inspect scientific text, table structure and
format-specific escaping, not merely the presence of a file. Mocked renderer
fault injection remains part of the offline suite and is labeled separately.

During the 2026-10-12 audit on Windows, the bundled Quarto CLI was discovered and the
public integration was attempted. `processx` failed to create its write pipe
(OS error 5); direct Quarto also failed to spawn Rscript (invalid handle).
These are recorded environmental integration failures, not successful renders.
See the report for any separately executed Pandoc/TeX conversions.

## Mutations and CmdStan

```powershell
New-Item -ItemType Directory -Path "$audit/mutations" -Force | Out-Null
# List selectable cases without running tests:
& $rscript --vanilla dev/validate-mutations.R --list
# Focused campaign for the two frequentist longitudinal functions:
& $rscript --vanilla dev/validate-mutations.R "$audit/mutations" 'arm-sign,missing-BH-family,truncate-scientific-CSV,false-render-success,parser-ragged,ignore-custom-chronology'
```

Run from the repository root into a fresh external directory. The script
copies the current working tree, including local/untracked files, excluding
top-level `.git`. It first runs an **unmutated copy** with the union of the
selected test-file filters and requires the designated relevant test labels
to exist. Failure, error, warning, or a missing designated test blocks the
campaign before any mutant is run. Each mutant then runs in a separate R
process. A silent marker outside its copy records execution of the altered
expression. Detection requires both that marker and a failed assertion in a
predeclared relevant test block. Failure CSVs retain the actual assertion
messages and must still be reviewed for their scientific meaning.

The summary distinguishes `detected`, `survived`, `non-reached`, `invalid`
(absent/non-unique mutation target or malformed mutated source), and `blocked`
(baseline, process/loading problem, or failures not establishing the intended
detection). Parse/load errors never count as detection. `invalid` or `blocked`
produces a nonzero exit status; survival/non-reach is a reported scientific
limit, not a false detection. Omitting the selection runs the same six cases
shown above. Historical cases remain explicitly selectable, including a Stan
case that is outside the current two-function campaign.

The chronology marker is inside the branch whose behavior changes. Historical
public calls did not exercise this internal fallback, so an unchanged test
result with no marker is classified `non-reached`, not equivalent. Public
custom chronology has separate tests. No arbitrary mutation score is a
completion criterion. Do not modify files while the copied baseline or any
mutation process is reading them; integrate all changes first.

The existing CmdStan smoke test covers all three supported likelihoods with
P=0 and P=1. It requires a configured existing CmdStan plus compiler and uses
short smoke chains, not scientific posterior validation. Run it explicitly:

```powershell
$env:MIRA_RUN_CMDSTAN_TESTS = 'true'
& $rscript --vanilla -e 'testthat::test_local(filter="dynamic-covariates_long", stop_on_failure=TRUE)'
```

CmdStan sampling was not run in this audit: the local directory was
inaccessible and no configured usable CmdStan was available. No installation
or permissions change was attempted for it.

## CI and interpretation

The workflow now configures release on Linux/Windows/macOS, Linux oldrel-1 and
Linux R 4.2.3 with compatible Matrix 1.6-5 and Deriv 4.2.0 pins. Its explicit
Depends/Imports/LinkingTo/Suggests dependency selection preserves MIRA's full
test dependencies without following Matrix's unrelated Enhances. A distinct Linux release job
explicitly installs Quarto and exercises installed-package HTML/DOCX. The
ordinary test jobs disable external integrations. Existing CmdStan workflow
remains separate. These CI jobs were inspected/configured, not dispatched.
The 2026-10-12 audit executed only local Windows R 4.4.3; R 4.5.3 lacked Imports.
Use the dated reports for the profiles actually executed in each later audit.

JOSS requires evidence beyond test outcomes. See `JOSS_READINESS.md` and
`AI_USAGE.md`; public history, research use, authorship and human review remain
the maintainer's responsibility. Neither passing tests nor these advanced
validation studies establish editorial eligibility.
