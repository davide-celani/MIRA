# Reproducing the longitudinal validation

This guide accompanies `VALIDATION_REPORT_2026-10-12.md`. Run from the repository
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

On the audited Windows host, the bundled Quarto CLI was discovered and the
public integration was attempted. `processx` failed to create its write pipe
(OS error 5); direct Quarto also failed to spawn Rscript (invalid handle).
These are recorded environmental integration failures, not successful renders.
See the report for any separately executed Pandoc/TeX conversions.

## Mutations and CmdStan

```powershell
New-Item -ItemType Directory -Path "$audit/mutations" -Force | Out-Null
& $rscript --vanilla dev/validate-mutations.R "$audit/mutations"
```

Run from the repository root into a fresh external directory. The script
copies the current working tree, including local/untracked files, excluding
top-level `.git`. Each targeted mutant runs in a separate R process. Only
relevant assertion failures count as detection; parse/load/infrastructure
errors do not. Inspect each failure CSV, not only the detection total. The
recorded chronology mutant survived because its internal fallback is not
exercised by the tested public path; public custom chronology has separate
tests. No arbitrary mutation score is a completion criterion.

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
Linux R 4.2.3 with a compatible Matrix pin. A distinct Linux release job
explicitly installs Quarto and exercises installed-package HTML/DOCX. The
ordinary test jobs disable external integrations. Existing CmdStan workflow
remains separate. These CI jobs were inspected/configured, not dispatched.
Only local Windows R 4.4.3 results were executed here; R 4.5.3 lacked Imports.

JOSS requires evidence beyond test outcomes. See `JOSS_READINESS.md` and
`AI_USAGE.md`; public history, research use, authorship and human review remain
the maintainer's responsibility. Neither passing tests nor these advanced
validation studies establish editorial eligibility.
