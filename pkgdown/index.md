# MIRA

**Bayesian Multilevel Inference for Longitudinal Data**

[<img src="https://img.shields.io/github/actions/workflow/status/davide-celani/MIRA/R-CMD-check.yaml?branch=main&amp;event=push&amp;style=flat&amp;label=R-CMD-check&amp;labelColor=30435A" alt="R-CMD-check" width="154" height="24">](https://github.com/davide-celani/MIRA/actions/workflows/R-CMD-check.yaml?query=branch%3Amain+event%3Apush)
[<img src="https://img.shields.io/github/actions/workflow/status/davide-celani/MIRA/cmdstan-smoke.yaml?branch=main&amp;event=push&amp;style=flat&amp;label=CmdStan%20smoke&amp;labelColor=30435A" alt="CmdStan smoke" width="183" height="24">](https://github.com/davide-celani/MIRA/actions/workflows/cmdstan-smoke.yaml?query=branch%3Amain+event%3Apush)
[<img src="https://img.shields.io/github/v/release/davide-celani/MIRA?style=flat&amp;label=GitHub%20release&amp;labelColor=30435A&amp;color=415979" alt="GitHub release" width="164" height="24">](https://github.com/davide-celani/MIRA/releases/latest)
[<img src="https://img.shields.io/badge/R-%3E%3D%204.2.0-415979?style=flat&amp;labelColor=30435A" alt="R &gt;= 4.2.0" width="94" height="24">](https://github.com/davide-celani/MIRA/blob/main/DESCRIPTION#L27)
[<img src="https://img.shields.io/badge/license-MIT-415979?style=flat&amp;labelColor=30435A" alt="License: MIT" width="94" height="24">](https://github.com/davide-celani/MIRA/blob/main/LICENSE.md)

[![Zenodo DOI for MIRA (10.5281/zenodo.23120183)](man/figures/badge-zenodo.svg)](https://doi.org/10.5281/zenodo.23120183)
[<img src="https://img.shields.io/badge/lifecycle-experimental-415979?style=flat&amp;labelColor=30435A" alt="Lifecycle: experimental" width="166" height="24">](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[<img src="https://img.shields.io/github/stars/davide-celani/MIRA?style=flat&amp;label=GitHub%20stars&amp;labelColor=30435A&amp;color=415979" alt="GitHub stars" width="123" height="24">](https://github.com/davide-celani/MIRA)

MIRA (**Multilevel Inference and Regression Analysis**) is an R package for
Bayesian modeling of continuous longitudinal outcomes. It connects repeated
measurements, treatment-arm contrasts, subject-level covariates, and
probabilities of meaningful change through an inspectable Stan model.
The package also provides frequentist exploration and Quarto report generation.

[Installation](#installation) · [Tutorials](#frequentist-longitudinal-tutorials) ·
[Reproducible example](https://github.com/davide-celani/MIRA#reproducible-example) ·
[Reference](reference/index.html) · [Changelog](news/index.html) ·
[Citation](authors.html#citation)

## Installation

MIRA requires R 4.2 or later and is installed from GitHub. Install CmdStanR
from the Stan R-universe repository before installing MIRA.

```r
install.packages(
  "cmdstanr",
  repos = c("https://stan-dev.r-universe.dev", getOption("repos"))
)
install.packages("remotes")
remotes::install_github("davide-celani/MIRA")
```

Bayesian fitting also requires a C++ toolchain and CmdStan. The
[installation guide](https://github.com/davide-celani/MIRA#installation)
covers CmdStan setup, Quarto reporting, and installation of archived releases.

## Frequentist longitudinal tutorials

The [frequentist handbook](articles/MIRA.html) maps the reading routes
and explains the ophthalmic examples. The tutorials use synthetic data to
illustrate MIRA's public API.

1. [Getting started](articles/frequentist-long-getting-started.html)
2. [Data preparation and detection](articles/frequentist-long-data-and-detection.html)
3. [Descriptive analysis](articles/frequentist-long-descriptive-analysis.html)
4. [Treatment groups](articles/frequentist-long-treatment-groups.html)
5. [Missing data](articles/frequentist-long-missing-data.html)
6. [Multiple outcomes](articles/frequentist-long-multiple-outcomes.html)
7. [A first longitudinal model](articles/frequentist-long-first-model.html)
8. [Advanced models](articles/frequentist-long-advanced-models.html)
9. [Reporting](articles/frequentist-long-reporting.html)
10. [Complete research workflow](articles/frequentist-long-complete-workflow.html)

## Citation

See [how to cite MIRA](authors.html#citation) for the citation and BibTeX entry.
The citation sources are
[`CITATION.cff`](https://github.com/davide-celani/MIRA/blob/main/CITATION.cff)
and [`inst/CITATION`](https://github.com/davide-celani/MIRA/blob/main/inst/CITATION).
