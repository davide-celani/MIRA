# MIRA

**Bayesian Multilevel Inference for Longitudinal Data**

MIRA (**Multilevel Inference and Regression Analysis**) is an R package for
Bayesian modeling of continuous longitudinal outcomes. It connects repeated
measurements, treatment-arm contrasts, subject-level covariates, and
probabilities of meaningful change through an inspectable Stan model.
The package also provides frequentist exploration and Quarto report generation.

Version `0.0.1` is being prepared for its first release. MIRA is under active
development; no Zenodo DOI has been assigned yet.

## What MIRA does

- Fits longitudinal multilevel models with correlated subject-specific
  random intercepts and slopes, a common time trajectory, and treatment-
  and covariate-specific trajectories.
- Supports two or more measurement occasions and two or more treatment arms.
- Offers Student-t and Gaussian identity-link models and a log-normal
  log-link model, with optional observation bounds handled by censoring.
- Encodes numeric, categorical, and logical subject-level covariates and
  records their scaling and reference levels.
- Exposes configurable priors, posterior contrasts, uncertainty intervals,
  directional probabilities, and clinically meaningful change thresholds.
- Summarizes existing-subject responses and new-subject latent and predictive
  responder probabilities separately.
- Provides posterior predictive checks, pointwise log-likelihood values,
  and MCMC diagnostics.
- Explores longitudinal datasets through `mira_info_long()` and writes
  frequentist report bundles through `mira_report_freq_long()`.

These are implemented features. See [NEWS.md](NEWS.md) for the release inventory.

## Installation

MIRA requires R 4.2 or later. It is installed from GitHub and is not currently
distributed through CRAN. Install CmdStanR from the Stan R-universe repository
before installing MIRA:

```r
install.packages(
  "cmdstanr",
  repos = c("https://stan-dev.r-universe.dev", getOption("repos"))
)
install.packages("remotes")
remotes::install_github("davide-celani/MIRA")
```

The development branch may change. After a tagged release is published,
install that release with `remotes::install_github("davide-celani/MIRA@v0.0.1")`.
This tag is a planned release identifier until it is published.

Bayesian fitting also requires a C++ toolchain and CmdStan. On Windows, install
the Rtools version that matches your R installation. Then run the following
one-time setup from R:

```r
cmdstanr::check_cmdstan_toolchain()
cmdstanr::install_cmdstan()
cmdstanr::cmdstan_version()
```

CmdStan is separate from the `cmdstanr` R package. Installing MIRA alone does
not install CmdStan. See the
[CmdStanR setup instructions](https://mc-stan.org/cmdstanr/articles/cmdstanr.html)
for operating-system requirements and existing installations.

HTML/PDF report rendering additionally requires the
[Quarto CLI](https://quarto.org/docs/get-started/). PDF rendering requires a
working TeX installation. The R package `quarto` alone does not install the
CLI. Report source files can be generated with `render = FALSE`.

## Reproducible example

This example uses simulated observations, not participant data. It runs data
preparation without a CmdStan installation:

```r
library(MIRA)

set.seed(20261003)
n <- 24L
arm <- rep(c("Control", "Treatment"), each = n / 2L)
baseline <- round(rnorm(n, mean = 60, sd = 7))
wide_data <- data.frame(
  patient = sprintf("P%02d", seq_len(n)),
  arm = arm,
  age = sample(45:80, n, replace = TRUE),
  BCVA_t0 = baseline,
  BCVA_t1 = pmin(100, pmax(0, round(
    baseline + 1 + 2 * (arm == "Treatment") + rnorm(n, 0, 2)
  ))),
  BCVA_t2 = pmin(100, pmax(0, round(
    baseline + 2 + 4 * (arm == "Treatment") + rnorm(n, 0, 2)
  )))
)

stan_data <- mira_data_long(
  data = wide_data,
  time_value = c(0, 3, 6),
  outcome = "BCVA",
  likelihood = "student_t",
  direction = "higher",
  meaningful_change = 5,
  meaningful_change_sd = 1,
  reference_arm = "Control",
  covariates = "age"
)

prior <- mira_prior_long(stan_data)
print(prior)
```

With CmdStan configured, fit the model and summarize the posterior:

```r
fit <- mira_fit_long(
  stan_data = stan_data,
  prior = prior,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 2000,
  iter_sampling = 3000,
  seed = 20261003
)

summary_mira <- mira_summary_long(
  fit = fit,
  stan_data = stan_data,
  credible_level = 0.95
)

summary_mira$change_from_baseline
summary_mira$treatment$from_baseline
summary_mira$covariates$from_baseline
summary_mira$individual_change
summary_mira$responders
summary_mira$population_clinical
summary_mira$diagnostics
summary_mira$quality_flags
```

The summary defaults to 90% credible intervals; the example explicitly requests
95%. Inspect convergence, divergences, effective sample sizes, and posterior
predictive adequacy before interpreting results. Adequate iteration counts and
prior sensitivity checks depend on the analysis.

For frequentist exploration of the same data:

```r
exploration <- mira_info_long(
  data = wide_data,
  id = "patient",
  time_vars = c("BCVA_t0", "BCVA_t1", "BCVA_t2"),
  arm = "arm",
  reference_arm = "Control",
  plots = FALSE,
  model = FALSE,
  verbose = FALSE
)

report <- mira_report_freq_long(
  x = exploration,
  output_dir = file.path(tempdir(), "mira-example-report"),
  render = FALSE,
  open = FALSE
)
```

## Data and model assumptions

The Bayesian interface expects wide data with one row and one unique `patient`
identifier per subject. Measurement columns use `<outcome>_t0`,
`<outcome>_t1`, and subsequent integer suffixes. Actual time values must be
finite and strictly increasing. Select one outcome per Bayesian fit and
explicitly choose the reference arm and covariates.

`likelihood = "auto"` selects Student-t for BCVA, log-normal for CMT, and a
Student-t fallback with a warning for other continuous outcomes. Check the
measurement scale and use an explicit family when appropriate. BCVA defaults
to bounds of 0 and 100; confirm these match the instrument used. Log-normal
models require strictly positive measurements.

The MCID has a Gamma prior informed by `meaningful_change` and
`meaningful_change_sd`, both in natural outcome units. Specify these from
external clinical justification. The present model does not learn clinical
importance from the outcome observations. For a fixed-threshold sensitivity
analysis, review the summary API and the underlying posterior draws carefully.

## Interpretation and current limits

Posterior probabilities are conditional on the data, likelihood, model
structure, and priors. Treatment contrasts alone do not establish a causal
effect. Small samples still require appropriate design and sensitivity checks.

The bundled Bayesian model has one subject grouping level. It does not
implement eyes nested within patients, an explicit crossover/carryover model,
or time-varying covariates. Modeling these designs requires additional model
development and validation. A fixed `study_eye` covariate does not create a
patient/eye random-effects hierarchy.

For bounded outcomes, quantities called `population_mean` in the Stan output
use a bounded latent-location transformation for identity-link models and a
bounded uncensored-mean transformation for log-normal models. They are not
exact marginal expectations of the censored observed distribution. Related
change and contrast summaries inherit this interpretation. Predictive draws
include observation noise and censoring. This distinction needs to be retained
in reporting and evaluated in future statistical validation.

The Student-t scale parameter is not its residual standard deviation. MIRA's
generated `residual_sd` uses the degrees-of-freedom correction. On the log-normal
branch, `sigma` is on the log scale; do not interpret it as an outcome-unit SD.

The current summary exposes pointwise log-likelihood information through
`$loo`; this does not automatically perform a complete PSIS-LOO analysis.
Choose a validation unit consistent with the scientific prediction task,
especially when observations are clustered within subjects.

## Documentation and tests

Function documentation is available in R, for example:

```r
?mira_data_long
?mira_prior_long
?mira_fit_long
?mira_summary_long
?mira_info_long
?mira_report_freq_long
```

From the repository root, run:

```r
devtools::test()
devtools::check(args = "--no-manual")
```

The ordinary suite tests data preparation, covariate encoding, exploration,
statistical reference calculations, summaries, and report bundles. Actual
CmdStan compilation and sampling are opt-in:

```r
Sys.setenv(MIRA_RUN_CMDSTAN_TESTS = "true")
devtools::test(filter = "dynamic-covariates_long")
```

The CmdStan smoke test checks software execution; its deliberately short chains
do not establish inferential accuracy or convergence. Quarto integration tests
require the CLI and are skipped when it is unavailable.

## Citation

The software citation title is **MIRA: Bayesian Multilevel Inference for
Longitudinal Data**. The author is Davide Celani. Retrieve the R citation with:

```r
citation("MIRA")
toBibtex(citation("MIRA"))
```

GitHub uses [CITATION.cff](CITATION.cff); R uses [inst/CITATION](inst/CITATION).
For an analysis that uses a specific archived release, cite its version DOI
once Zenodo assigns it. The general concept DOI identifies the project across
versions. Until the archive exists, the citation points to the source repository.

## Contributing and support

Open [an issue](https://github.com/davide-celani/MIRA/issues) for bugs or feature
requests. Read [CONTRIBUTING.md](CONTRIBUTING.md) for development and testing
guidance and [CODE_OF_CONDUCT.md](CODE_OF_CONDUCT.md) for community expectations.
MIRA is maintained by Davide Celani; response times depend on availability.

The [release guide](dev/RELEASE.md) describes preparation of a GitHub release
and Zenodo archive. Future JOSS work is tracked in
[dev/JOSS_READINESS.md](dev/JOSS_READINESS.md). AI assistance during this
preparation is recorded in [dev/AI_USAGE.md](dev/AI_USAGE.md).

## License

MIRA is distributed under the [MIT license](LICENSE.md).
Copyright 2026 Davide Celani.
