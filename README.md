# MIRA

<!-- badges: start -->
<!-- badges: end -->

# MIRA: Bayesian Modeling for Longitudinal Clinical Data

## Overview

**MIRA** is an R package for Bayesian analysis of longitudinal clinical data, with particular emphasis on studies characterized by:

- repeated measurements within the same subject;
- small to moderate sample sizes;
- multiple follow-up time points;
- clinically meaningful thresholds;
- individual-level and population-level inference;
- heterogeneous treatment responses.

MIRA was developed from a practical methodological problem frequently encountered in clinical research: longitudinal outcomes are often analyzed using repeated pairwise hypothesis tests, such as the Wilcoxon signed-rank test, followed by a frequentist effect size.

Although this approach can provide evidence regarding statistical differences between time points, it does not directly describe the full uncertainty surrounding the estimated effect, the probability and magnitude of clinically meaningful improvement, individual treatment trajectories, or heterogeneity between subjects.

MIRA provides a Bayesian hierarchical framework designed to extend longitudinal clinical inference beyond the question:

> "Is there a statistically significant difference?"

toward questions such as:

- How large is the estimated change?
- What is the uncertainty around that change?
- What is the posterior probability of improvement?
- What is the posterior probability that the change exceeds a clinically meaningful threshold?
- How heterogeneous are individual responses?
- How does each individual trajectory compare with the population trajectory?
- Is the observed effect clinically meaningful, rather than only statistically distinguishable from zero?

The package is intended as a practical interface for Bayesian longitudinal modeling using **Stan** and posterior simulation.

---

# Why MIRA?

A conventional longitudinal analysis based only on repeated Wilcoxon tests and effect sizes may provide a limited representation of the data.

For example, an analysis might report:

- a p-value for the comparison between baseline and follow-up;
- a second p-value for another follow-up comparison;
- an effect size such as `r`.

These quantities can be useful, but they do not directly provide a unified probabilistic description of the longitudinal process.

MIRA instead models the repeated observations jointly.

The resulting posterior distribution can be used to estimate:

1. population-level changes;
2. uncertainty around each estimated change;
3. standardized effect magnitude;
4. posterior probabilities of directional improvement;
5. posterior probabilities of clinically meaningful improvement;
6. individual-level changes;
7. subject-level response probabilities;
8. between-subject heterogeneity;
9. posterior predictive distributions;
10. model diagnostics;
11. predictive performance using approximate leave-one-out cross-validation.

The objective is not simply to replace a p-value with a Bayesian quantity. The objective is to extract substantially more clinically interpretable information from longitudinal data.

---

# Statistical framework

MIRA currently uses Bayesian hierarchical longitudinal models implemented in **Stan**.

A simplified version of the model can be represented as:

\[
y_{ij} \sim \text{Student-t}(\nu, \mu_{ij}, \sigma)
\]

where:

- \(y_{ij}\) is the observed outcome for subject \(i\) at measurement \(j\);
- \(\mu_{ij}\) is the expected value;
- \(\sigma\) represents residual variability;
- \(\nu\) controls the degrees of freedom of the Student-t distribution.

The use of a Student-t likelihood allows increased robustness to potentially influential observations compared with a purely Gaussian likelihood.

The longitudinal structure is represented using population-level effects together with subject-specific random effects.

A simplified formulation is:

\[
\mu_{ij} =
\alpha +
\beta_{\text{population}} \cdot t_j +
u_{0i} +
u_{1i} \cdot t_j
\]

where:

- \(\alpha\) is the population-level intercept;
- \(\beta_{\text{population}}\) is the population-level longitudinal slope;
- \(u_{0i}\) is the subject-specific random intercept;
- \(u_{1i}\) is the subject-specific random slope;
- \(t_j\) represents the measurement time.

This hierarchical structure allows MIRA to simultaneously estimate the overall population trajectory and individual deviations from that trajectory.

---

# Population-level inference

MIRA estimates posterior distributions for clinically relevant longitudinal contrasts.

Depending on the study design, these may include quantities such as:

- change from baseline to T1;
- change from baseline to T2;
- change from baseline to later follow-up times;
- change between consecutive follow-up visits;
- population-level longitudinal slope.

For each estimated parameter, MIRA can summarize:

- posterior mean;
- posterior median;
- posterior standard deviation;
- median absolute deviation;
- lower credible interval bound;
- upper credible interval bound;
- credible interval width;
- posterior probability of positive change;
- posterior probability of negative change.

For example:

```r
summary_mira$population
```

may return estimates such as:

| Parameter | Interpretation |
|---|---|
| `change_T1` | Population-level change at T1 |
| `change_T2` | Population-level change at T2 |
| `change_T2_vs_T1` | Additional change between T1 and T2 |
| `population_slope` | Overall longitudinal population slope |
| `sigma_intercept` | Between-subject intercept heterogeneity |
| `sigma_slope` | Between-subject slope heterogeneity |
| `rho_subject` | Correlation between individual intercepts and slopes |
| `sigma` | Residual variability |
| `nu` | Student-t degrees of freedom |

The primary interpretation is based on the estimated magnitude of change and its uncertainty rather than on dichotomous significance testing alone.

---

# Clinical meaningfulness and MCID

One of the central objectives of MIRA is to distinguish between:

1. evidence of any improvement;
2. evidence of clinically meaningful improvement.

A statistically detectable change is not necessarily clinically meaningful.

MIRA therefore allows the user to define a clinically meaningful change threshold, such as a **Minimal Clinically Important Difference (MCID)**.

If the clinically meaningful threshold is denoted by \(M\), MIRA evaluates:

\[
P(\Delta > M \mid \text{data})
\]

where \(\Delta\) represents the posterior distribution of the estimated change.

This allows direct estimation of quantities such as:

```text
P(population change > 0)
```

and:

```text
P(population change >= MCID)
```

For example:

```r
summary_mira$clinical
```

can provide:

```r
# A tibble containing:
# - MCID
# - P(population change T2 > 0)
# - P(population change T2 >= MCID)
# - P(population change T2 < MCID)
# - Mean population distance from MCID
# - P(population distance from MCID >= 0)
# - Standardized population change
```

This distinction is particularly important in clinical research.

A result may show overwhelming posterior evidence of improvement while simultaneously showing essentially no probability of reaching the predefined clinically meaningful threshold.

For example:

\[
P(\Delta > 0 \mid data) \approx 1
\]

does not imply:

\[
P(\Delta \geq MCID \mid data) \approx 1
\]

MIRA explicitly separates these two concepts.

---

# Individual-level inference

Population averages may conceal substantial differences between subjects.

A treatment may produce:

- substantial improvement in some patients;
- small improvement in others;
- no meaningful improvement in a subgroup;
- highly heterogeneous trajectories.

MIRA therefore estimates subject-specific posterior changes.

For each subject, the model can summarize quantities including:

- posterior mean change;
- posterior median change;
- posterior standard deviation;
- posterior credible interval;
- credible interval width;
- posterior change relative to the MCID;
- posterior probability of any improvement;
- posterior probability of clinically meaningful improvement.

These results are available through:

```r
summary_mira$individual
```

This allows the analysis to move beyond the statement:

> "The average patient improved."

toward:

> "How did individual patients respond, and with what posterior certainty?"

---

# Responder analysis

MIRA provides a posterior probability-based approach to responder classification.

A subject may be considered a responder when the posterior probability of clinically meaningful improvement exceeds a predefined threshold.

For example:

\[
P(\Delta_i \geq MCID \mid data) > 0.95
\]

could define a high-certainty responder.

The user can evaluate different probability thresholds, for example:

```text
0.50
0.80
0.95
```

through:

```r
summary_mira$responders
```

This approach avoids classifying subjects solely according to a single observed difference.

Instead, responder classification incorporates posterior uncertainty.

---

# Between-subject heterogeneity

MIRA explicitly models heterogeneity between individuals.

Important quantities include:

```r
summary_mira$heterogeneity
```

which may summarize:

- `sigma_intercept`;
- `sigma_slope`;
- `rho_subject`.

These parameters quantify whether subjects differ in:

1. their baseline levels;
2. their longitudinal trajectories;
3. the relationship between baseline level and longitudinal change.

For example, a larger `sigma_slope` suggests greater variability in treatment response trajectories between subjects.

This information is generally not available from a series of independent Wilcoxon tests.

---

# Posterior draws

All primary quantities are estimated as posterior distributions rather than as single deterministic values.

Posterior draws are accessible through:

```r
summary_mira$draws
```

These draws may include:

- population means at each time point;
- longitudinal contrasts;
- population slope;
- residual variability;
- random-effect parameters;
- Student-t degrees of freedom;
- population probabilities of improvement;
- population probabilities of meaningful improvement;
- individual-level posterior changes;
- individual-level responder indicators;
- pointwise log-likelihood values;
- posterior predictive replicated observations.

Access to posterior draws allows users to perform additional analyses beyond the standard MIRA summaries.

For example, users may calculate:

- custom posterior probabilities;
- alternative clinical thresholds;
- custom responder definitions;
- decision-oriented summaries;
- posterior risk differences;
- secondary standardized effects.

---

# Standardized population change

MIRA can quantify the magnitude of population-level change relative to model-estimated variability.

A standardized posterior change can be expressed conceptually as:

\[
d = \frac{\Delta}{\sigma}
\]

where:

- \(\Delta\) is the posterior longitudinal contrast;
- \(\sigma\) is the residual standard deviation.

Unlike a conventional single-value effect size, this quantity can itself be represented as a posterior distribution.

The resulting inference can therefore include:

- posterior mean or median standardized change;
- uncertainty around the standardized change;
- credible intervals;
- posterior directional probabilities.

This provides a probabilistic representation of effect magnitude rather than a single point estimate alone.

---

# Robust likelihood

Clinical datasets, particularly small datasets, may contain observations that have a disproportionate influence on Gaussian models.

MIRA uses a Student-t likelihood in the current hierarchical framework.

The degrees of freedom parameter:

\[
\nu
\]

is estimated from the data.

Lower values of \(\nu\) correspond to heavier tails, allowing increased robustness to potentially influential observations.

As \(\nu\) increases, the Student-t distribution increasingly resembles a Gaussian distribution.

This provides a flexible model that can adapt to departures from strict normality.

---

# Posterior predictive checks

MIRA includes posterior predictive replicated data.

These can be accessed through:

```r
summary_mira$ppc
```

and through posterior draws containing:

```text
y_rep
```

Posterior predictive checking evaluates whether observations simulated from the fitted model are compatible with important features of the observed data.

This is a central part of Bayesian workflow.

The objective is not only to estimate parameters but also to evaluate whether the fitted model can reasonably reproduce the observed data-generating structure.

---

# Predictive model assessment

MIRA provides approximate leave-one-out cross-validation information through:

```r
summary_mira$loo
```

This allows predictive assessment of the fitted model and can support comparisons between competing model specifications.

Possible future applications include comparison of:

- alternative random-effects structures;
- Gaussian versus robust likelihoods;
- models with and without specific covariates;
- different longitudinal formulations.

Predictive model comparison should complement, rather than replace, substantive and clinical model evaluation.

---

# MCMC diagnostics

Bayesian inference requires careful assessment of computational convergence.

MIRA provides diagnostic information through:

```r
summary_mira$diagnostics
```

Important diagnostics may include:

- R-hat;
- effective sample size;
- Monte Carlo error;
- divergent transitions;
- other sampling quality indicators.

A fitted model should not be interpreted solely because it returns numerical posterior estimates.

Reliable Bayesian inference requires adequate sampling diagnostics and convergence assessment.

---

# Quality flags

MIRA includes a dedicated output:

```r
summary_mira$quality_flags
```

designed to identify potential issues requiring attention before substantive interpretation.

Quality checks may include warnings related to:

- convergence;
- sampling diagnostics;
- effective sample size;
- posterior predictive performance;
- influential observations;
- approximate LOO diagnostics.

The purpose of this component is to encourage a workflow in which model adequacy is evaluated before final conclusions are drawn.

---

# Model information

General information about the fitted model is available through:

```r
summary_mira$model_information
```

This can be used to document the model specification and important analytical settings.

Such information is particularly useful for:

- reproducibility;
- reporting;
- supplementary materials;
- comparison of different analyses.

---

# Basic workflow

A typical MIRA workflow is:

```r
library(MIRA)
```

Prepare the longitudinal dataset:

```r
stan_data <- mira_prepare(
  data = data,
  id = "subject",
  time = "time",
  outcome = "outcome",
  meaningful_change = 5
)
```

Fit the Bayesian longitudinal model:

```r
fit <- mira_fit(
  data = stan_data
)
```

Create the comprehensive model summary:

```r
summary_mira <- mira_summary(
  fit = fit,
  meaningful_change = stan_data$meaningful_change,
  y = stan_data$y
)
```

Inspect population-level inference:

```r
summary_mira$population
```

Inspect clinical interpretation:

```r
summary_mira$clinical
```

Inspect individual trajectories:

```r
summary_mira$individual
```

Inspect responder probabilities:

```r
summary_mira$responders
```

Inspect between-subject heterogeneity:

```r
summary_mira$heterogeneity
```

Inspect posterior draws:

```r
summary_mira$draws
```

Inspect posterior predictive checks:

```r
summary_mira$ppc
```

Inspect predictive model assessment:

```r
summary_mira$loo
```

Inspect MCMC diagnostics:

```r
summary_mira$diagnostics
```

Inspect quality flags:

```r
summary_mira$quality_flags
```

Inspect model information:

```r
summary_mira$model_information
```

---

# Example of interpretation

Suppose a model estimates:

```text
Population change at T2 = 1.95
95% CrI = 1.77 to 2.13

P(change > 0) = 1.00

MCID = 5

P(change >= MCID) = 0.00
```

These results support two distinct conclusions.

First:

> There is overwhelming posterior evidence of improvement.

Second:

> The estimated improvement is substantially below the predefined threshold for clinical meaningfulness.

The appropriate interpretation is therefore not simply:

> "The treatment worked."

A more informative interpretation is:

> "The model estimated a highly probable improvement, but the magnitude of the improvement did not reach the predefined threshold for clinical meaningfulness."

This distinction represents one of the main motivations for MIRA.

---

# From statistical significance to clinical interpretation

MIRA is designed around the principle that statistical evidence and clinical importance are different concepts.

A conventional analysis may answer:

> Is the observed difference statistically distinguishable from zero?

MIRA additionally asks:

> What is the estimated magnitude of the change?

> How uncertain is this estimate?

> What is the probability that the change is beneficial?

> What is the probability that the change is clinically meaningful?

> How much do individual responses vary?

> Does the population average conceal clinically relevant heterogeneity?

The package is therefore intended to support richer interpretation of longitudinal clinical data.

---

# Current model capabilities

The current MIRA framework includes:

- Bayesian longitudinal hierarchical modeling;
- repeated measurements within subjects;
- population-level temporal effects;
- subject-specific random intercepts;
- subject-specific random slopes;
- correlation between random intercepts and slopes;
- robust Student-t likelihood;
- posterior longitudinal contrasts;
- posterior credible intervals;
- posterior directional probabilities;
- clinically meaningful change thresholds;
- MCID-based inference;
- population-level clinical interpretation;
- individual-level posterior estimates;
- posterior responder probabilities;
- heterogeneity assessment;
- standardized posterior effect estimates;
- posterior predictive replicated data;
- posterior predictive checks;
- pointwise log-likelihood extraction;
- approximate leave-one-out cross-validation;
- MCMC diagnostics;
- automated quality flags;
- comprehensive model summaries.

---

# Planned development

MIRA is under active development.

Future development is expected to focus on increasing model flexibility while maintaining a practical workflow for clinical researchers.

Potential areas of development include:

## Flexible number of time points

Support for longitudinal studies with different numbers of assessments, including:

- 2 time points;
- 3 time points;
- 4 time points;
- 5 time points;
- 6 or more time points.

The aim is to allow the model structure to adapt to the available longitudinal design rather than requiring a separate manually rewritten model for each number of assessments.

## Covariate extension

Future versions may allow incorporation of additional covariates, such as:

- age;
- sex;
- treatment;
- baseline clinical characteristics;
- other user-specified predictors.

These covariates may be incorporated into both population-level and longitudinal components of the model where appropriate.

## Personalized priors

Future development may allow greater control over prior specification.

Potential functionality may include:

- default weakly informative priors;
- user-defined priors;
- prior sensitivity analysis;
- standardized prior templates.

The goal is to allow users to incorporate justified prior information while maintaining transparent and reproducible model specifications.

## Expanded model comparison

Future versions may support structured comparison between alternative model formulations using:

- predictive performance;
- posterior predictive checks;
- LOO diagnostics;
- sensitivity analyses.

## Extended clinical decision metrics

Potential future development includes additional clinically oriented posterior summaries and responder definitions.

---

# Intended use

MIRA is particularly relevant for clinical and biomedical studies involving repeated measurements.

Potential applications include:

- ophthalmology;
- longitudinal treatment studies;
- pilot studies;
- rare diseases;
- retrospective cohorts;
- small clinical samples;
- early-phase clinical research;
- studies with heterogeneous individual responses.

The framework may be particularly useful when conventional statistical workflows are limited to repeated pairwise tests and single-number effect sizes.

MIRA does not assume that Bayesian analysis automatically makes a study more reliable.

The quality of inference still depends on:

- study design;
- data quality;
- sample size;
- measurement quality;
- model assumptions;
- prior specification;
- model diagnostics;
- posterior predictive adequacy.

The purpose of MIRA is to provide a richer analytical framework for these data, not to compensate for poor study design.

---

# Reproducibility

MIRA is built around reproducible Bayesian analysis.

Users should report:

- the R version;
- the MIRA version;
- the Stan/CmdStan configuration where applicable;
- model specification;
- prior distributions;
- number of chains;
- warm-up iterations;
- sampling iterations;
- convergence diagnostics;
- posterior predictive checks;
- clinical thresholds used;
- sensitivity analyses where relevant.

Clinical conclusions should be based on the complete posterior evidence rather than on a single threshold alone.

---

# Interpretation philosophy

MIRA emphasizes:

```text
Effect magnitude
        +
Uncertainty
        +
Posterior probability
        +
Clinical meaningfulness
        +
Individual heterogeneity
        +
Model adequacy
```

rather than relying exclusively on:

```text
p-value
        +
single effect size
```

The central question is therefore not only:

> "Is there evidence of a difference?"

but also:

> "How large is the change, how certain are we, is it clinically meaningful, and how consistently does it occur across individuals?"

---

# Limitations

MIRA is a statistical modeling framework and should not be interpreted as providing causal evidence by itself.

Posterior probabilities describe uncertainty conditional on:

- the observed data;
- the specified likelihood;
- the model structure;
- the prior distributions.

Results remain dependent on model assumptions and study design.

Small samples can benefit from hierarchical modeling and regularization, but Bayesian methods do not create information that is absent from the data.

Users should perform appropriate diagnostic and sensitivity analyses before drawing substantive conclusions.

---

# Citation

If you use MIRA in academic work, please cite the package according to the citation information provided by:

```r
citation("MIRA")
```

Additional citation information will be added as the package and associated methodological work develop.

---

# Contributing

Contributions, methodological discussion, bug reports, and feature suggestions are welcome.

Potential areas for contribution include:

- Bayesian model development;
- Stan optimization;
- diagnostic tools;
- visualization;
- simulation studies;
- documentation;
- clinical validation;
- reproducibility testing.

---

# License

License information will be provided with the package.

---

# Acknowledgments

MIRA was developed to address a practical problem in longitudinal clinical research: extracting clinically interpretable information from repeated measurements while explicitly representing uncertainty, individual heterogeneity, and the distinction between statistical change and clinically meaningful change.

The package is based on Bayesian hierarchical modeling and is implemented using R and Stan.
