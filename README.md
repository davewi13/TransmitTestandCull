# TransmitTestandCull

An R package for simulating stochastic infectious disease transmission in animal populations with periodic diagnostic testing and removal of positive animals.
The package was developed to support investigations of ovine pulmonary adenocarcinoma (OPA) control strategies but may also be applicable to other infectious diseases with latent and infectious stages, imperfect diagnostics, and test-and-cull interventions.

## Overview

`TransmitTestandCull` implements a stochastic individual-based transmission model with four disease states:

- Susceptible (S)
- Exposed / Latent (E)
- Infectious (I)
- Removed (R)

Transmission occurs from infectious individuals to susceptible individuals. Disease progression, transmission events, diagnostic outcomes, and removals are all modelled stochastically.

The package:

- Runs stochastic transmission simulations
- Incorporates diagnostic testing at specified time points
- Models removal of test-positive animals
- Incorporates periodic restocking
- Processes simulation outputs
- Visualises disease prevalence trajectories

## Model Assumptions

The current implementation assumes:

- Homogeneous mixing within the population.
- Frequency-dependent disease transmission.
- Individuals progress through latent and infectious stages before natural removal (removal due to testing can happen at any point in the disease progression).
- Latent and infectious periods follow Gamma distributions.
- Diagnostic sensitivity may differ between latent and infectious animals.
- Diagnostic specificity is constant for all uninfected animals.
- Positive test results lead to immediate removal from the flock.
- Restocking occurs through the addition of susceptible animals.

## Installation

### Install from GitHub

```r
remotes::install_github("davewi13/TransmitTestandCull")
```

### Install from a local source package

```r
install.packages(
  "TransmitTestandCull_1.0.0.tar.gz",
  repos = NULL,
  type = "source"
)
```

### Install from source directory

```bash
R CMD INSTALL .
```

## Main Functions

### `simulate_OPA()`

Runs the stochastic transmission model and returns individual-level infection histories.

#### Inputs

- Population size
- Initial prevalence
- Basic reproduction number (R0)
- Mean durations of latent and infectious periods
- Shape parameters for latent and infectious periods
- Diagnostic test sensitivity (separately for latent and infectious periods) and specificity
- Scanning schedule
- Restocking frequency

#### Output

A data frame containing infection histories for all individuals, including:

- Exposure time
- Infection time
- Removal time
- Indicator for disease-negative animals incorrectly culled

---

### `process_results()`

Takes the data frame of exposure, infection and removal times from simulate_OPA() and reformats it to give a single column for time with counts of exposures, infections and removals

---

### `plot_trajectory()`

Generates prevalence trajectories from one or more simulation runs.

The resulting plot includes:

- Individual simulation trajectories
- Median prevalence through time
- 95% prediction intervals

### `extend_data()`

Extends simulation outputs to create complete trajectories in the case where eradication occurred before the simulation would otherwise have ended.

---

### `run_parameter_combo()`

Explores the dynamics when two parameters are varied separately. Be very careful when using this! It is written quite specifically for the purposes of the OPA paper and doesn't necessarily generalise well depending on whether parameters are continuous or only take specific values.

---

## Example Workflow

```r
library(TransmitTestandCull)

params <- list(
  N0 = 1000,
  initial_prevalence = 0.03,
  sensitivity_infectious = 0.76,
  sensitivity_latent = 0.38,
  specificity = 0.998,
  prop_infectious = 0.33,
  scanning_timepoints = c(0, 365, 730, 1095),
  time_between_restocking = 365,
  R0 = 1.39,
  mu_L = 180,
  sh_L = 30,
  mu_I = 90,
  sh_I = 10
)

simulated_infections <- simulate_OPA(
  params = params,
  Tmax = 10 * 365
)

processed_results <- process_results(
  infections = simulated_infections,
  N0 = params$N0,
  time_between_restocking = params$time_between_restocking
)

plot_trajectory(
  temporal_infections = processed_results,
  Tmax = 10 * 365,
  time_between_scans = params$scanning_timepoints
)
```

## Applications

Potential applications include:

- Evaluating test-and-cull control strategies
- Assessing diagnostic test performance
- Estimating disease eradication probabilities
- Sensitivity analyses of transmission parameters
- Exploring long-term disease dynamics under different management scenarios

## Repository Structure

```text
TransmitTestandCull/
├── R/
├── man/
├── src/
├── analysis/
├── DESCRIPTION
├── NAMESPACE
└── README.md
```

The computationally intensive simulation code is implemented in C++ using Rcpp, while data processing and visualisation functions are implemented in R.
