# BayesMFRM

<!-- badges: start -->
[![Experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

## ⚠️ Experimental Early Version

**This is an experimental early version of BayesMFRM. The API may change, and some features may not be fully stable. Use with caution in production environments.**

## Overview

BayesMFRM is an R package for fitting Bayesian Many-Facet Rasch Models (MFRM) using Stan. It provides a user-friendly interface similar to brms but specialized for MFRM analysis, producing FACETS-style outputs including rater severity, item difficulty, fair averages, and infit/outfit statistics.

### Key Features

- 📊 **Simple formula interface** for specifying models (`score ~ person + item + rater`)
- 🔢 **RSM and PCM model families** (Rating Scale Model and Partial Credit Model)
- 🎯 **FACETS-style diagnostics** (infit/outfit, bias analysis, fair scores)
- 🔄 **Smart caching system** for fast model re-runs
- 📈 **Bayesian uncertainty quantification** for all parameters
- 🎨 **Flexible prior specification** system
- 📉 **Posterior predictive checks** and residual analysis
- ⚡ **Dynamic Stan code generation** with templates
- 🏗️ **Robust architecture** with comprehensive error checking

## Installation

### Prerequisites

This package requires CmdStan. Install it with:

```r
install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
library(cmdstanr)
install_cmdstan()
```

### Install BayesMFRM

**Using remotes:**
```r
# install.packages("remotes")
remotes::install_github("egumasa/BayesMFRM")
```

**Using devtools:**
```r
# install.packages("devtools") 
devtools::install_github("egumasa/BayesMFRM")
```

## Quick Start

```r
library(BayesMFRM)

# Load or create your rating data
# Data should have columns for person, item, rater, and score (1 to K)
data(ratings)  # example data

# Check data suitability
check_data_suitability(ratings, "score", c("person", "item", "rater"))

# Fit basic 3-facet model
fit1 <- bmfrm(score ~ person + item + rater, data = ratings, K = 6)

# View results
print(fit1)
summary(fit1)

# Model checking
pp_check(fit1)

# Residual analysis  
resid_tbl <- residuals(fit1)
rater_fit <- facet_fit(fit1, "rater")
item_fit <- facet_fit(fit1, "item")

# Fair scores (removing rater effects)
fair_scores_summary <- fair_scores(fit1)
```

### Model Families

```r
# Rating Scale Model (default) - shared thresholds across items
fit_rsm <- bmfrm(score ~ person + item + rater, data = ratings, K = 6, family = "rsm")

# Partial Credit Model - item-specific thresholds
fit_pcm <- bmfrm(score ~ person + item + rater, data = ratings, K = 6,
                 family = "pcm", threshold_facet = "item")
```

### Advanced Usage

```r
# Model with bias terms
fit2 <- bmfrm(
  score ~ person + item + rater + rater:item,  # rater×item bias
  data = ratings
)

# Custom priors
priors <- c(
  prior("normal(0, 2)", class = "theta"),                    # person abilities
  prior("normal(0, 1)", class = "item"),                     # item difficulties  
  prior("student_t(3, 0, 1)", class = "rater"),             # rater severities
  prior("normal(0, 0.5)", class = "bias", facet = "rater:item"), # bias effects
  prior("normal(0, 3)", class = "tau")                      # thresholds
)

fit3 <- bmfrm(
  score ~ person + item + rater + rater:item,
  data = ratings,
  priors = priors,
  chains = 4,
  iter = 2000
)

# Bias analysis
bias_summary <- summarise_bias(fit3)
print(bias_summary)

# Cross-criteria correlations
crit_cor <- residual_cor_criteria(resid_tbl)
print(crit_cor)
```

## Data Requirements

Your data must meet these requirements:

- **Score column**: Integer values from 1 to K (not 0-indexed)
- **Facet columns**: Any values convertible to factors (strings, numbers)
- **No missing values**: 
  - NA scores will be dropped with a message
  - NA in facets will cause an error (clean your data first)

```r
# Good data format:
ratings$score   # 1, 2, 3, 4, 5, 6
ratings$person  # "P001", "P002", ... or 1, 2, 3, ...
ratings$rater   # "R01", "R02", ... or 101, 102, 103, ...

# Will cause errors:
ratings$score   # 0, 1, 2, 3, 4, 5  (0-indexed - add 1 first)
ratings$rater   # "R01", NA, "R03"  (NA values not allowed)
```

## Current Status

**Version 0.1.0 (Experimental)**

### ✅ Implemented Features

- **Model Families**: Rating Scale Model (RSM) and Partial Credit Model (PCM)
- **Facets**: Main facet effects (person, item, rater, etc.)
- **Bias Terms**: 2-way bias facets (e.g., rater × item)
- **Interface**: Complete S3 method system
- **Diagnostics**: MCMC convergence checking, infit/outfit statistics
- **Analysis**: Fair score extraction, bias analysis, residual correlations
- **Caching**: Hash-based model and fit caching
- **Validation**: Comprehensive data and model validation

### 🚧 Planned Features

- Graded response models
- Multidimensional abilities
- Centrality/drift modeling
- More complex bias structures (rater×task, etc.)

## System Requirements

- **R**: >= 4.0.0
- **CmdStan**: Installed via cmdstanr  
- **C++ toolchain**: Required for Stan compilation
- **Memory**: Sufficient RAM for MCMC sampling (depends on model size)

## Dependencies

Core dependencies automatically installed:
- cmdstanr, posterior, dplyr, tidyr, tibble, rlang, glue, digest, stringr, purrr, magrittr, bayesplot

## Getting Help

- **Vignettes**: See `vignette/spec_v0.4.md` for detailed specification
- **Examples**: Check `example_usage.R` for workflow examples  
- **Issues**: Report bugs and feature requests on GitHub
- **Testing**: Run `test_bmfrm.R` to validate installation

## Citation

If you use BayesMFRM in your research, please cite:

```
[Citation information to be added when published]
```

## Contributing

This is an experimental package under active development. Contributions are welcome:

- **Bug reports**: Use GitHub issues
- **Feature requests**: Describe your use case
- **Code contributions**: Follow existing code style
- **Documentation**: Improvements always appreciated

## License

This work is licensed under a [Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License (CC BY-NC-SA 4.0)](https://creativecommons.org/licenses/by-nc-sa/4.0/).

---

**Note**: This package is in active development. The API and functionality may change in future versions. We recommend pinning to specific versions for reproducible research.
