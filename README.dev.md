# BayesMFRM Developer & Verification Guide (v0.1.0.9008)

The repository root is the **single canonical Git working tree** for package development, testing, simulation benchmarks, and `R CMD check`. Release history is maintained in [NEWS.md](NEWS.md), and user-facing documentation lives in [README.md](README.md), `vignette("getting-started", package = "BayesMFRM")` ([vignettes/getting-started.Rmd](vignettes/getting-started.Rmd)), and [docs/model-comparison-ratings1.md](docs/model-comparison-ratings1.md).

Use English for figure titles, axes, legends, annotations, tooltips, and captions in the package, examples, and interactive guide.

## Verification & Test Suite Commands

```sh
# 1. Fast unit & regression suite (no MCMC sampling; runs in ~5s)
Rscript --vanilla tools/test-local.R

# 2. Fixed-point Stan log-density & analytic HMC gradient parity checks against independent R math
Rscript --vanilla tools/check-rsm-stan.R
Rscript --vanilla tools/check-bias-stan.R
Rscript --vanilla tools/check-extensions-stan.R

# 3. Full integration test suite including live CmdStan MCMC runs
BAYESMFRM_RUN_STAN=true Rscript --vanilla tools/test-local.R

# 4. Clean CRAN-style build & check in tempdir() (keeps working tree clean)
Rscript --vanilla -e '
  td <- tempfile("bmfrm_check_"); dir.create(td)
  pkg <- system2("R", c("--vanilla", "CMD", "build", "."), stdout = TRUE)
  tarball <- list.files(".", "^BayesMFRM_.*\\.tar\\.gz$", full.names = TRUE)[1]
  file.rename(tarball, file.path(td, basename(tarball)))
  status <- system2("R", c("--vanilla", "CMD", "check", "--no-manual",
                           paste0("--output=", td), file.path(td, basename(tarball))))
  quit(status = status)
'
```

## Simulation & Calibration Benchmarks

All benchmark scripts load the local checkout via `pkgload::load_all(".")` and compile Stan binaries inside temporary directories (`PRECOMPILED_HEADERS = "false"`) so the repository stays clean:

- **Simulation-Based Calibration (Talts et al., 2018)**:
  `Rscript --vanilla simulations/sbc_calibration.R [S=40]`
  Validates rank uniformity, empirical ECDF bands, posterior contraction, and 50%/90% coverage across `Hier_RSM`, `Hier_PCM`, `Bias_StudentT`, and `Bias_Horseshoe` (`simulations/results/sbc/`).
- **Joint ML vs. Gauss-Hermite MML-EAP vs. BayesMFRM & Shrinkage Trade-Off**:
  `Rscript --vanilla simulations/bayesian_benefits.R [B=20]`
  Compares Joint ML, 31-node Gauss-Hermite MML-EAP, and BayesMFRM on examinee ability recovery, uncertainty coverage (including extreme response patterns), held-out log-loss, and null-vs-active cell RMSE in sparse interaction models (`simulations/results/`).
- **ACER ConQuest 5.47.5 Cross-Software Benchmark**:
  `Rscript --vanilla simulations/conquest_comparison.R [B=10]`
  Verifies parameter, EAP, standardized fair-score, and PV-averaged MNSQ parity across RSM, uniform PCM, heterogeneous PCM, and interaction models, plus sparse unbalanced replication benchmarks (`simulations/results/conquest/`).

The companion design study in `mfrm-replication-bayes` runs with
`Rscript --vanilla scripts/run-test-information.R` from that repository.
It reuses the archived four-rater, four-criterion interaction experiment
to compare rating counts, equal-budget rosters, and the allocation of 80
ratings across 20 examinees. Its figures, CSVs, and interpretation are in
`results/test-information/` and `reports/test-information.md`.

## Core Architecture & Design Contracts

1. **Role-Agnostic Formula & Prior Semantics (`R/formula_parsing.R`, `R/prior.R`, `R/model_spec.R`)**
   - `person_facet` (default `"person"`) identifies the unconstrained examinee ability facet (`theta`); all other main facets, step deviations (`tau`), and two-way interactions (`f1:f2`) are centered to sum to zero across declared levels (within each `threshold_facet` level for PCM).
   - `adjust_facets` identifies nuisance facets set to zero in `fair_scores()` and optional `generated quantities` (`mu_fair`); defaults to `intersect(c("rater", "interlocutor"), non_person_facets)` when `NULL`.
   - `prior()` normalizes canonical `class = "facet"` alongside legacy `"item"`, `"rater"`, and `"facet_<name>"` aliases via `prior_target_key()`, and validates distribution calls through `validate_stan_prior_spec()` (single-call AST parse + parameter domain whitelist, blocking semicolon/statement injection).

2. **Vectorized Stan Templates & Evaluation Contract (`inst/stan/*.stan`, `R/stan_generation.R`, `R/posterior_quantities.R`)**
   - Templates (`bmfrm_rsm.stan`, `bmfrm_pcm.stan`, `bmfrm_pcm_hetero.stan`) vectorize main-facet linear predictors, precompute cumulative step thresholds once per draw, and expand `// {{GQ_BLOCK}}` only when `save_gq = "full"`.
   - `resolve_eval_contract()` restricts R-side streaming recomputation to generated models and byte-verified bundled reference models (`inst/stan/models/*`), preventing silent linear-predictor mismatches on arbitrary custom `stan_file` inputs.

3. **Calibrated Fit, Multiplicity-Controlled Bias, & Portable Caching (`R/fit_stats.R`, `R/methods.R`, `R/caching.R`)**
   - `facet_fit()` streams matched realized $\text{MSQ}(y, \psi^{(d)})$ and replicated $\text{MSQ}(y_{\text{rep}}^{(d)}, \psi^{(d)})$ draws (with mid-$p$ tie handling for discrete small-$n$ cells) to report Sinharay (2006) `infit_ppp`, `outfit_ppp`, and predictive reference intervals as primary diagnostics.
   - `summarise_bias()` marks `n_observations == 0` cells as `"unobserved_cell"`, excludes both unobserved and likelihood-unidentified cells from tie-invariant `bfdr`, defaults to `flag_method = "bfdr"` (`fdr_thresh = 0.05`), and withholds flags for unchecked identification and reports single-cell threshold precision in `prob_threshold_stable`. BFDR decision stability is unassessed.
   - `fit_cmdstan_cached()` defaults to `cache_dir = tempdir()`, runs fresh chains when `seed = NULL` under `refit = "on_change"`, and stores relative `csv_rel` paths so moved cache directories remain loadable via `cmdstanr::as_cmdstan_fit()`.

4. **Test Organization (`tests/testthat/`)**
   - Tests are organized by module: `test-api-compatibility.R`, `test-rsm-pcm.R`, `test-rsm-priors.R`, `test-pcm-steps.R`, `test-heterogeneous-pcm.R`, `test-interactions.R`, `test-design-check.R`, `test-simulation.R`, `test-diagnostics-regression.R`, `test-bias-probability-diagnostics.R`, `test-predict-and-fair-scores.R`, `test-model-comparison.R`, `test-compile-and-cache.R`, `test-formula-and-contracts.R`, `test-residuals-and-fit.R`, `test-cv-and-streaming.R`, `test-methodology-and-roles.R`, `test-fixes-and-extensions.R`, `test-hardcoded-models.R`, and `test-tutorial-and-ux.R`.
