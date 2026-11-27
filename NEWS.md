# BayesMFRM Changelog

## 0.1.0.9008 (Development)

- Added `tab_bmfrm()` for compact console, knitted and standalone HTML tables.
  Named fits appear side by side with posterior means, credible intervals,
  rating counts and diagnostic status.

- `simulate_bmfrm()` reuses category probabilities across replications with
  fixed parameters. It and `pp_bmfrm()` compute cumulative probabilities by
  rubric in matrix operations. Posterior category predictions for RSM and
  uniform-category PCM reuse the category weights when returning probabilities.

- Added `rater_agreement()` and `plot_rater_agreement()` to compare observed
  exact agreement with model expectations and draw-matched predictive intervals
  on jointly scored performances. The writing example includes the results.

- Added `facet_separation()` for EAP or conventional Rasch separation, strata
  and reliability, and `facet_contrasts()` for within-draw main-effect
  differences. The writing vignette illustrates both with three raters.

- Added a writing-assessment example for HDIs, support intervals, ROPE decisions and Bayes factors. The vignette compares the evidence for one interaction with evidence for a prior scale and with held-out predictive accuracy. It also shows how the HDI changes on the odds scale and how to compare two raters' pass probabilities at three illustrative cut scores.

- Added `pp_bmfrm()` as a short name for `prior_predict_bmfrm()`. The writing vignette shows how to visualize criterion-by-rater draws with tidybayes and quantify their direction and practical equivalence with bayestestR.

- Interaction flags now require a rating and confirmed likelihood identification with either screening rule. An observed cell with an unchecked design has `flag_bias = "identification_unchecked"` and `bfdr = NA`. Cells without ratings are always labeled `"unobserved_cell"`.
- In `summarise_bias()`, `decision_stable` and `decision_status` refer to the selected rule. Under BFDR, set-level stability is unassessed. Code using the former fields for single-cell Monte Carlo precision should read `prob_threshold_stable` and `prob_threshold_status`.
- Added a worked simulated rating study, runnable help examples and `plot_bias()` for interaction means, credible intervals, practical thresholds and screening decisions.
- Observed ratings, predictions and standardized fair scores now share the chunked posterior category calculations. Standardized fair scores request only the expected score.

- Added `category_curves()` and `category_thresholds()` for probabilities at fixed abilities and adjacent-category crossings, including fitted interactions and item-specific PCM rubrics. The three ggplot2 displays are `plot_category_curves()`, the ridgeline `plot_wright_map()`, and `plot_pathway_map()` with separate intervals for facet measures and replicated fit.

- Added `category_contrast()` and `plot_category_contrast()` for paired, pointwise differences in score probabilities between two fitted profiles. Both show the number of observed ratings per profile and identify extrapolations. An optional margin in percentage points gives the posterior probability of practical equivalence at each score and ability.

- Added `test_information()` for ability-specific Fisher information under a planned set of ratings per examinee, with draw-matched posterior bounds and support for mixed-category PCM rubrics.

- Fitted models now support `posterior::as_draws()` and tidybayes extraction. `as_mirt()` exports a planned set of ratings at a chosen posterior draw for category traces and test information in mirt.

- Added a runnable priorsense example for the four-category writing assessment, including separate perturbations of the interaction prior and rating likelihood and a direct Gaussian-prior refit in the vignette.

### Breaking Default & API Improvements
- **Hierarchical ability by default**: `bmfrm()`, `parse_bmfrm_formula()`, `get_default_priors()`, and `simulate_bmfrm()` now default to `hierarchical = person_facet` (partial pooling of examinee abilities with posterior uncertainty in the population mean and SD). Pass `hierarchical = character(0)` (or supply `prior(..., class = "theta")`) for unpooled fixed-scale priors.
- **Compact CSV outputs by default**: `save_gq` now defaults to `"minimal"`, omitting observation-level Stan `generated quantities` and streaming residuals, fit statistics, fair scores, PSIS-LOO, and PPC on demand in R.
- **CRAN-compliant temporary caching & unseeded refits**: `cache_dir` now defaults to `tempdir()` rather than writing `"stan_cache"` into the working directory (use `tools::R_user_dir("BayesMFRM", "cache")` for persistent caching). Unseeded calls (`seed = NULL`) with `refit = "on_change"` now draw fresh Markov chains instead of reusing a cached fit, and cached `.rds` payloads store portable relative CSV paths (`csv_rel`) that relocate cleanly when a cache folder is moved.
- **Role-agnostic facet semantics**: Added `person_facet` (default `"person"`) and `adjust_facets` across `bmfrm()`, `parse_bmfrm_formula()`, `check_design_bmfrm()`, `simulate_bmfrm()`, `prior_predict_bmfrm()`, `fair_scores()`, `predict()`, `loo()`, and `cv_bmfrm()`, plus `prior(..., class = "facet", facet = "<name>")` so models no longer depend on column names `"person"`, `"item"`, or `"rater"` (legacy aliases `"item"`, `"rater"`, and `"facet_<name>"` remain supported).
- **Weakly-informative step threshold default**: `get_default_priors()` now sets `tau ~ normal(0, 2)` (previously `normal(0, 1)`) so multi-category rubrics spanning +/-3 logits are not overly compressed.

### Calibrated Fit & Bias Diagnostics
- **Sinharay (2006) posterior predictive p-values in `facet_fit()`**: `facet_fit()` now compares realized and replicated mean-square discrepancies (`MSQ(y, psi_d)` vs. `MSQ(y_rep_d, psi_d)`) on matched posterior draws and reports calibrated `infit_ppp`, `outfit_ppp`, and predictive reference bounds (`infit_rep_lo..hi`, `outfit_rep_lo..hi`) as primary outputs, eliminating small-`n_obs` posterior-variance inflation. Asymptotic Wilson-Hilferty ZSTD columns are opt-in via `include_zstd = TRUE`.
- **Zero-observation cells & BFDR multiplicity control in `summarise_bias()`**: Cells with `n_observations == 0` are now labeled `"unobserved_cell"` and excluded from Bayesian FDR (`bfdr`) and bias flagging even when algebraically determined by sum-to-zero constraints. `summarise_bias()` now defaults to `flag_method = "bfdr"` (`fdr_thresh = 0.05`) and adds configurable `ci_level` (`ci_lower`, `ci_upper`, also in `summarise_bias_ir()`).
- **Documented FACETS Fair Average vs. `fair_scores()`**: Documented the exact distinction between Linacre's single-origin Fair Average (`E[Y | theta_p, item = 0, rater = 0]`), `fair_scores(standardize = FALSE)`, and `fair_scores(standardize = TRUE)`, as well as double-centering spillover under latent local shrinkage (`local_shrinkage = "student_t"` / `"horseshoe"`).
- **Configurable horseshoe slab & convergence diagnostics**: Added `slab_df` and `slab_scale` to `prior()` for regularized horseshoe interactions; separated `max_treedepth` hits as `efficiency_warning` rather than `converged = FALSE` in `check_convergence()` while recording `rhat_max`, `ess_min`, and `ebfmi_min`.

### Performance, Security, & Maintainability
- **Vectorized Stan templates & AST prior validation**: Vectorized main-facet linear predictors, precomputed cumulative step thresholds outside observation loops, removed redundant `softmax`/`exp` passes in `generated quantities`, replaced greedy regex stripping with an explicit `// {{GQ_BLOCK}}` placeholder, removed unverified sibling path fallbacks, and enforced single-call AST/whitelist validation in `prior()` to prevent Stan code injection.
- **Simplified internals & module-aligned tests**: Extracted shared person-level fair-score and CV planning/refit helpers, reorganized the regression test suite by module, added Simulation-Based Calibration (`simulations/sbc_calibration.R`) and Gauss-Hermite MML-EAP benchmarks (`simulations/bayesian_benefits.R`, `simulations/conquest_comparison.R`), and unified documentation and vignettes in English with precomputed output examples.

## 0.1.0.9007 (Development)
- Hardened input validation, reserved Stan identifier checks, data suitability checks, signature-verified fit/CV caching (keeping chain CSV deletion opt-in in `clean_model_cache()`), column-wise unseen-level interaction streaming with bounded offset caches, and S3 generics (`stancode()`, `standata()`).

## 0.1.0.9006 (Development)
- Added heterogeneous-category Partial Credit Model support (`K` as a named integer vector per `threshold_facet` level) via `inst/stan/bmfrm_pcm_hetero.stan`, preserving unobserved rubric categories in softmax denominators while returning `NA` above `K_j` in `predict(..., type = "probs")`.

## 0.1.0.9005 (Development)
- Added `predict.bmfrm_fit()`, blueprint-weighted and decomposed `fair_scores()`, streaming `loo::loo.bmfrm_fit()`, target-aware `cv_bmfrm()` (rating-level and grouped unseen-level K-fold CV), and `compare_bmfrm()` with `obs_id` alignment.

## 0.1.0.9004
- Added `save_gq = "minimal"` R-side streaming evaluation, hierarchical main facets, hierarchical PCM (`hierarchical = "tau"`), Student-t and regularized horseshoe interaction priors, standardized `fair_scores()`, tie-invariant `bfdr`, autocorrelation-aware bias probability MCSE, and Schur-complement design rank checks.

## 0.1.0.9001 - 0.1.0.9003
- Introduced `hierarchical` and `priors` (`scale = "raw"` / `"marginal"`), within-level zero-sum PCM identification, hierarchical two-way interactions, and `summarise_bias(..., delta = ...)`.

## 0.1.0
- Initial experimental release supporting basic RSM/PCM Stan generation, bundled reference models, residuals, infit/outfit, fair scores, and posterior predictive checks.
