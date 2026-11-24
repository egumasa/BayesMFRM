# BayesMFRM Changelog

## 0.1.0.9005 (Development)

### Prediction, Fair-Score Extensions, & Target-Aware Model Comparison
- Added unified `predict.bmfrm_fit()` S3 method supporting `type = c("epred", "predict", "probs", "linpred")`, explicit prediction-target metadata (`bmfrm_prediction_meta`), population sampling for unseen hierarchical levels (`allow_new_levels = TRUE, sample_new_levels = "hierarchical"` vs `"zero"`), sum-to-zero constrained unseen-level interaction draws, and memory-bounded chunked streaming.
- Extended `fair_scores()` with `reference_weights` / `reference_grid$weight` blueprint weighting, exact draw-matched decomposition (`decompose = TRUE`), draw-matrix return (`return_draws = TRUE`) with preserved chain/draw IDs, and separate row-level vs. summary-level identifiability handling for zero-weight reference rows.
- Added target-aware Bayesian model comparison via `loo::loo.bmfrm_fit()` (true $O(D)$ streaming via `loo::loo.function()`, chain-aware `exp(log_lik)` `r_eff`, singleton-level detection, Pareto-$k$ / MCSE / MCMC convergence metadata), `cv_bmfrm()` (rating-level and grouped unseen-level K-fold CV with strict collateral-unseen target guardrails, dropleveled training folds, exact double-`logmeanexp` Monte Carlo marginalization, single-compile binary sharing, verified per-fold caching, and `refit_seed` for fixed-fold $M$/seed sensitivity audits), and `compare_bmfrm()` (strict compatibility and convergence gates, group-joint vs. group-clustered pointwise paired SEs, and explicit approximation caveats).
- Added pre-preprocessing observation identity tracking (`obs_id` argument in `bmfrm()` and `prepare_data_bmfrm()`, preserving explicit IDs, non-default unique `rownames`, or pre-NA row indices): `compare_bmfrm()` verifies observation correspondence via per-ID response and shared-facet/metadata hashes (catching same-score row swaps and aligning permuted datasets with preserved IDs), and `cv_bmfrm()$pointwise` outputs both the preserved `obs_id` and `row_index`.
- Added an empirical model-comparison guide and reproducible script on `sirt::data.ratings1` (`docs/model-comparison-ratings1.md`, `docs/compare-ratings1-cv.R`, `docs/ratings1_person_folds.csv`) covering PSIS-LOO, exact LOO refitting of high Pareto-$k$ observations, 5-fold person CV, multi-$M$/multi-seed/1D quadrature precision checks, and separation of predictive support from practical fair-score impact.
- Resolved `R CMD check` static-analysis notes and documentation warnings (`save_gq` in `parse_bmfrm_formula.Rd`, qualified `stats::sd()`, NSE global declarations in `R/globals.R`, and `.Rbuildignore` exclusions including `.Rprofile` so source builds do not trigger an empty `renv` bootstrap).

## 0.1.0.9004

### User Experience & Documentation
- Added the end-to-end introductory vignette `vignettes/getting-started.Rmd` covering model choice rationale, explicit `K` specification, design identifiability, prior predictive checks, diagnostic triage, sensitivity analysis, and manuscript reporting templates.
- Unified `README.md` across root and development directories with a self-contained runnable Quick Start (replacing the missing `data(ratings)` call with `simulate_bmfrm()`), consistent `bayesplot::pp_check()` usage, and explicit option/diagnostic decision tables.
- Separated development history into `NEWS.md` and updated `README.dev.md` to reflect current capabilities cleanly.
- `parse_bmfrm_formula()` now emits an informative message whenever `K` is inferred from `max(score)` rather than specified explicitly.
- `print.bmfrm_fit()` and `summary.bmfrm_fit()` now report **Design Identification** separately from **MCMC Convergence** and list prioritized next steps (`bayesplot::pp_check()`, `facet_fit()`, `fair_scores()`, `summarise_bias()`).

### Models, Diagnostics, & Scalability
- Added chain-aware probability MCSE and numerical decision stability to `summarise_bias()` and `summarise_bias_ir()`, using the shared internal indicator calculation. `diagnostics = TRUE` and `mcse_multiplier = 2` add annotations without changing flags or BFDR; `diagnostics = FALSE` retains the original columns. Missing chain information, failed/incomplete fit checks, and constant indicators cannot certify stable decisions. Numerical stability remains separate from likelihood identification and BFDR selection stability.
- Added `save_gq = "minimal"` with memory-bounded chunked streaming evaluation in R (`residuals()`, `facet_fit()`, `fair_scores()`, `bayesplot::pp_check()`).
- Added hierarchical main facets (`hierarchical = c("person", "rater", ...)`), hierarchical PCM (`family = "pcm", hierarchical = "tau"` with `tau_common` + level deviations), and cell-space local shrinkage for 2-way interactions (`local_shrinkage = "student_t"` with unit-variance normalization under `scale = "marginal"`, and `"horseshoe"` under `scale = "raw"`).
- Added `adjust_facets`, `standardize = TRUE`, `reference_grid`, `max_grid_cells`, `max_output_rows`, `likelihood_identified`, and `identification_basis` to `fair_scores()` (supporting person-level PCM, interactions involving `person`, and binary `K = 2` hierarchical PCM).
- Fixed Wright & Masters (1982) Wilson-Hilferty ZSTD formula, draw-level kurtosis $q_d$, multi-model `add_fit_cri_bayes()` joins, and bare-column NSE in `facet_infit_outfit()` and `add_fit_cri_bayes()`.
- Unified bias flagging across `summarise_bias()` and `summarise_bias_ir()` via a single directional `prob_thresh` criterion (`classify_bias_flags()`) and added tie-invariant Bayesian FDR (`bfdr`).
- Added an explicit internal evaluation contract (`resolve_eval_contract()`) blocking silent R-side recomputation for unverified custom Stan files while preserving verified bundled models.
- Added exact Schur-complement (`method = "auto"` / `"schur"`) design checking with within-person TSQR, factorized Kronecker operations (skipping zero-DOF 1-level bases and checking limits before allocation), and exact full-design condition numbers via Lanczos.
- Unified reserved Stan facet-name validation across non-hierarchical and hierarchical paths, and fixed `fit_cmdstan_cached()` when called with no `...` arguments.

## 0.1.0.9003
- Added hierarchical Gaussian two-way interactions via `hierarchical` and `priors` (`class = "sd"`), interaction-aware `check_design_bmfrm()`, and `summarise_bias(..., delta = ...)`.

## 0.1.0.9002
- Corrected PCM identification so step deviations sum to zero **within each `threshold_facet` level**, separating overall facet location from category step shape (breaking change relative to earlier across-level PCM centering).

## 0.1.0.9001
- Replaced the experimental `rsm_model()` constructor and `model` argument with `hierarchical` and `priors` (`scale = "raw"` / `"marginal"`), adding signature-verified fit caching and fixed-point Stan log-density/gradient parity checks.

## 0.1.0
- Initial experimental release supporting basic generated RSM/PCM models, bundled hardcoded models, residuals, infit/outfit, fair scores, and PPC.
