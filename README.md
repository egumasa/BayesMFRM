# BayesMFRM

<!-- badges: start -->
[![Experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

**BayesMFRM** (v0.1.0.9007, experimental) is an R package for Bayesian Many-Facet Rasch Models (MFRM) using Stan (`cmdstanr`). It supports Rating Scale Models (RSM), Partial Credit Models (PCM), hierarchical facets/steps, two-way interaction (bias) analysis with local shrinkage, structural design-identifiability checks, prior predictive checks, unified posterior prediction, target-aware model comparison (PSIS-LOO and K-fold CV), and scalable FACETS-style diagnostics.

> For a full step-by-step tutorial (research question -> model choices -> design check -> prior predictive check -> fitting -> diagnostics -> sensitivity analysis -> manuscript reporting), read the installed vignette:
> ```r
> vignette("getting-started", package = "BayesMFRM")
> ```
> Version history and breaking changes are maintained separately in [NEWS.md](NEWS.md).

---

## Installation

BayesMFRM requires [CmdStan](https://mc-stan.org/cmdstanr/) and a C++ toolchain.

```r
# 1. Install cmdstanr and CmdStan (once)
install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
cmdstanr::install_cmdstan()

# 2. Install BayesMFRM (build_vignettes = TRUE to install the tutorial vignette)
remotes::install_github("egumasa/BayesMFRM", build_vignettes = TRUE)
```

---

## Quick Start (Self-Contained & Runnable)

Every step below runs directly after `library(BayesMFRM)` using existing public functions—no external data file required. Always pass `K` explicitly as the **rubric maximum** (e.g., `K = 6` for a 1–6 scale), not just the highest score observed in a sample.

```r
library(BayesMFRM)

# 1. Create a small bridged rating plan (24 persons x 4 items x 4 raters, 8 ratings/person)
set.seed(42)
full_grid <- expand.grid(
  person = sprintf("P%02d", 1:24),
  item   = sprintf("I%02d", 1:4),
  rater  = sprintf("R%02d", 1:4),
  stringsAsFactors = FALSE
)
plan <- do.call(rbind, lapply(split(full_grid, full_grid$person), function(df) {
  df[sample(nrow(df), 8L), , drop = FALSE]
}))
rownames(plan) <- NULL

# 2. Check structural design identifiability BEFORE fitting
design <- check_design_bmfrm(score ~ person + item + rater, plan, K = 6)
print(design)

# 3. Prior predictive check: verify your priors imply plausible category proportions (1..6)
priors_rsm <- c(
  prior("normal(0, 0.5)", "mean", facet = "person"),
  prior("lognormal(-0.5, 0.3)", "sd", facet = "person"),
  prior("normal(0, 0.5)", "item", scale = "marginal"),
  prior("normal(0, 0.5)", "rater", scale = "marginal"),
  prior("normal(0, 0.6)", "tau", scale = "marginal")
)
ppc_prior <- prior_predict_bmfrm(
  score ~ person + item + rater, design = plan, K = 6,
  hierarchical = "person", priors = priors_rsm, ndraws = 100, seed = 43
)
round(prop.table(table(ppc_prior$y_rep)), 3)

# 4. Generate a sample dataset (replace with your own data.frame in real studies)
sim <- simulate_bmfrm(
  score ~ person + item + rater, design = plan, K = 6,
  hierarchical = "person", priors = priors_rsm, seed = 44
)
ratings <- sim$data
check_data_suitability(ratings, "score", c("person", "item", "rater"))

# 5. Fit the model (always specify K explicitly; use save_gq = "minimal" on large datasets)
fit <- bmfrm(
  score ~ person + item + rater,
  data = ratings,
  K = 6,
  family = "rsm",
  hierarchical = "person",
  priors = priors_rsm,
  design_check = "error",
  save_gq = "minimal",
  chains = 4,
  iter = 2000,
  seed = 45
)

# 6. Inspect in priority order: (1) Design & MCMC -> (2) Data fit -> (3) Substantive outputs
print(fit)
summary(fit)
bayesplot::pp_check(fit, group = "rater")          # Conditional posterior predictive check
rater_fit <- facet_fit(fit, "rater")               # Draw-consistent infit/outfit + ZSTD
fair_std  <- fair_scores(fit, standardize = TRUE)  # Common-item standardized fair scores
```

> **Note on `pp_check()`**: `BayesMFRM` registers the `pp_check.bmfrm_fit` S3 method for `bayesplot::pp_check()` without re-exporting the generic. Call **`bayesplot::pp_check(fit)`** (or load `library(bayesplot)` first).

---

## Choosing Model Options (Why & When)

| Option | Choices | When to use & what it assumes |
| --- | --- | --- |
| `K` | Integer rubric max (e.g., `K = 6`) | **Always specify explicitly.** Omitting `K` infers `max(score)` from the sample (with a message), which silently drops unobserved upper rubric categories. |
| `family` | `"rsm"` (default) vs `"pcm"` | Use `"rsm"` when all items/raters share one rating-scale step structure. Use `"pcm"` with `threshold_facet = "item"` (or `"rater"`) when category widths differ across levels. Add `hierarchical = c(..., "tau")` for partial pooling toward common steps (`tau_common`). |
| `hierarchical` | `character(0)` vs `"person"`, `"rater"`, `"item:rater"`, `"tau"` | Default `character(0)` uses fixed-scale priors. Include `"person"` for standard population ability pooling; include centered facets (e.g., `"rater"`) or interactions (`"item:rater"`) when partial pooling across levels/cells is theoretically appropriate and enough levels exist to inform the SD. |
| `scale` in `prior()` | `"raw"` (default) vs `"marginal"` | Under sum-to-zero constraints on $J$ levels, `"raw"` specifies the pre-centering SD (induced marginal SD is $\sigma\sqrt{1 - 1/J}$). Use `"marginal"` on Gaussian / unit-variance Student-t priors when you want the prior SD of each constrained level/cell to stay $\sigma$ regardless of $J$. |
| `local_shrinkage` | `"none"` (default), `"student_t"`, `"horseshoe"` | For hierarchical 2-way bias interactions (`item:rater`). `"none"` shrinks all cells isotropically; `"student_t"` (`scale = "marginal"`) allows heavier tails at matched cell SD; `"horseshoe"` (`scale = "raw"` only) concentrates mass near zero while leaving large deviations less shrunk. Always check induced prior spread with `prior_predict_bmfrm()`. |
| `delta` & `prob_thresh` in `summarise_bias()` | e.g., `delta = 0.3`, `prob_thresh = 0.95` | `delta` is the minimum **practically meaningful** bias on the adjacent-logit scale (avoid testing only `delta = 0`). `prob_thresh` is the directional posterior probability threshold (`0.975` corresponds to a two-sided 95% CrI excluding `[-delta, delta]`). Inspect `bfdr` for multiplicity. |
| `save_gq` | `"full"` (default) vs `"minimal"` | `"minimal"` omits observation-length `generated quantities` in Stan CSVs and computes `residuals()`, `facet_fit()`, `fair_scores()`, and `bayesplot::pp_check()` on demand in memory-bounded chunks—recommended whenever $N \ge 1{,}000$. |

---

## Diagnostic Triage & Interpretation Guardrails

Always inspect results in the following order—do **not** skip ahead to substantive claims when an earlier stage fails:

1. **Design Identification (`check_design_bmfrm()`, `fit$design`)**
   - *What it checks*: Algebraic rank of the conditional likelihood under sum-to-zero constraints (connectivity alone is insufficient).
   - *If it fails*: Fix the rating plan or restrict reporting strictly to identifiable contrasts (`likelihood_identified == TRUE` in `fair_scores()` and `summarise_bias()`). **MCMC convergence does not fix rank deficiency.**
2. **MCMC Convergence (`print(fit)`, `summary(fit)`)**
   - *What it checks*: Sampling stability ($\hat{R} < 1.01$, bulk/tail ESS $\ge 400$, 0 divergences, E-BFMI $\ge 0.3$, MCSE).
   - *If it fails*: Increase `adapt_delta` (e.g., `0.95`), `max_treedepth`, or iterations, or inspect overly diffuse priors before interpreting estimates. **Convergence means the sampler worked—not that the model fits your data.**
3. **Data Fit & Residual Structure (`bayesplot::pp_check()`, `facet_fit()`, `residual_cor_criteria()`)**
   - *What it checks*: Whether predicted category distributions (overall and by `group = "rater"` / `"item"`), infit/outfit MSQ, and cross-criterion residual correlations match observed ratings.
   - *Caution*: Wilson-Hilferty ZSTD is an approximation (especially for small `n_obs` or sparse designs); pair it with conditional `bayesplot::pp_check()` and MSQ credible intervals.
4. **Substantive Effects & Sensitivity (`fair_scores()`, `summarise_bias()`, alternative priors)**
   - *Fair scores*: Use `standardize = TRUE` to remove item-assignment imbalance across examinees; check `identification_basis`.
   - *Bias flags*: A flagged `item:rater` cell indicates a model-based relative interaction exceeding `delta`, **not** an automatic mandate to exclude a rater—verify cell observation count (`n_observations`), `likelihood_identified`, `bfdr`, and prior sensitivity first.

---

## Bias & Hierarchical PCM Example

```r
# Hierarchical person + rater + unit-variance Student-t bias interaction
priors_bias <- c(
  prior("normal(0, 1.5)", "mean", facet = "person"),
  prior("lognormal(-0.2, 0.4)", "sd", facet = "person"),
  prior("normal(0, 1)", "item", scale = "marginal"),
  prior("normal(0, 0.8)", "sd", facet = "rater", scale = "marginal"),
  prior("normal(0, 1)", "tau", scale = "marginal"),
  prior("normal(0, 0.5)", "sd", facet = "item:rater",
        scale = "marginal", local_shrinkage = "student_t")
)

fit_bias <- bmfrm(
  score ~ person + item + rater + item:rater,
  data = ratings,
  K = 6,
  hierarchical = c("person", "rater", "item:rater"),
  priors = priors_bias,
  save_gq = "minimal",
  chains = 4, iter = 2000, adapt_delta = 0.95, seed = 46
)

summarise_bias(fit_bias, "item:rater", prob_thresh = 0.95, delta = 0.3)
fair_scores(fit_bias, adjust_facets = "rater", standardize = TRUE)
```

Bias summaries include probability MCSE and `decision_stable` by default:

```r
bias <- summarise_bias(fit_bias, "item:rater", delta = 0.3,
                       prob_thresh = 0.95, mcse_multiplier = 2)
bias[, c("level1", "level2", "prob_gt_delta", "mcse_prob_gt_delta",
         "flag_bias", "decision_stable", "decision_status", "likelihood_identified")]
```

`mcse_prob_gt_delta` estimates simulation error in the probability of a positive
effect larger than `delta`, accounting for the original chains and their
autocorrelation. Negative and outside-`delta` probabilities have corresponding
MCSE columns. The approximate band `probability +/- 2 * MCSE` must lie wholly on
the same side of `prob_thresh` for a stable directional decision. A stable
`no_clear_bias` decision requires both directional upper bounds below the
threshold; it **does not prove that bias is absent**.

`decision_stable = FALSE` means the decision is sensitive to Monte Carlo error;
inspect mixing and consider longer sampling. `NA` means stability cannot be
assessed: check `decision_status` for unavailable chain information, failed or
incomplete sampler diagnostics, or unavailable MCSE. An indicator constant
within any chain has MCSE `NA`, including estimated probabilities of 0 or 1.
Numerical stability does not establish likelihood identification, validate the
BFDR selection, or establish model fit. Check `likelihood_identified`, `bfdr`,
PPC, and prior sensitivity separately. These bands are not credible intervals
for the effect. `diagnostics = FALSE` returns the original summary columns;
`summarise_bias_ir()` also supports these arguments for legacy raw fits.

---

## Data Requirements

- **Score column**: Finite integers `1..K` (1-indexed; add 1 if your rubric starts at 0).
- **Facet columns**: Factors, characters, or integers with no `NA` values (`NA` scores are dropped with a message; `NA` facets error).
- **Custom Stan files (`stan_file`)**: Supported for built-in Stan `generated quantities`. Unverified custom Stan files intentionally disable R-side recomputation (`standardize = TRUE` or custom `adjust_facets`) so linear-predictor differences cannot be miscomputed silently.

---

## Documentation & Development

- **Interactive English landing page (offline HTML, real `sirt::data.ratings1` results)**: [docs/index.html](docs/index.html)
- **Target-aware model comparison guide & reproducible script (`sirt::data.ratings1`, PSIS-LOO + 5-fold person CV)**: [docs/model-comparison-ratings1.md](docs/model-comparison-ratings1.md) / [docs/compare-ratings1-cv.R](docs/compare-ratings1-cv.R)
- **User tutorial (installed vignette)**: `vignette("getting-started", package = "BayesMFRM")` ([vignettes/getting-started.Rmd](vignettes/getting-started.Rmd))
- **Version history**: [NEWS.md](NEWS.md)
- **Developer & verification notes**: [README.dev.md](README.dev.md)
- **Local unit tests**: `Rscript --vanilla tools/test-local.R`

## License

Licensed under [CC BY-NC-SA 4.0](https://creativecommons.org/licenses/by-nc-sa/4.0/).
