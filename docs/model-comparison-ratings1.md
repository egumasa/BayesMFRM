# Target-Aware Bayesian Model Comparison on `sirt::data.ratings1`

This guide demonstrates how to run and interpret `BayesMFRM`'s model-comparison workflow (`loo::loo()`, `cv_bmfrm()`, and `compare_bmfrm()`) on the empirical writing-assessment dataset `sirt::data.ratings1` (135 students $\times$ 5 analytic criteria $\times$ 7 raters; 1,370 ratings on a 4-category rubric).

### Reproducibility Assets & Settings
- **Reproduction script**: [compare-ratings1-cv.R](compare-ratings1-cv.R)
- **Fixed 5-fold student assignment (135 examinees)**: [ratings1_person_folds.csv](ratings1_person_folds.csv)
- **Local checkout priority**: When run from the repository root, [compare-ratings1-cv.R](compare-ratings1-cv.R) loads the local source tree via `pkgload::load_all(".")` so results reflect the exact checked-out code.
- **Sampler settings**:
  - Full-data fits (for PSIS-LOO): `Main` (`iter = 1500, warmup = 500, adapt_delta = 0.95, seed = 202611`), `Student-t` (`iter = 2000, warmup = 1000, adapt_delta = 0.99, seed = 202612`), `Gaussian` (`iter = 2000, warmup = 1000, adapt_delta = 0.95, seed = 202613`), all with `chains = 4, max_treedepth = 12, obs_id = "obs_uid"`.
  - 5-fold person CV: `Main` (`iter = 1500, warmup = 500`), `Gaussian` (`iter = 2500, warmup = 1000, adapt_delta = 0.98`), `Student-t` (`iter = 2500, warmup = 1000, adapt_delta = 0.995`), with fold MCMC seeds locked to `refit_seed = 20261124` (`refit_seed + fold * 1000`).

---

## 0. Define the Prediction Target First

Before computing any information criterion, decide **which out-of-sample prediction task** matches your substantive question, and pass `obs_id = "obs_uid"` to `bmfrm()` so `compare_bmfrm()` verifies exact observation alignment across models.

| Prediction Task | Function & Specification | Train–Test Split Structure | Feasibility on `data.ratings1` |
|---|---|---|---|
| **New rating among known examinees, criteria, and raters** | `loo::loo(fit)` (PSIS-LOO) | Leave 1 rating out; all facet levels remain observed in training | **Valid** (all 1,370 observations leave every level connected) |
| **Generalization to new (unseen) examinees** | `cv_bmfrm(fit, target = "person", folds = ...)` | Hold out 5 groups of 27 students | **Valid** (all 5 training folds are connected and full-rank: `120/120` main, `144/144` interaction) |
| **Generalization strictly to new raters** | `cv_bmfrm(fit, target = "rater", ...)` | Hold out raters | **Invalid (correctly rejected)** |

### Why `target = "rater"` Must Stop

In `sirt::data.ratings1`, **89 of 135 students are rated by only a single rater**. Holding out any rater simultaneously removes every student rated exclusively by that rater from the training fold, turning the prediction task from *"new rater on known students"* into *"new rater $\times$ new students"*. With `strict_target = TRUE` (the default), `cv_bmfrm()` detects collateral unseen facets before fitting and halts with an informative error.

---

## 1. PSIS-LOO (Known Facet Levels) & Exact-LOO Verification of High Pareto-$k$

All three full-data fits converged cleanly ($\hat{R} \le 1.01$, 0 divergences, E-BFMI $> 0.61$). Streaming PSIS-LOO (`loo::loo(fit)`, `actual_target = "new_rating_known_levels"`) yields:

| Model | `elpd_loo` (SE) | `p_loo` (SE) | MCSE | Pareto-$k > 0.7$ | Min PSIS ESS | `elpd_diff` vs. Top (SE) | Diagnostic Note |
|---|---:|---:|---:|---:|---:|---:|---|
| **Gaussian Bias RSM** (`item:rater`, `none`) | `-1169.8` (`26.1`) | `128.1` (`5.4`) | `0.220` | `0 / 1370` (max `0.687`) | `376.4` | `0.0` (`0.0`) | Good |
| **Student-t Bias RSM** (`item:rater`, `student_t`) | `-1170.1` (`26.2`) | `128.5` (`5.4`) | `0.248` | `0 / 1370` (max `0.612`) | `565.1` | `-0.3` (`0.9`) | `\|elpd_diff\| < 4` |
| **Main-Effects RSM** (no interaction) | `-1200.7` (`26.1`) | `111.7` (`5.0`) | `NA` | `1 / 1370` (`k = 0.734`, obs 39) | `404.6` | `-30.9` (`8.7`) | `1 k_psis > 0.7` |

### Verifying High Pareto-$k$ (Observation 39) via Exact LOO Refit

Observation 39 in `Main` (student `200010210`, criterion `k4`, rater `db07`, score `3`) has `Pareto-k = 0.7338` because this student's five scores are `4, 4, 4, 3, 4`: removing the single `3` leaves an all-maximum pattern (`4, 4, 4, 4`), shifting the student's posterior substantially.

Never dismiss a high Pareto-$k$ warning merely because only one observation exceeded `0.7` (and note that 5-fold person CV evaluates a different target and does not validate rating-level PSIS). Refitting `Main` directly on `ratings[-39, ]` (`iter = 2000, warmup = 1000`, $\hat{R} = 1.0067$, `min_ess_bulk = 751.9`, 0 divergences) gives:

- PSIS estimate for observation 39: `-1.3597`
- Exact LOO refit value for observation 39: `-1.3438` (MCSE `0.0069`) — **discrepancy `+0.0158`**
- Exact-corrected total `Main` ELPD: `-1200.6728` $\to$ **`-1200.6569` (`SE = 26.0830`)**
- Paired difference `Main - Gaussian`: `-30.8521 (SE = 8.7066)` $\to$ **`-30.8362 (SE = 8.7064)`**

Direct refitting confirms the PSIS error on observation 39 is negligible (`0.016` nats) and that adding `item:rater` substantially improves rating-level prediction, whereas `Gaussian` and `Student-t` (`-0.3 ± 0.9`, `|elpd_diff| < 4`) show no clear predictive separation.

---

## 2. 5-Fold Person CV (Unseen Examinees) & Monte Carlo Integration Audit

Using the fixed 5-fold split in [ratings1_person_folds.csv](ratings1_person_folds.csv) (27 held-out students per fold), each fold refits the model on training students and marginalizes unseen abilities $\theta_g \sim \mathcal{N}(\mu_\theta^{(d)}, \sigma_\theta^{2(d)})$ across $M$ Monte Carlo draws per posterior draw via double `logmeanexp`.

### 2.1 Sensitivity to Inner Draws $M$, RNG Seed, and 1D Quadrature

On the identical cached fold posteriors (`refit_seed = 20261124`), varying `M` and the integration `seed` on `Main` gives:

| $M$ | Integration Seed | Joint `elpd_group` (SE) | Marginal `elpd_pointwise` | Student-Clustered SE | Naive Rating SE (Ref) | Min Median `mc_ess_m` | Min Draw `mc_ess_m` |
|---:|---:|---:|---:|---:|---:|---:|---:|
| `20` | `20261124` | `-1301.49` (`117.38`) | `-1724.17` | `161.23` | `17.80` | `1.01` (warns) | `1.00` |
| `100` | `20261124` | `-1301.41` (`117.37`) | `-1724.15` | `161.22` | `17.80` | `2.30` (warns) | `1.00` |
| `300` | `20261124` | `-1301.44` (`117.38`) | `-1724.28` | `161.27` | `17.80` | `6.20` (no warning) | `1.01` |
| `300` | `20261199` | `-1301.47` (`117.38`) | `-1724.26` | `161.26` | `17.80` | — | — |

Across all three models, `M = 300` agrees closely across independent seeds and deterministic 1D numerical quadrature:

| Model | `M = 100` (seed 1) | `M = 300` (seed 1) | `M = 300` (seed 2) | Max Per-Student Diff (`M = 300`) |
|---|---:|---:|---:|---:|
| `Main` | `-1301.408` | `-1301.441` | `-1301.474` | `0.0133` |
| `Gaussian` | `-1272.508` | `-1272.501` | `-1272.460` | `0.0160` |
| `Student-t` | `-1271.737` | `-1271.786` | `-1271.743` | `0.0133` |

- **Model contrast stability** (`M = 100`, `M = 300` seeds 1 & 2, and 61-node Gauss-Hermite quadrature):
  - `Gaussian - Main`: `+28.900` to `+29.014` (spread `< 0.12` vs. paired SE `7.36`)
  - `Student-t - Gaussian`: `+0.715` to `+0.771` (spread `< 0.06` vs. paired SE `0.74`)
- **Hardest individual students vs. 801-point deterministic quadrature** (`stats::integrate()` verified):
  - Lowest-ESS student `400130105` (30 ratings): exact `-18.5508` vs. `M = 300` `-18.5560` (error **`-0.0053`**)
  - Three 35-rating students (`700030217`, `300040123`, `400220106`): errors **`0.0019`–`0.0051`**

### Group-Joint Scores vs. Rating-Marginal Scores
1. **Use `score = "group"` as primary for unseen examinees**: Multiple ratings of the same unseen student share $\theta_g$, so joint group log predictive density (`-1301.4`) captures within-examinee dependence whereas marginal pointwise sums (`-1724.3`) do not.
2. **Clustered SEs vs. naive rating SEs**: Treating 1,370 ratings as independent (`SE = 17.8`) understates uncertainty by nearly an order of magnitude relative to student-clustered SEs (`161.3`), though clustered SEs remain finite-sample approximations because held-out folds still share raters, items, and overlapping training sets.

---

## 3. Comparing All Three Models on the Shared 5-Fold Split

| Model | Rank (All Folds) | Convergence (5 Folds) | `elpd_group` (SE, $G=135$) | `elpd_pointwise` (Clustered / Naive SE) |
|---|---|---|---|---|
| **Main-Effects RSM** | `120/120` | **5/5 converged** ($\hat{R} \le 1.008$, ESS $\ge 560$, 0 div) | `-1301.44` (`117.38`) | `-1724.28` (`161.27` / `17.80`) |
| **Gaussian Bias RSM** (`adapt_delta = 0.98`) | `144/144` | **5/5 converged** ($\hat{R} \le 1.008$, ESS $\ge 559$, 0 div) | `-1272.50` (`115.32`) | `-1711.31` (`160.48` / `18.29`) |
| **Student-t Bias RSM** (`adapt_delta = 0.995`) | `144/144` | **4/5 converged** (1/6,000 divergence on Fold 5; $\hat{R} \le 1.009$, ESS $\ge 707$) | `-1271.79` (`115.33`) | `-1711.02` (`160.50` / `18.32`) |

### Strict Comparison (`allow_unconverged = FALSE`, Default)

Comparing the fully converged `Main` and `Gaussian` models:
- **Group joint (`score = "group"`)**: `Main - Gaussian = -28.9 (SE = 7.4)` (~3.9 SE, `+0.215` nats/student).
- **Pointwise (`score = "pointwise"`)**: `Main - Gaussian = -13.0 (SE = 3.6)` (~3.6 SE, `+0.0095` nats/rating).
- Including `Student-t` under default `allow_unconverged = FALSE` raises an explicit error because Fold 5 recorded 1 divergent transition; passing `allow_unconverged = TRUE` for sensitivity inspection shows `Gaussian - Student-t = -0.7 (SE = 0.7)` (`|elpd_diff| < 4`).

### Separating Predictive Support from Practical Fair-Score Impact

A 3.9-SE predictive improvement confirms genuine `item:rater` structure, **not** that examinee fair scores change substantially:
- Across all 135 students on the 1–4 scale, `fair_scores(..., standardize = TRUE)` differs between `Gaussian` and `Main` by **RMS = `0.0090` score points** (median `|diff| = 0.0044`, max `|diff| = 0.0314`), with **Pearson $r = 0.99996$** and **Spearman $\rho = 0.99953$**.
- Substantive takeaway: modeling `item:rater` interactions is well-supported for **rater-by-criterion bias diagnostics and rater training**, while main-effects standardized fair scores already provide practically identical student rankings on this dataset.

---

## 4. Four Common Pitfalls to Avoid in Manuscripts

1. **Do not claim `|elpd_diff| < 4` proves model equivalence**: Small ELPD differences indicate weak evidence to prefer one specification over another (and noisy normal-theory SE approximations), not proof of zero difference.
2. **Do not report unconverged CV folds as winners**: Resolve sampler warnings first or report unconverged runs solely as labeled sensitivity checks.
3. **Never compare raw ELPD totals across PSIS-LOO (`-1170`) and 5-fold Person CV (`-1272`)**: They evaluate different conditioning sets and training sample sizes; compare models only *within* the same prediction target.
4. **Verify high Pareto-$k$ observations and inner Monte Carlo ESS directly**: Use exact LOO refits for $k > 0.7$ and check multiple `M` / `seed` settings on fixed `folds` and `refit_seed`.
