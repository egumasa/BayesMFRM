# Reproducible Target-Aware Model Comparison on sirt::data.ratings1
#
# Reproduces the exact PSIS-LOO, exact LOO refit for Observation 39, and 5-fold
# person-level CV analysis reported in docs/model-comparison-ratings1.md:
#   0. Checking why target = "rater" is invalid here (collateral unseen persons)
#   1. PSIS-LOO ("new_rating_known_levels") across 3 RSM specifications +
#      exact LOO refit verification for Observation 39 (Main RSM, Pareto-k > 0.7)
#   2. 5-fold person-level CV ("new_person") using the exact fixed person->fold
#      partition in docs/ratings1_person_folds.csv + M / seed sensitivity
#   3. 5-fold person-level CV on Gaussian & Student-t Bias RSMs using the EXACT
#      same person split, convergence gating, group-joint vs. rating-marginal
#      scoring, and group-clustered paired SEs
#   4. Practical assessment impact check (per-unit ELPD & fair_scores())
#
# Usage (from repository root):
#   Rscript --vanilla docs/compare-ratings1-cv.R [cache_dir] [output_rds] [folds_csv]
#
# Package loading:
#   When executed from a BayesMFRM source repository (with `pkgload` available),
#   this script ALWAYS loads the local repository source tree FIRST so a stale
#   installed development build cannot shadow the current implementation. Outside
#   a source tree, it falls back to an installed `BayesMFRM` library that provides
#   `cv_bmfrm(..., refit_seed = ...)` and `compare_bmfrm()`.

args <- commandArgs(trailingOnly = TRUE)
cache_dir <- if (length(args) >= 1L) args[1L] else file.path("stan_cache_ratings1")
out_rds   <- if (length(args) >= 2L) args[2L] else file.path(cache_dir, "ratings1_cv_results.rds")
folds_csv <- if (length(args) >= 3L) args[3L] else file.path("docs", "ratings1_person_folds.csv")

detect_repo_root <- function() {
  is_bmfrm_root <- function(dir_path) {
    desc <- file.path(dir_path, "DESCRIPTION")
    file.exists(desc) && any(grepl("^Package:\\s*BayesMFRM\\s*$", readLines(desc, warn = FALSE)))
  }
  if (is_bmfrm_root(".")) return(normalizePath(".", winslash = "/", mustWork = TRUE))
  cmd_all <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_all, value = TRUE)
  if (length(file_arg) >= 1L) {
    script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[1L]), winslash = "/", mustWork = FALSE))
    parent_dir <- dirname(script_dir)
    if (is_bmfrm_root(parent_dir)) return(parent_dir)
  }
  NULL
}

load_bayesmfrm_for_cv <- function() {
  repo_root <- detect_repo_root()
  if (!is.null(repo_root) && requireNamespace("pkgload", quietly = TRUE)) {
    message("Loading local BayesMFRM source tree via pkgload::load_all('", repo_root, "').")
    pkgload::load_all(repo_root, quiet = TRUE)
    return(list(mode = "pkgload_local_source", path = repo_root))
  }
  if (requireNamespace("BayesMFRM", quietly = TRUE) &&
      exists("cv_bmfrm", envir = asNamespace("BayesMFRM"), inherits = FALSE) &&
      "refit_seed" %in% names(formals(BayesMFRM::cv_bmfrm))) {
    suppressPackageStartupMessages(library(BayesMFRM))
    return(list(mode = "installed_library", path = find.package("BayesMFRM")))
  }
  stop(
    "This script requires the current BayesMFRM implementation (with obs_id alignment ",
    "and cv_bmfrm(..., refit_seed = ...)). Run from the repository root with pkgload ",
    "installed, or install the current repository via remotes::install_local('.')."
  )
}

pkg_load_info <- load_bayesmfrm_for_cv()
pkg_load_mode <- pkg_load_info$mode
suppressPackageStartupMessages({
  library(sirt)
  library(tidyr)
})

dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_rds), recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 0. Data preparation, explicit obs_id, & structural check for candidate targets
# -----------------------------------------------------------------------------
data(data.ratings1)
ratings <- tidyr::pivot_longer(
  data.ratings1,
  cols = k1:k5,
  names_to = "item",
  values_to = "score_raw"
)
ratings <- ratings[!is.na(ratings$score_raw), ]
ratings$person  <- as.character(ratings$idstud)
ratings$rater   <- as.character(ratings$rater)
ratings$score   <- as.integer(ratings$score_raw + 1L)
ratings$obs_uid <- sprintf("R%04d", seq_len(nrow(ratings)))
K <- 4L

# Load the EXACT fixed person -> fold mapping used in the published results
if (!file.exists(folds_csv)) {
  stop("Fixed fold mapping CSV not found: ", folds_csv)
}
pf_tab <- utils::read.csv(folds_csv, stringsAsFactors = FALSE)
pf_tab$person <- as.character(pf_tab$person)
pf_tab$fold   <- as.integer(pf_tab$fold)
shared_folds <- pf_tab$fold[match(ratings$person, pf_tab$person)]
if (anyNA(shared_folds)) {
  stop("Fold mapping CSV does not cover all examinees in sirt::data.ratings1")
}

raters_per_person <- tapply(ratings$rater, ratings$person, function(x) length(unique(x)))
cat(sprintf(
  "Dataset: %d ratings | %d persons | %d items | %d raters (loaded via %s)\n",
  nrow(ratings), length(unique(ratings$person)),
  length(unique(ratings$item)), length(unique(ratings$rater)), pkg_load_mode
))
cat(sprintf(
  "Single-rater examinees: %d / %d -> holding out any rater also makes persons unseen.\n\n",
  sum(raters_per_person == 1L), length(raters_per_person)
))

# Structural design checks BEFORE sampling
dc_main <- check_design_bmfrm(score ~ person + item + rater, ratings, K = K)
dc_bias <- check_design_bmfrm(score ~ person + item + rater + item:rater, ratings, K = K)
stopifnot(dc_main$identifiable, dc_bias$identifiable)

# -----------------------------------------------------------------------------
# 1. Fit (or load cached) 3 models on the full dataset & run PSIS-LOO
#    Sampling settings match the exact cached fits reported in the guide:
#      - Main RSM:      iter = 1500, warmup = 500,  adapt_delta = 0.95, seed = 202611
#      - Student-t RSM: iter = 2000, warmup = 1000, adapt_delta = 0.99, seed = 202612
#      - Gaussian RSM:  iter = 2000, warmup = 1000, adapt_delta = 0.95, seed = 202613
# -----------------------------------------------------------------------------
priors_main <- c(
  prior("normal(0, 1.0)", "mean", facet = "person"),
  prior("lognormal(-0.2, 0.4)", "sd", facet = "person"),
  prior("normal(0, 0.8)", "item", scale = "marginal"),
  prior("normal(0, 0.8)", "sd", facet = "rater", scale = "marginal"),
  prior("normal(0, 0.8)", "tau", scale = "marginal")
)
priors_st <- c(
  priors_main,
  list(prior("normal(0, 0.5)", "sd", facet = "item:rater",
             scale = "marginal", local_shrinkage = "student_t"))
)
priors_iso <- c(
  priors_main,
  list(prior("normal(0, 0.5)", "sd", facet = "item:rater",
             scale = "marginal", local_shrinkage = "none"))
)

fit_main <- bmfrm(
  score ~ person + item + rater, data = ratings, K = K,
  hierarchical = c("person", "rater"), priors = priors_main,
  save_gq = "minimal", cache_dir = cache_dir, model_name = "ratings1_main",
  design_check = "error", chains = 4, cores = 4, iter = 1500, warmup = 500,
  adapt_delta = 0.95, max_treedepth = 12, seed = 202611, obs_id = "obs_uid",
  cpp_options = list(PRECOMPILED_HEADERS = "false")
)
fit_st <- bmfrm(
  score ~ person + item + rater + item:rater, data = ratings, K = K,
  hierarchical = c("person", "rater", "item:rater"), priors = priors_st,
  save_gq = "minimal", cache_dir = cache_dir, model_name = "ratings1_bias_st",
  design_check = "error", chains = 4, cores = 4, iter = 2000, warmup = 1000,
  adapt_delta = 0.99, max_treedepth = 12, seed = 202612, obs_id = "obs_uid",
  cpp_options = list(PRECOMPILED_HEADERS = "false")
)
fit_iso <- bmfrm(
  score ~ person + item + rater + item:rater, data = ratings, K = K,
  hierarchical = c("person", "rater", "item:rater"), priors = priors_iso,
  save_gq = "minimal", cache_dir = cache_dir, model_name = "ratings1_bias_iso",
  design_check = "error", chains = 4, cores = 4, iter = 2000, warmup = 1000,
  adapt_delta = 0.95, max_treedepth = 12, seed = 202613, obs_id = "obs_uid",
  cpp_options = list(PRECOMPILED_HEADERS = "false")
)

cat("=== STEP 1: PSIS-LOO (target: new_rating_known_levels) ===\n")
loo_main <- loo::loo(fit_main, use_stan_log_lik = FALSE)
loo_st   <- loo::loo(fit_st,   use_stan_log_lik = FALSE)
loo_iso  <- loo::loo(fit_iso,  use_stan_log_lik = FALSE)

print(loo_main)
print(loo_st)
print(loo_iso)

cmp_loo <- compare_bmfrm(main = loo_main, student_t = loo_st, gaussian = loo_iso)
print(cmp_loo)

# Step 1b: Exact LOO refit for any high Pareto-k observation (Main RSM: obs 39)
# Keeps original `loo_main` untouched and saves exact-refit & hybrid-corrected
# totals, paired differences, SEs, and chain-aware MCSE separately.
bad_k_main <- loo_main$bmfrm_meta$bad_k_ids
exact_loo_checks <- list()
loo_exact_corrected <- NULL
if (length(bad_k_main) > 0L) {
  cat("\n--- Exact LOO refit verification for high Pareto-k observation(s) in Main RSM ---\n")
  loo_exact_dir <- file.path(cache_dir, "loo_exact_obs39")
  dir.create(loo_exact_dir, recursive = TRUE, showWarnings = FALSE)
  pt_main_corr <- as.numeric(loo_main$pointwise[, "elpd_loo"])
  for (obs_idx in bad_k_main) {
    fit_minus <- bmfrm(
      score ~ person + item + rater, data = ratings[-obs_idx, , drop = FALSE], K = K,
      hierarchical = c("person", "rater"), priors = priors_main,
      save_gq = "minimal", cache_dir = loo_exact_dir,
      model_name = paste0("ratings1_main_minus", obs_idx),
      design_check = "error", chains = 4, cores = 4, iter = 2000, warmup = 1000,
      adapt_delta = 0.95, max_treedepth = 12, seed = 202611,
      cpp_options = list(PRECOMPILED_HEADERS = "false")
    )
    params_m <- BayesMFRM:::extract_bmfrm_param_draws(fit_minus)
    ll_exact <- as.numeric(BayesMFRM:::eval_bmfrm_chunk(
      params_m,
      BayesMFRM:::map_eval_indices(ratings[obs_idx, , drop = FALSE], fit_minus$spec, params_m$labels),
      rows = 1L, resp_chunk = as.integer(ratings$score[obs_idx]),
      compute_mu = FALSE, compute_log_lik = TRUE
    )$log_lik)
    elpd_ex <- BayesMFRM:::logmeanexp_vec(ll_exact)
    cl_m <- BayesMFRM:::extract_verified_chain_layout(fit_minus$fit, params_m$D)
    ll_max <- max(ll_exact)
    p_mat <- matrix(exp(ll_exact - ll_max), nrow = cl_m$n_iter, ncol = cl_m$n_chains)
    mcse_ex <- posterior::mcse_mean(p_mat) / mean(p_mat)
    elpd_ps <- as.numeric(loo_main$pointwise[obs_idx, "elpd_loo"])
    pt_main_corr[obs_idx] <- elpd_ex
    cat(sprintf(
      "Obs %d (%s, Pareto-k = %.4f): PSIS = %.6f | Exact = %.6f (MCSE = %.6f) | Exact - PSIS = %+.6f | conv = %s (max_rhat = %.4f, min_ess = %.1f, div = %d)\n",
      obs_idx, ratings$obs_uid[obs_idx], loo_main$diagnostics$pareto_k[obs_idx],
      elpd_ps, elpd_ex, mcse_ex, elpd_ex - elpd_ps,
      fit_minus$convergence$converged, fit_minus$convergence$max_rhat,
      fit_minus$convergence$min_ess_bulk, fit_minus$convergence$n_divergent
    ))
    exact_loo_checks[[as.character(obs_idx)]] <- list(
      obs_idx = obs_idx,
      obs_uid = ratings$obs_uid[obs_idx],
      person = ratings$person[obs_idx],
      item = ratings$item[obs_idx],
      rater = ratings$rater[obs_idx],
      score = ratings$score[obs_idx],
      pareto_k = as.numeric(loo_main$diagnostics$pareto_k[obs_idx]),
      elpd_psis = elpd_ps,
      elpd_exact = elpd_ex,
      mcse_exact = mcse_ex,
      err_exact_minus_psis = elpd_ex - elpd_ps,
      converged = fit_minus$convergence$converged,
      max_rhat = fit_minus$convergence$max_rhat,
      min_ess_bulk = fit_minus$convergence$min_ess_bulk,
      n_divergent = fit_minus$convergence$n_divergent
    )
  }

  pt_iso <- as.numeric(loo_iso$pointwise[, "elpd_loo"])
  pt_st  <- as.numeric(loo_st$pointwise[, "elpd_loo"])
  N_loo <- length(pt_main_corr)
  elpd_main_hybrid <- sum(pt_main_corr)
  se_main_hybrid <- sqrt(N_loo * stats::var(pt_main_corr))
  diff_main_iso <- sum(pt_main_corr - pt_iso)
  se_diff_main_iso <- sqrt(N_loo * stats::var(pt_main_corr - pt_iso))
  diff_main_st <- sum(pt_main_corr - pt_st)
  se_diff_main_st <- sqrt(N_loo * stats::var(pt_main_corr - pt_st))
  loo_exact_corrected <- list(
    replaced_obs_idx = as.integer(bad_k_main),
    replaced_obs_uid = ratings$obs_uid[bad_k_main],
    pointwise_main_hybrid = pt_main_corr,
    elpd_main_hybrid = elpd_main_hybrid,
    se_main_hybrid = se_main_hybrid,
    diff_main_vs_gaussian = diff_main_iso,
    se_diff_main_vs_gaussian = se_diff_main_iso,
    diff_main_vs_student_t = diff_main_st,
    se_diff_main_vs_student_t = se_diff_main_st,
    per_rating_gain_gaussian_vs_main_nat = -diff_main_iso / N_loo,
    geom_prob_ratio_gaussian_vs_main = exp(-diff_main_iso / N_loo)
  )
  cat(sprintf(
    "Hybrid corrected Main ELPD: %.4f (SE = %.4f) | Main - Gaussian = %.4f (SE = %.4f) | +%.5f nat/rating (ratio %.4f)\n",
    elpd_main_hybrid, se_main_hybrid, diff_main_iso, se_diff_main_iso,
    loo_exact_corrected$per_rating_gain_gaussian_vs_main_nat,
    loo_exact_corrected$geom_prob_ratio_gaussian_vs_main
  ))
}

# Demonstrate the collateral-unseen guardrail for target = "rater":
cat("\n--- Collateral-unseen check for target = 'rater' ---\n")
rater_guard <- tryCatch(
  cv_bmfrm(fit_main, target = "rater", K_folds = 7L, seed = 20261124L),
  error = function(e) conditionMessage(e)
)
cat("Expected stop:", rater_guard, "\n\n")

# -----------------------------------------------------------------------------
# 2. Step 2: 5-fold person CV on Main-Effects RSM + M & seed sensitivity check
#    Inherits fit_main's iter = 1500, warmup = 500, adapt_delta = 0.95.
# -----------------------------------------------------------------------------
cat("=== STEP 2: 5-fold person CV (Main-Effects RSM) & M/seed sensitivity ===\n")
cv_main_20     <- suppressWarnings(cv_bmfrm(fit_main, target = "person", folds = shared_folds, M = 20L,  seed = 20261124L, refit_seed = 20261124L))
cv_main_100    <- suppressWarnings(cv_bmfrm(fit_main, target = "person", folds = shared_folds, M = 100L, seed = 20261124L, refit_seed = 20261124L))
cv_main_300    <- cv_bmfrm(fit_main, target = "person", folds = shared_folds, M = 300L, seed = 20261124L, refit_seed = 20261124L)
cv_main_300_s2 <- cv_bmfrm(fit_main, target = "person", folds = shared_folds, M = 300L, seed = 20261199L, refit_seed = 20261124L)

print(cv_main_300)
m_sens <- data.frame(
  M = c(20L, 100L, 300L, 300L),
  seed = c(20261124L, 20261124L, 20261124L, 20261199L),
  elpd_group = c(cv_main_20$elpd_group, cv_main_100$elpd_group, cv_main_300$elpd_group, cv_main_300_s2$elpd_group),
  se_group = c(cv_main_20$se_group, cv_main_100$se_group, cv_main_300$se_group, cv_main_300_s2$se_group),
  elpd_pointwise = c(cv_main_20$elpd_pointwise, cv_main_100$elpd_pointwise, cv_main_300$elpd_pointwise, cv_main_300_s2$elpd_pointwise),
  se_pointwise_clustered = c(cv_main_20$se_pointwise_clustered, cv_main_100$se_pointwise_clustered, cv_main_300$se_pointwise_clustered, cv_main_300_s2$se_pointwise_clustered),
  se_pointwise_naive = c(cv_main_20$se_pointwise, cv_main_100$se_pointwise, cv_main_300$se_pointwise, cv_main_300_s2$se_pointwise),
  min_med_ess_m = c(
    min(cv_main_20$fold_diagnostics$min_mc_ess_m, na.rm = TRUE),
    min(cv_main_100$fold_diagnostics$min_mc_ess_m, na.rm = TRUE),
    min(cv_main_300$fold_diagnostics$min_mc_ess_m, na.rm = TRUE),
    min(cv_main_300_s2$fold_diagnostics$min_mc_ess_m, na.rm = TRUE)
  ),
  min_draw_ess_m = c(
    min(cv_main_20$group_scores$mc_ess_m_min, na.rm = TRUE),
    min(cv_main_100$group_scores$mc_ess_m_min, na.rm = TRUE),
    min(cv_main_300$group_scores$mc_ess_m_min, na.rm = TRUE),
    min(cv_main_300_s2$group_scores$mc_ess_m_min, na.rm = TRUE)
  )
)
cat("\nM & seed sensitivity summary (Main-Effects RSM, identical cached folds):\n")
print(m_sens, row.names = FALSE)

# -----------------------------------------------------------------------------
# 3. Step 3: Same 5-fold person split on Gaussian & Student-t Bias RSMs
# -----------------------------------------------------------------------------
cat("\n=== STEP 3: 5-fold person CV (Gaussian & Student-t Bias RSMs) ===\n")
cv_iso_300 <- cv_bmfrm(
  fit_iso, target = "person", folds = shared_folds, M = 300L, seed = 20261124L,
  iter = 2500L, warmup = 1000L, adapt_delta = 0.98
)
print(cv_iso_300)

cv_st_300 <- cv_bmfrm(
  fit_st, target = "person", folds = shared_folds, M = 300L, seed = 20261124L,
  iter = 2500L, warmup = 1000L, adapt_delta = 0.995
)
print(cv_st_300)

cat("\n--- Strict converged comparison (allow_unconverged = FALSE) ---\n")
if (isTRUE(cv_main_300$all_folds_converged) && isTRUE(cv_iso_300$all_folds_converged)) {
  cmp_strict_grp <- compare_bmfrm(main = cv_main_300, gaussian = cv_iso_300, score = "group")
  cmp_strict_pt  <- compare_bmfrm(main = cv_main_300, gaussian = cv_iso_300, score = "pointwise")
  print(cmp_strict_grp)
  print(cmp_strict_pt)
}

strict_3_err <- tryCatch(
  {
    cmp_all3_strict <- compare_bmfrm(
      main = cv_main_300, student_t = cv_st_300, gaussian = cv_iso_300,
      score = "group"
    )
    print(cmp_all3_strict)
    NULL
  },
  error = function(e) conditionMessage(e)
)
if (!is.null(strict_3_err)) {
  cat("\nStrict 3-model guard triggered (as intended):\n ", strict_3_err, "\n")
  cat("\n--- Diagnostic / sensitivity 3-model table (allow_unconverged = TRUE) ---\n")
  cat("WARNING: Do NOT treat an unconverged model's tiny lead as a valid win;\n")
  cat("resolve fold diagnostics first or report the caveat explicitly.\n\n")
  cmp_3_grp <- compare_bmfrm(
    main = cv_main_300, student_t = cv_st_300, gaussian = cv_iso_300,
    score = "group", allow_unconverged = TRUE
  )
  cmp_3_pt <- compare_bmfrm(
    main = cv_main_300, student_t = cv_st_300, gaussian = cv_iso_300,
    score = "pointwise", allow_unconverged = TRUE
  )
  print(cmp_3_grp)
  print(cmp_3_pt)
} else {
  cmp_3_grp <- cmp_all3_strict
  cmp_3_pt <- compare_bmfrm(
    main = cv_main_300, student_t = cv_st_300, gaussian = cv_iso_300,
    score = "pointwise"
  )
}

# -----------------------------------------------------------------------------
# 4. Practical assessment impact check (Main vs. Gaussian fair_scores)
# -----------------------------------------------------------------------------
cat("\n=== STEP 4: PRACTICAL IMPACT CHECK (per-unit ELPD & fair_scores) ===\n")
fs_main <- fair_scores(fit_main, standardize = TRUE)
fs_iso  <- fair_scores(fit_iso,  standardize = TRUE)
m_fs <- match(fs_main$person, fs_iso$person)
fs_diff <- fs_iso$mean_fair_score[m_fs] - fs_main$mean_fair_score
practical_impact <- list(
  cv_group_nat_per_person = (cv_iso_300$elpd_group - cv_main_300$elpd_group) / nrow(cv_main_300$group_scores),
  cv_pointwise_nat_per_rating = (cv_iso_300$elpd_pointwise - cv_main_300$elpd_pointwise) / nrow(ratings),
  loo_hybrid_nat_per_rating = if (!is.null(loo_exact_corrected)) loo_exact_corrected$per_rating_gain_gaussian_vs_main_nat else NA_real_,
  fs_rms_diff = sqrt(mean(fs_diff^2)),
  fs_median_abs_diff = stats::median(abs(fs_diff)),
  fs_max_abs_diff = max(abs(fs_diff)),
  fs_pearson_r = stats::cor(fs_main$mean_fair_score, fs_iso$mean_fair_score[m_fs]),
  fs_spearman_rho = stats::cor(fs_main$mean_fair_score, fs_iso$mean_fair_score[m_fs], method = "spearman")
)
cat(sprintf(
  "Per-unit CV gain (Gaussian - Main): %+.4f nat/person | %+.5f nat/rating (geom prob ratio %.4f)\n",
  practical_impact$cv_group_nat_per_person,
  practical_impact$cv_pointwise_nat_per_rating,
  exp(practical_impact$cv_pointwise_nat_per_rating)
))
if (is.finite(practical_impact$loo_hybrid_nat_per_rating)) {
  cat(sprintf(
    "Per-unit hybrid LOO gain (Gaussian - Main): %+.5f nat/rating (geom prob ratio %.4f)\n",
    practical_impact$loo_hybrid_nat_per_rating,
    exp(practical_impact$loo_hybrid_nat_per_rating)
  ))
}
cat(sprintf(
  "Standardized fair_scores (1..4 scale): RMS diff = %.4f | median |diff| = %.4f | max |diff| = %.4f | Pearson r = %.5f | Spearman rho = %.5f\n",
  practical_impact$fs_rms_diff, practical_impact$fs_median_abs_diff, practical_impact$fs_max_abs_diff,
  practical_impact$fs_pearson_r, practical_impact$fs_spearman_rho
))

# -----------------------------------------------------------------------------
# 5. Save reproducible summary + session provenance
# -----------------------------------------------------------------------------
out_obj <- list(
  provenance = list(
    r_version = R.version.string,
    pkg_load_mode = pkg_load_mode,
    pkg_load_path = pkg_load_info$path,
    bayesmfrm_version = as.character(utils::packageVersion("BayesMFRM")),
    cmdstanr_version = as.character(utils::packageVersion("cmdstanr")),
    cmdstan_version = as.character(cmdstanr::cmdstan_version()),
    loo_version = as.character(utils::packageVersion("loo")),
    sirt_version = as.character(utils::packageVersion("sirt")),
    folds_csv = folds_csv,
    script = "docs/compare-ratings1-cv.R"
  ),
  session_info = utils::sessionInfo(),
  loo = list(main = loo_main, student_t = loo_st, gaussian = loo_iso,
             comparison = cmp_loo, exact_refit_checks = exact_loo_checks,
             exact_corrected = loo_exact_corrected),
  folds = shared_folds,
  m_sensitivity_main = m_sens,
  cv = list(main = cv_main_300, student_t = cv_st_300, gaussian = cv_iso_300),
  comparisons = list(group = cmp_3_grp, pointwise = cmp_3_pt),
  practical_impact = practical_impact
)
saveRDS(out_obj, file = out_rds)
cat("\nSaved reproducible comparison results to:", out_rds, "\n")
