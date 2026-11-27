# Reproducible 3-model Bayesian MFRM analysis of sirt::data.ratings1
# Generates the exact JSON summary embedded in docs/index.html.
#
# Usage:
#   From the repository root: Rscript --vanilla docs/reproduce-ratings1.R [cache_dir] [output_json]

args <- commandArgs(trailingOnly = TRUE)
cache_dir <- if (length(args) >= 1L) args[1L] else file.path(tempdir(), "stan_cache_ratings1")
out_json  <- if (length(args) >= 2L) args[2L] else file.path("docs", "ratings1_results.json")

# Match the implementation and templates in this checkout, not an older installed build.
suppressPackageStartupMessages(pkgload::load_all(".", quiet = TRUE,
                                                export_all = FALSE, export_imports = FALSE))
suppressPackageStartupMessages({
  library(sirt)
  library(tidyr)
  library(jsonlite)
})

dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_json), recursive = TRUE, showWarnings = FALSE)

# 1. Load empirical writing ratings (135 students x 5 criteria x 7 raters)
#    Raw scores are 0..3 -> shift +1 to 1..4 and specify K = 4 explicitly.
data(data.ratings1)
ratings <- tidyr::pivot_longer(data.ratings1, cols = k1:k5, names_to = "item", values_to = "score_raw")
ratings <- ratings[!is.na(ratings$score_raw), ]
ratings$person <- as.character(ratings$idstud)
ratings$rater  <- as.character(ratings$rater)
ratings$score  <- as.integer(ratings$score_raw + 1L)
K <- 4L

# 2. Structural design identifiability checks BEFORE sampling
dc_main <- check_design_bmfrm(score ~ person + item + rater, ratings, K = K)
dc_bias <- check_design_bmfrm(score ~ person + item + rater + item:rater, ratings, K = K)
stopifnot(dc_main$identifiable, dc_bias$identifiable)

# 3. Prior specifications (3-model comparison: Main RSM, Student-t Bias, Gaussian Bias)
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

# Real prior predictive check for the EXACT fitted priors_main vs. an overly diffuse baseline
priors_diffuse <- c(
  prior("normal(0, 3.0)", "mean", facet = "person"),
  prior("lognormal(1.0, 0.5)", "sd", facet = "person"),
  prior("normal(0, 3.0)", "item", scale = "marginal"),
  prior("normal(0, 3.0)", "sd", facet = "rater", scale = "marginal"),
  prior("normal(0, 3.0)", "tau", scale = "marginal")
)
summarise_prior_ppc <- function(priors_obj, seed) {
  pp <- prior_predict_bmfrm(
    score ~ person + item + rater, design = ratings, K = K,
    hierarchical = c("person", "rater"), priors = priors_obj,
    ndraws = 300L, seed = seed
  )
  props <- t(apply(pp$y_rep, 1L, function(y) as.numeric(prop.table(table(factor(y, 1:K))))))
  ext_draws <- props[, 1L] + props[, K]
  list(
    cat = seq_len(K),
    mean = as.numeric(colMeans(props)),
    lo95 = as.numeric(apply(props, 2L, stats::quantile, 0.025)),
    hi95 = as.numeric(apply(props, 2L, stats::quantile, 0.975)),
    extreme_mean = mean(ext_draws),
    extreme_lo95 = unname(stats::quantile(ext_draws, 0.025)),
    extreme_hi95 = unname(stats::quantile(ext_draws, 0.975))
  )
}
real_prior_ppc <- list(
  fitted_main = summarise_prior_ppc(priors_main, seed = 42L),
  diffuse_sd3 = summarise_prior_ppc(priors_diffuse, seed = 43L)
)

# 4. Fit all 3 models with CmdStan (save_gq = "minimal" keeps CSVs lightweight)
fit_main <- bmfrm(
  score ~ person + item + rater, data = ratings, K = K,
  hierarchical = c("person", "rater"), priors = priors_main,
  save_gq = "minimal", cache_dir = cache_dir, model_name = "ratings1_main",
  design_check = "error", chains = 4, cores = 4, iter = 1500, warmup = 500,
  adapt_delta = 0.95, max_treedepth = 12, seed = 202611,
  cpp_options = list(PRECOMPILED_HEADERS = "false")
)
fit_st <- bmfrm(
  score ~ person + item + rater + item:rater, data = ratings, K = K,
  hierarchical = c("person", "rater", "item:rater"), priors = priors_st,
  save_gq = "minimal", cache_dir = cache_dir, model_name = "ratings1_bias_st",
  design_check = "error", chains = 4, cores = 4, iter = 2000, warmup = 1000,
  adapt_delta = 0.99, max_treedepth = 12, seed = 202612,
  cpp_options = list(PRECOMPILED_HEADERS = "false")
)
fit_iso <- bmfrm(
  score ~ person + item + rater + item:rater, data = ratings, K = K,
  hierarchical = c("person", "rater", "item:rater"), priors = priors_iso,
  save_gq = "minimal", cache_dir = cache_dir, model_name = "ratings1_bias_iso",
  design_check = "error", chains = 4, cores = 4, iter = 1500, warmup = 500,
  adapt_delta = 0.95, max_treedepth = 12, seed = 202613,
  cpp_options = list(PRECOMPILED_HEADERS = "false")
)

# 5. Identical diagnostic battery across all 3 models (PPC + facet_fit + residual_cor)
extract_model_diagnostics <- function(fit_obj, label_str) {
  rt <- facet_summary(fit_obj, "rater")
  rf <- facet_fit(fit_obj, "rater", include_zstd = TRUE)
  it <- facet_summary(fit_obj, "item")
  if_fit <- facet_fit(fit_obj, "item", include_zstd = TRUE)
  res <- residuals(fit_obj)
  rc <- residual_cor_criteria(res, unit_facets = c("person", "rater"))
  cor_mat <- as.matrix(rc$cor_matrix)
  off_diag <- abs(cor_mat[row(cor_mat) != col(cor_mat)])

  # Conditional PPC (200 posterior predictive draws)
  q_ppc <- BayesMFRM:::compute_bmfrm_quantities(
    fit_obj, draw_idx = seq(1, 4000, by = 20), compute_y_rep = TRUE, seed = 42
  )
  obs_prop <- as.numeric(prop.table(table(factor(ratings$score, 1:K))))
  rep_props <- t(apply(q_ppc$y_rep, 1L, function(y) as.numeric(prop.table(table(factor(y, 1:K))))))
  lo_ov <- apply(rep_props, 2L, stats::quantile, 0.025)
  hi_ov <- apply(rep_props, 2L, stats::quantile, 0.975)
  ppc_overall <- data.frame(
    cat = 1:K, obs = obs_prop,
    pred_mean = as.numeric(colMeans(rep_props)),
    pred_lo = as.numeric(lo_ov),
    pred_hi = as.numeric(hi_ov),
    in_band = obs_prop >= lo_ov & obs_prop <= hi_ov
  )

  raters <- fit_obj$data_info$facet_labels$rater
  ppc_by_rater <- lapply(raters, function(r) {
    idx <- which(ratings$rater == r)
    op <- as.numeric(prop.table(table(factor(ratings$score[idx], 1:K))))
    rp <- t(apply(q_ppc$y_rep[, idx, drop = FALSE], 1L, function(y) as.numeric(prop.table(table(factor(y, 1:K))))))
    lo_r <- as.numeric(apply(rp, 2L, stats::quantile, 0.025))
    hi_r <- as.numeric(apply(rp, 2L, stats::quantile, 0.975))
    in_b <- op >= lo_r & op <= hi_r
    list(
      rater = r, n = length(idx), obs = op,
      pred = as.numeric(colMeans(rp)), lo = lo_r, hi = hi_r,
      in_band = in_b, n_in_band = sum(in_b)
    )
  })
  total_rater_cells_in_band <- sum(vapply(ppc_by_rater, function(x) x$n_in_band, integer(1)))
  db08_ppc_in_band <- ppc_by_rater[[match("db08", raters)]]$n_in_band

  rater_tbl <- as.data.frame(merge(
    rt,
    rf[, c("label", "n_obs", "infit_ppp", "outfit_ppp", "infit_msq", "infit_lo", "infit_hi", "infit_zstd",
           "outfit_msq", "outfit_lo", "outfit_hi", "outfit_zstd")],
    by.x = "level", by.y = "label", sort = FALSE
  ))
  item_tbl <- as.data.frame(merge(
    it,
    if_fit[, c("label", "n_obs", "infit_ppp", "outfit_ppp", "infit_msq", "infit_lo", "infit_hi", "infit_zstd",
               "outfit_msq", "outfit_lo", "outfit_hi", "outfit_zstd")],
    by.x = "level", by.y = "label", sort = FALSE
  ))
  db08_row <- rater_tbl[rater_tbl$level == "db08", , drop = FALSE]
  db07_row <- rater_tbl[rater_tbl$level == "db07", , drop = FALSE]

  list(
    label = label_str,
    conv = list(
      converged = fit_obj$convergence$converged,
      max_rhat = fit_obj$convergence$max_rhat,
      min_ess_bulk = fit_obj$convergence$min_ess_bulk,
      min_ess_tail = fit_obj$convergence$min_ess_tail,
      n_divergent = fit_obj$convergence$n_divergent,
      min_ebfmi = fit_obj$convergence$min_ebfmi
    ),
    summary_metrics = list(
      overall_ppc_in_band = sum(ppc_overall$in_band),
      rater_ppc_in_band = total_rater_cells_in_band,
      rater_ppc_total = length(raters) * K,
      db08_ppc_in_band = db08_ppc_in_band,
      db08_infit_msq = db08_row$infit_msq[1L],
      db08_infit_zstd = db08_row$infit_zstd[1L],
      db08_outfit_msq = db08_row$outfit_msq[1L],
      db07_infit_msq = db07_row$infit_msq[1L],
      db07_infit_zstd = db07_row$infit_zstd[1L],
      max_abs_resid_cor = max(off_diag, na.rm = TRUE)
    ),
    raters = rater_tbl,
    items = item_tbl,
    ppc_overall = ppc_overall,
    ppc_by_rater = ppc_by_rater,
    resid_cor = cor_mat
  )
}

diag_main <- extract_model_diagnostics(fit_main, "Main-Effects RSM")
diag_st   <- extract_model_diagnostics(fit_st,   "Bias RSM (Student-t)")
diag_iso  <- extract_model_diagnostics(fit_iso,  "Bias RSM (Isotropic Gaussian)")

# 6. Decompose Fair - Raw into:
#    (a) Pure Rater Adjustment = E[mu_fair - mu_actual] on the SAME posterior draws & criteria
#    (b) Shrinkage & Sampling Noise = E[mu_actual] - observed_raw_mean
fs_main <- fair_scores(fit_main, adjust_facets = "rater", standardize = TRUE)
th_main <- facet_summary(fit_main, "person")
tau_main <- as.data.frame(fit_main$fit$summary("tau", "mean", "sd", ~stats::quantile(.x, c(0.025, 0.975))))

st_actual <- BayesMFRM:::stream_bmfrm_fair_scores(fit_main, newdata = ratings, adjust_facets = character(0), summary = TRUE)
st_fair   <- BayesMFRM:::stream_bmfrm_fair_scores(fit_main, newdata = ratings, adjust_facets = "rater", summary = TRUE)
mu_act_draws  <- sweep(st_actual$sum_mu, 2L, st_actual$count_p, "/")
mu_fair_draws <- sweep(st_fair$sum_mu,   2L, st_fair$count_p,   "/")
rater_adj_draws <- mu_fair_draws - mu_act_draws

p_labels <- st_actual$labels
p_decomp <- data.frame(
  person = p_labels,
  mean_actual_expected = as.numeric(colMeans(mu_act_draws)),
  rater_adjustment = as.numeric(colMeans(rater_adj_draws)),
  rater_adj_lo = as.numeric(apply(rater_adj_draws, 2L, stats::quantile, 0.025)),
  rater_adj_hi = as.numeric(apply(rater_adj_draws, 2L, stats::quantile, 0.975)),
  stringsAsFactors = FALSE
)
p_raters <- tapply(ratings$rater, ratings$person, function(x) paste(sort(unique(x)), collapse = ", "))
fs_df <- as.data.frame(merge(fs_main, th_main[, c("level", "mean", "q2.5", "q97.5")], by.x = "person", by.y = "level", sort = FALSE))
fs_df <- as.data.frame(merge(fs_df, p_decomp, by = "person", sort = FALSE))
fs_df$raters_assigned <- as.character(unname(p_raters[fs_df$person]))
fs_df$fair_minus_raw <- as.numeric(fs_df$mean_fair_score - fs_df$mean_score)
fs_df$shrinkage_and_noise <- as.numeric(fs_df$mean_actual_expected - fs_df$mean_score)

# 7. Bias grids across delta = 0.0 .. 0.5 for both Student-t and Gaussian models
extract_bias_grid <- function(fit_obj) {
  deltas <- c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5)
  base <- summarise_bias(fit_obj, "item:rater", prob_thresh = 0.95, delta = 0.3)
  d_mat <- unclass(fit_obj$fit$draws("itemrater", format = "matrix"))
  d_mat <- d_mat[, base$variable, drop = FALSE]
  prob_gt <- lapply(deltas, function(dl) unname(colMeans(d_mat > dl)))
  prob_lt <- lapply(deltas, function(dl) unname(colMeans(d_mat < -dl)))
  names(prob_gt) <- names(prob_lt) <- sprintf("%.1f", deltas)
  list(
    cells = as.data.frame(base[, c("item", "rater", "bias_mean", "bias_sd", "lower_95", "upper_95",
                                   "n_observations", "likelihood_identified")]),
    prob_gt = prob_gt,
    prob_lt = prob_lt
  )
}

out <- list(
  provenance = list(
    r_version = R.version.string,
    bayesmfrm_version = as.character(utils::packageVersion("BayesMFRM")),
    cmdstanr_version = as.character(utils::packageVersion("cmdstanr")),
    cmdstan_version = as.character(cmdstanr::cmdstan_version()),
    sirt_version = as.character(utils::packageVersion("sirt")),
    posterior_version = as.character(utils::packageVersion("posterior")),
    script = "docs/reproduce-ratings1.R"
  ),
  meta = list(
    dataset = "sirt::data.ratings1",
    N = nrow(ratings), P = length(unique(ratings$person)),
    I = length(unique(ratings$item)), R = length(unique(ratings$rater)), K = K,
    design_main = list(connected = dc_main$connected, identifiable = dc_main$identifiable,
                       rank = dc_main$adjacent_logit_rank, params = dc_main$adjacent_logit_parameters),
    design_bias = list(connected = dc_bias$connected, identifiable = dc_bias$identifiable,
                       rank = dc_bias$adjacent_logit_rank, params = dc_bias$adjacent_logit_parameters),
    conv_main = diag_main$conv,
    conv_st   = diag_st$conv,
    conv_iso  = diag_iso$conv
  ),
  real_prior_ppc = real_prior_ppc,
  models = list(
    main = diag_main,
    st   = diag_st,
    iso  = diag_iso
  ),
  raters = diag_main$raters,
  items  = diag_main$items,
  tau    = tau_main,
  persons = fs_df,
  bias_st = extract_bias_grid(fit_st),
  bias_iso = extract_bias_grid(fit_iso)
)

jsonlite::write_json(out, out_json, auto_unbox = TRUE, digits = 4)

cat("\n=== 3-Model Diagnostic Comparison on sirt::data.ratings1 ===\n")
for (m in list(diag_main, diag_st, diag_iso)) {
  sm <- m$summary_metrics
  cat(sprintf("%-28s | Rhat=%.3f ESS=%4.0f Div=%d | Overall PPC=%d/4 | Rater PPC=%2d/%d (db08=%d/4) | db08 Infit=%.2f (z=%+.2f) | db07 Infit=%.2f (z=%+.2f) | Max|r|=%.3f\n",
              m$label, m$conv$max_rhat, m$conv$min_ess_bulk, m$conv$n_divergent,
              sm$overall_ppc_in_band, sm$rater_ppc_in_band, sm$rater_ppc_total, sm$db08_ppc_in_band,
              sm$db08_infit_msq, sm$db08_infit_zstd, sm$db07_infit_msq, sm$db07_infit_zstd, sm$max_abs_resid_cor))
}
cat("\nReal prior predictive check (fitted priors_main vs diffuse SD=3):\n")
cat("  fitted_main:", round(real_prior_ppc$fitted_main$mean, 3), "extreme=", round(real_prior_ppc$fitted_main$extreme_mean, 3), "\n")
cat("  diffuse_sd3:", round(real_prior_ppc$diffuse_sd3$mean, 3), "extreme=", round(real_prior_ppc$diffuse_sd3$extreme_mean, 3), "\n")

db08_only <- fs_df[fs_df$raters_assigned == "db08", c("person", "mean_score", "mean_actual_expected", "mean_fair_score", "fair_minus_raw", "rater_adjustment", "shrinkage_and_noise")]
if (nrow(db08_only)) {
  cat("\ndb08-only examinees (pure rater adjustment vs shrinkage+noise):\n")
  print(db08_only, digits = 3, row.names = FALSE)
}
cat("\nSaved:", out_json, "\n")
