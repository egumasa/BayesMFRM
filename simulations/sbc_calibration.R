# Simulation-Based Calibration (SBC; Talts et al., 2018) for BayesMFRM
#
# Validates posterior calibration across four core hierarchical architectures:
#   1. Hier_RSM       : Hierarchical examinees + hierarchical raters (RSM)
#   2. Hier_PCM       : Hierarchical examinees + partial-pooled PCM steps (tau_common + tau_dev)
#   3. Bias_StudentT  : Hierarchical RSM + two-way interaction (Student-t local shrinkage)
#   4. Bias_Horseshoe : Hierarchical RSM + two-way interaction (regularized horseshoe)
#
# Usage (from package root):
#   Rscript --vanilla simulations/sbc_calibration.R [S=40] [out_dir=simulations/results/sbc]

args <- commandArgs(trailingOnly = TRUE)
S <- if (length(args) >= 1L) as.integer(args[1L]) else 40L
out_dir <- if (length(args) >= 2L) args[2L] else file.path("simulations", "results", "sbc")
stopifnot(!is.na(S), S >= 10L)

req <- c("pkgload", "cmdstanr", "posterior", "digest")
stopifnot(all(vapply(req, requireNamespace, logical(1), quietly = TRUE)))
suppressPackageStartupMessages(
  pkgload::load_all(".", quiet = TRUE, export_all = FALSE, export_imports = FALSE)
)
if (!check_cmdstan()) stop("CmdStan is required to run simulations/sbc_calibration.R")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
cache_tmp <- tempfile("bmfrm_sbc_cache_")
dir.create(cache_tmp, recursive = TRUE, showWarnings = FALSE)
on.exit(unlink(cache_tmp, recursive = TRUE), add = TRUE)

# Connected bridged design: P = 24 examinees, I = 4 items, R = 6 raters, K = 4 categories
# Each examinee receives 8 ratings (N = 192), ensuring full rank for main + interaction models.
set.seed(20261127)
P <- 24L; I <- 4L; R <- 6L; K <- 4L
full_grid <- expand.grid(
  person = sprintf("P%02d", seq_len(P)),
  item   = sprintf("I%02d", seq_len(I)),
  rater  = sprintf("R%02d", seq_len(R)),
  stringsAsFactors = FALSE
)
# Balanced incomplete block: every item x rater cell appears 8 times (24 * 8 / 24 = 8)
cell_grid <- expand.grid(item = sprintf("I%02d", seq_len(I)), rater = sprintf("R%02d", seq_len(R)), stringsAsFactors = FALSE)
plan <- data.frame(
  person = rep(sprintf("P%02d", seq_len(P)), each = 8L),
  item   = character(P * 8L),
  rater  = character(P * 8L),
  stringsAsFactors = FALSE
)
for (p in seq_len(P)) {
  # Coprime stride 5 (gcd(5, 24) = 1) bridges all persons while keeping exact 8 obs/cell
  idx <- (((p - 1L) * 5L + seq_len(8L) - 1L) %% nrow(cell_grid)) + 1L
  rows <- ((p - 1L) * 8L + 1L):(p * 8L)
  plan$item[rows] <- cell_grid$item[idx]
  plan$rater[rows] <- cell_grid$rater[idx]
}

stopifnot(
  check_design_bmfrm(score ~ person + item + rater, plan, K = K)$identifiable,
  check_design_bmfrm(score ~ person + item + rater + item:rater, plan, K = K)$identifiable
)

L_ranks <- 39L # Thin posterior to L = 39 draws -> discrete ranks in 0..39 (40 equal-prob bins, 5 bins merged = 8 bins)

sbc_configs <- list(
  list(
    id = "Hier_RSM",
    formula = score ~ person + item + rater,
    family = "rsm",
    threshold_facet = NULL,
    hierarchical = c("person", "rater"),
    priors = c(
      prior("normal(0, 0.8)", "mean", facet = "person"),
      prior("lognormal(-0.2, 0.25)", "sd", facet = "person"),
      prior("normal(0, 0.7)", "facet", facet = "item", scale = "marginal"),
      prior("lognormal(-0.5, 0.3)", "sd", facet = "rater", scale = "marginal"),
      prior("normal(0, 1.0)", "tau", scale = "marginal")
    ),
    targets = list(
      mu_theta    = function(p) p$mu_theta,
      sigma_theta = function(p) p$sigma_theta,
      sigma_rater = function(p) p$sigma_rater,
      `theta[1]`  = function(p) p$theta[1L],
      `item[1]`   = function(p) p$item[1L],
      `rater[1]`  = function(p) p$rater[1L],
      `tau[1]`    = function(p) p$tau[1L]
    ),
    adapt_delta = 0.95
  ),
  list(
    id = "Hier_PCM",
    formula = score ~ person + item + rater,
    family = "pcm",
    threshold_facet = "item",
    hierarchical = c("person", "tau"),
    priors = c(
      prior("normal(0, 0.8)", "mean", facet = "person"),
      prior("lognormal(-0.2, 0.25)", "sd", facet = "person"),
      prior("normal(0, 0.7)", "facet", facet = "item", scale = "marginal"),
      prior("normal(0, 0.6)", "facet", facet = "rater", scale = "marginal"),
      prior("normal(0, 0.9)", "tau", scale = "marginal"),
      prior("lognormal(-0.8, 0.3)", "sd", facet = "tau", scale = "marginal")
    ),
    targets = list(
      mu_theta        = function(p) p$mu_theta,
      sigma_theta     = function(p) p$sigma_theta,
      sigma_tau       = function(p) p$sigma_tau,
      `tau_common[1]` = function(p) p$tau_common[1L],
      `tau[1,1]`      = function(p) p$tau[1L, 1L],
      `theta[1]`      = function(p) p$theta[1L],
      `item[1]`       = function(p) p$item[1L]
    ),
    adapt_delta = 0.96
  ),
  list(
    id = "Bias_StudentT",
    formula = score ~ person + item + rater + item:rater,
    family = "rsm",
    threshold_facet = NULL,
    hierarchical = c("person", "item:rater"),
    priors = c(
      prior("normal(0, 0.8)", "mean", facet = "person"),
      prior("lognormal(-0.2, 0.25)", "sd", facet = "person"),
      prior("normal(0, 0.7)", "facet", facet = "item", scale = "marginal"),
      prior("normal(0, 0.6)", "facet", facet = "rater", scale = "marginal"),
      prior("normal(0, 1.0)", "tau", scale = "marginal"),
      prior("lognormal(-0.9, 0.3)", "sd", facet = "item:rater",
            scale = "marginal", local_shrinkage = "student_t")
    ),
    targets = list(
      mu_theta          = function(p) p$mu_theta,
      sigma_theta       = function(p) p$sigma_theta,
      sigma_itemrater   = function(p) p$sigma_itemrater,
      `itemrater[1,1]`  = function(p) p$itemrater[1L, 1L],
      `itemrater[2,3]`  = function(p) p$itemrater[2L, 3L],
      `theta[1]`        = function(p) p$theta[1L],
      `tau[1]`          = function(p) p$tau[1L]
    ),
    adapt_delta = 0.98
  ),
  list(
    id = "Bias_Horseshoe",
    formula = score ~ person + item + rater + item:rater,
    family = "rsm",
    threshold_facet = NULL,
    hierarchical = c("person", "item:rater"),
    priors = c(
      prior("normal(0, 0.8)", "mean", facet = "person"),
      prior("lognormal(-0.2, 0.25)", "sd", facet = "person"),
      prior("normal(0, 0.7)", "facet", facet = "item", scale = "marginal"),
      prior("normal(0, 0.6)", "facet", facet = "rater", scale = "marginal"),
      prior("normal(0, 1.0)", "tau", scale = "marginal"),
      prior("lognormal(-1.1, 0.35)", "sd", facet = "item:rater",
            scale = "raw", local_shrinkage = "horseshoe", slab_df = 4, slab_scale = 1.5)
    ),
    targets = list(
      mu_theta          = function(p) p$mu_theta,
      sigma_theta       = function(p) p$sigma_theta,
      sigma_itemrater   = function(p) p$sigma_itemrater,
      `itemrater[1,1]`  = function(p) p$itemrater[1L, 1L],
      `itemrater[2,3]`  = function(p) p$itemrater[2L, 3L],
      `theta[1]`        = function(p) p$theta[1L],
      `tau[1]`          = function(p) p$tau[1L]
    ),
    adapt_delta = 0.98
  )
)

cat(sprintf("=== BayesMFRM Simulation-Based Calibration (S = %d replications/model, L = %d ranks) ===\n", S, L_ranks))

all_rank_rows <- list()
all_diag_rows <- list()

for (cfg in sbc_configs) {
  cat(sprintf("\n[%s] Compiling Stan model...\n", cfg$id))
  init_sim <- simulate_bmfrm(
    cfg$formula, plan, K = K, priors = cfg$priors,
    hierarchical = cfg$hierarchical, family = cfg$family,
    threshold_facet = cfg$threshold_facet, seed = 9000L
  )
  stan_code <- build_stan_code(init_sim$spec, cfg$priors, save_gq = "minimal")
  comp <- compile_bmfrm(
    stan_code, model_name = paste0("sbc_", tolower(cfg$id)),
    cache_dir = cache_tmp, cpp_options = list(PRECOMPILED_HEADERS = "false")
  )
  stan_vars <- unique(sub("\\[.*$", "", names(cfg$targets)))

  for (s in seq_len(S)) {
    sim_s <- simulate_bmfrm(
      cfg$formula, plan, K = K, priors = cfg$priors,
      hierarchical = cfg$hierarchical, family = cfg$family,
      threshold_facet = cfg$threshold_facet, seed = 100000L + match(cfg$id, vapply(sbc_configs, `[[`, "", "id")) * 1000L + s
    )
    di_s <- prepare_data_bmfrm(sim_s$spec, sim_s$data, cfg$priors)
    fit_s <- comp$model$sample(
      data = di_s$stan_data,
      seed = 200000L + s,
      chains = 2L, parallel_chains = 2L,
      iter_warmup = 500L, iter_sampling = 500L,
      adapt_delta = cfg$adapt_delta, max_treedepth = 12L,
      refresh = 0L, show_messages = FALSE, show_exceptions = FALSE
    )
    diag_s <- fit_s$diagnostic_summary(quiet = TRUE)
    sum_s <- fit_s$summary(stan_vars)
    draws_mat <- unclass(fit_s$draws(stan_vars, format = "matrix"))
    thin_idx <- round(seq(1, nrow(draws_mat), length.out = L_ranks))

    p_true <- sim_s$parameters[[1L]]
    for (p_nm in names(cfg$targets)) {
      truth_val <- cfg$targets[[p_nm]](p_true)
      post_full <- draws_mat[, p_nm]
      post_thin <- post_full[thin_idx]
      # Tie-broken rank in 0..L_ranks
      less <- sum(post_thin < truth_val)
      ties <- sum(post_thin == truth_val)
      rk <- less + if (ties > 0L) sample.int(ties + 1L, 1L) - 1L else 0L
      q_ci <- stats::quantile(post_full, probs = c(0.05, 0.25, 0.75, 0.95))

      all_rank_rows[[length(all_rank_rows) + 1L]] <- data.frame(
        model = cfg$id,
        replication = s,
        parameter = p_nm,
        truth = truth_val,
        post_mean = mean(post_full),
        post_sd = stats::sd(post_full),
        rank = rk,
        u = (rk + 0.5) / (L_ranks + 1L),
        in_50 = (truth_val >= q_ci[2L] && truth_val <= q_ci[3L]),
        in_90 = (truth_val >= q_ci[1L] && truth_val <= q_ci[4L]),
        stringsAsFactors = FALSE
      )
    }

    all_diag_rows[[length(all_diag_rows) + 1L]] <- data.frame(
      model = cfg$id,
      replication = s,
      max_rhat = max(sum_s$rhat, na.rm = TRUE),
      min_ess_bulk = min(sum_s$ess_bulk, na.rm = TRUE),
      n_divergent = sum(diag_s$num_divergent),
      n_max_treedepth = sum(diag_s$num_max_treedepth),
      stringsAsFactors = FALSE
    )
    unlink(fit_s$output_files())
    if (s %% 10L == 0L || s == S) {
      cat(sprintf("  [%s] Completed %d / %d replications\n", cfg$id, s, S))
    }
  }
}

ranks_df <- do.call(rbind, all_rank_rows)
diag_df <- do.call(rbind, all_diag_rows)

# Summarize calibration per (model, parameter):
# 8 coarse bins (5 adjacent ranks per bin out of 40) -> exact Pearson chi-square uniformity p-value
summary_list <- lapply(split(ranks_df, paste(ranks_df$model, ranks_df$parameter, sep = "||")), function(sub) {
  bins8 <- tabulate(floor(sub$rank / 5L) + 1L, nbins = 8L)
  chisq_p <- stats::chisq.test(bins8, p = rep(1 / 8, 8L))$p.value
  # Empirical ECDF max deviation vs Uniform(0, 1)
  grid_z <- (1:L_ranks) / (L_ranks + 1L)
  ecdf_diff <- vapply(grid_z, function(z) mean(sub$rank < z * (L_ranks + 1L)) - z, numeric(1))
  prior_var <- stats::var(sub$truth)
  contraction <- if (prior_var > 0) 1 - mean(sub$post_sd^2) / prior_var else NA_real_
  d_sub <- diag_df[diag_df$model == sub$model[1L], ]
  data.frame(
    model = sub$model[1L],
    parameter = sub$parameter[1L],
    S = nrow(sub),
    mean_norm_rank = mean(sub$u),
    sd_norm_rank = stats::sd(sub$u),
    max_ecdf_dev = max(abs(ecdf_diff)),
    chisq_8bin_p = chisq_p,
    coverage_50 = mean(sub$in_50),
    coverage_90 = mean(sub$in_90),
    contraction = contraction,
    div_free_rate = mean(d_sub$n_divergent == 0L),
    rhat_ok_rate = mean(d_sub$max_rhat < 1.02),
    stringsAsFactors = FALSE
  )
})
sbc_summary <- do.call(rbind, summary_list)
rownames(sbc_summary) <- NULL

write.csv(ranks_df, file.path(out_dir, "sbc_ranks.csv"), row.names = FALSE)
write.csv(diag_df, file.path(out_dir, "sbc_diagnostics.csv"), row.names = FALSE)
write.csv(sbc_summary, file.path(out_dir, "sbc_summary.csv"), row.names = FALSE)

# Plot ECDF difference envelopes across all 4 models
png(file.path(out_dir, "sbc_ecdf_ranks.png"), width = 1600, height = 1000, res = 150)
par(mfrow = c(2, 2), mar = c(4.2, 4.2, 3.0, 1.0))
grid_z <- (1:L_ranks) / (L_ranks + 1L)
# Pointwise 99% simultaneous-approximating binomial envelope (+/- 2.6 * sqrt(z*(1-z)/S))
band <- 2.6 * sqrt(grid_z * (1 - grid_z) / S)
cols <- c("#1b6ca8", "#d95f02", "#1b9e77", "#7570b3", "#e7298a", "#66a61e", "#a6761d")

for (cfg in sbc_configs) {
  sub_m <- ranks_df[ranks_df$model == cfg$id, ]
  d_m <- diag_df[diag_df$model == cfg$id, ]
  plot(grid_z, rep(0, length(grid_z)), type = "n", ylim = c(-max(band) * 1.4, max(band) * 1.4),
       xlab = "Fractional prior rank z", ylab = "ECDF(z) - z",
       main = sprintf("%s (S=%d, div-free=%.0f%%)", cfg$id, S, 100 * mean(d_m$n_divergent == 0L)))
  polygon(c(grid_z, rev(grid_z)), c(band, rev(-band)), col = "#eaeff5", border = NA)
  abline(h = 0, lty = 2, col = "gray40")
  p_names <- unique(sub_m$parameter)
  for (j in seq_along(p_names)) {
    r_j <- sub_m$rank[sub_m$parameter == p_names[j]]
    ecdf_d <- vapply(grid_z, function(z) mean(r_j < z * (L_ranks + 1L)) - z, numeric(1))
    lines(grid_z, ecdf_d, col = cols[((j - 1L) %% length(cols)) + 1L], lwd = 1.8)
  }
  legend("topleft", legend = p_names, col = cols[seq_along(p_names)], lwd = 1.6, cex = 0.68, bty = "n", ncol = 2)
}
dev.off()

cat("\n=== SBC Calibration Summary ===\n")
print(sbc_summary, digits = 3, row.names = FALSE)
cat(sprintf("\nSaved SBC tables and plot to: %s\n", out_dir))
