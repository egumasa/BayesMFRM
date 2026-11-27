# Empirical Benchmark: Joint ML vs. Gauss-Hermite MML-EAP vs. BayesMFRM
# and the Sparse Bias Shrinkage Trade-Off (Null-Cell Noise vs. Active-Cell Attenuation)
#
# Reproduce from the package root:
#   Rscript --vanilla simulations/bayesian_benefits.R [replications=20]
#
# Evaluates two honest empirical questions:
#   1. Ability & Held-Out Prediction (P = 40, n = 3 vs. n = 12 ratings/person):
#      Compares unpenalized Joint ML (diverges on extremes), Marginal ML + EAP
#      (31-node Gauss-Hermite quadrature; strong frequentist empirical-Bayes baseline),
#      and BayesMFRM hierarchical RSM (full parameter-uncertainty propagation).
#   2. Sparse Rater x Item Interaction Shrinkage Trade-Off (4 x 4 = 16 cells, 2 active biased cells):
#      Demonstrates why hierarchical shrinkage reduces overall & null-cell RMSE and
#      false discoveries while attenuating truly biased active cells toward zero.

args <- commandArgs(trailingOnly = TRUE)
B <- if (length(args)) as.integer(args[1L]) else 20L
stopifnot(!is.na(B), B > 0L)

required <- c("pkgload", "cmdstanr", "posterior", "digest")
stopifnot(all(vapply(required, requireNamespace, logical(1), quietly = TRUE)))
suppressPackageStartupMessages(
  pkgload::load_all(".", quiet = TRUE, export_all = FALSE, export_imports = FALSE)
)
if (!check_cmdstan()) stop("CmdStan is required to run simulations/bayesian_benefits.R")

out <- file.path("simulations", "results")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
scratch <- tempfile("bayesmfrm_sim_")
dir.create(scratch, recursive = TRUE, showWarnings = FALSE)
on.exit(unlink(scratch, recursive = TRUE), add = TRUE)

P <- 40L; I <- 4L; R <- 4L; K <- 6L
truth_item  <- c(-0.9, -0.3, 0.3, 0.9)
truth_rater <- c(-0.8, -0.2, 0.2, 0.8)
truth_tau   <- seq(-1.5, 1.5, length.out = K - 1L)
counts      <- rep(c(3L, 12L), each = P / 2L)

rsm_prob <- function(eta, tau) {
  logits <- sweep(outer(eta, 0:(K - 1L)), 2L, c(0, cumsum(tau)), "-")
  logits <- logits - apply(logits, 1L, max)
  z <- exp(logits)
  z / rowSums(z)
}

simulate_data <- function(seed) {
  set.seed(seed)
  theta <- stats::rnorm(P, 0, 1)
  d <- expand.grid(person = seq_len(P), item = seq_len(I), rater = seq_len(R))
  probs <- rsm_prob(theta[d$person] - truth_item[d$item] - truth_rater[d$rater], truth_tau)
  d$score <- vapply(seq_len(nrow(d)), function(j) sample.int(K, 1L, prob = probs[j, ]), integer(1))
  selected <- unlist(lapply(seq_len(P), function(p) sample(which(d$person == p), counts[p])))
  d$person <- sprintf("P%02d", d$person)
  d$item   <- sprintf("I%02d", d$item)
  d$rater  <- sprintf("R%02d", d$rater)
  list(train = d[selected, ], test = d[-selected, ], theta = theta)
}

# 1. Joint ML (unpenalized; infinite on all-min / all-max patterns)
unpack_jml <- function(x) {
  item  <- x[P + seq_len(I - 1L)]
  rater <- x[P + I - 1L + seq_len(R - 1L)]
  tau   <- x[P + I + R - 2L + seq_len(K - 2L)]
  list(theta = x[seq_len(P)], item = c(item, -sum(item)),
       rater = c(rater, -sum(rater)), tau = c(tau, -sum(tau)))
}

make_jml_objective <- function(d) {
  pi <- as.integer(factor(d$person)); ii <- as.integer(factor(d$item)); ri <- as.integer(factor(d$rater))
  indicators <- outer(d$score, seq_len(K - 1L), ">")
  group_sum <- function(x, id, n) as.numeric(rowsum(x, factor(id, levels = seq_len(n)), reorder = TRUE))
  function(x, gradient = FALSE) {
    a <- unpack_jml(x)
    eta <- a$theta[pi] - a$item[ii] - a$rater[ri]
    logits <- sweep(outer(eta, 0:(K - 1L)), 2L, c(0, cumsum(a$tau)), "-")
    shift <- apply(logits, 1L, max)
    z <- exp(logits - shift)
    if (!gradient) return(sum(log(rowSums(z)) + shift - logits[cbind(seq_len(nrow(d)), d$score)]))
    prob <- z / rowSums(z)
    g_eta <- as.vector(prob %*% seq_len(K)) - d$score
    g_theta <- group_sum(g_eta, pi, P)
    g_item  <- -group_sum(g_eta, ii, I)
    g_rater <- -group_sum(g_eta, ri, R)
    tails <- vapply(seq_len(K - 1L), function(k) rowSums(prob[, (k + 1L):K, drop = FALSE]), numeric(nrow(d)))
    g_tau <- colSums(indicators - tails)
    c(g_theta, g_item[-I] - g_item[I], g_rater[-R] - g_rater[R], g_tau[-(K - 1L)] - g_tau[K - 1L])
  }
}

# 2. Marginal ML + EAP via 31-node Gauss-Hermite quadrature (exact Golub-Welsch rule)
gauss_hermite_std_normal <- function(Q = 31L) {
  off <- sqrt(seq_len(Q - 1L))
  J_mat <- matrix(0, Q, Q)
  J_mat[cbind(seq_len(Q - 1L), 2:Q)] <- off
  J_mat[cbind(2:Q, seq_len(Q - 1L))] <- off
  eig <- eigen(J_mat, symmetric = TRUE)
  ord <- order(eig$values)
  list(nodes = eig$values[ord], weights = eig$vectors[1L, ord]^2)
}
gh_rule <- gauss_hermite_std_normal(31L)

fit_mml_eap <- function(d) {
  pi <- as.integer(factor(d$person, levels = sprintf("P%02d", seq_len(P))))
  ii <- as.integer(factor(d$item,   levels = sprintf("I%02d", seq_len(I))))
  ri <- as.integer(factor(d$rater,  levels = sprintf("R%02d", seq_len(R))))
  y  <- d$score
  Q  <- length(gh_rule$nodes)
  log_w <- log(gh_rule$weights)

  unpack_mml <- function(par) {
    mu_th <- par[1L]
    sd_th <- exp(par[2L])
    it <- par[2L + seq_len(I - 1L)]
    rt <- par[2L + (I - 1L) + seq_len(R - 1L)]
    tu <- par[2L + (I - 1L) + (R - 1L) + seq_len(K - 2L)]
    list(mu = mu_th, sd = sd_th,
         item = c(it, -sum(it)), rater = c(rt, -sum(rt)), tau = c(tu, -sum(tu)))
  }

  eval_person_grid <- function(a) {
    th_q <- a$mu + a$sd * gh_rule$nodes # length Q
    cum_tau <- c(0, cumsum(a$tau))
    offset_n <- -a$item[ii] - a$rater[ri]
    # N x Q log-likelihood matrix
    ll_nq <- matrix(0, nrow(d), Q)
    for (q in seq_len(Q)) {
      eta_q <- th_q[q] + offset_n
      logits_q <- sweep(outer(eta_q, 0:(K - 1L)), 2L, cum_tau, "-")
      m_q <- apply(logits_q, 1L, max)
      lse_q <- m_q + log(rowSums(exp(logits_q - m_q)))
      ll_nq[, q] <- logits_q[cbind(seq_len(nrow(d)), y)] - lse_q
    }
    # Sum within person -> P x Q
    ll_pq <- rowsum(ll_nq, factor(pi, levels = seq_len(P)), reorder = TRUE)
    log_post_pq <- sweep(ll_pq, 2L, log_w, "+")
    m_p <- apply(log_post_pq, 1L, max)
    log_marg_p <- m_p + log(rowSums(exp(log_post_pq - m_p)))
    list(th_q = th_q, log_post_pq = log_post_pq, log_marg_p = log_marg_p)
  }

  nll <- function(par) {
    a <- unpack_mml(par)
    -sum(eval_person_grid(a)$log_marg_p)
  }

  init_par <- c(0, 0, rep(0, (I - 1L) + (R - 1L) + (K - 2L)))
  opt <- stats::optim(init_par, nll, method = "BFGS", control = list(maxit = 1000L, reltol = 1e-10))
  a_hat <- unpack_mml(opt$par)
  eg <- eval_person_grid(a_hat)
  w_pq <- exp(sweep(eg$log_post_pq, 1L, eg$log_marg_p, "-"))
  eap_mean <- as.vector(w_pq %*% eg$th_q)
  eap_var  <- pmax(as.vector(w_pq %*% (eg$th_q^2)) - eap_mean^2, 1e-12)
  eap_sd   <- sqrt(eap_var)
  list(par = a_hat, eap = eap_mean, eap_sd = eap_sd,
       lower = eap_mean - 1.96 * eap_sd, upper = eap_mean + 1.96 * eap_sd,
       w_pq = w_pq, th_q = eg$th_q, convergence = opt$convergence)
}

# 3. Compile BayesMFRM hierarchical RSM (new package default: hierarchical = "person")
init_data <- simulate_data(20261124L)
spec <- parse_bmfrm_formula(score ~ person + item + rater, init_data$train, K = K, model_name = "simulation_rsm")
priors <- get_default_priors(spec$facets_main, spec$facets_bias, hierarchical = "person")
stan_code <- build_stan_code(spec, priors, save_gq = "minimal")
writeLines(stan_code, file.path(out, "simulation_rsm.stan"))
comp_rsm <- compile_bmfrm(stan_code, "simulation_rsm", cache_dir = scratch,
                          cpp_options = list(PRECOMPILED_HEADERS = "false"))

people <- vector("list", B); diagnostics <- vector("list", B); predictions <- vector("list", B)

cat("=== Part 1: Joint ML vs. MML-EAP vs. BayesMFRM (Ability & Prediction) ===\n")
for (b in seq_len(B)) {
  cat(sprintf("  Part 1 Replication %d / %d\n", b, B))
  dat <- simulate_data(20261123L + b)
  prepared <- prepare_data_bmfrm(spec, dat$train, priors)
  fit <- comp_rsm$model$sample(
    data = prepared$stan_data, seed = 74000L + b,
    chains = 4L, parallel_chains = 4L,
    iter_warmup = 750L, iter_sampling = 1000L,
    adapt_delta = 0.95, max_treedepth = 12L,
    refresh = 0L, show_messages = FALSE
  )
  vars <- c("theta", "item", "rater", "tau", "mu_theta", "sigma_theta")
  s <- fit$summary(vars)
  sampler <- fit$diagnostic_summary(quiet = TRUE)
  draws <- unclass(fit$draws(vars, format = "matrix"))
  theta_draws <- draws[, sprintf("theta[%d]", seq_len(P)), drop = FALSE]
  item_draws  <- draws[, sprintf("item[%d]", seq_len(I)), drop = FALSE]
  rater_draws <- draws[, sprintf("rater[%d]", seq_len(R)), drop = FALSE]
  tau_draws   <- draws[, sprintf("tau[%d]", seq_len(K - 1L)), drop = FALSE]

  # Joint ML
  obj_jml <- make_jml_objective(dat$train)
  start_jml <- c(colMeans(theta_draws), colMeans(item_draws)[-I],
                 colMeans(rater_draws)[-R], colMeans(tau_draws)[-(K - 1L)])
  jml <- stats::optim(start_jml, obj_jml, function(x) obj_jml(x, TRUE), method = "BFGS",
                      control = list(maxit = 2000L, reltol = 1e-10))
  jml_par <- unpack_jml(jml$par)
  extreme <- vapply(split(dat$train$score, dat$train$person),
                    function(y) all(y == 1L) || all(y == K), logical(1))

  # Gauss-Hermite MML + EAP
  mml <- fit_mml_eap(dat$train)

  ci_bayes <- apply(theta_draws, 2L, stats::quantile, probs = c(0.025, 0.975))
  people[[b]] <- data.frame(
    replication = b, person = seq_len(P), n = counts,
    truth = dat$theta,
    bayes = colMeans(theta_draws),
    mml_eap = mml$eap,
    jml = ifelse(extreme, NA_real_, jml_par$theta),
    extreme = extreme,
    bayes_lower = ci_bayes[1L, ], bayes_upper = ci_bayes[2L, ],
    mml_lower = mml$lower, mml_upper = mml$upper
  )
  diagnostics[[b]] <- data.frame(
    replication = b, max_rhat = max(s$rhat),
    min_ess_bulk = min(s$ess_bulk), min_ess_tail = min(s$ess_tail),
    divergences = sum(sampler$num_divergent),
    max_depth_hits = sum(sampler$num_max_treedepth),
    mml_convergence = mml$convergence, jml_convergence = jml$convergence
  )

  # Held-out predictive probabilities across ALL test observations (including extreme examinees!)
  cumulative_tau <- t(apply(tau_draws, 1L, cumsum))
  test <- dat$test
  pi_t <- as.integer(sub("P", "", test$person))
  ii_t <- as.integer(sub("I", "", test$item))
  ri_t <- as.integer(sub("R", "", test$rater))

  bayes_prob <- matrix(NA_real_, nrow(test), K)
  mml_prob   <- matrix(0, nrow(test), K)
  for (j in seq_len(nrow(test))) {
    eta_b <- theta_draws[, pi_t[j]] - item_draws[, ii_t[j]] - rater_draws[, ri_t[j]]
    logits_b <- outer(eta_b, 0:(K - 1L)) - cbind(0, cumulative_tau)
    z_b <- exp(logits_b - apply(logits_b, 1L, max))
    bayes_prob[j, ] <- colMeans(z_b / rowSums(z_b))

    # MML-EAP integrates over conditional ability quadrature weights given plug-in MML item/rater/tau
    eta_q <- mml$th_q - mml$par$item[ii_t[j]] - mml$par$rater[ri_t[j]]
    pr_q <- rsm_prob(eta_q, mml$par$tau)
    mml_prob[j, ] <- as.vector(mml$w_pq[pi_t[j], ] %*% pr_q)
  }

  jml_prob <- rsm_prob(jml_par$theta[pi_t] - jml_par$item[ii_t] - jml_par$rater[ri_t], jml_par$tau)
  for (p in which(extreme)) {
    rows <- which(pi_t == p)
    endpoint <- dat$train$score[match(sprintf("P%02d", p), dat$train$person)]
    jml_prob[rows, ] <- 0; jml_prob[rows, endpoint] <- 1
  }

  obs_idx <- cbind(seq_len(nrow(test)), test$score)
  comp_rows <- !extreme[pi_t]
  predictions[[b]] <- data.frame(
    replication = b,
    n_test_all = nrow(test),
    n_test_finite_jml = sum(comp_rows),
    bayes_logloss_all = mean(-log(bayes_prob[obs_idx])),
    mml_logloss_all   = mean(-log(mml_prob[obs_idx])),
    bayes_logloss_finite = mean(-log(bayes_prob[obs_idx][comp_rows])),
    mml_logloss_finite   = mean(-log(mml_prob[obs_idx][comp_rows])),
    jml_logloss_finite   = mean(-log(jml_prob[obs_idx][comp_rows]))
  )
  unlink(fit$output_files())
}

people <- do.call(rbind, people)
diagnostics <- do.call(rbind, diagnostics)
predictions <- do.call(rbind, predictions)

summary_ability <- do.call(rbind, lapply(sort(unique(counts)), function(n) {
  all_p <- people[people$n == n, ]
  fin_p <- all_p[!all_p$extreme, ]
  data.frame(
    n = n,
    persons_total = nrow(all_p),
    extreme_persons = sum(all_p$extreme),
    jml_rmse_finite   = sqrt(mean((fin_p$jml - fin_p$truth)^2)),
    mml_rmse_finite   = sqrt(mean((fin_p$mml_eap - fin_p$truth)^2)),
    bayes_rmse_finite = sqrt(mean((fin_p$bayes - fin_p$truth)^2)),
    mml_rmse_all      = sqrt(mean((all_p$mml_eap - all_p$truth)^2)),
    bayes_rmse_all    = sqrt(mean((all_p$bayes - all_p$truth)^2)),
    mml_cov95_all     = mean(all_p$mml_lower <= all_p$truth & all_p$truth <= all_p$mml_upper),
    bayes_cov95_all   = mean(all_p$bayes_lower <= all_p$truth & all_p$truth <= all_p$bayes_upper)
  )
}))

# -----------------------------------------------------------------------------
# Part 2: Sparse Interaction Bias Shrinkage Trade-Off (All vs. Null vs. Active Cells)
# -----------------------------------------------------------------------------
cat("\n=== Part 2: Sparse Bias Shrinkage Trade-Off (Null-Cell Noise vs. Active-Cell Attenuation) ===\n")
K_b <- 4L
tau_b <- c(-0.9, 0.0, 0.9)
raw_bias <- matrix(0, I, R)
raw_bias[1L, 1L] <-  0.80
raw_bias[4L, 4L] <- -0.80
truth_bias <- sweep(raw_bias, 1L, rowMeans(raw_bias), "-")
truth_bias <- sweep(truth_bias, 2L, colMeans(truth_bias), "-")
is_active_cell <- abs(truth_bias) >= 0.35 # The 2 genuinely biased cells (|B_ir| = 0.60 after centering)

# Sparse bridged plan: P = 48 examinees x 4 ratings = 192 ratings (12 ratings per item x rater cell, coprime stride 5)
grid_b <- expand.grid(item = seq_len(I), rater = seq_len(R))
idx_b <- unlist(lapply(seq_len(48L), function(p) (((p - 1L) * 5L + seq_len(4L) - 1L) %% 16L) + 1L))
plan_b <- data.frame(
  person = sprintf("P%02d", rep(seq_len(48L), each = 4L)),
  item   = sprintf("I%02d", grid_b$item[idx_b]),
  rater  = sprintf("R%02d", grid_b$rater[idx_b]),
  stringsAsFactors = FALSE
)
stopifnot(check_design_bmfrm(score ~ person + item + rater + item:rater, plan_b, K = K_b)$identifiable)

spec_unpen <- parse_bmfrm_formula(score ~ person + item + rater + item:rater,
                                  transform(plan_b, score = 1L), K = K_b,
                                  hierarchical = "person", model_name = "sim_bias_unpen")
priors_unpen <- c(
  prior("normal(0, 1.5)", "mean", facet = "person"),
  prior("lognormal(0, 0.4)", "sd", facet = "person"),
  prior("normal(0, 1.5)", "facet", facet = "item", scale = "marginal"),
  prior("normal(0, 1.5)", "facet", facet = "rater", scale = "marginal"),
  prior("normal(0, 2.0)", "tau", scale = "marginal"),
  prior("normal(0, 2.5)", "bias", facet = "item:rater", scale = "marginal")
)
comp_unpen <- compile_bmfrm(build_stan_code(spec_unpen, priors_unpen), "sim_bias_unpen",
                            cache_dir = scratch, cpp_options = list(PRECOMPILED_HEADERS = "false"))

spec_hs <- parse_bmfrm_formula(score ~ person + item + rater + item:rater,
                               transform(plan_b, score = 1L), K = K_b,
                               hierarchical = c("person", "item:rater"), model_name = "sim_bias_hs")
priors_hs <- c(
  prior("normal(0, 1.5)", "mean", facet = "person"),
  prior("lognormal(0, 0.4)", "sd", facet = "person"),
  prior("normal(0, 1.0)", "facet", facet = "item", scale = "marginal"),
  prior("normal(0, 1.0)", "facet", facet = "rater", scale = "marginal"),
  prior("normal(0, 2.0)", "tau", scale = "marginal"),
  prior("normal(0, 0.35)", "sd", facet = "item:rater", scale = "raw",
        local_shrinkage = "horseshoe", slab_df = 4, slab_scale = 1.5)
)
comp_hs <- compile_bmfrm(build_stan_code(spec_hs, priors_hs), "sim_bias_hs",
                         cache_dir = scratch, cpp_options = list(PRECOMPILED_HEADERS = "false"))

bias_rows <- vector("list", B * 2L)
cell_names <- sprintf("itemrater[%d,%d]", rep(seq_len(I), R), rep(seq_len(R), each = I))
truth_vec  <- as.vector(truth_bias)
active_vec <- as.vector(is_active_cell)

for (b in seq_len(B)) {
  set.seed(88000L + b)
  th_b <- stats::rnorm(48L, 0, 1)
  p_id <- as.integer(sub("^P", "", plan_b$person))
  i_id <- as.integer(sub("^I", "", plan_b$item))
  r_id <- as.integer(sub("^R", "", plan_b$rater))
  eta_b <- th_b[p_id] - truth_item[i_id] - truth_rater[r_id] + truth_bias[cbind(i_id, r_id)]
  logits <- sweep(outer(eta_b, 0:(K_b - 1L)), 2L, c(0, cumsum(tau_b)), "-")
  pr <- exp(logits - apply(logits, 1L, max)); pr <- pr / rowSums(pr)
  dat_b <- plan_b
  dat_b$score <- vapply(seq_len(nrow(dat_b)), function(j) sample.int(K_b, 1L, prob = pr[j, ]), integer(1))

  di_u <- prepare_data_bmfrm(spec_unpen, dat_b, priors_unpen)
  f_u <- comp_unpen$model$sample(di_u$stan_data, seed = 91000L + b, chains = 2L, parallel_chains = 2L,
                                 iter_warmup = 500L, iter_sampling = 600L, adapt_delta = 0.95,
                                 refresh = 0L, show_messages = FALSE, show_exceptions = FALSE)
  est_u <- colMeans(unclass(f_u$draws("itemrater", format = "matrix"))[, cell_names, drop = FALSE])
  unlink(f_u$output_files())

  di_h <- prepare_data_bmfrm(spec_hs, dat_b, priors_hs)
  f_h <- comp_hs$model$sample(di_h$stan_data, seed = 92000L + b, chains = 2L, parallel_chains = 2L,
                              iter_warmup = 600L, iter_sampling = 600L, adapt_delta = 0.98,
                              refresh = 0L, show_messages = FALSE, show_exceptions = FALSE)
  est_h <- colMeans(unclass(f_h$draws("itemrater", format = "matrix"))[, cell_names, drop = FALSE])
  unlink(f_h$output_files())

  bias_rows[[2L * b - 1L]] <- data.frame(
    replication = b, estimator = "Weakly-Unpenalized (SD=2.5)",
    rmse_all = sqrt(mean((est_u - truth_vec)^2)),
    rmse_null = sqrt(mean((est_u[!active_vec] - truth_vec[!active_vec])^2)),
    rmse_active = sqrt(mean((est_u[active_vec] - truth_vec[active_vec])^2)),
    stringsAsFactors = FALSE
  )
  bias_rows[[2L * b]] <- data.frame(
    replication = b, estimator = "Regularized Horseshoe",
    rmse_all = sqrt(mean((est_h - truth_vec)^2)),
    rmse_null = sqrt(mean((est_h[!active_vec] - truth_vec[!active_vec])^2)),
    rmse_active = sqrt(mean((est_h[active_vec] - truth_vec[active_vec])^2)),
    stringsAsFactors = FALSE
  )
}
bias_df <- do.call(rbind, bias_rows)
bias_summary <- stats::aggregate(cbind(rmse_all, rmse_null, rmse_active) ~ estimator, data = bias_df, FUN = mean)

write.csv(people, file.path(out, "person_results.csv"), row.names = FALSE)
write.csv(summary_ability, file.path(out, "summary.csv"), row.names = FALSE)
write.csv(diagnostics, file.path(out, "diagnostics.csv"), row.names = FALSE)
write.csv(predictions, file.path(out, "prediction_results.csv"), row.names = FALSE)
write.csv(bias_summary, file.path(out, "bias_tradeoff_summary.csv"), row.names = FALSE)
utils::capture.output(utils::sessionInfo(), file = file.path(out, "session_info.txt"))

png(file.path(out, "bayesian_benefits.png"), width = 1800, height = 600, res = 150)
par(mfrow = c(1, 3), mar = c(4.5, 4.3, 3.5, 1.0))
barplot(
  rbind(summary_ability$jml_rmse_finite, summary_ability$mml_rmse_finite, summary_ability$bayes_rmse_finite),
  beside = TRUE, names.arg = paste(summary_ability$n, "ratings"),
  col = c("#b0bec5", "#546e7a", "#197f9a"),
  ylab = "Ability RMSE (finite-JML subset)", main = "1. Ability: JML vs MML-EAP vs Bayes",
  ylim = c(0, max(summary_ability$jml_rmse_finite) * 1.28)
)
legend("topright", c("Joint ML", "MML-EAP", "BayesMFRM"), fill = c("#b0bec5", "#546e7a", "#197f9a"), bty = "n", cex = 0.82)

barplot(
  rbind(summary_ability$mml_cov95_all, summary_ability$bayes_cov95_all),
  beside = TRUE, names.arg = paste(summary_ability$n, "ratings"),
  col = c("#546e7a", "#197f9a"), ylim = c(0, 1.08),
  ylab = "Empirical 95% interval coverage (all persons)", main = "2. Uncertainty Calibration (incl. extremes)"
)
abline(h = 0.95, lty = 2, col = "#c62828", lwd = 1.5)
legend("bottomright", c("MML-EAP (plug-in)", "BayesMFRM (full)"), fill = c("#546e7a", "#197f9a"), bty = "n", cex = 0.82)

ord_b <- match(c("Weakly-Unpenalized (SD=2.5)", "Regularized Horseshoe"), bias_summary$estimator)
bs_ord <- bias_summary[ord_b, ]
barplot(
  as.matrix(bs_ord[, c("rmse_all", "rmse_null", "rmse_active")]),
  beside = TRUE, names.arg = c("All 16 cells", "Null/small cells", "Active biased cells"),
  col = c("#90a4ae", "#197f9a"),
  ylab = "Interaction RMSE (logits)", main = "3. Bias Shrinkage Trade-Off",
  ylim = c(0, max(bs_ord[, c("rmse_all", "rmse_null", "rmse_active")]) * 1.28)
)
legend("topright", c("Unpenalized", "Horseshoe"), fill = c("#90a4ae", "#197f9a"), bty = "n", cex = 0.82)
dev.off()

cat("\n=== Ability Summary (Joint ML vs. MML-EAP vs. BayesMFRM) ===\n")
print(summary_ability, digits = 4, row.names = FALSE)
cat("\n=== Held-Out Predictive Log-Loss ===\n")
print(colMeans(predictions[, -c(1:3)]), digits = 4)
cat("\n=== Sparse Bias Shrinkage Trade-Off ===\n")
print(bias_summary, digits = 4, row.names = FALSE)
