# Reproduce from the package root:
# Rscript --vanilla simulations/bayesian_benefits.R [replications=20]
# Uses this checkout's RSM generator and default priors, without installing it.
# Comparison: unpenalized joint ML with the identical likelihood/constraints.
# This is a small illustrative experiment, not a validation of every package API.
args <- commandArgs(trailingOnly = TRUE)
B <- if (length(args)) as.integer(args[1]) else 20L
stopifnot(!is.na(B), B > 0L)
required <- c("pkgload", "cmdstanr", "posterior", "digest")
stopifnot(all(vapply(required, requireNamespace, logical(1), quietly = TRUE)))
# Load the complete local package and its templates, even if another version is installed.
suppressPackageStartupMessages(pkgload::load_all(".", quiet = TRUE,
                                                export_all = FALSE, export_imports = FALSE))
out <- file.path("simulations", "results")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
scratch <- tempfile("bayesmfrm_sim_")
dir.create(scratch)
P <- 40L; I <- 4L; R <- 4L; K <- 6L
truth_item <- c(-0.9, -0.3, 0.3, 0.9)
truth_rater <- c(-0.8, -0.2, 0.2, 0.8)
truth_tau <- seq(-1.5, 1.5, length.out = K - 1L)
counts <- rep(c(3L, 12L), each = P / 2)

rsm_prob <- function(eta, tau) {
  logits <- sweep(outer(eta, 0:(K - 1L)), 2L, c(0, cumsum(tau)), "-")
  logits <- logits - apply(logits, 1L, max)
  z <- exp(logits)
  z / rowSums(z)
}

simulate_data <- function(seed) {
  set.seed(seed)
  theta <- rnorm(P, 0, 1) # Deliberately narrower than default prior N(0, 2^2).
  d <- expand.grid(person = seq_len(P), item = seq_len(I), rater = seq_len(R))
  probs <- rsm_prob(theta[d$person] - truth_item[d$item] - truth_rater[d$rater], truth_tau)
  d$score <- vapply(seq_len(nrow(d)), function(j) sample.int(K, 1L, prob = probs[j, ]), integer(1))
  selected <- unlist(lapply(seq_len(P), function(p) sample(which(d$person == p), counts[p])))
  d$person <- sprintf("P%02d", d$person)
  d$item <- sprintf("I%02d", d$item)
  d$rater <- sprintf("R%02d", d$rater)
  list(train = d[selected, ], test = d[-selected, ], theta = theta)
}

# Exact sum-to-zero contrasts remove the ML location indeterminacies.
unpack <- function(x) {
  item <- x[P + seq_len(I - 1L)]
  rater <- x[P + I - 1L + seq_len(R - 1L)]
  tau <- x[P + I + R - 2L + seq_len(K - 2L)]
  list(theta = x[seq_len(P)], item = c(item, -sum(item)),
       rater = c(rater, -sum(rater)), tau = c(tau, -sum(tau)))
}

make_objective <- function(d) {
  pi <- as.integer(factor(d$person)); ii <- as.integer(factor(d$item)); ri <- as.integer(factor(d$rater))
  indicators <- outer(d$score, seq_len(K - 1L), ">")
  group_sum <- function(x, id, n) as.numeric(rowsum(x, factor(id, levels = seq_len(n)), reorder = TRUE))
  function(x, gradient = FALSE) {
    a <- unpack(x)
    eta <- a$theta[pi] - a$item[ii] - a$rater[ri]
    logits <- sweep(outer(eta, 0:(K - 1L)), 2L, c(0, cumsum(a$tau)), "-")
    shift <- apply(logits, 1L, max)
    z <- exp(logits - shift)
    if (!gradient) return(sum(log(rowSums(z)) + shift - logits[cbind(seq_len(nrow(d)), d$score)]))
    prob <- z / rowSums(z)
    g_eta <- as.vector(prob %*% seq_len(K)) - d$score
    g_theta <- group_sum(g_eta, pi, P)
    g_item <- -group_sum(g_eta, ii, I)
    g_rater <- -group_sum(g_eta, ri, R)
    tails <- vapply(seq_len(K - 1L), function(k) rowSums(prob[, (k + 1L):K, drop = FALSE]), numeric(nrow(d)))
    g_tau <- colSums(indicators - tails)
    c(g_theta, g_item[-I] - g_item[I], g_rater[-R] - g_rater[R],
      g_tau[-(K - 1L)] - g_tau[K - 1L])
  }
}

init_data <- simulate_data(20261124L)
spec <- parse_bmfrm_formula(score ~ person + item + rater, init_data$train, K = K, model_name = "simulation_rsm")
priors <- get_default_priors(spec$facets_main, spec$facets_bias)
stan_code <- build_stan_code(spec, priors)
stan_path <- file.path(scratch, "simulation_rsm.stan")
writeLines(stan_code, stan_path)
writeLines(stan_code, file.path(out, "simulation_rsm.stan"))
# Disable shared precompiled headers: compilation stays in the temporary directory.
model <- cmdstanr::cmdstan_model(stan_path, cpp_options = list(PRECOMPILED_HEADERS = "false"))

# Validate the ML gradient independently by central differences, once.
objective <- make_objective(init_data$train)
set.seed(12)
probe <- rnorm(P + I + R + K - 4L, 0, 0.1)
h <- 1e-5
numeric_gradient <- vapply(seq_along(probe), function(j) {
  upper <- lower <- probe
  upper[j] <- upper[j] + h; lower[j] <- lower[j] - h
  (objective(upper) - objective(lower)) / (2 * h)
}, numeric(1))
gradient_error <- max(abs(numeric_gradient - objective(probe, TRUE)))
stopifnot(gradient_error < 1e-5)
cat("ML gradient maximum absolute error:", gradient_error, "\n")

people <- vector("list", B); diagnostics <- vector("list", B); predictions <- vector("list", B)
decision_examples <- NULL
for (b in seq_len(B)) {
  cat(sprintf("Replication %d / %d\n", b, B))
  dat <- simulate_data(20261123L + b)
  prepared <- prepare_data_bmfrm(spec, dat$train)
  fit <- model$sample(data = prepared$stan_data, seed = 74000L + b,
                     chains = 4L, parallel_chains = 4L,
                     iter_warmup = 750L, iter_sampling = 1000L,
                     adapt_delta = 0.95, max_treedepth = 12L,
                     refresh = 0L, show_messages = FALSE)
  vars <- c("theta", "item", "rater", "tau")
  s <- fit$summary(vars)
  sampler <- fit$diagnostic_summary()
  # Plain matrix indexing drops a single column to a vector; draws_matrix does not.
  draws <- unclass(fit$draws(vars, format = "matrix"))
  theta_draws <- draws[, sprintf("theta[%d]", seq_len(P)), drop = FALSE]
  item_draws <- draws[, sprintf("item[%d]", seq_len(I)), drop = FALSE]
  rater_draws <- draws[, sprintf("rater[%d]", seq_len(R)), drop = FALSE]
  tau_draws <- draws[, sprintf("tau[%d]", seq_len(K - 1L)), drop = FALSE]

  obj <- make_objective(dat$train)
  start <- c(colMeans(theta_draws), colMeans(item_draws)[-I],
             colMeans(rater_draws)[-R], colMeans(tau_draws)[-(K - 1L)])
  ml <- optim(start, obj, function(x) obj(x, TRUE), method = "BFGS",
              control = list(maxit = 2000L, reltol = 1e-10))
  ml_par <- unpack(ml$par)
  ml_gradient <- max(abs(obj(ml$par, TRUE)))
  stopifnot(ml$convergence == 0L, ml_gradient < 0.01)
  # All-minimum/all-maximum response patterns have infinite unpenalized ML ability.
  # Their arbitrary finite optimizer stopping values are NEVER used in RMSE.
  extreme <- vapply(split(dat$train$score, dat$train$person),
                    function(y) all(y == 1L) || all(y == K), logical(1))
  ci <- apply(theta_draws, 2L, quantile, probs = c(0.025, 0.975))
  people[[b]] <- data.frame(replication = b, person = seq_len(P), n = counts,
                            truth = dat$theta, bayes = colMeans(theta_draws),
                            ml = ifelse(extreme, NA_real_, ml_par$theta), extreme = extreme,
                            lower = ci[1L, ], upper = ci[2L, ],
                            prob_above_zero = colMeans(theta_draws > 0))
  diagnostics[[b]] <- data.frame(replication = b, max_rhat = max(s$rhat),
                                  min_ess_bulk = min(s$ess_bulk), min_ess_tail = min(s$ess_tail),
                                  divergences = sum(sampler$num_divergent),
                                  max_depth_hits = sum(sampler$num_max_treedepth),
                                  ml_convergence = ml$convergence, ml_max_gradient = ml_gradient)

  # Independent held-out responses from the same latent persons/items/raters.
  # Integrate category probabilities over joint posterior draws; retain their covariance.
  nd <- nrow(draws)
  cumulative_tau <- t(apply(tau_draws, 1L, cumsum))
  test <- dat$test
  pi <- as.integer(sub("P", "", test$person)); ii <- as.integer(sub("I", "", test$item)); ri <- as.integer(sub("R", "", test$rater))
  bayes_prob <- matrix(NA_real_, nrow(test), K)
  for (j in seq_len(nrow(test))) {
    eta <- theta_draws[, pi[j]] - item_draws[, ii[j]] - rater_draws[, ri[j]]
    logits <- outer(eta, 0:(K - 1L)) - cbind(0, cumulative_tau)
    z <- exp(logits - apply(logits, 1L, max))
    bayes_prob[j, ] <- colMeans(z / rowSums(z))
  }
  ml_prob <- rsm_prob(ml_par$theta[pi] - ml_par$item[ii] - ml_par$rater[ri], ml_par$tau)
  # Use the exact limiting ML predictive distributions for extreme persons.
  for (p in which(extreme)) {
    rows <- which(pi == p)
    endpoint <- dat$train$score[match(sprintf("P%02d", p), dat$train$person)]
    ml_prob[rows, ] <- 0
    ml_prob[rows, endpoint] <- 1
  }
  observed <- cbind(seq_len(nrow(test)), test$score)
  # Finite-ML group only for comparable log loss; no arbitrary probability clipping.
  comparable <- !extreme[pi]
  predictions[[b]] <- data.frame(replication = b, n_test = sum(comparable),
                                  bayes_logloss = mean(-log(bayes_prob[observed][comparable])),
                                  ml_logloss = mean(-log(ml_prob[observed][comparable])),
                                  bayes_brier_all = mean(rowSums((bayes_prob - outer(test$score, seq_len(K), "=="))^2)),
                                  ml_brier_all = mean(rowSums((ml_prob - outer(test$score, seq_len(K), "=="))^2)))
  if (b == 1L) {
    # Preselected example: R04, with true severity 0.8, threshold 0.5 logits.
    rd <- rater_draws[, R]
    decision_examples <- data.frame(parameter = "R04_severity", truth = truth_rater[R],
                                    posterior_mean = mean(rd), lower = unname(quantile(rd, 0.025)),
                                    upper = unname(quantile(rd, 0.975)),
                                    threshold = 0.5, probability_above = mean(rd > 0.5))
    # Standardize fair expected scores to the same four items, average rater = 0.
    fair <- matrix(0, nd, P)
    for (p in seq_len(P)) for (i in seq_len(I)) {
      logits <- outer(theta_draws[, p] - item_draws[, i], 0:(K - 1L)) - cbind(0, cumulative_tau)
      z <- exp(logits - apply(logits, 1L, max))
      fair[, p] <- fair[, p] + as.vector((z / rowSums(z)) %*% seq_len(K)) / I
    }
    fair_ci <- apply(fair, 2L, quantile, probs = c(0.025, 0.975))
    write.csv(data.frame(person = seq_len(P), n = counts, mean = colMeans(fair),
                         lower = fair_ci[1L, ], upper = fair_ci[2L, ],
                         probability_above_3.5 = colMeans(fair > 3.5)),
              file.path(out, "fair_scores_first_replication.csv"), row.names = FALSE)
    write.csv(dat$train, file.path(out, "train_first_replication.csv"), row.names = FALSE)
    write.csv(dat$test, file.path(out, "test_first_replication.csv"), row.names = FALSE)
    saveRDS(draws, file.path(out, "posterior_first_replication.rds"))
  }
  # Raw chain files are temporary; derived outputs and first posterior are retained.
  unlink(fit$output_files())
}

people <- do.call(rbind, people); diagnostics <- do.call(rbind, diagnostics)
predictions <- do.call(rbind, predictions)
summary_rows <- lapply(sort(unique(counts)), function(n) {
  all <- people[people$n == n, ]; paired <- all[!all$extreme, ]
  improvement <- tapply((paired$ml - paired$truth)^2 - (paired$bayes - paired$truth)^2,
                        paired$replication, mean)
  data.frame(n = n, persons_total = nrow(all), finite_ml_persons = nrow(paired),
             extreme_persons = sum(all$extreme),
             ml_rmse_finite = sqrt(mean((paired$ml - paired$truth)^2)),
             bayes_rmse_finite = sqrt(mean((paired$bayes - paired$truth)^2)),
             bayes_rmse_all = sqrt(mean((all$bayes - all$truth)^2)),
             mean_95_interval_width = mean(all$upper - all$lower),
             coverage_95 = mean(all$lower <= all$truth & all$truth <= all$upper),
             mse_improvement = mean(improvement),
             mse_improvement_mcse = sd(improvement) / sqrt(length(improvement)))
})
summary_table <- do.call(rbind, summary_rows)
write.csv(people, file.path(out, "person_results.csv"), row.names = FALSE)
write.csv(summary_table, file.path(out, "summary.csv"), row.names = FALSE)
write.csv(diagnostics, file.path(out, "diagnostics.csv"), row.names = FALSE)
write.csv(predictions, file.path(out, "prediction_results.csv"), row.names = FALSE)
write.csv(decision_examples, file.path(out, "decision_example.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(out, "session_info.txt"))
capture.output(list(replications = B, truth_sd_theta = 1, prior_sd_theta = 2,
                    truth_item = truth_item, truth_rater = truth_rater, truth_tau = truth_tau,
                    chains = 4, warmup = 750, sampling = 1000,
                    data_seed = "20261123 + replication", stan_seed = "74000 + replication",
                    stan_code_hash = digest::digest(stan_code),
                    cmdstan_version = as.character(cmdstanr::cmdstan_version()),
                    ml_gradient_validation = gradient_error), file = file.path(out, "settings.txt"))

png(file.path(out, "bayesian_benefits.png"), width = 1650, height = 600, res = 150)
par(mfrow = c(1, 3), mar = c(4.5, 4.2, 3.5, 1))
barplot(rbind(summary_table$ml_rmse_finite, summary_table$bayes_rmse_finite), beside = TRUE,
        names.arg = paste(counts[c(1, P)], "ratings"), col = c("#a9b7c5", "#197f9a"),
        ylab = "Ability RMSE (lower is better)", main = "1. More stable estimation", ylim = c(0, max(summary_table$ml_rmse_finite) * 1.3))
legend("topright", c("Joint ML", "Bayes"), fill = c("#a9b7c5", "#197f9a"), bty = "n", cex = 0.8)
plot(summary_table$n, summary_table$mean_95_interval_width, type = "b", pch = 19,
     col = "#197f9a", xlab = "Ratings per person", ylab = "Mean 95% credible interval width",
     main = "2. Quantified uncertainty", xaxt = "n", ylim = c(0, max(summary_table$mean_95_interval_width) * 1.2))
axis(1, at = summary_table$n)
losses <- colMeans(predictions[c("ml_logloss", "bayes_logloss")])
barplot(losses, names.arg = c("Joint ML", "Bayes"), col = c("#a9b7c5", "#197f9a"),
        ylab = "Held-out log loss (lower is better)", main = "3. Posterior prediction", ylim = c(0, max(losses) * 1.25))
mtext("RMSE and log loss exclude persons with infinite ML ability", side = 1, line = 3.5, cex = 0.55)
dev.off()
cat("\nAbility results (paired finite-ML cases for the RMSE comparison):\n")
print(summary_table, digits = 4, row.names = FALSE)
cat("\nPreselected first-replication decision example:\n")
print(decision_examples, digits = 4, row.names = FALSE)
cat("\nMean predictive scores:\n")
print(colMeans(predictions[-c(1, 2)]), digits = 4)
improvement <- predictions$ml_logloss - predictions$bayes_logloss
cat("Paired mean log-loss improvement:", mean(improvement), "MCSE:", sd(improvement) / sqrt(B), "\n")
cat("\nMCMC: max Rhat", max(diagnostics$max_rhat), "; min bulk ESS", min(diagnostics$min_ess_bulk),
    "; min tail ESS", min(diagnostics$min_ess_tail), "; divergences", sum(diagnostics$divergences),
    "; max depth hits", sum(diagnostics$max_depth_hits), "\n")
if (any(diagnostics$max_rhat > 1.01 | diagnostics$min_ess_bulk < 400 |
        diagnostics$min_ess_tail < 400 | diagnostics$divergences > 0 | diagnostics$max_depth_hits > 0)) {
  warning("Some MCMC diagnostics require review; see diagnostics.csv.")
}
