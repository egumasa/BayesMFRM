# Fixed-point log-density and finite-difference gradient parity check for
# extended BayesMFRM models:
#   1. Hierarchical main facets (person + rater + item)
#   2. Hierarchical PCM (common steps + level-specific deviations)
#   3. Sparse cell-space regularized horseshoe interaction prior
#   4. Cell-space Student-t local shrinkage interaction prior
#   5. Mixed non-Gaussian (student_t rater) + hierarchical person + marginal item
#   6. Minimal generated quantities (save_gq = "minimal") parity with full GQ
#
# Run: Rscript --vanilla tools/check-extensions-stan.R

script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1L])
root <- dirname(dirname(normalizePath(script)))
suppressPackageStartupMessages(pkgload::load_all(root, quiet = TRUE))
out <- file.path(root, ".cache", "extensions-parity")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

H <- function(n, scale = "raw") {
  b <- matrix(0, n, max(0L, n - 1L))
  if (n > 1L) for (j in seq_len(n - 1L)) {
    b[seq_len(j), j] <- -1 / sqrt(j * (j + 1))
    b[j + 1L, j] <- j / sqrt(j * (j + 1))
  }
  if (scale == "marginal" && n > 1L) b <- b * sqrt(n / (n - 1))
  b
}

double_center <- function(M) {
  M <- sweep(M, 1L, rowMeans(M), "-")
  sweep(M, 2L, colMeans(M), "-")
}

# Log-density for positive half-normal(0, s) on unconstrained h = log(sigma)
lp_half_normal <- function(h, s) {
  sig <- exp(h)
  -0.5 * (sig / s)^2 + h
}

# Log-density for inv_gamma(a, b) on unconstrained h = log(x)
lp_inv_gamma <- function(h, a, b) {
  x <- exp(h)
  -(a + 1) * h - b / x + h
}

# Log-density for half-cauchy(0, s) on unconstrained h = log(x)
lp_half_cauchy <- function(h, s) {
  x <- exp(h)
  stats::dcauchy(x, 0, s, log = TRUE) + h
}

d <- expand.grid(person = 1:5, item = 1:3, rater = 1:3)
d$score <- (seq_len(nrow(d)) - 1L) %% 4L + 1L

run_case <- function(name, spec, priors, ref_lp_fn, np, dat_in = d) {
  data <- prepare_data_bmfrm(spec, dat_in, priors)$stan_data
  code <- build_stan_code(spec, priors)
  compiled <- compile_bmfrm(code, paste0("ext_", name), out,
                            cpp_options = list(PRECOMPILED_HEADERS = "false"))
  dp <- file.path(out, paste0(name, "-data.json"))
  cmdstanr::write_stan_json(data, dp)

  points <- rbind(
    rep(0.11, np),
    0.25 * sin(seq_len(np)),
    0.45 * cos(seq_len(np)),
    seq(-0.5, 0.5, length.out = np) + 0.017
  )
  actual <- do.call(rbind, lapply(seq_len(nrow(points)), function(j) {
    qp <- file.path(out, paste0(name, "-", j, ".json"))
    op <- file.path(out, paste0(name, "-", j, ".csv"))
    cmdstanr::write_stan_json(list(params_r = as.numeric(points[j, ])), qp)
    rc <- system2(compiled$model$exe_file(),
                  c("log_prob", "jacobian=1", paste0("unconstrained_params=", shQuote(qp)),
                    "data", paste0("file=", shQuote(dp)),
                    "output", paste0("file=", shQuote(op)), "sig_figs=18"),
                  stdout = file.path(out, paste0(name, "-", j, ".log")), stderr = "")
    stopifnot(rc == 0L)
    utils::read.csv(op, comment.char = "#", check.names = FALSE)
  }))

  ref_lps <- vapply(seq_len(nrow(points)), function(j) ref_lp_fn(points[j, ], data), numeric(1))
  de <- max(abs((actual[[1]] - actual[[1]][1]) - (ref_lps - ref_lps[1])))

  # Central finite-difference check on gradient at point 2
  h_step <- 1e-6
  q0 <- points[2, ]
  num_g <- vapply(seq_len(np), function(k) {
    up <- dn <- q0; up[k] <- up[k] + h_step; dn[k] <- dn[k] - h_step
    (ref_lp_fn(up, data) - ref_lp_fn(dn, data)) / (2 * h_step)
  }, numeric(1))
  ge <- max(abs(as.numeric(actual[2, -1]) - num_g))
  cat(sprintf("%-24s dim=%2d  lp_diff_err=%.3e  grad_err=%.3e\n", name, np, de, ge))
  stopifnot(de < 1e-10, ge < 1e-6)
  data.frame(case = name, dim = np, lp_diff_err = de, grad_err = ge)
}

ll_ordinal <- function(eta, tau_mat, si, X, K) {
  logits <- outer(eta, 0:(K - 1L)) - t(vapply(si, function(j) c(0, cumsum(tau_mat[j, ])), numeric(K)))
  m <- apply(logits, 1L, max)
  sum(logits[cbind(seq_along(X), X)] - (m + log(rowSums(exp(logits - m)))))
}

rows <- list()

# Case 1: Hierarchical person + rater + item (minimal GQ)
ps1 <- c(
  prior("normal(0, 2)", "mean", facet = "person"),
  prior("lognormal(0, 0.4)", "sd", facet = "person"),
  prior("normal(0, 0.8)", "sd", facet = "item", scale = "marginal"),
  prior("normal(0, 0.7)", "sd", facet = "rater", scale = "marginal"),
  prior("normal(0, 1)", "tau", scale = "marginal")
)
s1 <- parse_bmfrm_formula(score ~ person + item + rater, d, K = 4,
                          hierarchical = c("person", "item", "rater"), save_gq = "minimal")
rows[[1]] <- run_case("hier_main_facets", s1, ps1, function(q, dat) {
  # Order in Stan parameters: theta_raw(5), mu_raw(1), lsig_raw(1), item_raw(2), sig_item(1), rater_raw(2), sig_rater(1), tau_raw(2)
  th_z <- q[1:5]; mu_z <- q[6]; ls_z <- q[7]
  it_z <- q[8:9]; hi_s <- q[10]; ra_z <- q[11:12]; hr_s <- q[13]; ta_z <- q[14:15]
  mu <- 2 * mu_z; sig <- exp(0.4 * ls_z); theta <- mu + sig * th_z
  item <- as.vector(exp(hi_s) * H(3, "marginal") %*% it_z)
  rater <- as.vector(exp(hr_s) * H(3, "marginal") %*% ra_z)
  tau <- matrix(as.vector(H(3, "marginal") %*% ta_z), 1, 3)
  eta <- theta[dat$person_id] - item[dat$item_id] - rater[dat$rater_id]
  lp <- -0.5 * sum(c(th_z, mu_z, ls_z, it_z, ra_z, ta_z)^2) +
    lp_half_normal(hi_s, 0.8) + lp_half_normal(hr_s, 0.7)
  ll_ordinal(eta, tau, rep(1L, dat$N), dat$X, dat$K) + lp
}, 15L)

# Case 2: Hierarchical PCM (common steps + level-specific deviations)
ps2 <- c(
  prior("normal(0, 2)", "mean", facet = "person"),
  prior("lognormal(0, 0.4)", "sd", facet = "person"),
  prior("normal(0, 1)", "item", scale = "marginal"),
  prior("normal(0, 1)", "rater", scale = "marginal"),
  prior("normal(0, 1)", "tau", scale = "marginal"),
  prior("normal(0, 0.5)", "sd", facet = "tau", scale = "marginal")
)
s2 <- parse_bmfrm_formula(score ~ person + item + rater, d, K = 4,
                          family = "pcm", threshold_facet = "item",
                          hierarchical = c("person", "tau"))
rows[[2]] <- run_case("hier_pcm_steps", s2, ps2, function(q, dat) {
  # Order: theta_raw(5), mu_z(1), ls_z(1), item_raw(2), rater_raw(2), tau_common_raw(2), tau_dev_raw(3*2=6), sig_tau(1)
  th_z <- q[1:5]; mu_z <- q[6]; ls_z <- q[7]
  it_z <- q[8:9]; ra_z <- q[10:11]
  tc_z <- q[12:13]; td_z <- matrix(q[14:19], 3, 2); ht_s <- q[20]
  theta <- 2 * mu_z + exp(0.4 * ls_z) * th_z
  item <- as.vector(H(3, "marginal") %*% it_z)
  rater <- as.vector(H(3, "marginal") %*% ra_z)
  Bt <- H(3, "marginal")
  tau_common <- as.vector(Bt %*% tc_z)
  tau_dev <- exp(ht_s) * td_z %*% t(Bt)
  tau_mat <- sweep(tau_dev, 2L, tau_common, "+")
  eta <- theta[dat$person_id] - item[dat$item_id] - rater[dat$rater_id]
  lp <- -0.5 * sum(c(th_z, mu_z, ls_z, it_z, ra_z, tc_z, as.vector(td_z))^2) +
    lp_half_normal(ht_s, 0.5)
  ll_ordinal(eta, tau_mat, dat$item_id, dat$X, dat$K) + lp
}, 20L)

# Case 3: Regularized horseshoe cell-space interaction (scale = "raw" required)
ps3 <- c(
  prior("normal(0, 2)", "mean", facet = "person"),
  prior("lognormal(0, 0.4)", "sd", facet = "person"),
  prior("normal(0, 1)", "item", scale = "marginal"),
  prior("normal(0, 1)", "rater", scale = "marginal"),
  prior("normal(0, 1)", "tau", scale = "marginal"),
  prior("normal(0, 0.5)", "sd", facet = "item:rater", scale = "raw", local_shrinkage = "horseshoe")
)
s3 <- parse_bmfrm_formula(score ~ person + item + rater + item:rater, d, K = 4,
                          hierarchical = c("person", "item:rater"))
rows[[3]] <- run_case("horseshoe_bias", s3, ps3, function(q, dat) {
  # Order: theta(5), mu(1), ls(1), item(2), rater(2), tau(2), b_raw(9), log_lam(9), log_c2(1), log_sig(1)
  th_z <- q[1:5]; mu_z <- q[6]; ls_z <- q[7]
  it_z <- q[8:9]; ra_z <- q[10:11]; ta_z <- q[12:13]
  b_z <- matrix(q[14:22], 3, 3); h_lam <- matrix(q[23:31], 3, 3); h_c2 <- q[32]; h_sig <- q[33]
  theta <- 2 * mu_z + exp(0.4 * ls_z) * th_z
  item <- as.vector(H(3, "marginal") %*% it_z)
  rater <- as.vector(H(3, "marginal") %*% ra_z)
  tau <- matrix(as.vector(H(3, "marginal") %*% ta_z), 1, 3)
  lam <- exp(h_lam); c2 <- exp(h_c2); sig <- exp(h_sig)
  lam_tilde <- sqrt(c2 * lam^2 / (c2 + sig^2 * lam^2))
  B <- double_center(sig * lam_tilde * b_z)
  eta <- theta[dat$person_id] - item[dat$item_id] - rater[dat$rater_id] + B[cbind(dat$item_id, dat$rater_id)]
  lp <- -0.5 * sum(c(th_z, mu_z, ls_z, it_z, ra_z, ta_z, as.vector(b_z))^2) +
    sum(lp_half_cauchy(as.vector(h_lam), 1)) + lp_inv_gamma(h_c2, 2, 8) + lp_half_normal(h_sig, 0.5)
  ll_ordinal(eta, tau, rep(1L, dat$N), dat$X, dat$K) + lp
}, 33L)

# Case 4: Unit-variance Student-t local shrinkage interaction (scale = "marginal" preserves SD = sigma)
ps4 <- ps3
ps4[[6L]] <- prior("normal(0, 0.5)", "sd", facet = "item:rater", scale = "marginal", local_shrinkage = "student_t")
rows[[4]] <- run_case("student_t_bias", s3, ps4, function(q, dat) {
  th_z <- q[1:5]; mu_z <- q[6]; ls_z <- q[7]
  it_z <- q[8:9]; ra_z <- q[10:11]; ta_z <- q[12:13]
  b_z <- matrix(q[14:22], 3, 3); h_lam <- matrix(q[23:31], 3, 3); h_sig <- q[32]
  theta <- 2 * mu_z + exp(0.4 * ls_z) * th_z
  item <- as.vector(H(3, "marginal") %*% it_z)
  rater <- as.vector(H(3, "marginal") %*% ra_z)
  tau <- matrix(as.vector(H(3, "marginal") %*% ta_z), 1, 3)
  lam <- exp(h_lam); sig <- exp(h_sig)
  sf <- sqrt((3 / 2) * (3 / 2))
  B <- double_center(sf * sig * sqrt(lam) * b_z)
  eta <- theta[dat$person_id] - item[dat$item_id] - rater[dat$rater_id] + B[cbind(dat$item_id, dat$rater_id)]
  lp <- -0.5 * sum(c(th_z, mu_z, ls_z, it_z, ra_z, ta_z, as.vector(b_z))^2) +
    sum(lp_inv_gamma(as.vector(h_lam), 1.5, 0.5)) + lp_half_normal(h_sig, 0.5)
  ll_ordinal(eta, tau, rep(1L, dat$N), dat$X, dat$K) + lp
}, 32L)

# Case 5: Mixed non-Gaussian rater (student_t) + hierarchical person + marginal item
ps5 <- c(
  prior("normal(0, 2)", "mean", facet = "person"),
  prior("lognormal(0, 0.4)", "sd", facet = "person"),
  prior("normal(0, 1)", "item", scale = "marginal"),
  prior("student_t(3, 0, 1)", "rater"),
  prior("normal(0, 1)", "tau", scale = "marginal")
)
s5 <- parse_bmfrm_formula(score ~ person + item + rater, d, K = 4, hierarchical = "person")
rows[[5]] <- run_case("mixed_student_t_hier", s5, ps5, function(q, dat) {
  # Order: theta(5), mu(1), ls(1), item(2), rater_raw(3), tau(2)
  th_z <- q[1:5]; mu_z <- q[6]; ls_z <- q[7]
  it_z <- q[8:9]; ra_raw <- q[10:12]; ta_z <- q[13:14]
  theta <- 2 * mu_z + exp(0.4 * ls_z) * th_z
  item <- as.vector(H(3, "marginal") %*% it_z)
  rater <- ra_raw - mean(ra_raw)
  tau <- matrix(as.vector(H(3, "marginal") %*% ta_z), 1, 3)
  eta <- theta[dat$person_id] - item[dat$item_id] - rater[dat$rater_id]
  lp <- -0.5 * sum(c(th_z, mu_z, ls_z, it_z, ta_z)^2) + sum(stats::dt(ra_raw, 3, log = TRUE))
  ll_ordinal(eta, tau, rep(1L, dat$N), dat$X, dat$K) + lp
}, 14L)

# Cases 6 & 7: Heterogeneous-category PCM (bmfrm_pcm_hetero.stan) with K = c(I1=2, I2=4, I3=6)
d_het <- expand.grid(person = 1:5, item = c("I1", "I2", "I3"), rater = 1:3, stringsAsFactors = FALSE)
K_het <- c(I1 = 2L, I2 = 4L, I3 = 6L)
d_het$score <- vapply(seq_len(nrow(d_het)), function(i) {
  ((i - 1L) %% K_het[[d_het$item[i]]]) + 1L
}, integer(1))

ll_hetero_pcm <- function(eta, tau_list, item_id, X, K_lvl) {
  sum(vapply(seq_along(X), function(n) {
    j <- item_id[n]
    Kj <- K_lvl[j]
    logits <- (0:(Kj - 1L)) * eta[n] - c(0, cumsum(tau_list[[j]]))
    m <- max(logits)
    logits[X[n]] - (m + log(sum(exp(logits - m))))
  }, numeric(1)))
}

# Case 6: Heterogeneous PCM with marginal Gaussian basis (B_tau_rows) + hierarchical rater
ps6 <- c(
  prior("normal(0, 2)", "mean", facet = "person"),
  prior("lognormal(0, 0.4)", "sd", facet = "person"),
  prior("normal(0, 0.9)", "item", scale = "marginal"),
  prior("normal(0, 0.7)", "sd", facet = "rater", scale = "marginal"),
  prior("normal(0, 0.8)", "tau", scale = "marginal")
)
s6 <- parse_bmfrm_formula(score ~ person + item + rater, d_het, K = K_het,
                          family = "pcm", threshold_facet = "item",
                          hierarchical = c("person", "rater"), save_gq = "minimal")
rows[[6]] <- run_case("hetero_pcm_marginal", s6, ps6, function(q, dat) {
  # Order: theta_raw(5), mu_z(1), ls_z(1), item_raw(2), rater_raw(2), sig_rater(1), tau_raw(0+2+4=6)
  th_z <- q[1:5]; mu_z <- q[6]; ls_z <- q[7]
  it_z <- q[8:9]; ra_z <- q[10:11]; hr_s <- q[12]; ta_z <- q[13:18]
  theta <- 2 * mu_z + exp(0.4 * ls_z) * th_z
  item <- as.vector(0.9 * H(3, "marginal") %*% it_z)
  rater <- as.vector(exp(hr_s) * H(3, "marginal") %*% ra_z)
  tau_list <- list(
    0,
    as.vector(0.8 * H(3, "marginal") %*% ta_z[1:2]),
    as.vector(0.8 * H(5, "marginal") %*% ta_z[3:6])
  )
  eta <- theta[dat$person_id] - item[dat$item_id] - rater[dat$rater_id]
  lp <- -0.5 * sum(c(th_z, mu_z, ls_z, it_z, ra_z, ta_z)^2) + lp_half_normal(hr_s, 0.7)
  ll_hetero_pcm(eta, tau_list, dat$item_id, dat$X, dat$K_level) + lp
}, 18L, dat_in = d_het)

# Case 7: Heterogeneous PCM with non-Gaussian raw Student-t steps (ragged centered tau_raw)
ps7 <- c(
  prior("normal(0, 2)", "mean", facet = "person"),
  prior("lognormal(0, 0.4)", "sd", facet = "person"),
  prior("normal(0, 1)", "item", scale = "marginal"),
  prior("normal(0, 1)", "rater", scale = "marginal"),
  prior("student_t(4, 0, 0.6)", "tau")
)
s7 <- parse_bmfrm_formula(score ~ person + item + rater, d_het, K = K_het,
                          family = "pcm", threshold_facet = "item",
                          hierarchical = "person", save_gq = "full")
rows[[7]] <- run_case("hetero_pcm_student_t", s7, ps7, function(q, dat) {
  # Order: theta_raw(5), mu_z(1), ls_z(1), item_raw(2), rater_raw(2), tau_raw(0+3+5=8)
  th_z <- q[1:5]; mu_z <- q[6]; ls_z <- q[7]
  it_z <- q[8:9]; ra_z <- q[10:11]; tr <- q[12:19]
  theta <- 2 * mu_z + exp(0.4 * ls_z) * th_z
  item <- as.vector(H(3, "marginal") %*% it_z)
  rater <- as.vector(H(3, "marginal") %*% ra_z)
  tau_list <- list(0, tr[1:3] - mean(tr[1:3]), tr[4:8] - mean(tr[4:8]))
  eta <- theta[dat$person_id] - item[dat$item_id] - rater[dat$rater_id]
  lp <- -0.5 * sum(c(th_z, mu_z, ls_z, it_z, ra_z)^2) +
    sum(stats::dt(tr / 0.6, df = 4, log = TRUE) - log(0.6))
  ll_hetero_pcm(eta, tau_list, dat$item_id, dat$X, dat$K_level) + lp
}, 19L, dat_in = d_het)

res_df <- do.call(rbind, rows)
utils::write.csv(res_df, file.path(out, "extensions-parity.csv"), row.names = FALSE)
cat("All extended model fixed-point parity checks passed!\n")
