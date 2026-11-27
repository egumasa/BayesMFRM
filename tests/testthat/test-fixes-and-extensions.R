# Tests for local bug fixes, API consistency, lightweight GQ, standardized fair
# scores, hierarchical main facets, sparse interaction priors, and hierarchical PCM.

test_that("add_zstd implements exact Wright & Masters Wilson-Hilferty formula", {
  df <- tibble::tibble(
    facet_value = 1:3,
    n_obs = c(10, 50, 100),
    infit_msq = c(1.0, 1.5, 0.7),
    outfit_msq = c(1.0, 1.5, 0.7)
  )
  res <- add_zstd(df)
  q_exp <- sqrt(2 / pmax(df$n_obs - 1, 1))
  expected <- (df$infit_msq^(1/3) - 1) * (3 / q_exp) + q_exp / 3
  expect_equal(res$infit_zstd, expected, tolerance = 1e-12)
  expect_equal(res$outfit_zstd, expected, tolerance = 1e-12)
  # Confirm clear misfit (MSQ = 1.5 at n = 50) is flagged (> 2), not shrunk by 1/3
  expect_gt(res$infit_zstd[2], 2.0)

  # When exact kurtosis-based q columns are supplied, add_zstd uses them
  df$infit_q <- c(0.4, 0.2, 0.15)
  df$outfit_q <- c(0.5, 0.25, 0.18)
  res_exact <- add_zstd(df)
  expect_equal(res_exact$infit_zstd, (df$infit_msq^(1/3) - 1) * (3 / df$infit_q) + df$infit_q / 3, tolerance = 1e-12)
  expect_equal(res_exact$outfit_zstd, (df$outfit_msq^(1/3) - 1) * (3 / df$outfit_q) + df$outfit_q / 3, tolerance = 1e-12)
})

test_that("facet_fit_from_matrices aligns point estimates with draw CrIs and preserves multi-model joins", {
  set.seed(101)
  D <- 200L; N <- 12L
  x <- rep(c(1L, 4L, 2L, 3L), 3L)
  fid <- rep(1:2, each = 6L)
  # Large posterior uncertainty in mu creates a gap between plug-in and draw mean
  mu_mat <- matrix(2.5 + stats::rnorm(D * N, 0, 0.8), D, N)
  sig_mat <- matrix(1.0, D, N)
  kurt_mat <- matrix(2.2, D, N)

  ff <- facet_fit_from_matrices(x, fid, mu_mat, sig_mat, kurt_mat = kurt_mat, model_name = "m1")
  # Draw-based point estimate is strictly inside its 95% CrI (no Jensen mismatch)
  expect_true(all(ff$infit_msq >= ff$infit_lo & ff$infit_msq <= ff$infit_hi))
  expect_true(all(ff$outfit_msq >= ff$outfit_lo & ff$outfit_msq <= ff$outfit_hi))
  expect_true(all(c("infit_msq_plugin", "outfit_msq_plugin", "infit_q", "outfit_q") %in% names(ff)))

  # Multi-model add_fit_cri_bayes must not duplicate rows
  res_m1 <- tibble::tibble(
    n = seq_len(N), x = x, rater = fid,
    mu_hat = colMeans(mu_mat), sigma2_hat = colMeans(sig_mat),
    resid = x - colMeans(mu_mat), weight = colMeans(sig_mat),
    z = (x - colMeans(mu_mat)) / sqrt(colMeans(sig_mat)),
    model = "m1",
    mu_draws = lapply(seq_len(N), function(j) as.numeric(mu_mat[, j])),
    sigma2_draws = lapply(seq_len(N), function(j) as.numeric(sig_mat[, j]))
  )
  res_m2 <- res_m1
  res_m2$model <- "m2"
  res_both <- dplyr::bind_rows(res_m1, res_m2)
  fs_both <- facet_infit_outfit(res_both, "rater")
  joined <- add_fit_cri_bayes(fs_both, res_both, "rater")
  expect_equal(nrow(joined), 4L)
  expect_true(all(c("infit_post_mean", "outfit_post_mean") %in% names(joined)))
})

test_that("summarise_bias respects prob_thresh directly without 95% CI contradiction", {
  d <- expand.grid(person = paste0("P", 1:4), item = paste0("I", 1:2), rater = c("R1", "R2"))
  d$score <- 1L
  spec <- parse_bmfrm_formula(score ~ person + item + rater + item:rater, d, K = 4)
  ch <- check_design_bmfrm(spec$formula, d, K = 4)

  # Construct 1000 draws where cell [1,1] has P(B > 0.2) = 0.960
  # (95% two-sided CI lower bound is < 0.2, while P(B > 0.2) >= 0.95)
  D <- 1000L
  c11 <- c(rep(0.15, 40L), rep(0.40, 960L))
  mat <- cbind(c11, -c11, -c11, c11)
  colnames(mat) <- c("itemrater[1,1]", "itemrater[2,1]", "itemrater[1,2]", "itemrater[2,2]")

  obj <- structure(list(
    spec = spec,
    design = ch,
    data_info = list(facet_labels = list(person = paste0("P", 1:4), item = paste0("I", 1:2), rater = c("R1", "R2")),
                     data_clean = d),
    fit = list(draws = function(...) mat)
  ), class = "bmfrm_fit")

  s95 <- summarise_bias(obj, "item:rater", prob_thresh = 0.95, delta = 0.2, flag_method = "prob")
  expect_equal(s95$flag_bias[1], "positive_bias")
  expect_equal(s95$flag_bias[2], "negative_bias")
  expect_true("bfdr" %in% names(s95))
  expect_equal(s95$bfdr[1], 0.04, tolerance = 1e-12)

  s975 <- summarise_bias(obj, "item:rater", prob_thresh = 0.975, delta = 0.2, flag_method = "prob")
  expect_true(all(s975$flag_bias == "no_clear_bias"))

  # Default flag_method = "bfdr" flags cells with bfdr <= 0.05 (here 0.04 <= 0.05)
  s_bfdr <- summarise_bias(obj, "item:rater", delta = 0.2)
  expect_equal(s_bfdr$flag_method[1], "bfdr")
  expect_equal(s_bfdr$flag_bias[1], "positive_bias")
})

test_that("facet_summary formats quantile column names without rounding collisions", {
  mat <- matrix(seq_len(200), 100, 2, dimnames = list(NULL, c("theta[1]", "theta[2]")))
  object <- structure(list(
    fit = list(summary = function(variables, ...) posterior::summarise_draws(posterior::as_draws_matrix(mat), ...)),
    spec = list(facets_main = "person"),
    data_info = list(facet_labels = list(person = c("p", "q")))
  ), class = "bmfrm_fit")
  x <- facet_summary(object, "person", probs = c(0.025, 0.5, 0.975))
  expect_true(all(c("q2.5", "q50", "q97.5") %in% names(x)))
  expect_false(any(c("q2", "q98") %in% names(x)))
})

test_that("print.residual_cor outputs real newlines and data_prep omits unused bias keys", {
  d <- expand.grid(person = 1:4, item = 1:3, rater = 1:2)
  d$resid <- sin(seq_len(nrow(d)))
  rc <- residual_cor_criteria(d, unit_facets = c("person", "rater"))
  txt <- paste(capture.output(print(rc)), collapse = "\n")
  expect_false(grepl("\\\\n", txt))

  d$score <- 1L
  spec <- parse_bmfrm_formula(score ~ person + item + rater + rater:item, d, K = 4)
  sd <- prepare_data_bmfrm(spec, d)$stan_data
  expect_false(any(grepl(":", names(sd), fixed = TRUE)))
  expect_false(any(grepl("^J_.*rateritem$", names(sd))))
})

test_that("save_gq = 'minimal' omits GQ and R on-demand engine reproduces exact moments and standardized fair scores", {
  # Unbalanced design: P1 takes easy item I1, P2 takes hard item I2, both have identical theta = 0
  d <- data.frame(
    person = c("P1", "P1", "P2", "P2"),
    item   = c("I1", "I1", "I2", "I2"),
    judge  = c("J1", "J2", "J1", "J2"),
    score  = c(3L, 4L, 1L, 2L),
    stringsAsFactors = FALSE
  )
  spec <- parse_bmfrm_formula(score ~ person + item + judge, d, K = 4, save_gq = "minimal")
  code <- build_stan_code(spec)
  expect_false(grepl("generated quantities", code, fixed = TRUE))

  # Construct deterministic 4-draw parameter posterior (no Stan GQ variables!)
  theta <- cbind(`theta[1]` = c(0, 0.2, -0.2, 0), `theta[2]` = c(0, 0.2, -0.2, 0))
  item  <- cbind(`item[1]` = rep(-1.0, 4), `item[2]` = rep(1.0, 4))
  judge <- cbind(`judge[1]` = rep(-0.5, 4), `judge[2]` = rep(0.5, 4))
  tau   <- cbind(`tau[1]` = rep(-0.8, 4), `tau[2]` = rep(0.0, 4), `tau[3]` = rep(0.8, 4))

  draws_map <- list(theta = theta, item = item, judge = judge, tau = tau)
  mock_fit <- list(
    metadata = function() list(stan_variables = names(draws_map)),
    draws = function(var, format = "matrix") draws_map[[var]]
  )
  di <- suppressWarnings(prepare_data_bmfrm(spec, d))
  obj <- structure(list(fit = mock_fit, spec = spec, data_info = di), class = "bmfrm_fit")

  # Warning when custom nuisance facet ('judge') is not matched by default
  expect_warning(fair_scores(obj), "No standard nuisance facet")

  # Explicit adjust_facets = "judge" with standardize = FALSE differs across persons due to item difficulty
  fs_obs <- fair_scores(obj, adjust_facets = "judge", standardize = FALSE)
  expect_gt(fs_obs$mean_fair_score[1], fs_obs$mean_fair_score[2] + 0.5)

  # With standardize = TRUE, both persons are evaluated on common reference grid {I1, I2}
  fs_std <- fair_scores(obj, adjust_facets = "judge", standardize = TRUE)
  expect_equal(fs_std$mean_fair_score[1], fs_std$mean_fair_score[2], tolerance = 1e-12)
  expect_equal(fs_std$sd_fair_score[1], fs_std$sd_fair_score[2], tolerance = 1e-12)
  expect_true(all(fs_std$standardized))

  # Verify on-demand residuals(), facet_fit(), and pp_check() all work with minimal GQ
  res <- residuals(obj)
  expect_equal(nrow(res), 4L)
  # Check exact first-observation expected score against hand calculation
  eta1 <- theta[, 1] - item[, 1] - judge[, 1]
  pr1 <- rsm_probabilities(eta1, tau[1, ])
  expect_equal(res$mu_hat[1], mean(as.vector(pr1 %*% 1:4)), tolerance = 1e-12)

  ff <- facet_fit(obj, "judge")
  expect_equal(nrow(ff), 2L)
  expect_true(all(is.finite(ff$infit_ppp) & ff$infit_ppp >= 0 & ff$infit_ppp <= 1))
  expect_false("infit_zstd" %in% names(ff))
  ff_z <- facet_fit(obj, "judge", include_zstd = TRUE)
  expect_true(all(is.finite(ff_z$infit_zstd)))
  expect_true(all(c("infit_q", "outfit_q") %in% names(ff_z)))

  ppc <- suppressMessages(pp_check(obj, type = "all", ndraws = 4, seed = 7))
  expect_equal(dim(ppc$yrep_subset), c(4L, 4L))
})

test_that("hierarchical main facets, sparse interaction priors, and hierarchical PCM simulate and generate valid Stan", {
  d <- expand.grid(person = paste0("P", 1:6), item = paste0("I", 1:4), rater = paste0("R", 1:3))
  d$score <- 1L

  # 1. Hierarchical person + rater + horseshoe interaction (horseshoe requires scale = "raw")
  f <- score ~ person + item + rater + item:rater
  h <- c("person", "rater", "item:rater")
  expect_error(
    prior("normal(0, 0.5)", "sd", facet = "item:rater", scale = "marginal", local_shrinkage = "horseshoe"),
    "global shrinkage scale"
  )
  ps <- c(
    prior("normal(0, 1)", "mean", facet = "person"),
    prior("lognormal(0, 0.3)", "sd", facet = "person"),
    prior("normal(0, 1)", "item", scale = "marginal"),
    prior("normal(0, 0.8)", "sd", facet = "rater", scale = "marginal"),
    prior("normal(0, 1)", "tau", scale = "marginal"),
    prior("normal(0, 0.5)", "sd", facet = "item:rater", scale = "raw", local_shrinkage = "horseshoe")
  )
  sim_hs <- simulate_bmfrm(f, d, K = 4, hierarchical = h, priors = ps, seed = 21)
  p1 <- sim_hs$parameters[[1L]]
  expect_true(all(c("sigma_rater", "sigma_itemrater", "itemrater") %in% names(p1)))
  expect_lt(abs(sum(p1$rater)), 1e-12)
  expect_lt(max(abs(c(rowSums(p1$itemrater), colSums(p1$itemrater)))), 1e-12)
  code_hs <- build_stan_code(sim_hs$spec)
  expect_match(code_hs, "sigma_rater_raw ~ normal(0, 0.8);", fixed = TRUE)
  expect_match(code_hs, "c2_itemrater ~ inv_gamma(hs_slab_shape_itemrater, hs_slab_rate_itemrater);", fixed = TRUE)
  expect_match(code_hs, "to_vector(lambda_itemrater) ~ cauchy(0, 1);", fixed = TRUE)

  # 2. Unit-variance Student-t local shrinkage interaction (scale = "marginal" matches SD = sigma, not sqrt(3)*sigma)
  ps_st <- ps
  ps_st[[6L]] <- prior(0.7, "sd", facet = "item:rater", scale = "marginal", local_shrinkage = "student_t")
  sim_st <- simulate_bmfrm(f, d, K = 4, hierarchical = h, priors = ps_st, seed = 22)
  expect_lt(max(abs(c(rowSums(sim_st$parameters[[1L]]$itemrater), colSums(sim_st$parameters[[1L]]$itemrater)))), 1e-12)
  expect_match(build_stan_code(sim_st$spec), "to_vector(lambda_itemrater) ~ inv_gamma(1.5, 0.5);", fixed = TRUE)
  pp_st <- prior_predict_bmfrm(f, d, K = 4, hierarchical = h, priors = ps_st, ndraws = 4000, seed = 24)
  b11 <- vapply(pp_st$parameters, function(p) p$itemrater[1, 1], numeric(1))
  # Conditional variance must match sigma^2 = 0.49 (not 3 * 0.49 = 1.47)
  expect_equal(stats::var(b11), 0.7^2, tolerance = 0.12)

  # 3. Hierarchical PCM (common steps + level-specific deviations)
  ps_pcm <- c(
    prior("normal(0, 1)", "mean", facet = "person"),
    prior("lognormal(0, 0.3)", "sd", facet = "person"),
    prior("normal(0, 1)", "item", scale = "marginal"),
    prior("normal(0, 0.7)", "sd", facet = "rater", scale = "marginal"),
    prior("normal(0, 1)", "tau", scale = "marginal"),
    prior("normal(0, 0.4)", "sd", facet = "tau", scale = "marginal")
  )
  sim_pcm <- simulate_bmfrm(
    score ~ person + item + rater, d, K = 4,
    family = "pcm", threshold_facet = "rater",
    hierarchical = c("person", "rater", "tau"),
    priors = ps_pcm, seed = 23
  )
  ppcm <- sim_pcm$parameters[[1L]]
  expect_true(all(c("tau_common", "sigma_tau", "tau") %in% names(ppcm)))
  expect_lt(abs(sum(ppcm$tau_common)), 1e-12)
  expect_lt(max(abs(rowSums(ppcm$tau))), 1e-12)
  code_pcm <- build_stan_code(sim_pcm$spec, save_gq = "full")
  expect_match(code_pcm, "tau_common", fixed = TRUE)
  expect_match(code_pcm, "sigma_tau_raw ~ normal(0, 0.4);", fixed = TRUE)
  # Stan full GQ must use cum_tau_fair (tau_common) when hierarchical PCM steps are on an adjusted facet
  expect_match(code_pcm, "cum_tau_fair = append_row(0.0, cumulative_sum(tau_common));", fixed = TRUE)

  # Non-adjusted threshold facet ("item") must keep item-specific cum_tau[threshold_n]
  sim_pcm_item <- simulate_bmfrm(
    score ~ person + item + rater, d, K = 4,
    family = "pcm", threshold_facet = "item",
    hierarchical = c("person", "tau"),
    priors = get_default_priors(c("person", "item", "rater"), family = "pcm", hierarchical = c("person", "tau")),
    seed = 25
  )
  expect_match(build_stan_code(sim_pcm_item$spec, save_gq = "full"), "eta_fair[n] - cum_tau[threshold_n]", fixed = TRUE)

  # Hierarchical PCM fair scores must be identical between full GQ, minimal GQ, and stale-Stan fallback
  di_pcm <- prepare_data_bmfrm(sim_pcm$spec, sim_pcm$data)
  tau_draws <- matrix(rep(as.vector(ppcm$tau), each = 4), 4)
  colnames(tau_draws) <- paste0("tau[", rep(1:3, 3), ",", rep(1:3, each = 3), "]")
  # Compute exact common-tau expected scores per observation
  obs <- di_pcm$data_clean
  eta_fair_true <- ppcm$theta[di_pcm$stan_data$person_id] - ppcm$item[di_pcm$stan_data$item_id]
  mu_fair_common <- as.vector(rsm_probabilities(eta_fair_true, ppcm$tau_common) %*% 1:4)
  mu_fair_wrong <- as.vector(pcm_probabilities(eta_fair_true, ppcm$tau[di_pcm$stan_data$rater_id, , drop = FALSE]) %*% 1:4)
  mu_fair_common_mat <- matrix(rep(mu_fair_common, each = 4), 4, dimnames = list(NULL, paste0("mu_fair[", seq_len(nrow(obs)), "]")))
  mu_fair_wrong_mat  <- matrix(rep(mu_fair_wrong,  each = 4), 4, dimnames = list(NULL, paste0("mu_fair[", seq_len(nrow(obs)), "]")))

  make_pcm_mock <- function(include_mu_fair, mu_fair_mat, stan_code_str) {
    vars <- c("theta", "item", "rater", "tau", "tau_common", if (include_mu_fair) "mu_fair")
    mock <- list(
      metadata = function() list(stan_variables = vars),
      draws = function(var, format = "matrix") {
        switch(var,
          theta = matrix(rep(ppcm$theta, each = 4), 4, dimnames = list(NULL, paste0("theta[", 1:6, "]"))),
          item  = matrix(rep(ppcm$item, each = 4), 4, dimnames = list(NULL, paste0("item[", 1:4, "]"))),
          rater = matrix(rep(ppcm$rater, each = 4), 4, dimnames = list(NULL, paste0("rater[", 1:3, "]"))),
          tau   = tau_draws,
          tau_common = matrix(rep(ppcm$tau_common, each = 4), 4, dimnames = list(NULL, paste0("tau_common[", 1:3, "]"))),
          mu_fair = mu_fair_mat
        )
      }
    )
    structure(list(fit = mock, spec = sim_pcm$spec, data_info = di_pcm, stan_code = stan_code_str), class = "bmfrm_fit")
  }

  fit_full    <- make_pcm_mock(TRUE,  mu_fair_common_mat, code_pcm)
  fit_minimal <- make_pcm_mock(FALSE, NULL, build_stan_code(sim_pcm$spec, save_gq = "minimal"))
  # Stale Stan code that used cum_tau[threshold_n] instead of cum_tau_fair must automatically fall back to R tau_common
  stale_code  <- gsub("cum_tau_fair", "cum_tau[threshold_n]", code_pcm, fixed = TRUE)
  stale_code  <- gsub("tau_common", "tau_unused", stale_code, fixed = TRUE)
  fit_stale   <- make_pcm_mock(TRUE,  mu_fair_wrong_mat, stale_code)

  expect_equal(suppressWarnings(fair_scores(fit_full, summary = FALSE)$mu_fair),
               suppressWarnings(fair_scores(fit_minimal, summary = FALSE)$mu_fair), tolerance = 1e-12)
  expect_equal(suppressWarnings(fair_scores(fit_full, summary = TRUE)$mean_fair_score),
               suppressWarnings(fair_scores(fit_minimal, summary = TRUE)$mean_fair_score), tolerance = 1e-12)
  expect_equal(suppressWarnings(fair_scores(fit_stale, summary = FALSE)$mu_fair),
               suppressWarnings(fair_scores(fit_minimal, summary = FALSE)$mu_fair), tolerance = 1e-12)
  expect_equal(suppressWarnings(fair_scores(fit_full, standardize = TRUE)$mean_fair_score),
               suppressWarnings(fair_scores(fit_minimal, standardize = TRUE)$mean_fair_score), tolerance = 1e-12)

  # Hierarchical PCM fair scores must NOT claim conditional-likelihood identification:
  # tau_common vs tau_dev split relies on the hierarchical step prior.
  fit_full$design <- check_design_bmfrm(sim_pcm$spec$formula, sim_pcm$data, K = 4,
                                        family = "pcm", threshold_facet = "rater")
  expect_true(fit_full$design$identifiable)
  expect_warning(fs_hpcm_obs <- fair_scores(fit_full, summary = FALSE), "hierarchical step prior")
  expect_true(all(!fs_hpcm_obs$likelihood_identified))
  expect_true(all(fs_hpcm_obs$identification_basis == "hierarchical_step_prior"))
  expect_warning(fs_hpcm_std <- fair_scores(fit_full, standardize = TRUE), "hierarchical step prior")
  expect_true(all(!fs_hpcm_std$likelihood_identified))
  expect_true(all(fs_hpcm_std$identification_basis == "hierarchical_step_prior"))
})

test_that("chunked streaming matches unchunked results identically and max_grid_cells guards auto-expansion", {
  d <- expand.grid(
    person = paste0("P", 1:4),
    item   = paste0("I", 1:3),
    task   = paste0("T", 1:2),
    rater  = paste0("R", 1:2),
    stringsAsFactors = FALSE
  )
  d$score <- (seq_len(nrow(d)) - 1L) %% 4L + 1L
  spec <- parse_bmfrm_formula(score ~ person + item + task + rater, d, K = 4, save_gq = "minimal")
  di <- prepare_data_bmfrm(spec, d)
  set.seed(77)
  D <- 12L
  draws_map <- list(
    theta = matrix(stats::rnorm(D * 4), D, 4, dimnames = list(NULL, paste0("theta[", 1:4, "]"))),
    item  = matrix(stats::rnorm(D * 3), D, 3, dimnames = list(NULL, paste0("item[", 1:3, "]"))),
    task  = matrix(stats::rnorm(D * 2), D, 2, dimnames = list(NULL, paste0("task[", 1:2, "]"))),
    rater = matrix(stats::rnorm(D * 2), D, 2, dimnames = list(NULL, paste0("rater[", 1:2, "]"))),
    tau   = matrix(stats::rnorm(D * 3), D, 3, dimnames = list(NULL, paste0("tau[", 1:3, "]")))
  )
  mock <- list(
    metadata = function() list(stan_variables = names(draws_map)),
    draws = function(var, format = "matrix") draws_map[[var]]
  )
  obj <- structure(list(fit = mock, spec = spec, data_info = di), class = "bmfrm_fit")

  # Chunked (chunk_size = 5) vs single-pass (chunk_size = 1000) must match to machine precision
  q_small <- compute_bmfrm_quantities(obj, compute_kurtosis = TRUE, compute_log_lik = TRUE, chunk_size = 5L)
  q_large <- compute_bmfrm_quantities(obj, compute_kurtosis = TRUE, compute_log_lik = TRUE, chunk_size = 1000L)
  expect_equal(q_small$mu, q_large$mu, tolerance = 1e-13)
  expect_equal(q_small$sigma2, q_large$sigma2, tolerance = 1e-13)
  expect_equal(q_small$kurt, q_large$kurt, tolerance = 1e-13)
  expect_equal(q_small$log_lik, q_large$log_lik, tolerance = 1e-13)

  ff_small <- facet_fit(obj, "rater", seed = 99L, include_zstd = TRUE, chunk_size = 5L)
  ff_large <- facet_fit(obj, "rater", seed = 99L, include_zstd = TRUE, chunk_size = 1000L)
  expect_equal(ff_small$infit_msq, ff_large$infit_msq, tolerance = 1e-13)
  expect_equal(ff_small$infit_ppp, ff_large$infit_ppp, tolerance = 1e-13)
  expect_equal(ff_small$outfit_zstd, ff_large$outfit_zstd, tolerance = 1e-13)

  fs_small <- fair_scores(obj, standardize = TRUE, chunk_size = 5L)
  fs_large <- fair_scores(obj, standardize = TRUE, chunk_size = 1000L)
  expect_equal(fs_small$mean_fair_score, fs_large$mean_fair_score, tolerance = 1e-13)
  expect_equal(fs_small$sd_fair_score, fs_large$sd_fair_score, tolerance = 1e-13)

  # Grid cap blocks combinatorial explosion (item x task = 6 cells > max_grid_cells = 4)
  expect_error(fair_scores(obj, standardize = TRUE, max_grid_cells = 4L), "max_grid_cells")
})

test_that("facet_fit falls back safely to saved mu/sigma2 for bundled hardcoded 4-facet models", {
  d <- expand.grid(person = paste0("P", 1:4), item = paste0("I", 1:2),
                   rater = c("R1", "R2"), scenario = c("S1", "S2"), stringsAsFactors = FALSE)
  d$score <- (seq_len(nrow(d)) - 1L) %% 4L + 1L
  spec <- parse_bmfrm_formula(score ~ person + item + rater + scenario, d, K = 4)
  spec$custom_stan <- TRUE
  di <- prepare_data_bmfrm(spec, d)
  N <- nrow(d); D <- 10L
  set.seed(88)
  mu_mat <- matrix(2.5 + stats::rnorm(D * N, 0, 0.1), D, N, dimnames = list(NULL, paste0("mu[", seq_len(N), "]")))
  sig_mat <- matrix(0.9 + stats::runif(D * N, 0, 0.1), D, N, dimnames = list(NULL, paste0("sigma2[", seq_len(N), "]")))
  # Mimic bundled four-facet/person-item-rater-scenario-main.stan parameter names (beta/gamma/delta, NOT item/rater/scenario)
  mock_bundled <- list(
    metadata = function() list(stan_variables = c("theta", "beta", "gamma", "delta", "tau", "mu", "sigma2", "mu_fair")),
    summary = function(var, ...) {
      m <- if (var == "mu") mu_mat else sig_mat
      tibble::tibble(variable = colnames(m), mean = colMeans(m))
    },
    draws = function(var, format = "matrix") {
      if (var == "mu") return(mu_mat)
      if (var == "sigma2") return(sig_mat)
      stop("Missing parameter: ", var)
    }
  )
  obj <- structure(list(fit = mock_bundled, spec = spec, data_info = di), class = "bmfrm_fit")
  expect_false(can_evaluate_in_r(obj))
  # Even if custom_stan were FALSE, missing item/rater/scenario must still prevent R-engine dispatch
  obj_flag_false <- obj; obj_flag_false$spec$custom_stan <- FALSE
  expect_false(can_evaluate_in_r(obj_flag_false))

  expect_equal(nrow(residuals(obj)), N)
  ff_def <- facet_fit(obj, "rater")
  expect_false("infit_zstd" %in% names(ff_def))
  ff <- facet_fit(obj, "rater", include_zstd = TRUE)
  expect_equal(nrow(ff), 2L)
  expect_true(all(is.finite(ff$infit_zstd)))
  expect_true(all(ff$zstd_method == "wilson_hilferty_chisq_approx"))
  expect_equal(nrow(facet_fit(obj_flag_false, "rater")), 2L)
})

test_that("summarise_bias BFDR is invariant to cell/label ordering under tied error probabilities", {
  d <- expand.grid(person = paste0("P", 1:6), item = paste0("I", 1:3), rater = c("R1", "R2"))
  d$score <- 1L
  spec <- parse_bmfrm_formula(score ~ person + item + rater + item:rater, d, K = 4)
  ch <- check_design_bmfrm(spec$formula, d, K = 4)

  # 6 cells: 2 cells with pep = 0.00, 4 tied cells with pep = 0.10
  # Expected tie-group BFDR for all 4 tied cells = (0 + 0 + 4 * 0.10) / 6 = 0.40 / 6
  D <- 100L
  zero_pep <- rep(0.5, D)
  tied_pep <- c(rep(0.0, 10L), rep(0.5, 90L))
  mat <- cbind(zero_pep, zero_pep, tied_pep, tied_pep, tied_pep, tied_pep)
  vars <- paste0("itemrater[", rep(1:3, 2), ",", rep(1:2, each = 3), "]")
  colnames(mat) <- vars

  make_obj <- function(m) {
    structure(list(
      spec = spec, design = ch,
      data_info = list(facet_labels = list(person = paste0("P", 1:6), item = paste0("I", 1:3), rater = c("R1", "R2")),
                       data_clean = d),
      fit = list(draws = function(...) m)
    ), class = "bmfrm_fit")
  }
  s1 <- summarise_bias(make_obj(mat), "item:rater", prob_thresh = 0.95, delta = 0.2)
  expect_equal(s1$bfdr[1:2], c(0, 0), tolerance = 1e-12)
  expect_equal(s1$bfdr[3:6], rep(0.4 / 6, 4), tolerance = 1e-12)

  # Permute columns so tied cells appear in reversed order: BFDR values must remain identical per cell
  perm <- c(6, 5, 4, 3, 2, 1)
  mat_perm <- mat[, perm]
  colnames(mat_perm) <- vars
  s2 <- summarise_bias(make_obj(mat_perm), "item:rater", prob_thresh = 0.95, delta = 0.2)
  expect_equal(s2$bfdr[1:4], rep(0.4 / 6, 4), tolerance = 1e-12)
  expect_equal(s2$bfdr[5:6], c(0, 0), tolerance = 1e-12)
})

test_that("fair_scores propagates design likelihood identifiability and warns on rank-deficient designs", {
  # Connected but rank-deficient design: item and rater are 1-to-1 confounded (I1<->R1, I2<->R2)
  d_bad <- data.frame(
    person = rep(c("P1", "P2"), each = 2),
    item   = rep(c("I1", "I2"), 2),
    rater  = rep(c("R1", "R2"), 2),
    score  = c(2L, 3L, 2L, 3L),
    stringsAsFactors = FALSE
  )
  spec_bad <- parse_bmfrm_formula(score ~ person + item + rater, d_bad, K = 4, save_gq = "minimal")
  di_bad <- suppressWarnings(prepare_data_bmfrm(spec_bad, d_bad))
  dc_bad <- check_design_bmfrm(spec_bad$formula, d_bad, K = 4)
  expect_true(dc_bad$connected)
  expect_false(dc_bad$identifiable)

  # Identical likelihood from two parameter sets shifted along the confounded item/rater null vector
  draws_bad <- list(
    theta = matrix(0, 4, 2, dimnames = list(NULL, c("theta[1]", "theta[2]"))),
    item  = matrix(rep(c(-0.5, 0.5), each = 4), 4, 2, dimnames = list(NULL, c("item[1]", "item[2]"))),
    rater = matrix(rep(c(0.5, -0.5), each = 4), 4, 2, dimnames = list(NULL, c("rater[1]", "rater[2]"))),
    tau   = matrix(rep(c(-0.5, 0, 0.5), each = 4), 4, 3, dimnames = list(NULL, paste0("tau[", 1:3, "]")))
  )
  mock_bad <- list(
    metadata = function() list(stan_variables = names(draws_bad)),
    draws = function(var, format = "matrix") draws_bad[[var]]
  )
  obj_bad <- structure(list(fit = mock_bad, spec = spec_bad, data_info = di_bad, design = dc_bad), class = "bmfrm_fit")

  expect_warning(fs_bad_obs <- fair_scores(obj_bad, summary = FALSE), "likelihood-unidentified")
  expect_true(all(!fs_bad_obs$likelihood_identified))
  expect_warning(fs_bad_std <- fair_scores(obj_bad, standardize = TRUE), "likelihood-unidentified")
  expect_true(all(!fs_bad_std$likelihood_identified))

  # Full-rank design must report likelihood_identified = TRUE without warning
  d_ok <- expand.grid(person = c("P1", "P2"), item = c("I1", "I2"), rater = c("R1", "R2"), stringsAsFactors = FALSE)
  d_ok$score <- c(2L, 3L, 2L, 3L, 3L, 2L, 3L, 2L)
  spec_ok <- parse_bmfrm_formula(score ~ person + item + rater, d_ok, K = 4, save_gq = "minimal")
  di_ok <- suppressWarnings(prepare_data_bmfrm(spec_ok, d_ok))
  dc_ok <- check_design_bmfrm(spec_ok$formula, d_ok, K = 4)
  expect_true(dc_ok$identifiable)
  obj_ok <- structure(list(fit = mock_bad, spec = spec_ok, data_info = di_ok, design = dc_ok), class = "bmfrm_fit")
  expect_no_warning(fs_ok <- fair_scores(obj_ok, standardize = TRUE))
  expect_true(all(fs_ok$likelihood_identified))
})

test_that("draw-level kurtosis q_d avoids swapping E[C]/E[W]^2 with E[C/W^2]", {
  # Two draws with markedly different variance W and fourth moment C
  x <- c(1L, 4L)
  fid <- c(1L, 1L)
  mu_mat <- rbind(c(2.5, 2.5), c(2.5, 2.5))
  sig_mat <- rbind(c(0.2, 0.2), c(2.0, 2.0))
  kurt_mat <- rbind(c(1.0, 1.0), c(6.0, 6.0))

  ff <- facet_fit_from_matrices(x, fid, mu_mat, sig_mat, kurt_mat = kurt_mat, model_name = "m", include_zstd = TRUE)
  # Per-draw q_outfit:
  # draw 1: C/W^2 - 1 = 1/0.04 - 1 = 24 -> sum=48 -> sqrt(48/4) = sqrt(12)
  # draw 2: C/W^2 - 1 = 6/4 - 1 = 0.5   -> sum=1  -> sqrt(1/4)  = 0.5
  q_d <- c(sqrt(12), 0.5)
  out_d <- c(mean((x - 2.5)^2 / 0.2), mean((x - 2.5)^2 / 2.0))
  z_d <- (out_d^(1/3) - 1) * (3 / q_d) + q_d / 3
  expect_equal(ff$outfit_q, mean(q_d), tolerance = 1e-12)
  expect_equal(ff$outfit_zstd, mean(z_d), tolerance = 1e-12)
  expect_equal(ff$zstd_method, "draw_wilson_hilferty_approx")
  # Confirm this differs markedly from the swapped E[C]/E[W]^2 formula
  q_swapped <- sqrt(2 * (3.5 / (1.1^2) - 1) / 4)
  expect_gt(abs(ff$outfit_q - q_swapped), 0.5)
})

test_that("bundled 3-facet and 4-facet models honor requested adjust_facets instead of misusing built-in mu_fair", {
  # Bundled 3-facet model removes BOTH rater and scenario in its built-in mu_fair
  m3 <- bmfrm_model("three-facet/person-rater-scenario-main")
  d3 <- data.frame(
    person   = c("P1", "P1"),
    rater    = c("R1", "R2"),
    scenario = c("S1", "S2"),
    score    = c(3L, 2L),
    stringsAsFactors = FALSE
  )
  spec3 <- parse_bmfrm_formula(score ~ person + rater + scenario, d3, K = 4)
  spec3$custom_stan <- TRUE
  di3 <- suppressWarnings(prepare_data_bmfrm(spec3, d3))

  theta_v <- 0.5
  rater_v <- c(-0.3, 0.3)
  scen_v  <- c(-0.8, 0.8)
  tau_v   <- c(-0.6, 0.0, 0.6)
  # Built-in Stan mu_fair sets rater = 0 AND scenario = 0 (equal for both rows)
  mu_fair_both <- as.vector(rsm_probabilities(c(theta_v, theta_v), tau_v) %*% 1:4)
  # Requested adjust_facets = "rater" keeps scenario in eta_fair (distinct across rows!)
  mu_fair_rater_only <- as.vector(rsm_probabilities(theta_v - scen_v, tau_v) %*% 1:4)
  expect_gt(abs(mu_fair_rater_only[1] - mu_fair_rater_only[2]), 0.5)

  D <- 4L
  draws3 <- list(
    theta    = matrix(theta_v, D, 1, dimnames = list(NULL, "theta[1]")),
    rater    = matrix(rep(rater_v, each = D), D, 2, dimnames = list(NULL, c("rater[1]", "rater[2]"))),
    scenario = matrix(rep(scen_v,  each = D), D, 2, dimnames = list(NULL, c("scenario[1]", "scenario[2]"))),
    tau      = matrix(rep(tau_v,   each = D), D, 3, dimnames = list(NULL, paste0("tau[", 1:3, "]"))),
    mu_fair  = matrix(rep(mu_fair_both, each = D), D, 2, dimnames = list(NULL, c("mu_fair[1]", "mu_fair[2]")))
  )
  mock3 <- list(
    metadata = function() list(stan_variables = names(draws3)),
    draws = function(var, format = "matrix") draws3[[var]]
  )
  obj3 <- structure(list(fit = mock3, spec = spec3, data_info = di3, stan_code = m3$code), class = "bmfrm_fit")

  expect_setequal(stan_mu_fair_adjusted_facets(obj3), c("rater", "scenario"))
  # Explicit adjust_facets = "rater" must NOT reuse Stan's rater+scenario mu_fair; it keeps scenario!
  fs_r <- fair_scores(obj3, adjust_facets = "rater", summary = FALSE)
  expect_equal(fs_r$mu_fair, mu_fair_rater_only, tolerance = 1e-12)
  # Default adjust_facets = NULL uses Stan's built-in mu_fair with an informative message
  expect_message(fs_def <- fair_scores(obj3, summary = FALSE), "Using custom Stan built-in mu_fair")
  expect_equal(fs_def$mu_fair, mu_fair_both, tolerance = 1e-12)

  # Bundled 4-facet model aliases item -> beta and removes rater + scenario in mu_fair:
  # requesting adjust_facets = "rater" cannot reuse mu_fair and cannot remap beta -> must error explicitly
  m4 <- bmfrm_model("four-facet/person-item-rater-scenario-main")
  d4 <- d3; d4$item <- c("I1", "I2")
  spec4 <- parse_bmfrm_formula(score ~ person + item + rater + scenario, d4, K = 4)
  spec4$custom_stan <- TRUE
  di4 <- suppressWarnings(prepare_data_bmfrm(spec4, d4))
  mock4 <- list(
    metadata = function() list(stan_variables = c("theta", "beta", "rater", "scenario", "tau", "mu_fair")),
    draws = function(var, format = "matrix") draws3$mu_fair
  )
  obj4 <- structure(list(fit = mock4, spec = spec4, data_info = di4, stan_code = m4$code), class = "bmfrm_fit")
  expect_setequal(stan_mu_fair_adjusted_facets(obj4), c("rater", "scenario"))
  expect_error(fair_scores(obj4, adjust_facets = "rater"), "does not match the requested adjust_facets")
})

test_that("Schur-complement design check scales to 1,500 persons (6,000 obs) and matches dense SVD", {
  # 1,500 persons x 2 items x 2 raters = 6,000 rows (n * p = 9,012,000 > default max_cells = 2e6)
  d_big <- expand.grid(person = seq_len(1500L), item = 1:2, rater = 1:2)
  dc_big <- check_design_bmfrm(score ~ person + item + rater, d_big, K = 4)
  expect_equal(dc_big$solver, "schur_complement")
  expect_true(dc_big$connected)
  expect_true(dc_big$identifiable)
  expect_equal(dc_big$rank, 1502L)
  expect_equal(dc_big$nullity, 0L)

  # Confounded 1,500-person design (item == rater): connected but rank-deficient (nullity = 1)
  d_bad <- d_big[d_big$item == d_big$rater, , drop = FALSE]
  dc_bad <- check_design_bmfrm(score ~ person + item + rater, d_bad, K = 4)
  expect_equal(dc_bad$solver, "schur_complement")
  expect_true(dc_bad$connected)
  expect_false(dc_bad$identifiable)
  expect_equal(dc_bad$rank, 1501L)
  expect_equal(dc_bad$nullity, 1L)
  expect_equal(ncol(dc_bad$effect_null_space), 1L)

  # Verify dense SVD and Schur solver agree exactly on an interaction design with unobserved cell
  d_int <- expand.grid(person = paste0("P", 1:6), item = paste0("I", 1:3), rater = c("A", "B"))
  d_int <- subset(d_int, !(item == "I3" & rater == "B"))
  ch_dense <- check_design_bmfrm(score ~ person + item + rater + item:rater, d_int, K = 4, method = "dense")
  ch_schur <- check_design_bmfrm(score ~ person + item + rater + item:rater, d_int, K = 4, method = "schur")
  expect_equal(ch_schur$rank, ch_dense$rank)
  expect_equal(ch_schur$nullity, ch_dense$nullity)
  expect_equal(ch_schur$identifiable, ch_dense$identifiable)
  # Null vector from Schur solver must satisfy X %*% v == 0 in dense design matrix
  v_dense_order <- ch_schur$null_space[colnames(ch_dense$design_matrix), , drop = FALSE]
  expect_lt(max(abs(ch_dense$design_matrix %*% v_dense_order)), 1e-11)

  # max_output_rows guards standardized row-level expansion (P * n_ref)
  d_small <- expand.grid(person = paste0("P", 1:4), item = paste0("I", 1:3), rater = c("R1", "R2"), stringsAsFactors = FALSE)
  d_small$score <- 1L
  spec_s <- parse_bmfrm_formula(score ~ person + item + rater, d_small, K = 4, save_gq = "minimal")
  di_s <- prepare_data_bmfrm(spec_s, d_small)
  mock_s <- list(
    metadata = function() list(stan_variables = c("theta", "item", "rater", "tau")),
    draws = function(var, format = "matrix") matrix(0, 2, switch(var, theta = 4, item = 3, rater = 2, tau = 3))
  )
  obj_s <- structure(list(fit = mock_s, spec = spec_s, data_info = di_s), class = "bmfrm_fit")
  expect_error(fair_scores(obj_s, standardize = TRUE, summary = FALSE, max_output_rows = 10L), "max_output_rows")

  # Condition number from Schur solver must match dense SVD exactly (not Schur-reduced singular values)
  set.seed(123)
  d_unbal <- rbind(
    expand.grid(person = paste0("P", 1:4), item = paste0("I", 1:3), rater = paste0("R", 1:3)),
    expand.grid(person = paste0("P", 5:10), item = paste0("I", 1:2), rater = paste0("R", 1:2))
  )
  cd <- check_design_bmfrm(score ~ person + item + rater, d_unbal, K = 4, method = "dense")
  cs <- check_design_bmfrm(score ~ person + item + rater, d_unbal, K = 4, method = "schur")
  expect_equal(cs$rank, cd$rank)
  expect_equal(cs$condition_number, cd$condition_number, tolerance = 1e-10)

  # Also check exact condition number match with an interaction term
  d_int_full <- expand.grid(person = paste0("P", 1:6), item = paste0("I", 1:3), rater = c("A", "B"))
  d_int_full <- d_int_full[-c(1L, 7L, 14L), , drop = FALSE] # mild imbalance, still full rank
  cd2 <- check_design_bmfrm(score ~ person + item + rater + item:rater, d_int_full, K = 4, method = "dense")
  cs2 <- check_design_bmfrm(score ~ person + item + rater + item:rater, d_int_full, K = 4, method = "schur")
  expect_true(cd2$identifiable && cs2$identifiable)
  expect_equal(cs2$condition_number, cd2$condition_number, tolerance = 1e-10)
})

test_that("person-level PCM (threshold_facet = 'person') supports standardized fair_scores", {
  d <- expand.grid(person = c("P1", "P2", "P3"), item = c("I1", "I2"), rater = c("R1", "R2"), stringsAsFactors = FALSE)
  d$score <- c(1L, 2L, 3L, 4L, 2L, 3L, 1L, 4L, 2L, 3L, 2L, 4L)
  spec <- parse_bmfrm_formula(score ~ person + item + rater, d, K = 4,
                              family = "pcm", threshold_facet = "person", save_gq = "minimal")
  di <- prepare_data_bmfrm(spec, d)
  dc <- check_design_bmfrm(spec$formula, d, K = 4, family = "pcm", threshold_facet = "person")

  theta_v <- c(-0.4, 0.1, 0.7)
  item_v  <- c(-0.6, 0.6)
  rater_v <- c(-0.3, 0.3)
  tau_m   <- rbind(c(-0.8, 0.0, 0.8), c(-0.3, 0.1, 0.2), c(-1.0, 0.2, 0.8)) # 3 persons x 3 steps
  D <- 4L
  tau_draws <- matrix(rep(as.vector(tau_m), each = D), D)
  colnames(tau_draws) <- paste0("tau[", rep(1:3, 3), ",", rep(1:3, each = 3), "]")
  draws_p <- list(
    theta = matrix(rep(theta_v, each = D), D, 3, dimnames = list(NULL, paste0("theta[", 1:3, "]"))),
    item  = matrix(rep(item_v,  each = D), D, 2, dimnames = list(NULL, paste0("item[", 1:2, "]"))),
    rater = matrix(rep(rater_v, each = D), D, 2, dimnames = list(NULL, paste0("rater[", 1:2, "]"))),
    tau   = tau_draws
  )
  mock_p <- list(
    metadata = function() list(stan_variables = names(draws_p)),
    draws = function(var, format = "matrix") draws_p[[var]]
  )
  obj <- structure(list(fit = mock_p, spec = spec, data_info = di, design = dc), class = "bmfrm_fit")

  fs_std <- fair_scores(obj, standardize = TRUE)
  expect_equal(nrow(fs_std), 3L)
  expect_true(all(fs_std$likelihood_identified))
  # Verify exact match with hand calculation using each person's own step profile tau_m[p, ]
  expected_fair <- vapply(1:3, function(p) {
    mean(as.vector(pcm_probabilities(theta_v[p] - item_v, tau_m[c(p, p), , drop = FALSE]) %*% 1:4))
  }, numeric(1))
  expect_equal(fs_std$mean_fair_score, expected_fair, tolerance = 1e-12)

  fs_rows <- fair_scores(obj, standardize = TRUE, summary = FALSE)
  expect_equal(nrow(fs_rows), 6L)
  expect_true(all(is.finite(fs_rows$mu_fair)))
})

test_that("standardized fair_scores supports interactions involving person in either facet order", {
  d <- expand.grid(person = c("P1", "P2", "P3"), item = c("I1", "I2"), rater = c("R1", "R2"), stringsAsFactors = FALSE)
  d$score <- c(1L, 2L, 3L, 4L, 2L, 3L, 1L, 4L, 2L, 3L, 2L, 4L)

  theta_v <- c(-0.5, 0.0, 0.6)
  item_v  <- c(-0.4, 0.4)
  rater_v <- c(-0.3, 0.3)
  tau_v   <- c(-0.7, 0.0, 0.7)
  # Double-centered 3 x 2 person:item interaction matrix
  B_pi <- rbind(c(0.5, -0.5), c(-0.2, 0.2), c(-0.3, 0.3))
  D <- 4L

  for (order_kind in c("person:item", "item:person")) {
    form <- if (order_kind == "person:item") {
      score ~ person + item + rater + person:item
    } else {
      # Putting item before person keeps formula interaction order as "item:person"
      score ~ item + person + rater + item:person
    }
    spec <- parse_bmfrm_formula(form, d, K = 4, save_gq = "minimal")
    expect_equal(spec$facets_bias, order_kind)
    di <- prepare_data_bmfrm(spec, d)
    dc <- check_design_bmfrm(spec$formula, d, K = 4)

    B_mat <- if (order_kind == "person:item") B_pi else t(B_pi)
    b_var <- gsub(":", "", order_kind)
    nr <- nrow(B_mat); nc <- ncol(B_mat)
    b_draws <- matrix(rep(as.vector(B_mat), each = D), D)
    colnames(b_draws) <- paste0(b_var, "[", rep(seq_len(nr), nc), ",", rep(seq_len(nc), each = nr), "]")

    draws_map <- list(
      theta = matrix(rep(theta_v, each = D), D, 3, dimnames = list(NULL, paste0("theta[", 1:3, "]"))),
      item  = matrix(rep(item_v,  each = D), D, 2, dimnames = list(NULL, paste0("item[", 1:2, "]"))),
      rater = matrix(rep(rater_v, each = D), D, 2, dimnames = list(NULL, paste0("rater[", 1:2, "]"))),
      tau   = matrix(rep(tau_v,   each = D), D, 3, dimnames = list(NULL, paste0("tau[", 1:3, "]")))
    )
    draws_map[[b_var]] <- b_draws
    mock <- list(
      metadata = function() list(stan_variables = names(draws_map)),
      draws = function(var, format = "matrix") draws_map[[var]]
    )
    obj <- structure(list(fit = mock, spec = spec, data_info = di, design = dc), class = "bmfrm_fit")

    fs_obs <- fair_scores(obj, standardize = FALSE)
    fs_std <- fair_scores(obj, standardize = TRUE)
    fs_std_chunked <- fair_scores(obj, standardize = TRUE, chunk_size = 2L)
    fs_rows <- fair_scores(obj, standardize = TRUE, summary = FALSE, chunk_size = 2L)

    expected_fair <- vapply(1:3, function(p) {
      mean(as.vector(rsm_probabilities(theta_v[p] - item_v + B_pi[p, ], tau_v) %*% 1:4))
    }, numeric(1))

    expect_equal(fs_std$mean_fair_score, expected_fair, tolerance = 1e-12)
    expect_equal(fs_std$mean_fair_score, fs_obs$mean_fair_score, tolerance = 1e-12)
    expect_equal(fs_std_chunked$mean_fair_score, expected_fair, tolerance = 1e-12)
    expect_true(all(fs_std$likelihood_identified))
    expect_equal(nrow(fs_rows), 6L)
    expect_equal(
      fs_rows$mu_fair,
      as.vector(vapply(1:3, function(p) as.vector(rsm_probabilities(theta_v[p] - item_v + B_pi[p, ], tau_v) %*% 1:4), numeric(2))),
      tolerance = 1e-12
    )
  }

  # Missing (P3, I2) cell creates a person:item null-space direction; check identifiability without crashing
  d_miss <- subset(d, !(person == "P3" & item == "I2"))
  spec_m <- parse_bmfrm_formula(score ~ person + item + rater + person:item, d_miss, K = 4, save_gq = "minimal")
  di_m <- prepare_data_bmfrm(spec_m, d_miss)
  dc_m <- check_design_bmfrm(spec_m$formula, d_miss, K = 4)
  expect_false(dc_m$identifiable)
  b_draws_m <- matrix(rep(as.vector(B_pi), each = D), D)
  colnames(b_draws_m) <- paste0("personitem[", rep(1:3, 2), ",", rep(1:2, each = 3), "]")
  draws_m <- list(
    theta = matrix(rep(theta_v, each = D), D, 3, dimnames = list(NULL, paste0("theta[", 1:3, "]"))),
    item  = matrix(rep(item_v,  each = D), D, 2, dimnames = list(NULL, paste0("item[", 1:2, "]"))),
    rater = matrix(rep(rater_v, each = D), D, 2, dimnames = list(NULL, paste0("rater[", 1:2, "]"))),
    tau   = matrix(rep(tau_v,   each = D), D, 3, dimnames = list(NULL, paste0("tau[", 1:3, "]"))),
    personitem = b_draws_m
  )
  obj_m <- structure(
    list(fit = list(metadata = function() list(stan_variables = names(draws_m)),
                    draws = function(var, format = "matrix") draws_m[[var]]),
         spec = spec_m, data_info = di_m, design = dc_m),
    class = "bmfrm_fit"
  )
  expect_warning(fs_miss <- fair_scores(obj_m, standardize = TRUE), "likelihood-unidentified")
  expect_false(fs_miss$likelihood_identified[fs_miss$person == "P3"])
})

test_that("binary hierarchical PCM (K = 2) is not falsely flagged as hierarchical_step_prior", {
  d <- expand.grid(person = paste0("P", 1:4), item = c("I1", "I2"), rater = c("R1", "R2"), stringsAsFactors = FALSE)
  d$score <- rep(c(1L, 2L), length.out = nrow(d))
  ps <- get_default_priors(c("person", "item", "rater"), family = "pcm", hierarchical = c("person", "tau"))
  sim_bin <- simulate_bmfrm(
    score ~ person + item + rater, d, K = 2,
    family = "pcm", threshold_facet = "rater",
    hierarchical = c("person", "tau"), priors = ps, seed = 41
  )
  di <- prepare_data_bmfrm(sim_bin$spec, sim_bin$data)
  dc <- check_design_bmfrm(sim_bin$spec$formula, sim_bin$data, K = 2,
                           family = "pcm", threshold_facet = "rater")
  expect_true(dc$identifiable)

  p1 <- sim_bin$parameters[[1L]]
  D <- 4L
  tau_draws <- matrix(rep(as.vector(p1$tau), each = D), D)
  colnames(tau_draws) <- paste0("tau[", 1:2, ",1]")
  draws_bin <- list(
    theta = matrix(rep(p1$theta, each = D), D, 4, dimnames = list(NULL, paste0("theta[", 1:4, "]"))),
    item  = matrix(rep(p1$item,  each = D), D, 2, dimnames = list(NULL, paste0("item[", 1:2, "]"))),
    rater = matrix(rep(p1$rater, each = D), D, 2, dimnames = list(NULL, paste0("rater[", 1:2, "]"))),
    tau   = tau_draws,
    tau_common = matrix(0, D, 1, dimnames = list(NULL, "tau_common[1]"))
  )
  obj_bin <- structure(
    list(fit = list(metadata = function() list(stan_variables = names(draws_bin)),
                    draws = function(var, format = "matrix") draws_bin[[var]]),
         spec = sim_bin$spec, data_info = di, design = dc),
    class = "bmfrm_fit"
  )

  expect_no_warning(fs_obs <- fair_scores(obj_bin, summary = FALSE))
  expect_true(all(fs_obs$likelihood_identified))
  expect_true(all(fs_obs$identification_basis == "conditional_likelihood"))

  expect_no_warning(fs_std <- fair_scores(obj_bin, standardize = TRUE))
  expect_true(all(fs_std$likelihood_identified))
  expect_true(all(fs_std$identification_basis == "conditional_likelihood"))
})

test_that("unverified custom Stan models cannot silently dispatch to R recomputation even with matching parameter names", {
  d <- expand.grid(person = c("P1", "P2"), item = c("I1", "I2"), rater = c("R1", "R2"), stringsAsFactors = FALSE)
  d$score <- c(1L, 2L, 3L, 4L, 2L, 3L, 1L, 4L)
  spec <- parse_bmfrm_formula(score ~ person + item + rater, d, K = 4)
  spec$custom_stan <- TRUE
  di <- suppressWarnings(prepare_data_bmfrm(spec, d))
  dc <- check_design_bmfrm(spec$formula, d, K = 4)

  # Custom Stan code with modified linear predictor (2 * theta) but identical variable names
  custom_code <- "
  model {
    for (n in 1:N) {
      real eta = 2 * theta[person_id[n]] - item[item_id[n]] - rater[rater_id[n]];
    }
  }
  generated quantities {
    for (n in 1:N) {
      real eta_fair = 2 * theta[person_id[n]] - item[item_id[n]];
    }
  }"

  theta_v <- c(0.8, -0.4)
  item_v  <- c(-0.3, 0.3)
  rater_v <- c(-0.2, 0.2)
  tau_v   <- c(-0.5, 0.0, 0.5)
  eta_fair_2theta <- 2 * theta_v[di$stan_data$person_id] - item_v[di$stan_data$item_id]
  mu_fair_2theta <- as.vector(rsm_probabilities(eta_fair_2theta, tau_v) %*% 1:4)
  D <- 4L
  draws_custom <- list(
    theta = matrix(rep(theta_v, each = D), D, 2, dimnames = list(NULL, c("theta[1]", "theta[2]"))),
    item  = matrix(rep(item_v,  each = D), D, 2, dimnames = list(NULL, c("item[1]", "item[2]"))),
    rater = matrix(rep(rater_v, each = D), D, 2, dimnames = list(NULL, c("rater[1]", "rater[2]"))),
    tau   = matrix(rep(tau_v,   each = D), D, 3, dimnames = list(NULL, paste0("tau[", 1:3, "]"))),
    mu_fair = matrix(rep(mu_fair_2theta, each = D), D, nrow(d), dimnames = list(NULL, paste0("mu_fair[", seq_len(nrow(d)), "]")))
  )
  mock_custom <- list(
    metadata = function() list(stan_variables = names(draws_custom)),
    draws = function(var, format = "matrix") draws_custom[[var]]
  )
  obj <- structure(
    list(fit = mock_custom, spec = spec, data_info = di, design = dc, stan_code = custom_code),
    class = "bmfrm_fit"
  )

  expect_equal(get_eval_contract(obj)$source, "custom_unverified")
  expect_false(can_evaluate_in_r(obj))
  expect_null(stan_mu_fair_adjusted_facets(obj))

  # Default fair_scores() uses Stan's own mu_fair (~3.666 for P1) rather than recomputing 1 * theta
  expect_message(fs_stan <- fair_scores(obj, summary = FALSE), "Using custom Stan built-in mu_fair")
  expect_equal(fs_stan$mu_fair, mu_fair_2theta, tolerance = 1e-12)
  expect_true(all(is.na(fs_stan$likelihood_identified)))

  # Standardized fair_scores(), explicit adjust_facets (even empty), and missing-GQ residuals/pp_check/facet_fit must stop
  expect_error(fair_scores(obj, standardize = TRUE), "Unverified custom Stan model")
  expect_error(fair_scores(obj, adjust_facets = "rater"), "Unverified custom Stan model")
  expect_error(fair_scores(obj, adjust_facets = character(0)), "Unverified custom Stan model")
  expect_error(residuals(obj), "Unverified custom Stan model")
  expect_error(pp_check(obj), "Unverified custom Stan model")
  expect_error(facet_fit(obj, "rater"), "Unverified custom Stan model")

  # Modifying a bundled Stan file's linear predictor must also downgrade it to custom_unverified
  m3 <- bmfrm_model("three-facet/person-rater-scenario-main")
  tampered_code <- sub("real eta = theta\\[j\\]", "real eta = 2 * theta[j]", m3$code)
  d3 <- data.frame(person = c("P1", "P1"), rater = c("R1", "R2"), scenario = c("S1", "S2"), score = c(2L, 3L))
  spec3 <- parse_bmfrm_formula(score ~ person + rater + scenario, d3, K = 4)
  spec3$custom_stan <- TRUE
  expect_equal(resolve_eval_contract(spec3, m3$code)$source, "bundled_verified")
  expect_true(resolve_eval_contract(spec3, m3$code)$r_recompute)
  expect_equal(resolve_eval_contract(spec3, tampered_code)$source, "custom_unverified")
  expect_false(resolve_eval_contract(spec3, tampered_code)$r_recompute)
})

test_that("summarise_bias_ir and summarise_bias apply the exact same directional probability rule", {
  # 1000 draws where cell [1,1] has P(B > 0) = 0.960 (95% two-sided CrI includes 0, P > 0 >= 0.95)
  D <- 1000L
  c11 <- c(rep(-0.05, 40L), rep(0.25, 960L))
  mat <- cbind(c11, -c11, -c11, c11)
  colnames(mat) <- c("bias[1,1]", "bias[2,1]", "bias[1,2]", "bias[2,2]")
  arr <- array(mat, dim = c(D, 1L, 4L), dimnames = list(NULL, NULL, colnames(mat)))

  mock_legacy <- list(draws = function(var, format = "draws_array") arr)
  ir95  <- summarise_bias_ir(mock_legacy, item_labels = c("I1", "I2"), rater_labels = c("R1", "R2"), prob_thresh = 0.95)
  ir975 <- summarise_bias_ir(mock_legacy, item_labels = c("I1", "I2"), rater_labels = c("R1", "R2"), prob_thresh = 0.975)

  d <- expand.grid(person = paste0("P", 1:4), item = c("I1", "I2"), rater = c("R1", "R2"))
  d$score <- 1L
  spec <- parse_bmfrm_formula(score ~ person + item + rater + item:rater, d, K = 4)
  ch <- check_design_bmfrm(spec$formula, d, K = 4)
  mat_new <- mat
  colnames(mat_new) <- sub("^bias", "itemrater", colnames(mat))
  obj_new <- structure(list(
    spec = spec, design = ch,
    data_info = list(facet_labels = list(person = paste0("P", 1:4), item = c("I1", "I2"), rater = c("R1", "R2")),
                     data_clean = d),
    fit = list(draws = function(...) mat_new)
  ), class = "bmfrm_fit")
  new95  <- summarise_bias(obj_new, "item:rater", prob_thresh = 0.95, delta = 0, flag_method = "prob")
  new975 <- summarise_bias(obj_new, "item:rater", prob_thresh = 0.975, delta = 0, flag_method = "prob")

  # Align by (item, rater) and confirm identical flags at both 0.95 and 0.975
  ord_ir  <- order(ir95$i, ir95$r)
  ord_new <- order(new95$i, new95$r)
  expect_equal(ir95$flag_bias[ord_ir], new95$flag_bias[ord_new])
  expect_equal(ir95$flag_bias[ord_ir][1], "positive_bias")
  expect_equal(ir975$flag_bias[ord_ir], new975$flag_bias[ord_new])
  expect_true(all(ir975$flag_bias == "no_clear_bias"))
})

test_that("facet_infit_outfit and add_fit_cri_bayes support bare names, strings, and string variables", {
  N <- 8L; D <- 20L
  set.seed(501)
  mu_m <- matrix(2.5, D, N)
  sg_m <- matrix(1.0, D, N)
  df <- tibble::tibble(
    n = seq_len(N),
    x = c(1L, 2L, 3L, 4L, 2L, 3L, 2L, 3L),
    rater = rep(1:2, each = 4L),
    mu_hat = colMeans(mu_m),
    sigma2_hat = colMeans(sg_m),
    resid = x - mu_hat,
    weight = sigma2_hat,
    z = resid / sqrt(sigma2_hat),
    model = "m1",
    mu_draws = lapply(seq_len(N), function(j) as.numeric(mu_m[, j])),
    sigma2_draws = lapply(seq_len(N), function(j) as.numeric(sg_m[, j]))
  )

  # Bare unquoted column name must NOT fail with "object 'rater' not found"
  f_bare <- facet_infit_outfit(df, rater)
  f_str  <- facet_infit_outfit(df, "rater")
  col_var <- "rater"
  f_var  <- facet_infit_outfit(df, col_var)
  expect_equal(f_bare, f_str)
  expect_equal(f_var, f_str)

  c_bare <- add_fit_cri_bayes(f_bare, df, rater)
  c_str  <- add_fit_cri_bayes(f_str, df, "rater")
  c_var  <- add_fit_cri_bayes(f_var, df, col_var)
  expect_equal(c_bare, c_str)
  expect_equal(c_var, c_str)
})

test_that("reserved Stan facet names are rejected consistently across non-hierarchical and hierarchical paths", {
  d <- expand.grid(person = paste0("P", 1:4), theta = c("A", "B"), tau = c("T1", "T2"),
                   item_raw = c("X1", "X2"), model = c("M1", "M2"), stringsAsFactors = FALSE)
  d$score <- rep(1:4, length.out = nrow(d))

  # Default non-hierarchical RSM must reject reserved names just like hierarchical models
  expect_error(parse_bmfrm_formula(score ~ person + theta, d, K = 4), "nonreserved Stan identifiers")
  expect_error(parse_bmfrm_formula(score ~ person + tau, d, K = 4), "nonreserved Stan identifiers")
  expect_error(parse_bmfrm_formula(score ~ person + item_raw, d, K = 4), "nonreserved Stan identifiers")
  expect_error(parse_bmfrm_formula(score ~ person + model, d, K = 4), "nonreserved Stan identifiers")

  # Even bypassing parse_bmfrm_formula(), resolve_bmfrm_spec() must reject before early return
  bad_spec <- list(
    response = "score", facets_main = c("person", "theta"), facets_bias = character(),
    K = 4L, family = "rsm", hierarchical = character(), save_gq = "full"
  )
  expect_error(resolve_bmfrm_spec(bad_spec), "nonreserved Stan identifiers")
})

test_that("check_design_bmfrm guards oversized interactions before allocating Kronecker bases and matches dense SVD in factored Schur mode", {
  # 200 persons x 10 items with person:item: p_other = 9 + 199*9 = 1800, p_other^2 = 3.24e6 > 2e6
  # Must fail fast BEFORE allocating 28.7 MB Kronecker arrays (and also for 10,000 x 10 -> 72 GB!)
  d_200 <- expand.grid(person = seq_len(200L), item = seq_len(10L))
  expect_error(check_design_bmfrm(score ~ person + item + person:item, d_200, K = 4), "exceeds max_cells limit")

  d_10k <- data.frame(
    person = factor(rep(seq_len(10000L), each = 2L)),
    item   = factor(rep(seq_len(10L), length.out = 20000L))
  )
  expect_error(check_design_bmfrm(score ~ person + item + person:item, d_10k, K = 4), "exceeds max_cells limit")

  # Factorized Schur mode on person:item (both full-rank and rank-deficient missing cell) must match dense SVD exactly
  d_pi <- expand.grid(person = paste0("P", 1:10), item = paste0("I", 1:3), rater = c("R1", "R2"), stringsAsFactors = FALSE)
  ch_d1 <- check_design_bmfrm(score ~ person + item + rater + person:item, d_pi, K = 4, method = "dense")
  ch_s1 <- check_design_bmfrm(score ~ person + item + rater + person:item, d_pi, K = 4, method = "schur")
  expect_equal(ch_s1$rank, ch_d1$rank)
  expect_equal(ch_s1$nullity, ch_d1$nullity)
  expect_equal(ch_s1$condition_number, ch_d1$condition_number, tolerance = 1e-10)

  d_pi_miss <- subset(d_pi, !(person == "P10" & item == "I3"))
  ch_d2 <- check_design_bmfrm(score ~ person + item + rater + person:item, d_pi_miss, K = 4, method = "dense")
  ch_s2 <- check_design_bmfrm(score ~ person + item + rater + person:item, d_pi_miss, K = 4, method = "schur")
  expect_equal(ch_s2$rank, ch_d2$rank)
  expect_equal(ch_s2$nullity, ch_d2$nullity)
  # Verify factored apply_kron_null produces the exact same effect_null_space subspace (up to sign/rotation)
  Q_d <- qr.Q(qr(ch_d2$effect_null_space))
  Q_s <- qr.Q(qr(ch_s2$effect_null_space))
  expect_lt(max(abs(Q_d %*% t(Q_d) - Q_s %*% t(Q_s))), 1e-11)

  # 1-level facet in person:item has 0 interaction DOF: must NOT allocate a P x (P-1) person basis
  d_1item <- data.frame(
    person = factor(seq_len(1600L)),
    item   = factor(rep("I1", 1600L))
  )
  ch_1item <- check_design_bmfrm(score ~ person + item + person:item, d_1item, K = 4)
  expect_true(ch_1item$identifiable)
  expect_equal(ch_1item$rank, 1600L)
  expect_equal(ch_1item$nullity, 0L)
  expect_equal(ch_1item$solver, "schur_complement")

  # Dense path with 1-level facet in person:item must also work cleanly with 0 interaction DOF
  d_1item_small <- data.frame(person = factor(1:5), item = factor(rep("I1", 5L)))
  ch_1item_d <- check_design_bmfrm(score ~ person + item + person:item, d_1item_small, K = 4, method = "dense")
  expect_true(ch_1item_d$identifiable)
  expect_equal(ch_1item_d$rank, 5L)
  expect_equal(ch_1item_d$nullity, 0L)
})
