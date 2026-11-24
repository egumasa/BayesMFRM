# Tests for target-aware Bayesian model comparison:
# loo.bmfrm_fit(), cv_bmfrm() (rating & grouped marginalized CV), and compare_bmfrm()

make_mock_cmp_fit <- function(d = NULL,
                              K = 4L,
                              D_iter = 40L,
                              n_chains = 4L,
                              family = "rsm",
                              threshold_facet = NULL,
                              hierarchical = c("person", "rater"),
                              include_bias = FALSE,
                              bias_shrinkage = "none",
                              custom_stan = FALSE,
                              save_stan_ll = TRUE,
                              include_chains = TRUE,
                              converged = TRUE,
                              obs_id = NULL,
                              seed = 20261124L) {
  set.seed(seed)
  if (is.null(d)) {
    d <- expand.grid(
      person = paste0("P", 1:6),
      item = c("I1", "I2", "I3"),
      rater = c("R1", "R2", "R3", "R4"),
      stringsAsFactors = FALSE
    )
    # Realistic scores loosely aligned with person/item/rater ordering (prevents extreme Pareto-k on mock data)
    p_eff <- seq(-1, 1, length.out = 6L)[match(d$person, paste0("P", 1:6))]
    i_eff <- seq(-0.5, 0.5, length.out = 3L)[match(d$item, c("I1", "I2", "I3"))]
    r_eff_sim <- seq(-0.7, 0.7, length.out = 4L)[match(d$rater, c("R1", "R2", "R3", "R4"))]
    eta_sim <- p_eff - i_eff - r_eff_sim
    d$score <- pmax(1L, pmin(K, as.integer(round(2.5 + 0.6 * eta_sim + stats::rnorm(nrow(d), 0, 0.35)))))
  }
  D <- as.integer(D_iter * n_chains)

  # Helper for mildly autocorrelated mock chains (avoids negative-lag white-noise ESS capping warning)
  ar_draws <- function(ncol_m, means, sds) {
    out <- matrix(0, nrow = D, ncol = ncol_m)
    for (j in seq_len(ncol_m)) {
      mu_j <- if (length(means) == 1L) means else means[j]
      sd_j <- if (length(sds) == 1L) sds else sds[j]
      for (c_i in seq_len(n_chains)) {
        v <- numeric(D_iter)
        v[1L] <- stats::rnorm(1L, 0, sd_j)
        for (t in 2:D_iter) v[t] <- 0.35 * v[t - 1L] + sqrt(1 - 0.35^2) * stats::rnorm(1L, 0, sd_j)
        rows_c <- ((c_i - 1L) * D_iter + 1L):(c_i * D_iter)
        out[rows_c, j] <- mu_j + v
      }
    }
    out
  }

  fml <- if (include_bias) {
    score ~ person + item + rater + item:rater
  } else {
    score ~ person + item + rater
  }
  spec <- parse_bmfrm_formula(
    fml, d, K = K, family = family,
    threshold_facet = threshold_facet,
    hierarchical = hierarchical
  )
  priors <- NULL
  if (include_bias && "item:rater" %in% hierarchical && bias_shrinkage != "none") {
    sc <- if (bias_shrinkage == "student_t") "marginal" else "raw"
    dp <- get_default_priors(spec$facets_main, spec$facets_bias, spec$family, spec$hierarchical)
    dp <- Filter(function(p) !(identical(p$class, "sd") && identical(p$facet, "item:rater")), dp)
    priors <- c(dp, list(
      prior("normal(0, 0.5)", class = "sd", facet = "item:rater",
            scale = sc, local_shrinkage = bias_shrinkage)
    ))
  }
  spec <- resolve_bmfrm_spec(spec, priors)
  if (custom_stan) spec$custom_stan <- TRUE
  spec$eval_contract <- resolve_eval_contract(
    spec,
    stan_code = if (custom_stan) "parameters { vector[6] theta; }" else NULL
  )

  data_info <- prepare_data_bmfrm(spec, d, obs_id = obs_id)
  dc <- check_design_bmfrm(fml, data_info$data_clean, K = K, family = family, threshold_facet = threshold_facet)

  J_p <- length(data_info$facet_labels$person)
  J_i <- length(data_info$facet_labels$item)
  J_r <- length(data_info$facet_labels$rater)

  theta <- ar_draws(J_p, seq(-1, 1, length.out = J_p), 0.12)
  colnames(theta) <- paste0("theta[", seq_len(J_p), "]")

  item_raw <- ar_draws(J_i, seq(-0.5, 0.5, length.out = J_i), 0.08)
  item_draws <- if (J_i > 1L) sweep(item_raw, 1L, rowMeans(item_raw), "-") else matrix(0, D, 1L)
  colnames(item_draws) <- paste0("item[", seq_len(J_i), "]")

  rater_raw <- ar_draws(J_r, seq(-0.7, 0.7, length.out = J_r), 0.10)
  rater_draws <- if (J_r > 1L) sweep(rater_raw, 1L, rowMeans(rater_raw), "-") else matrix(0, D, 1L)
  colnames(rater_draws) <- paste0("rater[", seq_len(J_r), "]")

  draw_store <- list(theta = theta, item = item_draws, rater = rater_draws)

  if ("person" %in% hierarchical) {
    draw_store$mu_theta <- ar_draws(1L, 0.1, 0.06)
    colnames(draw_store$mu_theta) <- "mu_theta"
    draw_store$sigma_theta <- matrix(pmax(0.2, ar_draws(1L, 0.8, 0.06)), D, 1L, dimnames = list(NULL, "sigma_theta"))
  }
  if ("rater" %in% hierarchical) {
    draw_store$sigma_rater <- matrix(pmax(0.2, ar_draws(1L, 0.65, 0.06)), D, 1L, dimnames = list(NULL, "sigma_rater"))
  }
  if ("item" %in% hierarchical) {
    draw_store$sigma_item <- matrix(pmax(0.2, ar_draws(1L, 0.5, 0.05)), D, 1L, dimnames = list(NULL, "sigma_item"))
  }

  if (identical(family, "rsm")) {
    tau_raw <- ar_draws(K - 1L, seq(-1, 1, length.out = K - 1L), 0.05)
    tau_draws <- sweep(tau_raw, 1L, rowMeans(tau_raw), "-")
    colnames(tau_draws) <- paste0("tau[", seq_len(K - 1L), "]")
    draw_store$tau <- tau_draws
  } else {
    J_tf <- length(data_info$facet_labels[[threshold_facet]])
    tau_common <- sweep(
      matrix(stats::rnorm(D * (K - 1L), rep(seq(-1, 1, length.out = K - 1L), each = D), 0.06), D, K - 1L),
      1L, rep(0, D), "-"
    )
    tau_common <- sweep(tau_common, 1L, rowMeans(tau_common), "-")
    colnames(tau_common) <- paste0("tau_common[", seq_len(K - 1L), "]")
    tau_mat <- matrix(0, nrow = D, ncol = J_tf * (K - 1L))
    cn <- character(J_tf * (K - 1L))
    for (k in seq_len(K - 1L)) {
      for (j in seq_len(J_tf)) {
        col_i <- j + (k - 1L) * J_tf
        cn[col_i] <- paste0("tau[", j, ",", k, "]")
        tau_mat[, col_i] <- tau_common[, k] + (j - mean(seq_len(J_tf))) * 0.1 * (k - (K / 2))
      }
    }
    colnames(tau_mat) <- cn
    draw_store$tau <- tau_mat
    if ("tau" %in% hierarchical) {
      draw_store$tau_common <- tau_common
      draw_store$sigma_tau <- matrix(stats::runif(D, 0.2, 0.3), D, 1L, dimnames = list(NULL, "sigma_tau"))
    }
  }

  if (include_bias) {
    b_dc <- matrix(0, nrow = D, ncol = J_i * J_r)
    if (J_i > 1L && J_r > 1L) {
      b_raw <- ar_draws(J_i * J_r, 0, 0.12)
      for (dd in seq_len(D)) {
        m <- matrix(b_raw[dd, ], nrow = J_i, ncol = J_r)
        m <- sweep(m, 1L, rowMeans(m), "-")
        m <- sweep(m, 2L, colMeans(m), "-")
        b_dc[dd, ] <- as.vector(m)
      }
    }
    colnames(b_dc) <- as.vector(outer(seq_len(J_i), seq_len(J_r), function(i, j) paste0("itemrater[", i, ",", j, "]")))
    draw_store$itemrater <- b_dc
    if ("item:rater" %in% hierarchical) {
      draw_store$sigma_itemrater <- matrix(pmax(0.1, ar_draws(1L, 0.35, 0.04)), D, 1L, dimnames = list(NULL, "sigma_itemrater"))
      if (bias_shrinkage == "horseshoe") {
        draw_store$c2_itemrater <- matrix(pmax(0.5, ar_draws(1L, 2.0, 0.2)), D, 1L, dimnames = list(NULL, "c2_itemrater"))
      }
    }
  }

  make_obj <- function(store, ll_arr = NULL) {
    conv <- list(
      converged = converged,
      status = if (converged) "ok" else "warning",
      n_chains = if (include_chains) n_chains else 1L,
      max_rhat = if (converged) 1.003 else 1.08,
      min_ess_bulk = if (converged) 850 else 120,
      n_divergent = if (converged) 0L else 5L
    )
    fit_mock <- list(
      num_chains = if (include_chains) function() n_chains else NULL,
      metadata = function() {
        md <- list(stan_variables = c(names(store), if (!is.null(ll_arr)) "log_lik"))
        if (include_chains) md$iter_sampling <- D_iter
        md
      },
      draws = function(var, format = "matrix") {
        if (identical(var, "log_lik") && !is.null(ll_arr)) {
          if (identical(format, "draws_array")) return(ll_arr)
          mat <- matrix(ll_arr, nrow = D, ncol = dim(ll_arr)[3L])
          colnames(mat) <- dimnames(ll_arr)[[3L]]
          return(mat)
        }
        if (!var %in% names(store)) stop("Missing variable: ", var)
        m <- store[[var]]
        if (identical(format, "draws_array")) {
          if (!include_chains) stop("No chain structure")
          arr <- array(m, dim = c(D_iter, n_chains, ncol(m)),
                       dimnames = list(NULL, NULL, colnames(m)))
          return(arr)
        }
        m
      }
    )
    structure(
      list(
        fit = fit_mock, spec = spec, data_info = data_info,
        design = dc, priors = spec$priors, stan_code = NULL,
        convergence = conv
      ),
      class = "bmfrm_fit"
    )
  }

  obj <- make_obj(draw_store)
  if (save_stan_ll && !custom_stan) {
    ll_mat <- compute_bmfrm_quantities(obj, compute_log_lik = TRUE)$log_lik
    N_obs <- nrow(data_info$data_clean)
    # Deliberately reverse column order in Stan array to test strict 1..N ordering!
    rev_ord <- rev(seq_len(N_obs))
    ll_arr <- array(ll_mat[, rev_ord, drop = FALSE], dim = c(D_iter, n_chains, N_obs),
                    dimnames = list(NULL, NULL, paste0("log_lik[", rev_ord, "]")))
    obj <- make_obj(draw_store, ll_arr = ll_arr)
  }
  obj
}

test_that("loo.bmfrm_fit matches Stan log_lik vs low-memory R streaming, chunk invariance, and enforces chain & contract guards", {
  obj <- make_mock_cmp_fit(D_iter = 40L, n_chains = 4L, include_bias = TRUE,
                           hierarchical = c("person", "rater", "item:rater"), save_stan_ll = TRUE)

  l_stan <- loo::loo(obj, use_stan_log_lik = TRUE)
  l_r_c1 <- loo::loo(obj, use_stan_log_lik = FALSE, chunk_size = 1L)
  l_r_c7 <- loo::loo(obj, use_stan_log_lik = FALSE, chunk_size = 7L)

  expect_s3_class(l_stan, "bmfrm_loo")
  expect_s3_class(l_r_c1, "bmfrm_loo")
  expect_equal(l_stan$bmfrm_meta$source, "stan_log_lik")
  expect_equal(l_r_c1$bmfrm_meta$source, "r_streaming")
  expect_equal(l_r_c1$bmfrm_meta$actual_target, "new_rating_known_levels")

  # Exact agreement between Stan full GQ and R minimal streaming across chunk sizes
  expect_equal(l_r_c1$estimates, l_stan$estimates, tolerance = 1e-12)
  expect_equal(l_r_c7$estimates, l_stan$estimates, tolerance = 1e-12)
  expect_equal(l_r_c1$pointwise, l_stan$pointwise, tolerance = 1e-12)
  expect_equal(l_r_c1$diagnostics$pareto_k, l_stan$diagnostics$pareto_k, tolerance = 1e-12)
  expect_equal(l_r_c1$diagnostics$r_eff, l_stan$diagnostics$r_eff, tolerance = 1e-12)
  expect_equal(l_r_c1$bmfrm_meta$mcse_elpd_loo, l_stan$bmfrm_meta$mcse_elpd_loo, tolerance = 1e-12)

  # Verify r_eff was computed on exp(log_lik), NOT raw log_lik
  ll_arr_raw <- obj$fit$draws("log_lik", format = "draws_array")
  idx <- as.integer(sub(".*\\[([0-9]+)\\]", "\\1", dimnames(ll_arr_raw)[[3L]]))
  ll_arr_ord <- ll_arr_raw[, , order(idx), drop = FALSE]
  ref_reff <- as.numeric(loo::relative_eff(exp(ll_arr_ord)))
  expect_equal(as.numeric(l_r_c1$diagnostics$r_eff), ref_reff, tolerance = 1e-12)

  # Missing chain layout must NOT silently assume r_eff = 1
  obj_noch <- make_mock_cmp_fit(save_stan_ll = FALSE, include_chains = FALSE)
  expect_error(loo::loo(obj_noch), "refusing to silently assume r_eff = 1")
  # Explicit r_eff = NA succeeds deliberately
  l_na <- loo::loo(obj_noch, r_eff = NA)
  expect_s3_class(l_na, "bmfrm_loo")

  # Custom unverified Stan model is blocked from R-side log_lik recomputation
  obj_cust <- make_mock_cmp_fit(custom_stan = TRUE, save_stan_ll = FALSE)
  expect_error(loo::loo(obj_cust), "Unverified custom Stan model")
})

test_that("loo.bmfrm_fit flags singleton observations and cv_bmfrm enforces strict target, hierarchy, degeneracy, and droplevels", {
  # 1. Singleton examinee in PSIS-LOO triggers warning and explicit actual_target
  d_sing <- rbind(
    expand.grid(person = paste0("P", 1:4), item = c("I1", "I2"), rater = c("R1", "R2"), stringsAsFactors = FALSE),
    data.frame(person = "P_SINGLE", item = "I1", rater = "R1", stringsAsFactors = FALSE)
  )
  d_sing$score <- rep(1:4, length.out = nrow(d_sing))
  obj_sing <- make_mock_cmp_fit(d = d_sing, save_stan_ll = FALSE)
  expect_warning(l_sing <- loo::loo(obj_sing), "sole rating")
  expect_equal(l_sing$bmfrm_meta$actual_target, "new_rating_with_singleton_levels")
  expect_equal(l_sing$bmfrm_meta$singleton_obs, nrow(d_sing))

  # 2. Collateral unseen examinee when holding out a rater MUST stop cv_bmfrm()
  d_collat <- rbind(
    expand.grid(person = paste0("P", 1:4), item = c("I1", "I2"), rater = c("R1", "R2", "R3"), stringsAsFactors = FALSE),
    data.frame(person = "P_ONLY_R4", item = c("I1", "I2"), rater = "R4", stringsAsFactors = FALSE)
  )
  d_collat$score <- rep(1:4, length.out = nrow(d_collat))
  obj_collat <- make_mock_cmp_fit(d = d_collat, hierarchical = c("person", "rater"), save_stan_ll = FALSE)
  expect_error(
    cv_bmfrm(obj_collat, target = "rater", K_folds = 4L, seed = 1L),
    "also makes other facet\\(s\\) unseen in training: person \\(P_ONLY_R4\\)"
  )

  # 3. Rating CV fold dropping all ratings of a level stops under strict_target = TRUE
  obj_std <- make_mock_cmp_fit(save_stan_ll = FALSE)
  bad_folds <- ifelse(obj_std$data_info$data_clean$person == "P1", 1L, 2L)
  expect_error(
    cv_bmfrm(obj_std, target = "rating", folds = bad_folds),
    "drops all training observations for facet\\(s\\): person"
  )

  # 4. Hierarchy guards for main effect, interaction, and PCM tau
  obj_no_h_r <- make_mock_cmp_fit(hierarchical = "person", save_stan_ll = FALSE)
  expect_error(cv_bmfrm(obj_no_h_r, target = "rater", K_folds = 2L), "Missing from hierarchical: rater")

  obj_no_h_b <- make_mock_cmp_fit(include_bias = TRUE, hierarchical = c("person", "rater"), save_stan_ll = FALSE)
  expect_error(cv_bmfrm(obj_no_h_b, target = "rater", K_folds = 2L), "Missing from hierarchical: item:rater")

  obj_no_h_tau <- make_mock_cmp_fit(family = "pcm", threshold_facet = "rater",
                                    hierarchical = c("person", "rater"), save_stan_ll = FALSE)
  expect_error(cv_bmfrm(obj_no_h_tau, target = "rater", K_folds = 2L), "Missing from hierarchical: tau")

  # 5. 1-level degeneracy stop when 2 raters with item:rater are split into 2 folds
  d_2r <- expand.grid(person = paste0("P", 1:4), item = c("I1", "I2"), rater = c("R1", "R2"), stringsAsFactors = FALSE)
  d_2r$score <- rep(1:4, length.out = nrow(d_2r))
  obj_2r <- make_mock_cmp_fit(d = d_2r, include_bias = TRUE,
                              hierarchical = c("person", "rater", "item:rater"), save_stan_ll = FALSE)
  expect_error(
    cv_bmfrm(obj_2r, target = "rater", K_folds = 2L, seed = 1L),
    "<= 1 observed level, collapsing interaction 'item:rater'"
  )

  # 6. Training fold preprocessing drops unused factor levels AND locks K!
  d_fac <- expand.grid(person = paste0("P", 1:4), item = c("I1", "I2", "I3"), rater = c("R1", "R2", "R3", "R4"))
  d_fac$score <- rep(1:3, length.out = nrow(d_fac)) # max score 3 on K = 4 rubric!
  obj_fac <- make_mock_cmp_fit(d = d_fac, K = 4L, hierarchical = c("person", "rater"), save_stan_ll = FALSE)
  seen_levels <- list()
  seen_K <- integer(0)
  mock_refit <- function(train_df, fold, spec, fit_args) {
    seen_levels[[fold]] <<- levels(train_df$rater)
    seen_K <<- c(seen_K, spec$K)
    make_mock_cmp_fit(d = train_df, K = spec$K, D_iter = 20L, n_chains = 4L,
                      hierarchical = spec$hierarchical, save_stan_ll = FALSE, seed = 100L + fold)
  }
  cv_r <- suppressWarnings(
    cv_bmfrm(obj_fac, target = "rater", K_folds = 2L, M = 10L, seed = 5L, refit_fn = mock_refit)
  )
  expect_true(all(lengths(seen_levels) == 2L)) # held-out raters dropped cleanly from factor levels!
  expect_identical(seen_K, c(4L, 4L))          # K locked to 4, never re-inferred as 3!
})

test_that("cv_bmfrm computes exact double-logmeanexp group/pointwise scores, M=1, chunk invariance, and guarded compare_bmfrm()", {
  obj_main <- make_mock_cmp_fit(D_iter = 40L, n_chains = 4L, include_bias = FALSE,
                                hierarchical = c("person", "rater"), save_stan_ll = FALSE, seed = 11L)
  shared_d <- obj_main$data_info$data_clean
  obj_bias <- make_mock_cmp_fit(d = shared_d, D_iter = 40L, n_chains = 4L, include_bias = TRUE,
                                hierarchical = c("person", "rater", "item:rater"),
                                bias_shrinkage = "student_t", save_stan_ll = FALSE, seed = 22L)
  obj_diff_y <- make_mock_cmp_fit(D_iter = 40L, n_chains = 4L, include_bias = FALSE,
                                  hierarchical = c("person", "rater"), save_stan_ll = FALSE, seed = 999L)

  refit_for <- function(include_b, sh = "none", converged = TRUE, fail_fold = NULL) {
    function(train_df, fold, spec, fit_args) {
      if (!is.null(fail_fold) && fold == fail_fold) stop("simulated Stan failure")
      make_mock_cmp_fit(
        d = train_df, K = spec$K, D_iter = 40L, n_chains = 4L,
        include_bias = include_b, bias_shrinkage = sh,
        hierarchical = spec$hierarchical, save_stan_ll = FALSE,
        converged = converged, seed = 500L + fold
      )
    }
  }

  # 1. Rating-level CV
  cv_rat_1 <- cv_bmfrm(obj_main, target = "rating", K_folds = 3L, seed = 42L,
                       refit_fn = refit_for(FALSE))
  cv_rat_2 <- cv_bmfrm(obj_bias, target = "rating", folds = cv_rat_1$folds, seed = 42L,
                       refit_fn = refit_for(TRUE, "student_t"))
  expect_equal(cv_rat_1$actual_target, "new_rating_known_levels")
  expect_true(cv_rat_1$all_folds_converged)
  expect_true(cv_rat_1$all_folds_identifiable)
  expect_true(is.finite(cv_rat_1$elpd_pointwise))
  expect_true(is.na(cv_rat_1$elpd_group))

  cmp_rat <- compare_bmfrm(main = cv_rat_1, bias = cv_rat_2)
  expect_s3_class(cmp_rat, "bmfrm_model_comparison")
  expect_equal(attr(cmp_rat, "comparison_meta")$uncertainty_unit, "rating")

  # 2. Grouped rater CV: exact double-logmeanexp, chunk_size invariance, M = 1 valid
  cv_grp_1 <- cv_bmfrm(obj_main, target = "rater", K_folds = 2L, M = 25L, seed = 77L,
                       chunk_size = 50L, refit_fn = refit_for(FALSE))
  cv_grp_1_c1 <- cv_bmfrm(obj_main, target = "rater", folds = cv_grp_1$folds, M = 25L, seed = 77L,
                          chunk_size = 1L, refit_fn = refit_for(FALSE))
  cv_grp_2 <- cv_bmfrm(obj_bias, target = "rater", folds = cv_grp_1$folds, M = 25L, seed = 77L,
                       refit_fn = refit_for(TRUE, "student_t"))
  cv_grp_m1 <- cv_bmfrm(obj_main, target = "rater", folds = cv_grp_1$folds, M = 1L, seed = 77L,
                        refit_fn = refit_for(FALSE))

  # Bit-for-bit chunk invariance
  expect_equal(cv_grp_1$elpd_group, cv_grp_1_c1$elpd_group, tolerance = 1e-12)
  expect_equal(cv_grp_1$elpd_pointwise, cv_grp_1_c1$elpd_pointwise, tolerance = 1e-12)
  expect_equal(cv_grp_1$group_scores$elpd, cv_grp_1_c1$group_scores$elpd, tolerance = 1e-12)

  # Group joint density != sum of independent rating marginal densities (shared u_{gdm}!)
  expect_false(isTRUE(all.equal(cv_grp_1$elpd_group, cv_grp_1$elpd_pointwise)))
  expect_equal(cv_grp_m1$integration_diagnostics$M, 1L)
  expect_true(all(cv_grp_m1$group_scores$mc_ess_m_median == 1))
  expect_true(all(cv_grp_1$group_scores$mc_ess_m_median >= 1 & cv_grp_1$group_scores$mc_ess_m_median <= 25))

  # 3. Guarded compare_bmfrm(): group SE & group-clustered pointwise SE
  cmp_grp <- compare_bmfrm(main = cv_grp_1, bias = cv_grp_2, score = "auto")
  expect_equal(attr(cmp_grp, "comparison_meta")$score_unit, "group")
  expect_equal(attr(cmp_grp, "comparison_meta")$uncertainty_unit, "group")

  cmp_grp_pt <- compare_bmfrm(main = cv_grp_1, bias = cv_grp_2, score = "pointwise")
  expect_equal(attr(cmp_grp_pt, "comparison_meta")$uncertainty_unit, "group_clustered_pointwise")
  # Verify SE matches group-clustered paired difference exactly (NOT naive N-rating SE!)
  canon_g <- cv_grp_1$group_scores$group
  d_g <- as.numeric(tapply(cv_grp_1$pointwise$elpd - cv_grp_2$pointwise$elpd,
                           cv_grp_1$pointwise$group, sum)[canon_g])
  expected_se <- sqrt(length(canon_g) * stats::var(d_g))
  expect_equal(as.numeric(cmp_grp_pt[2L, "se_diff"]), expected_se, tolerance = 1e-12)

  # 4. Refusal on incompatible targets, methods, splits, response data, and unconverged/failed folds
  l_main <- loo::loo(obj_main)
  l_diff_y <- loo::loo(obj_diff_y)
  expect_error(compare_bmfrm(l_main, l_diff_y), "different response values")
  expect_error(compare_bmfrm(l_main, cv_rat_1), "different evaluation methods")
  expect_error(compare_bmfrm(cv_rat_1, cv_grp_1), "different prediction targets")

  cv_diff_folds <- cv_bmfrm(obj_main, target = "rating", K_folds = 2L, seed = 99L,
                            refit_fn = refit_for(FALSE))
  expect_error(compare_bmfrm(cv_rat_1, cv_diff_folds), "different fold splits")

  cv_unconv <- suppressWarnings(
    cv_bmfrm(obj_main, target = "rating", folds = cv_rat_1$folds, seed = 42L,
             refit_fn = refit_for(FALSE, converged = FALSE))
  )
  expect_false(cv_unconv$all_folds_converged)
  expect_error(compare_bmfrm(cv_rat_1, cv_unconv), "did not converge")
  expect_s3_class(compare_bmfrm(cv_rat_1, cv_unconv, allow_unconverged = TRUE), "bmfrm_model_comparison")

  # Failed fold never silently dropped -> aggregate ELPD is NA and comparison blocked
  cv_fail <- suppressWarnings(
    cv_bmfrm(obj_main, target = "rating", folds = cv_rat_1$folds, seed = 42L,
             refit_fn = refit_for(FALSE, fail_fold = 2L))
  )
  expect_true(cv_fail$any_refit_failed)
  expect_true(is.na(cv_fail$elpd_pointwise))
  expect_error(compare_bmfrm(cv_rat_1, cv_fail), "failed refitting")
})

test_that("compare_bmfrm verifies pre-preprocessing obs_id and catches same-score row swaps", {
  obj_init <- make_mock_cmp_fit(include_bias = FALSE, seed = 20261124L)
  d_base <- obj_init$data_info$data_clean
  obj_base <- make_mock_cmp_fit(d = d_base, include_bias = FALSE, seed = 20261124L)

  # Find two observations with the EXACT same score but different person/item/rater
  same_score_pair <- NULL
  for (s_val in unique(d_base$score)) {
    idx_s <- which(d_base$score == s_val)
    if (length(idx_s) >= 2L) {
      for (a in idx_s) {
        for (b in idx_s) {
          if (a < b && (d_base$person[a] != d_base$person[b] || d_base$item[a] != d_base$item[b])) {
            same_score_pair <- c(a, b)
            break
          }
        }
        if (!is.null(same_score_pair)) break
      }
    }
    if (!is.null(same_score_pair)) break
  }
  expect_true(!is.null(same_score_pair))

  # 1. Exact reproduction of the silent same-score row-swap bug (reset rownames, no explicit obs_id):
  #    Both models have identical y vector and 1:N rownames, but swapped rows -> MUST ERROR!
  d_swap_reset <- d_base
  d_swap_reset[same_score_pair, ] <- d_swap_reset[rev(same_score_pair), ]
  rownames(d_swap_reset) <- NULL
  obj_swap_reset <- make_mock_cmp_fit(d = d_swap_reset, include_bias = FALSE, seed = 20261124L)
  loo_base <- loo::loo(obj_base, use_stan_log_lik = FALSE)
  loo_swap_reset <- loo::loo(obj_swap_reset, use_stan_log_lik = FALSE)
  expect_error(
    compare_bmfrm(a = loo_base, b = loo_swap_reset),
    "observation alignment mismatch on shared facet"
  )

  # 2. Row permutation WITH preserved rownames -> aligns by obs_id and yields exact 0 diff & 0 SE!
  perm <- seq_len(nrow(d_base))
  perm[same_score_pair] <- rev(same_score_pair)
  d_perm_rn <- d_base[perm, , drop = FALSE]
  obj_perm_rn <- make_mock_cmp_fit(d = d_perm_rn, include_bias = FALSE, seed = 20261124L)
  loo_perm_rn <- loo::loo(obj_perm_rn, use_stan_log_lik = FALSE)
  cmp_perm_rn <- compare_bmfrm(a = loo_base, b = loo_perm_rn)
  expect_equal(as.numeric(cmp_perm_rn[2L, "elpd_diff"]), 0, tolerance = 1e-12)
  expect_equal(as.numeric(cmp_perm_rn[2L, "se_diff"]), 0, tolerance = 1e-12)

  # 3. Explicit obs_id column survives pre-NA filtering and row permutation cleanly
  d_with_id <- d_base
  d_with_id$rating_uid <- sprintf("OBS_%03d", seq_len(nrow(d_with_id)))
  # Insert an NA row that will be dropped before modeling
  d_with_na <- rbind(
    data.frame(person = "P1", item = "I1", rater = "R1", score = NA_integer_, rating_uid = "OBS_NA", stringsAsFactors = FALSE),
    d_with_id
  )
  d_with_id_perm <- d_with_id[rev(seq_len(nrow(d_with_id))), , drop = FALSE]
  rownames(d_with_id_perm) <- NULL

  obj_id_1 <- suppressMessages(make_mock_cmp_fit(d = d_with_na, obs_id = "rating_uid", seed = 20261124L))
  obj_id_2 <- make_mock_cmp_fit(d = d_with_id_perm, obs_id = "rating_uid", seed = 20261124L)
  expect_equal(obj_id_1$data_info$obs_id_source, "explicit_column")
  expect_equal(obj_id_1$data_info$orig_row_index, seq_len(nrow(d_with_id)) + 1L)

  loo_id_1 <- loo::loo(obj_id_1, use_stan_log_lik = FALSE)
  loo_id_2 <- loo::loo(obj_id_2, use_stan_log_lik = FALSE)
  cmp_id <- compare_bmfrm(m1 = loo_id_1, m2 = loo_id_2)
  expect_equal(as.numeric(cmp_id[2L, "elpd_diff"]), 0, tolerance = 1e-12)
  expect_equal(as.numeric(cmp_id[2L, "se_diff"]), 0, tolerance = 1e-12)

  # Also check grouped CV alignment under row permutation with explicit obs_id
  refit_fixed <- function(train_data, fold, fit_orig, ...) {
    make_mock_cmp_fit(d = train_data, include_bias = FALSE, seed = 900L + fold)
  }
  folds_orig <- ifelse(d_with_id$person %in% c("P1", "P2", "P3"), 1L, 2L)
  folds_rev <- rev(folds_orig)
  cv_1 <- suppressWarnings(cv_bmfrm(obj_id_1, target = "person", folds = folds_orig, M = 15L, seed = 55L, refit_fn = refit_fixed))
  cv_2 <- suppressWarnings(cv_bmfrm(obj_id_2, target = "person", folds = folds_rev, M = 15L, seed = 55L, refit_fn = refit_fixed))
  # pointwise$obs_id must expose the actual preserved obs_id (not 1:N) with row_index separate
  expect_identical(cv_1$pointwise$obs_id, d_with_id$rating_uid)
  expect_identical(cv_1$pointwise$row_index, seq_len(nrow(d_with_id)))
  expect_identical(cv_2$pointwise$obs_id, rev(d_with_id$rating_uid))
  expect_identical(cv_2$pointwise$row_index, seq_len(nrow(d_with_id)))
  cmp_cv_grp <- compare_bmfrm(m1 = cv_1, m2 = cv_2, score = "group")
  cmp_cv_pt <- compare_bmfrm(m1 = cv_1, m2 = cv_2, score = "pointwise")
  expect_s3_class(cmp_cv_grp, "bmfrm_model_comparison")
  expect_s3_class(cmp_cv_pt, "bmfrm_model_comparison")
  # Wrong fold assignment for the same obs_ids must be rejected
  expect_error(
    compare_bmfrm(m1 = cv_1, m2 = suppressWarnings(cv_bmfrm(obj_id_2, target = "person", folds = folds_orig, M = 15L, seed = 55L, refit_fn = refit_fixed))),
    "different fold splits"
  )

  # 4. Different pre-NA row indices (different NA rows dropped) must be rejected
  d_na_a <- rbind(data.frame(person = "P1", item = "I1", rater = "R1", score = NA_integer_, stringsAsFactors = FALSE), d_base)
  d_na_b <- rbind(d_base, data.frame(person = "P1", item = "I1", rater = "R1", score = NA_integer_, stringsAsFactors = FALSE))
  obj_na_a <- suppressMessages(make_mock_cmp_fit(d = d_na_a, seed = 20261124L))
  obj_na_b <- suppressMessages(make_mock_cmp_fit(d = d_na_b, seed = 20261124L))
  expect_error(
    compare_bmfrm(a = loo::loo(obj_na_a, use_stan_log_lik = FALSE),
                  b = loo::loo(obj_na_b, use_stan_log_lik = FALSE)),
    "different observation sets"
  )
})
