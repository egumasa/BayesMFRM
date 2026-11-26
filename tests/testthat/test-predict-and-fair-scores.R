# Tests for unified predict.bmfrm_fit() API and enhanced fair_scores()
# (weighted standardization, draw-matched decomposition, posterior draws)

make_mock_predict_fit <- function(D = 120L,
                                  family = "rsm",
                                  threshold_facet = NULL,
                                  hierarchical = c("person", "rater"),
                                  include_bias = FALSE,
                                  bias_shrinkage = "none",
                                  unbalanced = FALSE,
                                  custom_stan = FALSE,
                                  save_stan_gq = FALSE,
                                  seed = 20261124L) {
  set.seed(seed)
  if (!unbalanced) {
    d <- expand.grid(
      person = paste0("P", 1:4),
      item = c("I1", "I2", "I3"),
      rater = c("R1", "R2"),
      stringsAsFactors = FALSE
    )
  } else {
    # Unbalanced design: P1 rated by R1 on I1,I2; P2 by R2 on I2,I3; P3,P4 bridged
    d <- rbind(
      data.frame(person = "P1", item = c("I1", "I2"), rater = "R1", stringsAsFactors = FALSE),
      data.frame(person = "P2", item = c("I2", "I3"), rater = "R2", stringsAsFactors = FALSE),
      expand.grid(person = c("P3", "P4"), item = c("I1", "I2", "I3"), rater = c("R1", "R2"), stringsAsFactors = FALSE)
    )
  }
  d$score <- rep(c(1L, 2L, 3L, 4L), length.out = nrow(d))
  K <- 4L

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
  spec$eval_contract <- resolve_eval_contract(spec, stan_code = if (custom_stan) "parameters { vector[4] theta; }" else NULL)

  data_info <- prepare_data_bmfrm(spec, d)
  dc <- check_design_bmfrm(fml, d, K = K, family = family, threshold_facet = threshold_facet)

  # Construct deterministic-structured posterior draws
  theta <- cbind(
    stats::rnorm(D, -1.2, 0.15),
    stats::rnorm(D,  0.1, 0.15),
    stats::rnorm(D,  0.8, 0.15),
    stats::rnorm(D,  1.6, 0.15)
  )
  colnames(theta) <- paste0("theta[", 1:4, "]")

  item_raw <- cbind(
    stats::rnorm(D, -0.7, 0.08),
    stats::rnorm(D,  0.1, 0.08),
    stats::rnorm(D,  0.6, 0.08)
  )
  item_draws <- sweep(item_raw, 1L, rowMeans(item_raw), "-")
  colnames(item_draws) <- paste0("item[", 1:3, "]")

  r1 <- stats::rnorm(D, 0.9, 0.10)
  rater_draws <- cbind(r1, -r1)
  colnames(rater_draws) <- paste0("rater[", 1:2, "]")

  draw_store <- list(
    theta = theta,
    item = item_draws,
    rater = rater_draws
  )

  if ("person" %in% hierarchical) {
    draw_store$mu_theta <- matrix(stats::rnorm(D, 0.3, 0.08), nrow = D, ncol = 1L, dimnames = list(NULL, "mu_theta"))
    draw_store$sigma_theta <- matrix(stats::runif(D, 0.9, 1.1), nrow = D, ncol = 1L, dimnames = list(NULL, "sigma_theta"))
  }
  if ("rater" %in% hierarchical) {
    draw_store$sigma_rater <- matrix(stats::runif(D, 0.8, 1.0), nrow = D, ncol = 1L, dimnames = list(NULL, "sigma_rater"))
  }
  if ("item" %in% hierarchical) {
    draw_store$sigma_item <- matrix(stats::runif(D, 0.5, 0.7), nrow = D, ncol = 1L, dimnames = list(NULL, "sigma_item"))
  }

  if (identical(family, "rsm")) {
    tau_draws <- cbind(
      stats::rnorm(D, -1.0, 0.06),
      stats::rnorm(D,  0.0, 0.06),
      stats::rnorm(D,  1.0, 0.06)
    )
    tau_draws <- sweep(tau_draws, 1L, rowMeans(tau_draws), "-")
    colnames(tau_draws) <- paste0("tau[", 1:3, "]")
    draw_store$tau <- tau_draws
  } else {
    J_tf <- length(data_info$facet_labels[[threshold_facet]])
    tau_common <- cbind(
      stats::rnorm(D, -1.0, 0.05),
      stats::rnorm(D,  0.0, 0.05),
      stats::rnorm(D,  1.0, 0.05)
    )
    tau_common <- sweep(tau_common, 1L, rowMeans(tau_common), "-")
    colnames(tau_common) <- paste0("tau_common[", 1:3, "]")
    tau_mat <- matrix(0, nrow = D, ncol = J_tf * (K - 1L))
    cn <- character(J_tf * (K - 1L))
    for (k in 1:(K - 1L)) {
      for (j in seq_len(J_tf)) {
        col_i <- j + (k - 1L) * J_tf
        cn[col_i] <- paste0("tau[", j, ",", k, "]")
        tau_mat[, col_i] <- tau_common[, k] + (j - mean(seq_len(J_tf))) * 0.15 * (k - 2)
      }
    }
    colnames(tau_mat) <- cn
    draw_store$tau <- tau_mat
    if ("tau" %in% hierarchical) {
      draw_store$tau_common <- tau_common
      draw_store$sigma_tau <- matrix(stats::runif(D, 0.25, 0.35), nrow = D, ncol = 1L, dimnames = list(NULL, "sigma_tau"))
    }
  }

  if (include_bias) {
    # 3 x 2 double-centered bias matrix
    b_raw <- matrix(stats::rnorm(D * 6L, 0, 0.3), nrow = D, ncol = 6L)
    b_dc <- matrix(0, nrow = D, ncol = 6L)
    for (dd in seq_len(D)) {
      m <- matrix(b_raw[dd, ], nrow = 3L, ncol = 2L)
      m <- sweep(m, 1L, rowMeans(m), "-")
      m <- sweep(m, 2L, colMeans(m), "-")
      b_dc[dd, ] <- as.vector(m)
    }
    cn <- as.vector(outer(1:3, 1:2, function(i, j) paste0("itemrater[", i, ",", j, "]")))
    colnames(b_dc) <- cn
    draw_store$itemrater <- b_dc
    if ("item:rater" %in% hierarchical) {
      draw_store$sigma_itemrater <- matrix(stats::runif(D, 0.35, 0.45), nrow = D, ncol = 1L, dimnames = list(NULL, "sigma_itemrater"))
      if (bias_shrinkage == "horseshoe") {
        draw_store$c2_itemrater <- matrix(stats::runif(D, 1.5, 2.5), nrow = D, ncol = 1L, dimnames = list(NULL, "c2_itemrater"))
      }
    }
  }

  make_fit_obj <- function(store) {
    mock_cmdstan <- list(
      num_chains = function() 4L,
      metadata = function() list(stan_variables = names(store), iter_sampling = as.integer(D / 4L)),
      draws = function(var, format = "matrix") {
        if (!var %in% names(store)) stop("Missing variable: ", var)
        store[[var]]
      }
    )
    structure(
      list(fit = mock_cmdstan, spec = spec, data_info = data_info, design = dc, stan_code = NULL),
      class = "bmfrm_fit"
    )
  }

  obj <- make_fit_obj(draw_store)
  if (save_stan_gq) {
    q_full <- compute_bmfrm_quantities(obj, adjust_facets = character(0))
    q_fair <- compute_bmfrm_quantities(obj, adjust_facets = "rater")
    colnames(q_full$mu) <- paste0("mu[", seq_len(nrow(d)), "]")
    colnames(q_fair$mu) <- paste0("mu_fair[", seq_len(nrow(d)), "]")
    draw_store$mu <- q_full$mu
    draw_store$mu_fair <- q_fair$mu
    obj <- make_fit_obj(draw_store)
  }
  obj
}

test_that("predict.bmfrm_fit computes consistent epred, probs, predict, and linpred across chunks", {
  obj <- make_mock_predict_fit(D = 80L, include_bias = TRUE, hierarchical = c("person", "rater", "item:rater"))
  nd <- obj$data_info$data_clean[1:6, c("person", "item", "rater")]

  # Raw draws
  ep_draws <- predict(obj, newdata = nd, type = "epred", summary = FALSE)
  pr_draws <- predict(obj, newdata = nd, type = "probs", summary = FALSE)
  lp_draws <- predict(obj, newdata = nd, type = "linpred", summary = FALSE)
  yp_draws <- predict(obj, newdata = nd, type = "predict", summary = FALSE, seed = 42L)
  yp_draws2 <- predict(obj, newdata = nd, type = "predict", summary = FALSE, seed = 42L)

  expect_equal(dim(ep_draws), c(80L, 6L))
  expect_equal(dim(pr_draws), c(80L, 6L, 4L))
  expect_equal(dim(lp_draws), c(80L, 6L))
  expect_equal(dim(yp_draws), c(80L, 6L))
  expect_true(is.integer(yp_draws))
  expect_true(all(yp_draws >= 1L & yp_draws <= 4L))
  expect_identical(yp_draws, yp_draws2)

  # Exact mathematical identities: sum_k p_k = 1 and sum_k k * p_k = epred
  prob_sums <- unname(apply(pr_draws, c(1L, 2L), sum))
  expect_equal(prob_sums, matrix(1, 80L, 6L), tolerance = 1e-12, check.attributes = FALSE)
  ep_from_probs <- pr_draws[, , 1L] + 2 * pr_draws[, , 2L] + 3 * pr_draws[, , 3L] + 4 * pr_draws[, , 4L]
  expect_equal(unname(ep_from_probs), as.matrix(ep_draws), tolerance = 1e-12, check.attributes = FALSE)

  # Match internal compute_bmfrm_quantities()
  q_ref <- compute_bmfrm_quantities(obj, newdata = nd)
  expect_equal(as.matrix(ep_draws), q_ref$mu, tolerance = 1e-12, check.attributes = FALSE)

  # Summary mode & streaming chunk_size invariance
  ep_sum_full <- predict(obj, newdata = nd, type = "epred", summary = TRUE)
  ep_sum_chunk <- predict(obj, newdata = nd, type = "epred", summary = TRUE, chunk_size = 2L)
  expect_equal(ep_sum_chunk, ep_sum_full, tolerance = 1e-12)
  expect_equal(ep_sum_full$estimate, unname(colMeans(ep_draws)), tolerance = 1e-12)
  expect_equal(ep_sum_full$sd, unname(apply(ep_draws, 2L, stats::sd)), tolerance = 1e-12)
  expect_true(all(c("estimate", "sd", "q2.5", "q97.5") %in% names(ep_sum_full)))

  pr_sum_full <- predict(obj, newdata = nd, type = "probs", summary = TRUE)
  pr_sum_chunk <- predict(obj, newdata = nd, type = "probs", summary = TRUE, chunk_size = 2L)
  expect_equal(pr_sum_chunk, pr_sum_full, tolerance = 1e-12)
  expect_equal(names(pr_sum_full), paste0("prob_", 1:4))
  expect_equal(rowSums(as.matrix(pr_sum_full)), rep(1, 6L), tolerance = 1e-12)

  # ndraws / draw_ids and max_cells memory guard
  ep_sub <- predict(obj, newdata = nd, type = "epred", summary = FALSE, draw_ids = c(1L, 5L, 10L))
  expect_equal(as.matrix(ep_sub), ep_draws[c(1L, 5L, 10L), , drop = FALSE], check.attributes = FALSE)
  ep_nd <- predict(obj, newdata = nd, type = "epred", summary = FALSE, ndraws = 15L, seed = 7L)
  expect_equal(nrow(ep_nd), 15L)
  expect_error(predict(obj, newdata = nd, type = "probs", summary = FALSE, max_cells = 100), "exceeding max_cells")
  expect_error(predict(obj, newdata = nd, ndraws = 10L, draw_ids = 1:10), "not both")

  # Metadata contract
  meta <- attr(ep_sum_full, "bmfrm_prediction_meta")
  expect_equal(meta$prediction_target, "conditional_known_levels")
  expect_equal(meta$type, "epred")
})

test_that("predict.bmfrm_fit distinguishes zero-adjusted vs population-sampled new raters and guards non-hierarchical facets", {
  obj <- make_mock_predict_fit(D = 400L, hierarchical = c("person", "rater"))

  # 1. Default allow_new_levels = FALSE errors clearly
  nd_new_rater <- data.frame(person = c("P4", "P4"), item = c("I1", "I3"), rater = c("R_NEW", "R_NEW"), stringsAsFactors = FALSE)
  expect_error(predict(obj, newdata = nd_new_rater), "Unknown level\\(s\\) in facet 'rater'")

  # 2. Target C: adjust_facets = "rater" sets rater severity = 0 (rater column not even required!)
  nd_adj <- data.frame(person = c("P4", "P4"), item = c("I1", "I3"), stringsAsFactors = FALSE)
  pred_zero_adj <- predict(obj, newdata = nd_adj, type = "epred", adjust_facets = "rater", summary = FALSE)
  expect_equal(attr(pred_zero_adj, "bmfrm_prediction_meta")$prediction_target, "adjusted_reference_zero")

  # 3. Target B: population new rater drawn from N(0, sigma_rater)
  pred_pop_rater <- predict(
    obj, newdata = nd_new_rater, type = "epred",
    allow_new_levels = TRUE, sample_new_levels = "hierarchical",
    summary = FALSE, seed = 123L
  )
  meta_pop <- attr(pred_pop_rater, "bmfrm_prediction_meta")
  expect_equal(meta_pop$prediction_target, "population_new_levels")
  expect_equal(meta_pop$new_levels$rater, "R_NEW")

  # Strictly larger predictive SD under population rater sampling than under fixed rater = 0
  expect_gt(stats::sd(pred_pop_rater[, 1L]), stats::sd(pred_zero_adj[, 1L]) * 2)

  # Repeated rows with the same new rater ID share the EXACT same drawn rater severity per draw!
  lp_pop <- predict(
    obj, newdata = nd_new_rater, type = "linpred",
    allow_new_levels = TRUE, sample_new_levels = "hierarchical",
    summary = FALSE, seed = 123L
  )
  item_mat <- obj$fit$draws("item", format = "matrix")
  # eta(I1) - eta(I3) = -(item[I1] - item[I3]) exactly because theta[P4] and rater[R_NEW] cancel!
  expect_equal(unname(lp_pop[, 1L] - lp_pop[, 2L]), unname(-(item_mat[, 1L] - item_mat[, 3L])), tolerance = 1e-12)

  # Jensen's inequality: for high-ability P4 (+1.6 logits) on easy I1 (-0.7 logits),
  # E[mu(eta - rater_new)] < mu(eta - 0) because mu() is concave near upper bound 4!
  expect_lt(mean(pred_pop_rater[, 1L]), mean(pred_zero_adj[, 1L]) - 0.03)

  # 4. Non-hierarchical facet ('item') MUST error under sample_new_levels = "hierarchical"
  nd_new_item <- data.frame(person = "P1", item = "I_NEW", rater = "R1", stringsAsFactors = FALSE)
  expect_error(
    predict(obj, newdata = nd_new_item, allow_new_levels = TRUE, sample_new_levels = "hierarchical"),
    "non-hierarchical facet 'item'"
  )
  # Explicit sample_new_levels = "zero" succeeds and records zero_new_levels
  pred_zero_new <- predict(obj, newdata = nd_new_item, allow_new_levels = TRUE, sample_new_levels = "zero")
  expect_equal(attr(pred_zero_new, "bmfrm_prediction_meta")$prediction_target, "zero_new_levels")
})

test_that("predict.bmfrm_fit supports hierarchical interactions (Gaussian, Student-t, horseshoe), hierarchical PCM, and blocks custom_unverified", {
  # Hierarchical interactions with unseen rater
  for (sh in c("none", "student_t", "horseshoe")) {
    obj_b <- make_mock_predict_fit(
      D = 60L, include_bias = TRUE,
      hierarchical = c("person", "rater", "item:rater"),
      bias_shrinkage = sh
    )
    nd <- data.frame(
      person = c("P_NEW", "P_NEW"),
      item = c("I1", "I1"),
      rater = c("R_NEW", "R_NEW"),
      stringsAsFactors = FALSE
    )
    lp <- predict(obj_b, newdata = nd, type = "linpred",
                  allow_new_levels = TRUE, sample_new_levels = "hierarchical",
                  summary = FALSE, seed = 99L)
    # Identical unseen (person, item, rater) rows get identical draw within each iteration
    expect_equal(lp[, 1L], lp[, 2L], tolerance = 1e-12)
    expect_true(all(is.finite(lp)))
  }

  # Non-hierarchical interaction raises error if unseen rater is sampled hierarchically
  obj_nonhier_b <- make_mock_predict_fit(D = 40L, include_bias = TRUE, hierarchical = c("person", "rater"))
  nd_r <- data.frame(person = "P1", item = "I1", rater = "R_NEW", stringsAsFactors = FALSE)
  expect_error(
    predict(obj_nonhier_b, newdata = nd_r, allow_new_levels = TRUE, sample_new_levels = "hierarchical"),
    "non-hierarchical term 'item:rater'"
  )

  # Hierarchical PCM with unseen threshold_facet level
  obj_pcm <- make_mock_predict_fit(
    D = 60L, family = "pcm", threshold_facet = "rater",
    hierarchical = c("person", "rater", "tau")
  )
  pr_pcm <- predict(obj_pcm, newdata = nd_r, type = "probs",
                    allow_new_levels = TRUE, sample_new_levels = "hierarchical",
                    summary = TRUE, seed = 11L)
  expect_equal(sum(as.numeric(pr_pcm[1L, ])), 1, tolerance = 1e-12)

  # Non-hierarchical PCM errors on unseen threshold_facet under hierarchical sampling
  obj_pcm_fixed <- make_mock_predict_fit(
    D = 40L, family = "pcm", threshold_facet = "rater",
    hierarchical = c("person", "rater")
  )
  expect_error(
    predict(obj_pcm_fixed, newdata = nd_r, allow_new_levels = TRUE, sample_new_levels = "hierarchical"),
    "'tau' was not fit hierarchically"
  )

  # Unverified custom Stan model is blocked
  obj_cust <- make_mock_predict_fit(D = 40L, custom_stan = TRUE, hierarchical = character())
  expect_error(predict(obj_cust), "Unverified custom Stan model")
})

test_that("fair_scores supports weighted reference_grid, draw-matched decomposition, and return_draws contrasts", {
  obj <- make_mock_predict_fit(D = 120L, unbalanced = TRUE, hierarchical = c("person", "rater"))

  # 1. Three equivalent ways to specify 40/30/30 item blueprint weights
  grid_items <- data.frame(item = c("I1", "I2", "I3"), stringsAsFactors = FALSE)
  grid_with_w <- data.frame(item = c("I1", "I2", "I3"), weight = c(40, 30, 30), stringsAsFactors = FALSE)

  fs_vec <- fair_scores(obj, standardize = TRUE, reference_grid = grid_items,
                        reference_weights = c(0.4, 0.3, 0.3), return_draws = TRUE)
  fs_list <- fair_scores(obj, standardize = TRUE,
                         reference_weights = list(item = c(I1 = 0.4, I2 = 0.3, I3 = 0.3)),
                         return_draws = TRUE)
  fs_col <- fair_scores(obj, reference_grid = grid_with_w, return_draws = TRUE)
  fs_chunk <- fair_scores(obj, standardize = TRUE, reference_grid = grid_items,
                          reference_weights = c(0.4, 0.3, 0.3), chunk_size = 1L)

  expect_equal(fs_vec$mean_fair_score, fs_list$mean_fair_score, tolerance = 1e-12)
  expect_equal(fs_vec$mean_fair_score, fs_col$mean_fair_score, tolerance = 1e-12)
  expect_equal(fs_vec$mean_fair_score, fs_chunk$mean_fair_score, tolerance = 1e-12)
  expect_equal(attr(fs_vec, "draws")$fair, attr(fs_list, "draws")$fair, tolerance = 1e-12)

  # Verify weighted fair score matches manual weighted predict(..., adjust_facets = "rater")
  nd_p1 <- data.frame(person = "P1", item = c("I1", "I2", "I3"), stringsAsFactors = FALSE)
  ep_p1 <- predict(obj, newdata = nd_p1, type = "epred", adjust_facets = "rater", summary = FALSE)
  manual_p1 <- as.numeric(ep_p1 %*% c(0.4, 0.3, 0.3))
  expect_equal(unname(attr(fs_vec, "draws")$fair[, "P1"]), manual_p1, tolerance = 1e-12)

  # Weight validation errors
  expect_error(fair_scores(obj, standardize = TRUE, reference_weights = c(-0.1, 0.5, 0.6)), "non-negative")
  expect_error(fair_scores(obj, standardize = TRUE, reference_weights = list(item = c(I1 = 0.5, I99 = 0.5))), "unknown level")
  expect_error(fair_scores(obj, standardize = TRUE, reference_weights = list(item = c(I1 = 0.5, I2 = 0.5))), "missing weight")

  # 2. Exact draw-matched decomposition under observed design (standardize = FALSE)
  fs_obs_dec <- fair_scores(obj, standardize = FALSE, decompose = TRUE, return_draws = TRUE)
  expect_equal(fs_obs_dec$adjustment + fs_obs_dec$shrinkage_and_noise, fs_obs_dec$fair_minus_raw, tolerance = 1e-12)
  d_obs <- attr(fs_obs_dec, "draws")
  expect_equal(d_obs$adjustment, d_obs$fair - d_obs$actual, tolerance = 1e-12)
  expect_equal(fs_obs_dec$adjustment_sd, unname(apply(d_obs$adjustment, 2L, stats::sd)), tolerance = 1e-12)

  # Also works identically via saved-Stan mu + mu_fair path!
  obj_gq <- make_mock_predict_fit(D = 120L, unbalanced = TRUE, hierarchical = c("person", "rater"), save_stan_gq = TRUE)
  fs_gq_dec <- fair_scores(obj_gq, standardize = FALSE, decompose = TRUE, return_draws = TRUE)
  expect_equal(fs_gq_dec$adjustment, fs_obs_dec$adjustment, tolerance = 1e-12)
  expect_equal(fs_gq_dec$adjustment_sd, fs_obs_dec$adjustment_sd, tolerance = 1e-12)

  # 3. Exact 3-part draw-matched decomposition under weighted standardized design
  fs_std_dec <- fair_scores(
    obj, standardize = TRUE,
    reference_weights = list(item = c(I1 = 0.4, I2 = 0.3, I3 = 0.3)),
    decompose = TRUE, return_draws = TRUE
  )
  expect_equal(
    fs_std_dec$standardization_shift + fs_std_dec$adjustment + fs_std_dec$shrinkage_and_noise,
    fs_std_dec$fair_minus_raw,
    tolerance = 1e-12
  )
  # Unbalanced P1 (only rated on easier I1, I2) has non-zero standardization shift to I1/I2/I3
  expect_true(abs(fs_std_dec$standardization_shift[fs_std_dec$person == "P1"]) > 0.01)

  # Balanced equal-weight design has zero standardization_shift identically
  obj_bal <- make_mock_predict_fit(D = 60L, unbalanced = FALSE, hierarchical = c("person", "rater"))
  fs_bal_dec <- fair_scores(obj_bal, standardize = TRUE, decompose = TRUE)
  expect_equal(fs_bal_dec$standardization_shift, rep(0, 4L), tolerance = 1e-12)

  # 4. Draw-matched counterfactual contrast: same 40/30/30 items, R1 (harsh) vs R2 (lenient)
  g_r1 <- data.frame(item = c("I1", "I2", "I3"), rater = "R1", weight = c(0.4, 0.3, 0.3), stringsAsFactors = FALSE)
  g_r2 <- data.frame(item = c("I1", "I2", "I3"), rater = "R2", weight = c(0.4, 0.3, 0.3), stringsAsFactors = FALSE)
  fs_r1 <- fair_scores(obj, adjust_facets = character(0), reference_grid = g_r1, return_draws = TRUE)
  fs_r2 <- fair_scores(obj, adjust_facets = character(0), reference_grid = g_r2, return_draws = TRUE)
  diff_draws <- attr(fs_r2, "draws")$fair - attr(fs_r1, "draws")$fair
  # Since R1 is harsh (+0.9) and R2 is lenient (-0.9), R2 gives strictly higher expected scores on every draw
  expect_true(all(diff_draws > 0))
  # Guard errors when summary = FALSE
  expect_error(fair_scores(obj, summary = FALSE, decompose = TRUE), "requires summary = TRUE")
  expect_error(fair_scores(obj, summary = FALSE, return_draws = TRUE), "requires summary = TRUE")
})

test_that("unseen-rater interactions preserve sum-to-zero constraints, 1-level zero-df boundary, and existing cells", {
  for (sh in c("none", "student_t", "horseshoe")) {
    obj_b <- make_mock_predict_fit(
      D = 300L, include_bias = TRUE,
      hierarchical = c("person", "rater", "item:rater"),
      bias_shrinkage = sh, seed = 2026L
    )
    # Build matching object with sigma_itemrater = 0 to isolate pure interaction contribution
    obj_zero_b <- obj_b
    orig_draws <- obj_b$fit$draws
    obj_zero_b$fit$draws <- function(var, format = "matrix") {
      m <- orig_draws(var, format = format)
      if (identical(var, "sigma_itemrater")) m[] <- 0
      m
    }

    nd_new <- data.frame(
      person = rep("P1", 3L),
      item = c("I1", "I2", "I3"),
      rater = rep("R_NEW", 3L),
      stringsAsFactors = FALSE
    )
    lp_b <- predict(obj_b, newdata = nd_new, type = "linpred",
                    allow_new_levels = TRUE, sample_new_levels = "hierarchical",
                    summary = FALSE, seed = 55L)
    lp_0 <- predict(obj_zero_b, newdata = nd_new, type = "linpred",
                    allow_new_levels = TRUE, sample_new_levels = "hierarchical",
                    summary = FALSE, seed = 55L)
    b_col <- as.matrix(lp_b) - as.matrix(lp_0)

    # Exact sum-to-zero along known item axis for ALL 3 interaction modes!
    expect_lt(max(abs(rowSums(b_col))), 1e-13)
    # Non-trivial variance & negative within-column correlation
    expect_gt(stats::sd(b_col[, 1L]), 0.1)
    expect_lt(stats::cor(b_col[, 1L], b_col[, 2L]), -0.25)

    # Known rater R1 predictions are completely unchanged when R_NEW is also in newdata
    nd_mix <- rbind(
      data.frame(person = "P1", item = c("I1", "I2", "I3"), rater = "R1", stringsAsFactors = FALSE),
      nd_new
    )
    lp_known_only <- predict(obj_b, newdata = nd_mix[1:3, ], type = "linpred", summary = FALSE)
    lp_mix <- predict(obj_b, newdata = nd_mix, type = "linpred",
                      allow_new_levels = TRUE, sample_new_levels = "hierarchical",
                      summary = FALSE, seed = 55L)
    expect_equal(as.matrix(lp_mix)[, 1:3], as.matrix(lp_known_only), tolerance = 1e-13, check.attributes = FALSE)
  }

  # 1-level known item (0 df for item:rater): interaction for unseen rater MUST be identically 0!
  d1 <- expand.grid(person = paste0("P", 1:6), item = "I1", rater = c("R1", "R2"), stringsAsFactors = FALSE)
  d1$score <- rep(1:4, length.out = nrow(d1))
  spec1 <- parse_bmfrm_formula(score ~ person + item + rater + item:rater, d1, K = 4L,
                               hierarchical = c("person", "rater", "item:rater"))
  spec1 <- resolve_bmfrm_spec(spec1)
  spec1$eval_contract <- resolve_eval_contract(spec1)
  di1 <- prepare_data_bmfrm(spec1, d1)
  D1 <- 40L
  store1 <- list(
    theta = matrix(stats::rnorm(D1 * 6L), D1, 6L, dimnames = list(NULL, paste0("theta[", 1:6, "]"))),
    item = matrix(0, D1, 1L, dimnames = list(NULL, "item[1]")),
    rater = cbind(stats::rnorm(D1, 0.5, 0.1), stats::rnorm(D1, -0.5, 0.1)),
    itemrater = matrix(0, D1, 2L, dimnames = list(NULL, c("itemrater[1,1]", "itemrater[1,2]"))),
    tau = matrix(rep(c(-1, 0, 1), each = D1), D1, 3L, dimnames = list(NULL, paste0("tau[", 1:3, "]"))),
    mu_theta = matrix(0, D1, 1L, dimnames = list(NULL, "mu_theta")),
    sigma_theta = matrix(1, D1, 1L, dimnames = list(NULL, "sigma_theta")),
    sigma_rater = matrix(0.8, D1, 1L, dimnames = list(NULL, "sigma_rater")),
    sigma_itemrater = matrix(0.5, D1, 1L, dimnames = list(NULL, "sigma_itemrater"))
  )
  colnames(store1$rater) <- c("rater[1]", "rater[2]")
  obj1 <- structure(
    list(
      fit = list(metadata = function() list(stan_variables = names(store1)),
                 draws = function(var, format = "matrix") store1[[var]]),
      spec = spec1, data_info = di1, design = NULL, stan_code = NULL
    ),
    class = "bmfrm_fit"
  )
  obj1_zero <- obj1
  obj1_zero$fit$draws <- function(var, format = "matrix") {
    m <- store1[[var]]
    if (identical(var, "sigma_itemrater")) m[] <- 0
    m
  }
  nd1 <- data.frame(person = "P1", item = "I1", rater = "R_NEW", stringsAsFactors = FALSE)
  lp1_a <- predict(obj1, newdata = nd1, type = "linpred", allow_new_levels = TRUE,
                   sample_new_levels = "hierarchical", summary = FALSE, seed = 10L)
  lp1_b <- predict(obj1_zero, newdata = nd1, type = "linpred", allow_new_levels = TRUE,
                   sample_new_levels = "hierarchical", summary = FALSE, seed = 10L)
  expect_equal(as.matrix(lp1_a), as.matrix(lp1_b), tolerance = 0, check.attributes = FALSE)
})

test_that("zero-weight unidentified rows are flagged FALSE under summary = FALSE while summary = TRUE stays identified", {
  # Build a design where I1 and I2 are fully observed and identified, while factor level I3 is unobserved
  d <- expand.grid(
    person = paste0("P", 1:4),
    item = c("I1", "I2"),
    rater = c("R1", "R2"),
    stringsAsFactors = FALSE
  )
  d$item <- factor(d$item, levels = c("I1", "I2", "I3"))
  d$score <- rep(1:4, length.out = nrow(d))
  spec <- parse_bmfrm_formula(score ~ person + item + rater, d, K = 4L)
  spec <- resolve_bmfrm_spec(spec)
  spec$eval_contract <- resolve_eval_contract(spec)
  di <- prepare_data_bmfrm(spec, d)
  dc <- check_design_bmfrm(score ~ person + item + rater, d, K = 4L)

  D <- 40L
  store <- list(
    theta = matrix(0, D, 4L, dimnames = list(NULL, paste0("theta[", 1:4, "]"))),
    item = matrix(0, D, 3L, dimnames = list(NULL, paste0("item[", 1:3, "]"))),
    rater = matrix(0, D, 2L, dimnames = list(NULL, paste0("rater[", 1:2, "]"))),
    tau = matrix(rep(c(-1, 0, 1), each = D), D, 3L, dimnames = list(NULL, paste0("tau[", 1:3, "]")))
  )
  obj <- structure(
    list(
      fit = list(metadata = function() list(stan_variables = names(store)),
                 draws = function(var, format = "matrix") store[[var]]),
      spec = spec, data_info = di, design = dc, stan_code = NULL
    ),
    class = "bmfrm_fit"
  )

  ref_g <- data.frame(item = c("I1", "I2", "I3"), stringsAsFactors = FALSE)
  # Summary = TRUE uses only active positive-weight rows (I1, I2), so all persons are identified!
  fs_sum <- fair_scores(obj, standardize = TRUE, reference_grid = ref_g,
                        reference_weights = c(0.5, 0.5, 0), summary = TRUE)
  expect_true(all(fs_sum$likelihood_identified))

  # Summary = FALSE returns row-level identifiability: I3 row MUST be FALSE even though weight == 0!
  expect_warning(
    fs_row <- fair_scores(obj, standardize = TRUE, reference_grid = ref_g,
                          reference_weights = c(0.5, 0.5, 0), summary = FALSE),
    "likelihood-unidentified"
  )
  p1_rows <- fs_row[fs_row$person == "P1", ]
  expect_identical(p1_rows$likelihood_identified, c(TRUE, TRUE, FALSE))
})

test_that("predict and fair_scores preserve original draw_ids, rownames, chain, and iteration metadata", {
  obj <- make_mock_predict_fit(D = 120L, hierarchical = c("person", "rater"))
  nd <- obj$data_info$data_clean[1:2, c("person", "item", "rater")]

  sel <- c(5L, 2L, 31L, 95L)
  ep <- predict(obj, newdata = nd, type = "epred", draw_ids = sel, summary = FALSE)
  pr <- predict(obj, newdata = nd, type = "probs", draw_ids = sel, summary = FALSE)
  ep_sum <- predict(obj, newdata = nd, type = "epred", draw_ids = sel, summary = TRUE)

  expect_identical(rownames(ep), paste0("draw_", sel))
  expect_identical(dimnames(pr)[[1L]], paste0("draw_", sel))

  meta <- attr(ep, "bmfrm_prediction_meta")
  expect_identical(meta$draw_ids, sel)
  expect_identical(meta$draw_info$draw_id, sel)
  expect_identical(meta$draw_info$chain, c(1L, 1L, 2L, 4L))
  expect_identical(meta$draw_info$iteration, c(5L, 2L, 1L, 5L))
  expect_identical(attr(ep_sum, "bmfrm_prediction_meta")$draw_ids, sel)

  fs <- fair_scores(obj, standardize = TRUE, return_draws = TRUE)
  expect_identical(rownames(attr(fs, "draws")$fair)[1:3], c("draw_1", "draw_2", "draw_3"))
  expect_equal(nrow(attr(fs, "reference_meta")$draw_info), 120L)
  expect_identical(attr(fs, "reference_meta")$draw_info$chain[31L], 2L)
  expect_identical(attr(fs, "reference_meta")$draw_info$iteration[31L], 1L)
})

test_that("chunked new-level generation avoids upfront D x n_new allocation and is chunk_size invariant", {
  obj <- make_mock_predict_fit(D = 120L, hierarchical = c("person", "rater"))
  n_new <- 200L
  nd_many <- data.frame(
    person = paste0("P_NEW_", seq_len(n_new)),
    item = rep(c("I1", "I2", "I3"), length.out = n_new),
    rater = rep(c("R1", "R2"), length.out = n_new),
    stringsAsFactors = FALSE
  )

  # prepare_eval_plan must NOT expand theta to D x (4 + n_new)
  plan <- BayesMFRM:::prepare_eval_plan(
    obj, newdata = nd_many, adjust_facets = character(0), draw_idx = NULL,
    allow_new_levels = TRUE, sample_new_levels = "hierarchical"
  )
  expect_equal(ncol(plan$params$theta), 4L)
  expect_equal(length(plan$main_new$person$seeds), n_new)
  expect_lt(as.numeric(utils::object.size(plan$main_new)), 10000)

  # chunk_size = 1 vs chunk_size = 50 produce bit-for-bit identical predictions!
  p_c1 <- predict(obj, newdata = nd_many, type = "epred",
                  allow_new_levels = TRUE, sample_new_levels = "hierarchical",
                  summary = TRUE, chunk_size = 1L, seed = 77L)
  p_c50 <- predict(obj, newdata = nd_many, type = "epred",
                   allow_new_levels = TRUE, sample_new_levels = "hierarchical",
                   summary = TRUE, chunk_size = 50L, seed = 77L)
  expect_equal(p_c1, p_c50, tolerance = 1e-12)
})
