# 4-stage verification suite for Heterogeneous-Category PCM (v1):
#   Stage 1: Exact uniform-PCM parity (deterministic < 1e-12 on identical params/draws)
#   Stage 2: Mixed 2/4/6-category simulation, Stan codegen, posterior engine, predict, fair_scores, LOO/CV
#   Stage 3: Boundary, sparse/unbalanced, unobserved rubric categories, CV fold preservation, error guards

make_mock_hetero_pcm_fit <- function(d = NULL,
                                     K = c(I1 = 2L, I2 = 4L, I3 = 6L),
                                     D_iter = 30L,
                                     n_chains = 4L,
                                     hierarchical = c("person", "rater"),
                                     include_bias = FALSE,
                                     obs_id = NULL,
                                     seed = 20261124L) {
  set.seed(seed)
  if (is.null(d)) {
    d <- expand.grid(
      person = paste0("P", 1:6),
      item = names(K),
      rater = c("R1", "R2", "R3"),
      stringsAsFactors = FALSE
    )
    p_eff <- seq(-1, 1, length.out = 6L)[match(d$person, paste0("P", 1:6))]
    r_eff <- seq(-0.5, 0.5, length.out = 3L)[match(d$rater, c("R1", "R2", "R3"))]
    d$score <- vapply(seq_len(nrow(d)), function(i) {
      Ki <- as.integer(K[[d$item[i]]])
      val <- round(1 + (Ki - 1) * stats::plogis(p_eff[i] - r_eff[i] + stats::rnorm(1, 0, 0.4)))
      as.integer(max(1L, min(Ki, val)))
    }, integer(1))
  }
  D <- as.integer(D_iter * n_chains)

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
    fml, d, K = K, family = "pcm", threshold_facet = "item",
    hierarchical = hierarchical, save_gq = "minimal"
  )
  spec <- resolve_bmfrm_spec(spec, NULL)
  spec$eval_contract <- resolve_eval_contract(spec, stan_code = NULL)
  data_info <- prepare_data_bmfrm(spec, d, obs_id = obs_id)
  dc <- check_design_bmfrm(fml, data_info$data_clean, K = spec$K, family = "pcm", threshold_facet = "item")

  labels <- data_info$facet_labels
  J_p <- length(labels$person)
  J_i <- length(labels$item)
  J_r <- length(labels$rater)

  theta <- ar_draws(J_p, seq(-1, 1, length.out = J_p), 0.12)
  colnames(theta) <- paste0("theta[", seq_len(J_p), "]")

  item_raw <- ar_draws(J_i, seq(-0.4, 0.4, length.out = J_i), 0.08)
  item_draws <- if (J_i > 1L) sweep(item_raw, 1L, rowMeans(item_raw), "-") else matrix(0, D, 1L)
  colnames(item_draws) <- paste0("item[", seq_len(J_i), "]")

  rater_raw <- ar_draws(J_r, seq(-0.5, 0.5, length.out = J_r), 0.09)
  rater_draws <- if (J_r > 1L) sweep(rater_raw, 1L, rowMeans(rater_raw), "-") else matrix(0, D, 1L)
  colnames(rater_draws) <- paste0("rater[", seq_len(J_r), "]")

  if (isTRUE(spec$K_heterogeneous)) {
    K_lvl <- as.integer(unname(data_info$K_by_level[labels$item]))
    tau_blocks <- lapply(seq_len(J_i), function(j) {
      Kj <- K_lvl[j]
      if (Kj == 2L) {
        matrix(0, nrow = D, ncol = 1L)
      } else {
        tr <- ar_draws(Kj - 1L, seq(-0.6, 0.6, length.out = Kj - 1L), 0.07)
        sweep(tr, 1L, rowMeans(tr), "-")
      }
    })
    tau_draws <- do.call(cbind, tau_blocks)
    colnames(tau_draws) <- paste0("tau[", seq_len(ncol(tau_draws)), "]")
  } else {
    K_u <- get_spec_K_max(spec)
    tau_cols <- list()
    for (j in seq_len(J_i)) {
      tr <- ar_draws(K_u - 1L, seq(-0.5, 0.5, length.out = K_u - 1L), 0.07)
      tc <- sweep(tr, 1L, rowMeans(tr), "-")
      tau_cols[[j]] <- tc
    }
    # Column-major [j, k] order matching Stan matrix[J_item, K - 1]
    tau_draws <- matrix(0, nrow = D, ncol = J_i * (K_u - 1L))
    cn <- character(J_i * (K_u - 1L))
    pos <- 1L
    for (k in seq_len(K_u - 1L)) {
      for (j in seq_len(J_i)) {
        tau_draws[, pos] <- tau_cols[[j]][, k]
        cn[pos] <- paste0("tau[", j, ",", k, "]")
        pos <- pos + 1L
      }
    }
    colnames(tau_draws) <- cn
  }

  mu_theta <- matrix(stats::rnorm(D, 0, 0.1), ncol = 1L, dimnames = list(NULL, "mu_theta"))
  sigma_theta <- matrix(stats::rlnorm(D, 0, 0.1), ncol = 1L, dimnames = list(NULL, "sigma_theta"))
  sigma_rater <- matrix(stats::rlnorm(D, -0.5, 0.1), ncol = 1L, dimnames = list(NULL, "sigma_rater"))

  all_cols <- list(theta, item_draws, rater_draws, tau_draws, mu_theta, sigma_theta, sigma_rater)
  if (include_bias) {
    br <- array(ar_draws(J_i * J_r, 0, 0.15), dim = c(D, J_i, J_r))
    for (dd in seq_len(D)) {
      m <- br[dd, , ]
      m <- sweep(m, 1L, rowMeans(m), "-")
      m <- sweep(m, 2L, colMeans(m), "-")
      br[dd, , ] <- m
    }
    b_mat <- matrix(0, nrow = D, ncol = J_i * J_r)
    b_cn <- character(J_i * J_r)
    pos <- 1L
    for (r_i in seq_len(J_r)) {
      for (i_i in seq_len(J_i)) {
        b_mat[, pos] <- br[, i_i, r_i]
        b_cn[pos] <- paste0("itemrater[", i_i, ",", r_i, "]")
        pos <- pos + 1L
      }
    }
    colnames(b_mat) <- b_cn
    all_cols[[length(all_cols) + 1L]] <- b_mat
    if ("item:rater" %in% hierarchical) {
      all_cols[[length(all_cols) + 1L]] <- matrix(
        stats::rlnorm(D, -1.0, 0.1), ncol = 1L, dimnames = list(NULL, "sigma_itemrater")
      )
    }
  }

  mat <- do.call(cbind, all_cols)
  chain_vec <- rep(seq_len(n_chains), each = D_iter)
  iter_vec <- rep(seq_len(D_iter), times = n_chains)

  mock_fit <- list(
    draws = function(variables = NULL, format = "matrix") {
      if (is.null(variables)) {
        sel <- colnames(mat)
      } else {
        sel <- unlist(lapply(variables, function(v) {
          colnames(mat)[colnames(mat) == v | startsWith(colnames(mat), paste0(v, "["))]
        }), use.names = FALSE)
      }
      sub_m <- mat[, sel, drop = FALSE]
      if (identical(format, "df")) {
        df_out <- as.data.frame(sub_m)
        df_out$.chain <- chain_vec
        df_out$.iteration <- iter_vec
        df_out$.draw <- seq_len(D)
        return(df_out)
      }
      sub_m
    },
    metadata = function() list(
      stan_variables = unique(sub("\\[.*$", "", colnames(mat))),
      iter_sampling = D_iter
    ),
    num_chains = function() n_chains,
    summary = function(variables = NULL, ...) {
      sub_m <- mock_fit$draws(variables = variables, format = "matrix")
      dots <- list(...)
      want_sd <- any(vapply(dots, function(x) identical(x, "sd"), logical(1)))
      out <- tibble::tibble(
        variable = colnames(sub_m),
        mean = unname(colMeans(sub_m))
      )
      if (want_sd) {
        out$sd <- unname(apply(sub_m, 2L, stats::sd))
      }
      out$q2.5 <- unname(apply(sub_m, 2L, stats::quantile, probs = 0.025))
      out$q97.5 <- unname(apply(sub_m, 2L, stats::quantile, probs = 0.975))
      out
    }
  )

  structure(
    list(
      fit = mock_fit,
      spec = spec,
      data_info = data_info,
      design = dc,
      priors = spec$priors,
      hierarchical = spec$hierarchical,
      stan_code = build_stan_code(spec, spec$priors),
      convergence = list(
        converged = TRUE, status = "ok", max_rhat = 1.002,
        min_ess_bulk = 800, min_ess_tail = 750, min_ebfmi = 0.85,
        max_mcse_mean = 0.005, n_divergent = 0L, n_max_treedepth = 0L
      ),
      fit_args = list(
        formula = fml, K = spec$K, priors = spec$priors, family = "pcm",
        threshold_facet = "item", hierarchical = spec$hierarchical,
        save_gq = "minimal", chains = n_chains, iter = 1000L, warmup = 500L,
        seed = seed, obs_id = obs_id, sampler_args = list()
      )
    ),
    class = "bmfrm_fit"
  )
}

test_that("Stage 1: Uniform named-K PCM matches scalar-K PCM deterministically (< 1e-12) and keeps hierarchical tau", {
  set.seed(101)
  plan <- expand.grid(
    person = paste0("P", 1:6),
    item = c("I1", "I2", "I3"),
    rater = c("R1", "R2", "R3"),
    stringsAsFactors = FALSE
  )
  sim <- simulate_bmfrm(
    score ~ person + item + rater, design = plan, K = 4L,
    family = "pcm", threshold_facet = "item", hierarchical = c("person", "rater"), seed = 102L
  )
  d <- sim$data
  K_named <- c(I1 = 4L, I2 = 4L, I3 = 4L)

  # Design check exact parity
  dc_sc <- check_design_bmfrm(score ~ person + item + rater, d, K = 4L, family = "pcm", threshold_facet = "item")
  dc_nm <- check_design_bmfrm(score ~ person + item + rater, d, K = K_named, family = "pcm", threshold_facet = "item")
  expect_false(dc_nm$K_heterogeneous)
  expect_identical(dc_sc$adjacent_logit_rank, dc_nm$adjacent_logit_rank)
  expect_identical(dc_sc$step_parameters, dc_nm$step_parameters)
  expect_equal(dc_sc$singular_values, dc_nm$singular_values, tolerance = 1e-12)

  # Uniform named K still supports hierarchical = "tau"
  sp_hier <- parse_bmfrm_formula(
    score ~ person + item + rater, d, K = K_named,
    family = "pcm", threshold_facet = "item", hierarchical = c("person", "tau")
  )
  expect_false(sp_hier$K_heterogeneous)
  expect_true("tau" %in% sp_hier$hierarchical)
  expect_true(grepl("tau_common", build_stan_code(sp_hier), fixed = TRUE))

  # Identical draws: uniform named K vs scalar K -> exact < 1e-12 across predict, fair_scores, loo
  fit_sc <- make_mock_hetero_pcm_fit(d = d, K = 4L, seed = 303L)
  fit_nm <- make_mock_hetero_pcm_fit(d = d, K = K_named, seed = 303L)
  expect_identical(fit_sc$stan_code, fit_nm$stan_code)
  expect_equal(predict(fit_sc, type = "epred")$estimate, predict(fit_nm, type = "epred")$estimate, tolerance = 1e-12)
  expect_equal(as.matrix(predict(fit_sc, type = "probs")), as.matrix(predict(fit_nm, type = "probs")), tolerance = 1e-12)
  expect_equal(fair_scores(fit_sc, standardize = TRUE)$mean_fair_score,
               fair_scores(fit_nm, standardize = TRUE)$mean_fair_score, tolerance = 1e-12)

  loo_sc <- suppressWarnings(loo::loo(fit_sc, use_stan_log_lik = FALSE))
  loo_nm <- suppressWarnings(loo::loo(fit_nm, use_stan_log_lik = FALSE))
  expect_equal(loo_sc$pointwise[, "elpd_loo"], loo_nm$pointwise[, "elpd_loo"], tolerance = 1e-12)
  cmp_u <- compare_bmfrm(scalar_pcm = loo_sc, named_pcm = loo_nm)
  expect_equal(as.numeric(cmp_u[2L, "elpd_diff"]), 0, tolerance = 1e-12)
  expect_equal(as.numeric(cmp_u[2L, "se_diff"]), 0, tolerance = 1e-12)

  # Also verify eval_softmax_ragged directly against uniform eval_bmfrm_chunk on identical draws (< 1e-12)
  params_u <- BayesMFRM:::extract_bmfrm_param_draws(fit_sc)
  ids_u <- BayesMFRM:::map_eval_indices(d, fit_sc$spec, params_u$labels)
  ch_u <- BayesMFRM:::eval_bmfrm_chunk(
    params_u, ids_u, seq_len(nrow(d)), resp_chunk = d$score,
    compute_mu = TRUE, compute_linpred = TRUE, compute_probs = TRUE,
    compute_kurtosis = TRUE, compute_log_lik = TRUE
  )
  # Build matching pcm_ragged tau_info from the exact same cum_tau_list
  J_i <- length(params_u$labels$item)
  ct_by_lv <- lapply(seq_len(J_i), function(j) {
    do.call(cbind, lapply(seq_len(4L), function(k) params_u$tau_info$cum_tau_list[[k]][, j]))
  })
  tau_rag <- list(mode = "pcm_ragged", threshold_facet = "item", K_level = rep(4L, J_i), K_max = 4L, cum_tau_by_level = ct_by_lv)
  ch_r <- BayesMFRM:::eval_softmax_ragged(
    eta = ch_u$eta, tau_info = tau_rag, t_idx = ids_u$item, resp_chunk = d$score,
    compute_mu = TRUE, compute_probs = TRUE, compute_kurtosis = TRUE, compute_log_lik = TRUE
  )
  expect_equal(unname(ch_r$mu), unname(ch_u$mu), tolerance = 1e-12)
  expect_equal(unname(ch_r$sigma2), unname(ch_u$sigma2), tolerance = 1e-12)
  expect_equal(unname(ch_r$probs), unname(ch_u$probs), tolerance = 1e-12)
  expect_equal(unname(ch_r$kurt), unname(ch_u$kurt), tolerance = 1e-12)
  expect_equal(unname(ch_r$log_lik), unname(ch_u$log_lik), tolerance = 1e-12)
})

test_that("Stage 2: Mixed 2/4/6-category PCM supports simulation, marginal scaling, Stan codegen, predict, fair_scores, LOO, and CV", {
  K_mix <- c(I1 = 2L, I2 = 4L, I3 = 6L)
  plan <- expand.grid(
    person = paste0("P", 1:8),
    item = c("I1", "I2", "I3"),
    rater = c("R1", "R2", "R3"),
    stringsAsFactors = FALSE
  )

  # 1. Simulation & prior predictive check on mixed 2/4/6 rubrics
  priors_mix <- c(
    prior("normal(0, 1)", "mean", facet = "person"),
    prior("lognormal(0, 0.5)", "sd", facet = "person"),
    prior("normal(0, 0.8)", "item", scale = "marginal"),
    prior("normal(0, 0.5)", "sd", facet = "rater", scale = "marginal"),
    prior("normal(0, 0.8)", "tau", scale = "marginal")
  )
  sim <- simulate_bmfrm(
    score ~ person + item + rater, design = plan, K = K_mix,
    family = "pcm", threshold_facet = "item",
    hierarchical = c("person", "rater"),
    priors = priors_mix,
    nsim = 1L, seed = 404L
  )
  d <- sim$data
  expect_true(all(d$score[d$item == "I1"] %in% 1:2))
  expect_true(all(d$score[d$item == "I2"] %in% 1:4))
  expect_true(all(d$score[d$item == "I3"] %in% 1:6))
  tau_sim <- sim$parameters[[1L]]$tau
  expect_identical(unname(tau_sim$I1), 0)
  expect_equal(length(tau_sim$I2), 3L)
  expect_equal(sum(tau_sim$I2), 0, tolerance = 1e-12)
  expect_equal(length(tau_sim$I3), 5L)
  expect_equal(sum(tau_sim$I3), 0, tolerance = 1e-12)

  # Verify exact marginal prior SD = 1 across both K=4 and K=6 rows of B_tau_rows
  spec_m <- resolve_bmfrm_spec(
    parse_bmfrm_formula(
      score ~ person + item + rater, d, K = K_mix, family = "pcm",
      threshold_facet = "item", hierarchical = c("person", "rater")
    ),
    priors = priors_mix
  )
  di_m <- prepare_data_bmfrm(spec_m, d, obs_id = sprintf("OBS_%03d", seq_len(nrow(d))))
  B_rows <- di_m$stan_data$B_tau_rows
  # Row 1 (I1, K=2) -> all 0; Rows 2..4 (I2, K=4) & 5..9 (I3, K=6) -> row squared norms == 1.0!
  expect_equal(sum(B_rows[1L, ]^2), 0, tolerance = 1e-12)
  expect_equal(rowSums(B_rows[2:4, , drop = FALSE]^2), rep(1, 3L), tolerance = 1e-12)
  expect_equal(rowSums(B_rows[5:9, , drop = FALSE]^2), rep(1, 5L), tolerance = 1e-12)
  expect_identical(di_m$stan_data$N_tau_free, 0L + 2L + 4L)
  expect_identical(di_m$stan_data$N_tau_steps, 1L + 3L + 5L)

  # Also check binary level LAST (K = c(I1=6, I2=4, I3=2)) so 0-length end slice indices stay in bounds
  d_rev <- d
  d_rev$score <- ifelse(d_rev$item == "I1", pmin(d_rev$score, 6L),
                        ifelse(d_rev$item == "I2", pmin(d_rev$score, 4L), pmin(d_rev$score, 2L)))
  spec_rev <- resolve_bmfrm_spec(
    parse_bmfrm_formula(
      score ~ person + item + rater, d_rev, K = c(I1 = 6L, I2 = 4L, I3 = 2L),
      family = "pcm", threshold_facet = "item", hierarchical = c("person", "rater")
    ),
    priors = priors_mix
  )
  di_rev <- prepare_data_bmfrm(spec_rev, d_rev)
  expect_true(all(di_rev$stan_data$tau_free_start >= 1L & di_rev$stan_data$tau_free_start <= di_rev$stan_data$N_tau_free))
  expect_identical(di_rev$stan_data$tau_free_len, c(4L, 2L, 0L))

  # Stan codegen for both gaussian_basis & non-Gaussian student_t tau, and hierarchical interaction
  code_full <- build_stan_code(spec_m, save_gq = "full")
  code_min <- build_stan_code(spec_m, save_gq = "minimal")
  expect_true(grepl("B_tau_rows", code_full, fixed = TRUE))
  expect_true(grepl("generated quantities", code_full, fixed = TRUE))
  expect_false(grepl("generated quantities", code_min, fixed = TRUE))

  priors_ng <- c(
    prior("normal(0, 1)", "mean", facet = "person"),
    prior("lognormal(0, 0.5)", "sd", facet = "person"),
    prior("normal(0, 0.8)", "item", scale = "raw"),
    prior("normal(0, 0.5)", "sd", facet = "rater", scale = "raw"),
    prior("normal(0, 0.3)", "sd", facet = "item:rater", scale = "raw"),
    prior("student_t(3, 0, 1)", "tau", scale = "raw")
  )
  spec_ng <- resolve_bmfrm_spec(
    parse_bmfrm_formula(
      score ~ person + item + rater + item:rater, d, K = K_mix, family = "pcm",
      threshold_facet = "item", hierarchical = c("person", "rater", "item:rater")
    ),
    priors = priors_ng
  )
  code_ng <- build_stan_code(spec_ng, save_gq = "minimal")
  expect_true(grepl("tau_raw_len", code_ng, fixed = TRUE))
  expect_true(grepl("tau_raw ~ student_t(3, 0, 1);", code_ng, fixed = TRUE))

  # Diagnostic variable monitor excludes constant K_j == 2 tau[1]
  dvars <- BayesMFRM:::model_diagnostic_variables(spec_m, di_m)
  expect_false("tau[1]" %in% dvars)
  expect_true(all(paste0("tau[", 2:9, "]") %in% dvars))

  # 2. Posterior engine: predict(), residuals(), facet_fit(), fair_scores(), pp_check(), loo(), cv_bmfrm()
  d$oid <- sprintf("H%03d", seq_len(nrow(d)))
  fit_h <- make_mock_hetero_pcm_fit(d = d, K = K_mix, obs_id = "oid", seed = 505L)
  fit_hb <- make_mock_hetero_pcm_fit(d = d, K = K_mix, include_bias = TRUE,
                                     hierarchical = c("person", "rater", "item:rater"), obs_id = "oid", seed = 506L)

  # predict(type = "probs") -> K_obs column, NA_real_ for non-existent categories, row sum == 1 over valid categories
  p_sum <- predict(fit_h, type = "probs", summary = TRUE)
  expect_true("K_obs" %in% names(p_sum))
  expect_identical(p_sum$K_obs, as.integer(unname(K_mix[d$item])))
  i1_rows <- which(d$item == "I1")
  i2_rows <- which(d$item == "I2")
  i3_rows <- which(d$item == "I3")
  expect_true(all(is.na(as.matrix(p_sum[i1_rows, paste0("prob_", 3:6)]))))
  expect_equal(rowSums(as.matrix(p_sum[i1_rows, c("prob_1", "prob_2")])), rep(1, length(i1_rows)), tolerance = 1e-12)
  expect_true(all(is.na(as.matrix(p_sum[i2_rows, c("prob_5", "prob_6")]))))
  expect_equal(rowSums(as.matrix(p_sum[i2_rows, paste0("prob_", 1:4)])), rep(1, length(i2_rows)), tolerance = 1e-12)
  expect_false(anyNA(as.matrix(p_sum[i3_rows, paste0("prob_", 1:6)])))
  expect_equal(rowSums(as.matrix(p_sum[i3_rows, paste0("prob_", 1:6)])), rep(1, length(i3_rows)), tolerance = 1e-12)

  # predict(type = "probs", summary = FALSE) -> D x N x 6 with NA for c > K_j
  p_raw <- predict(fit_h, type = "probs", summary = FALSE, ndraws = 10L, seed = 1L)
  expect_identical(dim(p_raw), c(10L, nrow(d), 6L))
  expect_true(all(is.na(p_raw[, i1_rows, 3:6])))
  expect_false(anyNA(p_raw[, i1_rows, 1:2]))
  expect_identical(attr(p_raw, "bmfrm_prediction_meta")$K_obs, as.integer(unname(K_mix[d$item])))

  # predict(type = "predict") stays strictly within 1..K_j
  yrep_raw <- predict(fit_h, type = "predict", summary = FALSE, ndraws = 20L, seed = 2L)
  expect_true(all(yrep_raw[, i1_rows] %in% 1:2))
  expect_true(all(yrep_raw[, i2_rows] %in% 1:4))
  expect_true(all(yrep_raw[, i3_rows] %in% 1:6))

  # residuals() & facet_fit()
  res_h <- residuals(fit_h)
  expect_equal(nrow(res_h), nrow(d))
  expect_true(all(is.finite(res_h$z)))
  ff_i <- facet_fit(fit_h, "item")
  expect_equal(nrow(ff_i), 3L)
  expect_true(all(is.finite(ff_i$infit_msq)))

  # fair_scores(): reference_meta ALWAYS attached even with return_draws = FALSE!
  fs_unstd <- fair_scores(fit_h, standardize = FALSE, return_draws = FALSE)
  rm_u <- attr(fs_unstd, "reference_meta")
  expect_false(is.null(rm_u))
  expect_identical(rm_u$score_unit, "native_weighted_mean")
  expect_false(rm_u$standardized)
  expect_true(isTRUE(rm_u$K_heterogeneous))
  expect_identical(as.integer(rm_u$K_by_level[c("I1", "I2", "I3")]), c(2L, 4L, 6L))

  fs_std <- fair_scores(fit_h, standardize = TRUE, decompose = TRUE, return_draws = FALSE)
  rm_s <- attr(fs_std, "reference_meta")
  expect_false(is.null(rm_s))
  expect_true(rm_s$standardized)
  expect_identical(rm_s$score_unit, "native_weighted_mean")
  expect_equal(sum(rm_s$reference_weights), 1, tolerance = 1e-12)
  expect_false("expected_total_score" %in% names(fs_std))

  # pp_check() emits grouping recommendation message when ungrouped
  expect_message(
    suppressWarnings(bayesplot::pp_check(fit_h, type = "bars", ndraws = 10L, seed = 3L)),
    "Heterogeneous category counts"
  )

  # LOO & compare_bmfrm() between two heterogeneous PCMs with same K_by_level
  loo_h <- suppressWarnings(loo::loo(fit_h, use_stan_log_lik = FALSE))
  loo_hb <- suppressWarnings(loo::loo(fit_hb, use_stan_log_lik = FALSE))
  cmp_h <- compare_bmfrm(main = loo_h, bias = loo_hb)
  expect_s3_class(cmp_h, "bmfrm_model_comparison")

  # Grouped person CV preserves K_by_level & obs_id cleanly
  refit_h <- function(train_df, fold, spec, fit_args) {
    make_mock_hetero_pcm_fit(d = train_df, K = K_mix, seed = 700L + fold)
  }
  cv_h <- suppressWarnings(cv_bmfrm(fit_h, target = "person", K_folds = 2L, M = 15L, seed = 99L, refit_fn = refit_h))
  expect_true(cv_h$all_folds_identifiable)
  expect_identical(cv_h$pointwise$obs_id, d$oid)
  expect_identical(cv_h$pointwise$row_index, seq_len(nrow(d)))
  expect_true(is.finite(cv_h$elpd_group))
})

test_that("Stage 3: Boundary cases, unobserved rubric categories, CV guards, and v1 error rules", {
  K_mix <- c(I1 = 2L, I2 = 4L, I3 = 6L)
  plan <- expand.grid(
    person = paste0("P", 1:6),
    item = c("I1", "I2", "I3"),
    rater = c("R1", "R2", "R3"),
    stringsAsFactors = FALSE
  )
  set.seed(808)
  d <- plan
  # Intentionally keep I3 observed scores <= 4 even though rubric K["I3"] = 6!
  d$score <- ifelse(d$item == "I1", sample(1:2, nrow(d), replace = TRUE),
                    ifelse(d$item == "I2", sample(1:4, nrow(d), replace = TRUE),
                           sample(1:4, nrow(d), replace = TRUE)))

  # 1. Unobserved upper rubric categories (5, 6) preserved in design check, predict(), and CV
  dc <- check_design_bmfrm(score ~ person + item + rater, d, K = K_mix, family = "pcm", threshold_facet = "item")
  expect_identical(dc$step_parameters, 0L + 2L + 4L)
  expect_identical(dc$unobserved_categories$I3, c(5L, 6L))

  fit_unobs <- make_mock_hetero_pcm_fit(d = d, K = K_mix, seed = 809L)
  p_u <- predict(fit_unobs, type = "probs", summary = TRUE)
  i3_rows <- which(d$item == "I3")
  i1_rows <- which(d$item == "I1")
  # Existing-but-unobserved categories 5 & 6 have strictly positive probability on I3
  expect_true(all(p_u$prob_5[i3_rows] > 0 & p_u$prob_6[i3_rows] > 0))
  # Non-existent categories 3..6 remain NA_real_ on I1
  expect_true(all(is.na(as.matrix(p_u[i1_rows, paste0("prob_", 3:6)]))))

  # 2. All v1 validation & safety guards must error cleanly
  # Score out of range on K=2 item
  d_bad <- d
  d_bad$score[which(d_bad$item == "I1")[1L]] <- 3L
  expect_error(
    parse_bmfrm_formula(score ~ person + item + rater, d_bad, K = K_mix, family = "pcm", threshold_facet = "item"),
    "exceed valid category range"
  )

  # Vector K rejected for RSM
  expect_error(
    parse_bmfrm_formula(score ~ person + item + rater, d, K = K_mix, family = "rsm"),
    "only supported for family = 'pcm'"
  )

  # Missing / extra / unnamed K rejected
  expect_error(
    parse_bmfrm_formula(score ~ person + item + rater, d, K = c(2L, 4L, 6L), family = "pcm", threshold_facet = "item"),
    "named integer vector"
  )
  expect_error(
    parse_bmfrm_formula(score ~ person + item + rater, d, K = c(I1 = 2L, I2 = 4L), family = "pcm", threshold_facet = "item"),
    "missing level"
  )
  expect_error(
    parse_bmfrm_formula(score ~ person + item + rater, d, K = c(I1 = 2L, I2 = 4L, I3 = 6L, I9 = 3L), family = "pcm", threshold_facet = "item"),
    "unknown level"
  )

  # hierarchical = "tau" rejected when K is actually heterogeneous
  expect_error(
    parse_bmfrm_formula(score ~ person + item + rater, d, K = K_mix, family = "pcm",
                        threshold_facet = "item", hierarchical = c("person", "tau")),
    "Hierarchical step sharing"
  )

  # Adjusting threshold_facet in heterogeneous PCM rejected in fair_scores() & predict()
  expect_error(fair_scores(fit_unobs, adjust_facets = "item"), "heterogeneous-category PCM")
  expect_error(predict(fit_unobs, adjust_facets = "item"), "heterogeneous-category PCM")

  # Unseen threshold_facet level rejected even under allow_new_levels = TRUE AND sample_new_levels = "zero"
  new_item_df <- data.frame(person = "P1", item = "I_NEW", rater = "R1", stringsAsFactors = FALSE)
  expect_error(
    predict(fit_unobs, newdata = new_item_df, allow_new_levels = TRUE, sample_new_levels = "hierarchical"),
    "not supported in v1 heterogeneous-category PCM"
  )
  expect_error(
    predict(fit_unobs, newdata = new_item_df, allow_new_levels = TRUE, sample_new_levels = "zero"),
    "not supported in v1 heterogeneous-category PCM"
  )

  # Grouped CV holding out threshold_facet rejected in v1 heterogeneous PCM
  expect_error(
    cv_bmfrm(fit_unobs, target = "facet", facet = "item", K_folds = 2L),
    "not supported in v1 heterogeneous-category PCM"
  )

  # Pre-fit fold check stops if a person split drops a threshold_facet level from training!
  d_drop <- d
  # Make I3 appear ONLY for P1
  d_drop <- d_drop[!(d_drop$item == "I3" & d_drop$person != "P1"), , drop = FALSE]
  fit_drop <- make_mock_hetero_pcm_fit(d = d_drop, K = K_mix, seed = 810L)
  bad_folds <- ifelse(d_drop$person == "P1", 1L, 2L)
  expect_error(
    cv_bmfrm(fit_drop, target = "person", folds = bad_folds, M = 5L),
    "also makes other facet\\(s\\) unseen in training: item"
  )

  # compare_bmfrm() rejects K_obs mismatch (heterogeneous 2/4/6 PCM vs uniform K=6 RSM)
  d_rsm <- d
  spec_r6 <- resolve_bmfrm_spec(parse_bmfrm_formula(score ~ person + item + rater, d_rsm, K = 6L, family = "rsm"))
  di_r6 <- prepare_data_bmfrm(spec_r6, d_rsm)
  # Build a uniform K=6 PCM fit on the exact same data & compare against heterogeneous 2/4/6 PCM
  fit_u6 <- make_mock_hetero_pcm_fit(d = d, K = 6L, seed = 811L)
  loo_mix <- suppressWarnings(loo::loo(fit_unobs, use_stan_log_lik = FALSE))
  loo_u6 <- suppressWarnings(loo::loo(fit_u6, use_stan_log_lik = FALSE))
  expect_error(
    compare_bmfrm(hetero = loo_mix, uniform6 = loo_u6),
    "K_obs mismatch"
  )
})

test_that("Stage 3b: Regression tests for review boundary cases (CV K-mismatch, all-NA item, summary tau labels, strict_target=FALSE, fair_scores summary=FALSE)", {
  K_mix <- c(I1 = 2L, I2 = 4L, I3 = 6L)
  plan <- expand.grid(
    person = paste0("P", 1:6),
    item = c("I1", "I2", "I3"),
    rater = c("R1", "R2", "R3"),
    stringsAsFactors = FALSE
  )
  set.seed(909)
  d <- plan
  d$score <- ifelse(d$item == "I1", sample(1:2, nrow(d), replace = TRUE),
                    ifelse(d$item == "I2", sample(1:4, nrow(d), replace = TRUE),
                           sample(1:6, nrow(d), replace = TRUE)))

  # Regression 1: Uniform K=6 original model must reject custom refit_fn returning hetero K=c(2,4,6)
  fit_u6 <- make_mock_hetero_pcm_fit(d = d, K = 6L, seed = 910L)
  bad_refit_hetero <- function(train_df, fold, spec, fit_args) {
    make_mock_hetero_pcm_fit(d = train_df, K = K_mix, seed = 920L + fold)
  }
  expect_error(
    cv_bmfrm(fit_u6, target = "person", K_folds = 2L, M = 10L, seed = 911L, refit_fn = bad_refit_hetero),
    "changed per-level/per-observation category counts K"
  )
  # Also verify uniform RSM original model rejects custom refit_fn returning hetero PCM with same K_max=6
  fit_rsm6 <- fit_u6
  fit_rsm6$spec$family <- "rsm"
  fit_rsm6$spec$threshold_facet <- NULL
  fit_rsm6$spec$K_by_level <- NULL
  fit_rsm6$data_info$K_by_level <- NULL
  expect_error(
    cv_bmfrm(fit_rsm6, target = "person", K_folds = 2L, M = 10L, seed = 911L, refit_fn = bad_refit_hetero),
    "changed per-level/per-observation category counts K"
  )

  # Regression 2: All-NA level in a character column must NOT break scalar-K PCM or valid named-K PCM
  d_allna <- d
  d_allna$item <- as.character(d_allna$item)
  d_allna$score[d_allna$item == "I3"] <- NA_integer_

  # 2a. Scalar-K PCM (backward compatibility: retains all 36 non-NA rows)
  spec_sc_na <- parse_bmfrm_formula(
    score ~ person + item + rater, d_allna, K = 4, family = "pcm", threshold_facet = "item"
  )
  di_sc_na <- suppressMessages(prepare_data_bmfrm(spec_sc_na, d_allna))
  expect_identical(nrow(di_sc_na$data_clean), 36L)
  expect_identical(di_sc_na$facet_labels$item, c("I1", "I2"))
  expect_identical(di_sc_na$K_by_level, c(I1 = 4L, I2 = 4L))
  expect_identical(di_sc_na$stan_data$K, 4)

  # 2b. Heterogeneous named-K PCM where I3 has all NA scores (valid pre-NA rubric, mapped to I1 & I2 post-NA)
  spec_het_na <- parse_bmfrm_formula(
    score ~ person + item + rater, d_allna, K = K_mix, family = "pcm", threshold_facet = "item"
  )
  di_het_na <- suppressMessages(prepare_data_bmfrm(spec_het_na, d_allna))
  expect_identical(nrow(di_het_na$data_clean), 36L)
  expect_identical(di_het_na$facet_labels$item, c("I1", "I2"))
  expect_identical(di_het_na$K_by_level, c(I1 = 2L, I2 = 4L))
  expect_identical(di_het_na$stan_data$K_level, c(2L, 4L))
  expect_identical(di_het_na$stan_data$tau_start, c(1L, 2L))
  expect_identical(di_het_na$stan_data$tau_len, c(1L, 3L))

  # 2c. Genuine unknown level typo still errors at parse time
  expect_error(
    parse_bmfrm_formula(
      score ~ person + item + rater, d_allna,
      K = c(I1 = 2L, I2 = 4L, I3 = 6L, I_TYPO = 3L),
      family = "pcm", threshold_facet = "item"
    ),
    "unknown level"
  )

  # 2d. K = c(I1=2, I2=4, I3=6) with I2 & I3 all-NA normalizes to uniform binary K=2 PCM
  #     matching scalar K=2 byte-for-byte (both via parse_bmfrm_formula and via pre-built hetero spec)
  d_only_i1 <- d
  d_only_i1$item <- as.character(d_only_i1$item)
  d_only_i1$score[d_only_i1$item %in% c("I2", "I3")] <- NA_integer_

  spec_k2_sc <- resolve_bmfrm_spec(
    parse_bmfrm_formula(score ~ person + item + rater, d_only_i1, K = 2L, family = "pcm", threshold_facet = "item")
  )
  di_k2_sc <- suppressMessages(prepare_data_bmfrm(spec_k2_sc, d_only_i1))

  spec_k2_nm <- resolve_bmfrm_spec(
    parse_bmfrm_formula(score ~ person + item + rater, d_only_i1, K = K_mix, family = "pcm", threshold_facet = "item")
  )
  di_k2_nm <- suppressMessages(prepare_data_bmfrm(spec_k2_nm, d_only_i1))

  expect_false(isTRUE(spec_k2_nm$K_heterogeneous))
  expect_identical(spec_k2_nm$K, 2L)
  expect_identical(nrow(di_k2_nm$data_clean), sum(d$item == "I1"))
  expect_identical(di_k2_nm$K_by_level, c(I1 = 2L))
  expect_false(isTRUE(di_k2_nm$K_heterogeneous))
  expect_identical(di_k2_nm$stan_data, di_k2_sc$stan_data)
  expect_identical(build_stan_code(spec_k2_nm), build_stan_code(spec_k2_sc))
  expect_identical(dim(di_k2_nm$stan_data$B_tau), c(1L, 0L))
  expect_false("tau" %in% BayesMFRM:::model_diagnostic_variables(spec_k2_nm, di_k2_nm))

  # Also when a pre-built heterogeneous spec (from full data without NAs) is passed to prepare_data_bmfrm()
  spec_full_het <- resolve_bmfrm_spec(
    parse_bmfrm_formula(score ~ person + item + rater, d, K = K_mix, family = "pcm", threshold_facet = "item")
  )
  expect_true(isTRUE(spec_full_het$K_heterogeneous))
  di_k2_from_full <- suppressMessages(prepare_data_bmfrm(spec_full_het, d_only_i1))
  expect_false(isTRUE(di_k2_from_full$K_heterogeneous))
  expect_false(isTRUE(di_k2_from_full$spec$K_heterogeneous))
  expect_identical(di_k2_from_full$K_by_level, c(I1 = 2L))
  expect_identical(di_k2_from_full$stan_data, di_k2_sc$stan_data)
  expect_identical(build_stan_code(di_k2_from_full$spec), build_stan_code(spec_k2_sc))

  # 2e. Multi-level post-NA uniform case: K = c(I1=4, I2=4, I3=6) with I3 all-NA -> uniform K=4 PCM
  d_446 <- d
  d_446$item <- as.character(d_446$item)
  d_446$score[d_446$item == "I3"] <- NA_integer_
  spec_44_sc <- resolve_bmfrm_spec(
    parse_bmfrm_formula(score ~ person + item + rater, d_446, K = 4L, family = "pcm", threshold_facet = "item")
  )
  di_44_sc <- suppressMessages(prepare_data_bmfrm(spec_44_sc, d_446))
  spec_44_nm <- resolve_bmfrm_spec(
    parse_bmfrm_formula(score ~ person + item + rater, d_446, K = c(I1 = 4L, I2 = 4L, I3 = 6L),
                        family = "pcm", threshold_facet = "item")
  )
  di_44_nm <- suppressMessages(prepare_data_bmfrm(spec_44_nm, d_446))
  expect_false(isTRUE(spec_44_nm$K_heterogeneous))
  expect_identical(di_44_nm$stan_data, di_44_sc$stan_data)
  expect_identical(build_stan_code(spec_44_nm), build_stan_code(spec_44_sc))

  # Regression 3: summary() decodes packed tau[1..9] cleanly into (item=..., step=...) with zero NAs
  fit_h <- make_mock_hetero_pcm_fit(d = d, K = K_mix, seed = 930L)
  summ_lines <- utils::capture.output(summary(fit_h))
  summ_txt <- paste(summ_lines, collapse = "\n")
  expect_true(grepl("tau\\[1\\] \\(item=I1, step=1\\)", summ_txt))
  expect_true(grepl("tau\\[2\\] \\(item=I2, step=1\\)", summ_txt))
  expect_true(grepl("tau\\[3\\] \\(item=I2, step=2\\)", summ_txt))
  expect_true(grepl("tau\\[4\\] \\(item=I2, step=3\\)", summ_txt))
  expect_true(grepl("tau\\[5\\] \\(item=I3, step=1\\)", summ_txt))
  expect_true(grepl("tau\\[9\\] \\(item=I3, step=5\\)", summ_txt))
  expect_false(grepl("item=NA|step=NA", summ_txt))

  # Regression 4: strict_target = FALSE must NOT bypass pre-fit unseen threshold_facet guard (0 refit calls!)
  d_drop <- d[!(d$item == "I3" & d$person != "P1"), , drop = FALSE]
  fit_drop <- make_mock_hetero_pcm_fit(d = d_drop, K = K_mix, seed = 940L)
  bad_folds <- ifelse(d_drop$person == "P1", 1L, 2L)
  refit_calls <- 0L
  spy_refit <- function(train_df, fold, spec, fit_args) {
    refit_calls <<- refit_calls + 1L
    stop("spy_refit should never have been invoked")
  }
  expect_error(
    cv_bmfrm(fit_drop, target = "person", folds = bad_folds, strict_target = FALSE, M = 5L, refit_fn = spy_refit),
    "not supported in v1 heterogeneous-category PCM"
  )
  expect_identical(refit_calls, 0L)

  # Regression 5: fair_scores(summary = FALSE) preserves reference_meta on ALL 3 return paths
  # 5a. Standardized summary = FALSE
  fs_std_row <- fair_scores(fit_h, standardize = TRUE, summary = FALSE)
  rm_std_row <- attr(fs_std_row, "reference_meta")
  expect_false(is.null(rm_std_row))
  expect_true(isTRUE(rm_std_row$standardized))
  expect_false(isTRUE(rm_std_row$summary))
  expect_identical(rm_std_row$score_unit, "native_row_expected_score")
  expect_true(isTRUE(rm_std_row$K_heterogeneous))
  expect_identical(as.integer(rm_std_row$K_by_level[c("I1", "I2", "I3")]), c(2L, 4L, 6L))
  expect_equal(sum(rm_std_row$reference_weights), 1, tolerance = 1e-12)

  # 5b. Unstandardized R-streaming summary = FALSE
  fs_unstd_row <- fair_scores(fit_h, standardize = FALSE, summary = FALSE)
  rm_unstd_row <- attr(fs_unstd_row, "reference_meta")
  expect_false(is.null(rm_unstd_row))
  expect_false(isTRUE(rm_unstd_row$standardized))
  expect_false(isTRUE(rm_unstd_row$summary))
  expect_identical(rm_unstd_row$score_unit, "native_row_expected_score")
  expect_true(isTRUE(rm_unstd_row$K_heterogeneous))
  expect_identical(as.integer(rm_unstd_row$K_by_level[c("I1", "I2", "I3")]), c(2L, 4L, 6L))

  # 5c. Unstandardized Stan-fast (saved mu_fair) summary = FALSE
  fit_h_stan <- fit_h
  q_fair <- BayesMFRM:::compute_bmfrm_quantities(fit_h, adjust_facets = "rater")
  mu_fair_mat <- q_fair$mu
  colnames(mu_fair_mat) <- paste0("mu_fair[", seq_len(ncol(mu_fair_mat)), "]")
  orig_draws_fn <- fit_h$fit$draws
  orig_meta_fn <- fit_h$fit$metadata
  fit_h_stan$fit$draws <- function(variables = NULL, format = "matrix") {
    if (identical(variables, "mu_fair")) return(mu_fair_mat)
    orig_draws_fn(variables = variables, format = format)
  }
  fit_h_stan$fit$metadata <- function() {
    md <- orig_meta_fn()
    md$stan_variables <- unique(c(md$stan_variables, "mu_fair"))
    md
  }
  fs_stan_row <- fair_scores(fit_h_stan, standardize = FALSE, summary = FALSE)
  rm_stan_row <- attr(fs_stan_row, "reference_meta")
  expect_false(is.null(rm_stan_row))
  expect_false(isTRUE(rm_stan_row$standardized))
  expect_false(isTRUE(rm_stan_row$summary))
  expect_identical(rm_stan_row$score_unit, "native_row_expected_score")
  expect_true(isTRUE(rm_stan_row$K_heterogeneous))
  expect_identical(as.integer(rm_stan_row$K_by_level[c("I1", "I2", "I3")]), c(2L, 4L, 6L))
})

test_that("mixed-rubric curves and ridges agree with independently indexed adjacent odds", {
  K <- c(I1 = 2L, I2 = 4L, I3 = 6L)
  obj <- make_mock_hetero_pcm_fit(K = K, include_bias = TRUE)
  profiles <- data.frame(item = c("I3", "I1", "I2"), rater = c("R1", "R2", "R3"))
  grid <- seq(-4, 4, length.out = 33)
  curves <- category_curves(obj, profiles, theta = grid)
  cuts <- category_thresholds(obj, profiles)
  expect_equal(nrow(cuts), 9L)
  all_draws <- obj$fit$draws()
  tau_start <- c(1L, 2L, 5L)
  for (i in 1:3) {
    j <- match(profiles$item[i], names(K))
    r <- match(profiles$rater[i], obj$data_info$facet_labels$rater)
    Kj <- K[[j]]
    s <- seq_len(Kj - 1L)
    pos <- all_draws[, sprintf("item[%d]", j)] + all_draws[, sprintf("rater[%d]", r)] -
      all_draws[, sprintf("itemrater[%d,%d]", j, r)]
    steps <- all_draws[, sprintf("tau[%d]", tau_start[j] + s - 1L), drop = FALSE] + pos
    got <- curves[curves$profile == i, ]
    expect_equal(sort(unique(got$category)), seq_len(Kj))
    expect_equal(cuts$median[cuts$profile == i], unname(apply(steps, 2, median)), tolerance = 1e-12)
    # The probability ratio between neighbors is exp(theta - threshold).
    ref_prob <- function(ability) {
      w <- matrix(1, nrow(steps), Kj)
      for (k in s) w[, k + 1L] <- w[, k] * exp(ability - steps[, k])
      w / rowSums(w)
    }
    expected <- t(vapply(grid, function(x) colMeans(ref_prob(x)), numeric(Kj)))
    expect_equal(got$probability, as.vector(expected), tolerance = 1e-12)
    expect_equal(rowSums(expected), rep(1, length(grid)), tolerance = 1e-12)
    # Compare the draw-level crossings used for the ridges, including every binary step.
    ct <- bmfrm_profile_cumulative(bmfrm_curve_profiles(obj, profiles, 0.95), i)
    actual_steps <- bmfrm_step_draws(ct)
    for (k in s) {
      at_crossing <- ref_prob(actual_steps[, k])
      expect_lt(max(abs(at_crossing[, k] - at_crossing[, k + 1L])), 1e-12)
    }
  }
  expect_no_warning(ggplot2::ggplot_build(plot_wright_map(obj, profiles)))
  for (stat in c("infit", "outfit")) {
    f <- facet_fit(obj, "rater", seed = 17)
    m <- facet_summary(obj, "rater", probs = c(0.025, 0.975))
    p <- plot_pathway_map(obj, "rater", statistic = stat, seed = 17)
    expect_equal(p$data$level, m$level)
    expect_equal(p$data$lower, m$q2.5)
    expect_equal(p$data$msq, f[[paste0(stat, "_msq")]])
    expect_equal(p$data$rep_lower, f[[paste0(stat, "_rep_lo")]])
    expect_equal(p$data$rep_upper, f[[paste0(stat, "_rep_hi")]])
    expect_equal(p$data$ppp, f[[paste0(stat, "_ppp")]])
    expect_no_warning(ggplot2::ggplot_build(p))
  }
})
