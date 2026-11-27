# Target-aware Bayesian model comparison for BayesMFRM:
#   1. loo.bmfrm_fit(): low-memory PSIS-LOO with chain-aware exp(log_lik) r_eff,
#      singleton-level verification, native loo diagnostics, and eval-contract guard.
#   2. cv_bmfrm(): target-first K-fold CV supporting rating-level refits and
#      grouped unseen-level Monte Carlo marginalization (exact double-logmeanexp),
#      dropleveled training prep, per-fold design + convergence diagnostics.
#   3. compare_bmfrm(): guarded model comparison verifying prediction target,
#      evaluation data/encoding, fold split, score unit, and cluster-aware SEs.

#' Compute canonical evaluation fingerprint for comparison compatibility checks
#' @keywords internal
build_comparison_key <- function(method, target, target_facet = NULL, actual_target,
                                 known_facets, unseen_facets,
                                 data_clean, response_col, K,
                                 K_folds = NULL, fold_ids = NULL,
                                 group_ids = NULL, aggregation = "sum",
                                 data_info = NULL, facets_main = NULL,
                                 K_by_level = NULL, threshold_facet = NULL) {
  N <- nrow(data_clean)
  if (!is.null(data_info) && !is.null(data_info$obs_id) && length(data_info$obs_id) == N) {
    obs_ids <- as.character(data_info$obs_id)
    obs_id_source <- data_info$obs_id_source %||% "data_info"
  } else {
    rn <- rownames(data_clean)
    has_custom_rn <- !is.null(rn) &&
      length(rn) == N &&
      !anyNA(rn) &&
      !anyDuplicated(rn) &&
      !identical(rn, as.character(seq_len(N)))
    if (has_custom_rn) {
      obs_ids <- as.character(rn)
      obs_id_source <- "rownames"
    } else {
      obs_ids <- as.character(seq_len(N))
      obs_id_source <- "fallback_seq"
    }
  }
  if (anyNA(obs_ids) || anyDuplicated(obs_ids)) {
    stop("Observation IDs used for model comparison must be unique and non-missing")
  }

  ord <- order(obs_ids)
  obs_ids_sorted <- obs_ids[ord]
  resp_int <- as.integer(data_clean[[response_col]])

  # Per-observation rubric category count vector K_obs (uniform models -> rep(K, N))
  eff_K_by_level <- K_by_level %||% if (!is.null(data_info)) data_info$K_by_level else NULL
  if (!is.null(eff_K_by_level) && !is.null(threshold_facet) && threshold_facet %in% names(data_clean)) {
    K_obs <- as.integer(unname(eff_K_by_level[as.character(data_clean[[threshold_facet]])]))
  } else {
    K_obs <- rep(as.integer(max(K)), N)
  }
  K_max_int <- as.integer(max(K_obs))

  if (is.null(facets_main)) {
    facets_main <- sort(unique(c(as.character(known_facets), as.character(unseen_facets))))
  } else {
    facets_main <- sort(unique(as.character(facets_main)))
  }
  facets_in_data <- intersect(facets_main, names(data_clean))
  facet_hashes_by_id <- stats::setNames(
    lapply(facets_in_data, function(f) {
      digest::digest(list(obs_ids = obs_ids_sorted, values = as.character(data_clean[[f]])[ord]))
    }),
    facets_in_data
  )

  # Also hash any extra atomic/factor non-response columns present in data_clean
  # (excluding facet columns already hashed above) to catch row swaps in shared metadata.
  extra_cols <- setdiff(names(data_clean), c(response_col, facets_in_data))
  extra_cols <- extra_cols[vapply(extra_cols, function(cn) {
    v <- data_clean[[cn]]
    is.atomic(v) || is.factor(v)
  }, logical(1))]
  extra_col_hashes_by_id <- stats::setNames(
    lapply(sort(extra_cols), function(cn) {
      digest::digest(list(obs_ids = obs_ids_sorted, values = as.character(data_clean[[cn]])[ord]))
    }),
    sort(extra_cols)
  )

  fold_ids_int <- if (is.null(fold_ids)) NULL else as.integer(fold_ids)
  group_ids_chr <- if (is.null(group_ids)) NULL else as.character(group_ids)

  list(
    method = method,
    target = target,
    target_facet = target_facet,
    actual_target = actual_target,
    known_facets = sort(as.character(known_facets)),
    unseen_facets = sort(as.character(unseen_facets)),
    facets_main = facets_main,
    n_obs = N,
    obs_ids = obs_ids,
    obs_id_source = obs_id_source,
    obs_id_set_hash = digest::digest(obs_ids_sorted),
    obs_id_order_hash = digest::digest(obs_ids),
    K = K_max_int,
    K_obs_hash_by_id = digest::digest(list(obs_ids = obs_ids_sorted, K_obs = K_obs[ord])),
    K_heterogeneous = length(unique(K_obs)) > 1L,
    response_col = as.character(response_col),
    response_hash = digest::digest(list(obs_ids = obs_ids_sorted, response = resp_int[ord], K_obs = K_obs[ord])),
    facet_hashes_by_id = facet_hashes_by_id,
    extra_col_hashes_by_id = extra_col_hashes_by_id,
    K_folds = if (is.null(K_folds)) NULL else as.integer(K_folds),
    fold_ids = fold_ids_int,
    fold_hash_by_id = if (is.null(fold_ids_int)) NULL else digest::digest(list(obs_ids = obs_ids_sorted, fold = fold_ids_int[ord])),
    group_ids = group_ids_chr,
    group_hash_by_id = if (is.null(group_ids_chr)) NULL else digest::digest(list(obs_ids = obs_ids_sorted, group = group_ids_chr[ord])),
    aggregation = aggregation
  )
}

#' Detect observations that are the sole observation of any main facet level
#' @keywords internal
detect_singleton_observations <- function(data_clean, facets_main) {
  singleton_facets <- list()
  singleton_mask <- logical(nrow(data_clean))
  for (f in facets_main) {
    vals <- as.character(data_clean[[f]])
    tab <- table(vals)
    sing_lv <- names(tab)[tab == 1L]
    if (length(sing_lv)) {
      singleton_facets[[f]] <- sing_lv
      singleton_mask <- singleton_mask | (vals %in% sing_lv)
    }
  }
  list(
    singleton_obs = which(singleton_mask),
    singleton_facets = singleton_facets
  )
}

#' Drop unused factor levels on all model facet columns so training folds rebuild clean indices
#' @keywords internal
drop_unused_facet_levels <- function(df, facets_main) {
  out <- df
  for (f in facets_main) {
    if (f %in% names(out)) {
      vals <- as.character(out[[f]])
      out[[f]] <- factor(vals, levels = unique(vals))
    }
  }
  rownames(out) <- NULL
  out
}

#' Compute PSIS-LOO cross-validation for a fitted BayesMFRM model
#'
#' Evaluates pointwise leave-one-rating-out expected log predictive density
#' (`elpd_loo`) via [loo::loo()]. Uses Stan's saved `log_lik` array when
#' available (`save_gq = "full"`), or streams pointwise log-likelihoods on
#' demand via `loo::loo.function()` when `save_gq = "minimal"` without
#' allocating a full `D x N` matrix.
#'
#' The rating-level prediction task is `"new_rating_known_levels"`: leaving out
#' one rating while all main facet levels remain observed in training. Sole
#' observations of any main facet level are flagged with
#' `actual_target = "new_rating_with_singleton_levels"`, because leaving them
#' out removes that facet level from training.
#'
#' Relative effective sample sizes (`r_eff`) are computed on likelihoods
#' `exp(log_lik)` preserving the original iteration-by-chain layout; missing
#' chain structure raises an error rather than silently defaulting to `r_eff = 1`.
#' When Pareto-k exceeds the diagnostic threshold (k > 0.7), verify the
#' contribution by exact leave-one-out refitting or rating-level [cv_bmfrm()]
#' rather than assuming its error is negligible.
#'
#' @param x A `bmfrm_fit` object.
#' @param r_eff Optional numeric vector of relative effective sample sizes
#'   (length `N`) or `NA` to deliberately omit `r_eff` adjustment. When `NULL`
#'   (default), computed automatically from `exp(log_lik)` using the verified
#'   iteration-by-chain layout.
#' @param save_psis Logical; passed to [loo::loo()]. Default `FALSE`.
#' @param cores Number of cores for `loo` computations.
#' @param chunk_size Optional observation chunk size when computing `r_eff` via
#'   the streaming R engine (auto-sized to ~16 MB when `NULL`).
#' @param use_stan_log_lik Optional logical; `TRUE` forces Stan's saved
#'   `log_lik`, `FALSE` forces the verified R streaming engine, `NULL`
#'   (default) uses Stan `log_lik` when saved and falls back to R streaming
#'   otherwise.
#' @param ... Additional arguments passed to [loo::loo()].
#' @return A `bmfrm_loo` object (inheriting from `psis_loo` and `loo`) with
#'   `bmfrm_meta` (target verification, Pareto-k diagnostics, PSIS ESS, MCSE,
#'   and original MCMC convergence) and `comparison_key`.
#' @importFrom loo loo
#' @exportS3Method loo::loo bmfrm_fit
loo.bmfrm_fit <- function(x,
                          r_eff = NULL,
                          save_psis = FALSE,
                          cores = getOption("mc.cores", 1L),
                          chunk_size = NULL,
                          use_stan_log_lik = NULL,
                          ...) {
  if (!inherits(x, "bmfrm_fit")) stop("x must be a bmfrm_fit object")
  obs_data <- x$data_info$data_clean
  spec <- x$spec
  N <- nrow(obs_data)

  # 1. Verify whether leaving out any single rating creates an unseen facet level
  sing <- detect_singleton_observations(obs_data, spec$facets_main)
  actual_target <- if (length(sing$singleton_obs)) {
    "new_rating_with_singleton_levels"
  } else {
    "new_rating_known_levels"
  }
  if (length(sing$singleton_obs)) {
    warning(
      length(sing$singleton_obs), " observation(s) are the sole rating for a level of facet(s) (",
      paste(names(sing$singleton_facets), collapse = ", "),
      "). Leave-one-out for those rows leaves 0 training observations for that level; ",
      "actual_target recorded as 'new_rating_with_singleton_levels'.",
      call. = FALSE
    )
  }

  # 2. Check original MCMC convergence
  conv <- x$convergence
  if (!is.null(conv) && !isTRUE(conv$converged)) {
    warning(
      "Original fit did not satisfy MCMC convergence checks (status = '",
      conv$status %||% "unknown", "'); PSIS-LOO estimates may be unreliable.",
      call. = FALSE
    )
  }

  has_stan_ll <- fit_has_variable(x$fit, "log_lik")
  if (is.null(use_stan_log_lik)) {
    use_stan <- isTRUE(has_stan_ll)
  } else {
    use_stan <- isTRUE(use_stan_log_lik)
    if (use_stan && !has_stan_ll) {
      stop("use_stan_log_lik = TRUE requested, but fitted object does not contain Stan 'log_lik'")
    }
  }

  if (use_stan) {
    ll_arr <- tryCatch(x$fit$draws("log_lik", format = "draws_array"), error = function(e) NULL)
    if (is.null(ll_arr) || length(dim(ll_arr)) != 3L) {
      if (is.null(r_eff)) {
        stop(
          "Stan 'log_lik' draws could not be extracted in iteration x chain x observation layout. ",
          "Cannot compute chain-aware relative_eff(); supply explicit r_eff (or r_eff = NA deliberately)."
        )
      }
      ll_mat <- extract_ordered_draws(x$fit, "log_lik", N)
      loo_obj <- loo::loo(ll_mat, r_eff = r_eff, save_psis = save_psis, cores = cores, ...)
      used_source <- "stan_log_lik_matrix"
    } else {
      vn <- dimnames(ll_arr)[[3L]]
      idx <- as.integer(sub(".*\\[([0-9]+)\\]", "\\1", vn))
      if (length(idx) != N || anyNA(idx)) stop("Stan log_lik dimensions do not match N = ", N)
      ll_arr <- ll_arr[, , order(idx), drop = FALSE]
      if (is.null(r_eff)) {
        if (dim(ll_arr)[2L] < 1L || dim(ll_arr)[1L] < 4L) {
          stop("Insufficient iteration/chain dimensions to compute relative_eff(); supply explicit r_eff or r_eff = NA")
        }
        r_eff <- unname(as.numeric(loo::relative_eff(exp(ll_arr), cores = cores)))
      }
      loo_obj <- loo::loo(ll_arr, r_eff = r_eff, save_psis = save_psis, cores = cores, ...)
      used_source <- "stan_log_lik"
    }
  } else {
    if (!can_evaluate_in_r(x)) {
      stop(
        "Cannot compute log-likelihood in R for PSIS-LOO: ", get_eval_contract(x)$reason,
        if (has_stan_ll) " Set use_stan_log_lik = TRUE to use saved Stan log_lik." else ""
      )
    }
    params <- extract_bmfrm_param_draws(x)
    ids <- map_eval_indices(obs_data, spec, params$labels)
    resp <- as.integer(obs_data[[spec$response]])
    ll_data <- as.data.frame(ids, stringsAsFactors = FALSE)
    # Facet names cannot contain dots, so this cannot overwrite a facet index.
    ll_data$.bmfrm_response <- resp

    if (is.null(r_eff)) {
      cl <- extract_verified_chain_layout(x$fit, params$D)
      if (!cl$valid || cl$n_iter < 4L) {
        stop(
          "Verified iteration-by-chain layout is unavailable; refusing to silently assume r_eff = 1. ",
          "Ensure chain metadata is present or pass explicit r_eff (or r_eff = NA deliberately)."
        )
      }
      # Stream r_eff in chunks of size chunk_size: O(D * chunk_size) memory, zero D x N retention
      cs <- default_bmfrm_chunk_size(params$D, N, chunk_size)
      r_eff <- numeric(N)
      starts <- seq.int(1L, N, by = cs)
      for (s in starts) {
        e <- min(N, s + cs - 1L)
        rows <- s:e
        ll_ch <- eval_bmfrm_chunk(
          params, ids, rows, resp_chunk = resp[rows],
          compute_mu = FALSE, compute_log_lik = TRUE
        )$log_lik
        lik_arr <- array(exp(ll_ch), dim = c(cl$n_iter, cl$n_chains, length(rows)))
        r_eff[rows] <- as.numeric(loo::relative_eff(lik_arr, cores = cores))
      }
    }

    # True O(D) per-observation streaming PSIS-LOO via loo::loo.function()
    loo_obj <- loo::loo(
      eval_bmfrm_log_lik_obs,
      data = ll_data,
      draws = params,
      r_eff = r_eff,
      save_psis = save_psis,
      cores = cores,
      ...
    )
    used_source <- "r_streaming"
  }

  pk_tab <- loo::pareto_k_table(loo_obj)
  bad_k <- loo::pareto_k_ids(loo_obj)
  mcse_val <- tryCatch(loo::mcse_loo(loo_obj), error = function(e) NA_real_)
  n_eff_min <- if (!is.null(loo_obj$diagnostics$n_eff)) min(loo_obj$diagnostics$n_eff, na.rm = TRUE) else NA_real_

  rec <- if (length(bad_k) > 0L) {
    paste0(length(bad_k), " observation(s) exceed loo's Pareto-k diagnostic threshold; ",
           "verify with exact refit cross-validation via cv_bmfrm(fit, target = 'rating').")
  } else {
    "All observations satisfy loo's Pareto-k diagnostic threshold."
  }

  loo_obj$bmfrm_meta <- list(
    source = used_source,
    prediction_target = "rating",
    actual_target = actual_target,
    singleton_obs = sing$singleton_obs,
    singleton_facets = sing$singleton_facets,
    pareto_k_table = pk_tab,
    bad_k_ids = bad_k,
    mcse_elpd_loo = mcse_val,
    psis_n_eff_min = n_eff_min,
    mcmc_convergence = conv,
    recommendation = rec
  )
  loo_obj$comparison_key <- build_comparison_key(
    method = "psis_loo",
    target = "rating",
    target_facet = NULL,
    actual_target = actual_target,
    known_facets = if (length(sing$singleton_obs)) setdiff(spec$facets_main, names(sing$singleton_facets)) else spec$facets_main,
    unseen_facets = if (length(sing$singleton_obs)) names(sing$singleton_facets) else character(0),
    data_clean = obs_data,
    response_col = spec$response,
    K = get_spec_K_max(spec),
    aggregation = "sum",
    data_info = x$data_info,
    facets_main = spec$facets_main,
    K_by_level = x$data_info$K_by_level %||% spec$K_by_level,
    threshold_facet = spec$threshold_facet
  )
  class(loo_obj) <- unique(c("bmfrm_loo", class(loo_obj)))
  loo_obj
}

#' Print method for bmfrm_loo objects
#' @param x A `bmfrm_loo` object.
#' @param ... Passed to next `print` method.
#' @export
print.bmfrm_loo <- function(x, ...) {
  meta <- x$bmfrm_meta
  cat("BayesMFRM PSIS-LOO\n")
  cat("  Prediction target :", meta$actual_target, "\n")
  cat("  Log-lik source    :", meta$source, "\n")
  if (is.finite(meta$mcse_elpd_loo)) {
    cat("  MCSE(elpd_loo)    :", sprintf("%.4f", meta$mcse_elpd_loo), "\n")
  }
  if (is.finite(meta$psis_n_eff_min)) {
    cat("  Min PSIS ESS      :", sprintf("%.1f", meta$psis_n_eff_min), "\n")
  }
  if (!is.null(meta$mcmc_convergence)) {
    cat("  Original MCMC     :", if (isTRUE(meta$mcmc_convergence$converged)) "converged" else paste0("NOT converged (", meta$mcmc_convergence$status, ")"), "\n")
  }
  if (length(meta$bad_k_ids)) {
    cat("  Note              :", meta$recommendation, "\n")
  }
  cat("\n")
  NextMethod("print")
}

#' Build balanced rating-level or grouped folds
#' @keywords internal
build_cv_folds <- function(data_clean, facets_main, target, target_facet, K_folds, seed = NULL) {
  N <- nrow(data_clean)
  if (!is.numeric(K_folds) || length(K_folds) != 1L || !is.finite(K_folds) ||
      K_folds != floor(K_folds) || K_folds < 2L) {
    stop("K_folds must be an integer >= 2")
  }
  K_folds <- as.integer(K_folds)

  make_fn <- function() {
    if (identical(target, "rating")) {
      if (K_folds > N) stop("K_folds (", K_folds, ") exceeds number of observations (", N, ")")
      # Stratify round-robin by the facet with the smallest minimum level count to avoid dropping levels
      min_counts <- vapply(facets_main, function(f) min(table(as.character(data_clean[[f]]))), numeric(1))
      strat_f <- facets_main[which.min(min_counts)]
      vals <- as.character(data_clean[[strat_f]])
      folds <- integer(N)
      offset <- 0L
      for (lv in sample(unique(vals))) {
        idx <- which(vals == lv)
        idx <- idx[sample.int(length(idx))]
        folds[idx] <- ((seq_along(idx) - 1L + offset) %% K_folds) + 1L
        offset <- (offset + length(idx)) %% K_folds
      }
      return(folds)
    }
    g_vals <- as.character(data_clean[[target_facet]])
    u_groups <- unique(g_vals)
    G <- length(u_groups)
    if (K_folds > G) {
      stop("K_folds (", K_folds, ") exceeds number of unique levels in '", target_facet, "' (", G, ")")
    }
    perm_groups <- sample(u_groups)
    g_fold <- ((seq_len(G) - 1L) %% K_folds) + 1L
    g_fold[match(g_vals, perm_groups)]
  }
  if (!is.null(seed)) with_bmfrm_seed(seed, make_fn()) else make_fn()
}

#' Validate model hierarchy requirements for unseen-level CV
#' @keywords internal
validate_cv_hierarchy <- function(spec, target_facet) {
  if (isTRUE(spec$K_heterogeneous) && identical(spec$threshold_facet, target_facet)) {
    stop(
      "Grouped CV holding out threshold_facet '", target_facet,
      "' is not supported in v1 heterogeneous-category PCM because unseen levels have ",
      "undefined rubric category counts and step distributions."
    )
  }
  hier <- spec$hierarchical %||% character(0)
  missing_terms <- character(0)
  if (!(target_facet %in% hier)) {
    missing_terms <- c(missing_terms, target_facet)
  }
  for (term in spec$facets_bias) {
    comps <- strsplit(term, ":", fixed = TRUE)[[1L]]
    if (target_facet %in% comps && !(term %in% hier)) {
      missing_terms <- c(missing_terms, term)
    }
  }
  if (identical(spec$family, "pcm") && identical(spec$threshold_facet, target_facet) && !("tau" %in% hier)) {
    missing_terms <- c(missing_terms, "tau")
  }
  if (length(missing_terms)) {
    stop(
      "Grouped CV for unseen '", target_facet, "' requires hierarchical distributions for all model ",
      "components involving '", target_facet, "'. Missing from hierarchical: ",
      paste(missing_terms, collapse = ", "), "."
    )
  }
  invisible(TRUE)
}

#' Score held-out observations for one rating-level CV fold (known levels)
#' @keywords internal
score_rating_fold <- function(fold_fit, test_df, chunk_size = NULL) {
  spec <- fold_fit$spec
  params <- extract_bmfrm_param_draws(fold_fit)
  ids <- map_eval_indices(test_df, spec, params$labels)
  resp <- as.integer(test_df[[spec$response]])
  N_test <- nrow(test_df)
  D <- params$D
  cl <- extract_verified_chain_layout(fold_fit$fit, D)
  cs <- default_bmfrm_chunk_size(D, N_test, chunk_size)

  lpd_vec <- numeric(N_test)
  mcse_vec <- rep(NA_real_, N_test)
  starts <- seq.int(1L, N_test, by = cs)
  for (s in starts) {
    e <- min(N_test, s + cs - 1L)
    rows <- s:e
    ll_ch <- eval_bmfrm_chunk(
      params, ids, rows, resp_chunk = resp[rows],
      compute_mu = FALSE, compute_log_lik = TRUE
    )$log_lik
    for (loc in seq_along(rows)) {
      ll_i <- ll_ch[, loc]
      lpd_vec[rows[loc]] <- logmeanexp_vec(ll_i)
      if (cl$valid && cl$n_chains >= 2L && cl$n_iter >= 4L) {
        m_i <- max(ll_i)
        p_mat <- matrix(exp(ll_i - m_i), nrow = cl$n_iter, ncol = cl$n_chains)
        p_bar <- mean(p_mat)
        if (is.finite(p_bar) && p_bar > 0) {
          se_p <- tryCatch(posterior::mcse_mean(p_mat), error = function(e) NA_real_)
          if (is.finite(se_p) && se_p >= 0) mcse_vec[rows[loc]] <- se_p / p_bar
        }
      }
    }
  }
  list(lpd_pointwise = lpd_vec, mcse_pointwise = mcse_vec)
}

#' Score held-out groups for one grouped CV fold via exact double-logmeanexp marginalization
#'
#' Evaluates, for each held-out group `g` of `target_facet`:
#'   - Group joint lpd: `logmeanexp_d(logmeanexp_m(sum_{i in g} l_i(phi_d, u_{gdm})))`
#'   - Rating-level lpd: `logmeanexp_d(logmeanexp_m(l_i(phi_d, u_{gdm})))`
#' sharing the complete joint latent vector `u_{gdm}` (main effect, sum-to-zero
#' constrained interactions anchored to training-fold partner levels, and
#' hierarchical PCM steps when applicable) across all observations i in g.
#'
#' @keywords internal
score_grouped_fold <- function(fold_fit, test_df, target_facet, M = 100L,
                               seed = NULL, chunk_size = NULL) {
  spec <- fold_fit$spec
  params <- extract_bmfrm_param_draws(fold_fit)
  D <- params$D
  K <- params$K
  labels <- params$labels
  cl <- extract_verified_chain_layout(fold_fit$fit, D)

  # Pre-extract compact hyperprior vectors for target_facet (anchored to training fold!)
  if (identical(target_facet, "person")) {
    mu_main <- extract_scalar_draws(fold_fit$fit, "mu_theta")
    sd_main <- extract_scalar_draws(fold_fit$fit, "sigma_theta")
  } else {
    mu_main <- rep(0, D)
    sd_main <- extract_scalar_draws(fold_fit$fit, paste0("sigma_", target_facet))
  }

  # Interactions involving target_facet: anchor constrained axis strictly to training-fold partner levels
  target_biases <- list()
  known_biases <- list()
  for (term in names(params$bias_draws)) {
    b <- params$bias_draws[[term]]
    if (target_facet %in% c(b$f1, b$f2)) {
      partner <- if (identical(b$f1, target_facet)) b$f2 else b$f1
      n_partner <- length(labels[[partner]])
      b_name <- gsub(":", "", term)
      sd_b <- extract_scalar_draws(fold_fit$fit, paste0("sigma_", b_name))
      b_spec <- spec$gaussian_model$biases[[term]]
      b_mode <- if (!is.null(b_spec)) b_spec$mode else "gaussian_basis"
      b_scale <- if (!is.null(b_spec)) (b_spec$scale %||% "raw") else "raw"
      c2_b <- if (identical(b_mode, "shrinkage_horseshoe")) {
        extract_scalar_draws(fold_fit$fit, paste0("c2_", b_name))
      } else NULL
      target_biases[[term]] <- list(
        partner = partner,
        n_known = n_partner,
        info = list(mode = b_mode, scale = b_scale, sd = sd_b, c2 = c2_b)
      )
    } else {
      known_biases[[term]] <- b
    }
  }

  # PCM step thresholds if threshold_facet == target_facet
  is_target_pcm <- identical(params$tau_info$mode, "pcm") &&
    identical(params$tau_info$threshold_facet, target_facet)
  pcm_hyper <- NULL
  if (is_target_pcm) {
    tc <- extract_ordered_draws(fold_fit$fit, "tau_common", K - 1L)
    sd_tau <- extract_scalar_draws(fold_fit$fit, "sigma_tau")
    dev_scale <- spec$gaussian_model$step$dev_scale %||% "raw"
    B_dev <- if (K > 2L) rsm_prior_basis(K - 1L, dev_scale) else NULL
    pcm_hyper <- list(tau_common = tc, sd = sd_tau, B_dev = B_dev)
  }

  # Map known facets (all facets except target_facet) strictly against training-fold labels
  known_facets <- setdiff(spec$facets_main, target_facet)
  known_ids <- map_eval_indices(test_df, spec, labels, adjust_facets = target_facet)
  resp_all <- as.integer(test_df[[spec$response]])
  g_vals <- as.character(test_df[[target_facet]])
  u_groups <- unique(g_vals)
  G <- length(u_groups)

  lpd_pointwise <- numeric(nrow(test_df))
  mcse_pointwise <- rep(NA_real_, nrow(test_df))
  lpd_group <- numeric(G)
  mcse_group <- rep(NA_real_, G)
  mc_ess_m_med <- numeric(G)
  mc_ess_m_min <- numeric(G)
  n_obs_group <- integer(G)

  run_scoring <- function() {
    group_seeds <- sample.int(.Machine$integer.max, G, replace = TRUE)

    for (gi in seq_len(G)) {
      g_id <- u_groups[gi]
      g_rows <- which(g_vals == g_id)
      n_g <- length(g_rows)
      n_obs_group[gi] <<- n_g
      y_g <- resp_all[g_rows]

      cs_g <- default_bmfrm_chunk_size(D, n_g, chunk_size)

      # Fixed cum_tau indices when tau does not depend on target_facet
      fixed_t_idx <- if (!is_target_pcm && (params$tau_info$mode %in% c("pcm", "pcm_ragged"))) {
        known_ids[[params$tau_info$threshold_facet]][g_rows]
      } else NULL

      # Accumulate exact group joint log-lik (D x M) while streaming rating chunks in O(D * cs_g)
      ll_joint_dm <- matrix(0, nrow = D, ncol = M)
      m_seeds <- with_bmfrm_seed(group_seeds[gi], sample.int(.Machine$integer.max, M, replace = TRUE))
      g_starts <- seq.int(1L, n_g, by = cs_g)

      p_idx_by_term <- lapply(target_biases, function(tb) known_ids[[tb$partner]][g_rows])
      scope_cols_by_term <- lapply(p_idx_by_term, function(v) sort(unique(v)))
      n_active_terms <- sum(vapply(target_biases, function(tb) tb$n_known > 1L, logical(1)))
      g_cache <- if (length(g_starts) > 1L && n_active_terms > 0L) {
        new_interaction_offset_cache(D, cs_g, n_concurrent_offsets = as.numeric(M) * n_active_terms)
      } else NULL

      for (gs in g_starts) {
        ge <- min(n_g, gs + cs_g - 1L)
        loc <- gs:ge
        rows_c <- g_rows[loc]
        nc <- length(loc)
        need_after <- (ge < n_g)

        # Known-facet base linear predictor strictly for this observation chunk (D x nc)
        eta_known_c <- matrix(0, nrow = D, ncol = nc)
        if (!identical(target_facet, "person")) {
          eta_known_c <- eta_known_c + params$theta[, known_ids$person[rows_c], drop = FALSE]
        }
        for (f in intersect(names(params$facet_draws), known_facets)) {
          eta_known_c <- eta_known_c - params$facet_draws[[f]][, known_ids[[f]][rows_c], drop = FALSE]
        }
        for (b in known_biases) {
          cidx <- known_ids[[b$f1]][rows_c] + (known_ids[[b$f2]][rows_c] - 1L) * b$n1
          eta_known_c <- eta_known_c + b$draws[, cidx, drop = FALSE]
        }

        max_c <- matrix(-Inf, nrow = D, ncol = nc)
        sum_exp_c <- matrix(0, nrow = D, ncol = nc)

        for (m in seq_len(M)) {
          u_draw <- with_bmfrm_seed(m_seeds[m], {
            u_main <- mu_main + sd_main * stats::rnorm(D)
            sub_seeds <- stats::setNames(integer(length(target_biases)), names(target_biases))
            for (term in names(target_biases)) {
              if (target_biases[[term]]$n_known > 1L) {
                sub_seeds[[term]] <- sample.int(.Machine$integer.max, 1L)
              }
            }
            cum_tau_m <- NULL
            if (is_target_pcm) {
              if (K > 2L) {
                z_dev <- matrix(stats::rnorm(D * (K - 2L)), nrow = D, ncol = K - 2L)
                tau_m <- pcm_hyper$tau_common + (z_dev %*% t(pcm_hyper$B_dev)) * pcm_hyper$sd
              } else {
                tau_m <- pcm_hyper$tau_common
              }
              cum_tau_m <- matrix(0, nrow = D, ncol = K)
              if (K > 1L) {
                for (k in 2:K) cum_tau_m[, k] <- cum_tau_m[, k - 1L] + tau_m[, k - 1L]
              }
            }
            list(u_main = u_main, sub_seeds = sub_seeds, cum_tau = cum_tau_m)
          })

          eta_c <- eta_known_c
          if (identical(target_facet, "person")) {
            eta_c <- eta_c + u_draw$u_main
          } else {
            eta_c <- eta_c - u_draw$u_main
          }
          for (term in names(target_biases)) {
            tb <- target_biases[[term]]
            if (tb$n_known > 1L) {
              eta_c <- eta_c + eval_constrained_interaction_chunk(
                D = D,
                n_known = tb$n_known,
                b_info = tb$info,
                slice_seed = u_draw$sub_seeds[[term]],
                chunk_k_idx = p_idx_by_term[[term]][loc],
                scope_cols = scope_cols_by_term[[term]],
                need_after_chunk = need_after,
                cache_key = paste0(term, ":", m),
                offset_cache = g_cache
              )
            }
          }

          if (!is_target_pcm && identical(params$tau_info$mode, "pcm_ragged")) {
            ll_c <- eval_softmax_ragged(
              eta = eta_c, tau_info = params$tau_info, t_idx = fixed_t_idx[loc],
              resp_chunk = y_g[loc], compute_mu = FALSE, compute_probs = FALSE,
              compute_kurtosis = FALSE, compute_log_lik = TRUE, compute_y_rep = FALSE
            )$log_lik
          } else {
            get_ct <- function(k) {
              if (is_target_pcm) {
                matrix(u_draw$cum_tau[, k], nrow = D, ncol = nc)
              } else if (identical(params$tau_info$mode, "shared")) {
                matrix(params$tau_info$cum_tau[, k], nrow = D, ncol = nc)
              } else {
                params$tau_info$cum_tau_list[[k]][, fixed_t_idx[loc], drop = FALSE]
              }
            }

            max_l <- matrix(0, nrow = D, ncol = nc)
            if (K >= 2L) {
              for (k in 2:K) max_l <- pmax(max_l, (k - 1) * eta_c - get_ct(k))
            }
            denom <- exp(-max_l)
            if (K >= 2L) {
              for (k in 2:K) denom <- denom + exp((k - 1) * eta_c - get_ct(k) - max_l)
            }
            log_denom <- log(denom) + max_l
            ll_c <- matrix(NA_real_, nrow = D, ncol = nc)
            y_loc <- y_g[loc]
            for (k in seq_len(K)) {
              cols <- which(y_loc == k)
              if (length(cols)) {
                lk <- if (k == 1L) matrix(0, nrow = D, ncol = length(cols)) else {
                  ct_cols <- if (is_target_pcm) {
                    matrix(u_draw$cum_tau[, k], nrow = D, ncol = length(cols))
                  } else if (identical(params$tau_info$mode, "shared")) {
                    matrix(params$tau_info$cum_tau[, k], nrow = D, ncol = length(cols))
                  } else {
                    params$tau_info$cum_tau_list[[k]][, fixed_t_idx[loc[cols]], drop = FALSE]
                  }
                  (k - 1) * eta_c[, cols, drop = FALSE] - ct_cols
                }
                ll_c[, cols] <- lk - log_denom[, cols, drop = FALSE]
              }
            }
          }

          ll_joint_dm[, m] <- ll_joint_dm[, m] + rowSums(ll_c)
          prev_max <- max_c
          new_max <- pmax(prev_max, ll_c)
          prev_scale <- ifelse(is.finite(prev_max), exp(prev_max - new_max), 0)
          sum_exp_c <- sum_exp_c * prev_scale + exp(ll_c - new_max)
          max_c <- new_max
        }

        # Exact rating-level marginal lpd summarized immediately for this chunk (D x nc)
        lpd_c_mat <- max_c + log(sum_exp_c / M)
        for (loc_i in seq_len(nc)) {
          v_d <- lpd_c_mat[, loc_i]
          lpd_pointwise[rows_c[loc_i]] <<- logmeanexp_vec(v_d)
          if (cl$valid && cl$n_chains >= 2L && cl$n_iter >= 4L) {
            mi <- max(v_d)
            p_mat <- matrix(exp(v_d - mi), nrow = cl$n_iter, ncol = cl$n_chains)
            p_bar <- mean(p_mat)
            if (is.finite(p_bar) && p_bar > 0) {
              se_p <- tryCatch(posterior::mcse_mean(p_mat), error = function(e) NA_real_)
              if (is.finite(se_p) && se_p >= 0) mcse_pointwise[rows_c[loc_i]] <<- se_p / p_bar
            }
          }
        }
      }

      # Exact group joint lpd: logmeanexp_d( logmeanexp_m( sum_{i in g} l_i ) )
      lpd_g_d <- logmeanexp_rows(ll_joint_dm)
      lpd_group[gi] <<- logmeanexp_vec(lpd_g_d)

      if (M > 1L) {
        m_row <- apply(ll_joint_dm, 1L, max)
        w_dm <- exp(ll_joint_dm - m_row)
        w_dm <- w_dm / rowSums(w_dm)
        ess_d <- 1 / rowSums(w_dm * w_dm)
        mc_ess_m_med[gi] <<- stats::median(ess_d)
        mc_ess_m_min[gi] <<- min(ess_d)
      } else {
        mc_ess_m_med[gi] <<- 1
        mc_ess_m_min[gi] <<- 1
      }

      if (cl$valid && cl$n_chains >= 2L && cl$n_iter >= 4L) {
        mg <- max(lpd_g_d)
        p_mat <- matrix(exp(lpd_g_d - mg), nrow = cl$n_iter, ncol = cl$n_chains)
        p_bar <- mean(p_mat)
        if (is.finite(p_bar) && p_bar > 0) {
          se_p <- tryCatch(posterior::mcse_mean(p_mat), error = function(e) NA_real_)
          if (is.finite(se_p) && se_p >= 0) mcse_group[gi] <<- se_p / p_bar
        }
      }
    }
  }

  if (!is.null(seed)) with_bmfrm_seed(seed, run_scoring()) else run_scoring()

  list(
    lpd_pointwise = lpd_pointwise,
    mcse_pointwise = mcse_pointwise,
    group_scores = tibble::tibble(
      group = u_groups,
      n_obs = n_obs_group,
      elpd = lpd_group,
      mcse = mcse_group,
      mc_ess_m_median = mc_ess_m_med,
      mc_ess_m_min = mc_ess_m_min
    )
  )
}

#' Offset a non-negative 32-bit RNG seed without integer overflow
#' @keywords internal
offset_bmfrm_seed <- function(base_seed, step = 0L) {
  if (is.null(base_seed)) return(NULL)
  s <- as.numeric(base_seed)
  modulus <- as.numeric(.Machine$integer.max) + 1
  if (length(s) != 1L || !is.finite(s) || s != floor(s) || s < 0 || s > .Machine$integer.max) {
    stop("seed must be a non-negative 32-bit integer (0..", .Machine$integer.max, ")")
  }
  st <- as.numeric(step)
  if (length(st) != 1L || !is.finite(st) || st != floor(st)) {
    stop("step must be a finite integer")
  }
  as.integer((s + st) %% modulus)
}

#' Build per-group known-partner scope and last-row liveness metadata
#' @keywords internal
build_group_bias_scope <- function(plan, g_rows) {
  if (length(plan$bias_new) == 0L) return(NULL)
  out <- list()
  for (term in names(plan$bias_new)) {
    bn <- plan$bias_new[[term]]
    if (identical(bn$mode, "zero")) next
    i1_g <- plan$ids[[bn$f1]][g_rows]
    i2_g <- plan$ids[[bn$f2]][g_rows]
    n2_new <- length(bn$col_seeds)
    c_scope <- vector("list", n2_new)
    c_last <- integer(n2_new)
    if (n2_new > 0L && bn$n1_known > 1L) {
      kn <- which(i1_g <= bn$n1_known & i2_g > bn$n2_known)
      if (length(kn)) {
        u2 <- i2_g[kn] - bn$n2_known
        sp1 <- split.default(i1_g[kn], u2)
        for (nm in names(sp1)) c_scope[[as.integer(nm)]] <- sort(unique(sp1[[nm]]))
        rk <- g_rows[kn]
        for (k_i in seq_along(kn)) {
          u_i <- u2[k_i]
          if (rk[k_i] > c_last[u_i]) c_last[u_i] <- rk[k_i]
        }
      }
    }
    n1_new <- length(bn$row_seeds)
    r_scope <- vector("list", n1_new)
    r_last <- integer(n1_new)
    if (n1_new > 0L && bn$n2_known > 1L) {
      nk <- which(i1_g > bn$n1_known & i2_g <= bn$n2_known)
      if (length(nk)) {
        u1 <- i1_g[nk] - bn$n1_known
        sp2 <- split.default(i2_g[nk], u1)
        for (nm in names(sp2)) r_scope[[as.integer(nm)]] <- sort(unique(sp2[[nm]]))
        rk <- g_rows[nk]
        for (k_i in seq_along(nk)) {
          u_i <- u1[k_i]
          if (rk[k_i] > r_last[u_i]) r_last[u_i] <- rk[k_i]
        }
      }
    }
    out[[term]] <- list(
      col_scope = c_scope, col_last_row = c_last,
      row_scope = r_scope, row_last_row = r_last
    )
  }
  out
}

#' Score a CV fold with arbitrary hierarchical unseen facets (strict_target = FALSE)
#' @keywords internal
score_general_unseen_fold <- function(fold_fit, test_df, target_facet = NULL,
                                      M = 100L, seed = NULL, chunk_size = NULL) {
  run_all <- function() {
    spec <- fold_fit$spec
    base_plan <- prepare_eval_plan(
      fold_fit, test_df,
      allow_new_levels = TRUE,
      sample_new_levels = "hierarchical"
    )
    base_plan$offset_cache <- NULL
    params <- base_plan$params
    D <- params$D
    cl <- extract_verified_chain_layout(fold_fit$fit, D)
    resp_all <- as.integer(test_df[[spec$response]])
    N_test <- nrow(test_df)

    lpd_pointwise <- numeric(N_test)
    mcse_pointwise <- rep(NA_real_, N_test)

    summarize_chunk_lpd <- function(lpd_c_mat, rows_c) {
      for (loc_i in seq_along(rows_c)) {
        v_d <- lpd_c_mat[, loc_i]
        lpd_pointwise[rows_c[loc_i]] <<- logmeanexp_vec(v_d)
        if (cl$valid && cl$n_chains >= 2L && cl$n_iter >= 4L) {
          mi <- max(v_d)
          p_mat <- matrix(exp(v_d - mi), nrow = cl$n_iter, ncol = cl$n_chains)
          p_bar <- mean(p_mat)
          if (is.finite(p_bar) && p_bar > 0) {
            se_p <- tryCatch(posterior::mcse_mean(p_mat), error = function(e) NA_real_)
            if (is.finite(se_p) && se_p >= 0) mcse_pointwise[rows_c[loc_i]] <<- se_p / p_bar
          }
        }
      }
    }

    if (is.null(target_facet)) {
      is_unseen_row <- rep(FALSE, N_test)
      for (f in names(base_plan$main_new)) {
        is_unseen_row <- is_unseen_row | (base_plan$ids[[f]] > base_plan$main_new[[f]]$n_known)
      }
      known_rows <- which(!is_unseen_row)
      unseen_rows <- which(is_unseen_row)

      if (length(known_rows) > 0L) {
        sc_k <- score_rating_fold(fold_fit, test_df[known_rows, , drop = FALSE], chunk_size = chunk_size)
        lpd_pointwise[known_rows] <- sc_k$lpd_pointwise
        mcse_pointwise[known_rows] <- sc_k$mcse_pointwise
      }

      if (length(unseen_rows) > 0L) {
        n_u <- length(unseen_rows)
        cs_u <- default_bmfrm_chunk_size(D, n_u, chunk_size)
        m_seeds <- sample.int(.Machine$integer.max, M, replace = M > .Machine$integer.max)
        u_starts <- seq.int(1L, n_u, by = cs_u)
        n_conc_u <- if (length(u_starts) > 1L && length(base_plan$bias_new) > 0L) {
          count_concurrent_interaction_slices(base_plan$ids, base_plan$bias_new, unseen_rows, cs_u)
        } else 0L
        r_cache <- if (n_conc_u > 0L) {
          new_interaction_offset_cache(D, cs_u, n_concurrent_offsets = as.numeric(M) * n_conc_u)
        } else NULL
        for (us in u_starts) {
          ue <- min(n_u, us + cs_u - 1L)
          rows_c <- unseen_rows[us:ue]
          nc <- length(rows_c)
          max_c <- matrix(-Inf, nrow = D, ncol = nc)
          sum_exp_c <- matrix(0, nrow = D, ncol = nc)
          for (m in seq_len(M)) {
            ll_c <- eval_bmfrm_chunk_from_plan(
              base_plan, rows_c, resp_chunk = resp_all[rows_c],
              compute_mu = FALSE, compute_log_lik = TRUE,
              seed_salt = m_seeds[m], offset_cache = r_cache
            )$log_lik
            prev_max <- max_c
            new_max <- pmax(prev_max, ll_c)
            prev_scale <- ifelse(is.finite(prev_max), exp(prev_max - new_max), 0)
            sum_exp_c <- sum_exp_c * prev_scale + exp(ll_c - new_max)
            max_c <- new_max
          }
          summarize_chunk_lpd(max_c + log(sum_exp_c / M), rows_c)
        }
      }

      return(list(lpd_pointwise = lpd_pointwise, mcse_pointwise = mcse_pointwise))
    }

    g_vals <- as.character(test_df[[target_facet]])
    u_groups <- unique(g_vals)
    G <- length(u_groups)

    lpd_group <- numeric(G)
    mcse_group <- rep(NA_real_, G)
    mc_ess_m_med <- numeric(G)
    mc_ess_m_min <- numeric(G)
    n_obs_group <- integer(G)

    group_seeds <- sample.int(.Machine$integer.max, G, replace = TRUE)
    for (gi in seq_len(G)) {
      g_rows <- which(g_vals == u_groups[gi])
      n_g <- length(g_rows)
      n_obs_group[gi] <- n_g
      cs_g <- default_bmfrm_chunk_size(D, n_g, chunk_size)
      ll_joint_dm <- matrix(0, nrow = D, ncol = M)

      m_seeds <- with_bmfrm_seed(group_seeds[gi], sample.int(.Machine$integer.max, M, replace = M > .Machine$integer.max))
      g_starts <- seq.int(1L, n_g, by = cs_g)
      g_scope <- build_group_bias_scope(base_plan, g_rows)
      n_conc_g <- if (length(g_starts) > 1L && length(base_plan$bias_new) > 0L) {
        count_concurrent_interaction_slices(base_plan$ids, base_plan$bias_new, g_rows, cs_g)
      } else 0L
      g_cache <- if (n_conc_g > 0L) {
        new_interaction_offset_cache(D, cs_g, n_concurrent_offsets = as.numeric(M) * n_conc_g)
      } else NULL

      for (gs in g_starts) {
        ge <- min(n_g, gs + cs_g - 1L)
        rows_c <- g_rows[gs:ge]
        nc <- length(rows_c)
        max_c <- matrix(-Inf, nrow = D, ncol = nc)
        sum_exp_c <- matrix(0, nrow = D, ncol = nc)
        for (m in seq_len(M)) {
          ll_c <- eval_bmfrm_chunk_from_plan(
            base_plan, rows_c, resp_chunk = resp_all[rows_c],
            compute_mu = FALSE, compute_log_lik = TRUE,
            seed_salt = m_seeds[m], offset_cache = g_cache,
            bias_scope_override = g_scope
          )$log_lik
          ll_joint_dm[, m] <- ll_joint_dm[, m] + rowSums(ll_c)
          prev_max <- max_c
          new_max <- pmax(prev_max, ll_c)
          prev_scale <- ifelse(is.finite(prev_max), exp(prev_max - new_max), 0)
          sum_exp_c <- sum_exp_c * prev_scale + exp(ll_c - new_max)
          max_c <- new_max
        }
        summarize_chunk_lpd(max_c + log(sum_exp_c / M), rows_c)
      }

      lpd_g_d <- logmeanexp_rows(ll_joint_dm)
      lpd_group[gi] <- logmeanexp_vec(lpd_g_d)

      if (M > 1L) {
        m_row <- apply(ll_joint_dm, 1L, max)
        w_dm <- exp(ll_joint_dm - m_row)
        w_dm <- w_dm / rowSums(w_dm)
        ess_d <- 1 / rowSums(w_dm * w_dm)
        mc_ess_m_med[gi] <- stats::median(ess_d)
        mc_ess_m_min[gi] <- min(ess_d)
      } else {
        mc_ess_m_med[gi] <- 1
        mc_ess_m_min[gi] <- 1
      }

      if (cl$valid && cl$n_chains >= 2L && cl$n_iter >= 4L) {
        mg <- max(lpd_g_d)
        p_mat <- matrix(exp(lpd_g_d - mg), nrow = cl$n_iter, ncol = cl$n_chains)
        p_bar <- mean(p_mat)
        if (is.finite(p_bar) && p_bar > 0) {
          se_p <- tryCatch(posterior::mcse_mean(p_mat), error = function(e) NA_real_)
          if (is.finite(se_p) && se_p >= 0) mcse_group[gi] <- se_p / p_bar
        }
      }
    }

    list(
      lpd_pointwise = lpd_pointwise,
      mcse_pointwise = mcse_pointwise,
      group_scores = tibble::tibble(
        group = u_groups,
        n_obs = n_obs_group,
        elpd = lpd_group,
        mcse = mcse_group,
        mc_ess_m_median = mc_ess_m_med,
        mc_ess_m_min = mc_ess_m_min
      )
    )
  }

  if (!is.null(seed)) with_bmfrm_seed(seed, run_all()) else run_all()
}

#' Target-aware K-fold cross-validation for BayesMFRM models
#'
#' Performs K-fold cross-validation with explicit prediction-target
#' verification, training-fold factor level dropping, per-fold structural design
#' and MCMC convergence checks, and exact Monte Carlo marginalization over
#' unseen hierarchical groups when evaluating generalization to new examinees,
#' raters, or other facets.
#'
#' `target = "rating"` evaluates ratings from known facet levels. Grouped
#' targets hold out entire levels and integrate their effects over the learned
#' hierarchical population with `M` draws per posterior draw, using double
#' log-mean-exp. With `strict_target = TRUE`, collateral unseen facets cause an
#' error rather than silently changing the prediction target.
#'
#' For heterogeneous-category PCM, training folds preserve `K_by_level` even
#' when upper categories are unobserved. Holding out `threshold_facet` itself is
#' unsupported in v1. Passing `mc_ess_m_median >= 5` alone does not guarantee
#' integration accuracy: check larger `M` and independent `seed` values while
#' keeping both `folds` (e.g., `folds = cv_init$folds`) and `refit_seed` fixed
#' so the train/test split and cached fold posterior draws remain identical.
#'
#' @param object A `bmfrm_fit` object.
#' @param target Prediction target: `"rating"` (leave ratings out, known facet
#'   levels), `"person"` (leave examinees out), `"rater"` (leave raters out),
#'   or `"facet"` (leave levels of `facet` out).
#' @param facet Character name of the held-out facet when `target = "facet"`.
#' @param K_folds Integer number of folds (`>= 2`, default `5L`). Ignored if
#'   `folds` is supplied.
#' @param folds Optional integer vector of length `N` assigning each
#'   observation in `object$data_info$data_clean` to a fold `1..K_folds`.
#' @param M Positive integer number of Monte Carlo population draws per
#'   posterior draw for unseen-group marginalization (default `100L`; `M = 1L`
#'   is also valid). Recorded alongside integration ESS diagnostics; `D * M` is
#'   never treated as `D * M` independent posterior draws.
#' @param seed Optional integer RNG seed for fold construction, fold refits
#'   (unless `refit_seed` is specified), and Monte Carlo marginalization.
#' @param refit_seed Optional integer RNG seed used exclusively for fold MCMC
#'   refits. Fix both `folds = cv_init$folds` and `refit_seed` while varying
#'   `seed` and `M` to assess Monte Carlo integration sensitivity on identical
#'   cached fold posterior draws (varying `seed` with `folds = NULL` also
#'   changes the fold split).
#' @param strict_target Logical; if `TRUE` (default), stops whenever any fold
#'   causes an unintended facet to become unseen (e.g., holding out a rater
#'   also drops an examinee rated only by that rater, or a rating fold drops
#'   all ratings of some level).
#' @param allow_degenerate_interaction Logical; if `FALSE` (default), stops if
#'   any training fold reduces a facet or interaction axis to `<= 1` observed
#'   level (collapsing the interaction to 0 df).
#' @param chunk_size Optional observation chunk size for streaming evaluation.
#'   When an unseen level's observations span multiple chunks, constrained
#'   two-way interaction centering offsets (`D * 8` bytes each) are cached
#'   across chunks using a byte budget estimated from
#'   `M * n_concurrent_offsets * D * 8` up to a 64 MiB cap (releasing each
#'   offset immediately after its final chunk and retaining 0 bytes for
#'   single-chunk levels). If `M * n_concurrent_offsets * D * 8` exceeds the
#'   64 MiB cap, uncached offsets are deterministically recomputed per chunk as
#'   an intentional memory-vs-compute trade-off with identical numerical results.
#' @param refit_fn Optional custom refit function `function(train_df, fold, spec, fit_args)`
#'   returning a `bmfrm_fit` object (useful for custom clusters or fast unit
#'   tests). When `NULL`, refits via [bmfrm()] using the locked model
#'   specification (`K`, `priors`, `family`, `threshold_facet`, `hierarchical`,
#'   and sampler settings from `object$fit_args`).
#' @param ... Optional overrides passed to [bmfrm()] when refitting folds.
#' @return A `bmfrm_cv` object containing:
#'   \itemize{
#'     \item `target`, `target_facet`, `actual_target`: requested and verified
#'       prediction task (`"new_rating_known_levels"` or `"new_<facet>_known_other_facets"`).
#'     \item `elpd_group`, `se_group`: held-out group joint ELPD and group-level
#'       SE (primary for grouped CV; `NA` for rating CV).
#'     \item `elpd_pointwise`, `se_pointwise_clustered`, `se_pointwise`: marginal
#'       rating ELPD, group-clustered SE, and naive rating-level SE for reference.
#'     \item `pointwise`: tibble with preserved `obs_id`, cleaned-data `row_index`,
#'       `fold`, `group`, `elpd`, and chain-aware `mcse`.
#'     \item `group_scores`: joint scores with `group`, `n_obs`, `elpd`, `mcse`,
#'       `mc_ess_m_median`, and `mc_ess_m_min` (`NULL` for rating CV).
#'     \item `fold_diagnostics`: per-fold design connectivity/identifiability,
#'       rank, `max_rhat`, `min_ess_bulk`, `n_divergent`, ELPD, `min_mc_ess_m`, and `status`.
#'     \item `integration_diagnostics`, `comparison_key`: Monte Carlo metadata
#'       and the observation-aligned compatibility key used by [compare_bmfrm()].
#'   }
#' @export
cv_bmfrm <- function(object,
                     target = c("rating", "person", "rater", "facet"),
                     facet = NULL,
                     K_folds = 5L,
                     folds = NULL,
                     M = 100L,
                     seed = NULL,
                     refit_seed = NULL,
                     strict_target = TRUE,
                     allow_degenerate_interaction = FALSE,
                     chunk_size = NULL,
                     refit_fn = NULL,
                     ...) {
  if (!inherits(object, "bmfrm_fit")) stop("object must be a bmfrm_fit object")
  target <- match.arg(target)
  spec <- object$spec
  obs_data <- object$data_info$data_clean
  N <- nrow(obs_data)

  if (!can_evaluate_in_r(object)) {
    stop("cv_bmfrm() requires a verified R-engine model: ", get_eval_contract(object)$reason)
  }
  if (!is.numeric(M) || length(M) != 1L || !is.finite(M) || M != floor(M) || M < 1) {
    stop("M must be a positive integer")
  }
  M <- as.integer(M)
  seed <- offset_bmfrm_seed(seed, 0L)
  refit_seed <- offset_bmfrm_seed(refit_seed, 0L)

  overrides <- list(...)
  if (length(overrides) > 0L) {
    if (is.null(names(overrides)) || any(!nzchar(names(overrides)))) {
      stop("Additional refit arguments in '...' must be named")
    }
    locked_bmfrm <- intersect(
      names(overrides),
      c("formula", "data", "K", "priors", "stan_file", "family",
        "threshold_facet", "hierarchical", "save_csvs", "save_gq", "obs_id")
    )
    if (length(locked_bmfrm) > 0L) {
      stop("Cannot override locked model/data specification argument(s) in cv_bmfrm(): ",
           paste(locked_bmfrm, collapse = ", "))
    }
    reserved_dots <- intersect(
      names(overrides),
      c("seed", "parallel_chains", "output_dir", "iter_sampling", "iter_warmup", "save_warmup")
    )
    if (length(reserved_dots) > 0L) {
      stop("Do not pass reserved cmdstanr argument(s) (",
           paste(reserved_dots, collapse = ", "),
           ") in '...'; use iter, warmup, chains, cores, cache_dir, seed, and refit_seed instead")
    }
  }
  base_fit_args <- object$fit_args %||% list(
    formula = spec$formula,
    K = spec$K,
    priors = object$priors,
    stan_file = NULL,
    family = spec$family,
    threshold_facet = spec$threshold_facet,
    hierarchical = spec$hierarchical,
    save_gq = "minimal",
    cache_dir = spec$cache_dir %||% "stan_cache",
    cpp_options = list(),
    iter = 2000L,
    warmup = 1000L,
    chains = 4L,
    cores = getOption("mc.cores", 1L),
    seed = seed,
    sampler_args = list()
  )
  iter_cv <- overrides$iter %||% base_fit_args$iter %||% 2000L
  if (!is.numeric(iter_cv) || length(iter_cv) != 1L || !is.finite(iter_cv) ||
      iter_cv != floor(iter_cv) || iter_cv < 2L) {
    stop("iter must be a finite integer >= 2")
  }
  warmup_cv <- if (!is.null(overrides$warmup)) {
    overrides$warmup
  } else if (!is.null(overrides$iter)) {
    floor(iter_cv / 2)
  } else {
    base_fit_args$warmup %||% floor(iter_cv / 2)
  }
  if (!is.numeric(warmup_cv) || length(warmup_cv) != 1L || !is.finite(warmup_cv) ||
      warmup_cv != floor(warmup_cv) || warmup_cv < 0L || warmup_cv >= iter_cv) {
    stop("warmup must be a non-negative integer strictly less than iter")
  }
  chains_cv <- overrides$chains %||% base_fit_args$chains %||% 4L
  if (!is.numeric(chains_cv) || length(chains_cv) != 1L || !is.finite(chains_cv) ||
      chains_cv != floor(chains_cv) || chains_cv < 1L) {
    stop("chains must be a positive integer")
  }
  cores_cv <- overrides$cores %||% base_fit_args$cores %||% getOption("mc.cores", 1L)
  if (is.null(cores_cv) || (length(cores_cv) == 1L && is.na(cores_cv) && !("cores" %in% names(overrides)))) {
    cores_cv <- 1L
  }
  if (!is.numeric(cores_cv) || length(cores_cv) != 1L || !is.finite(cores_cv) ||
      cores_cv != floor(cores_cv) || cores_cv < 1L) {
    stop("cores must be a positive integer")
  }
  refit_mode_cv <- match.arg(overrides$refit %||% "on_change", c("on_change", "never", "always"))
  cache_dir_cv <- overrides$cache_dir %||% base_fit_args$cache_dir %||% spec$cache_dir %||% "stan_cache"
  if (!is.character(cache_dir_cv) || length(cache_dir_cv) != 1L || is.na(cache_dir_cv) || !nzchar(cache_dir_cv)) {
    stop("cache_dir must be a non-empty character path")
  }

  # 1. Resolve target facet and validate full model hierarchy for unseen levels
  target_facet <- switch(
    target,
    rating = NULL,
    person = "person",
    rater = {
      if (!("rater" %in% spec$facets_main)) {
        stop("target = 'rater' requested, but 'rater' is not in model facets (",
             paste(spec$facets_main, collapse = ", "), ")")
      }
      "rater"
    },
    facet = {
      if (is.null(facet) || !is.character(facet) || length(facet) != 1L || !(facet %in% spec$facets_main)) {
        stop("For target = 'facet', supply one valid model facet in 'facet' (available: ",
             paste(spec$facets_main, collapse = ", "), ")")
      }
      facet
    }
  )

  if (!is.null(target_facet)) {
    validate_cv_hierarchy(spec, target_facet)
  }

  # 2. Build or validate fold assignment aligned strictly to obs_id = 1..N
  if (is.null(folds)) {
    folds <- build_cv_folds(obs_data, spec$facets_main, target, target_facet, K_folds, seed = seed)
    K_folds <- max(folds)
  } else {
    if (!is.numeric(folds) || length(folds) != N || any(!is.finite(folds)) ||
        any(folds != floor(folds)) || any(folds < 1L)) {
      stop("folds must be an integer vector of length nrow(data_clean) (", N, ") with positive fold IDs")
    }
    folds <- as.integer(folds)
    u_f <- sort(unique(folds))
    if (length(u_f) < 2L || !identical(u_f, seq_len(length(u_f)))) {
      stop("folds must contain consecutive integers 1..K_folds with K_folds >= 2")
    }
    K_folds <- length(u_f)
    if (!is.null(target_facet)) {
      g_vals <- as.character(obs_data[[target_facet]])
      n_folds_per_group <- tapply(folds, g_vals, function(z) length(unique(z)))
      split_groups <- names(n_folds_per_group)[n_folds_per_group > 1L]
      if (length(split_groups)) {
        stop(
          "For grouped CV on '", target_facet, "', each level must belong to a single fold; ",
          "split level(s): ", paste(utils::head(split_groups, 5L), collapse = ", ")
        )
      }
    }
  }

  # 3. Pre-fit audit across ALL folds: dropleveled training prep, actual unseen facets,
  #    1-level degeneracy check, and per-fold structural design identifiability!
  fold_plans <- vector("list", K_folds)
  all_unseen_facets <- character(0)
  for (k in seq_len(K_folds)) {
    train_idx <- which(folds != k)
    test_idx <- which(folds == k)
    if (!length(train_idx) || !length(test_idx)) {
      stop("Fold ", k, " has empty training or test split")
    }
    train_df <- drop_unused_facet_levels(obs_data[train_idx, , drop = FALSE], spec$facets_main)
    test_df <- obs_data[test_idx, , drop = FALSE]

    # Actual unseen facets on this fold
    unseen_k <- list()
    for (f in spec$facets_main) {
      tr_lv <- levels(train_df[[f]])
      te_lv <- unique(as.character(test_df[[f]]))
      u_lv <- setdiff(te_lv, tr_lv)
      if (length(u_lv)) unseen_k[[f]] <- u_lv
    }
    all_unseen_facets <- union(all_unseen_facets, names(unseen_k))

    if (!is.null(target_facet)) {
      tr_t <- levels(train_df[[target_facet]])
      te_t <- unique(as.character(test_df[[target_facet]]))
      overlap_t <- intersect(tr_t, te_t)
      if (length(overlap_t)) {
        stop("Fold ", k, " has overlapping '", target_facet, "' levels between train and test: ",
             paste(utils::head(overlap_t, 5L), collapse = ", "))
      }
    }

    if (isTRUE(strict_target)) {
      if (is.null(target_facet)) {
        if (length(unseen_k)) {
          details <- paste(vapply(names(unseen_k), function(f) {
            paste0(f, " (", paste(utils::head(unseen_k[[f]], 3L), collapse = ", "), ")")
          }, character(1)), collapse = "; ")
          stop(
            "Fold ", k, " in target = 'rating' drops all training observations for facet(s): ",
            details, ". Supply custom folds keeping all facet levels in training or choose a grouped target."
          )
        }
      } else {
        collateral <- setdiff(names(unseen_k), target_facet)
        if (length(collateral)) {
          details <- paste(vapply(collateral, function(f) {
            paste0(f, " (", paste(utils::head(unseen_k[[f]], 3L), collapse = ", "), ")")
          }, character(1)), collapse = "; ")
          stop(
            "Fold ", k, " holding out '", target_facet, "' also makes other facet(s) unseen in training: ",
            details, ". This changes the prediction target; adjust the rating plan or fold split, ",
            "or handle joint unseen facets explicitly."
          )
        }
      }
    }

    # Mandatory pre-fit guard for v1 heterogeneous-category PCM: unseen threshold_facet
    # levels are NEVER supported regardless of strict_target (must stop before any refit).
    if (isTRUE(spec$K_heterogeneous) &&
        identical(spec$family, "pcm") &&
        !is.null(spec$threshold_facet) &&
        (spec$threshold_facet %in% names(unseen_k))) {
      u_tf <- unseen_k[[spec$threshold_facet]]
      stop(
        "Fold ", k, " leaves threshold_facet '", spec$threshold_facet, "' level(s) (",
        paste(utils::head(u_tf, 5L), collapse = ", "),
        ") unseen in training, which is not supported in v1 heterogeneous-category PCM ",
        "(even when strict_target = FALSE), because their rubric category count and step distribution are undefined."
      )
    }

    # Validate full model hierarchy for every unseen facet before running any refits
    for (uf in names(unseen_k)) {
      validate_cv_hierarchy(spec, uf)
    }

    # Check training facet sizes (1-level degeneracy for interactions)
    if (!isTRUE(allow_degenerate_interaction) && length(spec$facets_bias) > 0L) {
      for (term in spec$facets_bias) {
        comps <- strsplit(term, ":", fixed = TRUE)[[1L]]
        degen <- comps[vapply(comps, function(f) nlevels(train_df[[f]]) <= 1L, logical(1))]
        if (length(degen)) {
          stop(
            "Fold ", k, " reduces training facet(s) ", paste(degen, collapse = ", "),
            " to <= 1 observed level, collapsing interaction '", term,
            "' to 0 degrees of freedom. Use fewer folds or set allow_degenerate_interaction = TRUE deliberately."
          )
        }
      }
    }

    K_fold_dc <- if (isTRUE(spec$K_heterogeneous) && !is.null(spec$K_by_level)) {
      spec$K_by_level[levels(train_df[[spec$threshold_facet]])]
    } else {
      spec$K
    }
    dc_k <- check_design_bmfrm(
      spec$formula, train_df, K = K_fold_dc,
      family = spec$family, threshold_facet = spec$threshold_facet
    )
    fold_plans[[k]] <- list(
      fold = k,
      train_idx = train_idx,
      test_idx = test_idx,
      train_df = train_df,
      test_df = test_df,
      unseen_k = unseen_k,
      design = dc_k
    )
  }

  # 4. Locked refit specification (never re-infer K or alter priors/sampler settings)
  base_fit_args$iter <- as.integer(iter_cv)
  base_fit_args$warmup <- as.integer(warmup_cv)
  base_fit_args$chains <- as.integer(chains_cv)
  base_fit_args$cores <- as.integer(cores_cv)
  base_fit_args$cache_dir <- cache_dir_cv
  compiled_cv_model <- NULL
  default_refit_fn <- function(train_df, fold, spec, fit_args) {
    fold_spec <- spec
    if (!isTRUE(fold_spec$custom_stan)) {
      fold_spec$save_gq <- "minimal"
    }
    if (!is.null(spec$K_by_level) && !is.null(spec$threshold_facet)) {
      tr_tf <- levels(train_df[[spec$threshold_facet]])
      fold_spec$K_by_level <- spec$K_by_level[tr_tf]
      fold_spec$K_max <- as.integer(max(fold_spec$K_by_level))
      fold_spec$K_heterogeneous <- length(unique(fold_spec$K_by_level)) > 1L
      fold_spec$K <- if (fold_spec$K_heterogeneous) fold_spec$K_by_level else as.integer(fold_spec$K_by_level[1L])
    }
    fold_spec <- resolve_bmfrm_spec(fold_spec, object$priors)
    stan_code <- if (isTRUE(fold_spec$custom_stan)) {
      if (!is.null(fit_args$stan_file) && file.exists(fit_args$stan_file)) {
        paste(readLines(fit_args$stan_file, warn = FALSE), collapse = "\n")
      } else if (!is.null(object$stan_code) && is.character(object$stan_code) && nzchar(object$stan_code)) {
        object$stan_code
      } else {
        stop("Custom Stan code is unavailable for CV refit (missing both stan_file on disk and object$stan_code)")
      }
    } else if (!is.null(fit_args$stan_file)) {
      paste(readLines(fit_args$stan_file, warn = FALSE), collapse = "\n")
    } else {
      build_stan_code(fold_spec, object$priors, save_gq = "minimal")
    }
    fold_spec$eval_contract <- resolve_eval_contract(fold_spec, stan_code)

    # Compile once and share the compiled CmdStan binary across all folds
    if (is.null(compiled_cv_model)) {
      compiled_cv_model <<- compile_bmfrm(
        stan_code = stan_code,
        model_name = overrides$model_name %||% fold_spec$model_name,
        cache_dir = cache_dir_cv,
        cpp_options = overrides$cpp_options %||% fit_args$cpp_options %||% list()
      )
    }

    data_info_k <- prepare_data_bmfrm(fold_spec, train_df, priors = object$priors)
    base_refit_seed <- if (!is.null(refit_seed)) refit_seed else seed
    fold_seed <- if (!is.null(base_refit_seed)) {
      offset_bmfrm_seed(base_refit_seed, fold * 1000L)
    } else if (!is.null(fit_args$seed)) {
      offset_bmfrm_seed(fit_args$seed, fold * 1000L)
    } else {
      NULL
    }
    iter_k <- as.integer(iter_cv)
    warmup_k <- as.integer(warmup_cv)
    chains_k <- as.integer(chains_cv)
    cores_k <- as.integer(cores_cv)
    refit_mode <- refit_mode_cv
    fold_cache_dir <- file.path(cache_dir_cv, "cv_folds")
    base_fold_rds <- file.path(
      fold_cache_dir,
      paste0(overrides$model_name %||% fold_spec$model_name, "_fold", fold, "_", compiled_cv_model$hash, ".rds")
    )

    sampler_args <- fit_args$sampler_args %||% list()
    reserved_ov <- c(
      "cpp_options", "iter", "warmup", "chains", "cores", "refit",
      "cache_dir", "design_check", "model_name"
    )
    for (nm in setdiff(names(overrides), reserved_ov)) sampler_args[[nm]] <- overrides[[nm]]

    extra_sampling_k <- c(
      list(
        iter_sampling = iter_k - warmup_k,
        iter_warmup = warmup_k,
        save_warmup = FALSE
      ),
      sampler_args
    )
    fold_rds <- resolve_fit_cache_file(
      base_cache_file = base_fold_rds,
      code = compiled_cv_model$model$code(),
      data = data_info_k$stan_data,
      seed = fold_seed,
      chains = chains_k,
      parallel_chains = min(chains_k, cores_k),
      extra_sampling = extra_sampling_k
    )

    fit_k <- do.call(
      fit_cmdstan_cached,
      c(
        list(
          model = compiled_cv_model$model,
          data = data_info_k$stan_data,
          file = fold_rds,
          refit = refit_mode,
          seed = fold_seed,
          chains = chains_k,
          parallel_chains = min(chains_k, cores_k)
        ),
        extra_sampling_k
      )
    )
    out_k <- postprocess_bmfrm(
      fit = fit_k,
      spec = fold_spec,
      data_info = data_info_k,
      priors = object$priors,
      stan_code = stan_code,
      stan_code_hash = compiled_cv_model$hash
    )
    out_k$design <- fold_plans[[fold]]$design
    out_k
  }
  runner <- if (!is.null(refit_fn)) refit_fn else default_refit_fn

  # 5. Fit & score every fold (never silently dropping failed folds!)
  lpd_pointwise <- rep(NA_real_, N)
  mcse_pointwise <- rep(NA_real_, N)
  group_tables <- list()
  fold_diag_rows <- vector("list", K_folds)
  n_draws_per_fold <- rep(NA_integer_, K_folds)

  for (k in seq_len(K_folds)) {
    fp <- fold_plans[[k]]
    dc_k <- fp$design
    fit_err <- NULL
    fold_fit <- tryCatch(
      runner(fp$train_df, fold = k, spec = spec, fit_args = base_fit_args),
      error = function(e) { fit_err <<- conditionMessage(e); NULL }
    )

    if (is.null(fold_fit)) {
      fold_diag_rows[[k]] <- tibble::tibble(
        fold = k, n_train = length(fp$train_idx), n_test = length(fp$test_idx),
        n_groups = if (is.null(target_facet)) NA_integer_ else length(unique(as.character(fp$test_df[[target_facet]]))),
        design_connected = isTRUE(dc_k$connected),
        design_identifiable = isTRUE(dc_k$identifiable),
        design_rank = dc_k$adjacent_logit_rank,
        design_params = dc_k$adjacent_logit_parameters,
        converged = FALSE, max_rhat = NA_real_, min_ess_bulk = NA_real_,
        n_divergent = NA_integer_, elpd_pointwise = NA_real_, elpd_group = NA_real_,
        min_mc_ess_m = NA_real_, status = paste0("refit_failed: ", fit_err)
      )
      next
    }

    # Verify unused factor levels did NOT leak into the training fit!
    for (f in spec$facets_main) {
      if (!setequal(fold_fit$data_info$facet_labels[[f]], levels(fp$train_df[[f]]))) {
        stop(
          "Fold ", k, " refit retained facet labels not observed in training for '", f,
          "'. Training fold preprocessing must drop unused factor levels."
        )
      }
    }
    # Verify category definitions match the original fit per observation / per level,
    # regardless of whether the original fit or fold refit uses uniform or heterogeneous K.
    resolve_fold_K_obs <- function(f_spec, f_di, df_check) {
      tf_f <- f_spec$threshold_facet
      k_map_f <- f_di$K_by_level %||% f_spec$K_by_level
      if (identical(f_spec$family, "pcm") && !is.null(tf_f) && !is.null(k_map_f)) {
        out_k <- as.integer(unname(k_map_f[as.character(df_check[[tf_f]])]))
        if (anyNA(out_k)) {
          stop("Fold ", k, " refit is missing per-level category count(s) for '", tf_f, "'")
        }
        return(out_k)
      }
      rep(as.integer(get_spec_K_max(f_spec)), nrow(df_check))
    }
    exp_K_tr <- resolve_fold_K_obs(spec, object$data_info, fp$train_df)
    got_K_tr <- resolve_fold_K_obs(fold_fit$spec, fold_fit$data_info, fp$train_df)
    if (!identical(got_K_tr, exp_K_tr)) {
      stop("Fold ", k, " refit changed per-level/per-observation category counts K relative to the original fit")
    }
    if (isTRUE(fold_fit$spec$K_heterogeneous) && !isTRUE(spec$K_heterogeneous)) {
      stop("Fold ", k, " refit changed category definition from uniform K to heterogeneous K_by_level")
    }
    if (!isTRUE(fold_fit$spec$K_heterogeneous) && isTRUE(spec$K_heterogeneous)) {
      tr_tf <- levels(fp$train_df[[spec$threshold_facet]])
      exp_tr_u <- unique(as.integer(unname((object$data_info$K_by_level %||% spec$K_by_level)[tr_tf])))
      if (length(exp_tr_u) > 1L) {
        stop("Fold ", k, " refit changed category definition from heterogeneous K_by_level to uniform K")
      }
    }
    if (!isTRUE(spec$K_heterogeneous) &&
        !identical(as.integer(get_spec_K_max(fold_fit$spec)), as.integer(get_spec_K_max(spec)))) {
      stop("Fold ", k, " refit changed K from ", get_spec_K_max(spec), " to ", get_spec_K_max(fold_fit$spec))
    }

    n_draws_per_fold[k] <- nrow(unclass(fold_fit$fit$draws("theta", format = "matrix")))
    conv_k <- fold_fit$convergence
    conv_ok <- if (is.null(conv_k) || is.null(conv_k$converged)) NA else isTRUE(conv_k$converged)

    fold_score_seed <- offset_bmfrm_seed(seed, k * 10007L)
    score_err <- NULL
    sc <- tryCatch({
      if (is.null(target_facet)) {
        if (length(fp$unseen_k) == 0L) {
          score_rating_fold(fold_fit, fp$test_df, chunk_size = chunk_size)
        } else {
          score_general_unseen_fold(
            fold_fit, fp$test_df, target_facet = NULL,
            M = M, seed = fold_score_seed, chunk_size = chunk_size
          )
        }
      } else {
        if (length(setdiff(names(fp$unseen_k), target_facet)) == 0L) {
          score_grouped_fold(
            fold_fit, fp$test_df, target_facet = target_facet,
            M = M, seed = fold_score_seed, chunk_size = chunk_size
          )
        } else {
          score_general_unseen_fold(
            fold_fit, fp$test_df, target_facet = target_facet,
            M = M, seed = fold_score_seed, chunk_size = chunk_size
          )
        }
      }
    }, error = function(e) { score_err <<- conditionMessage(e); NULL })

    if (is.null(sc)) {
      fold_diag_rows[[k]] <- tibble::tibble(
        fold = k, n_train = length(fp$train_idx), n_test = length(fp$test_idx),
        n_groups = if (is.null(target_facet)) NA_integer_ else length(unique(as.character(fp$test_df[[target_facet]]))),
        design_connected = isTRUE(dc_k$connected),
        design_identifiable = isTRUE(dc_k$identifiable),
        design_rank = dc_k$adjacent_logit_rank,
        design_params = dc_k$adjacent_logit_parameters,
        converged = conv_ok,
        max_rhat = if (!is.null(conv_k$max_rhat)) conv_k$max_rhat else NA_real_,
        min_ess_bulk = if (!is.null(conv_k$min_ess_bulk)) conv_k$min_ess_bulk else NA_real_,
        n_divergent = if (!is.null(conv_k$n_divergent)) as.integer(conv_k$n_divergent) else NA_integer_,
        elpd_pointwise = NA_real_, elpd_group = NA_real_,
        min_mc_ess_m = NA_real_, status = paste0("scoring_failed: ", score_err)
      )
      next
    }

    if (is.null(target_facet)) {
      lpd_pointwise[fp$test_idx] <- sc$lpd_pointwise
      mcse_pointwise[fp$test_idx] <- sc$mcse_pointwise
      fold_elpd_pt <- sum(sc$lpd_pointwise)
      fold_elpd_grp <- NA_real_
      fold_min_ess_m <- NA_real_
    } else {
      lpd_pointwise[fp$test_idx] <- sc$lpd_pointwise
      mcse_pointwise[fp$test_idx] <- sc$mcse_pointwise
      gs <- sc$group_scores
      gs$fold <- k
      group_tables[[length(group_tables) + 1L]] <- gs
      fold_elpd_pt <- sum(sc$lpd_pointwise)
      fold_elpd_grp <- sum(gs$elpd)
      fold_min_ess_m <- min(gs$mc_ess_m_median, na.rm = TRUE)
    }

    st_label <- if (isTRUE(conv_ok) && isTRUE(dc_k$identifiable)) {
      "ok"
    } else if (identical(conv_ok, FALSE)) {
      "convergence_warning"
    } else if (!isTRUE(dc_k$identifiable)) {
      "rank_deficient_train_fold"
    } else {
      "convergence_unverified"
    }

    fold_diag_rows[[k]] <- tibble::tibble(
      fold = k, n_train = length(fp$train_idx), n_test = length(fp$test_idx),
      n_groups = if (is.null(target_facet)) NA_integer_ else length(unique(as.character(fp$test_df[[target_facet]]))),
      design_connected = isTRUE(dc_k$connected),
      design_identifiable = isTRUE(dc_k$identifiable),
      design_rank = dc_k$adjacent_logit_rank,
      design_params = dc_k$adjacent_logit_parameters,
      converged = conv_ok,
      max_rhat = if (!is.null(conv_k$max_rhat)) conv_k$max_rhat else NA_real_,
      min_ess_bulk = if (!is.null(conv_k$min_ess_bulk)) conv_k$min_ess_bulk else NA_real_,
      n_divergent = if (!is.null(conv_k$n_divergent)) as.integer(conv_k$n_divergent) else NA_integer_,
      elpd_pointwise = fold_elpd_pt,
      elpd_group = fold_elpd_grp,
      min_mc_ess_m = fold_min_ess_m,
      status = st_label
    )
  }

  fold_diag <- dplyr::bind_rows(fold_diag_rows)
  any_refit_failed <- any(grepl("^(refit_failed|scoring_failed):", fold_diag$status))
  all_converged <- !any_refit_failed && all(isTRUE_vec(fold_diag$converged))
  all_identifiable <- all(fold_diag$design_identifiable)

  if (any_refit_failed) {
    warning("One or more CV folds failed refitting or scoring; aggregate ELPD set to NA (failed folds are never silently excluded).", call. = FALSE)
  }
  if (!any_refit_failed && !all_converged) {
    warning("One or more CV folds did not pass convergence checks; inspect $fold_diagnostics before comparing models.", call. = FALSE)
  }
  if (!all_identifiable) {
    warning("One or more training folds were rank-deficient; inspect $fold_diagnostics.", call. = FALSE)
  }

  collateral_unseen <- sort(setdiff(all_unseen_facets, target_facet))
  actual_target <- if (is.null(target_facet)) {
    if (length(all_unseen_facets) == 0L) "new_rating_known_levels" else "new_rating_with_unseen_levels"
  } else {
    if (length(collateral_unseen) == 0L) {
      paste0("new_", target_facet, "_known_other_facets")
    } else {
      paste0("new_", target_facet, "_with_unseen_", paste(collateral_unseen, collapse = "_"))
    }
  }
  known_f <- setdiff(spec$facets_main, all_unseen_facets)
  comp_key <- build_comparison_key(
    method = "kfold_cv",
    target = target,
    target_facet = target_facet,
    actual_target = actual_target,
    known_facets = known_f,
    unseen_facets = all_unseen_facets,
    data_clean = obs_data,
    response_col = spec$response,
    K = get_spec_K_max(spec),
    K_folds = K_folds,
    fold_ids = folds,
    group_ids = if (is.null(target_facet)) NULL else as.character(obs_data[[target_facet]]),
    aggregation = "sum",
    data_info = object$data_info,
    facets_main = spec$facets_main,
    K_by_level = object$data_info$K_by_level %||% spec$K_by_level,
    threshold_facet = spec$threshold_facet
  )

  pointwise_tbl <- tibble::tibble(
    obs_id = comp_key$obs_ids,
    row_index = seq_len(N),
    fold = folds,
    group = if (is.null(target_facet)) NA_character_ else as.character(obs_data[[target_facet]]),
    elpd = lpd_pointwise,
    mcse = mcse_pointwise
  )

  if (any_refit_failed) {
    elpd_pt <- NA_real_
    se_pt <- NA_real_
    se_pt_clust <- NA_real_
    group_scores_tbl <- NULL
    elpd_grp <- NA_real_
    se_grp <- NA_real_
  } else {
    elpd_pt <- sum(lpd_pointwise)
    se_pt <- sqrt(N * stats::var(lpd_pointwise))
    if (is.null(target_facet)) {
      se_pt_clust <- se_pt
      group_scores_tbl <- NULL
      elpd_grp <- NA_real_
      se_grp <- NA_real_
    } else {
      group_scores_tbl <- dplyr::bind_rows(group_tables)
      # Preserve canonical order matching unique(obs_data[[target_facet]])
      canon_g <- unique(as.character(obs_data[[target_facet]]))
      group_scores_tbl <- group_scores_tbl[match(canon_g, group_scores_tbl$group), , drop = FALSE]
      G <- nrow(group_scores_tbl)
      elpd_grp <- sum(group_scores_tbl$elpd)
      se_grp <- if (G >= 2L) sqrt(G * stats::var(group_scores_tbl$elpd)) else NA_real_
      # Cluster-robust SE for rating-level pointwise sum (summed within each held-out group g)
      pt_by_g <- as.numeric(tapply(lpd_pointwise, as.character(obs_data[[target_facet]]), sum)[canon_g])
      se_pt_clust <- if (G >= 2L) sqrt(G * stats::var(pt_by_g)) else NA_real_
      if (M > 1L && any(group_scores_tbl$mc_ess_m_median < 5, na.rm = TRUE)) {
        warning(
          "Low Monte Carlo integration ESS (group median < 5 out of M = ", M,
          ") detected in one or more held-out groups; increase M and verify stability across M and seeds ",
          "(note: passing this threshold alone does not guarantee integration accuracy).",
          call. = FALSE
        )
      }
    }
  }

  structure(
    list(
      target = target,
      target_facet = target_facet,
      actual_target = actual_target,
      K_folds = K_folds,
      folds = folds,
      elpd_pointwise = elpd_pt,
      se_pointwise = se_pt,
      se_pointwise_clustered = se_pt_clust,
      elpd_group = elpd_grp,
      se_group = se_grp,
      pointwise = pointwise_tbl,
      group_scores = group_scores_tbl,
      fold_diagnostics = fold_diag,
      all_folds_converged = all_converged,
      all_folds_identifiable = all_identifiable,
      any_refit_failed = any_refit_failed,
      integration_diagnostics = list(
        M = if (is.null(target_facet) && length(all_unseen_facets) == 0L) NA_integer_ else M,
        n_draws_per_fold = n_draws_per_fold,
        seed = seed,
        note = if (is.null(target_facet) && length(all_unseen_facets) == 0L) {
          "Rating-level CV evaluates conditional log-likelihood directly (no population marginalization)."
        } else {
          "Marginalized via M population draws per posterior draw d; D*M is NOT D*M independent posterior draws."
        }
      ),
      comparison_key = comp_key
    ),
    class = "bmfrm_cv"
  )
}

#' Logical helper treating NA as FALSE
#' @keywords internal
isTRUE_vec <- function(x) {
  !is.na(x) & as.logical(x)
}

#' Print method for bmfrm_cv objects
#' @param x A `bmfrm_cv` object.
#' @param ... Unused.
#' @export
print.bmfrm_cv <- function(x, ...) {
  cat("BayesMFRM Target-Aware K-Fold Cross-Validation\n")
  cat("  Target            :", x$target,
      if (!is.null(x$target_facet)) paste0(" (held-out facet: ", x$target_facet, ")") else "", "\n")
  cat("  Verified task     :", x$actual_target, "\n")
  cat("  Folds             :", x$K_folds, "| All converged:", x$all_folds_converged,
      "| All identifiable:", x$all_folds_identifiable, "\n")
  if (!is.null(x$target_facet)) {
    cat(sprintf("  Group joint ELPD  : %.2f (SE = %.2f across G = %d groups)\n",
                x$elpd_group, x$se_group, if (is.null(x$group_scores)) 0L else nrow(x$group_scores)))
    cat(sprintf("  Rating-level ELPD : %.2f (group-clustered SE = %.2f; naive rating SE = %.2f)\n",
                x$elpd_pointwise, x$se_pointwise_clustered, x$se_pointwise))
    cat("  MC marginalization: M =", x$integration_diagnostics$M,
        "| min median ESS_m =",
        if (!is.null(x$group_scores)) sprintf("%.1f", min(x$group_scores$mc_ess_m_median, na.rm = TRUE)) else "NA", "\n")
  } else {
    cat(sprintf("  Rating-level ELPD : %.2f (SE = %.2f across N = %d ratings)\n",
                x$elpd_pointwise, x$se_pointwise, nrow(x$pointwise)))
  }
  invisible(x)
}

#' Build a minimal loo-compatible wrapper for paired loo::loo_compare() on arbitrary units
#' @keywords internal
make_unit_elpd_loo_object <- function(elpd_vec) {
  n <- length(elpd_vec)
  elpd_sum <- sum(elpd_vec)
  se_val <- if (n >= 2L) sqrt(n * stats::var(elpd_vec)) else NA_real_
  est <- matrix(c(elpd_sum, se_val), nrow = 1L, ncol = 2L,
                dimnames = list("elpd_cv", c("Estimate", "SE")))
  pt <- matrix(elpd_vec, ncol = 1L, dimnames = list(NULL, "elpd_cv"))
  structure(
    list(estimates = est, pointwise = pt),
    class = c("kfold", "loo")
  )
}

#' Guarded model comparison for BayesMFRM PSIS-LOO and CV objects
#'
#' Validates that all candidate models were evaluated on the exact same
#' prediction task (target + verified `actual_target` known/unseen facets),
#' observation set (preserved `obs_id`, response values, category encoding,
#' shared facet/metadata alignment), fold split, held-out groups, and aggregation
#' before delegating to [loo::loo_compare()].
#'
#' Observation correspondence uses preserved pre-preprocessing `obs_id` (an
#' explicit ID column/vector in [bmfrm()], non-default unique row names, or
#' pre-NA row indices), plus per-ID hashes of responses, per-observation rubric
#' category counts (`K_obs`), and shared facet and metadata columns. Permuted
#' observations are aligned by preserved IDs; permutations without preserved
#' IDs (including swaps of equal-score rows) are rejected when shared facets
#' disagree. Compared models must have identical `K_obs` definitions; uniform
#' named-vector PCM remains comparable to scalar-K RSM/PCM.
#'
#' For grouped CV (`target != "rating"`), paired differences are computed across
#' held-out groups (`score = "group"` for joint ELPD, `score = "pointwise"` for
#' rating-level ELPD clustered by group). Multiple ratings from one held-out
#' examinee/rater are not treated as independent. These SEs remain approximate:
#' shared known items/raters and overlapping training folds induce dependence.
#'
#' Do not compare raw ELPD across different prediction targets or methods
#' (e.g., PSIS-LOO and person CV). Inspect `diag_diff` and `diag_elpd`:
#' `abs(elpd_diff) < 4` flags unreliable normal-approximation SEs, not equivalence;
#' `p_worse` is a rounded normal-approximation probability.
#'
#' @param ... Two or more `bmfrm_loo` or `bmfrm_cv` objects (optionally named).
#' @param score Which score unit to compare on `bmfrm_cv` objects: `"auto"`
#'   (default: `"group"` for grouped CV, `"pointwise"` for rating-level CV /
#'   PSIS-LOO), `"group"` (held-out group joint ELPD), or `"pointwise"`
#'   (rating-level ELPD; uses group-clustered paired SE for grouped CV).
#' @param allow_unconverged Logical; if `FALSE` (default), stops when any
#'   candidate CV object has failed/unconverged folds or a `bmfrm_loo` object
#'   failed original MCMC convergence checks.
#' @return A `bmfrm_model_comparison` object (inheriting from `compare.loo`)
#'   with `attr(res, "comparison_meta")` documenting the verified comparison
#'   contract, uncertainty unit, and approximation caveats.
#' @export
compare_bmfrm <- function(...,
                          score = c("auto", "pointwise", "group"),
                          allow_unconverged = FALSE) {
  score <- match.arg(score)
  dots <- list(...)
  if (length(dots) == 1L && is.list(dots[[1L]]) &&
      !inherits(dots[[1L]], c("bmfrm_loo", "bmfrm_cv"))) {
    dots <- dots[[1L]]
  }
  if (length(dots) < 2L) {
    stop("compare_bmfrm() requires at least two bmfrm_loo or bmfrm_cv objects")
  }
  nm <- names(dots)
  if (is.null(nm) || any(!nzchar(nm))) {
    exprs <- as.list(substitute(list(...)))[-1L]
    fallback <- if (length(exprs) == length(dots)) {
      vapply(exprs, function(e) paste(deparse(e), collapse = ""), character(1))
    } else {
      paste0("model", seq_along(dots))
    }
    if (is.null(nm)) nm <- fallback else nm[!nzchar(nm)] <- fallback[!nzchar(nm)]
  }
  names(dots) <- make.unique(nm)

  for (i in seq_along(dots)) {
    if (!inherits(dots[[i]], c("bmfrm_loo", "bmfrm_cv")) || is.null(dots[[i]]$comparison_key)) {
      stop("Argument '", names(dots)[i], "' is not a valid bmfrm_loo or bmfrm_cv object with comparison_key")
    }
  }

  ref_name <- names(dots)[1L]
  ref_key <- dots[[1L]]$comparison_key
  ref_obs_ids <- as.character(ref_key$obs_ids)
  explicit_sources <- c("explicit_column", "explicit_vector", "attribute")

  # 1. Enforce ALL 5 comparison compatibility conditions + obs_id / shared-facet alignment
  perm_indices <- vector("list", length(dots))
  names(perm_indices) <- names(dots)
  for (i in seq_along(dots)) {
    obj <- dots[[i]]
    oname <- names(dots)[i]
    k <- obj$comparison_key

    if (!identical(k$method, ref_key$method)) {
      stop("Cannot compare different evaluation methods: '", ref_name, "' (", ref_key$method,
           ") vs '", oname, "' (", k$method, ").")
    }
    if (!identical(k$target, ref_key$target) ||
        !identical(k$target_facet, ref_key$target_facet) ||
        !identical(k$actual_target, ref_key$actual_target) ||
        !identical(k$unseen_facets, ref_key$unseen_facets)) {
      stop("Cannot compare different prediction targets: '", ref_name, "' (", ref_key$actual_target,
           ") vs '", oname, "' (", k$actual_target, ").")
    }
    if (!identical(k$n_obs, ref_key$n_obs)) {
      stop("Cannot compare models evaluated on different sample sizes ('",
           ref_name, "' vs '", oname, "').")
    }

    k_obs_ids <- as.character(k$obs_ids)
    if (length(k_obs_ids) != ref_key$n_obs || anyNA(k_obs_ids) || anyDuplicated(k_obs_ids)) {
      stop("Cannot compare '", oname, "': observation IDs are missing, wrong length, or duplicated.")
    }
    if (!is.null(ref_key$obs_id_set_hash) && !is.null(k$obs_id_set_hash)) {
      if (!identical(k$obs_id_set_hash, ref_key$obs_id_set_hash)) {
        stop("Cannot compare models evaluated on different observation sets / pre-preprocessing obs_id ('",
             ref_name, "' vs '", oname, "').")
      }
    } else if (!setequal(k_obs_ids, ref_obs_ids)) {
      stop("Cannot compare models evaluated on different observation sets ('",
           ref_name, "' vs '", oname, "').")
    }
    if (!is.null(ref_key$K_obs_hash_by_id) && !is.null(k$K_obs_hash_by_id)) {
      if (!identical(k$K_obs_hash_by_id, ref_key$K_obs_hash_by_id)) {
        stop(
          "Cannot compare models with different per-observation category definitions (K_obs mismatch between '",
          ref_name, "' and '", oname,
          "'). v1 compare_bmfrm() requires identical rubric category definitions across compared models."
        )
      }
    } else if (!identical(k$K, ref_key$K)) {
      stop("Cannot compare models evaluated on different K encoding ('",
           ref_name, "' vs '", oname, "').")
    }
    if (!identical(k$response_hash, ref_key$response_hash)) {
      stop("Cannot compare models evaluated on different response values, obs_id-to-response mapping, or K encoding ('",
           ref_name, "' vs '", oname, "').")
    }

    # Verify observation correspondence across shared facet columns (and shared extra columns)
    if (!is.null(ref_key$facet_hashes_by_id) && !is.null(k$facet_hashes_by_id)) {
      shared_facets <- intersect(names(ref_key$facet_hashes_by_id), names(k$facet_hashes_by_id))
      both_explicit <- (ref_key$obs_id_source %in% explicit_sources) &&
        (k$obs_id_source %in% explicit_sources)
      if (length(shared_facets) == 0L && !both_explicit) {
        stop("Cannot verify observation correspondence between '", ref_name, "' and '", oname,
             "': models share no facet names and do not provide an explicit obs_id.")
      }
      bad_facets <- shared_facets[!vapply(shared_facets, function(f) {
        identical(k$facet_hashes_by_id[[f]], ref_key$facet_hashes_by_id[[f]])
      }, logical(1))]
      if (length(bad_facets) > 0L) {
        stop("Cannot compare '", ref_name, "' and '", oname,
             "': observation alignment mismatch on shared facet(s) [",
             paste(bad_facets, collapse = ", "),
             "] for the same obs_id. Ensure both models use identically ordered rows or pass an explicit unique obs_id.")
      }
    }
    all_ref_cols <- c(ref_key$facet_hashes_by_id, ref_key$extra_col_hashes_by_id)
    all_k_cols <- c(k$facet_hashes_by_id, k$extra_col_hashes_by_id)
    if (length(all_ref_cols) > 0L && length(all_k_cols) > 0L) {
      checked_f <- if (!is.null(ref_key$facet_hashes_by_id) && !is.null(k$facet_hashes_by_id)) {
        intersect(names(ref_key$facet_hashes_by_id), names(k$facet_hashes_by_id))
      } else character(0)
      shared_extra <- setdiff(intersect(names(all_ref_cols), names(all_k_cols)), checked_f)
      bad_extra <- shared_extra[!vapply(shared_extra, function(cn) {
        identical(all_k_cols[[cn]], all_ref_cols[[cn]])
      }, logical(1))]
      if (length(bad_extra) > 0L) {
        stop("Cannot compare '", ref_name, "' and '", oname,
             "': observation alignment mismatch on shared column(s) [",
             paste(bad_extra, collapse = ", "),
             "] for the same obs_id. Ensure both models use identically ordered rows or pass an explicit unique obs_id.")
      }
    }

    if (!identical(k$K_folds, ref_key$K_folds)) {
      stop("Cannot compare CV objects with different fold splits (different K_folds: '", ref_name, "' vs '", oname, "').")
    }
    if (!is.null(ref_key$fold_hash_by_id) && !is.null(k$fold_hash_by_id)) {
      if (!identical(k$fold_hash_by_id, ref_key$fold_hash_by_id)) {
        stop("Cannot compare CV objects with different fold splits ('", ref_name, "' vs '", oname, "').")
      }
    } else if (!identical(k$fold_ids, ref_key$fold_ids)) {
      stop("Cannot compare CV objects with different fold splits ('", ref_name, "' vs '", oname, "').")
    }
    if (!is.null(ref_key$group_hash_by_id) && !is.null(k$group_hash_by_id)) {
      if (!identical(k$group_hash_by_id, ref_key$group_hash_by_id)) {
        stop("Cannot compare CV objects with different held-out group assignments ('", ref_name, "' vs '", oname, "').")
      }
    } else if (!identical(k$group_ids, ref_key$group_ids)) {
      stop("Cannot compare CV objects with different held-out group assignments ('", ref_name, "' vs '", oname, "').")
    }
    if (!identical(k$aggregation, ref_key$aggregation)) {
      stop("Cannot compare objects with different aggregation rules ('", ref_name, "' vs '", oname, "').")
    }

    perm_idx <- match(ref_obs_ids, k_obs_ids)
    if (anyNA(perm_idx)) {
      stop("Cannot align observations between '", ref_name, "' and '", oname, "' by obs_id.")
    }
    perm_indices[[i]] <- perm_idx

    # Diagnostic / failure gates
    if (inherits(obj, "bmfrm_cv")) {
      if (isTRUE(obj$any_refit_failed)) {
        stop("Cannot compare '", oname, "': one or more CV folds failed refitting.")
      }
      if (!isTRUE(allow_unconverged) && !isTRUE(obj$all_folds_converged)) {
        stop("Cannot compare '", oname, "': one or more CV folds did not converge. ",
             "Inspect $fold_diagnostics or pass allow_unconverged = TRUE deliberately.")
      }
    } else if (inherits(obj, "bmfrm_loo")) {
      conv <- obj$bmfrm_meta$mcmc_convergence
      if (!isTRUE(allow_unconverged) && !is.null(conv) && !isTRUE(conv$converged)) {
        stop("Cannot compare '", oname, "': original MCMC fit did not converge. ",
             "Inspect diagnostics or pass allow_unconverged = TRUE deliberately.")
      }
    }
  }

  # 2. Resolve score unit & proper uncertainty unit (with explicit obs_id / group alignment)
  is_grouped_cv <- identical(ref_key$method, "kfold_cv") && !is.null(ref_key$target_facet)
  if (identical(score, "auto")) {
    eff_score <- if (is_grouped_cv) "group" else "pointwise"
  } else {
    eff_score <- score
  }
  if (identical(eff_score, "group") && !is_grouped_cv) {
    stop("score = 'group' is only available for grouped CV objects (target != 'rating')")
  }

  limitations <- character(0)
  if (identical(ref_key$method, "psis_loo")) {
    aligned_dots <- lapply(seq_along(dots), function(i) {
      o <- dots[[i]]
      p <- perm_indices[[i]]
      if (!identical(p, seq_len(ref_key$n_obs))) {
        o$pointwise <- o$pointwise[p, , drop = FALSE]
      }
      o
    })
    names(aligned_dots) <- names(dots)
    comp_mat <- loo::loo_compare(x = aligned_dots)
    uncertainty_unit <- "rating"
    n_units <- ref_key$n_obs
    any_bad_k <- any(vapply(dots, function(o) length(o$bmfrm_meta$bad_k_ids) > 0L, logical(1)))
    if (any_bad_k) {
      limitations <- c(limitations, "High Pareto-k observations present; confirm ranking with exact refit LOO or cv_bmfrm(..., target = 'rating').")
    }
  } else if (!is_grouped_cv) {
    unit_objs <- lapply(seq_along(dots), function(i) {
      o <- dots[[i]]
      make_unit_elpd_loo_object(o$pointwise$elpd[perm_indices[[i]]])
    })
    names(unit_objs) <- names(dots)
    comp_mat <- loo::loo_compare(x = unit_objs)
    uncertainty_unit <- "rating"
    n_units <- ref_key$n_obs
  } else if (identical(eff_score, "group")) {
    canon_g <- as.character(dots[[1L]]$group_scores$group)
    unit_objs <- lapply(seq_along(dots), function(i) {
      o <- dots[[i]]
      g_idx <- match(canon_g, as.character(o$group_scores$group))
      if (anyNA(g_idx)) {
        stop("Cannot align held-out groups between '", ref_name, "' and '", names(dots)[i], "'.")
      }
      make_unit_elpd_loo_object(o$group_scores$elpd[g_idx])
    })
    names(unit_objs) <- names(dots)
    comp_mat <- loo::loo_compare(x = unit_objs)
    uncertainty_unit <- "group"
    n_units <- length(canon_g)
    if (n_units < 10L) {
      limitations <- c(limitations, paste0("Small number of held-out groups (G = ", n_units, "); normal-approximation SEs are rough."))
    }
    limitations <- c(limitations, "Crossed known facets (e.g., shared items/examinees across groups) and overlapping CV training folds induce residual dependence not captured by group-level i.i.d. SEs.")
  } else {
    # Grouped CV with score = "pointwise": align by obs_id, then cluster by held-out group BEFORE paired SE!
    canon_g <- as.character(dots[[1L]]$group_scores$group)
    ref_pt_groups <- as.character(dots[[1L]]$pointwise$group)
    unit_objs <- lapply(seq_along(dots), function(i) {
      o <- dots[[i]]
      aligned_elpd <- o$pointwise$elpd[perm_indices[[i]]]
      by_g <- as.numeric(tapply(aligned_elpd, ref_pt_groups, sum)[canon_g])
      make_unit_elpd_loo_object(by_g)
    })
    names(unit_objs) <- names(dots)
    comp_mat <- loo::loo_compare(x = unit_objs)
    uncertainty_unit <- "group_clustered_pointwise"
    n_units <- length(canon_g)
    if (n_units < 10L) {
      limitations <- c(limitations, paste0("Small number of held-out groups (G = ", n_units, "); cluster-robust SEs are rough."))
    }
    limitations <- c(limitations, "Pointwise scores were clustered within held-out groups for SE estimation; crossed-facet dependence and overlapping CV training folds remain approximation limits.")
  }

  attr(comp_mat, "comparison_meta") <- list(
    method = ref_key$method,
    target = ref_key$target,
    actual_target = ref_key$actual_target,
    known_facets = ref_key$known_facets,
    unseen_facets = ref_key$unseen_facets,
    score_unit = eff_score,
    uncertainty_unit = uncertainty_unit,
    n_units = n_units,
    limitations = limitations
  )
  class(comp_mat) <- unique(c("bmfrm_model_comparison", class(comp_mat)))
  comp_mat
}

#' Print method for bmfrm_model_comparison objects
#' @param x A `bmfrm_model_comparison` object.
#' @param ... Passed to next `print` method.
#' @export
print.bmfrm_model_comparison <- function(x, ...) {
  meta <- attr(x, "comparison_meta")
  cat("BayesMFRM Model Comparison\n")
  cat("  Method / Target   :", meta$method, "|", meta$actual_target, "\n")
  cat("  Score / SE unit   :", meta$score_unit, "| SE over", meta$n_units, meta$uncertainty_unit, "unit(s)\n")
  if (length(meta$limitations)) {
    for (lim in meta$limitations) cat("  Note              :", lim, "\n")
  }
  cat("\n")
  NextMethod("print")
}
