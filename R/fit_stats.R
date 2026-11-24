#' Extract residual info from a cmdstanr fit (with optional posterior draws)
#'
#' @param fit         cmdstanr fit object (must have mu[n], sigma2[n] in GQ)
#' @param stan_data   list with X and one or more facet columns
#' @param facet_cols  which columns from stan_data are facets
#' @param save_draws  if TRUE, save posterior draws of mu and sigma2 per obs
#' @param model_name  optional model label
#'
#' @return tibble with:
#'   n, x, mu_hat, sigma2_hat, resid, z, weight, model,
#'   facet columns,
#'   and (if save_draws=TRUE):
#'       mu_draws[n]      list-column of length-D vectors
#'       sigma2_draws[n]  list-column of length-D vectors
#' @export
extract_residuals <- function(
    fit,
    stan_data,
    facet_cols = NULL,
    save_draws = TRUE,
    model_name = deparse(substitute(fit))
) {
  # ───────────────────────────────────────────────────────────────
  # 1. Posterior summaries for mu and sigma2
  # ───────────────────────────────────────────────────────────────
  mu_summ  <- fit$summary("mu")
  sig_summ <- fit$summary("sigma2")
  
  get_index <- function(x) as.integer(stringr::str_extract(x, "(?<=\\[)\\d+(?=\\])"))
  
  mu_df <- mu_summ %>%
    dplyr::transmute(
      n      = get_index(variable),
      mu_hat = mean
    )
  
  sig_df <- sig_summ %>%
    dplyr::transmute(
      n          = get_index(variable),
      sigma2_hat = mean
    )
  
  # ───────────────────────────────────────────────────────────────
  # 2. Determine facets
  # ───────────────────────────────────────────────────────────────
  N <- length(stan_data$X)
  
  if (is.null(facet_cols)) {
    lens <- vapply(stan_data, length, integer(1))
    facet_cols <- names(stan_data)[lens == N & names(stan_data) != "X"]
  }
  
  # ───────────────────────────────────────────────────────────────
  # 3. Build observation-level table
  # ───────────────────────────────────────────────────────────────
  obs_df <- tibble::tibble(
    n = seq_len(N),
    x = stan_data$X
  )
  
  for (nm in facet_cols) {
    obs_df[[nm]] <- stan_data[[nm]]
  }
  
  # ───────────────────────────────────────────────────────────────
  # 4. Add point summaries
  # ───────────────────────────────────────────────────────────────
  out <- obs_df %>%
    dplyr::left_join(mu_df,  by = "n") %>%
    dplyr::left_join(sig_df, by = "n") %>%
    dplyr::mutate(
      resid  = x - mu_hat,
      weight = sigma2_hat,
      z      = resid / sqrt(sigma2_hat),
      model  = model_name
    )
  
  # ───────────────────────────────────────────────────────────────
  # 5. Optional: Attach posterior draws (compact list-column)
  # ───────────────────────────────────────────────────────────────
  if (save_draws) {
    # raw draws (iterations x N)
    mu_mat  <- unclass(fit$draws("mu",     format = "matrix"))
    sig_mat <- unclass(fit$draws("sigma2", format = "matrix"))
    
    # sort columns by [n]
    mu_mat  <- mu_mat[,  order(get_index(colnames(mu_mat))), drop = FALSE]
    sig_mat <- sig_mat[, order(get_index(colnames(sig_mat))), drop = FALSE]
    
    # attach list columns without transposing to an intermediate data.frame
    out$mu_draws     <- lapply(seq_len(N), function(j) as.numeric(mu_mat[, j]))
    out$sigma2_draws <- lapply(seq_len(N), function(j) as.numeric(sig_mat[, j]))
  }
  
  # attach metadata: facet columns
  attr(out, "facet_cols") <- facet_cols
  
  out
}

#' Resolve a facet column argument passed either as a bare name or a string
#' @keywords internal
resolve_facet_name_nse <- function(facet_quo, col_names = character()) {
  expr <- rlang::quo_get_expr(facet_quo)
  if (rlang::is_string(expr)) {
    return(expr)
  }
  if (rlang::is_symbol(expr)) {
    sym_name <- rlang::as_name(expr)
    if (sym_name %in% col_names) {
      return(sym_name)
    }
    val <- tryCatch(rlang::eval_tidy(facet_quo), error = function(e) NULL)
    if (is.character(val) && length(val) == 1L && !is.na(val)) {
      return(val)
    }
    return(sym_name)
  }
  val <- rlang::eval_tidy(facet_quo)
  if (is.character(val) && length(val) == 1L && !is.na(val)) {
    return(val)
  }
  rlang::as_name(facet_quo)
}

#' Compute infit and outfit mean square for any facet
#'
#' @param res_df tibble from extract_residuals()
#' @param facet  facet column to group by; can be a bare name or a string
#'               e.g. RaterID, "RaterID", ExamineeID, "TaskID", "CriterionID", ...
#'
#' @return tibble with facet_value, model, n_obs, infit_msq, outfit_msq
#' @export
facet_infit_outfit <- function(res_df, facet) {
  # Capture before evaluation so both bare column names (rater) and strings ("rater") work
  facet_name <- resolve_facet_name_nse(rlang::enquo(facet), names(res_df))
  facet_sym <- rlang::sym(facet_name)
  
  res_df %>%
    dplyr::group_by(!!facet_sym, model) %>%
    dplyr::summarise(
      n_obs      = dplyr::n(),
      outfit_msq = mean(z^2, na.rm = TRUE),
      infit_msq  = sum(weight * z^2, na.rm = TRUE) / sum(weight, na.rm = TRUE),
      .groups    = "drop"
    ) %>%
    dplyr::rename(facet_value = !!facet_sym)
}

#' Bayesian CrIs and posterior means for infit & outfit using draws stored in residuals
#'
#' Computes draw-level infit and outfit mean squares for each facet level,
#' returning posterior means alongside credible intervals so point estimates
#' and intervals come from the same posterior distribution (avoiding Jensen's
#' inequality discrepancies between plug-in means and draw intervals).
#'
#' @param fit_stats tibble from facet_infit_outfit()
#' @param res_df tibble from extract_residuals() with draws
#' @param facet facet name (bare column name or string)
#' @param prob probability for credible intervals
#' @return fit_stats tibble with added posterior summary and CrI columns
#' @export
add_fit_cri_bayes <- function(fit_stats,
                              res_df,
                              facet,
                              prob = 0.95) {
  facet_name <- resolve_facet_name_nse(rlang::enquo(facet), names(res_df))

  prob_lo <- (1 - prob) / 2
  prob_hi <- 1 - prob_lo

  has_model <- "model" %in% names(res_df) && "model" %in% names(fit_stats)
  models <- if (has_model) unique(res_df$model) else NA_character_

  cri_df <- purrr::map_dfr(models, function(mod) {
    sub_df <- if (has_model) res_df[res_df$model == mod, , drop = FALSE] else res_df
    facet_levels <- sort(unique(sub_df[[facet_name]]))
    purrr::map_dfr(
      facet_levels,
      function(f_val) {
        idx <- which(sub_df[[facet_name]] == f_val)

        mu_mat  <- do.call(cbind, sub_df$mu_draws[idx])
        sig_mat <- do.call(cbind, sub_df$sigma2_draws[idx])

        D <- nrow(mu_mat)
        n_obs <- length(idx)

        X_sub <- matrix(sub_df$x[idx], nrow = D, ncol = n_obs, byrow = TRUE)

        z2 <- (X_sub - mu_mat)^2 / sig_mat
        w  <- sig_mat

        outfit <- rowMeans(z2)
        infit  <- rowSums(w * z2) / rowSums(w)

        row <- tibble::tibble(
          facet_value = f_val,
          outfit_post_mean = mean(outfit),
          outfit_lo = unname(stats::quantile(outfit, prob_lo)),
          outfit_hi = unname(stats::quantile(outfit, prob_hi)),
          infit_post_mean = mean(infit),
          infit_lo  = unname(stats::quantile(infit,  prob_lo)),
          infit_hi  = unname(stats::quantile(infit,  prob_hi))
        )
        if (has_model) row$model <- mod
        row
      }
    )
  })

  join_by <- if (has_model) c("facet_value", "model") else "facet_value"
  fit_stats %>%
    dplyr::left_join(cri_df, by = join_by)
}

#' Compute facet-level infit/outfit directly from draw matrices without list-columns
#'
#' @param x Integer vector of observed ratings (length N).
#' @param facet_id Integer vector of facet indices (length N).
#' @param mu_mat Numeric matrix (D x N) of expected scores per draw.
#' @param sig_mat Numeric matrix (D x N) of model variances per draw.
#' @param kurt_mat Optional numeric matrix (D x N) of fourth central moments
#'   `E[(Y - mu)^4]` per draw for exact Wright & Masters ZSTD variance `q`.
#' @param model_name Character model label.
#' @param prob Credible interval probability (default 0.95).
#' @return Tibble with draw-consistent posterior mean infit/outfit, credible
#'   intervals, plug-in values, and optional exact `outfit_q` / `infit_q`.
#' @keywords internal
facet_fit_from_matrices <- function(x, facet_id, mu_mat, sig_mat,
                                    kurt_mat = NULL, model_name = "model",
                                    prob = 0.95) {
  prob_lo <- (1 - prob) / 2
  prob_hi <- 1 - prob_lo
  D <- nrow(mu_mat)
  mu_hat <- colMeans(mu_mat)
  sig_hat <- colMeans(sig_mat)
  kurt_hat <- if (!is.null(kurt_mat)) colMeans(kurt_mat) else NULL

  levels_obs <- sort(unique(facet_id))
  rows <- lapply(levels_obs, function(f_val) {
    idx <- which(facet_id == f_val)
    n_obs <- length(idx)
    x_sub <- x[idx]
    mu_sub <- mu_mat[, idx, drop = FALSE]
    sig_sub <- sig_mat[, idx, drop = FALSE]

    # Plug-in mean squares
    z2_plug <- (x_sub - mu_hat[idx])^2 / sig_hat[idx]
    w_plug <- sig_hat[idx]
    outfit_plug <- mean(z2_plug, na.rm = TRUE)
    infit_plug <- sum(w_plug * z2_plug, na.rm = TRUE) / sum(w_plug, na.rm = TRUE)

    # Posterior draw-level mean squares (consistent with credible intervals)
    X_mat <- matrix(x_sub, nrow = D, ncol = n_obs, byrow = TRUE)
    z2_draw <- (X_mat - mu_sub)^2 / sig_sub
    outfit_draws <- rowMeans(z2_draw)
    infit_draws <- rowSums(sig_sub * z2_draw) / rowSums(sig_sub)

    out <- tibble::tibble(
      facet_value = f_val,
      model = model_name,
      n_obs = n_obs,
      infit_msq = mean(infit_draws),
      outfit_msq = mean(outfit_draws),
      infit_lo = unname(stats::quantile(infit_draws, prob_lo)),
      infit_hi = unname(stats::quantile(infit_draws, prob_hi)),
      outfit_lo = unname(stats::quantile(outfit_draws, prob_lo)),
      outfit_hi = unname(stats::quantile(outfit_draws, prob_hi)),
      infit_msq_plugin = infit_plug,
      outfit_msq_plugin = outfit_plug
    )
    if (!is.null(kurt_mat)) {
      kurt_sub <- kurt_mat[, idx, drop = FALSE]
      w2_draw <- sig_sub * sig_sub
      # Compute q_d and Wilson-Hilferty ZSTD inside each posterior draw
      # (avoids swapping expectations: E[C]/E[W]^2 != E[C/W^2])
      q_out_draws <- sqrt(pmax(rowSums(pmax(kurt_sub / w2_draw - 1, 0)) / (n_obs^2), 1e-12))
      q_in_draws <- sqrt(pmax(rowSums(pmax(kurt_sub - w2_draw, 0)) / (rowSums(sig_sub)^2), 1e-12))
      out_z_draws <- (outfit_draws^(1/3) - 1) * (3 / q_out_draws) + q_out_draws / 3
      in_z_draws <- (infit_draws^(1/3) - 1) * (3 / q_in_draws) + q_in_draws / 3
      out$outfit_q <- mean(q_out_draws)
      out$infit_q <- mean(q_in_draws)
      out$infit_zstd <- mean(in_z_draws)
      out$outfit_zstd <- mean(out_z_draws)
      out$infit_zstd_lo <- unname(stats::quantile(in_z_draws, prob_lo))
      out$infit_zstd_hi <- unname(stats::quantile(in_z_draws, prob_hi))
      out$outfit_zstd_lo <- unname(stats::quantile(out_z_draws, prob_lo))
      out$outfit_zstd_hi <- unname(stats::quantile(out_z_draws, prob_hi))
      out$zstd_method <- "draw_wilson_hilferty_approx"
    }
    out
  })
  dplyr::bind_rows(rows)
}

#' Add approximate Wilson-Hilferty / Wright & Masters z-standardized infit/outfit (ZSTD)
#'
#' Applies the Wright & Masters (1982) Wilson-Hilferty cube-root approximation:
#' `ZSTD = (MSQ^(1/3) - 1) * (3 / q) + q / 3`.
#' When `outfit_q` / `infit_q` are already present, uses them; otherwise falls
#' back to the chi-square approximation `q = sqrt(2 / df)` (`df = max(n_obs - 1, 1)`).
#'
#' Note: Wilson-Hilferty ZSTD is an asymptotic approximation whose calibration
#' degrades for small `n_obs` or sparse designs with high posterior uncertainty;
#' it should be treated as a descriptive flag alongside conditional posterior
#' predictive checks (`bayesplot::pp_check(..., group = ...)`) rather than an exact test.
#'
#' @param fit_stats tibble from facet_infit_outfit() or facet_fit_from_matrices()
#' @return fit_stats tibble with added approximate z-standardized columns
#' @export
add_zstd <- function(fit_stats) {
  df <- pmax(fit_stats$n_obs - 1, 1)
  q_approx <- sqrt(2 / df)
  q_outfit <- if ("outfit_q" %in% names(fit_stats)) fit_stats$outfit_q else q_approx
  q_infit <- if ("infit_q" %in% names(fit_stats)) fit_stats$infit_q else q_approx
  q_outfit <- pmax(q_outfit, 1e-12)
  q_infit <- pmax(q_infit, 1e-12)
  method_label <- if ("outfit_q" %in% names(fit_stats)) "wilson_hilferty_kurtosis_approx" else "wilson_hilferty_chisq_approx"

  fit_stats %>%
    dplyr::mutate(
      df          = df,
      outfit_zstd = (outfit_msq^(1/3) - 1) * (3 / q_outfit) + q_outfit / 3,
      infit_zstd  = (infit_msq^(1/3)  - 1) * (3 / q_infit)  + q_infit / 3,
      zstd_method = method_label
    )
}

#' Compute infit/outfit for all facets in residual table
#'
#' @param res_df tibble from extract_residuals()
#' @return tibble with all facets' infit/outfit statistics
#' @export
facet_infit_outfit_all <- function(res_df) {
  facet_cols <- attr(res_df, "facet_cols")
  
  purrr::map_dfr(
    facet_cols,
    ~ facet_infit_outfit(res_df, .x) %>%
      dplyr::mutate(facet = .x),
    .id = NULL
  ) %>%
    dplyr::relocate(facet, facet_value)
}
