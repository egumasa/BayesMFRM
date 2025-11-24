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
    mu_mat  <- fit$draws("mu",     format = "matrix")
    sig_mat <- fit$draws("sigma2", format = "matrix")
    
    # sort columns by [n]
    mu_mat  <- mu_mat[,  order(get_index(colnames(mu_mat))), drop=FALSE]
    sig_mat <- sig_mat[, order(get_index(colnames(sig_mat))), drop=FALSE]
    
    # attach list columns
    out$mu_draws    <- split(as.data.frame(t(mu_mat)),    seq_len(N))
    out$sigma2_draws <- split(as.data.frame(t(sig_mat)), seq_len(N))
    
    # each list element will be a numeric vector of draws
    out <- out %>%
      dplyr::mutate(
        mu_draws    = purrr::map(mu_draws,    unlist),
        sigma2_draws = purrr::map(sigma2_draws, unlist)
      )
  }
  
  # attach metadata: facet columns
  attr(out, "facet_cols") <- facet_cols
  
  out
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
  # allow both strings and bare names
  facet_sym <- if (is.character(facet)) rlang::sym(facet) else rlang::enquo(facet)
  
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

#' Bayesian CrIs for infit & outfit using draws stored in residuals
#' 
#' @param fit_stats tibble from facet_infit_outfit()
#' @param res_df tibble from extract_residuals() with draws
#' @param facet facet name
#' @param prob probability for credible intervals
#' @return fit_stats tibble with added CrI columns
#' @export
add_fit_cri_bayes <- function(fit_stats,
                              res_df,
                              facet,
                              prob = 0.95) {
  facet_name <- if (is.character(facet)) facet else rlang::as_name(rlang::enquo(facet))
  
  prob_lo <- (1 - prob) / 2
  prob_hi <- 1 - prob_lo
  
  facet_levels <- sort(unique(res_df[[facet_name]]))
  
  cri_df <- purrr::map_dfr(
    facet_levels,
    function(f_val) {
      idx <- which(res_df[[facet_name]] == f_val)
      
      mu_list  <- res_df$mu_draws[idx]      # list of numeric vectors
      sig_list <- res_df$sigma2_draws[idx]
      
      # matrix: draws x obs_in_facet
      mu_mat  <- do.call(cbind, mu_list)
      sig_mat <- do.call(cbind, sig_list)
      
      D <- nrow(mu_mat)
      n_obs <- length(idx)
      
      X_sub <- matrix(res_df$x[idx], nrow=D, ncol=n_obs, byrow=TRUE)
      
      z2 <- (X_sub - mu_mat)^2 / sig_mat
      w  <- sig_mat
      
      outfit <- rowMeans(z2)
      infit  <- rowSums(w * z2) / rowSums(w)
      
      tibble::tibble(
        facet_value = f_val,
        outfit_lo = stats::quantile(outfit, prob_lo),
        outfit_hi = stats::quantile(outfit, prob_hi),
        infit_lo  = stats::quantile(infit,  prob_lo),
        infit_hi  = stats::quantile(infit,  prob_hi)
      )
    }
  )
  
  fit_stats %>%
    dplyr::left_join(cri_df, by = "facet_value")
}

#' Add approximate z-standardized infit/outfit
#' (Wright & Masters cube-root transformation, approximate df = n_obs - 1)
#' 
#' @param fit_stats tibble from facet_infit_outfit()
#' @return fit_stats tibble with added z-standardized columns
#' @export
add_zstd <- function(fit_stats) {
  fit_stats %>%
    dplyr::mutate(
      df          = pmax(n_obs - 1, 1),
      # cube-root transform
      outfit_zstd = ( (outfit_msq^(1/3) - 1) * sqrt(df) ) / 1.4142,
      infit_zstd  = ( (infit_msq^(1/3)  - 1) * sqrt(df) ) / 1.4142
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
