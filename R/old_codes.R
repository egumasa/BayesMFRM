##### HELPER to LOAD model from disk ######

library(cmdstanr)
library(digest)

fit_cmdstan_cached <- function(
  model,           # compiled cmdstan_model
  data,            # list passed to model$sample()
  file,            # path to .rds cache file
  refit = c("never", "always", "on_change"),
  seed = NULL,
  chains = 4,
  parallel_chains = chains,
  ...              # additional cmdstanr::sample() args
) {
  refit <- match.arg(refit)

  # ---- helpers ----
  compute_hash <- function(model, data) {
    digest::digest(list(code = model$code(), data = data))
  }

  get_model_name_safe <- function(model) {
    md <- try(model$metadata(), silent = TRUE)
    if (inherits(md, "try-error") || is.null(md$model_name)) {
      return(NA_character_)
    } else {
      return(md$model_name)
    }
  }

  # ---- ensure parent folder exists ----
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)

  # ---- persistent CSV directory ----
  csv_dir <- file.path(dirname(file), "csv")
  dir.create(csv_dir, recursive = TRUE, showWarnings = FALSE)

  # ---- current model signature ----
  current_hash       <- compute_hash(model, data)
  current_model_name <- get_model_name_safe(model)

  # -------------------------------------------------------------------
  # Case 1: Explicit refit = "always"
  # -------------------------------------------------------------------
  if (refit == "always") {
    message("Refitting model (refit = 'always')...")
    fit <- model$sample(
      data = data,
      seed = seed,
      output_dir = csv_dir,
      chains = chains,
      parallel_chains = parallel_chains,
      ...
    )

    meta <- list(
      hash       = current_hash,
      model_name = current_model_name
    )

    saveRDS(list(fit = fit, meta = meta), file = file)
    return(fit)
  }

  # -------------------------------------------------------------------
  # Case 2: Cached file exists → validate before reuse
  # -------------------------------------------------------------------
  if (file.exists(file)) {
    obj <- readRDS(file)

    # New-style object with meta
    if (is.list(obj) && !is.null(obj$fit) && !is.null(obj$meta$hash)) {

      cached_hash       <- obj$meta$hash
      cached_model_name <- obj$meta$model_name %||% NA_character_

      # --- safety check: model name mismatch ---
      if (!is.na(current_model_name) &&
          !is.na(cached_model_name) &&
          !identical(current_model_name, cached_model_name)) {

        stop(
          "Cached file '", file, "' was created for model '", cached_model_name,
          "', but current model is '", current_model_name, "'.\n",
          "To avoid contamination, please use a different 'file' path or delete the old cache."
        )
      }

      # --- safety check: hash mismatch ---
      if (!identical(cached_hash, current_hash)) {
        if (refit == "never") {
          stop(
            "Cached file '", file, "' does not match current model+data hash, ",
            "but refit = 'never'.\n",
            "Refusing to reuse potentially contaminated cache. ",
            "Use refit = 'on_change' or 'always', or delete the file."
          )
        }

        # refit == "on_change"
        message("Model or data changed → refitting (refit = 'on_change').")
        fit <- model$sample(
          data = data,
          seed = seed,
          output_dir = csv_dir,
          chains = chains,
          parallel_chains = parallel_chains,
          ...
        )

        meta <- list(
          hash       = current_hash,
          model_name = current_model_name
        )

        saveRDS(list(fit = fit, meta = meta), file = file)
        return(fit)
      }

      # --- hashes match → safe to reuse ---
      message("Loading cached fit (model & data unchanged).")
      return(obj$fit)
    }

    # -----------------------------------------------------------------
    # Old-style RDS: just a CmdStanMCMC object, no meta
    # -----------------------------------------------------------------
    if (refit == "never") {
      message("Loading old-style cached fit (no meta; refit = 'never'). ",
              "Cannot fully verify contamination, use with care.")
      return(obj)
    }

    if (refit == "on_change") {
      message("Migrating old cached fit (no meta) → wrapping & adding hash (no refit).")

      # Ensure CSVs exist
      missing <- !file.exists(obj$output_files())
      if (any(missing)) {
        stop("Old cached fit refers to missing CSV files. You must refit the model.")
      }

      meta <- list(
        hash       = current_hash,
        model_name = current_model_name
      )

      obj_new <- list(fit = obj, meta = meta)
      saveRDS(obj_new, file = file)
      return(obj_new$fit)
    }
  } else {
    message("No cached file found → fitting model.")
  }

  # -------------------------------------------------------------------
  # Case 3: Must refit (no cache, or migrated/invalid)
  # -------------------------------------------------------------------
  message("Refitting model...")
  fit <- model$sample(
    data = data,
    seed = seed,
    output_dir = csv_dir,
    chains = chains,
    parallel_chains = parallel_chains,
    ...
  )

  meta <- list(
    hash       = current_hash,
    model_name = current_model_name
  )

  saveRDS(list(fit = fit, meta = meta), file = file)
  return(fit)
}


## ============================================================
## Running ppc
## ============================================================

run_ppc <- function(fit, y_obs, ndraws = 200, seed = 123) {
  # --- 1. Extract posterior draws for y_rep ---
  y_rep_mat <- fit$draws("y_rep", format = "matrix")
  
  # --- 2. Subset draws ---
  set.seed(seed)
  idx <- sample(seq_len(nrow(y_rep_mat)), ndraws)
  yrep_subset <- y_rep_mat[idx, , drop = FALSE]
  
  # --- 3. PPC plots ---
  message("Plotting PPC bars...")
  p1 <- ppc_bars(y = y_obs, yrep = yrep_subset)
  
  message("Plotting PPC mean distribution...")
  p2 <- ppc_stat(y = y_obs, yrep = yrep_subset, stat = "mean")
  
  message("Plotting PPC variance distribution...")
  p3 <- ppc_stat(y = y_obs, yrep = yrep_subset, stat = function(x) var(x))
  
  # --- 4. Return ggplot objects for flexible usage ---
  return(list(
    ppc_bars = p1,
    ppc_mean = p2,
    ppc_var  = p3,
    yrep_subset = yrep_subset
  ))
}

library(dplyr)
library(tidyr)
library(stringr)
library(purrr)


################## FIT STATS #####################

library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(rlang)

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
    transmute(
      n      = get_index(variable),
      mu_hat = mean
    )
  
  sig_df <- sig_summ %>%
    transmute(
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
  obs_df <- tibble(
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
    left_join(mu_df,  by = "n") %>%
    left_join(sig_df, by = "n") %>%
    mutate(
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
      mutate(
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
facet_infit_outfit <- function(res_df, facet) {
  # allow both strings and bare names
  facet_sym <- if (is.character(facet)) sym(facet) else enquo(facet)
  
  res_df %>%
    group_by(!!facet_sym, model) %>%
    summarise(
      n_obs      = n(),
      outfit_msq = mean(z^2, na.rm = TRUE),
      infit_msq  = sum(weight * z^2, na.rm = TRUE) / sum(weight, na.rm = TRUE),
      .groups    = "drop"
    ) %>%
    rename(facet_value = !!facet_sym)
}

#' Bayesian CrIs for infit & outfit using draws stored in residuals
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
      
      tibble(
        facet_value = f_val,
        outfit_lo = quantile(outfit, prob_lo),
        outfit_hi = quantile(outfit, prob_hi),
        infit_lo  = quantile(infit,  prob_lo),
        infit_hi  = quantile(infit,  prob_hi)
      )
    }
  )
  
  fit_stats %>%
    left_join(cri_df, by = "facet_value")
}

#' Add approximate z-standardized infit/outfit
#' (Wright & Masters cube-root transformation, approximate df = n_obs - 1)
add_zstd <- function(fit_stats) {
  fit_stats %>%
    mutate(
      df          = pmax(n_obs - 1, 1),
      # cube-root transform
      outfit_zstd = ( (outfit_msq^(1/3) - 1) * sqrt(df) ) / 1.4142,
      infit_zstd  = ( (infit_msq^(1/3)  - 1) * sqrt(df) ) / 1.4142
    )
}


facet_infit_outfit_all <- function(res_df) {
  facet_cols <- attr(res_df, "facet_cols")
  
  map_dfr(
    facet_cols,
    ~ facet_infit_outfit(res_df, .x) %>%
      mutate(facet = .x),
    .id = NULL
  ) %>%
    relocate(facet, facet_value)
}


library(posterior)

#' Summarize rater × criterion bias effects
#'
#' @param fit          cmdstanr fit with parameter bias[I,R]
#' @param item_labels  optional character vector length I
#' @param rater_labels optional character vector length R
#' @param prob_thresh  posterior probability threshold for flagging (e.g., .95)
#'
#' @return tibble with item, rater, bias_mean, sd, CrI, prob_pos, prob_neg, flag
summarise_bias_ir <- function(fit,
                              item_labels  = NULL,
                              rater_labels = NULL,
                              prob_thresh  = 0.95) {
  # Draws for bias as a matrix
  draws_bias <- fit$draws("bias", format = "draws_array")
  # draws_bias: iterations x chains x (I*R) -> convert to draws_df
  draws_df <- as_draws_df(draws_bias)
  
  # Long format: variable = "bias[i,r]"
  bias_long <- draws_df %>%
    select(starts_with("bias[")) %>%
    pivot_longer(
      cols = everything(),
      names_to = "param",
      values_to = "value"
    )
  
  # Extract i and r indices
  bias_long <- bias_long %>%
    mutate(
      idx = str_match(param, "bias\\[(\\d+),(\\d+)\\]"),
      i   = as.integer(idx[, 2]),
      r   = as.integer(idx[, 3])
    ) %>%
    select(-idx)
  
  # Summaries per (i,r)
  bias_summ <- bias_long %>%
    group_by(i, r) %>%
    summarise(
      bias_mean = mean(value),
      bias_sd   = sd(value),
      lower_95  = quantile(value, 0.025),
      upper_95  = quantile(value, 0.975),
      prob_pos  = mean(value > 0),
      prob_neg  = mean(value < 0),
      .groups   = "drop"
    ) %>%
    mutate(
      # label if CrI excludes 0 and posterior prob is high
      flag_bias = case_when(
        lower_95 > 0 & prob_pos >= prob_thresh ~ "positive_bias",
        upper_95 < 0 & prob_neg >= prob_thresh ~ "negative_bias",
        TRUE                                   ~ "no_clear_bias"
      )
    )
  
  # Attach labels if provided
  if (!is.null(item_labels)) {
    stopifnot(length(item_labels) >= max(bias_summ$i))
    bias_summ <- bias_summ %>%
      mutate(item = item_labels[i])
  } else {
    bias_summ <- bias_summ %>%
      mutate(item = paste0("Item_", i))
  }
  
  if (!is.null(rater_labels)) {
    stopifnot(length(rater_labels) >= max(bias_summ$r))
    bias_summ <- bias_summ %>%
      mutate(rater = rater_labels[r])
  } else {
    bias_summ <- bias_summ %>%
      mutate(rater = paste0("Rater_", r))
  }
  
  bias_summ %>%
    select(item, rater, i, r,
           bias_mean, bias_sd, lower_95, upper_95,
           prob_pos, prob_neg, flag_bias) %>%
    arrange(item, rater)
}

