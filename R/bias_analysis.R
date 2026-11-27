#' Classify directional bias flags using a single consistent probability threshold
#' @keywords internal
classify_bias_flags <- function(prob_pos, prob_neg, prob_thresh = 0.95) {
  if (!is.numeric(prob_thresh) || length(prob_thresh) != 1L || !is.finite(prob_thresh) ||
      prob_thresh <= 0.5 || prob_thresh > 1) {
    stop("prob_thresh must be in (0.5, 1]")
  }
  unname(ifelse(
    prob_pos >= prob_thresh, "positive_bias",
    ifelse(prob_neg >= prob_thresh, "negative_bias", "no_clear_bias")
  ))
}

# Internal shared diagnostics; keep the iteration x chain layout throughout.
validate_bias_diagnostics <- function(diagnostics, mcse_multiplier) {
  if (!is.logical(diagnostics) || length(diagnostics) != 1L || is.na(diagnostics)) {
    stop("diagnostics must be TRUE or FALSE")
  }
  if (!is.numeric(mcse_multiplier) || length(mcse_multiplier) != 1L ||
      !is.finite(mcse_multiplier) || mcse_multiplier <= 0) {
    stop("mcse_multiplier must be one finite positive number")
  }
}

bias_probability_diagnostics <- function(draws, variables, prob_pos, prob_neg,
                                         delta, prob_thresh, mcse_multiplier,
                                         convergence = NULL) {
  n <- length(variables)
  mcse <- matrix(NA_real_, n, 3L)
  dims <- dim(draws)
  has_chains <- length(dims) == 3L && dims[1L] > 0L && dims[2L] > 0L &&
    !is.null(dimnames(draws)[[3L]]) &&
    !anyDuplicated(dimnames(draws)[[3L]]) &&
    all(variables %in% dimnames(draws)[[3L]])
  if (has_chains) {
    indicator_mcse <- function(x) {
      # Constant chains do not establish zero Monte Carlo error, even when
      # different chains are stuck on different sides of the threshold.
      if (nrow(x) < 4L || any(!is.finite(x)) ||
          any(vapply(seq_len(ncol(x)), function(j) all(x[, j] == x[1L, j]), logical(1)))) {
        return(NA_real_)
      }
      ans <- posterior::mcse_mean(x * 1)
      if (is.finite(ans) && ans >= 0) ans else NA_real_
    }
    for (j in seq_len(n)) {
      x <- matrix(draws[, , variables[j]], nrow = dims[1L], ncol = dims[2L])
      mcse[j, ] <- c(indicator_mcse(x > delta), indicator_mcse(x < -delta),
                     indicator_mcse(abs(x) > delta))
    }
  }

  flags <- classify_bias_flags(prob_pos, prob_neg, prob_thresh)
  status <- rep("mcse_unavailable", n)
  for (j in seq_len(n)) {
    p <- c(prob_pos[j], prob_neg[j])
    se <- mcse[j, 1:2]
    # A positive/negative decision depends only on its own directional
    # probability because prob_thresh > 0.5 and delta >= 0.
    needed <- switch(flags[j], positive_bias = 1L, negative_bias = 2L, 1:2)
    if (!all(is.finite(se[needed]))) next
    lower <- p - mcse_multiplier * se
    upper <- p + mcse_multiplier * se
    stable <- switch(flags[j],
      positive_bias = lower[1L] > prob_thresh,
      negative_bias = lower[2L] > prob_thresh,
      all(upper < prob_thresh))
    status[j] <- if (stable) paste0("stable_", flags[j]) else "boundary_sensitive"
  }
  # Numerical stability is conditional on the fit's existing sampler checks.
  # Identification remains a separate property in summarise_bias()'s output.
  if (!has_chains) {
    status[] <- "chain_structure_unavailable"
  } else if (dims[2L] < 2L || is.null(convergence$converged) ||
             is.na(convergence$converged) || identical(convergence$status, "incomplete")) {
    status[] <- "fit_diagnostics_unavailable"
  } else if (!isTRUE(convergence$converged)) {
    status[] <- "fit_diagnostics_failed"
  }
  stable <- rep(NA, n)
  stable[grepl("^stable_", status)] <- TRUE
  stable[status == "boundary_sensitive"] <- FALSE
  tibble::tibble(
    mcse_prob_gt_delta = mcse[, 1L], mcse_prob_lt_neg_delta = mcse[, 2L],
    mcse_prob_outside_delta = mcse[, 3L], prob_thresh = prob_thresh,
    mcse_multiplier = mcse_multiplier, decision_stable = stable,
    decision_status = status
  )
}

#' Summarize rater x criterion bias effects (legacy raw `bias[I,R]` fits)
#'
#' Applies the same single directional posterior probability rule as
#' [summarise_bias()] (`prob_pos >= prob_thresh` for `"positive_bias"`,
#' `prob_neg >= prob_thresh` for `"negative_bias"`). Pass `prob_thresh = 0.975`
#' to match a two-sided 95% credible interval excluding 0.
#'
#' @param fit          cmdstanr fit with parameter `bias[I,R]`
#' @param item_labels  optional character vector length I
#' @param rater_labels optional character vector length R
#' @param prob_thresh  posterior directional probability threshold in (0.5, 1]
#'   for flagging (default 0.95; e.g., use 0.975 to match a two-sided 95% CrI
#'   excluding 0)
#' @inheritParams summarise_bias
#'
#' @return tibble with item, rater, bias_mean, sd, CrI, prob_pos, prob_neg, flag.
#'   With `diagnostics = TRUE`, also returns `mcse_prob_pos`, `mcse_prob_neg`,
#'   `prob_thresh`, `mcse_multiplier`, `decision_stable`, and `decision_status`.
#' @details Uses the same probability MCSE and numerical stability rules as
#'   [summarise_bias()], with `delta = 0`. Sampler checks are obtained from the
#'   raw fit; if unavailable, stability is `NA`. This legacy interface does not
#'   check likelihood identification.
#' @export
summarise_bias_ir <- function(fit,
                              item_labels  = NULL,
                              rater_labels = NULL,
                              prob_thresh  = 0.95,
                              diagnostics = TRUE,
                              mcse_multiplier = 2) {
  validate_bias_diagnostics(diagnostics, mcse_multiplier)
  if (!is.numeric(prob_thresh) || length(prob_thresh) != 1L || !is.finite(prob_thresh) ||
      prob_thresh <= 0.5 || prob_thresh > 1) {
    stop("prob_thresh must be in (0.5, 1]")
  }
  # Draws for bias as a matrix
  draws_bias <- fit$draws("bias", format = "draws_array")
  # draws_bias: iterations x chains x (I*R) -> convert to draws_df
  draws_df <- posterior::as_draws_df(draws_bias)
  
  # Long format: variable = "bias[i,r]"
  bias_long <- tibble::as_tibble(draws_df) %>%
    dplyr::select(dplyr::starts_with("bias[")) %>%
    tidyr::pivot_longer(
      cols = dplyr::everything(),
      names_to = "param",
      values_to = "value"
    )
  
  # Extract i and r indices
  bias_long <- bias_long %>%
    dplyr::mutate(
      idx = stringr::str_match(param, "bias\\[(\\d+),(\\d+)\\]"),
      i   = as.integer(idx[, 2]),
      r   = as.integer(idx[, 3])
    ) %>%
    dplyr::select(-idx)
  
  # Summaries per (i,r)
  bias_summ <- bias_long %>%
    dplyr::group_by(i, r) %>%
    dplyr::summarise(
      bias_mean = mean(value),
      bias_sd   = stats::sd(value),
      lower_95  = stats::quantile(value, 0.025),
      upper_95  = stats::quantile(value, 0.975),
      prob_pos  = mean(value > 0),
      prob_neg  = mean(value < 0),
      .groups   = "drop"
    ) %>%
    dplyr::mutate(
      # Unified single directional posterior probability criterion (matches summarise_bias())
      flag_bias = classify_bias_flags(prob_pos, prob_neg, prob_thresh = prob_thresh)
    )
  
  # Attach labels if provided
  if (!is.null(item_labels)) {
    stopifnot(length(item_labels) >= max(bias_summ$i))
    bias_summ <- bias_summ %>%
      dplyr::mutate(item = item_labels[i])
  } else {
    bias_summ <- bias_summ %>%
      dplyr::mutate(item = paste0("Item_", i))
  }
  
  if (!is.null(rater_labels)) {
    stopifnot(length(rater_labels) >= max(bias_summ$r))
    bias_summ <- bias_summ %>%
      dplyr::mutate(rater = rater_labels[r])
  } else {
    bias_summ <- bias_summ %>%
      dplyr::mutate(rater = paste0("Rater_", r))
  }
  
  bias_summ <- bias_summ %>%
    dplyr::select(item, rater, i, r,
           bias_mean, bias_sd, lower_95, upper_95,
           prob_pos, prob_neg, flag_bias)
  ord_i <- if (is.null(item_labels)) bias_summ$i else bias_summ$item
  ord_r <- if (is.null(rater_labels)) bias_summ$r else bias_summ$rater
  bias_summ <- bias_summ[order(ord_i, ord_r), , drop = FALSE]

  if (diagnostics) {
    convergence <- tryCatch(check_convergence(fit), error = function(e) NULL)
    diagnostic <- bias_probability_diagnostics(
      draws_bias, paste0("bias[", bias_summ$i, ",", bias_summ$r, "]"),
      bias_summ$prob_pos, bias_summ$prob_neg, delta = 0,
      prob_thresh = prob_thresh, mcse_multiplier = mcse_multiplier,
      convergence = convergence
    )
    diagnostic$mcse_prob_outside_delta <- NULL
    names(diagnostic)[1:2] <- c("mcse_prob_pos", "mcse_prob_neg")
    bias_summ <- dplyr::bind_cols(bias_summ, diagnostic)
  }
  bias_summ
}
