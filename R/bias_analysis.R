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
#' Applies the single-cell rule from [summarise_bias()] with
#' `flag_method = "prob"` and `delta = 0` (`prob_pos >= prob_thresh` for `"positive_bias"`,
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
                              ci_level     = 0.95,
                              diagnostics = TRUE,
                              mcse_multiplier = 2) {
  validate_bias_diagnostics(diagnostics, mcse_multiplier)
  if (!is.numeric(prob_thresh) || length(prob_thresh) != 1L || !is.finite(prob_thresh) ||
      prob_thresh <= 0.5 || prob_thresh > 1) {
    stop("prob_thresh must be in (0.5, 1]")
  }
  if (!is.numeric(ci_level) || length(ci_level) != 1L || !is.finite(ci_level) ||
      ci_level <= 0 || ci_level >= 1) {
    stop("ci_level must be in (0, 1)")
  }
  ci_alpha <- (1 - ci_level) / 2
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
      ci_lower  = stats::quantile(value, ci_alpha),
      ci_upper  = stats::quantile(value, 1 - ci_alpha),
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
           bias_mean, bias_sd, ci_lower, ci_upper, lower_95, upper_95,
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

#' Plot posterior intervals and screening decisions for interaction cells
#'
#' Each row is one pair of facet levels. The point is the posterior mean and
#' the line is the credible interval selected by `ci_level`. Dashed vertical
#' lines show the practical threshold on the logit scale. A positive interaction
#' raises the adjacent-category log odds relative to the additive prediction.
#' Filled points mark cells selected by the rule in `flag_method`. Open points
#' mark eligible cells without a flag. Crosses show estimates for which a
#' decision was withheld. Interval endpoints and screening use different rules.
#'
#' @param x A table from [summarise_bias()] for one interaction and one set of
#'   thresholds. Subset rows to focus on particular raters or criteria.
#' @return Invisibly returns the plotted rows.
#' @examples
#' data(bmfrm_example)
#' b <- bmfrm_example$bias
#' plot_bias(b)
#' b[b$flag_bias == "positive_bias", c("level1", "level2", "bfdr")]
#' @export
plot_bias <- function(x) {
  required <- c("facet1", "facet2", "level1", "level2", "bias_mean", "ci_lower",
                "ci_upper", "ci_level", "delta", "flag_method", "flag_bias",
                "n_observations", "likelihood_identified")
  if (!is.data.frame(x) || !nrow(x) || !all(required %in% names(x))) {
    stop("Supply one or more rows from summarise_bias()")
  }
  settings <- intersect(c("facet1", "facet2", "delta", "ci_level", "flag_method",
                          "fdr_thresh", "prob_thresh"), names(x))
  if (nrow(unique(x[settings])) != 1L) stop("Plot one interaction and one set of screening settings at a time")
  eligible <- x$likelihood_identified %in% TRUE & x$n_observations > 0L
  flagged <- eligible & x$flag_bias %in% c("positive_bias", "negative_bias")
  symbols <- ifelse(flagged, 19L, ifelse(eligible, 1L, 4L))
  colors <- ifelse(flagged, "#125B72", ifelse(eligible, "#333333", "#747474"))
  y <- rev(seq_len(nrow(x)))
  labels <- paste(x$level1, x$level2, sep = " / ")
  delta <- x$delta[1L]
  limits <- range(x$ci_lower, x$ci_upper, x$bias_mean, -delta, delta, finite = TRUE)
  bfdr_rule <- x$flag_method[1L] == "bfdr"
  rule <- if (bfdr_rule) "BFDR flag" else "Probability flag"
  cutoff <- if (bfdr_rule) x[["fdr_thresh"]] else x[["prob_thresh"]]
  if (length(cutoff)) rule <- sprintf("%s (%g%%)", rule, 100 * cutoff[1L])
  xlab <- "Interaction effect (logits)"
  if (delta > 0) xlab <- sprintf("%s   (dashed lines at +/-%g)", xlab, delta)
  old <- graphics::par(mar = c(4.3, 2, 3.6, 1))
  on.exit(graphics::par(old))
  left <- max(graphics::strwidth(labels, units = "inches", cex = 0.85)) + 0.35
  margins <- graphics::par("mai")
  margins[2L] <- min(left, graphics::par("fin")[1L] * 0.45)
  graphics::par(mai = margins)
  graphics::plot(NA, xlim = limits, ylim = c(0.5, length(y) + 0.5), yaxt = "n",
                 ylab = "", xlab = xlab, bty = "n")
  graphics::title(main = sprintf("%s by %s (%g%% credible intervals)",
                                 x$facet1[1L], x$facet2[1L], 100 * x$ci_level[1L]), line = 2)
  bounds <- graphics::par("usr")
  graphics::legend(mean(bounds[1:2]), bounds[4L], xpd = NA, xjust = 0.5, yjust = 0,
                   legend = c(rule, "Unflagged", "Withheld"),
                   pch = c(19, 1, 4), col = c("#125B72", "#333333", "#747474"),
                   horiz = TRUE, bty = "n", cex = 0.8)
  graphics::abline(v = 0, col = "#A6A6A6")
  if (delta > 0) graphics::abline(v = c(-delta, delta), lty = 2, col = "#93733B")
  graphics::axis(2, at = y, labels = labels, las = 1, tick = FALSE, cex.axis = 0.85)
  graphics::segments(x$ci_lower, y, x$ci_upper, y, col = colors)
  graphics::points(x$bias_mean, y, pch = symbols, col = colors)
  invisible(x)
}
