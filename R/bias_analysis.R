#' Summarize rater × criterion bias effects
#'
#' @param fit          cmdstanr fit with parameter bias[I,R]
#' @param item_labels  optional character vector length I
#' @param rater_labels optional character vector length R
#' @param prob_thresh  posterior probability threshold for flagging (e.g., .95)
#'
#' @return tibble with item, rater, bias_mean, sd, CrI, prob_pos, prob_neg, flag
#' @export
summarise_bias_ir <- function(fit,
                              item_labels  = NULL,
                              rater_labels = NULL,
                              prob_thresh  = 0.95) {
  # Draws for bias as a matrix
  draws_bias <- fit$draws("bias", format = "draws_array")
  # draws_bias: iterations x chains x (I*R) -> convert to draws_df
  draws_df <- posterior::as_draws_df(draws_bias)
  
  # Long format: variable = "bias[i,r]"
  bias_long <- draws_df %>%
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
      bias_sd   = sd(value),
      lower_95  = stats::quantile(value, 0.025),
      upper_95  = stats::quantile(value, 0.975),
      prob_pos  = mean(value > 0),
      prob_neg  = mean(value < 0),
      .groups   = "drop"
    ) %>%
    dplyr::mutate(
      # label if CrI excludes 0 and posterior prob is high
      flag_bias = dplyr::case_when(
        lower_95 > 0 & prob_pos >= prob_thresh ~ "positive_bias",
        upper_95 < 0 & prob_neg >= prob_thresh ~ "negative_bias",
        TRUE                                   ~ "no_clear_bias"
      )
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
  
  bias_summ %>%
    dplyr::select(item, rater, i, r,
           bias_mean, bias_sd, lower_95, upper_95,
           prob_pos, prob_neg, flag_bias) %>%
    dplyr::arrange(item, rater)
}
