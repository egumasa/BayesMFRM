#' Compare category probabilities at the same ability
#'
#' For each score and ability, subtract the second profile's probability from
#' the first profile's probability within each posterior draw. Positive values
#' mean that the first profile is more likely to yield that score.
#'
#' @inheritParams category_curves
#' @param newdata Exactly two rows with the conditioning facets required by
#'   [category_curves()]. Row order sets the direction (first minus second).
#'   Compare scores with the same meaning in both profiles. The function checks
#'   that the declared numbers of categories agree. Rubric meanings must be
#'   established by the user.
#' @param margin_pp Optional practical equivalence margin in percentage points,
#'   greater than 0 and at most 100. For example, 5 asks how often the absolute
#'   difference is at most five percentage points. Choose it for the intended
#'   scoring decision. The default leaves practical equivalence unassessed.
#' @details For a rater comparison, hold the criterion and other facets fixed.
#'   Differences include rater severity and any fitted interactions. The counts
#'   are ratings used in the fit matching all facets of each profile, across
#'   abilities and scores. A zero count marks a model-based extrapolation.
#'
#'   The rating scale model assumes common step parameters, while a partial
#'   credit model allows steps to vary by threshold-facet level. Common steps
#'   can coexist with different score probabilities at the same ability.
#'   To examine practical homogeneity, choose a meaningful margin and read
#'   the posterior probability within that margin for a specified score and
#'   ability. A wide interval spanning zero can also include sizeable differences.
#'   An assessment of homogeneity across an entire rubric or ability range
#'   requires a joint criterion across those scores and abilities.
#' @return `category_contrast()` returns a tibble. `difference_pp` is the
#'   posterior mean of the paired probability differences multiplied by 100.
#'   `lower_pp` and `upper_pp` bound their equal-tailed interval at `ci_level`.
#'   Each row identifies `theta`, `category`, `profile_1`, `profile_2`, their
#'   rating counts (`n_1`, `n_2`), and whether either is `extrapolated`.
#'   `prob_within_margin` is the fraction of draws with an absolute difference
#'   at most `margin_pp`. Both fields are missing when no margin is supplied.
#'   `plot_category_contrast()` returns a ggplot with this table in `$data`.
#' @seealso [category_curves()], [summarise_bias()]
#' @examples
#' \donttest{
#' if (check_cmdstan()) {
#'   source(system.file("examples", "fit-writing.R", package = "BayesMFRM"),
#'          local = TRUE)
#'   ab <- data.frame(item = "Content", rater = c("A", "B"))
#'   d <- category_contrast(fit, ab, theta = 0, margin_pp = 5)
#'   d[d$category %in% 3:4, c("category", "difference_pp", "lower_pp", "upper_pp",
#'                        "n_1", "n_2", "prob_within_margin")]
#'   # At zero logits, a 4 is about 56 percentage points more likely with A.
#'   # For a 3 the interval spans zero, and only 16.7% of draws fall within +/-5 points.
#'   plot_category_contrast(fit, ab, margin_pp = 5)
#' }
#' }
#' @export
category_contrast <- function(object, newdata, theta = seq(-4, 4, length.out = 81),
                              ci_level = 0.95, margin_pp = NULL) {
  if (!is.data.frame(newdata) || nrow(newdata) != 2L) {
    stop("Supply exactly two profile rows in newdata (first minus second)")
  }
  if (!is.null(margin_pp) && (!is.numeric(margin_pp) || !is.null(dim(margin_pp)) ||
      length(margin_pp) != 1L || !is.finite(margin_pp) || margin_pp <= 0 || margin_pp > 100)) {
    stop("margin_pp must be a single number in (0, 100] percentage points")
  }
  theta <- bmfrm_curve_grid(theta)
  profiles <- bmfrm_curve_profiles(object, newdata, ci_level)
  first <- bmfrm_profile_cumulative(profiles, 1L)
  second <- bmfrm_profile_cumulative(profiles, 2L)
  K <- ncol(first)
  if (K != ncol(second)) {
    stop("Profiles must use the same number of rubric categories with comparable meanings (",
         K, " and ", ncol(second), " supplied)")
  }

  observed <- object$data_info$stan_data
  counts <- vapply(1:2, function(i) {
    keep <- rep(TRUE, observed$N)
    for (f in names(profiles$ids)) {
      keep <- keep & observed[[paste0(f, "_id")]] == profiles$ids[[f]][i]
    }
    sum(keep)
  }, integer(1))
  labels <- bmfrm_profile_labels(newdata)
  size <- default_bmfrm_chunk_size(2 * nrow(first) * K, length(theta))
  chunks <- split(theta, ceiling(seq_along(theta) / size))
  result <- dplyr::bind_rows(lapply(chunks, function(grid) {
    delta <- 100 * (bmfrm_curve_draws(first, grid) - bmfrm_curve_draws(second, grid))
    bounds <- apply(delta, 2L, stats::quantile, probs = profiles$probs, names = FALSE)
    within <- if (is.null(margin_pp)) NA_real_ else colMeans(abs(delta) <= margin_pp)
    tibble::tibble(theta = rep(grid, K), category = rep(seq_len(K), each = length(grid)),
                   difference_pp = colMeans(delta), lower_pp = bounds[1L, ], upper_pp = bounds[2L, ],
                   ci_level = ci_level, margin_pp = margin_pp %||% NA_real_, prob_within_margin = within,
                   profile_1 = labels[1L], profile_2 = labels[2L], n_1 = counts[1L], n_2 = counts[2L],
                   extrapolated = any(counts == 0L))
  }))
  result[order(result$category, result$theta), ]
}

#' @rdname category_contrast
#' @export
plot_category_contrast <- function(object, newdata, theta = seq(-4, 4, length.out = 81),
                                   ci_level = 0.95, margin_pp = NULL) {
  dat <- category_contrast(object, newdata, theta, ci_level, margin_pp)
  counts <- c(dat$n_1[1L], dat$n_2[1L])
  label <- sprintf("%s (n = %d%s)", c(dat$profile_1[1L], dat$profile_2[1L]), counts,
                   ifelse(counts == 0L, ", model extrapolation", ""))
  single <- length(unique(dat$theta)) == 1L
  caption <- sprintf("Posterior means and %g%% pointwise credible intervals.", 100 * ci_level)
  p <- ggplot2::ggplot(dat, ggplot2::aes(.data$theta, .data$difference_pp))
  if (!is.null(margin_pp)) {
    p <- p + ggplot2::geom_hline(yintercept = c(-margin_pp, margin_pp),
                                 linetype = "dotted", colour = "grey45")
    caption <- paste(caption, sprintf("Dotted lines mark the chosen margin (+/-%g points).", margin_pp))
  }
  p <- p + ggplot2::geom_hline(yintercept = 0, colour = "grey40")
  if (single) {
    p <- p + ggplot2::aes(x = factor(.data$category)) +
      ggplot2::geom_linerange(ggplot2::aes(ymin = .data$lower_pp, ymax = .data$upper_pp),
                              colour = "#205E72", linewidth = 0.8) +
      ggplot2::geom_point(colour = "#205E72", size = 2)
  } else {
    p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lower_pp, ymax = .data$upper_pp),
                                  fill = "#6BA8AE", alpha = 0.35) +
      ggplot2::geom_line(colour = "#205E72", linewidth = 0.8) +
      ggplot2::facet_wrap(~category, labeller = ggplot2::labeller(category = function(k) paste("Score", k)))
  }
  p + ggplot2::labs(title = "Difference in category probabilities",
                    subtitle = paste(strwrap(paste(label, collapse = "  minus  "), width = 90),
                                     collapse = "\n"),
                    x = if (single) sprintf("Score category at ability %g logits", dat$theta[1L]) else
                      "Ability (logits)",
                    y = "Probability difference (percentage points)", caption = caption) +
    ggplot2::theme_minimal()
}
