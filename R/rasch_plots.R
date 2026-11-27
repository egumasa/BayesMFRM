#' Plot a fitted many-facet Rasch model
#'
#' The default interval plot shows the same main effects as [tab_bmfrm()].
#' Points are posterior means and bars are equal-tailed credible intervals.
#' The zero line is the logit reference. A level with no ratings or a known
#' identification problem has a distinct symbol. For interaction screening
#' with a practical threshold, use [plot_bias()].
#' @param x A fitted `bmfrm_fit` model.
#' @param type `"intervals"` (default), `"wright"` for the ridgeline Wright
#'   map, `"category"` for category characteristic curves, `"surface"` for
#'   coloured category surfaces over ability and severity shift, or `"pathway"`
#'   for measures against infit or outfit.
#' @param ... For `"intervals"`, supply `facet` (one or more main facets or
#'   interactions) and optionally `ci_level`. For `"wright"`, `"category"`, and `"surface"`,
#'   supply `newdata` with the facet profiles to display. For `"pathway"`,
#'   specify one main `facet`. Further arguments are passed to
#'   [plot_wright_map()], [plot_category_curves()], [plot_category_surfaces()] or [plot_pathway_map()].
#' @return Draws the plot and invisibly returns it for reuse. Surface plots are
#'   lattice `trellis` objects. Other types are ggplot objects.
#' @examples
#' \donttest{
#' if (check_cmdstan()) {
#'   source(system.file("examples", "fit-writing.R", package = "BayesMFRM"), local = TRUE)
#'   plot(fit, facet = "rater")
#'   profiles <- data.frame(item = "Content", rater = c("A", "B"))
#'   plot(fit, type = "wright", newdata = profiles)
#'   plot(fit, type = "category", newdata = profiles)
#'   plot(fit, type = "pathway", facet = "rater")
#' }
#' }
#' @export
plot.bmfrm_fit <- function(x, type = c("intervals", "wright", "category", "pathway", "surface"), ...) {
  p <- switch(match.arg(type),
              intervals = bmfrm_interval_plot(x, ...),
              wright = plot_wright_map(x, ...),
              category = plot_category_curves(x, ...),
              surface = plot_category_surfaces(x, ...),
              pathway = plot_pathway_map(x, ...))
  print(p)
  invisible(p)
}

bmfrm_interval_plot <- function(object, facet = NULL, ci_level = 0.95) {
  dat <- tab_bmfrm(object, facets = facet, ci_level = ci_level)$estimates
  label <- ifelse(is.na(dat$level2), dat$level, paste(dat$level, dat$level2, sep = " / "))
  dat$facet <- factor(dat$facet, levels = unique(dat$facet))
  ids <- seq_len(nrow(dat))
  dat$row <- factor(ids, levels = rev(ids))
  dat$status <- ifelse(dat$n_obs == 0L, "No ratings",
                       ifelse(is.na(dat$likelihood_identified), "Unchecked",
                              ifelse(dat$likelihood_identified, "Identified", "Unidentified")))
  ggplot2::ggplot(dat, ggplot2::aes(.data$mean, .data$row)) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey60", linetype = "dashed") +
    ggplot2::geom_segment(ggplot2::aes(x = .data$lower, xend = .data$upper, yend = .data$row),
                          linewidth = 0.8, colour = "#246780") +
    ggplot2::geom_point(ggplot2::aes(shape = .data$status), size = 2.7, colour = "#164755") +
    ggplot2::scale_shape_manual(values = c("Identified" = 16, "Unidentified" = 1,
                                            "No ratings" = 4, "Unchecked" = 2),
                                guide = if (all(dat$status == "Identified")) "none" else "legend") +
    ggplot2::scale_y_discrete(labels = stats::setNames(label, ids)) +
    ggplot2::facet_grid(rows = ggplot2::vars(facet), scales = "free_y", space = "free_y") +
    ggplot2::labs(title = "Facet estimates", x = "Measure (logits)", y = NULL,
                  shape = "Design status",
                  caption = sprintf("Posterior means and %g%% equal-tailed credible intervals.", 100 * ci_level)) +
    ggplot2::theme_minimal() +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   strip.text.y = ggplot2::element_text(angle = 0, face = "bold"),
                   legend.position = "bottom")
}

#' Plot conditional category characteristic curves
#'
#' @inheritParams category_curves
#' @return A ggplot with the [category_curves()] table in `$data`. Each panel
#'   is a facet profile. Lines show posterior mean response probabilities and
#'   bands show pointwise credible intervals. For a single ability, the plot
#'   shows points and vertical intervals.
#' @examples
#' \donttest{
#' if (check_cmdstan()) {
#'   source(system.file("examples", "fit-writing.R", package = "BayesMFRM"),
#'          local = TRUE)
#'   profiles <- data.frame(item = "Content", rater = c("A", "B"))
#'   plot_category_curves(fit, profiles)
#'   # Compare the height of the category 4 curve across panels at ability zero.
#' }
#' }
#' @importFrom rlang .data
#' @export
plot_category_curves <- function(object, newdata, theta = seq(-4, 4, length.out = 81),
                                 ci_level = 0.95) {
  dat <- category_curves(object, newdata, theta, ci_level)
  profile_labels <- stats::setNames(bmfrm_profile_labels(newdata), seq_len(nrow(newdata)))
  p <- ggplot2::ggplot(dat, ggplot2::aes(.data$theta, .data$probability,
                                         colour = factor(.data$category),
                                         group = .data$category))
  single <- length(unique(dat$theta)) == 1L
  if (single) {
    p <- p + ggplot2::aes(x = factor(.data$category)) +
      ggplot2::geom_linerange(ggplot2::aes(ymin = .data$lower, ymax = .data$upper)) +
      ggplot2::geom_point(size = 2) +
      ggplot2::labs(subtitle = sprintf("Ability = %g logits", dat$theta[1L]))
  } else {
    p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lower, ymax = .data$upper,
                                               fill = factor(.data$category)),
                                  colour = NA, alpha = 0.16, show.legend = FALSE) +
      ggplot2::geom_line(linewidth = 0.8)
  }
  p + ggplot2::facet_wrap(~profile, labeller = ggplot2::labeller(profile = profile_labels)) +
    ggplot2::scale_colour_viridis_d(end = 0.85, aesthetics = c("colour", "fill")) +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::labs(title = "Category characteristic curves",
                  x = if (single) "Score category" else "Ability (logits)",
                  y = "Category probability", colour = "Score", fill = "Score",
                  caption = sprintf("Posterior means and %g%% pointwise credible intervals.",
                                    100 * ci_level)) +
    ggplot2::theme_minimal()
}

#' Draw a ridgeline Wright map on the ability scale
#'
#' The top ridge describes variation among the fitted persons' posterior mean
#' abilities. Each tick is one person. Each remaining ridge is the posterior
#' distribution of one adjacent-category threshold for the indicated profile.
#' The dot marks its median and the segment its credible interval. Farther
#' right means higher ability is needed to make the two adjacent scores equally
#' likely. Each density is scaled to its own peak. A constant sample appears as
#' a tick or point. Ridge heights do not encode sample sizes.
#'
#' @inheritParams category_curves
#' @return A ggplot with the [category_thresholds()] table in `$data`.
#' @examples
#' \donttest{
#' if (check_cmdstan()) {
#'   source(system.file("examples", "fit-writing.R", package = "BayesMFRM"),
#'          local = TRUE)
#'   profiles <- data.frame(item = c("Content", "Evidence", "Structure"), rater = "A")
#'   plot_wright_map(fit, profiles)
#'   # Compare the 3 to 4 thresholds with the range covered by the writers' ticks.
#' }
#' }
#' @export
plot_wright_map <- function(object, newdata, ci_level = 0.95) {
  profiles <- bmfrm_curve_profiles(object, newdata, ci_level)
  abilities <- colMeans(profiles$params$theta)
  ticks <- data.frame(x = abilities, row = 1)
  ridge_density <- function(x, row, kind) {
    if (length(x) < 2L || diff(range(x)) == 0) return(NULL)
    d <- stats::density(x)
    data.frame(x = d$x, row = row, top = row + 0.8 * d$y / max(d$y), kind = kind)
  }
  ridges <- list(ridge_density(abilities, 1, "Persons"))
  tables <- vector("list", nrow(newdata))
  baseline <- 0
  for (i in seq_len(nrow(newdata))) {
    steps <- bmfrm_step_draws(bmfrm_profile_cumulative(profiles, i))
    k <- seq_len(ncol(steps))
    rows <- baseline - (k - 1L)
    tables[[i]] <- bmfrm_crossings(steps, i, profiles$probs, ci_level)
    tables[[i]]$row <- rows
    ridges[[i + 1L]] <- dplyr::bind_rows(lapply(k, function(j) {
      ridge_density(steps[, j], rows[j], "Thresholds")
    }))
    baseline <- min(rows) - 1.35
  }
  ridges <- dplyr::bind_rows(ridges)
  markers <- dplyr::bind_rows(tables)
  cuts <- add_curve_profile_levels(markers[names(markers) != "row"], newdata)
  fixed <- vapply(newdata, function(x) length(unique(x)) == 1L, logical(1))
  subtitle <- if (any(fixed)) bmfrm_profile_labels(newdata[1L, fixed, drop = FALSE])
  labels <- if (all(fixed)) rep("", nrow(newdata)) else bmfrm_profile_labels(newdata[, !fixed, drop = FALSE])
  cut_labels <- trimws(sprintf("%s   %d to %d", labels[cuts$profile], cuts$from, cuts$to))
  locations <- c(1, markers$row)
  p <- ggplot2::ggplot(cuts) +
    ggplot2::geom_hline(yintercept = locations, colour = "grey90", linewidth = 0.3)
  if (nrow(ridges)) {
    p <- p + ggplot2::geom_ribbon(data = ridges, ggplot2::aes(x = .data$x, ymin = .data$row,
                               ymax = .data$top, group = .data$row, fill = .data$kind), alpha = 0.45) +
      ggplot2::scale_fill_manual(values = c(Persons = "#546875", Thresholds = "#24817F"),
                                 guide = "none")
  }
  p + ggplot2::geom_segment(data = ticks, ggplot2::aes(x = .data$x, xend = .data$x,
                          y = .data$row - 0.12, yend = .data$row), colour = "#394C58") +
    ggplot2::geom_segment(data = markers, ggplot2::aes(x = .data$lower, xend = .data$upper,
                         y = .data$row, yend = .data$row), colour = "#145D61", linewidth = 0.9) +
    ggplot2::geom_point(data = markers, ggplot2::aes(.data$median, .data$row),
                        colour = "#145D61", size = 1.9) +
    ggplot2::scale_y_continuous(breaks = locations, labels = c("Person abilities", cut_labels)) +
    ggplot2::labs(title = "Wright map", subtitle = subtitle,
                  x = "Ability and adjacent-category thresholds (logits)",
                  y = NULL, caption = sprintf(paste0("Ticks: person posterior means.  ",
                     "Points and bars: threshold medians and %g%% credible intervals."), 100 * ci_level)) +
    ggplot2::theme_minimal() +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(),
                   panel.grid.minor.y = ggplot2::element_blank())
}

#' Plot facet measures against posterior predictive fit
#'
#' @param object A fitted `bmfrm_fit` that supports [facet_summary()] and
#'   [facet_fit()] with replicated fit statistics.
#' @param facet Name of a fitted main facet, e.g. `"rater"`.
#' @param statistic Realized mean square on the horizontal axis, either
#'   `"infit"` (information weighted) or `"outfit"` (unweighted).
#' @inheritParams category_curves
#' @param seed Random seed passed to [facet_fit()] for replicated ratings.
#' @return A ggplot. Its `$data` includes `level`, `n_obs`, the facet measure's
#'   posterior `mean`, `lower`, `upper`, the realized posterior mean `msq`,
#'   replicated bounds `rep_lower`, `rep_upper`, and the matched upper-tail
#'   posterior predictive probability `ppp`.
#' @details For a rater, higher on the plot means more severe. For the ability
#'   facet it means greater ability. A value to the right of one indicates more
#'   residual variation than expected. The vertical interval expresses
#'   uncertainty about the measure. The horizontal interval shows the range
#'   expected from replicated ratings of that level under the fitted model.
#'   Use the matched posterior predictive probability to investigate misfit.
#' @examples
#' \donttest{
#' if (check_cmdstan()) {
#'   source(system.file("examples", "fit-writing.R", package = "BayesMFRM"),
#'          local = TRUE)
#'   p <- plot_pathway_map(fit, "rater")
#'   print(p)
#'   p$data[, c("level", "n_obs", "msq", "rep_lower", "rep_upper", "ppp")]
#'   # Compare each observed infit with the replicated range for that rater.
#' }
#' }
#' @export
plot_pathway_map <- function(object, facet, statistic = c("infit", "outfit"),
                             ci_level = 0.95, seed = 123L) {
  statistic <- match.arg(statistic)
  probs <- bmfrm_interval_probs(ci_level)
  measures <- facet_summary(object, facet, probs = probs)
  fit <- facet_fit(object, facet, prob = ci_level, seed = seed)
  i <- match(fit$label, measures$level)
  dat <- tibble::tibble(level = fit$label, n_obs = fit$n_obs, mean = measures$mean[i],
                        lower = measures[[fmt_prob_col(probs[1L])]][i],
                        upper = measures[[fmt_prob_col(probs[2L])]][i],
                        msq = fit[[paste0(statistic, "_msq")]],
                        rep_lower = fit[[paste0(statistic, "_rep_lo")]],
                        rep_upper = fit[[paste0(statistic, "_rep_hi")]],
                        ppp = fit[[paste0(statistic, "_ppp")]])
  if (all(is.na(dat$rep_lower))) {
    stop("The pathway map needs replicated fit statistics. Save y_rep or fit a generated model.")
  }
  ggplot2::ggplot(dat, ggplot2::aes(.data$msq, .data$mean)) +
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed", colour = "grey55") +
    ggplot2::geom_segment(ggplot2::aes(x = .data$rep_lower, xend = .data$rep_upper,
                          yend = .data$mean), colour = "#A65E26", linewidth = 0.9) +
    ggplot2::geom_linerange(ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
                            colour = "#246780", linewidth = 0.9) +
    ggplot2::geom_point(size = 2.4) +
    ggplot2::geom_text(ggplot2::aes(label = .data$level), hjust = -0.45, vjust = -0.55) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.07, 0.15))) +
    ggplot2::labs(title = "Pathway map", subtitle = paste("Facet:", facet),
                  x = paste(tools::toTitleCase(statistic), "mean square"),
                  y = if (facet == spec_person_facet(object$spec)) "Ability (logits)" else
                    if (facet == "rater") "Rater severity (logits)" else
                    paste0(tools::toTitleCase(facet), " measure (logits)"),
                  caption = sprintf(paste0("Vertical: %g%% credible interval for the measure.\n",
                                    "Horizontal: %g%% reference interval for replicated fit."),
                                    100 * ci_level, 100 * ci_level)) +
    ggplot2::theme_minimal()
}
