#' Exact agreement between raters on the same performances
#'
#' Compare the observed proportion of identical scores with its model expectation
#' and a reference interval from replicated ratings on the observed assignments.
#'
#' @param object A fitted `bmfrm_fit` supporting category probabilities in R.
#' @param facet Name of the rater main facet. All other main facets, including
#'   examinee, criterion and any task facet, are held fixed when matching.
#' @param occasion Optional additional column name(s) in the fitted data that
#'   identify a performance or scoring occasion. Each rater must have at most
#'   one observed score per matching combination.
#' @param ci_level Probability of the central posterior predictive interval.
#' @param seed Seed for replicated scores. A fixed seed preserves the caller's
#'   random-number state. Use `NULL` to draw from the current random stream.
#' @param x A result table from `rater_agreement()` or a subset of its rows.
#' @details Each unordered pair of different raters scoring the same performance
#'   on the same rubric contributes one opportunity. Three raters contribute
#'   three opportunities. A score is matched on its numeric category. Paired
#'   ratings must have equal numbers of categories with comparable meanings.
#'   Missing ratings and unmatched performances contribute no opportunities.
#'
#'   Given a posterior draw and conditional independence of the two ratings,
#'   the chance of exact agreement is the sum over categories of the products
#'   of the two probabilities. The expectation averages these chances over
#'   opportunities and draws. To obtain the predictive interval, one replicate
#'   score is generated per rating per draw and reused in every pair involving
#'   that rating. The observed design, fitted abilities, and any interactions
#'   enter both calculations.
#' @return A tibble with one row per pair with shared ratings and an
#'   opportunity-weighted overall row. `scope`, `rater1`, `rater2` and `facet`
#'   identify the comparison. `n_opportunities` is the denominator and
#'   `n_exact` is the observed count. `observed` and `expected` are proportions
#'   from 0 to 1. `rep_lower` and `rep_upper` are empirical quantiles of
#'   replicated proportions at `ci_level` (pointwise, without interpolation).
#'   With no shared ratings, the overall row gives zero counts and missing
#'   rates and bounds. `plot_rater_agreement()` returns a ggplot.
#' @references [FACETS Table 7 agreement statistics](https://www.winsteps.com/facetman64/table7agreementstatistics.htm).
#' @seealso [facet_separation()], [facet_contrasts()], [facet_fit()]
#' @examples
#' data(bmfrm_example)
#' bmfrm_example$agreement
#' plot_rater_agreement(bmfrm_example$agreement)
#' # A and B gave identical scores on 8 of 46 shared writer-criterion combinations.
#' \donttest{
#' if (check_cmdstan()) {
#'   source(system.file("examples", "fit-writing.R", package = "BayesMFRM"), local = TRUE)
#'   rater_agreement(fit, "rater", seed = 123)
#' }
#' }
#' @export
rater_agreement <- function(object, facet = "rater", occasion = NULL,
                            ci_level = 0.95, seed = 123L) {
  if (!inherits(object, "bmfrm_fit")) stop("object must be a bmfrm_fit")
  spec <- object$spec
  if (length(facet) != 1L || !is.character(facet) ||
      !facet %in% setdiff(spec$facets_main, spec_person_facet(spec))) {
    stop("facet must name a non-ability main facet identifying the raters")
  }
  probs <- bmfrm_interval_probs(ci_level)
  df <- object$data_info$data_clean
  if (!is.null(occasion) && (!is.character(occasion) ||
      !all(occasion %in% setdiff(names(df), c(spec$facets_main, spec$response))))) {
    stop("occasion must name additional performance identifiers in the fitted data")
  }
  keys <- c(setdiff(spec$facets_main, facet), occasion)
  if (anyNA(df[keys])) stop("Matching columns must have nonmissing values")
  if (anyDuplicated(df[c(keys, facet)])) {
    stop("Repeated ratings by one rater for the same matching combination. ",
         "Supply an occasion column identifying the shared performance.")
  }
  groups <- dplyr::group_rows(dplyr::group_by(df, dplyr::across(dplyr::all_of(keys))))
  groups <- groups[lengths(groups) >= 2L]
  overall <- tibble::tibble(facet = facet, scope = "overall", rater1 = NA_character_,
                            rater2 = NA_character_, n_opportunities = 0L, n_exact = 0L)
  if (!length(groups)) {
    return(dplyr::mutate(overall, observed = NA_real_, expected = NA_real_,
                         rep_lower = NA_real_, rep_upper = NA_real_, ci_level = ci_level))
  }
  if (!can_evaluate_in_r(object)) stop("Agreement predictions require a verified R prediction model")

  labels <- object$data_info$facet_labels
  ids <- map_eval_indices(df, spec, labels)
  rater <- ids[[facet]]
  groups <- lapply(groups, function(rows) rows[order(rater[rows])])
  pairs <- lapply(groups, utils::combn, m = 2L)
  all_pairs <- do.call(cbind, pairs)
  a <- all_pairs[1L, ]
  b <- all_pairs[2L, ]
  by_pair <- dplyr::group_by(tibble::tibble(i = rater[a], j = rater[b]), .data$i, .data$j)
  pair_id <- dplyr::group_indices(by_pair)
  raters <- dplyr::group_keys(by_pair)
  n_pairs <- nrow(raters)
  scores <- df[[spec$response]]
  rater_labels <- labels[[facet]]
  tab <- tibble::tibble(facet = facet, scope = "pair", rater1 = rater_labels[raters$i],
                        rater2 = rater_labels[raters$j],
                        n_opportunities = tabulate(pair_id, n_pairs),
                        n_exact = tabulate(pair_id[scores[a] == scores[b]], n_pairs))
  overall$n_opportunities <- length(a)
  overall$n_exact <- sum(tab$n_exact)
  tab <- dplyr::bind_rows(tab, overall)

  params <- extract_bmfrm_param_draws(object)
  K <- rep.int(params$K, nrow(df))
  if (identical(params$tau_info$mode, "pcm_ragged")) {
    K <- params$tau_info$K_level[ids[[spec$threshold_facet]]]
    if (any(K[a] != K[b])) stop("Paired ratings must have the same number of rubric categories")
  }
  expected <- numeric(n_pairs)
  replicated <- matrix(0L, nrow = params$D, ncol = n_pairs)
  offset <- 0L
  with_bmfrm_seed(seed, {
    for (g in seq_along(groups)) {
      rows <- groups[[g]]
      ch <- eval_bmfrm_chunk(params, ids, rows, compute_mu = FALSE,
                             compute_probs = TRUE, compute_y_rep = TRUE)
      # Each column of the pair matrix points to two of the shared rating rows.
      within_group <- matrix(match(pairs[[g]], rows), nrow = 2L)
      categories <- seq_len(K[rows[1L]])
      for (j in seq_len(ncol(within_group))) {
        p <- pair_id[offset + j]
        i1 <- within_group[1L, j]
        i2 <- within_group[2L, j]
        joint <- ch$probs[, i1, categories] * ch$probs[, i2, categories]
        expected[p] <- expected[p] + sum(joint) / params$D
        replicated[, p] <- replicated[, p] + (ch$y_rep[, i1] == ch$y_rep[, i2])
      }
      offset <- offset + ncol(within_group)
    }
  })
  replicated <- cbind(replicated, rowSums(replicated))
  bounds <- apply(replicated, 2L, stats::quantile, probs = probs, type = 1, names = FALSE)
  tab$observed <- tab$n_exact / tab$n_opportunities
  tab$expected <- c(expected, sum(expected)) / tab$n_opportunities
  tab$rep_lower <- bounds[1L, ] / tab$n_opportunities
  tab$rep_upper <- bounds[2L, ] / tab$n_opportunities
  tab$ci_level <- ci_level
  tab
}

#' @rdname rater_agreement
#' @export
plot_rater_agreement <- function(x) {
  required <- c("scope", "rater1", "rater2", "n_opportunities", "observed", "expected",
                "rep_lower", "rep_upper", "ci_level")
  if (!is.data.frame(x) || !nrow(x) || !all(required %in% names(x)) ||
      length(unique(x$ci_level)) != 1L) stop("Supply rows from one rater_agreement() result")
  if (all(x$n_opportunities == 0)) stop("No jointly rated performances to plot")
  x$position <- rev(seq_len(nrow(x)))
  pair <- ifelse(x$scope == "overall", "Overall", paste(x$rater1, x$rater2, sep = " / "))
  ggplot2::ggplot(x, ggplot2::aes(y = .data$position)) +
    ggplot2::geom_segment(ggplot2::aes(x = .data$rep_lower, xend = .data$rep_upper,
                                       yend = .data$position), colour = "#7099A6", linewidth = 1.2) +
    ggplot2::geom_point(ggplot2::aes(x = .data$observed, shape = "Observed"), size = 2.7) +
    ggplot2::geom_point(ggplot2::aes(x = .data$expected, shape = "Model expectation"),
                        size = 4, colour = "#A35021") +
    ggplot2::scale_shape_manual(NULL, values = c("Observed" = 16, "Model expectation" = 5),
                                breaks = c("Observed", "Model expectation")) +
    ggplot2::scale_y_continuous(breaks = x$position, labels = sprintf("%s  (n = %d)", pair, x$n_opportunities)) +
    ggplot2::scale_x_continuous(labels = function(x) paste0(round(100 * x), "%")) +
    ggplot2::labs(title = "Exact agreement on shared ratings",
                  subtitle = sprintf("Lines show %g%% posterior predictive intervals", 100 * x$ci_level[1]),
                  x = "Proportion of identical scores", y = NULL) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                   legend.position = "bottom")
}
