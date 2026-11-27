facet_parameter_name <- function(object, facet) {
  if (!is.character(facet) || length(facet) != 1L || !facet %in% object$spec$facets_main) {
    stop("Facet '", paste(facet, collapse = ", "), "' not found in model. Available facets: ",
         paste(object$spec$facets_main, collapse = ", "))
  }
  mapped <- get_eval_contract(object)$param_map
  if (facet %in% names(mapped)) return(unname(mapped[[facet]]))
  if (identical(facet, spec_person_facet(object$spec))) "theta" else facet
}

facet_observation_counts <- function(object, facet) {
  ids <- object$data_info$stan_data[[paste0(facet, "_id")]]
  if (is.null(ids)) stop("Observed rating indices are missing for facet '", facet, "'")
  tabulate(ids, length(object$data_info$facet_labels[[facet]]))
}

#' How clearly are the levels of a facet separated?
#'
#' Compares the spread of main-effect estimates across observed levels with
#' their typical uncertainty. Choose a facet to distinguish writer ability,
#' criterion difficulty or rater severity. High rater reliability here means
#' severity differences are well resolved relative to marginal uncertainty.
#' Inter-rater agreement concerns the scores raters assign to the same work.
#'
#' @param object A fitted `bmfrm_fit` (EAP method) or a data frame with one row
#'   per level. Supply `level`, `mean` and `sd` for EAP (posterior mean and
#'   posterior SD), or `level`, `measure` and `se` for conventional Rasch
#'   estimates and their standard errors. An optional `facet` column selects
#'   matching rows. Fits exclude levels with no observed rating. For a table,
#'   include only eligible rows or supply `n_obs` to exclude zero counts.
#' @param facet Name of one main facet.
#' @param method `"eap"` (default) or `"rasch"`. The Rasch option takes
#'   externally calibrated measures and SEs. Match the intended handling of
#'   extreme scores and choice of SE in the input.
#' @details Each included level has equal weight. Write V for the sample
#'   variance of estimates (denominator n - 1), and U for the average squared
#'   posterior SD or SE. EAP uses B = V, as in
#'   `mirt::empirical_rxx(..., T_as_X = FALSE)`. Rasch uses B = max(0, V - U),
#'   the sample-variance version of the conventional error correction.
#'   Separation G = sqrt(B/U), reliability = B/(B+U), and
#'   strata H = (4G+1)/3. For EAP, G and H are descriptive transformations
#'   of the empirical reliability. H is fractional and can exceed the number
#'   of observed levels. It does not assign levels to groups.
#'
#'   These plug-in summaries depend on the sampled levels, design and, for a
#'   Bayesian fit, the priors. Sparse or disconnected assignments can give
#'   prior-sensitive results. Use [facet_contrasts()] for uncertainty about
#'   individual differences, especially with few raters.
#' @return One row giving `facet`, `method`, the included `n_levels` and
#'   excluded `n_excluded`, `estimate_sd` (sqrt(V)), `rms_uncertainty`
#'   (sqrt(U)), `signal_sd` (sqrt(B)), `separation`, `strata` and
#'   `reliability`. At least two eligible levels are required. At zero
#'   uncertainty with positive spread, separation and strata are infinite and
#'   reliability is one. With zero spread and zero uncertainty all three
#'   indices are undefined (`NA`).
#' @references Wright, B. D. & Masters, G. N. (2002).
#'   [Number of person or item strata](https://www.rasch.org/rmt/rmt163f.htm).
#'   See also the [FACETS Table 7 manual](https://www.winsteps.com/facetman64/table7summarystatistics.htm)
#'   and [mirt empirical reliability](https://philchalmers.github.io/mirt/reference/empirical_rxx.html).
#' @seealso [rater_agreement()] for identical scores on shared performances,
#'   [facet_fit()], [facet_summary()], [plot_pathway_map()]
#' @examples
#' data(bmfrm_example)
#' facet_separation(bmfrm_example$rater, "rater")
#' # Reliability is about 0.83 for severity differences among these three raters.
#' calibration <- data.frame(level = c("A", "B", "C"),
#'                           measure = c(-1, 0, 1), se = c(0.4, 0.5, 0.6))
#' facet_separation(calibration, "rater", method = "rasch")
#' @export
facet_separation <- function(object, facet, method = c("eap", "rasch")) {
  method <- match.arg(method)
  if (!is.character(facet) || length(facet) != 1L || is.na(facet) || !nzchar(facet)) {
    stop("facet must name one main facet")
  }
  if (inherits(object, "bmfrm_fit")) {
    if (method == "rasch") stop("Supply a table of conventional Rasch measures and SEs for method = 'rasch'")
    tab <- facet_summary(object, facet, probs = NULL)
    tab$n_obs <- facet_observation_counts(object, facet)
  } else if (is.data.frame(object)) {
    tab <- object
    if ("facet" %in% names(tab)) tab <- tab[which(tab$facet == facet), , drop = FALSE]
  } else {
    stop("object must be a bmfrm_fit or a data frame")
  }
  cols <- if (method == "eap") c("mean", "sd") else c("measure", "se")
  if (!all(c("level", cols) %in% names(tab))) {
    stop("Required columns for ", method, ": level, ", paste(cols, collapse = ", "))
  }
  if (anyNA(tab$level) || anyDuplicated(tab$level)) stop("level labels must be unique and nonmissing")
  n_excluded <- 0L
  if ("n_obs" %in% names(tab)) {
    n <- tab$n_obs
    if (!is.numeric(n) || any(!is.finite(n) | n < 0 | n != floor(n))) {
      stop("n_obs must contain nonnegative integer counts")
    }
    n_excluded <- sum(n == 0)
    tab <- tab[n > 0, , drop = FALSE]
  }
  if (nrow(tab) < 2L) stop("At least two observed levels are required for facet '", facet, "'")
  m <- tab[[cols[1L]]]
  s <- tab[[cols[2L]]]
  if (!is.numeric(m) || !is.numeric(s) || any(!is.finite(m)) || any(!is.finite(s) | s < 0)) {
    stop("Estimates must be finite and uncertainty values must be finite and nonnegative")
  }
  v <- stats::var(m)
  u <- mean(s^2)
  b <- if (method == "rasch") max(0, v - u) else v
  g <- if (b + u == 0) NA_real_ else sqrt(b / u)
  tibble::tibble(facet = facet, method = method, n_levels = nrow(tab), n_excluded = n_excluded,
                 estimate_sd = sqrt(v), rms_uncertainty = sqrt(u), signal_sd = sqrt(b),
                 separation = g, strata = (4 * g + 1) / 3,
                 reliability = if (b + u == 0) NA_real_ else b / (b + u))
}

#' Compare selected facet levels within posterior draws
#'
#' Computes main-effect differences for all pairs of the explicitly selected
#' levels. Within-draw subtraction preserves the posterior dependence between
#' estimates, including dependence induced by centering.
#' @param object A fitted `bmfrm_fit`.
#' @param facet Name of one main facet.
#' @param levels Character vector of two or more distinct observed labels.
#'   For each pair the later label is the comparison and the earlier one is
#'   the reference. Pass just two labels for a single comparison.
#' @param ci_level Probability in (0, 1) for each equal-tailed credible interval.
#' @return One row per pair with the `facet`, `reference`, `comparison`,
#'   posterior mean `difference` (comparison minus reference, in logits),
#'   `lower`, `upper`, `ci_level`, and posterior probability `prob_positive`
#'   of a strictly positive difference. The intervals are pointwise.
#'   `likelihood_identified` reports whether the design check identifies the
#'   contrast (`NA` if unchecked). A structurally unidentified contrast can
#'   be sensitive to priors even when both levels have ratings.
#' @details For rater severity a positive difference means the comparison rater
#'   tends to score more strictly. If the model includes interactions, a
#'   specific criterion can show a different pattern. Use [category_contrast()]
#'   for a comparison of predicted scores at specified profiles and abilities.
#' @seealso [facet_separation()], [facet_summary()]
#' @examples
#' \donttest{
#' if (check_cmdstan()) {
#'   source(system.file("examples", "fit-writing.R", package = "BayesMFRM"), local = TRUE)
#'   facet_contrasts(fit, "rater", c("A", "B", "C"))
#'   # An interval spanning zero leaves the direction of that difference uncertain.
#' }
#' }
#' @export
facet_contrasts <- function(object, facet, levels, ci_level = 0.95) {
  if (!inherits(object, "bmfrm_fit")) stop("object must be a bmfrm_fit")
  param <- facet_parameter_name(object, facet)
  probs <- bmfrm_interval_probs(ci_level)
  if (!is.character(levels) || length(levels) < 2L || anyNA(levels) || anyDuplicated(levels)) {
    stop("levels must contain at least two distinct, nonmissing labels")
  }
  ix <- match(levels, object$data_info$facet_labels[[facet]])
  if (anyNA(ix)) stop("Unknown level for facet '", facet, "': ", paste(levels[is.na(ix)], collapse = ", "))
  if (any(facet_observation_counts(object, facet)[ix] == 0L)) stop("Each selected level must have observed ratings")
  vars <- sprintf("%s[%d]", param, ix)
  draws <- object$fit$draws(vars, format = "matrix")[, vars, drop = FALSE]
  pairs <- utils::combn(seq_along(ix), 2L)
  ref <- pairs[1L, ]
  cmp <- pairs[2L, ]
  stats <- vapply(seq_along(ref), function(i) {
    d <- draws[, cmp[i]] - draws[, ref[i]]
    c(mean(d), stats::quantile(d, probs, names = FALSE), mean(d > 0))
  }, numeric(4L))

  dc <- object$design
  null <- dc$effect_null_space
  identified <- rep(NA, length(ref))
  rows <- match(paste(facet, levels, sep = ":"), rownames(null))
  if (is.matrix(null) && !anyNA(rows) && length(dc$tolerance) == 1L) {
    delta <- null[rows[cmp], , drop = FALSE] - null[rows[ref], , drop = FALSE]
    identified <- unname(rowSums(delta^2)) < dc$tolerance^2
  }
  tibble::tibble(facet = facet, reference = levels[ref], comparison = levels[cmp],
                 difference = stats[1L, ], lower = stats[2L, ], upper = stats[3L, ],
                 ci_level = ci_level, prob_positive = stats[4L, ], likelihood_identified = identified)
}
