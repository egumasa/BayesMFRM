#' Posterior means and credible intervals for model parameters
#'
#' `coef()` reports posterior means. `confint()` reports marginal equal-tailed
#' Bayesian credible intervals from the same retained draws. Both use the
#' original parameter scale. For example, positive rater main effects indicate
#' greater severity; `tau` represents centered step deviations. Use
#' [category_thresholds()] for adjacent-category crossings on the ability scale.
#'
#' @param object A fitted `bmfrm_fit` model.
#' @param parm An optional vector of fitted facet names (e.g. `"rater"`),
#'   two-way interactions in either order (e.g. `"item:rater"`), saved Stan
#'   variables (e.g. `"sigma_theta"`), or indexed names (e.g. `"rater[3]"`).
#'   Positive integer positions select entries in the order of `coef(object)`.
#'   The default includes the main facet measures, centered steps (`tau`), and
#'   fitted interactions. For a custom Stan model with its own parameter names,
#'   supply those names explicitly.
#' @param level Probability mass of the equal-tailed credible interval.
#' @param ... Reserved for other methods.
#' @return `coef()` returns a named numeric vector. `confint()` returns a
#'   numeric matrix with one row per parameter and columns for the lower
#'   and upper posterior quantiles. Row names use Stan indices in fitted order,
#'   including for reversed interaction names. [get_facet_labels()] gives their
#'   original level labels. For scores on the rating scale use
#'   [predict.bmfrm_fit()]; for a difference between levels use [facet_contrasts()].
#' @examples
#' \donttest{
#' if (check_cmdstan()) {
#'   source(system.file("examples", "fit-writing.R", package = "BayesMFRM"), local = TRUE)
#'   get_facet_labels(fit$data_info, "rater")
#'   coef(fit, parm = "rater")
#'   confint(fit, parm = "rater[3]")  # Rater C: about 0.07 to 0.57 logits
#'   confint(fit, parm = "tau", level = 0.90)
#' }
#' }
#' @name bmfrm_coefficients
NULL

#' @rdname bmfrm_coefficients
#' @importFrom stats coef
#' @export
coef.bmfrm_fit <- function(object, parm = NULL, ...) {
  rlang::check_dots_empty()
  s <- bmfrm_parameter_summary(object, parm, "mean")
  stats::setNames(s$mean, s$variable)
}

#' @rdname bmfrm_coefficients
#' @importFrom stats confint
#' @export
confint.bmfrm_fit <- function(object, parm = NULL, level = 0.95, ...) {
  rlang::check_dots_empty()
  probs <- bmfrm_interval_probs(level)
  s <- bmfrm_parameter_summary(object, parm, ~stats::quantile(.x, probs, names = FALSE))
  bounds <- as.matrix(s[, -1L, drop = FALSE])
  dimnames(bounds) <- list(s$variable, paste0(format(100 * probs, trim = TRUE, digits = 15L), "%"))
  bounds
}

bmfrm_parameter_summary <- function(object, parm, ...) {
  facets <- object$spec$facets_main
  interactions <- object$spec$facets_bias
  aliases <- c(stats::setNames(vapply(facets, function(f) facet_parameter_name(object, f), character(1)), facets),
               stats::setNames(gsub(":", "", interactions), interactions))
  defaults <- unique(c(unname(aliases[facets]), "tau", unname(aliases[interactions])))
  if (is.null(parm) || is.numeric(parm)) {
    variables <- defaults
    if (isTRUE(object$spec$custom_stan)) {
      variables <- Filter(function(v) fit_has_variable(object$fit, v), variables)
      if (!length(variables)) stop("Specify a saved Stan variable in parm for this custom model.")
    }
    if (is.numeric(parm)) {
      if (!length(parm) || !is.null(dim(parm)) || any(!is.finite(parm)) ||
          any(parm != floor(parm)) || any(parm < 1)) stop("parm positions must be positive integers")
      # CmdStan's expanded names allow selection before summarising any draws.
      saved <- if (is.function(object$fit$metadata)) object$fit$metadata()$model_params
      if (length(saved) && any(grepl("[", saved, fixed = TRUE))) {
        blocks <- split(saved, sub("\\[.*$", "", saved))
        index <- unlist(blocks[variables], use.names = FALSE)
        if (any(parm > length(index))) stop("parm position exceeds the number of coefficients")
        selected <- index[parm]
        s <- object$fit$summary(unique(selected), ...)
        return(s[match(selected, s$variable), , drop = FALSE])
      }
    }
  } else {
    if (!is.character(parm) || !length(parm) || !is.null(dim(parm)) ||
        anyNA(parm) || any(!nzchar(parm))) stop("parm must contain parameter names or positive positions")
    variables <- match_interaction_names(parm, interactions)
    mapped <- variables %in% names(aliases)
    variables[mapped] <- aliases[variables[mapped]]
  }
  s <- object$fit$summary(variables, ...)
  if (is.numeric(parm)) {
    if (any(parm > nrow(s))) stop("parm position exceeds the number of coefficients")
    s <- s[parm, , drop = FALSE]
  }
  s
}
