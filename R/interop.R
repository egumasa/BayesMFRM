#' Use a fitted model with posterior and tidybayes
#'
#' Extract joint draws from the underlying CmdStan fit with their original
#' chain and iteration identifiers. Use the facet labels to interpret indices.
#'
#' @param x A fitted `bmfrm_fit`.
#' @param ... Arguments to the CmdStan fit's `$draws()` method. For example,
#'   pass `variables = "rater"` to read only rater severity draws.
#' @return A posterior draws object, by default a draws array. The
#'   `.draw` column of its data frame representation identifies the joint
#'   draw used by [as_mirt()].
#' @examples
#' \donttest{
#' if (check_cmdstan()) {
#'   source(system.file("examples", "fit-writing.R", package = "BayesMFRM"),
#'          local = TRUE)
#'   draws <- posterior::as_draws(fit, variables = "rater")
#'   posterior::summarise_draws(draws)
#'   if (requireNamespace("tidybayes", quietly = TRUE)) {
#'     by_rater <- tidybayes::spread_draws(draws, rater[r])
#'     by_rater$label <- get_facet_labels(fit$data_info, "rater")[by_rater$r]
#'     head(by_rater) # One severity estimate per draw and rater
#'   }
#' }
#' }
#' @exportS3Method posterior::as_draws bmfrm_fit
as_draws.bmfrm_fit <- function(x, ...) x$fit$draws(...)

#' Export a fixed rating plan to mirt
#'
#' Use a single joint posterior draw to evaluate a planned set of ratings on
#' the fitted ability scale with the optional mirt package. Each row becomes
#' one unit-slope generalized partial credit item, including its fitted facet
#' effects, interactions and rubric steps.
#'
#' @inheritParams category_curves
#' @param newdata A data frame as in [test_information()], in the desired
#'   item order. Repeat a row to schedule a second independent rating under
#'   the same conditions.
#' @param draw One integer identifying a retained posterior draw (the
#'   `.draw` in `posterior::as_draws_df(object)`).
#' @return A mirt `SingleGroupClass` object with fixed item parameters.
#'   Item `rating_1` corresponds to input row 1. In an item's probability
#'   trace, column position k corresponds to rubric score k. The response
#'   matrix in this object is a synthetic template that declares each rubric's
#'   categories. Use the original fitted model for observed-data fit statistics.
#' @details For adjacent thresholds b1, b2, ... in the selected draw, the
#'   mirt category intercepts are 0, -b1, -(b1 + b2), ... . Conditional test
#'   information assumes local independence. Use [test_information()] to
#'   summarize the information across joint draws, including credible limits.
#' @examples
#' \donttest{
#' if (check_cmdstan() && requireNamespace("mirt", quietly = TRUE)) {
#'   source(system.file("examples", "fit-writing.R", package = "BayesMFRM"),
#'          local = TRUE)
#'   plan <- expand.grid(item = c("Content", "Evidence"), rater = c("A", "B"))
#'   bank <- as_mirt(fit, plan, draw = 1L)
#'   grid <- matrix(seq(-3, 3, by = 0.25), ncol = 1)
#'   info <- mirt::testinfo(bank, Theta = grid)
#'   plot(grid[, 1], info, type = "l", xlab = "Ability (logits)",
#'        ylab = "Conditional test information")
#'   # Each value describes the four planned ratings at that ability in draw 1.
#'   mirt::probtrace(mirt::extract.item(bank, 1L), Theta = matrix(0))
#'   # Probabilities of scores 1 through 4 for Content scored by A at ability 0.
#' }
#' }
#' @export
as_mirt <- function(object, newdata, draw) {
  if (!requireNamespace("mirt", quietly = TRUE)) stop("Install the mirt package to use as_mirt()")
  if (missing(draw) || length(draw) != 1L || !is.null(dim(draw))) {
    stop("Supply one posterior draw number")
  }
  profiles <- bmfrm_curve_profiles(object, newdata, draw_idx = draw)
  intercepts <- lapply(seq_len(nrow(newdata)), function(i) {
    -as.numeric(bmfrm_profile_cumulative(profiles, i))
  })
  K <- lengths(intercepts)
  scores <- seq_len(max(K)) - 1L
  template <- vapply(K, function(k) 1L + scores %% k, integer(length(scores)))
  colnames(template) <- paste0("rating_", seq_along(K))
  pars <- mirt::mirt(template, 1L, "gpcm", pars = "values", verbose = FALSE)
  pars$value[pars$name == "a1"] <- 1
  for (i in seq_along(K)) {
    rows <- which(pars$item == colnames(template)[i])
    ds <- paste0("d", seq_len(K[i]) - 1L)
    pars$value[rows[match(ds, pars$name[rows])]] <- intercepts[[i]]
  }
  pars$est[] <- FALSE
  mirt::mirt(template, 1L, "gpcm", pars = pars, TOL = NA, verbose = FALSE)
}
