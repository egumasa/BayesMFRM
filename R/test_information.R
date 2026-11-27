#' Conditional test information for a planned set of ratings
#'
#' At each fixed ability, sum the Fisher information contributed by the
#' ratings scheduled for one examinee. For a unit-slope Rasch rating, the
#' information is the conditional variance of its rubric score. The sum assumes
#' the ratings are independent given ability and the fitted model parameters.
#'
#' @inheritParams category_curves
#' @param newdata Data frame with the conditioning facets described in
#'   [category_curves()]. Supply one row per planned rating of the same
#'   examinee. Each row uses its own fitted rubric. Repeating a row schedules
#'   another conditionally independent rating under those conditions.
#' @details The function conditions on the selected fitted levels. It computes
#'   the information sum separately in each joint posterior draw, then
#'   summarizes the sums. The intervals reflect calibration uncertainty.
#'   At fixed calibration the local likelihood standard error is the reciprocal
#'   square root of information. A posterior standard deviation for ability
#'   also depends on the prior, observed responses and calibration uncertainty.
#'   Dependence between ratings of the same performance beyond the fitted
#'   facets can reduce the precision gained from additional ratings.
#' @return A tibble with `theta`, the posterior mean `information`, equal-tailed
#'   credible limits `lower` and `upper`, `ci_level`, and `n_ratings` per examinee.
#' @seealso [category_curves()], [as_mirt()]
#' @examples
#' \donttest{
#' if (check_cmdstan()) {
#'   source(system.file("examples", "fit-writing.R", package = "BayesMFRM"),
#'          local = TRUE)
#'   plan <- expand.grid(item = c("Content", "Evidence"), rater = c("A", "B"))
#'   grid <- seq(-3, 3, by = 0.25)
#'   one <- test_information(fit, subset(plan, rater == "A"), grid)
#'   both <- test_information(fit, plan, grid)
#'   matplot(grid, cbind(one$information, both$information), type = "l", lty = 1:2,
#'           xlab = "Ability (logits)", ylab = "Conditional test information")
#'   legend("topright", c("Rater A", "Raters A and B"), lty = 1:2, col = 1:2, bty = "n")
#'   # With the same two criteria, the second rater adds two ratings per examinee.
#' }
#' }
#' @export
test_information <- function(object, newdata, theta = seq(-4, 4, length.out = 81), ci_level = 0.95) {
  theta <- bmfrm_curve_grid(theta)
  profiles <- bmfrm_curve_profiles(object, newdata, ci_level)
  D <- profiles$params$D
  size <- default_bmfrm_chunk_size(D * as.double(profiles$params$K), length(theta))
  chunks <- split(theta, ceiling(seq_along(theta) / size))
  dplyr::bind_rows(lapply(chunks, function(grid) {
    total <- matrix(0, nrow = D, ncol = length(grid))
    for (i in seq_len(nrow(newdata))) {
      ct <- bmfrm_profile_cumulative(profiles, i)
      p <- bmfrm_curve_draws(ct, grid)
      dim(p) <- c(D * length(grid), ncol(ct))
      scores <- seq_len(ncol(ct))
      mu <- drop(p %*% scores)
      for (k in scores) {
        total <- total + p[, k] * (k - mu)^2
      }
    }
    limits <- apply(total, 2L, stats::quantile, probs = profiles$probs, names = FALSE)
    tibble::tibble(theta = grid, information = colMeans(total), lower = limits[1L, ],
                   upper = limits[2L, ], ci_level = ci_level, n_ratings = nrow(newdata))
  }))
}
