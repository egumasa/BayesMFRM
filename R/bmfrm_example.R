#' Worked four-category writing assessment
#'
#' A simulated rating study of 36 writers, three criteria and three raters.
#' Each writer received six ratings. Raters A and B were simulated with
#' opposing deviations on Content and Evidence. The summaries come from one
#' four-chain RSM fit with an item-by-rater interaction and seed 2409.
#'
#' @format A list with six tables.
#' * `ratings` contains one observed score (1 to 4) per person, item and rater.
#' * `rater` gives posterior rater severity in logits.
#' * `rater_fit` gives matched observed and replicated fit statistics by rater.
#' * `agreement` compares exact agreement on shared ratings with its model
#'   expectation and 95% posterior predictive interval for each pair and overall.
#' * `fair` gives expected scores on the common three-criterion grid.
#' * `bias` gives interaction summaries at delta = 0.3 and BFDR = 0.05.
#' @source Generated with `data-raw/bmfrm_example.R` in the source repository.
#' @seealso [bmfrm()], [summarise_bias()], [plot_bias()], [fair_scores()]
#' @examples
#' data(bmfrm_example)
#' head(bmfrm_example$ratings)
#' plot_bias(bmfrm_example$bias)
"bmfrm_example"
