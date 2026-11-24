#' Run posterior predictive checks
#'
#' @param fit cmdstanr fit object with y_rep in generated quantities (or unused
#'   when `y_rep` is supplied directly)
#' @param y_obs observed response data
#' @param ndraws number of draws to use for PPC
#' @param seed random seed
#' @param group Optional grouping vector for conditional checks.
#' @param y_rep Optional precomputed replicated response matrix (`D x N`),
#'   e.g. from on-demand computation when `save_gq = "minimal"`.
#' 
#' @return list with PPC plots and y_rep subset
#' @export
run_ppc <- function(fit, y_obs, ndraws = 200, seed = 123, group = NULL, y_rep = NULL) {
  # --- 1. Extract posterior draws for y_rep ---
  if (is.null(y_rep)) {
    y_rep_mat <- fit$draws("y_rep", format = "matrix")
  } else {
    y_rep_mat <- y_rep
  }
  
  # --- 2. Subset draws ---
  if (length(ndraws) != 1L || !is.finite(ndraws) || ndraws < 1 || ndraws != floor(ndraws)) stop("ndraws must be a positive integer")
  if (!is.null(group) && (length(group) != length(y_obs) || anyNA(group))) stop("group must match all observed ratings without missing values")
  if (is.null(y_rep) || nrow(y_rep_mat) > ndraws) {
    idx <- with_bmfrm_seed(seed, sample.int(nrow(y_rep_mat), min(ndraws, nrow(y_rep_mat))))
    yrep_subset <- y_rep_mat[idx, , drop = FALSE]
  } else {
    yrep_subset <- y_rep_mat
  }
  
  # --- 3. PPC plots ---
  message("Plotting PPC bars...")
  p1 <- if (is.null(group)) bayesplot::ppc_bars(y = y_obs, yrep = yrep_subset) else
    bayesplot::ppc_bars_grouped(y = y_obs, yrep = yrep_subset, group = group)
  
  message("Plotting PPC mean distribution...")
  p2 <- if (is.null(group)) bayesplot::ppc_stat(y = y_obs, yrep = yrep_subset, stat = "mean") else
    bayesplot::ppc_stat_grouped(y = y_obs, yrep = yrep_subset, group = group, stat = "mean")
  
  message("Plotting PPC variance distribution...")
  p3 <- if (is.null(group)) bayesplot::ppc_stat(y = y_obs, yrep = yrep_subset, stat = stats::var) else
    bayesplot::ppc_stat_grouped(y = y_obs, yrep = yrep_subset, group = group, stat = stats::var)
  
  # --- 4. Return ggplot objects for flexible usage ---
  return(list(
    ppc_bars = p1,
    ppc_mean = p2,
    ppc_var  = p3,
    yrep_subset = yrep_subset
  ))
}
