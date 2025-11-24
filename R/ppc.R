#' Run posterior predictive checks
#'
#' @param fit cmdstanr fit object with y_rep in generated quantities
#' @param y_obs observed response data
#' @param ndraws number of draws to use for PPC
#' @param seed random seed
#' 
#' @return list with PPC plots and y_rep subset
#' @export
run_ppc <- function(fit, y_obs, ndraws = 200, seed = 123) {
  # --- 1. Extract posterior draws for y_rep ---
  y_rep_mat <- fit$draws("y_rep", format = "matrix")
  
  # --- 2. Subset draws ---
  set.seed(seed)
  idx <- sample(seq_len(nrow(y_rep_mat)), ndraws)
  yrep_subset <- y_rep_mat[idx, , drop = FALSE]
  
  # --- 3. PPC plots ---
  message("Plotting PPC bars...")
  p1 <- bayesplot::ppc_bars(y = y_obs, yrep = yrep_subset)
  
  message("Plotting PPC mean distribution...")
  p2 <- bayesplot::ppc_stat(y = y_obs, yrep = yrep_subset, stat = "mean")
  
  message("Plotting PPC variance distribution...")
  p3 <- bayesplot::ppc_stat(y = y_obs, yrep = yrep_subset, stat = function(x) var(x))
  
  # --- 4. Return ggplot objects for flexible usage ---
  return(list(
    ppc_bars = p1,
    ppc_mean = p2,
    ppc_var  = p3,
    yrep_subset = yrep_subset
  ))
}
