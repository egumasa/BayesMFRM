#' Print method for bmfrm_fit objects
#' @param x bmfrm_fit object
#' @param ... unused
#' @export
print.bmfrm_fit <- function(x, ...) {
  cat("Bayesian Many-Facet Rasch Model\n")
  cat("================================\n\n")
  
  cat("Formula:", deparse(x$spec$formula), "\n")
  cat("Family: ", x$spec$family, "\n")
  cat("Observations:", x$data_info$stan_data$N, "\n")
  cat("Rating categories: 1 to", x$spec$K, "\n\n")
  
  cat("Facets:\n")
  cat("  Main facets:", paste(x$spec$facets_main, collapse = ", "), "\n")
  if (length(x$spec$facets_bias) > 0) {
    cat("  Bias facets:", paste(x$spec$facets_bias, collapse = ", "), "\n")
  }
  
  cat("\nFacet levels:\n")
  for (facet in x$spec$facets_main) {
    n_levels <- x$data_info$stan_data[[paste0("J_", facet)]]
    cat("  ", facet, ": ", n_levels, " levels\n", sep = "")
  }
  
  cat("\nMCMC sampling:\n")
  cat("  Chains:", x$fit$num_chains(), "\n")
  cat("  Iterations per chain:", x$fit$metadata()$iter_sampling, "\n")
  cat("  Warmup per chain:", x$fit$metadata()$iter_warmup, "\n")
  
  cat("\nConvergence:\n")
  if (x$convergence$converged) {
    cat("  ✓ Model appears to have converged\n")
  } else {
    cat("  ✗ Potential convergence issues detected\n")
    cat("    Max R-hat:", round(x$convergence$max_rhat, 3), "\n")
    if (x$convergence$n_divergent > 0) {
      cat("    Divergent transitions:", x$convergence$n_divergent, "\n")
    }
  }
  
  cat("\nUse summary() for detailed results\n")
  cat("Use pp_check() for posterior predictive checks\n")
  cat("Use residuals() for residual analysis\n")
}

#' Summary method for bmfrm_fit objects
#' @param object bmfrm_fit object
#' @param facets which facets to summarize (default: all main facets)
#' @param probs quantiles for credible intervals
#' @param ... unused
#' @export
summary.bmfrm_fit <- function(object, facets = NULL, probs = c(0.025, 0.5, 0.975), ...) {
  
  if (is.null(facets)) {
    facets <- object$spec$facets_main
  }
  
  cat("Bayesian Many-Facet Rasch Model Summary\n")
  cat("=======================================\n\n")
  
  # Model info
  cat("Formula:", deparse(object$spec$formula), "\n")
  cat("Family: ", object$spec$family, "\n") 
  cat("Observations:", object$data_info$stan_data$N, "\n")
  cat("Rating categories: 1 to", object$spec$K, "\n\n")
  
  # MCMC diagnostics
  cat("MCMC Diagnostics:\n")
  cat("  Chains:", object$fit$num_chains(), "\n")
  cat("  Total iterations:", object$fit$num_chains() * object$fit$metadata()$iter_sampling, "\n")
  cat("  Max R-hat:", round(object$convergence$max_rhat, 3), "\n")
  cat("  Min ESS (bulk):", object$convergence$min_ess_bulk, "\n")
  cat("  Min ESS (tail):", object$convergence$min_ess_tail, "\n")
  
  if (object$convergence$n_divergent > 0) {
    cat("  Divergent transitions:", object$convergence$n_divergent, "\n")
  }
  
  if (!object$convergence$converged) {
    cat("  ⚠ Potential convergence issues detected!\n")
  }
  cat("\n")
  
  # Threshold parameters
  cat("Threshold Parameters:\n")
  tau_summary <- object$fit$summary("tau")
  for (i in seq_len(nrow(tau_summary))) {
    cat(sprintf("  tau[%d]: %6.2f [%6.2f, %6.2f]\n",
                i, 
                tau_summary$mean[i],
                tau_summary$q5[i],
                tau_summary$q95[i]))
  }
  cat("\n")
  
  # Facet summaries
  for (facet in facets) {
    cat(paste0(toupper(substring(facet, 1, 1)), substring(facet, 2)), " Parameters:\n")
    facet_summ <- facet_summary(object, facet, probs = probs)
    print(facet_summ, row.names = FALSE)
    cat("\n")
  }
  
  # Bias summaries if present
  if (length(object$spec$facets_bias) > 0) {
    cat("Bias Effects:\n")
    for (bias_facet in object$spec$facets_bias) {
      cat("  ", bias_facet, ":\n")
      bias_summ <- summarise_bias(object, bias_facet)
      n_flagged <- sum(bias_summ$flag_bias != "no_clear_bias")
      cat("    ", n_flagged, " bias effects flagged (out of ", nrow(bias_summ), ")\n")
    }
    cat("    Use summarise_bias(fit) for details\n\n")
  }
  
  invisible(object)
}

#' Residuals method for bmfrm_fit objects
#' @param object bmfrm_fit object
#' @param save_draws whether to include posterior draws
#' @param ... unused
#' @export
residuals.bmfrm_fit <- function(object, save_draws = FALSE, ...) {
  extract_residuals(
    fit = object$fit,
    stan_data = object$data_info$stan_data,
    facet_cols = object$spec$facets_main,
    save_draws = save_draws,
    model_name = object$spec$model_name
  )
}

#' Posterior predictive check method for bmfrm_fit objects
#' @param object bmfrm_fit object
#' @param type type of PPC ("bars", "stat", or "all")
#' @param ndraws number of posterior draws
#' @param seed random seed
#' @param ... unused
#' @export
pp_check.bmfrm_fit <- function(object, type = c("bars", "stat", "all"), ndraws = 200, seed = 123, ...) {
  
  type <- match.arg(type)
  
  y_obs <- object$data_info$stan_data$X
  
  if (type == "all") {
    return(run_ppc(object$fit, y_obs, ndraws = ndraws, seed = seed))
  } else if (type == "bars") {
    ppc_result <- run_ppc(object$fit, y_obs, ndraws = ndraws, seed = seed)
    return(ppc_result$ppc_bars)
  } else if (type == "stat") {
    ppc_result <- run_ppc(object$fit, y_obs, ndraws = ndraws, seed = seed)
    return(list(
      ppc_mean = ppc_result$ppc_mean,
      ppc_var = ppc_result$ppc_var
    ))
  }
}

#' Get Stan code from bmfrm_fit object
#' @param object bmfrm_fit object
#' @param ... unused
#' @export
stancode.bmfrm_fit <- function(object, ...) {
  cat(object$stan_code)
  invisible(object$stan_code)
}

#' Get Stan data from bmfrm_fit object
#' @param object bmfrm_fit object
#' @param ... unused
#' @export
standata.bmfrm_fit <- function(object, ...) {
  object$data_info$stan_data
}

#' Extract facet-level summaries
#' @param object bmfrm_fit object
#' @param facet name of facet to summarize
#' @param probs quantiles for credible intervals
#' @param ... unused
#' @export
facet_summary <- function(object, facet, probs = c(0.025, 0.5, 0.975), ...) {
  UseMethod("facet_summary")
}

#' @export
facet_summary.bmfrm_fit <- function(object, facet, probs = c(0.025, 0.5, 0.975), ...) {
  
  # Check facet exists
  if (!facet %in% object$spec$facets_main) {
    stop("Facet '", facet, "' not found in model. Available facets: ", 
         paste(object$spec$facets_main, collapse = ", "))
  }
  
  # Map facet name to parameter name
  if (facet == "person") {
    param_name <- "theta"
  } else {
    param_name <- facet
  }
  
  # Extract parameter summary
  param_summary <- object$fit$summary(param_name, 
                                     ~quantile(.x, probs = probs))
  
  # Get original labels
  labels <- object$data_info$facet_labels[[facet]]
  
  # Build result
  result <- tibble::tibble(
    facet = facet,
    level = labels,
    index = seq_along(labels)
  )
  
  # Add parameter summaries
  result$mean <- param_summary$mean
  result$sd <- param_summary$sd
  
  # Add quantiles
  for (i in seq_along(probs)) {
    prob_name <- paste0("q", round(probs[i] * 100))
    result[[prob_name]] <- param_summary[[ncol(param_summary) - length(probs) + i]]
  }
  
  result
}

#' Extract fair scores
#' @param object bmfrm_fit object
#' @param summary whether to return summary or observation-level scores
#' @param ... unused
#' @export
fair_scores <- function(object, summary = TRUE, ...) {
  UseMethod("fair_scores")
}

#' @export
fair_scores.bmfrm_fit <- function(object, summary = TRUE, ...) {
  
  # Extract mu_fair from generated quantities
  mu_fair_summary <- object$fit$summary("mu_fair")
  
  # Build observation-level data
  obs_data <- object$data_info$data_clean
  
  # Add fair scores
  obs_data$mu_fair <- mu_fair_summary$mean
  obs_data$mu_fair_sd <- mu_fair_summary$sd
  
  if (!summary) {
    return(obs_data)
  }
  
  # Person-level summary
  person_summary <- obs_data %>%
    dplyr::group_by(.data[[object$spec$facets_main[1]]]) %>%  # assumes first facet is person
    dplyr::summarise(
      n_obs = dplyr::n(),
      mean_score = mean(.data[[object$spec$response]], na.rm = TRUE),
      mean_fair_score = mean(mu_fair, na.rm = TRUE),
      sd_fair_score = sqrt(mean(mu_fair_sd^2)),  # approximate
      .groups = "drop"
    )
  
  person_summary
}

#' Extract facet fit statistics
#' @param object bmfrm_fit object
#' @param facet name of facet
#' @param ... unused
#' @export
facet_fit <- function(object, facet, ...) {
  UseMethod("facet_fit")
}

#' @export
facet_fit.bmfrm_fit <- function(object, facet, ...) {
  
  # Get residuals with draws for Bayesian fit stats
  res_df <- extract_residuals(
    fit = object$fit,
    stan_data = object$data_info$stan_data,
    facet_cols = object$spec$facets_main,
    save_draws = TRUE,
    model_name = object$spec$model_name
  )
  
  # Basic fit statistics
  fit_stats <- facet_infit_outfit(res_df, facet)
  
  # Add Bayesian credible intervals
  fit_stats <- add_fit_cri_bayes(fit_stats, res_df, facet)
  
  # Add z-standardized statistics
  fit_stats <- add_zstd(fit_stats)
  
  # Add original labels
  labels <- object$data_info$facet_labels[[facet]]
  fit_stats$label <- labels[fit_stats$facet_value]
  
  # Reorder columns
  fit_stats <- fit_stats %>%
    dplyr::select(facet_value, label, n_obs, infit_msq, outfit_msq,
                  infit_lo, infit_hi, outfit_lo, outfit_hi,
                  infit_zstd, outfit_zstd, everything())
  
  fit_stats
}

#' Summarise bias effects (wrapper)
#' @param object bmfrm_fit object  
#' @param facet bias facet name (if NULL, uses first bias facet)
#' @param prob_thresh probability threshold for flagging bias
#' @param ... unused
#' @export
summarise_bias <- function(object, facet = NULL, prob_thresh = 0.95, ...) {
  UseMethod("summarise_bias")
}

#' @export
summarise_bias.bmfrm_fit <- function(object, facet = NULL, prob_thresh = 0.95, ...) {
  
  if (length(object$spec$facets_bias) == 0) {
    stop("No bias facets in model")
  }
  
  if (is.null(facet)) {
    facet <- object$spec$facets_bias[1]
    message("Using bias facet: ", facet)
  }
  
  if (!facet %in% object$spec$facets_bias) {
    stop("Bias facet '", facet, "' not found. Available: ", 
         paste(object$spec$facets_bias, collapse = ", "))
  }
  
  # Currently only support rater:item bias
  if (facet == "rater:item") {
    item_labels <- object$data_info$facet_labels$item
    rater_labels <- object$data_info$facet_labels$rater
    
    return(summarise_bias_ir(
      fit = object$fit,
      item_labels = item_labels,
      rater_labels = rater_labels,
      prob_thresh = prob_thresh
    ))
  } else {
    stop("Bias facet '", facet, "' not yet supported. Currently only 'rater:item' is implemented.")
  }
}
