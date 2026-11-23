#' Print method for bmfrm_fit objects
#' @param x bmfrm_fit object
#' @param ... unused
#' @export
print.bmfrm_fit <- function(x, ...) {
  cat("Bayesian Many-Facet Rasch Model\n")
  cat("================================\n\n")
  
  cat("Formula:", deparse(x$spec$formula), "\n")
  cat("Family: ", x$spec$family, "\n")
  cat("Hierarchical facets:", paste(x$hierarchical %||% if (!is.null(x$model)) "person" else character(), collapse = ", "), "\n")
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
    cat("  OK: Model appears to have converged\n")
  } else {
    cat("  CHECK: Potential convergence issues detected\n")
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
  cat("Hierarchical facets:", paste(object$hierarchical %||% if (!is.null(object$model)) "person" else character(), collapse = ", "), "\n")
  cat("Observations:", object$data_info$stan_data$N, "\n")
  cat("Rating categories: 1 to", object$spec$K, "\n\n")
  
  # MCMC diagnostics
  cat("MCMC Diagnostics:\n")
  cat("  Chains:", object$fit$num_chains(), "\n")
  cat("  Total iterations:", object$fit$num_chains() * object$fit$metadata()$iter_sampling, "\n")
  cat("  Max R-hat:", round(object$convergence$max_rhat, 3), "\n")
  cat("  Min ESS (bulk):", object$convergence$min_ess_bulk, "\n")
  cat("  Min ESS (tail):", object$convergence$min_ess_tail, "\n")
  cat("  Min E-BFMI:", object$convergence$min_ebfmi, "\n")
  cat("  Max mean MCSE:", object$convergence$max_mcse_mean, "\n")
  cat("  Maximum-depth hits:", object$convergence$n_max_treedepth, "\n")
  
  if (object$convergence$n_divergent > 0) {
    cat("  Divergent transitions:", object$convergence$n_divergent, "\n")
  }
  
  if (!object$convergence$converged) {
    cat("  Warning: Potential convergence issues detected!\n")
  }
  cat("\n")
  
  # Threshold parameters
  cat("Centered step deviations (tau):\n")
  tau_summary <- object$fit$summary("tau", "mean", ~stats::quantile(.x, probs = range(probs)))
  for (i in seq_len(nrow(tau_summary))) {
    label <- tau_summary$variable[i]
    if (identical(object$spec$family, "pcm")) {
      idx <- as.integer(strsplit(gsub("^tau\\[|\\]$", "", label), ",")[[1L]])
      label <- paste0(label, " (", object$spec$threshold_facet, "=",
        object$data_info$facet_labels[[object$spec$threshold_facet]][idx[1L]], ", step=", idx[2L], ")")
    }
    cat(sprintf("  %s: %6.2f [%6.2f, %6.2f]\n",
                label,
                tau_summary$mean[i],
                tau_summary[[3L]][i],
                tau_summary[[ncol(tau_summary)]][i]))
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
      n_flagged <- sum(bias_summ$flag_bias %in% c("positive_bias", "negative_bias"))
      cat("    ", n_flagged, " bias effects flagged (out of ", nrow(bias_summ), ")\n")
      cat("    ", sum(bias_summ$flag_bias == "not_likelihood_identified"), " effects lack likelihood identification\n")
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

  # Column names in stan_data have "_id" suffix (e.g., "person_id", "rater_id")
  facet_cols_stan <- paste0(object$spec$facets_main, "_id")

  res_df <- extract_residuals(
    fit = object$fit,
    stan_data = object$data_info$stan_data,
    facet_cols = facet_cols_stan,
    save_draws = save_draws,
    model_name = object$spec$model_name
  )

  # Rename columns from "_id" suffix back to original names
  for (i in seq_along(object$spec$facets_main)) {
    orig_name <- object$spec$facets_main[i]
    stan_name <- facet_cols_stan[i]
    if (stan_name %in% names(res_df)) {
      names(res_df)[names(res_df) == stan_name] <- orig_name
    }
  }

  attr(res_df, "facet_cols") <- object$spec$facets_main
  res_df
}

#' Posterior predictive check method for bmfrm_fit objects
#' @param object bmfrm_fit object
#' @param group Optional column name in the cleaned data for conditional PPC.
#' @param type type of PPC ("bars", "stat", or "all")
#' @param ndraws number of posterior draws
#' @param seed random seed
#' @param ... unused
#' @importFrom bayesplot pp_check
#' @export pp_check.bmfrm_fit
#' @exportS3Method bayesplot::pp_check bmfrm_fit
pp_check.bmfrm_fit <- function(object, type = c("bars", "stat", "all"), ndraws = 200, seed = 123, group = NULL, ...) {
  
  type <- match.arg(type)
  
  y_obs <- object$data_info$stan_data$X
  if (!is.null(group)) {
    if (!is.character(group) || length(group) != 1L || !group %in% names(object$data_info$data_clean)) stop("group must name a column in the fitted data")
    group <- object$data_info$data_clean[[group]]
  }
  
  if (type == "all") {
    return(run_ppc(object$fit, y_obs, ndraws = ndraws, seed = seed, group = group))
  } else if (type == "bars") {
    ppc_result <- run_ppc(object$fit, y_obs, ndraws = ndraws, seed = seed, group = group)
    return(ppc_result$ppc_bars)
  } else if (type == "stat") {
    ppc_result <- run_ppc(object$fit, y_obs, ndraws = ndraws, seed = seed, group = group)
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
  param_summary <- object$fit$summary(param_name, "mean", "sd",
                                     ~stats::quantile(.x, probs = probs))
  
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
#'
#' Person summaries average joint posterior draws over that person's observed
#' item mix. They do not standardize different persons to a common item design.
#' PCM steps on an adjusted facet (rater or interlocutor) require a reference
#' step profile, which is not implemented; this case raises an error.
#' Interactions involving rater or interlocutor are removed; other interactions
#' remain in the prediction.
#' @param object bmfrm_fit object
#' @param summary whether to return summary or observation-level scores
#' @param ... unused
#' @export
fair_scores <- function(object, summary = TRUE, ...) {
  UseMethod("fair_scores")
}

#' @export
fair_scores.bmfrm_fit <- function(object, summary = TRUE, ...) {
  if (identical(object$spec$family, "pcm") && object$spec$threshold_facet %in% c("rater", "interlocutor")) {
    stop("Fair scores are not defined for PCM steps on an adjusted facet; a reference step profile must be specified first")
  }
  
  obs_data <- object$data_info$data_clean
  draws <- unclass(object$fit$draws("mu_fair", format = "matrix"))
  index <- as.integer(sub(".*\\[([0-9]+)\\]", "\\1", colnames(draws)))
  draws <- draws[, order(index), drop = FALSE]
  if (ncol(draws) != nrow(obs_data)) stop("mu_fair draws do not match fitted observations")
  obs_data$mu_fair <- colMeans(draws)
  obs_data$mu_fair_sd <- apply(draws, 2L, stats::sd)
  if (!summary) return(obs_data)
  labels <- object$data_info$facet_labels$person
  rows <- lapply(labels, function(label) {
    idx <- which(as.character(obs_data$person) == label)
    if (!length(idx)) return(NULL)
    aggregate_draws <- rowMeans(draws[, idx, drop = FALSE])
    data.frame(person = label, n_obs = length(idx),
               mean_score = mean(obs_data[[object$spec$response]][idx]),
               mean_fair_score = mean(aggregate_draws), sd_fair_score = stats::sd(aggregate_draws),
               lower_95 = unname(stats::quantile(aggregate_draws, .025)),
               upper_95 = unname(stats::quantile(aggregate_draws, .975)))
  })
  tibble::as_tibble(do.call(rbind, rows))

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

  # Column names in stan_data have "_id" suffix (e.g., "person_id", "rater_id")
  facet_cols_stan <- paste0(object$spec$facets_main, "_id")

  # Get residuals with draws for Bayesian fit stats
  res_df <- extract_residuals(
    fit = object$fit,
    stan_data = object$data_info$stan_data,
    facet_cols = facet_cols_stan,
    save_draws = TRUE,
    model_name = object$spec$model_name
  )

  # Rename columns from "_id" suffix back to original names for downstream processing
  for (i in seq_along(object$spec$facets_main)) {
    orig_name <- object$spec$facets_main[i]
    stan_name <- facet_cols_stan[i]
    if (stan_name %in% names(res_df)) {
      names(res_df)[names(res_df) == stan_name] <- orig_name
    }
  }

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
                  infit_zstd, outfit_zstd, dplyr::everything())
  
  fit_stats
}

#' Summarise bias effects (wrapper)
#' @param object bmfrm_fit object  
#' @param facet Two-way interaction name in either component order (if NULL,
#'   uses the first interaction). The generated parameter name and original
#'   level labels are used. For legacy raw CmdStan fits named `bias[I,R]`, use
#'   summarise_bias_ir().
#' @param prob_thresh probability threshold for flagging bias
#' @param delta Nonnegative practically meaningful effect on the adjacent-logit
#'   scale. Report probabilities beyond +/-delta; flags also require the 95%
#'   interval to exclude that region. This is not a multiplicity correction.
#' @return A tibble with original labels, posterior summaries, directional and
#'   practical-effect probabilities, cell observation counts and a likelihood
#'   identification indicator when a design check is available. Unidentified
#'   cells retain posterior summaries but receive flag not_likelihood_identified.
#' @param ... unused
#' @export
summarise_bias <- function(object, facet = NULL, prob_thresh = 0.95, delta = 0, ...) {
  UseMethod("summarise_bias")
}

#' @export
summarise_bias.bmfrm_fit <- function(object, facet = NULL, prob_thresh = 0.95, delta = 0, ...) {
  if (!is.numeric(delta) || length(delta) != 1L || !is.finite(delta) || delta < 0) stop("delta must be one finite nonnegative number")
  if (!is.numeric(prob_thresh) || length(prob_thresh) != 1L || !is.finite(prob_thresh) || prob_thresh < .5 || prob_thresh > 1) stop("prob_thresh must be between .5 and 1")
  
  if (length(object$spec$facets_bias) == 0) {
    stop("No bias facets in model")
  }
  
  if (is.null(facet)) {
    facet <- object$spec$facets_bias[1]
    message("Using bias facet: ", facet)
  }
  
  if (!is.character(facet) || length(facet) != 1L || is.na(facet)) stop("facet must identify one interaction")
  facet <- match_interaction_names(facet, object$spec$facets_bias)
  if (!facet %in% object$spec$facets_bias) {
    stop("Bias facet '", facet, "' not found. Available: ", 
         paste(object$spec$facets_bias, collapse = ", "))
  }
  
  f <- strsplit(facet, ":", fixed = TRUE)[[1L]]
  labels <- object$data_info$facet_labels[f]
  n <- lengths(labels); name <- gsub(":", "", facet)
  draws <- unclass(object$fit$draws(name, format = "matrix"))
  cells <- expand.grid(index1 = seq_len(n[1L]), index2 = seq_len(n[2L]))
  variable <- paste0(name, "[", cells$index1, ",", cells$index2, "]")
  if (!all(variable %in% colnames(draws))) stop("Generated interaction draws are missing; summarise_bias_ir() supports legacy raw fits with a parameter named bias")
  draws <- draws[, variable, drop = FALSE]
  ci <- apply(draws, 2L, stats::quantile, probs = c(.025, .975))
  result <- tibble::tibble(facet = facet, facet1 = f[1L], facet2 = f[2L],
    level1 = labels[[1L]][cells$index1], level2 = labels[[2L]][cells$index2],
    index1 = cells$index1, index2 = cells$index2, variable = variable,
    bias_mean = colMeans(draws), bias_sd = apply(draws, 2L, stats::sd), lower_95 = ci[1L, ], upper_95 = ci[2L, ],
    prob_pos = colMeans(draws > 0), prob_neg = colMeans(draws < 0), delta = delta,
    prob_gt_delta = colMeans(draws > delta), prob_lt_neg_delta = colMeans(draws < -delta),
    prob_outside_delta = colMeans(abs(draws) > delta))
  result$flag_bias <- ifelse(result$lower_95 > delta & result$prob_gt_delta >= prob_thresh, "positive_bias",
    ifelse(result$upper_95 < -delta & result$prob_lt_neg_delta >= prob_thresh, "negative_bias", "no_clear_bias"))
  dat <- object$data_info$data_clean
  idx1 <- match(as.character(dat[[f[1L]]]), labels[[1L]])
  idx2 <- match(as.character(dat[[f[2L]]]), labels[[2L]])
  result$n_observations <- tabulate(idx1 + (idx2 - 1L) * n[1L], prod(n))
  dc <- object$design
  result$likelihood_identified <- NA
  if (!is.null(dc$interaction_levels[[facet]])) {
    rows <- dc$interaction_levels[[facet]]$map_row
    result$likelihood_identified <- rowSums(dc$effect_null_space[rows, , drop = FALSE]^2) < dc$tolerance^2
    result$flag_bias[!result$likelihood_identified] <- "not_likelihood_identified"
  }
  # Preserve familiar column names for the item/rater use case in either order.
  if (setequal(f, c("item", "rater"))) {
    result$item <- if (f[1L] == "item") result$level1 else result$level2
    result$rater <- if (f[1L] == "rater") result$level1 else result$level2
    result$i <- if (f[1L] == "item") result$index1 else result$index2
    result$r <- if (f[1L] == "rater") result$index1 else result$index2
  }
  result
}
