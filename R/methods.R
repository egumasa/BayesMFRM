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
  if (isTRUE(x$spec$K_heterogeneous) && !is.null(x$spec$K_by_level)) {
    cat("Rating categories: heterogeneous by ", x$spec$threshold_facet, " (",
        paste(sprintf("%s=1..%d", names(x$spec$K_by_level), as.integer(x$spec$K_by_level)), collapse = ", "),
        ")\n\n", sep = "")
  } else {
    cat("Rating categories: 1 to", get_spec_K_max(x$spec), "\n\n")
  }
  
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
  
  cat("\nDesign Identification:\n")
  if (is.null(x$design)) {
    cat("  NOT CHECKED: design_check = 'none' or custom Stan model (run check_design_bmfrm() to inspect)\n")
  } else if (isTRUE(x$design$identifiable)) {
    cat("  OK: Conditional likelihood is structurally identifiable (rank ",
        x$design$adjacent_logit_rank, "/", x$design$adjacent_logit_parameters, ")\n", sep = "")
  } else {
    n_unid <- x$design$adjacent_logit_parameters - x$design$adjacent_logit_rank
    cat("  WARNING: Rank-deficient design (", n_unid, " unidentified direction(s); connected = ",
        x$design$connected, ")\n", sep = "")
    if (length(x$design$unobserved_step_levels)) {
      cat("    Unobserved PCM step levels:", paste(x$design$unobserved_step_levels, collapse = ", "), "\n")
    }
    cat("    Proper priors yield a posterior even without likelihood identification.\n")
    cat("    Inspect fit$design and check 'likelihood_identified' in fair_scores() / summarise_bias().\n")
  }

  cat("\nMCMC Convergence:\n")
  if (x$convergence$converged) {
    cat("  OK: MCMC diagnostics passed (R-hat < 1.01, ESS >= 400, 0 divergences, E-BFMI >= 0.3)\n")
    if (!is.null(x$design) && !isTRUE(x$design$identifiable)) {
      cat("    Note: Sampling stability does NOT resolve the design rank deficiency above.\n")
    } else {
      cat("    Note: Convergence checks sampling stability only, not model fit to data.\n")
    }
  } else {
    cat("  CHECK: Potential MCMC convergence issues detected\n")
    cat("    Max R-hat:", round(x$convergence$max_rhat, 3), "\n")
    if (!is.null(x$convergence$min_ess_bulk) && is.finite(x$convergence$min_ess_bulk)) {
      cat("    Min ESS (bulk / tail):", round(x$convergence$min_ess_bulk, 1), "/",
          round(x$convergence$min_ess_tail, 1), "\n")
    }
    if (x$convergence$n_divergent > 0) {
      cat("    Divergent transitions:", x$convergence$n_divergent, "\n")
    }
    if (!is.null(x$convergence$n_max_treedepth) && x$convergence$n_max_treedepth > 0) {
      cat("    Max treedepth hits:", x$convergence$n_max_treedepth, "\n")
    }
    cat("    Resolve sampling issues before interpreting parameter estimates.\n")
  }
  
  cat("\nNext steps (in priority order):\n")
  cat("  1. Design & convergence : summary(fit)\n")
  cat("  2. Data fit (PPC)       : bayesplot::pp_check(fit, group = \"rater\")\n")
  cat("  3. Residuals & misfit   : residuals(fit), facet_fit(fit, \"rater\")\n")
  cat("  4. Substantive outputs  : fair_scores(fit, standardize = TRUE), summarise_bias(fit, delta = 0.3)\n")
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
  if (isTRUE(object$spec$K_heterogeneous) && !is.null(object$spec$K_by_level)) {
    cat("Rating categories: heterogeneous by ", object$spec$threshold_facet, " (",
        paste(sprintf("%s=1..%d", names(object$spec$K_by_level), as.integer(object$spec$K_by_level)), collapse = ", "),
        ")\n\n", sep = "")
  } else {
    cat("Rating categories: 1 to", get_spec_K_max(object$spec), "\n\n")
  }
  
  # Design identification (always inspect before MCMC convergence and estimates)
  cat("Design Identification:\n")
  if (is.null(object$design)) {
    cat("  Status: Not checked (design_check = 'none' or custom Stan model)\n\n")
  } else if (isTRUE(object$design$identifiable)) {
    cat("  Status: OK (structurally identifiable; rank ",
        object$design$adjacent_logit_rank, "/", object$design$adjacent_logit_parameters,
        ", connected = ", object$design$connected, ")\n\n", sep = "")
  } else {
    n_unid <- object$design$adjacent_logit_parameters - object$design$adjacent_logit_rank
    cat("  Status: WARNING - rank-deficient design (", n_unid,
        " unidentified direction(s); connected = ", object$design$connected, ")\n", sep = "")
    cat("  Note: Unidentified contrasts depend on the prior even when MCMC converges.\n")
    cat("        Inspect fit$design and check 'likelihood_identified' before reporting.\n\n")
  }

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
  tau_tf <- object$spec$threshold_facet
  tau_tf_labels <- if (!is.null(tau_tf)) object$data_info$facet_labels[[tau_tf]] else NULL
  tau_start_vec <- object$data_info$stan_data$tau_start
  tau_len_vec <- object$data_info$stan_data$tau_len
  if (isTRUE(object$spec$K_heterogeneous) &&
      (is.null(tau_start_vec) || is.null(tau_len_vec)) &&
      !is.null(tau_tf_labels)) {
    k_map_s <- object$data_info$K_by_level %||% object$spec$K_by_level
    if (!is.null(k_map_s)) {
      tau_len_vec <- as.integer(unname(k_map_s[tau_tf_labels]) - 1L)
      tau_start_vec <- as.integer(c(1L, utils::head(cumsum(tau_len_vec), -1L) + 1L))
    }
  }
  for (i in seq_len(nrow(tau_summary))) {
    label <- tau_summary$variable[i]
    if (identical(object$spec$family, "pcm")) {
      idx <- as.integer(strsplit(gsub("^tau\\[|\\]$", "", label), ",")[[1L]])
      if (length(idx) == 1L && !is.null(tau_start_vec) && !is.null(tau_len_vec) && !is.null(tau_tf_labels)) {
        p <- idx[1L]
        lv_i <- which(p >= tau_start_vec & p < (tau_start_vec + tau_len_vec))
        if (length(lv_i) == 1L) {
          step_m <- p - tau_start_vec[lv_i] + 1L
          label <- paste0(label, " (", tau_tf, "=", tau_tf_labels[lv_i], ", step=", step_m, ")")
        }
      } else if (length(idx) >= 2L && !is.null(tau_tf_labels)) {
        label <- paste0(label, " (", tau_tf, "=", tau_tf_labels[idx[1L]], ", step=", idx[2L], ")")
      }
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
      cat("    Numerical decision stability: ", sum(bias_summ$decision_stable %in% TRUE),
          " stable, ", sum(bias_summ$decision_stable %in% FALSE), " boundary-sensitive, ",
          sum(is.na(bias_summ$decision_stable)), " unavailable\n", sep = "")
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

  if (fit_has_variable(object$fit, "mu") && fit_has_variable(object$fit, "sigma2")) {
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
    return(res_df)
  }

  # On-demand computation for save_gq = "minimal" (requires verified R-engine contract)
  if (!can_evaluate_in_r(object)) {
    stop(
      "Cannot recompute residuals in R: ", get_eval_contract(object)$reason,
      " Save 'mu' and 'sigma2' in generated quantities or fit a generated model."
    )
  }
  if (!save_draws) {
    st <- stream_bmfrm_residuals(object)
    mu_hat <- st$mu_hat
    sig_hat <- st$sigma2_hat
    N <- length(mu_hat)
    q <- NULL
  } else {
    q <- compute_bmfrm_quantities(object)
    N <- ncol(q$mu)
    mu_hat <- unname(colMeans(q$mu))
    sig_hat <- unname(colMeans(q$sigma2))
  }
  x <- object$data_info$stan_data$X
  res_df <- tibble::tibble(n = seq_len(N), x = x)
  for (f in object$spec$facets_main) {
    res_df[[f]] <- object$data_info$stan_data[[paste0(f, "_id")]]
  }
  res_df <- res_df %>%
    dplyr::mutate(
      mu_hat = mu_hat,
      sigma2_hat = sig_hat,
      resid = x - mu_hat,
      weight = sigma2_hat,
      z = resid / sqrt(sigma2_hat),
      model = object$spec$model_name
    )
  if (save_draws) {
    res_df$mu_draws <- lapply(seq_len(N), function(j) as.numeric(q$mu[, j]))
    res_df$sigma2_draws <- lapply(seq_len(N), function(j) as.numeric(q$sigma2[, j]))
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
  if (is.null(group) && isTRUE(object$spec$K_heterogeneous)) {
    message(
      "Note: Heterogeneous category counts across '", object$spec$threshold_facet,
      "' are pooled in ungrouped PPC; pass group = '", object$spec$threshold_facet,
      "' to inspect category distributions on each valid rubric scale."
    )
  }
  if (!is.null(group)) {
    if (!is.character(group) || length(group) != 1L || !group %in% names(object$data_info$data_clean)) stop("group must name a column in the fitted data")
    group <- object$data_info$data_clean[[group]]
  }

  y_rep_override <- NULL
  if (!fit_has_variable(object$fit, "y_rep")) {
    if (!can_evaluate_in_r(object)) {
      stop(
        "Cannot recompute posterior predictive draws in R: ", get_eval_contract(object)$reason,
        " Save 'y_rep' in generated quantities or fit a generated model."
      )
    }
    theta_all <- extract_ordered_draws(object$fit, "theta")
    total_draws <- nrow(theta_all)
    if (length(ndraws) != 1L || !is.finite(ndraws) || ndraws < 1 || ndraws != floor(ndraws)) stop("ndraws must be a positive integer")
    idx <- with_bmfrm_seed(seed, sample.int(total_draws, min(ndraws, total_draws)))
    y_rep_override <- compute_bmfrm_quantities(
      object, draw_idx = idx, compute_y_rep = TRUE, seed = seed
    )$y_rep
  }
  
  if (type == "all") {
    return(run_ppc(object$fit, y_obs, ndraws = ndraws, seed = seed, group = group, y_rep = y_rep_override))
  } else if (type == "bars") {
    ppc_result <- run_ppc(object$fit, y_obs, ndraws = ndraws, seed = seed, group = group, y_rep = y_rep_override)
    return(ppc_result$ppc_bars)
  } else if (type == "stat") {
    ppc_result <- run_ppc(object$fit, y_obs, ndraws = ndraws, seed = seed, group = group, y_rep = y_rep_override)
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
  
  # Add quantiles with exact non-colliding names (e.g., q2.5, q50, q97.5)
  fmt_q <- function(p) {
    paste0("q", sub("\\.?0+$", "", formatC(p * 100, format = "f", digits = 4)))
  }
  for (i in seq_along(probs)) {
    prob_name <- fmt_q(probs[i])
    result[[prob_name]] <- param_summary[[ncol(param_summary) - length(probs) + i]]
  }
  
  result
}

#' Predict from a fitted BayesMFRM model
#'
#' Computes posterior expectations (`type = "epred"`), replicated discrete
#' ratings (`type = "predict"`), category probabilities (`type = "probs"`), or
#' logit-scale linear predictors (`type = "linpred"`) on observed or new data
#' using the shared memory-bounded posterior evaluation engine.
#'
#' Explicitly distinguishes three prediction targets:
#' \itemize{
#'   \item **Conditional on known levels** (default `adjust_facets = character(0)`,
#'     `allow_new_levels = FALSE`): predicts scores for observed examinees,
#'     criteria, and raters.
#'   \item **Population new levels** (`allow_new_levels = TRUE`,
#'     `sample_new_levels = "hierarchical"`): draws unseen facet levels (e.g., a
#'     new examinee or a new rater) once per posterior iteration from their
#'     learned hierarchical population distribution, sharing the same draw across
#'     repeated rows of the same unseen level ID. Raises an error if an unseen
#'     level belongs to a non-hierarchical facet so unseen levels are never
#'     silently replaced with 0.
#'   \item **Adjusted / reference-zero predictions** (e.g.,
#'     `adjust_facets = "rater"`): fixes adjusted nuisance facets (and any
#'     interactions involving them) to 0 on the logit scale. Note that because
#'     the RSM/PCM expectation is nonlinear, fixing a rater effect to 0 is not
#'     identical to averaging over a new rater drawn from the rater population.
#' }
#'
#' Metadata describing the prediction target, adjusted facets, and any unseen
#' levels is attached as `attr(res, "bmfrm_prediction_meta")`.
#'
#' @param object A `bmfrm_fit` object.
#' @param newdata Optional data.frame of facet combinations to predict for.
#'   Defaults to the cleaned fitted dataset (`object$data_info$data_clean`).
#' @param type Quantity to compute: `"epred"` (expected score in `1..K`,
#'   default), `"predict"` (replicated discrete rating in `1..K`), `"probs"`
#'   (category probabilities `1..K`), or `"linpred"` (logit linear predictor).
#' @param adjust_facets Character vector of non-person main facets to set to 0
#'   (and exclude from `newdata` requirements). Default `character(0)` predicts
#'   conditional on all model facets.
#' @param allow_new_levels Logical; allow facet levels in `newdata` that were
#'   not present in the fitted data? Default `FALSE` (raises an error).
#' @param sample_new_levels How to handle unseen levels when
#'   `allow_new_levels = TRUE`: `"hierarchical"` (default; sample from learned
#'   hierarchical population distributions and error on non-hierarchical facets)
#'   or `"zero"` (explicitly set unseen effects to 0 on the logit scale).
#' @param summary Logical; if `TRUE` (default), return a row-aligned tibble of
#'   posterior summaries (`estimate`, `sd`, and `probs` quantiles, or
#'   `prob_1..prob_K` for `type = "probs"`). If `FALSE`, return raw posterior
#'   draws (`D x N` matrix, or `D x N x K` array for `type = "probs"`).
#' @param probs Numeric vector of quantile probabilities for `summary = TRUE`
#'   (ignored for `type = "probs"`). Default `c(0.025, 0.975)`.
#' @param ndraws Optional positive integer specifying a random subset of
#'   posterior draws to evaluate. Mutually exclusive with `draw_ids`.
#' @param draw_ids Optional integer vector of specific posterior draw indices.
#'   Mutually exclusive with `ndraws`.
#' @param seed Optional integer RNG seed for reproducible new-level sampling,
#'   `ndraws` subsampling, and `type = "predict"` replicated draws.
#' @param chunk_size Optional observation chunk size for memory-bounded
#'   streaming evaluation.
#' @param max_cells Maximum total numeric entries allowed when `summary = FALSE`
#'   (default `5e7`, ~400 MB) to prevent accidental memory exhaustion on large
#'   `D x N x K` probability arrays.
#' @param ... Unused.
#' @return A tibble (when `summary = TRUE`) or numeric/integer matrix or 3D
#'   array (when `summary = FALSE`), with attribute `"bmfrm_prediction_meta"`.
#' @export
predict.bmfrm_fit <- function(object,
                              newdata = NULL,
                              type = c("epred", "predict", "probs", "linpred"),
                              adjust_facets = character(0),
                              allow_new_levels = FALSE,
                              sample_new_levels = c("hierarchical", "zero"),
                              summary = TRUE,
                              probs = c(0.025, 0.975),
                              ndraws = NULL,
                              draw_ids = NULL,
                              seed = NULL,
                              chunk_size = NULL,
                              max_cells = 5e7,
                              ...) {
  type <- match.arg(type)
  sample_new_levels <- match.arg(sample_new_levels)
  if (is.null(adjust_facets)) adjust_facets <- character(0)
  resolved_adjust <- resolve_adjust_facets(
    object$spec,
    adjust_facets = adjust_facets,
    warn_if_empty = FALSE
  )
  stream_bmfrm_predictions(
    object = object,
    newdata = newdata,
    type = type,
    adjust_facets = resolved_adjust,
    allow_new_levels = allow_new_levels,
    sample_new_levels = sample_new_levels,
    summary = summary,
    probs = probs,
    ndraws = ndraws,
    draw_ids = draw_ids,
    seed = seed,
    chunk_size = chunk_size,
    max_cells = max_cells
  )
}

#' Extract fair scores
#'
#' Computes expected scores after setting adjusted nuisance facets (and any
#' interactions involving them) to zero (average severity). By default,
#' `"rater"` and `"interlocutor"` are detected automatically; if neither is
#' present and other non-person facets exist, a warning is issued so custom
#' facet names (e.g., `"judge"`, `"examiner"`) are not silently left unadjusted.
#'
#' When `standardize = TRUE` (or `reference_grid` / `reference_weights` is
#' supplied), person summaries are evaluated over a common reference condition
#' spanning the unadjusted facets, optionally weighted by `reference_weights`
#' (e.g., rubric blueprint weights such as `list(item = c(content = 0.4, org = 0.3, lang = 0.3))`
#' or a numeric row-weight vector / `reference_grid$weight` column).
#'
#' When `decompose = TRUE` (with `summary = TRUE`), computes exact draw-matched
#' decompositions of `Fair - Raw` on the same posterior draws:
#' \itemize{
#'   \item Under `standardize = FALSE`: `fair_minus_raw = adjustment + shrinkage_and_noise`,
#'     separating pure nuisance-facet adjustment (`E[mu_fair - mu_actual]`) from
#'     hierarchical shrinkage and sampling noise (`E[mu_actual] - mean_score`).
#'   \item Under `standardize = TRUE`: `fair_minus_raw = standardization_shift + adjustment + shrinkage_and_noise`,
#'     further separating the design/criterion weighting shift from pure rater adjustment.
#' }
#'
#' When `return_draws = TRUE` (with `summary = TRUE`), attaches person-level
#' posterior draw matrices (`D x P`) as `attr(res, "draws")` so users can
#' compute draw-matched score contrasts across reference conditions without
#' subtracting marginal credible intervals.
#'
#' For non-hierarchical PCM with steps on an adjusted facet, no reference step
#' profile exists and an error is raised. For hierarchical PCM
#' (`hierarchical = "tau"`), the shared `tau_common` profile is used.
#'
#' @param object bmfrm_fit object
#' @param summary whether to return person-level summaries or row-level scores
#' @param adjust_facets Optional character vector of nuisance facet names to
#'   remove (default NULL auto-detects `"rater"` and `"interlocutor"`).
#' @param standardize Logical; if TRUE, standardize every person to the common
#'   reference grid of unadjusted non-person facets (or `reference_grid`).
#' @param reference_grid Optional data.frame of unadjusted non-person facet
#'   combinations (and optional `weight` column) used for standardized
#'   evaluation.
#' @param reference_weights Optional weights over `reference_grid`: either a
#'   non-negative numeric vector of length `nrow(reference_grid)` or a named
#'   list of named facet-level weight vectors (e.g.,
#'   `list(item = c(k1 = 0.4, k2 = 0.3, k3 = 0.3))`). Normalized internally to
#'   sum to 1.
#' @param decompose Logical; if TRUE (requires `summary = TRUE`), add exact
#'   draw-matched decomposition columns separating pure facet adjustment,
#'   optional standardization shift, and shrinkage/sampling noise with 95% CrIs.
#' @param return_draws Logical; if TRUE (requires `summary = TRUE`), attach
#'   `D x P` posterior draw matrices as `attr(res, "draws")` and reference
#'   metadata as `attr(res, "reference_meta")`.
#' @param max_grid_cells Maximum number of cells allowed when auto-expanding the
#'   common reference grid for `standardize = TRUE` (default `10000L`). Prevents
#'   combinatorial explosion with many unadjusted facets; supply an explicit
#'   `reference_grid` or raise `max_grid_cells` deliberately for larger grids.
#' @param max_output_rows Maximum total person x reference-grid rows allowed
#'   when `standardize = TRUE` and `summary = FALSE` (default `1000000L`).
#' @param chunk_size Optional observation chunk size for streaming posterior
#'   evaluation (auto-sized to ~16 MB per chunk when NULL).
#' @param ... unused
#' @export
fair_scores <- function(object, summary = TRUE, adjust_facets = NULL,
                        standardize = FALSE, reference_grid = NULL,
                        reference_weights = NULL,
                        decompose = FALSE,
                        return_draws = FALSE,
                        max_grid_cells = 10000L, max_output_rows = 1e6L,
                        chunk_size = NULL, ...) {
  UseMethod("fair_scores")
}

#' @export
fair_scores.bmfrm_fit <- function(object, summary = TRUE, adjust_facets = NULL,
                                  standardize = FALSE, reference_grid = NULL,
                                  reference_weights = NULL,
                                  decompose = FALSE,
                                  return_draws = FALSE,
                                  max_grid_cells = 10000L, max_output_rows = 1e6L,
                                  chunk_size = NULL, ...) {
  if (isTRUE(decompose) && !isTRUE(summary)) {
    stop("decompose = TRUE requires summary = TRUE")
  }
  if (isTRUE(return_draws) && !isTRUE(summary)) {
    stop("return_draws = TRUE requires summary = TRUE; use predict(..., summary = FALSE) for row-level draw matrices")
  }

  has_stan_mu_fair <- fit_has_variable(object$fit, "mu_fair")
  use_default_adjust <- is.null(adjust_facets)
  stan_adjusted <- stan_mu_fair_adjusted_facets(object)
  use_reference <- isTRUE(standardize) || !is.null(reference_grid) || !is.null(reference_weights)
  resolved_adjust <- resolve_adjust_facets(
    object$spec,
    adjust_facets = adjust_facets,
    warn_if_empty = !(has_stan_mu_fair && use_default_adjust &&
                      !use_reference &&
                      length(setdiff(object$spec$facets_main, "person")) <= 1L)
  )

  is_hetero_pcm <- isTRUE(object$spec$K_heterogeneous)
  K_by_level_meta <- object$data_info$K_by_level %||% object$spec$K_by_level

  if (identical(object$spec$family, "pcm") && object$spec$threshold_facet %in% resolved_adjust) {
    if (is_hetero_pcm) {
      stop(
        "Fair scores are not defined when adjusting threshold_facet '",
        object$spec$threshold_facet, "' in heterogeneous-category PCM because no common ",
        "reference step profile exists across differing category counts."
      )
    }
    if (!fit_has_variable(object$fit, "tau_common")) {
      stop("Fair scores are not defined for PCM steps on an adjusted facet; a reference step profile must be specified first")
    }
  }

  warn_identifiability <- function(id_info) {
    basis <- id_info$identification_basis
    if (any(basis == "unidentified_predictor_or_step", na.rm = TRUE)) {
      warning(
        "One or more fair scores depend on likelihood-unidentified predictor/step directions and rely on the prior; ",
        "inspect 'likelihood_identified' and 'identification_basis'.",
        call. = FALSE
      )
    } else if (any(basis == "hierarchical_step_prior", na.rm = TRUE)) {
      warning(
        "Hierarchical PCM fair scores rely on the hierarchical step prior to separate 'tau_common' from 'tau_dev' ",
        "(not identified by the conditional likelihood alone; see 'likelihood_identified' and 'identification_basis').",
        call. = FALSE
      )
    }
  }

  obs_data <- object$data_info$data_clean
  labels <- object$data_info$facet_labels$person

  # Standardized / weighted common-reference fair scores across unadjusted facets (zero P x n_ref expansion)
  if (use_reference) {
    if (!can_evaluate_in_r(object, adjust_facets = resolved_adjust)) {
      stop("Cannot compute standardized fair scores in R: ", get_eval_contract(object)$reason)
    }
    keep_facets <- setdiff(object$spec$facets_main, c("person", resolved_adjust))
    if (is.null(reference_grid)) {
      if (length(keep_facets)) {
        ref_list <- object$data_info$facet_labels[keep_facets]
        n_cells <- prod(as.numeric(lengths(ref_list)))
        if (!is.numeric(max_grid_cells) || length(max_grid_cells) != 1L || !is.finite(max_grid_cells) || max_grid_cells < 1) {
          stop("max_grid_cells must be a positive number")
        }
        if (n_cells > max_grid_cells) {
          stop(
            "Auto-expanded reference grid has ", format(n_cells, big.mark = ","),
            " cells across unadjusted facets (", paste(keep_facets, collapse = ", "),
            "), exceeding max_grid_cells = ", format(max_grid_cells, big.mark = ","),
            ". Supply an explicit 'reference_grid' or raise 'max_grid_cells' deliberately."
          )
        }
        reference_grid <- expand.grid(ref_list, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
      } else {
        reference_grid <- data.frame(.row = 1L)
      }
    } else {
      if (!is.data.frame(reference_grid) || !nrow(reference_grid)) {
        stop("reference_grid must be a non-empty data.frame")
      }
      missing_cols <- setdiff(keep_facets, names(reference_grid))
      if (length(missing_cols)) {
        stop("reference_grid is missing unadjusted facet column(s): ", paste(missing_cols, collapse = ", "))
      }
    }

    ref_w <- resolve_reference_weights(
      reference_grid = reference_grid,
      reference_weights = reference_weights,
      facet_labels = object$data_info$facet_labels
    )
    has_explicit_weights <- !is.null(reference_weights) ||
      ("weight" %in% names(reference_grid) && !("weight" %in% names(object$data_info$facet_labels)))

    n_ref <- nrow(reference_grid)
    P <- length(labels)
    if (!summary) {
      total_rows <- as.numeric(P) * as.numeric(n_ref)
      if (!is.numeric(max_output_rows) || length(max_output_rows) != 1L || !is.finite(max_output_rows) || max_output_rows < 1) {
        stop("max_output_rows must be a positive number")
      }
      if (total_rows > max_output_rows) {
        stop(
          "Returning row-level standardized fair scores requires ", format(total_rows, big.mark = ","),
          " rows (P * nrow(reference_grid)), exceeding max_output_rows = ", format(max_output_rows, big.mark = ","),
          ". Use summary = TRUE or raise 'max_output_rows' deliberately."
        )
      }
    }

    st <- stream_bmfrm_standardized_fair_scores(
      object, reference_grid = reference_grid, adjust_facets = resolved_adjust,
      weights = ref_w, summary = summary, chunk_size = chunk_size
    )
    if (!summary) {
      out_grid <- reference_grid[rep(seq_len(n_ref), times = P), , drop = FALSE]
      out_grid$person <- rep(labels, each = n_ref)
      if (".row" %in% names(out_grid)) out_grid$.row <- NULL
      if (has_explicit_weights) out_grid$weight <- rep(ref_w, times = P)
      rownames(out_grid) <- NULL
      # Evaluate row-level identifiability per (person, reference row) rather than copying summary identifiability
      id_rows <- check_fair_score_identifiability(object, out_grid, resolved_adjust)
      warn_identifiability(id_rows)
      out_grid$mu_fair <- st$mu_fair
      out_grid$mu_fair_sd <- st$mu_fair_sd
      out_grid$likelihood_identified <- id_rows$likelihood_identified
      out_grid$identification_basis <- id_rows$identification_basis
      res_row <- tibble::as_tibble(out_grid)
      if (has_explicit_weights || is_hetero_pcm) {
        ref_grid_meta <- reference_grid
        if (".row" %in% names(ref_grid_meta)) ref_grid_meta$.row <- NULL
        rm_list <- list(
          adjust_facets = resolved_adjust,
          standardized = TRUE,
          summary = FALSE,
          score_unit = "native_row_expected_score",
          K_heterogeneous = is_hetero_pcm,
          K_by_level = K_by_level_meta,
          reference_grid = ref_grid_meta,
          reference_weights = ref_w
        )
        if (is_hetero_pcm) {
          rm_list$comparability_note <- paste0(
            "Heterogeneous category counts across '", object$spec$threshold_facet,
            "': each row's mu_fair is an expected score in the native rubric unit (1..K_j) of its ",
            object$spec$threshold_facet, " level (not rescaled to 0-1); downstream aggregation across ",
            "levels should use the common reference grid and reference_weights."
          )
        }
        attr(res_row, "reference_meta") <- rm_list
      }
      return(res_row)
    }

    active_grid <- reference_grid[ref_w > 0, , drop = FALSE]
    id_std <- check_standardized_fair_identifiability(object, labels, active_grid, resolved_adjust)
    warn_identifiability(id_std)

    D <- nrow(st$sum_mu)
    draw_meta <- extract_draw_index_info(object$fit, D)
    fair_draws_mat <- st$sum_mu
    dimnames(fair_draws_mat) <- list(draw_meta$draw_names, labels)

    if (isTRUE(decompose)) {
      actual_st <- stream_bmfrm_fair_scores(
        object, newdata = obs_data, adjust_facets = character(0),
        summary = TRUE, chunk_size = chunk_size
      )
      obs_fair_st <- stream_bmfrm_fair_scores(
        object, newdata = obs_data, adjust_facets = resolved_adjust,
        summary = TRUE, chunk_size = chunk_size
      )
      actual_draws_mat <- matrix(NA_real_, nrow = D, ncol = P, dimnames = list(draw_meta$draw_names, labels))
      obs_fair_draws_mat <- matrix(NA_real_, nrow = D, ncol = P, dimnames = list(draw_meta$draw_names, labels))
      adj_draws_mat <- matrix(NA_real_, nrow = D, ncol = P, dimnames = list(draw_meta$draw_names, labels))
      std_shift_mat <- matrix(NA_real_, nrow = D, ncol = P, dimnames = list(draw_meta$draw_names, labels))
      shrink_mat <- matrix(NA_real_, nrow = D, ncol = P, dimnames = list(draw_meta$draw_names, labels))
    }

    rows <- lapply(seq_along(labels), function(j) {
      label <- labels[j]
      aggregate_draws <- fair_draws_mat[, j]
      obs_idx <- which(as.character(obs_data$person) == label)
      raw_m <- if (length(obs_idx)) mean(obs_data[[object$spec$response]][obs_idx]) else NA_real_
      row_df <- data.frame(
        person = label,
        n_obs = length(obs_idx),
        mean_score = raw_m,
        mean_fair_score = mean(aggregate_draws),
        sd_fair_score = stats::sd(aggregate_draws),
        lower_95 = unname(stats::quantile(aggregate_draws, .025)),
        upper_95 = unname(stats::quantile(aggregate_draws, .975)),
        likelihood_identified = id_std$likelihood_identified[j],
        identification_basis = id_std$identification_basis[j],
        standardized = TRUE,
        stringsAsFactors = FALSE
      )
      if (isTRUE(decompose)) {
        if (actual_st$count_p[j] > 0L && obs_fair_st$count_p[j] > 0L) {
          act_d <- actual_st$sum_mu[, j] / actual_st$count_p[j]
          obs_f_d <- obs_fair_st$sum_mu[, j] / obs_fair_st$count_p[j]
          adj_d <- obs_f_d - act_d
          std_d <- aggregate_draws - obs_f_d
          shr_d <- act_d - raw_m
          actual_draws_mat[, j] <<- act_d
          obs_fair_draws_mat[, j] <<- obs_f_d
          adj_draws_mat[, j] <<- adj_d
          std_shift_mat[, j] <<- std_d
          shrink_mat[, j] <<- shr_d
          row_df$mean_actual_expected <- mean(act_d)
          row_df$sd_actual_expected <- stats::sd(act_d)
          row_df$mean_obs_fair_score <- mean(obs_f_d)
          row_df$adjustment <- mean(adj_d)
          row_df$adjustment_sd <- stats::sd(adj_d)
          row_df$adjustment_lower_95 <- unname(stats::quantile(adj_d, .025))
          row_df$adjustment_upper_95 <- unname(stats::quantile(adj_d, .975))
          row_df$standardization_shift <- mean(std_d)
          row_df$standardization_sd <- stats::sd(std_d)
          row_df$standardization_lower_95 <- unname(stats::quantile(std_d, .025))
          row_df$standardization_upper_95 <- unname(stats::quantile(std_d, .975))
          row_df$shrinkage_and_noise <- mean(shr_d)
          row_df$shrinkage_lower_95 <- unname(stats::quantile(shr_d, .025))
          row_df$shrinkage_upper_95 <- unname(stats::quantile(shr_d, .975))
          row_df$fair_minus_raw <- mean(aggregate_draws) - raw_m
        } else {
          row_df$mean_actual_expected <- NA_real_
          row_df$sd_actual_expected <- NA_real_
          row_df$mean_obs_fair_score <- NA_real_
          row_df$adjustment <- NA_real_
          row_df$adjustment_sd <- NA_real_
          row_df$adjustment_lower_95 <- NA_real_
          row_df$adjustment_upper_95 <- NA_real_
          row_df$standardization_shift <- NA_real_
          row_df$standardization_sd <- NA_real_
          row_df$standardization_lower_95 <- NA_real_
          row_df$standardization_upper_95 <- NA_real_
          row_df$shrinkage_and_noise <- NA_real_
          row_df$shrinkage_lower_95 <- NA_real_
          row_df$shrinkage_upper_95 <- NA_real_
          row_df$fair_minus_raw <- NA_real_
        }
      }
      row_df
    })
    res <- tibble::as_tibble(do.call(rbind, rows))
    if (isTRUE(return_draws)) {
      draw_list <- list(fair = fair_draws_mat)
      if (isTRUE(decompose)) {
        draw_list$actual <- actual_draws_mat
        draw_list$obs_fair <- obs_fair_draws_mat
        draw_list$adjustment <- adj_draws_mat
        draw_list$standardization_shift <- std_shift_mat
        draw_list$shrinkage_and_noise <- shrink_mat
      }
      attr(res, "draws") <- draw_list
    }
    if (isTRUE(return_draws) || has_explicit_weights || is_hetero_pcm) {
      ref_grid_meta <- reference_grid
      if (".row" %in% names(ref_grid_meta)) ref_grid_meta$.row <- NULL
      rm_list <- list(
        adjust_facets = resolved_adjust,
        standardized = TRUE,
        score_unit = "native_weighted_mean",
        K_heterogeneous = is_hetero_pcm,
        K_by_level = K_by_level_meta,
        reference_grid = ref_grid_meta,
        reference_weights = ref_w,
        draw_ids = draw_meta$draw_ids,
        draw_info = draw_meta$draw_info
      )
      if (is_hetero_pcm) {
        rm_list$comparability_note <- paste0(
          "Heterogeneous category counts across '", object$spec$threshold_facet,
          "': scores are weighted means in native rubric units on the common reference configuration ",
          "(not rescaled to 0-1)."
        )
      }
      attr(res, "reference_meta") <- rm_list
    }
    return(res)
  }

  # Observation-level fair scores:
  # Saved Stan mu_fair may ONLY be used when its actual adjustment specification
  # matches the requested adjust_facets (or adjust_facets = NULL on a custom model).
  hier_pcm_on_adjusted <- identical(object$spec$family, "pcm") &&
    (object$spec$threshold_facet %in% resolved_adjust) &&
    fit_has_variable(object$fit, "tau_common")
  stan_uses_common_tau <- !is.null(object$stan_code) &&
    grepl("eta_fair\\s*-\\s*tau_common\\[k-1\\]", object$stan_code)

  if (isTRUE(object$spec$custom_stan)) {
    contract <- get_eval_contract(object)
    if (use_default_adjust && has_stan_mu_fair) {
      effective_adjust <- if (!is.null(stan_adjusted)) stan_adjusted else resolved_adjust
      can_use_stan <- TRUE
      if (!is.null(stan_adjusted) && !setequal(stan_adjusted, resolved_adjust)) {
        message("Using custom Stan built-in mu_fair (adjusts: ",
                paste(stan_adjusted, collapse = ", "), ").")
      } else if (identical(contract$source, "custom_unverified")) {
        message("Using custom Stan built-in mu_fair (adjustment semantics are defined by the custom Stan model).")
      }
    } else {
      effective_adjust <- resolved_adjust
      can_use_stan <- has_stan_mu_fair && !is.null(stan_adjusted) && setequal(stan_adjusted, resolved_adjust)
    }
  } else {
    effective_adjust <- resolved_adjust
    can_use_stan <- has_stan_mu_fair &&
      !is.null(stan_adjusted) &&
      setequal(stan_adjusted, resolved_adjust) &&
      (!hier_pcm_on_adjusted || stan_uses_common_tau)
  }

  id_obs <- check_fair_score_identifiability(object, obs_data, effective_adjust)
  warn_identifiability(id_obs)

  collapse_person_basis <- function(bases) {
    if (all(is.na(bases))) return(NA_character_)
    if (any(bases == "unidentified_predictor_or_step", na.rm = TRUE)) return("unidentified_predictor_or_step")
    if (any(bases == "hierarchical_step_prior", na.rm = TRUE)) return("hierarchical_step_prior")
    "conditional_likelihood"
  }

  if (can_use_stan) {
    draws <- unclass(object$fit$draws("mu_fair", format = "matrix"))
    index <- as.integer(sub(".*\\[([0-9]+)\\]", "\\1", colnames(draws)))
    draws <- draws[, order(index), drop = FALSE]
    if (ncol(draws) != nrow(obs_data)) stop("mu_fair draws do not match fitted observations")
    obs_data$mu_fair <- colMeans(draws)
    obs_data$mu_fair_sd <- apply(draws, 2L, stats::sd)
    obs_data$likelihood_identified <- id_obs$likelihood_identified
    obs_data$identification_basis <- id_obs$identification_basis
    if (!summary) {
      if (is_hetero_pcm) {
        attr(obs_data, "reference_meta") <- list(
          adjust_facets = effective_adjust,
          standardized = FALSE,
          summary = FALSE,
          score_unit = "native_row_expected_score",
          K_heterogeneous = TRUE,
          K_by_level = K_by_level_meta,
          reference_grid = NULL,
          reference_weights = NULL,
          comparability_note = paste0(
            "Heterogeneous category counts across '", object$spec$threshold_facet,
            "': each row's mu_fair is an expected score in the native rubric unit (1..K_j) of its ",
            object$spec$threshold_facet, " level (not rescaled to 0-1); unstandardized person-level ",
            "aggregation depends on each examinee's assigned levels."
          )
        )
      }
      return(obs_data)
    }

    mu_draws_stan <- NULL
    if (isTRUE(decompose)) {
      if (fit_has_variable(object$fit, "mu")) {
        mu_draws_stan <- unclass(object$fit$draws("mu", format = "matrix"))
        mu_idx <- as.integer(sub(".*\\[([0-9]+)\\]", "\\1", colnames(mu_draws_stan)))
        mu_draws_stan <- mu_draws_stan[, order(mu_idx), drop = FALSE]
        if (ncol(mu_draws_stan) != nrow(obs_data)) stop("mu draws do not match fitted observations")
      } else if (can_evaluate_in_r(object, adjust_facets = character(0))) {
        actual_st <- stream_bmfrm_fair_scores(
          object, newdata = obs_data, adjust_facets = character(0),
          summary = TRUE, chunk_size = chunk_size
        )
      } else {
        stop("decompose = TRUE requires saved Stan 'mu' draws or a verified R-engine model")
      }
    }

    D <- nrow(draws)
    kept_labels <- character(0)
    fair_cols <- list()
    actual_cols <- list()
    adj_cols <- list()
    shrink_cols <- list()

    rows <- lapply(seq_along(labels), function(j) {
      label <- labels[j]
      idx <- which(as.character(obs_data$person) == label)
      if (!length(idx)) return(NULL)
      person_id_ok <- if (all(is.na(id_obs$likelihood_identified[idx]))) NA else all(id_obs$likelihood_identified[idx])
      aggregate_draws <- rowMeans(draws[, idx, drop = FALSE])
      raw_m <- mean(obs_data[[object$spec$response]][idx])
      kept_labels <<- c(kept_labels, label)
      fair_cols[[label]] <<- aggregate_draws
      row_df <- data.frame(
        person = label, n_obs = length(idx),
        mean_score = raw_m,
        mean_fair_score = mean(aggregate_draws), sd_fair_score = stats::sd(aggregate_draws),
        lower_95 = unname(stats::quantile(aggregate_draws, .025)),
        upper_95 = unname(stats::quantile(aggregate_draws, .975)),
        likelihood_identified = person_id_ok,
        identification_basis = collapse_person_basis(id_obs$identification_basis[idx]),
        stringsAsFactors = FALSE
      )
      if (isTRUE(decompose)) {
        act_d <- if (!is.null(mu_draws_stan)) {
          rowMeans(mu_draws_stan[, idx, drop = FALSE])
        } else {
          actual_st$sum_mu[, j] / actual_st$count_p[j]
        }
        adj_d <- aggregate_draws - act_d
        shr_d <- act_d - raw_m
        actual_cols[[label]] <<- act_d
        adj_cols[[label]] <<- adj_d
        shrink_cols[[label]] <<- shr_d
        row_df$mean_actual_expected <- mean(act_d)
        row_df$sd_actual_expected <- stats::sd(act_d)
        row_df$adjustment <- mean(adj_d)
        row_df$adjustment_sd <- stats::sd(adj_d)
        row_df$adjustment_lower_95 <- unname(stats::quantile(adj_d, .025))
        row_df$adjustment_upper_95 <- unname(stats::quantile(adj_d, .975))
        row_df$shrinkage_and_noise <- mean(shr_d)
        row_df$shrinkage_lower_95 <- unname(stats::quantile(shr_d, .025))
        row_df$shrinkage_upper_95 <- unname(stats::quantile(shr_d, .975))
        row_df$fair_minus_raw <- mean(aggregate_draws) - raw_m
      }
      row_df
    })
    res <- tibble::as_tibble(do.call(rbind, rows))
    if (isTRUE(return_draws)) {
      draw_meta_stan <- extract_draw_index_info(object$fit, D)
      set_dn <- function(m) { rownames(m) <- draw_meta_stan$draw_names; m }
      draw_list <- list(fair = set_dn(do.call(cbind, fair_cols)))
      if (isTRUE(decompose)) {
        draw_list$actual <- set_dn(do.call(cbind, actual_cols))
        draw_list$adjustment <- set_dn(do.call(cbind, adj_cols))
        draw_list$shrinkage_and_noise <- set_dn(do.call(cbind, shrink_cols))
      }
      attr(res, "draws") <- draw_list
    }
    if (isTRUE(return_draws) || is_hetero_pcm) {
      draw_meta_stan <- extract_draw_index_info(object$fit, D)
      rm_list <- list(
        adjust_facets = effective_adjust,
        standardized = FALSE,
        score_unit = "native_weighted_mean",
        K_heterogeneous = is_hetero_pcm,
        K_by_level = K_by_level_meta,
        reference_grid = NULL,
        reference_weights = NULL,
        draw_ids = draw_meta_stan$draw_ids,
        draw_info = draw_meta_stan$draw_info
      )
      if (is_hetero_pcm) {
        rm_list$comparability_note <- paste0(
          "Heterogeneous category counts across '", object$spec$threshold_facet,
          "': unstandardized averages depend on each examinee's assigned levels; ",
          "use standardize = TRUE or reference_grid for cross-examinee comparisons."
        )
      }
      attr(res, "reference_meta") <- rm_list
    }
    return(res)
  }

  if (!can_evaluate_in_r(object, adjust_facets = resolved_adjust)) {
    stop(
      "Saved Stan 'mu_fair' does not match the requested adjust_facets (",
      paste(resolved_adjust, collapse = ", "),
      "), and R-side recomputation is not available: ",
      get_eval_contract(object)$reason,
      if (has_stan_mu_fair) " Pass adjust_facets = NULL with standardize = FALSE to use the model's built-in mu_fair, or fit a generated model." else " Fit a generated model."
    )
  }

  # Memory-bounded streaming path when computing fair scores in R
  st <- stream_bmfrm_fair_scores(
    object, newdata = obs_data, adjust_facets = resolved_adjust,
    summary = summary, chunk_size = chunk_size
  )
  if (!summary) {
    obs_data$mu_fair <- st$mu_fair
    obs_data$mu_fair_sd <- st$mu_fair_sd
    obs_data$likelihood_identified <- id_obs$likelihood_identified
    obs_data$identification_basis <- id_obs$identification_basis
    if (is_hetero_pcm) {
      attr(obs_data, "reference_meta") <- list(
        adjust_facets = resolved_adjust,
        standardized = FALSE,
        summary = FALSE,
        score_unit = "native_row_expected_score",
        K_heterogeneous = TRUE,
        K_by_level = K_by_level_meta,
        reference_grid = NULL,
        reference_weights = NULL,
        comparability_note = paste0(
          "Heterogeneous category counts across '", object$spec$threshold_facet,
          "': each row's mu_fair is an expected score in the native rubric unit (1..K_j) of its ",
          object$spec$threshold_facet, " level (not rescaled to 0-1); unstandardized person-level ",
          "aggregation depends on each examinee's assigned levels."
        )
      )
    }
    return(obs_data)
  }

  actual_st <- if (isTRUE(decompose)) {
    stream_bmfrm_fair_scores(
      object, newdata = obs_data, adjust_facets = character(0),
      summary = TRUE, chunk_size = chunk_size
    )
  } else {
    NULL
  }

  fair_cols <- list()
  actual_cols <- list()
  adj_cols <- list()
  shrink_cols <- list()

  rows <- lapply(seq_along(labels), function(j) {
    if (st$count_p[j] == 0L) return(NULL)
    label <- labels[j]
    idx <- which(as.character(obs_data$person) == label)
    person_id_ok <- if (all(is.na(id_obs$likelihood_identified[idx]))) NA else all(id_obs$likelihood_identified[idx])
    aggregate_draws <- st$sum_mu[, j] / st$count_p[j]
    raw_m <- mean(obs_data[[object$spec$response]][idx])
    fair_cols[[label]] <<- aggregate_draws
    row_df <- data.frame(
      person = label, n_obs = length(idx),
      mean_score = raw_m,
      mean_fair_score = mean(aggregate_draws), sd_fair_score = stats::sd(aggregate_draws),
      lower_95 = unname(stats::quantile(aggregate_draws, .025)),
      upper_95 = unname(stats::quantile(aggregate_draws, .975)),
      likelihood_identified = person_id_ok,
      identification_basis = collapse_person_basis(id_obs$identification_basis[idx]),
      stringsAsFactors = FALSE
    )
    if (isTRUE(decompose)) {
      act_d <- actual_st$sum_mu[, j] / actual_st$count_p[j]
      adj_d <- aggregate_draws - act_d
      shr_d <- act_d - raw_m
      actual_cols[[label]] <<- act_d
      adj_cols[[label]] <<- adj_d
      shrink_cols[[label]] <<- shr_d
      row_df$mean_actual_expected <- mean(act_d)
      row_df$sd_actual_expected <- stats::sd(act_d)
      row_df$adjustment <- mean(adj_d)
      row_df$adjustment_sd <- stats::sd(adj_d)
      row_df$adjustment_lower_95 <- unname(stats::quantile(adj_d, .025))
      row_df$adjustment_upper_95 <- unname(stats::quantile(adj_d, .975))
      row_df$shrinkage_and_noise <- mean(shr_d)
      row_df$shrinkage_lower_95 <- unname(stats::quantile(shr_d, .025))
      row_df$shrinkage_upper_95 <- unname(stats::quantile(shr_d, .975))
      row_df$fair_minus_raw <- mean(aggregate_draws) - raw_m
    }
    row_df
  })
  res <- tibble::as_tibble(do.call(rbind, rows))
  if (isTRUE(return_draws)) {
    draw_meta_stream <- extract_draw_index_info(object$fit, nrow(st$sum_mu))
    set_dn <- function(m) { rownames(m) <- draw_meta_stream$draw_names; m }
    draw_list <- list(fair = set_dn(do.call(cbind, fair_cols)))
    if (isTRUE(decompose)) {
      draw_list$actual <- set_dn(do.call(cbind, actual_cols))
      draw_list$adjustment <- set_dn(do.call(cbind, adj_cols))
      draw_list$shrinkage_and_noise <- set_dn(do.call(cbind, shrink_cols))
    }
    attr(res, "draws") <- draw_list
  }
  if (isTRUE(return_draws) || is_hetero_pcm) {
    draw_meta_stream <- extract_draw_index_info(object$fit, nrow(st$sum_mu))
    rm_list <- list(
      adjust_facets = resolved_adjust,
      standardized = FALSE,
      score_unit = "native_weighted_mean",
      K_heterogeneous = is_hetero_pcm,
      K_by_level = K_by_level_meta,
      reference_grid = NULL,
      reference_weights = NULL,
      draw_ids = draw_meta_stream$draw_ids,
      draw_info = draw_meta_stream$draw_info
    )
    if (is_hetero_pcm) {
      rm_list$comparability_note <- paste0(
        "Heterogeneous category counts across '", object$spec$threshold_facet,
        "': unstandardized averages depend on each examinee's assigned levels; ",
        "use standardize = TRUE or reference_grid for cross-examinee comparisons."
      )
    }
    attr(res, "reference_meta") <- rm_list
  }
  res
}

#' Extract facet fit statistics
#' @param object bmfrm_fit object
#' @param facet name of facet
#' @param prob credible interval probability (default 0.95)
#' @param chunk_size Optional observation chunk size for streaming evaluation.
#' @param ... unused
#' @export
facet_fit <- function(object, facet, prob = 0.95, chunk_size = NULL, ...) {
  UseMethod("facet_fit")
}

#' @export
facet_fit.bmfrm_fit <- function(object, facet, prob = 0.95, chunk_size = NULL, ...) {
  if (!facet %in% object$spec$facets_main) {
    stop("Facet '", facet, "' not found in model. Available facets: ",
         paste(object$spec$facets_main, collapse = ", "))
  }

  stan_data <- object$data_info$stan_data
  facet_id <- stan_data[[paste0(facet, "_id")]]
  x <- stan_data$X

  # Only dispatch to the R recomputation engine when the explicit evaluation
  # contract verifies R-engine compatibility AND all required parameters exist.
  can_compute_exact <- can_evaluate_in_r(object)
  if (can_compute_exact) {
    fit_stats <- stream_bmfrm_facet_fit(object, facet = facet, prob = prob, chunk_size = chunk_size)
    fit_stats$df <- pmax(fit_stats$n_obs - 1, 1)
  } else {
    if (!fit_has_variable(object$fit, "mu") || !fit_has_variable(object$fit, "sigma2")) {
      stop(
        "Cannot compute facet fit statistics: ", get_eval_contract(object)$reason,
        " Save 'mu' and 'sigma2' in generated quantities or fit a generated model."
      )
    }
    mu_mat <- extract_ordered_draws(object$fit, "mu", length(x))
    sig_mat <- extract_ordered_draws(object$fit, "sigma2", length(x))
    fit_stats <- facet_fit_from_matrices(
      x = x,
      facet_id = facet_id,
      mu_mat = mu_mat,
      sig_mat = sig_mat,
      kurt_mat = NULL,
      model_name = object$spec$model_name,
      prob = prob
    )
    # Add approximate chi-square Wilson-Hilferty ZSTD when kurtosis is unavailable
    fit_stats <- add_zstd(fit_stats)
  }

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
#' @param prob_thresh Posterior directional probability threshold for flagging
#'   bias beyond `+/-delta` (default 0.95; e.g., use 0.975 to match a two-sided
#'   95% credible interval excluding `[-delta, delta]`).
#' @param delta Nonnegative practically meaningful effect on the adjacent-logit
#'   scale. Reports `prob_gt_delta`, `prob_lt_neg_delta`, `prob_outside_delta`,
#'   and per-cell Bayesian expected false discovery rate (`bfdr`) among
#'   likelihood-identified cells.
#' @param diagnostics Logical; include probability Monte Carlo standard errors
#'   (MCSE) and numerical decision stability (default `TRUE`). `FALSE` retains
#'   the original summary columns and skips these calculations.
#' @param mcse_multiplier Positive multiplier for the approximate Monte Carlo
#'   error band `probability +/- mcse_multiplier * MCSE` (default 2). This does
#'   not change `prob_thresh`, `delta`, or `flag_bias`.
#' @return A tibble with original labels, posterior summaries, directional and
#'   practical-effect probabilities, `bfdr`, cell observation counts, and a
#'   likelihood identification indicator when a design check is available.
#'   Unidentified cells retain posterior summaries but receive flag
#'   `not_likelihood_identified`.
#'   With diagnostics, additional columns are `mcse_prob_gt_delta`,
#'   `mcse_prob_lt_neg_delta`, `mcse_prob_outside_delta`, `prob_thresh`,
#'   `mcse_multiplier`, `decision_stable`, and `decision_status`.
#' @details Probability MCSE uses [posterior::mcse_mean()] on binary indicators
#'   in the original iteration-by-chain layout, accounting for autocorrelation.
#'   Missing chain structure, fewer than four iterations, or an indicator
#'   constant within any chain gives unavailable MCSE (`NA`), never zero by
#'   assumption. At least two chains and passing fit-level sampler checks are
#'   required to assess decision stability.
#'
#'   `decision_status` is `stable_positive_bias` or `stable_negative_bias` when
#'   the corresponding probability's lower Monte Carlo bound is strictly above
#'   `prob_thresh`; `stable_no_clear_bias` when both upper bounds are strictly
#'   below it; and `boundary_sensitive` otherwise. Only the relevant direction
#'   needs finite MCSE for a positive/negative decision; both are required for
#'   no clear bias. Unavailable diagnostics are reported as `mcse_unavailable`,
#'   `chain_structure_unavailable`, `fit_diagnostics_unavailable`, or
#'   `fit_diagnostics_failed`. `decision_stable` is `TRUE` for stable statuses,
#'   `FALSE` for boundary sensitivity, and `NA` when assessment is unavailable.
#'
#'   Stability concerns numerical error in the directional probability rule;
#'   inspect `likelihood_identified` separately even when stability is `TRUE`.
#'   These approximate error bands are not credible intervals for the bias and
#'   do not certify BFDR ranking/selection, equivalence, or model adequacy.
#'   `stable_no_clear_bias` does not establish absence of a meaningful effect.
#'   For boundary sensitivity, inspect mixing and consider longer sampling.
#'   For unavailable MCSE, inspect chain traces and threshold crossings; do not
#'   treat an estimated probability of 0 or 1 as exact certainty.
#' @param ... unused
#' @export
summarise_bias <- function(object, facet = NULL, prob_thresh = 0.95, delta = 0,
                           diagnostics = TRUE, mcse_multiplier = 2, ...) {
  UseMethod("summarise_bias")
}

#' @export
summarise_bias.bmfrm_fit <- function(object, facet = NULL, prob_thresh = 0.95, delta = 0,
                                    diagnostics = TRUE, mcse_multiplier = 2, ...) {
  validate_bias_diagnostics(diagnostics, mcse_multiplier)
  if (!is.numeric(delta) || length(delta) != 1L || !is.finite(delta) || delta < 0) stop("delta must be one finite nonnegative number")
  if (!is.numeric(prob_thresh) || length(prob_thresh) != 1L || !is.finite(prob_thresh) || prob_thresh <= .5 || prob_thresh > 1) stop("prob_thresh must be in (0.5, 1]")
  
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
  chained <- if (diagnostics) {
    tryCatch(object$fit$draws(name, format = "draws_array"), error = function(e) NULL)
  } else NULL
  if (length(dim(chained)) == 3L) {
    draws <- matrix(chained, nrow = prod(dim(chained)[1:2]),
                    dimnames = list(NULL, dimnames(chained)[[3L]]))
  } else {
    # Older adapters may return only matrices. Never invent their chains.
    chained <- NULL
    draws <- unclass(object$fit$draws(name, format = "matrix"))
  }
  cells <- expand.grid(index1 = seq_len(n[1L]), index2 = seq_len(n[2L]))
  variable <- paste0(name, "[", cells$index1, ",", cells$index2, "]")
  if (!all(variable %in% colnames(draws))) stop("Generated interaction draws are missing; summarise_bias_ir() supports legacy raw fits with a parameter named bias")
  draws <- draws[, variable, drop = FALSE]
  ci <- apply(draws, 2L, stats::quantile, probs = c(.025, .975))
  result <- tibble::tibble(facet = facet, facet1 = f[1L], facet2 = f[2L],
    level1 = labels[[1L]][cells$index1], level2 = labels[[2L]][cells$index2],
    index1 = cells$index1, index2 = cells$index2, variable = variable,
    bias_mean = unname(colMeans(draws)), bias_sd = unname(apply(draws, 2L, stats::sd)),
    lower_95 = unname(ci[1L, ]), upper_95 = unname(ci[2L, ]),
    prob_pos = unname(colMeans(draws > 0)), prob_neg = unname(colMeans(draws < 0)), delta = delta,
    prob_gt_delta = unname(colMeans(draws > delta)), prob_lt_neg_delta = unname(colMeans(draws < -delta)),
    prob_outside_delta = unname(colMeans(abs(draws) > delta)))
  # Single consistent probability criterion shared with summarise_bias_ir()
  result$flag_bias <- classify_bias_flags(result$prob_gt_delta, result$prob_lt_neg_delta, prob_thresh = prob_thresh)
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
  # Tie-invariant Bayesian FDR across likelihood-identified cells:
  # All cells sharing the same posterior error probability (within 1e-12) receive
  # the cumulative mean evaluated at the END of their tie group, so selection
  # never depends on arbitrary cell/label ordering.
  pep <- 1 - pmax(result$prob_gt_delta, result$prob_lt_neg_delta)
  result$bfdr <- NA_real_
  elig <- which(is.na(result$likelihood_identified) | result$likelihood_identified)
  if (length(elig)) {
    ord <- elig[order(pep[elig])]
    pep_ord <- pep[ord]
    cum_bfdr <- cumsum(pep_ord) / seq_along(pep_ord)
    tie_key <- round(pep_ord, 12)
    # Last index of each tie group gives the exact set-level BFDR when including all tied cells
    end_idx <- stats::ave(seq_along(tie_key), tie_key, FUN = max)
    result$bfdr[ord] <- cum_bfdr[end_idx]
  }
  # Preserve familiar column names for the item/rater use case in either order.
  if (setequal(f, c("item", "rater"))) {
    result$item <- if (f[1L] == "item") result$level1 else result$level2
    result$rater <- if (f[1L] == "rater") result$level1 else result$level2
    result$i <- if (f[1L] == "item") result$index1 else result$index2
    result$r <- if (f[1L] == "rater") result$index1 else result$index2
  }
  if (diagnostics) {
    result <- dplyr::bind_cols(result, bias_probability_diagnostics(
      chained, variable, result$prob_gt_delta, result$prob_lt_neg_delta,
      delta = delta, prob_thresh = prob_thresh, mcse_multiplier = mcse_multiplier,
      convergence = object$convergence
    ))
  }
  result
}
