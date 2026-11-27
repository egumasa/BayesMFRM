#' Inspect a fitted Bayesian many-facet Rasch model
#'
#' Typing a fit in the console shows its formula, rating scale, number of
#' retained ratings, facet sizes, and separate design and MCMC checks.
#' This overview uses stored metadata. Use [tab_bmfrm()] for estimates and
#' credible intervals or [summary.bmfrm_fit()] for parameter diagnostics.
#' @param x A `bmfrm_fit` object.
#' @param ... Reserved for other methods.
#' @param details If `TRUE`, also print sampling settings, the design rank
#'   and guidance for following up on diagnostics.
#' @return The fit, invisibly.
#' @examples
#' \donttest{
#' if (check_cmdstan()) {
#'   source(system.file("examples", "fit-writing.R", package = "BayesMFRM"), local = TRUE)
#'   print(fit)
#'   print(fit, details = TRUE)
#' }
#' }
#' @export
print.bmfrm_fit <- function(x, ..., details = FALSE) {
  pf <- spec_person_facet(x$spec)
  if (!details) {
    m <- bmfrm_model_overview(x, "Model")
    cat("Bayesian Many-Facet Rasch Model\n", m$formula, "\nFamily: ", m$family,
        " | Categories: ", m$categories, " | Ratings: ", m$n_obs,
        "\nFacet levels: ", m$facet_levels, "\n", sep = "")
    if (!identical(pf, "person")) cat("Ability facet: ", pf, "\n", sep = "")
    pooled <- x$hierarchical %||% x$spec$hierarchical
    if (length(pooled)) cat("Hierarchical facets: ", paste(pooled, collapse = ", "), "\n", sep = "")
    cat("Design: ", m$design, "\n", sep = "")
    cv <- x$convergence
    cat("MCMC:", m$mcmc)
    chains <- cv$n_chains %||% x$fit_args$chains
    if (length(chains)) cat(" | chains:", chains)
    if (length(cv$max_rhat) && is.finite(cv$max_rhat)) cat(" | max R-hat:", format(round(cv$max_rhat, 3)))
    if (length(cv$n_divergent)) cat(" | divergences:", cv$n_divergent)
    cat("\n")
    return(invisible(x))
  }
  cat("Bayesian Many-Facet Rasch Model\n")
  cat("================================\n\n")

  cat("Formula:", deparse(x$spec$formula), "\n")
  cat("Family: ", x$spec$family, "\n")
  if (!identical(pf, "person")) cat("Ability facet:", pf, "\n")
  cat("Hierarchical facets:", paste(x$hierarchical %||% if (!is.null(x$model)) pf else character(), collapse = ", "), "\n")
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
  if (is.null(x$design$identifiable) || is.na(x$design$identifiable)) {
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
  if (identical(x$convergence$status, "incomplete")) {
    cat("  Incomplete: too few chains or missing diagnostic values. See summary(fit).\n")
  } else if (is.null(x$convergence$converged) || is.na(x$convergence$converged)) {
    cat("  Unchecked: no convergence status recorded.\n")
  } else if (isTRUE(x$convergence$converged)) {
    cat("  OK: MCMC diagnostics passed (0 divergences, R-hat & ESS within thresholds)\n")
    if (isTRUE(x$convergence$efficiency_warning) || isTRUE(x$convergence$n_max_treedepth > 0)) {
      cat("    Efficiency note:", x$convergence$n_max_treedepth,
          "transition(s) hit max_treedepth (consider increasing max_treedepth).\n")
    }
    if (identical(x$design$identifiable, FALSE)) {
      cat("    Note: Sampling stability does NOT resolve the design rank deficiency above.\n")
    } else {
      cat("    Note: Convergence checks sampling stability only, not model fit to data.\n")
    }
  } else {
    cat("  CHECK: Potential MCMC convergence issues detected\n")
    if (!is.null(x$convergence$max_rhat) && is.finite(x$convergence$max_rhat)) {
      cat("    Max R-hat:", round(x$convergence$max_rhat, 3), "\n")
    }
    if (!is.null(x$convergence$min_ess_bulk) && is.finite(x$convergence$min_ess_bulk)) {
      cat("    Min ESS (bulk / tail):", round(x$convergence$min_ess_bulk, 1), "/",
          if (!is.null(x$convergence$min_ess_tail) && is.finite(x$convergence$min_ess_tail)) round(x$convergence$min_ess_tail, 1) else NA, "\n")
    }
    if (isTRUE(x$convergence$n_divergent > 0)) {
      cat("    Divergent transitions:", x$convergence$n_divergent, "\n")
    }
    if (isTRUE(x$convergence$n_max_treedepth > 0)) {
      cat("    Max treedepth hits (efficiency):", x$convergence$n_max_treedepth, "\n")
    }
    cat("    Resolve sampling issues before interpreting parameter estimates.\n")
  }

  non_p <- setdiff(x$spec$facets_main, pf)
  ex_facet <- if (length(non_p)) non_p[1L] else pf
  cat("\nNext steps (in priority order):\n")
  cat("  1. Design & convergence : summary(fit)\n")
  cat("  2. Calibrated fit (PPP) : facet_fit(fit, \"", ex_facet, "\"), bayesplot::pp_check(fit, group = \"", ex_facet, "\")\n", sep = "")
  cat("  3. Substantive outputs  : fair_scores(fit, standardize = TRUE)")
  if (length(x$spec$facets_bias)) cat(", summarise_bias(fit, delta = 0.25)")
  cat("\n")
  invisible(x)
}

#' Extract a fitted model's formula or observation count
#'
#' @param x,object A `bmfrm_fit` object.
#' @param ... Reserved for other methods.
#' @return `formula()` returns the fitted formula with its environment.
#'   `nobs()` returns the number of individual ratings retained for fitting.
#' @name bmfrm_fit_accessors
NULL

#' @rdname bmfrm_fit_accessors
#' @importFrom stats formula
#' @export
formula.bmfrm_fit <- function(x, ...) x$spec$formula

#' @rdname bmfrm_fit_accessors
#' @importFrom stats nobs
#' @export
nobs.bmfrm_fit <- function(object, ...) object$data_info$stan_data$N

#' Summary method for bmfrm_fit objects
#' @param object bmfrm_fit object
#' @param facets which facets to summarize (default: all main facets)
#' @param probs quantiles for credible intervals
#' @param ... unused
#' @export
summary.bmfrm_fit <- function(object, facets = NULL, probs = c(0.025, 0.5, 0.975), ...) {
  if (!is.null(probs) && (!is.numeric(probs) || any(!is.finite(probs)) || any(probs < 0 | probs > 1))) {
    stop("probs must be numeric probabilities in [0, 1]")
  }
  if (is.null(facets)) {
    facets <- object$spec$facets_main
  }
  pf <- spec_person_facet(object$spec)

  cat("Bayesian Many-Facet Rasch Model Summary\n")
  cat("=======================================\n\n")

  cat("Formula:", deparse(object$spec$formula), "\n")
  cat("Family: ", object$spec$family, "\n")
  if (!identical(pf, "person")) cat("Ability facet:", pf, "\n")
  cat("Hierarchical facets:", paste(object$hierarchical %||% if (!is.null(object$model)) pf else character(), collapse = ", "), "\n")
  cat("Observations:", object$data_info$stan_data$N, "\n")
  if (isTRUE(object$spec$K_heterogeneous) && !is.null(object$spec$K_by_level)) {
    cat("Rating categories: heterogeneous by ", object$spec$threshold_facet, " (",
        paste(sprintf("%s=1..%d", names(object$spec$K_by_level), as.integer(object$spec$K_by_level)), collapse = ", "),
        ")\n\n", sep = "")
  } else {
    cat("Rating categories: 1 to", get_spec_K_max(object$spec), "\n\n")
  }

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

  fmt_conv <- function(v, digits = NULL) {
    if (is.null(v) || length(v) == 0L || !is.numeric(v) || is.na(v)) return("NA")
    if (!is.null(digits)) as.character(round(v, digits)) else as.character(v)
  }
  cat("MCMC Diagnostics:\n")
  cat("  Chains:", object$fit$num_chains(), "\n")
  cat("  Total iterations:", object$fit$num_chains() * object$fit$metadata()$iter_sampling, "\n")
  cat("  Max R-hat:", fmt_conv(object$convergence$max_rhat, 3), "\n")
  cat("  Min ESS (bulk):", fmt_conv(object$convergence$min_ess_bulk), "\n")
  cat("  Min ESS (tail):", fmt_conv(object$convergence$min_ess_tail), "\n")
  cat("  Min E-BFMI:", fmt_conv(object$convergence$min_ebfmi), "\n")
  cat("  Max mean MCSE:", fmt_conv(object$convergence$max_mcse_mean), "\n")
  cat("  Maximum-depth hits (efficiency):", fmt_conv(object$convergence$n_max_treedepth), "\n")
  if (isTRUE(object$convergence$n_divergent > 0)) {
    cat("  Divergent transitions:", object$convergence$n_divergent, "\n")
  }
  if (!isTRUE(object$convergence$converged)) {
    cat("  Warning: Potential convergence issues detected!\n")
  }
  cat("\n")

  cat("Centered step deviations (tau):\n")
  tau_probs <- if (length(probs) == 0L) c(0.025, 0.975) else sort(unique(range(probs)))
  tau_summary <- if (!isTRUE(object$spec$custom_stan) || fit_has_variable(object$fit, "tau")) {
    tryCatch(
      object$fit$summary("tau", "mean", ~stats::quantile(.x, probs = tau_probs)),
      error = function(e) if (isTRUE(object$spec$custom_stan)) NULL else stop(e)
    )
  } else {
    NULL
  }
  if (is.null(tau_summary)) {
    cat("  (Parameter 'tau' not found in custom Stan model)\n\n")
  } else {
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
      lo_col <- min(3L, ncol(tau_summary))
      hi_col <- ncol(tau_summary)
      cat(sprintf("  %s: %6.2f [%6.2f, %6.2f]\n",
                  label, tau_summary$mean[i], tau_summary[[lo_col]][i], tau_summary[[hi_col]][i]))
    }
    cat("\n")
  }

  for (facet in facets) {
    cat(paste0(toupper(substring(facet, 1, 1)), substring(facet, 2)), " Parameters:\n")
    facet_summ <- tryCatch(
      facet_summary(object, facet, probs = probs),
      error = function(e) if (isTRUE(object$spec$custom_stan)) NULL else stop(e)
    )
    if (is.null(facet_summ)) {
      cat("  (Parameter for facet '", facet, "' not found in custom Stan model)\n\n", sep = "")
    } else {
      print(facet_summ, row.names = FALSE)
      cat("\n")
    }
  }

  if (length(object$spec$facets_bias) > 0) {
    cat("Bias Effects (observed, likelihood-identified cells):\n")
    for (bias_facet in object$spec$facets_bias) {
      cat("  ", bias_facet, ":\n")
      bias_summ <- tryCatch(
        summarise_bias(object, bias_facet),
        error = function(e) if (isTRUE(object$spec$custom_stan)) NULL else stop(e)
      )
      if (is.null(bias_summ)) {
        cat("    (Parameter for interaction '", bias_facet,
            "' not found in custom Stan model; use summarise_bias_ir(object$fit) for legacy 'bias' parameters)\n",
            sep = "")
      } else {
        n_flagged <- sum(bias_summ$flag_bias %in% c("positive_bias", "negative_bias"))
        eligible <- !is.na(bias_summ$bfdr)
        cat("    ", n_flagged, " flagged among ", sum(eligible), " eligible cells",
            " (BFDR <= ", bias_summ$fdr_thresh[1L], ", delta = ", bias_summ$delta[1L], " logits)\n", sep = "")
        excluded <- table(bias_summ$flag_bias[!eligible])
        if (length(excluded)) {
          cat("    Withheld: ", paste(names(excluded), excluded, sep = " = ", collapse = ", "), "\n", sep = "")
        }
        stable <- bias_summ$prob_threshold_stable[eligible]
        if (length(stable)) {
          cat("    Single-cell probability threshold (", bias_summ$prob_thresh[1L], "): ", sum(stable %in% TRUE),
              " resolved, ", sum(stable %in% FALSE), " boundary-sensitive, ",
              sum(is.na(stable)), " unavailable. BFDR stability is unassessed.\n", sep = "")
        }
      }
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
  facet_cols_stan <- paste0(object$spec$facets_main, "_id")
  res_reserved <- c(
    "n", "x", "mu_hat", "sigma2_hat", "resid", "weight", "z", "model",
    "mu_draws", "sigma2_draws"
  )
  resolve_res_facet_col <- function(nm, taken = character()) {
    cand <- if (nm %in% res_reserved) paste0("facet_", nm) else nm
    while (cand %in% taken || cand %in% res_reserved) cand <- paste0("facet_", cand)
    cand
  }

  if (fit_has_variable(object$fit, "mu") && fit_has_variable(object$fit, "sigma2")) {
    res_df <- extract_residuals(
      fit = object$fit,
      stan_data = object$data_info$stan_data,
      facet_cols = facet_cols_stan,
      save_draws = save_draws,
      model_name = object$spec$model_name
    )
    out_facet_cols <- character(length(object$spec$facets_main))
    for (i in seq_along(object$spec$facets_main)) {
      orig_name <- object$spec$facets_main[i]
      stan_name <- facet_cols_stan[i]
      col_pos <- match(stan_name, names(res_df))
      target_name <- resolve_res_facet_col(orig_name, setdiff(names(res_df)[-col_pos], res_reserved))
      out_facet_cols[i] <- target_name
      if (!is.na(col_pos)) {
        flab <- object$data_info$facet_labels[[orig_name]]
        if (!is.null(flab) && is.numeric(res_df[[col_pos]])) {
          res_df[[col_pos]] <- flab[as.integer(res_df[[col_pos]])]
        }
        names(res_df)[col_pos] <- target_name
      }
    }
    attr(res_df, "facet_cols") <- out_facet_cols
    return(res_df)
  }

  if (!can_evaluate_in_r(object)) {
    stop(
      "Cannot recompute residuals in R: ", get_eval_contract(object)$reason,
      " Save 'mu' and 'sigma2' in generated quantities or fit a generated model."
    )
  }
  if (!save_draws) {
    st <- stream_bmfrm_residuals(object)
    mu_hat <- st$mu_hat
    sig_hat <- pmax(st$sigma2_hat, 1e-12)
    N <- length(mu_hat)
    q <- NULL
  } else {
    q <- compute_bmfrm_quantities(object)
    N <- ncol(q$mu)
    mu_hat <- unname(colMeans(q$mu))
    sig_hat <- pmax(unname(colMeans(q$sigma2)), 1e-12)
  }
  x_obs <- object$data_info$stan_data$X
  resid_vec <- x_obs - mu_hat
  z_vec <- resid_vec / sqrt(sig_hat)
  res_df <- tibble::tibble(n = seq_len(N), x = x_obs)
  out_facet_cols <- character(length(object$spec$facets_main))
  for (i in seq_along(object$spec$facets_main)) {
    f <- object$spec$facets_main[i]
    idx_f <- object$data_info$stan_data[[paste0(f, "_id")]]
    flab <- object$data_info$facet_labels[[f]]
    target_name <- resolve_res_facet_col(f, names(res_df))
    out_facet_cols[i] <- target_name
    res_df[[target_name]] <- if (!is.null(flab)) flab[as.integer(idx_f)] else idx_f
  }
  res_df$mu_hat <- mu_hat
  res_df$sigma2_hat <- sig_hat
  res_df$resid <- resid_vec
  res_df$weight <- sig_hat
  res_df$z <- z_vec
  res_df$model <- object$spec$model_name %||% "model"
  if (save_draws) {
    res_df$mu_draws <- lapply(seq_len(N), function(j) as.numeric(q$mu[, j]))
    res_df$sigma2_draws <- lapply(seq_len(N), function(j) as.numeric(q$sigma2[, j]))
  }
  attr(res_df, "facet_cols") <- out_facet_cols
  res_df
}

#' Posterior predictive check method for bmfrm_fit objects
#' @seealso [pp_bmfrm()] for prior predictive simulation on a planned design.
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

  ppc_result <- run_ppc(object$fit, y_obs, ndraws = ndraws, seed = seed, group = group, y_rep = y_rep_override)
  if (type == "all") return(ppc_result)
  if (type == "bars") return(ppc_result$ppc_bars)
  list(ppc_mean = ppc_result$ppc_mean, ppc_var = ppc_result$ppc_var)
}

#' Get Stan code from a fitted model object
#' @param object A fitted model object (e.g., `bmfrm_fit`).
#' @param ... Unused.
#' @export
stancode <- function(object, ...) UseMethod("stancode")

#' @rdname stancode
#' @export stancode.bmfrm_fit
#' @export
stancode.bmfrm_fit <- function(object, ...) {
  cat(object$stan_code)
  invisible(object$stan_code)
}

#' Get Stan data from a fitted model object
#' @param object A fitted model object (e.g., `bmfrm_fit`).
#' @param ... Unused.
#' @export
standata <- function(object, ...) UseMethod("standata")

#' @rdname standata
#' @export standata.bmfrm_fit
#' @export
standata.bmfrm_fit <- function(object, ...) object$data_info$stan_data

#' Extract facet-level summaries
#' @param object bmfrm_fit object
#' @param facet name of facet to summarize
#' @param probs quantiles for credible intervals
#' @param ... unused
#' @return A table with facet labels, posterior means, SDs and requested quantiles.
#'   For a rater main effect, positive values indicate greater severity.
#' @examples
#' data(bmfrm_example)
#' # From facet_summary(fit, "rater") in the worked example.
#' bmfrm_example$rater[, c("level", "mean", "q2.5", "q97.5")]
#' \donttest{
#' if (check_cmdstan()) {
#'   source(system.file("examples", "fit-writing.R", package = "BayesMFRM"), local = TRUE)
#'   facet_summary(fit, "rater")
#' }
#' }
#' @seealso [bmfrm()] for the complete fit.
#' @export
facet_summary <- function(object, facet, probs = c(0.025, 0.5, 0.975), ...) {
  UseMethod("facet_summary")
}

#' @export
facet_summary.bmfrm_fit <- function(object, facet, probs = c(0.025, 0.5, 0.975), ...) {
  if (!is.null(probs) && (!is.numeric(probs) || any(!is.finite(probs)) || any(probs < 0 | probs > 1))) {
    stop("probs must be numeric probabilities in [0, 1]")
  }
  probs <- unique(as.numeric(probs))

  param_name <- facet_parameter_name(object, facet)
  param_summary <- if (length(probs) > 0L) {
    object$fit$summary(param_name, "mean", "sd", ~stats::quantile(.x, probs = probs))
  } else {
    object$fit$summary(param_name, "mean", "sd")
  }
  if ("variable" %in% names(param_summary) && nrow(param_summary) > 1L &&
      all(grepl("\\[\\d+\\]$", param_summary$variable))) {
    p_idx <- as.integer(stringr::str_extract(param_summary$variable, "(?<=\\[)\\d+(?=\\])"))
    param_summary <- param_summary[order(p_idx), , drop = FALSE]
  }

  labels <- object$data_info$facet_labels[[facet]]
  if (nrow(param_summary) != length(labels)) {
    stop("Parameter '", param_name, "' has ", nrow(param_summary),
         " entries, expected ", length(labels), " for facet '", facet, "'")
  }

  result <- tibble::tibble(
    facet = facet,
    level = labels,
    index = seq_along(labels),
    mean = unname(param_summary$mean),
    sd = unname(param_summary$sd)
  )
  fmt_q <- function(p) paste0("q", sub("\\.?0+$", "", formatC(p * 100, format = "f", digits = 4)))
  for (i in seq_along(probs)) {
    result[[fmt_q(probs[i])]] <- param_summary[[ncol(param_summary) - length(probs) + i]]
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
#' @param object A `bmfrm_fit` object.
#' @param newdata Optional data.frame of facet combinations to predict for.
#' @param type `"epred"` (default) or `"response"` for expected rubric scores,
#'   `"predict"` for replicated ratings, `"probs"` for category probabilities,
#'   or `"linpred"` / `"link"` for the adjacent-logit predictor before steps.
#' @param adjust_facets Character vector of non-ability main facets to set to 0.
#' @param allow_new_levels Logical; allow unseen facet levels in `newdata`?
#' @param sample_new_levels `"hierarchical"` (default) or `"zero"`.
#' @param summary Logical; return summary tibble (`TRUE`) or draw array (`FALSE`).
#' @param probs Quantile probabilities for `summary = TRUE` (default `c(0.025, 0.975)`).
#' @param ndraws Optional random subset size of posterior draws.
#' @param draw_ids Optional specific posterior draw indices.
#' @param seed Optional integer RNG seed.
#' @param chunk_size Optional observation chunk size.
#' @param max_cells Maximum numeric entries allowed when `summary = FALSE`.
#' @param ... Reserved. Select uncertainty bounds with `probs` and select the
#'   prediction target with `type`.
#' @return A tibble (when `summary = TRUE`) or matrix/3D array (`summary = FALSE`).
#'   An expected rating is on the original rubric scale. A replicated rating
#'   also includes variation from assigning an observed category.
#'   Scalar summaries contain `estimate`, `sd` and the requested quantiles in
#'   input row order. Category probability summaries give one posterior mean
#'   per supported category, with `NA` padding for heterogeneous rubrics.
#'   For `type = "probs", summary = FALSE`, the dimensions are posterior draw,
#'   prediction row and score category. To compare pass probabilities, sum the
#'   passing categories within each draw, then subtract the two rows in that draw.
#' @examples
#' \donttest{
#' if (check_cmdstan()) {
#'   source(system.file("examples", "fit-writing.R", package = "BayesMFRM"), local = TRUE)
#'   # One column per supported category for this single design row.
#'   p <- predict(fit, ratings[1L, ], type = "probs", summary = FALSE)
#'   graphics::barplot(colMeans(p[, 1L, ]), names.arg = 1:4, ylim = c(0, 1),
#'                     xlab = "Rubric category", ylab = "Posterior probability")
#'   profiles <- data.frame(person = "P03", item = "Structure", rater = c("A", "B"))
#'   p <- predict(fit, profiles, type = "probs", summary = FALSE)
#'   pass <- p[, , 3] + p[, , 4]  # Suppose 3 or 4 on Structure counts as a pass.
#'   gap <- 100 * (pass[, 2] - pass[, 1])
#'   stats::quantile(gap, c(0.025, 0.975))  # Either rater could yield the higher probability.
#'   mean(abs(gap) <= 5)  # Probability the difference is within five percentage points.
#' }
#' }
#' @seealso [bmfrm()], [fitted.bmfrm_fit()], [fair_scores()]
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
  rlang::check_dots_empty()
  if (identical(type, "response")) type <- "epred"
  if (identical(type, "link")) type <- "linpred"
  type <- match.arg(type)
  sample_new_levels <- match.arg(sample_new_levels)
  if (is.null(adjust_facets)) adjust_facets <- character(0)
  resolved_adjust <- resolve_adjust_facets(object$spec, adjust_facets = adjust_facets, warn_if_empty = FALSE)
  stream_bmfrm_predictions(
    object = object, newdata = newdata, type = type, adjust_facets = resolved_adjust,
    allow_new_levels = allow_new_levels, sample_new_levels = sample_new_levels,
    summary = summary, probs = probs, ndraws = ndraws, draw_ids = draw_ids,
    seed = seed, chunk_size = chunk_size, max_cells = max_cells
  )
}

#' Fitted expected ratings on the original rubric
#'
#' @param object A fitted `bmfrm_fit` model that supports [predict.bmfrm_fit()].
#' @param ... Reserved for other methods.
#' @return A numeric vector of posterior mean expected scores, one per
#'   retained rating in fitted data order. Observed facet levels are used.
#'   `predict(object)` also reports uncertainty around each expected score.
#' @importFrom stats fitted
#' @export
fitted.bmfrm_fit <- function(object, ...) {
  rlang::check_dots_empty()
  stats::predict(object, type = "epred", probs = numeric())$estimate
}

# Shared internal metadata builder for fair_scores()
build_fair_reference_meta <- function(object, adjust_facets, standardized, summary,
                                      reference_grid = NULL, reference_weights = NULL,
                                      draw_meta = NULL) {
  is_hetero <- isTRUE(object$spec$K_heterogeneous)
  rm_list <- list(
    adjust_facets = adjust_facets,
    standardized = standardized,
    summary = summary,
    score_unit = if (summary) "native_weighted_mean" else "native_row_expected_score",
    K_heterogeneous = is_hetero,
    K_by_level = object$data_info$K_by_level %||% object$spec$K_by_level,
    reference_grid = reference_grid,
    reference_weights = reference_weights
  )
  if (!is.null(draw_meta)) {
    rm_list$draw_ids <- draw_meta$draw_ids
    rm_list$draw_info <- draw_meta$draw_info
  }
  if (is_hetero) {
    tf <- object$spec$threshold_facet
    rm_list$comparability_note <- if (standardized && summary) {
      paste0("Heterogeneous category counts across '", tf,
             "': scores are weighted means in native rubric units on the common reference configuration (not rescaled to 0-1).")
    } else if (standardized) {
      paste0("Heterogeneous category counts across '", tf,
             "': each row's mu_fair is an expected score in the native rubric unit (1..K_j) of its ",
             tf, " level (not rescaled to 0-1); downstream aggregation across levels should use the common reference grid and reference_weights.")
    } else {
      paste0("Heterogeneous category counts across '", tf,
             "': unstandardized averages depend on each examinee's assigned levels; use standardize = TRUE or reference_grid for cross-examinee comparisons.")
    }
  }
  rm_list
}

# Shared internal person-level summary & decomposition builder for fair_scores()
summarise_person_fair_draws <- function(object, pf, labels, fair_draws_by_j, keep_j,
                                        obs_data, id_per_j, basis_per_j,
                                        standardized, decompose, return_draws,
                                        actual_draws_by_j = NULL, obs_fair_draws_by_j = NULL) {
  D <- length(fair_draws_by_j[[keep_j[1L]]])
  draw_meta <- extract_draw_index_info(object$fit, D)
  kept_labels <- labels[keep_j]
  n_kept <- length(keep_j)

  if (isTRUE(return_draws)) {
    alloc_dm <- function() matrix(NA_real_, nrow = D, ncol = n_kept, dimnames = list(draw_meta$draw_names, kept_labels))
    fair_mat <- alloc_dm()
    if (isTRUE(decompose)) {
      act_mat <- alloc_dm()
      adj_mat <- alloc_dm()
      shr_mat <- alloc_dm()
      if (standardized) {
        obs_f_mat <- alloc_dm()
        std_mat <- alloc_dm()
      }
    }
  }

  resp_vec <- obs_data[[object$spec$response]]
  obs_pf <- as.character(obs_data[[pf]])

  rows <- lapply(seq_along(keep_j), function(pos) {
    j <- keep_j[pos]
    label <- labels[j]
    f_d <- fair_draws_by_j[[j]]
    obs_idx <- which(obs_pf == label)
    raw_m <- if (length(obs_idx)) mean(resp_vec[obs_idx]) else NA_real_
    if (isTRUE(return_draws)) fair_mat[, pos] <<- f_d

    row_df <- data.frame(
      n_obs = length(obs_idx),
      mean_score = raw_m,
      mean_fair_score = mean(f_d),
      sd_fair_score = stats::sd(f_d),
      lower_95 = unname(stats::quantile(f_d, 0.025)),
      upper_95 = unname(stats::quantile(f_d, 0.975)),
      likelihood_identified = id_per_j[j],
      identification_basis = basis_per_j[j],
      stringsAsFactors = FALSE
    )
    row_df <- cbind(stats::setNames(data.frame(label, stringsAsFactors = FALSE), pf), row_df)
    if (standardized) row_df$standardized <- TRUE

    if (isTRUE(decompose)) {
      act_d <- actual_draws_by_j[[j]]
      has_obs <- !is.null(act_d) && (!standardized || !is.null(obs_fair_draws_by_j[[j]]))
      if (has_obs) {
        if (standardized) {
          obs_f_d <- obs_fair_draws_by_j[[j]]
          adj_d <- obs_f_d - act_d
          std_d <- f_d - obs_f_d
        } else {
          adj_d <- f_d - act_d
        }
        shr_d <- act_d - raw_m
        if (isTRUE(return_draws)) {
          act_mat[, pos] <<- act_d
          adj_mat[, pos] <<- adj_d
          shr_mat[, pos] <<- shr_d
          if (standardized) {
            obs_f_mat[, pos] <<- obs_f_d
            std_mat[, pos] <<- std_d
          }
        }
        row_df$mean_actual_expected <- mean(act_d)
        row_df$sd_actual_expected <- stats::sd(act_d)
        if (standardized) row_df$mean_obs_fair_score <- mean(obs_f_d)
        row_df$adjustment <- mean(adj_d)
        row_df$adjustment_sd <- stats::sd(adj_d)
        row_df$adjustment_lower_95 <- unname(stats::quantile(adj_d, 0.025))
        row_df$adjustment_upper_95 <- unname(stats::quantile(adj_d, 0.975))
        if (standardized) {
          row_df$standardization_shift <- mean(std_d)
          row_df$standardization_sd <- stats::sd(std_d)
          row_df$standardization_lower_95 <- unname(stats::quantile(std_d, 0.025))
          row_df$standardization_upper_95 <- unname(stats::quantile(std_d, 0.975))
        }
        row_df$shrinkage_and_noise <- mean(shr_d)
        row_df$shrinkage_lower_95 <- unname(stats::quantile(shr_d, 0.025))
        row_df$shrinkage_upper_95 <- unname(stats::quantile(shr_d, 0.975))
        row_df$fair_minus_raw <- mean(f_d) - raw_m
      } else {
        row_df$mean_actual_expected <- NA_real_
        row_df$sd_actual_expected <- NA_real_
        if (standardized) row_df$mean_obs_fair_score <- NA_real_
        row_df$adjustment <- NA_real_
        row_df$adjustment_sd <- NA_real_
        row_df$adjustment_lower_95 <- NA_real_
        row_df$adjustment_upper_95 <- NA_real_
        if (standardized) {
          row_df$standardization_shift <- NA_real_
          row_df$standardization_sd <- NA_real_
          row_df$standardization_lower_95 <- NA_real_
          row_df$standardization_upper_95 <- NA_real_
        }
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
    dl <- list(fair = fair_mat)
    if (isTRUE(decompose)) {
      dl$actual <- act_mat
      if (standardized) dl$obs_fair <- obs_f_mat
      dl$adjustment <- adj_mat
      if (standardized) dl$standardization_shift <- std_mat
      dl$shrinkage_and_noise <- shr_mat
    }
    attr(res, "draws") <- dl
  }
  list(result = res, draw_meta = draw_meta)
}

#' Extract fair scores
#'
#' Computes expected scores on the rating rubric scale after setting adjusted
#' nuisance facets (and interactions involving them) to zero.
#'
#' @section Relationship to FACETS "Fair Average":
#' Distinguish three quantities when comparing across software:
#' \enumerate{
#'   \item **FACETS "Fair Average"** (`E[Y | theta_p, item = 0, rater = 0]`):
#'     evaluates the nonlinear category expectation at the single reference
#'     origin `0` for all non-person facets.
#'   \item **`fair_scores(fit, standardize = FALSE)`**: averages `E[Y_n | theta_p, item_{i(n)}, rater = 0]`
#'     over each examinee's *observed* tasks/items (retaining task difficulty
#'     differences if task assignment is unbalanced).
#'   \item **`fair_scores(fit, standardize = TRUE)`** (recommended for cross-examinee
#'     comparison): averages `E[Y | theta_p, item_g, rater = 0]` across every
#'     cell `g` of the common reference task/item grid (`reference_grid`,
#'     optionally weighted by `reference_weights`). By Jensen's inequality on
#'     the nonlinear ordinal expectation, averaging expected scores across real
#'     item difficulties differs slightly from evaluating a single hypothetical
#'     item at difficulty `0` (FACETS Fair Average).
#' }
#'
#' @param object `bmfrm_fit` object
#' @param summary whether to return examinee-level summaries (`TRUE`) or row-level scores (`FALSE`)
#' @param adjust_facets Optional character vector of nuisance facet names to set
#'   to zero (defaults to `adjust_facets` set in [bmfrm()], or auto-detects
#'   `"rater"` and `"interlocutor"`).
#' @param standardize Logical; evaluate every examinee over the common reference
#'   grid of unadjusted facets?
#' @param reference_grid Optional data.frame of unadjusted facet combinations.
#' @param reference_weights Optional numeric or named-list weights over `reference_grid`.
#' @param decompose Logical; include draw-matched `Fair - Raw` decomposition?
#' @param return_draws Logical; attach `D x P` posterior draw matrices as `attr(res, "draws")`?
#' @param max_grid_cells Maximum auto-expanded reference grid cells (default `10000L`).
#' @param max_output_rows Maximum row-level standardized output rows (default `1e6L`).
#' @param chunk_size Optional observation chunk size for streaming evaluation.
#' @param ... unused
#' @return By default, one row per person. `mean_score` is the observed average,
#'   `mean_fair_score` the posterior mean adjusted score and `lower_95` and
#'   `upper_95` its credible limits. Check `likelihood_identified` before
#'   interpreting an adjustment. The limits describe an expected score.
#' @examples
#' data(bmfrm_example)
#' # These rows were computed with fair_scores(fit, standardize = TRUE).
#' s <- bmfrm_example$fair
#' head(s[, c("person", "mean_score", "mean_fair_score", "lower_95", "upper_95")])
#' with(s, {
#'   plot(mean_score, mean_fair_score, xlim = c(1, 4), ylim = c(1, 4), asp = 1,
#'        xlab = "Observed average", ylab = "Common-grid expected score")
#'   segments(mean_score, lower_95, mean_score, upper_95, col = "grey60")
#'   abline(0, 1, lty = 2)
#'   points(mean_score, mean_fair_score, pch = 16)
#' })
#' # Points above the diagonal have a higher expected score on the common grid.
#' \donttest{
#' if (check_cmdstan()) {
#'   source(system.file("examples", "fit-writing.R", package = "BayesMFRM"), local = TRUE)
#'   head(fair_scores(fit, standardize = TRUE))
#' }
#' }
#' @seealso [bmfrm()] for fitting and scoring the example data.
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
  if (isTRUE(decompose) && !isTRUE(summary)) stop("decompose = TRUE requires summary = TRUE")
  if (isTRUE(return_draws) && !isTRUE(summary)) {
    stop("return_draws = TRUE requires summary = TRUE; use predict(..., summary = FALSE) for row-level draw matrices")
  }

  pf <- spec_person_facet(object$spec)
  has_stan_mu_fair <- fit_has_variable(object$fit, "mu_fair")
  use_default_adjust <- is.null(adjust_facets)
  stan_adjusted <- stan_mu_fair_adjusted_facets(object)
  use_reference <- isTRUE(standardize) || !is.null(reference_grid) || !is.null(reference_weights)
  generic_adjust <- resolve_adjust_facets(
    object$spec,
    adjust_facets = adjust_facets,
    warn_if_empty = !(has_stan_mu_fair && use_default_adjust && !use_reference &&
                      length(setdiff(object$spec$facets_main, pf)) <= 1L) &&
                    !(isTRUE(object$spec$custom_stan) && use_default_adjust && !is.null(stan_adjusted))
  )
  if (isTRUE(object$spec$custom_stan) && use_default_adjust && !is.null(stan_adjusted)) {
    resolved_adjust <- intersect(stan_adjusted, setdiff(object$spec$facets_main, pf))
    if (use_reference && !setequal(resolved_adjust, generic_adjust)) {
      message("Using custom Stan built-in mu_fair (adjusts: ", paste(resolved_adjust, collapse = ", "), ").")
    }
  } else {
    resolved_adjust <- generic_adjust
  }

  fs_row_reserved <- c("mu_fair", "mu_fair_sd", "likelihood_identified", "identification_basis")
  preserve_colliding_facets <- function(df_in, facet_names) {
    for (fn in intersect(facet_names, intersect(names(df_in), fs_row_reserved))) {
      alt <- paste0("facet_", fn)
      while (alt %in% names(df_in) || alt %in% fs_row_reserved) alt <- paste0("facet_", alt)
      df_in[[alt]] <- df_in[[fn]]
    }
    df_in
  }

  is_hetero_pcm <- isTRUE(object$spec$K_heterogeneous)
  if (identical(object$spec$family, "pcm") && object$spec$threshold_facet %in% resolved_adjust) {
    if (is_hetero_pcm) {
      stop("Fair scores are not defined when adjusting threshold_facet '", object$spec$threshold_facet,
           "' in heterogeneous-category PCM because no common reference step profile exists across differing category counts.")
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
        "inspect 'likelihood_identified' and 'identification_basis'.", call. = FALSE
      )
    } else if (any(basis == "hierarchical_step_prior", na.rm = TRUE)) {
      warning(
        "Hierarchical PCM fair scores rely on the hierarchical step prior to separate 'tau_common' from 'tau_dev' ",
        "(not identified by the conditional likelihood alone; see 'likelihood_identified' and 'identification_basis').",
        call. = FALSE
      )
    }
  }

  collapse_person_basis <- function(bases) {
    if (all(is.na(bases))) return(NA_character_)
    if (any(bases == "unidentified_predictor_or_step", na.rm = TRUE)) return("unidentified_predictor_or_step")
    if (any(bases == "hierarchical_step_prior", na.rm = TRUE)) return("hierarchical_step_prior")
    "conditional_likelihood"
  }

  obs_data <- object$data_info$data_clean
  labels <- object$data_info$facet_labels[[pf]]
  P <- length(labels)

  # 1. Standardized / weighted common-reference fair scores
  if (use_reference) {
    if (!can_evaluate_in_r(object, adjust_facets = resolved_adjust)) {
      stop("Cannot compute standardized fair scores in R: ", get_eval_contract(object)$reason)
    }
    keep_facets <- setdiff(object$spec$facets_main, c(pf, resolved_adjust))
    if (is.null(reference_grid)) {
      if (length(keep_facets)) {
        ref_list <- object$data_info$facet_labels[keep_facets]
        n_cells <- prod(as.numeric(lengths(ref_list)))
        if (!is.numeric(max_grid_cells) || length(max_grid_cells) != 1L || !is.finite(max_grid_cells) || max_grid_cells < 1) {
          stop("max_grid_cells must be a positive number")
        }
        if (n_cells > max_grid_cells) {
          stop("Auto-expanded reference grid has ", format(n_cells, big.mark = ","),
               " cells across unadjusted facets (", paste(keep_facets, collapse = ", "),
               "), exceeding max_grid_cells = ", format(max_grid_cells, big.mark = ","),
               ". Supply an explicit 'reference_grid' or raise 'max_grid_cells' deliberately.")
        }
        reference_grid <- expand.grid(ref_list, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
      } else {
        reference_grid <- data.frame(.row = 1L)
      }
    } else {
      if (!is.data.frame(reference_grid) || !nrow(reference_grid)) stop("reference_grid must be a non-empty data.frame")
      missing_cols <- setdiff(keep_facets, names(reference_grid))
      if (length(missing_cols)) {
        stop("reference_grid is missing unadjusted facet column(s): ", paste(missing_cols, collapse = ", "))
      }
    }

    ref_w <- resolve_reference_weights(reference_grid, reference_weights, object$data_info$facet_labels)
    has_explicit_weights <- !is.null(reference_weights) ||
      ("weight" %in% names(reference_grid) && !("weight" %in% names(object$data_info$facet_labels)))

    n_ref <- nrow(reference_grid)
    if (!summary) {
      total_rows <- as.numeric(P) * as.numeric(n_ref)
      if (!is.numeric(max_output_rows) || length(max_output_rows) != 1L || !is.finite(max_output_rows) || max_output_rows < 1) {
        stop("max_output_rows must be a positive number")
      }
      if (total_rows > max_output_rows) {
        stop("Returning row-level standardized fair scores requires ", format(total_rows, big.mark = ","),
             " rows (P * nrow(reference_grid)), exceeding max_output_rows = ", format(max_output_rows, big.mark = ","),
             ". Use summary = TRUE or raise 'max_output_rows' deliberately.")
      }
    }

    st <- stream_bmfrm_standardized_fair_scores(
      object, reference_grid = reference_grid, adjust_facets = resolved_adjust,
      weights = ref_w, summary = summary, chunk_size = chunk_size
    )
    ref_grid_meta <- reference_grid
    if (".row" %in% names(ref_grid_meta)) ref_grid_meta$.row <- NULL

    if (!summary) {
      out_grid <- reference_grid[rep(seq_len(n_ref), times = P), , drop = FALSE]
      out_grid[[pf]] <- rep(labels, each = n_ref)
      if (".row" %in% names(out_grid)) out_grid$.row <- NULL
      rownames(out_grid) <- NULL
      id_rows <- check_fair_score_identifiability(object, out_grid, resolved_adjust)
      warn_identifiability(id_rows)
      out_grid <- preserve_colliding_facets(out_grid, keep_facets)
      if (has_explicit_weights) {
        w_col <- if ("weight" %in% keep_facets) "reference_weight" else "weight"
        out_grid[[w_col]] <- rep(ref_w, times = P)
      }
      out_grid$mu_fair <- st$mu_fair
      out_grid$mu_fair_sd <- st$mu_fair_sd
      out_grid$likelihood_identified <- id_rows$likelihood_identified
      out_grid$identification_basis <- id_rows$identification_basis
      res_row <- tibble::as_tibble(out_grid)
      attr(res_row, "reference_meta") <- build_fair_reference_meta(
        object, resolved_adjust, standardized = TRUE, summary = FALSE,
        reference_grid = ref_grid_meta, reference_weights = ref_w
      )
      return(res_row)
    }

    active_grid <- reference_grid[ref_w > 0, , drop = FALSE]
    id_std <- check_standardized_fair_identifiability(object, labels, active_grid, resolved_adjust)
    warn_identifiability(id_std)

    fair_by_j <- lapply(seq_len(P), function(j) st$sum_mu[, j])
    act_by_j <- obs_f_by_j <- NULL
    if (isTRUE(decompose)) {
      actual_st <- stream_bmfrm_fair_scores(object, obs_data, character(0), summary = TRUE, chunk_size = chunk_size)
      obs_fair_st <- stream_bmfrm_fair_scores(object, obs_data, resolved_adjust, summary = TRUE, chunk_size = chunk_size)
      act_by_j <- lapply(seq_len(P), function(j) if (actual_st$count_p[j] > 0L) actual_st$sum_mu[, j] / actual_st$count_p[j] else NULL)
      obs_f_by_j <- lapply(seq_len(P), function(j) if (obs_fair_st$count_p[j] > 0L) obs_fair_st$sum_mu[, j] / obs_fair_st$count_p[j] else NULL)
    }

    built <- summarise_person_fair_draws(
      object = object, pf = pf, labels = labels, fair_draws_by_j = fair_by_j, keep_j = seq_len(P),
      obs_data = obs_data, id_per_j = id_std$likelihood_identified, basis_per_j = id_std$identification_basis,
      standardized = TRUE, decompose = decompose, return_draws = return_draws,
      actual_draws_by_j = act_by_j, obs_fair_draws_by_j = obs_f_by_j
    )
    res <- built$result
    attr(res, "reference_meta") <- build_fair_reference_meta(
      object, resolved_adjust, standardized = TRUE, summary = TRUE,
      reference_grid = ref_grid_meta, reference_weights = ref_w, draw_meta = built$draw_meta
    )
    return(res)
  }

  # 2. Unstandardized observation-level / observed-mean fair scores
  hier_pcm_on_adjusted <- identical(object$spec$family, "pcm") &&
    (object$spec$threshold_facet %in% resolved_adjust) &&
    fit_has_variable(object$fit, "tau_common")
  stan_uses_common_tau <- !is.null(object$stan_code) &&
    grepl("cum_tau_fair|tau_common\\[k-1\\]", object$stan_code)

  if (isTRUE(object$spec$custom_stan)) {
    contract <- get_eval_contract(object)
    if (use_default_adjust && has_stan_mu_fair) {
      effective_adjust <- if (!is.null(stan_adjusted)) stan_adjusted else resolved_adjust
      can_use_stan <- TRUE
      if (!is.null(stan_adjusted) && !setequal(stan_adjusted, generic_adjust)) {
        message("Using custom Stan built-in mu_fair (adjusts: ", paste(stan_adjusted, collapse = ", "), ").")
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

  obs_pf <- as.character(obs_data[[pf]])
  id_per_j <- vapply(labels, function(lb) {
    idx <- which(obs_pf == lb)
    if (!length(idx) || all(is.na(id_obs$likelihood_identified[idx]))) NA else all(id_obs$likelihood_identified[idx])
  }, logical(1))
  basis_per_j <- vapply(labels, function(lb) {
    idx <- which(obs_pf == lb)
    if (!length(idx)) NA_character_ else collapse_person_basis(id_obs$identification_basis[idx])
  }, character(1))

  if (can_use_stan) {
    draws <- unclass(object$fit$draws("mu_fair", format = "matrix"))
    index <- as.integer(sub(".*\\[([0-9]+)\\]", "\\1", colnames(draws)))
    draws <- draws[, order(index), drop = FALSE]
    if (ncol(draws) != nrow(obs_data)) stop("mu_fair draws do not match fitted observations")
    if (!summary) {
      obs_out <- preserve_colliding_facets(obs_data, object$spec$facets_main)
      obs_out$mu_fair <- colMeans(draws)
      obs_out$mu_fair_sd <- apply(draws, 2L, stats::sd)
      obs_out$likelihood_identified <- id_obs$likelihood_identified
      obs_out$identification_basis <- id_obs$identification_basis
      attr(obs_out, "reference_meta") <- build_fair_reference_meta(
        object, effective_adjust, standardized = FALSE, summary = FALSE
      )
      return(obs_out)
    }

    act_by_j <- NULL
    if (isTRUE(decompose)) {
      if (fit_has_variable(object$fit, "mu")) {
        mu_d <- unclass(object$fit$draws("mu", format = "matrix"))
        mu_idx <- as.integer(sub(".*\\[([0-9]+)\\]", "\\1", colnames(mu_d)))
        mu_d <- mu_d[, order(mu_idx), drop = FALSE]
        if (ncol(mu_d) != nrow(obs_data)) stop("mu draws do not match fitted observations")
        act_by_j <- lapply(labels, function(lb) {
          idx <- which(obs_pf == lb)
          if (length(idx)) rowMeans(mu_d[, idx, drop = FALSE]) else NULL
        })
      } else if (can_evaluate_in_r(object, adjust_facets = character(0))) {
        actual_st <- stream_bmfrm_fair_scores(object, obs_data, character(0), summary = TRUE, chunk_size = chunk_size)
        act_by_j <- lapply(seq_len(P), function(j) if (actual_st$count_p[j] > 0L) actual_st$sum_mu[, j] / actual_st$count_p[j] else NULL)
      } else {
        stop("decompose = TRUE requires saved Stan 'mu' draws or a verified R-engine model")
      }
    }

    fair_by_j <- lapply(labels, function(lb) {
      idx <- which(obs_pf == lb)
      if (length(idx)) rowMeans(draws[, idx, drop = FALSE]) else NULL
    })
    keep_j <- which(!vapply(fair_by_j, is.null, logical(1)))
    built <- summarise_person_fair_draws(
      object = object, pf = pf, labels = labels, fair_draws_by_j = fair_by_j, keep_j = keep_j,
      obs_data = obs_data, id_per_j = id_per_j, basis_per_j = basis_per_j,
      standardized = FALSE, decompose = decompose, return_draws = return_draws,
      actual_draws_by_j = act_by_j
    )
    res <- built$result
    attr(res, "reference_meta") <- build_fair_reference_meta(
      object, effective_adjust, standardized = FALSE, summary = TRUE, draw_meta = built$draw_meta
    )
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

  st <- stream_bmfrm_fair_scores(object, obs_data, resolved_adjust, summary = summary, chunk_size = chunk_size)
  if (!summary) {
    obs_out <- preserve_colliding_facets(obs_data, object$spec$facets_main)
    obs_out$mu_fair <- st$mu_fair
    obs_out$mu_fair_sd <- st$mu_fair_sd
    obs_out$likelihood_identified <- id_obs$likelihood_identified
    obs_out$identification_basis <- id_obs$identification_basis
    attr(obs_out, "reference_meta") <- build_fair_reference_meta(
      object, resolved_adjust, standardized = FALSE, summary = FALSE
    )
    return(obs_out)
  }

  act_by_j <- NULL
  if (isTRUE(decompose)) {
    actual_st <- stream_bmfrm_fair_scores(object, obs_data, character(0), summary = TRUE, chunk_size = chunk_size)
    act_by_j <- lapply(seq_len(P), function(j) if (actual_st$count_p[j] > 0L) actual_st$sum_mu[, j] / actual_st$count_p[j] else NULL)
  }
  fair_by_j <- lapply(seq_len(P), function(j) if (st$count_p[j] > 0L) st$sum_mu[, j] / st$count_p[j] else NULL)
  keep_j <- which(st$count_p > 0L)

  built <- summarise_person_fair_draws(
    object = object, pf = pf, labels = labels, fair_draws_by_j = fair_by_j, keep_j = keep_j,
    obs_data = obs_data, id_per_j = id_per_j, basis_per_j = basis_per_j,
    standardized = FALSE, decompose = decompose, return_draws = return_draws,
    actual_draws_by_j = act_by_j
  )
  res <- built$result
  attr(res, "reference_meta") <- build_fair_reference_meta(
    object, resolved_adjust, standardized = FALSE, summary = TRUE, draw_meta = built$draw_meta
  )
  res
}

#' Extract calibrated facet fit statistics (Sinharay 2006 posterior predictive p-values)
#'
#' Computes draw-matched realized and replicated mean-square discrepancies
#' `MSQ(y, psi_d)` and `MSQ(y_rep_d, psi_d)` across posterior draws, returning
#' calibrated posterior predictive p-values (`infit_ppp`, `outfit_ppp`) and
#' posterior predictive reference bounds (`infit_rep_lo..hi`, `outfit_rep_lo..hi`)
#' as primary fit diagnostics.
#'
#' Why PPP is primary: posterior parameter uncertainty inflates realized
#' `(y - mu_d)^2 / sigma2_d` for sparse levels or small `n_obs`, causing raw MSQ
#' rules of thumb and asymptotic Wilson-Hilferty `|ZSTD| > 2` thresholds to emit
#' systematic false alarms. Comparing `MSQ(y, psi_d)` against `MSQ(y_rep_d, psi_d)`
#' on the same posterior draws accounts for exact discrete category distributions,
#' sample size `n_obs`, and parameter uncertainty (`PPP` near `0` indicates
#' underfit/erratic noise; `PPP` near `1` indicates overfit/muted dispersion).
#'
#' @param object `bmfrm_fit` object
#' @param facet name of main facet
#' @param prob credible and predictive reference interval probability (default `0.95`)
#' @param seed Optional integer RNG seed for reproducible posterior predictive
#'   draws when `y_rep` is streamed in R (default `123L`).
#' @param include_zstd Logical; include uncalibrated asymptotic Wilson-Hilferty
#'   `infit_zstd` / `outfit_zstd` columns for legacy descriptive comparison?
#'   Default `FALSE`.
#' @param chunk_size Optional observation chunk size for streaming evaluation.
#' @param ... unused
#' @return One row per facet level. `infit_ppp` and `outfit_ppp` compare
#'   observed dispersion with replicated dispersion. Values near zero suggest
#'   unexpectedly variable ratings. Values near one suggest unusually predictable
#'   ratings. Read these alongside `n_obs` and the replicated reference limits.
#' @examples
#' data(bmfrm_example)
#' # From facet_fit(fit, "rater"). Compare the observed MSQ to its reference limits.
#' bmfrm_example$rater_fit[, c("label", "n_obs", "infit_msq", "infit_rep_lo",
#'                             "infit_rep_hi", "infit_ppp")]
#' \donttest{
#' if (check_cmdstan()) {
#'   source(system.file("examples", "fit-writing.R", package = "BayesMFRM"), local = TRUE)
#'   facet_fit(fit, "rater")[, c("label", "infit_msq", "infit_rep_lo", "infit_rep_hi")]
#' }
#' }
#' @seealso [plot_pathway_map()] for a display of fit and severity,
#'   [facet_separation()] and [facet_contrasts()] for differences among levels.
#' @export
facet_fit <- function(object, facet, prob = 0.95, seed = 123L, include_zstd = FALSE, chunk_size = NULL, ...) {
  UseMethod("facet_fit")
}

#' @export
facet_fit.bmfrm_fit <- function(object, facet, prob = 0.95, seed = 123L, include_zstd = FALSE, chunk_size = NULL, ...) {
  if (!facet %in% object$spec$facets_main) {
    stop("Facet '", facet, "' not found in model. Available facets: ",
         paste(object$spec$facets_main, collapse = ", "))
  }
  if (!is.numeric(prob) || length(prob) != 1L || !is.finite(prob) || prob <= 0 || prob >= 1) {
    stop("prob must be in (0, 1)")
  }

  stan_data <- object$data_info$stan_data
  facet_id <- stan_data[[paste0(facet, "_id")]]
  x <- stan_data$X

  if (can_evaluate_in_r(object)) {
    fit_stats <- stream_bmfrm_facet_fit(
      object, facet = facet, prob = prob, seed = seed,
      include_zstd = include_zstd, chunk_size = chunk_size
    )
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
    y_rep_mat <- if (fit_has_variable(object$fit, "y_rep")) {
      extract_ordered_draws(object$fit, "y_rep", length(x))
    } else {
      NULL
    }
    fit_stats <- facet_fit_from_matrices(
      x = x, facet_id = facet_id, mu_mat = mu_mat, sig_mat = sig_mat,
      kurt_mat = NULL, y_rep_mat = y_rep_mat, model_name = object$spec$model_name,
      prob = prob, include_zstd = FALSE
    )
    fit_stats$df <- pmax(fit_stats$n_obs - 1, 1)
    if (isTRUE(include_zstd)) {
      fit_stats <- add_zstd(fit_stats)
    }
  }

  labels <- object$data_info$facet_labels[[facet]]
  fit_stats$label <- labels[fit_stats$facet_value]

  front_cols <- c(
    "facet_value", "label", "n_obs",
    "infit_ppp", "outfit_ppp",
    "infit_msq", "infit_rep_mean", "outfit_msq", "outfit_rep_mean",
    "infit_lo", "infit_hi", "outfit_lo", "outfit_hi",
    "infit_rep_lo", "infit_rep_hi", "outfit_rep_lo", "outfit_rep_hi"
  )
  front_cols <- intersect(front_cols, names(fit_stats))
  fit_stats[, c(front_cols, setdiff(names(fit_stats), front_cols)), drop = FALSE]
}

#' Summarise two-way interaction (bias) effects with BFDR multiplicity control
#'
#' Summarises posterior interaction cells while accounting for structural
#' likelihood identification, zero-observation cells, and multiple testing.
#'
#' \itemize{
#'   \item **Zero-observation cells (`n_observations == 0`)**: even when sum-to-zero
#'     constraints algebraically determine an unobserved cell from other cells,
#'     no ratings occurred in that cell; it is labeled `"unobserved_cell"` and
#'     excluded from `bfdr` and bias flagging.
#'   \item **Multiplicity control (`flag_method = "bfdr"`, default)**: flags
#'     observed, likelihood-identified cells whose set-level Bayesian expected
#'     false discovery rate (`bfdr`, evaluated tie-invariantly over directional
#'     exceedance beyond `+/-delta`) satisfies `bfdr <= fdr_thresh` (default `0.05`).
#'     Pass `flag_method = "prob"` to use per-cell unadjusted directional
#'     probabilities in either direction beyond `delta`.
#' }
#'
#' @param object `bmfrm_fit` object
#' @param facet Two-way interaction name in either component order (defaults to first interaction).
#' @param prob_thresh Per-cell directional probability threshold in `(0.5, 1]`
#'   (used for `flag_method = "prob"` and directional MCSE stability; default `0.95`).
#' @param delta Nonnegative practical effect threshold on the logit scale (default `0`).
#' @param flag_method Multiplicity rule for `flag_bias`: `"bfdr"` (default; flags
#'   cells with `bfdr <= fdr_thresh`) or `"prob"` (per-cell `prob_thresh`).
#' @param fdr_thresh Target Bayesian false discovery rate in `(0, 0.5)` when
#'   `flag_method = "bfdr"` (default `0.05`).
#' @param ci_level Credible interval probability in `(0, 1)` for `ci_lower` and
#'   `ci_upper` (default `0.95`; e.g., use `2 * prob_thresh - 1` to align two-sided
#'   interval bounds with one-sided `prob_thresh`). `lower_95` and `upper_95` are
#'   also retained for backward compatibility.
#' @param diagnostics Logical; include autocorrelation-aware probability MCSE
#'   and decision stability columns (default `TRUE`).
#' @param mcse_multiplier Positive multiplier for MCSE error bands (default `2`).
#' @param ... unused
#' @return One row per pair of facet levels. `bias_mean`, `ci_lower` and
#'   `ci_upper` describe the posterior interaction in logits. A positive value
#'   raises expected ratings relative to the main effects. Both screening rules
#'   require observations and confirmed likelihood identification. Inspect
#'   `flag_bias` for the selected rule. Withheld cells have `bfdr = NA` and a
#'   reason in `flag_bias`: `unobserved_cell`, `not_likelihood_identified`
#'   or `identification_unchecked`. Their posterior probabilities remain available.
#'   With diagnostics enabled, `prob_threshold_stable` and `prob_threshold_status`
#'   describe Monte Carlo precision relative to the single-cell `prob_thresh`.
#'   `decision_stable` and `decision_status` refer to the selected rule.
#'   For eligible BFDR cells these are `NA` and `bfdr_stability_not_assessed`.
#' @examples
#' data(bmfrm_example)
#' # From summarise_bias(fit, "item:rater", delta = 0.3).
#' b <- bmfrm_example$bias
#' b[, c("level1", "level2", "n_observations", "likelihood_identified",
#'       "bias_mean", "bfdr", "flag_bias", "prob_threshold_status")]
#' plot_bias(b)
#' \donttest{
#' if (check_cmdstan()) {
#'   source(system.file("examples", "fit-writing.R", package = "BayesMFRM"), local = TRUE)
#'   plot_bias(summarise_bias(fit, "item:rater", delta = 0.3))
#' }
#' }
#' @seealso [bmfrm()] for a reproducible fit and [plot_bias()] for the symbols.
#' @export
summarise_bias <- function(object, facet = NULL, prob_thresh = 0.95, delta = 0,
                           flag_method = c("bfdr", "prob"), fdr_thresh = 0.05,
                           ci_level = 0.95, diagnostics = TRUE, mcse_multiplier = 2, ...) {
  UseMethod("summarise_bias")
}

#' @export
summarise_bias.bmfrm_fit <- function(object, facet = NULL, prob_thresh = 0.95, delta = 0,
                                    flag_method = c("bfdr", "prob"), fdr_thresh = 0.05,
                                    ci_level = 0.95, diagnostics = TRUE, mcse_multiplier = 2, ...) {
  flag_method <- match.arg(flag_method)
  validate_bias_diagnostics(diagnostics, mcse_multiplier)
  if (!is.numeric(delta) || length(delta) != 1L || !is.finite(delta) || delta < 0) {
    stop("delta must be one finite nonnegative number")
  }
  if (!is.numeric(prob_thresh) || length(prob_thresh) != 1L || !is.finite(prob_thresh) || prob_thresh <= 0.5 || prob_thresh > 1) {
    stop("prob_thresh must be in (0.5, 1]")
  }
  if (!is.numeric(fdr_thresh) || length(fdr_thresh) != 1L || !is.finite(fdr_thresh) || fdr_thresh <= 0 || fdr_thresh >= 0.5) {
    stop("fdr_thresh must be in (0, 0.5)")
  }
  if (!is.numeric(ci_level) || length(ci_level) != 1L || !is.finite(ci_level) || ci_level <= 0 || ci_level >= 1) {
    stop("ci_level must be in (0, 1)")
  }

  if (length(object$spec$facets_bias) == 0) stop("No bias facets in model")
  if (is.null(facet)) {
    facet <- object$spec$facets_bias[1L]
    message("Using bias facet: ", facet)
  }
  if (!is.character(facet) || length(facet) != 1L || is.na(facet)) stop("facet must identify one interaction")
  facet <- match_interaction_names(facet, object$spec$facets_bias)
  if (!facet %in% object$spec$facets_bias) {
    stop("Bias facet '", facet, "' not found. Available: ", paste(object$spec$facets_bias, collapse = ", "))
  }

  f <- strsplit(facet, ":", fixed = TRUE)[[1L]]
  labels <- object$data_info$facet_labels[f]
  n <- lengths(labels)
  name <- gsub(":", "", facet)
  chained <- if (diagnostics) {
    tryCatch(object$fit$draws(name, format = "draws_array"), error = function(e) NULL)
  } else NULL
  if (length(dim(chained)) == 3L) {
    draws <- matrix(chained, nrow = prod(dim(chained)[1:2]), dimnames = list(NULL, dimnames(chained)[[3L]]))
  } else {
    chained <- NULL
    draws <- unclass(object$fit$draws(name, format = "matrix"))
  }

  cells <- expand.grid(index1 = seq_len(n[1L]), index2 = seq_len(n[2L]))
  variable <- paste0(name, "[", cells$index1, ",", cells$index2, "]")
  if (!all(variable %in% colnames(draws))) {
    stop("Generated interaction draws are missing; summarise_bias_ir() supports legacy raw fits with a parameter named bias")
  }
  draws <- draws[, variable, drop = FALSE]

  ci_alpha <- (1 - ci_level) / 2
  ci_custom <- apply(draws, 2L, stats::quantile, probs = c(ci_alpha, 1 - ci_alpha))
  ci_95 <- if (abs(ci_level - 0.95) < 1e-12) ci_custom else apply(draws, 2L, stats::quantile, probs = c(0.025, 0.975))

  dat <- object$data_info$data_clean
  idx1 <- match(as.character(dat[[f[1L]]]), labels[[1L]])
  idx2 <- match(as.character(dat[[f[2L]]]), labels[[2L]])
  n_obs_cell <- tabulate(idx1 + (idx2 - 1L) * n[1L], prod(n))

  dc <- object$design
  lik_id <- rep(NA, nrow(cells))
  if (is.matrix(dc$effect_null_space) && length(dc$tolerance) == 1L &&
      length(dc$interaction_levels[[facet]]$map_row) == nrow(cells)) {
    rows <- dc$interaction_levels[[facet]]$map_row
    lik_id <- rowSums(dc$effect_null_space[rows, , drop = FALSE]^2) < dc$tolerance^2
  }

  p_gt <- unname(colMeans(draws > delta))
  p_lt <- unname(colMeans(draws < -delta))

  # Tie-invariant Bayesian FDR strictly across observed, likelihood-identified cells
  pep <- 1 - pmax(p_gt, p_lt)
  bfdr_vec <- rep(NA_real_, nrow(cells))
  eligible <- (lik_id %in% TRUE) & n_obs_cell > 0L
  elig <- which(eligible)
  if (length(elig)) {
    ord <- elig[order(pep[elig])]
    pep_ord <- pep[ord]
    cum_bfdr <- cumsum(pep_ord) / seq_along(pep_ord)
    tie_key <- round(pep_ord, 12)
    end_idx <- stats::ave(seq_along(tie_key), tie_key, FUN = max)
    bfdr_vec[ord] <- cum_bfdr[end_idx]
  }

  if (identical(flag_method, "bfdr")) {
    flag_vec <- rep("no_clear_bias", nrow(cells))
    sig <- !is.na(bfdr_vec) & (bfdr_vec <= fdr_thresh + 1e-12) & (pmax(p_gt, p_lt) > 0.5)
    flag_vec[sig & (p_gt >= p_lt)] <- "positive_bias"
    flag_vec[sig & (p_lt > p_gt)] <- "negative_bias"
  } else {
    flag_vec <- classify_bias_flags(p_gt, p_lt, prob_thresh = prob_thresh)
  }
  flag_vec[is.na(lik_id)] <- "identification_unchecked"
  flag_vec[lik_id %in% FALSE] <- "not_likelihood_identified"
  flag_vec[n_obs_cell == 0L] <- "unobserved_cell"

  result <- tibble::tibble(
    facet = facet, facet1 = f[1L], facet2 = f[2L],
    level1 = labels[[1L]][cells$index1], level2 = labels[[2L]][cells$index2],
    index1 = cells$index1, index2 = cells$index2, variable = variable,
    bias_mean = unname(colMeans(draws)), bias_sd = unname(apply(draws, 2L, stats::sd)),
    ci_level = ci_level,
    ci_lower = unname(ci_custom[1L, ]), ci_upper = unname(ci_custom[2L, ]),
    lower_95 = unname(ci_95[1L, ]), upper_95 = unname(ci_95[2L, ]),
    prob_pos = unname(colMeans(draws > 0)), prob_neg = unname(colMeans(draws < 0)), delta = delta,
    prob_gt_delta = p_gt, prob_lt_neg_delta = p_lt,
    prob_outside_delta = unname(colMeans(abs(draws) > delta)),
    flag_method = flag_method, fdr_thresh = fdr_thresh, prob_thresh = prob_thresh,
    flag_bias = flag_vec,
    n_observations = n_obs_cell,
    likelihood_identified = lik_id,
    bfdr = bfdr_vec
  )

  if (setequal(f, c("item", "rater"))) {
    result$item <- if (f[1L] == "item") result$level1 else result$level2
    result$rater <- if (f[1L] == "rater") result$level1 else result$level2
    result$i <- if (f[1L] == "item") result$index1 else result$index2
    result$r <- if (f[1L] == "rater") result$index1 else result$index2
  }
  if (diagnostics) {
    diag_tbl <- bias_probability_diagnostics(
      chained, variable, result$prob_gt_delta, result$prob_lt_neg_delta,
      delta = delta, prob_thresh = prob_thresh, mcse_multiplier = mcse_multiplier,
      convergence = object$convergence
    )
    diag_tbl$prob_thresh <- NULL  # Stored with the rule even when diagnostics are disabled.
    diag_tbl$prob_threshold_stable <- diag_tbl$decision_stable
    diag_tbl$prob_threshold_status <- sub("^stable_.*", "resolved", diag_tbl$decision_status)
    if (flag_method == "bfdr") {
      diag_tbl$decision_stable[] <- NA
      diag_tbl$decision_status[] <- "bfdr_stability_not_assessed"
    }
    diag_tbl$decision_stable[!eligible] <- NA
    diag_tbl$decision_status[!eligible] <- flag_vec[!eligible]
    result <- dplyr::bind_cols(result, diag_tbl)
  }
  result
}
