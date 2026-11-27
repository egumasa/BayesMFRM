# Chunked, memory-bounded on-demand posterior computation engine for BayesMFRM.
# Avoids storing K x (D x N) intermediate matrices and streams observation chunks
# directly into facet-level / person-level accumulators so minimal GQ remains
# lightweight all the way through residuals(), facet_fit(), fair_scores(), and pp_check().

#' Check whether a fitted cmdstanr object contains a given Stan variable
#' @keywords internal
fit_has_variable <- function(fit, var_name) {
  md <- tryCatch(fit$metadata(), error = function(e) NULL)
  if (!is.null(md)) {
    vars <- md$stan_variables %||% md$model_params
    if (!is.null(vars)) {
      base <- unique(sub("\\[.*$", "", vars))
      return(var_name %in% base)
    }
  }
  res <- tryCatch(fit$draws(var_name, format = "matrix"), error = function(e) NULL)
  !is.null(res)
}

#' Resolve the internal post-processing evaluation contract for a model specification
#'
#' Explicitly declares whether the R posterior-quantity engine is semantically
#' verified to match the fitted Stan model's linear predictor and link function,
#' and which facets any built-in Stan `mu_fair` adjusts. Generated models and
#' byte-for-byte verified bundled models under `inst/stan/models/` receive an
#' explicit capability record; arbitrary custom Stan files are marked
#' `"custom_unverified"` so R-side recomputation cannot silently evaluate a
#' different likelihood or linear predictor (e.g., `2 * theta`).
#'
#' @keywords internal
resolve_eval_contract <- function(spec, stan_code = NULL) {
  if (!isTRUE(spec$custom_stan)) {
    return(list(
      source = "generated",
      model_id = NULL,
      r_recompute = TRUE,
      stan_mu_fair_adjust = spec_adjust_facets(spec),
      supported_ops = c("residuals", "facet_fit", "fair_scores", "standardized_fair_scores", "pp_check", "predict", "loo"),
      reason = "Generated standard MFRM (RSM/PCM) with verified R-engine semantics."
    ))
  }

  matched <- match_bundled_stan_model(stan_code)
  if (!is.null(matched)) {
    spec_family <- spec$family %||% "rsm"
    spec_bias <- spec$facets_bias %||% character(0)
    structure_ok <- identical(spec_family, matched$family) &&
      !is.null(spec$facets_main) &&
      setequal(spec$facets_main, matched$facets_main) &&
      length(spec_bias) == length(matched$facets_bias) &&
      setequal(spec_bias, matched$facets_bias)
    if (structure_ok) {
      return(list(
        source = "bundled_verified",
        model_id = matched$model_id,
        param_map = matched$param_map,
        r_recompute = isTRUE(matched$r_recompute),
        stan_mu_fair_adjust = matched$stan_mu_fair_adjust,
        supported_ops = matched$supported_ops,
        reason = if (isTRUE(matched$r_recompute)) {
          paste0("Verified bundled Stan model '", matched$model_id, "' matches the standard additive RSM R engine.")
        } else {
          paste0("Verified bundled Stan model '", matched$model_id,
                 "' uses non-generated parameter names; only Stan generated quantities are supported.")
        }
      ))
    }
  }

  list(
    source = "custom_unverified",
    model_id = NULL,
    r_recompute = FALSE,
    stan_mu_fair_adjust = NULL,
    supported_ops = "stan_gq_only",
    reason = "Unverified custom Stan model: R-side recomputation is disabled because linear predictor and link semantics cannot be verified."
  )
}

#' Retrieve (or lazily resolve) the internal evaluation contract on a bmfrm_fit
#' @keywords internal
get_eval_contract <- function(object) {
  ec <- object$spec$eval_contract
  if (is.list(ec) && !is.null(ec$source)) {
    is_custom <- isTRUE(object$spec$custom_stan)
    consistent <- (!is_custom && identical(ec$source, "generated")) ||
      (is_custom && !identical(ec$source, "generated"))
    if (consistent) {
      if (identical(ec$source, "bundled_verified") && is.null(ec$param_map)) {
        reg <- if (!is.null(ec$model_id)) bundled_stan_model_registry()[[ec$model_id]] else NULL
        if (is.null(reg) && !is.null(object$stan_code)) {
          reg <- match_bundled_stan_model(object$stan_code)
        }
        if (!is.null(reg$param_map)) ec$param_map <- reg$param_map
      }
      return(ec)
    }
  }
  resolve_eval_contract(object$spec, object$stan_code)
}

#' Check whether a bmfrm_fit is explicitly verified for R-engine recomputation
#'
#' Requires BOTH an explicit internal evaluation contract certifying semantic
#' compatibility (`r_recompute = TRUE`, granted only to generated models and
#' byte-for-byte verified bundled models whose parameter conventions match) AND
#' presence of every required Stan parameter in the fit object. Arbitrary custom
#' Stan models always return `FALSE` even if their variable names coincide with
#' generated models.
#'
#' @keywords internal
can_evaluate_in_r <- function(object, adjust_facets = character()) {
  contract <- get_eval_contract(object)
  if (!isTRUE(contract$r_recompute)) return(FALSE)

  spec <- object$spec
  if (is.null(spec$facets_main) || is.null(spec$K)) return(FALSE)
  if (!fit_has_variable(object$fit, "theta")) return(FALSE)
  pf <- spec_person_facet(spec)

  other_facets <- setdiff(spec$facets_main, c(pf, adjust_facets))
  for (f in other_facets) {
    if (!fit_has_variable(object$fit, f)) return(FALSE)
  }
  for (bias_facet in spec$facets_bias) {
    comps <- strsplit(bias_facet, ":", fixed = TRUE)[[1L]]
    if (any(comps %in% adjust_facets)) next
    if (!fit_has_variable(object$fit, gsub(":", "", bias_facet))) return(FALSE)
  }

  is_pcm <- identical(spec$family, "pcm")
  if (is_pcm && (spec$threshold_facet %in% adjust_facets)) {
    return(fit_has_variable(object$fit, "tau_common"))
  }
  fit_has_variable(object$fit, "tau")
}

#' Inspect which non-person main facets are removed by a Stan model's built-in `mu_fair`
#'
#' Reads the verified adjustment specification from the model's internal
#' evaluation contract. Returns `NULL` for unverified custom Stan models rather
#' than guessing from regex matching.
#'
#' @keywords internal
stan_mu_fair_adjusted_facets <- function(object) {
  if (is.null(object$spec$facets_main)) return(NULL)
  get_eval_contract(object)$stan_mu_fair_adjust
}

#' Evaluate conditional-likelihood identifiability of fair-score rows and distinguish
#' hierarchical step-prior assumptions
#'
#' Returns a list with logical vector `likelihood_identified` and character
#' vector `identification_basis`:
#'   - `"conditional_likelihood"`: identified algebraically from the conditional likelihood
#'   - `"hierarchical_step_prior"`: predictor identified, but `tau_common` vs `tau_dev`
#'     decomposition relies on the hierarchical step prior rather than the conditional likelihood
#'   - `"unidentified_predictor_or_step"`: depends on a predictor null-space direction
#'     or an unobserved PCM step level
#'   - `NA_character_`: no design check attached to the fit object
#'
#' @keywords internal
check_fair_score_identifiability <- function(object, eval_df, adjust_facets = character()) {
  n <- nrow(eval_df)
  if (identical(get_eval_contract(object)$source, "custom_unverified")) {
    return(list(
      likelihood_identified = rep(NA, n),
      identification_basis = rep(NA_character_, n)
    ))
  }
  dc <- object$design
  K_max_val <- if (!is.null(object$spec$K)) get_spec_K_max(object$spec) else NA_integer_
  n_free_steps <- if (!is.na(K_max_val)) max(K_max_val - 2L, 0L) else NA_integer_
  has_free_steps <- is.na(n_free_steps) || n_free_steps > 0L
  uses_hier_common_tau <- identical(object$spec$family, "pcm") &&
    has_free_steps &&
    (object$spec$threshold_facet %in% adjust_facets) &&
    fit_has_variable(object$fit, "tau_common")

  if (is.null(dc) || is.null(dc$effect_null_space) || is.null(dc$levels)) {
    if (uses_hier_common_tau) {
      return(list(
        likelihood_identified = rep(FALSE, n),
        identification_basis = rep("hierarchical_step_prior", n)
      ))
    }
    return(list(
      likelihood_identified = rep(NA, n),
      identification_basis = rep(NA_character_, n)
    ))
  }

  null_mat <- dc$effect_null_space
  if (!ncol(null_mat)) {
    pred_ok <- rep(TRUE, n)
  } else {
    tol2 <- (dc$tolerance %||% 1e-10)^2
    row_names <- rownames(null_mat)
    proj <- matrix(0, nrow = n, ncol = ncol(null_mat))
    pf <- spec_person_facet(object$spec)

    p_idx <- match(paste0(pf, ":", as.character(eval_df[[pf]])), row_names)
    if (anyNA(p_idx)) {
      return(list(likelihood_identified = rep(NA, n), identification_basis = rep(NA_character_, n)))
    }
    proj <- proj + null_mat[p_idx, , drop = FALSE]

    keep_facets <- setdiff(object$spec$facets_main, c(pf, adjust_facets))
    for (f in keep_facets) {
      f_idx <- match(paste0(f, ":", as.character(eval_df[[f]])), row_names)
      if (anyNA(f_idx)) {
        return(list(likelihood_identified = rep(NA, n), identification_basis = rep(NA_character_, n)))
      }
      proj <- proj - null_mat[f_idx, , drop = FALSE]
    }

    for (term in object$spec$facets_bias) {
      comps <- strsplit(term, ":", fixed = TRUE)[[1L]]
      if (any(comps %in% adjust_facets)) next
      i1 <- match(as.character(eval_df[[comps[1L]]]), object$data_info$facet_labels[[comps[1L]]])
      i2 <- match(as.character(eval_df[[comps[2L]]]), object$data_info$facet_labels[[comps[2L]]])
      b_idx <- match(paste0(term, "[", i1, ",", i2, "]"), row_names)
      if (anyNA(b_idx)) {
        return(list(likelihood_identified = rep(NA, n), identification_basis = rep(NA_character_, n)))
      }
      proj <- proj + null_mat[b_idx, , drop = FALSE]
    }
    pred_ok <- rowSums(proj^2) < tol2
  }

  # Unobserved PCM step levels on an unadjusted threshold facet are unidentified
  if (identical(object$spec$family, "pcm") &&
      has_free_steps &&
      !(object$spec$threshold_facet %in% adjust_facets) &&
      length(dc$unobserved_step_levels)) {
    tf <- object$spec$threshold_facet
    pred_ok <- pred_ok & !(as.character(eval_df[[tf]]) %in% dc$unobserved_step_levels)
  }

  # For hierarchical PCM with free steps (K > 2) on an adjusted facet, total steps
  # tau = tau_common + tau_dev may be identified, but the split into tau_common (used
  # for fair scores) is NOT identified by the conditional likelihood alone (shifting
  # tau_common by +c and tau_dev by -c leaves the conditional likelihood invariant).
  if (uses_hier_common_tau) {
    basis <- ifelse(pred_ok, "hierarchical_step_prior", "unidentified_predictor_or_step")
    return(list(
      likelihood_identified = rep(FALSE, n),
      identification_basis = basis
    ))
  }

  basis <- ifelse(pred_ok, "conditional_likelihood", "unidentified_predictor_or_step")
  list(likelihood_identified = pred_ok, identification_basis = basis)
}

#' Check standardized fair-score identifiability in O((P + n_ref) * q) without expanding P x n_ref
#' @keywords internal
check_standardized_fair_identifiability <- function(object, person_labels, reference_grid, adjust_facets = character()) {
  P <- length(person_labels)
  n_ref <- nrow(reference_grid)
  if (identical(get_eval_contract(object)$source, "custom_unverified")) {
    return(list(likelihood_identified = rep(NA, P), identification_basis = rep(NA_character_, P)))
  }
  dc <- object$design
  K_max_val <- if (!is.null(object$spec$K)) get_spec_K_max(object$spec) else NA_integer_
  n_free_steps <- if (!is.na(K_max_val)) max(K_max_val - 2L, 0L) else NA_integer_
  has_free_steps <- is.na(n_free_steps) || n_free_steps > 0L
  uses_hier_common_tau <- identical(object$spec$family, "pcm") &&
    has_free_steps &&
    (object$spec$threshold_facet %in% adjust_facets) &&
    fit_has_variable(object$fit, "tau_common")

  if (is.null(dc) || is.null(dc$effect_null_space) || is.null(dc$levels)) {
    if (uses_hier_common_tau) {
      return(list(likelihood_identified = rep(FALSE, P), identification_basis = rep("hierarchical_step_prior", P)))
    }
    return(list(likelihood_identified = rep(NA, P), identification_basis = rep(NA_character_, P)))
  }

  pf <- spec_person_facet(object$spec)
  null_mat <- dc$effect_null_space
  q <- ncol(null_mat)
  if (!q) {
    pred_ok <- rep(TRUE, P)
  } else {
    tol2 <- (dc$tolerance %||% 1e-10)^2
    row_names <- rownames(null_mat)
    p_idx <- match(paste0(pf, ":", as.character(person_labels)), row_names)
    if (anyNA(p_idx)) {
      return(list(likelihood_identified = rep(NA, P), identification_basis = rep(NA_character_, P)))
    }
    V_p <- null_mat[p_idx, , drop = FALSE] # P x q

    # Non-person null-space projection U_r (n_ref x q)
    U_r <- matrix(0, nrow = n_ref, ncol = q)
    keep_facets <- setdiff(object$spec$facets_main, c(pf, adjust_facets))
    for (f in keep_facets) {
      f_idx <- match(paste0(f, ":", as.character(reference_grid[[f]])), row_names)
      if (anyNA(f_idx)) {
        return(list(likelihood_identified = rep(NA, P), identification_basis = rep(NA_character_, P)))
      }
      U_r <- U_r - null_mat[f_idx, , drop = FALSE]
    }
    person_terms <- list()
    for (term in object$spec$facets_bias) {
      comps <- strsplit(term, ":", fixed = TRUE)[[1L]]
      if (any(comps %in% adjust_facets)) next
      if (pf %in% comps) {
        other <- setdiff(comps, pf)
        other_idx <- match(as.character(reference_grid[[other]]), object$data_info$facet_labels[[other]])
        if (anyNA(other_idx)) {
          return(list(likelihood_identified = rep(NA, P), identification_basis = rep(NA_character_, P)))
        }
        person_terms[[length(person_terms) + 1L]] <- list(
          term = term,
          person_first = identical(comps[1L], pf),
          other_idx = other_idx
        )
      } else {
        i1 <- match(as.character(reference_grid[[comps[1L]]]), object$data_info$facet_labels[[comps[1L]]])
        i2 <- match(as.character(reference_grid[[comps[2L]]]), object$data_info$facet_labels[[comps[2L]]])
        b_idx <- match(paste0(term, "[", i1, ",", i2, "]"), row_names)
        if (anyNA(b_idx)) {
          return(list(likelihood_identified = rep(NA, P), identification_basis = rep(NA_character_, P)))
        }
        U_r <- U_r + null_mat[b_idx, , drop = FALSE]
      }
    }

    if (!length(person_terms)) {
      u0 <- U_r[1L, , drop = FALSE]
      ref_spread <- if (n_ref > 1L) max(rowSums(sweep(U_r, 2L, u0, "-")^2)) else 0
      if (ref_spread >= tol2) {
        pred_ok <- rep(FALSE, P)
      } else {
        pred_ok <- rowSums(sweep(V_p, 2L, u0, "+")^2) < tol2
      }
    } else {
      p_label_idx <- match(as.character(person_labels), object$data_info$facet_labels[[pf]])
      if (anyNA(p_label_idx)) {
        return(list(likelihood_identified = rep(NA, P), identification_basis = rep(NA_character_, P)))
      }
      pred_ok <- logical(P)
      for (p in seq_len(P)) {
        W_p <- sweep(U_r, 2L, V_p[p, ], "+")
        for (pt in person_terms) {
          i1 <- if (pt$person_first) rep(p_label_idx[p], n_ref) else pt$other_idx
          i2 <- if (pt$person_first) pt$other_idx else rep(p_label_idx[p], n_ref)
          b_idx <- match(paste0(pt$term, "[", i1, ",", i2, "]"), row_names)
          if (anyNA(b_idx)) {
            return(list(likelihood_identified = rep(NA, P), identification_basis = rep(NA_character_, P)))
          }
          W_p <- W_p + null_mat[b_idx, , drop = FALSE]
        }
        pred_ok[p] <- all(rowSums(W_p^2) < tol2)
      }
    }
  }

  if (identical(object$spec$family, "pcm") &&
      has_free_steps &&
      !(object$spec$threshold_facet %in% adjust_facets) &&
      length(dc$unobserved_step_levels)) {
    tf <- object$spec$threshold_facet
    if (identical(tf, pf)) {
      pred_ok <- pred_ok & !(as.character(person_labels) %in% dc$unobserved_step_levels)
    } else if (any(as.character(reference_grid[[tf]]) %in% dc$unobserved_step_levels)) {
      pred_ok <- rep(FALSE, P)
    }
  }

  if (uses_hier_common_tau) {
    basis <- ifelse(pred_ok, "hierarchical_step_prior", "unidentified_predictor_or_step")
    return(list(likelihood_identified = rep(FALSE, P), identification_basis = basis))
  }
  basis <- ifelse(pred_ok, "conditional_likelihood", "unidentified_predictor_or_step")
  list(likelihood_identified = pred_ok, identification_basis = basis)
}

#' Extract indexed parameter draws ordered strictly by 1..J (or column-major matrix)
#' @keywords internal
extract_ordered_draws <- function(fit, var_name, expected_cols = NULL) {
  mat <- unclass(fit$draws(var_name, format = "matrix"))
  cn <- colnames(mat)
  if (!is.null(cn) && length(cn) > 1L) {
    if (all(grepl("\\[\\d+,\\d+\\]$", cn))) {
      m <- stringr::str_match(cn, "\\[(\\d+),(\\d+)\\]$")
      r <- as.integer(m[, 2L])
      c <- as.integer(m[, 3L])
      ord <- order(c, r)
      mat <- mat[, ord, drop = FALSE]
    } else if (all(grepl("\\[\\d+\\]$", cn))) {
      idx <- as.integer(stringr::str_extract(cn, "(?<=\\[)\\d+(?=\\])"))
      mat <- mat[, order(idx), drop = FALSE]
    }
  }
  if (!is.null(expected_cols) && ncol(mat) != expected_cols) {
    stop("Parameter '", var_name, "' has ", ncol(mat), " columns, expected ", expected_cols)
  }
  mat
}

#' Resolve default nuisance facets for fair score adjustment
#' @keywords internal
resolve_adjust_facets <- function(spec, adjust_facets = NULL, warn_if_empty = TRUE) {
  facets_main <- spec$facets_main
  pf <- spec_person_facet(spec)
  if (is.null(adjust_facets)) {
    if (is.null(facets_main)) {
      return(c("rater", "interlocutor"))
    }
    non_person <- setdiff(facets_main, pf)
    if (isTRUE(spec$adjust_facets_explicit)) {
      return(intersect(spec$adjust_facets %||% character(0), non_person))
    }
    found <- intersect(c("rater", "interlocutor"), non_person)
    if (!length(found) && length(non_person) >= 1L && warn_if_empty) {
      warning(
        "No standard nuisance facet ('rater' or 'interlocutor') found in model facets (",
        paste(non_person, collapse = ", "),
        "). No facet was adjusted; specify 'adjust_facets' explicitly (e.g., adjust_facets = 'judge') ",
        "or pass adjust_facets = character(0) to silence this warning.",
        call. = FALSE
      )
    }
    return(found)
  }
  non_person <- setdiff(facets_main, pf)
  if (!is.character(adjust_facets) || anyNA(adjust_facets)) {
    stop("adjust_facets must be a character vector of facet names")
  }
  if (pf %in% adjust_facets) {
    stop("'", pf, "' cannot be included in adjust_facets")
  }
  unknown <- setdiff(adjust_facets, facets_main)
  if (length(unknown)) {
    stop("Unknown facet(s) in adjust_facets: ", paste(unknown, collapse = ", "),
         ". Available main facets: ", paste(non_person, collapse = ", "))
  }
  unique(adjust_facets)
}

#' Extract compact parameter draw matrices once for chunked evaluation
#' @keywords internal
extract_bmfrm_param_draws <- function(object, adjust_facets = character(), draw_idx = NULL) {
  spec <- object$spec
  pf <- spec_person_facet(spec)
  labels <- object$data_info$facet_labels
  K <- get_spec_K_max(spec)
  family <- spec$family %||% "rsm"
  is_hetero <- isTRUE(spec$K_heterogeneous)

  theta <- extract_ordered_draws(object$fit, "theta", length(labels[[pf]]))
  if (!is.null(draw_idx)) theta <- theta[draw_idx, , drop = FALSE]
  D <- nrow(theta)

  other_facets <- setdiff(spec$facets_main, c(pf, adjust_facets))
  facet_draws <- list()
  for (f in other_facets) {
    fd <- extract_ordered_draws(object$fit, f, length(labels[[f]]))
    if (!is.null(draw_idx)) fd <- fd[draw_idx, , drop = FALSE]
    facet_draws[[f]] <- fd
  }

  bias_draws <- list()
  for (bias_facet in spec$facets_bias) {
    comps <- strsplit(bias_facet, ":", fixed = TRUE)[[1L]]
    if (any(comps %in% adjust_facets)) next
    b_name <- gsub(":", "", bias_facet)
    n1 <- length(labels[[comps[1L]]])
    n2 <- length(labels[[comps[2L]]])
    bd <- extract_ordered_draws(object$fit, b_name, n1 * n2)
    if (!is.null(draw_idx)) bd <- bd[draw_idx, , drop = FALSE]
    bias_draws[[bias_facet]] <- list(f1 = comps[1L], f2 = comps[2L], n1 = n1, draws = bd)
  }

  is_pcm <- identical(family, "pcm")
  use_common_tau <- FALSE
  if (is_pcm && spec$threshold_facet %in% adjust_facets) {
    if (is_hetero) {
      stop(
        "Fair scores / adjusted predictions cannot adjust threshold_facet '",
        spec$threshold_facet, "' in heterogeneous-category PCM because no common ",
        "reference step profile exists across differing category counts."
      )
    }
    if (fit_has_variable(object$fit, "tau_common")) {
      use_common_tau <- TRUE
    } else {
      stop("Fair scores are not defined for non-hierarchical PCM steps on an adjusted facet; ",
           "use hierarchical = 'tau' to define a common reference step profile")
    }
  }

  if (!is_pcm || use_common_tau) {
    tau_var <- if (use_common_tau) "tau_common" else "tau"
    tau_draws <- extract_ordered_draws(object$fit, tau_var, K - 1L)
    if (!is.null(draw_idx)) tau_draws <- tau_draws[draw_idx, , drop = FALSE]
    cum_tau <- matrix(0, nrow = D, ncol = K)
    if (K > 1L) {
      for (k in 2:K) cum_tau[, k] <- cum_tau[, k - 1L] + tau_draws[, k - 1L]
    }
    tau_info <- list(mode = "shared", cum_tau = cum_tau)
  } else if (is_hetero) {
    tf <- spec$threshold_facet
    tf_labels <- labels[[tf]]
    J_tf <- length(tf_labels)
    K_lvl <- as.integer(unname((object$data_info$K_by_level %||% spec$K_by_level)[tf_labels]))
    tau_len <- K_lvl - 1L
    tau_start <- as.integer(c(1L, utils::head(cumsum(tau_len), -1L) + 1L))
    tau_draws <- extract_ordered_draws(object$fit, "tau", sum(tau_len))
    if (!is.null(draw_idx)) tau_draws <- tau_draws[draw_idx, , drop = FALSE]
    cum_tau_by_level <- vector("list", J_tf)
    for (j in seq_len(J_tf)) {
      Kj <- K_lvl[j]
      ct_j <- matrix(0, nrow = D, ncol = Kj)
      if (Kj >= 2L) {
        s <- tau_start[j]
        for (k in 2:Kj) {
          ct_j[, k] <- ct_j[, k - 1L] + tau_draws[, s + k - 2L]
        }
      }
      cum_tau_by_level[[j]] <- ct_j
    }
    tau_info <- list(
      mode = "pcm_ragged",
      threshold_facet = tf,
      K_level = K_lvl,
      K_max = K,
      cum_tau_by_level = cum_tau_by_level
    )
  } else {
    tf <- spec$threshold_facet
    J_tf <- length(labels[[tf]])
    tau_draws <- extract_ordered_draws(object$fit, "tau", J_tf * (K - 1L))
    if (!is.null(draw_idx)) tau_draws <- tau_draws[draw_idx, , drop = FALSE]
    cum_tau_list <- vector("list", K)
    cum_tau_list[[1L]] <- matrix(0, nrow = D, ncol = J_tf)
    for (k in 2:K) {
      step_cols <- seq_len(J_tf) + (k - 2L) * J_tf
      cum_tau_list[[k]] <- cum_tau_list[[k - 1L]] + tau_draws[, step_cols, drop = FALSE]
    }
    tau_info <- list(mode = "pcm", threshold_facet = tf, cum_tau_list = cum_tau_list)
  }

  list(
    D = D, K = K, person_facet = pf, labels = labels,
    theta = theta, facet_draws = facet_draws,
    bias_draws = bias_draws, tau_info = tau_info
  )
}

#' Map data frame facet columns to integer indices using stored facet_labels
#' @keywords internal
map_eval_indices <- function(df, spec, labels, adjust_facets = character()) {
  ids <- list()
  for (f in spec$facets_main) {
    if (f %in% adjust_facets) next
    if (!f %in% names(df)) stop("Column '", f, "' not found in evaluation data")
    vals <- as.character(df[[f]])
    idx <- match(vals, labels[[f]])
    if (anyNA(idx)) stop("Unknown level(s) in facet '", f, "' for posterior quantity evaluation")
    ids[[f]] <- idx
  }
  ids
}

#' Default memory-bounded chunk size (target ~16 MB per D x chunk_size matrix)
#' @keywords internal
default_bmfrm_chunk_size <- function(D, N, chunk_size = NULL) {
  if (!is.null(chunk_size)) {
    if (!is.numeric(chunk_size) || length(chunk_size) != 1L || !is.finite(chunk_size) || chunk_size < 1) {
      stop("chunk_size must be a positive integer")
    }
    return(as.integer(min(N, floor(chunk_size))))
  }
  # 2e6 doubles ~= 16 MB per temporary matrix
  as.integer(max(1L, min(N, floor(2e6 / max(D, 1L)))))
}

#' Vectorized two-pass softmax over valid per-level categories 1..K_j for pcm_ragged
#' @keywords internal
eval_softmax_ragged <- function(eta, tau_info, t_idx, resp_chunk = NULL,
                                compute_mu = TRUE,
                                compute_probs = FALSE,
                                compute_kurtosis = FALSE,
                                compute_log_lik = FALSE,
                                compute_y_rep = FALSE, compute_variance = compute_mu) {
  D <- nrow(eta)
  nc <- ncol(eta)
  K_max <- tau_info$K_max
  K_rows <- tau_info$K_level[t_idx]

  out <- list()
  mu_mat <- if (compute_mu || compute_kurtosis) matrix(NA_real_, nrow = D, ncol = nc) else NULL
  s2_mat <- if (compute_variance) matrix(NA_real_, nrow = D, ncol = nc) else NULL
  probs_arr <- if (compute_probs) array(NA_real_, dim = c(D, nc, K_max)) else NULL
  kurt_mat <- if (compute_kurtosis) matrix(0, nrow = D, ncol = nc) else NULL
  ll_mat <- if (compute_log_lik) matrix(NA_real_, nrow = D, ncol = nc) else NULL
  y_rep_mat <- if (compute_y_rep) matrix(1L, nrow = D, ncol = nc) else NULL
  u_all <- if (compute_y_rep) matrix(stats::runif(D * nc), nrow = D, ncol = nc) else NULL

  for (Kj in unique(K_rows)) {
    sub <- which(K_rows == Kj)
    nc_s <- length(sub)
    eta_s <- eta[, sub, drop = FALSE]
    t_sub <- t_idx[sub]
    u_lv <- unique(t_sub)

    ct_list <- vector("list", Kj)
    ct_list[[1L]] <- matrix(0, nrow = D, ncol = nc_s)
    if (Kj >= 2L) {
      for (k in 2:Kj) {
        m_k <- matrix(0, nrow = D, ncol = nc_s)
        for (lv in u_lv) {
          pos <- which(t_sub == lv)
          m_k[, pos] <- tau_info$cum_tau_by_level[[lv]][, k]
        }
        ct_list[[k]] <- m_k
      }
    }

    max_logit <- matrix(0, nrow = D, ncol = nc_s)
    if (Kj >= 2L) {
      for (k in 2:Kj) {
        lk <- (k - 1) * eta_s - ct_list[[k]]
        max_logit <- pmax(max_logit, lk)
      }
    }

    e1 <- exp(-max_logit)
    denom <- e1
    num_mu <- if (compute_mu || compute_kurtosis) e1 else NULL
    num_ex2 <- if (compute_variance) e1 else NULL
    if (Kj >= 2L) {
      for (k in 2:Kj) {
        ek <- exp((k - 1) * eta_s - ct_list[[k]] - max_logit)
        denom <- denom + ek
        if (compute_mu || compute_kurtosis) num_mu <- num_mu + k * ek
        if (compute_variance) num_ex2 <- num_ex2 + (k * k) * ek
      }
    }
    mu_s <- if (compute_mu || compute_kurtosis) num_mu / denom else NULL
    if (compute_mu) mu_mat[, sub] <- mu_s
    if (compute_variance) s2_mat[, sub] <- pmax(num_ex2 / denom - mu_s * mu_s, 1e-12)

    if (compute_probs || compute_kurtosis || compute_log_lik || compute_y_rep) {
      log_denom <- if (compute_log_lik) log(denom) + max_logit else NULL
      u_s <- if (compute_y_rep) u_all[, sub, drop = FALSE] else NULL
      yr_s <- if (compute_y_rep) matrix(1L, nrow = D, ncol = nc_s) else NULL
      cum_p <- if (compute_y_rep) matrix(0, nrow = D, ncol = nc_s) else NULL
      kurt_s <- if (compute_kurtosis) matrix(0, nrow = D, ncol = nc_s) else NULL
      resp_s <- if (compute_log_lik) resp_chunk[sub] else NULL

      for (k in seq_len(Kj)) {
        lk <- if (k == 1L) matrix(0, nrow = D, ncol = nc_s) else (k - 1) * eta_s - ct_list[[k]]
        if (compute_probs || compute_kurtosis || compute_y_rep) {
          pk <- exp(lk - max_logit) / denom
          if (compute_probs) probs_arr[, sub, k] <- pk
          if (compute_kurtosis) kurt_s <- kurt_s + ((k - mu_s)^4) * pk
          if (compute_y_rep && k < Kj) {
            cum_p <- cum_p + pk
            yr_s <- yr_s + (u_s > cum_p)
          }
        }
        if (compute_log_lik) {
          cols <- which(resp_s == k)
          if (length(cols)) {
            ll_mat[, sub[cols]] <- lk[, cols, drop = FALSE] - log_denom[, cols, drop = FALSE]
          }
        }
      }
      if (compute_kurtosis) kurt_mat[, sub] <- kurt_s
      if (compute_y_rep) y_rep_mat[, sub] <- yr_s
    }
  }

  if (compute_mu) out$mu <- mu_mat
  if (compute_variance) out$sigma2 <- s2_mat
  if (compute_probs) out$probs <- probs_arr
  if (compute_kurtosis) out$kurt <- kurt_mat
  if (compute_log_lik) out$log_lik <- ll_mat
  if (compute_y_rep) {
    storage.mode(y_rep_mat) <- "integer"
    out$y_rep <- y_rep_mat
  }
  out
}

#' Evaluate posterior moments, linear predictors, probabilities, or draws for one chunk
#' @keywords internal
eval_bmfrm_chunk <- function(params, ids, rows, resp_chunk = NULL,
                             compute_mu = TRUE,
                             compute_linpred = FALSE,
                             compute_probs = FALSE,
                             compute_kurtosis = FALSE,
                             compute_log_lik = FALSE,
                             compute_y_rep = FALSE) {
  pf <- params$person_facet %||% "person"

  eta <- params$theta[, ids[[pf]][rows], drop = FALSE]
  for (f in names(params$facet_draws)) {
    eta <- eta - params$facet_draws[[f]][, ids[[f]][rows], drop = FALSE]
  }
  for (b in params$bias_draws) {
    cell_idx <- if (!is.null(b$row_cell_idx)) {
      b$row_cell_idx[rows]
    } else {
      ids[[b$f1]][rows] + (ids[[b$f2]][rows] - 1L) * b$n1
    }
    eta <- eta + b$draws[, cell_idx, drop = FALSE]
  }

  tau_info <- params$tau_info
  t_idx <- if (tau_info$mode != "shared") ids[[tau_info$threshold_facet]][rows]
  eval_bmfrm_categories(eta, params$K, tau_info, t_idx, resp_chunk,
                        compute_mu, compute_linpred, compute_probs, compute_kurtosis,
                        compute_log_lik, compute_y_rep)
}

# Per-chunk thresholds for population predictions of unseen uniform-PCM levels.
new_pcm_chunk_thresholds <- function(tau_info, tau_new, t_idx, D, K, seed_salt) {
  nc <- length(t_idx)
  n_k_tf <- tau_new$n_known
  u_pos <- which(t_idx > n_k_tf)
  if (!length(u_pos)) return(NULL)
  k_pos <- which(t_idx <= n_k_tf)
  chunk_cum_tau_pcm <- lapply(seq_len(K), function(k) {
    mat_k <- matrix(0, nrow = D, ncol = nc)
    if (length(k_pos)) mat_k[, k_pos] <- tau_info$cum_tau_list[[k]][, t_idx[k_pos], drop = FALSE]
    mat_k
  })
  u_local <- t_idx[u_pos] - n_k_tf
  if (identical(tau_new$mode, "zero")) {
    for (k in seq_len(K)) {
      chunk_cum_tau_pcm[[k]][, u_pos] <- tau_new$cum_zero[, k]
    }
  } else {
    for (m in unique(u_local)) {
      s_t <- mix_bmfrm_seed(tau_new$seeds[m], seed_salt)
      tau_m <- with_bmfrm_seed(s_t, {
        if (K > 2L) {
          z_dev <- matrix(stats::rnorm(D * (K - 2L)), nrow = D, ncol = K - 2L)
          tau_new$tau_common + (z_dev %*% t(tau_new$B_dev)) * tau_new$sd
        } else {
          tau_new$tau_common
        }
      })
      cum_m <- matrix(0, nrow = D, ncol = K)
      if (K > 1L) {
        for (k in 2:K) cum_m[, k] <- cum_m[, k - 1L] + tau_m[, k - 1L]
      }
      pos_m <- u_pos[u_local == m]
      for (k in seq_len(K)) {
        chunk_cum_tau_pcm[[k]][, pos_m] <- cum_m[, k]
      }
    }
  }
  chunk_cum_tau_pcm
}

# The ordinary and prediction paths build eta differently, then use the same
# category probabilities, moments, log likelihood and replicated-rating draws.
eval_bmfrm_categories <- function(eta, K, tau_info, t_idx, resp_chunk = NULL,
                                  compute_mu = TRUE, compute_linpred = FALSE,
                                  compute_probs = FALSE, compute_kurtosis = FALSE,
                                  compute_log_lik = FALSE, compute_y_rep = FALSE,
                                  tau_new = NULL, seed_salt = NULL,
                                  compute_variance = compute_mu) {
  D <- nrow(eta)
  nc <- ncol(eta)
  out <- list()
  if (compute_linpred) out$eta <- eta

  need_softmax <- compute_mu || compute_probs || compute_kurtosis || compute_log_lik || compute_y_rep
  if (!need_softmax) return(out)

  if (identical(tau_info$mode, "pcm_ragged")) {
    rag <- eval_softmax_ragged(
      eta = eta, tau_info = tau_info, t_idx = t_idx, resp_chunk = resp_chunk,
      compute_mu = compute_mu, compute_probs = compute_probs,
      compute_kurtosis = compute_kurtosis, compute_log_lik = compute_log_lik,
      compute_y_rep = compute_y_rep, compute_variance = compute_variance
    )
    for (nm in names(rag)) out[[nm]] <- rag[[nm]]
    return(out)
  }
  chunk_cum_tau_pcm <- if (!is.null(tau_new) && tau_info$mode == "pcm") {
    new_pcm_chunk_thresholds(tau_info, tau_new, t_idx, D, K, seed_salt)
  }
  get_cum_tau_k <- function(k) {
    if (tau_info$mode == "shared") {
      matrix(tau_info$cum_tau[, k], nrow = D, ncol = nc)
    } else if (!is.null(chunk_cum_tau_pcm)) {
      chunk_cum_tau_pcm[[k]]
    } else {
      tau_info$cum_tau_list[[k]][, t_idx, drop = FALSE]
    }
  }

  # Pass 1: max_logit only (zero K-length matrix lists retained)
  max_logit <- matrix(0, nrow = D, ncol = nc)
  if (K >= 2L) {
    for (k in 2:K) {
      lk <- (k - 1) * eta - get_cum_tau_k(k)
      max_logit <- pmax(max_logit, lk)
    }
  }

  # Pass 2: denom, E[Y], E[Y^2] in one sweep without storing prob_list
  e1 <- exp(-max_logit)
  denom <- e1
  num_mu <- if (compute_mu || compute_kurtosis) e1 else NULL
  num_ex2 <- if (compute_variance) e1 else NULL
  if (K >= 2L) {
    for (k in 2:K) {
      ek <- exp((k - 1) * eta - get_cum_tau_k(k) - max_logit)
      denom <- denom + ek
      if (compute_mu || compute_kurtosis) num_mu <- num_mu + k * ek
      if (compute_variance) num_ex2 <- num_ex2 + (k * k) * ek
    }
  }
  mu <- if (compute_mu || compute_kurtosis) num_mu / denom else NULL
  if (compute_mu) out$mu <- mu
  if (compute_variance) out$sigma2 <- pmax(num_ex2 / denom - mu * mu, 1e-12)

  # Pass 3 (only if probs, kurtosis, log_lik, or y_rep requested)
  if (compute_probs || compute_kurtosis || compute_log_lik || compute_y_rep) {
    probs_arr <- if (compute_probs) array(NA_real_, dim = c(D, nc, K)) else NULL
    kurt <- if (compute_kurtosis) matrix(0, nrow = D, ncol = nc) else NULL
    ll <- if (compute_log_lik) matrix(NA_real_, nrow = D, ncol = nc) else NULL
    log_denom <- if (compute_log_lik) log(denom) + max_logit else NULL
    u <- if (compute_y_rep) matrix(stats::runif(D * nc), nrow = D, ncol = nc) else NULL
    y_rep <- if (compute_y_rep) matrix(1L, nrow = D, ncol = nc) else NULL
    cum_p <- if (compute_y_rep) matrix(0, nrow = D, ncol = nc) else NULL

    for (k in seq_len(K)) {
      lk <- if (k == 1L) matrix(0, nrow = D, ncol = nc) else (k - 1) * eta - get_cum_tau_k(k)
      if (compute_probs || compute_kurtosis || compute_y_rep) {
        pk <- exp(lk - max_logit) / denom
        if (compute_probs) probs_arr[, , k] <- pk
        if (compute_kurtosis) kurt <- kurt + ((k - mu)^4) * pk
        if (compute_y_rep && k < K) {
          cum_p <- cum_p + pk
          y_rep <- y_rep + (u > cum_p)
        }
      }
      if (compute_log_lik) {
        cols <- which(resp_chunk == k)
        if (length(cols)) {
          ll[, cols] <- lk[, cols, drop = FALSE] - log_denom[, cols, drop = FALSE]
        }
      }
    }
    if (compute_probs) out$probs <- probs_arr
    if (compute_kurtosis) out$kurt <- kurt
    if (compute_log_lik) out$log_lik <- ll
    if (compute_y_rep) {
      storage.mode(y_rep) <- "integer"
      out$y_rep <- y_rep
    }
  }

  out
}

#' Compute observation- or grid-level posterior quantities from parameter draws
#'
#' Evaluates `mu`, `sigma2`, and optional `kurt`, `log_lik`, and `y_rep` in
#' memory-bounded observation chunks using two-pass softmax without retaining
#' K-fold intermediate logit/probability lists.
#'
#' @param object A `bmfrm_fit` object.
#' @param newdata Optional data.frame specifying the evaluation grid.
#' @param adjust_facets Character vector of nuisance facets set to 0.
#' @param draw_idx Optional integer vector of draw indices to subset.
#' @param compute_kurtosis Logical; compute fourth central moments `E[(Y-mu)^4]`.
#' @param compute_log_lik Logical; compute pointwise log-likelihood.
#' @param compute_y_rep Logical; sample replicated ratings `y_rep` (`1..K`).
#' @param seed Optional RNG seed for `y_rep` generation.
#' @param chunk_size Optional observation chunk size (auto-sized ~16 MB if NULL).
#' @return A list with matrices (`D x N`) for requested quantities.
#' @keywords internal
compute_bmfrm_quantities <- function(object,
                                     newdata = NULL,
                                     adjust_facets = character(),
                                     draw_idx = NULL,
                                     compute_kurtosis = FALSE,
                                     compute_log_lik = FALSE,
                                     compute_y_rep = FALSE,
                                     seed = NULL,
                                     chunk_size = NULL) {
  spec <- object$spec
  df <- if (is.null(newdata)) object$data_info$data_clean else newdata
  N <- nrow(df)
  params <- extract_bmfrm_param_draws(object, adjust_facets = adjust_facets, draw_idx = draw_idx)
  ids <- map_eval_indices(df, spec, params$labels, adjust_facets = adjust_facets)
  D <- params$D
  cs <- default_bmfrm_chunk_size(D, N, chunk_size)
  resp <- if (compute_log_lik) as.integer(df[[spec$response]]) else NULL

  mu_all <- matrix(NA_real_, nrow = D, ncol = N)
  sig_all <- matrix(NA_real_, nrow = D, ncol = N)
  kurt_all <- if (compute_kurtosis) matrix(NA_real_, nrow = D, ncol = N) else NULL
  ll_all <- if (compute_log_lik) matrix(NA_real_, nrow = D, ncol = N) else NULL
  yrep_all <- if (compute_y_rep) matrix(NA_integer_, nrow = D, ncol = N) else NULL

  run_chunks <- function() {
    starts <- seq.int(1L, N, by = cs)
    for (s in starts) {
      e <- min(N, s + cs - 1L)
      rows <- s:e
      ch <- eval_bmfrm_chunk(
        params, ids, rows,
        resp_chunk = if (!is.null(resp)) resp[rows] else NULL,
        compute_kurtosis = compute_kurtosis,
        compute_log_lik = compute_log_lik,
        compute_y_rep = compute_y_rep
      )
      mu_all[, rows] <<- ch$mu
      sig_all[, rows] <<- ch$sigma2
      if (compute_kurtosis) kurt_all[, rows] <<- ch$kurt
      if (compute_log_lik) ll_all[, rows] <<- ch$log_lik
      if (compute_y_rep) yrep_all[, rows] <<- ch$y_rep
    }
  }
  if (compute_y_rep && !is.null(seed)) with_bmfrm_seed(seed, run_chunks()) else run_chunks()

  out <- list(mu = mu_all, sigma2 = sig_all)
  if (compute_kurtosis) out$kurt <- kurt_all
  if (compute_log_lik) out$log_lik <- ll_all
  if (compute_y_rep) out$y_rep <- yrep_all
  out
}

#' Streaming observation-level residual point summaries (O(N) memory, no D x N retention)
#' @keywords internal
stream_bmfrm_residuals <- function(object, chunk_size = NULL) {
  df <- object$data_info$data_clean
  N <- nrow(df)
  params <- extract_bmfrm_param_draws(object)
  ids <- map_eval_indices(df, object$spec, params$labels)
  cs <- default_bmfrm_chunk_size(params$D, N, chunk_size)

  mu_hat <- numeric(N)
  sig_hat <- numeric(N)
  starts <- seq.int(1L, N, by = cs)
  for (s in starts) {
    e <- min(N, s + cs - 1L)
    rows <- s:e
    ch <- eval_bmfrm_chunk(params, ids, rows)
    mu_hat[rows] <- unname(colMeans(ch$mu))
    sig_hat[rows] <- unname(colMeans(ch$sigma2))
  }
  list(mu_hat = mu_hat, sigma2_hat = sig_hat)
}

#' Streaming facet-level infit/outfit with Sinharay (2006) posterior predictive p-values
#'
#' Evaluates realized and replicated discrepancies `MSQ(y, psi_d)` vs.
#' `MSQ(y_rep_d, psi_d)` on identical posterior draws in `O(D * J_facet)` memory,
#' returning calibrated posterior predictive p-values (`infit_ppp`, `outfit_ppp`)
#' and posterior predictive reference bounds (`infit_rep_lo..hi`, `outfit_rep_lo..hi`).
#'
#' @keywords internal
stream_bmfrm_facet_fit <- function(object, facet, prob = 0.95, seed = 123L,
                                   include_zstd = FALSE, chunk_size = NULL) {
  df <- object$data_info$data_clean
  N <- nrow(df)
  x <- object$data_info$stan_data$X
  facet_id <- object$data_info$stan_data[[paste0(facet, "_id")]]
  J <- length(object$data_info$facet_labels[[facet]])

  params <- extract_bmfrm_param_draws(object)
  ids <- map_eval_indices(df, object$spec, params$labels)
  D <- params$D
  cs <- default_bmfrm_chunk_size(D, N, chunk_size)

  sum_z2 <- matrix(0, nrow = D, ncol = J)
  sum_wz2 <- matrix(0, nrow = D, ncol = J)
  sum_z2_rep <- matrix(0, nrow = D, ncol = J)
  sum_wz2_rep <- matrix(0, nrow = D, ncol = J)
  sum_w <- matrix(0, nrow = D, ncol = J)
  sum_c_over_w2_m1 <- if (isTRUE(include_zstd)) matrix(0, nrow = D, ncol = J) else NULL
  sum_c_minus_w2 <- if (isTRUE(include_zstd)) matrix(0, nrow = D, ncol = J) else NULL
  n_obs_vec <- integer(J)

  mu_hat <- numeric(N)
  sig_hat <- numeric(N)

  run_fit_stream <- function() {
    starts <- seq.int(1L, N, by = cs)
    for (s in starts) {
      e <- min(N, s + cs - 1L)
      rows <- s:e
      ch <- eval_bmfrm_chunk(
        params, ids, rows,
        compute_kurtosis = isTRUE(include_zstd),
        compute_y_rep = TRUE
      )
      mu_hat[rows] <<- unname(colMeans(ch$mu))
      sig_hat[rows] <<- unname(colMeans(ch$sigma2))

      X_mat <- matrix(x[rows], nrow = D, ncol = length(rows), byrow = TRUE)
      resid2 <- (X_mat - ch$mu)^2
      z2_chunk <- resid2 / ch$sigma2

      resid2_rep <- (ch$y_rep - ch$mu)^2
      z2_rep_chunk <- resid2_rep / ch$sigma2

      if (isTRUE(include_zstd)) {
        w2_chunk <- ch$sigma2 * ch$sigma2
        c_over_w2_m1 <- pmax(ch$kurt / w2_chunk - 1, 0)
        c_minus_w2 <- pmax(ch$kurt - w2_chunk, 0)
      }

      f_chunk <- facet_id[rows]
      for (lv in unique(f_chunk)) {
        loc <- which(f_chunk == lv)
        n_obs_vec[lv] <<- n_obs_vec[lv] + length(loc)
        sum_z2[, lv] <<- sum_z2[, lv] + rowSums(z2_chunk[, loc, drop = FALSE])
        sum_wz2[, lv] <<- sum_wz2[, lv] + rowSums(resid2[, loc, drop = FALSE])
        sum_z2_rep[, lv] <<- sum_z2_rep[, lv] + rowSums(z2_rep_chunk[, loc, drop = FALSE])
        sum_wz2_rep[, lv] <<- sum_wz2_rep[, lv] + rowSums(resid2_rep[, loc, drop = FALSE])
        sum_w[, lv] <<- sum_w[, lv] + rowSums(ch$sigma2[, loc, drop = FALSE])
        if (isTRUE(include_zstd)) {
          sum_c_over_w2_m1[, lv] <<- sum_c_over_w2_m1[, lv] + rowSums(c_over_w2_m1[, loc, drop = FALSE])
          sum_c_minus_w2[, lv] <<- sum_c_minus_w2[, lv] + rowSums(c_minus_w2[, loc, drop = FALSE])
        }
      }
    }
  }
  if (!is.null(seed)) with_bmfrm_seed(seed, run_fit_stream()) else run_fit_stream()

  prob_lo <- (1 - prob) / 2
  prob_hi <- 1 - prob_lo
  tol_tie <- 1e-12
  obs_levels <- which(n_obs_vec > 0L)
  out_rows <- lapply(obs_levels, function(lv) {
    n_obs <- n_obs_vec[lv]
    idx <- which(facet_id == lv)
    outfit_draws <- sum_z2[, lv] / n_obs
    infit_draws <- sum_wz2[, lv] / sum_w[, lv]
    outfit_rep_draws <- sum_z2_rep[, lv] / n_obs
    infit_rep_draws <- sum_wz2_rep[, lv] / sum_w[, lv]

    in_diff <- infit_rep_draws - infit_draws
    out_diff <- outfit_rep_draws - outfit_draws
    infit_ppp <- mean((in_diff > tol_tie) + 0.5 * (abs(in_diff) <= tol_tie))
    outfit_ppp <- mean((out_diff > tol_tie) + 0.5 * (abs(out_diff) <= tol_tie))

    z2_plug <- (x[idx] - mu_hat[idx])^2 / sig_hat[idx]
    w_plug <- sig_hat[idx]

    out <- tibble::tibble(
      facet_value = lv,
      model = object$spec$model_name,
      n_obs = n_obs,
      infit_ppp = infit_ppp,
      outfit_ppp = outfit_ppp,
      infit_msq = mean(infit_draws),
      outfit_msq = mean(outfit_draws),
      infit_lo = unname(stats::quantile(infit_draws, prob_lo)),
      infit_hi = unname(stats::quantile(infit_draws, prob_hi)),
      outfit_lo = unname(stats::quantile(outfit_draws, prob_lo)),
      outfit_hi = unname(stats::quantile(outfit_draws, prob_hi)),
      infit_rep_mean = mean(infit_rep_draws),
      outfit_rep_mean = mean(outfit_rep_draws),
      infit_rep_lo = unname(stats::quantile(infit_rep_draws, prob_lo)),
      infit_rep_hi = unname(stats::quantile(infit_rep_draws, prob_hi)),
      outfit_rep_lo = unname(stats::quantile(outfit_rep_draws, prob_lo)),
      outfit_rep_hi = unname(stats::quantile(outfit_rep_draws, prob_hi)),
      infit_msq_plugin = sum(w_plug * z2_plug, na.rm = TRUE) / sum(w_plug, na.rm = TRUE),
      outfit_msq_plugin = mean(z2_plug, na.rm = TRUE)
    )
    if (isTRUE(include_zstd)) {
      q_out_draws <- sqrt(pmax(sum_c_over_w2_m1[, lv] / (n_obs^2), 1e-12))
      q_in_draws <- sqrt(pmax(sum_c_minus_w2[, lv] / (sum_w[, lv]^2), 1e-12))
      outfit_zstd_draws <- (outfit_draws^(1/3) - 1) * (3 / q_out_draws) + q_out_draws / 3
      infit_zstd_draws <- (infit_draws^(1/3) - 1) * (3 / q_in_draws) + q_in_draws / 3
      out$outfit_q <- mean(q_out_draws)
      out$infit_q <- mean(q_in_draws)
      out$infit_zstd <- mean(infit_zstd_draws)
      out$outfit_zstd <- mean(outfit_zstd_draws)
      out$infit_zstd_lo <- unname(stats::quantile(infit_zstd_draws, prob_lo))
      out$infit_zstd_hi <- unname(stats::quantile(infit_zstd_draws, prob_hi))
      out$outfit_zstd_lo <- unname(stats::quantile(outfit_zstd_draws, prob_lo))
      out$outfit_zstd_hi <- unname(stats::quantile(outfit_zstd_draws, prob_hi))
      out$zstd_method <- "uncalibrated_draw_wilson_hilferty_approx"
    }
    out
  })
  dplyr::bind_rows(out_rows)
}

#' Streaming fair scores (O(D * J_person) memory when summary = TRUE)
#' @keywords internal
stream_bmfrm_fair_scores <- function(object, newdata, adjust_facets,
                                     summary = TRUE, chunk_size = NULL) {
  spec <- object$spec
  params <- extract_bmfrm_param_draws(object, adjust_facets = adjust_facets)
  pf <- params$person_facet %||% "person"
  ids <- map_eval_indices(newdata, spec, params$labels, adjust_facets = adjust_facets)
  D <- params$D
  N <- nrow(newdata)
  P <- length(params$labels[[pf]])
  cs <- default_bmfrm_chunk_size(D, N, chunk_size)

  if (!summary) {
    mu_mean <- numeric(N)
    mu_sd <- numeric(N)
    starts <- seq.int(1L, N, by = cs)
    for (s in starts) {
      e <- min(N, s + cs - 1L)
      rows <- s:e
      ch <- eval_bmfrm_chunk(params, ids, rows)
      mu_mean[rows] <- unname(colMeans(ch$mu))
      mu_sd[rows] <- unname(apply(ch$mu, 2L, stats::sd))
    }
    return(list(mu_fair = mu_mean, mu_fair_sd = mu_sd))
  }

  sum_mu <- matrix(0, nrow = D, ncol = P)
  count_p <- integer(P)
  p_ids <- ids[[pf]]
  starts <- seq.int(1L, N, by = cs)
  for (s in starts) {
    e <- min(N, s + cs - 1L)
    rows <- s:e
    ch <- eval_bmfrm_chunk(params, ids, rows)
    pc <- p_ids[rows]
    for (pv in unique(pc)) {
      loc <- which(pc == pv)
      count_p[pv] <- count_p[pv] + length(loc)
      sum_mu[, pv] <- sum_mu[, pv] + rowSums(ch$mu[, loc, drop = FALSE])
    }
  }

  list(sum_mu = sum_mu, count_p = count_p, labels = params$labels[[pf]])
}

#' Streaming standardized fair scores without expanding a P x n_ref data.frame
#' @keywords internal
stream_bmfrm_standardized_fair_scores <- function(object, reference_grid, adjust_facets,
                                                  weights = NULL,
                                                  summary = TRUE, chunk_size = NULL) {
  spec <- object$spec
  params <- extract_bmfrm_param_draws(object, adjust_facets = adjust_facets)
  pf <- params$person_facet %||% "person"
  D <- params$D
  K <- params$K
  P <- length(params$labels[[pf]])
  n_ref <- nrow(reference_grid)
  if (is.null(weights)) {
    weights <- rep(1 / n_ref, n_ref)
  } else {
    weights <- as.numeric(weights)
    if (length(weights) != n_ref || any(!is.finite(weights)) || any(weights < 0) || sum(weights) <= 0) {
      stop("weights must be a finite non-negative numeric vector of length nrow(reference_grid) with positive sum")
    }
    weights <- weights / sum(weights)
  }
  ref_ids <- map_eval_indices(reference_grid, spec, params$labels, adjust_facets = c(pf, adjust_facets))

  budget <- default_bmfrm_chunk_size(D, P * n_ref, chunk_size)
  ref_chunk_sz <- max(1L, min(n_ref, budget))
  p_chunk_sz <- max(1L, min(P, floor(budget / ref_chunk_sz)))

  sum_mu <- if (summary) matrix(0, nrow = D, ncol = P) else NULL
  mu_mean <- if (!summary) numeric(P * n_ref) else NULL
  mu_sd <- if (!summary) numeric(P * n_ref) else NULL

  tau_info <- params$tau_info
  ref_biases <- Filter(function(b) !identical(b$f1, pf) && !identical(b$f2, pf), params$bias_draws)
  person_biases <- Filter(function(b) identical(b$f1, pf) || identical(b$f2, pf), params$bias_draws)
  ref_starts <- seq.int(1L, n_ref, by = ref_chunk_sz)
  p_starts <- seq.int(1L, P, by = p_chunk_sz)

  for (rs in ref_starts) {
    re <- min(n_ref, rs + ref_chunk_sz - 1L)
    r_rows <- rs:re
    nr <- length(r_rows)
    w_chunk <- weights[r_rows]

    base_eta <- matrix(0, nrow = D, ncol = nr)
    for (f in names(params$facet_draws)) {
      base_eta <- base_eta - params$facet_draws[[f]][, ref_ids[[f]][r_rows], drop = FALSE]
    }
    for (b in ref_biases) {
      cell_idx <- ref_ids[[b$f1]][r_rows] + (ref_ids[[b$f2]][r_rows] - 1L) * b$n1
      base_eta <- base_eta + b$draws[, cell_idx, drop = FALSE]
    }
    is_pcm_mode <- tau_info$mode %in% c("pcm", "pcm_ragged")
    is_person_pcm <- is_pcm_mode && identical(tau_info$threshold_facet, pf)
    t_idx_r <- if (is_pcm_mode && !is_person_pcm) {
      ref_ids[[tau_info$threshold_facet]][r_rows]
    } else {
      NULL
    }

    for (ps in p_starts) {
      pe <- min(P, ps + p_chunk_sz - 1L)
      p_rows <- ps:pe
      np <- length(p_rows)

      p_rep <- rep(p_rows, each = nr)
      r_rep <- rep(seq_len(nr), times = np)
      eta <- params$theta[, p_rep, drop = FALSE] + base_eta[, r_rep, drop = FALSE]
      for (b in person_biases) {
        i1 <- if (identical(b$f1, pf)) p_rep else ref_ids[[b$f1]][r_rows][r_rep]
        i2 <- if (identical(b$f2, pf)) p_rep else ref_ids[[b$f2]][r_rows][r_rep]
        cell_idx <- i1 + (i2 - 1L) * b$n1
        eta <- eta + b$draws[, cell_idx, drop = FALSE]
      }

      t_rep <- if (is_person_pcm) p_rep else t_idx_r[r_rep]
      mu_block <- eval_bmfrm_categories(eta, K, tau_info, t_rep, compute_variance = FALSE)$mu

      if (summary) {
        for (j_loc in seq_len(np)) {
          cols <- ((j_loc - 1L) * nr + 1L):(j_loc * nr)
          sum_mu[, p_rows[j_loc]] <- sum_mu[, p_rows[j_loc]] + as.numeric(mu_block[, cols, drop = FALSE] %*% w_chunk)
        }
      } else {
        for (j_loc in seq_len(np)) {
          cols <- ((j_loc - 1L) * nr + 1L):(j_loc * nr)
          out_idx <- (p_rows[j_loc] - 1L) * n_ref + r_rows
          mu_mean[out_idx] <- unname(colMeans(mu_block[, cols, drop = FALSE]))
          mu_sd[out_idx] <- unname(apply(mu_block[, cols, drop = FALSE], 2L, stats::sd))
        }
      }
    }
  }

  if (summary) {
    return(list(sum_mu = sum_mu, count_p = rep(1, P), labels = params$labels[[pf]], weights = weights))
  }
  list(mu_fair = mu_mean, mu_fair_sd = mu_sd, weights = weights)
}

#' Resolve and normalize reference-grid weights for standardized fair scores
#' @keywords internal
resolve_reference_weights <- function(reference_grid, reference_weights = NULL,
                                      facet_labels = list()) {
  n_ref <- nrow(reference_grid)
  if (is.null(reference_weights)) {
    if ("weight" %in% names(reference_grid) && !("weight" %in% names(facet_labels))) {
      w <- reference_grid$weight
      if (!is.numeric(w) || any(!is.finite(w)) || any(w < 0) || sum(w) <= 0) {
        stop("reference_grid$weight must contain finite non-negative numbers with a positive sum")
      }
      return(as.numeric(w) / sum(w))
    }
    return(rep(1 / n_ref, n_ref))
  }

  if (is.numeric(reference_weights)) {
    if (length(reference_weights) != n_ref) {
      stop("Numeric reference_weights must have length nrow(reference_grid) (", n_ref, "), got ", length(reference_weights))
    }
    if (any(!is.finite(reference_weights)) || any(reference_weights < 0) || sum(reference_weights) <= 0) {
      stop("reference_weights must contain finite non-negative numbers with a positive sum")
    }
    return(as.numeric(reference_weights) / sum(reference_weights))
  }

  if (is.list(reference_weights)) {
    w_names <- names(reference_weights)
    if (is.null(w_names) || !length(w_names) || any(!nzchar(w_names)) || anyDuplicated(w_names)) {
      stop("List reference_weights must be a non-empty named list keyed by facet names in reference_grid")
    }
    missing_f <- setdiff(w_names, names(reference_grid))
    if (length(missing_f)) {
      stop("reference_weights contains facet(s) absent from reference_grid: ", paste(missing_f, collapse = ", "))
    }
    w <- rep(1, n_ref)
    for (f in w_names) {
      wf <- reference_weights[[f]]
      if (!is.numeric(wf) || is.null(names(wf)) || any(!nzchar(names(wf))) || anyDuplicated(names(wf))) {
        stop("reference_weights[['", f, "']] must be a uniquely named numeric vector")
      }
      if (any(!is.finite(wf)) || any(wf < 0) || sum(wf) <= 0) {
        stop("reference_weights[['", f, "']] must contain finite non-negative numbers with a positive sum")
      }
      if (!is.null(facet_labels[[f]])) {
        bad_lv <- setdiff(names(wf), facet_labels[[f]])
        if (length(bad_lv)) {
          stop("reference_weights[['", f, "']] contains unknown level(s): ", paste(bad_lv, collapse = ", "))
        }
      }
      grid_vals <- as.character(reference_grid[[f]])
      missing_lv <- setdiff(unique(grid_vals), names(wf))
      if (length(missing_lv)) {
        stop("reference_weights[['", f, "']] is missing weight(s) for reference_grid level(s): ",
             paste(missing_lv, collapse = ", "))
      }
      w <- w * unname(wf[grid_vals])
    }
    if (sum(w) <= 0) {
      stop("Combined reference_weights have zero total weight across reference_grid")
    }
    return(as.numeric(w) / sum(w))
  }

  stop("reference_weights must be NULL, a numeric vector of length nrow(reference_grid), or a named list of facet weights")
}

#' Extract a scalar Stan parameter vector of length D
#' @keywords internal
extract_scalar_draws <- function(fit, var_name, draw_idx = NULL) {
  mat <- unclass(fit$draws(var_name, format = "matrix"))
  if (ncol(mat) != 1L) {
    stop("Expected scalar parameter '", var_name, "', found ", ncol(mat), " columns")
  }
  v <- as.numeric(mat[, 1L])
  if (!is.null(draw_idx)) v <- v[draw_idx]
  v
}

#' Resolve draw indices from ndraws / draw_ids
#' @keywords internal
resolve_draw_idx <- function(D_total, ndraws = NULL, draw_ids = NULL, seed = NULL) {
  if (!is.null(ndraws) && !is.null(draw_ids)) {
    stop("Specify either 'ndraws' or 'draw_ids', not both")
  }
  if (!is.null(draw_ids)) {
    if (!is.numeric(draw_ids) || !length(draw_ids) || any(!is.finite(draw_ids)) ||
        any(draw_ids != floor(draw_ids)) || any(draw_ids < 1L) || any(draw_ids > D_total)) {
      stop("draw_ids must be integers between 1 and ", D_total)
    }
    return(as.integer(draw_ids))
  }
  if (!is.null(ndraws)) {
    if (!is.numeric(ndraws) || length(ndraws) != 1L || !is.finite(ndraws) || ndraws < 1) {
      stop("ndraws must be a positive integer")
    }
    nd <- as.integer(min(D_total, floor(ndraws)))
    if (nd < D_total) {
      draw_fn <- function() sort(sample.int(D_total, size = nd, replace = FALSE))
      return(if (!is.null(seed)) with_bmfrm_seed(seed, draw_fn()) else draw_fn())
    }
  }
  NULL
}

#' Extract draw identifiers (draw_id, chain, iteration) for posterior draws
#' @keywords internal
extract_draw_index_info <- function(fit, D_total, draw_idx = NULL) {
  ids <- if (is.null(draw_idx)) seq_len(D_total) else as.integer(draw_idx)
  n_chains <- tryCatch({
    if (is.function(fit$num_chains)) fit$num_chains() else NULL
  }, error = function(e) NULL)
  if (is.null(n_chains)) {
    md <- tryCatch(fit$metadata(), error = function(e) NULL)
    if (!is.null(md) && is.numeric(md$iter_sampling) && length(md$iter_sampling) == 1L && md$iter_sampling >= 1) {
      thin <- if (is.numeric(md$thin) && length(md$thin) == 1L && md$thin >= 1) as.integer(md$thin) else 1L
      per_chain <- as.integer(floor(md$iter_sampling / thin))
      if (per_chain >= 1L && (D_total %% per_chain == 0L)) {
        n_chains <- as.integer(D_total / per_chain)
      }
    }
  }
  if (is.numeric(n_chains) && length(n_chains) == 1L && is.finite(n_chains) &&
      n_chains >= 1 && (D_total %% as.integer(n_chains) == 0L)) {
    C <- as.integer(n_chains)
    S <- as.integer(D_total / C)
    chain_vec <- as.integer(((ids - 1L) %/% S) + 1L)
    iter_vec <- as.integer(((ids - 1L) %% S) + 1L)
  } else {
    chain_vec <- rep(NA_integer_, length(ids))
    iter_vec <- rep(NA_integer_, length(ids))
  }
  list(
    draw_ids = ids,
    draw_names = paste0("draw_", ids),
    draw_info = tibble::tibble(
      draw_id = ids,
      chain = chain_vec,
      iteration = iter_vec
    )
  )
}

#' Draw a single sum-to-zero constrained 1D interaction slice along a known facet axis
#'
#' Generates a `D x n_known` interaction slice for one unseen level of the partner
#' facet, enforcing exact sum-to-zero along the `n_known` observed levels of the
#' known facet, exact 0 when `n_known <= 1L`, and matching the fitted model's
#' prior mode (`gaussian_basis`, `shrinkage_student_t`, or `shrinkage_horseshoe`)
#' and scale convention (`raw` vs `marginal`) without altering existing cells.
#'
#' @keywords internal
draw_constrained_interaction_slice <- function(D, n_known, b_info, basis, seed) {
  if (n_known <= 1L) {
    return(matrix(0, nrow = D, ncol = n_known))
  }
  with_bmfrm_seed(seed, {
    b_mode <- b_info$mode
    sd_b <- b_info$sd
    if (identical(b_mode, "shrinkage_student_t")) {
      lam <- matrix(1 / stats::rgamma(D * n_known, shape = 1.5, rate = 0.5),
                    nrow = D, ncol = n_known)
      z <- matrix(stats::rnorm(D * n_known), nrow = D, ncol = n_known)
      s_axis <- if (identical(b_info$scale, "marginal")) sqrt(n_known / (n_known - 1L)) else 1
      u <- (sqrt(lam) * z) * (s_axis * sd_b)
      u - rowMeans(u)
    } else if (identical(b_mode, "shrinkage_horseshoe")) {
      lam_raw <- matrix(abs(stats::rcauchy(D * n_known, location = 0, scale = 1)),
                        nrow = D, ncol = n_known)
      sig2 <- sd_b * sd_b
      c2 <- b_info$c2
      lam_tilde <- sqrt((c2 * (lam_raw * lam_raw)) / (c2 + sig2 * (lam_raw * lam_raw)))
      z <- matrix(stats::rnorm(D * n_known), nrow = D, ncol = n_known)
      u <- (lam_tilde * z) * sd_b
      u - rowMeans(u)
    } else {
      z <- matrix(stats::rnorm(D * (n_known - 1L)), nrow = D, ncol = n_known - 1L)
      (z %*% t(basis)) * sd_b
    }
  })
}

#' Draw an unconstrained population interaction cell when both facets are unseen
#' @keywords internal
draw_unconstrained_interaction_cell <- function(D, b_info, seed) {
  with_bmfrm_seed(seed, {
    b_mode <- b_info$mode
    sd_b <- b_info$sd
    if (identical(b_mode, "shrinkage_student_t")) {
      lam <- 1 / stats::rgamma(D, shape = 1.5, rate = 0.5)
      sqrt(lam) * stats::rnorm(D) * sd_b
    } else if (identical(b_mode, "shrinkage_horseshoe")) {
      lam_raw <- abs(stats::rcauchy(D, location = 0, scale = 1))
      sig2 <- sd_b * sd_b
      c2 <- b_info$c2
      lam_tilde <- sqrt((c2 * (lam_raw * lam_raw)) / (c2 + sig2 * (lam_raw * lam_raw)))
      lam_tilde * stats::rnorm(D) * sd_b
    } else {
      stats::rnorm(D) * sd_b
    }
  })
}

#' Prepare lightweight evaluation plan with chunk-scoped new-level generation
#'
#' Validates all columns, unseen levels, and hierarchical requirements up front
#' WITHOUT allocating `D x n_new` matrices for unseen levels. Instead stores
#' compact length-D hyperprior vectors and one 32-bit seed per unique unseen
#' level so each chunk generates only the `<= chunk_size` unseen levels it uses
#' while preserving exact cross-chunk consistency and sum-to-zero interaction
#' constraints.
#'
#' @keywords internal
prepare_eval_plan <- function(object, newdata,
                              adjust_facets = character(),
                              draw_idx = NULL,
                              allow_new_levels = FALSE,
                              sample_new_levels = c("hierarchical", "zero")) {
  sample_new_levels <- match.arg(sample_new_levels)
  spec <- object$spec
  params <- extract_bmfrm_param_draws(object, adjust_facets = adjust_facets, draw_idx = draw_idx)
  D <- params$D
  K <- params$K
  labels <- params$labels
  hier <- spec$hierarchical %||% character()

  ids <- list()
  new_levels <- list()
  main_new <- list()

  # 1. Validate main facets and capture compact hyperprior vectors + per-level seeds
  for (f in spec$facets_main) {
    if (f %in% adjust_facets) next
    if (!f %in% names(newdata)) {
      stop("Column '", f, "' not found in newdata")
    }
    vals <- newdata[[f]]
    if (anyNA(vals)) {
      stop("Column '", f, "' in newdata contains NA values")
    }
    vals_chr <- as.character(vals)
    known <- labels[[f]]
    idx <- match(vals_chr, known)
    unseen <- unique(vals_chr[is.na(idx)])

    if (length(unseen)) {
      if (identical(params$tau_info$mode, "pcm_ragged") &&
          identical(f, params$tau_info$threshold_facet)) {
        stop(
          "Unseen level(s) of threshold_facet '", f, "' (",
          paste(utils::head(unseen, 5L), collapse = ", "),
          ") are not supported in v1 heterogeneous-category PCM (even with allow_new_levels = TRUE ",
          "or sample_new_levels = 'zero'), because their rubric category count and step distribution are undefined."
        )
      }
      if (!isTRUE(allow_new_levels)) {
        preview <- paste(utils::head(unseen, 5L), collapse = ", ")
        if (length(unseen) > 5L) preview <- paste0(preview, ", ...")
        stop(
          "Unknown level(s) in facet '", f, "': ", preview,
          ". Set allow_new_levels = TRUE (with sample_new_levels = 'hierarchical' or 'zero') ",
          "or adjust '", f, "' explicitly via adjust_facets."
        )
      }
      new_levels[[f]] <- unseen
      n_new <- length(unseen)
      idx[is.na(idx)] <- length(known) + match(vals_chr[is.na(idx)], unseen)

      if (identical(sample_new_levels, "zero")) {
        main_new[[f]] <- list(mode = "zero", n_known = length(known), seeds = integer(0))
      } else if (identical(f, spec_person_facet(spec))) {
        pf_cur <- spec_person_facet(spec)
        if (!(pf_cur %in% hier) ||
            !fit_has_variable(object$fit, "mu_theta") ||
            !fit_has_variable(object$fit, "sigma_theta")) {
          stop(
            "Cannot sample new levels for non-hierarchical facet '", pf_cur, "' under ",
            "sample_new_levels = 'hierarchical'. Refit with hierarchical = c(..., '", pf_cur, "') ",
            "or explicitly set sample_new_levels = 'zero'."
          )
        }
        main_new[[f]] <- list(
          mode = "hierarchical_person",
          n_known = length(known),
          mu = extract_scalar_draws(object$fit, "mu_theta", draw_idx),
          sd = extract_scalar_draws(object$fit, "sigma_theta", draw_idx),
          seeds = sample.int(.Machine$integer.max, n_new, replace = TRUE)
        )
      } else {
        sig_var <- paste0("sigma_", f)
        if (!(f %in% hier) || !fit_has_variable(object$fit, sig_var)) {
          stop(
            "Cannot sample new levels for non-hierarchical facet '", f, "' under ",
            "sample_new_levels = 'hierarchical'. Refit with hierarchical = c(..., '", f, "'), ",
            "set adjust_facets = '", f, "' to fix its effect at 0, or explicitly pass ",
            "sample_new_levels = 'zero'."
          )
        }
        main_new[[f]] <- list(
          mode = "hierarchical_facet",
          n_known = length(known),
          sd = extract_scalar_draws(object$fit, sig_var, draw_idx),
          seeds = sample.int(.Machine$integer.max, n_new, replace = TRUE)
        )
      }
    }
    ids[[f]] <- idx
  }

  # 2. PCM step thresholds for unseen threshold_facet levels (compact hyperprior + seeds only)
  tau_new <- NULL
  if (identical(params$tau_info$mode, "pcm_ragged")) {
    tf <- params$tau_info$threshold_facet
    if (length(new_levels[[tf]]) > 0L) {
      stop(
        "Unseen level(s) of threshold_facet '", tf, "' (",
        paste(utils::head(new_levels[[tf]], 5L), collapse = ", "),
        ") are not supported in v1 heterogeneous-category PCM (even with allow_new_levels = TRUE ",
        "or sample_new_levels = 'zero'), because their rubric category count and step distribution are undefined."
      )
    }
  } else if (identical(params$tau_info$mode, "pcm")) {
    tf <- params$tau_info$threshold_facet
    n_new_tf <- length(new_levels[[tf]])
    if (n_new_tf > 0L) {
      n_known_tf <- length(labels[[tf]])
      if (identical(sample_new_levels, "zero")) {
        if (fit_has_variable(object$fit, "tau_common")) {
          tc <- extract_ordered_draws(object$fit, "tau_common", K - 1L)
          if (!is.null(draw_idx)) tc <- tc[draw_idx, , drop = FALSE]
        } else {
          tc <- matrix(0, nrow = D, ncol = K - 1L)
        }
        cum_zero <- matrix(0, nrow = D, ncol = K)
        if (K > 1L) {
          for (k in 2:K) cum_zero[, k] <- cum_zero[, k - 1L] + tc[, k - 1L]
        }
        tau_new <- list(mode = "zero", n_known = n_known_tf, cum_zero = cum_zero)
      } else {
        if (!("tau" %in% hier) ||
            !fit_has_variable(object$fit, "tau_common") ||
            !fit_has_variable(object$fit, "sigma_tau")) {
          stop(
            "Cannot sample new PCM step thresholds for unseen level(s) of '", tf,
            "' because 'tau' was not fit hierarchically. Refit with hierarchical = c(..., 'tau') ",
            "or explicitly set sample_new_levels = 'zero'."
          )
        }
        tc <- extract_ordered_draws(object$fit, "tau_common", K - 1L)
        if (!is.null(draw_idx)) tc <- tc[draw_idx, , drop = FALSE]
        sd_tau <- extract_scalar_draws(object$fit, "sigma_tau", draw_idx)
        dev_scale <- spec$gaussian_model$step$dev_scale %||% "raw"
        B_dev <- if (K > 2L) rsm_prior_basis(K - 1L, dev_scale) else NULL
        tau_new <- list(
          mode = "hierarchical",
          n_known = n_known_tf,
          tau_common = tc,
          sd = sd_tau,
          B_dev = B_dev,
          seeds = sample.int(.Machine$integer.max, n_new_tf, replace = TRUE)
        )
      }
    }
  }

  # 3. Constrained two-way interactions involving unseen levels (compact hyperprior + seeds only)
  bias_new <- list()
  for (term in names(params$bias_draws)) {
    b <- params$bias_draws[[term]]
    f1 <- b$f1
    f2 <- b$f2
    n1_known <- b$n1
    n2_known <- length(labels[[f2]])
    i1 <- ids[[f1]]
    i2 <- ids[[f2]]
    is_new_pair <- (i1 > n1_known) | (i2 > n2_known)
    if (!any(is_new_pair)) next

    n1_new <- length(new_levels[[f1]])
    n2_new <- length(new_levels[[f2]])

    if (identical(sample_new_levels, "zero")) {
      bias_new[[term]] <- list(
        mode = "zero", f1 = f1, f2 = f2,
        n1_known = n1_known, n2_known = n2_known
      )
      next
    }

    b_name <- gsub(":", "", term)
    sig_var <- paste0("sigma_", b_name)
    if (!(term %in% hier) || !fit_has_variable(object$fit, sig_var)) {
      stop(
        "Cannot sample new interaction cells for non-hierarchical term '", term,
        "' under sample_new_levels = 'hierarchical'. Refit with hierarchical = c(..., '", term,
        "'), adjust one of its facets via adjust_facets, or explicitly pass sample_new_levels = 'zero'."
      )
    }
    sd_b <- extract_scalar_draws(object$fit, sig_var, draw_idx)
    b_spec <- spec$gaussian_model$biases[[term]]
    b_mode <- if (!is.null(b_spec)) b_spec$mode else "gaussian_basis"
    b_scale <- if (!is.null(b_spec)) (b_spec$scale %||% "raw") else "raw"
    c2_b <- NULL
    if (identical(b_mode, "shrinkage_horseshoe")) {
      c2_var <- paste0("c2_", b_name)
      if (!fit_has_variable(object$fit, c2_var)) {
        stop("Cannot sample regularized horseshoe interaction '", term, "': missing parameter '", c2_var, "'")
      }
      c2_b <- extract_scalar_draws(object$fit, c2_var, draw_idx)
    }

    both_new_rows <- which((i1 > n1_known) & (i2 > n2_known))
    both_new_map <- integer(nrow(newdata))
    both_new_seeds <- integer(0)
    if (length(both_new_rows)) {
      pair_keys <- paste(i1[both_new_rows], i2[both_new_rows], sep = "\r")
      uniq_keys <- unique(pair_keys)
      both_new_map[both_new_rows] <- match(pair_keys, uniq_keys)
      both_new_seeds <- sample.int(.Machine$integer.max, length(uniq_keys), replace = TRUE)
    }

    col_scope <- vector("list", n2_new)
    col_last_row <- integer(n2_new)
    if (n2_new > 0L && n1_known > 1L) {
      kn_rows <- which(i1 <= n1_known & i2 > n2_known)
      if (length(kn_rows)) {
        u2_all <- i2[kn_rows] - n2_known
        sp_i1 <- split.default(i1[kn_rows], u2_all)
        for (nm_u in names(sp_i1)) {
          col_scope[[as.integer(nm_u)]] <- sort(unique(sp_i1[[nm_u]]))
        }
        for (k_i in seq_along(kn_rows)) {
          u_i <- u2_all[k_i]
          r_i <- kn_rows[k_i]
          if (r_i > col_last_row[u_i]) col_last_row[u_i] <- r_i
        }
      }
    }

    row_scope <- vector("list", n1_new)
    row_last_row <- integer(n1_new)
    if (n1_new > 0L && n2_known > 1L) {
      nk_rows <- which(i1 > n1_known & i2 <= n2_known)
      if (length(nk_rows)) {
        u1_all <- i1[nk_rows] - n1_known
        sp_i2 <- split.default(i2[nk_rows], u1_all)
        for (nm_u in names(sp_i2)) {
          row_scope[[as.integer(nm_u)]] <- sort(unique(sp_i2[[nm_u]]))
        }
        for (k_i in seq_along(nk_rows)) {
          u_i <- u1_all[k_i]
          r_i <- nk_rows[k_i]
          if (r_i > row_last_row[u_i]) row_last_row[u_i] <- r_i
        }
      }
    }

    bias_new[[term]] <- list(
      mode = b_mode,
      scale = b_scale,
      f1 = f1, f2 = f2,
      n1_known = n1_known, n2_known = n2_known,
      sd = sd_b,
      c2 = c2_b,
      col_seeds = if (n2_new > 0L && n1_known > 1L) sample.int(.Machine$integer.max, n2_new, replace = TRUE) else integer(0),
      col_scope = col_scope,
      col_last_row = col_last_row,
      row_seeds = if (n1_new > 0L && n2_known > 1L) sample.int(.Machine$integer.max, n1_new, replace = TRUE) else integer(0),
      row_scope = row_scope,
      row_last_row = row_last_row,
      both_new_map = both_new_map,
      both_new_seeds = both_new_seeds
    )
  }

  n_conc_init <- if (length(bias_new) > 0L && !identical(sample_new_levels, "zero")) {
    count_concurrent_interaction_slices(ids, bias_new, seq_len(nrow(newdata)), 1L)
  } else 0L

  list(
    params = params,
    ids = ids,
    new_levels = new_levels,
    main_new = main_new,
    tau_new = tau_new,
    bias_new = bias_new,
    sample_new_levels = sample_new_levels,
    offset_cache = if (length(bias_new) > 0L && !identical(sample_new_levels, "zero")) {
      new_interaction_offset_cache(D, 1L, n_concurrent_offsets = n_conc_init)
    } else NULL
  )
}

#' Mix a per-entity 32-bit seed with an optional Monte Carlo iteration salt
#' @keywords internal
mix_bmfrm_seed <- function(base_seed, salt = NULL) {
  if (is.null(salt)) return(base_seed)
  as.integer(((as.numeric(base_seed) * 1048583 + as.numeric(salt)) %% 2147483647) + 1)
}

#' Count peak simultaneously live constrained interaction slices across chunks
#' @keywords internal
count_concurrent_interaction_slices <- function(ids, bias_new, rows_seq, cs) {
  n_seq <- length(rows_seq)
  cs <- as.integer(max(1L, cs))
  n_chunks <- as.integer((n_seq + cs - 1L) %/% cs)
  if (n_chunks <= 1L || length(bias_new) == 0L) return(0L)

  diff_active <- integer(n_chunks + 2L)
  add_axis_spans <- function(u_vec, pos_in_seq, n_u) {
    if (n_u <= 0L || length(pos_in_seq) <= 1L) return(invisible(NULL))
    ch_vec <- (pos_in_seq - 1L) %/% cs + 1L
    first_ch <- integer(n_u)
    last_ch <- integer(n_u)
    for (idx in seq_along(u_vec)) {
      u <- u_vec[idx]
      ch <- ch_vec[idx]
      if (first_ch[u] == 0L) first_ch[u] <- ch
      if (ch > last_ch[u]) last_ch[u] <- ch
    }
    multi <- which(first_ch > 0L & last_ch > first_ch)
    for (u in multi) {
      diff_active[first_ch[u]] <<- diff_active[first_ch[u]] + 1L
      diff_active[last_ch[u] + 1L] <<- diff_active[last_ch[u] + 1L] - 1L
    }
    invisible(NULL)
  }

  for (term in names(bias_new)) {
    bn <- bias_new[[term]]
    if (identical(bn$mode, "zero")) next
    i1 <- ids[[bn$f1]][rows_seq]
    i2 <- ids[[bn$f2]][rows_seq]
    if (bn$n1_known > 1L && length(bn$col_seeds) > 0L) {
      kn <- which(i1 <= bn$n1_known & i2 > bn$n2_known)
      if (length(kn) > 1L) {
        add_axis_spans(i2[kn] - bn$n2_known, kn, length(bn$col_seeds))
      }
    }
    if (bn$n2_known > 1L && length(bn$row_seeds) > 0L) {
      nk <- which(i1 > bn$n1_known & i2 <= bn$n2_known)
      if (length(nk) > 1L) {
        add_axis_spans(i1[nk] - bn$n1_known, nk, length(bn$row_seeds))
      }
    }
  }

  as.integer(max(0L, cumsum(diff_active)))
}

#' Create a byte-bounded cache for length-D interaction centering offsets
#'
#' Sizes the default byte budget from the estimated simultaneously live offsets
#' (`n_concurrent_offsets * D * 8`, e.g., `M * n_active_slices * D * 8` in CV)
#' and current chunk footprint, capped at `hard_max_bytes` (default 64 MiB).
#' When demand exceeds `max_bytes`, cached entries remain pinned until their
#' final chunk and remaining offsets are deterministically recomputed per chunk.
#' @keywords internal
new_interaction_offset_cache <- function(D, eff_chunk_size = 1L,
                                         n_concurrent_offsets = 0L,
                                         max_bytes = NULL,
                                         hard_max_bytes = 64 * 1024^2) {
  if (is.null(max_bytes)) {
    entry_bytes <- as.numeric(D) * 8
    est_bytes <- as.numeric(max(0, n_concurrent_offsets)) * entry_bytes
    chunk_bytes <- 16 * as.numeric(D) * as.numeric(max(1L, eff_chunk_size)) * 8
    target <- max(1024^2, chunk_bytes, est_bytes)
    max_bytes <- min(as.numeric(hard_max_bytes), target)
  }
  env <- new.env(hash = TRUE, parent = emptyenv())
  env$data <- new.env(hash = TRUE, parent = emptyenv())
  env$bytes <- 0
  env$max_bytes <- as.numeric(max_bytes)
  env
}

#' Generate one raw length-D interaction column inside an active RNG guard
#' @keywords internal
draw_raw_interaction_col_unseeded <- function(D, n_known, b_info, col_seed) {
  set.seed(col_seed)
  b_mode <- b_info$mode
  sd_b <- b_info$sd
  s_axis <- if (identical(b_info$scale, "marginal")) sqrt(n_known / (n_known - 1L)) else 1
  if (identical(b_mode, "shrinkage_student_t")) {
    lam <- 1 / stats::rgamma(D, shape = 1.5, rate = 0.5)
    (sqrt(lam) * stats::rnorm(D)) * (s_axis * sd_b)
  } else if (identical(b_mode, "shrinkage_horseshoe")) {
    lam_raw <- abs(stats::rcauchy(D, location = 0, scale = 1))
    sig2 <- sd_b * sd_b
    c2 <- b_info$c2
    lam_tilde <- sqrt((c2 * (lam_raw * lam_raw)) / (c2 + sig2 * (lam_raw * lam_raw)))
    (lam_tilde * stats::rnorm(D)) * sd_b
  } else {
    stats::rnorm(D) * (s_axis * sd_b)
  }
}

#' Compute length-D sum-to-zero centering offset and requested chunk columns
#' @keywords internal
compute_interaction_center_unseeded <- function(D, n_known, b_info, slice_seed,
                                                scope_cols, u_c = integer(0)) {
  is_gauss <- !identical(b_info$mode, "shrinkage_student_t") &&
    !identical(b_info$mode, "shrinkage_horseshoe")
  iter_cols <- if (is_gauss) scope_cols else seq_len(n_known)
  k_iter <- length(iter_cols)
  sum_raw <- numeric(D)
  k_c <- length(u_c)
  raw_c <- if (k_c > 0L) matrix(0, nrow = D, ncol = k_c) else NULL
  c_map <- if (k_c > 0L) match(iter_cols, u_c, nomatch = 0L) else integer(k_iter)

  for (idx in seq_len(k_iter)) {
    ell <- iter_cols[idx]
    v <- draw_raw_interaction_col_unseeded(
      D, n_known, b_info, mix_bmfrm_seed(slice_seed, ell)
    )
    sum_raw <- sum_raw + v
    pos_c <- c_map[idx]
    if (pos_c > 0L) raw_c[, pos_c] <- v
  }

  alpha <- if (is_gauss && k_iter < n_known) {
    1 - sqrt((n_known - k_iter) / n_known)
  } else {
    1
  }
  list(center = (alpha / k_iter) * sum_raw, raw_c = raw_c)
}

#' Evaluate constrained 1D interaction values for only the columns used in a chunk
#' @keywords internal
eval_constrained_interaction_chunk <- function(D, n_known, b_info, slice_seed,
                                               chunk_k_idx, scope_cols,
                                               need_after_chunk = FALSE,
                                               cache_key = NULL,
                                               offset_cache = NULL) {
  nc <- length(chunk_k_idx)
  if (n_known <= 1L || nc == 0L) {
    return(matrix(0, nrow = D, ncol = nc))
  }
  u_c <- unique(chunk_k_idx)
  k_c <- length(u_c)
  entry_bytes <- as.numeric(D) * 8

  has_hit <- !is.null(offset_cache) && !is.null(cache_key) &&
    exists(cache_key, envir = offset_cache$data, inherits = FALSE)

  if (has_hit) {
    center_d <- offset_cache$data[[cache_key]]
    if (!isTRUE(need_after_chunk)) {
      rm(list = cache_key, envir = offset_cache$data)
      offset_cache$bytes <- max(0, offset_cache$bytes - entry_bytes)
    }
    raw_c <- with_bmfrm_seed(1L, {
      out_c <- matrix(0, nrow = D, ncol = k_c)
      for (j in seq_len(k_c)) {
        out_c[, j] <- draw_raw_interaction_col_unseeded(
          D, n_known, b_info, mix_bmfrm_seed(slice_seed, u_c[j])
        )
      }
      out_c
    })
    sub_mat <- raw_c - center_d
    return(sub_mat[, match(chunk_k_idx, u_c), drop = FALSE])
  }

  res <- with_bmfrm_seed(1L, {
    compute_interaction_center_unseeded(
      D = D, n_known = n_known, b_info = b_info,
      slice_seed = slice_seed, scope_cols = scope_cols, u_c = u_c
    )
  })
  if (isTRUE(need_after_chunk) && !is.null(offset_cache) && !is.null(cache_key)) {
    if (offset_cache$bytes + entry_bytes <= offset_cache$max_bytes) {
      offset_cache$data[[cache_key]] <- res$center
      offset_cache$bytes <- offset_cache$bytes + entry_bytes
    }
  }
  sub_mat <- res$raw_c - res$center
  sub_mat[, match(chunk_k_idx, u_c), drop = FALSE]
}

#' Evaluate a single observation chunk using chunk-scoped new-level generation
#' @keywords internal
eval_bmfrm_chunk_from_plan <- function(plan, rows, resp_chunk = NULL,
                                       compute_mu = TRUE,
                                       compute_linpred = FALSE,
                                       compute_probs = FALSE,
                                       compute_kurtosis = FALSE,
                                       compute_log_lik = FALSE,
                                       compute_y_rep = FALSE,
                                       seed_salt = NULL,
                                       offset_cache = plan$offset_cache,
                                       bias_scope_override = NULL) {
  params <- plan$params
  ids <- plan$ids
  D <- params$D
  nc <- length(rows)

  # 1. Person / ability linear predictor contribution
  pf <- plan$params$person_facet %||% "person"
  p_id <- ids[[pf]][rows]
  p_new_info <- plan$main_new[[pf]]
  if (is.null(p_new_info)) {
    eta <- params$theta[, p_id, drop = FALSE]
  } else {
    eta <- matrix(0, nrow = D, ncol = nc)
    n_k <- p_new_info$n_known
    k_pos <- which(p_id <= n_k)
    if (length(k_pos)) eta[, k_pos] <- params$theta[, p_id[k_pos], drop = FALSE]
    u_pos <- which(p_id > n_k)
    if (length(u_pos) && !identical(p_new_info$mode, "zero")) {
      u_local <- p_id[u_pos] - n_k
      for (m in unique(u_local)) {
        s_m <- mix_bmfrm_seed(p_new_info$seeds[m], seed_salt)
        v <- with_bmfrm_seed(s_m, p_new_info$mu + p_new_info$sd * stats::rnorm(D))
        eta[, u_pos[u_local == m]] <- v
      }
    }
  }

  # 2. Non-person main facet contributions
  for (f in names(params$facet_draws)) {
    f_id <- ids[[f]][rows]
    f_new_info <- plan$main_new[[f]]
    if (is.null(f_new_info)) {
      eta <- eta - params$facet_draws[[f]][, f_id, drop = FALSE]
    } else {
      n_k <- f_new_info$n_known
      k_pos <- which(f_id <= n_k)
      if (length(k_pos)) {
        eta[, k_pos] <- eta[, k_pos, drop = FALSE] - params$facet_draws[[f]][, f_id[k_pos], drop = FALSE]
      }
      u_pos <- which(f_id > n_k)
      if (length(u_pos) && !identical(f_new_info$mode, "zero")) {
        u_local <- f_id[u_pos] - n_k
        for (m in unique(u_local)) {
          s_m <- mix_bmfrm_seed(f_new_info$seeds[m], seed_salt)
          v <- with_bmfrm_seed(s_m, f_new_info$sd * stats::rnorm(D))
          pos_m <- u_pos[u_local == m]
          eta[, pos_m] <- eta[, pos_m, drop = FALSE] - v
        }
      }
    }
  }

  # 3. Two-way interaction contributions (constrained along known axes!)
  for (term in names(params$bias_draws)) {
    b <- params$bias_draws[[term]]
    i1 <- ids[[b$f1]][rows]
    i2 <- ids[[b$f2]][rows]
    b_info <- plan$bias_new[[term]]
    if (is.null(b_info)) {
      cell_idx <- i1 + (i2 - 1L) * b$n1
      eta <- eta + b$draws[, cell_idx, drop = FALSE]
    } else {
      n1 <- b_info$n1_known
      n2 <- b_info$n2_known
      kk <- which(i1 <= n1 & i2 <= n2)
      if (length(kk)) {
        cell_idx <- i1[kk] + (i2[kk] - 1L) * n1
        eta[, kk] <- eta[, kk, drop = FALSE] + b$draws[, cell_idx, drop = FALSE]
      }
      if (!identical(b_info$mode, "zero")) {
        max_row_c <- max(rows)
        ov <- if (!is.null(bias_scope_override)) bias_scope_override[[term]] else NULL
        c_scope_list <- if (!is.null(ov)) ov$col_scope else b_info$col_scope
        c_last_vec <- if (!is.null(ov)) ov$col_last_row else b_info$col_last_row
        r_scope_list <- if (!is.null(ov)) ov$row_scope else b_info$row_scope
        r_last_vec <- if (!is.null(ov)) ov$row_last_row else b_info$row_last_row

        # Known f1 x New f2: sum-to-zero along known f1 (identically 0 when n1 <= 1L)
        kn <- which(i1 <= n1 & i2 > n2)
        if (n1 > 1L && length(kn)) {
          u2 <- i2[kn] - n2
          for (m2 in unique(u2)) {
            s_c <- mix_bmfrm_seed(b_info$col_seeds[m2], seed_salt)
            pos <- kn[u2 == m2]
            sc_m2 <- c_scope_list[[m2]]
            if (is.null(sc_m2)) sc_m2 <- sort(unique(i1[pos]))
            eta[, pos] <- eta[, pos, drop = FALSE] + eval_constrained_interaction_chunk(
              D = D, n_known = n1, b_info = b_info, slice_seed = s_c,
              chunk_k_idx = i1[pos], scope_cols = sc_m2,
              need_after_chunk = (c_last_vec[m2] > max_row_c),
              cache_key = paste0(term, ":c:", s_c),
              offset_cache = offset_cache
            )
          }
        }
        # New f1 x Known f2: sum-to-zero along known f2 (identically 0 when n2 <= 1L)
        nk <- which(i1 > n1 & i2 <= n2)
        if (n2 > 1L && length(nk)) {
          u1 <- i1[nk] - n1
          for (m1 in unique(u1)) {
            s_r <- mix_bmfrm_seed(b_info$row_seeds[m1], seed_salt)
            pos <- nk[u1 == m1]
            sc_m1 <- r_scope_list[[m1]]
            if (is.null(sc_m1)) sc_m1 <- sort(unique(i2[pos]))
            eta[, pos] <- eta[, pos, drop = FALSE] + eval_constrained_interaction_chunk(
              D = D, n_known = n2, b_info = b_info, slice_seed = s_r,
              chunk_k_idx = i2[pos], scope_cols = sc_m1,
              need_after_chunk = (r_last_vec[m1] > max_row_c),
              cache_key = paste0(term, ":r:", s_r),
              offset_cache = offset_cache
            )
          }
        }
        # New f1 x New f2: unconstrained population cell
        nn <- which(i1 > n1 & i2 > n2)
        if (length(nn)) {
          p_map <- b_info$both_new_map[rows[nn]]
          for (pm in unique(p_map)) {
            s_nn <- mix_bmfrm_seed(b_info$both_new_seeds[pm], seed_salt)
            v <- draw_unconstrained_interaction_cell(D, b_info, s_nn)
            pos <- nn[p_map == pm]
            eta[, pos] <- eta[, pos, drop = FALSE] + v
          }
        }
      }
    }
  }

  tau_info <- params$tau_info
  t_idx <- if (tau_info$mode != "shared") ids[[tau_info$threshold_facet]][rows]
  eval_bmfrm_categories(eta, params$K, tau_info, t_idx, resp_chunk,
                        compute_mu, compute_linpred, compute_probs, compute_kurtosis,
                        compute_log_lik, compute_y_rep, tau_new = plan$tau_new, seed_salt = seed_salt)
}

#' Format quantile column names matching facet_summary() (e.g., q2.5, q97.5)
#' @keywords internal
fmt_prob_col <- function(p) {
  paste0("q", sub("\\.?0+$", "", formatC(p * 100, format = "f", digits = 4)))
}

#' Classify the semantic prediction target for metadata clarity
#' @keywords internal
classify_prediction_target <- function(adjust_facets, new_levels, sample_new_levels) {
  has_adj <- length(adjust_facets) > 0L
  has_new <- length(new_levels) > 0L
  if (!has_adj && !has_new) return("conditional_known_levels")
  if (has_adj && !has_new) return("adjusted_reference_zero")
  if (!has_adj && has_new) {
    return(if (identical(sample_new_levels, "hierarchical")) "population_new_levels" else "zero_new_levels")
  }
  if (identical(sample_new_levels, "hierarchical")) "adjusted_and_population_new_levels" else "adjusted_and_zero_new_levels"
}

#' Unified streaming prediction engine powering predict.bmfrm_fit()
#' @keywords internal
stream_bmfrm_predictions <- function(object,
                                     newdata = NULL,
                                     type = c("epred", "predict", "probs", "linpred"),
                                     adjust_facets = character(),
                                     allow_new_levels = FALSE,
                                     sample_new_levels = c("hierarchical", "zero"),
                                     summary = TRUE,
                                     probs = c(0.025, 0.975),
                                     ndraws = NULL,
                                     draw_ids = NULL,
                                     seed = NULL,
                                     chunk_size = NULL,
                                     max_cells = 5e7) {
  type <- match.arg(type)
  sample_new_levels <- match.arg(sample_new_levels)
  if (isTRUE(object$spec$K_heterogeneous) &&
      identical(object$spec$family, "pcm") &&
      !is.null(object$spec$threshold_facet) &&
      (object$spec$threshold_facet %in% adjust_facets)) {
    stop(
      "Adjusting threshold_facet '", object$spec$threshold_facet,
      "' is not supported in v1 heterogeneous-category PCM because no shared reference ",
      "category count or step profile exists across levels."
    )
  }
  if (!can_evaluate_in_r(object, adjust_facets = adjust_facets)) {
    stop("Cannot compute predictions in R for this fit: ", get_eval_contract(object)$reason)
  }

  df <- if (is.null(newdata)) object$data_info$data_clean else newdata
  if (!is.data.frame(df) || !nrow(df)) {
    stop("newdata must be a non-empty data.frame")
  }
  if (isTRUE(summary) && (!is.numeric(probs) || any(!is.finite(probs)) || any(probs < 0 | probs > 1))) {
    stop("probs must be numeric probabilities in [0, 1]")
  }
  if (!is.numeric(max_cells) || length(max_cells) != 1L || !is.finite(max_cells) || max_cells < 1) {
    stop("max_cells must be a positive number")
  }

  D_total <- nrow(unclass(object$fit$draws("theta", format = "matrix")))
  draw_idx <- resolve_draw_idx(D_total, ndraws = ndraws, draw_ids = draw_ids, seed = seed)
  draw_meta <- extract_draw_index_info(object$fit, D_total, draw_idx)

  N <- nrow(df)
  K <- get_spec_K_max(object$spec)
  is_hetero <- isTRUE(object$spec$K_heterogeneous)
  K_by_level_obj <- object$data_info$K_by_level %||% object$spec$K_by_level
  D <- length(draw_meta$draw_ids)

  if (!isTRUE(summary)) {
    req_cells <- as.numeric(D) * as.numeric(N) * if (identical(type, "probs")) as.numeric(K) else 1
    if (req_cells > max_cells) {
      stop(
        "Returning raw posterior ", type, " draws requires ", format(req_cells, big.mark = ","),
        " numeric entries (", D, " draws x ", N, " rows",
        if (identical(type, "probs")) paste0(" x ", K, " categories") else "",
        "), exceeding max_cells = ", format(max_cells, big.mark = ","),
        ". Use summary = TRUE, subset with ndraws/draw_ids, or raise max_cells deliberately."
      )
    }
  }

  run_eval <- function() {
    plan <- prepare_eval_plan(
      object, newdata = df, adjust_facets = adjust_facets, draw_idx = draw_idx,
      allow_new_levels = allow_new_levels, sample_new_levels = sample_new_levels
    )
    eff_cs <- default_bmfrm_chunk_size(
      D * if (identical(type, "probs")) K else 1L, N, chunk_size
    )
    starts <- seq.int(1L, N, by = eff_cs)
    n_conc_p <- if (length(starts) > 1L && length(plan$bias_new) > 0L && !identical(sample_new_levels, "zero")) {
      count_concurrent_interaction_slices(plan$ids, plan$bias_new, seq_len(N), eff_cs)
    } else 0L
    plan$offset_cache <- if (n_conc_p > 0L) {
      new_interaction_offset_cache(D, eff_cs, n_concurrent_offsets = n_conc_p)
    } else NULL

    if (!isTRUE(summary)) {
      if (identical(type, "probs")) {
        out_arr <- array(NA_real_, dim = c(D, N, K),
                         dimnames = list(draw_meta$draw_names, NULL, paste0("cat_", seq_len(K))))
        for (s in starts) {
          e <- min(N, s + eff_cs - 1L)
          rows <- s:e
          ch <- eval_bmfrm_chunk_from_plan(plan, rows, compute_mu = FALSE, compute_probs = TRUE)
          out_arr[, rows, ] <- ch$probs
        }
        return(list(value = out_arr, new_levels = plan$new_levels))
      }

      out_mat <- if (identical(type, "predict")) {
        matrix(NA_integer_, nrow = D, ncol = N, dimnames = list(draw_meta$draw_names, NULL))
      } else {
        matrix(NA_real_, nrow = D, ncol = N, dimnames = list(draw_meta$draw_names, NULL))
      }
      for (s in starts) {
        e <- min(N, s + eff_cs - 1L)
        rows <- s:e
        ch <- eval_bmfrm_chunk_from_plan(
          plan, rows,
          compute_mu = identical(type, "epred"),
          compute_linpred = identical(type, "linpred"),
          compute_y_rep = identical(type, "predict")
        )
        out_mat[, rows] <- switch(
          type,
          epred = ch$mu,
          linpred = ch$eta,
          predict = ch$y_rep
        )
      }
      return(list(value = out_mat, new_levels = plan$new_levels))
    }

    # Streaming summary path (O(D * chunk_size) memory including unseen levels!)
    if (identical(type, "probs")) {
      prob_means <- matrix(NA_real_, nrow = N, ncol = K)
      colnames(prob_means) <- paste0("prob_", seq_len(K))
      for (s in starts) {
        e <- min(N, s + eff_cs - 1L)
        rows <- s:e
        ch <- eval_bmfrm_chunk_from_plan(plan, rows, compute_mu = FALSE, compute_probs = TRUE)
        for (k in seq_len(K)) {
          pk_mat <- matrix(ch$probs[, , k], nrow = D, ncol = length(rows))
          prob_means[rows, k] <- unname(colMeans(pk_mat))
        }
      }
      prob_tbl <- tibble::as_tibble(prob_means)
      if (is_hetero) {
        K_obs_vec <- as.integer(unname(K_by_level_obj[as.character(df[[object$spec$threshold_facet]])]))
        prob_tbl <- tibble::add_column(prob_tbl, K_obs = K_obs_vec, .before = 1L)
      }
      return(list(value = prob_tbl, new_levels = plan$new_levels))
    }

    est <- numeric(N)
    sd_v <- numeric(N)
    q_mat <- if (length(probs)) matrix(NA_real_, nrow = N, ncol = length(probs)) else NULL
    for (s in starts) {
      e <- min(N, s + eff_cs - 1L)
      rows <- s:e
      ch <- eval_bmfrm_chunk_from_plan(
        plan, rows,
        compute_mu = identical(type, "epred"),
        compute_linpred = identical(type, "linpred"),
        compute_y_rep = identical(type, "predict")
      )
      mat <- switch(type, epred = ch$mu, linpred = ch$eta, predict = ch$y_rep)
      est[rows] <- unname(colMeans(mat))
      sd_v[rows] <- unname(apply(mat, 2L, stats::sd))
      if (length(probs)) {
        q_chunk <- apply(mat, 2L, stats::quantile, probs = probs, names = FALSE)
        if (length(probs) == 1L) {
          q_mat[rows, 1L] <- as.numeric(q_chunk)
        } else {
          q_mat[rows, ] <- t(q_chunk)
        }
      }
    }
    res_df <- tibble::tibble(estimate = est, sd = sd_v)
    if (length(probs)) {
      q_names <- vapply(probs, fmt_prob_col, character(1))
      for (j in seq_along(probs)) res_df[[q_names[j]]] <- q_mat[, j]
    }
    list(value = res_df, new_levels = plan$new_levels)
  }

  eval_res <- if (!is.null(seed)) with_bmfrm_seed(seed, run_eval()) else run_eval()
  out <- eval_res$value
  K_obs_meta <- if (is_hetero) {
    as.integer(unname(K_by_level_obj[as.character(df[[object$spec$threshold_facet]])]))
  } else {
    NULL
  }
  attr(out, "bmfrm_prediction_meta") <- list(
    type = type,
    prediction_target = classify_prediction_target(adjust_facets, eval_res$new_levels, sample_new_levels),
    adjust_facets = adjust_facets,
    allow_new_levels = isTRUE(allow_new_levels),
    sample_new_levels = sample_new_levels,
    new_levels = eval_res$new_levels,
    K_heterogeneous = is_hetero,
    K_by_level = K_by_level_obj,
    K_obs = K_obs_meta,
    n_draws = D,
    draw_ids = draw_meta$draw_ids,
    draw_info = draw_meta$draw_info,
    seed = seed
  )
  out
}

#' Numerically stable log(mean(exp(x))) for a numeric vector
#' @keywords internal
logmeanexp_vec <- function(x) {
  m <- max(x)
  if (!is.finite(m)) return(m)
  m + log(mean(exp(x - m)))
}

#' Numerically stable row-wise log(mean(exp(X))) for a matrix
#' @keywords internal
logmeanexp_rows <- function(X) {
  m <- apply(X, 1L, max)
  out <- m + log(rowMeans(exp(X - m)))
  out[!is.finite(m)] <- m[!is.finite(m)]
  out
}

#' Single-observation log-likelihood evaluator for loo::loo.function() (O(D) memory)
#'
#' Expects `data_i` to be a 1-row data.frame containing integer facet indices
#' (`person`, plus each unadjusted non-person facet) and integer `.bmfrm_response`,
#' and `draws` to be the compact parameter list from `extract_bmfrm_param_draws()`.
#' Returns a length-D numeric vector of pointwise log-likelihoods without
#' materializing any D x N intermediate matrices.
#'
#' @keywords internal
eval_bmfrm_log_lik_obs <- function(data_i, draws) {
  D <- draws$D
  K <- draws$K
  y <- as.integer(data_i$.bmfrm_response[[1L]])

  pf <- draws$person_facet %||% "person"
  eta <- draws$theta[, as.integer(data_i[[pf]][[1L]])]
  for (f in names(draws$facet_draws)) {
    eta <- eta - draws$facet_draws[[f]][, as.integer(data_i[[f]][[1L]])]
  }
  for (b in draws$bias_draws) {
    i1 <- as.integer(data_i[[b$f1]][[1L]])
    i2 <- as.integer(data_i[[b$f2]][[1L]])
    cell_idx <- i1 + (i2 - 1L) * b$n1
    eta <- eta + b$draws[, cell_idx]
  }

  tau_info <- draws$tau_info
  if (identical(tau_info$mode, "pcm_ragged")) {
    t_idx <- as.integer(data_i[[tau_info$threshold_facet]][[1L]])
    K_i <- tau_info$K_level[t_idx]
    ct_mat <- tau_info$cum_tau_by_level[[t_idx]]
    max_logit <- numeric(D)
    if (K_i >= 2L) {
      for (k in 2:K_i) {
        lk <- (k - 1) * eta - ct_mat[, k]
        max_logit <- pmax(max_logit, lk)
      }
    }
    denom <- exp(-max_logit)
    if (K_i >= 2L) {
      for (k in 2:K_i) {
        denom <- denom + exp((k - 1) * eta - ct_mat[, k] - max_logit)
      }
    }
    ly <- if (y == 1L) numeric(D) else (y - 1) * eta - ct_mat[, y]
    return(as.numeric(ly - (log(denom) + max_logit)))
  }

  t_idx <- if (tau_info$mode == "pcm") as.integer(data_i[[tau_info$threshold_facet]][[1L]]) else NULL
  get_ct_k <- function(k) {
    if (tau_info$mode == "shared") {
      tau_info$cum_tau[, k]
    } else {
      tau_info$cum_tau_list[[k]][, t_idx]
    }
  }

  max_logit <- numeric(D)
  if (K >= 2L) {
    for (k in 2:K) {
      lk <- (k - 1) * eta - get_ct_k(k)
      max_logit <- pmax(max_logit, lk)
    }
  }
  denom <- exp(-max_logit)
  if (K >= 2L) {
    for (k in 2:K) {
      denom <- denom + exp((k - 1) * eta - get_ct_k(k) - max_logit)
    }
  }
  ly <- if (y == 1L) numeric(D) else (y - 1) * eta - get_ct_k(y)
  as.numeric(ly - (log(denom) + max_logit))
}

#' Extract verified iteration-by-chain layout without inventing chains
#' @keywords internal
extract_verified_chain_layout <- function(fit, D_total) {
  di <- extract_draw_index_info(fit, D_total, draw_idx = NULL)
  ch <- di$draw_info$chain
  it <- di$draw_info$iteration
  if (anyNA(ch) || anyNA(it)) {
    return(list(valid = FALSE, n_iter = NA_integer_, n_chains = NA_integer_, chain_id = NULL))
  }
  n_chains <- max(ch)
  n_iter <- max(it)
  if (n_chains < 1L || n_iter < 1L || (n_chains * n_iter != D_total)) {
    return(list(valid = FALSE, n_iter = NA_integer_, n_chains = NA_integer_, chain_id = NULL))
  }
  list(valid = TRUE, n_iter = as.integer(n_iter), n_chains = as.integer(n_chains), chain_id = ch)
}
