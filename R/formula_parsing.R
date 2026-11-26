#' Parse bmfrm formula to extract facets and model structure
#'
#' @param formula R formula (e.g., score ~ person + item + rater + (rater:item)).
#'   The response must be a single column name; transform scores in `data` first.
#' @param data data.frame containing variables
#' @param K Number of rating categories: a scalar integer `>= 2`, a named
#'   integer vector matching `threshold_facet` levels for heterogeneous PCM,
#'   or `NULL` (inferred from `max(score)`).
#' @param family model family: "rsm" (Rating Scale Model) or "pcm" (Partial Credit Model)
#' @param threshold_facet for PCM only: which facet gets item-specific thresholds (e.g., "item")
#' @param model_name optional model name for caching
#' @param cache_dir cache directory
#' @param hierarchical Names of facets or terms with hierarchical distributions
#'   (e.g., `"person"`, centered main facets such as `"rater"` / `"item"`,
#'   two-way interactions such as `"item:rater"`, or `"tau"` for hierarchical
#'   PCM). Component order in interaction names does not affect prior matching.
#' @param save_gq Either `"full"` (default; generate observation-level
#'   `log_lik`, `y_rep`, `mu`, `sigma2`, `mu_fair` in Stan `generated quantities`)
#'   or `"minimal"` (omit observation-level `generated quantities` and compute
#'   downstream quantities in memory-bounded R chunks).
#'
#' @return bmfrm_spec object with model specification
#' @export
parse_bmfrm_formula <- function(formula,
                               data,
                               K = NULL,
                               family = "rsm",
                               threshold_facet = NULL,
                               model_name = NULL,
                               cache_dir = "stan_cache",
                               hierarchical = character(),
                               save_gq = c("full", "minimal")) {

  # Validate family
  allowed_families <- c("rsm", "pcm")
  family <- match.arg(family, allowed_families)
  save_gq <- match.arg(save_gq)
  
  # Parse formula terms
  response_name <- bmfrm_response_name(formula)
  terms_obj <- stats::terms(formula)
  
  # Check response exists in data
  if (!response_name %in% names(data)) {
    stop("Response variable '", response_name, "' not found in data")
  }
  
  # Extract RHS variables
  rhs_vars <- attr(terms_obj, "term.labels")
  
  # Separate main facets from bias facets (interactions)
  is_interaction <- grepl(":", rhs_vars)
  facets_main <- rhs_vars[!is_interaction]
  facets_bias <- rhs_vars[is_interaction]
  
  # Remove parentheses from bias terms if present
  facets_bias <- gsub("[()]", "", facets_bias)
  validate_facet_names(facets_main)
  validate_interaction_terms(facets_main, facets_bias)
  
  # Check all variables exist in data
  all_vars <- unique(c(facets_main, unlist(strsplit(facets_bias, ":"))))
  missing_vars <- setdiff(all_vars, names(data))
  if (length(missing_vars) > 0) {
    stop("Variables not found in data: ", paste(missing_vars, collapse = ", "))
  }

  # Validate threshold_facet for PCM

  if (family == "pcm") {
    if (is.null(threshold_facet)) {
      stop("PCM family requires 'threshold_facet' to specify which facet gets ",
           "item-specific thresholds (e.g., threshold_facet = 'item')")
    }
    if (!threshold_facet %in% facets_main) {
      stop("threshold_facet '", threshold_facet, "' must be one of the main facets: ",
           paste(facets_main, collapse = ", "))
    }
  } else if (family == "rsm") {
    if (!is.null(threshold_facet)) {
      message("Note: 'threshold_facet' is ignored for RSM family (thresholds are shared)")
      threshold_facet <- NULL
    }
  }

  hierarchical <- validate_hierarchical(hierarchical, facets_main, facets_bias, family)

  # Infer K if not provided
  if (is.null(K)) {
    resp <- data[[response_name]]
    if (!all(resp == floor(resp), na.rm = TRUE)) {
      stop("Response must be integer for rating-scale MFRM.")
    }
    K <- as.integer(max(resp, na.rm = TRUE))
    if (min(resp, na.rm = TRUE) < 1) {
      stop("Response values must start from 1, not 0. Add 1 to your scores if needed.")
    }
    message(
      "K was not specified; inferred K = ", K, " from max(", response_name, "). ",
      "If your rating rubric defines higher categories that did not occur in this sample ",
      "(or uses different category counts across ", threshold_facet %||% "items", "), ",
      "pass K explicitly (e.g., K = 6 or a named vector for PCM) so the intended scale is modeled."
    )
  }

  k_info <- resolve_category_spec(
    K = K,
    data = data,
    family = family,
    threshold_facet = threshold_facet,
    response_name = response_name,
    check_response = TRUE,
    drop_na_response_levels = TRUE
  )

  if (isTRUE(k_info$K_heterogeneous) && ("tau" %in% hierarchical)) {
    stop(
      "Hierarchical step sharing (hierarchical = 'tau') is not supported in v1 when PCM category ",
      "counts K vary across '", threshold_facet, "' levels. Step indices do not correspond across ",
      "different rubric lengths; omit 'tau' from hierarchical (other facets and interactions may ",
      "still be hierarchical) or use uniform K across levels."
    )
  }

  # Generate model name if not provided
  if (is.null(model_name)) {
    # Create a hash based on formula and facets
    facet_string <- paste(c(facets_main, facets_bias), collapse = "_")
    model_name <- paste0("bmfrm_", substr(digest::digest(facet_string), 1, 8))
  }
  
  # Create spec object
  spec <- list(
    response        = response_name,
    facets_main     = facets_main,
    facets_bias     = facets_bias,
    K               = k_info$K,
    K_by_level      = k_info$K_by_level,
    K_max           = k_info$K_max,
    K_heterogeneous = k_info$K_heterogeneous,
    family          = family,
    hierarchical    = hierarchical,
    threshold_facet = threshold_facet,
    save_gq         = save_gq,
    formula         = formula,
    model_name      = model_name,
    cache_dir       = cache_dir
  )
  
  class(spec) <- "bmfrm_spec"
  spec
}

#' Resolve and validate scalar or per-level category count specification K
#' @keywords internal
resolve_category_spec <- function(K,
                                  data,
                                  family,
                                  threshold_facet = NULL,
                                  response_name = NULL,
                                  check_response = FALSE,
                                  drop_na_response_levels = FALSE) {
  is_named_or_vec <- length(K) > 1L || (!is.null(names(K)) && length(names(K)) > 0L)

  if (!is_named_or_vec) {
    validate_category_count(K)
    K_int <- as.integer(K)
    if (identical(family, "pcm")) {
      tf_vec <- data[[threshold_facet]]
      if (isTRUE(drop_na_response_levels) && !is.null(response_name) &&
          response_name %in% names(data) && !is.factor(tf_vec)) {
        non_na_resp <- !is.na(data[[response_name]])
        tf_fac <- if (any(non_na_resp)) factor(tf_vec[non_na_resp]) else factor(tf_vec)
      } else {
        tf_fac <- if (is.factor(tf_vec)) tf_vec else factor(tf_vec)
      }
      tf_levels <- levels(tf_fac)
      K_by_level <- stats::setNames(rep(K_int, length(tf_levels)), tf_levels)
    } else {
      K_by_level <- NULL
    }
    K_max <- K_int
    K_heterogeneous <- FALSE
    K_out <- K
  } else {
    if (!identical(family, "pcm")) {
      stop(
        "Category-count vector K is only supported for family = 'pcm' ",
        "(family = 'rsm' requires a single shared category count K)."
      )
    }
    if (is.null(threshold_facet) || !threshold_facet %in% names(data)) {
      stop("Named vector K requires a valid threshold_facet column in data")
    }
    if (!is.numeric(K) || any(!is.finite(K)) || any(K < 2) || any(K != floor(K))) {
      stop("All elements of K must be finite integers of at least 2")
    }
    nm <- names(K)
    if (is.null(nm) || anyNA(nm) || any(!nzchar(nm))) {
      stop("Vector K for PCM must be a named integer vector matching threshold_facet levels")
    }
    if (anyDuplicated(nm)) {
      stop("Vector K contains duplicate level names: ", paste(unique(nm[duplicated(nm)]), collapse = ", "))
    }

    tf_vec <- data[[threshold_facet]]
    tf_fac <- if (is.factor(tf_vec)) tf_vec else factor(tf_vec)
    tf_levels <- levels(tf_fac)

    missing_lv <- setdiff(tf_levels, nm)
    if (length(missing_lv) > 0L) {
      stop(
        "Named vector K is missing level(s) of threshold_facet '", threshold_facet, "': ",
        paste(missing_lv, collapse = ", ")
      )
    }
    extra_lv <- setdiff(nm, tf_levels)
    if (length(extra_lv) > 0L) {
      stop(
        "Named vector K contains unknown level(s) not present in threshold_facet '",
        threshold_facet, "': ", paste(extra_lv, collapse = ", ")
      )
    }

    # After validating input rubric names against full data, project onto levels
    # that survive missing-score filtering (for non-factor columns) so uniform
    # remaining levels (including binary K = 2) normalize to the uniform PCM path.
    active_tf_levels <- tf_levels
    if (isTRUE(drop_na_response_levels) && !is.null(response_name) &&
        response_name %in% names(data) && !is.factor(tf_vec)) {
      non_na_resp <- !is.na(data[[response_name]])
      if (any(non_na_resp)) {
        active_tf_levels <- levels(factor(tf_vec[non_na_resp]))
      }
    }

    K_by_level <- stats::setNames(as.integer(K[active_tf_levels]), active_tf_levels)
    K_max <- max(K_by_level)
    K_heterogeneous <- length(unique(K_by_level)) > 1L
    K_out <- if (K_heterogeneous) K_by_level else as.integer(K_by_level[1L])
  }

  if (isTRUE(check_response) && !is.null(response_name) && response_name %in% names(data)) {
    resp <- data[[response_name]]
    non_na <- !is.na(resp)
    if (!is.numeric(resp) || any(!is.finite(resp[non_na])) ||
        any(resp[non_na] != floor(resp[non_na]))) {
      stop("Response must contain finite integer scores")
    }
    if (identical(family, "pcm") && !is.null(K_by_level)) {
      obs_tf <- as.character(data[[threshold_facet]][non_na])
      if (anyNA(obs_tf)) {
        stop("Missing values found in threshold_facet '", threshold_facet, "'. Please clean your data first.")
      }
      obs_K <- unname(K_by_level[obs_tf])
      bad <- resp[non_na] < 1 | resp[non_na] > obs_K
      if (any(bad)) {
        bad_lv <- unique(obs_tf[bad])
        stop(
          "Response values exceed valid category range 1..K_j for threshold_facet '",
          threshold_facet, "' level(s): ",
          paste(sprintf("%s (1..%d)", bad_lv, K_by_level[bad_lv]), collapse = ", ")
        )
      }
    } else {
      if (any(resp[non_na] < 1 | resp[non_na] > K_max)) {
        stop("Response values must be between 1 and K")
      }
    }
  }

  list(
    K = K_out,
    K_by_level = K_by_level,
    K_max = as.integer(K_max),
    K_heterogeneous = isTRUE(K_heterogeneous)
  )
}

bmfrm_response_name <- function(formula) {
  if (!inherits(formula, "formula") || length(formula) != 3L || !is.symbol(formula[[2L]])) {
    stop("The formula response must be a single column name; transform scores in data before fitting")
  }
  as.character(formula[[2L]])
}

#' Get maximum category count from a bmfrm_spec
#' @keywords internal
get_spec_K_max <- function(spec) {
  if (!is.null(spec$K_max)) return(as.integer(spec$K_max))
  as.integer(max(spec$K))
}

#' Print method for bmfrm_spec
#' @param x bmfrm_spec object
#' @param ... unused
#' @export
print.bmfrm_spec <- function(x, ...) {
  cat("BMFRM Model Specification\n")
  cat("Formula:", deparse(x$formula), "\n")
  cat("Response:", x$response, "\n")
  if (isTRUE(x$K_heterogeneous) && !is.null(x$K_by_level)) {
    cat("Rating categories: heterogeneous by ", x$threshold_facet, " (",
        paste(sprintf("%s=1..%d", names(x$K_by_level), as.integer(x$K_by_level)), collapse = ", "),
        ")\n", sep = "")
  } else {
    cat("Rating categories: 1 to", get_spec_K_max(x), "\n")
  }
  family_desc <- switch(x$family,
    "rsm" = "RSM (Rating Scale Model - shared thresholds)",
    "pcm" = "PCM (Partial Credit Model - facet-specific thresholds)",
    x$family
  )
  cat("Family:", family_desc, "\n")
  if (x$family == "pcm" && !is.null(x$threshold_facet)) {
    cat("Threshold facet:", x$threshold_facet, "\n")
  }
  cat("Main facets:", paste(x$facets_main, collapse = ", "), "\n")
  if (length(x$facets_bias) > 0) {
    cat("Bias facets:", paste(x$facets_bias, collapse = ", "), "\n")
  }
  cat("Model name:", x$model_name, "\n")
}

#' Validate bmfrm model specification
#'
#' @param spec bmfrm_spec object
#' @param data data.frame
#'
#' @return TRUE if valid (invisibly), otherwise stops with error
#' @keywords internal
validate_spec <- function(spec, data) {
  # Check required facets
  if (!"person" %in% spec$facets_main) {
    stop("Model must include 'person' as a main facet")
  }
  
  # Check for at least one item-like facet
  item_facets <- intersect(c("item", "criterion", "task"), spec$facets_main)
  if (length(item_facets) == 0) {
    message("Note: No standard item facet found. Make sure you have an appropriate difficulty facet.")
  }
  
  # Validate bias facets reference existing main facets
  if (length(spec$facets_bias) > 0) {
    all_bias_components <- unique(unlist(strsplit(spec$facets_bias, ":")))
    missing_components <- setdiff(all_bias_components, spec$facets_main)
    if (length(missing_components) > 0) {
      stop("Bias facets reference undefined main facets: ", 
           paste(missing_components, collapse = ", "))
    }
  }
  
  # Check data has sufficient observations
  if (nrow(data) < 10) {
    warning("Very few observations (", nrow(data), "). Results may be unreliable.")
  }
  
  invisible(TRUE)
}
