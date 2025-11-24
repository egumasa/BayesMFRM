#' Parse bmfrm formula to extract facets and model structure
#'
#' @param formula R formula (e.g., score ~ person + item + rater + (rater:item))
#' @param data data.frame containing variables
#' @param K number of rating categories (if NULL, inferred from data)
#' @param family model family (currently only "rating_scale")
#' @param model_name optional model name for caching
#' @param cache_dir cache directory
#'
#' @return bmfrm_spec object with model specification
#' @export
parse_bmfrm_formula <- function(formula, 
                               data, 
                               K = NULL, 
                               family = "rating_scale",
                               model_name = NULL,
                               cache_dir = "stan_cache") {
  
  # Validate family
  allowed_families <- c("rating_scale")
  family <- match.arg(family, allowed_families)
  
  # Parse formula terms
  terms_obj <- stats::terms(formula)
  response_name <- all.vars(formula)[1]
  
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
  
  # Check all variables exist in data
  all_vars <- unique(c(facets_main, unlist(strsplit(facets_bias, ":"))))
  missing_vars <- setdiff(all_vars, names(data))
  if (length(missing_vars) > 0) {
    stop("Variables not found in data: ", paste(missing_vars, collapse = ", "))
  }
  
  # Infer K if not provided
  if (is.null(K)) {
    resp <- data[[response_name]]
    if (!all(resp == floor(resp), na.rm = TRUE)) {
      stop("Response must be integer for rating-scale MFRM.")
    }
    K <- max(resp, na.rm = TRUE)
    if (min(resp, na.rm = TRUE) < 1) {
      stop("Response values must start from 1, not 0. Add 1 to your scores if needed.")
    }
  }
  
  # Generate model name if not provided
  if (is.null(model_name)) {
    # Create a hash based on formula and facets
    facet_string <- paste(c(facets_main, facets_bias), collapse = "_")
    model_name <- paste0("bmfrm_", substr(digest::digest(facet_string), 1, 8))
  }
  
  # Create spec object
  spec <- list(
    response     = response_name,
    facets_main  = facets_main,
    facets_bias  = facets_bias,
    K            = K,
    family       = family,
    formula      = formula,
    model_name   = model_name,
    cache_dir    = cache_dir
  )
  
  class(spec) <- "bmfrm_spec"
  spec
}

#' Print method for bmfrm_spec
#' @param x bmfrm_spec object
#' @param ... unused
#' @export
print.bmfrm_spec <- function(x, ...) {
  cat("BMFRM Model Specification\n")
  cat("Formula:", deparse(x$formula), "\n")
  cat("Response:", x$response, "\n")
  cat("Rating categories: 1 to", x$K, "\n")
  cat("Family:", x$family, "\n")
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
