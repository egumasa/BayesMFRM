#' Create a prior specification for bmfrm models
#'
#' @param spec prior specification string (e.g., "normal(0, 2)")
#' @param class parameter class ("theta", "item", "rater", "bias", "tau", or "facet_<name>")
#' @param facet interaction facet name (required for class = "bias")
#'
#' @return bmfrm_prior object
#' @export
#'
#' @examples
#' # Basic priors
#' prior("normal(0, 2)", class = "theta")
#' prior("normal(0, 1)", class = "item")
#' prior("student_t(3, 0, 1)", class = "rater")
#' prior("normal(0, 3)", class = "tau")
#'
#' # Bias prior (requires facet specification)
#' prior("normal(0, 0.5)", class = "bias", facet = "rater:item")
prior <- function(spec, class, facet = NULL) {
  # Validate class
  allowed_classes <- c("theta", "item", "rater", "bias", "tau")
  if (!class %in% allowed_classes && !grepl("^facet_", class)) {
    stop("Invalid class '", class, "'. Must be one of: ", 
         paste(allowed_classes, collapse = ", "), 
         " or 'facet_<name>'")
  }
  
  # Validate bias class requires facet
  if (class == "bias" && is.null(facet)) {
    stop("Class 'bias' requires 'facet' specification (e.g., facet = 'rater:item')")
  }
  
  # Validate facet format for bias
  if (class == "bias" && !is.null(facet)) {
    if (!grepl(":", facet)) {
      stop("Bias facet must be an interaction term (e.g., 'rater:item')")
    }
  }
  
  # Validate spec format - basic check for distribution name and parentheses
  if (!grepl("^[a-zA-Z_]+\\s*\\(.*\\)$", spec)) {
    stop("Invalid prior specification '", spec, "'. Must be in format 'distribution(parameters)'")
  }
  
  structure(
    list(
      spec  = spec,    # e.g. "normal(0,2)"
      class = class,   # e.g. "theta", "item", "rater", "bias", "tau"
      facet = facet    # e.g. "rater:item" (only for bias)
    ),
    class = "bmfrm_prior"
  )
}

#' Print method for bmfrm_prior
#' @param x bmfrm_prior object
#' @param ... unused
#' @export
print.bmfrm_prior <- function(x, ...) {
  if (x$class == "bias") {
    cat("Prior: ", x$spec, " for ", x$class, " (", x$facet, ")\n", sep = "")
  } else {
    cat("Prior: ", x$spec, " for ", x$class, "\n", sep = "")
  }
}

#' Get default priors for bmfrm models
#'
#' @param facets_main vector of main facet names
#' @param facets_bias vector of bias facet names
#'
#' @return list of bmfrm_prior objects
#' @export
get_default_priors <- function(facets_main = c("person", "item", "rater"), 
                               facets_bias = NULL) {
  priors <- list()
  
  # Default priors for main facets
  for (facet in facets_main) {
    if (facet == "person") {
      priors[[length(priors) + 1]] <- prior("normal(0, 2)", class = "theta")
    } else if (facet == "item") {
      priors[[length(priors) + 1]] <- prior("normal(0, 1)", class = "item")
    } else if (facet == "rater") {
      priors[[length(priors) + 1]] <- prior("normal(0, 1)", class = "rater")
    } else {
      # Other main facets get a generic name
      class_name <- paste0("facet_", facet)
      priors[[length(priors) + 1]] <- prior("normal(0, 1)", class = class_name)
    }
  }
  
  # Default priors for bias facets
  if (!is.null(facets_bias)) {
    for (bias_facet in facets_bias) {
      priors[[length(priors) + 1]] <- prior("normal(0, 0.5)", 
                                             class = "bias", 
                                             facet = bias_facet)
    }
  }
  
  # Default prior for thresholds
  priors[[length(priors) + 1]] <- prior("normal(0, 3)", class = "tau")
  
  priors
}

#' Validate a list of priors
#'
#' @param priors list of bmfrm_prior objects or NULL
#' @param facets_main vector of main facet names
#' @param facets_bias vector of bias facet names
#'
#' @return validated list of priors
#' @keywords internal
validate_priors <- function(priors, facets_main, facets_bias = NULL) {
  if (is.null(priors)) {
    return(get_default_priors(facets_main, facets_bias))
  }
  
  if (!is.list(priors)) {
    stop("priors must be a list of prior() objects")
  }
  
  # Check all elements are bmfrm_prior objects
  is_prior <- sapply(priors, inherits, "bmfrm_prior")
  if (!all(is_prior)) {
    stop("All elements of priors must be prior() objects")
  }
  
  # Extract classes and facets from priors
  prior_classes <- sapply(priors, function(p) p$class)
  prior_facets <- sapply(priors, function(p) p$facet %||% NA_character_)
  
  # Check for required priors
  required_classes <- c("theta", "tau")
  for (facet in facets_main) {
    if (facet == "person") {
      # theta already in required_classes
    } else if (facet == "item") {
      required_classes <- c(required_classes, "item")
    } else if (facet == "rater") {
      required_classes <- c(required_classes, "rater") 
    } else {
      required_classes <- c(required_classes, paste0("facet_", facet))
    }
  }
  
  missing_classes <- setdiff(required_classes, prior_classes)
  if (length(missing_classes) > 0) {
    stop("Missing priors for classes: ", paste(missing_classes, collapse = ", "))
  }
  
  # Check bias priors have valid facet specifications
  bias_priors <- priors[prior_classes == "bias"]
  if (length(bias_priors) > 0) {
    bias_facet_specs <- sapply(bias_priors, function(p) p$facet)
    if (!is.null(facets_bias)) {
      missing_bias <- setdiff(facets_bias, bias_facet_specs)
      if (length(missing_bias) > 0) {
        stop("Missing bias priors for facets: ", paste(missing_bias, collapse = ", "))
      }
    }
  }
  
  priors
}
