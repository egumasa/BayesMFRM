#' Create a prior specification for bmfrm models
#'
#' @param spec Prior distribution string. For a population mean or SD, a numeric
#'   value fixes that parameter (SD must be positive).
#' @param class Parameter class: theta, item, rater, bias, tau, facet_<name>,
#'   mean, or sd. The last two specify a hierarchical population.
#' @param facet Target facet for mean/sd, or interaction name for bias.
#' @param scale Meaning of the SD for a constrained Gaussian facet or step prior:
#'   "raw" retains the historical pre-centering SD; "marginal" specifies the
#'   constrained level SD. Marginal scaling requires a normal(0, sd) prior and
#'   generated additive Gaussian RSM. Singleton blocks remain identically zero.
#'
#' @return bmfrm_prior object
#' @export
#'
#' @examples
#' # Basic priors
#' prior("normal(0, 2)", class = "theta")
#' prior("normal(0, 1)", class = "item")
#' prior("student_t(3, 0, 1)", class = "rater")
#' prior("normal(0, 1)", class = "tau")
#'
#' # Bias prior (requires facet specification)
#' prior("normal(0, 0.5)", class = "bias", facet = "rater:item")
prior <- function(spec, class, facet = NULL, scale = c("raw", "marginal")) {
  scale <- match.arg(scale)
  # Validate class
  allowed_classes <- c("theta", "item", "rater", "bias", "tau", "mean", "sd")
  if (!is.character(class) || length(class) != 1L || is.na(class)) stop("class must be one parameter class")
  if (!class %in% allowed_classes && !grepl("^facet_", class)) {
    stop("Invalid class '", class, "'. Must be one of: ", 
         paste(allowed_classes, collapse = ", "), 
         " or 'facet_<name>'")
  }
  
  # Validate bias class requires facet
  if (class == "bias" && is.null(facet)) {
    stop("Class 'bias' requires 'facet' specification (e.g., facet = 'rater:item')")
  }
  if (!class %in% c("bias", "mean", "sd") && !is.null(facet)) stop("facet is only used for bias and population priors")
  
  # Validate facet format for bias
  if (class == "bias" && !is.null(facet)) {
    if (!is.character(facet) || length(facet) != 1L || is.na(facet) || !grepl(":", facet)) {
      stop("Bias facet must be an interaction term (e.g., 'rater:item')")
    }
  }
  
  # Validate spec format - basic check for distribution name and parentheses
  population <- class %in% c("mean", "sd")
  if (population && (!is.character(facet) || length(facet) != 1L || is.na(facet) || !nzchar(facet))) stop("Population priors require a target facet")
  if (population) {
    rsm_population_prior(spec, if (class == "mean") "normal" else "lognormal", class == "sd")
  } else if (!is.character(spec) || length(spec) != 1L || is.na(spec) ||
             !grepl("^[a-zA-Z_]+\\s*\\(.*\\)$", spec)) {
    stop("Invalid prior specification '", spec, "'. Must be in format 'distribution(parameters)'")
  }
  if (scale == "marginal") {
    if (class %in% c("theta", "bias", "mean", "sd")) stop("Marginal scaling applies only to constrained main facets and tau")
    p <- parse_bmfrm_prior_literal(spec)
    if (p$distribution != "normal" || length(p$args) != 2L || p$args[1L] != 0 || p$args[2L] <= 0) stop("Marginal scaling requires normal(0, positive_sd)")
  }
  
  structure(
    list(
      spec  = spec,    # e.g. "normal(0,2)"
      class = class,   # e.g. "theta", "item", "rater", "bias", "tau"
      facet = facet,   # e.g. "rater:item" (only for bias)
      scale = scale
    ),
    class = "bmfrm_prior"
  )
}

#' Print method for bmfrm_prior
#' @param x bmfrm_prior object
#' @param ... unused
#' @export
print.bmfrm_prior <- function(x, ...) {
  if (x$class %in% c("bias", "mean", "sd")) {
    cat("Prior: ", x$spec, " for ", x$class, " (", x$facet, ")\n", sep = "")
  } else {
    cat("Prior: ", x$spec, " for ", x$class, "\n", sep = "")
  }
  if (identical(x$scale, "marginal")) cat("  SD refers to each constrained level (singleton blocks are zero).\n")
}

#' Combine prior specifications without flattening their fields
#' @param ... Prior objects or lists of prior objects.
#' @param recursive Unused; recursive concatenation is not supported.
#' @return A list of bmfrm_prior objects.
#' @export
c.bmfrm_prior <- function(..., recursive = FALSE) {
  if (recursive) stop("Recursive prior concatenation is not supported")
  pieces <- lapply(list(...), function(x) {
    if (inherits(x, "bmfrm_prior")) return(list(x))
    if (is.list(x) && all(vapply(x, inherits, logical(1), "bmfrm_prior"))) return(x)
    stop("Combine only prior objects or lists of priors")
  })
  unlist(pieces, recursive = FALSE, use.names = FALSE)
}

#' Get default priors for bmfrm models
#'
#' @param facets_main vector of main facet names
#' @param facets_bias vector of bias facet names
#' @param family model family: "rsm" or "pcm"
#' @param hierarchical Names of facets with hierarchical distributions. Currently
#'   only "person" is supported; character(0) retains independent priors.
#'
#' @return list of bmfrm_prior objects
#' @export
get_default_priors <- function(facets_main = c("person", "item", "rater"),
                               facets_bias = NULL,
                               family = "rsm", hierarchical = character()) {
  hierarchical <- validate_hierarchical(hierarchical, facets_main, facets_bias, family)
  priors <- list()
  
  # Default priors for main facets
  for (facet in facets_main) {
    if (facet == "person") {
      if ("person" %in% hierarchical) {
        priors <- c(priors, list(prior("normal(0, 2.5)", "mean", facet = "person"), prior("lognormal(0, 0.5)", "sd", facet = "person")))
      } else priors[[length(priors) + 1]] <- prior("normal(0, 2)", class = "theta")
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
  priors[[length(priors) + 1]] <- prior("normal(0, 1)", class = "tau")
  
  priors
}

#' Validate a list of priors
#'
#' @param priors list of bmfrm_prior objects or NULL
#' @param facets_main vector of main facet names
#' @param facets_bias vector of bias facet names
#' @param family model family: "rsm" or "pcm"
#' @param hierarchical Names of facets with hierarchical distributions. Currently
#'   only "person" is supported; character(0) retains independent priors.
#'
#' @return validated list of priors
#' @keywords internal
validate_priors <- function(priors, facets_main, facets_bias = NULL, family = "rsm", hierarchical = character()) {
  defaults <- get_default_priors(facets_main, facets_bias, family, hierarchical)
  if (is.null(priors)) priors <- defaults
  if (!is.list(priors) || !length(priors)) stop("priors must be a list of prior() objects")
  if (!all(vapply(priors, inherits, logical(1), "bmfrm_prior"))) stop("All elements of priors must be prior() objects")
  # Revalidate edited objects and normalize saved priors without the scale field.
  priors <- lapply(priors, function(p) prior(p$spec, p$class, p$facet, p$scale %||% "raw"))
  key <- function(p) paste(p$class, p$facet %||% "", sep = ":")
  keys <- vapply(priors, key, character(1))
  required <- vapply(defaults, key, character(1))
  if (anyDuplicated(keys)) stop("Duplicate prior for the same parameter class/facet")
  extra <- setdiff(keys, required)
  if (length(extra)) stop("Prior classes/facets not used by this model: ", paste(extra, collapse = ", "))
  missing <- setdiff(required, keys)
  if (length(missing)) stop("Missing priors for classes/facets: ", paste(missing, collapse = ", "))
  priors
}
