#' Validate a Stan univariate continuous prior specification (AST-guarded)
#' @keywords internal
validate_stan_prior_spec <- function(spec) {
  if (!is.character(spec) || length(spec) != 1L || is.na(spec) || !nzchar(trimws(spec))) {
    stop("Invalid prior specification. Must be a single non-empty distribution string such as 'normal(0, 1)'")
  }
  if (grepl(";", spec, fixed = TRUE)) {
    stop("Invalid prior specification '", spec, "': semicolons and multiple statements are not allowed")
  }
  parsed <- tryCatch(parse(text = spec, keep.source = FALSE), error = function(e) NULL)
  if (is.null(parsed) || length(parsed) != 1L) {
    stop("Invalid prior specification '", spec, "'. Must be a single call in format 'distribution(parameters)'")
  }
  p <- tryCatch(
    parse_bmfrm_prior_literal(spec),
    error = function(e) stop("Invalid prior specification '", spec, "': ", conditionMessage(e), call. = FALSE)
  )
  dist <- p$distribution
  args <- p$args
  n_a <- length(args)

  ok <- switch(
    dist,
    std_normal         = (n_a == 0L),
    normal             = (n_a == 2L && args[2L] > 0),
    logistic           = (n_a == 2L && args[2L] > 0),
    cauchy             = (n_a == 2L && args[2L] > 0),
    double_exponential = (n_a == 2L && args[2L] > 0),
    gumbel             = (n_a == 2L && args[2L] > 0),
    lognormal          = (n_a == 2L && args[2L] > 0),
    student_t          = (n_a == 3L && args[1L] > 0 && args[3L] > 0),
    gamma              = (n_a == 2L && args[1L] > 0 && args[2L] > 0),
    inv_gamma          = (n_a == 2L && args[1L] > 0 && args[2L] > 0),
    weibull            = (n_a == 2L && args[1L] > 0 && args[2L] > 0),
    frechet            = (n_a == 2L && args[1L] > 0 && args[2L] > 0),
    exponential        = (n_a == 1L && args[1L] > 0),
    rayleigh           = (n_a == 1L && args[1L] > 0),
    uniform            = (n_a == 2L && args[1L] < args[2L]),
    FALSE
  )
  if (!isTRUE(ok)) {
    stop(
      "Invalid or unsupported Stan prior '", spec,
      "'. Supported distributions require literal numeric arguments with valid positive scale/shape parameters."
    )
  }
  p
}

#' Normalize a prior object's lookup key across canonical and legacy class aliases
#' @keywords internal
prior_target_key <- function(p, person_facet = "person", facets_bias = character()) {
  cl <- p$class
  fc <- p$facet %||% ""
  if (nzchar(fc) && grepl(":", fc, fixed = TRUE) && length(facets_bias)) {
    fc <- match_interaction_names(fc, facets_bias)
  }
  if (cl %in% c("item", "rater")) {
    return(paste0("facet:", cl))
  }
  if (grepl("^facet_", cl)) {
    return(paste0("facet:", sub("^facet_", "", cl)))
  }
  if (identical(cl, "facet")) {
    return(paste0("facet:", fc))
  }
  if (cl %in% c("mean", "sd") && (fc == "" || identical(fc, "person"))) {
    return(paste(cl, person_facet, sep = ":"))
  }
  paste(cl, fc, sep = ":")
}

#' Create a prior specification for bmfrm models
#'
#' @param spec Prior distribution string (e.g., `"normal(0, 1)"`), validated as
#'   a single literal distribution call. For a population mean or SD, a numeric
#'   scalar fixes that parameter (SD must be positive).
#' @param class Parameter class:
#'   \itemize{
#'     \item `"facet"` (with `facet = "<name>"`): fixed-scale prior for any
#'       sum-to-zero main facet (recommended role-agnostic form; legacy aliases
#'       `"item"`, `"rater"`, and `"facet_<name>"` remain supported).
#'     \item `"theta"`: non-hierarchical ability prior when `person_facet` is omitted
#'       from `hierarchical`.
#'     \item `"tau"`: category step threshold prior (or common step prior in
#'       hierarchical PCM).
#'     \item `"bias"` (with `facet = "f1:f2"`): non-hierarchical two-way interaction prior.
#'     \item `"mean"`, `"sd"`: population mean (`person_facet` only) or hierarchical
#'       scale (`person_facet`, any centered main facet, `"tau"`, or `"f1:f2"`).
#'   }
#' @param facet Target facet name for `class = "facet"`, `"mean"`, `"sd"`, or
#'   two-way interaction term for `class = "bias"` / `"sd"`.
#' @param scale Meaning of the SD for a constrained Gaussian facet, step, or
#'   interaction prior: `"raw"` (pre-centering SD) or `"marginal"` (post-centering
#'   per-level / per-cell SD; requires `normal(0, sd)` or hierarchical `class = "sd"`).
#' @param local_shrinkage Optional local-global shrinkage structure for
#'   hierarchical two-way interactions (`class = "sd"`, `facet = "f1:f2"`):
#'   `"none"`, `"student_t"` (unit-variance Student-t mixture, df = 3), or
#'   `"horseshoe"` (regularized horseshoe in latent cell space; requires `scale = "raw"`).
#'   Note: local shrinkage acts on latent cells prior to double centering
#'   `B = P1 U P2`, so sparse latent deviations redistribute small centering
#'   offsets across rows/columns of `B`; combine `delta > 0` and BFDR control in
#'   [summarise_bias()] when screening for substantive bias.
#' @param slab_df Positive degrees of freedom for the regularized horseshoe slab
#'   `c^2 ~ inv_gamma(slab_df / 2, slab_df * slab_scale^2 / 2)` (default `4`).
#' @param slab_scale Positive scale for the regularized horseshoe slab (default
#'   `2`, yielding `c^2 ~ inv_gamma(2, 8)` with `slab_df = 4`).
#' @param person_facet Name of the ability facet when specifying `class = "mean"`
#'   with a custom examinee column name (default `"person"`).
#'
#' @return A `bmfrm_prior` object.
#' @export
prior <- function(spec, class, facet = NULL, scale = c("raw", "marginal"),
                  local_shrinkage = c("none", "horseshoe", "student_t"),
                  slab_df = 4, slab_scale = 2, person_facet = "person") {
  scale <- match.arg(scale)
  local_shrinkage <- match.arg(local_shrinkage)

  allowed_classes <- c("theta", "facet", "item", "rater", "bias", "tau", "mean", "sd")
  if (!is.character(class) || length(class) != 1L || is.na(class)) stop("class must be one parameter class")
  if (!class %in% allowed_classes && !grepl("^facet_", class)) {
    stop("Invalid class '", class, "'. Must be one of: ",
         paste(allowed_classes, collapse = ", "),
         " or 'facet_<name>'")
  }
  if (!is.character(person_facet) || length(person_facet) != 1L || is.na(person_facet) || !nzchar(person_facet)) {
    stop("person_facet must be a single non-empty facet name")
  }

  if (identical(class, "facet")) {
    if (!is.character(facet) || length(facet) != 1L || is.na(facet) || !nzchar(facet) || grepl(":", facet, fixed = TRUE)) {
      stop("Class 'facet' requires a single main facet name in 'facet' (e.g., facet = 'judge')")
    }
  } else if (class == "bias") {
    if (is.null(facet)) {
      stop("Class 'bias' requires 'facet' specification (e.g., facet = 'rater:item')")
    }
    if (!is.character(facet) || length(facet) != 1L || is.na(facet) || !grepl(":", facet, fixed = TRUE)) {
      stop("Bias facet must be an interaction term (e.g., 'rater:item')")
    }
  } else if (class %in% c("mean", "sd")) {
    if (is.null(facet) && identical(class, "mean")) {
      facet <- person_facet
    }
    if (!is.character(facet) || length(facet) != 1L || is.na(facet) || !nzchar(facet)) {
      stop("Population priors require a target facet")
    }
  } else if (!is.null(facet)) {
    stop("facet is only used for facet, bias, and population (mean/sd) priors")
  }

  if (!is.numeric(slab_df) || length(slab_df) != 1L || !is.finite(slab_df) || slab_df <= 0) {
    stop("slab_df must be a finite positive number")
  }
  if (!is.numeric(slab_scale) || length(slab_scale) != 1L || !is.finite(slab_scale) || slab_scale <= 0) {
    stop("slab_scale must be a finite positive number")
  }

  population <- class %in% c("mean", "sd")
  interaction_sd <- population && grepl(":", facet, fixed = TRUE)
  person_pop <- population && (facet %in% c("person", person_facet))
  centered_sd <- population && !person_pop

  if (local_shrinkage != "none" && !interaction_sd) {
    stop("local_shrinkage is only supported for hierarchical two-way interaction priors (class = 'sd', facet = 'f1:f2')")
  }
  if (local_shrinkage == "horseshoe" && scale == "marginal") {
    stop("local_shrinkage = 'horseshoe' uses a distribution-specific global shrinkage scale rather than a marginal cell SD; use scale = 'raw'")
  }

  if (population) {
    if (centered_sd) {
      if (class != "sd") {
        stop("Centered facets, steps, and interactions have no population mean; specify class = 'sd'")
      }
      positive_sd_prior(spec)
    } else {
      rsm_population_prior(spec, if (class == "mean") "normal" else "lognormal", class == "sd")
    }
  } else {
    validate_stan_prior_spec(spec)
  }

  if (scale == "marginal") {
    if (class %in% c("theta", "mean") || (class == "sd" && person_pop)) {
      stop("Marginal scaling applies only to constrained main facets, interactions and tau")
    }
    if (!centered_sd) {
      p <- parse_bmfrm_prior_literal(spec)
      if (p$distribution != "normal" || length(p$args) != 2L || p$args[1L] != 0 || p$args[2L] <= 0) {
        stop("Marginal scaling requires normal(0, positive_sd)")
      }
    }
  }

  obj <- list(
    spec  = spec,
    class = class,
    facet = facet,
    scale = scale
  )
  if (local_shrinkage != "none") {
    obj$local_shrinkage <- local_shrinkage
    if (local_shrinkage == "horseshoe") {
      obj$slab_df <- as.numeric(slab_df)
      obj$slab_scale <- as.numeric(slab_scale)
    }
  }
  structure(obj, class = "bmfrm_prior")
}

#' Print method for bmfrm_prior
#' @param x bmfrm_prior object
#' @param ... unused
#' @export
print.bmfrm_prior <- function(x, ...) {
  if (x$class %in% c("facet", "bias", "mean", "sd") && !is.null(x$facet)) {
    cat("Prior: ", x$spec, " for ", x$class, " (", x$facet, ")\n", sep = "")
  } else {
    cat("Prior: ", x$spec, " for ", x$class, "\n", sep = "")
  }
  if (identical(x$scale, "marginal")) cat("  SD refers to each constrained level (singleton blocks are zero).\n")
  if (!is.null(x$local_shrinkage) && x$local_shrinkage != "none") {
    desc <- if (x$local_shrinkage == "student_t") {
      "student_t (df = 3, unit-variance normalized)"
    } else {
      sprintf("horseshoe (slab_df = %g, slab_scale = %g)", x$slab_df %||% 4, x$slab_scale %||% 2)
    }
    cat("  Local shrinkage:", desc, "\n")
  }
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
#' By default, `person_facet` is partially pooled (`hierarchical = person_facet`)
#' with `mean ~ normal(0, 2.5)` and `sd ~ lognormal(0, 0.5)`, centered main
#' facets receive `normal(0, 1)`, two-way interactions receive `normal(0, 0.5)`,
#' and category step thresholds `tau` receive weakly-informative `normal(0, 2)`
#' so wide multi-category rubrics (e.g. `K = 6` spanning +/-3 logits) are not
#' overly compressed. Inspect induced response distributions before fitting via
#' [pp_bmfrm()].
#'
#' @param facets_main vector of main facet names
#' @param facets_bias vector of bias facet names
#' @param family model family: `"rsm"` or `"pcm"`
#' @param hierarchical Names of facets with hierarchical distributions (default
#'   `person_facet`; pass `character(0)` for fixed-scale priors).
#' @param person_facet Name of the unconstrained ability facet (default `"person"`).
#'
#' @return list of `bmfrm_prior` objects
#' @export
get_default_priors <- function(facets_main = c("person", "item", "rater"),
                               facets_bias = NULL,
                               family = "rsm",
                               hierarchical = person_facet,
                               person_facet = "person") {
  hierarchical <- validate_hierarchical(hierarchical, facets_main, facets_bias, family)
  priors <- list()

  for (facet in facets_main) {
    if (identical(facet, person_facet)) {
      if (person_facet %in% hierarchical) {
        priors <- c(priors, list(
          prior("normal(0, 2.5)", "mean", facet = person_facet, person_facet = person_facet),
          prior("lognormal(0, 0.5)", "sd", facet = person_facet, person_facet = person_facet)
        ))
      } else {
        priors[[length(priors) + 1L]] <- prior("normal(0, 2)", class = "theta")
      }
    } else if (facet %in% hierarchical) {
      priors[[length(priors) + 1L]] <- prior("normal(0, 1)", class = "sd", facet = facet, person_facet = person_facet)
    } else {
      priors[[length(priors) + 1L]] <- prior("normal(0, 1)", class = "facet", facet = facet, person_facet = person_facet)
    }
  }

  if (!is.null(facets_bias)) {
    for (bias_facet in facets_bias) {
      priors[[length(priors) + 1L]] <- prior(
        "normal(0, 0.5)",
        class = if (bias_facet %in% hierarchical) "sd" else "bias",
        facet = bias_facet,
        person_facet = person_facet
      )
    }
  }

  priors[[length(priors) + 1L]] <- prior("normal(0, 2)", class = "tau", person_facet = person_facet)
  if ("tau" %in% hierarchical) {
    priors[[length(priors) + 1L]] <- prior("normal(0, 0.5)", class = "sd", facet = "tau", person_facet = person_facet)
  }

  priors
}

#' Validate and complete a list of priors against model facets
#' @keywords internal
validate_priors <- function(priors, facets_main, facets_bias = NULL, family = "rsm",
                            hierarchical = person_facet, person_facet = "person") {
  if (missing(hierarchical) && !is.null(priors) && is.list(priors) && length(priors) > 0L &&
      all(vapply(priors, inherits, logical(1), "bmfrm_prior"))) {
    has_theta <- any(vapply(priors, function(p) identical(p$class, "theta"), logical(1)))
    has_person_pop <- any(vapply(priors, function(p) {
      p$class %in% c("mean", "sd") && (is.null(p$facet) || p$facet %in% c("person", person_facet))
    }, logical(1)))
    if (has_theta && !has_person_pop) {
      hierarchical <- setdiff(hierarchical, person_facet)
    }
  }
  defaults <- get_default_priors(
    facets_main = facets_main,
    facets_bias = facets_bias,
    family = family,
    hierarchical = hierarchical,
    person_facet = person_facet
  )
  if (is.null(priors)) return(defaults)
  if (!is.list(priors) || !length(priors)) stop("priors must be a list of prior() objects")
  if (!all(vapply(priors, inherits, logical(1), "bmfrm_prior"))) stop("All elements of priors must be prior() objects")

  priors <- lapply(priors, function(p) {
    fc <- if (is.null(p$facet)) NULL else match_interaction_names(p$facet, facets_bias)
    if (p$class %in% c("mean", "sd") && identical(fc, "person") && !(fc %in% facets_main) && (person_facet %in% facets_main)) {
      fc <- person_facet
    }
    if (identical(p$class, "mean") && !identical(fc, person_facet)) {
      stop("Centered facets, steps, and interactions have no population mean; specify class = 'sd'")
    }
    prior(
      spec = p$spec,
      class = p$class,
      facet = fc,
      scale = p$scale %||% "raw",
      local_shrinkage = p$local_shrinkage %||% "none",
      slab_df = p$slab_df %||% 4,
      slab_scale = p$slab_scale %||% 2,
      person_facet = person_facet
    )
  })

  keys <- vapply(priors, prior_target_key, character(1), person_facet = person_facet, facets_bias = facets_bias)
  req_keys <- vapply(defaults, prior_target_key, character(1), person_facet = person_facet, facets_bias = facets_bias)
  if (anyDuplicated(keys)) stop("Duplicate prior for the same parameter class/facet")
  extra <- setdiff(keys, req_keys)
  if (length(extra)) {
    stop("Prior classes/facets not used by this model: ", paste(extra, collapse = ", "))
  }
  missing_keys <- setdiff(req_keys, keys)
  if (length(missing_keys)) {
    stop("Missing priors for classes/facets: ", paste(missing_keys, collapse = ", "))
  }
  priors
}
