#' Create a prior specification for bmfrm models
#'
#' @param spec Prior distribution string. For a population mean or SD, a numeric
#'   value fixes that parameter (SD must be positive).
#' @param class Parameter class: theta, item, rater, bias, tau, facet_<name>,
#'   mean, or sd. The last two specify a hierarchical population or scale.
#' @param facet Target facet for mean/sd, or interaction name for bias.
#'   For a hierarchical main facet (`"rater"`, `"item"`, etc.), hierarchical PCM
#'   steps (`"tau"`), or a hierarchical two-way interaction (`"item:rater"`),
#'   use `class = "sd"` and its facet/term name. Supported positive-scale priors
#'   are `normal(0,s)`, `student_t(df,0,s)`, and `cauchy(0,s)` (restricted to
#'   positive values), or `lognormal(m,s)`. A positive numeric value fixes the
#'   SD. Population means (`class = "mean"`) apply only to `"person"`.
#' @param scale Meaning of the SD for a constrained Gaussian facet, step, or
#'   interaction prior: `"raw"` retains the historical pre-centering SD;
#'   `"marginal"` specifies the constrained level SD. Marginal scaling requires
#'   a `normal(0, sd)` prior (or a hierarchical `class = "sd"` prior). Singleton
#'   blocks remain identically zero.
#' @param local_shrinkage Optional local-global shrinkage structure for
#'   hierarchical two-way interactions (`class = "sd"` with an interaction
#'   `facet`): `"none"` (isotropic Gaussian conditional prior, default),
#'   `"student_t"` (unit-variance Student-t scale mixture with df = 3,
#'   `lambda ~ inv_gamma(1.5, 0.5)`, so conditional cell variance given `sigma`
#'   is exactly `sigma^2` under `scale = "marginal"` and
#'   `sigma^2 * (1 - 1/J1) * (1 - 1/J2)` under `scale = "raw"`), or
#'   `"horseshoe"` (regularized horseshoe in latent cell space followed by
#'   double centering; because its global scale parameter is a
#'   distribution-specific shrinkage scale rather than a marginal cell SD,
#'   `"horseshoe"` requires `scale = "raw"`). Note that local shrinkage acts on
#'   unconstrained latent cells `U` prior to double centering `B = P1 U P2`:
#'   double centering preserves label-permutation invariance and exact zero
#'   row/column sums, but redistributes a single sparse latent deviation across
#'   other cells of `B`. Thus latent-cell sparsity is distinct from exact zeros
#'   in the reported relative bias matrix `B`; use a practical threshold
#'   `delta > 0` in `summarise_bias()` to separate large deviations from small
#'   centering spillover.
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
#' # Bias and hierarchical scale priors
#' prior("normal(0, 0.5)", class = "bias", facet = "rater:item")
#' prior("normal(0, 0.5)", class = "sd", facet = "item:rater", scale = "marginal")
#' prior("normal(0, 0.5)", class = "sd", facet = "item:rater", local_shrinkage = "horseshoe")
#' prior("normal(0, 0.5)", class = "sd", facet = "item:rater", scale = "marginal", local_shrinkage = "student_t")
#' prior("normal(0, 1)", class = "sd", facet = "rater", scale = "marginal")
prior <- function(spec, class, facet = NULL, scale = c("raw", "marginal"),
                  local_shrinkage = c("none", "horseshoe", "student_t")) {
  scale <- match.arg(scale)
  local_shrinkage <- match.arg(local_shrinkage)

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

  population <- class %in% c("mean", "sd")
  if (population && (!is.character(facet) || length(facet) != 1L || is.na(facet) || !nzchar(facet))) {
    stop("Population priors require a target facet")
  }
  interaction_sd <- population && grepl(":", facet, fixed = TRUE)
  person_pop <- population && identical(facet, "person")
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
      # person population mean (normal) or sd (lognormal or fixed positive)
      rsm_population_prior(spec, if (class == "mean") "normal" else "lognormal", class == "sd")
    }
  } else if (!is.character(spec) || length(spec) != 1L || is.na(spec) ||
             !grepl("^[a-zA-Z_]+\\s*\\(.*\\)$", spec)) {
    stop("Invalid prior specification '", spec, "'. Must be in format 'distribution(parameters)'")
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
    spec  = spec,    # e.g. "normal(0,2)"
    class = class,   # e.g. "theta", "item", "rater", "bias", "tau", "mean", "sd"
    facet = facet,   # e.g. "rater:item" or "person" or "rater" or "tau"
    scale = scale
  )
  if (local_shrinkage != "none") {
    obj$local_shrinkage <- local_shrinkage
  }
  structure(obj, class = "bmfrm_prior")
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
  if (!is.null(x$local_shrinkage) && x$local_shrinkage != "none") {
    desc <- if (x$local_shrinkage == "student_t") {
      "student_t (df = 3, unit-variance normalized)"
    } else {
      "horseshoe (global shrinkage scale; not a marginal cell SD)"
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
#' @param facets_main vector of main facet names
#' @param facets_bias vector of bias facet names
#' @param family model family: "rsm" or "pcm"
#' @param hierarchical Names of facets with hierarchical distributions. Supports
#'   any main facet (e.g. `"person"`, `"rater"`, `"item"`), two-way interaction
#'   terms (e.g. `"item:rater"`), and `"tau"` for hierarchical PCM partial
#'   pooling toward common steps. `character(0)` retains fixed-scale priors.
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
        priors <- c(priors, list(
          prior("normal(0, 2.5)", "mean", facet = "person"),
          prior("lognormal(0, 0.5)", "sd", facet = "person")
        ))
      } else {
        priors[[length(priors) + 1L]] <- prior("normal(0, 2)", class = "theta")
      }
    } else if (facet %in% hierarchical) {
      priors[[length(priors) + 1L]] <- prior("normal(0, 1)", class = "sd", facet = facet)
    } else if (facet == "item") {
      priors[[length(priors) + 1L]] <- prior("normal(0, 1)", class = "item")
    } else if (facet == "rater") {
      priors[[length(priors) + 1L]] <- prior("normal(0, 1)", class = "rater")
    } else {
      class_name <- paste0("facet_", facet)
      priors[[length(priors) + 1L]] <- prior("normal(0, 1)", class = class_name)
    }
  }

  # Default priors for bias facets
  if (!is.null(facets_bias)) {
    for (bias_facet in facets_bias) {
      priors[[length(priors) + 1L]] <- prior(
        "normal(0, 0.5)",
        class = if (bias_facet %in% hierarchical) "sd" else "bias",
        facet = bias_facet
      )
    }
  }

  # Default prior(s) for thresholds
  priors[[length(priors) + 1L]] <- prior("normal(0, 1)", class = "tau")
  if ("tau" %in% hierarchical) {
    priors[[length(priors) + 1L]] <- prior("normal(0, 0.5)", class = "sd", facet = "tau")
  }

  priors
}

#' Validate a list of priors
#'
#' @param priors list of bmfrm_prior objects or NULL
#' @param facets_main vector of main facet names
#' @param facets_bias vector of bias facet names
#' @param family model family: "rsm" or "pcm"
#' @param hierarchical Names of facets with hierarchical distributions.
#'
#' @return validated list of priors
#' @keywords internal
validate_priors <- function(priors, facets_main, facets_bias = NULL, family = "rsm", hierarchical = character()) {
  defaults <- get_default_priors(facets_main, facets_bias, family, hierarchical)
  if (is.null(priors)) priors <- defaults
  if (!is.list(priors) || !length(priors)) stop("priors must be a list of prior() objects")
  if (!all(vapply(priors, inherits, logical(1), "bmfrm_prior"))) stop("All elements of priors must be prior() objects")
  # Revalidate edited objects and normalize saved priors without scale/local_shrinkage fields.
  priors <- lapply(priors, function(p) {
    prior(
      p$spec,
      p$class,
      if (is.null(p$facet)) NULL else match_interaction_names(p$facet, facets_bias),
      p$scale %||% "raw",
      p$local_shrinkage %||% "none"
    )
  })
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
