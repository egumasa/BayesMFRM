# A facet-indexed extension point: no separate argument for every future facet.
# Other hierarchical facets need their own identified parameterization first.
validate_hierarchical <- function(hierarchical, facets_main, facets_bias, family) {
  if (!is.character(hierarchical) || anyNA(hierarchical) || any(!nzchar(hierarchical)) || anyDuplicated(hierarchical)) {
    stop("hierarchical must be a character vector of unique facet names")
  }
  hierarchical <- match_interaction_names(hierarchical, facets_bias)
  if (anyDuplicated(hierarchical)) stop("hierarchical must contain unique terms, including reversed interaction names")
  if (any(!hierarchical %in% c(facets_main, facets_bias))) stop("hierarchical contains a facet absent from the formula")
  if (length(hierarchical) && !family %in% c("rsm", "pcm")) stop("Hierarchical priors support RSM/PCM only")
  if (any(!hierarchical %in% c("person", facets_bias))) stop("Hierarchical priors currently support only person and two-way interactions; other main facets are not implemented")
  unname(hierarchical)
}

# Internal normalization shared by generation, Stan data, and Stan code.
# Public configuration lives in hierarchical and priors; this list is not an input API.
resolve_bmfrm_spec <- function(spec, priors = NULL) {
  validate_interaction_terms(spec$facets_main, spec$facets_bias)
  spec$hierarchical <- validate_hierarchical(spec$hierarchical %||% character(), spec$facets_main, spec$facets_bias, spec$family)
  priors <- validate_priors(priors %||% spec$priors, spec$facets_main, spec$facets_bias, spec$family, spec$hierarchical)
  spec$priors <- priors
  spec$gaussian_model <- NULL
  marginal <- any(vapply(priors, function(p) p$scale == "marginal", logical(1)))
  gaussian_pcm <- spec$family == "pcm" && !isTRUE(spec$custom_stan) && all(vapply(priors, function(p) {
    x <- tryCatch(parse_bmfrm_prior_literal(p$spec), error = function(e) NULL)
    !is.null(x) && x$distribution == "normal" && length(x$args) == 2L && x$args[2L] > 0
  }, logical(1)))
  if (!length(spec$hierarchical) && !marginal && !gaussian_pcm) return(spec)
  if (!spec$family %in% c("rsm", "pcm") || !"person" %in% spec$facets_main) {
    stop("Hierarchical priors and marginal scaling support generated RSM/PCM only")
  }
  facets <- setdiff(spec$facets_main, "person")
  reserved <- c("theta", "tau", "mu_theta", "sigma_theta", "N", "K", "X", "log_lik", "y_rep", "mu", "sigma2", "mu_fair")
  if (any(facets %in% reserved) || any(grepl("_raw$|^(B_|prior_sd_|ability_|estimate_)", facets)) ||
      any(!grepl("^[A-Za-z][A-Za-z0-9_]*$", facets))) stop("Facet names must be nonreserved Stan identifiers")
  keys <- vapply(priors, function(p) paste(p$class, p$facet %||% "", sep = ":"), character(1))
  get_prior <- function(cl, facet = "") priors[[match(paste(cl, facet, sep = ":"), keys)]]
  normal <- function(p) {
    x <- parse_bmfrm_prior_literal(p$spec)
    if (x$distribution != "normal" || length(x$args) != 2L || x$args[2L] <= 0) {
      stop("Hierarchical priors and marginal scaling currently require Gaussian facet/step and conditional ability priors")
    }
    x$args
  }
  if ("person" %in% spec$hierarchical) {
    mu <- rsm_population_prior(get_prior("mean", "person")$spec, "normal", FALSE)
    sd <- rsm_population_prior(get_prior("sd", "person")$spec, "lognormal", TRUE)
  } else {
    theta <- normal(get_prior("theta"))
    mu <- rsm_population_prior(theta[1L], "normal", FALSE)
    sd <- rsm_population_prior(theta[2L], "lognormal", TRUE)
  }
  blocks <- lapply(facets, function(f) get_prior(if (f %in% c("item", "rater")) f else paste0("facet_", f)))
  names(blocks) <- facets
  tau <- get_prior("tau")
  biases <- lapply(spec$facets_bias, function(term) {
    hierarchical <- term %in% spec$hierarchical
    p <- get_prior(if (hierarchical) "sd" else "bias", term)
    sd <- if (hierarchical) interaction_sd_prior(p$spec) else list(estimated = FALSE, value = normal(p)[2L])
    list(term = term, name = gsub(":", "", term), facets = strsplit(term, ":", fixed = TRUE)[[1L]],
      hierarchical = hierarchical, sd = sd, scale = p$scale)
  })
  names(biases) <- spec$facets_bias
  spec$gaussian_model <- list(schema = "bayesmfrm.gaussian_ordinal.v4", family = spec$family,
    threshold_facet = spec$threshold_facet, hierarchical = spec$hierarchical,
    ability_mean = mu, ability_sd = sd,
    facet_sd = vapply(blocks, function(p) normal(p)[2L], numeric(1)),
    facet_scale = vapply(blocks, `[[`, character(1), "scale"),
    step_sd = normal(tau)[2L], step_scale = tau$scale, biases = biases)
  spec
}

rsm_population_prior <- function(x, distribution, positive) {
  if (is.numeric(x) && length(x) == 1L && is.finite(x) && (!positive || x > 0)) {
    return(list(estimated = FALSE, location = unname(x), scale = 0))
  }
  if (!is.character(x) || length(x) != 1L || is.na(x)) stop("Invalid population prior or fixed value")
  p <- parse_bmfrm_prior_literal(x)
  if (p$distribution != distribution || length(p$args) != 2L || p$args[2L] <= 0) {
    stop("Population prior requires ", distribution, "(location, positive_scale)")
  }
  list(estimated = TRUE, location = unname(p$args[1L]), scale = unname(p$args[2L]))
}

isometric_rsm_basis <- function(n) {
  if (n == 1L) return(matrix(numeric(), 1L, 0L))
  unname(sum_zero_basis(n) * sqrt(n / (n - 1L)))
}

rsm_prior_basis <- function(n, scale) {
  if (scale == "marginal") isometric_rsm_basis(n) else unname(sum_zero_basis(n))
}

gaussian_model_data <- function(model, labels, K) {
  mu <- model$ability_mean
  sd <- model$ability_sd
  out <- list(estimate_mu = as.integer(mu$estimated), estimate_sigma = as.integer(sd$estimated),
    ability_mean_location = mu$location, ability_mean_scale = mu$scale,
    ability_sd_location = sd$location, ability_sd_scale = sd$scale,
    B_tau = rsm_prior_basis(K - 1L, model$step_scale), prior_sd_tau = model$step_sd)
  for (f in names(model$facet_sd)) {
    out[[paste0("B_", f)]] <- rsm_prior_basis(length(labels[[f]]), model$facet_scale[[f]])
    out[[paste0("prior_sd_", f)]] <- unname(model$facet_sd[f])
  }
  for (b in model$biases) {
    out[[paste0("B_", b$name, "_1")]] <- rsm_prior_basis(length(labels[[b$facets[1L]]]), b$scale)
    out[[paste0("B_", b$name, "_2")]] <- rsm_prior_basis(length(labels[[b$facets[2L]]]), b$scale)
    out[[paste0("estimate_sd_", b$name)]] <- as.integer(b$sd$estimated)
    out[[paste0("fixed_sd_", b$name)]] <- if (b$sd$estimated) 1 else b$sd$value
  }
  out
}

sample_gaussian_model <- function(model, labels, K) {
  mu <- model$ability_mean
  sd <- model$ability_sd
  mu <- if (mu$estimated) stats::rnorm(1L, mu$location, mu$scale) else mu$location
  sd <- if (sd$estimated) stats::rlnorm(1L, sd$location, sd$scale) else sd$location
  out <- list(theta = stats::rnorm(length(labels$person), mu, sd))
  for (f in names(model$facet_sd)) {
    n <- length(labels[[f]])
    out[[f]] <- as.vector(model$facet_sd[f] * rsm_prior_basis(n, model$facet_scale[[f]]) %*% stats::rnorm(n - 1L))
  }
  B <- rsm_prior_basis(K - 1L, model$step_scale)
  if (model$family == "pcm") {
    J <- length(labels[[model$threshold_facet]])
    z <- matrix(stats::rnorm(J * (K - 2L)), J, K - 2L)
    out$tau <- unname(model$step_sd * z %*% t(B))
  } else out$tau <- as.vector(model$step_sd * B %*% stats::rnorm(K - 2L))
  if ("person" %in% model$hierarchical) {
    out$mu_theta <- mu
    out$sigma_theta <- sd
  }
  for (b in model$biases) {
    n <- lengths(labels[b$facets])
    sd <- sample_interaction_sd(b$sd)
    z <- matrix(stats::rnorm(prod(n - 1L)), n[1L] - 1L, n[2L] - 1L)
    out[[b$name]] <- unname(sd * rsm_prior_basis(n[1L], b$scale) %*% z %*% t(rsm_prior_basis(n[2L], b$scale)))
    if (b$hierarchical) out[[paste0("sigma_", b$name)]] <- sd
  }
  out
}

model_diagnostic_variables <- function(spec, data_info) {
  model <- spec$gaussian_model
  if (is.null(model)) return(NULL)
  facets <- names(model$facet_sd)[lengths(data_info$facet_labels[names(model$facet_sd)]) > 1L]
  variable <- c("theta", facets, if (spec$K > 2L) "tau",
    if (model$ability_mean$estimated) "mu_theta",
    if (model$ability_sd$estimated) "sigma_theta")
  raw <- c("theta_raw", if (length(facets)) paste0(facets, "_raw"), if (spec$K > 2L) "tau_raw",
    if (model$ability_mean$estimated) "mu_theta_raw",
    if (model$ability_sd$estimated) "log_sigma_theta_raw")
  for (b in model$biases) {
    if (all(lengths(data_info$facet_labels[b$facets]) > 1L)) {
      variable <- c(variable, b$name); raw <- c(raw, paste0(b$name, "_raw"))
    }
    if (b$sd$estimated) {
      variable <- c(variable, paste0("sigma_", b$name)); raw <- c(raw, paste0("sigma_", b$name, "_raw"))
    }
  }
  unique(c(raw, variable))
}

interaction_key <- function(x) {
  if (!length(x)) return(character())
  vapply(strsplit(x, ":", fixed = TRUE), function(parts) paste(sort(trimws(parts)), collapse = ":"), character(1))
}

match_interaction_names <- function(x, terms) {
  at <- match(interaction_key(x), interaction_key(terms))
  x[!is.na(at)] <- terms[at[!is.na(at)]]
  x
}

validate_interaction_terms <- function(facets, terms) {
  if (!length(terms)) return(invisible(TRUE))
  parts <- strsplit(terms, ":", fixed = TRUE)
  if (any(lengths(parts) != 2L) || any(vapply(parts, anyDuplicated, integer(1)) > 0L)) stop("Only two-way interactions between distinct main facets are supported")
  if (any(!unlist(parts) %in% facets)) stop("Interaction components must be present as main facets")
  if (anyDuplicated(interaction_key(terms))) stop("Duplicate interaction term")
  names <- gsub(":", "", terms)
  generated <- c(facets, names, paste0("sigma_", names))
  if (anyDuplicated(generated) || any(!grepl("^[A-Za-z][A-Za-z0-9_]*$", generated)) ||
      any(names %in% c("theta", "tau", "mu_theta", "sigma_theta", "mu", "mu_fair", "sigma2", "log_lik", "y_rep"))) {
    stop("Facet and interaction names collide or are invalid Stan identifiers")
  }
  invisible(TRUE)
}

# The positive-support prior is shared by simulation and generated Stan. A normal,
# Student-t, or Cauchy literal is restricted to zero location, then folded/truncated.
interaction_sd_prior <- function(x) {
  if (is.numeric(x) && length(x) == 1L && is.finite(x) && x > 0) return(list(estimated = FALSE, value = unname(x)))
  if (!is.character(x) || length(x) != 1L || is.na(x)) stop("Interaction SD must be positive or a supported positive-scale prior")
  p <- parse_bmfrm_prior_literal(x)
  a <- p$args; dist <- p$distribution
  valid <- (dist %in% c("normal", "cauchy") && length(a) == 2L && a[1L] == 0 && a[2L] > 0) ||
    (dist == "student_t" && length(a) == 3L && a[1L] > 0 && a[2L] == 0 && a[3L] > 0) ||
    (dist == "lognormal" && length(a) == 2L && a[2L] > 0)
  if (!valid) stop("Interaction SD prior requires normal(0,s), student_t(df,0,s), cauchy(0,s), or lognormal(m,s), with positive scales")
  list(estimated = TRUE, spec = x, distribution = dist, args = a)
}

sample_interaction_sd <- function(p) {
  if (!p$estimated) return(p$value)
  a <- p$args
  switch(p$distribution, normal = abs(stats::rnorm(1L, 0, a[2L])),
    student_t = abs(a[3L] * stats::rt(1L, a[1L])), cauchy = abs(stats::rcauchy(1L, 0, a[2L])),
    lognormal = stats::rlnorm(1L, a[1L], a[2L]))
}
