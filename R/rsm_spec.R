# A facet-indexed extension point: no separate argument for every future facet.
# Other hierarchical facets need their own identified parameterization first.
validate_hierarchical <- function(hierarchical, facets_main, facets_bias, family) {
  if (!is.character(hierarchical) || anyNA(hierarchical) || any(!nzchar(hierarchical)) || anyDuplicated(hierarchical)) {
    stop("hierarchical must be a character vector of unique facet names")
  }
  if (any(!hierarchical %in% facets_main)) stop("hierarchical contains a facet absent from the formula")
  if (length(hierarchical) && (family != "rsm" || length(facets_bias))) stop("Hierarchical priors support additive RSM only")
  if (any(hierarchical != "person")) stop("Hierarchical priors currently support only person; other facets are not implemented")
  unname(hierarchical)
}

# Internal normalization shared by generation, Stan data, and Stan code.
# Public configuration lives in hierarchical and priors; this list is not an input API.
resolve_bmfrm_spec <- function(spec, priors = NULL) {
  spec$hierarchical <- validate_hierarchical(spec$hierarchical %||% character(), spec$facets_main, spec$facets_bias, spec$family)
  priors <- validate_priors(priors %||% spec$priors, spec$facets_main, spec$facets_bias, spec$family, spec$hierarchical)
  spec$priors <- priors
  spec$rsm_model <- NULL
  marginal <- any(vapply(priors, function(p) p$scale == "marginal", logical(1)))
  if (!length(spec$hierarchical) && !marginal) return(spec)
  if (spec$family != "rsm" || length(spec$facets_bias) || !"person" %in% spec$facets_main) {
    stop("Hierarchical priors and marginal scaling support generated additive RSM only")
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
  spec$rsm_model <- list(schema = "bayesmfrm.gaussian_rsm.v2", hierarchical = spec$hierarchical,
    ability_mean = mu, ability_sd = sd,
    facet_sd = vapply(blocks, function(p) normal(p)[2L], numeric(1)),
    facet_scale = vapply(blocks, `[[`, character(1), "scale"),
    step_sd = normal(tau)[2L], step_scale = tau$scale)
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

rsm_model_data <- function(model, labels, K) {
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
  out
}

sample_rsm_model <- function(model, labels, K) {
  mu <- model$ability_mean
  sd <- model$ability_sd
  mu <- if (mu$estimated) stats::rnorm(1L, mu$location, mu$scale) else mu$location
  sd <- if (sd$estimated) stats::rlnorm(1L, sd$location, sd$scale) else sd$location
  out <- list(theta = stats::rnorm(length(labels$person), mu, sd))
  for (f in names(model$facet_sd)) {
    n <- length(labels[[f]])
    out[[f]] <- as.vector(model$facet_sd[f] * rsm_prior_basis(n, model$facet_scale[[f]]) %*% stats::rnorm(n - 1L))
  }
  out$tau <- as.vector(model$step_sd * rsm_prior_basis(K - 1L, model$step_scale) %*% stats::rnorm(K - 2L))
  if ("person" %in% model$hierarchical) {
    out$mu_theta <- mu
    out$sigma_theta <- sd
  }
  out
}

rsm_diagnostic_variables <- function(spec, data_info) {
  model <- spec$rsm_model
  if (is.null(model)) return(NULL)
  facets <- names(model$facet_sd)[lengths(data_info$facet_labels[names(model$facet_sd)]) > 1L]
  variable <- c("theta", facets, if (spec$K > 2L) "tau",
    if (model$ability_mean$estimated) "mu_theta",
    if (model$ability_sd$estimated) "sigma_theta")
  raw <- c("theta_raw", if (length(facets)) paste0(facets, "_raw"), if (spec$K > 2L) "tau_raw",
    if (model$ability_mean$estimated) "mu_theta_raw",
    if (model$ability_sd$estimated) "log_sigma_theta_raw")
  unique(c(raw, variable))
}
