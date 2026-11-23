#' Specify an isometric Gaussian RSM and its population priors
#'
#' An explicit alternative to the historical raw-centered model. Every non-person
#' facet and the shared step vector use J - 1 orthonormal zero-sum coordinates,
#' scaled so each nonconstant level has the requested marginal prior SD. Singleton
#' facets and the binary step are identically zero. Ability is mu + sigma * z,
#' with independent standard-normal z. The population mean and SD may each be
#' fixed or estimated. This specification supports additive RSM only.
#'
#' @param ability_mean A finite fixed mean or literal normal(mean, sd) expression.
#' @param ability_sd A positive fixed SD or literal lognormal(log_mean, log_sd)
#'   expression. The latter describes the SD on its positive scale.
#' @param facet_sd One positive marginal SD for all non-person facets, or a named
#'   numeric vector specifying every non-person facet exactly.
#' @param step_sd Positive marginal SD of each shared step when K > 2.
#' @return A bmfrm_rsm_model for the model argument of bmfrm(), simulate_bmfrm(),
#'   and prior_predict_bmfrm(). Do not combine it with the historical priors argument.
#' @export
rsm_model <- function(ability_mean = "normal(0, 2.5)",
                      ability_sd = "lognormal(0, 0.5)", facet_sd = 1, step_sd = 1) {
  positive <- function(x) is.numeric(x) && length(x) > 0L && all(is.finite(x)) && all(x > 0)
  if (!positive(facet_sd) || (length(facet_sd) > 1L && is.null(names(facet_sd))) ||
      (!is.null(names(facet_sd)) && (anyNA(names(facet_sd)) || any(!nzchar(names(facet_sd))) || anyDuplicated(names(facet_sd))))) {
    stop("facet_sd must be positive: an unnamed scalar or a uniquely named vector")
  }
  if (!positive(step_sd) || length(step_sd) != 1L) stop("step_sd must be one positive finite SD")
  structure(list(schema = "bayesmfrm.isometric_rsm.v1",
    ability_mean = rsm_population_prior(ability_mean, "normal", FALSE),
    ability_sd = rsm_population_prior(ability_sd, "lognormal", TRUE),
    facet_sd = facet_sd, step_sd = unname(step_sd)), class = "bmfrm_rsm_model")
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

bind_rsm_model <- function(model, spec) {
  if (is.null(model)) return(NULL)
  if (!inherits(model, "bmfrm_rsm_model") || !identical(model$schema, "bayesmfrm.isometric_rsm.v1")) stop("Use rsm_model() to specify model")
  if (spec$family != "rsm" || length(spec$facets_bias) || !"person" %in% spec$facets_main) stop("rsm_model() supports additive RSM only")
  facets <- setdiff(spec$facets_main, "person")
  # Revalidate serializable specifications, including user-edited list fields.
  restore <- function(p, distribution) {
    if (!is.list(p) || !is.logical(p$estimated) || length(p$estimated) != 1L || is.na(p$estimated)) stop("Invalid population specification")
    if (!p$estimated) return(p$location)
    if (!is.numeric(p$location) || length(p$location) != 1L || !is.finite(p$location) ||
        !is.numeric(p$scale) || length(p$scale) != 1L || !is.finite(p$scale) || p$scale <= 0) stop("Invalid population specification")
    paste0(distribution, "(", sprintf("%.17g", p$location), ",", sprintf("%.17g", p$scale), ")")
  }
  facet_sd <- if (!length(facets) && !length(model$facet_sd)) 1 else model$facet_sd
  model <- rsm_model(restore(model$ability_mean, "normal"), restore(model$ability_sd, "lognormal"), facet_sd, model$step_sd)
  reserved <- c("theta", "tau", "mu_theta", "sigma_theta", "N", "K", "X", "log_lik", "y_rep", "mu", "sigma2", "mu_fair")
  if (any(facets %in% reserved) || any(grepl("_raw$|^(B_|prior_sd_|ability_|estimate_)", facets)) ||
      any(!grepl("^[A-Za-z][A-Za-z0-9_]*$", facets))) stop("Facet names must be nonreserved Stan identifiers")
  if (is.null(names(model$facet_sd))) {
    model$facet_sd <- stats::setNames(rep(model$facet_sd, length(facets)), facets)
  } else {
    if (!setequal(names(model$facet_sd), facets)) stop("Named facet_sd must specify every non-person facet exactly")
    model$facet_sd <- model$facet_sd[facets]
  }
  model
}

isometric_rsm_basis <- function(n) {
  if (n == 1L) return(matrix(numeric(), 1L, 0L))
  unname(sum_zero_basis(n) * sqrt(n / (n - 1L)))
}

rsm_model_data <- function(model, labels, K) {
  mu <- model$ability_mean
  sd <- model$ability_sd
  out <- list(estimate_mu = as.integer(mu$estimated), estimate_sigma = as.integer(sd$estimated),
    ability_mean_location = mu$location, ability_mean_scale = mu$scale,
    ability_sd_location = sd$location, ability_sd_scale = sd$scale,
    B_tau = isometric_rsm_basis(K - 1L), prior_sd_tau = model$step_sd)
  for (f in names(model$facet_sd)) {
    out[[paste0("B_", f)]] <- isometric_rsm_basis(length(labels[[f]]))
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
    out[[f]] <- as.vector(model$facet_sd[f] * isometric_rsm_basis(n) %*% stats::rnorm(n - 1L))
  }
  out$tau <- as.vector(model$step_sd * isometric_rsm_basis(K - 1L) %*% stats::rnorm(K - 2L))
  out$mu_theta <- mu
  out$sigma_theta <- sd
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

#' @export
print.bmfrm_rsm_model <- function(x, ...) {
  cat("Additive isometric RSM (explicit Gaussian prior specification)\n")
  cat("Ability mean:", if (x$ability_mean$estimated) "estimated (normal prior)" else x$ability_mean$location, "\n")
  cat("Ability SD:", if (x$ability_sd$estimated) "estimated (lognormal prior)" else x$ability_sd$location, "\n")
  cat("Facet and step SDs refer to constrained marginal levels; singleton blocks are zero.\n")
  invisible(x)
}
