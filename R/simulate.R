#' Simulate ratings from the generated additive RSM
#'
#' With parameters = NULL, sample each raw effect from its specified prior and
#' apply the same centering as Stan. This supports normal and student_t priors
#' with literal numeric arguments. Other Stan prior expressions are rejected.
#' Fixed parameters must already satisfy the constraints. Unnamed vectors follow
#' facet_labels order; named vectors are matched to the original facet labels.
#'
#' @param formula Model formula with person and additive facet columns.
#' @param design Planned rating rows. Existing responses are ignored.
#' @param K Number of categories.
#' @param priors List of prior specifications; NULL uses package defaults.
#' @param parameters NULL or a list with theta, each other facet, and tau.
#'   With an explicit model, also supply mu_theta and sigma_theta.
#' @param nsim Number of simulated datasets.
#' @param seed Optional seed, restored on exit.
#' @param model Optional explicit rsm_model(); NULL retains the historical priors.
#' @return A bmfrm_simulation with y_rep (simulations by ratings), generating
#'   parameters, facet_labels, specification, and data for nsim = 1.
#' @export
simulate_bmfrm <- function(formula, design, K, priors = NULL, parameters = NULL,
                           nsim = 1L, seed = NULL, model = NULL) {
  validate_category_count(K)
  if (length(nsim) != 1L || !is.finite(nsim) || nsim < 1 || nsim != floor(nsim)) stop("nsim must be a positive integer")
  if (!is.data.frame(design) || !nrow(design)) stop("design must contain rating rows")
  response <- all.vars(formula)[1L]
  fake <- design
  fake[[response]] <- 1L
  spec <- parse_bmfrm_formula(formula, fake, K = K)
  model <- bind_rsm_model(model, spec)
  spec$rsm_model <- model
  if (!is.null(model) && !is.null(priors)) stop("Use model or priors, not both")
  if (length(spec$facets_bias) || !"person" %in% spec$facets_main ||
      any(!spec$facets_main %in% names(design))) stop("Simulation currently supports additive RSM models")
  fs <- lapply(design[spec$facets_main], function(x) if (is.factor(x)) x else factor(x))
  if (any(vapply(fs, anyNA, logical(1)))) stop("Missing values in design facets")
  labels <- lapply(fs, levels)
  ids <- lapply(fs, as.integer)
  if (is.null(model)) priors <- validate_priors(priors, spec$facets_main, spec$facets_bias)
  param_names <- c(ifelse(spec$facets_main == "person", "theta", spec$facets_main), "tau")
  sizes <- c(lengths(labels), K - 1L)
  names(sizes) <- param_names
  if (!is.null(model)) {
    param_names <- c(param_names, "mu_theta", "sigma_theta")
    sizes <- c(sizes, mu_theta = 1L, sigma_theta = 1L)
  }
  if (!is.null(parameters)) {
    if (!is.list(parameters) || anyDuplicated(names(parameters)) || !setequal(names(parameters), param_names)) stop("parameters must contain exactly: ", paste(param_names, collapse = ", "))
    parameters <- parameters[param_names]
    for (j in seq_along(parameters)) {
      x <- parameters[[j]]
      if (!is.numeric(x) || length(x) != sizes[j] || any(!is.finite(x))) stop("Invalid parameter vector: ", param_names[j])
      if (j <= length(labels) && !is.null(names(x))) {
        if (anyDuplicated(names(x)) || !setequal(names(x), labels[[j]])) stop("Parameter labels do not match: ", param_names[j])
        x <- x[labels[[j]]]
      }
      if (!param_names[j] %in% c("theta", "mu_theta", "sigma_theta") && abs(sum(x)) > 1e-8) stop("Fixed parameter must sum to zero: ", param_names[j])
      parameters[[j]] <- unname(x)
    }
    if (!is.null(model)) {
      if (parameters$sigma_theta <= 0) stop("sigma_theta must be positive")
      if ((!model$ability_mean$estimated && parameters$mu_theta != model$ability_mean$location) ||
          (!model$ability_sd$estimated && parameters$sigma_theta != model$ability_sd$location)) stop("Fixed population parameters must match model")
    }
  }
  classes <- c(ifelse(spec$facets_main == "person", "theta",
                     ifelse(spec$facets_main %in% c("item", "rater"), spec$facets_main,
                            paste0("facet_", spec$facets_main))), "tau")
  prior_classes <- vapply(priors, function(p) p$class, character(1))
  if (anyDuplicated(prior_classes)) stop("Use one prior for each parameter class")
  sample_one <- function() {
    pars <- parameters
    if (is.null(pars) && !is.null(model)) pars <- sample_rsm_model(model, labels, K)
    if (is.null(pars)) {
      pars <- lapply(seq_along(param_names), function(j) {
        x <- sample_bmfrm_prior(priors[[match(classes[j], prior_classes)]]$spec, sizes[j])
        if (param_names[j] != "theta") x <- x - mean(x)
        x
      })
      names(pars) <- param_names
    }
    eta <- pars$theta[ids$person]
    for (f in setdiff(spec$facets_main, "person")) eta <- eta - pars[[f]][ids[[f]]]
    prob <- rsm_probabilities(eta, pars$tau)
    cum <- t(apply(prob, 1L, cumsum))
    y <- 1L + rowSums(stats::runif(nrow(design)) > cum[, -K, drop = FALSE])
    list(parameters = pars, y = as.integer(y))
  }
  generated <- with_bmfrm_seed(seed, replicate(nsim, sample_one(), simplify = FALSE))
  yrep <- do.call(rbind, lapply(generated, `[[`, "y"))
  dat <- design
  dat[[response]] <- yrep[1L, ]
  structure(list(y_rep = yrep, parameters = lapply(generated, `[[`, "parameters"),
                 data = if (nsim == 1L) dat else NULL, design = design,
                 facet_labels = labels, spec = spec, priors = priors, model = model,
                 source = if (is.null(parameters)) "prior" else "fixed", seed = seed),
            class = "bmfrm_simulation")
}

#' Generate prior predictive ratings on a planned RSM design
#' @inheritParams simulate_bmfrm
#' @param ndraws Number of prior predictive datasets.
#' @return A bmfrm_simulation. Inspect its y_rep before fitting observed data.
#' @export
prior_predict_bmfrm <- function(formula, design, K, priors = NULL, ndraws = 200L, seed = NULL, model = NULL) {
  simulate_bmfrm(formula, design, K, priors, nsim = ndraws, seed = seed, model = model)
}

with_bmfrm_seed <- function(seed, code) {
  if (is.null(seed)) return(force(code))
  existed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (existed) old <- get(".Random.seed", envir = .GlobalEnv)
  on.exit(if (existed) assign(".Random.seed", old, envir = .GlobalEnv) else
    if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv))
  set.seed(seed)
  force(code)
}

rsm_probabilities <- function(eta, tau) {
  logits <- sweep(outer(eta, 0:length(tau)), 2L, c(0, cumsum(tau)), "-")
  z <- exp(logits - apply(logits, 1L, max))
  z / rowSums(z)
}

sample_bmfrm_prior <- function(spec, n) {
  p <- parse_bmfrm_prior_literal(spec)
  distribution <- p$distribution
  args <- p$args
  if (distribution == "normal" && length(args) == 2L && args[2L] > 0) return(stats::rnorm(n, args[1L], args[2L]))
  if (distribution == "student_t" && length(args) == 3L && args[1L] > 0 && args[3L] > 0) return(args[2L] + args[3L] * stats::rt(n, args[1L]))
  stop("Simulation supports normal(mean, sd) and student_t(df, location, scale) with positive scales")
}

parse_bmfrm_prior_literal <- function(spec) {
  expr <- tryCatch(str2lang(spec), error = function(e) NULL)
  if (!is.call(expr) || !is.symbol(expr[[1L]])) stop("Unsupported simulation prior: ", spec)
  distribution <- as.character(expr[[1L]])
  if (any(nzchar(names(as.list(expr)[-1L])))) stop("Prior arguments must be positional numeric literals")
  literal <- function(x) {
    if (is.numeric(x) && length(x) == 1L && is.finite(x)) return(x)
    if (is.call(x) && length(x) == 2L && is.symbol(x[[1L]]) &&
        as.character(x[[1L]]) %in% c("+", "-") && is.numeric(x[[2L]]) &&
        length(x[[2L]]) == 1L && is.finite(x[[2L]])) {
      return(if (as.character(x[[1L]]) == "-") -x[[2L]] else x[[2L]])
    }
    stop("Simulation priors require literal numeric arguments")
  }
  args <- vapply(as.list(expr)[-1L], literal, numeric(1))
  list(distribution = distribution, args = args)
}
