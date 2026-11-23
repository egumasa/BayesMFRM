#' Simulate ratings from generated additive RSM or PCM
#'
#' With parameters = NULL, generate parameters and ratings from the same priors
#' and constraints used by Stan. Independent raw-scale priors support normal and
#' student_t numeric literals. Hierarchical ability and marginal scaling use
#' Gaussian zero-sum blocks and the declared population priors.
#' Fixed parameters must already satisfy the constraints. Unnamed vectors follow
#' facet_labels order; named vectors are matched to the original facet labels.
#'
#' @param formula Model formula with person and additive facet columns.
#' @param design Planned rating rows. Existing responses are ignored.
#' @param K Number of categories.
#' @param priors List of prior specifications; NULL uses package defaults.
#' @param parameters NULL or a list with theta, each other facet, and tau.
#'   With hierarchical ability, also supply mu_theta and sigma_theta.
#' @param nsim Number of simulated datasets.
#' @param seed Optional seed, restored on exit.
#' @param hierarchical Names of facets with hierarchical distributions; currently
#'   only "person" is implemented. Its population uses prior classes mean and sd
#'   with facet = "person". The default character(0) retains independent priors.
#' @param family "rsm" for common steps or "pcm" for steps within threshold_facet.
#' @param threshold_facet Required for PCM. Its step deviations sum to zero within
#'   each level. Fixed tau must be a levels-by-(K-1) matrix; optional row names
#'   match facet labels and column names match step numbers 1,...,K-1.
#' @return A bmfrm_simulation with y_rep (simulations by ratings), generating
#'   parameters, facet_labels, specification, and data for nsim = 1.
#' @export
simulate_bmfrm <- function(formula, design, K, priors = NULL, parameters = NULL,
                           nsim = 1L, seed = NULL, hierarchical = character(),
                           family = c("rsm", "pcm"), threshold_facet = NULL) {
  family <- match.arg(family)
  validate_category_count(K)
  if (length(nsim) != 1L || !is.finite(nsim) || nsim < 1 || nsim != floor(nsim)) stop("nsim must be a positive integer")
  if (!is.data.frame(design) || !nrow(design)) stop("design must contain rating rows")
  response <- all.vars(formula)[1L]
  fake <- design
  fake[[response]] <- 1L
  spec <- resolve_bmfrm_spec(parse_bmfrm_formula(formula, fake, K = K, hierarchical = hierarchical,
    family = family, threshold_facet = threshold_facet), priors)
  model <- spec$gaussian_model
  priors <- spec$priors
  if (length(spec$facets_bias) || !"person" %in% spec$facets_main ||
      any(!spec$facets_main %in% names(design))) stop("Simulation currently supports additive RSM/PCM models")
  fs <- lapply(design[spec$facets_main], function(x) if (is.factor(x)) x else factor(x))
  if (any(vapply(fs, anyNA, logical(1)))) stop("Missing values in design facets")
  labels <- lapply(fs, levels)
  ids <- lapply(fs, as.integer)
  param_names <- c(ifelse(spec$facets_main == "person", "theta", spec$facets_main), "tau")
  sizes <- c(lengths(labels), K - 1L)
  names(sizes) <- param_names
  J_step <- if (family == "pcm") length(labels[[spec$threshold_facet]]) else 1L
  sizes["tau"] <- J_step * (K - 1L)
  if ("person" %in% hierarchical) {
    param_names <- c(param_names, "mu_theta", "sigma_theta")
    sizes <- c(sizes, mu_theta = 1L, sigma_theta = 1L)
  }
  if (!is.null(parameters)) {
    if (!is.list(parameters) || anyDuplicated(names(parameters)) || !setequal(names(parameters), param_names)) stop("parameters must contain exactly: ", paste(param_names, collapse = ", "))
    parameters <- parameters[param_names]
    for (j in seq_along(parameters)) {
      x <- parameters[[j]]
      if (!is.numeric(x) || length(x) != sizes[j] || any(!is.finite(x))) stop("Invalid parameter vector: ", param_names[j])
      if (param_names[j] == "tau" && family == "pcm") {
        if (!is.matrix(x) || !identical(dim(x), c(J_step, as.integer(K - 1L)))) stop("PCM tau must have threshold-facet levels as rows and K-1 steps as columns")
        if (!is.null(rownames(x))) {
          if (anyDuplicated(rownames(x)) || !setequal(rownames(x), labels[[spec$threshold_facet]])) stop("PCM tau row labels do not match threshold_facet")
          x <- x[labels[[spec$threshold_facet]], , drop = FALSE]
        }
        if (!is.null(colnames(x))) {
          if (anyDuplicated(colnames(x)) || !setequal(colnames(x), as.character(seq_len(K - 1L)))) stop("PCM tau column names must be step numbers")
          x <- x[, as.character(seq_len(K - 1L)), drop = FALSE]
        }
        if (any(abs(rowSums(x)) > 1e-8)) stop("PCM tau must sum to zero within every threshold-facet level")
        parameters[[j]] <- unname(x)
        next
      }
      if (!is.null(dim(x))) stop("Non-PCM parameter must be a vector: ", param_names[j])
      if (j <= length(labels) && !is.null(names(x))) {
        if (anyDuplicated(names(x)) || !setequal(names(x), labels[[j]])) stop("Parameter labels do not match: ", param_names[j])
        x <- x[labels[[j]]]
      }
      if (!param_names[j] %in% c("theta", "mu_theta", "sigma_theta") && abs(sum(x)) > 1e-8) stop("Fixed parameter must sum to zero: ", param_names[j])
      parameters[[j]] <- unname(x)
    }
    if ("person" %in% hierarchical) {
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
    if (is.null(pars) && !is.null(model)) pars <- sample_gaussian_model(model, labels, K)
    if (is.null(pars)) {
      pars <- lapply(seq_along(param_names), function(j) {
        x <- sample_bmfrm_prior(priors[[match(classes[j], prior_classes)]]$spec, sizes[j])
        if (param_names[j] == "tau" && family == "pcm") {
          x <- matrix(x, J_step, K - 1L)
          x <- x - rowMeans(x)
        } else if (param_names[j] != "theta") x <- x - mean(x)
        x
      })
      names(pars) <- param_names
    }
    eta <- pars$theta[ids$person]
    for (f in setdiff(spec$facets_main, "person")) eta <- eta - pars[[f]][ids[[f]]]
    prob <- if (family == "pcm") pcm_probabilities(eta, pars$tau[ids[[spec$threshold_facet]], , drop = FALSE]) else rsm_probabilities(eta, pars$tau)
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
                 facet_labels = labels, spec = spec, priors = priors, hierarchical = hierarchical,
                 source = if (is.null(parameters)) "prior" else "fixed", seed = seed),
            class = "bmfrm_simulation")
}

#' Generate prior predictive ratings on a planned RSM/PCM design
#' @inheritParams simulate_bmfrm
#' @param ndraws Number of prior predictive datasets.
#' @return A bmfrm_simulation. Inspect its y_rep before fitting observed data.
#' @export
prior_predict_bmfrm <- function(formula, design, K, priors = NULL, ndraws = 200L, seed = NULL,
                               hierarchical = character(), family = c("rsm", "pcm"), threshold_facet = NULL) {
  simulate_bmfrm(formula, design, K, priors, nsim = ndraws, seed = seed, hierarchical = hierarchical,
    family = match.arg(family), threshold_facet = threshold_facet)
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

pcm_probabilities <- function(eta, tau) {
  steps <- tau
  if (ncol(steps) > 1L) for (j in 2:ncol(steps)) steps[, j] <- steps[, j] + steps[, j - 1L]
  logits <- outer(eta, 0:ncol(steps)) - cbind(0, steps)
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
