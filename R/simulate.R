#' Simulate ratings from generated RSM or PCM, including two-way interactions
#'
#' With parameters = NULL, generate parameters and ratings from the same priors
#' and constraints used by Stan. Independent raw-scale priors support normal and
#' student_t and double_exponential numeric literals. Hierarchical priors and marginal scaling use
#' Gaussian zero-sum blocks and the declared population priors.
#' Fixed parameters must already satisfy the constraints. Unnamed vectors follow
#' facet_labels order; named vectors are matched to the original facet labels.
#'
#' @param formula Model formula with person, main facets and optional two-way interactions.
#' @param design Planned rating rows. Existing responses are ignored.
#' @param K Number of categories (`>= 2`): scalar integer for RSM or uniform
#'   PCM, or a named integer vector matching `threshold_facet` levels for
#'   heterogeneous-category PCM.
#' @param priors List of prior specifications; NULL uses package defaults.
#' @param parameters NULL or a list with theta, each other facet, and tau.
#'   With hierarchical ability, also supply mu_theta and sigma_theta.
#'   Interactions use a matrix named by concatenating their formula-order facet
#'   names (e.g. itemrater), with zero row and column sums. Dimnames match facet
#'   labels. Hierarchical interactions also require sigma_<interaction_name>.
#' @param nsim Number of simulated datasets.
#' @param seed Optional seed, restored on exit.
#' @param hierarchical Names of facets or terms with hierarchical distributions:
#'   `"person"` (classes `"mean"` and `"sd"` with `facet = "person"`), centered
#'   main facets (e.g., `"rater"`, `"item"` with `class = "sd"`), two-way
#'   interaction terms (e.g., `"item:rater"` with `class = "sd"`), and `"tau"`
#'   for uniform hierarchical PCM step sharing. The default `character(0)`
#'   retains fixed-scale priors.
#' @param family "rsm" for common steps or "pcm" for steps within threshold_facet.
#' @param threshold_facet Required for PCM. Its step deviations sum to zero within
#'   each level. For uniform `K`, fixed `tau` must be a levels-by-(K-1) matrix;
#'   optional row names match facet labels and column names match step numbers
#'   1,...,K-1. For heterogeneous `K`, fixed `tau` must be a named list of numeric
#'   vectors of length `K_j - 1`, each summing to zero (`0` for binary `K_j = 2`).
#' @return A bmfrm_simulation with y_rep (simulations by ratings), generating
#'   parameters, facet_labels, specification, and data for nsim = 1.
#' @export
simulate_bmfrm <- function(formula, design, K, priors = NULL, parameters = NULL,
                           nsim = 1L, seed = NULL, hierarchical = character(),
                           family = c("rsm", "pcm"), threshold_facet = NULL) {
  family <- match.arg(family)
  if (length(nsim) != 1L || !is.finite(nsim) || nsim < 1 || nsim != floor(nsim)) stop("nsim must be a positive integer")
  if (!is.data.frame(design) || !nrow(design)) stop("design must contain rating rows")
  response <- bmfrm_response_name(formula)
  fake <- design
  fake[[response]] <- 1L
  spec <- resolve_bmfrm_spec(parse_bmfrm_formula(formula, fake, K = K, hierarchical = hierarchical,
    family = family, threshold_facet = threshold_facet), priors)
  model <- spec$gaussian_model
  priors <- spec$priors
  hierarchical <- spec$hierarchical
  K_max <- get_spec_K_max(spec)
  is_hetero <- isTRUE(spec$K_heterogeneous)
  K_by_level <- spec$K_by_level
  if (!"person" %in% spec$facets_main ||
      any(!spec$facets_main %in% names(design))) stop("Simulation requires person and named main facets")
  fs <- lapply(design[spec$facets_main], function(x) if (is.factor(x)) x else factor(x))
  if (any(vapply(fs, anyNA, logical(1)))) stop("Missing values in design facets")
  labels <- lapply(fs, levels)
  ids <- lapply(fs, as.integer)
  param_names <- c(ifelse(spec$facets_main == "person", "theta", spec$facets_main), "tau")
  sizes <- c(lengths(labels), K_max - 1L)
  names(sizes) <- param_names
  J_step <- if (family == "pcm") length(labels[[spec$threshold_facet]]) else 1L
  sizes["tau"] <- if (is_hetero) sum(as.integer(K_by_level[labels[[spec$threshold_facet]]]) - 1L) else J_step * (K_max - 1L)
  if ("person" %in% hierarchical) {
    param_names <- c(param_names, "mu_theta", "sigma_theta")
    sizes <- c(sizes, mu_theta = 1L, sigma_theta = 1L)
  }
  hier_main <- setdiff(intersect(hierarchical, spec$facets_main), "person")
  for (hf in hier_main) {
    nm <- paste0("sigma_", hf)
    param_names <- c(param_names, nm)
    sizes <- c(sizes, stats::setNames(1L, nm))
  }
  if ("tau" %in% hierarchical) {
    param_names <- c(param_names, "tau_common", "sigma_tau")
    sizes <- c(sizes, tau_common = K_max - 1L, sigma_tau = 1L)
  }
  bias_names <- gsub(":", "", spec$facets_bias)
  bias_dims <- lapply(strsplit(spec$facets_bias, ":", fixed = TRUE), function(f) lengths(labels[f]))
  names(bias_dims) <- bias_names
  for (j in seq_along(bias_names)) {
    param_names <- c(param_names, bias_names[j]); sizes <- c(sizes, stats::setNames(prod(bias_dims[[j]]), bias_names[j]))
    if (spec$facets_bias[j] %in% hierarchical) {
      nm <- paste0("sigma_", bias_names[j]); param_names <- c(param_names, nm); sizes <- c(sizes, stats::setNames(1L, nm))
    }
  }
  if (!is.null(parameters)) {
    if (!is.list(parameters) || anyDuplicated(names(parameters)) || !setequal(names(parameters), param_names)) stop("parameters must contain exactly: ", paste(param_names, collapse = ", "))
    parameters <- parameters[param_names]
    for (j in seq_along(parameters)) {
      x <- parameters[[j]]
      if (param_names[j] == "tau" && is_hetero) {
        tf_lv <- labels[[spec$threshold_facet]]
        if (!is.list(x) || length(x) != length(tf_lv)) {
          stop("Heterogeneous PCM tau must be a list of step vectors, one per level of ", spec$threshold_facet)
        }
        if (!is.null(names(x))) {
          if (anyDuplicated(names(x)) || !setequal(names(x), tf_lv)) {
            stop("Heterogeneous PCM tau list names do not match threshold_facet levels")
          }
          x <- x[tf_lv]
        } else {
          names(x) <- tf_lv
        }
        for (lv in tf_lv) {
          Kj <- as.integer(K_by_level[[lv]])
          tv <- x[[lv]]
          if (!is.numeric(tv) || length(tv) != (Kj - 1L) || any(!is.finite(tv))) {
            stop("PCM tau for level '", lv, "' must be a finite numeric vector of length ", Kj - 1L)
          }
          if (abs(sum(tv)) > 1e-8) {
            stop("PCM tau for level '", lv, "' must sum to zero")
          }
          x[[lv]] <- as.numeric(unname(tv))
        }
        parameters[[j]] <- x
        next
      }
      if (!is.numeric(x) || length(x) != sizes[j] || any(!is.finite(x))) stop("Invalid parameter vector: ", param_names[j])
      if (param_names[j] %in% bias_names) {
        f <- strsplit(spec$facets_bias[match(param_names[j], bias_names)], ":", fixed = TRUE)[[1L]]
        if (!is.matrix(x) || !identical(dim(x), unname(bias_dims[[param_names[j]]]))) stop("Interaction must be a matrix with one dimension per facet: ", param_names[j])
        for (axis in 1:2) if (!is.null(dimnames(x)[[axis]])) {
          nm <- dimnames(x)[[axis]]
          if (anyDuplicated(nm) || !setequal(nm, labels[[f[axis]]])) stop("Interaction matrix labels do not match facets")
          x <- if (axis == 1L) x[labels[[f[axis]]], , drop = FALSE] else x[, labels[[f[axis]]], drop = FALSE]
        }
        if (any(abs(c(rowSums(x), colSums(x))) > 1e-8)) stop("Interaction rows and columns must each sum to zero")
        parameters[[j]] <- unname(x); next
      }
      if (param_names[j] %in% paste0("sigma_", bias_names)) {
        b <- model$biases[[match(sub("^sigma_", "", param_names[j]), bias_names)]]
        if (x <= 0 || (!b$sd$estimated && x != b$sd$value)) stop("Interaction SD must be positive and match any fixed prior value")
        parameters[[j]] <- unname(x); next
      }
      if (param_names[j] %in% paste0("sigma_", hier_main)) {
        fs <- model$facet_specs[[sub("^sigma_", "", param_names[j])]]
        if (x <= 0 || (!fs$sd$estimated && x != fs$sd$value)) stop("Hierarchical facet SD must be positive and match any fixed prior value")
        parameters[[j]] <- unname(x); next
      }
      if (param_names[j] == "sigma_tau") {
        if (x <= 0 || (!model$step$dev_sd$estimated && x != model$step$dev_sd$value)) stop("Hierarchical step SD must be positive and match any fixed prior value")
        parameters[[j]] <- unname(x); next
      }
      if (param_names[j] == "tau_common") {
        if (!is.null(dim(x))) stop("tau_common must be a vector")
        if (abs(sum(x)) > 1e-8) stop("Fixed tau_common must sum to zero")
        parameters[[j]] <- unname(x); next
      }
      if (param_names[j] == "tau" && family == "pcm") {
        if (!is.matrix(x) || !identical(dim(x), c(J_step, as.integer(K_max - 1L)))) stop("PCM tau must have threshold-facet levels as rows and K-1 steps as columns")
        if (!is.null(rownames(x))) {
          if (anyDuplicated(rownames(x)) || !setequal(rownames(x), labels[[spec$threshold_facet]])) stop("PCM tau row labels do not match threshold_facet")
          x <- x[labels[[spec$threshold_facet]], , drop = FALSE]
        }
        if (!is.null(colnames(x))) {
          if (anyDuplicated(colnames(x)) || !setequal(colnames(x), as.character(seq_len(K_max - 1L)))) stop("PCM tau column names must be step numbers")
          x <- x[, as.character(seq_len(K_max - 1L)), drop = FALSE]
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
  sample_one <- function() {
    pars <- parameters
    if (is.null(pars) && !is.null(model)) {
      pars <- sample_gaussian_model(model, labels, K_max, K_by_level = K_by_level)
    }
    if (is.null(pars)) {
      pars <- lapply(seq_along(classes), function(j) {
        if (param_names[j] == "tau" && is_hetero) {
          tf_lv <- labels[[spec$threshold_facet]]
          p_spec <- priors[[match(classes[j], prior_classes)]]$spec
          tl <- lapply(tf_lv, function(lv) {
            Kj <- as.integer(K_by_level[[lv]])
            if (Kj == 2L) return(0)
            v <- sample_bmfrm_prior(p_spec, Kj - 1L)
            v - mean(v)
          })
          names(tl) <- tf_lv
          return(tl)
        }
        x <- sample_bmfrm_prior(priors[[match(classes[j], prior_classes)]]$spec, sizes[j])
        if (param_names[j] == "tau" && family == "pcm") {
          x <- matrix(x, J_step, K_max - 1L)
          x <- x - rowMeans(x)
        } else if (param_names[j] != "theta") x <- x - mean(x)
        x
      })
      names(pars) <- param_names[seq_along(classes)]
      for (j in seq_along(bias_names)) {
        p <- priors[[which(prior_classes == "bias" & vapply(priors, function(p) identical(p$facet, spec$facets_bias[j]), logical(1)))]]
        n <- bias_dims[[j]]
        x <- matrix(sample_bmfrm_prior(p$spec, prod(n)), n[1L], n[2L])
        x <- sweep(x, 1L, rowMeans(x), "-")
        pars[[bias_names[j]]] <- unname(sweep(x, 2L, colMeans(x), "-"))
      }
    }
    eta <- pars$theta[ids$person]
    for (f in setdiff(spec$facets_main, "person")) eta <- eta - pars[[f]][ids[[f]]]
    for (j in seq_along(bias_names)) {
      f <- strsplit(spec$facets_bias[j], ":", fixed = TRUE)[[1L]]
      eta <- eta + pars[[bias_names[j]]][cbind(ids[[f[1L]]], ids[[f[2L]]])]
    }
    if (is_hetero) {
      tf_ids <- ids[[spec$threshold_facet]]
      tf_lv <- labels[[spec$threshold_facet]]
      u <- stats::runif(nrow(design))
      y <- integer(nrow(design))
      for (j in seq_along(tf_lv)) {
        loc <- which(tf_ids == j)
        if (!length(loc)) next
        Kj <- as.integer(K_by_level[[tf_lv[j]]])
        prob_j <- rsm_probabilities(eta[loc], pars$tau[[j]])
        cum_j <- t(apply(prob_j, 1L, cumsum))
        y[loc] <- 1L + rowSums(u[loc] > cum_j[, -Kj, drop = FALSE])
      }
    } else {
      prob <- if (family == "pcm") pcm_probabilities(eta, pars$tau[ids[[spec$threshold_facet]], , drop = FALSE]) else rsm_probabilities(eta, pars$tau)
      cum <- t(apply(prob, 1L, cumsum))
      y <- 1L + rowSums(stats::runif(nrow(design)) > cum[, -K_max, drop = FALSE])
    }
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
  if (distribution == "double_exponential" && length(args) == 2L && args[2L] > 0) {
    return(args[1L] + stats::rexp(n, rate = 1 / args[2L]) * sample(c(-1, 1), n, replace = TRUE))
  }
  stop("Simulation supports normal, student_t and double_exponential numeric priors with positive scales")
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
