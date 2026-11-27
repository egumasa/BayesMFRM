# Unified model specification shared by simulation, Stan data preparation,
# Stan code generation, and MCMC diagnostics.
# Supports hierarchical main facets (person + any centered facet), hierarchical
# PCM steps (common steps + level deviations), and hierarchical two-way
# interactions with optional cell-space local shrinkage (horseshoe / student_t),
# as well as mixed Gaussian + non-Gaussian priors without silent path switches.

validate_hierarchical <- function(hierarchical, facets_main, facets_bias, family) {
  if (!is.character(hierarchical) || anyNA(hierarchical) || any(!nzchar(hierarchical)) || anyDuplicated(hierarchical)) {
    stop("hierarchical must be a character vector of unique facet names")
  }
  hierarchical <- match_interaction_names(hierarchical, facets_bias)
  if (anyDuplicated(hierarchical)) stop("hierarchical must contain unique terms, including reversed interaction names")
  allowed_targets <- c(facets_main, facets_bias, if (identical(family, "pcm")) "tau")
  if ("tau" %in% hierarchical && !identical(family, "pcm")) {
    stop("Hierarchical step prior ('tau') is supported for PCM family only")
  }
  if (any(!hierarchical %in% allowed_targets)) {
    stop("hierarchical contains a facet absent from the formula: ",
         paste(setdiff(hierarchical, allowed_targets), collapse = ", "))
  }
  if (length(hierarchical) && !family %in% c("rsm", "pcm")) {
    stop("Hierarchical priors support RSM/PCM only")
  }
  unname(hierarchical)
}

is_normal_prior_literal <- function(spec) {
  x <- tryCatch(parse_bmfrm_prior_literal(spec), error = function(e) NULL)
  !is.null(x) && x$distribution == "normal" && length(x$args) == 2L && x$args[2L] > 0
}

is_reserved_stan_identifier <- function(ids) {
  if (!length(ids)) return(logical(0))
  reserved <- c(
    "theta", "tau", "tau_common", "tau_dev", "mu_theta", "sigma_theta",
    "N", "K", "X", "log_lik", "y_rep", "mu", "sigma2", "mu_fair",
    "eta", "eta_full", "eta_fair", "logits", "logits_full", "logits_fair",
    "p_full", "p_fair", "ex", "ex2", "ex_fair",
    "data", "model", "parameters", "transformed", "generated", "quantities",
    "functions", "real", "int", "vector", "matrix", "array", "for", "in",
    "if", "else", "while", "return", "target",
    # https://mc-stan.org/docs/reference-manual/arithmetic-expressions.html
    "repeat", "until", "then", "true", "false", "struct", "typedef", "export",
    "auto", "extern", "var", "static", "lower", "upper", "offset", "multiplier",
    "complex", "simplex", "unit_vector", "sum_to_zero_vector", "sum_to_zero_matrix",
    "ordered", "positive_ordered", "row_vector", "cholesky_factor_corr",
    "column_stochastic_matrix", "row_stochastic_matrix", "cholesky_factor_cov",
    "corr_matrix", "cov_matrix", "print", "reject", "profile", "fatal_error", "jacobian",
    "fvar", "STAN_MAJOR", "STAN_MINOR", "STAN_PATCH",
    "STAN_MATH_MAJOR", "STAN_MATH_MINOR", "STAN_MATH_PATCH",
    "alignas", "alignof", "and", "and_eq", "asm", "bitand", "bitor", "bool",
    "break", "case", "catch", "char", "char16_t", "char32_t", "class", "compl",
    "const", "constexpr", "const_cast", "continue", "decltype", "default", "delete",
    "do", "double", "dynamic_cast", "enum", "explicit", "float", "friend", "goto",
    "inline", "long", "mutable", "namespace", "new", "noexcept", "not", "not_eq",
    "nullptr", "operator", "or", "or_eq", "private", "protected", "public", "register",
    "reinterpret_cast", "short", "signed", "sizeof", "static_assert", "static_cast",
    "switch", "template", "this", "thread_local", "throw", "try", "typeid", "typename",
    "union", "unsigned", "using", "virtual", "void", "volatile", "wchar_t", "xor", "xor_eq",
    # Local identifiers in the generated model, including ragged PCM blocks.
    "threshold", "n", "k", "i", "j", "mj", "rj", "dj", "s", "rs", "fs", "raw_j", "Kj", "t_start",
    "tau_step", "centered", "raw_eff", "lam_tilde"
  )
  (ids %in% reserved) |
    grepl("__$|_raw$|_id$|_n$|^(B_|J_|K_|N_tau_|tau_|mean_|fixed_sd_|scale_factor_|prior_sd_|ability_|estimate_|sigma_|lambda_|c2_)", ids) |
    !grepl("^[A-Za-z][A-Za-z0-9_]*$", ids)
}

validate_facet_names <- function(facets_main) {
  facets <- setdiff(facets_main, "person")
  if (any(is_reserved_stan_identifier(facets))) {
    stop("Facet names must be nonreserved Stan identifiers")
  }
  invisible(TRUE)
}

# Internal normalization shared by generation, Stan data, and Stan code.
# Public configuration lives in hierarchical, priors, and save_gq.
resolve_bmfrm_spec <- function(spec, priors = NULL) {
  validate_facet_names(spec$facets_main)
  validate_interaction_terms(spec$facets_main, spec$facets_bias)
  spec$hierarchical <- validate_hierarchical(spec$hierarchical %||% character(), spec$facets_main, spec$facets_bias, spec$family)
  spec$save_gq <- match.arg(spec$save_gq %||% "full", c("full", "minimal"))
  if (isTRUE(spec$custom_stan)) {
    if (!is.null(priors)) stop("priors cannot be supplied with a custom Stan model; specify priors inside the Stan file")
    if (length(spec$hierarchical)) stop("Hierarchical ability and marginal scaling require generated Stan code")
    # The Stan file owns its priors; package defaults do not describe that model.
    spec$priors <- NULL
    spec$gaussian_model <- NULL
    return(spec)
  }
  priors <- validate_priors(priors %||% spec$priors, spec$facets_main, spec$facets_bias, spec$family, spec$hierarchical)
  spec$priors <- priors
  spec$gaussian_model <- NULL

  marginal <- any(vapply(priors, function(p) identical(p$scale, "marginal"), logical(1)))
  gaussian_pcm <- identical(spec$family, "pcm") && !isTRUE(spec$custom_stan) &&
    all(vapply(priors, function(p) is_normal_prior_literal(p$spec), logical(1)))
  has_shrinkage <- any(vapply(priors, function(p) !is.null(p$local_shrinkage) && p$local_shrinkage != "none", logical(1)))

  # Preserve exact legacy code path when no hierarchical/marginal/PCM-gaussian/shrinkage features are used
  if (!length(spec$hierarchical) && !marginal && !gaussian_pcm && !has_shrinkage) {
    return(spec)
  }

  if (!spec$family %in% c("rsm", "pcm") || !"person" %in% spec$facets_main) {
    stop("Hierarchical priors and marginal scaling support generated RSM/PCM only")
  }

  facets <- setdiff(spec$facets_main, "person")

  keys <- vapply(priors, function(p) paste(p$class, p$facet %||% "", sep = ":"), character(1))
  get_prior <- function(cl, facet = "") priors[[match(paste(cl, facet, sep = ":"), keys)]]

  # 1. Ability block (person)
  if ("person" %in% spec$hierarchical) {
    mu <- rsm_population_prior(get_prior("mean", "person")$spec, "normal", FALSE)
    sd <- rsm_population_prior(get_prior("sd", "person")$spec, "lognormal", TRUE)
    ability_spec <- list(mode = "hierarchical", mean = mu, sd = sd)
  } else {
    tp <- get_prior("theta")
    if (is_normal_prior_literal(tp$spec)) {
      args <- parse_bmfrm_prior_literal(tp$spec)$args
      mu <- rsm_population_prior(args[1L], "normal", FALSE)
      sd <- rsm_population_prior(args[2L], "lognormal", TRUE)
      ability_spec <- list(mode = "gaussian", mean = mu, sd = sd)
    } else {
      # Non-Gaussian person prior (e.g. student_t, double_exponential) allowed alongside other hierarchical/marginal terms
      parse_bmfrm_prior_literal(tp$spec)
      ability_spec <- list(mode = "non_gaussian", spec = tp$spec,
                           mean = rsm_population_prior(0, "normal", FALSE),
                           sd = rsm_population_prior(1, "lognormal", TRUE))
    }
  }

  # 2. Other main facets
  facet_specs <- lapply(facets, function(f) {
    if (f %in% spec$hierarchical) {
      p <- get_prior("sd", f)
      list(facet = f, mode = "hierarchical_sd", sd = positive_sd_prior(p$spec), scale = p$scale)
    } else {
      cl <- if (f %in% c("item", "rater")) f else paste0("facet_", f)
      p <- get_prior(cl)
      if (is_normal_prior_literal(p$spec)) {
        args <- parse_bmfrm_prior_literal(p$spec)$args
        list(facet = f, mode = "gaussian_basis", sd = list(estimated = FALSE, value = unname(args[2L])), scale = p$scale)
      } else {
        if (identical(p$scale, "marginal")) {
          stop("Marginal scaling requires normal(0, positive_sd)")
        }
        parse_bmfrm_prior_literal(p$spec)
        list(facet = f, mode = "non_gaussian_raw", spec = p$spec, scale = "raw")
      }
    }
  })
  names(facet_specs) <- facets

  # 3. Step (tau) block
  tau_p <- get_prior("tau")
  if ("tau" %in% spec$hierarchical) {
    if (isTRUE(spec$K_heterogeneous)) {
      stop(
        "Hierarchical step sharing (hierarchical = 'tau') is not supported in v1 when PCM category ",
        "counts K vary across '", spec$threshold_facet, "' levels."
      )
    }
    # Hierarchical PCM: tau[j, k] = tau_common[k] + tau_dev[j, k], each row sums to zero
    if (!is_normal_prior_literal(tau_p$spec)) {
      stop("Hierarchical PCM common step prior ('tau') requires a normal(0, positive_sd) distribution")
    }
    tau_dev_p <- get_prior("sd", "tau")
    step_spec <- list(
      mode = "hierarchical_pcm",
      common_sd = unname(parse_bmfrm_prior_literal(tau_p$spec)$args[2L]),
      common_scale = tau_p$scale,
      dev_sd = positive_sd_prior(tau_dev_p$spec),
      dev_scale = tau_dev_p$scale
    )
  } else if (is_normal_prior_literal(tau_p$spec)) {
    step_spec <- list(
      mode = "gaussian_basis",
      sd = unname(parse_bmfrm_prior_literal(tau_p$spec)$args[2L]),
      scale = tau_p$scale
    )
  } else {
    if (identical(tau_p$scale, "marginal")) {
      stop("Marginal scaling requires normal(0, positive_sd)")
    }
    parse_bmfrm_prior_literal(tau_p$spec)
    step_spec <- list(mode = "non_gaussian_raw", spec = tau_p$spec, scale = "raw")
  }

  # 4. Interaction (bias) blocks
  biases <- lapply(spec$facets_bias, function(term) {
    is_hier <- term %in% spec$hierarchical
    p <- get_prior(if (is_hier) "sd" else "bias", term)
    ls <- p$local_shrinkage %||% "none"
    bname <- gsub(":", "", term)
    fpair <- strsplit(term, ":", fixed = TRUE)[[1L]]
    if (is_hier) {
      sd_info <- positive_sd_prior(p$spec)
      mode <- if (ls == "none") "gaussian_basis" else paste0("shrinkage_", ls)
      list(term = term, name = bname, facets = fpair,
           hierarchical = TRUE, mode = mode, local_shrinkage = ls,
           sd = sd_info, scale = p$scale)
    } else if (is_normal_prior_literal(p$spec)) {
      args <- parse_bmfrm_prior_literal(p$spec)$args
      list(term = term, name = bname, facets = fpair,
           hierarchical = FALSE, mode = "gaussian_basis", local_shrinkage = "none",
           sd = list(estimated = FALSE, value = unname(args[2L])), scale = p$scale)
    } else {
      if (identical(p$scale, "marginal")) {
        stop("Marginal scaling requires normal(0, positive_sd)")
      }
      parse_bmfrm_prior_literal(p$spec)
      list(term = term, name = bname, facets = fpair,
           hierarchical = FALSE, mode = "non_gaussian_raw", local_shrinkage = "none",
           spec = p$spec, scale = "raw")
    }
  })
  names(biases) <- spec$facets_bias

  # Backward-compatible fields alongside rich per-block specs
  gaussian_facets <- facets[vapply(facet_specs, function(x) x$mode == "gaussian_basis", logical(1))]
  facet_sd_vec <- if (length(gaussian_facets)) {
    vapply(facet_specs[gaussian_facets], function(x) x$sd$value, numeric(1))
  } else {
    stats::setNames(numeric(), character())
  }
  facet_scale_vec <- if (length(gaussian_facets)) {
    vapply(facet_specs[gaussian_facets], function(x) x$scale, character(1))
  } else {
    stats::setNames(character(), character())
  }

  spec$gaussian_model <- list(
    schema = "bayesmfrm.unified_ordinal.v5",
    family = spec$family,
    threshold_facet = spec$threshold_facet,
    K_heterogeneous = isTRUE(spec$K_heterogeneous),
    K_by_level = spec$K_by_level,
    hierarchical = spec$hierarchical,
    ability = ability_spec,
    ability_mean = ability_spec$mean,
    ability_sd = ability_spec$sd,
    facet_specs = facet_specs,
    facet_sd = facet_sd_vec,
    facet_scale = facet_scale_vec,
    step = step_spec,
    step_sd = if (step_spec$mode == "gaussian_basis") step_spec$sd else NA_real_,
    step_scale = step_spec$scale %||% step_spec$common_scale,
    biases = biases
  )
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

# Marginal scale correction factor for double-centered J1 x J2 matrices
double_center_scale_factor <- function(n1, n2, scale) {
  if (!identical(scale, "marginal") || n1 <= 1L || n2 <= 1L) return(1)
  sqrt((n1 / (n1 - 1L)) * (n2 / (n2 - 1L)))
}

gaussian_model_data <- function(model, labels, K, K_by_level = NULL) {
  out <- list()

  # Ability data (preserved verbatim for Gaussian / hierarchical person compatibility)
  if (model$ability$mode %in% c("hierarchical", "gaussian")) {
    mu <- model$ability$mean
    sd <- model$ability$sd
    out <- c(out, list(
      estimate_mu = as.integer(mu$estimated),
      estimate_sigma = as.integer(sd$estimated),
      ability_mean_location = mu$location,
      ability_mean_scale = mu$scale,
      ability_sd_location = sd$location,
      ability_sd_scale = sd$scale
    ))
  }

  # Step (tau) data
  if (isTRUE(model$K_heterogeneous)) {
    tf_labels <- labels[[model$threshold_facet]]
    K_lvl <- as.integer(unname((K_by_level %||% model$K_by_level)[tf_labels]))
    M_lvl <- K_lvl - 1L
    K_max <- max(K_lvl)
    N_tau_steps <- sum(M_lvl)
    tau_start <- as.integer(c(1L, utils::head(cumsum(M_lvl), -1L) + 1L))
    if (model$step$mode == "gaussian_basis") {
      d_lvl <- pmax(K_lvl - 2L, 0L)
      tau_free_start <- as.integer(c(1L, utils::head(cumsum(d_lvl), -1L) + 1L))
      tau_free_start[d_lvl == 0L] <- 1L
      B_tau_rows <- matrix(0, nrow = N_tau_steps, ncol = max(K_max - 2L, 1L))
      for (j in seq_along(K_lvl)) {
        if (d_lvl[j] > 0L) {
          rows_j <- tau_start[j]:(tau_start[j] + M_lvl[j] - 1L)
          B_tau_rows[rows_j, seq_len(d_lvl[j])] <- rsm_prior_basis(M_lvl[j], model$step$scale)
        }
      }
      out$N_tau_free <- as.integer(sum(d_lvl))
      out$tau_free_start <- tau_free_start
      out$tau_free_len <- as.integer(d_lvl)
      out$B_tau_rows <- B_tau_rows
      out$prior_sd_tau <- model$step$sd
    } else {
      tau_raw_len <- ifelse(K_lvl >= 3L, M_lvl, 0L)
      tau_raw_start <- as.integer(c(1L, utils::head(cumsum(tau_raw_len), -1L) + 1L))
      tau_raw_start[tau_raw_len == 0L] <- 1L
      out$N_tau_raw <- as.integer(sum(tau_raw_len))
      out$tau_raw_start <- tau_raw_start
      out$tau_raw_len <- as.integer(tau_raw_len)
    }
  } else if (model$step$mode == "gaussian_basis") {
    out$B_tau <- rsm_prior_basis(K - 1L, model$step$scale)
    out$prior_sd_tau <- model$step$sd
  } else if (model$step$mode == "hierarchical_pcm") {
    out$B_tau_common <- rsm_prior_basis(K - 1L, model$step$common_scale)
    out$prior_sd_tau_common <- model$step$common_sd
    out$B_tau_dev <- rsm_prior_basis(K - 1L, model$step$dev_scale)
    out$estimate_sd_tau <- as.integer(model$step$dev_sd$estimated)
    out$fixed_sd_tau <- if (model$step$dev_sd$estimated) 1 else model$step$dev_sd$value
  }

  # Main facet data
  for (f in names(model$facet_specs)) {
    fs <- model$facet_specs[[f]]
    if (fs$mode == "gaussian_basis") {
      out[[paste0("B_", f)]] <- rsm_prior_basis(length(labels[[f]]), fs$scale)
      out[[paste0("prior_sd_", f)]] <- unname(fs$sd$value)
    } else if (fs$mode == "hierarchical_sd") {
      out[[paste0("B_", f)]] <- rsm_prior_basis(length(labels[[f]]), fs$scale)
      out[[paste0("estimate_sd_", f)]] <- as.integer(fs$sd$estimated)
      out[[paste0("fixed_sd_", f)]] <- if (fs$sd$estimated) 1 else fs$sd$value
    }
  }

  # Interaction data
  for (b in model$biases) {
    n1 <- length(labels[[b$facets[1L]]])
    n2 <- length(labels[[b$facets[2L]]])
    if (b$mode == "gaussian_basis") {
      out[[paste0("B_", b$name, "_1")]] <- rsm_prior_basis(n1, b$scale)
      out[[paste0("B_", b$name, "_2")]] <- rsm_prior_basis(n2, b$scale)
      out[[paste0("estimate_sd_", b$name)]] <- as.integer(b$sd$estimated)
      out[[paste0("fixed_sd_", b$name)]] <- if (b$sd$estimated) 1 else b$sd$value
    } else if (b$mode %in% c("shrinkage_horseshoe", "shrinkage_student_t")) {
      out[[paste0("estimate_sd_", b$name)]] <- as.integer(b$sd$estimated)
      out[[paste0("fixed_sd_", b$name)]] <- if (b$sd$estimated) 1 else b$sd$value
      out[[paste0("scale_factor_", b$name)]] <- double_center_scale_factor(n1, n2, b$scale)
    }
  }

  out
}

sample_gaussian_model <- function(model, labels, K, K_by_level = NULL) {
  # 1. Person ability
  if (model$ability$mode %in% c("hierarchical", "gaussian")) {
    mu <- model$ability$mean
    sd <- model$ability$sd
    mu <- if (mu$estimated) stats::rnorm(1L, mu$location, mu$scale) else mu$location
    sd <- if (sd$estimated) stats::rlnorm(1L, sd$location, sd$scale) else sd$location
    out <- list(theta = stats::rnorm(length(labels$person), mu, sd))
  } else {
    out <- list(theta = sample_bmfrm_prior(model$ability$spec, length(labels$person)))
    mu <- NA_real_; sd <- NA_real_
  }

  # 2. Other main facets
  for (f in names(model$facet_specs)) {
    fs <- model$facet_specs[[f]]
    n <- length(labels[[f]])
    if (fs$mode == "gaussian_basis") {
      out[[f]] <- as.vector(fs$sd$value * rsm_prior_basis(n, fs$scale) %*% stats::rnorm(n - 1L))
    } else if (fs$mode == "hierarchical_sd") {
      s_val <- sample_positive_sd(fs$sd)
      out[[f]] <- as.vector(s_val * rsm_prior_basis(n, fs$scale) %*% stats::rnorm(n - 1L))
      out[[paste0("sigma_", f)]] <- s_val
    } else {
      x <- sample_bmfrm_prior(fs$spec, n)
      out[[f]] <- x - mean(x)
    }
  }

  # 3. Steps (tau)
  if (isTRUE(model$K_heterogeneous)) {
    tf_labels <- labels[[model$threshold_facet]]
    K_lvl <- as.integer(unname((K_by_level %||% model$K_by_level)[tf_labels]))
    tau_list <- vector("list", length(tf_labels))
    names(tau_list) <- tf_labels
    for (j in seq_along(tf_labels)) {
      Kj <- K_lvl[j]
      if (Kj == 2L) {
        tau_list[[j]] <- 0
      } else if (model$step$mode == "gaussian_basis") {
        Bj <- rsm_prior_basis(Kj - 1L, model$step$scale)
        tau_list[[j]] <- as.vector(model$step$sd * Bj %*% stats::rnorm(Kj - 2L))
      } else {
        xj <- sample_bmfrm_prior(model$step$spec, Kj - 1L)
        tau_list[[j]] <- xj - mean(xj)
      }
    }
    out$tau <- tau_list
  } else if (model$step$mode == "gaussian_basis") {
    B <- rsm_prior_basis(K - 1L, model$step$scale)
    if (model$family == "pcm") {
      J <- length(labels[[model$threshold_facet]])
      z <- matrix(stats::rnorm(J * (K - 2L)), J, K - 2L)
      out$tau <- unname(model$step$sd * z %*% t(B))
    } else {
      out$tau <- as.vector(model$step$sd * B %*% stats::rnorm(K - 2L))
    }
  } else if (model$step$mode == "hierarchical_pcm") {
    J <- length(labels[[model$threshold_facet]])
    B_c <- rsm_prior_basis(K - 1L, model$step$common_scale)
    B_d <- rsm_prior_basis(K - 1L, model$step$dev_scale)
    tau_common <- as.vector(model$step$common_sd * B_c %*% stats::rnorm(K - 2L))
    sigma_tau <- sample_positive_sd(model$step$dev_sd)
    z_dev <- matrix(stats::rnorm(J * (K - 2L)), J, K - 2L)
    tau_dev <- unname(sigma_tau * z_dev %*% t(B_d))
    out$tau_common <- tau_common
    out$sigma_tau <- sigma_tau
    out$tau <- unname(sweep(tau_dev, 2L, tau_common, "+"))
  } else {
    if (model$family == "pcm") {
      J <- length(labels[[model$threshold_facet]])
      x <- matrix(sample_bmfrm_prior(model$step$spec, J * (K - 1L)), J, K - 1L)
      out$tau <- unname(x - rowMeans(x))
    } else {
      x <- sample_bmfrm_prior(model$step$spec, K - 1L)
      out$tau <- x - mean(x)
    }
  }

  # Preserve exact RNG order for person population metadata
  if ("person" %in% model$hierarchical) {
    out$mu_theta <- mu
    out$sigma_theta <- sd
  }

  # 4. Interactions
  for (b in model$biases) {
    n <- lengths(labels[b$facets])
    if (b$mode == "gaussian_basis") {
      s_val <- sample_positive_sd(b$sd)
      z <- matrix(stats::rnorm(prod(n - 1L)), n[1L] - 1L, n[2L] - 1L)
      out[[b$name]] <- unname(s_val * rsm_prior_basis(n[1L], b$scale) %*% z %*% t(rsm_prior_basis(n[2L], b$scale)))
      if (b$hierarchical) out[[paste0("sigma_", b$name)]] <- s_val
    } else if (b$mode == "shrinkage_horseshoe") {
      s_val <- sample_positive_sd(b$sd)
      z <- matrix(stats::rnorm(prod(n)), n[1L], n[2L])
      lam <- matrix(abs(stats::rcauchy(prod(n), 0, 1)), n[1L], n[2L])
      # c2 ~ inv_gamma(2, 8) <=> 1 / rgamma(1, shape = 2, rate = 8)
      c2 <- 1 / stats::rgamma(1L, shape = 2, rate = 8)
      lam_tilde <- sqrt((c2 * lam^2) / (c2 + (s_val^2) * (lam^2)))
      sf <- double_center_scale_factor(n[1L], n[2L], b$scale)
      u <- sf * s_val * lam_tilde * z
      u <- sweep(u, 1L, rowMeans(u), "-")
      out[[b$name]] <- unname(sweep(u, 2L, colMeans(u), "-"))
      out[[paste0("sigma_", b$name)]] <- s_val
    } else if (b$mode == "shrinkage_student_t") {
      s_val <- sample_positive_sd(b$sd)
      z <- matrix(stats::rnorm(prod(n)), n[1L], n[2L])
      # lam ~ inv_gamma(1.5, 0.5) => E[lam] = 1, so sqrt(lam) * z is unit-variance
      # Student-t (df = 3), preserving exact conditional cell SD = s_val under scale = "marginal".
      lam <- matrix(1 / stats::rgamma(prod(n), shape = 1.5, rate = 0.5), n[1L], n[2L])
      sf <- double_center_scale_factor(n[1L], n[2L], b$scale)
      u <- sf * s_val * sqrt(lam) * z
      u <- sweep(u, 1L, rowMeans(u), "-")
      out[[b$name]] <- unname(sweep(u, 2L, colMeans(u), "-"))
      out[[paste0("sigma_", b$name)]] <- s_val
    } else {
      x <- matrix(sample_bmfrm_prior(b$spec, prod(n)), n[1L], n[2L])
      x <- sweep(x, 1L, rowMeans(x), "-")
      out[[b$name]] <- unname(sweep(x, 2L, colMeans(x), "-"))
    }
  }

  out
}

model_diagnostic_variables <- function(spec, data_info) {
  model <- spec$gaussian_model
  if (is.null(model)) return(NULL)

  all_other <- names(model$facet_specs)
  active_facets <- all_other[lengths(data_info$facet_labels[all_other]) > 1L]

  K_max <- get_spec_K_max(spec)
  is_hetero <- isTRUE(spec$K_heterogeneous)
  tau_vars <- character(0)
  if (is_hetero) {
    tf_labels <- data_info$facet_labels[[spec$threshold_facet]]
    K_lvl <- as.integer(unname((data_info$K_by_level %||% spec$K_by_level)[tf_labels]))
    tau_len <- K_lvl - 1L
    tau_start <- as.integer(c(1L, utils::head(cumsum(tau_len), -1L) + 1L))
    active_idx <- unlist(lapply(seq_along(K_lvl), function(j) {
      if (K_lvl[j] >= 3L) seq.int(tau_start[j], tau_start[j] + tau_len[j] - 1L) else integer(0)
    }), use.names = FALSE)
    if (length(active_idx) > 0L) tau_vars <- paste0("tau[", active_idx, "]")
  } else if (K_max > 2L) {
    tau_vars <- "tau"
  }

  variable <- c(
    "theta",
    active_facets,
    tau_vars,
    if (model$ability$mode == "hierarchical" && model$ability$mean$estimated) "mu_theta",
    if (model$ability$mode == "hierarchical" && model$ability$sd$estimated) "sigma_theta"
  )
  raw <- c(
    "theta_raw",
    if (length(active_facets)) paste0(active_facets, "_raw"),
    if (K_max > 2L && model$step$mode != "hierarchical_pcm") "tau_raw",
    if (model$ability$mode == "hierarchical" && model$ability$mean$estimated) "mu_theta_raw",
    if (model$ability$mode == "hierarchical" && model$ability$sd$estimated) "log_sigma_theta_raw"
  )

  # Hierarchical main facet SDs (included even when J_f == 1 because unseen-level prediction uses sigma_f)
  for (f in all_other) {
    fs <- model$facet_specs[[f]]
    if (fs$mode == "hierarchical_sd" && fs$sd$estimated) {
      variable <- c(variable, paste0("sigma_", f))
      raw <- c(raw, paste0("sigma_", f, "_raw"))
    }
  }

  # Hierarchical PCM step parameters
  if (model$step$mode == "hierarchical_pcm" && K_max > 2L) {
    variable <- c(variable, "tau_common")
    raw <- c(raw, "tau_common_raw", "tau_dev_raw")
    if (model$step$dev_sd$estimated) {
      variable <- c(variable, "sigma_tau")
      raw <- c(raw, "sigma_tau_raw")
    }
  }

  # Interactions
  for (b in model$biases) {
    active_b <- all(lengths(data_info$facet_labels[b$facets]) > 1L)
    if (active_b) {
      variable <- c(variable, b$name)
      raw <- c(raw, paste0(b$name, "_raw"))
      if (b$mode %in% c("shrinkage_horseshoe", "shrinkage_student_t")) {
        raw <- c(raw, paste0("lambda_", b$name))
      }
    }
    if (b$mode == "shrinkage_horseshoe") {
      raw <- c(raw, paste0("c2_", b$name))
    }
    if (b$hierarchical && b$sd$estimated) {
      variable <- c(variable, paste0("sigma_", b$name))
      raw <- c(raw, paste0("sigma_", b$name, "_raw"))
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
  if (any(lengths(parts) != 2L) || any(vapply(parts, anyDuplicated, integer(1)) > 0L)) {
    stop("Only two-way interactions between distinct main facets are supported")
  }
  if (any(!unlist(parts) %in% facets)) stop("Interaction components must be present as main facets")
  if (anyDuplicated(interaction_key(terms))) stop("Duplicate interaction term")
  names <- gsub(":", "", terms)
  generated <- c(facets, names, paste0("sigma_", names))
  if (anyDuplicated(generated) || any(is_reserved_stan_identifier(names)) ||
      any(!grepl("^[A-Za-z][A-Za-z0-9_]*$", generated))) {
    stop("Facet and interaction names collide or are invalid Stan identifiers")
  }
  invisible(TRUE)
}

# Positive-support SD prior shared by hierarchical main facets, PCM step
# deviations, and two-way interactions.
positive_sd_prior <- function(x) {
  if (is.numeric(x) && length(x) == 1L && is.finite(x) && x > 0) return(list(estimated = FALSE, value = unname(x)))
  if (!is.character(x) || length(x) != 1L || is.na(x)) stop("Hierarchical SD must be positive or a supported positive-scale prior")
  p <- parse_bmfrm_prior_literal(x)
  a <- p$args; dist <- p$distribution
  valid <- (dist %in% c("normal", "cauchy") && length(a) == 2L && a[1L] == 0 && a[2L] > 0) ||
    (dist == "student_t" && length(a) == 3L && a[1L] > 0 && a[2L] == 0 && a[3L] > 0) ||
    (dist == "lognormal" && length(a) == 2L && a[2L] > 0)
  if (!valid) stop("Hierarchical SD prior requires normal(0,s), student_t(df,0,s), cauchy(0,s), or lognormal(m,s), with positive scales")
  list(estimated = TRUE, spec = x, distribution = dist, args = a)
}

# Backward-compatible alias
interaction_sd_prior <- positive_sd_prior

sample_positive_sd <- function(p) {
  if (!p$estimated) return(p$value)
  a <- p$args
  switch(p$distribution,
    normal = abs(stats::rnorm(1L, 0, a[2L])),
    student_t = abs(a[3L] * stats::rt(1L, a[1L])),
    cauchy = abs(stats::rcauchy(1L, 0, a[2L])),
    lognormal = stats::rlnorm(1L, a[1L], a[2L])
  )
}

sample_interaction_sd <- sample_positive_sd
