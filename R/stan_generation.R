#' Build Stan code from specification and priors
#'
#' @param spec `bmfrm_spec` object from [parse_bmfrm_formula()]
#' @param priors list of `bmfrm_prior` objects
#' @param save_gq Either `"minimal"` (default; omits observation-level
#'   `generated quantities` so downstream S3 methods stream quantities on demand
#'   in R) or `"full"` (saves `log_lik`, `y_rep`, `mu`, `sigma2`, `mu_fair` in
#'   Stan `generated quantities`).
#'
#' @return character string with complete Stan code
#' @export
build_stan_code <- function(spec, priors = NULL, save_gq = NULL) {
  if (isTRUE(spec$custom_stan)) {
    stop(
      "Cannot generate Stan code from a custom_stan specification because its priors and model structure ",
      "are defined by the external Stan file. Use stancode(fit) or fit$stan_code instead."
    )
  }
  if (!is.null(save_gq)) {
    spec$save_gq <- match.arg(save_gq, c("minimal", "full"))
  }

  spec <- resolve_bmfrm_spec(spec, priors)
  priors <- spec$priors
  pf <- spec_person_facet(spec)

  template_path <- get_stan_template_path(spec$family, heterogeneous = isTRUE(spec$K_heterogeneous))
  stan_code <- paste(readLines(template_path, warn = FALSE), collapse = "\n")

  stan_code <- gsub("\\{\\{PERSON_FACET\\}\\}", pf, stan_code)
  if (spec$family == "pcm" && !is.null(spec$threshold_facet)) {
    stan_code <- gsub("\\{\\{THRESHOLD_FACET\\}\\}", spec$threshold_facet, stan_code)
  }

  replacements <- ordinal_stan_replacements(spec)
  replacements$FACET_DATA_DECL <- generate_facet_data_decl(spec)
  replacements$FACET_PARAM_DECL <- generate_facet_param_decl(spec)
  replacements$BIAS_PARAM_DECL <- generate_bias_param_decl(spec)
  replacements$FACET_TRANSFORM <- generate_facet_transform(spec)
  replacements$BIAS_TRANSFORM <- generate_bias_transform(spec)
  replacements$PRIOR_BLOCK <- generate_prior_block(spec, priors)
  replacements$ETA_VEC_MAIN <- generate_eta_vec_main(spec, fair = FALSE)
  replacements$ETA_VEC_BIAS <- generate_eta_vec_bias(spec, target_var = "eta", fair = FALSE)
  replacements$GQ_BLOCK <- if (identical(spec$save_gq, "full")) generate_gq_block(spec) else ""

  for (placeholder in names(replacements)) {
    pattern <- paste0("// \\{\\{", placeholder, "\\}\\}")
    stan_code <- gsub(pattern, replacements[[placeholder]], stan_code, fixed = FALSE)
  }

  sub("\n{3,}$", "\n", stan_code)
}

#' Generate data block declarations for facets
#' @keywords internal
generate_facet_data_decl <- function(spec) {
  int_lines <- vapply(spec$facets_main, function(f) paste0("  int<lower=1> J_", f, ";"), character(1))
  arr_lines <- vapply(spec$facets_main, function(f) paste0("  array[N] int<lower=1, upper=J_", f, "> ", f, "_id;"), character(1))
  paste(c(int_lines, "", arr_lines), collapse = "\n")
}

#' Generate parameter declarations for main facets
#' @keywords internal
generate_facet_param_decl <- function(spec) {
  pf <- spec_person_facet(spec)
  other_facets <- setdiff(spec$facets_main, pf)
  lines <- character(0)

  for (facet in other_facets) {
    if (is.null(spec$gaussian_model)) {
      lines <- c(lines, paste0("  vector[J_", facet, "] ", facet, "_raw;"))
      next
    }
    fs <- spec$gaussian_model$facet_specs[[facet]]
    if (fs$mode == "gaussian_basis") {
      lines <- c(lines, paste0("  vector[J_", facet, " - 1] ", facet, "_raw;"))
    } else if (fs$mode == "hierarchical_sd") {
      lines <- c(
        lines,
        paste0("  vector[J_", facet, " - 1] ", facet, "_raw;"),
        paste0("  vector<lower=0>[estimate_sd_", facet, "] sigma_", facet, "_raw;")
      )
    } else {
      lines <- c(lines, paste0("  vector[J_", facet, "] ", facet, "_raw;"))
    }
  }

  paste(lines, collapse = "\n")
}

#' Generate parameter declarations for bias facets
#' @keywords internal
generate_bias_param_decl <- function(spec) {
  lines <- character(0)
  for (bias_facet in spec$facets_bias) {
    components <- strsplit(bias_facet, ":", fixed = TRUE)[[1L]]
    f1 <- components[1L]
    f2 <- components[2L]
    bias_name <- gsub(":", "", bias_facet)
    if (!is.null(spec$gaussian_model)) {
      b <- spec$gaussian_model$biases[[bias_facet]]
      if (b$mode == "gaussian_basis") {
        lines <- c(
          lines,
          paste0("  matrix[J_", f1, " - 1, J_", f2, " - 1] ", bias_name, "_raw;"),
          paste0("  vector<lower=0>[estimate_sd_", bias_name, "] sigma_", bias_name, "_raw;")
        )
      } else if (b$mode == "shrinkage_horseshoe") {
        lines <- c(
          lines,
          paste0("  matrix[J_", f1, ", J_", f2, "] ", bias_name, "_raw;"),
          paste0("  matrix<lower=0>[J_", f1, ", J_", f2, "] lambda_", bias_name, ";"),
          paste0("  real<lower=0> c2_", bias_name, ";"),
          paste0("  vector<lower=0>[estimate_sd_", bias_name, "] sigma_", bias_name, "_raw;")
        )
      } else if (b$mode == "shrinkage_student_t") {
        lines <- c(
          lines,
          paste0("  matrix[J_", f1, ", J_", f2, "] ", bias_name, "_raw;"),
          paste0("  matrix<lower=0>[J_", f1, ", J_", f2, "] lambda_", bias_name, ";"),
          paste0("  vector<lower=0>[estimate_sd_", bias_name, "] sigma_", bias_name, "_raw;")
        )
      } else {
        lines <- c(lines, paste0("  matrix[J_", f1, ", J_", f2, "] ", bias_name, "_raw;"))
      }
    } else {
      lines <- c(lines, paste0("  matrix[J_", f1, ", J_", f2, "] ", bias_name, "_raw;"))
    }
  }
  paste(lines, collapse = "\n")
}

#' Generate transformed parameters for main facets
#' @keywords internal
generate_facet_transform <- function(spec) {
  pf <- spec_person_facet(spec)
  other_facets <- setdiff(spec$facets_main, pf)
  vector_lines <- character(0)
  centering_lines <- character(0)

  for (facet in other_facets) {
    if (!is.null(spec$gaussian_model)) {
      fs <- spec$gaussian_model$facet_specs[[facet]]
      if (fs$mode == "gaussian_basis") {
        vector_lines <- c(
          vector_lines,
          paste0("  vector[J_", facet, "] ", facet, " = prior_sd_", facet, " * B_", facet, " * ", facet, "_raw;")
        )
        next
      } else if (fs$mode == "hierarchical_sd") {
        vector_lines <- c(
          vector_lines,
          paste0("  real<lower=0> sigma_", facet, " = estimate_sd_", facet,
                 " ? sigma_", facet, "_raw[1] : fixed_sd_", facet, ";"),
          paste0("  vector[J_", facet, "] ", facet, " = sigma_", facet, " * B_", facet, " * ", facet, "_raw;")
        )
        next
      }
    }
    vector_lines <- c(vector_lines, paste0("  vector[J_", facet, "] ", facet, ";"))
    centering_lines <- c(
      centering_lines,
      paste0("  ", facet, " = ", facet, "_raw - mean(", facet, "_raw);")
    )
  }

  paste(c(vector_lines, centering_lines), collapse = "\n")
}

#' Generate transformed parameters for bias facets
#' @keywords internal
generate_bias_transform <- function(spec) {
  lines <- character(0)
  for (bias_facet in spec$facets_bias) {
    bias_name <- gsub(":", "", bias_facet)
    components <- strsplit(bias_facet, ":", fixed = TRUE)[[1L]]
    dims <- paste0("J_", components, collapse = ", ")
    if (!is.null(spec$gaussian_model)) {
      b <- spec$gaussian_model$biases[[bias_facet]]
      if (b$mode == "gaussian_basis") {
        lines <- c(
          lines,
          paste0("  real<lower=0> sigma_", bias_name, " = estimate_sd_", bias_name,
                 " ? sigma_", bias_name, "_raw[1] : fixed_sd_", bias_name, ";"),
          paste0("  matrix[", dims, "] ", bias_name, " = sigma_", bias_name,
                 " * B_", bias_name, "_1 * ", bias_name, "_raw * B_", bias_name, "_2';")
        )
        next
      } else if (b$mode == "shrinkage_horseshoe") {
        lines <- c(
          lines,
          paste0("  real<lower=0> sigma_", bias_name, " = estimate_sd_", bias_name,
                 " ? sigma_", bias_name, "_raw[1] : fixed_sd_", bias_name, ";"),
          paste0("  matrix[", dims, "] ", bias_name, ";"),
          "  {",
          paste0("    matrix[", dims, "] lam_tilde = sqrt(c2_", bias_name,
                 " * square(lambda_", bias_name, ") ./ (c2_", bias_name,
                 " + square(sigma_", bias_name, ") * square(lambda_", bias_name, ")));"),
          paste0("    matrix[", dims, "] raw_eff = scale_factor_", bias_name,
                 " * sigma_", bias_name, " * (lam_tilde .* ", bias_name, "_raw);"),
          paste0("    matrix[", dims, "] centered;"),
          "    for (i in 1:rows(raw_eff)) centered[i, ] = raw_eff[i, ] - mean(raw_eff[i, ]);",
          paste0("    for (j in 1:cols(centered)) ", bias_name, "[:, j] = centered[:, j] - mean(centered[:, j]);"),
          "  }"
        )
        next
      } else if (b$mode == "shrinkage_student_t") {
        lines <- c(
          lines,
          paste0("  real<lower=0> sigma_", bias_name, " = estimate_sd_", bias_name,
                 " ? sigma_", bias_name, "_raw[1] : fixed_sd_", bias_name, ";"),
          paste0("  matrix[", dims, "] ", bias_name, ";"),
          "  {",
          paste0("    matrix[", dims, "] raw_eff = scale_factor_", bias_name,
                 " * sigma_", bias_name, " * (sqrt(lambda_", bias_name, ") .* ", bias_name, "_raw);"),
          paste0("    matrix[", dims, "] centered;"),
          "    for (i in 1:rows(raw_eff)) centered[i, ] = raw_eff[i, ] - mean(raw_eff[i, ]);",
          paste0("    for (j in 1:cols(centered)) ", bias_name, "[:, j] = centered[:, j] - mean(centered[:, j]);"),
          "  }"
        )
        next
      }
    }
    lines <- c(
      lines,
      paste0("  matrix[", dims, "] ", bias_name, ";"),
      "  {",
      paste0("    matrix[rows(", bias_name, "_raw), cols(", bias_name, "_raw)] centered;"),
      paste0("    for (i in 1:rows(", bias_name, "_raw)) centered[i, ] = ", bias_name, "_raw[i, ] - mean(", bias_name, "_raw[i, ]);"),
      paste0("    for (j in 1:cols(centered)) ", bias_name, "[:, j] = centered[:, j] - mean(centered[:, j]);"),
      "  }"
    )
  }
  paste(lines, collapse = "\n")
}

#' Generate prior statements
#' @keywords internal
generate_prior_block <- function(spec, priors) {
  if (!is.null(spec$gaussian_model)) {
    gm <- spec$gaussian_model
    lines <- character(0)

    if (gm$ability$mode %in% c("hierarchical", "gaussian")) {
      lines <- c(lines, "  theta_raw ~ std_normal();", "  mu_theta_raw ~ std_normal();", "  log_sigma_theta_raw ~ std_normal();")
    } else {
      lines <- c(lines, paste0("  theta_raw ~ ", gm$ability$spec, ";"))
    }

    if (gm$step$mode == "gaussian_basis") {
      tau_lhs <- if (spec$family == "pcm" && !isTRUE(spec$K_heterogeneous)) "to_vector(tau_raw)" else "tau_raw"
      lines <- c(lines, paste0("  ", tau_lhs, " ~ std_normal();"))
    } else if (gm$step$mode == "hierarchical_pcm") {
      lines <- c(lines, "  tau_common_raw ~ std_normal();", "  to_vector(tau_dev_raw) ~ std_normal();")
      if (gm$step$dev_sd$estimated) {
        lines <- c(lines, paste0("  sigma_tau_raw ~ ", gm$step$dev_sd$spec, ";"))
      }
    } else {
      tau_lhs <- if (spec$family == "pcm" && !isTRUE(spec$K_heterogeneous)) "to_vector(tau_raw)" else "tau_raw"
      lines <- c(lines, paste0("  ", tau_lhs, " ~ ", gm$step$spec, ";"))
    }

    for (f in names(gm$facet_specs)) {
      fs <- gm$facet_specs[[f]]
      if (fs$mode == "gaussian_basis") {
        lines <- c(lines, paste0("  ", f, "_raw ~ std_normal();"))
      } else if (fs$mode == "hierarchical_sd") {
        lines <- c(lines, paste0("  ", f, "_raw ~ std_normal();"))
        if (fs$sd$estimated) {
          lines <- c(lines, paste0("  sigma_", f, "_raw ~ ", fs$sd$spec, ";"))
        }
      } else {
        lines <- c(lines, paste0("  ", f, "_raw ~ ", fs$spec, ";"))
      }
    }

    for (b in gm$biases) {
      if (b$mode == "gaussian_basis") {
        lines <- c(lines, paste0("  to_vector(", b$name, "_raw) ~ std_normal();"))
        if (b$sd$estimated) {
          lines <- c(lines, paste0("  sigma_", b$name, "_raw ~ ", b$sd$spec, ";"))
        }
      } else if (b$mode == "shrinkage_horseshoe") {
        lines <- c(
          lines,
          paste0("  to_vector(", b$name, "_raw) ~ std_normal();"),
          paste0("  to_vector(lambda_", b$name, ") ~ cauchy(0, 1);"),
          paste0("  c2_", b$name, " ~ inv_gamma(hs_slab_shape_", b$name, ", hs_slab_rate_", b$name, ");")
        )
        if (b$sd$estimated) {
          lines <- c(lines, paste0("  sigma_", b$name, "_raw ~ ", b$sd$spec, ";"))
        }
      } else if (b$mode == "shrinkage_student_t") {
        lines <- c(
          lines,
          paste0("  to_vector(", b$name, "_raw) ~ std_normal();"),
          paste0("  to_vector(lambda_", b$name, ") ~ inv_gamma(1.5, 0.5);")
        )
        if (b$sd$estimated) {
          lines <- c(lines, paste0("  sigma_", b$name, "_raw ~ ", b$sd$spec, ";"))
        }
      } else {
        lines <- c(lines, paste0("  to_vector(", b$name, "_raw) ~ ", b$spec, ";"))
      }
    }

    return(paste(lines, collapse = "\n"))
  }

  lines <- character(0)
  for (prior_obj in priors) {
    if (prior_obj$class == "theta") {
      lines <- c(lines, paste0("  theta_raw ~ ", prior_obj$spec, ";"))
    } else if (prior_obj$class == "tau") {
      if (spec$family == "pcm" && !isTRUE(spec$K_heterogeneous)) {
        lines <- c(lines, paste0("  to_vector(tau_raw) ~ ", prior_obj$spec, ";"))
      } else {
        lines <- c(lines, paste0("  tau_raw ~ ", prior_obj$spec, ";"))
      }
    } else if (prior_obj$class == "bias") {
      bias_name <- gsub(":", "", prior_obj$facet)
      lines <- c(lines, paste0("  to_vector(", bias_name, "_raw) ~ ", prior_obj$spec, ";"))
    } else {
      key <- prior_target_key(prior_obj, person_facet = spec_person_facet(spec), facets_bias = spec$facets_bias)
      if (startsWith(key, "facet:")) {
        facet_name <- sub("^facet:", "", key)
        lines <- c(lines, paste0("  ", facet_name, "_raw ~ ", prior_obj$spec, ";"))
      }
    }
  }
  paste(lines, collapse = "\n")
}

ordinal_stan_replacements <- function(spec) {
  pf <- spec_person_facet(spec)
  Jp <- paste0("J_", pf)
  J <- paste0("J_", spec$threshold_facet)
  is_hetero <- isTRUE(spec$K_heterogeneous)
  if (is.null(spec$gaussian_model)) {
    if (is_hetero) {
      return(list(
        MODEL_DATA = paste(c(
          "  int<lower=1> N_tau_raw;",
          paste0("  array[", J, "] int<lower=1> tau_raw_start;"),
          paste0("  array[", J, "] int<lower=0> tau_raw_len;")
        ), collapse = "\n"),
        ABILITY_PARAM_DECL = paste0("vector[", Jp, "] theta_raw;"),
        ABILITY_TRANSFORM = paste0("vector[", Jp, "] theta = theta_raw;"),
        THRESHOLD_PARAM_DECL = "vector[N_tau_raw] tau_raw;",
        THRESHOLD_TRANSFORM = paste(c(
          "vector[N_tau_steps] tau = rep_vector(0.0, N_tau_steps);",
          paste0("  for (j in 1:", J, ") {"),
          "    int mj = tau_len[j];",
          "    int rj = tau_raw_len[j];",
          "    if (rj > 0) {",
          "      int s = tau_start[j];",
          "      int rs = tau_raw_start[j];",
          "      vector[mj] raw_j = tau_raw[rs:(rs + mj - 1)];",
          "      tau[s:(s + mj - 1)] = raw_j - mean(raw_j);",
          "    }",
          "  }"
        ), collapse = "\n")
      ))
    }
    return(list(
      MODEL_DATA = "",
      ABILITY_PARAM_DECL = paste0("vector[", Jp, "] theta_raw;"),
      ABILITY_TRANSFORM = paste0("vector[", Jp, "] theta = theta_raw;"),
      THRESHOLD_PARAM_DECL = if (spec$family == "pcm") paste0("matrix[", J, ", K-1] tau_raw;") else "vector[K - 1] tau_raw;",
      THRESHOLD_TRANSFORM = if (spec$family == "pcm") {
        paste0("matrix[", J, ", K - 1] tau;\n  for (j in 1:", J, ") tau[j] = tau_raw[j] - mean(tau_raw[j]);")
      } else {
        "vector[K - 1] tau = tau_raw - mean(tau_raw);"
      }
    ))
  }

  gm <- spec$gaussian_model
  fields <- character(0)

  if (gm$ability$mode %in% c("hierarchical", "gaussian")) {
    fields <- c(
      fields,
      "int<lower=0, upper=1> estimate_mu;",
      "int<lower=0, upper=1> estimate_sigma;",
      "real ability_mean_location;",
      "real<lower=0> ability_mean_scale;",
      "real ability_sd_location;",
      "real<lower=0> ability_sd_scale;"
    )
    ability_param <- paste(c(
      paste0("vector[", Jp, "] theta_raw;"),
      "vector[estimate_mu] mu_theta_raw;",
      "vector[estimate_sigma] log_sigma_theta_raw;"
    ), collapse = "\n  ")
    ability_transform <- paste(c(
      "real mu_theta = estimate_mu ? ability_mean_location + ability_mean_scale * mu_theta_raw[1] : ability_mean_location;",
      "real<lower=0> sigma_theta = estimate_sigma ? exp(ability_sd_location + ability_sd_scale * log_sigma_theta_raw[1]) : ability_sd_location;",
      paste0("vector[", Jp, "] theta = mu_theta + sigma_theta * theta_raw;")
    ), collapse = "\n  ")
  } else {
    ability_param <- paste0("vector[", Jp, "] theta_raw;")
    ability_transform <- paste0("vector[", Jp, "] theta = theta_raw;")
  }

  if (is_hetero) {
    if (gm$step$mode == "gaussian_basis") {
      fields <- c(
        fields,
        "int<lower=1> N_tau_free;",
        paste0("array[", J, "] int<lower=1> tau_free_start;"),
        paste0("array[", J, "] int<lower=0> tau_free_len;"),
        "matrix[N_tau_steps, K_max - 2] B_tau_rows;",
        "real<lower=0> prior_sd_tau;"
      )
      thresh_param <- "vector[N_tau_free] tau_raw;"
      thresh_transform <- paste(c(
        "vector[N_tau_steps] tau = rep_vector(0.0, N_tau_steps);",
        paste0("  for (j in 1:", J, ") {"),
        "    int mj = tau_len[j];",
        "    int dj = tau_free_len[j];",
        "    if (dj > 0) {",
        "      int s = tau_start[j];",
        "      int fs = tau_free_start[j];",
        "      tau[s:(s + mj - 1)] = prior_sd_tau * (B_tau_rows[s:(s + mj - 1), 1:dj] * tau_raw[fs:(fs + dj - 1)]);",
        "    }",
        "  }"
      ), collapse = "\n")
    } else {
      fields <- c(
        fields,
        "int<lower=1> N_tau_raw;",
        paste0("array[", J, "] int<lower=1> tau_raw_start;"),
        paste0("array[", J, "] int<lower=0> tau_raw_len;")
      )
      thresh_param <- "vector[N_tau_raw] tau_raw;"
      thresh_transform <- paste(c(
        "vector[N_tau_steps] tau = rep_vector(0.0, N_tau_steps);",
        paste0("  for (j in 1:", J, ") {"),
        "    int mj = tau_len[j];",
        "    int rj = tau_raw_len[j];",
        "    if (rj > 0) {",
        "      int s = tau_start[j];",
        "      int rs = tau_raw_start[j];",
        "      vector[mj] raw_j = tau_raw[rs:(rs + mj - 1)];",
        "      tau[s:(s + mj - 1)] = raw_j - mean(raw_j);",
        "    }",
        "  }"
      ), collapse = "\n")
    }
  } else if (gm$step$mode == "gaussian_basis") {
    fields <- c(fields, "matrix[K - 1, K - 2] B_tau;", "real<lower=0> prior_sd_tau;")
    thresh_param <- if (spec$family == "pcm") paste0("matrix[", J, ", K - 2] tau_raw;") else "vector[K - 2] tau_raw;"
    thresh_transform <- if (spec$family == "pcm") {
      paste0("matrix[", J, ", K - 1] tau = prior_sd_tau * tau_raw * B_tau';")
    } else {
      "vector[K - 1] tau = prior_sd_tau * B_tau * tau_raw;"
    }
  } else if (gm$step$mode == "hierarchical_pcm") {
    fields <- c(
      fields,
      "matrix[K - 1, K - 2] B_tau_common;",
      "real<lower=0> prior_sd_tau_common;",
      "matrix[K - 1, K - 2] B_tau_dev;",
      "int<lower=0, upper=1> estimate_sd_tau;",
      "real<lower=0> fixed_sd_tau;"
    )
    thresh_param <- paste(c(
      "vector[K - 2] tau_common_raw;",
      paste0("matrix[", J, ", K - 2] tau_dev_raw;"),
      "vector<lower=0>[estimate_sd_tau] sigma_tau_raw;"
    ), collapse = "\n  ")
    thresh_transform <- paste(c(
      "vector[K - 1] tau_common = prior_sd_tau_common * B_tau_common * tau_common_raw;",
      "real<lower=0> sigma_tau = estimate_sd_tau ? sigma_tau_raw[1] : fixed_sd_tau;",
      paste0("matrix[", J, ", K - 1] tau_dev = sigma_tau * tau_dev_raw * B_tau_dev';"),
      paste0("matrix[", J, ", K - 1] tau;"),
      paste0("for (j in 1:", J, ") tau[j] = tau_common' + tau_dev[j];")
    ), collapse = "\n  ")
  } else {
    thresh_param <- if (spec$family == "pcm") paste0("matrix[", J, ", K - 1] tau_raw;") else "vector[K - 1] tau_raw;"
    thresh_transform <- if (spec$family == "pcm") {
      paste0("matrix[", J, ", K - 1] tau;\n  for (j in 1:", J, ") tau[j] = tau_raw[j] - mean(tau_raw[j]);")
    } else {
      "vector[K - 1] tau = tau_raw - mean(tau_raw);"
    }
  }

  for (f in names(gm$facet_specs)) {
    fs <- gm$facet_specs[[f]]
    if (fs$mode == "gaussian_basis") {
      fields <- c(
        fields,
        paste0("matrix[J_", f, ", J_", f, " - 1] B_", f, ";"),
        paste0("real<lower=0> prior_sd_", f, ";")
      )
    } else if (fs$mode == "hierarchical_sd") {
      fields <- c(
        fields,
        paste0("matrix[J_", f, ", J_", f, " - 1] B_", f, ";"),
        paste0("int<lower=0, upper=1> estimate_sd_", f, ";"),
        paste0("real<lower=0> fixed_sd_", f, ";")
      )
    }
  }

  for (b in gm$biases) {
    if (b$mode == "gaussian_basis") {
      fields <- c(
        fields,
        paste0("matrix[J_", b$facets[1L], ", J_", b$facets[1L], " - 1] B_", b$name, "_1;"),
        paste0("matrix[J_", b$facets[2L], ", J_", b$facets[2L], " - 1] B_", b$name, "_2;"),
        paste0("int<lower=0, upper=1> estimate_sd_", b$name, ";"),
        paste0("real<lower=0> fixed_sd_", b$name, ";")
      )
    } else if (b$mode %in% c("shrinkage_horseshoe", "shrinkage_student_t")) {
      fields <- c(
        fields,
        paste0("int<lower=0, upper=1> estimate_sd_", b$name, ";"),
        paste0("real<lower=0> fixed_sd_", b$name, ";"),
        paste0("real<lower=0> scale_factor_", b$name, ";")
      )
      if (b$mode == "shrinkage_horseshoe") {
        fields <- c(
          fields,
          paste0("real<lower=0> hs_slab_shape_", b$name, ";"),
          paste0("real<lower=0> hs_slab_rate_", b$name, ";")
        )
      }
    }
  }

  list(
    MODEL_DATA = paste(fields, collapse = "\n  "),
    ABILITY_PARAM_DECL = ability_param,
    ABILITY_TRANSFORM = ability_transform,
    THRESHOLD_PARAM_DECL = thresh_param,
    THRESHOLD_TRANSFORM = thresh_transform
  )
}

#' Generate vectorized main-facet subtraction expression
#' @keywords internal
generate_eta_vec_main <- function(spec, fair = FALSE) {
  pf <- spec_person_facet(spec)
  adj <- if (fair) spec_adjust_facets(spec) else character(0)
  keep <- setdiff(spec$facets_main, c(pf, adj))
  if (!length(keep)) return("")
  paste(paste0("- ", keep, "[", keep, "_id]"), collapse = "\n                    ")
}

#' Generate 2D interaction accumulation loop on a vectorized linear predictor
#' @keywords internal
generate_eta_vec_bias <- function(spec, target_var = "eta", fair = FALSE) {
  adj <- if (fair) spec_adjust_facets(spec) else character(0)
  terms <- character(0)
  for (bias_facet in spec$facets_bias) {
    comps <- strsplit(bias_facet, ":", fixed = TRUE)[[1L]]
    if (fair && any(comps %in% adj)) next
    bname <- gsub(":", "", bias_facet)
    terms <- c(terms, paste0(bname, "[", comps[1L], "_id[n], ", comps[2L], "_id[n]]"))
  }
  if (!length(terms)) return("")
  paste(c(
    "for (n in 1:N) {",
    paste0("      ", target_var, "[n] += ", paste(terms, collapse = " + "), ";"),
    "    }"
  ), collapse = "\n    ")
}

#' Generate explicit generated quantities block (omitted cleanly when save_gq == "minimal")
#' @keywords internal
generate_gq_block <- function(spec) {
  pf <- spec_person_facet(spec)
  adj <- spec_adjust_facets(spec)
  has_fair_diff <- length(adj) > 0L
  is_hetero <- isTRUE(spec$K_heterogeneous)
  tf <- spec$threshold_facet

  eta_full_main <- generate_eta_vec_main(spec, fair = FALSE)
  eta_full_bias <- generate_eta_vec_bias(spec, target_var = "eta_full", fair = FALSE)
  eta_fair_main <- if (has_fair_diff) generate_eta_vec_main(spec, fair = TRUE) else ""
  eta_fair_bias <- if (has_fair_diff) generate_eta_vec_bias(spec, target_var = "eta_fair", fair = TRUE) else ""

  use_common_fair_tau <- identical(spec$family, "pcm") &&
    !is_hetero &&
    !is.null(spec$gaussian_model) &&
    identical(spec$gaussian_model$step$mode, "hierarchical_pcm") &&
    (tf %in% adj)

  k_upper <- if (is_hetero) "K_max" else "K"
  lines <- c(
    "generated quantities {",
    "  vector[N] log_lik;",
    paste0("  array[N] int<lower=1, upper=", k_upper, "> y_rep;"),
    "  vector[N] mu;",
    "  vector[N] sigma2;",
    "  vector[N] mu_fair;",
    "  {",
    paste0("    vector[N] eta_full = theta[", pf, "_id]",
           if (nzchar(eta_full_main)) paste0(" ", gsub("\n\\s*", " ", eta_full_main)) else "", ";")
  )
  if (nzchar(eta_full_bias)) lines <- c(lines, paste0("    ", eta_full_bias))
  if (has_fair_diff) {
    lines <- c(
      lines,
      paste0("    vector[N] eta_fair = theta[", pf, "_id]",
             if (nzchar(eta_fair_main)) paste0(" ", gsub("\n\\s*", " ", eta_fair_main)) else "", ";")
    )
    if (nzchar(eta_fair_bias)) lines <- c(lines, paste0("    ", eta_fair_bias))
  }

  if (identical(spec$family, "rsm")) {
    lines <- c(
      lines,
      "    vector[K] cum_tau = append_row(0.0, cumulative_sum(tau));",
      "    for (n in 1:N) {",
      "      vector[K] logits_full = cat_seq * eta_full[n] - cum_tau;",
      "      real lse_full = log_sum_exp(logits_full);",
      "      vector[K] p_full = softmax(logits_full);",
      "      real ex = dot_product(cat_val, p_full);",
      "      mu[n] = ex;",
      "      sigma2[n] = fmax(dot_product(cat_val2, p_full) - ex * ex, 0.0);",
      "      log_lik[n] = logits_full[X[n]] - lse_full;",
      "      y_rep[n] = categorical_rng(p_full);",
      if (has_fair_diff) {
        c("      vector[K] p_fair = softmax(cat_seq * eta_fair[n] - cum_tau);",
          "      mu_fair[n] = dot_product(cat_val, p_fair);")
      } else {
        "      mu_fair[n] = ex;"
      },
      "    }"
    )
  } else if (!is_hetero) {
    lines <- c(
      lines,
      paste0("    array[J_", tf, "] vector[K] cum_tau;"),
      paste0("    for (j in 1:J_", tf, ") cum_tau[j] = append_row(0.0, cumulative_sum(tau[j]'));"),
      if (use_common_fair_tau) "    vector[K] cum_tau_fair = append_row(0.0, cumulative_sum(tau_common));",
      "    for (n in 1:N) {",
      paste0("      int threshold_n = ", tf, "_id[n];"),
      "      vector[K] logits_full = cat_seq * eta_full[n] - cum_tau[threshold_n];",
      "      real lse_full = log_sum_exp(logits_full);",
      "      vector[K] p_full = softmax(logits_full);",
      "      real ex = dot_product(cat_val, p_full);",
      "      mu[n] = ex;",
      "      sigma2[n] = fmax(dot_product(cat_val2, p_full) - ex * ex, 0.0);",
      "      log_lik[n] = logits_full[X[n]] - lse_full;",
      "      y_rep[n] = categorical_rng(p_full);",
      if (has_fair_diff) {
        c(paste0("      vector[K] p_fair = softmax(cat_seq * eta_fair[n] - ",
                 if (use_common_fair_tau) "cum_tau_fair" else "cum_tau[threshold_n]", ");"),
          "      mu_fair[n] = dot_product(cat_val, p_fair);")
      } else {
        "      mu_fair[n] = ex;"
      },
      "    }"
    )
  } else {
    lines <- c(
      lines,
      "    vector[N_tau_steps] cum_tau;",
      paste0("    for (j in 1:J_", tf, ") {"),
      "      int s = tau_start[j];",
      "      int mj = tau_len[j];",
      "      cum_tau[s:(s + mj - 1)] = cumulative_sum(tau[s:(s + mj - 1)]);",
      "    }",
      "    for (n in 1:N) {",
      paste0("      int threshold_n = ", tf, "_id[n];"),
      "      int Kj = K_level[threshold_n];",
      "      int t_start = tau_start[threshold_n];",
      "      vector[Kj] logits_full = head(cat_seq, Kj) * eta_full[n];",
      "      logits_full[2:Kj] -= cum_tau[t_start:(t_start + Kj - 2)];",
      "      real lse_full = log_sum_exp(logits_full);",
      "      vector[Kj] p_full = softmax(logits_full);",
      "      real ex = dot_product(head(cat_val, Kj), p_full);",
      "      mu[n] = ex;",
      "      sigma2[n] = fmax(dot_product(head(cat_val2, Kj), p_full) - ex * ex, 0.0);",
      "      log_lik[n] = logits_full[X[n]] - lse_full;",
      "      y_rep[n] = categorical_rng(p_full);",
      if (has_fair_diff) {
        c("      vector[Kj] logits_fair = head(cat_seq, Kj) * eta_fair[n];",
          "      logits_fair[2:Kj] -= cum_tau[t_start:(t_start + Kj - 2)];",
          "      mu_fair[n] = dot_product(head(cat_val, Kj), softmax(logits_fair));")
      } else {
        "      mu_fair[n] = ex;"
      },
      "    }"
    )
  }

  lines <- c(lines, "  }", "}")
  paste(lines, collapse = "\n")
}

#' Verify whether the current working directory is the BayesMFRM package source root
#' @keywords internal
is_bayesmfrm_source_root <- function(dir = ".") {
  desc <- file.path(dir, "DESCRIPTION")
  if (!file.exists(desc)) return(FALSE)
  lines <- readLines(desc, n = 10L, warn = FALSE)
  any(grepl("^Package:\\s*BayesMFRM\\s*$", lines))
}

#' Get Stan template path without unverified relative sibling fallbacks
#' @keywords internal
get_stan_template_path <- function(family = "rsm", heterogeneous = FALSE) {
  template_name <- switch(family,
    "rsm" = "bmfrm_rsm.stan",
    "pcm" = if (isTRUE(heterogeneous)) "bmfrm_pcm_hetero.stan" else "bmfrm_pcm.stan",
    "rating_scale" = "bmfrm_rsm.stan",
    "partial_credit" = if (isTRUE(heterogeneous)) "bmfrm_pcm_hetero.stan" else "bmfrm_pcm.stan",
    stop("Unknown family: ", family, ". Use 'rsm' or 'pcm'.")
  )

  pkg_path <- system.file("stan", template_name, package = "BayesMFRM")
  if (nzchar(pkg_path) && file.exists(pkg_path)) {
    return(pkg_path)
  }

  inst_path <- file.path("inst", "stan", template_name)
  if (is_bayesmfrm_source_root(".") && file.exists(inst_path)) {
    return(normalizePath(inst_path, mustWork = TRUE))
  }

  stop("Stan template '", template_name, "' not found. Please check BayesMFRM installation.")
}
