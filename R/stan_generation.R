#' Build Stan code from specification and priors
#'
#' @param spec bmfrm_spec object from parse_bmfrm_formula()
#' @param priors list of bmfrm_prior objects
#' @param save_gq Either `"full"` (default; saves observation-level `log_lik`,
#'   `y_rep`, `mu`, `sigma2`, `mu_fair` in Stan `generated quantities`) or
#'   `"minimal"` (omits observation-level `generated quantities` to save disk
#'   and memory on large datasets; downstream S3 methods compute quantities on
#'   demand in R).
#'
#' @return character string with complete Stan code
#' @export
build_stan_code <- function(spec, priors = NULL, save_gq = NULL) {
  if (!is.null(save_gq)) {
    spec$save_gq <- match.arg(save_gq, c("full", "minimal"))
  }

  spec <- resolve_bmfrm_spec(spec, priors)
  priors <- spec$priors

  # Load appropriate template based on family (and heterogeneous PCM flag)
  template_path <- get_stan_template_path(spec$family, heterogeneous = isTRUE(spec$K_heterogeneous))

  stan_code <- readLines(template_path, warn = FALSE)
  stan_code <- paste(stan_code, collapse = "\n")

  # For PCM, replace threshold facet placeholder first
  if (spec$family == "pcm" && !is.null(spec$threshold_facet)) {
    stan_code <- gsub("\\{\\{THRESHOLD_FACET\\}\\}", spec$threshold_facet, stan_code)
  }

  # Generate replacements for each placeholder
  replacements <- ordinal_stan_replacements(spec)

  # 1. FACET_DATA_DECL
  replacements$FACET_DATA_DECL <- generate_facet_data_decl(spec)

  # 2. FACET_PARAM_DECL
  replacements$FACET_PARAM_DECL <- generate_facet_param_decl(spec)

  # 3. BIAS_PARAM_DECL
  replacements$BIAS_PARAM_DECL <- generate_bias_param_decl(spec)

  # 4. FACET_TRANSFORM
  replacements$FACET_TRANSFORM <- generate_facet_transform(spec)

  # 5. BIAS_TRANSFORM
  replacements$BIAS_TRANSFORM <- generate_bias_transform(spec)

  # 6. PRIOR_BLOCK
  replacements$PRIOR_BLOCK <- generate_prior_block(spec, priors)

  # 7. INDEX_EXTRACT
  replacements$INDEX_EXTRACT <- generate_index_extract(spec)

  # 8. ETA_SUM
  replacements$ETA_SUM <- generate_eta_sum(spec, fair = FALSE)

  # 9. ETA_SUM_FAIR
  replacements$ETA_SUM_FAIR <- generate_eta_sum(spec, fair = TRUE)

  # 10. FAIR_TAU_EXPR (PCM only: use tau_common when hierarchical PCM steps lie on an adjusted nuisance facet)
  use_common_fair_tau <- identical(spec$family, "pcm") &&
    !is.null(spec$gaussian_model) &&
    identical(spec$gaussian_model$step$mode, "hierarchical_pcm") &&
    (spec$threshold_facet %in% c("rater", "interlocutor"))
  replacements$FAIR_TAU_EXPR <- if (use_common_fair_tau) "tau_common[k-1]" else "tau[threshold_n, k-1]"

  # Replace all placeholders
  for (placeholder in names(replacements)) {
    pattern <- paste0("// \\{\\{", placeholder, "\\}\\}")
    replacement <- replacements[[placeholder]]
    stan_code <- gsub(pattern, replacement, stan_code, fixed = FALSE)
  }

  # Strip generated quantities block for minimal output mode
  if (identical(spec$save_gq, "minimal")) {
    stan_code <- sub("\ngenerated quantities \\{[\\s\\S]*$", "", stan_code, perl = TRUE)
  }

  stan_code
}

#' Generate data block declarations for facets
#' @keywords internal
generate_facet_data_decl <- function(spec) {
  int_lines <- ""
  for (facet in spec$facets_main) {
    int_lines <- c(int_lines,
      paste0("  int<lower=1> J_", facet, ";")
    )
  }

  array_lines <- ""
  for (facet in spec$facets_main) {
    array_lines <- c(
      array_lines,
      paste0("  array[N] int<lower=1, upper=J_", facet, "> ", facet, "_id", ";")
    )
  }

  paste(
    c(int_lines, "", array_lines),
    collapse = "\n"
  )
}

#' Generate parameter declarations for main facets
#' @keywords internal
generate_facet_param_decl <- function(spec) {
  lines <- ""

  # Skip person (theta) as it's always present in template
  other_facets <- setdiff(spec$facets_main, "person")

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

  if (length(spec$facets_bias) > 0) {
    for (bias_facet in spec$facets_bias) {
      components <- strsplit(bias_facet, ":")[[1]]
      f1 <- components[1]
      f2 <- components[2]
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
        next
      }
      lines <- c(lines, paste0("  matrix[J_", f1, ", J_", f2, "] ", bias_name, "_raw;"))
    }
  }

  paste(lines, collapse = "\n")
}

#' Generate transformed parameters for main facets
#' @keywords internal
generate_facet_transform <- function(spec) {
  other_facets <- setdiff(spec$facets_main, "person")

  vector_lines <- ""
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
      # non_gaussian_raw falls through to standard centering
    }
    vector_lines <- c(
      vector_lines,
      paste0("  vector[J_", facet, "] ", facet, ";")
    )
    centering_lines <- c(
      centering_lines,
      paste0("  // ", facet, " (sum-to-zero constraint)"),
      "  {",
      paste0("    real mean_", facet, " = mean(", facet, "_raw);"),
      paste0("    ", facet, " = ", facet, "_raw - mean_", facet, ";"),
      "  }"
    )
  }

  paste(
    c(vector_lines, centering_lines),
    collapse = "\n"
  )
}

#' Generate transformed parameters for bias facets
#' @keywords internal
generate_bias_transform <- function(spec) {
  lines <- ""

  if (length(spec$facets_bias) > 0) {
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
            paste0("  // ", bias_facet, " regularized horseshoe bias (double-centered in cell space)"),
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
            paste0("    for (i in 1:rows(raw_eff)) centered[i, ] = raw_eff[i, ] - mean(raw_eff[i, ]);"),
            paste0("    for (j in 1:cols(centered)) ", bias_name, "[:, j] = centered[:, j] - mean(centered[:, j]);"),
            "  }"
          )
          next
        } else if (b$mode == "shrinkage_student_t") {
          lines <- c(
            lines,
            paste0("  // ", bias_facet, " Student-t local shrinkage bias (double-centered in cell space)"),
            paste0("  real<lower=0> sigma_", bias_name, " = estimate_sd_", bias_name,
                   " ? sigma_", bias_name, "_raw[1] : fixed_sd_", bias_name, ";"),
            paste0("  matrix[", dims, "] ", bias_name, ";"),
            "  {",
            paste0("    matrix[", dims, "] raw_eff = scale_factor_", bias_name,
                   " * sigma_", bias_name, " * (sqrt(lambda_", bias_name, ") .* ", bias_name, "_raw);"),
            paste0("    matrix[", dims, "] centered;"),
            paste0("    for (i in 1:rows(raw_eff)) centered[i, ] = raw_eff[i, ] - mean(raw_eff[i, ]);"),
            paste0("    for (j in 1:cols(centered)) ", bias_name, "[:, j] = centered[:, j] - mean(centered[:, j]);"),
            "  }"
          )
          next
        }
        # non_gaussian_raw falls through to standard double centering
      }
      lines <- c(
        lines,
        paste0("  // ", bias_facet, " bias (sum-to-zero constraints)"),
        paste0("  matrix[", dims, "] ", bias_name, ";"),
        "  {",
        paste0("    matrix[rows(", bias_name, "_raw), cols(", bias_name, "_raw)] centered;"),
        paste0("    // Row-wise centering"),
        paste0("    for (i in 1:rows(", bias_name, "_raw)) {"),
        paste0("      centered[i, ] = ", bias_name, "_raw[i, ] - mean(", bias_name, "_raw[i, ]);"),
        "    }",
        paste0("    // Column-wise centering"),
        paste0("    for (j in 1:cols(centered)) {"),
        paste0("      ", bias_name, "[:, j] = centered[:, j] - mean(centered[:, j]);"),
        "    }",
        "  }"
      )
    }
  }

  paste(lines, collapse = "\n")
}

#' Generate prior statements
#' @keywords internal
generate_prior_block <- function(spec, priors) {
  if (!is.null(spec$gaussian_model)) {
    gm <- spec$gaussian_model
    lines <- character(0)

    # 1. Ability priors (preserve exact ordering for Gaussian/hierarchical compatibility)
    if (gm$ability$mode %in% c("hierarchical", "gaussian")) {
      lines <- c(lines, "  theta_raw ~ std_normal();", "  mu_theta_raw ~ std_normal();", "  log_sigma_theta_raw ~ std_normal();")
    } else {
      lines <- c(lines, paste0("  theta_raw ~ ", gm$ability$spec, ";"))
    }

    # 2. Step (tau) priors
    if (gm$step$mode == "gaussian_basis") {
      tau_lhs <- if (spec$family == "pcm" && !isTRUE(spec$K_heterogeneous)) "to_vector(tau_raw)" else "tau_raw"
      lines <- c(lines, paste0("  ", tau_lhs, " ~ std_normal();"))
    } else if (gm$step$mode == "hierarchical_pcm") {
      lines <- c(
        lines,
        "  tau_common_raw ~ std_normal();",
        "  to_vector(tau_dev_raw) ~ std_normal();"
      )
      if (gm$step$dev_sd$estimated) {
        lines <- c(lines, paste0("  sigma_tau_raw ~ ", gm$step$dev_sd$spec, ";"))
      }
    } else {
      tau_lhs <- if (spec$family == "pcm" && !isTRUE(spec$K_heterogeneous)) "to_vector(tau_raw)" else "tau_raw"
      lines <- c(lines, paste0("  ", tau_lhs, " ~ ", gm$step$spec, ";"))
    }

    # 3. Main facet priors
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

    # 4. Interaction priors
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
          paste0("  c2_", b$name, " ~ inv_gamma(2, 8);")
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

  lines <- ""
  for (prior_obj in priors) {
    if (prior_obj$class == "theta") {
      lines <- c(lines, paste0("  theta_raw ~ ", prior_obj$spec, ";"))
    } else if (prior_obj$class == "item") {
      lines <- c(lines, paste0("  item_raw ~ ", prior_obj$spec, ";"))
    } else if (prior_obj$class == "rater") {
      lines <- c(lines, paste0("  rater_raw ~ ", prior_obj$spec, ";"))
    } else if (prior_obj$class == "tau") {
      if (spec$family == "pcm" && !isTRUE(spec$K_heterogeneous)) {
        lines <- c(lines, paste0("  to_vector(tau_raw) ~ ", prior_obj$spec, ";"))
      } else {
        lines <- c(lines, paste0("  tau_raw ~ ", prior_obj$spec, ";"))
      }
    } else if (prior_obj$class == "bias") {
      bias_name <- gsub(":", "", prior_obj$facet)
      lines <- c(lines, paste0("  to_vector(", bias_name, "_raw) ~ ", prior_obj$spec, ";"))
    } else if (grepl("^facet_", prior_obj$class)) {
      facet_name <- gsub("^facet_", "", prior_obj$class)
      lines <- c(lines, paste0("  ", facet_name, "_raw ~ ", prior_obj$spec, ";"))
    }
  }

  paste(lines, collapse = "\n")
}

ordinal_stan_replacements <- function(spec) {
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
        ABILITY_PARAM_DECL = "vector[J_person] theta_raw;",
        ABILITY_TRANSFORM = "vector[J_person] theta = theta_raw;",
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
      ABILITY_PARAM_DECL = "vector[J_person] theta_raw;",
      ABILITY_TRANSFORM = "vector[J_person] theta = theta_raw;",
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

  # 1. Ability fields/decls/transforms
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
      "vector[J_person] theta_raw;",
      "vector[estimate_mu] mu_theta_raw;",
      "vector[estimate_sigma] log_sigma_theta_raw;"
    ), collapse = "\n  ")
    ability_transform <- paste(c(
      "real mu_theta = estimate_mu ? ability_mean_location + ability_mean_scale * mu_theta_raw[1] : ability_mean_location;",
      "real<lower=0> sigma_theta = estimate_sigma ? exp(ability_sd_location + ability_sd_scale * log_sigma_theta_raw[1]) : ability_sd_location;",
      "vector[J_person] theta = mu_theta + sigma_theta * theta_raw;"
    ), collapse = "\n  ")
  } else {
    ability_param <- "vector[J_person] theta_raw;"
    ability_transform <- "vector[J_person] theta = theta_raw;"
  }

  # 2. Step (tau) fields/decls/transforms
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

  # 3. Main facet data fields
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

  # 4. Interaction data fields
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

#' Generate index extraction code
#' @keywords internal
generate_index_extract <- function(spec) {
  lines <- ""

  other_facets <- setdiff(spec$facets_main, "person")

  for (facet in other_facets) {
    lines <- c(lines, paste0("    int ", facet, "_n = ", facet, "_id[n];"))
  }

  paste(lines, collapse = "\n")
}

#' Generate eta sum (linear predictor)
#' @keywords internal
generate_eta_sum <- function(spec, fair = FALSE) {
  lines <- ""

  other_facets <- setdiff(spec$facets_main, "person")

  for (facet in other_facets) {
    if (fair && facet %in% c("rater", "interlocutor")) {
      next
    }
    lines <- c(lines, paste0("+ (- ", facet, "[", facet, "_n])"))
  }

  if (length(spec$facets_bias) > 0) {
    for (bias_facet in spec$facets_bias) {
      components <- strsplit(bias_facet, ":")[[1]]
      f1 <- components[1]
      f2 <- components[2]
      bias_clean <- gsub(":", "", bias_facet)

      if (fair && (f1 %in% c("rater", "interlocutor") || f2 %in% c("rater", "interlocutor"))) {
        next
      }

      lines <- c(lines, paste0("+ ", bias_clean, "[", f1, "_n, ", f2, "_n]"))
    }
  }

  if (length(lines) == 0) {
    return("")
  }

  paste(lines, collapse = "\n               ")
}

#' Get Stan template path
#'
#' @param family Model family: "rsm" (Rating Scale Model) or "pcm" (Partial Credit Model)
#' @param heterogeneous Logical; if `TRUE` and `family = "pcm"`, return the
#'   packed ragged heterogeneous-category PCM template (`bmfrm_pcm_hetero.stan`).
#'
#' @return Path to the Stan template file
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
  if (file.exists(inst_path)) {
    return(inst_path)
  }

  workspace_path <- file.path("BayesMFRM", "inst", "stan", template_name)
  if (file.exists(workspace_path)) {
    return(workspace_path)
  }

  stop("Stan template '", template_name, "' not found. Please check BayesMFRM installation.")
}
