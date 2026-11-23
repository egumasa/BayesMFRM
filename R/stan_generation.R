#' Build Stan code from specification and priors
#'
#' @param spec bmfrm_spec object from parse_bmfrm_formula()
#' @param priors list of bmfrm_prior objects
#'
#' @return character string with complete Stan code
#' @export
build_stan_code <- function(spec, priors = NULL) {

  spec <- resolve_bmfrm_spec(spec, priors)
  priors <- spec$priors

  # Load appropriate template based on family
  template_path <- get_stan_template_path(spec$family)

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

  # Replace all placeholders
  for (placeholder in names(replacements)) {
    pattern <- paste0("// \\{\\{", placeholder, "\\}\\}")
    replacement <- replacements[[placeholder]]
    stan_code <- gsub(pattern, replacement, stan_code, fixed = FALSE)
  }

  stan_code
}

#' Generate data block declarations for facets
#' @keywords internal
generate_facet_data_decl <- function(spec) {
  empty_line <- character(0) 
  # First collect all int declarations
  int_lines <- ""
  for (facet in spec$facets_main) {
    int_lines <- c(int_lines,
      paste0("  int<lower=1> J_", facet, ";")
    )
  }
  
  # Then collect all array declarations
  array_lines <- ""
  for (facet in spec$facets_main) {
    array_lines <- c(
      array_lines,
      paste0("  array[N] int<lower=1, upper=J_", facet, "> ", facet, "_id", ";")
    )
  }
  
  # Combine both blocks
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
    dim <- if (is.null(spec$gaussian_model)) paste0("J_", facet) else paste0("J_", facet, " - 1")
    lines <- c(lines, paste0("  vector[", dim, "] ", facet, "_raw;"))
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
      # Create matrix parameter for bias
      bias_name <- gsub(":", "", bias_facet)
      lines <- c(lines, paste0("  matrix[J_", f1, ", J_", f2, "] ", bias_name, "_raw;"))
    }
  }
  
  paste(lines, collapse = "\n")
}

#' Generate transformed parameters for main facets
#' @keywords internal
generate_facet_transform <- function(spec) {
  
  # Skip person (theta) as it's handled in template
  other_facets <- setdiff(spec$facets_main, "person")
  
  vector_lines <- ""
  centering_lines <- character(0)

  for (facet in other_facets) {
    if (!is.null(spec$gaussian_model)) {
      vector_lines <- c(vector_lines, paste0("  vector[J_", facet, "] ", facet,
        " = prior_sd_", facet, " * B_", facet, " * ", facet, "_raw;"))
      next
    }
    vector_lines <- c(
               vector_lines,
               paste0("  vector[J_", facet, "] ", facet, ";"))

    centering_lines <- c(
               centering_lines,
               paste0("  // ", facet, " (sum-to-zero constraint)"),
               "  {",
               paste0("    real mean_", facet, " = mean(", facet, "_raw);"),
               paste0("    ", facet, " = ", facet, "_raw - mean_", facet, ";"),
               "  }")
  }

  # Combine both blocks
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
      lines <- c(lines,
                 paste0("  // ", bias_facet, " bias (sum-to-zero constraints)"),
                 paste0("  matrix[rows(", bias_name, "_raw), cols(", bias_name, "_raw)] ", bias_name, ";"),
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
                 "  }")
    }
  }
  
  paste(lines, collapse = "\n")
}

#' Generate prior statements
#' @keywords internal
generate_prior_block <- function(spec, priors) {
  if (!is.null(spec$gaussian_model)) {
    variables <- c("theta_raw", "mu_theta_raw", "log_sigma_theta_raw",
                   if (spec$family == "pcm") "to_vector(tau_raw)" else "tau_raw",
                   if (length(spec$gaussian_model$facet_sd)) paste0(names(spec$gaussian_model$facet_sd), "_raw"))
    return(paste0("  ", variables, " ~ std_normal();", collapse = "\n"))
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
      # For PCM, tau is a matrix so we need to_vector()
      if (spec$family == "pcm") {
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
  if (is.null(spec$gaussian_model)) return(list(
    MODEL_DATA = "",
    ABILITY_PARAM_DECL = "vector[J_person] theta_raw;",
    ABILITY_TRANSFORM = "vector[J_person] theta = theta_raw;",
    THRESHOLD_PARAM_DECL = if (spec$family == "pcm") paste0("matrix[", J, ", K-1] tau_raw;") else "vector[K - 1] tau_raw;",
    THRESHOLD_TRANSFORM = if (spec$family == "pcm") paste0("matrix[", J, ", K - 1] tau;\n  for (j in 1:", J,
      ") tau[j] = tau_raw[j] - mean(tau_raw[j]);") else "vector[K - 1] tau = tau_raw - mean(tau_raw);"))
  fields <- c("int<lower=0, upper=1> estimate_mu;", "int<lower=0, upper=1> estimate_sigma;",
    "real ability_mean_location;", "real<lower=0> ability_mean_scale;",
    "real ability_sd_location;", "real<lower=0> ability_sd_scale;",
    "matrix[K - 1, K - 2] B_tau;", "real<lower=0> prior_sd_tau;")
  for (f in names(spec$gaussian_model$facet_sd)) fields <- c(fields,
    paste0("matrix[J_", f, ", J_", f, " - 1] B_", f, ";"), paste0("real<lower=0> prior_sd_", f, ";"))
  list(MODEL_DATA = paste(fields, collapse = "\n  "),
    ABILITY_PARAM_DECL = paste(c("vector[J_person] theta_raw;", "vector[estimate_mu] mu_theta_raw;",
      "vector[estimate_sigma] log_sigma_theta_raw;"), collapse = "\n  "),
    ABILITY_TRANSFORM = paste(c(
      "real mu_theta = estimate_mu ? ability_mean_location + ability_mean_scale * mu_theta_raw[1] : ability_mean_location;",
      "real<lower=0> sigma_theta = estimate_sigma ? exp(ability_sd_location + ability_sd_scale * log_sigma_theta_raw[1]) : ability_sd_location;",
      "vector[J_person] theta = mu_theta + sigma_theta * theta_raw;"), collapse = "\n  "),
    THRESHOLD_PARAM_DECL = if (spec$family == "pcm") paste0("matrix[", J, ", K - 2] tau_raw;") else "vector[K - 2] tau_raw;",
    THRESHOLD_TRANSFORM = if (spec$family == "pcm") paste0("matrix[", J,
      ", K - 1] tau = prior_sd_tau * tau_raw * B_tau';") else "vector[K - 1] tau = prior_sd_tau * B_tau * tau_raw;")
}

#' Generate index extraction code
#' @keywords internal
generate_index_extract <- function(spec) {
  lines <- ""
  
  # Skip person as it's in template
  other_facets <- setdiff(spec$facets_main, "person")
  
  for (facet in other_facets) {
    lines <- c(lines, paste0("    int ", facet, "_n = ", facet, "_id[n];"))
  }
  
  # Note: Bias facets reuse the main facet indices (e.g., rater_n, item_n)
  # No additional index extraction needed for bias terms
  
  paste(lines, collapse = "\n")
}

#' Generate eta sum (linear predictor)
#' @keywords internal  
generate_eta_sum <- function(spec, fair = FALSE) {
  lines <- ""
  
  # Always subtract other main facets (except person which is added)
  other_facets <- setdiff(spec$facets_main, "person")
  
  for (facet in other_facets) {
    # For fair scores, exclude nuisance facets like rater
    if (fair && facet %in% c("rater", "interlocutor")) {
      next  # skip nuisance facets in fair scores
    }
    lines <- c(lines, paste0("+ (- ", facet, "[", facet, "_n])"))
  }
  
  # Add bias terms (only for non-fair scores)
  if (!fair && length(spec$facets_bias) > 0) {
    for (bias_facet in spec$facets_bias) {
      components <- strsplit(bias_facet, ":")[[1]]
      f1 <- components[1]
      f2 <- components[2] 
      bias_clean <- gsub(":", "", bias_facet)
      
      # Skip bias terms involving nuisance facets in fair scores
      if (fair && (f1 %in% c("rater", "interlocutor") || f2 %in% c("rater", "interlocutor"))) {
        next
      }
      
      # Use the main facet indices directly (e.g., rater_n, item_n)
      lines <- c(lines, paste0("+ ", bias_clean, "[", f1, "_n, ", f2, "_n]"))
    }
  }
  
  if (length(lines) == 0) {
    return("")
  }
  
  # Join lines with proper spacing
  paste(lines, collapse = "\n               ")
}

#' Get Stan template path
#'
#' @param family Model family: "rsm" (Rating Scale Model) or "pcm" (Partial Credit Model)
#'
#' @return Path to the Stan template file
#' @keywords internal
get_stan_template_path <- function(family = "rsm") {
  # Determine template filename based on family
  template_name <- switch(family,
    "rsm" = "bmfrm_rsm.stan",
    "pcm" = "bmfrm_pcm.stan",
    # Backward compatibility: treat old family names
    "rating_scale" = "bmfrm_rsm.stan",
    "partial_credit" = "bmfrm_pcm.stan",
    stop("Unknown family: ", family, ". Use 'rsm' or 'pcm'.")
  )

  # Try package installation location first
  pkg_path <- system.file("stan", template_name, package = "BayesMFRM")
  if (nzchar(pkg_path) && file.exists(pkg_path)) {
    return(pkg_path)
  }

  # Try inst/ directory for development (from package root)
  inst_path <- file.path("inst", "stan", template_name)
  if (file.exists(inst_path)) {
    return(inst_path)
  }

  # Try BayesMFRM/inst/ directory for development (from workspace root)
  workspace_path <- file.path("BayesMFRM", "inst", "stan", template_name)
  if (file.exists(workspace_path)) {
    return(workspace_path)
  }

  stop("Stan template '", template_name, "' not found. Please check BayesMFRM installation.")
}
