#' Build Stan code from specification and priors
#'
#' @param spec bmfrm_spec object from parse_bmfrm_formula()
#' @param priors list of bmfrm_prior objects
#'
#' @return character string with complete Stan code
#' @export
build_stan_code <- function(spec, priors) {
  
  # Validate priors
  priors <- validate_priors(priors, spec$facets_main, spec$facets_bias)
  
  # Load base template
  template_path <- get_stan_template_path()
  
  stan_code <- readLines(template_path, warn = FALSE)
  stan_code <- paste(stan_code, collapse = "\n")
  
  # Generate replacements for each placeholder
  replacements <- list()
  
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
  lines <- character(0)
  
  # Main facets - declare dimensions first, then arrays
  for (facet in spec$facets_main) {
    lines <- c(lines, 
               paste0("  int<lower=1> J_", facet, ";"),
               paste0("  array[N] int<lower=1, upper=J_", facet, "> ", facet, ";"))
  }
  
  # Note: Bias facets don't need separate index arrays - they reuse main facet indices
  # The bias parameters are matrices indexed by the existing main facet arrays
  
  paste(lines, collapse = "\n")
}

#' Generate parameter declarations for main facets
#' @keywords internal  
generate_facet_param_decl <- function(spec) {
  lines <- character(0)
  
  # Skip person (theta) as it's always present in template
  other_facets <- setdiff(spec$facets_main, "person")
  
  for (facet in other_facets) {
    lines <- c(lines, paste0("  vector[J_", facet, "] ", facet, "_raw;"))
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
  lines <- character(0)
  
  # Skip person (theta) as it's handled in template
  other_facets <- setdiff(spec$facets_main, "person")
  
  for (facet in other_facets) {
    lines <- c(lines,
               paste0("  // ", facet, " (sum-to-zero constraint)"),
               paste0("  vector[J_", facet, "] ", facet, ";"),
               "  {",
               paste0("    real mean_", facet, " = mean(", facet, "_raw);"),
               paste0("    ", facet, " = ", facet, "_raw - mean_", facet, ";"),
               "  }")
  }
  
  paste(lines, collapse = "\n")
}

#' Generate transformed parameters for bias facets
#' @keywords internal 
generate_bias_transform <- function(spec) {
  lines <- character(0)
  
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
  lines <- character(0)
  
  for (prior_obj in priors) {
    if (prior_obj$class == "theta") {
      lines <- c(lines, paste0("  theta_raw ~ ", prior_obj$spec, ";"))
    } else if (prior_obj$class == "item") {
      lines <- c(lines, paste0("  item_raw ~ ", prior_obj$spec, ";"))
    } else if (prior_obj$class == "rater") {
      lines <- c(lines, paste0("  rater_raw ~ ", prior_obj$spec, ";"))
    } else if (prior_obj$class == "tau") {
      lines <- c(lines, paste0("  tau_raw ~ ", prior_obj$spec, ";"))
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

#' Generate index extraction code
#' @keywords internal
generate_index_extract <- function(spec) {
  lines <- character(0)
  
  # Skip person as it's in template
  other_facets <- setdiff(spec$facets_main, "person")
  
  for (facet in other_facets) {
    lines <- c(lines, paste0("    int ", facet, "_n = ", facet, "[n];"))
  }
  
  # Note: Bias facets reuse the main facet indices (e.g., rater_n, item_n)
  # No additional index extraction needed for bias terms
  
  paste(lines, collapse = "\n")
}

#' Generate eta sum (linear predictor)
#' @keywords internal  
generate_eta_sum <- function(spec, fair = FALSE) {
  lines <- character(0)
  
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
#' @keywords internal
get_stan_template_path <- function() {
  # Try package installation location first
  pkg_path <- system.file("stan", "bmfrm_template.stan", package = "BayesMFRM")
  if (file.exists(pkg_path)) {
    return(pkg_path)
  }
  
  # Try inst/ directory for development
  inst_path <- file.path("inst", "stan", "bmfrm_template.stan")
  if (file.exists(inst_path)) {
    return(inst_path)
  }
  
  stop("Stan template not found. Please check BayesMFRM installation.")
}
