#' Prepare data for bmfrm Stan model
#'
#' @param spec bmfrm_spec object from parse_bmfrm_formula()
#' @param data data.frame with original data
#'
#' @return list with stan_data, data_clean, and facet_labels
#' @export
prepare_data_bmfrm <- function(spec, data) {
  
  # Validate spec
  validate_spec(spec, data)
  
  # 1. Drop missing scores with message
  n_before <- nrow(data)
  df <- data[!is.na(data[[spec$response]]), , drop = FALSE]
  n_after <- nrow(df)
  
  if (n_after < n_before) {
    message("Dropped ", n_before - n_after, " observations with missing scores")
  }
  
  if (n_after == 0) {
    stop("No valid observations remaining after removing missing scores")
  }
  
  # 2. Check for missing values in facet columns
  all_facet_vars <- unique(c(spec$facets_main, 
                            unlist(strsplit(spec$facets_bias, ":"))))
  
  for (var in all_facet_vars) {
    if (any(is.na(df[[var]]))) {
      stop("Missing values found in facet '", var, "'. Please clean your data first.")
    }
  }
  
  # 3. Create integer indices for each main facet
  idx_list   <- list()
  J_list     <- list()
  label_map  <- list()
  
  for (f in spec$facets_main) {
    fac_vec    <- df[[f]]
    fac_factor <- factor(fac_vec)
    idx        <- as.integer(fac_factor)          # 1...J_f
    idx_list[[f]]   <- idx
    J_list[[f]]     <- nlevels(fac_factor)
    label_map[[f]]  <- levels(fac_factor)         # store original labels
  }
  
  # 4. Handle bias facets - create combined indices
  bias_idx_list <- list()
  bias_J_list <- list()
  bias_label_map <- list()
  
  if (length(spec$facets_bias) > 0) {
    for (bias_facet in spec$facets_bias) {
      # Parse interaction (e.g., "rater:item" -> c("rater", "item"))
      components <- strsplit(bias_facet, ":")[[1]]
      if (length(components) != 2) {
        stop("Only 2-way interactions supported for bias facets: ", bias_facet)
      }
      
      facet1 <- components[1]
      facet2 <- components[2]
      
      # Get indices for each component
      idx1 <- idx_list[[facet1]]
      idx2 <- idx_list[[facet2]]
      
      bias_idx_list[[paste0(facet1, "_", bias_facet)]] <- idx1
      bias_idx_list[[paste0(facet2, "_", bias_facet)]] <- idx2
      
      bias_J_list[[paste0("J_", facet1, "_", gsub(":", "", bias_facet))]] <- J_list[[facet1]]
      bias_J_list[[paste0("J_", facet2, "_", gsub(":", "", bias_facet))]] <- J_list[[facet2]]
      
      bias_label_map[[bias_facet]] <- list(
        facet1 = list(name = facet1, labels = label_map[[facet1]]),
        facet2 = list(name = facet2, labels = label_map[[facet2]])
      )
    }
  }
  
  # 5. Build stan_data list
  stan_data <- list(
    N = nrow(df),
    K = spec$K,
    X = as.integer(df[[spec$response]])
  )
  
  # Add main facet data
  for (f in spec$facets_main) {
    stan_data[[f]]               <- idx_list[[f]]
    stan_data[[paste0("J_", f)]] <- J_list[[f]]
  }
  
  # Add bias facet data
  for (name in names(bias_idx_list)) {
    stan_data[[name]] <- bias_idx_list[[name]]
  }
  for (name in names(bias_J_list)) {
    stan_data[[name]] <- bias_J_list[[name]]
  }
  
  # 6. Validate response values are in expected range
  score_range <- range(stan_data$X)
  if (score_range[1] < 1 || score_range[2] > spec$K) {
    stop("Response values must be between 1 and ", spec$K, 
         ", but found range ", score_range[1], " to ", score_range[2])
  }
  
  # 7. Return structured result
  list(
    stan_data    = stan_data,
    data_clean   = df,
    facet_labels = label_map,
    bias_labels  = bias_label_map
  )
}

#' Get facet labels from data_info
#'
#' @param data_info output from prepare_data_bmfrm()
#' @param facet facet name
#'
#' @return character vector of labels
#' @export
get_facet_labels <- function(data_info, facet) {
  if (facet %in% names(data_info$facet_labels)) {
    return(data_info$facet_labels[[facet]])
  } else if (facet %in% names(data_info$bias_labels)) {
    # For bias facets, return combined labels
    bias_info <- data_info$bias_labels[[facet]]
    return(paste(bias_info$facet1$name, bias_info$facet1$labels, 
                 bias_info$facet2$name, bias_info$facet2$labels, 
                 sep = ":"))
  } else {
    stop("Facet '", facet, "' not found in data_info")
  }
}

#' Check data suitability for MFRM
#'
#' @param data data.frame
#' @param response response column name
#' @param facets facet column names
#'
#' @return list of diagnostic information
#' @export
check_data_suitability <- function(data, response, facets) {
  results <- list()
  
  # Basic data info
  results$n_obs <- nrow(data)
  results$response_range <- range(data[[response]], na.rm = TRUE)
  results$missing_response <- sum(is.na(data[[response]]))
  
  # Facet information
  facet_info <- list()
  for (facet in facets) {
    facet_info[[facet]] <- list(
      n_levels = length(unique(data[[facet]])),
      missing = sum(is.na(data[[facet]])),
      min_obs_per_level = min(table(data[[facet]]))
    )
  }
  results$facet_info <- facet_info
  
  # Design matrix sparsity
  design_complete <- nrow(expand.grid(lapply(facets, function(f) unique(data[[f]]))))
  design_observed <- nrow(unique(data[facets]))
  results$design_sparsity <- 1 - (design_observed / design_complete)
  
  # Warnings
  warnings <- character(0)
  if (results$missing_response > 0) {
    warnings <- c(warnings, paste("Missing responses:", results$missing_response))
  }
  
  for (facet in facets) {
    if (facet_info[[facet]]$missing > 0) {
      warnings <- c(warnings, paste("Missing values in", facet, ":", facet_info[[facet]]$missing))
    }
    if (facet_info[[facet]]$min_obs_per_level < 3) {
      warnings <- c(warnings, paste("Some", facet, "levels have very few observations"))
    }
  }
  
  if (results$design_sparsity > 0.8) {
    warnings <- c(warnings, "Very sparse design - many facet combinations unobserved")
  }
  
  results$warnings <- warnings
  
  # Overall recommendation
  if (length(warnings) == 0) {
    results$recommendation <- "Data looks suitable for MFRM analysis"
  } else if (length(warnings) <= 2) {
    results$recommendation <- "Data mostly suitable, but check warnings"
  } else {
    results$recommendation <- "Data may not be suitable for MFRM - address warnings first"
  }
  
  class(results) <- "bmfrm_data_check"
  results
}

#' Print method for data check results
#' @param x bmfrm_data_check object
#' @param ... unused
#' @export
print.bmfrm_data_check <- function(x, ...) {
  cat("BMFRM Data Suitability Check\n")
  cat("============================\n\n")
  
  cat("Observations:", x$n_obs, "\n")
  cat("Response range:", x$response_range[1], "to", x$response_range[2], "\n")
  cat("Missing responses:", x$missing_response, "\n\n")
  
  cat("Facet Information:\n")
  for (facet in names(x$facet_info)) {
    info <- x$facet_info[[facet]]
    cat("  ", facet, ": ", info$n_levels, " levels, ", 
        info$missing, " missing, min obs per level: ", info$min_obs_per_level, "\n", sep = "")
  }
  
  cat("\nDesign sparsity:", round(x$design_sparsity, 3), "\n")
  
  if (length(x$warnings) > 0) {
    cat("\nWarnings:\n")
    for (w in x$warnings) {
      cat("  *", w, "\n")
    }
  }
  
  cat("\nRecommendation:", x$recommendation, "\n")
}
