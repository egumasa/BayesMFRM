#' Prepare data for bmfrm Stan model
#'
#' @param spec bmfrm_spec object from parse_bmfrm_formula()
#' @param data data.frame with original data
#' @param priors Optional prior specifications; otherwise use priors stored in
#'   spec, or the defaults. Use the same priors for build_stan_code().
#' @param obs_id Optional observation identifier preserved from before missing-score
#'   filtering for model-comparison alignment: either a single column name in
#'   `data` or a character/integer vector of length `nrow(data)` with unique
#'   non-missing values on retained rows. When `NULL` (default), preserves
#'   non-default unique `rownames(data)` or original pre-filtering row indices.
#'
#' @return list with stan_data, data_clean, facet_labels, bias_labels, and
#'   pre-preprocessing observation identity metadata (`obs_id`, `obs_id_source`,
#'   `orig_row_index`)
#' @export
prepare_data_bmfrm <- function(spec, data, priors = NULL, obs_id = NULL) {
  spec <- resolve_bmfrm_spec(spec, priors)
  
  # Validate spec
  validate_spec(spec, data)
  
  # 1. Capture pre-preprocessing observation identity BEFORE dropping missing scores
  n_before <- nrow(data)
  keep_mask <- !is.na(data[[spec$response]])
  orig_row_index <- which(keep_mask)

  if (!is.null(obs_id)) {
    if (is.character(obs_id) && length(obs_id) == 1L && obs_id %in% names(data)) {
      raw_obs_id <- as.character(data[[obs_id]])
      obs_id_source <- "explicit_column"
    } else if (length(obs_id) == n_before) {
      raw_obs_id <- as.character(obs_id)
      obs_id_source <- "explicit_vector"
    } else {
      stop("obs_id must be either a column name in data or a vector of length nrow(data)")
    }
    obs_id_clean <- raw_obs_id[keep_mask]
    if (anyNA(obs_id_clean) || any(!nzchar(obs_id_clean))) {
      stop("obs_id cannot contain NA or empty values on retained observations")
    }
    if (anyDuplicated(obs_id_clean)) {
      stop("obs_id must be unique across retained observations")
    }
  } else if (is.character(attr(data, "bmfrm_obs_id")) &&
             length(attr(data, "bmfrm_obs_id")) == n_before) {
    obs_id_clean <- as.character(attr(data, "bmfrm_obs_id"))[keep_mask]
    if (anyNA(obs_id_clean) || anyDuplicated(obs_id_clean)) {
      stop("attr(data, 'bmfrm_obs_id') must be unique and non-missing on retained observations")
    }
    obs_id_source <- "attribute"
  } else {
    rn <- rownames(data)
    has_custom_rn <- !is.null(rn) &&
      length(rn) == n_before &&
      !anyNA(rn) &&
      !anyDuplicated(rn) &&
      !identical(rn, as.character(seq_len(n_before)))
    if (has_custom_rn) {
      obs_id_clean <- rn[keep_mask]
      obs_id_source <- "rownames"
    } else {
      obs_id_clean <- as.character(orig_row_index)
      obs_id_source <- "pre_na_index"
    }
  }

  # Drop missing scores with message
  df <- data[keep_mask, , drop = FALSE]
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
    fac_factor <- if (is.factor(fac_vec)) fac_vec else factor(fac_vec)
    idx        <- as.integer(fac_factor)          # 1...J_f
    idx_list[[f]]   <- idx
    J_list[[f]]     <- nlevels(fac_factor)
    label_map[[f]]  <- levels(fac_factor)         # store original labels
  }
  
  # 4. Handle bias facets - store label metadata (Stan reuses main facet IDs and J_*)
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
      
      bias_label_map[[bias_facet]] <- list(
        facet1 = list(name = facet1, labels = label_map[[facet1]]),
        facet2 = list(name = facet2, labels = label_map[[facet2]])
      )
    }
  }
  
  response <- df[[spec$response]]
  validate_category_count(spec$K)
  if (!is.numeric(response) || any(!is.finite(response)) || any(response != floor(response))) {
    stop("Response must contain finite integer scores")
  }

  # 5. Build stan_data list
  stan_data <- list(
    N = nrow(df),
    K = spec$K,
    X = as.integer(df[[spec$response]])
  )
  
  # Add main facet data
  for (f in spec$facets_main) {
    stan_data[[paste0(f, "_id")]] <- idx_list[[f]]
    stan_data[[paste0("J_", f)]] <- J_list[[f]]
  }
  if (!is.null(spec$gaussian_model)) {
    stan_data <- c(stan_data, gaussian_model_data(spec$gaussian_model, label_map, spec$K))
  }
  
  # 6. Validate response values are in expected range
  score_range <- range(stan_data$X)
  if (score_range[1] < 1 || score_range[2] > spec$K) {
    stop("Response values must be between 1 and ", spec$K, 
         ", but found range ", score_range[1], " to ", score_range[2])
  }
  
  # 7. Return structured result
  list(
    stan_data      = stan_data,
    data_clean     = df,
    facet_labels   = label_map,
    bias_labels    = bias_label_map,
    obs_id         = obs_id_clean,
    obs_id_source  = obs_id_source,
    orig_row_index = orig_row_index
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
  facet <- match_interaction_names(facet, names(data_info$bias_labels))
  if (facet %in% names(data_info$facet_labels)) {
    return(data_info$facet_labels[[facet]])
  } else if (facet %in% names(data_info$bias_labels)) {
    # For bias facets, return combined labels
    bias_info <- data_info$bias_labels[[facet]]
    cells <- expand.grid(first = bias_info$facet1$labels, second = bias_info$facet2$labels, stringsAsFactors = FALSE)
    return(paste(bias_info$facet1$name, cells$first,
                 bias_info$facet2$name, cells$second,
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
  design_complete <- prod(vapply(facets, function(f) length(unique(data[[f]])), integer(1)))
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
