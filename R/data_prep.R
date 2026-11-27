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
  if (!is.numeric(response) || any(!is.finite(response)) || any(response != floor(response))) {
    stop("Response must contain finite integer scores")
  }
  spec_eff <- spec
  is_hetero <- isTRUE(spec_eff$K_heterogeneous)

  # 5. Build stan_data list
  # Separate input-level rubric validation from post-NA-drop level mapping so
  # dropping an all-NA level from a character column never rejects valid scalar K
  # or valid pre-NA named K definitions. If all retained levels share a single K
  # (including K = 2), normalize smoothly to the uniform PCM path.
  if (is_hetero) {
    tf <- spec_eff$threshold_facet
    tf_labels <- label_map[[tf]]
    raw_K_by_level <- spec_eff$K_by_level
    if (is.null(raw_K_by_level) || is.null(names(raw_K_by_level))) {
      stop("Heterogeneous PCM specification requires named K_by_level")
    }
    missing_lv <- setdiff(tf_labels, names(raw_K_by_level))
    if (length(missing_lv) > 0L) {
      stop(
        "Named vector K is missing level(s) of threshold_facet '", tf, "': ",
        paste(missing_lv, collapse = ", ")
      )
    }
    input_tf_vec <- data[[tf]]
    input_tf_levels <- if (is.factor(input_tf_vec)) {
      levels(input_tf_vec)
    } else {
      unique(as.character(input_tf_vec[!is.na(input_tf_vec)]))
    }
    unknown_lv <- setdiff(names(raw_K_by_level), union(tf_labels, input_tf_levels))
    if (length(unknown_lv) > 0L) {
      stop(
        "Named vector K contains unknown level(s) not present in threshold_facet '", tf, "': ",
        paste(unknown_lv, collapse = ", ")
      )
    }
    K_level_named <- stats::setNames(as.integer(raw_K_by_level[tf_labels]), tf_labels)
    obs_tf <- as.character(df[[tf]])
    obs_K <- unname(K_level_named[obs_tf])
    bad <- response < 1 | response > obs_K
    if (any(bad)) {
      bad_lv <- unique(obs_tf[bad])
      stop(
        "Response values exceed valid category range 1..K_j for threshold_facet '",
        tf, "' level(s): ",
        paste(sprintf("%s (1..%d)", bad_lv, K_level_named[bad_lv]), collapse = ", ")
      )
    }
    K_level_int <- as.integer(unname(K_level_named))
    K_max_int <- as.integer(max(K_level_int))
    if (length(unique(K_level_int)) == 1L) {
      # All surviving levels share a single category count (e.g., K = 2 only):
      # normalize internal specification and Stan data to the uniform PCM path.
      is_hetero <- FALSE
      spec_eff$K <- K_max_int
      spec_eff$K_by_level <- K_level_named
      spec_eff$K_max <- K_max_int
      spec_eff$K_heterogeneous <- FALSE
      if (!is.null(spec_eff$gaussian_model)) {
        spec_eff <- resolve_bmfrm_spec(spec_eff, spec_eff$priors)
      }
      stan_data <- list(
        N = nrow(df),
        K = K_max_int,
        X = as.integer(response)
      )
    } else {
      tau_len_int <- as.integer(K_level_int - 1L)
      tau_start_int <- as.integer(c(1L, utils::head(cumsum(tau_len_int), -1L) + 1L))
      stan_data <- list(
        N = nrow(df),
        K_max = K_max_int,
        X = as.integer(response),
        K_level = K_level_int,
        N_tau_steps = as.integer(sum(tau_len_int)),
        tau_start = tau_start_int,
        tau_len = tau_len_int
      )
    }
  } else {
    K_scalar <- if (length(spec_eff$K) == 1L) spec_eff$K else as.integer(get_spec_K_max(spec_eff))
    validate_category_count(K_scalar)
    K_max_int <- as.integer(K_scalar)
    stan_data <- list(
      N = nrow(df),
      K = K_scalar,
      X = as.integer(response)
    )
    K_level_named <- if (identical(spec_eff$family, "pcm") && !is.null(spec_eff$threshold_facet)) {
      tf_labels <- label_map[[spec_eff$threshold_facet]]
      stats::setNames(rep(K_max_int, length(tf_labels)), tf_labels)
    } else {
      NULL
    }
  }
  
  # Add main facet data
  for (f in spec_eff$facets_main) {
    stan_data[[paste0(f, "_id")]] <- idx_list[[f]]
    stan_data[[paste0("J_", f)]] <- J_list[[f]]
  }
  if (!is.null(spec_eff$gaussian_model)) {
    stan_data <- c(
      stan_data,
      gaussian_model_data(spec_eff$gaussian_model, label_map, K_max_int, K_by_level = K_level_named)
    )
  } else if (is_hetero) {
    # Legacy/no-gaussian_model raw ragged step offsets
    K_level_int <- as.integer(unname(K_level_named))
    tau_raw_len <- ifelse(K_level_int >= 3L, K_level_int - 1L, 0L)
    tau_raw_start <- as.integer(cumsum(c(1L, utils::head(tau_raw_len, -1L))))
    tau_raw_start[tau_raw_len == 0L] <- 1L
    stan_data$N_tau_raw <- as.integer(sum(tau_raw_len))
    stan_data$tau_raw_start <- as.integer(tau_raw_start)
    stan_data$tau_raw_len <- as.integer(tau_raw_len)
  }
  
  # 6. Validate response values are in expected range (already checked per level above for hetero)
  if (!is_hetero) {
    score_range <- range(stan_data$X)
    if (score_range[1] < 1 || score_range[2] > K_max_int) {
      stop("Response values must be between 1 and ", K_max_int, 
           ", but found range ", score_range[1], " to ", score_range[2])
    }
  }
  
  # 7. Return structured result
  list(
    stan_data       = stan_data,
    data_clean      = df,
    facet_labels    = label_map,
    bias_labels     = bias_label_map,
    spec            = spec_eff,
    K_by_level      = K_level_named,
    K_max           = K_max_int,
    K_heterogeneous = is_hetero,
    obs_id          = obs_id_clean,
    obs_id_source   = obs_id_source,
    orig_row_index  = orig_row_index
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
  if (!is.data.frame(data) || nrow(data) == 0L) {
    stop("data must be a non-empty data.frame")
  }
  if (!is.character(response) || length(response) != 1L || is.na(response) ||
      !(response %in% names(data))) {
    stop("Response column '", response, "' not found in data")
  }
  if (!is.character(facets) || length(facets) == 0L || anyNA(facets)) {
    stop("facets must be a non-empty character vector of column names")
  }
  missing_cols <- setdiff(facets, names(data))
  if (length(missing_cols) > 0L) {
    stop("Facet column(s) not found in data: ", paste(missing_cols, collapse = ", "))
  }

  results <- list()
  resp_vec <- data[[response]]
  non_na_resp <- !is.na(resp_vec)
  unsuitable <- FALSE
  warnings <- character(0)

  # Basic data info and response validity checks
  results$n_obs <- nrow(data)
  results$missing_response <- sum(!non_na_resp)
  if (!is.numeric(resp_vec)) {
    unsuitable <- TRUE
    results$response_range <- c(NA_real_, NA_real_)
    warnings <- c(warnings, paste0("Response '", response, "' must be numeric integer category scores (>= 1)"))
  } else if (!any(non_na_resp)) {
    unsuitable <- TRUE
    results$response_range <- c(NA_real_, NA_real_)
    warnings <- c(warnings, paste0("Response '", response, "' has no non-missing observations"))
  } else {
    obs_r <- resp_vec[non_na_resp]
    fin_r <- obs_r[is.finite(obs_r)]
    results$response_range <- if (length(fin_r) > 0L) range(fin_r) else c(NA_real_, NA_real_)
    if (any(!is.finite(obs_r))) {
      unsuitable <- TRUE
      warnings <- c(warnings, paste0("Response '", response, "' contains non-finite values (Inf/-Inf/NaN)"))
    }
    if (length(fin_r) > 0L && any(fin_r != floor(fin_r))) {
      unsuitable <- TRUE
      warnings <- c(warnings, paste0("Response '", response, "' contains non-integer values; scores must be integers 1..K"))
    }
    if (length(fin_r) > 0L && any(fin_r < 1)) {
      unsuitable <- TRUE
      warnings <- c(warnings, paste0("Response '", response, "' contains values < 1; category scores must start at 1"))
    }
  }
  if (results$missing_response > 0) {
    warnings <- c(warnings, paste("Missing responses:", results$missing_response))
  }
  
  # Facet information (exclude NA from level count; keep declared factor levels consistent with table())
  facet_info <- list()
  for (facet in facets) {
    v <- data[[facet]]
    lv <- if (is.factor(v)) levels(v) else unique(v[!is.na(v)])
    n_lv <- length(lv)
    tab <- if (n_lv > 0L) table(factor(v[!is.na(v)], levels = lv)) else integer(0)
    facet_info[[facet]] <- list(
      n_levels = n_lv,
      missing = sum(is.na(v)),
      min_obs_per_level = if (length(tab) > 0L) as.integer(min(tab)) else 0L
    )
  }
  results$facet_info <- facet_info
  
  # Design matrix sparsity (evaluated on non-missing facet & response rows)
  design_complete <- prod(vapply(facet_info, function(fi) as.numeric(fi$n_levels), numeric(1)))
  valid_rows <- stats::complete.cases(data[facets]) & non_na_resp
  if (!any(valid_rows)) unsuitable <- TRUE
  design_observed <- if (any(valid_rows)) nrow(unique(data[valid_rows, facets, drop = FALSE])) else 0L
  results$design_sparsity <- if (design_complete > 0) 1 - (design_observed / design_complete) else NA_real_
  
  for (facet in facets) {
    if (facet_info[[facet]]$missing > 0) {
      warnings <- c(warnings, paste("Missing values in", facet, ":", facet_info[[facet]]$missing))
    }
    if (facet_info[[facet]]$n_levels == 0L) {
      unsuitable <- TRUE
      warnings <- c(warnings, paste("Facet", facet, "has no valid levels"))
    } else if (facet_info[[facet]]$min_obs_per_level < 3) {
      warnings <- c(warnings, paste("Some", facet, "levels have very few observations"))
    }
  }
  
  if (is.na(results$design_sparsity)) {
    warnings <- c(warnings, "Design sparsity cannot be computed because at least one facet has no valid levels")
  } else if (results$design_sparsity > 0.8) {
    warnings <- c(warnings, "Very sparse design - many facet combinations unobserved")
  }
  
  results$warnings <- warnings
  
  # Overall recommendation
  if (unsuitable) {
    results$recommendation <- "Data is not suitable for MFRM - resolve invalid or missing data before fitting"
  } else if (length(warnings) == 0) {
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
