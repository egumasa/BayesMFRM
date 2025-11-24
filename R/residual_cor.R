#' Compute residual correlations across criteria/items
#'
#' This function analyzes correlations of residuals across criteria/items,
#' which can help identify multidimensionality or other model violations.
#'
#' @param resid_tbl tibble from extract_residuals() containing residuals
#' @param person_facet name of person facet (default: "person")
#' @param item_facet name of item/criterion facet (default: "item")
#'
#' @return correlation matrix or tibble with residual correlations
#' @export
residual_cor_criteria <- function(resid_tbl, 
                                  person_facet = "person", 
                                  item_facet = "item") {
  
  # Check required columns exist
  required_cols <- c(person_facet, item_facet, "resid")
  missing_cols <- setdiff(required_cols, names(resid_tbl))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Create wide format: one row per person, one column per item
  crit_wide <- resid_tbl %>%
    dplyr::select(dplyr::all_of(c(person_facet, item_facet, "resid"))) %>%
    dplyr::distinct() %>%
    tidyr::pivot_wider(
      names_from = !!rlang::sym(item_facet),
      values_from = resid,
      names_prefix = "item_"
    )
  
  # Remove person column for correlation
  cor_data <- crit_wide %>%
    dplyr::select(-!!rlang::sym(person_facet))
  
  # Compute correlation matrix
  crit_cor <- stats::cor(cor_data, use = "pairwise.complete.obs")
  
  # Add item labels as row/column names
  item_names <- gsub("item_", "", colnames(cor_data))
  rownames(crit_cor) <- item_names
  colnames(crit_cor) <- item_names
  
  # Return both matrix and tidy format
  structure(
    list(
      cor_matrix = crit_cor,
      cor_tidy = as.data.frame(crit_cor) %>%
        tibble::rownames_to_column("item1") %>%
        tidyr::pivot_longer(
          cols = -item1,
          names_to = "item2", 
          values_to = "correlation"
        ) %>%
        dplyr::filter(item1 != item2),  # Remove diagonal
      n_pairs = nrow(cor_data),
      method = "Pearson correlation of residuals"
    ),
    class = "residual_cor"
  )
}

#' Print method for residual correlation results
#' @param x residual_cor object
#' @param ... unused
#' @export
print.residual_cor <- function(x, ...) {
  cat("Residual Correlations Across Criteria\\n")
  cat("====================================\\n\\n")
  
  cat("Method:", x$method, "\\n")
  cat("Number of person-level observations:", x$n_pairs, "\\n\\n")
  
  cat("Correlation Matrix:\\n")
  print(round(x$cor_matrix, 3))
  
  cat("\\nSummary of correlations:\\n")
  cors <- x$cor_tidy$correlation
  cors <- cors[!is.na(cors)]
  
  if (length(cors) > 0) {
    cat("  Mean absolute correlation:", round(mean(abs(cors)), 3), "\\n")
    cat("  Range:", round(range(cors), 3), "\\n")
    cat("  SD:", round(sd(cors), 3), "\\n")
    
    # Flag high correlations
    high_cors <- cors[abs(cors) > 0.3]
    if (length(high_cors) > 0) {
      cat("\\n  Note:", length(high_cors), "correlations > |0.3| detected\\n")
      cat("  This may indicate multidimensionality or local dependence\\n")
    }
  }
}

#' Plot method for residual correlations
#' @param x residual_cor object
#' @param ... unused
#' @export
plot.residual_cor <- function(x, ...) {
  if (!requireNamespace("corrplot", quietly = TRUE)) {
    stop("Package 'corrplot' is required for plotting. Install with install.packages('corrplot')")
  }
  
  corrplot::corrplot(
    x$cor_matrix,
    method = "color",
    type = "upper",
    order = "hclust",
    tl.cex = 0.8,
    tl.col = "black",
    tl.srt = 45,
    title = "Residual Correlations Across Items",
    mar = c(0, 0, 1, 0)
  )
}

#' Get residual correlation diagnostics
#'
#' @param resid_cor residual_cor object
#' @param threshold threshold for flagging high correlations
#'
#' @return list with diagnostic information
#' @export
diagnose_residual_cors <- function(resid_cor, threshold = 0.3) {
  
  cors <- resid_cor$cor_tidy$correlation
  cors <- cors[!is.na(cors)]
  
  if (length(cors) == 0) {
    return(list(
      message = "No correlations computed",
      high_cors = data.frame(),
      summary = list(n_cors = 0)
    ))
  }
  
  # Find high correlations
  high_cors <- resid_cor$cor_tidy %>%
    dplyr::filter(abs(correlation) > threshold) %>%
    dplyr::arrange(dplyr::desc(abs(correlation)))
  
  # Overall summary
  summary_stats <- list(
    n_cors = length(cors),
    mean_abs_cor = mean(abs(cors)),
    sd_cor = sd(cors),
    max_abs_cor = max(abs(cors)),
    n_high_cors = nrow(high_cors),
    prop_high_cors = nrow(high_cors) / length(cors)
  )
  
  # Interpretation
  if (summary_stats$prop_high_cors > 0.1) {
    message <- "Many high residual correlations detected. Consider multidimensional model."
  } else if (summary_stats$prop_high_cors > 0.05) {
    message <- "Some high residual correlations detected. Check for local dependence."
  } else {
    message <- "Residual correlations appear acceptable for unidimensional model."
  }
  
  list(
    message = message,
    high_cors = high_cors,
    summary = summary_stats
  )
}
