#' Compute residual correlations across items with explicit observation pairing
#'
#' These are descriptive Pearson correlations, not calibrated tests of local
#' independence. Repeated ratings require an explicit pairing or averaging rule.
#' @param resid_tbl Residual table.
#' @param person_facet Person column.
#' @param item_facet Item column.
#' @param unit_facets Columns identifying a comparable unit, e.g. person and rater.
#' @param aggregate How to handle repeated unit-item cells: error or mean.
#' @return A residual_cor object including pair-specific observation counts.
#' @export
residual_cor_criteria <- function(resid_tbl, person_facet = "person", item_facet = "item",
                                  unit_facets = person_facet, aggregate = c("error", "mean")) {
  aggregate <- match.arg(aggregate)
  needed <- unique(c(unit_facets, item_facet, "resid"))
  if (!length(unit_facets) || item_facet %in% unit_facets) stop("unit_facets must identify units separately from items")
  if (!all(needed %in% names(resid_tbl))) stop("Missing required columns: ", paste(setdiff(needed, names(resid_tbl)), collapse = ", "))
  dat <- resid_tbl[, needed, drop = FALSE]
  keys <- c(unit_facets, item_facet)
  if (anyNA(dat[keys])) stop("Missing values in residual pairing columns")
  if (!is.numeric(dat$resid)) stop("resid must be numeric")
  if (anyDuplicated(dat[keys])) {
    if (aggregate == "error") stop("Repeated unit-item ratings: specify unit_facets (e.g. person and rater) or aggregate = 'mean'")
    dat <- dat %>% dplyr::group_by(dplyr::across(dplyr::all_of(keys))) %>%
      dplyr::summarise(resid = mean(resid, na.rm = TRUE), .groups = "drop")
  }
  wide <- tidyr::pivot_wider(dat, id_cols = dplyr::all_of(unit_facets),
                            names_from = dplyr::all_of(item_facet), values_from = "resid",
                            names_prefix = "item_")
  values <- as.matrix(wide[, setdiff(names(wide), unit_facets), drop = FALSE])
  if (ncol(values) < 2L || nrow(values) < 2L) stop("Residual correlations require at least two items and two units")
  cors <- stats::cor(values, use = "pairwise.complete.obs")
  counts <- crossprod(!is.na(values))
  labels <- sub("^item_", "", colnames(values))
  dimnames(cors) <- dimnames(counts) <- list(labels, labels)
  tidy <- as.data.frame(cors) %>% tibble::rownames_to_column("item1") %>%
    tidyr::pivot_longer(cols = -item1, names_to = "item2", values_to = "correlation") %>%
    dplyr::filter(item1 != item2)
  tidy$n_pairs <- counts[cbind(match(tidy$item1, labels), match(tidy$item2, labels))]
  structure(list(cor_matrix = cors, cor_tidy = tidy, n_pairs = nrow(values),
                 pair_counts = counts, unit_facets = unit_facets, aggregation = aggregate,
                 method = "Descriptive Pearson correlation of residuals"), class = "residual_cor")
}

#' Print method for residual correlation results
#' @param x residual_cor object
#' @param ... unused
#' @export
print.residual_cor <- function(x, ...) {
  cat("Residual Correlations Across Criteria\n")
  cat("====================================\n\n")
  
  cat("Method:", x$method, "\n")
  cat("Number of unit-level observations:", x$n_pairs, "\n\n")
  
  cat("Correlation Matrix:\n")
  print(round(x$cor_matrix, 3))
  
  cat("\nSummary of correlations:\n")
  cors <- x$cor_tidy$correlation
  cors <- cors[!is.na(cors)]
  
  if (length(cors) > 0) {
    cat("  Mean absolute correlation:", round(mean(abs(cors)), 3), "\n")
    cat("  Range:", round(range(cors), 3), "\n")
    cat("  SD:", round(stats::sd(cors), 3), "\n")
    
    # Flag high correlations
    high_cors <- cors[abs(cors) > 0.3]
    if (length(high_cors) > 0) {
      cat("\n  Note:", length(high_cors), "correlations > |0.3| detected\n")
      cat("  This may indicate multidimensionality or local dependence\n")
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
    sd_cor = stats::sd(cors),
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
