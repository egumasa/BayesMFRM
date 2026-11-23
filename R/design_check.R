#' Check a planned or observed RSM/PCM design, including two-way interactions
#'
#' Checks the facet incidence graph and the rank of the linear predictor after
#' imposing the same sum-to-zero constraints as the generated model. Connectivity
#' alone does not establish identification. This is a numerical structural rank
#' check, not a test for separation, practical precision, or adequate categories.
#' Unused levels of factor columns are retained as planned levels.
#'
#' @param formula Model formula, e.g. score ~ person + item + rater.
#' @param data Planned rating rows; the response column is optional.
#' @param K Number of categories, at least two.
#' @param family "rsm" or "pcm".
#' @param observed_only Drop rows with missing responses before checking.
#' @param tol Relative singular-value tolerance.
#' @param max_cells Limit for the dense matrix used by this initial implementation.
#' @param threshold_facet Required for PCM: the facet with within-level centered steps.
#' @return A bmfrm_design_check with component membership, rank, null spaces,
#'   the constraint map, and observation counts. rank/null_space describe the
#'   category-independent predictor. Adjacent-logit rank also includes step
#'   contrasts: K-2 for RSM, or K-2 for each observed threshold-facet level in PCM.
#'   Unobserved PCM step levels are reported separately; no large Kronecker matrix
#'   is constructed. This is a property of all adjacent logits, not just observed
#'   response categories. condition_number concerns the predictor only.
#'   Interaction contrasts are included in rank/null_space. interaction_levels
#'   retains every planned cell, its observation count and constraint-map row.
#'   Sum-to-zero interaction constraints are unweighted across declared levels.
#' @export
check_design_bmfrm <- function(formula, data, K, family = "rsm",
                               observed_only = FALSE, tol = 1e-10,
                               max_cells = 2e6, threshold_facet = NULL) {
  family <- match.arg(family, c("rsm", "pcm"))
  validate_category_count(K)
  if (!is.data.frame(data) || !nrow(data)) stop("data must contain planned rating rows")
  if (length(tol) != 1L || !is.finite(tol) || tol <= 0) stop("tol must be positive")
  if (length(max_cells) != 1L || !is.finite(max_cells) || max_cells <= 0) stop("max_cells must be positive")
  vars <- all.vars(formula)
  response <- vars[1L]
  terms <- attr(stats::terms(formula), "term.labels")
  interactions <- terms[grepl(":", terms, fixed = TRUE)]
  facets <- terms[!grepl(":", terms, fixed = TRUE)]
  validate_interaction_terms(facets, interactions)
  if (any(!facets %in% names(data)) || !"person" %in% facets) {
    stop("Design checking requires person and named main facet columns")
  }
  if (family == "pcm" && (!is.character(threshold_facet) || length(threshold_facet) != 1L ||
      is.na(threshold_facet) || !threshold_facet %in% facets)) stop("PCM requires threshold_facet among the main facets")
  n_planned <- nrow(data)
  # Keep explicit factor levels even when observed-only filtering removes a level.
  fs <- lapply(data[facets], function(x) if (is.factor(x)) x else factor(x))
  if (any(vapply(fs, anyNA, logical(1)))) stop("Missing values in design facets")
  keep <- rep(TRUE, n_planned)
  if (observed_only) {
    if (!response %in% names(data)) stop("observed_only requires a response column")
    keep <- !is.na(data[[response]])
  }
  if (!any(keep)) stop("No observed rating rows remain")
  labels <- lapply(fs, levels)
  ids <- lapply(fs, function(x) as.integer(x)[keep])
  sizes <- lengths(labels)
  n <- sum(keep)
  interaction_dims <- lapply(strsplit(interactions, ":", fixed = TRUE), function(f) sizes[f])
  main_p <- sum(sizes - (facets != "person"))
  p <- main_p + sum(vapply(interaction_dims, function(n) prod(n - 1L), numeric(1)))
  if (n * p > max_cells || p * p > max_cells) {
    stop("Design exceeds the dense rank-check limit; increase max_cells deliberately or use a smaller design")
  }
  bases <- lapply(seq_along(facets), function(j) {
    if (facets[j] == "person") return(diag(sizes[j]))
    sum_zero_basis(sizes[j])
  })
  blocks <- lapply(seq_along(facets), function(j) {
    b <- bases[[j]][ids[[j]], , drop = FALSE]
    if (facets[j] != "person") b <- -b
    colnames(b) <- if (!ncol(b)) NULL else if (facets[j] == "person") paste0("theta[", labels[[j]], "]") else
      paste0(facets[j], "_contrast[", seq_len(ncol(b)), "]")
    b
  })
  X <- do.call(cbind, blocks)
  interaction_maps <- interaction_levels <- list()
  for (j in seq_along(interactions)) {
    term <- interactions[j]; f <- strsplit(term, ":", fixed = TRUE)[[1L]]
    pair_n <- interaction_dims[[j]]
    # Column-major order agrees with R/Stan matrix vectorization.
    map <- kronecker(sum_zero_basis(pair_n[2L]), sum_zero_basis(pair_n[1L]))
    index <- ids[[f[1L]]] + (ids[[f[2L]]] - 1L) * pair_n[1L]
    block <- map[index, , drop = FALSE]
    colnames(block) <- if (ncol(block)) paste0(term, "_contrast[", seq_len(ncol(block)), "]") else NULL
    X <- cbind(X, block)
    interaction_maps[[term]] <- map
    cells <- expand.grid(index1 = seq_len(pair_n[1L]), index2 = seq_len(pair_n[2L]))
    cells$level1 <- labels[[f[1L]]][cells$index1]; cells$level2 <- labels[[f[2L]]][cells$index2]
    cells$n_observations <- tabulate(index, prod(pair_n))
    interaction_levels[[term]] <- cells
  }
  offset <- c(0L, cumsum(sizes))
  nodes <- data.frame(facet = rep(facets, sizes), label = unlist(labels, use.names = FALSE),
                      n_observations = unlist(lapply(seq_along(facets), function(j) tabulate(ids[[j]], sizes[j]))),
                      stringsAsFactors = FALSE)
  parent <- seq_len(nrow(nodes))
  root <- function(k) {
    while (parent[k] != k) k <- parent[k]
    k
  }
  for (row in seq_len(n)) {
    anchor <- root(offset[1L] + ids[[1L]][row])
    for (j in seq_along(facets)[-1L]) {
      other <- root(offset[j] + ids[[j]][row])
      if (anchor != other) parent[other] <- anchor
    }
  }
  roots <- vapply(seq_along(parent), root, integer(1))
  nodes$component <- match(roots, unique(roots))
  s <- svd(X, nu = 0L, nv = p)
  rank <- sum(s$d > tol * max(s$d))
  null <- s$v[, if (rank < p) (rank + 1L):p else integer(), drop = FALSE]
  rownames(null) <- colnames(X)
  if (ncol(null)) colnames(null) <- paste0("direction_", seq_len(ncol(null)))
  effect_names <- c(paste(nodes$facet, nodes$label, sep = ":"), unlist(lapply(names(interaction_levels), function(term) {
    x <- interaction_levels[[term]]; paste0(term, "[", x$index1, ",", x$index2, "]")
  }), use.names = FALSE))
  if (length(effect_names) * p > max_cells) stop("Constraint map exceeds the dense rank-check limit")
  effect_map <- matrix(0, length(effect_names), p, dimnames = list(effect_names, colnames(X)))
  start <- 0L
  for (j in seq_along(facets)) {
    q <- ncol(bases[[j]])
    if (q) effect_map[offset[j] + seq_len(sizes[j]), start + seq_len(q)] <- bases[[j]]
    start <- start + q
  }
  row_start <- nrow(nodes)
  for (term in interactions) {
    map <- interaction_maps[[term]]
    rows <- row_start + seq_len(nrow(map)); cols <- start + seq_len(ncol(map))
    effect_map[rows, cols] <- map
    interaction_levels[[term]]$map_row <- rows
    row_start <- row_start + nrow(map); start <- start + ncol(map)
  }
  step_levels <- if (family == "pcm") labels[[threshold_facet]] else "shared"
  step_observed <- if (family == "pcm") tabulate(ids[[threshold_facet]], length(step_levels)) > 0L else TRUE
  step_parameters <- length(step_levels) * (K - 2L)
  step_rank <- sum(step_observed) * (K - 2L)
  structure(list(
    formula = formula, family = family, K = K, threshold_facet = threshold_facet, observed_only = observed_only,
    n_observations = n, n_dropped = n_planned - n,
    levels = nodes, interaction_levels = interaction_levels,
    main_parameters = main_p, interaction_parameters = p - main_p,
    n_components = length(unique(roots)), connected = length(unique(roots)) == 1L,
    rank = rank, n_parameters = p, nullity = p - rank, identifiable = rank == p && step_rank == step_parameters,
    step_rank = step_rank, step_parameters = step_parameters, step_nullity = step_parameters - step_rank,
    unobserved_step_levels = step_levels[!step_observed],
    adjacent_logit_rank = rank + step_rank, adjacent_logit_parameters = p + step_parameters,
    condition_number = if (rank == p) max(s$d) / min(s$d) else Inf,
    singular_values = s$d, tolerance = tol, design_matrix = X,
    constraint_map = effect_map, null_space = null, effect_null_space = effect_map %*% null
  ), class = "bmfrm_design_check")
}

sum_zero_basis <- function(n) {
  if (n == 1L) return(matrix(numeric(), 1L, 0L))
  b <- stats::contr.helmert(n)
  sweep(b, 2L, sqrt(colSums(b^2)), "/")
}

validate_category_count <- function(K) {
  if (!is.numeric(K) || length(K) != 1L || !is.finite(K) || K < 2 || K != floor(K)) {
    stop("K must be one integer of at least 2")
  }
  invisible(TRUE)
}

#' @export
print.bmfrm_design_check <- function(x, ...) {
  cat(toupper(x$family), "design check\n")
  cat("Rows:", x$n_observations, "  Dropped:", x$n_dropped, "\n")
  cat("Incidence components:", x$n_components, "\n")
  cat("Constrained predictor rank:", x$rank, "/", x$n_parameters, "\n")
  cat("Unidentified directions:", x$nullity, "\n")
  cat("Adjacent-logit rank including steps:", x$adjacent_logit_rank, "/", x$adjacent_logit_parameters, "\n")
  cat("Likelihood identification:", if (x$identifiable) "full structural rank" else "rank deficient", "\n")
  cat("This does not certify practical precision or finite maximum-likelihood estimates.\n")
  invisible(x)
}
