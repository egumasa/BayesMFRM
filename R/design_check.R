#' Check a planned or observed additive RSM design
#'
#' Checks the facet incidence graph and the rank of the linear predictor after
#' imposing the same sum-to-zero constraints as the generated RSM. Connectivity
#' alone does not establish identification. This is a numerical structural rank
#' check, not a test for separation, practical precision, or adequate categories.
#' Unused levels of factor columns are retained as planned levels.
#'
#' @param formula Model formula, e.g. score ~ person + item + rater.
#' @param data Planned rating rows; the response column is optional.
#' @param K Number of categories, at least two.
#' @param family Currently only additive "rsm" is supported.
#' @param observed_only Drop rows with missing responses before checking.
#' @param tol Relative singular-value tolerance.
#' @param max_cells Limit for the dense matrix used by this initial implementation.
#' @return A bmfrm_design_check with component membership, rank, null spaces,
#'   the constraint map, and observation counts. Adjacent-logit rank includes
#'   the K - 2 free common threshold coordinates.
#' @export
check_design_bmfrm <- function(formula, data, K, family = "rsm",
                               observed_only = FALSE, tol = 1e-10,
                               max_cells = 2e6) {
  validate_category_count(K)
  if (!is.data.frame(data) || !nrow(data)) stop("data must contain planned rating rows")
  if (length(tol) != 1L || !is.finite(tol) || tol <= 0) stop("tol must be positive")
  if (length(max_cells) != 1L || !is.finite(max_cells) || max_cells <= 0) stop("max_cells must be positive")
  vars <- all.vars(formula)
  response <- vars[1L]
  facets <- attr(stats::terms(formula), "term.labels")
  if (family != "rsm" || any(!facets %in% names(data)) ||
      any(grepl(":", facets)) || !"person" %in% facets) {
    stop("Design checking currently supports additive RSM formulas with person and named facet columns")
  }
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
  p <- sum(sizes - (facets != "person"))
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
  effect_map <- matrix(0, nrow(nodes), p,
                       dimnames = list(paste(nodes$facet, nodes$label, sep = ":"), colnames(X)))
  start <- 0L
  for (j in seq_along(facets)) {
    q <- ncol(bases[[j]])
    if (q) effect_map[offset[j] + seq_len(sizes[j]), start + seq_len(q)] <- bases[[j]]
    start <- start + q
  }
  structure(list(
    formula = formula, family = family, K = K, observed_only = observed_only,
    n_observations = n, n_dropped = n_planned - n,
    levels = nodes, n_components = length(unique(roots)), connected = length(unique(roots)) == 1L,
    rank = rank, n_parameters = p, nullity = p - rank, identifiable = rank == p,
    adjacent_logit_rank = rank + K - 2L, adjacent_logit_parameters = p + K - 2L,
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
  cat("Additive RSM design check\n")
  cat("Rows:", x$n_observations, "  Dropped:", x$n_dropped, "\n")
  cat("Incidence components:", x$n_components, "\n")
  cat("Constrained predictor rank:", x$rank, "/", x$n_parameters, "\n")
  cat("Unidentified directions:", x$nullity, "\n")
  cat("Likelihood identification:", if (x$identifiable) "full structural rank" else "rank deficient", "\n")
  cat("This does not certify practical precision or finite maximum-likelihood estimates.\n")
  invisible(x)
}
