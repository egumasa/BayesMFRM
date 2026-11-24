#' Check a planned or observed RSM/PCM design, including two-way interactions
#'
#' Checks the facet incidence graph and the rank of the linear predictor after
#' imposing the same sum-to-zero constraints as the generated model. Connectivity
#' alone does not establish identification. This is a numerical structural rank
#' check, not a test for separation, practical precision, or adequate categories.
#' Unused levels of factor columns are retained as planned levels.
#'
#' For small/medium designs (`n * p <= max_cells` and `p * p <= max_cells`),
#' uses dense SVD. For large designs (e.g., thousands of persons), `method = "auto"`
#' automatically switches to an exact Schur-complement (within-person Gram
#' matrix) solver that scales in `O(n * p_other + p_other^3)` without ever
#' forming an `n x P` dense matrix, while still computing the exact structural
#' rank and null space.
#'
#' @param formula Model formula, e.g. score ~ person + item + rater.
#' @param data Planned rating rows; the response column is optional.
#' @param K Number of categories, at least two.
#' @param family "rsm" or "pcm".
#' @param observed_only Drop rows with missing responses before checking.
#' @param tol Relative singular-value tolerance.
#' @param max_cells Cell limit for dense SVD (and for `p_other^2` in Schur mode).
#' @param threshold_facet Required for PCM: the facet with within-level centered steps.
#' @param method Solver mode: `"auto"` (dense SVD when within `max_cells`,
#'   Schur-complement solver otherwise), `"dense"` (strict dense SVD limit), or
#'   `"schur"` (always use within-person Schur-complement solver).
#' @return A bmfrm_design_check with component membership, rank, null spaces,
#'   the constraint map, and observation counts.
#' @export
check_design_bmfrm <- function(formula, data, K, family = "rsm",
                               observed_only = FALSE, tol = 1e-10,
                               max_cells = 2e6, threshold_facet = NULL,
                               method = c("auto", "dense", "schur")) {
  family <- match.arg(family, c("rsm", "pcm"))
  method <- match.arg(method)
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
      is.na(threshold_facet) || !threshold_facet %in% facets)) {
    stop("PCM requires threshold_facet among the main facets")
  }
  n_planned <- nrow(data)
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
  P <- sizes[["person"]]
  interaction_dims <- lapply(strsplit(interactions, ":", fixed = TRUE), function(f) sizes[f])
  main_p <- sum(sizes - (facets != "person"))
  p <- main_p + sum(vapply(interaction_dims, function(sz) prod(sz - 1L), numeric(1)))
  p_other <- p - P
  n_effects <- sum(sizes) + sum(vapply(interaction_dims, prod, numeric(1)))

  # Decide dense vs Schur-complement solver and enforce ALL memory limits BEFORE
  # allocating any contrast bases, Kronecker products, or cell tables.
  exceeds_dense <- (as.numeric(n) * as.numeric(p) > max_cells) ||
    (as.numeric(p) * as.numeric(p) > max_cells) ||
    (as.numeric(n_effects) * as.numeric(p) > max_cells)
  if (method == "dense" && exceeds_dense) {
    stop("Design exceeds the dense rank-check limit; use method = 'auto'/'schur' or increase max_cells deliberately")
  }
  use_schur <- (method == "schur") || (method == "auto" && exceeds_dense)
  if (use_schur && as.numeric(p_other) * as.numeric(p_other) > max_cells) {
    stop("Non-person parameter block (p_other = ", p_other,
         ") exceeds max_cells limit; increase max_cells deliberately")
  }

  # Incidence graph connected components (O(n * num_facets), independent of matrix size)
  offset <- c(0L, cumsum(sizes))
  nodes <- data.frame(
    facet = rep(facets, sizes),
    label = unlist(labels, use.names = FALSE),
    n_observations = unlist(lapply(seq_along(facets), function(j) tabulate(ids[[j]], sizes[j]))),
    stringsAsFactors = FALSE
  )
  parent <- seq_len(nrow(nodes))
  root <- function(k) {
    while (parent[k] != k) {
      parent[k] <<- parent[parent[k]]
      k <- parent[k]
    }
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

  # Factorized row-wise Kronecker product: evaluates rows of kronecker(B2, B1)
  # for cell indices (i1, i2) without materializing the (J1*J2) x (q1*q2) matrix.
  kron_rows <- function(B1, B2, i1, i2) {
    q1 <- ncol(B1)
    q2 <- ncol(B2)
    if (q1 == 0L || q2 == 0L) return(matrix(numeric(), nrow = length(i1), ncol = 0L))
    A1 <- B1[i1, , drop = FALSE]
    A2 <- B2[i2, , drop = FALSE]
    A1[, rep(seq_len(q1), times = q2), drop = FALSE] *
      A2[, rep(seq_len(q2), each = q1), drop = FALSE]
  }

  # Factorized Kronecker matmul: computes kronecker(B2, B1) %*% V in O(k * J1 * J2)
  # memory using vec(B1 %*% mat(v, q1, q2) %*% t(B2)) without forming kronecker(B2, B1).
  apply_kron_null <- function(B1, B2, V) {
    q1 <- ncol(B1)
    q2 <- ncol(B2)
    J1 <- nrow(B1)
    J2 <- nrow(B2)
    k_cols <- ncol(V)
    if (k_cols == 0L) return(matrix(numeric(), nrow = J1 * J2, ncol = 0L))
    out <- matrix(0, nrow = J1 * J2, ncol = k_cols)
    if (q1 == 0L || q2 == 0L) return(out)
    tB2 <- t(B2)
    for (col_j in seq_len(k_cols)) {
      out[, col_j] <- as.vector(B1 %*% matrix(V[, col_j], nrow = q1, ncol = q2) %*% tB2)
    }
    out
  }

  # Prepare compact 1D factor bases and interaction cell metadata (zero Kronecker allocation)
  bases_other <- lapply(setdiff(facets, "person"), function(f) sum_zero_basis(sizes[[f]]))
  names(bases_other) <- setdiff(facets, "person")
  person_sz_basis <- NULL
  get_sz_basis <- function(f) {
    if (f != "person") return(bases_other[[f]])
    if (is.null(person_sz_basis)) person_sz_basis <<- sum_zero_basis(P)
    person_sz_basis
  }

  interaction_factors <- interaction_levels <- list()
  for (j in seq_along(interactions)) {
    term <- interactions[j]
    f <- strsplit(term, ":", fixed = TRUE)[[1L]]
    pair_n <- interaction_dims[[j]]
    q_term <- as.integer(prod(pair_n - 1L))
    if (q_term == 0L) {
      # Skip basis allocation when either facet has <= 1 level (0 interaction DOF)
      B1 <- NULL
      B2 <- NULL
    } else {
      B1 <- get_sz_basis(f[1L])
      B2 <- get_sz_basis(f[2L])
    }
    interaction_factors[[term]] <- list(
      f1 = f[1L], f2 = f[2L], B1 = B1, B2 = B2,
      nr = as.integer(prod(pair_n)), q = q_term
    )
    index <- ids[[f[1L]]] + (ids[[f[2L]]] - 1L) * pair_n[1L]
    cells <- expand.grid(index1 = seq_len(pair_n[1L]), index2 = seq_len(pair_n[2L]))
    cells$level1 <- labels[[f[1L]]][cells$index1]
    cells$level2 <- labels[[f[2L]]][cells$index2]
    cells$n_observations <- tabulate(index, prod(pair_n))
    interaction_levels[[term]] <- cells
  }

  effect_names <- c(paste(nodes$facet, nodes$label, sep = ":"), unlist(lapply(names(interaction_levels), function(term) {
    x <- interaction_levels[[term]]
    paste0(term, "[", x$index1, ",", x$index2, "]")
  }), use.names = FALSE))

  # Assign map_row on interaction_levels consistently in both dense and schur modes
  row_start <- nrow(nodes)
  for (term in interactions) {
    nr <- interaction_factors[[term]]$nr
    interaction_levels[[term]]$map_row <- row_start + seq_len(nr)
    row_start <- row_start + nr
  }

  if (!use_schur) {
    # Dense SVD path (exact backward compatibility)
    bases <- lapply(seq_along(facets), function(j) {
      if (facets[j] == "person") return(diag(sizes[j]))
      bases_other[[facets[j]]]
    })
    blocks <- lapply(seq_along(facets), function(j) {
      b <- bases[[j]][ids[[j]], , drop = FALSE]
      if (facets[j] != "person") b <- -b
      colnames(b) <- if (!ncol(b)) NULL else if (facets[j] == "person") paste0("theta[", labels[[j]], "]") else
        paste0(facets[j], "_contrast[", seq_len(ncol(b)), "]")
      b
    })
    X <- do.call(cbind, blocks)
    for (term in interactions) {
      fac <- interaction_factors[[term]]
      if (fac$q) {
        block <- kron_rows(fac$B1, fac$B2, ids[[fac$f1]], ids[[fac$f2]])
        colnames(block) <- paste0(term, "_contrast[", seq_len(ncol(block)), "]")
        X <- cbind(X, block)
      }
    }

    s <- svd(X, nu = 0L, nv = p)
    rank <- sum(s$d > tol * max(s$d))
    null <- s$v[, if (rank < p) (rank + 1L):p else integer(), drop = FALSE]
    rownames(null) <- colnames(X)
    if (ncol(null)) colnames(null) <- paste0("direction_", seq_len(ncol(null)))

    effect_map <- matrix(0, length(effect_names), p, dimnames = list(effect_names, colnames(X)))
    start <- 0L
    for (j in seq_along(facets)) {
      q <- ncol(bases[[j]])
      if (q) effect_map[offset[j] + seq_len(sizes[j]), start + seq_len(q)] <- bases[[j]]
      start <- start + q
    }
    for (term in interactions) {
      fac <- interaction_factors[[term]]
      rows_t <- interaction_levels[[term]]$map_row
      if (fac$q) {
        cols_t <- start + seq_len(fac$q)
        effect_map[rows_t, cols_t] <- kronecker(fac$B2, fac$B1)
        start <- start + fac$q
      }
    }
    effect_null <- effect_map %*% null
    cond_num <- if (rank == p) max(s$d) / min(s$d) else Inf
    sing_vals <- s$d
    solver_used <- "dense_svd"
  } else {
    # Exact Schur-complement (within-person TSQR) solver:
    # Eliminates person indicators analytically: X = [X_P, X_O] where X_P^T X_P = diag(n_p).
    # Never materializes Kronecker interaction bases; evaluates rows and null maps in factored form.
    other_facets <- setdiff(facets, "person")
    p_id <- ids[["person"]]
    n_p <- tabulate(p_id, P)
    obs_p <- which(n_p > 0L)
    unobs_p <- which(n_p == 0L)
    rank_p <- length(obs_p)

    # Column names matching dense ordering (person first in Schur coordinate, then mapped cleanly)
    col_names_p <- paste0("theta[", labels[["person"]], "]")
    col_names_o <- unlist(c(
      lapply(other_facets, function(f) {
        q <- ncol(bases_other[[f]])
        if (q) paste0(f, "_contrast[", seq_len(q), "]") else character()
      }),
      lapply(interactions, function(term) {
        q <- interaction_factors[[term]]$q
        if (q) paste0(term, "_contrast[", seq_len(q), "]") else character()
      })
    ), use.names = FALSE)

    if (p_other == 0L) {
      rank_o <- 0L
      null_o <- matrix(numeric(), 0L, 0L)
      S_po <- matrix(0, P, 0L)
      sing_vals <- sqrt(sort(n_p, decreasing = TRUE))
    } else {
      # Numerically stable within-person TSQR (avoids forming normal equations X^T X - S^T D^-1 S
      # which would square the condition number and suffer sqrt(eps) cancellation):
      #   Pass 1: compute per-person sums S_po (P x p_other) and means M_po = diag(1/n_p) S_po
      #   Pass 2: stream centered chunks Xo_tilde = Xo - M_po[p_id, ] through incremental QR
      #   Final : SVD of the compact p_other x p_other upper-triangular factor R_acc
      build_xo_chunk <- function(r_idx) {
        Xo_blocks <- list()
        for (f in other_facets) {
          b <- bases_other[[f]]
          if (ncol(b)) Xo_blocks[[length(Xo_blocks) + 1L]] <- -b[ids[[f]][r_idx], , drop = FALSE]
        }
        for (term in interactions) {
          fac <- interaction_factors[[term]]
          if (fac$q) {
            Xo_blocks[[length(Xo_blocks) + 1L]] <- kron_rows(
              fac$B1, fac$B2, ids[[fac$f1]][r_idx], ids[[fac$f2]][r_idx]
            )
          }
        }
        do.call(cbind, Xo_blocks)
      }

      S_po <- matrix(0, P, p_other)
      G_xx <- matrix(0, p_other, p_other)
      chunk_sz <- as.integer(max(1L, min(n, floor(max_cells / max(p_other, 1L)))))
      starts <- seq.int(1L, n, by = chunk_sz)
      for (st_i in starts) {
        en_i <- min(n, st_i + chunk_sz - 1L)
        r_idx <- st_i:en_i
        Xo_c <- build_xo_chunk(r_idx)
        G_xx <- G_xx + crossprod(Xo_c)
        rs <- rowsum(Xo_c, p_id[r_idx], reorder = FALSE)
        S_po[as.integer(rownames(rs)), ] <- S_po[as.integer(rownames(rs)), , drop = FALSE] + rs
      }
      inv_np <- ifelse(n_p > 0L, 1 / n_p, 0)
      M_po <- S_po * inv_np

      # Pass 2: incremental QR reduction of within-person demeaned X_O
      R_acc <- matrix(0, 0L, p_other)
      for (st_i in starts) {
        en_i <- min(n, st_i + chunk_sz - 1L)
        r_idx <- st_i:en_i
        Xo_tilde <- build_xo_chunk(r_idx) - M_po[p_id[r_idx], , drop = FALSE]
        stacked <- rbind(R_acc, Xo_tilde)
        q_step <- qr(stacked, LAPACK = TRUE)
        R_raw <- qr.R(q_step)
        # Undo column pivoting so columns remain in original parameter order
        piv_inv <- order(q_step$pivot)
        R_acc <- R_raw[, piv_inv, drop = FALSE]
      }

      # Exact SVD of the compact R_acc (same singular values and right singular vectors as M_P X_O!)
      s_o <- svd(R_acc, nu = 0L, nv = p_other)
      scale_ref <- max(c(s_o$d, sqrt(n_p), 1))
      rank_o <- sum(s_o$d > tol * scale_ref)
      null_o <- s_o$v[, if (rank_o < p_other) (rank_o + 1L):p_other else integer(), drop = FALSE]
      sing_vals <- s_o$d
    }

    rank <- rank_p + rank_o
    n_null <- p - rank
    # Reconstruct exact null space in parameter coordinates [person, other] (and effect coordinates)
    if (n_null == 0L) {
      null <- matrix(numeric(), p, 0L, dimnames = list(c(col_names_p, col_names_o), NULL))
      effect_null <- matrix(numeric(), length(effect_names), 0L, dimnames = list(effect_names, NULL))
    } else {
      null_cols <- list()
      # 1. Unobserved persons: unit vector e_j on person j, zero elsewhere
      for (uj in unobs_p) {
        v <- numeric(p)
        v[uj] <- 1
        null_cols[[length(null_cols) + 1L]] <- v
      }
      # 2. Lifted Schur null vectors: v_P[obs] = -(1 / n_p) * S_po[obs, ] %*% v_O
      if (ncol(null_o) > 0L) {
        inv_np <- ifelse(n_p > 0L, 1 / n_p, 0)
        for (k in seq_len(ncol(null_o))) {
          v_o <- null_o[, k]
          v_p <- -inv_np * as.vector(S_po %*% v_o)
          v <- c(v_p, v_o)
          v <- v / sqrt(sum(v^2))
          null_cols[[length(null_cols) + 1L]] <- v
        }
      }
      null <- do.call(cbind, null_cols)
      rownames(null) <- c(col_names_p, col_names_o)
      colnames(null) <- paste0("direction_", seq_len(ncol(null)))

      # Map directly into effect_null_space (nrow = length(effect_names), ncol = n_null) without forming p x effect_names
      effect_null <- matrix(0, length(effect_names), n_null,
                            dimnames = list(effect_names, colnames(null)))
      p_facet_idx <- match("person", facets)
      effect_null[offset[p_facet_idx] + seq_len(P), ] <- null[seq_len(P), , drop = FALSE]
      o_pos <- P
      for (f in other_facets) {
        fi <- match(f, facets)
        b <- bases_other[[f]]
        q <- ncol(b)
        if (q) {
          effect_null[offset[fi] + seq_len(sizes[[f]]), ] <- b %*% null[o_pos + seq_len(q), , drop = FALSE]
          o_pos <- o_pos + q
        }
      }
      for (term in interactions) {
        fac <- interaction_factors[[term]]
        q <- fac$q
        if (q) {
          effect_null[interaction_levels[[term]]$map_row, ] <- apply_kron_null(
            fac$B1, fac$B2, null[o_pos + seq_len(q), , drop = FALSE]
          )
          o_pos <- o_pos + q
        }
      }
    }

    X <- NULL
    effect_map <- NULL
    # Compute the EXACT condition number kappa(X) = sigma_max(X) / sigma_min(X) of the
    # full design matrix X = [X_P, X_O] via symmetric Lanczos on G = X^T X and G^-1
    # (using exact Schur-complement matvec in O(P * p_other + p_other^2), matching dense SVD).
    if (rank < p) {
      cond_num <- Inf
    } else if (p_other == 0L) {
      cond_num <- sqrt(max(n_p) / min(n_p))
    } else {
      lanczos_max_eig <- function(matvec, dim_n, max_iter = 60L, conv_tol = 1e-12) {
        m <- min(dim_n, max_iter)
        # Deterministic starting vector covering both person and non-person blocks
        v <- sin(seq_len(dim_n) * 1.61803398875) + 0.5
        v <- v / sqrt(sum(v^2))
        Q <- matrix(0, dim_n, m)
        alpha <- numeric(m)
        beta <- numeric(m)
        prev_max <- -Inf
        k_used <- 0L
        for (k in seq_len(m)) {
          Q[, k] <- v
          w <- matvec(v)
          if (k > 1L) w <- w - beta[k - 1L] * Q[, k - 1L]
          a_k <- sum(v * w)
          alpha[k] <- a_k
          w <- w - a_k * v
          # Full reorthogonalization (two passes) for machine-precision stability
          if (k >= 1L) {
            Qk <- Q[, seq_len(k), drop = FALSE]
            w <- w - as.vector(Qk %*% crossprod(Qk, w))
            w <- w - as.vector(Qk %*% crossprod(Qk, w))
          }
          b_k <- sqrt(sum(w^2))
          beta[k] <- b_k
          k_used <- k
          if (k >= 4L || b_k < 1e-14) {
            T_mat <- diag(alpha[seq_len(k)], k, k)
            if (k > 1L) {
              for (j in seq_len(k - 1L)) {
                T_mat[j, j + 1L] <- beta[j]
                T_mat[j + 1L, j] <- beta[j]
              }
            }
            cur_max <- max(eigen(T_mat, symmetric = TRUE, only.values = TRUE)$values)
            if (b_k < 1e-14 || abs(cur_max - prev_max) <= conv_tol * max(abs(cur_max), 1)) {
              return(cur_max)
            }
            prev_max <- cur_max
          }
          v <- w / b_k
        }
        T_mat <- diag(alpha[seq_len(k_used)], k_used, k_used)
        if (k_used > 1L) {
          for (j in seq_len(k_used - 1L)) {
            T_mat[j, j + 1L] <- beta[j]
            T_mat[j + 1L, j] <- beta[j]
          }
        }
        max(eigen(T_mat, symmetric = TRUE, only.values = TRUE)$values)
      }

      inv_np_full <- 1 / n_p
      # G_O^{-1} = V diag(1 / d^2) V^T from exact TSQR SVD (s_o)
      V_o <- s_o$v
      inv_d2 <- 1 / (s_o$d^2)
      apply_G <- function(u) {
        u_p <- u[seq_len(P)]
        u_o <- u[(P + 1L):p]
        c(n_p * u_p + as.vector(S_po %*% u_o),
          as.vector(crossprod(S_po, u_p)) + as.vector(G_xx %*% u_o))
      }
      apply_G_inv <- function(u) {
        u_p <- u[seq_len(P)]
        u_o <- u[(P + 1L):p]
        rhs_o <- u_o - as.vector(crossprod(S_po, inv_np_full * u_p))
        w_o <- as.vector(V_o %*% (inv_d2 * as.vector(crossprod(V_o, rhs_o))))
        w_p <- inv_np_full * (u_p - as.vector(S_po %*% w_o))
        c(w_p, w_o)
      }
      lam_max <- lanczos_max_eig(apply_G, p)
      lam_inv_max <- lanczos_max_eig(apply_G_inv, p)
      cond_num <- sqrt(lam_max * lam_inv_max)
    }
    solver_used <- "schur_complement"
  }

  step_levels <- if (family == "pcm") labels[[threshold_facet]] else "shared"
  step_observed <- if (family == "pcm") tabulate(ids[[threshold_facet]], length(step_levels)) > 0L else TRUE
  step_parameters <- length(step_levels) * (K - 2L)
  step_rank <- sum(step_observed) * (K - 2L)
  structure(list(
    formula = formula, family = family, K = K, threshold_facet = threshold_facet, observed_only = observed_only,
    solver = solver_used,
    n_observations = n, n_dropped = n_planned - n,
    levels = nodes, interaction_levels = interaction_levels,
    main_parameters = main_p, interaction_parameters = p - main_p,
    n_components = length(unique(roots)), connected = length(unique(roots)) == 1L,
    rank = rank, n_parameters = p, nullity = p - rank, identifiable = rank == p && step_rank == step_parameters,
    step_rank = step_rank, step_parameters = step_parameters, step_nullity = step_parameters - step_rank,
    unobserved_step_levels = step_levels[!step_observed],
    adjacent_logit_rank = rank + step_rank, adjacent_logit_parameters = p + step_parameters,
    condition_number = cond_num,
    singular_values = sing_vals, tolerance = tol, design_matrix = X,
    constraint_map = effect_map, null_space = null, effect_null_space = effect_null
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
  cat(toupper(x$family), "design check (", x$solver %||% "dense_svd", ")\n", sep = "")
  cat("Rows:", x$n_observations, "  Dropped:", x$n_dropped, "\n")
  cat("Incidence components:", x$n_components, "\n")
  cat("Constrained predictor rank:", x$rank, "/", x$n_parameters, "\n")
  cat("Unidentified directions:", x$nullity, "\n")
  cat("Adjacent-logit rank including steps:", x$adjacent_logit_rank, "/", x$adjacent_logit_parameters, "\n")
  cat("Likelihood identification:", if (x$identifiable) "full structural rank" else "rank deficient", "\n")
  cat("This does not certify practical precision or finite maximum-likelihood estimates.\n")
  invisible(x)
}
