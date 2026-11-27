#' Conditional category probabilities and adjacent-category thresholds
#'
#' For each combination of facet levels, evaluate the fitted rating scale over
#' a grid of fixed abilities. The probability curves and the Wright map use
#' the same draw-by-draw adjacent-category thresholds.
#'
#' @param object A fitted `bmfrm_fit` with a verified R evaluation contract.
#' @param newdata Data frame with one row per profile. Supply the fitted level
#'   of every non-ability facet, e.g. `data.frame(item = "Content", rater = "A")`.
#'   Include the ability facet only when it also indexes an interaction or PCM
#'   thresholds. Its label then selects those parameters while `theta` varies.
#' @param theta Finite abilities in logits. The default grid spans -4 to 4.
#' @param ci_level Probability mass for equal-tailed posterior intervals.
#' @details Within each posterior draw the log odds of category `k + 1`
#'   relative to `k` are the ability minus their threshold. Thresholds include
#'   all fitted main and interaction effects for the chosen profile. An
#'   adjacent pair is equally likely at its threshold even if a third category
#'   is more probable. Step order is preserved, including disordered steps.
#'   The median threshold can differ from the intersection of the posterior
#'   mean curves. In a heterogeneous PCM each profile uses its own rubric.
#'   A profile with no observed ratings is a model-based extrapolation.
#' @return `category_curves()` returns a tibble with the input row number
#'   (`profile`), its levels (`facet_<name>`), ability (`theta`), integer
#'   `category`, posterior mean `probability`, and pointwise `lower` and
#'   `upper` bounds at `ci_level`. The means sum to one within each profile
#'   and ability. Only categories on that profile's rubric appear.
#'   `category_thresholds()` returns one row per adjacent pair, with
#'   `from`, `to`, and posterior `median`, `lower`, `upper` in logits.
#'   Both tables can be used directly with ggplot2 or exported with
#'   `utils::write.csv()`.
#' @seealso [plot_category_curves()], [category_contrast()], [plot_wright_map()],
#'   [predict.bmfrm_fit()]
#' @examples
#' \donttest{
#' if (check_cmdstan()) {
#'   source(system.file("examples", "fit-writing.R", package = "BayesMFRM"),
#'          local = TRUE)
#'   profiles <- data.frame(item = "Content", rater = c("A", "B"))
#'   p <- category_curves(fit, profiles, theta = 0)
#'   p[p$category == 4, c("facet_rater", "probability", "lower", "upper")]
#'   # At ability zero the top score is more probable with A than with B.
#'   category_thresholds(fit, profiles)
#'
#'   # Fit different rubrics to simulated ratings on the same assignment design.
#'   rubric <- c(Content = 2L, Evidence = 4L, Structure = 6L)
#'   f <- score ~ person + item + rater
#'   sim <- simulate_bmfrm(f, ratings, K = rubric, family = "pcm",
#'                         threshold_facet = "item", seed = 412)
#'   pcm <- bmfrm(f, sim$data, K = rubric, family = "pcm", threshold_facet = "item",
#'                seed = 413, chains = 4, cores = 2, iter = 2000, adapt_delta = 0.95,
#'                cpp_options = list(PRECOMPILED_HEADERS = "false"))
#'   mixed <- data.frame(item = names(rubric), rater = "A")
#'   table(category_curves(pcm, mixed, theta = 0)$facet_item) # 2, 4, 6 categories
#'   plot_wright_map(pcm, mixed) # 1, 3, 5 adjacent thresholds
#' }
#' }
#' @export
category_curves <- function(object, newdata, theta = seq(-4, 4, length.out = 81),
                            ci_level = 0.95) {
  theta <- bmfrm_curve_grid(theta)
  profiles <- bmfrm_curve_profiles(object, newdata, ci_level)
  tables <- lapply(seq_len(nrow(newdata)), function(i) {
    ct <- bmfrm_profile_cumulative(profiles, i)
    D <- nrow(ct)
    K <- ncol(ct)
    size <- default_bmfrm_chunk_size(D * as.double(K), length(theta))
    chunks <- split(theta, ceiling(seq_along(theta) / size))
    dplyr::bind_rows(lapply(chunks, function(grid) {
      p <- bmfrm_curve_draws(ct, grid)
      bounds <- apply(p, 2L, stats::quantile, probs = profiles$probs, names = FALSE)
      tibble::tibble(profile = i, theta = rep(grid, K),
                     category = rep(seq_len(K), each = length(grid)),
                     probability = colMeans(p), lower = bounds[1L, ], upper = bounds[2L, ],
                     ci_level = ci_level)
    }))
  })
  add_curve_profile_levels(dplyr::bind_rows(tables), newdata)
}

bmfrm_curve_grid <- function(theta) {
  if (!is.numeric(theta) || !is.null(dim(theta)) || !length(theta) || any(!is.finite(theta))) {
    stop("theta must be a nonempty vector of finite abilities")
  }
  sort(unique(theta))
}

# Columns run through the ability grid for each successive score category.
bmfrm_curve_draws <- function(ct, grid) {
  D <- nrow(ct)
  K <- ncol(ct)
  eta <- matrix(rep(grid, each = D), nrow = D)
  p <- eval_bmfrm_categories(eta, K, list(mode = "shared", cum_tau = ct), NULL,
                             compute_mu = FALSE, compute_probs = TRUE)$probs
  dim(p) <- c(D, length(grid) * K)
  p
}

#' @rdname category_curves
#' @export
category_thresholds <- function(object, newdata, ci_level = 0.95) {
  profiles <- bmfrm_curve_profiles(object, newdata, ci_level)
  cuts <- lapply(seq_len(nrow(newdata)), function(i) {
    steps <- bmfrm_step_draws(bmfrm_profile_cumulative(profiles, i))
    bmfrm_crossings(steps, i, profiles$probs, ci_level)
  })
  add_curve_profile_levels(dplyr::bind_rows(cuts), newdata)
}

bmfrm_curve_profiles <- function(object, newdata, ci_level) {
  probs <- bmfrm_interval_probs(ci_level)
  if (!inherits(object, "bmfrm_fit") || !can_evaluate_in_r(object)) {
    stop("Category diagnostics require a fitted model with a verified R evaluation contract")
  }
  spec <- object$spec
  pf <- spec_person_facet(spec)
  bias_facets <- unlist(strsplit(as.character(spec$facets_bias), ":", fixed = TRUE))
  need_person <- pf %in% bias_facets || identical(spec$threshold_facet, pf)
  omit <- if (!need_person) pf else character()
  required <- setdiff(spec$facets_main, omit)
  if (!is.data.frame(newdata) || !nrow(newdata) || anyDuplicated(names(newdata)) ||
      !setequal(names(newdata), required)) {
    stop("Supply one row per profile in newdata with exactly these facets: ",
         paste(required, collapse = ", "))
  }
  ids <- map_eval_indices(newdata, spec, object$data_info$facet_labels, adjust_facets = omit)
  list(params = extract_bmfrm_param_draws(object), ids = ids, probs = probs)
}

# Adding the fitted profile offset to each adjacent logit is equivalent to
# subtracting k copies from the cumulative threshold for category k + 1.
bmfrm_profile_cumulative <- function(profiles, i) {
  params <- profiles$params
  ti <- params$tau_info
  offset <- bmfrm_add_facet_effects(matrix(0, params$D, 1L), params, profiles$ids, i)[, 1L]
  if (ti$mode == "shared") {
    ct <- ti$cum_tau
  } else {
    level <- profiles$ids[[ti$threshold_facet]][i]
    ct <- if (ti$mode == "pcm_ragged") ti$cum_tau_by_level[[level]] else
      do.call(cbind, lapply(ti$cum_tau_list, function(m) m[, level]))
  }
  ct - outer(offset, seq_len(ncol(ct)) - 1L)
}

bmfrm_interval_probs <- function(ci_level) {
  if (!is.numeric(ci_level) || length(ci_level) != 1L || !is.finite(ci_level) ||
      ci_level <= 0 || ci_level >= 1) stop("ci_level must be a probability in (0, 1)")
  c((1 - ci_level) / 2, (1 + ci_level) / 2)
}

add_curve_profile_levels <- function(table, newdata) {
  levels <- newdata[table$profile, , drop = FALSE]
  names(levels) <- paste0("facet_", names(levels), recycle0 = TRUE)
  dplyr::bind_cols(table, levels)
}

bmfrm_step_draws <- function(ct) ct[, -1L, drop = FALSE] - ct[, -ncol(ct), drop = FALSE]

bmfrm_crossings <- function(steps, i, probs, ci_level) {
  qs <- apply(steps, 2L, stats::quantile, probs = c(probs[1L], 0.5, probs[2L]), names = FALSE)
  from <- seq_len(ncol(steps))
  tibble::tibble(profile = i, from = from, to = from + 1L,
                 median = qs[2L, ], lower = qs[1L, ], upper = qs[3L, ], ci_level = ci_level)
}

bmfrm_profile_labels <- function(newdata) {
  if (!ncol(newdata)) return(rep("Overall", nrow(newdata)))
  values <- Map(function(nm, x) paste0(nm, " = ", x), names(newdata), newdata)
  do.call(paste, c(unname(values), sep = ", "))
}
