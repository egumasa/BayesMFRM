#' A compact table of Bayesian many-facet Rasch estimates
#'
#' Display posterior means and equal-tailed credible intervals in logits,
#' with one column per fitted model. Named fits supply column headings.
#' The table prints in the console and renders in a knitted report.
#'
#' @param ... One or more `bmfrm_fit` objects.
#' @param facets Character vector of main facets or fitted two-way interactions
#'   (either order). By default each model shows its non-ability main facets.
#'   An ability-only model shows its ability facet. Request the ability facet
#'   explicitly to include examinees in other models.
#' @param ci_level Probability for the equal-tailed credible intervals.
#' @param digits Integer from 0 to 15 for display. Numeric results retain
#'   their original precision.
#' @param file Optional path for a standalone UTF-8 HTML table. HTML export
#'   uses the suggested package knitr. An existing file at this path is replaced.
#' @details Rows align by facet and original level labels, including both
#'   components of an interaction. A dash denotes a level or term absent from
#'   that model's displayed estimates. A level with zero ratings is marked
#'   "no ratings". Effects known to be structurally unidentified are marked
#'   "unidentified". The footer reports the design and MCMC checks for each fit.
#'   Compare models with comparable facet meanings, orientation and scale.
#'   Intervals describe individual parameters. Use [facet_contrasts()] to
#'   estimate differences within a fit and [summarise_bias()] for bias screening.
#' @return A `bmfrm_table` list. `estimates` is a long tibble containing the
#'   model label, facet, role, first and second level, rating count, posterior mean,
#'   posterior SD, interval bounds and likelihood-identification status
#'   (`NA` if unchecked). `models` contains model metadata. `ci_level` and
#'   `digits` record the display settings. The object can be saved with [saveRDS()].
#' @seealso [facet_summary()], [summarise_bias()], [rater_agreement()]
#' @examples
#' data(bmfrm_example)
#' bmfrm_example$model_table
#' subset(bmfrm_example$model_table$estimates, facet == "rater")
#' \donttest{
#' if (check_cmdstan()) {
#'   source(system.file("examples", "fit-writing.R", package = "BayesMFRM"), local = TRUE)
#'   tab_bmfrm(Writing = fit, facets = c("rater", "item:rater"), ci_level = 0.90)
#'   if (requireNamespace("knitr", quietly = TRUE)) {
#'     tab_bmfrm(Writing = fit, file = tempfile(fileext = ".html"))
#'   }
#' }
#' }
#' @export
tab_bmfrm <- function(..., facets = NULL, ci_level = 0.95, digits = 2L, file = NULL) {
  fits <- list(...)
  if (!length(fits) || !all(vapply(fits, inherits, logical(1), "bmfrm_fit"))) {
    stop("Supply one or more bmfrm_fit objects")
  }
  probs <- bmfrm_interval_probs(ci_level)
  if (!is.numeric(digits) || length(digits) != 1L || !is.finite(digits) ||
      digits != floor(digits) || digits < 0 || digits > 15) stop("digits must be an integer in 0..15")
  if (!is.null(file) && (!is.character(file) || length(file) != 1L || is.na(file) || !nzchar(file))) {
    stop("file must be one HTML output path")
  }
  if (!is.null(file) && !requireNamespace("knitr", quietly = TRUE)) stop("Install knitr for HTML export")

  labels <- names(fits) %||% rep("", length(fits))
  unnamed <- !nzchar(labels)
  labels[unnamed] <- paste("Model", which(unnamed))
  if (anyDuplicated(labels)) stop("Use distinct model names")
  available <- unique(unlist(lapply(fits, function(x) c(x$spec$facets_main, x$spec$facets_bias))))
  interactions <- grep(":", available, value = TRUE)
  reverse <- vapply(strsplit(interactions, ":", fixed = TRUE),
                    function(z) paste(rev(z), collapse = ":"), character(1L))
  if (!is.null(facets) && (!is.character(facets) || !length(facets) || anyNA(facets) ||
      anyDuplicated(interaction_key(facets)) || !all(facets %in% c(available, reverse)))) {
    stop("facets must name distinct fitted main facets or interactions")
  }

  estimates <- lapply(seq_along(fits), function(i) {
    fit <- fits[[i]]
    main <- fit$spec$facets_main
    selected <- facets
    if (is.null(selected)) {
      selected <- setdiff(main, spec_person_facet(fit$spec))
      if (!length(selected)) selected <- main
    }
    fitted <- match_interaction_names(selected, fit$spec$facets_bias)
    keep <- which(fitted %in% c(main, fit$spec$facets_bias))
    dplyr::bind_rows(lapply(keep, function(j) {
      bmfrm_table_term(fit, labels[i], fitted[j], selected[j], probs, ci_level)
    }))
  })
  models <- dplyr::bind_rows(Map(bmfrm_model_overview, fits, labels))
  x <- structure(list(estimates = dplyr::bind_rows(estimates), models = models,
                      ci_level = ci_level, digits = digits), class = "bmfrm_table")
  if (!is.null(file)) {
    html <- knitr::kable(bmfrm_display_table(x), format = "html", row.names = FALSE,
                         caption = bmfrm_table_caption(x), escape = TRUE)
    writeLines(enc2utf8(c('<!doctype html><html lang="en"><head><meta charset="utf-8">',
      '<meta name="viewport" content="width=device-width, initial-scale=1"><title>Bayesian many-facet Rasch estimates</title>',
      '<style>body{font:16px/1.5 system-ui,sans-serif;color:#182d39;margin:2rem}main{overflow-x:auto}',
      'table{border-collapse:collapse}caption{text-align:left;font-weight:650;padding-bottom:1em}',
      'th,td{padding:.35em .8em;border-bottom:1px solid #d2dadd;vertical-align:top}thead{background:#edf3f4}</style>',
      '</head><body><main>', as.character(html), '</main></body></html>')), file, useBytes = TRUE)
  }
  x
}

bmfrm_table_term <- function(fit, label, term, requested, probs, ci_level) {
  if (term %in% fit$spec$facets_main) {
    s <- facet_summary(fit, term, probs = probs)
    lo <- fmt_prob_col(probs[1L])
    hi <- fmt_prob_col(probs[2L])
    lower <- if (lo == hi && probs[1L] != probs[2L]) facet_summary(fit, term, probs[1L])[[lo]] else s[[lo]]
    dc <- fit$design
    null <- dc$effect_null_space
    rows <- match(paste(term, s$level, sep = ":"), rownames(null))
    identified <- rep(NA, nrow(s))
    if (is.matrix(null) && length(dc$tolerance) == 1L) {
      found <- !is.na(rows)
      identified[found] <- rowSums(null[rows[found], , drop = FALSE]^2) < dc$tolerance^2
    }
    role <- if (term == spec_person_facet(fit$spec)) "ability" else "main"
    tibble::tibble(model = label, facet = term, role = role, level = s$level, level2 = NA_character_,
                   n_obs = facet_observation_counts(fit, term), mean = s$mean, sd = s$sd,
                   lower = lower, upper = s[[hi]], likelihood_identified = identified)
  } else {
    s <- summarise_bias(fit, term, ci_level = ci_level, diagnostics = FALSE)
    reversed <- requested != term
    tibble::tibble(model = label, facet = requested, role = "interaction",
                   level = if (reversed) s$level2 else s$level1,
                   level2 = if (reversed) s$level1 else s$level2, n_obs = s$n_observations,
                   mean = s$bias_mean, sd = s$bias_sd, lower = s$ci_lower, upper = s$ci_upper,
                   likelihood_identified = s$likelihood_identified)
  }
}

bmfrm_model_overview <- function(fit, label) {
  spec <- fit$spec
  dc <- fit$design$identifiable
  cv <- fit$convergence
  mcmc <- if (identical(cv$status, "incomplete")) "Incomplete" else if (isTRUE(cv$converged)) {
    if (isTRUE(cv$efficiency_warning) || isTRUE(cv$n_max_treedepth > 0)) "Passed with treedepth hits" else "Passed"
  } else if (identical(cv$converged, FALSE)) "Check diagnostics" else "Unchecked"
  rubric <- if (isTRUE(spec$K_heterogeneous)) {
    paste(paste(sort(unique(spec$K_by_level)), collapse = "/"), "by", spec$threshold_facet)
  } else as.character(get_spec_K_max(spec))
  counts <- lengths(fit$data_info$facet_labels[spec$facets_main])
  tibble::tibble(model = label, formula = paste(deparse(stats::formula(fit)), collapse = " "),
                 family = toupper(spec$family %||% "rsm"), categories = rubric, n_obs = stats::nobs(fit),
                 facet_levels = paste(names(counts), counts, collapse = ", "),
                 design = if (isTRUE(dc)) "Identified" else if (identical(dc, FALSE)) "Rank deficient" else "Unchecked",
                 mcmc = mcmc)
}

bmfrm_table_caption <- function(x) {
  sprintf("Posterior means and %g%% equal-tailed credible intervals (logits)", 100 * x$ci_level)
}

bmfrm_display_table <- function(x) {
  est <- x$estimates
  models <- x$models
  keys <- c("facet", "role", "level", "level2")
  rows <- unique(est[keys])
  fmt <- function(z) formatC(round(z, x$digits) + 0, digits = x$digits, format = "f")
  est$cell <- paste0(fmt(est$mean), " [", fmt(est$lower), ", ", fmt(est$upper), "]")
  est$cell <- paste0(est$cell, ifelse(est$n_obs == 0L, " (no ratings)", ""),
                     ifelse(est$likelihood_identified %in% FALSE, " (unidentified)", ""))
  level <- ifelse(is.na(rows$level2), rows$level, paste(rows$level, rows$level2, sep = " / "))
  footer <- c("Formula", "Family (categories)", "Ratings", "Facet levels", "Design", "MCMC")
  roles <- unique(rows[c("facet", "role")])
  mixed <- roles$facet[duplicated(roles$facet)]
  group <- ifelse(rows$facet %in% mixed, paste0(rows$facet, " (", rows$role, ")"), rows$facet)
  display <- data.frame(Facet = c(group, "Model", rep("", length(footer) - 1L)),
                        Level = c(level, footer), stringsAsFactors = FALSE)
  for (i in seq_len(nrow(models))) {
    m <- models[i, ]
    cells <- dplyr::left_join(rows, est[est$model == m$model, c(keys, "cell")], by = keys)$cell
    cells[is.na(cells)] <- "-"
    display[[i + 2L]] <- c(cells, m$formula, paste0(m$family, " (", m$categories, ")"),
                           as.character(m$n_obs), m$facet_levels, m$design, m$mcmc)
  }
  names(display)[-(1:2)] <- models$model
  display
}

#' @export
print.bmfrm_table <- function(x, ...) {
  cat(bmfrm_table_caption(x), "\n", sep = "")
  tab <- bmfrm_display_table(x)
  columns <- Map(function(header, cells) format(c(header, cells), justify = "left"), names(tab), tab)
  writeLines(trimws(do.call(paste, c(unname(columns), sep = "  ")), which = "right"))
  invisible(x)
}

#' @exportS3Method knitr::knit_print
knit_print.bmfrm_table <- function(x, ...) {
  tab <- bmfrm_display_table(x)
  fmt <- if (knitr::is_html_output()) "html" else NULL
  knitr::knit_print(knitr::kable(tab, format = fmt, row.names = FALSE,
                                caption = bmfrm_table_caption(x), escape = TRUE), ...)
}
