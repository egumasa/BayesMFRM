#' Category probability surfaces over ability and a severity shift
#'
#' How would the chance of each rubric score change if a given rating profile
#' became more severe? The horizontal axes are ability and a hypothetical
#' additive shift in severity, both in logits. A positive shift raises every
#' adjacent-category threshold for the profile by the same amount. The fitted
#' facet effects, interactions and step spacings retain their joint posterior.
#' Thus a surface at ability \eqn{\theta} and shift \eqn{s} equals the profile's
#' ordinary category curve at \eqn{\theta - s}. Zero shift recovers the fitted
#' curve. The model has one ability dimension. To compare actual raters,
#' supply their profiles at zero shift or use [category_contrast()].
#'
#' @inheritParams category_curves
#' @param severity_shift Finite changes in severity relative to each fitted
#'   profile, in logits. The default spans -2 to 2 in steps of 0.25.
#' @param category Optional vector of integer rubric scores to plot. The
#'   default shows every supported score. For different rubrics each panel
#'   uses a category supported by its profile.
#' @return `category_surfaces()` returns a tibble with the same profile,
#'   category, posterior mean probability, and pointwise credible bounds as
#'   [category_curves()], together with `theta` and `severity_shift`.
#'   `plot_category_surfaces()` returns a lattice `trellis` object. Each panel
#'   shows one category for one profile. Height and colour both encode posterior
#'   mean probability on the same 0 to 1 scale across panels. Inspect `lower`
#'   and `upper` in the table for uncertainty at a chosen coordinate.
#' @examples
#' \donttest{
#' if (check_cmdstan()) {
#'   source(system.file("examples", "fit-writing.R", package = "BayesMFRM"), local = TRUE)
#'   profile <- data.frame(item = "Content", rater = "A")
#'   p <- category_surfaces(fit, profile, theta = 0, severity_shift = c(0, 1))
#'   subset(p, category == 4, select = c(severity_shift, probability, lower, upper))
#'   plot(fit, type = "surface", newdata = profile)
#' }
#' }
#' @export
category_surfaces <- function(object, newdata, theta = seq(-4, 4, by = 0.25),
                              severity_shift = seq(-2, 2, by = 0.25), ci_level = 0.95) {
  grid <- expand.grid(theta = bmfrm_curve_grid(theta),
                      severity_shift = bmfrm_curve_grid(severity_shift, "severity_shift"),
                      KEEP.OUT.ATTRS = FALSE)
  grid$effective_theta <- grid$theta - grid$severity_shift
  curves <- category_curves(object, newdata, theta = unique(grid$effective_theta), ci_level = ci_level)
  names(curves)[names(curves) == "theta"] <- "effective_theta"
  dat <- merge(grid, curves, by = "effective_theta", sort = FALSE)
  dat$effective_theta <- NULL
  tibble::as_tibble(dat[order(dat$profile, dat$category, dat$severity_shift, dat$theta), ])
}

#' @rdname category_surfaces
#' @export
plot_category_surfaces <- function(object, newdata, theta = seq(-4, 4, by = 0.25),
                                   severity_shift = seq(-2, 2, by = 0.25), category = NULL,
                                   ci_level = 0.95) {
  if (length(bmfrm_curve_grid(theta)) < 2L ||
      length(bmfrm_curve_grid(severity_shift, "severity_shift")) < 2L) {
    stop("A surface plot needs at least two distinct values on each horizontal axis")
  }
  dat <- category_surfaces(object, newdata, theta, severity_shift, ci_level)
  if (!is.null(category)) {
    if (!is.numeric(category) || !is.null(dim(category)) || !length(category) ||
        !all(category %in% dat$category)) stop("category must contain supported rubric scores")
    dat <- dat[dat$category %in% category, ]
  }
  panels <- unique(dat[c("profile", "category")])
  key <- paste(dat$profile, dat$category)
  dat$panel <- factor(key, levels = unique(key))
  profiles <- bmfrm_profile_labels(newdata)
  titles <- paste("Score", panels$category)
  if (nrow(newdata) > 1L) titles <- paste(titles, profiles[panels$profile], sep = "  |  ")
  caption <- "Ability and severity shift in logits. Height and colour show posterior mean probability."
  if (nrow(newdata) == 1L) caption <- paste(profiles[1L], caption, sep = "\n")
  lattice::wireframe(probability ~ theta * severity_shift | panel, data = dat,
    as.table = TRUE, strip = lattice::strip.custom(factor.levels = titles),
    main = "Category probability surfaces", sub = list(caption, cex = 0.8, font = 1),
    xlab = list("Ability", rot = -24, cex = 0.78), ylab = list("Severity shift", rot = 30, cex = 0.78),
    zlab = list("Probability", rot = 90, cex = 0.75),
    zlim = c(0, 1), drape = TRUE, at = seq(0, 1, length.out = 65L),
    col.regions = grDevices::hcl.colors(64L, "viridis"), col = "#394855", lwd = 0.3,
    colorkey = list(space = "right", labels = list(at = seq(0, 1, by = 0.25))),
    scales = list(arrows = FALSE, cex = 0.65, x = list(distance = 1.65), y = list(distance = 1.65),
                  z = list(distance = 1, at = c(0, 0.5, 1))),
    screen = list(z = -40, x = -62), aspect = c(1, 0.7),
    par.settings = list(axis.line = list(col = "grey55"), strip.background = list(col = "#EBF1F3")))
}
