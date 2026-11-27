# Run with library(BayesMFRM) and CmdStan. The density estimates use 40,000 draws.
if (!requireNamespace("bayestestR", quietly = TRUE) ||
    !requireNamespace("logspline", quietly = TRUE)) {
  stop("Install bayestestR and logspline to run this example.")
}
data(bmfrm_example, package = "BayesMFRM")
fit_evidence <- bmfrm(score ~ person + item + rater + item:rater, bmfrm_example$ratings,
                      K = 4, design_check = "error", seed = 2411, chains = 4, cores = 2,
                      warmup = 1000, iter = 11000, adapt_delta = 0.95,
                      cpp_options = list(PRECOMPILED_HEADERS = "false"))
stopifnot(isTRUE(fit_evidence$convergence$converged))
i <- match("Structure", get_facet_labels(fit_evidence$data_info, "item"))
r <- match("B", get_facet_labels(fit_evidence$data_info, "rater"))
cell <- sprintf("itemrater[%d,%d]", i, r)
draws <- posterior::as_draws(fit_evidence, variables = c(cell, "itemrater_raw"))
b <- posterior::as_draws_df(draws)[[cell]]

# The fitted Gaussian prior acts on independent contrasts. Project it onto this cell.
dat <- standata(fit_evidence)
raw_sd <- dat$fixed_sd_itemrater
cell_sd <- raw_sd * sqrt(sum(dat$B_itemrater_1[i, ]^2) * sum(dat$B_itemrater_2[r, ]^2))
stopifnot(dat$estimate_sd_itemrater == 0L, isTRUE(all.equal(cell_sd, 1 / 3)))
cell_prior <- stats::qnorm(stats::ppoints(length(b)), sd = cell_sd)
margin <- c(-0.3, 0.3)
hdi <- bayestestR::hdi(b, ci = 0.95)
support <- bayestestR::si(b, prior = cell_prior, BF = 1)
in_rope <- bayestestR::rope(b, range = margin, ci = 1)$ROPE_Percentage
decision <- bayestestR::equivalence_test(b, range = margin, ci = 0.95, ci_method = "HDI")
point_bf <- bayestestR::bayesfactor_parameters(b, prior = cell_prior, null = 0)
region_bf <- bayestestR::bayesfactor_parameters(b, prior = cell_prior, null = margin)
print(hdi)
print(support)
print(data.frame(prior_in_rope = diff(stats::pnorm(margin, sd = cell_sd)),
                 posterior_in_rope = in_rope, decision = decision$ROPE_Equivalence,
                 BF_unrestricted_vs_zero = exp(point_bf$log_BF),
                 BF_outside_vs_inside = exp(region_bf$log_BF)), row.names = FALSE)

# Both models keep the same likelihood. Only the prior on the four contrasts changes.
z <- posterior::as_draws_matrix(posterior::subset_draws(draws, variable = "itemrater_raw"))
interaction_contrasts <- raw_sd * z
tighter_sd <- raw_sd / sqrt(1.25)
log_ratio <- rowSums(stats::dnorm(interaction_contrasts, sd = tighter_sd, log = TRUE) -
                     stats::dnorm(interaction_contrasts, sd = raw_sd, log = TRUE))
ratio <- exp(log_ratio)
print(c(BF_tighter_vs_original = mean(ratio),
        MCSE = posterior::mcse_mean(matrix(ratio, ncol = posterior::nchains(draws)))))

grid <- seq(-1, 1.1, length.out = 512)
curves <- data.frame(
  logit = rep(grid, 2), distribution = rep(c("Prior", "Posterior"), each = length(grid)),
  density = c(stats::dnorm(grid, sd = cell_sd), logspline::dlogspline(grid, logspline::logspline(b))))
evidence_plot <- ggplot2::ggplot(curves, ggplot2::aes(logit, density)) +
  ggplot2::annotate("rect", xmin = margin[1], xmax = margin[2], ymin = -Inf, ymax = Inf,
                    fill = "#E6F0EB") +
  ggplot2::geom_vline(xintercept = 0, linetype = "dotted", colour = "grey45") +
  ggplot2::geom_line(ggplot2::aes(colour = distribution, linetype = distribution), linewidth = 1) +
  ggplot2::annotate("segment", x = hdi$CI_low, xend = hdi$CI_high, y = 0, yend = 0,
                    linewidth = 2, colour = "#1B5C79") +
  ggplot2::scale_colour_manual(values = c(Posterior = "#1B5C79", Prior = "#9B4B2C")) +
  ggplot2::scale_linetype_manual(values = c(Posterior = "solid", Prior = "dashed")) +
  ggplot2::labs(title = "How did the ratings change support for rater B on Structure?",
                subtitle = "Shading marks the practical margin of ±0.3 logits.",
                x = "Criterion-by-rater interaction (logits)", y = "Density",
                colour = NULL, linetype = NULL,
                caption = paste("The blue bar spans the 95% HDI.",
                                "Values gain support where the posterior rises above the prior.")) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(legend.position = "top", plot.caption = ggplot2::element_text(hjust = 0))
print(evidence_plot)
