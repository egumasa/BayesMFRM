# Run with library(BayesMFRM) and the optional priorsense package.
# All calculations below refer to the bundled four-category writing RSM.
if (!requireNamespace("priorsense", quietly = TRUE)) stop("Install priorsense to run this example.")
source(system.file("examples", "fit-writing.R", package = "BayesMFRM"), local = TRUE)
stopifnot(isTRUE(fit$convergence$converged), standata(fit)$estimate_sd_itemrater == 0L,
          standata(fit)$fixed_sd_itemrater == 0.5,
          identical(get_facet_labels(fit$data_info, "item"), c("Content", "Evidence", "Structure")),
          identical(get_facet_labels(fit$data_info, "rater"), c("A", "B", "C")))
draws <- posterior::as_draws_df(posterior::as_draws(fit,
  variables = c("itemrater_raw", "itemrater", "rater", "item", "tau")))

# One factor only: the four independent Gaussian interaction coordinates.
# The nine double-centered cells are dependent. Their joint prior is determined by these four.
z <- posterior::as_draws_matrix(posterior::subset_draws(draws, variable = "itemrater_raw"))
stopifnot(ncol(z) == 4L)
draws$lprior_interaction <- -0.5 * rowSums(z^2)

# Conditional likelihood of the original 216 scores, retaining all fitted effects.
y <- fit$data_info$data_clean$score
p <- predict(fit, type = "probs", summary = FALSE)
stopifnot(identical(draws$.draw, attr(p, "bmfrm_prediction_meta")$draw_ids))
draws$log_lik <- Reduce(`+`, lapply(seq_along(y), function(n) log(p[, n, y[n]])))
stopifnot(all(is.finite(draws$log_lik)))

# Unit-slope information at ability zero for A and B rating Content and Evidence.
# The variance within each rating is conditional on the same posterior draw.
plan <- expand.grid(item = c("Content", "Evidence"), rater = c("A", "B"))
info <- numeric(nrow(draws))
for (j in 1:2) for (r in 1:2) {
  eta <- -draws[[sprintf("item[%d]", j)]] - draws[[sprintf("rater[%d]", r)]] +
    draws[[sprintf("itemrater[%d,%d]", j, r)]]
  logits <- matrix(0, nrow(draws), 4L)
  for (k in 1:3) logits[, k + 1L] <- logits[, k] + eta - draws[[sprintf("tau[%d]", k)]]
  probs <- exp(logits - apply(logits, 1L, max))
  probs <- probs / rowSums(probs)
  mu <- drop(probs %*% (1:4))
  info <- info + rowSums(probs * outer(mu, 1:4, "-")^2)
}
draws$information <- info
draws$content_A <- draws[["itemrater[1,1]"]]
draws$rater_A <- draws[["rater[1]"]]
stopifnot(abs(mean(info) - test_information(fit, plan, theta = 0)$information) < 1e-10)

targets <- c("content_A", "rater_A", "information")
ps <- priorsense::create_priorsense_data(
  posterior::subset_draws(draws, variable = c(targets, "lprior_interaction", "log_lik")),
  log_prior_name = "lprior_interaction")
sensitivity <- priorsense::powerscale_sensitivity(ps, variables = targets)
print(as.data.frame(sensitivity), row.names = FALSE)
# Each row is one posterior quantity. The prior column refers only to the interaction prior.

sequence <- priorsense::powerscale_sequence(ps, variables = targets,
  lower_alpha = 0.8, upper_alpha = 1.25, length = 5L)
scaled_means <- posterior::summarise_draws(sequence, "mean")[[1L]]
print(scaled_means[, c("component", ".powerscale_alpha", "variable", "mean", "pareto_k")],
      row.names = FALSE)
# Pareto k diagnoses the importance approximation. Inspect it before interpreting the curves.
labels <- c(content_A = "Content × A interaction (logits)", rater_A = "Rater A severity (logits)",
            information = "Four-rating information at ability 0")
power_plot <- ggplot2::ggplot(scaled_means,
    ggplot2::aes(.powerscale_alpha, mean, colour = component, linetype = component)) +
  ggplot2::geom_vline(xintercept = 1, colour = "grey75", linetype = "dotted") +
  ggplot2::geom_line(linewidth = 0.9) + ggplot2::geom_point(size = 1.8) +
  ggplot2::facet_wrap(~variable, scales = "free_y", labeller = ggplot2::as_labeller(labels)) +
  ggplot2::scale_colour_manual(values = c(prior = "#256368", likelihood = "#A55230"),
                                labels = c("Interaction prior", "Rating likelihood")) +
  ggplot2::scale_linetype_discrete(labels = c("Interaction prior", "Rating likelihood")) +
  ggplot2::labs(title = "How strongly does the writing assessment depend on the interaction prior?",
                x = "Power (1 = fitted model)", y = "Posterior mean", colour = NULL, linetype = NULL) +
  ggplot2::theme_minimal(base_size = 11) + ggplot2::theme(legend.position = "bottom")
print(power_plot)
