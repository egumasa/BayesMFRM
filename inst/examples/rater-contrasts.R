# Run after fit-writing.R. Positive differences indicate greater severity.
pairs <- facet_contrasts(fit, "rater", c("A", "B", "C"))
labels <- paste(pairs$comparison, "minus", pairs$reference)
pairs$pair <- factor(labels, levels = rev(labels))
contrast_plot <- ggplot2::ggplot(pairs, ggplot2::aes(x = difference, y = pair)) +
  ggplot2::geom_vline(xintercept = 0, colour = "grey55", linetype = "dashed") +
  ggplot2::geom_segment(ggplot2::aes(x = lower, xend = upper, yend = pair),
                        linewidth = 1, colour = "#216782") +
  ggplot2::geom_point(size = 3, colour = "#154760") +
  ggplot2::scale_x_continuous(breaks = seq(-0.25, 1, by = 0.25)) +
  ggplot2::labs(title = "How much stricter is one rater than another?",
                subtitle = "Posterior mean and 95% credible interval for each pair",
                x = "Difference in rater severity (logits)", y = NULL) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(),
                 panel.grid.minor = ggplot2::element_blank())
