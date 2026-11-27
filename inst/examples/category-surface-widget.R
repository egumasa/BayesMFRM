# Run after fit-writing.R. Rotate the surface and select a score with the menu.
profile <- data.frame(item = "Content", rater = "A")
surface <- category_surfaces(fit, profile)
abilities <- sort(unique(surface$theta))
shifts <- sort(unique(surface$severity_shift))
scores <- sort(unique(surface$category))
widget <- plotly::plot_ly()
for (k in scores) {
  rows <- surface[surface$category == k, ]
  rows <- rows[order(rows$severity_shift, rows$theta), ]
  z <- matrix(rows$probability, ncol = length(abilities), byrow = TRUE)
  bounds <- sprintf("%.1f%% to %.1f%%", 100 * rows$lower, 100 * rows$upper)
  # Plotly's surface hover text is indexed by x then y, while z is indexed by y then x.
  widget <- plotly::add_surface(widget, x = abilities, y = shifts, z = z,
    text = matrix(bounds, nrow = length(abilities)),
    hovertemplate = paste0("Ability: %{x:.2f}<br>Severity shift: %{y:.2f}",
      "<br>Probability: %{z:.1%}<br>95% credible interval: %{text}<extra></extra>"),
    name = paste("Score", k), visible = k == max(scores), colorscale = "Viridis",
    cmin = 0, cmax = 1, colorbar = list(title = "Probability", tickformat = ".0%"))
}
buttons <- lapply(scores, function(k) list(label = paste("Score", k), method = "restyle",
                                          args = list("visible", as.list(scores == k))))
widget <- plotly::layout(widget,
  title = list(text = "Score probabilities for Content rated by A", y = 0.97),
  margin = list(t = 65, b = 35),
  scene = list(xaxis = list(title = "Ability (logits)"),
               yaxis = list(title = "Severity shift (logits)"),
               zaxis = list(title = "Probability", range = c(0, 1), tickformat = ".0%")),
  updatemenus = list(list(buttons = buttons, active = length(scores) - 1L, x = 0, y = 1)),
  uirevision = "preserve-camera")
widget
