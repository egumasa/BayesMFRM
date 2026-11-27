# Rebuild the worked example from the repository root with Rscript --vanilla.
pkgload::load_all(".", quiet = TRUE)
model <- score ~ person + item + rater + item:rater
persons <- sprintf("P%02d", 1:36)
items <- c("Content", "Evidence", "Structure")
raters <- c("A", "B", "C")
grid <- expand.grid(item = items, rater = raters, stringsAsFactors = FALSE)
set.seed(2407)
plan <- do.call(rbind, lapply(persons, function(id) {
  data.frame(person = id, grid[sample.int(nrow(grid), 6L), ], row.names = NULL)
}))
# A rates Content more generously and Evidence more severely; B does the reverse.
parameters <- list(theta = stats::rnorm(length(persons)), item = c(-0.4, 0, 0.4),
                   rater = c(-0.5, 0, 0.5), tau = c(-1, 0, 1),
                   itemrater = matrix(c(1.5, -1.5, 0, -1.5, 1.5, 0, 0, 0, 0), 3))
ratings <- simulate_bmfrm(model, plan, K = 4, parameters = parameters, seed = 2408)$data
fit <- bmfrm(model, ratings, K = 4, design_check = "error", seed = 2409,
             chains = 4, cores = 2, iter = 2000, adapt_delta = 0.95,
             cpp_options = list(PRECOMPILED_HEADERS = "false"))
stopifnot(isTRUE(fit$convergence$converged))
bmfrm_example <- list(ratings = ratings,
                      bias = summarise_bias(fit, "item:rater", delta = 0.3),
                      fair = fair_scores(fit, standardize = TRUE),
                      rater_fit = facet_fit(fit, "rater"),
                      rater = facet_summary(fit, "rater"))
save(bmfrm_example, file = "data/bmfrm_example.rda", version = 2, compress = "gzip")
