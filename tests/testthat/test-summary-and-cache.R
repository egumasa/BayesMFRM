test_that("fair-score uncertainty preserves dependence and groups by person", {
  mat <- cbind(`mu_fair[1]` = c(1, 2, 3, 4), `mu_fair[2]` = c(4, 3, 2, 1),
                `mu_fair[3]` = c(2, 3, 4, 5), `mu_fair[4]` = c(2, 3, 4, 5))
  object <- structure(list(fit = list(draws = function(...) mat),
                           spec = list(facets_main = c("item", "person"), response = "score"),
                           data_info = list(facet_labels = list(person = c("p", "q")),
                                            data_clean = data.frame(person = c("p", "p", "q", "q"),
                                                                    item = c(1, 2, 1, 2), score = c(2, 3, 3, 4)))),
                      class = "bmfrm_fit")
  x <- fair_scores(object)
  expect_equal(x$person, c("p", "q"))
  expect_equal(x$mean_fair_score, c(2.5, 3.5))
  expect_equal(x$sd_fair_score, c(0, stats::sd(2:5)))
  expect_equal(x$lower_95[1], x$upper_95[1])
})

test_that("facet summaries retain means and standard deviations with custom quantiles", {
  mat <- matrix(seq_len(20), 10, 2, dimnames = list(NULL, c("theta[1]", "theta[2]")))
  object <- structure(list(fit = list(summary = function(variables, ...) posterior::summarise_draws(posterior::as_draws_matrix(mat), ...)),
                           spec = list(facets_main = "person"),
                           data_info = list(facet_labels = list(person = c("p", "q")))), class = "bmfrm_fit")
  x <- facet_summary(object, "person")
  expect_equal(unname(x$mean), unname(colMeans(mat)))
  expect_equal(unname(x$sd), unname(apply(mat, 2, stats::sd)))
  expect_equal(x$q50, c(5.5, 15.5))
})

test_that("cache keys include settings and missing CSVs invalidate a fit", {
  directory <- tempfile("cache_test_")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE))
  counter <- 0L
  model <- list(code = function() "model {}", sample = function(output_dir, ...) {
    counter <<- counter + 1L
    path <- file.path(output_dir, paste0("chain", counter, ".csv"))
    writeLines("test", path)
    list(id = counter, output_files = function() path)
  })
  file <- file.path(directory, "fit.rds")
  x <- fit_cmdstan_cached(model, list(N = 1), file, "on_change", seed = 1, iter_sampling = 100)
  y <- fit_cmdstan_cached(model, list(N = 1), file, "on_change", seed = 1, iter_sampling = 100)
  expect_equal(x$id, y$id)
  z <- fit_cmdstan_cached(model, list(N = 1), file, "on_change", seed = 1, iter_sampling = 200)
  expect_equal(z$id, 2L)
  expect_error(fit_cmdstan_cached(model, list(N = 1), file, "never", seed = 2, iter_sampling = 200), "cannot be verified")
  unlink(z$output_files())
  expect_error(fit_cmdstan_cached(model, list(N = 1), file, "never", seed = 1, iter_sampling = 200), "missing")
  w <- fit_cmdstan_cached(model, list(N = 1), file, "on_change", seed = 1, iter_sampling = 200)
  expect_equal(w$id, 3L)

  # Unseeded calls (seed = NULL) run fresh chains under refit = "on_change", whereas seeded calls reuse cache
  file_min <- file.path(directory, "fit_minimal_dots.rds")
  m1 <- fit_cmdstan_cached(model, list(N = 1), file_min, "on_change")
  m2 <- fit_cmdstan_cached(model, list(N = 1), file_min, "on_change")
  expect_equal(m1$id, 4L)
  expect_equal(m2$id, 5L)
  s1 <- fit_cmdstan_cached(model, list(N = 1), file_min, "on_change", seed = 7)
  s2 <- fit_cmdstan_cached(model, list(N = 1), file_min, "on_change", seed = 7)
  expect_equal(s1$id, 6L)
  expect_equal(s2$id, 6L)
})

test_that("PPC supports groups, small posterior samples, and reproducible selection", {
  mat <- matrix(rep(1:3, 8), 4, 6)
  fit <- list(draws = function(...) mat)
  before <- .Random.seed
  x <- suppressMessages(run_ppc(fit, c(1, 2, 3, 1, 2, 3), ndraws = 20, seed = 1, group = rep(1:2, each = 3)))
  expect_identical(.Random.seed, before)
  expect_equal(nrow(x$yrep_subset), 4L)
  expect_true(inherits(x$ppc_bars, "ggplot"))
})
