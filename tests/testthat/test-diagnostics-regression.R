mock_diagnostic_fit <- function(rhat = 1.001, depth = c(0, 0), ebfmi = c(.8, .9), chains = 2L) {
  list(
    metadata = function() list(model_params = c("theta_raw[1]", "mu[1]")),
    summary = function(variables, ...) {
      stopifnot(identical(variables, "theta_raw"))
      data.frame(variable = "theta_raw[1]", mean = 0, sd = 1,
                 rhat = rhat, ess_bulk = 1000, ess_tail = 900, mcse_mean = .03)
    },
    diagnostic_summary = function(...) list(num_divergent = c(0, 0), num_max_treedepth = depth, ebfmi = ebfmi),
    num_chains = function() chains
  )
}

test_that("diagnostics use the sampler's configured depth limit and include E-BFMI", {
  good <- check_convergence(mock_diagnostic_fit())
  expect_true(good$converged)
  expect_equal(good$status, "ok")
  expect_equal(good$max_mcse_mean, .03)
  expect_false(check_convergence(mock_diagnostic_fit(depth = c(1, 0)))$converged)
  expect_false(check_convergence(mock_diagnostic_fit(ebfmi = c(.2, .9)))$converged)
})

test_that("missing Rhat and single chains cannot pass the diagnostic gate", {
  expect_false(check_convergence(mock_diagnostic_fit(rhat = NA_real_))$converged)
  expect_equal(check_convergence(mock_diagnostic_fit(rhat = NA_real_))$status, "incomplete")
  expect_false(check_convergence(mock_diagnostic_fit(chains = 1L))$converged)
})

test_that("prior concatenation preserves prior objects as documented", {
  ps <- c(prior("normal(0, 2)", "theta"), prior("normal(0, 1)", "item"),
          prior("normal(0, 1)", "tau"))
  expect_length(ps, 3L)
  expect_true(all(vapply(ps, inherits, logical(1), "bmfrm_prior")))
  expect_no_error(validate_priors(ps, c("person", "item")))
})

test_that("integer validation is applied even with an explicit K", {
  d <- data.frame(person = rep(1:5, 2), item = rep(1:2, each = 5), score = rep(2.7, 10))
  expect_error(parse_bmfrm_formula(score ~ person + item, d, K = 6), "integer")
})

test_that("repeated ratings require an explicit residual pairing rule", {
  d <- expand.grid(person = 1:5, item = 1:3, rater = 1:2)
  d$resid <- sin(seq_len(nrow(d)))
  expect_error(residual_cor_criteria(d), "unit_facets|aggregate")
  x <- residual_cor_criteria(d, unit_facets = c("person", "rater"))
  expect_equal(dim(x$cor_matrix), c(3L, 3L))
  expect_equal(x$n_pairs, 10L)
  expect_equal(unname(x$pair_counts[1, 2]), 10L)
  y <- residual_cor_criteria(d, aggregate = "mean")
  expect_equal(y$n_pairs, 5L)
  expect_equal(y$aggregation, "mean")
})
