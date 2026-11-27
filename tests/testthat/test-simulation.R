test_that("binary RSM agrees with the Bernoulli logistic special case", {
  prob <- rsm_probabilities(c(-2, 0, 2), 0)
  expect_equal(prob[, 2], stats::plogis(c(-2, 0, 2)))
  expect_equal(rowSums(prob), rep(1, 3))
  expect_equal(rsm_probabilities(0, c(-1, 0, 1))[1, ],
               c(1, exp(1), exp(1), 1) / (2 + 2 * exp(1)))
})

test_that("simulated scores obey known probabilities and restore seeded RNG state", {
  d <- data.frame(person = "p", item = "i", rater = "r")
  pars <- list(theta = log(3), item = 0, rater = 0, tau = 0)
  set.seed(100)
  before <- .Random.seed
  x <- simulate_bmfrm(score ~ person + item + rater, d, 2, parameters = pars, nsim = 4000, seed = 99)
  expect_identical(.Random.seed, before)
  expect_equal(dim(x$y_rep), c(4000L, 1L))
  expect_lt(abs(mean(x$y_rep == 2) - .75), .03)
  expect_identical(x, simulate_bmfrm(score ~ person + item + rater, d, 2, parameters = pars, nsim = 4000, seed = 99))
})

test_that("prior generation applies centering and retains the induced prior variance", {
  expect_identical(pp_bmfrm, prior_predict_bmfrm)
  d <- expand.grid(person = 1:2, item = 1:4, rater = 1:3)
  x <- pp_bmfrm(score ~ person + item + rater, d, 6, ndraws = 1500, seed = 1)
  expect_true(all(vapply(x$parameters, function(p) abs(sum(p$item)) < 1e-10 && abs(sum(p$rater)) < 1e-10 && abs(sum(p$tau)) < 1e-10, logical(1))))
  item1 <- vapply(x$parameters, function(p) p$item[1L], numeric(1))
  expect_equal(stats::var(item1), 1 - 1/4, tolerance = .12)
  expect_true(all(x$y_rep >= 1 & x$y_rep <= 6))
})

test_that("fixed mixed-rubric replications draw fresh scores in design order", {
  d <- data.frame(person = "P", item = factor(c("long", "binary", "long", "binary"),
                  levels = c("binary", "unused", "long")))
  K <- c(binary = 2L, unused = 3L, long = 4L)
  pars <- list(theta = 0, item = c(0, 0, 0), tau = list(binary = 0, unused = c(-1, 1),
                                                         long = c(-1, 0, 1)))
  x <- simulate_bmfrm(score ~ person + item, d, K, family = "pcm", threshold_facet = "item",
                      parameters = pars, nsim = 16L, seed = 42)
  u <- with_bmfrm_seed(42, t(matrix(runif(16L * nrow(d)), nrow(d))))
  cuts <- cumsum(c(1, exp(1), exp(1)) / (2 + 2 * exp(1)))
  expected <- vapply(seq_len(nrow(d)), function(j) {
    thresholds <- if (d$item[j] == "binary") 0.5 else cuts
    as.integer(1L + rowSums(outer(u[, j], thresholds, ">")))
  }, integer(16L))
  expect_identical(x$y_rep, expected)
})

test_that("invalid priors and fixed constraints are rejected rather than silently changed", {
  expect_error(sample_bmfrm_prior("normal(0, sqrt(2))", 2), "literal")
  expect_error(sample_bmfrm_prior("normal(0, -1)", 2), "positive")
  d <- expand.grid(person = 1:2, item = 1:2)
  expect_error(simulate_bmfrm(score ~ person + item, d, 2,
                              parameters = list(theta = c(0, 0), item = c(1, 0), tau = 0)), "sum to zero")
})
