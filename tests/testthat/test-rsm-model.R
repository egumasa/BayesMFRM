test_that("isometric blocks have the declared constrained covariance", {
  for (n in c(2L, 3L, 7L)) {
    B <- isometric_rsm_basis(n)
    expect_equal(colSums(B), rep(0, n - 1L), tolerance = 1e-14)
    V <- 1.7^2 * tcrossprod(B)
    expect_equal(diag(V), rep(1.7^2, n), tolerance = 1e-14)
    expect_equal(V[1, 2], -1.7^2 / (n - 1L), tolerance = 1e-14)
    expect_equal(qr(B)$rank, n - 1L)
  }
  expect_equal(dim(isometric_rsm_basis(1)), c(1L, 0L))
})

test_that("normal prior and likelihood are invariant to a facet-label permutation", {
  # Transport an independently chosen point between the two level orderings.
  B <- isometric_rsm_basis(4)
  permutation <- c(3L, 1L, 4L, 2L)
  z <- c(.3, -1.2, .8)
  effect <- as.vector(B %*% z)
  transported <- as.vector(crossprod(B, effect[permutation]) * 3 / 4)
  expect_equal(sum(stats::dnorm(z, log = TRUE)), sum(stats::dnorm(transported, log = TRUE)), tolerance = 1e-14)
  reconstructed <- as.vector(B %*% transported)
  expect_equal(reconstructed, effect[permutation], tolerance = 1e-14)
  theta <- c(-1, .4, 1.2, -.2)
  old <- rsm_probabilities(theta - effect, c(-.6, .1, .5))
  new <- rsm_probabilities(theta - reconstructed[match(seq_len(4), permutation)], c(-.6, .1, .5))
  expect_equal(new, old, tolerance = 1e-14)
})

test_that("population priors and facet declarations are explicit and validated", {
  d <- expand.grid(person = 1:3, item = 1:2, rater = 1:2)
  d$score <- 1L
  s <- parse_bmfrm_formula(score ~ person + item + rater, d, K = 4)
  m <- bind_rsm_model(rsm_model(facet_sd = c(rater = .6, item = 1.2)), s)
  expect_identical(names(m$facet_sd), c("item", "rater"))
  expect_true(m$ability_mean$estimated && m$ability_sd$estimated)
  expect_false(rsm_model(ability_mean = 2, ability_sd = .7)$ability_mean$estimated)
  expect_error(rsm_model(ability_sd = 0), "Invalid")
  expect_error(rsm_model(ability_sd = "normal(0, 1)"), "lognormal")
  expect_error(rsm_model(ability_mean = "normal(0, sqrt(2))"), "literal")
  expect_error(rsm_model(ability_mean = "normal(sd=2, mean=0)"), "positional")
  expect_error(rsm_model(facet_sd = c(1, 2)), "named")
  expect_error(bind_rsm_model(rsm_model(facet_sd = c(item = 1)), s), "every")
  pcm <- s; pcm$family <- "pcm"
  expect_error(bind_rsm_model(rsm_model(), pcm), "additive RSM")
  expect_error(simulate_bmfrm(s$formula, d, 4, priors = list(), model = m), "not both")
})

test_that("prior generation uses the same hyperparameters and marginal level SDs", {
  d <- expand.grid(person = 1:2, item = 1:3, rater = 1:2)
  m <- rsm_model(ability_mean = "normal(1, 0.4)", ability_sd = "lognormal(-0.2, 0.3)",
                 facet_sd = c(item = .8, rater = 1.2), step_sd = .6)
  s <- prior_predict_bmfrm(score ~ person + item + rater, d, 4, model = m, ndraws = 3000, seed = 14)
  p <- s$parameters
  mu <- vapply(p, `[[`, numeric(1), "mu_theta")
  sigma <- vapply(p, `[[`, numeric(1), "sigma_theta")
  z <- vapply(p, function(x) (x$theta[1] - x$mu_theta) / x$sigma_theta, numeric(1))
  expect_equal(mean(mu), 1, tolerance = .04)
  expect_equal(stats::sd(mu), .4, tolerance = .04)
  expect_equal(mean(log(sigma)), -.2, tolerance = .04)
  expect_equal(stats::sd(log(sigma)), .3, tolerance = .04)
  expect_equal(mean(z), 0, tolerance = .07)
  expect_equal(stats::var(z), 1, tolerance = .12)
  expect_equal(stats::var(vapply(p, function(x) x$item[1], numeric(1))), .8^2, tolerance = .08)
  expect_equal(stats::var(vapply(p, function(x) x$rater[1], numeric(1))), 1.2^2, tolerance = .12)
  expect_equal(stats::var(vapply(p, function(x) x$tau[1], numeric(1))), .6^2, tolerance = .06)
  expect_true(all(vapply(p, function(x) max(abs(c(sum(x$item), sum(x$rater), sum(x$tau)))) < 1e-12, logical(1))))
})

test_that("binary and singleton blocks stay fixed and are excluded from MCMC diagnostics", {
  d <- data.frame(person = rep(1:5, 2), item = "only", rater = "only")
  s <- simulate_bmfrm(score ~ person + item + rater, d, 2,
                       model = rsm_model(ability_mean = 1, ability_sd = .7), seed = 5)
  expect_equal(s$parameters[[1]]$item, 0)
  expect_equal(s$parameters[[1]]$rater, 0)
  expect_equal(s$parameters[[1]]$tau, 0)
  expect_equal(s$parameters[[1]]$mu_theta, 1)
  expect_equal(s$parameters[[1]]$sigma_theta, .7)
  di <- prepare_data_bmfrm(s$spec, s$data)
  expect_identical(rsm_diagnostic_variables(s$spec, di), c("theta_raw", "theta"))
  expect_equal(dim(di$stan_data$B_tau), c(1L, 0L))
  expect_error(simulate_bmfrm(s$spec$formula, d, 2, model = s$model,
    parameters = within(s$parameters[[1]], sigma_theta <- -1)), "positive")
})

test_that("explicit unused levels survive generation, design checks, and Stan preparation", {
  d <- expand.grid(person = 1:5, item = 1:2, rater = 1:2)
  d$rater <- factor(d$rater, levels = c(1, 2, 3))
  s <- simulate_bmfrm(score ~ person + item + rater, d, 3, model = rsm_model(), seed = 8)
  di <- prepare_data_bmfrm(s$spec, s$data)
  dc <- check_design_bmfrm(s$spec$formula, s$data, 3)
  expect_identical(di$facet_labels, s$facet_labels)
  expect_equal(di$stan_data$J_rater, 3)
  expect_equal(length(s$parameters[[1]]$rater), 3)
  expect_false(dc$identifiable)
})

test_that("a bound person-only model can be reused without spurious raw parameters", {
  d <- data.frame(person = rep(1:5, 2))
  s <- simulate_bmfrm(score ~ person, d, 2, model = rsm_model(), seed = 1)
  again <- simulate_bmfrm(score ~ person, d, 2, model = s$model, seed = 1)
  expect_identical(s$y_rep, again$y_rep)
  expect_false(grepl("  _raw ~", build_stan_code(s$spec), fixed = TRUE))
})
