# Test fixture for complete prior lists, not a public configuration constructor.
rsm_test_priors <- function(mean = "normal(0, 2.5)", sd = "lognormal(0, 0.5)",
                            item = 1, rater = 1, step = 1, scale = "marginal") {
  c(prior(mean, "mean", facet = "person"), prior(sd, "sd", facet = "person"),
    prior(sprintf("normal(0, %.17g)", item), "item", scale = scale),
    prior(sprintf("normal(0, %.17g)", rater), "rater", scale = scale),
    prior(sprintf("normal(0, %.17g)", step), "tau", scale = scale))
}

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

test_that("population priors and facet declarations share the existing prior API", {
  d <- expand.grid(person = 1:3, item = 1:2, rater = 1:2)
  d$score <- 1L
  s <- parse_bmfrm_formula(score ~ person + item + rater, d, K = 4, hierarchical = "person")
  resolved <- resolve_bmfrm_spec(s, rsm_test_priors(item = 1.2, rater = .6))
  m <- resolved$gaussian_model
  expect_identical(names(m$facet_sd), c("item", "rater"))
  expect_true(m$ability_mean$estimated && m$ability_sd$estimated)
  expect_false(resolve_bmfrm_spec(s, rsm_test_priors(mean = 2, sd = .7))$gaussian_model$ability_mean$estimated)
  expect_error(prior(0, "sd", facet = "person"), "Invalid")
  expect_error(prior("normal(0, 1)", "sd", facet = "person"), "lognormal")
  expect_error(prior("normal(0, sqrt(2))", "mean", facet = "person"), "literal")
  expect_error(prior("normal(sd=2, mean=0)", "mean", facet = "person"), "positional")
  expect_error(prior("student_t(3, 0, 1)", "item", scale = "marginal"), "normal")
  expect_error(prior("normal(1, 1)", "item", scale = "marginal"), "normal")
  expect_error(prior("normal(0, 1)", "theta", scale = "marginal"), "constrained")
  expect_error(resolve_bmfrm_spec(s, rsm_test_priors()[-3]), "Missing")
  expect_error(resolve_bmfrm_spec(s, c(rsm_test_priors(), list(prior(0, "mean", facet = "person")))), "Duplicate")
  expect_error(resolve_bmfrm_spec(s, c(rsm_test_priors(), list(prior("normal(0, 2)", "theta")))), "not used")
  pcm <- s; pcm$family <- "pcm"; pcm$threshold_facet <- "item"
  expect_no_error(resolve_bmfrm_spec(pcm, rsm_test_priors()))
  indep <- s; indep$hierarchical <- character()
  expect_error(resolve_bmfrm_spec(indep, rsm_test_priors()), "not used")
})

test_that("prior generation uses the same hyperparameters and marginal level SDs", {
  d <- expand.grid(person = 1:2, item = 1:3, rater = 1:2)
  ps <- rsm_test_priors(mean = "normal(1, 0.4)", sd = "lognormal(-0.2, 0.3)",
                        item = .8, rater = 1.2, step = .6)
  s <- prior_predict_bmfrm(score ~ person + item + rater, d, 4, hierarchical = "person", priors = ps, ndraws = 3000, seed = 14)
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
                       hierarchical = "person", priors = rsm_test_priors(mean = 1, sd = .7), seed = 5)
  expect_equal(s$parameters[[1]]$item, 0)
  expect_equal(s$parameters[[1]]$rater, 0)
  expect_equal(s$parameters[[1]]$tau, 0)
  expect_equal(s$parameters[[1]]$mu_theta, 1)
  expect_equal(s$parameters[[1]]$sigma_theta, .7)
  di <- prepare_data_bmfrm(s$spec, s$data)
  expect_identical(model_diagnostic_variables(s$spec, di), c("theta_raw", "theta"))
  expect_equal(dim(di$stan_data$B_tau), c(1L, 0L))
  expect_error(simulate_bmfrm(s$spec$formula, d, 2, hierarchical = s$hierarchical, priors = s$priors,
    parameters = within(s$parameters[[1]], sigma_theta <- -1)), "positive")
})

test_that("explicit unused levels survive generation, design checks, and Stan preparation", {
  d <- expand.grid(person = 1:5, item = 1:2, rater = 1:2)
  d$rater <- factor(d$rater, levels = c(1, 2, 3))
  s <- simulate_bmfrm(score ~ person + item + rater, d, 3, hierarchical = "person", seed = 8)
  di <- prepare_data_bmfrm(s$spec, s$data)
  dc <- check_design_bmfrm(s$spec$formula, s$data, 3)
  expect_identical(di$facet_labels, s$facet_labels)
  expect_equal(di$stan_data$J_rater, 3)
  expect_equal(length(s$parameters[[1]]$rater), 3)
  expect_false(dc$identifiable)
})

test_that("person-only priors can be reused without spurious raw parameters", {
  d <- data.frame(person = rep(1:5, 2))
  s <- simulate_bmfrm(score ~ person, d, 2, hierarchical = "person", seed = 1)
  again <- simulate_bmfrm(score ~ person, d, 2, hierarchical = s$hierarchical, priors = s$priors, seed = 1)
  expect_identical(s$y_rep, again$y_rep)
  expect_false(grepl("  _raw ~", build_stan_code(s$spec), fixed = TRUE))
})
