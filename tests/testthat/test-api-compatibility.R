test_that("Gaussian raw and marginal scales keep their distinct covariance semantics", {
  for (n in c(1L, 2L, 4L, 7L)) {
    projection <- diag(n) - matrix(1 / n, n, n)
    raw <- rsm_prior_basis(n, "raw")
    expect_equal(tcrossprod(raw), projection, tolerance = 1e-14)
    marginal <- rsm_prior_basis(n, "marginal")
    expected <- if (n == 1L) projection else projection * n / (n - 1L)
    expect_equal(tcrossprod(marginal), expected, tolerance = 1e-14)
  }
})

test_that("marginal Gaussian priors do not require hierarchical abilities", {
  d <- expand.grid(person = 1:5, item = 1:4, rater = 1:3)
  ps <- c(prior("normal(1, 2)", "theta"), prior("normal(0, .8)", "item", scale = "marginal"),
    prior("normal(3, .6)", "rater"), prior("normal(0, 1)", "tau"))
  sim <- simulate_bmfrm(score ~ person + item + rater, d, 4, priors = ps, seed = 7)
  expect_identical(sim$hierarchical, character())
  expect_setequal(names(sim$parameters[[1L]]), c("theta", "item", "rater", "tau"))
  spec <- parse_bmfrm_formula(sim$spec$formula, sim$data, 4)
  # Both low-level entry points resolve the same public configuration.
  expect_identical(prepare_data_bmfrm(spec, sim$data, ps)$stan_data,
                   prepare_data_bmfrm(sim$spec, sim$data)$stan_data)
  expect_identical(build_stan_code(spec, ps), build_stan_code(sim$spec))
  std <- prepare_data_bmfrm(spec, sim$data, ps)$stan_data
  expect_equal(diag(tcrossprod(std$B_item) * std$prior_sd_item^2), rep(.8^2, 4))
  expect_equal(diag(tcrossprod(std$B_rater) * std$prior_sd_rater^2), rep(.6^2 * 2/3, 3))
  pcm <- spec; pcm$family <- "pcm"; pcm$threshold_facet <- "item"
  expect_no_error(build_stan_code(pcm, ps))
  expect_no_error(prepare_data_bmfrm(pcm, sim$data, ps))
})

test_that("hierarchical defaults retain pre-centering facet and step SDs", {
  d <- expand.grid(person = 1:5, item = 1:4, rater = 1:3)
  sim <- simulate_bmfrm(score ~ person + item + rater, d, 4, hierarchical = "person", seed = 8)
  std <- prepare_data_bmfrm(sim$spec, sim$data)$stan_data
  expect_equal(diag(tcrossprod(std$B_item)), rep(3/4, 4))
  expect_equal(diag(tcrossprod(std$B_tau)), rep(2/3, 3))
  expect_identical(sim$priors, get_default_priors(hierarchical = "person"))
  expect_identical(sim$y_rep, simulate_bmfrm(sim$spec$formula, d, 4,
    priors = sim$priors, hierarchical = sim$hierarchical, seed = 8)$y_rep)
})

test_that("unused prior settings fail instead of silently changing the model", {
  d <- expand.grid(person = 1:5, item = 1:2, rater = 1:2)
  ps <- get_default_priors(hierarchical = character())
  expect_error(simulate_bmfrm(score ~ person + item + rater, d, 4,
    priors = ps, hierarchical = "person"), "not used")
  ps[[1]] <- prior("student_t(4, 0, 1)", "theta")
  ps[[2]] <- prior("normal(0, 1)", "item", scale = "marginal")
  # Mixed non-Gaussian + marginal priors are now supported in the unified pipeline
  mixed <- simulate_bmfrm(score ~ person + item + rater, d, 4, priors = ps, seed = 9)
  expect_match(build_stan_code(mixed$spec), "theta_raw ~ student_t(4, 0, 1);", fixed = TRUE)
  expect_error(bmfrm(score ~ person + item + rater, d, K = 4, model = list()), "removed")
  expect_false("rsm_model" %in% getNamespaceExports("BayesMFRM"))
})

test_that("API migration preserves frozen legacy and explicit marginal targets", {
  # Captured before this change from b0d059c, seed 37, without new code paths.
  old <- dget(testthat::test_path("fixtures", "rsm-api-b0d059c.R"))
  d <- expand.grid(person = 1:4, item = 1:3, rater = 1:2)
  legacy_default_ps <- c(prior("normal(0, 2)", "theta"), prior("normal(0, 1)", "item"),
    prior("normal(0, 1)", "rater"), prior("normal(0, 1)", "tau"))
  ps <- c(prior("normal(1, 1.2)", "theta"), prior("normal(2, 0.8)", "item"),
    prior("student_t(4, 1, 0.6)", "rater"), prior("normal(3, 0.7)", "tau"))
  iso <- c(prior("normal(0, 2.5)", "mean", facet = "person"), prior("lognormal(0, 0.5)", "sd", facet = "person"),
    prior("normal(0, .8)", "item", scale = "marginal"),
    prior("normal(0, 1.2)", "rater", scale = "marginal"),
    prior("normal(0, .6)", "tau", scale = "marginal"))
  generated <- list(
    default = simulate_bmfrm(score ~ person + item + rater, d, 4,
      priors = legacy_default_ps, hierarchical = character(), seed = 37, nsim = 3),
    custom = simulate_bmfrm(score ~ person + item + rater, d, 4, priors = ps, seed = 37, nsim = 3),
    isometric = simulate_bmfrm(score ~ person + item + rater, d, 4, priors = iso,
      hierarchical = "person", seed = 37, nsim = 3))
  for (name in names(generated)) {
    s <- generated[[name]]; dat <- d; dat$score <- s$y_rep[1L, ]
    expect_identical(s$y_rep, old[[name]]$y)
    # Continuous draws can differ by machine roundoff across math libraries.
    expect_equal(s$parameters, old[[name]]$parameters, tolerance = 1e-14)
    expect_identical(prepare_data_bmfrm(s$spec, dat)$stan_data, old[[name]]$prepared$data)
  }
  # New defaults use hierarchical = person_facet and weakly-informative tau ~ normal(0, 2)
  new_def <- simulate_bmfrm(score ~ person + item + rater, d, 4, seed = 37, nsim = 1)
  expect_identical(new_def$hierarchical, "person")
  expect_equal(prepare_data_bmfrm(new_def$spec, new_def$data)$stan_data$prior_sd_tau, 2)
})

test_that("ten facets use names without adding per-facet function arguments", {
  facets <- c("person", "item", "rater", paste0("f", 4:10))
  design <- expand.grid(stats::setNames(rep(list(1:2), 10), facets))
  formula <- stats::reformulate(facets, response = "score")
  s <- simulate_bmfrm(formula, design, 3, hierarchical = "person", seed = 11)
  expect_length(s$facet_labels, 10L)
  expect_equal(length(s$parameters[[1L]]$f10), 2L)
  expect_true(check_design_bmfrm(formula, s$data, 3)$identifiable)
  std <- prepare_data_bmfrm(s$spec, s$data)$stan_data
  expect_equal(std$J_f10, 2)
  expect_equal(dim(std$B_f10), c(2L, 1L))
  expect_true(grepl("vector[J_f10 - 1] f10_raw;", build_stan_code(s$spec), fixed = TRUE))
  expect_error(prior("normal(0, 2.5)", "mean", facet = "f10"), "no population mean")
  wrong_target <- c(s$priors, list(prior("normal(0, 1)", "sd", facet = "f10")))
  expect_error(simulate_bmfrm(formula, design, 3, hierarchical = "person", priors = wrong_target), "not used")
  s_hf <- simulate_bmfrm(formula, design, 3, hierarchical = c("person", "f10"), seed = 12)
  expect_true("sigma_f10" %in% names(s_hf$parameters[[1L]]))
  expect_match(build_stan_code(s_hf$spec), "sigma_f10_raw", fixed = TRUE)
  expect_error(simulate_bmfrm(formula, design, 3, hierarchical = "unknown"), "absent")
  expect_error(simulate_bmfrm(formula, design, 3, hierarchical = c("person", "person")), "unique")
})
