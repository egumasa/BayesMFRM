test_that("the previous PCM constraints allow an item-step alias", {
  beta <- c(-.5, .2, .3)
  tau <- matrix(c(-.2, .1, .1, .3, -.1, -.2, -.1, 0, .1), 3)
  shift <- c(-.4, .1, .3)
  shifted <- sweep(tau, 1, shift, "-")
  expect_equal(sum(beta + shift), 0)
  expect_equal(colSums(shifted), rep(0, 3), tolerance = 1e-14)
  eta <- .4 - beta
  expect_equal(pcm_probabilities(eta, tau), pcm_probabilities(eta - shift, shifted), tolerance = 1e-14)
  # That shift is excluded by the corrected within-level step constraint.
  centered <- tau - rowMeans(tau)
  expect_gt(max(abs(rowSums(sweep(centered, 1, shift, "-")))), .5)
})

test_that("PCM adjacent-logit rank agrees with an independently assembled matrix", {
  d <- expand.grid(person = 1:4, item = 1:3, rater = 1:2)
  x <- check_design_bmfrm(score ~ person + item + rater, d, 4, family = "pcm", threshold_facet = "item")
  S <- 3L
  incidence <- diag(3)[d$item, ]
  step_basis <- stats::contr.helmert(S)
  full <- cbind(x$design_matrix[rep(seq_len(nrow(d)), each = S), ], -kronecker(incidence, step_basis))
  expect_equal(qr(full)$rank, ncol(full))
  expect_equal(x$adjacent_logit_rank, qr(full)$rank)
  expect_equal(x$adjacent_logit_parameters, ncol(full))
  old <- cbind(x$design_matrix[rep(seq_len(nrow(d)), each = S), ],
               -kronecker(incidence %*% stats::contr.helmert(3), diag(S)))
  expect_equal(ncol(old) - qr(old)$rank, 2L)
  # Preserve an unobserved planned level and account for its step contrasts.
  d$item <- factor(d$item, levels = 1:4)
  y <- check_design_bmfrm(score ~ person + item + rater, d, 4, family = "pcm", threshold_facet = "item")
  expect_equal(y$step_rank, 6)
  expect_equal(y$step_parameters, 8)
  expect_identical(y$unobserved_step_levels, "4")
  expect_false(y$identifiable)
  full <- cbind(y$design_matrix[rep(seq_len(nrow(d)), each = S), ],
                -kronecker(diag(4)[as.integer(d$item), ], step_basis))
  expect_equal(y$adjacent_logit_rank, qr(full)$rank)
})

test_that("PCM step SD refers to within-level constraints and is independent of item count", {
  for (J in c(1L, 3L, 7L)) {
    d <- expand.grid(person = 1:12, task = seq_len(J))
    for (scale in c("raw", "marginal")) {
      ps <- c(prior("normal(0, 2)", "theta"), prior("normal(0, 1)", "facet_task"),
               prior("normal(0, .7)", "tau", scale = scale))
      s <- simulate_bmfrm(score ~ person + task, d, 4, family = "pcm", threshold_facet = "task", priors = ps, seed = 14)
      expect_equal(dim(s$parameters[[1]]$tau), c(J, 3L))
      expect_equal(rowSums(s$parameters[[1]]$tau), rep(0, J), tolerance = 1e-14)
      info <- prepare_data_bmfrm(s$spec, s$data)
      V <- info$stan_data$prior_sd_tau^2 * tcrossprod(info$stan_data$B_tau)
      expected <- if (scale == "raw") .7^2 * 2/3 else .7^2
      expect_equal(diag(V), rep(expected, 3))
      expect_equal(V[1, 2], -expected/2)
    }
  }
  d <- expand.grid(person = 1:4, item = 1:3)
  ps <- c(prior("normal(0, 2)", "theta"), prior("normal(0, 1)", "item"),
           prior("normal(0, .7)", "tau", scale = "marginal"))
  s <- prior_predict_bmfrm(score ~ person + item, d, 4, priors = ps,
    family = "pcm", threshold_facet = "item", ndraws = 2500, seed = 72)
  x <- vapply(s$parameters, function(p) p$tau[1, 1], numeric(1))
  y <- vapply(s$parameters, function(p) p$tau[1, 2], numeric(1))
  z <- vapply(s$parameters, function(p) p$tau[2, 1], numeric(1))
  expect_equal(stats::var(x), .49, tolerance = .05)
  expect_equal(stats::cov(x, y), -.245, tolerance = .04)
  expect_lt(abs(stats::cov(x, z)), .04)
})

test_that("binary PCM is the binary RSM and singleton items retain nonbinary steps", {
  d <- expand.grid(person = 1:4, item = 1:3)
  pars <- list(theta = c(-1, -.2, .3, .9), item = c(-.4, .1, .3), tau = matrix(0, 3, 1))
  pcm <- simulate_bmfrm(score ~ person + item, d, 2, parameters = pars, family = "pcm", threshold_facet = "item", seed = 91)
  pars$tau <- 0
  rsm <- simulate_bmfrm(score ~ person + item, d, 2, parameters = pars, seed = 91)
  expect_identical(pcm$y_rep, rsm$y_rep)
  expect_false("tau_raw" %in% model_diagnostic_variables(pcm$spec, prepare_data_bmfrm(pcm$spec, pcm$data)))
  d$item <- "one"
  single <- simulate_bmfrm(score ~ person + item, d, 4, family = "pcm", threshold_facet = "item", seed = 3)
  expect_gt(max(abs(single$parameters[[1]]$tau)), .1)
})

test_that("fixed PCM step matrices validate labels and every row constraint", {
  d <- expand.grid(person = 1:4, task = c("A", "B", "C"))
  pars <- list(theta = c(-1, -.2, .3, .9), task = c(-.4, .1, .3),
    tau = matrix(c(-1, -.5, -.7, 0, .2, .1, 1, .3, .6), 3, dimnames = list(c("A", "B", "C"), 1:3)))
  a <- simulate_bmfrm(score ~ person + task, d, 4, parameters = pars, family = "pcm", threshold_facet = "task", seed = 5)
  pars$tau <- pars$tau[3:1, 3:1]
  b <- simulate_bmfrm(score ~ person + task, d, 4, parameters = pars, family = "pcm", threshold_facet = "task", seed = 5)
  expect_identical(a$y_rep, b$y_rep)
  pars$tau <- matrix(c(1, -1, 0), 3, 3)
  expect_error(simulate_bmfrm(score ~ person + task, d, 4, parameters = pars, family = "pcm", threshold_facet = "task"), "within every")
})

test_that("raw Student-t PCM priors generate the same row-centered step structure", {
  d <- expand.grid(person = 1:4, item = 1:3)
  ps <- c(prior("student_t(4, 0, 2)", "theta"), prior("normal(0, 1)", "item"), prior("student_t(4, 0, .5)", "tau"))
  s <- simulate_bmfrm(score ~ person + item, d, 4, family = "pcm", threshold_facet = "item", priors = ps, seed = 17)
  expect_null(s$spec$gaussian_model)
  expect_equal(rowSums(s$parameters[[1]]$tau), rep(0, 3), tolerance = 1e-14)
  expect_match(build_stan_code(s$spec), "to_vector(tau_raw) ~ student_t(4, 0, .5)", fixed = TRUE)
})

test_that("fair scores cannot silently leave adjusted-facet PCM steps in place", {
  x <- structure(list(spec = list(family = "pcm", threshold_facet = "rater")), class = "bmfrm_fit")
  expect_error(fair_scores(x), "reference step profile")
})
