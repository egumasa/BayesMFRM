# Tests for user-facing UX, diagnostic clarity, K inference messaging, and
# end-to-end verification of the Quick Start and getting-started vignette workflow.

test_that("parse_bmfrm_formula notifies explicitly when K is inferred from max(score)", {
  d <- expand.grid(person = paste0("P", 1:6), item = c("I1", "I2"), rater = c("R1", "R2"), stringsAsFactors = FALSE)
  d$score <- rep(1:4, length.out = nrow(d))

  # Omitting K must inform the user that K was inferred from max(score)
  expect_message(
    spec_inf <- parse_bmfrm_formula(score ~ person + item + rater, d),
    "inferred K = 4"
  )
  expect_equal(spec_inf$K, 4L)

  # Passing explicit rubric upper bound K = 6 keeps K = 6 without an inference message
  expect_no_message(
    spec_exp <- parse_bmfrm_formula(score ~ person + item + rater, d, K = 6)
  )
  expect_equal(spec_exp$K, 6L)
})

test_that("print() and summary() separate design identification from MCMC convergence and cite bayesplot::pp_check()", {
  # Connected but rank-deficient design (item <-> rater confounded)
  d_bad <- data.frame(
    person = rep(c("P1", "P2", "P3", "P4"), each = 2),
    item   = rep(c("I1", "I2"), 4),
    rater  = rep(c("R1", "R2"), 4),
    score  = c(2L, 3L, 2L, 3L, 3L, 2L, 3L, 2L),
    stringsAsFactors = FALSE
  )
  spec <- parse_bmfrm_formula(score ~ person + item + rater, d_bad, K = 4, save_gq = "minimal")
  di <- suppressWarnings(prepare_data_bmfrm(spec, d_bad))
  dc <- check_design_bmfrm(spec$formula, d_bad, K = 4)
  expect_true(dc$connected)
  expect_false(dc$identifiable)

  D <- 10L
  draws_map <- list(
    theta = matrix(0, D, 4, dimnames = list(NULL, paste0("theta[", 1:4, "]"))),
    item  = matrix(0, D, 2, dimnames = list(NULL, paste0("item[", 1:2, "]"))),
    rater = matrix(0, D, 2, dimnames = list(NULL, paste0("rater[", 1:2, "]"))),
    tau   = matrix(0, D, 3, dimnames = list(NULL, paste0("tau[", 1:3, "]")))
  )
  mock_fit <- list(
    num_chains = function() 4L,
    metadata = function() list(iter_sampling = 1000L, iter_warmup = 1000L, stan_variables = names(draws_map)),
    draws = function(var, format = "matrix") draws_map[[var]],
    summary = function(var, ...) {
      m <- draws_map[[var]]
      posterior::summarise_draws(posterior::as_draws_matrix(m), ...)
    }
  )
  conv_ok <- list(
    converged = TRUE, max_rhat = 1.001, min_ess_bulk = 1200, min_ess_tail = 900,
    min_ebfmi = 0.9, max_mcse_mean = 0.01, n_divergent = 0L, n_max_treedepth = 0L
  )
  obj_bad <- structure(
    list(fit = mock_fit, spec = spec, data_info = di, design = dc, convergence = conv_ok),
    class = "bmfrm_fit"
  )

  # Both views flag rank deficiency even with stable chains. The overview uses no draws.
  stored <- obj_bad
  stored$fit <- NULL
  compact <- paste(capture.output(returned <- withVisible(print(stored))), collapse = "\n")
  expect_identical(returned, list(value = stored, visible = FALSE))
  expect_match(compact, "Design: Rank deficient", fixed = TRUE)
  expect_match(compact, "MCMC: Passed", fixed = TRUE)
  stored$convergence <- NULL
  expect_true(any(grepl("MCMC: Unchecked", capture.output(print(stored)), fixed = TRUE)))
  stored$convergence <- list(status = "incomplete", converged = FALSE)
  expect_true(any(grepl("MCMC: Incomplete", capture.output(print(stored)), fixed = TRUE)))
  txt_print <- paste(capture.output(print(obj_bad, details = TRUE)), collapse = "\n")
  expect_match(txt_print, "Design Identification:", fixed = TRUE)
  expect_match(txt_print, "WARNING: Rank-deficient design", fixed = TRUE)
  expect_match(txt_print, "does NOT resolve the design rank deficiency", fixed = TRUE)
  expect_match(txt_print, "bayesplot::pp_check(", fixed = TRUE)

  txt_sum <- paste(capture.output(summary(obj_bad)), collapse = "\n")
  expect_match(txt_sum, "Design Identification:", fixed = TRUE)
  expect_match(txt_sum, "WARNING - rank-deficient design", fixed = TRUE)

  # Identifiable design must display OK for identification + scope note for convergence
  d_ok <- expand.grid(person = paste0("P", 1:4), item = c("I1", "I2"), rater = c("R1", "R2"), stringsAsFactors = FALSE)
  d_ok$score <- rep(1:4, length.out = nrow(d_ok))
  dc_ok <- check_design_bmfrm(score ~ person + item + rater, d_ok, K = 4)
  obj_ok <- obj_bad
  obj_ok$design <- dc_ok
  txt_ok <- paste(capture.output(print(obj_ok, details = TRUE)), collapse = "\n")
  expect_match(txt_ok, "OK: Conditional likelihood is structurally identifiable", fixed = TRUE)
  expect_match(txt_ok, "sampling stability only, not model fit to data", fixed = TRUE)
})

test_that("Quick Start and getting-started tutorial workflow runs cleanly end-to-end via public APIs", {
  # Exact small bridged plan from README Quick Start & vignettes/getting-started.Rmd
  set.seed(2026)
  full_grid <- expand.grid(
    person = sprintf("P%02d", 1:24),
    item   = sprintf("I%02d", 1:4),
    rater  = sprintf("R%02d", 1:4),
    stringsAsFactors = FALSE
  )
  plan <- do.call(rbind, lapply(split(full_grid, full_grid$person), function(df) {
    df[sample(nrow(df), 8L), , drop = FALSE]
  }))
  rownames(plan) <- NULL

  # 1. Design check
  dc <- check_design_bmfrm(score ~ person + item + rater, plan, K = 6)
  expect_true(dc$connected)
  expect_true(dc$identifiable)

  # 2. Prior predictive check: moderate vs wide SD must demonstrate extreme-category concentration
  priors_mod <- c(
    prior("normal(0, 0.5)", "mean", facet = "person"),
    prior("lognormal(-0.5, 0.3)", "sd", facet = "person"),
    prior("normal(0, 0.5)", "item", scale = "marginal"),
    prior("normal(0, 0.5)", "rater", scale = "marginal"),
    prior("normal(0, 0.6)", "tau", scale = "marginal")
  )
  priors_wide <- c(
    prior("normal(0, 3.0)", "mean", facet = "person"),
    prior("lognormal(1.0, 0.5)", "sd", facet = "person"),
    prior("normal(0, 3.0)", "item", scale = "marginal"),
    prior("normal(0, 3.0)", "rater", scale = "marginal"),
    prior("normal(0, 3.0)", "tau", scale = "marginal")
  )
  pp_mod  <- prior_predict_bmfrm(score ~ person + item + rater, plan, K = 6,
                                 hierarchical = "person", priors = priors_mod, ndraws = 60, seed = 101)
  pp_wide <- prior_predict_bmfrm(score ~ person + item + rater, plan, K = 6,
                                 hierarchical = "person", priors = priors_wide, ndraws = 60, seed = 101)
  ext_mod  <- mean(pp_mod$y_rep %in% c(1L, 6L))
  ext_wide <- mean(pp_wide$y_rep %in% c(1L, 6L))
  expect_gt(ext_wide, ext_mod + 0.15)

  # 3. Simulated dataset + suitability check
  sim <- simulate_bmfrm(score ~ person + item + rater, design = plan, K = 6,
                        hierarchical = "person", priors = priors_mod, seed = 202)
  ratings <- sim$data
  suit <- suppressMessages(check_data_suitability(ratings, "score", c("person", "item", "rater")))
  expect_s3_class(suit, "bmfrm_data_check")
  expect_equal(suit$n_obs, nrow(ratings))
  expect_equal(length(suit$warnings), 0L)

  # 4. Verify all downstream public S3 post-processing (including bayesplot::pp_check)
  spec <- parse_bmfrm_formula(score ~ person + item + rater, ratings, K = 6,
                              hierarchical = "person", save_gq = "minimal")
  spec <- resolve_bmfrm_spec(spec, priors_mod)
  di <- prepare_data_bmfrm(spec, ratings, priors_mod)
  p1 <- sim$parameters[[1L]]
  D <- 8L
  draws_map <- list(
    theta = matrix(rep(p1$theta, each = D), D, 24, dimnames = list(NULL, paste0("theta[", 1:24, "]"))),
    item  = matrix(rep(p1$item,  each = D), D,  4, dimnames = list(NULL, paste0("item[",  1:4,  "]"))),
    rater = matrix(rep(p1$rater, each = D), D,  4, dimnames = list(NULL, paste0("rater[", 1:4,  "]"))),
    tau   = matrix(rep(p1$tau,   each = D), D,  5, dimnames = list(NULL, paste0("tau[",   1:5,  "]")))
  )
  mock_fit <- list(
    num_chains = function() 4L,
    metadata = function() list(iter_sampling = 1000L, iter_warmup = 1000L, stan_variables = names(draws_map)),
    draws = function(var, format = "matrix") draws_map[[var]]
  )
  obj <- structure(list(fit = mock_fit, spec = spec, data_info = di, design = dc), class = "bmfrm_fit")

  ppc_plot <- suppressMessages(bayesplot::pp_check(obj, group = "rater", ndraws = 8, seed = 1))
  expect_true(inherits(ppc_plot, "ggplot"))

  res_tbl <- residuals(obj)
  expect_equal(nrow(res_tbl), nrow(ratings))

  rf <- facet_fit(obj, "rater")
  expect_equal(nrow(rf), 4L)

  rc <- residual_cor_criteria(res_tbl, unit_facets = c("person", "rater"))
  expect_true(inherits(rc, "residual_cor"))

  fs <- fair_scores(obj, adjust_facets = "rater", standardize = TRUE)
  expect_equal(nrow(fs), 24L)
  expect_true(all(fs$likelihood_identified))
})
