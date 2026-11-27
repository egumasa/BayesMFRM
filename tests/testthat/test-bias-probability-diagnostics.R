# Fixed probabilities, with genuinely retained iterations and chains.
make_bias_probability_fit <- function(p = .98, iterations = 1000L, chains = 4L,
                                      ordered = FALSE, delta = 0) {
  set.seed(614)
  x <- replicate(chains, {
    v <- c(rep(1 + delta, round(p * iterations)),
           rep(-1 - delta, iterations - round(p * iterations)))
    if (ordered) v else sample(v)
  })
  variables <- c("itemrater[1,1]", "itemrater[2,1]", "itemrater[1,2]", "itemrater[2,2]")
  a <- array(c(x, -x, -x, x), c(iterations, chains, 4L),
             dimnames = list(NULL, NULL, variables))
  # Deliberately permute variable storage order and label sort order.
  a <- a[, , c(3L, 1L, 4L, 2L), drop = FALSE]
  dat <- expand.grid(person = paste0("P", 1:3), item = c("Z", "A"), rater = c("Y", "B"))
  dat$score <- 1L
  spec <- parse_bmfrm_formula(score ~ person + item + rater + item:rater, dat, K = 3)
  structure(list(spec = spec,
    design = check_design_bmfrm(spec$formula, dat, K = 3),
    data_info = list(data_clean = dat,
      facet_labels = list(person = paste0("P", 1:3), item = c("Z", "A"), rater = c("Y", "B"))),
    convergence = list(converged = TRUE, status = "ok"),
    fit = list(draws = function(var, format = "matrix") {
      if (format == "draws_array") a else
        matrix(a, nrow = iterations * chains, dimnames = list(NULL, dimnames(a)[[3L]]))
    })), class = "bmfrm_fit")
}

test_that("probability MCSE uses original chains and preserves original bias summaries", {
  obj <- make_bias_probability_fit(delta = .3)
  out <- summarise_bias(obj, "item:rater", delta = .3)
  old <- summarise_bias(obj, "item:rater", delta = .3, diagnostics = FALSE)
  expect_equal(out[, names(old)], old)
  expect_equal(out$prob_thresh, rep(.95, 4))
  expect_equal(out$mcse_multiplier, rep(2, 4))
  a <- obj$fit$draws("itemrater", format = "draws_array")
  for (j in seq_len(nrow(out))) {
    x <- a[, , out$variable[j]]
    expect_equal(out$mcse_prob_gt_delta[j], posterior::mcse_mean((x > .3) * 1))
    expect_equal(out$mcse_prob_lt_neg_delta[j], posterior::mcse_mean((x < -.3) * 1))
  }
  expect_true(all(is.na(out$mcse_prob_outside_delta))) # all |bias| > delta
  expect_true(all(out$prob_threshold_stable))
  expect_true(all(out$prob_threshold_status == "resolved"))
  expect_true(all(is.na(out$decision_stable)))
  expect_true(all(out$decision_status == "bfdr_stability_not_assessed"))
  expect_equal(summarise_bias(obj, "rater:item", delta = .3), out)
})

test_that("autocorrelation changes MCSE even when posterior probabilities are identical", {
  independent <- summarise_bias(make_bias_probability_fit(p = .96), "item:rater")
  sticky <- summarise_bias(make_bias_probability_fit(p = .96, ordered = TRUE), "item:rater")
  expect_equal(sticky$prob_gt_delta, independent$prob_gt_delta)
  expect_gt(sticky$mcse_prob_gt_delta[1], 3 * independent$mcse_prob_gt_delta[1])
  expect_true(independent$prob_threshold_stable[1])
  expect_false(sticky$prob_threshold_stable[1])
  expect_identical(sticky$prob_threshold_status[1], "boundary_sensitive")
})

test_that("posterior draws_array subclasses retain chain layout during MCSE extraction", {
  obj <- make_bias_probability_fit()
  expected <- summarise_bias(obj, "item:rater")
  original <- obj$fit$draws
  obj$fit$draws <- function(var, format = "matrix") {
    if (format == "draws_array") posterior::as_draws_array(original(var, format)) else
      original(var, format)
  }
  out <- summarise_bias(obj, "item:rater")
  expect_equal(out, expected)
  x <- posterior::extract_variable_matrix(obj$fit$draws("itemrater", "draws_array"),
                                          out$variable[1L])
  expect_equal(out$mcse_prob_gt_delta[1L], posterior::mcse_mean((x > 0) * 1))
})

test_that("boundary equality is unresolved and diagnostic multiplier never changes flags", {
  obj <- make_bias_probability_fit(p = .95)
  out <- summarise_bias(obj, "item:rater")
  expect_true(all(out$prob_threshold_status == "boundary_sensitive"))
  expect_true(all(!out$prob_threshold_stable))
  expect_equal(out$flag_bias[1], "positive_bias")
  obj <- make_bias_probability_fit()
  ordinary <- summarise_bias(obj, "item:rater")
  conservative <- summarise_bias(obj, "item:rater", mcse_multiplier = 20)
  expect_equal(ordinary$flag_bias, conservative$flag_bias)
  expect_equal(ordinary$bfdr, conservative$bfdr)
  expect_true(all(ordinary$prob_threshold_stable))
  expect_true(all(!conservative$prob_threshold_stable))
  no_clear <- summarise_bias(obj, "item:rater", prob_thresh = .999)
  expect_true(all(no_clear$flag_bias %in% c("positive_bias", "negative_bias")))
  expect_true(all(no_clear$prob_threshold_status == "resolved"))
  expect_true(all(no_clear$decision_status == "bfdr_stability_not_assessed"))
  single_rule <- summarise_bias(obj, "item:rater", flag_method = "prob", prob_thresh = .999)
  expect_true(all(single_rule$flag_bias == "no_clear_bias"))
  expect_true(all(single_rule$decision_status == "stable_no_clear_bias"))
})

test_that("ROPE probabilities use delta and have their own MCSE", {
  set.seed(951)
  x <- matrix(sample(c(-1, 0, 1), 4000L, replace = TRUE), 1000L, 4L)
  a <- array(x, c(1000L, 4L, 1L), dimnames = list(NULL, NULL, "b"))
  ans <- bias_probability_diagnostics(a, "b", mean(x > .3), mean(x < -.3),
                                      .3, .95, 2, list(converged = TRUE))
  expect_equal(ans$mcse_prob_outside_delta, posterior::mcse_mean((abs(x) > .3) * 1))
  expect_equal(ans$mcse_prob_gt_delta, posterior::mcse_mean((x > .3) * 1))
  expect_equal(ans$mcse_prob_lt_neg_delta, posterior::mcse_mean((x < -.3) * 1))
  expect_true(ans$decision_stable)
  expect_identical(ans$decision_status, "stable_no_clear_bias")
})

test_that("constant indicators and stuck chains never imply zero MC error", {
  constant <- summarise_bias(make_bias_probability_fit(p = 1), "item:rater")
  expect_true(all(is.na(constant$mcse_prob_gt_delta)))
  expect_true(all(is.na(constant$prob_threshold_stable)))
  expect_true(all(constant$prob_threshold_status == "mcse_unavailable"))
  x <- cbind(rep(1, 100), rep(-1, 100))
  a <- array(x, c(100L, 2L, 1L), dimnames = list(NULL, NULL, "b"))
  ans <- bias_probability_diagnostics(a, "b", .5, .5, 0, .95, 2, list(converged = TRUE))
  expect_true(is.na(ans$mcse_prob_gt_delta))
  expect_true(is.na(ans$decision_stable))
  expect_identical(ans$decision_status, "mcse_unavailable")
  # A single constant chain still makes this event's MCSE unavailable.
  x[, 2L] <- rep(c(-1, 1), 50)
  a[] <- x
  ans <- bias_probability_diagnostics(a, "b", .75, .25, 0, .95, 2, list(converged = TRUE))
  expect_true(is.na(ans$mcse_prob_gt_delta))
  expect_true(is.na(ans$decision_stable))
})

test_that("sampler checks and chain availability gate stability, not substantive flags", {
  obj <- make_bias_probability_fit()
  expected <- summarise_bias(obj, "item:rater")
  for (cv in list(list(converged = FALSE, status = "warning"),
                 list(converged = FALSE, status = "incomplete"), NULL)) {
    obj$convergence <- cv
    out <- summarise_bias(obj, "item:rater")
    expect_true(all(is.na(out$prob_threshold_stable)))
    expect_equal(out$flag_bias, expected$flag_bias)
    expect_equal(out$mcse_prob_gt_delta, expected$mcse_prob_gt_delta)
    expect_true(all(out$prob_threshold_status == if (identical(cv$status, "warning"))
      "fit_diagnostics_failed" else "fit_diagnostics_unavailable"))
  }
  single <- summarise_bias(make_bias_probability_fit(chains = 1L), "item:rater")
  expect_true(all(single$prob_threshold_status == "fit_diagnostics_unavailable"))
  expect_true(all(is.na(single$prob_threshold_stable)))
  mat <- obj$fit$draws("itemrater", "matrix")
  obj$fit$draws <- function(...) mat
  unknown <- summarise_bias(obj, "item:rater")
  expect_true(all(unknown$prob_threshold_status == "chain_structure_unavailable"))
  expect_true(all(is.na(unknown$prob_threshold_stable)))
})

test_that("numerical stability remains distinct from likelihood identification", {
  obj <- make_bias_probability_fit()
  rows <- obj$design$interaction_levels[["item:rater"]]$map_row
  # Simulate the design check marking these contrasts unidentifiable.
  obj$design$effect_null_space <- matrix(0, nrow(obj$design$effect_null_space), 1L)
  obj$design$effect_null_space[rows, 1L] <- 1
  out <- summarise_bias(obj, "item:rater")
  expect_true(all(!out$likelihood_identified))
  expect_true(all(out$flag_bias == "not_likelihood_identified"))
  expect_true(all(is.na(out$bfdr)))
  expect_true(all(out$prob_threshold_stable))
  obj$design <- NULL
  for (method in c("bfdr", "prob")) {
    unchecked <- summarise_bias(obj, "item:rater", flag_method = method)
    expect_true(all(is.na(unchecked$likelihood_identified)))
    expect_true(all(is.na(unchecked$bfdr)))
    expect_true(all(unchecked$flag_bias == "identification_unchecked"))
    expect_true(all(unchecked$decision_status == "identification_unchecked"))
    expect_true(all(is.na(unchecked$decision_stable)))
    expect_equal(unchecked$mcse_prob_gt_delta, out$mcse_prob_gt_delta)
    expect_true(all(unchecked$prob_threshold_stable))
  }
})

test_that("legacy and generated summaries share MCSE and maintain cell alignment", {
  obj <- make_bias_probability_fit()
  a <- obj$fit$draws("itemrater", "draws_array")
  dimnames(a)[[3L]] <- sub("itemrater", "bias", dimnames(a)[[3L]])
  legacy <- list(draws = function(...) posterior::as_draws_array(a),
                 metadata = function() list(stan_variables = "bias"),
                 summary = function(...) data.frame(rhat = 1, ess_bulk = 2000,
                   ess_tail = 2000, mcse_mean = .01),
                 diagnostic_summary = function(...) list(ebfmi = rep(.8, 4),
                   num_divergent = rep(0, 4), num_max_treedepth = rep(0, 4)),
                 num_chains = function() 4L)
  out <- summarise_bias_ir(legacy, c("Z", "A"), c("Y", "B"))
  modern <- summarise_bias(obj, "item:rater", flag_method = "prob")
  modern <- modern[order(modern$item, modern$rater), ]
  expect_equal(out$flag_bias, modern$flag_bias)
  expect_equal(out$mcse_prob_pos, modern$mcse_prob_gt_delta)
  expect_equal(out$mcse_prob_neg, modern$mcse_prob_lt_neg_delta)
  expect_equal(out$decision_status, modern$decision_status)
  old <- summarise_bias_ir(legacy, c("Z", "A"), c("Y", "B"), diagnostics = FALSE)
  expect_equal(out[, names(old)], old)
  legacy$metadata <- NULL
  expect_true(all(summarise_bias_ir(legacy)$decision_status == "fit_diagnostics_unavailable"))
})

test_that("diagnostic arguments fail early and diagnostics can be skipped", {
  obj <- make_bias_probability_fit()
  for (bad in list(NA, 1, "yes", c(TRUE, FALSE))) {
    expect_error(summarise_bias(obj, "item:rater", diagnostics = bad), "diagnostics")
    expect_error(summarise_bias_ir(NULL, diagnostics = bad), "diagnostics")
  }
  for (bad in list(0, -1, Inf, NA, "2", c(1, 2))) {
    expect_error(summarise_bias(obj, "item:rater", mcse_multiplier = bad), "mcse_multiplier")
    expect_error(summarise_bias_ir(NULL, mcse_multiplier = bad), "mcse_multiplier")
  }
  draw_fn <- obj$fit$draws
  obj$fit$draws <- function(var, format = "matrix") {
    if (format == "draws_array") stop("Chain array must not be requested")
    draw_fn(var, format)
  }
  out <- summarise_bias(obj, "item:rater", diagnostics = FALSE)
  expect_false("decision_stable" %in% names(out))
})

test_that("bias interval plots use distinct symbols for screened and withheld estimates", {
  obj <- make_bias_probability_fit()
  tab <- summarise_bias(obj, "item:rater", diagnostics = FALSE)
  tab$flag_bias[2L] <- "no_clear_bias"
  # Old saved tables can contain a signed flag despite unknown identification.
  tab$likelihood_identified[3L] <- NA
  tab$n_observations[4L] <- 0L
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  symbols <- NULL
  testthat::with_mocked_bindings(
    points.default = function(x, y, pch, ...) symbols <<- pch,
    expect_identical(plot_bias(tab), tab),
    .package = "graphics"
  )
  expect_equal(symbols, c(19L, 1L, 4L, 4L))
  tab$delta[1L] <- 0.5
  expect_error(plot_bias(tab), "one set of screening settings")
})
