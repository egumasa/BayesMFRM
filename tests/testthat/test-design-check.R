test_that("connectedness does not conceal item-rater aliasing", {
  d <- expand.grid(person = 1:6, item = 1:3)
  d$rater <- d$item
  x <- check_design_bmfrm(score ~ person + item + rater, d, K = 6)
  expect_true(x$connected)
  expect_false(x$identifiable)
  expect_equal(x$rank, 8L)
  expect_equal(x$n_parameters, 10L)
  expect_equal(x$nullity, 2L)
  expect_lt(max(abs(x$design_matrix %*% x$null_space)), 1e-10)
})

test_that("complete RSM designs are identified and row order is irrelevant", {
  d <- expand.grid(person = 1:6, item = 1:3, rater = 1:2)
  x <- check_design_bmfrm(score ~ person + item + rater, d, K = 6)
  y <- check_design_bmfrm(score ~ person + item + rater, d[nrow(d):1, ], K = 6)
  expect_true(x$connected)
  expect_true(x$identifiable)
  expect_equal(x$rank, 9L)
  expect_equal(x$adjacent_logit_rank, 13L)
  expect_equal(x$rank, y$rank)
  expect_equal(x$condition_number, y$condition_number, tolerance = 1e-10)
})

test_that("planned and observed designs can be checked separately", {
  d <- expand.grid(person = 1:6, item = 1:3, rater = 1:3)
  d$score <- ifelse(d$rater == d$item, 3L, NA_integer_)
  planned <- check_design_bmfrm(score ~ person + item + rater, d, K = 6)
  observed <- check_design_bmfrm(score ~ person + item + rater, d, K = 6, observed_only = TRUE)
  expect_true(planned$identifiable)
  expect_false(observed$identifiable)
  expect_equal(observed$n_observations, 18L)
  expect_equal(observed$n_dropped, 36L)
})

test_that("disconnected components and unused declared levels are reported", {
  d <- rbind(expand.grid(person = 1:3, item = 1:2, rater = 1:2),
             expand.grid(person = 4:6, item = 3:4, rater = 3:4))
  x <- check_design_bmfrm(score ~ person + item + rater, d, K = 6)
  expect_equal(x$n_components, 2L)
  expect_false(x$identifiable)
  d$person <- factor(d$person, levels = 1:7)
  y <- check_design_bmfrm(score ~ person + item + rater, d, K = 6)
  expect_equal(y$n_components, 3L)
  expect_true(any(y$levels$n_observations == 0))
})

test_that("unsupported models and invalid designs do not receive an identification verdict", {
  d <- expand.grid(person = 1:4, item = 1:2, rater = 1:2)
  expect_error(check_design_bmfrm(score ~ person + item + rater:item, d, K = 6), "additive RSM")
  expect_error(check_design_bmfrm(score ~ person + item, d, K = 6, family = "pcm"), "threshold_facet")
  expect_error(check_design_bmfrm(score ~ person + item, d, K = 1), "K")
  d$person[1] <- NA
  expect_error(check_design_bmfrm(score ~ person + item, d, K = 6), "Missing")
})

test_that("singleton centered facets and binary thresholds add no free coordinates", {
  d <- data.frame(person = 1:4, item = 1, rater = 1)
  x <- check_design_bmfrm(score ~ person + item + rater, d, K = 2)
  expect_true(x$identifiable)
  expect_equal(x$n_parameters, 4L)
  expect_equal(x$adjacent_logit_rank, 4L)
})
