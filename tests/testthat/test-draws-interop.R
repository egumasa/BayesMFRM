test_that("posterior and tidybayes preserve rater, iteration and chain indices", {
  raw <- posterior::as_draws_array(array(seq_len(12), dim = c(3, 2, 2),
    dimnames = list(NULL, NULL, c("rater[1]", "rater[2]"))))
  fit <- structure(list(fit = list(draws = function(variables = NULL) {
    posterior::subset_draws(raw, variable = variables)
  })), class = "bmfrm_fit")
  expect_identical(posterior::as_draws(fit, variables = "rater[2]"),
                   posterior::subset_draws(raw, variable = "rater[2]"))
  expect_identical(posterior::as_draws_df(fit), posterior::as_draws_df(raw))

  skip_if_not_installed("tidybayes")
  expect_identical(tidybayes::spread_draws(fit, rater[r]), tidybayes::spread_draws(raw, rater[r]))
})

test_that("bayestestR uses the selected parameter draws for direction and full ROPE mass", {
  skip_if_not_installed("bayestestR")
  b <- c(-0.3, 0.2, 0.4, 0.8, -0.5, 0.3)
  raw <- posterior::as_draws_array(array(c(b, -b), dim = c(3, 2, 2),
    dimnames = list(NULL, NULL, c("rater[1]", "rater[2]"))))
  fit <- structure(list(fit = list(draws = function(variables) {
    posterior::subset_draws(raw, variable = variables)
  })), class = "bmfrm_fit")
  selected <- posterior::as_draws(fit, variables = "rater[1]")
  direction <- bayestestR::p_direction(selected)
  equiv <- bayestestR::rope(selected, range = c(-0.3, 0.3), ci = 1)
  expect_identical(as.character(direction$Parameter), "rater[1]")
  expect_equal(direction$pd, max(mean(b > 0), mean(b < 0)))
  expect_equal(equiv$ROPE_Percentage, mean(abs(b) <= 0.3))
})
