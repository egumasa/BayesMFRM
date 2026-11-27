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
