agreement_fit <- function(d, D = 200L) {
  spec <- parse_bmfrm_formula(score ~ writer + judge, d, K = 2, person_facet = "writer")
  expect_warning(info <- prepare_data_bmfrm(spec, d), "Very few observations")
  # Shared ability draws make the mean of products differ from the product of means.
  eta <- rep(c(-log(9), log(9)), length.out = D)
  mat <- list(theta = matrix(eta, D, length(info$facet_labels$writer)),
              judge = matrix(0, D, length(info$facet_labels$judge)), tau = matrix(0, D, 1L))
  fit <- list(metadata = function() list(stan_variables = names(mat)),
              draws = function(variables, format) mat[[variables]])
  structure(list(spec = spec, fit = fit, data_info = info), class = "bmfrm_fit")
}

test_that("exact agreement uses matched probabilities and shares replicated scores across pairs", {
  d <- data.frame(writer = "one", judge = c("B", "C", "A"), score = c(2, 1, 2))
  fit <- agreement_fit(d)
  set.seed(987)
  before <- .Random.seed
  x <- rater_agreement(fit, "judge")
  expect_identical(before, .Random.seed)
  expect_equal(x, rater_agreement(fit, "judge"))
  expect_equal(x$observed, c(1, 0, 0, 1/3))
  expect_equal(x$n_opportunities, c(1, 1, 1, 3))
  expect_equal(x$expected, rep(.82, 4))
  # With three binary scores only one or three pairs can agree in a replicate.
  expect_equal(x$rep_lower[4], 1/3)
  expect_equal(x$rep_upper[4], 1)
  expect_equal(x$expected, rater_agreement(fit, "judge", seed = 7)$expected)
  expect_s3_class(plot_rater_agreement(x), "ggplot")
  expect_error(rater_agreement(fit, "writer"), "non-ability")
  expect_error(rater_agreement(fit, "judge", ci_level = 0), "ci_level")
  fit$spec$custom_stan <- TRUE
  expect_error(rater_agreement(fit, "judge"), "verified")
})

test_that("repeated performances, missing ratings and unmatched raters keep correct denominators", {
  d <- data.frame(writer = "one", judge = c("A", "B", "A", "B", "C", "D"),
                  session = c(1, 1, 2, 2, 2, 3), score = c(1, 1, 2, NA, 1, 2))
  fit <- suppressMessages(agreement_fit(d))
  expect_error(rater_agreement(fit, "judge"), "occasion")
  expect_error(rater_agreement(fit, "judge", occasion = "score"), "identifier")
  x <- rater_agreement(fit, "judge", occasion = "session")
  expect_equal(x$rater1, c("A", "A", NA))
  expect_equal(x$rater2, c("B", "C", NA))
  expect_equal(x$n_opportunities, c(1, 1, 2))
  expect_equal(x$n_exact, c(1, 0, 1))
  fit$data_info$data_clean$session[1] <- NA
  expect_error(rater_agreement(fit, "judge", occasion = "session"), "nonmissing")

  d$session <- seq_len(nrow(d))
  none <- rater_agreement(suppressMessages(agreement_fit(d)), "judge", occasion = "session")
  expect_equal(none$scope, "overall")
  expect_equal(none$n_opportunities, 0L)
  expect_equal(none$n_exact, 0L)
  expect_true(all(is.na(none[c("observed", "expected", "rep_lower", "rep_upper")])))
  expect_error(plot_rater_agreement(none), "No jointly rated")
})
