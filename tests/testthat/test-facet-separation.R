test_that("EAP and conventional separation use distinct variance conventions", {
  tab <- data.frame(level = letters[1:3], mean = -1:1, sd = c(.4, .5, .6))
  u <- mean(tab$sd^2)
  eap <- facet_separation(tab, "rater")
  expect_equal(eap$reliability, 1 / (1 + u))
  expect_equal(eap$separation, sqrt(1 / u))
  expect_equal(eap$strata, (4 * sqrt(1 / u) + 1) / 3)
  names(tab)[2:3] <- c("measure", "se")
  rasch <- facet_separation(tab, "rater", method = "rasch")
  expect_equal(rasch$reliability, 1 - u)
  expect_equal(rasch$separation, sqrt((1 - u) / u))
  expect_equal(rasch$rms_uncertainty, sqrt(u))
  expect_identical(rasch$method, "rasch")
  tab$se <- 2
  zero <- facet_separation(tab, "rater", method = "rasch")
  expect_equal(c(zero$separation, zero$strata, zero$reliability), c(0, 1/3, 0))
})

test_that("supplied summaries validate eligible levels and zero uncertainty", {
  tab <- data.frame(level = c("B", "A", "unused", "A"), facet = c(rep("rater", 3), "item"),
                    mean = c(1, 0, NA, 2), sd = c(0, 0, NA, 1), n_obs = c(1, 4, 0, 2))
  x <- facet_separation(tab, "rater")
  expect_equal(c(x$n_levels, x$n_excluded, x$separation, x$strata, x$reliability),
               c(2, 1, Inf, Inf, 1))
  tab$mean[1] <- 0
  x <- facet_separation(tab, "rater")
  expect_true(all(is.na(c(x$separation, x$strata, x$reliability))))
  expect_error(facet_separation(tab, "item"), "two observed")
  expect_error(facet_separation(tab, "rater", "rasch"), "measure, se")
  tab$sd[2] <- -1
  expect_error(facet_separation(tab, "rater"), "nonnegative")
  tab$n_obs[2] <- NA
  expect_error(facet_separation(tab, "rater"), "integer counts")
  expect_error(facet_separation(rbind(tab, tab[1, ]), "rater"), "unique")
})

test_that("the EAP reliability matches mirt on the bundled writing example", {
  skip_if_not_installed("mirt")
  data("bmfrm_example", package = "BayesMFRM", envir = environment())
  tab <- bmfrm_example$rater
  expected <- mirt::empirical_rxx(cbind(F1 = tab$mean, SE_F1 = tab$sd))
  expect_equal(facet_separation(tab, "rater")$reliability, unname(expected))
})

test_that("fitted summaries exclude unused levels and preserve joint draws and design information", {
  dat <- data.frame(writer = c("w1", "w1", "w2"), rater = c("A", "B", "C"), score = c(1, 2, 2))
  dc <- check_design_bmfrm(score ~ writer + rater, dat, K = 2, person_facet = "writer")
  mat <- cbind(2:6, 1:5, 4:8, 11:15)
  colnames(mat) <- sprintf("rater[%d]", c(2, 1, 3, 4))
  object <- structure(list(
    spec = list(facets_main = c("writer", "rater"), person_facet = "writer"),
    fit = list(draws = function(variables, format) mat[, colnames(mat) %in% variables, drop = FALSE],
               summary = function(variables, ...) posterior::summarise_draws(posterior::as_draws_matrix(mat), ...)),
    data_info = list(facet_labels = list(rater = c("A", "B", "C", "unused")),
                     stan_data = list(rater_id = 1:3)), design = dc), class = "bmfrm_fit")

  x <- facet_separation(object, "rater")
  expect_equal(c(x$n_levels, x$n_excluded, x$estimate_sd), c(3, 1, sd(c(3, 4, 6))))
  expect_error(facet_separation(object, "rater", method = "rasch"), "conventional Rasch")
  y <- facet_contrasts(object, "rater", c("C", "B", "A"), ci_level = .8)
  expect_equal(y$comparison, c("B", "A", "A"))
  expect_equal(y$reference, c("C", "C", "B"))
  expect_equal(y$difference, c(-2, -3, -1))
  expect_equal(y$upper, y$lower)
  expect_equal(y$lower, y$difference)
  expect_equal(y$prob_positive, c(0, 0, 0))
  expect_equal(y$likelihood_identified, c(FALSE, FALSE, TRUE))
  expect_error(facet_contrasts(object, "rater", c("A", "unused")), "observed")
  expect_error(facet_contrasts(object, "rater", c("A", "Z")), "Unknown level")
  expect_error(facet_contrasts(object, "rater", c("A", "A")), "distinct")
  expect_error(facet_contrasts(object, "rater", c("A", "B"), ci_level = 1), "ci_level")
  object$design <- NULL
  expect_true(is.na(facet_contrasts(object, "rater", c("A", "B"))$likelihood_identified))
})
