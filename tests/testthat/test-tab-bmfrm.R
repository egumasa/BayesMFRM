table_test_fit <- function(judges) {
  draws <- cbind(`judge[2]` = c(2, 4, 6, 7), `judge[1]` = c(-2, -4, -6, -7),
                 `theta[1]` = c(1, 2, 3, 4))
  fit <- list(summary = function(variables, ...) {
    mat <- draws[, startsWith(colnames(draws), paste0(variables, "[")), drop = FALSE]
    posterior::summarise_draws(posterior::as_draws_matrix(mat), ...)
  })
  structure(list(fit = fit, spec = list(formula = score ~ writer + judge, family = "rsm",
    person_facet = "writer", facets_main = c("writer", "judge"), facets_bias = character(), K = 4L),
    data_info = list(facet_labels = list(writer = "P", judge = judges),
                     stan_data = list(N = 2L, judge_id = c(1L, 1L), writer_id = c(1L, 1L)))),
    class = "bmfrm_fit")
}

test_that("model tables align original levels and keep unrounded marginal summaries", {
  a <- table_test_fit(c("A", "B"))
  b <- table_test_fit(c("C", "A"))
  b$design <- list(identifiable = FALSE, tolerance = 1e-8, effect_null_space =
    matrix(c(1, 0), 2L, 1L, dimnames = list(c("judge:C", "judge:A"), NULL)))
  b$convergence <- list(status = "incomplete", converged = FALSE)
  x <- tab_bmfrm(First = a, Second = b, ci_level = 0.8, digits = 1)
  expect_equal(unique(x$estimates$facet), "judge")
  expect_equal(x$estimates$mean, rep(c(-4.75, 4.75), 2))
  expect_equal(x$estimates$lower[1L], unname(quantile(c(-2,-4,-6,-7), 0.1)))
  expect_equal(x$models$design, c("Unchecked", "Rank deficient"))
  expect_equal(x$models$mcmc, c("Unchecked", "Incomplete"))
  shown <- bmfrm_display_table(x)
  expect_equal(shown$Level[1:3], c("A", "B", "C"))
  expect_match(shown$Second[1L], "4.8 .*no ratings")
  expect_match(shown$Second[3L], "-4.8 .*unidentified")
  expect_equal(c(shown$First[3L], shown$Second[2L]), c("-", "-"))
  expect_match(capture.output(print(x))[2L], "First.*Second")
  expect_equal(tab_bmfrm(a, facets = "writer")$estimates$role, "ability")
  a$spec$facets_main <- "writer"
  expect_equal(tab_bmfrm(a)$estimates$facet, "writer")
  expect_error(tab_bmfrm(), "bmfrm_fit")
  expect_error(tab_bmfrm(a, facets = "judge"), "facets")
  expect_error(tab_bmfrm(M = a, M = b), "distinct model names")
  expect_error(tab_bmfrm(a, digits = 0.5), "digits")
  expect_error(tab_bmfrm(a, ci_level = 1), "ci_level")
})

test_that("HTML export escapes model and level labels", {
  skip_if_not_installed("knitr")
  path <- tempfile(fileext = ".html")
  a <- table_test_fit(c("<A>", "B"))
  x <- tab_bmfrm(`A & B` = a, file = path)
  html <- paste(readLines(path), collapse = "\n")
  expect_match(html, "&lt;A&gt;", fixed = TRUE)
  expect_match(html, "A &amp; B", fixed = TRUE)
  expect_match(html, "<table>", fixed = TRUE)
  expect_equal(x$estimates$level, c("<A>", "B"))
})
