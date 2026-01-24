# Tests for RSM/PCM family support

test_that("bmfrm() accepts rsm and pcm family values", {
  expect_true("family" %in% names(formals(bmfrm)))
  expect_true("threshold_facet" %in% names(formals(bmfrm)))
})

test_that("parse_bmfrm_formula validates family parameter", {
  test_data <- data.frame(
    score = 1:10,
    person = rep(1:5, 2),
    item = rep(1:2, 5)
  )

  # RSM should work
 spec_rsm <- parse_bmfrm_formula(
    score ~ person + item,
    data = test_data,
    family = "rsm"
  )
  expect_equal(spec_rsm$family, "rsm")
  expect_null(spec_rsm$threshold_facet)

  # PCM requires threshold_facet
  expect_error(
    parse_bmfrm_formula(
      score ~ person + item,
      data = test_data,
      family = "pcm"
    ),
    "threshold_facet"
  )

  # PCM with valid threshold_facet
  spec_pcm <- parse_bmfrm_formula(
    score ~ person + item,
    data = test_data,
    family = "pcm",
    threshold_facet = "item"
  )
  expect_equal(spec_pcm$family, "pcm")
  expect_equal(spec_pcm$threshold_facet, "item")
})

test_that("parse_bmfrm_formula validates threshold_facet exists in main facets", {
  test_data <- data.frame(
    score = 1:10,
    person = rep(1:5, 2),
    item = rep(1:2, 5)
  )

  expect_error(
    parse_bmfrm_formula(
      score ~ person + item,
      data = test_data,
      family = "pcm",
      threshold_facet = "nonexistent"
    ),
    "must be one of the main facets"
  )
})

test_that("get_stan_template_path returns correct paths", {
  # RSM template
  rsm_path <- get_stan_template_path("rsm")
  expect_true(file.exists(rsm_path))
  expect_true(grepl("bmfrm_rsm\\.stan$", rsm_path))

  # PCM template
  pcm_path <- get_stan_template_path("pcm")
  expect_true(file.exists(pcm_path))
  expect_true(grepl("bmfrm_pcm\\.stan$", pcm_path))

  # Invalid family should error
  expect_error(
    get_stan_template_path("invalid"),
    "Unknown family"
  )
})

test_that("RSM template has correct structure", {
  rsm_path <- get_stan_template_path("rsm")
  rsm_code <- paste(readLines(rsm_path), collapse = "\n")

  # Should have vector tau (not ordered, not matrix)
  expect_true(grepl("vector\\[K-1\\] tau_raw", rsm_code))
  # Should NOT have ordered tau
  expect_false(grepl("ordered\\[K-1\\] tau_raw", rsm_code))
  # Should have shared tau indexing
  expect_true(grepl("tau\\[k-1\\]", rsm_code))
})

test_that("PCM template has correct structure", {
  pcm_path <- get_stan_template_path("pcm")
  pcm_code <- paste(readLines(pcm_path), collapse = "\n")

  # Should have matrix tau with threshold facet placeholder
  expect_true(grepl("matrix\\[J_\\{\\{THRESHOLD_FACET\\}\\}, K-1\\] tau_raw", pcm_code))
  # Should have facet-specific tau indexing
  expect_true(grepl("tau\\[threshold_n, k-1\\]", pcm_code))
})

test_that("print.bmfrm_spec shows family information", {
  test_data <- data.frame(
    score = 1:10,
    person = rep(1:5, 2),
    item = rep(1:2, 5)
  )

  spec_rsm <- parse_bmfrm_formula(
    score ~ person + item,
    data = test_data,
    family = "rsm"
  )

  output_rsm <- capture.output(print(spec_rsm))
  expect_true(any(grepl("RSM", output_rsm)))

  spec_pcm <- parse_bmfrm_formula(
    score ~ person + item,
    data = test_data,
    family = "pcm",
    threshold_facet = "item"
  )

  output_pcm <- capture.output(print(spec_pcm))
  expect_true(any(grepl("PCM", output_pcm)))
  expect_true(any(grepl("item", output_pcm)))
})

test_that("build_stan_code generates correct RSM code", {
  test_data <- data.frame(
    score = 1:10,
    person = rep(1:5, 2),
    item = rep(1:2, 5)
  )

  spec <- parse_bmfrm_formula(
    score ~ person + item,
    data = test_data,
    family = "rsm"
  )

  priors <- get_default_priors(spec$facets_main, spec$facets_bias, spec$family)
  stan_code <- build_stan_code(spec, priors)

  # Should have vector tau
  expect_true(grepl("vector\\[K-1\\] tau_raw", stan_code))
  # Should have shared tau indexing (no matrix access)
  expect_true(grepl("tau\\[k-1\\]", stan_code))
  expect_false(grepl("tau\\[.*_n, k-1\\]", stan_code))
})

test_that("build_stan_code generates correct PCM code", {
  test_data <- data.frame(
    score = 1:10,
    person = rep(1:5, 2),
    item = rep(1:2, 5)
  )

  spec <- parse_bmfrm_formula(
    score ~ person + item,
    data = test_data,
    family = "pcm",
    threshold_facet = "item"
  )

  priors <- get_default_priors(spec$facets_main, spec$facets_bias, spec$family)
  stan_code <- build_stan_code(spec, priors)

  # Should have matrix tau with item dimension
  expect_true(grepl("matrix\\[J_item, K-1\\] tau_raw", stan_code))
  # Should have facet-specific tau indexing
  expect_true(grepl("tau\\[threshold_n, k-1\\]", stan_code))
  # Should have threshold_n extraction
  expect_true(grepl("int threshold_n = item_id\\[n\\]", stan_code))
})

# Integration tests - only run if cmdstan is available
test_that("bmfrm() works with RSM family (integration)", {
  skip_if_not(check_cmdstan(), "CmdStan not available")
  skip_on_cran()

  set.seed(123)
  test_data <- data.frame(
    score = sample(1:6, 100, replace = TRUE),
    person = rep(1:10, each = 10),
    item = rep(1:5, 20),
    rater = sample(1:3, 100, replace = TRUE)
  )

  fit <- suppressWarnings(
    bmfrm(
      score ~ person + item + rater,
      data = test_data,
      K = 6,
      family = "rsm",
      chains = 1,
      iter = 100,
      warmup = 50
    )
  )

  expect_s3_class(fit, "bmfrm_fit")
  expect_equal(fit$spec$family, "rsm")
})

test_that("bmfrm() works with PCM family (integration)", {
  skip_if_not(check_cmdstan(), "CmdStan not available")
  skip_on_cran()

  set.seed(123)
  test_data <- data.frame(
    score = sample(1:6, 100, replace = TRUE),
    person = rep(1:10, each = 10),
    item = rep(1:5, 20),
    rater = sample(1:3, 100, replace = TRUE)
  )

  fit <- suppressWarnings(
    bmfrm(
      score ~ person + item + rater,
      data = test_data,
      K = 6,
      family = "pcm",
      threshold_facet = "item",
      chains = 1,
      iter = 100,
      warmup = 50
    )
  )

  expect_s3_class(fit, "bmfrm_fit")
  expect_equal(fit$spec$family, "pcm")
  expect_equal(fit$spec$threshold_facet, "item")
})