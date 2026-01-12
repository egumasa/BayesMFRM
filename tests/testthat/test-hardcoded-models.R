# Tests for hardcoded Stan model support

test_that("bmfrm_models() lists available models", {
  models <- bmfrm_models()


  expect_type(models, "character")
  expect_true(length(models) >= 2)
  expect_true(any(grepl("four-facet", models)))
  expect_true(any(grepl("three-facet", models)))
})

test_that("bmfrm_model() returns correct structure", {
  model <- bmfrm_model("four-facet/person-item-rater-scenario-main")

  expect_s3_class(model, "bmfrm_hardcoded_model")
  expect_true(file.exists(model$path))
  expect_type(model$code, "character")
  expect_true(nchar(model$code) > 100)
  expect_null(model$stanmodel)  # compile = FALSE by default
})

test_that("bmfrm_model() fails gracefully for non-existent model", {
  expect_error(
    bmfrm_model("non-existent/model"),
    "not found"
  )
})

test_that("hardcoded Stan files have correct variable names", {
  # Four-facet model
  model_4f <- bmfrm_model("four-facet/person-item-rater-scenario-main")

  # Check for prepare_data_bmfrm() style variable names

  expect_true(grepl("J_person", model_4f$code))
  expect_true(grepl("J_item", model_4f$code))
  expect_true(grepl("J_rater", model_4f$code))
  expect_true(grepl("J_scenario", model_4f$code))
  expect_true(grepl("person_id", model_4f$code))
  expect_true(grepl("item_id", model_4f$code))
  expect_true(grepl("rater_id", model_4f$code))
  expect_true(grepl("scenario_id", model_4f$code))

  # Should NOT have old-style names
 expect_false(grepl("ExamineeID", model_4f$code))
  expect_false(grepl("ItemID", model_4f$code))
  expect_false(grepl("RaterID", model_4f$code))
  expect_false(grepl("ScenarioID", model_4f$code))

  # Three-facet model
  model_3f <- bmfrm_model("three-facet/person-rater-scenario-main")

  expect_true(grepl("J_person", model_3f$code))
  expect_true(grepl("J_rater", model_3f$code))
  expect_true(grepl("J_scenario", model_3f$code))
  expect_true(grepl("person_id", model_3f$code))
  expect_true(grepl("rater_id", model_3f$code))
  expect_true(grepl("scenario_id", model_3f$code))
})

test_that("print method works for bmfrm_hardcoded_model", {
  model <- bmfrm_model("four-facet/person-item-rater-scenario-main")

  output <- capture.output(print(model))
  expect_true(any(grepl("BayesMFRM Hardcoded Model", output)))
  expect_true(any(grepl("Name:", output)))
  expect_true(any(grepl("Path:", output)))
})

test_that("bmfrm() accepts stan_file argument", {
  # Just check the function signature accepts the argument
  # Full integration test requires Stan compilation
  expect_true("stan_file" %in% names(formals(bmfrm)))
})

# Integration test - only run if cmdstan is available
test_that("bmfrm() works with hardcoded model (integration)", {
  skip_if_not(check_cmdstan(), "CmdStan not available")
  skip_on_cran()

  # Create minimal test data
  set.seed(123)
  test_data <- data.frame(
    score = sample(1:6, 100, replace = TRUE),
    person = rep(1:10, each = 10),
    item = rep(1:5, 20),
    rater = sample(1:3, 100, replace = TRUE),
    scenario = sample(1:2, 100, replace = TRUE)
  )

  # Test with hardcoded model
  # Suppress ESS warnings (expected with minimal iterations)
  fit <- suppressWarnings(
    bmfrm(
      score ~ person + item + rater + scenario,
      data = test_data,
      K = 6,
      stan_file = bmfrm_model("four-facet/person-item-rater-scenario-main")$path,
      chains = 1,
      iter = 100,
      warmup = 50
    )
  )

  expect_s3_class(fit, "bmfrm_fit")
  expect_true(!is.null(fit$fit))
})
