test_that("streaming LOO preserves a facet named y", {
  d <- expand.grid(person = c("P1", "P2"), y = c("Y1", "Y2"),
                   replicate = 1:3, stringsAsFactors = FALSE)
  d$score <- ifelse(d$y == "Y1", 2L, 1L)
  D <- 80L
  theta <- cbind(seq(-0.2, 0.2, length.out = D), seq(0.3, 0.7, length.out = D))
  y <- matrix(rep(c(-1, 1), each = D), D, 2L)
  colnames(theta) <- paste0("theta[", 1:2, "]")
  colnames(y) <- paste0("y[", 1:2, "]")
  # Independent binary likelihood oracle, using the original labels and scores.
  eta <- theta[, match(d$person, c("P1", "P2"))] - y[, match(d$y, c("Y1", "Y2"))]
  ll <- log(plogis(sweep(eta, 2L, 2 * d$score - 3, "*")))
  colnames(ll) <- paste0("log_lik[", seq_len(nrow(d)), "]")
  for (family in c("rsm", "pcm")) {
    spec <- resolve_bmfrm_spec(parse_bmfrm_formula(score ~ person + y, d, K = 2,
      family = family, threshold_facet = if (family == "pcm") "y" else NULL))
    di <- prepare_data_bmfrm(spec, d)
    tau <- if (family == "rsm") matrix(0, D, 1L, dimnames = list(NULL, "tau[1]")) else
      matrix(0, D, 2L, dimnames = list(NULL, c("tau[1,1]", "tau[2,1]")))
    store <- list(theta = theta, y = y, tau = tau, log_lik = ll)
    fit <- list(
      metadata = function() list(stan_variables = names(store), iter_sampling = 40L),
      num_chains = function() 2L,
      draws = function(variables, format = "matrix") {
        out <- store[[variables]]
        if (format == "draws_array") return(array(out, c(40L, 2L, ncol(out)),
                                                  dimnames = list(NULL, NULL, colnames(out))))
        out
      }
    )
    obj <- structure(list(fit = fit, spec = spec, data_info = di,
                          convergence = list(converged = TRUE)), class = "bmfrm_fit")
    expected <- suppressWarnings(loo::loo(ll, r_eff = NA, cores = 1))
    saved <- suppressWarnings(loo::loo(obj, r_eff = NA, cores = 1, use_stan_log_lik = TRUE))
    streamed <- suppressWarnings(loo::loo(obj, r_eff = NA, cores = 1, use_stan_log_lik = FALSE))
    expect_equal(unname(streamed$pointwise), unname(expected$pointwise), tolerance = 1e-12)
    expect_equal(streamed$estimates, saved$estimates, tolerance = 1e-12)
  }
})

test_that("response expressions are rejected consistently before data processing", {
  d <- data.frame(person = rep(1:5, 2), score = rep(1:2, each = 5))
  for (f in list(I(score + 1) ~ person, log(score) ~ person, cbind(score, score) ~ person, ~ person)) {
    expect_error(parse_bmfrm_formula(f, d, K = 3), "single column name")
    expect_error(check_design_bmfrm(f, d, K = 3), "single column name")
    expect_error(simulate_bmfrm(f, d, K = 3), "single column name")
  }
  d[["shifted score"]] <- d$score + 1L
  spec <- parse_bmfrm_formula(`shifted score` ~ person, d, K = 3)
  expect_identical(prepare_data_bmfrm(spec, d)$stan_data$X, d[["shifted score"]])
})

test_that("observed PCM designs retain the rubric for all planned levels", {
  d <- expand.grid(person = letters[1:3], item = c("a", "b"), stringsAsFactors = FALSE)
  d$score <- c(1, 2, 1, NA, NA, NA)
  for (K in list(4L, c(a = 4L, b = 2L))) {
    for (method in c("dense", "schur")) {
      char <- check_design_bmfrm(score ~ person + item, d, K = K, family = "pcm",
                                threshold_facet = "item", observed_only = TRUE, method = method)
      f <- d
      f$item <- factor(f$item)
      fact <- check_design_bmfrm(score ~ person + item, f, K = K, family = "pcm",
                                threshold_facet = "item", observed_only = TRUE, method = method)
      expect_equal(char, fact)
      expect_false(char$identifiable)
      expect_equal(char$n_observations, 3L)
      expect_equal(char$n_dropped, 3L)
      expect_equal(char$adjacent_logit_rank, 5L)
      expect_equal(char$adjacent_logit_parameters, if (length(K) == 1L) 8L else 6L)
    }
  }
})

test_that("Stan keywords and generated identifiers are rejected before compilation", {
  for (name in c("lower", "upper", "offset", "multiplier", "threshold", "n", "k", "j",
                 "K_max", "K_level", "N_tau_steps", "tau_start", "tau_step", "raw_j",
                 "simplex", "row_vector", "sum_to_zero_vector", "profile", "foo__",
                 "class", "template", "fixed_sd_item", "scale_factor_itemrater", "mean_item", "centered")) {
    d <- data.frame(person = 1:10, score = rep(1:2, 5))
    d[[name]] <- rep(c("a", "b"), 5)
    expect_error(parse_bmfrm_formula(stats::reformulate(c("person", name), "score"),
                                    d, K = 3), "nonreserved Stan identifiers")
  }
  expect_silent(validate_facet_names(c("person", "item", "rater", "y", "task2")))
})

test_that("syntax validation invokes Stan and returns compiler diagnostics", {
  skip_if_not(check_cmdstan())
  expect_true(validate_stan_code("parameters { real z; } model { z ~ normal(0, 1); }")$is_valid)
  bad <- suppressMessages(validate_stan_code("this is not valid Stan"))
  expect_false(bad$is_valid)
  expect_match(bad$error, "[Ss]yntax")
  missing <- suppressMessages(validate_stan_code("model { target += not_declared; }"))
  expect_false(missing$is_valid)
  expect_match(missing$error, "not_declared")
  d <- expand.grid(person = 1:3, y = c("a", "b"))
  d$score <- rep(1:2, 3)
  for (K in list(3L, c(a = 2L, b = 3L))) {
    spec <- parse_bmfrm_formula(score ~ person + y, d, K = K, family = "pcm", threshold_facet = "y")
    checked <- validate_stan_code(build_stan_code(spec))
    expect_true(checked$is_valid, info = checked$error)
  }
})

test_that("residual pairing is independent of item label spelling", {
  d <- expand.grid(person = 1:4, item = c("item1", "item2"), stringsAsFactors = FALSE)
  d$resid <- c(-2, -1, 1, 2, 1, -2, -1, 2)
  a <- residual_cor_criteria(d)
  d$item <- sub("item", "label", d$item)
  b <- residual_cor_criteria(d)
  expect_equal(unname(a$cor_matrix), unname(b$cor_matrix))
  expect_equal(a$cor_tidy$correlation, rep(0.3, 2))
  expect_equal(a$cor_tidy$n_pairs, rep(4, 2))
  expect_setequal(a$cor_tidy$item1, c("item1", "item2"))
})

test_that("singleton strata do not overwrite other CV fold assignments", {
  d <- data.frame(person = c("a", "a", "b", "b", "c"), item = "i")
  for (seed in 1:10) {
    folds <- build_cv_folds(d, c("person", "item"), "rating", NULL, 2L, seed = seed)
    expect_setequal(folds[1:2], 1:2)
    expect_setequal(folds[3:4], 1:2)
    expect_true(folds[5] %in% 1:2)
    expect_equal(sort(as.integer(table(folds))), c(2L, 3L))
  }
})

test_that("ignored fitting options fail before compilation and custom priors remain unknown", {
  calls <- 0L
  testthat::local_mocked_bindings(
    check_cmdstan = function() TRUE,
    compile_bmfrm = function(stan_code, ...) {
      calls <<- calls + 1L
      list(model = NULL, cache_file = "unused", hash = digest::digest(stan_code))
    },
    fit_cmdstan_cached = function(...) list(),
    check_convergence = function(...) list(converged = TRUE),
    .package = "BayesMFRM"
  )
  d <- data.frame(person = rep(1:5, 2), rater = rep(1:2, each = 5), score = rep(1:2, each = 5))
  file <- bmfrm_model("two-facet/person-rater-main")$path
  expect_error(bmfrm(score ~ person + rater, d, K = 2, stan_file = file,
                     priors = prior("normal(0, 0.01)", "theta")), "priors cannot be supplied")
  expect_error(bmfrm(score ~ person + rater, d, K = 2, save_csvs = FALSE), "save_csvs must be TRUE")
  expect_equal(calls, 0L)
  x <- suppressMessages(bmfrm(score ~ person + rater, d, K = 2, stan_file = file))
  expect_null(x$priors)
  expect_null(x$spec$priors)
  expect_null(x$data_info$spec$priors)
  expect_null(x$fit_args$priors)
  expect_match(x$stan_code, "theta_raw +~ student_t\\(3, 0, 4\\)")
  expect_equal(calls, 1L)
})

test_that("cache cleanup removes old Unix binaries and preserves unrelated files", {
  cache <- tempfile("audit-cache-")
  dir.create(cache)
  on.exit(unlink(cache, recursive = TRUE))
  old <- file.path(cache, c("model.stan", "model", "model.exe", "model.rds"))
  keep <- file.path(cache, c("recent.stan", "recent", "notes", "draws.csv"))
  invisible(file.create(c(old, keep)))
  Sys.setFileTime(old, Sys.time() - 60 * 86400)
  Sys.setFileTime(keep[3:4], Sys.time() - 60 * 86400)
  expect_setequal(suppressMessages(clean_model_cache(cache, dry_run = TRUE)), old)
  expect_true(all(file.exists(c(old, keep))))
  expect_setequal(suppressMessages(clean_model_cache(cache)), old)
  expect_false(any(file.exists(old)))
  expect_true(all(file.exists(keep)))
})
