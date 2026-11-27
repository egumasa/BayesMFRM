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

test_that("formula parsing and design check handle K=NULL errors, offset(), dot formulas, and reserved names", {
  d_bad_chr <- data.frame(person = c("P1", "P2"), rater = c("R1", "R2"), score = c("a", "b"))
  expect_error(parse_bmfrm_formula(score ~ person + rater, d_bad_chr, K = NULL), "finite integer scores")

  d_bad_inf <- data.frame(person = c("P1", "P2"), rater = c("R1", "R2"), score = c(1, Inf))
  expect_error(parse_bmfrm_formula(score ~ person + rater, d_bad_inf, K = NULL), "finite integer scores")

  d_bad_na <- data.frame(person = c("P1", "P2"), rater = c("R1", "R2"), score = c(NA_real_, NA_real_))
  expect_error(parse_bmfrm_formula(score ~ person + rater, d_bad_na, K = NULL), "No valid observations")

  d_ok <- expand.grid(person = c("P1", "P2", "P3"), rater = c("R1", "R2"), stringsAsFactors = FALSE)
  d_ok$score <- c(1L, 2L, 3L, 2L, 3L, 1L)
  spec_dot <- parse_bmfrm_formula(score ~ ., d_ok, K = 3L)
  expect_identical(spec_dot$facets_main, c("person", "rater"))
  dc_dot <- check_design_bmfrm(score ~ ., d_ok, K = 3L)
  expect_true(dc_dot$connected)

  d_off <- d_ok
  d_off$w <- 1
  expect_error(parse_bmfrm_formula(score ~ person + rater + offset(w), d_off, K = 3L), "Formula offsets are not supported")
  expect_error(check_design_bmfrm(score ~ person + rater + offset(w), d_off, K = 3L), "Formula offsets are not supported")

  for (bad_nm in c("theta", "tau", "model", "item_raw")) {
    d_res <- data.frame(person = c("P1", "P2"), score = c(1L, 2L))
    d_res[[bad_nm]] <- c("A", "B")
    fml_bad <- stats::reformulate(c("person", bad_nm), "score")
    expect_error(check_design_bmfrm(fml_bad, d_res, K = 2L), "nonreserved Stan identifiers")
  }
})

test_that("check_data_suitability validates columns, ignores NA in level counts, and aligns unused factor levels", {
  d <- data.frame(
    person = c("P1", "P1", "P2", "P2", NA_character_),
    rater = factor(c("R1", "R2", "R1", "R2", "R1"), levels = c("R1", "R2", "R_UNUSED")),
    score = c(1L, 2L, 2L, 3L, 2L),
    stringsAsFactors = FALSE
  )
  expect_error(check_data_suitability(d, "missing_col", c("person", "rater")), "not found in data")
  expect_error(check_data_suitability(d, "score", c("person", "nope")), "not found in data")

  suit <- check_data_suitability(d, "score", c("person", "rater"))
  expect_equal(suit$facet_info$person$n_levels, 2L)
  expect_equal(suit$facet_info$rater$n_levels, 3L)
  expect_equal(suit$facet_info$rater$min_obs_per_level, 0L)
  expect_gt(length(suit$warnings), 0L)
})

test_that("fit caching prevents signature collisions, recovers from corrupt RDS, and cleans CSV subdirectories safely", {
  cache <- tempfile("fit-cache-audit-")
  dir.create(cache, recursive = TRUE)
  on.exit(unlink(cache, recursive = TRUE))

  base_rds <- file.path(cache, "bmfrm_rsm_abcd.rds")
  p1 <- BayesMFRM:::resolve_fit_cache_file(base_rds, "code1", list(N = 4L), seed = 1L, chains = 2L, parallel_chains = 2L)
  p2 <- BayesMFRM:::resolve_fit_cache_file(base_rds, "code1", list(N = 5L), seed = 1L, chains = 2L, parallel_chains = 2L)
  expect_false(identical(p1, p2))

  csv_dir <- file.path(cache, "csv")
  dir.create(csv_dir, recursive = TRUE)
  active_csv <- file.path(csv_dir, "active-chain1.csv")
  orphan_csv <- file.path(csv_dir, "orphan-chain1.csv")
  top_csv <- file.path(cache, "draws.csv")
  cv_dir <- file.path(cache, "cv_folds")
  cv_csv_dir <- file.path(cv_dir, "csv")
  dir.create(cv_csv_dir, recursive = TRUE)
  old_fold_rds <- file.path(cv_dir, "fold1.rds")
  old_fold_csv <- file.path(cv_csv_dir, "fold1-chain1.csv")
  invisible(file.create(c(active_csv, orphan_csv, top_csv, old_fold_rds, old_fold_csv)))

  # Mock cmdstanr model for corrupt RDS recovery & active CSV tracking
  sampled <- 0L
  mock_fit <- list(output_files = function() active_csv)
  mock_model <- list(
    code = function() "parameters { real z; }",
    sample = function(data, output_dir, ...) {
      sampled <<- sampled + 1L
      mock_fit
    }
  )
  writeLines("not a valid RDS binary", p1)
  res_fit <- suppressMessages(fit_cmdstan_cached(mock_model, list(N = 4L), file = p1, refit = "on_change", seed = 1L, chains = 2L))
  expect_equal(sampled, 1L)
  expect_identical(res_fit$output_files(), active_csv)

  # Legacy base_rds match check
  file.copy(p1, base_rds, overwrite = TRUE)
  unlink(p1)
  expect_identical(
    BayesMFRM:::resolve_fit_cache_file(base_rds, "parameters { real z; }", list(N = 4L), seed = 1L, chains = 2L, parallel_chains = 2L),
    base_rds
  )

  # Age files: base_rds is recent (references active_csv); orphan_csv, top_csv, old_fold_rds, old_fold_csv are 60 days old
  Sys.setFileTime(c(active_csv, orphan_csv, top_csv, old_fold_rds, old_fold_csv), Sys.time() - 60 * 86400)
  Sys.setFileTime(base_rds, Sys.time())
  # Default include_csv = FALSE: CSVs are untouched so externally saved/in-memory fits remain usable
  removed_def <- suppressMessages(clean_model_cache(cache, max_age_days = 30))
  expect_setequal(removed_def, old_fold_rds)
  expect_true(all(file.exists(c(base_rds, active_csv, orphan_csv, old_fold_csv, top_csv))))

  # Explicit include_csv = TRUE: removes unreferenced expired CSVs while preserving active_csv and top_csv
  removed_csv <- suppressMessages(clean_model_cache(cache, max_age_days = 30, include_csv = TRUE))
  expect_setequal(removed_csv, c(orphan_csv, old_fold_csv))
  expect_true(all(file.exists(c(base_rds, active_csv, top_csv))))

  # refit = "never" MUST stop before calling model$sample() when target cache file does not exist
  sampled_before <- sampled
  missing_rds <- file.path(cache, "nonexistent_fit.rds")
  expect_error(
    fit_cmdstan_cached(mock_model, list(N = 4L), file = missing_rds, refit = "never", seed = 1L, chains = 2L),
    "cannot be verified"
  )
  expect_equal(sampled, sampled_before)
})

test_that("bmfrm pre-validates sampler arguments and stan_file + minimal GQ before compilation", {
  calls <- 0L
  testthat::local_mocked_bindings(
    check_cmdstan = function() TRUE,
    compile_bmfrm = function(stan_code, ...) {
      calls <<- calls + 1L
      list(model = NULL, cache_file = "unused", hash = digest::digest(stan_code))
    },
    .package = "BayesMFRM"
  )
  d <- data.frame(person = rep(1:4, 2), rater = rep(1:2, each = 4), score = rep(1:2, each = 4))
  f_stan <- bmfrm_model("two-facet/person-rater-main")$path

  expect_error(bmfrm(score ~ person + rater, d, K = 2, iter = 500, warmup = 500), "warmup must be a non-negative integer strictly less than iter")
  expect_error(bmfrm(score ~ person + rater, d, K = 2, chains = 0), "chains must be a positive integer")
  expect_error(bmfrm(score ~ person + rater, d, K = 2, iter_sampling = 200), "reserved cmdstanr argument")
  expect_error(bmfrm(score ~ person + rater, d, K = 2, stan_file = f_stan, save_gq = "minimal"), "save_gq = 'minimal' requires generated Stan code")
  expect_equal(calls, 0L)
})

test_that("bundled 4-facet beta mapping, stancode/standata generics, facet_summary index sorting, and fair_scores reference_meta work", {
  d4 <- expand.grid(
    person = c("P1", "P2", "P3"),
    item = paste0("Item_", 1:10),
    rater = c("R1", "R2"),
    scenario = c("S1", "S2"),
    stringsAsFactors = FALSE
  )
  d4$score <- rep(1:3, length.out = nrow(d4))
  f4_path <- bmfrm_model("four-facet/person-item-rater-scenario-main")$path
  f4_code <- paste(readLines(f4_path), collapse = "\n")
  spec4 <- parse_bmfrm_formula(score ~ person + item + rater + scenario, d4, K = 3L)
  spec4$custom_stan <- TRUE
  spec4$eval_contract <- BayesMFRM:::resolve_eval_contract(spec4, f4_code)
  di4 <- prepare_data_bmfrm(spec4, d4)

  D <- 16L
  mk_mat <- function(prefix, vals) {
    m <- matrix(rep(vals, each = D), nrow = D, ncol = length(vals))
    colnames(m) <- paste0(prefix, "[", seq_along(vals), "]")
    m
  }
  beta_true <- seq(-0.9, 0.9, length.out = 10L)
  store4 <- list(
    theta = mk_mat("theta", c(-0.5, 0.0, 0.5)),
    beta = mk_mat("beta", beta_true),
    rater = mk_mat("rater", c(-0.2, 0.2)),
    scenario = mk_mat("scenario", c(-0.1, 0.1)),
    tau = mk_mat("tau", c(-0.4, 0.4))
  )
  # Deliberately return lexicographically ordered [1], [10], [2], ... in summary()
  lex_ord <- order(paste0("beta[", 1:10, "]"))
  fit_mock4 <- list(
    num_chains = function() 2L,
    metadata = function() list(stan_variables = names(store4), iter_sampling = 8L),
    draws = function(var = NULL, format = "matrix") {
      if (is.null(var)) return(do.call(cbind, store4))
      if (!var %in% names(store4)) stop("Missing ", var)
      m <- store4[[var]]
      if (format == "draws_array") return(array(m, c(8L, 2L, ncol(m)), dimnames = list(NULL, NULL, colnames(m))))
      m
    },
    summary = function(variables = NULL, ...) {
      m <- store4[[variables]]
      if (identical(variables, "beta")) m <- m[, lex_ord, drop = FALSE]
      s_df <- tibble::tibble(
        variable = colnames(m),
        mean = colMeans(m),
        sd = apply(m, 2L, stats::sd)
      )
      dots <- list(...)
      fml_idx <- which(vapply(dots, inherits, logical(1), "formula"))
      if (length(fml_idx)) {
        q_fn <- rlang::as_function(dots[[fml_idx[1L]]])
        q_mat <- t(apply(m, 2L, q_fn))
        if (nrow(q_mat) == 1L && ncol(m) > 1L) q_mat <- t(q_mat)
        q_df <- as.data.frame(q_mat)
        s_df <- dplyr::bind_cols(s_df, q_df)
      }
      s_df
    }
  )
  obj4 <- structure(
    list(fit = fit_mock4, spec = spec4, data_info = di4, stan_code = f4_code,
         convergence = list(converged = TRUE, status = "ok")),
    class = "bmfrm_fit"
  )

  invisible(utils::capture.output(sc_ret <- stancode(obj4)))
  expect_identical(sc_ret, f4_code)
  expect_identical(standata(obj4), di4$stan_data)

  fs_item <- facet_summary(obj4, "item", probs = c(0.1, 0.9))
  expect_equal(nrow(fs_item), 10L)
  expect_identical(fs_item$level, di4$facet_labels$item)
  expect_equal(fs_item$mean, beta_true, tolerance = 1e-12)

  expect_no_warning(sum_txt <- utils::capture.output(summary(obj4, probs = 0.5)))
  expect_match(paste(sum_txt, collapse = "\n"), "Item +Parameters:")

  # Verify R-streaming & standardized fair_scores() + saved Stan mu_fair fast-path attach reference_meta unconditionally
  d2 <- expand.grid(person = paste0("Person_", 1:5), rater = c("Rater_A", "Rater_B"), stringsAsFactors = FALSE)
  d2$score <- rep(1:3, length.out = nrow(d2))
  spec2 <- resolve_bmfrm_spec(parse_bmfrm_formula(score ~ person + rater, d2, K = 3L))
  di2 <- prepare_data_bmfrm(spec2, d2)
  dc2 <- check_design_bmfrm(score ~ person + rater, d2, K = 3L)

  store2_min <- list(
    theta = mk_mat("theta", seq(-0.4, 0.4, length.out = 5L)),
    rater = mk_mat("rater", c(-0.2, 0.2)),
    tau = mk_mat("tau", c(-0.5, 0.5))
  )
  fit_mock2_min <- list(
    num_chains = function() 2L,
    metadata = function() list(stan_variables = names(store2_min), iter_sampling = 8L),
    draws = function(var = NULL, format = "matrix") store2_min[[var]]
  )
  obj2_min <- structure(list(fit = fit_mock2_min, spec = spec2, data_info = di2,
                             design = dc2, convergence = list(converged = TRUE)), class = "bmfrm_fit")

  for (std_flag in c(FALSE, TRUE)) {
    fs_sum <- fair_scores(obj2_min, summary = TRUE, standardize = std_flag, return_draws = FALSE)
    fs_drw <- fair_scores(obj2_min, summary = TRUE, standardize = std_flag, return_draws = TRUE)
    fs_row <- fair_scores(obj2_min, summary = FALSE, standardize = std_flag)
    expect_type(attr(fs_sum, "reference_meta"), "list")
    expect_type(attr(fs_drw, "draws"), "list")
    expect_type(attr(fs_drw, "reference_meta"), "list")
    expect_type(attr(fs_row, "reference_meta"), "list")
  }

  store2 <- c(store2_min, list(mu_fair = mk_mat("mu_fair", rep(seq(1.6, 2.4, length.out = 5L), times = 2L))))
  fit_mock2 <- list(
    num_chains = function() 2L,
    metadata = function() list(stan_variables = names(store2), iter_sampling = 8L),
    draws = function(var = NULL, format = "matrix") store2[[var]]
  )
  obj2 <- structure(list(fit = fit_mock2, spec = spec2, data_info = di2,
                         design = dc2, convergence = list(converged = TRUE)), class = "bmfrm_fit")
  fs_fast_sum <- fair_scores(obj2, summary = TRUE, return_draws = FALSE)
  fs_fast_row <- fair_scores(obj2, summary = FALSE)
  expect_type(attr(fs_fast_sum, "reference_meta"), "list")
  expect_type(attr(fs_fast_row, "reference_meta"), "list")

  # Verify residuals.bmfrm_fit restores original string labels (not 1..J integers)
  res2 <- residuals(obj2)
  expect_setequal(unique(res2$person), paste0("Person_", 1:5))
  expect_setequal(unique(res2$rater), c("Rater_A", "Rater_B"))
})

test_that("residual correlations deduplicate counts, handle NA plots, and fit/bias/ppc helpers handle edge cases", {
  # 1. Unique unordered pairs in print/diagnose while preserving symmetric cor_tidy
  d_rc <- expand.grid(person = paste0("P", 1:6), item = c("I1", "I2", "I3"), stringsAsFactors = FALSE)
  d_rc$resid <- c(
    -2, -1, 0, 1, 2, 3,
    -2, -1, 0, 1, 2, 3,
    rep(0, 6) # zero variance -> NA correlations with I3!
  )
  rc <- suppressWarnings(residual_cor_criteria(d_rc))
  expect_equal(nrow(rc$cor_tidy), 6L) # all ordered off-diagonal pairs kept
  diag_rc <- diagnose_residual_cors(rc, threshold = 0.2)
  expect_equal(diag_rc$summary$n_cors, 1L)
  expect_equal(diag_rc$summary$n_high_cors, 1L)
  expect_equal(nrow(diag_rc$high_cors), 1L)
  out_txt <- utils::capture.output(print(rc))
  expect_match(paste(out_txt, collapse = "\n"), "1 correlations > \\|0\\.3\\| detected")
  # Plot with NA entries and hclust must fall back cleanly instead of crashing
  tf_pdf <- tempfile(fileext = ".pdf")
  grDevices::pdf(tf_pdf)
  on.exit(unlink(tf_pdf), add = TRUE)
  expect_no_error(plot(rc, order = "hclust"))
  grDevices::dev.off()

  # 2. extract_residuals() excludes 2D basis matrices even when length(B_item) == N, and add_fit_cri_bayes checks draws
  stan_dat_2 <- list(
    N = 2L,
    X = c(1L, 2L),
    person_id = c(1L, 2L),
    item_id = c(1L, 2L),
    B_item = matrix(c(-0.7071068, 0.7071068), nrow = 2L, ncol = 1L)
  )
  mock_gq_fit <- list(
    summary = function(variables) {
      if (variables == "mu") return(tibble::tibble(variable = c("mu[1]", "mu[2]"), mean = c(1.25, 1.75)))
      if (variables == "sigma2") return(tibble::tibble(variable = c("sigma2[1]", "sigma2[2]"), mean = c(0.2, 0.2)))
      stop("unexpected")
    }
  )
  res_no_draws <- extract_residuals(mock_gq_fit, stan_dat_2, save_draws = FALSE)
  expect_false("B_item" %in% names(res_no_draws))
  expect_true(all(c("person_id", "item_id") %in% names(res_no_draws)))
  fs_tbl <- facet_infit_outfit(res_no_draws, person_id)
  expect_error(add_fit_cri_bayes(fs_tbl, res_no_draws, person_id), "save_draws = TRUE")

  # 3. summarise_bias_ir numeric order for default labels vs alphabetical for explicit labels
  D_b <- 20L
  cn_b <- as.vector(outer(1:10, 1:2, function(i, r) paste0("bias[", i, ",", r, "]")))
  mat_b <- matrix(0.1, nrow = D_b, ncol = length(cn_b), dimnames = list(NULL, cn_b))
  mock_b <- list(
    metadata = function() list(stan_variables = "bias"),
    draws = function(variables, format = "df") as.data.frame(mat_b)
  )
  sb_def <- summarise_bias_ir(mock_b, diagnostics = FALSE)
  expect_identical(unique(sb_def$item), paste0("Item_", 1:10))
  sb_exp <- summarise_bias_ir(mock_b, item_labels = sprintf("Lbl_%02d", 10:1), rater_labels = c("RB", "RA"), diagnostics = FALSE)
  expect_identical(unique(sb_exp$item), sprintf("Lbl_%02d", 1:10))

  # 4. run_ppc handles singleton group in variance PPC without ggplot missing-value warning
  y_ppc <- c(1L, 2L, 3L, 2L, 1L)
  yrep_ppc <- matrix(rep(y_ppc, each = 10L), nrow = 10L, ncol = 5L)
  colnames(yrep_ppc) <- paste0("y_rep[", 1:5, "]")
  mock_ppc <- list(draws = function(variables, format = "matrix") yrep_ppc)
  grp_ppc <- c("G1", "G1", "G2", "G2", "G_SINGLE")
  expect_no_warning(suppressMessages(ppc_out <- run_ppc(mock_ppc, y_ppc, group = grp_ppc, ndraws = 5L)))
  expect_s3_class(ppc_out$ppc_var, "ggplot")
})

test_that("predict rejects non-hierarchical 1-known-level interactions and cv_bmfrm handles large seeds and strict_target = FALSE", {
  # 1. Non-hierarchical interaction with 1 known level on f1 must NOT silently zero out under sample_new_levels = "hierarchical"
  d_1i <- expand.grid(person = paste0("P", 1:5), item = "I1", rater = c("R1", "R2"), stringsAsFactors = FALSE)
  d_1i$score <- rep(1:3, length.out = nrow(d_1i))
  spec_1i <- resolve_bmfrm_spec(parse_bmfrm_formula(
    score ~ person + item + rater + item:rater, d_1i, K = 3L,
    hierarchical = c("person", "rater") # item:rater intentionally NOT hierarchical!
  ))
  di_1i <- prepare_data_bmfrm(spec_1i, d_1i)
  D <- 16L
  mk_m <- function(nm, v) {
    m <- matrix(rep(v, each = D), nrow = D, ncol = length(v))
    colnames(m) <- paste0(nm, "[", seq_along(v), "]")
    m
  }
  store_1i <- list(
    theta = mk_m("theta", seq(-0.4, 0.4, length.out = 5L)),
    mu_theta = matrix(0, D, 1L, dimnames = list(NULL, "mu_theta")),
    sigma_theta = matrix(0.5, D, 1L, dimnames = list(NULL, "sigma_theta")),
    item = mk_m("item", 0),
    rater = mk_m("rater", c(-0.2, 0.2)),
    sigma_rater = matrix(0.4, D, 1L, dimnames = list(NULL, "sigma_rater")),
    itemrater = matrix(0, D, 2L, dimnames = list(NULL, c("itemrater[1,1]", "itemrater[1,2]"))),
    tau = mk_m("tau", c(-0.5, 0.5))
  )
  fit_1i <- structure(list(
    fit = list(
      num_chains = function() 2L,
      metadata = function() list(stan_variables = names(store_1i), iter_sampling = 8L),
      draws = function(var, format = "matrix") store_1i[[var]]
    ),
    spec = spec_1i, data_info = di_1i, convergence = list(converged = TRUE)
  ), class = "bmfrm_fit")
  nd_new_r <- data.frame(person = "P1", item = "I1", rater = "R_NEW", stringsAsFactors = FALSE)
  expect_error(
    predict(fit_1i, newdata = nd_new_r, allow_new_levels = TRUE, sample_new_levels = "hierarchical"),
    "non-hierarchical term 'item:rater'"
  )

  # 2. Near-INT_MAX seed safety & strict_target = FALSE marginalization in cv_bmfrm()
  expect_true(is.finite(BayesMFRM:::offset_bmfrm_seed(2147480000L, 10007L)))
  expect_gt(BayesMFRM:::offset_bmfrm_seed(2147480000L, 10007L), 0L)

  make_cv_mock <- function(df_in, hier = c("person", "rater"), seed_val = 101L) {
    set.seed(seed_val)
    sp <- resolve_bmfrm_spec(parse_bmfrm_formula(
      score ~ person + item + rater, df_in, K = 3L, hierarchical = hier
    ))
    di <- suppressWarnings(prepare_data_bmfrm(sp, df_in))
    dc <- check_design_bmfrm(score ~ person + item + rater, di$data_clean, K = 3L)
    Jp <- length(di$facet_labels$person)
    Ji <- length(di$facet_labels$item)
    Jr <- length(di$facet_labels$rater)
    D_cv <- 32L
    st <- list(
      theta = matrix(stats::rnorm(D_cv * Jp, 0, 0.3), D_cv, Jp, dimnames = list(NULL, paste0("theta[", seq_len(Jp), "]"))),
      item = {
        m <- matrix(stats::rnorm(D_cv * Ji, 0, 0.2), D_cv, Ji)
        if (Ji > 1L) m <- sweep(m, 1L, rowMeans(m), "-") else m[] <- 0
        colnames(m) <- paste0("item[", seq_len(Ji), "]")
        m
      },
      rater = {
        m <- matrix(stats::rnorm(D_cv * Jr, 0, 0.25), D_cv, Jr)
        if (Jr > 1L) m <- sweep(m, 1L, rowMeans(m), "-") else m[] <- 0
        colnames(m) <- paste0("rater[", seq_len(Jr), "]")
        m
      },
      tau = {
        m <- matrix(rep(c(-0.5, 0.5), each = D_cv), D_cv, 2L)
        colnames(m) <- c("tau[1]", "tau[2]")
        m
      }
    )
    if ("person" %in% hier) {
      st$mu_theta <- matrix(stats::rnorm(D_cv, 0, 0.1), D_cv, 1L, dimnames = list(NULL, "mu_theta"))
      st$sigma_theta <- matrix(rep(0.6, D_cv), D_cv, 1L, dimnames = list(NULL, "sigma_theta"))
    }
    if ("rater" %in% hier) {
      st$sigma_rater <- matrix(rep(0.5, D_cv), D_cv, 1L, dimnames = list(NULL, "sigma_rater"))
    }
    structure(list(
      fit = list(
        num_chains = function() 4L,
        metadata = function() list(stan_variables = names(st), iter_sampling = 8L),
        draws = function(var, format = "matrix") {
          m <- st[[var]]
          if (format == "draws_array") return(array(m, c(8L, 4L, ncol(m)), dimnames = list(NULL, NULL, colnames(m))))
          m
        }
      ),
      spec = sp, data_info = di, design = dc, convergence = list(converged = TRUE)
    ), class = "bmfrm_fit")
  }

  # R4 connects to P1..P2 AND is the sole rater of P_ONLY_R4:
  # all training folds stay connected & identifiable, while holding out R4 makes P_ONLY_R4 collateral unseen!
  d_collat <- rbind(
    expand.grid(person = paste0("P", 1:4), item = c("I1", "I2"), rater = c("R1", "R2", "R3"), stringsAsFactors = FALSE),
    expand.grid(person = c("P1", "P2", "P_ONLY_R4"), item = c("I1", "I2"), rater = "R4", stringsAsFactors = FALSE)
  )
  d_collat$score <- rep(1:3, length.out = nrow(d_collat))

  # Collateral unseen person when person is NOT hierarchical MUST fail pre-fit (0 refit calls)
  refit_calls <- 0L
  obj_nonhier_p <- make_cv_mock(d_collat, hier = "rater")
  expect_error(
    cv_bmfrm(
      obj_nonhier_p, target = "rater", K_folds = 2L, strict_target = FALSE, seed = 2147480000L,
      refit_fn = function(train_df, fold, spec, fit_args) {
        refit_calls <<- refit_calls + 1L
        make_cv_mock(train_df, hier = spec$hierarchical, seed_val = fold)
      }
    ),
    "Missing from hierarchical: person"
  )
  expect_equal(refit_calls, 0L)

  # When both person and rater ARE hierarchical, strict_target = FALSE marginalizes cleanly for both rating and grouped CV!
  obj_hier_both <- make_cv_mock(d_collat, hier = c("person", "rater"), seed_val = 10L)
  obj_hier_alt <- make_cv_mock(d_collat, hier = c("person", "rater"), seed_val = 20L)
  refit_ok <- function(train_df, fold, spec, fit_args) {
    make_cv_mock(train_df, hier = spec$hierarchical, seed_val = 300L + fold)
  }

  folds_rat_unseen <- ifelse(d_collat$person == "P_ONLY_R4", 1L, 2L)
  cv_rat_u1 <- cv_bmfrm(obj_hier_both, target = "rating", folds = folds_rat_unseen, M = 15L,
                        seed = 2147480000L, strict_target = FALSE, chunk_size = 50L, refit_fn = refit_ok)
  cv_rat_u2 <- cv_bmfrm(obj_hier_both, target = "rating", folds = folds_rat_unseen, M = 15L,
                        seed = 2147480000L, strict_target = FALSE, chunk_size = 1L, refit_fn = refit_ok)
  expect_equal(cv_rat_u1$actual_target, "new_rating_with_unseen_levels")
  expect_true(is.finite(cv_rat_u1$elpd_pointwise))
  expect_equal(cv_rat_u1$elpd_pointwise, cv_rat_u2$elpd_pointwise, tolerance = 1e-12)

  folds_grp_collat <- ifelse(d_collat$rater %in% c("R1", "R4"), 1L, 2L)
  cv_grp_c1 <- cv_bmfrm(obj_hier_both, target = "rater", folds = folds_grp_collat, M = 20L,
                        seed = 2147480000L, strict_target = FALSE, chunk_size = 50L, refit_fn = refit_ok)
  cv_grp_c1_ch1 <- cv_bmfrm(obj_hier_both, target = "rater", folds = folds_grp_collat, M = 20L,
                            seed = 2147480000L, strict_target = FALSE, chunk_size = 1L, refit_fn = refit_ok)
  cv_grp_c2 <- cv_bmfrm(obj_hier_alt, target = "rater", folds = folds_grp_collat, M = 20L,
                        seed = 2147480000L, strict_target = FALSE, chunk_size = 50L, refit_fn = refit_ok)
  expect_equal(cv_grp_c1$actual_target, "new_rater_with_unseen_person")
  expect_true(is.finite(cv_grp_c1$elpd_group))
  expect_true(is.finite(cv_grp_c1$elpd_pointwise))
  expect_equal(cv_grp_c1$elpd_group, cv_grp_c1_ch1$elpd_group, tolerance = 1e-12)
  expect_equal(cv_grp_c1$elpd_pointwise, cv_grp_c1_ch1$elpd_pointwise, tolerance = 1e-12)
  expect_s3_class(compare_bmfrm(m1 = cv_grp_c1, m2 = cv_grp_c2), "bmfrm_model_comparison")
})

test_that("review regressions: bmfrm refit=never, all-NA suitability, legacy 4-facet param_map, seed bounds, RNG preservation, chunked memory, and 1-level partner interactions", {
  # 1. Public bmfrm(..., refit = "never") stops before sampling when data changes (missing new signature cache)
  sampler_calls <- 0L
  tmp_cache <- tempfile("bmfrm-never-")
  dir.create(tmp_cache)
  on.exit(unlink(tmp_cache, recursive = TRUE), add = TRUE)
  mock_stan_model <- list(
    code = function() "parameters { real z; }",
    sample = function(data, output_dir, ...) {
      sampler_calls <<- sampler_calls + 1L
      csv_p <- file.path(output_dir, paste0("c", sampler_calls, ".csv"))
      writeLines("stub", csv_p)
      list(
        output_files = function() csv_p,
        num_chains = function() 1L,
        metadata = function() list(stan_variables = c("theta", "rater", "tau"), iter_sampling = 4L),
        draws = function(var = NULL, format = "matrix") {
          st <- list(
            theta = matrix(0, 4L, 2L, dimnames = list(NULL, c("theta[1]", "theta[2]"))),
            rater = matrix(0, 4L, 2L, dimnames = list(NULL, c("rater[1]", "rater[2]"))),
            tau = matrix(0, 4L, 1L, dimnames = list(NULL, "tau[1]"))
          )
          if (is.null(var)) return(do.call(cbind, st))
          st[[var]]
        }
      )
    }
  )
  testthat::local_mocked_bindings(
    check_cmdstan = function() TRUE,
    compile_bmfrm = function(stan_code, model_name, cache_dir, ...) {
      list(
        model = mock_stan_model,
        hash = "h1",
        stan_file = file.path(cache_dir, "m.stan"),
        cache_file = file.path(cache_dir, "m.rds")
      )
    },
    check_convergence = function(...) list(converged = TRUE, status = "ok"),
    .package = "BayesMFRM"
  )
  d_a <- data.frame(person = c("P1", "P1", "P2", "P2"), rater = c("R1", "R2", "R1", "R2"), score = c(1L, 2L, 1L, 2L))
  d_b <- d_a
  d_b$score <- c(2L, 1L, 2L, 1L)
  fit_a <- suppressWarnings(suppressMessages(bmfrm(score ~ person + rater, d_a, K = 2L, cache_dir = tmp_cache, chains = 1L, iter = 10L, warmup = 5L, seed = 1L)))
  expect_equal(sampler_calls, 1L)
  expect_error(
    suppressWarnings(suppressMessages(bmfrm(score ~ person + rater, d_b, K = 2L, cache_dir = tmp_cache, refit = "never", chains = 1L, iter = 10L, warmup = 5L, seed = 1L))),
    "cannot be verified"
  )
  expect_equal(sampler_calls, 1L)

  # 2. check_data_suitability() handles an all-NA facet without crashing on design_sparsity
  d_all_na <- data.frame(
    person = rep(NA_character_, 4L),
    rater = c("R1", "R2", "R1", "R2"),
    score = c(1L, 2L, 1L, 2L),
    stringsAsFactors = FALSE
  )
  expect_no_error(suit_na <- check_data_suitability(d_all_na, "score", c("person", "rater")))
  expect_equal(suit_na$facet_info$person$n_levels, 0L)
  expect_true(is.na(suit_na$design_sparsity))
  expect_true(any(grepl("Facet person has no valid levels", suit_na$warnings)))
  expect_no_error(invisible(utils::capture.output(print(suit_na))))

  # 3. Legacy saved 4-facet fit lacking param_map in spec$eval_contract lazily backfills item -> beta
  f4_code <- paste(readLines(bmfrm_model("four-facet/person-item-rater-scenario-main")$path), collapse = "\n")
  d4_leg <- expand.grid(person = c("P1", "P2"), item = c("I1", "I2"), rater = c("R1", "R2"), scenario = c("S1", "S2"), stringsAsFactors = FALSE)
  d4_leg$score <- rep(1:2, length.out = nrow(d4_leg))
  spec4_leg <- parse_bmfrm_formula(score ~ person + item + rater + scenario, d4_leg, K = 2L)
  spec4_leg$custom_stan <- TRUE
  # Simulate pre-patch saved eval_contract WITHOUT param_map
  spec4_leg$eval_contract <- list(
    source = "bundled_verified",
    model_id = "four-facet/person-item-rater-scenario-main",
    r_recompute = FALSE,
    stan_mu_fair_adjust = c("rater", "scenario"),
    supported_ops = c("residuals", "facet_fit", "fair_scores", "pp_check", "loo"),
    reason = "legacy"
  )
  di4_leg <- prepare_data_bmfrm(spec4_leg, d4_leg)
  obj4_leg <- structure(list(
    fit = list(
      num_chains = function() 1L,
      metadata = function() list(stan_variables = c("theta", "beta", "rater", "scenario", "tau"), iter_sampling = 4L),
      summary = function(variables, ...) {
        if (!identical(variables, "beta")) stop("Expected beta, got ", variables)
        tibble::tibble(variable = c("beta[1]", "beta[2]"), mean = c(-0.3, 0.3), sd = c(0.1, 0.1), `5%` = c(-0.4, 0.2), `95%` = c(-0.2, 0.4))
      }
    ),
    spec = spec4_leg, data_info = di4_leg, stan_code = f4_code, convergence = list(converged = TRUE)
  ), class = "bmfrm_fit")
  fs_leg <- facet_summary(obj4_leg, "item", probs = c(0.05, 0.95))
  expect_equal(fs_leg$mean, c(-0.3, 0.3))

  # 4. offset_bmfrm_seed() boundary invariants (0, 1, INT_MAX, wraparound)
  imax <- .Machine$integer.max
  expect_identical(BayesMFRM:::offset_bmfrm_seed(0L, 0L), 0L)
  expect_identical(BayesMFRM:::offset_bmfrm_seed(1L, 0L), 1L)
  expect_identical(BayesMFRM:::offset_bmfrm_seed(imax, 0L), imax)
  expect_identical(BayesMFRM:::offset_bmfrm_seed(imax, 1L), 0L)
  expect_identical(BayesMFRM:::offset_bmfrm_seed(imax, 2L), 1L)
  expect_error(BayesMFRM:::offset_bmfrm_seed(-1L, 0L), "non-negative 32-bit integer")
  expect_error(BayesMFRM:::offset_bmfrm_seed(as.numeric(imax) + 1, 0L), "non-negative 32-bit integer")

  # 5. Degenerate 1-level training partner interaction (2 known items x 1 known rater, new rater R2):
  #    Known f1 (2 items) x New f2 (R2) MUST preserve sum-to-zero uncertainty across items, matching exact quadrature!
  d_deg_tr <- expand.grid(person = c("P1", "P2", "P3"), item = c("I1", "I2"), rater = "R1", stringsAsFactors = FALSE)
  d_deg_tr$score <- rep(1:2, length.out = nrow(d_deg_tr))
  sp_deg <- resolve_bmfrm_spec(parse_bmfrm_formula(
    score ~ person + item + rater + item:rater, d_deg_tr, K = 2L,
    hierarchical = c("person", "rater", "item:rater")
  ))
  di_deg <- suppressWarnings(prepare_data_bmfrm(sp_deg, d_deg_tr))
  D_deg <- 2000L
  st_deg <- list(
    theta = matrix(0, D_deg, 3L, dimnames = list(NULL, paste0("theta[", 1:3, "]"))),
    mu_theta = matrix(0, D_deg, 1L, dimnames = list(NULL, "mu_theta")),
    sigma_theta = matrix(0, D_deg, 1L, dimnames = list(NULL, "sigma_theta")),
    item = matrix(0, D_deg, 2L, dimnames = list(NULL, c("item[1]", "item[2]"))),
    rater = matrix(0, D_deg, 1L, dimnames = list(NULL, "rater[1]")),
    sigma_rater = matrix(0, D_deg, 1L, dimnames = list(NULL, "sigma_rater")),
    itemrater = matrix(0, D_deg, 2L, dimnames = list(NULL, c("itemrater[1,1]", "itemrater[2,1]"))),
    sigma_itemrater = matrix(2.0, D_deg, 1L, dimnames = list(NULL, "sigma_itemrater")),
    tau = matrix(0, D_deg, 1L, dimnames = list(NULL, "tau[1]"))
  )
  fit_deg <- structure(list(
    fit = list(
      num_chains = function() 4L,
      metadata = function() list(stan_variables = names(st_deg), iter_sampling = D_deg %/% 4L),
      draws = function(var, format = "matrix") {
        m <- st_deg[[var]]
        if (format == "draws_array") return(array(m, c(D_deg %/% 4L, 4L, ncol(m)), dimnames = list(NULL, NULL, colnames(m))))
        m
      }
    ),
    spec = sp_deg, data_info = di_deg, convergence = list(converged = TRUE)
  ), class = "bmfrm_fit")

  nd_r2 <- data.frame(person = c("P1", "P1"), item = c("I1", "I2"), rater = c("R2", "R2"), score = c(1L, 1L), stringsAsFactors = FALSE)
  set.seed(777L)
  rng_before <- .Random.seed
  lp_r2 <- predict(fit_deg, newdata = nd_r2, type = "linpred", summary = FALSE,
                   allow_new_levels = TRUE, sample_new_levels = "hierarchical", seed = 42L)
  expect_identical(.Random.seed, rng_before)
  # Sum across the 2 known items for new rater R2 is identically 0; item contrast SD is 2 * sqrt(2) ~ 2.828427
  expect_equal(max(abs(lp_r2[, 1L] + lp_r2[, 2L])), 0, tolerance = 1e-12)
  expect_equal(stats::sd(lp_r2[, 1L] - lp_r2[, 2L]), 2 * sqrt(2), tolerance = 0.12)

  # Compare score_grouped_fold and score_general_unseen_fold against exact 1D numerical quadrature (-1.706054)
  exact_quad <- log(stats::integrate(
    function(z) stats::plogis(sqrt(2) * z) * stats::plogis(-sqrt(2) * z) * stats::dnorm(z),
    lower = -Inf, upper = Inf
  )$value)
  sc_grp_deg <- BayesMFRM:::score_grouped_fold(fit_deg, nd_r2, target_facet = "rater", M = 1500L, seed = 123L)
  expect_identical(.Random.seed, rng_before)
  sc_gen_deg <- BayesMFRM:::score_general_unseen_fold(fit_deg, nd_r2, target_facet = "rater", M = 1500L, seed = 123L)
  expect_identical(.Random.seed, rng_before)
  expect_equal(sc_grp_deg$group_scores$elpd, exact_quad, tolerance = 0.02)
  expect_equal(sc_gen_deg$group_scores$elpd, exact_quad, tolerance = 0.02)

  # 6. Rating-level unseen CV respects chunk_size = 1L strictly (1 row per chunk) and preserves caller .Random.seed
  seen_chunk_sizes <- integer(0)
  orig_eval_chunk <- BayesMFRM:::eval_bmfrm_chunk_from_plan
  testthat::local_mocked_bindings(
    eval_bmfrm_chunk_from_plan = function(plan, rows, ...) {
      seen_chunk_sizes <<- c(seen_chunk_sizes, length(rows))
      orig_eval_chunk(plan, rows, ...)
    },
    .package = "BayesMFRM"
  )
  nd_many <- data.frame(
    person = rep("P_UNSEEN", 6L),
    item = rep(c("I1", "I2"), 3L),
    rater = rep("R1", 6L),
    score = rep(1:2, 3L),
    stringsAsFactors = FALSE
  )
  sc_rat_ch1 <- BayesMFRM:::score_general_unseen_fold(fit_deg, nd_many, target_facet = NULL, M = 5L, seed = 0L, chunk_size = 1L)
  expect_identical(.Random.seed, rng_before)
  expect_true(length(seen_chunk_sizes) == 6L * 5L && all(seen_chunk_sizes == 1L))
  sc_rat_ch6 <- BayesMFRM:::score_general_unseen_fold(fit_deg, nd_many, target_facet = NULL, M = 5L, seed = 0L, chunk_size = 6L)
  expect_equal(sc_rat_ch1$lpd_pointwise, sc_rat_ch6$lpd_pointwise, tolerance = 1e-12)

  # Public cv_bmfrm() with seed = 0L and seed = INT_MAX preserves caller .Random.seed
  d_full_deg <- rbind(
    d_deg_tr,
    data.frame(person = c("P1", "P2", "P_R2"), item = c("I1", "I2", "I1"), rater = "R2", score = c(1L, 2L, 1L), stringsAsFactors = FALSE)
  )
  sp_full <- resolve_bmfrm_spec(parse_bmfrm_formula(
    score ~ person + item + rater + item:rater, d_full_deg, K = 2L,
    hierarchical = c("person", "rater", "item:rater")
  ))
  di_full <- suppressWarnings(prepare_data_bmfrm(sp_full, d_full_deg))
  D_sm <- 16L
  mk_zero_fit <- function(df_sub) {
    sp_s <- resolve_bmfrm_spec(parse_bmfrm_formula(
      score ~ person + item + rater + item:rater, df_sub, K = 2L,
      hierarchical = c("person", "rater", "item:rater")
    ))
    di_s <- suppressWarnings(prepare_data_bmfrm(sp_s, df_sub))
    dc_s <- check_design_bmfrm(score ~ person + item + rater + item:rater, di_s$data_clean, K = 2L)
    Jp <- length(di_s$facet_labels$person)
    Ji <- length(di_s$facet_labels$item)
    Jr <- length(di_s$facet_labels$rater)
    st_s <- list(
      theta = matrix(0, D_sm, Jp, dimnames = list(NULL, paste0("theta[", seq_len(Jp), "]"))),
      mu_theta = matrix(0, D_sm, 1L, dimnames = list(NULL, "mu_theta")),
      sigma_theta = matrix(0.4, D_sm, 1L, dimnames = list(NULL, "sigma_theta")),
      item = matrix(0, D_sm, Ji, dimnames = list(NULL, paste0("item[", seq_len(Ji), "]"))),
      rater = matrix(0, D_sm, Jr, dimnames = list(NULL, paste0("rater[", seq_len(Jr), "]"))),
      sigma_rater = matrix(0.4, D_sm, 1L, dimnames = list(NULL, "sigma_rater")),
      itemrater = matrix(0, D_sm, Ji * Jr, dimnames = list(NULL, as.vector(outer(seq_len(Ji), seq_len(Jr), function(a, b) paste0("itemrater[", a, ",", b, "]"))))),
      sigma_itemrater = matrix(1.0, D_sm, 1L, dimnames = list(NULL, "sigma_itemrater")),
      tau = matrix(0, D_sm, 1L, dimnames = list(NULL, "tau[1]"))
    )
    structure(list(
      fit = list(
        num_chains = function() 4L,
        metadata = function() list(stan_variables = names(st_s), iter_sampling = 4L),
        draws = function(var, format = "matrix") {
          m <- st_s[[var]]
          if (format == "draws_array") return(array(m, c(4L, 4L, ncol(m)), dimnames = list(NULL, NULL, colnames(m))))
          m
        }
      ),
      spec = sp_s, data_info = di_s, design = dc_s, convergence = list(converged = TRUE)
    ), class = "bmfrm_fit")
  }
  obj_full_deg <- mk_zero_fit(d_full_deg)
  for (s_edge in c(0L, imax)) {
    cv_edge <- suppressWarnings(cv_bmfrm(
      obj_full_deg, target = "rater", K_folds = 2L, M = 10L, seed = s_edge,
      strict_target = FALSE, allow_degenerate_interaction = TRUE,
      refit_fn = function(train_df, fold, spec, fit_args) mk_zero_fit(train_df)
    ))
    expect_identical(.Random.seed, rng_before)
    expect_true(is.finite(cv_edge$elpd_group))
  }
})

test_that("round-2 review regressions: 1-level hierarchical convergence vars, linear CV seed work, chunked score_grouped_fold, and response suitability", {
  # 1. Convergence diagnostics include sigma_rater and c2_itemrater even when rater has 1 training level
  d_1r <- expand.grid(person = paste0("P", 1:4), item = c("I1", "I2"), rater = "R1", stringsAsFactors = FALSE)
  d_1r$score <- rep(1:3, length.out = nrow(d_1r))
  sp_1r <- parse_bmfrm_formula(
    score ~ person + item + rater + item:rater, d_1r, K = 3L,
    hierarchical = c("person", "rater", "item:rater")
  )
  dp_1r <- get_default_priors(sp_1r$facets_main, sp_1r$facets_bias, sp_1r$family, sp_1r$hierarchical)
  dp_1r <- Filter(function(p) !(identical(p$class, "sd") && identical(p$facet, "item:rater")), dp_1r)
  pr_hs <- c(dp_1r, list(prior("normal(0, 0.5)", class = "sd", facet = "item:rater", scale = "raw", local_shrinkage = "horseshoe")))
  sp_1r <- resolve_bmfrm_spec(sp_1r, pr_hs)
  di_1r <- suppressWarnings(prepare_data_bmfrm(sp_1r, d_1r, priors = pr_hs))

  dvars <- BayesMFRM:::model_diagnostic_variables(sp_1r, di_1r)
  expect_true(all(c("sigma_rater", "sigma_rater_raw", "c2_itemrater", "sigma_itemrater", "sigma_itemrater_raw") %in% dvars))
  expect_false(any(c("rater", "rater_raw", "itemrater", "itemrater_raw", "lambda_itemrater") %in% dvars))

  mock_conv_fit <- list(
    num_chains = function() 4L,
    diagnostic_summary = function(quiet = TRUE) list(
      ebfmi = rep(0.9, 4L), num_divergent = rep(0L, 4L), num_max_treedepth = rep(0L, 4L)
    ),
    summary = function(variables, ...) {
      bad <- variables %in% c("sigma_rater", "sigma_rater_raw", "c2_itemrater")
      tibble::tibble(
        variable = variables,
        mean = 0.5,
        sd = 0.1,
        rhat = ifelse(bad, 1.5, 1.001),
        ess_bulk = ifelse(bad, 5, 800),
        ess_tail = ifelse(bad, 5, 800),
        mcse_mean = 0.01
      )
    }
  )
  conv_res <- BayesMFRM:::check_convergence(mock_conv_fit, variables = dvars)
  expect_false(conv_res$converged)
  expect_equal(conv_res$max_rhat, 1.5)
  expect_equal(conv_res$min_ess_bulk, 5)

  # 2. Linear (non-quadratic) seed-mixing work in score_general_unseen_fold() + cross-chunk level consistency
  d_base_fit <- expand.grid(person = paste0("P", 1:4), item = c("I1", "I2"), rater = c("R1", "R2"), stringsAsFactors = FALSE)
  d_base_fit$score <- rep(1:3, length.out = nrow(d_base_fit))
  sp_b <- resolve_bmfrm_spec(parse_bmfrm_formula(score ~ person + item + rater, d_base_fit, K = 3L, hierarchical = c("person", "rater")))
  di_b <- prepare_data_bmfrm(sp_b, d_base_fit)
  D_b <- 32L
  st_b <- list(
    theta = matrix(0, D_b, 4L, dimnames = list(NULL, paste0("theta[", 1:4, "]"))),
    mu_theta = matrix(0, D_b, 1L, dimnames = list(NULL, "mu_theta")),
    sigma_theta = matrix(0.7, D_b, 1L, dimnames = list(NULL, "sigma_theta")),
    item = matrix(rep(c(-0.2, 0.2), each = D_b), D_b, 2L, dimnames = list(NULL, c("item[1]", "item[2]"))),
    rater = matrix(rep(c(-0.1, 0.1), each = D_b), D_b, 2L, dimnames = list(NULL, c("rater[1]", "rater[2]"))),
    sigma_rater = matrix(0.5, D_b, 1L, dimnames = list(NULL, "sigma_rater")),
    tau = matrix(rep(c(-0.4, 0.4), each = D_b), D_b, 2L, dimnames = list(NULL, c("tau[1]", "tau[2]")))
  )
  fit_b <- structure(list(
    fit = list(
      num_chains = function() 4L,
      metadata = function() list(stan_variables = names(st_b), iter_sampling = 8L),
      draws = function(var, format = "matrix") {
        m <- st_b[[var]]
        if (format == "draws_array") return(array(m, c(8L, 4L, ncol(m)), dimnames = list(NULL, NULL, colnames(m))))
        m
      }
    ),
    spec = sp_b, data_info = di_b, convergence = list(converged = TRUE)
  ), class = "bmfrm_fit")

  mix_calls <- 0L
  orig_mix <- BayesMFRM:::mix_bmfrm_seed
  count_mix_for_N <- function(n_unseen) {
    mix_calls <<- 0L
    nd_u <- data.frame(
      person = paste0("NEW_P_", seq_len(n_unseen)),
      item = rep(c("I1", "I2"), length.out = n_unseen),
      rater = rep(c("R1", "R2"), length.out = n_unseen),
      score = rep(1:3, length.out = n_unseen),
      stringsAsFactors = FALSE
    )
    testthat::with_mocked_bindings(
      mix_bmfrm_seed = function(base_seed, salt = NULL) {
        mix_calls <<- mix_calls + length(base_seed)
        orig_mix(base_seed, salt)
      },
      BayesMFRM:::score_general_unseen_fold(fit_b, nd_u, target_facet = NULL, M = 5L, seed = 11L, chunk_size = 1L),
      .package = "BayesMFRM"
    )
    mix_calls
  }
  c100 <- count_mix_for_N(100L)
  c200 <- count_mix_for_N(200L)
  expect_equal(c100, 500L)
  expect_equal(c200, 1000L) # Exact 2x linear scaling (never 4x quadratic!)

  # Cross-chunk consistency when the same unseen person spans multiple chunks
  nd_shared_p <- data.frame(
    person = c("P_U1", "P_U2", "P_U1", "P_U2"),
    item = c("I1", "I1", "I2", "I2"),
    rater = c("R1", "R2", "R1", "R2"),
    score = c(1L, 3L, 2L, 2L),
    stringsAsFactors = FALSE
  )
  sc_u_c1 <- BayesMFRM:::score_general_unseen_fold(fit_b, nd_shared_p, target_facet = NULL, M = 25L, seed = 77L, chunk_size = 1L)
  sc_u_c4 <- BayesMFRM:::score_general_unseen_fold(fit_b, nd_shared_p, target_facet = NULL, M = 25L, seed = 77L, chunk_size = 4L)
  expect_equal(sc_u_c1$lpd_pointwise, sc_u_c4$lpd_pointwise, tolerance = 1e-12)

  # 3. Standard score_grouped_fold() is strictly chunked in observation memory (0 D x N_group allocations at chunk_size = 1)
  D_mem <- 256L
  st_mem <- list(
    theta = matrix(0.1, D_mem, 4L, dimnames = list(NULL, paste0("theta[", 1:4, "]"))),
    mu_theta = matrix(0, D_mem, 1L, dimnames = list(NULL, "mu_theta")),
    sigma_theta = matrix(0.6, D_mem, 1L, dimnames = list(NULL, "sigma_theta")),
    item = matrix(rep(c(-0.2, 0.2), each = D_mem), D_mem, 2L, dimnames = list(NULL, c("item[1]", "item[2]"))),
    rater = matrix(rep(c(-0.1, 0.1), each = D_mem), D_mem, 2L, dimnames = list(NULL, c("rater[1]", "rater[2]"))),
    sigma_rater = matrix(0.5, D_mem, 1L, dimnames = list(NULL, "sigma_rater")),
    tau = matrix(rep(c(-0.4, 0.4), each = D_mem), D_mem, 2L, dimnames = list(NULL, c("tau[1]", "tau[2]")))
  )
  fit_mem <- structure(list(
    fit = list(
      num_chains = function() 4L,
      metadata = function() list(stan_variables = names(st_mem), iter_sampling = 64L),
      draws = function(var, format = "matrix") {
        m <- st_mem[[var]]
        if (format == "draws_array") return(array(m, c(64L, 4L, ncol(m)), dimnames = list(NULL, NULL, colnames(m))))
        m
      }
    ),
    spec = sp_b, data_info = di_b, convergence = list(converged = TRUE)
  ), class = "bmfrm_fit")
  N_grp <- 500L
  nd_grp_big <- data.frame(
    person = rep(paste0("P", 1:4), length.out = N_grp),
    item = rep(c("I1", "I2"), length.out = N_grp),
    rater = rep("R_HELD", N_grp),
    score = rep(1:3, length.out = N_grp),
    stringsAsFactors = FALSE
  )
  grp_c1 <- BayesMFRM:::score_grouped_fold(fit_mem, nd_grp_big, target_facet = "rater", M = 3L, seed = 55L, chunk_size = 1L)
  grp_all <- BayesMFRM:::score_grouped_fold(fit_mem, nd_grp_big, target_facet = "rater", M = 3L, seed = 55L, chunk_size = N_grp)
  expect_equal(grp_c1$group_scores$elpd, grp_all$group_scores$elpd, tolerance = 1e-12)
  expect_equal(grp_c1$lpd_pointwise, grp_all$lpd_pointwise, tolerance = 1e-12)

  if (isTRUE(capabilities("profmem"))) {
    tf_pm <- tempfile(fileext = ".out")
    on.exit(unlink(tf_pm), add = TRUE)
    thresh_bytes <- D_mem * N_grp * 8
    utils::Rprofmem(tf_pm, threshold = thresh_bytes)
    invisible(BayesMFRM:::score_grouped_fold(fit_mem, nd_grp_big, target_facet = "rater", M = 1L, seed = 55L, chunk_size = 1L))
    utils::Rprofmem(NULL)
    pm_lines <- if (file.exists(tf_pm)) readLines(tf_pm, warn = FALSE) else character(0)
    pm_lines <- pm_lines[nzchar(trimws(pm_lines)) & !grepl("^new page:", pm_lines)]
    expect_equal(length(pm_lines), 0L)
  }

  # 4. check_data_suitability() flags character, Inf, all-NA, fractional, and <= 0 responses as unsuitable
  d_bal <- expand.grid(person = paste0("P", 1:4), rater = paste0("R", 1:4), stringsAsFactors = FALSE)
  d_bal$score <- rep(1:4, length.out = nrow(d_bal))
  suit_good <- check_data_suitability(d_bal, "score", c("person", "rater"))
  expect_equal(length(suit_good$warnings), 0L)
  expect_match(suit_good$recommendation, "^Data looks suitable")

  bad_scores <- list(
    chr = rep(c("a", "b"), length.out = 16L),
    inf = c(rep(1, 15L), Inf),
    all_na = rep(NA_real_, 16L),
    frac = rep(1.5, 16L),
    zero = rep(0L, 16L)
  )
  for (nm_bs in names(bad_scores)) {
    d_bad <- d_bal
    d_bad$score <- bad_scores[[nm_bs]]
    s_bad <- check_data_suitability(d_bad, "score", c("person", "rater"))
    expect_gt(length(s_bad$warnings), 0L)
    expect_match(s_bad$recommendation, "^Data is not suitable", label = nm_bs)
  }

  # 5. Fixing both folds = cv_init$folds and refit_seed keeps CV splits comparable across seed/M sensitivity runs
  refit_det <- function(train_df, fold, spec, fit_args) {
    di_tr <- suppressWarnings(prepare_data_bmfrm(spec, train_df))
    dc_tr <- check_design_bmfrm(spec$formula, di_tr$data_clean, K = 3L)
    Jp_tr <- length(di_tr$facet_labels$person)
    Ji_tr <- length(di_tr$facet_labels$item)
    Jr_tr <- length(di_tr$facet_labels$rater)
    st_tr <- list(
      theta = matrix(0, D_b, Jp_tr, dimnames = list(NULL, paste0("theta[", seq_len(Jp_tr), "]"))),
      mu_theta = matrix(0, D_b, 1L, dimnames = list(NULL, "mu_theta")),
      sigma_theta = matrix(0.7, D_b, 1L, dimnames = list(NULL, "sigma_theta")),
      item = matrix(0, D_b, Ji_tr, dimnames = list(NULL, paste0("item[", seq_len(Ji_tr), "]"))),
      rater = matrix(0, D_b, Jr_tr, dimnames = list(NULL, paste0("rater[", seq_len(Jr_tr), "]"))),
      sigma_rater = matrix(0.5, D_b, 1L, dimnames = list(NULL, "sigma_rater")),
      tau = matrix(rep(c(-0.4, 0.4), each = D_b), D_b, 2L, dimnames = list(NULL, c("tau[1]", "tau[2]")))
    )
    structure(list(
      fit = list(
        num_chains = function() 4L,
        metadata = function() list(stan_variables = names(st_tr), iter_sampling = 8L),
        draws = function(var, format = "matrix") {
          m <- st_tr[[var]]
          if (format == "draws_array") return(array(m, c(8L, 4L, ncol(m)), dimnames = list(NULL, NULL, colnames(m))))
          m
        }
      ),
      spec = spec, data_info = di_tr, design = dc_tr, convergence = list(converged = TRUE)
    ), class = "bmfrm_fit")
  }
  cv_s1 <- suppressWarnings(cv_bmfrm(fit_b, target = "rater", K_folds = 2L, M = 10L, seed = 101L, refit_seed = 777L, refit_fn = refit_det))
  cv_s2 <- suppressWarnings(cv_bmfrm(fit_b, target = "rater", folds = cv_s1$folds, M = 20L, seed = 202L, refit_seed = 777L, refit_fn = refit_det))
  expect_identical(cv_s1$folds, cv_s2$folds)
  expect_s3_class(compare_bmfrm(s1 = cv_s1, s2 = cv_s2), "bmfrm_model_comparison")
})

test_that("round-3 and round-4 review regressions: column-wise reproducible interaction evaluation, 512/513 boundary invariance, zero single-use retention, and group cache release", {
  orig_center <- BayesMFRM:::compute_interaction_center_unseeded
  orig_col <- BayesMFRM:::draw_raw_interaction_col_unseeded
  make_bias_grp_fit <- function(J_items, D_draws = 32L, K_cat = 2L, shrinkage = "none") {
    d <- expand.grid(
      person = c("P1", "P2", "P3"),
      item = paste0("I", seq_len(J_items)),
      rater = c("R1", "R2"),
      stringsAsFactors = FALSE
    )
    d$score <- rep(seq_len(K_cat), length.out = nrow(d))
    sp <- parse_bmfrm_formula(
      score ~ person + item + rater + item:rater, d, K = K_cat,
      hierarchical = c("person", "rater", "item:rater")
    )
    pr <- NULL
    if (!identical(shrinkage, "none")) {
      dp <- get_default_priors(sp$facets_main, sp$facets_bias, sp$family, sp$hierarchical)
      dp <- Filter(function(p) !(identical(p$class, "sd") && identical(p$facet, "item:rater")), dp)
      sc <- if (identical(shrinkage, "student_t")) "marginal" else "raw"
      pr <- c(dp, list(prior("normal(0, 0.35)", class = "sd", facet = "item:rater", scale = sc, local_shrinkage = shrinkage)))
    }
    sp <- resolve_bmfrm_spec(sp, pr)
    di <- prepare_data_bmfrm(sp, d, priors = pr)
    item_mat <- sweep(matrix(seq(-0.3, 0.3, length.out = J_items), nrow = D_draws, ncol = J_items, byrow = TRUE), 1L, 0, "-")
    item_mat <- sweep(item_mat, 1L, rowMeans(item_mat), "-")
    colnames(item_mat) <- paste0("item[", seq_len(J_items), "]")
    b_mat <- matrix(0, nrow = D_draws, ncol = J_items * 2L)
    colnames(b_mat) <- as.vector(outer(seq_len(J_items), 1:2, function(i, r) paste0("itemrater[", i, ",", r, "]")))
    tau_mat <- matrix(0, nrow = D_draws, ncol = K_cat - 1L, dimnames = list(NULL, paste0("tau[", seq_len(K_cat - 1L), "]")))
    st <- list(
      theta = matrix(0, D_draws, 3L, dimnames = list(NULL, paste0("theta[", 1:3, "]"))),
      mu_theta = matrix(0, D_draws, 1L, dimnames = list(NULL, "mu_theta")),
      sigma_theta = matrix(0.5, D_draws, 1L, dimnames = list(NULL, "sigma_theta")),
      item = item_mat,
      rater = matrix(rep(c(-0.1, 0.1), each = D_draws), D_draws, 2L, dimnames = list(NULL, c("rater[1]", "rater[2]"))),
      sigma_rater = matrix(0.4, D_draws, 1L, dimnames = list(NULL, "sigma_rater")),
      itemrater = b_mat,
      sigma_itemrater = matrix(0.35, D_draws, 1L, dimnames = list(NULL, "sigma_itemrater")),
      tau = tau_mat
    )
    if (identical(shrinkage, "horseshoe")) {
      st$c2_itemrater <- matrix(2.0, D_draws, 1L, dimnames = list(NULL, "c2_itemrater"))
    }
    structure(list(
      fit = list(
        num_chains = function() 4L,
        metadata = function() list(stan_variables = names(st), iter_sampling = D_draws %/% 4L),
        draws = function(var, format = "matrix") {
          m <- st[[var]]
          if (format == "draws_array") return(array(m, c(D_draws %/% 4L, 4L, ncol(m)), dimnames = list(NULL, NULL, colnames(m))))
          m
        }
      ),
      spec = sp, data_info = di, convergence = list(converged = TRUE)
    ), class = "bmfrm_fit")
  }

  measure_work <- function(expr_fn) {
    center_calls <- 0L
    col_draws <- 0L
    ans <- testthat::with_mocked_bindings(
      compute_interaction_center_unseeded = function(D, n_known, b_info, slice_seed, scope_cols, u_c = integer(0)) {
        center_calls <<- center_calls + 1L
        orig_center(D, n_known, b_info, slice_seed, scope_cols, u_c)
      },
      draw_raw_interaction_col_unseeded = function(D, n_known, b_info, col_seed) {
        col_draws <<- col_draws + 1L
        orig_col(D, n_known, b_info, col_seed)
      },
      expr_fn(),
      .package = "BayesMFRM"
    )
    list(ans = ans, center_calls = center_calls, col_draws = col_draws)
  }

  # 1. Linear scaling across J = 30, 60, 120 at chunk_size = 1 vs J (both grouped and general paths)
  for (J in c(30L, 60L, 120L)) {
    obj_j <- make_bias_grp_fit(J, D_draws = 32L, K_cat = 2L)
    nd_j <- data.frame(
      person = "P1", item = paste0("I", seq_len(J)), rater = "R_NEW",
      score = rep(1:2, length.out = J), stringsAsFactors = FALSE
    )
    m_grp_1 <- measure_work(function() BayesMFRM:::score_grouped_fold(obj_j, nd_j, "rater", M = 3L, seed = 91L, chunk_size = 1L))
    m_grp_w <- measure_work(function() BayesMFRM:::score_grouped_fold(obj_j, nd_j, "rater", M = 3L, seed = 91L, chunk_size = J))
    expect_equal(m_grp_1$center_calls, 3L)
    expect_equal(m_grp_1$col_draws, 3L * (2L * J - 1L))
    expect_equal(m_grp_1$ans, m_grp_w$ans, tolerance = 1e-12)

    m_gen_1 <- measure_work(function() BayesMFRM:::score_general_unseen_fold(obj_j, nd_j, "rater", M = 3L, seed = 91L, chunk_size = 1L))
    m_gen_w <- measure_work(function() BayesMFRM:::score_general_unseen_fold(obj_j, nd_j, "rater", M = 3L, seed = 91L, chunk_size = J))
    expect_equal(m_gen_1$center_calls, 3L)
    expect_equal(m_gen_1$col_draws, 3L * (2L * J - 1L))
    expect_equal(m_gen_1$ans, m_gen_w$ans, tolerance = 1e-12)
  }

  # 2. Exact 512 / 513 / 1000 boundary check: NO cache-wipe thrashing in CV or predict()!
  obj_30 <- make_bias_grp_fit(30L, D_draws = 32L, K_cat = 2L)
  nd_both <- data.frame(
    person = "P_NEW", item = paste0("I", 1:30), rater = "R_NEW",
    score = rep(1:2, 15L), stringsAsFactors = FALSE
  )
  for (M_val in c(100L, 512L, 513L, 1000L)) {
    w1 <- measure_work(function() BayesMFRM:::score_general_unseen_fold(obj_30, nd_both, "rater", M = M_val, seed = 91L, chunk_size = 1L))
    w30 <- measure_work(function() BayesMFRM:::score_general_unseen_fold(obj_30, nd_both, "rater", M = M_val, seed = 91L, chunk_size = 30L))
    expect_equal(w1$center_calls, M_val)
    expect_equal(w30$center_calls, M_val)
    expect_equal(w1$ans$lpd_pointwise, w30$ans$lpd_pointwise, tolerance = 1e-12)
  }

  obj_p4 <- make_bias_grp_fit(4L, D_draws = 32L, K_cat = 2L)
  for (G_val in c(512L, 513L)) {
    np_g <- data.frame(
      person = "P1", item = rep(c("I1", "I2"), each = G_val),
      rater = rep(paste0("NEW_R", seq_len(G_val)), 2L),
      stringsAsFactors = FALSE
    )
    p1 <- measure_work(function() predict(obj_p4, np_g, type = "linpred", allow_new_levels = TRUE, summary = FALSE, seed = 27L, chunk_size = 1L))
    pw <- measure_work(function() predict(obj_p4, np_g, type = "linpred", allow_new_levels = TRUE, summary = FALSE, seed = 27L, chunk_size = nrow(np_g)))
    expect_equal(p1$center_calls, G_val)
    expect_equal(p1$ans, pw$ans, tolerance = 1e-12)
  }

  # 3. Zero cache retention for single-use unseen raters (100 raters x 1 item out of J = 100, D = 256, chunk_size = 1)
  obj_100 <- make_bias_grp_fit(100L, D_draws = 256L, K_cat = 2L)
  nd_once <- data.frame(person = "P1", item = "I1", rater = paste0("NEW_R", 1:100), stringsAsFactors = FALSE)
  orig_eval_ic <- BayesMFRM:::eval_constrained_interaction_chunk
  peak_entries <- 0L
  peak_bytes <- 0
  w_once <- measure_work(function() {
    testthat::with_mocked_bindings(
      eval_constrained_interaction_chunk = function(..., offset_cache = NULL) {
        out <- orig_eval_ic(..., offset_cache = offset_cache)
        if (!is.null(offset_cache)) {
          peak_entries <<- max(peak_entries, length(offset_cache$data))
          peak_bytes <<- max(peak_bytes, offset_cache$bytes)
        }
        out
      },
      predict(obj_100, nd_once, type = "epred", summary = TRUE, allow_new_levels = TRUE, seed = 91L, chunk_size = 1L),
      .package = "BayesMFRM"
    )
  })
  expect_equal(peak_entries, 0L)
  expect_equal(peak_bytes, 0)
  expect_equal(w_once$col_draws, 100L) # 1 column per single-use rater (never 100 * 100 = 10,000!)

  # Also verify direct single-row plan evaluation retains 0 cache entries across M = 128
  np_one <- data.frame(person = "P_NEW", item = "I1", rater = "R_NEW", score = 1L, stringsAsFactors = FALSE)
  plan_one <- BayesMFRM:::with_bmfrm_seed(10L, BayesMFRM:::prepare_eval_plan(obj_100, np_one, allow_new_levels = TRUE))
  for (m in 1:128) invisible(BayesMFRM:::eval_bmfrm_chunk_from_plan(plan_one, 1L, compute_linpred = TRUE, compute_mu = FALSE, seed_salt = m))
  expect_equal(length(plan_one$offset_cache$data), 0L)
  expect_equal(plan_one$offset_cache$bytes, 0)

  # 4. Completed groups in general CV never leave stale cache entries for subsequent groups,
  #    and byte cap + Student-t / horseshoe chunk invariance hold strictly
  small_4 <- make_bias_grp_fit(4L, D_draws = 32L, K_cat = 2L)
  ng_6 <- expand.grid(item = paste0("I", 1:4), rater = paste0("RN", 1:6), stringsAsFactors = FALSE)
  ng_6$person <- "P_NEW"
  ng_6$score <- 1L
  seen_g <- character()
  entries_at_group_start <- integer()
  orig_eval_plan <- BayesMFRM:::eval_bmfrm_chunk_from_plan
  testthat::with_mocked_bindings(
    eval_bmfrm_chunk_from_plan = function(plan, rows, ..., offset_cache = plan$offset_cache) {
      gid <- unique(as.character(ng_6$rater[rows]))[1L]
      if (!gid %in% seen_g) {
        seen_g <<- c(seen_g, gid)
        entries_at_group_start <<- c(entries_at_group_start, if (is.null(offset_cache)) 0L else length(offset_cache$data))
      }
      orig_eval_plan(plan, rows, ..., offset_cache = offset_cache)
    },
    invisible(BayesMFRM:::score_general_unseen_fold(small_4, ng_6, "rater", M = 100L, seed = 551L, chunk_size = 1L)),
    .package = "BayesMFRM"
  )
  expect_identical(entries_at_group_start, rep(0L, 6L))

  # Strict byte-budget cap check + Student-t / horseshoe chunk invariance
  tiny_cache <- BayesMFRM:::new_interaction_offset_cache(D = 32L, eff_chunk_size = 1L, max_bytes = 1024)
  for (m_i in 1:20) {
    invisible(BayesMFRM:::eval_constrained_interaction_chunk(
      D = 32L, n_known = 10L, b_info = list(mode = "gaussian_basis", scale = "raw", sd = rep(0.5, 32L)),
      slice_seed = m_i, chunk_k_idx = 1L, scope_cols = 1:10,
      need_after_chunk = TRUE, cache_key = paste0("k:", m_i), offset_cache = tiny_cache
    ))
    expect_lte(tiny_cache$bytes, 1024)
  }
  expect_equal(length(tiny_cache$data), 4L) # 4 * (32 * 8) = 1024 bytes max, pinned without wipe-all!

  for (sh_mode in c("student_t", "horseshoe")) {
    obj_sh <- make_bias_grp_fit(6L, D_draws = 32L, K_cat = 2L, shrinkage = sh_mode)
    nd_sh <- data.frame(person = "P_NEW", item = paste0("I", 1:6), rater = "R_NEW", score = rep(1:2, 3L), stringsAsFactors = FALSE)
    s_sh_1 <- BayesMFRM:::score_general_unseen_fold(obj_sh, nd_sh, "rater", M = 10L, seed = 44L, chunk_size = 1L)
    s_sh_6 <- BayesMFRM:::score_general_unseen_fold(obj_sh, nd_sh, "rater", M = 10L, seed = 44L, chunk_size = 6L)
    expect_equal(s_sh_1, s_sh_6, tolerance = 1e-12)
  }

  # 5. Live cumulative memory & Rprofmem check on score_grouped_fold WITH interaction
  #    (no u_draws retention, peak offset cache = D * M * 8 = 32 KB << 1.14 MB, and 0 bytes left at group end)
  obj_m30 <- make_bias_grp_fit(30L, D_draws = 32L, K_cat = 2L)
  nm_500 <- data.frame(person = "P1", item = rep(paste0("I", 1:30), length.out = 500L), rater = "R_NEW", score = 1L, stringsAsFactors = FALSE)
  orig_lme <- BayesMFRM:::logmeanexp_vec
  peak_cache_bytes <- 0
  end_cache_bytes <- -1
  has_u_draws <- FALSE
  testthat::with_mocked_bindings(
    logmeanexp_vec = function(x) {
      fr <- parent.frame()
      if (exists("u_draws", envir = fr, inherits = FALSE)) has_u_draws <<- TRUE
      if (exists("g_cache", envir = fr, inherits = FALSE) && !is.null(fr$g_cache)) {
        peak_cache_bytes <<- max(peak_cache_bytes, fr$g_cache$bytes)
        end_cache_bytes <<- fr$g_cache$bytes
      }
      orig_lme(x)
    },
    invisible(BayesMFRM:::score_grouped_fold(obj_m30, nm_500, "rater", M = 128L, seed = 55L, chunk_size = 1L)),
    .package = "BayesMFRM"
  )
  expect_false(has_u_draws)
  expect_equal(peak_cache_bytes, 128 * 32 * 8)
  expect_lt(peak_cache_bytes, 32 * 500 * 8)
  expect_equal(end_cache_bytes, 0)

  # 6. Practical MCMC draw count (D = 4000, M = 100) in both general & ordinary CV:
  #    cache budget accounts for M * n_concurrent_offsets * D * 8 (~3.2 MB <= 64 MiB),
  #    so chunk_size = 1 computes centering exactly M = 100 times (not 712 / 1,392)!
  for (J_prac in c(10L, 20L)) {
    obj_4k <- make_bias_grp_fit(J_prac, D_draws = 4000L, K_cat = 2L)
    nd_gen_4k <- data.frame(
      person = "P_NEW",
      item = paste0("I", seq_len(J_prac)),
      rater = "R_NEW",
      score = rep(1:2, length.out = J_prac),
      stringsAsFactors = FALSE
    )
    gen_c1 <- measure_work(function() BayesMFRM:::score_general_unseen_fold(
      obj_4k, nd_gen_4k, "rater", M = 100L, seed = 901L, chunk_size = 1L
    ))
    gen_all <- measure_work(function() BayesMFRM:::score_general_unseen_fold(
      obj_4k, nd_gen_4k, "rater", M = 100L, seed = 901L, chunk_size = J_prac
    ))
    expect_equal(gen_c1$center_calls, 100L)
    expect_equal(gen_all$center_calls, 100L)
    expect_equal(gen_c1$col_draws, 100L * (2L * J_prac - 1L))
    expect_equal(gen_all$col_draws, 100L * J_prac)
    expect_equal(gen_c1$ans, gen_all$ans, tolerance = 1e-12)
  }

  obj_4k_20 <- make_bias_grp_fit(20L, D_draws = 4000L, K_cat = 2L)
  nd_ord_4k <- data.frame(
    person = "P1",
    item = paste0("I", 1:20),
    rater = "R_NEW",
    score = rep(1:2, 10L),
    stringsAsFactors = FALSE
  )
  ord_c1 <- measure_work(function() BayesMFRM:::score_grouped_fold(
    obj_4k_20, nd_ord_4k, "rater", M = 100L, seed = 902L, chunk_size = 1L
  ))
  ord_all <- measure_work(function() BayesMFRM:::score_grouped_fold(
    obj_4k_20, nd_ord_4k, "rater", M = 100L, seed = 902L, chunk_size = 20L
  ))
  expect_equal(ord_c1$center_calls, 100L)
  expect_equal(ord_all$center_calls, 100L)
  expect_equal(ord_c1$col_draws, 100L * (2L * 20L - 1L))
  expect_equal(ord_all$col_draws, 100L * 20L)
  expect_equal(ord_c1$ans, ord_all$ans, tolerance = 1e-12)

  # Verify hard 64 MiB cap when estimated concurrent offsets exceed hard_max_bytes
  capped_cache <- BayesMFRM:::new_interaction_offset_cache(
    D = 4000L, eff_chunk_size = 1L, n_concurrent_offsets = 3000L
  )
  expect_equal(capped_cache$max_bytes, 64 * 1024^2)
})

test_that("deep-audit round regressions: facet column collisions, CV overrides, custom_stan guards, interaction keywords, bundled fair_scores, compare_bmfrm cross-columns, NSE shadowing, and orphan binaries", {
  # 1. Facets named 'x' or 'weight' never clobber observed ratings, weights, or identifiability
  d_fx <- expand.grid(person = c("P1", "P2", "P3"), x = c("X1", "X2"), rater = c("R1", "R2"), stringsAsFactors = FALSE)
  d_fx$score <- rep(1:4, length.out = nrow(d_fx))
  sp_fx_min <- resolve_bmfrm_spec(parse_bmfrm_formula(score ~ person + x + rater, d_fx, K = 4L, save_gq = "minimal"))
  sp_fx_min$eval_contract <- resolve_eval_contract(sp_fx_min, stan_code = NULL)
  di_fx <- prepare_data_bmfrm(sp_fx_min, d_fx)
  dc_fx <- check_design_bmfrm(score ~ person + x + rater, di_fx$data_clean, K = 4L)
  D_fx <- 8L
  draws_fx <- cbind(
    matrix(rep(c(-0.5, 0, 0.5), each = D_fx), D_fx, 3L, dimnames = list(NULL, paste0("theta[", 1:3, "]"))),
    matrix(rep(c(-0.3, 0.3), each = D_fx), D_fx, 2L, dimnames = list(NULL, paste0("x[", 1:2, "]"))),
    matrix(rep(c(-0.2, 0.2), each = D_fx), D_fx, 2L, dimnames = list(NULL, paste0("rater[", 1:2, "]"))),
    matrix(rep(c(-0.6, 0, 0.6), each = D_fx), D_fx, 3L, dimnames = list(NULL, paste0("tau[", 1:3, "]")))
  )
  mock_fx_min <- list(
    draws = function(variables = NULL, format = "matrix") {
      if (is.null(variables)) return(draws_fx)
      sel <- unlist(lapply(variables, function(v) colnames(draws_fx)[colnames(draws_fx) == v | startsWith(colnames(draws_fx), paste0(v, "["))]))
      draws_fx[, sel, drop = FALSE]
    },
    metadata = function() list(stan_variables = c("theta", "x", "rater", "tau"), iter_sampling = 4L),
    num_chains = function() 2L
  )
  fit_fx_min <- structure(list(fit = mock_fx_min, spec = sp_fx_min, data_info = di_fx, design = dc_fx), class = "bmfrm_fit")
  res_fx_min <- residuals(fit_fx_min)
  expect_false(anyDuplicated(names(res_fx_min)) > 0L)
  expect_identical(as.integer(res_fx_min$x), as.integer(d_fx$score))
  expect_identical(as.character(res_fx_min$facet_x), as.character(d_fx$x))
  expect_equal(nrow(facet_infit_outfit(res_fx_min, "x")), 2L)
  expect_setequal(facet_infit_outfit(res_fx_min, "x")$facet_value, c("X1", "X2"))

  # Also verify save_gq = "full" path with mu/sigma2 present
  q_fx <- BayesMFRM:::compute_bmfrm_quantities(fit_fx_min)
  mu_fx <- q_fx$mu; colnames(mu_fx) <- paste0("mu[", seq_len(ncol(mu_fx)), "]")
  s2_fx <- q_fx$sigma2; colnames(s2_fx) <- paste0("sigma2[", seq_len(ncol(s2_fx)), "]")
  draws_fx_full <- cbind(draws_fx, mu_fx, s2_fx)
  mock_fx_full <- list(
    draws = function(variables = NULL, format = "matrix") {
      sel <- unlist(lapply(variables, function(v) colnames(draws_fx_full)[colnames(draws_fx_full) == v | startsWith(colnames(draws_fx_full), paste0(v, "["))]))
      draws_fx_full[, sel, drop = FALSE]
    },
    metadata = function() list(stan_variables = c("theta", "x", "rater", "tau", "mu", "sigma2")),
    summary = function(v) tibble::tibble(
      variable = paste0(v, "[", seq_len(ncol(mu_fx)), "]"),
      mean = if (v == "mu") colMeans(mu_fx) else colMeans(s2_fx)
    )
  )
  fit_fx_full <- fit_fx_min; fit_fx_full$fit <- mock_fx_full
  res_fx_full <- residuals(fit_fx_full)
  expect_false(anyDuplicated(names(res_fx_full)) > 0L)
  expect_identical(as.integer(res_fx_full$x), as.integer(d_fx$score))
  expect_equal(res_fx_full$resid, res_fx_min$resid, tolerance = 1e-12)
  rc_fx <- residual_cor_criteria(res_fx_full, item_facet = "x", unit_facets = c("person", "rater"))
  expect_setequal(rownames(rc_fx$cor_matrix), c("X1", "X2"))

  # Facet named 'weight' under fair_scores(standardize = TRUE, summary = FALSE) with explicit reference_weights
  d_fw <- expand.grid(person = c("P1", "P2", "P3"), weight = c("heavy", "light"), rater = c("R1", "R2"), stringsAsFactors = FALSE)
  d_fw$score <- rep(1:4, length.out = nrow(d_fw))
  sp_fw <- resolve_bmfrm_spec(parse_bmfrm_formula(score ~ person + weight + rater, d_fw, K = 4L, save_gq = "minimal"))
  sp_fw$eval_contract <- resolve_eval_contract(sp_fw, stan_code = NULL)
  di_fw <- prepare_data_bmfrm(sp_fw, d_fw)
  dc_fw <- check_design_bmfrm(score ~ person + weight + rater, di_fw$data_clean, K = 4L)
  draws_fw <- draws_fx
  colnames(draws_fw)[colnames(draws_fw) %in% c("x[1]", "x[2]")] <- c("weight[1]", "weight[2]")
  mock_fw <- list(
    draws = function(variables = NULL, format = "matrix") {
      sel <- unlist(lapply(variables, function(v) colnames(draws_fw)[colnames(draws_fw) == v | startsWith(colnames(draws_fw), paste0(v, "["))]))
      draws_fw[, sel, drop = FALSE]
    },
    metadata = function() list(stan_variables = c("theta", "weight", "rater", "tau"), iter_sampling = 4L),
    num_chains = function() 2L
  )
  fit_fw <- structure(list(fit = mock_fw, spec = sp_fw, data_info = di_fw, design = dc_fw), class = "bmfrm_fit")
  fs_fw <- fair_scores(fit_fw, standardize = TRUE, summary = FALSE,
                       reference_weights = list(weight = c(heavy = 0.7, light = 0.3)))
  expect_setequal(unique(as.character(fs_fw$weight)), c("heavy", "light"))
  expect_equal(fs_fw$reference_weight[1:2], c(0.7, 0.3), tolerance = 1e-12)
  expect_true(all(fs_fw$likelihood_identified))

  # 2. cv_bmfrm() pre-validates overrides, adjusts default warmup when iter is overridden, and rejects locked/reserved args
  fit_fx_min$fit_args <- list(iter = 2000L, warmup = 1000L, chains = 2L, cores = 1L, cache_dir = tempdir())
  seen_fit_args <- list()
  invisible(suppressWarnings(cv_bmfrm(
    fit_fx_min, target = "rating", K_folds = 2L, seed = 11L, iter = 800L,
    refit_fn = function(train_df, fold, spec, fit_args) {
      seen_fit_args[[fold]] <<- fit_args
      fit_fx_min
    }
  )))
  expect_identical(seen_fit_args[[1L]]$iter, 800L)
  expect_identical(seen_fit_args[[1L]]$warmup, 400L)
  expect_error(cv_bmfrm(fit_fx_min, target = "rating", K_folds = 2L, iter = 500L, warmup = 500L), "warmup must be a non-negative integer strictly less than iter")
  expect_error(cv_bmfrm(fit_fx_min, target = "rating", K_folds = 2L, iter_sampling = 100L), "reserved cmdstanr argument")
  expect_error(cv_bmfrm(fit_fx_min, target = "rating", K_folds = 2L, K = 5L), "locked model/data specification")
  expect_error(cv_bmfrm(fit_fx_min, target = "rating", K_folds = 2L, cache_dir = ""), "cache_dir must be a non-empty")

  # 3. build_stan_code() rejects custom_stan specs instead of emitting improper zero-prior code
  sp_cust <- parse_bmfrm_formula(score ~ person + rater, d_fx[, c("person", "rater", "score")], K = 4L)
  sp_cust$custom_stan <- TRUE
  expect_error(build_stan_code(sp_cust), "custom_stan")

  # 4. Concatenated interaction names colliding with Stan keywords/locals are rejected pre-compile
  for (bad_pair in list(c("e", "x"), c("e", "ta"), c("fo", "r"), c("m", "odel"))) {
    d_bi <- data.frame(person = paste0("P", 1:6), score = rep(1:3, 2L), stringsAsFactors = FALSE)
    d_bi[[bad_pair[1L]]] <- rep(c("A", "B"), 3L)
    d_bi[[bad_pair[2L]]] <- rep(c("C", "D"), each = 3L)
    fml_bi <- stats::reformulate(c("person", bad_pair, paste(bad_pair, collapse = ":")), "score")
    expect_error(parse_bmfrm_formula(fml_bi, d_bi, K = 3L), "collide or are invalid Stan identifiers")
    expect_error(check_design_bmfrm(fml_bi, d_bi, K = 3L), "collide or are invalid Stan identifiers")
  }

  # 5. Bundled 3-facet fair_scores() uses consistent default adjust_facets across standardize = FALSE and TRUE
  m3_b <- bmfrm_model("three-facet/person-rater-scenario-main")
  d3_b <- expand.grid(person = c("P1", "P2"), rater = c("R1", "R2"), scenario = c("S1", "S2"), stringsAsFactors = FALSE)
  d3_b$score <- rep(1:4, 2L)
  sp3_b <- parse_bmfrm_formula(score ~ person + rater + scenario, d3_b, K = 4L)
  sp3_b$custom_stan <- TRUE
  di3_b <- suppressWarnings(prepare_data_bmfrm(sp3_b, d3_b))
  dc3_b <- check_design_bmfrm(score ~ person + rater + scenario, di3_b$data_clean, K = 4L)
  draws3_b <- list(
    theta = matrix(rep(c(-0.4, 0.5), each = 4L), 4L, 2L, dimnames = list(NULL, c("theta[1]", "theta[2]"))),
    rater = matrix(rep(c(-0.3, 0.3), each = 4L), 4L, 2L, dimnames = list(NULL, c("rater[1]", "rater[2]"))),
    scenario = matrix(rep(c(-0.7, 0.7), each = 4L), 4L, 2L, dimnames = list(NULL, c("scenario[1]", "scenario[2]"))),
    tau = matrix(rep(c(-0.6, 0, 0.6), each = 4L), 4L, 3L, dimnames = list(NULL, paste0("tau[", 1:3, "]")))
  )
  eta_f3 <- c(-0.4, 0.5)[di3_b$stan_data$person_id]
  mu_f3 <- as.vector(BayesMFRM:::rsm_probabilities(eta_f3, c(-0.6, 0, 0.6)) %*% 1:4)
  draws3_b$mu_fair <- matrix(rep(mu_f3, each = 4L), 4L, nrow(d3_b), dimnames = list(NULL, paste0("mu_fair[", seq_len(nrow(d3_b)), "]")))
  mock3_b <- list(
    metadata = function() list(stan_variables = names(draws3_b), iter_sampling = 2L),
    num_chains = function() 2L,
    draws = function(var = NULL, format = "matrix") {
      if (format == "df") {
        out <- as.data.frame(draws3_b$theta)
        out$.chain <- c(1L, 1L, 2L, 2L); out$.iteration <- c(1L, 2L, 1L, 2L); out$.draw <- 1:4
        return(out)
      }
      draws3_b[[var]]
    }
  )
  obj3_b <- structure(list(fit = mock3_b, spec = sp3_b, data_info = di3_b, design = dc3_b, stan_code = m3_b$code), class = "bmfrm_fit")
  expect_message(fs3_unstd <- fair_scores(obj3_b, standardize = FALSE), "Using custom Stan built-in mu_fair")
  expect_message(fs3_std <- fair_scores(obj3_b, standardize = TRUE), "Using custom Stan built-in mu_fair")
  expect_setequal(attr(fs3_unstd, "reference_meta")$adjust_facets, c("rater", "scenario"))
  expect_setequal(attr(fs3_std, "reference_meta")$adjust_facets, c("rater", "scenario"))
  expect_equal(fs3_unstd$mean_fair_score, fs3_std$mean_fair_score, tolerance = 1e-12)

  # 6 & 10. compare_bmfrm() detects row swaps on facet-vs-extra columns and deduplicates explicit names
  d_cmp1 <- data.frame(
    person = c("P1", "P1", "P2", "P2"),
    item = c("I1", "I2", "I1", "I2"),
    rater = c("R1", "R1", "R2", "R2"),
    score = c(2L, 2L, 3L, 3L),
    stringsAsFactors = FALSE
  )
  d_cmp2 <- d_cmp1[c(2L, 1L, 3L, 4L), ]
  rownames(d_cmp2) <- NULL
  sp_c2f <- resolve_bmfrm_spec(parse_bmfrm_formula(score ~ person + rater, d_cmp1, K = 4L, save_gq = "minimal"))
  sp_c2f$eval_contract <- resolve_eval_contract(sp_c2f, NULL)
  di_c2f <- suppressWarnings(prepare_data_bmfrm(sp_c2f, d_cmp1))
  sp_c3f <- resolve_bmfrm_spec(parse_bmfrm_formula(score ~ person + item + rater, d_cmp2, K = 4L, save_gq = "minimal"))
  sp_c3f$eval_contract <- resolve_eval_contract(sp_c3f, NULL)
  di_c3f <- suppressWarnings(prepare_data_bmfrm(sp_c3f, d_cmp2))
  mk_cmp_fit <- function(sp, di, extra_item = FALSE) {
    dm <- cbind(
      matrix(rep(c(-0.3, 0.3), each = 8L), 8L, 2L, dimnames = list(NULL, c("theta[1]", "theta[2]"))),
      matrix(rep(c(-0.2, 0.2), each = 8L), 8L, 2L, dimnames = list(NULL, c("rater[1]", "rater[2]"))),
      matrix(rep(c(-0.5, 0, 0.5), each = 8L), 8L, 3L, dimnames = list(NULL, paste0("tau[", 1:3, "]")))
    )
    if (extra_item) dm <- cbind(dm, matrix(rep(c(-0.4, 0.4), each = 8L), 8L, 2L, dimnames = list(NULL, c("item[1]", "item[2]"))))
    mf <- list(
      draws = function(variables = NULL, format = "matrix") {
        if (format == "df") {
          out <- as.data.frame(dm[, 1:2]); out$.chain <- rep(1:2, each = 4L); out$.iteration <- rep(1:4, 2L); out$.draw <- 1:8
          return(out)
        }
        sel <- unlist(lapply(variables, function(v) colnames(dm)[colnames(dm) == v | startsWith(colnames(dm), paste0(v, "["))]))
        dm[, sel, drop = FALSE]
      },
      metadata = function() list(stan_variables = unique(sub("\\[.*$", "", colnames(dm))), iter_sampling = 4L),
      num_chains = function() 2L
    )
    structure(list(fit = mf, spec = sp, data_info = di, convergence = list(converged = TRUE, status = "ok")), class = "bmfrm_fit")
  }
  loo_2f <- suppressWarnings(loo::loo(mk_cmp_fit(sp_c2f, di_c2f, FALSE), use_stan_log_lik = FALSE))
  loo_3f_swap <- suppressWarnings(loo::loo(mk_cmp_fit(sp_c3f, di_c3f, TRUE), use_stan_log_lik = FALSE))
  expect_error(compare_bmfrm(two = loo_2f, three = loo_3f_swap), "observation alignment mismatch on shared column\\(s\\) \\[item\\]")
  cmp_dup <- compare_bmfrm(m = loo_2f, m = loo_2f)
  expect_false(anyDuplicated(cmp_dup$table$model) > 0L)

  # 7, 8, 11, 12. NSE caller variable 'x', custom_stan summary() with legacy bias, stripped facet_cols, and orphan binary cleanup
  x <- "rater"
  expect_setequal(facet_infit_outfit(res_fx_min, x)$facet_value, c("R1", "R2"))
  res_stripped <- as.data.frame(res_fx_min)
  attr(res_stripped, "facet_cols") <- NULL
  expect_setequal(unique(facet_infit_outfit_all(res_stripped)$facet), c("person", "facet_x", "rater"))

  sp_cb <- sp_cust; sp_cb$facets_bias <- "person:rater"
  fit_cb <- fit_fx_min; fit_cb$spec <- sp_cb; fit_cb$convergence <- list(converged = TRUE, status = "ok", max_rhat = 1.0, min_ess_bulk = 500, min_ess_tail = 500, n_divergent = 0L, n_max_treedepth = 0L, min_ebfmi = 0.9)
  fit_cb$fit$summary <- function(v, ...) stop("parameter not found")
  summ_cb <- paste(utils::capture.output(summary(fit_cb)), collapse = "\n")
  expect_match(summ_cb, "not found in custom Stan model")

  tmp_oc <- tempfile("orphan_cache_")
  dir.create(tmp_oc)
  on.exit(unlink(tmp_oc, recursive = TRUE), add = TRUE)
  orphan_bin <- file.path(tmp_oc, paste0("mdl_", paste(rep("a", 32L), collapse = "")))
  user_note <- file.path(tmp_oc, paste0("note_", paste(rep("b", 32L), collapse = "")))
  writeBin(as.raw(c(0xcf, 0xfa, 0xed, 0xfe, 0x00, 0x00)), orphan_bin)
  writeLines("plain text user note", user_note)
  Sys.setFileTime(orphan_bin, Sys.time() - 40 * 86400)
  Sys.setFileTime(user_note, Sys.time() - 40 * 86400)
  deleted_oc <- suppressMessages(clean_model_cache(tmp_oc, max_age_days = 30))
  expect_true(orphan_bin %in% deleted_oc)
  expect_false(file.exists(orphan_bin))
  expect_true(file.exists(user_note))
})
