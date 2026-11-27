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


test_that("ignored fitting options fail before compilation and custom priors remain unknown", {
  calls <- 0L
  testthat::local_mocked_bindings(cmdstan_version = function() numeric_version("2.39.0"),
                                 .package = "cmdstanr")
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


