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
  testthat::local_mocked_bindings(cmdstan_version = function() numeric_version("2.39.0"),
                                 .package = "cmdstanr")
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


