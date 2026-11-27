test_that("Sinharay (2006) posterior predictive p-values remain calibrated on sparse designs and detect noisy raters", {
  # Sparse design: 40 examinees x 2 ratings each (n_obs = 2 per person!), 4 items, 4 raters
  set.seed(2026)
  P <- 40L; I <- 4L; R <- 4L; K <- 4L
  design <- data.frame(
    person = paste0("P", rep(seq_len(P), each = 2L)),
    item   = paste0("I", ((seq_len(2L * P) - 1L) %% I) + 1L),
    rater  = paste0("R", ((rep(seq_len(P), each = 2L) + rep(0:1, P) - 1L) %% R) + 1L),
    stringsAsFactors = FALSE
  )
  theta_true <- seq(-1.5, 1.5, length.out = P)
  item_true  <- c(-0.6, -0.2, 0.2, 0.6)
  rater_true <- c(-0.4, -0.1, 0.1, 0.4)
  tau_true   <- c(-0.9, 0.0, 0.9)

  sim <- simulate_bmfrm(
    score ~ person + item + rater, design, K = K,
    parameters = list(theta = theta_true, item = item_true, rater = rater_true, tau = tau_true),
    seed = 101
  )
  dat <- sim$data
  spec <- parse_bmfrm_formula(score ~ person + item + rater, dat, K = K, save_gq = "minimal")
  di <- prepare_data_bmfrm(spec, dat)

  # Emulate posterior draws with realistic small-n examinee posterior uncertainty (SD = 0.65)
  D <- 120L
  theta_draws <- sweep(matrix(stats::rnorm(D * P, 0, 0.65), D, P), 2L, theta_true, "+")
  colnames(theta_draws) <- paste0("theta[", seq_len(P), "]")
  item_draws <- matrix(rep(item_true, each = D), D, I, dimnames = list(NULL, paste0("item[", seq_len(I), "]")))
  rater_draws <- matrix(rep(rater_true, each = D), D, R, dimnames = list(NULL, paste0("rater[", seq_len(R), "]")))
  tau_draws <- matrix(rep(tau_true, each = D), D, K - 1L, dimnames = list(NULL, paste0("tau[", seq_len(K - 1L), "]")))

  draws_store <- list(theta = theta_draws, item = item_draws, rater = rater_draws, tau = tau_draws)
  mock_fit <- list(
    metadata = function() list(stan_variables = names(draws_store), iter_sampling = D / 2L),
    num_chains = function() 2L,
    draws = function(var, format = "matrix") draws_store[[var]]
  )
  obj <- structure(list(fit = mock_fit, spec = spec, data_info = di), class = "bmfrm_fit")

  # Realized MSQ is systematically > 1 due to theta posterior variance on n_obs = 2,
  # whereas replicated MSQ matches it so infit_ppp centers near 0.5 without mass false alarms!
  ff_p <- facet_fit(obj, "person", seed = 42)
  expect_gt(mean(ff_p$infit_msq), 1.05)
  expect_equal(mean(ff_p$infit_msq), mean(ff_p$infit_rep_mean), tolerance = 0.15)
  expect_gt(mean(ff_p$infit_ppp), 0.30)
  expect_lt(mean(ff_p$infit_ppp), 0.70)
  expect_false("infit_zstd" %in% names(ff_p))

  # Inject erratic extreme noise into Rater R1 only: infit_ppp drops near 0
  dat_noisy <- dat
  r1_rows <- which(dat_noisy$rater == "R1")
  dat_noisy$score[r1_rows] <- ifelse(seq_along(r1_rows) %% 2L == 1L, 1L, K)
  di_noisy <- prepare_data_bmfrm(spec, dat_noisy)
  obj_noisy <- structure(list(fit = mock_fit, spec = spec, data_info = di_noisy), class = "bmfrm_fit")
  ff_r <- facet_fit(obj_noisy, "rater", seed = 42)
  expect_lt(ff_r$infit_ppp[ff_r$label == "R1"], 0.05)
})

test_that("custom role names (person_facet, adjust_facets, class = 'facet') work end-to-end without 'person'/'item'/'rater'", {
  d <- expand.grid(
    student = paste0("S", 1:6),
    task    = paste0("T", 1:3),
    judge   = paste0("J", 1:3),
    stringsAsFactors = FALSE
  )
  ps <- c(
    prior("normal(0, 1.5)", "mean", facet = "student", person_facet = "student"),
    prior("lognormal(0, 0.4)", "sd", facet = "student", person_facet = "student"),
    prior("normal(0, 0.8)", "facet", facet = "task", scale = "marginal", person_facet = "student"),
    prior("normal(0, 0.6)", "sd", facet = "judge", scale = "marginal", person_facet = "student"),
    prior("normal(0, 1.5)", "tau", scale = "marginal", person_facet = "student")
  )
  fml <- rating ~ student + task + judge
  sim <- simulate_bmfrm(
    fml, d, K = 4, priors = ps,
    hierarchical = c("student", "judge"),
    person_facet = "student", seed = 303
  )
  expect_identical(sim$spec$person_facet, "student")
  expect_setequal(sim$hierarchical, c("student", "judge"))

  spec <- parse_bmfrm_formula(
    fml, sim$data, K = 4,
    person_facet = "student", adjust_facets = "judge",
    hierarchical = c("student", "judge")
  )
  dc <- check_design_bmfrm(fml, sim$data, K = 4, person_facet = "student")
  expect_true(dc$identifiable)

  code <- build_stan_code(spec, ps, save_gq = "full")
  expect_match(code, "theta[student_id]", fixed = TRUE)
  expect_false(grepl("person_id", code, fixed = TRUE))
  expect_match(code, "eta_fair = theta[student_id] - task[task_id];", fixed = TRUE)

  di <- prepare_data_bmfrm(spec, sim$data, ps)
  D <- 16L
  p1 <- sim$parameters[[1L]]
  mk <- function(nm, v) matrix(rep(v, each = D), D, length(v), dimnames = list(NULL, paste0(nm, "[", seq_along(v), "]")))
  store <- list(
    theta = mk("theta", p1$theta),
    task  = mk("task", p1$task),
    judge = mk("judge", p1$judge),
    tau   = mk("tau", p1$tau),
    mu_theta = matrix(p1$mu_theta, D, 1, dimnames = list(NULL, "mu_theta")),
    sigma_theta = matrix(p1$sigma_theta, D, 1, dimnames = list(NULL, "sigma_theta")),
    sigma_judge = matrix(p1$sigma_judge, D, 1, dimnames = list(NULL, "sigma_judge"))
  )
  mock_fit <- list(
    metadata = function() list(stan_variables = names(store), iter_sampling = D / 2L, iter_warmup = 10L),
    num_chains = function() 2L,
    draws = function(var = NULL, format = "matrix") {
      if (format == "df") {
        df <- as.data.frame(store$theta)
        df$.chain <- rep(1:2, each = D / 2L)
        df$.iteration <- rep(seq_len(D / 2L), 2L)
        df$.draw <- seq_len(D)
        return(df)
      }
      store[[var]]
    },
    summary = function(variables, ...) posterior::summarise_draws(posterior::as_draws_matrix(store[[variables]]), ...)
  )
  obj <- structure(
    list(fit = mock_fit, spec = resolve_bmfrm_spec(spec, ps), data_info = di,
         priors = ps, design = dc, hierarchical = c("student", "judge"),
         convergence = list(converged = TRUE, status = "ok", max_rhat = 1.0, min_ess_bulk = 400)),
    class = "bmfrm_fit"
  )

  fs <- fair_scores(obj, standardize = TRUE, decompose = TRUE)
  expect_true("student" %in% names(fs))
  expect_equal(attr(fs, "reference_meta")$adjust_facets, "judge")

  fsum <- facet_summary(obj, "student")
  expect_equal(fsum$level, paste0("S", 1:6))

  # Unseen student prediction works via person_facet = "student"
  new_df <- data.frame(student = "S_new", task = "T1", judge = "J1", stringsAsFactors = FALSE)
  pr <- predict(obj, newdata = new_df, allow_new_levels = TRUE, sample_new_levels = "hierarchical", seed = 11)
  expect_equal(nrow(pr), 1L)
  expect_true(is.finite(pr$estimate))
})

test_that("summarise_bias excludes zero-observation cells from BFDR/flags and honors ci_level", {
  # 3x3 item x rater grid where cell (I3, R3) has 0 observations
  grid_full <- expand.grid(person = paste0("P", 1:6), item = paste0("I", 1:3), rater = paste0("R", 1:3), stringsAsFactors = FALSE)
  d <- subset(grid_full, !(item == "I3" & rater == "R3"))
  d$score <- ((seq_len(nrow(d)) - 1L) %% 4L) + 1L

  spec <- parse_bmfrm_formula(score ~ person + item + rater + item:rater, d, K = 4)
  dc <- check_design_bmfrm(spec$formula, d, K = 4)

  # Construct zero-sum 3x3 interaction draws where unobserved (I3, R3) is algebraically non-zero (+0.4)
  B <- matrix(c(0.4, -0.2, -0.2,
                -0.2, 0.4, -0.2,
                -0.2, -0.2, 0.4), nrow = 3L, byrow = TRUE)
  D <- 200L
  set.seed(55)
  mat <- do.call(rbind, lapply(seq_len(D), function(d_i) as.vector(B + stats::rnorm(1, 0, 0.05) * B)))
  colnames(mat) <- paste0("itemrater[", rep(1:3, 3), ",", rep(1:3, each = 3), "]")

  di <- prepare_data_bmfrm(spec, d)
  # Case A: design = NULL (e.g. design_check = "none"): unobserved cell must still get "unobserved_cell" and NA bfdr!
  obj_nodc <- structure(
    list(spec = spec, design = NULL, data_info = di,
         convergence = list(converged = TRUE, status = "ok"),
         fit = list(draws = function(...) mat)),
    class = "bmfrm_fit"
  )
  sb <- summarise_bias(obj_nodc, "item:rater", ci_level = 0.80)
  cell_33 <- sb[sb$item == "I3" & sb$rater == "R3", ]
  expect_equal(cell_33$n_observations, 0L)
  expect_equal(cell_33$flag_bias, "unobserved_cell")
  expect_true(is.na(cell_33$bfdr))
  expect_equal(cell_33$decision_status, "unobserved_cell")
  expect_true(all(sb$flag_bias[sb$n_observations > 0] == "identification_unchecked"))
  expect_true(all(is.na(sb$bfdr)))

  # Case B: with design check attached, rank-deficient cells get "not_likelihood_identified" and NA bfdr
  obj_dc <- obj_nodc; obj_dc$design <- dc
  sb_dc <- summarise_bias(obj_dc, "item:rater", ci_level = 0.80)
  expect_true(is.na(sb_dc$bfdr[sb_dc$item == "I3" & sb_dc$rater == "R3"]))

  # ci_level = 0.80 sets ci_lower/ci_upper to 10th and 90th percentiles while preserving lower_95/upper_95
  q80 <- stats::quantile(mat[, "itemrater[1,1]"], probs = c(0.10, 0.90))
  expect_equal(sb$ci_lower[1], unname(q80[1]), tolerance = 1e-12)
  expect_equal(sb$ci_upper[1], unname(q80[2]), tolerance = 1e-12)
  expect_lt(sb$lower_95[1], sb$ci_lower[1])
})

test_that("prior() blocks Stan code injection, configures horseshoe slab, and stores portable csv_rel cache metadata", {
  # L10: AST guard rejects semicolon / multi-statement Stan injection
  expect_error(prior("normal(0, 1); target += 100", "item"), "semicolons|Invalid")
  expect_error(prior("normal(0, 1) + 1", "item"), "Invalid|Unsupported")
  expect_error(prior("unknown_rng(0, 1)", "item"), "Invalid|unsupported")

  # L6: Custom horseshoe slab parameters propagate into Stan data
  d <- expand.grid(person = paste0("P", 1:4), item = paste0("I", 1:2), rater = paste0("R", 1:2))
  d$score <- 1L
  ps_hs <- c(
    prior("normal(0, 2)", "mean", facet = "person"),
    prior("lognormal(0, 0.5)", "sd", facet = "person"),
    prior("normal(0, 1)", "facet", facet = "item"),
    prior("normal(0, 1)", "facet", facet = "rater"),
    prior("normal(0, 2)", "tau"),
    prior("normal(0, 0.5)", "sd", facet = "item:rater", local_shrinkage = "horseshoe",
          slab_df = 6, slab_scale = 1.5)
  )
  sp_hs <- parse_bmfrm_formula(score ~ person + item + rater + item:rater, d, K = 3,
                               hierarchical = c("person", "item:rater"))
  sd_hs <- prepare_data_bmfrm(sp_hs, d, ps_hs)$stan_data
  expect_equal(sd_hs$hs_slab_shape_itemrater, 3.0)
  expect_equal(sd_hs$hs_slab_rate_itemrater, 0.5 * 6 * 1.5^2)

  # G6: Saved cache RDS stores portable relative csv_rel path that resolves after directory move
  dir_a <- tempfile("cache_orig_")
  dir_b <- tempfile("cache_moved_")
  dir.create(dir_a, recursive = TRUE)
  on.exit(unlink(c(dir_a, dir_b), recursive = TRUE), add = TRUE)
  model_mock <- list(
    code = function() "model {}",
    sample = function(output_dir, ...) {
      p <- file.path(output_dir, "chain1.csv")
      writeLines("# dummy", p)
      list(id = 1L, output_files = function() p)
    }
  )
  rds_a <- file.path(dir_a, "fit.rds")
  fit_cmdstan_cached(model_mock, list(N = 1L), rds_a, refit = "on_change", seed = 11L)
  payload_a <- readRDS(rds_a)
  expect_equal(payload_a$meta$csv_rel, file.path("csv", "chain1.csv"))
  expect_true(file.rename(dir_a, dir_b))
  rds_b <- file.path(dir_b, "fit.rds")
  payload_b <- readRDS(rds_b)
  expect_false(all(file.exists(payload_b$fit$output_files())))
  expect_true(all(file.exists(file.path(dirname(rds_b), payload_b$meta$csv_rel))))
})
