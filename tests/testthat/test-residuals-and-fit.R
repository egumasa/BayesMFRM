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


