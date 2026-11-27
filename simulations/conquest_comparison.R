# Reproducible Cross-Software Validation & Sparse-Design Benchmark:
# ACER ConQuest 5.47.5 (MML + EAP) vs. BayesMFRM (Stan HMC/NUTS)
#
# Usage (from repository root):
#   Rscript --vanilla simulations/conquest_comparison.R [B=10] [P_exp1=200] [P_exp2=120] [out_dir=simulations/results/conquest]
#
# Structure:
#   - Experiment 1 (Cross-software parity under well-connected balanced data, P = 200):
#       M1_RSM            : Rating Scale Model (K = 5)
#       M2_PCM_Uniform    : Uniform Partial Credit Model (K = 4, item-specific steps)
#       M3_PCM_Hetero     : Heterogeneous-Category PCM (K = c(2, 3, 4, 5, 4))
#       M4_RSM_Bias       : RSM with two-way rater x item interaction (K = 4)
#     Evaluates both a weakly-informative marginal prior (normal(0, 2)) and the
#     recommended hierarchical marginal prior (normal(0, 1)), comparing item,
#     rater, step, interaction, population SD, person EAP, fair scores, and
#     Infit/Outfit MSQ against ACER ConQuest.
#
#   - Experiment 2 (Sparse & unbalanced design + localized rater x item bias, B = 10):
#       Compares unpenalized ConQuest MML against BayesMFRM hierarchical shrinkage
#       models (Gaussian, Student-t, and prior-RMS-calibrated Horseshoe) across
#       replications on RMSE (all cells vs. low-count cells), examinee ability &
#       standardized fair-score recovery (10-rating vs. 3-rating & extreme examinees),
#       95% interval coverage, and empirical Hessian stability.

args <- commandArgs(trailingOnly = TRUE)
B_reps <- if (length(args) >= 1L) as.integer(args[1L]) else 10L
P_exp1 <- if (length(args) >= 2L) as.integer(args[2L]) else 200L
P_exp2 <- if (length(args) >= 3L) as.integer(args[3L]) else 120L
out_dir <- if (length(args) >= 4L) args[4L] else file.path("simulations", "results", "conquest")
conquest_bin <- Sys.getenv("CONQUEST_BIN", "/Applications/ConQuest/ConQuest")

stopifnot(!is.na(B_reps), B_reps >= 1L, !is.na(P_exp1), P_exp1 >= 40L, !is.na(P_exp2), P_exp2 >= 30L)
if (!file.exists(conquest_bin)) {
  stop("ACER ConQuest executable not found at: ", conquest_bin)
}

req_pkgs <- c("pkgload", "cmdstanr", "posterior", "digest", "jsonlite")
missing_pkgs <- req_pkgs[!vapply(req_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0L) {
  stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "))
}

suppressPackageStartupMessages(
  pkgload::load_all(".", quiet = TRUE, export_all = FALSE, export_imports = FALSE)
)
if (!check_cmdstan()) stop("CmdStan is required to run simulations/conquest_comparison.R")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
stan_cache_dir <- tempfile("bmfrm_cq_stan_cache_")
dir.create(stan_cache_dir, recursive = TRUE, showWarnings = FALSE)
on.exit(unlink(stan_cache_dir, recursive = TRUE), add = TRUE)

# -----------------------------------------------------------------------------
# 0. Provenance & Reproducibility Helpers
# -----------------------------------------------------------------------------
git_info <- function(...) {
  tryCatch(
    suppressWarnings(system2("git", c(...), stdout = TRUE, stderr = FALSE)),
    error = function(e) character(0)
  )
}

probe_conquest_version <- function(bin) {
  tmp <- tempfile("cq_ver_")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  writeLines("quit;", file.path(tmp, "q.cqc"))
  out <- system2(bin, stdin = file.path(tmp, "q.cqc"), stdout = TRUE, stderr = TRUE)
  ver_line <- grep("ConQuest version:", out, value = TRUE)
  if (length(ver_line)) trimws(ver_line[1L]) else "unknown"
}

cq_version_str <- probe_conquest_version(conquest_bin)
cat(sprintf("=== ACER ConQuest vs. BayesMFRM Simulation Benchmark ===\n"))
cat(sprintf("ConQuest : %s (%s)\n", cq_version_str, conquest_bin))
cat(sprintf("CmdStan  : %s\n", as.character(cmdstanr::cmdstan_version())))
cat(sprintf("Settings : Exp1 P=%d | Exp2 P=%d, B=%d replications | 4 chains x (1000 warmup + 1000 sampling)\n\n",
            P_exp1, P_exp2, B_reps))

r_src_files <- sort(c(
  "DESCRIPTION", "NAMESPACE",
  file.path("R", list.files("R", "\\.R$")),
  file.path("inst", list.files("inst", recursive = TRUE))
))
r_src_files <- r_src_files[file.exists(r_src_files)]

manifest <- list(
  timestamp_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  git_branch = paste(git_info("branch", "--show-current"), collapse = ""),
  git_commit = paste(git_info("rev-parse", "HEAD"), collapse = ""),
  git_dirty = length(git_info("status", "--porcelain")) > 0L,
  r_version = R.version.string,
  bayesmfrm_version = as.character(utils::packageVersion("BayesMFRM")),
  cmdstan_version = as.character(cmdstanr::cmdstan_version()),
  conquest_version = cq_version_str,
  conquest_binary = conquest_bin,
  config = list(
    B_replications_exp2 = B_reps,
    P_exp1 = P_exp1,
    P_exp2 = P_exp2,
    mcmc = list(
      chains = 4L, cores = 4L, warmup = 1000L, iter = 2000L,
      total_post_draws = 4000L, adapt_delta_main = 0.95,
      adapt_delta_shrinkage = 0.98, max_treedepth = 12L,
      save_gq = "minimal"
    ),
    conquest = list(
      method = "gauss", nodes = 25L, converge = 0.0001,
      deviancechange = 0.0001, max_iter = 2000L,
      n_plausible = 50L, fitdraws = 50L,
      constraints = "items", stderr_primary = "empirical",
      stderr_fallback = "quick"
    ),
    seeds = list(
      exp1_data_base = 20260100L,
      exp1_stan_base = 71000L,
      exp1_cq_base = 101L,
      exp2_data_base = 20261127L,
      exp2_stan_base = 81000L,
      exp2_cq_base = 200L,
      horseshoe_calibration_seed = 91001L
    )
  ),
  source_md5 = as.list(tools::md5sum(r_src_files))
)
jsonlite::write_json(manifest, file.path(out_dir, "manifest.json"), pretty = TRUE, auto_unbox = TRUE)
utils::capture.output(utils::sessionInfo(), file = file.path(out_dir, "session_info.txt"))

# -----------------------------------------------------------------------------
# 1. Ordinal Probability, Fair-Score, & Simulation Helpers
# -----------------------------------------------------------------------------
double_center <- function(M) {
  M <- sweep(M, 1L, rowMeans(M), "-")
  sweep(M, 2L, colMeans(M), "-")
}

ordinal_probs <- function(eta, tau_vec) {
  K <- length(tau_vec) + 1L
  logits <- sweep(outer(eta, 0:(K - 1L)), 2L, c(0, cumsum(tau_vec)), "-")
  logits <- logits - apply(logits, 1L, max)
  z <- exp(logits)
  z / rowSums(z)
}

compute_point_fair_scores <- function(theta_vec, item_vec, tau_list) {
  # Standardized expected score across all items 1..I with rater = 0 and bias = 0
  P <- length(theta_vec)
  I <- length(item_vec)
  acc <- numeric(P)
  for (i in seq_len(I)) {
    pr <- ordinal_probs(theta_vec - item_vec[i], tau_list[[i]])
    Ki <- ncol(pr)
    acc <- acc + as.vector(pr %*% seq_len(Ki)) / I
  }
  acc
}

simulate_mfrm_dataset <- function(design_df, theta, item, rater, tau_list, bias_mat = NULL, seed = NULL) {
  run_sim <- function() {
    N <- nrow(design_df)
    p_idx <- as.integer(sub("^P", "", design_df$person))
    i_idx <- as.integer(sub("^I", "", design_df$item))
    r_idx <- as.integer(sub("^R", "", design_df$rater))
    eta <- theta[p_idx] - item[i_idx] - rater[r_idx]
    if (!is.null(bias_mat)) {
      eta <- eta + bias_mat[cbind(i_idx, r_idx)]
    }
    y <- integer(N)
    u <- stats::runif(N)
    for (i in unique(i_idx)) {
      loc <- which(i_idx == i)
      pr <- ordinal_probs(eta[loc], tau_list[[i]])
      Ki <- ncol(pr)
      cum_p <- t(apply(pr, 1L, cumsum))
      y[loc] <- 1L + rowSums(u[loc] > cum_p[, -Ki, drop = FALSE])
    }
    out <- design_df
    out$score <- as.integer(y)
    out
  }
  if (is.null(seed)) return(run_sim())
  existed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (existed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (existed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(seed)
  run_sim()
}

# -----------------------------------------------------------------------------
# 2. ACER ConQuest 5.47.5 Batch Runner & Self-Checking Output Parser
# -----------------------------------------------------------------------------
parse_conquest_term_block <- function(lines, header_regex) {
  h_idx <- grep(header_regex, lines)
  if (!length(h_idx)) return(character(0))
  start <- h_idx[1L]
  sep_eq <- grep("^={20,}", lines)
  end <- sep_eq[sep_eq > start][1L]
  if (is.na(end)) end <- length(lines)
  lines[start:end]
}

parse_se_token <- function(tok_vec) {
  clean <- ifelse(nzchar(tok_vec) & tok_vec != ".", tok_vec, NA_character_)
  as.numeric(clean)
}

parse_1d_shw_table <- function(block_lines, expected_n) {
  # Matches e.g.: "1   1   -0.791   0.037   1.16 (0.80,1.20)  1.6   1.16 (0.80,1.20)  1.5"
  # or constrained '*' / missing empirical SE '.'
  rx <- "^\\s*(\\d+)\\s+(\\S+)\\s+([-+]?\\d+\\.\\d+)(\\*?)(?:\\s+(\\d+\\.\\d+|\\.))?\\s+(\\d+\\.\\d+)\\s+\\([^)]+\\)\\s+[-+]?\\d+\\.\\d+\\s+(\\d+\\.\\d+)"
  hits <- regmatches(block_lines, regexec(rx, block_lines, perl = TRUE))
  hits <- Filter(function(m) length(m) > 0L, hits)
  if (length(hits) != expected_n) {
    stop(sprintf("Expected %d rows in ConQuest 1D table, parsed %d", expected_n, length(hits)))
  }
  data.frame(
    idx = as.integer(vapply(hits, `[`, character(1), 2L)),
    label = vapply(hits, `[`, character(1), 3L),
    shw_est = as.numeric(vapply(hits, `[`, character(1), 4L)),
    constrained = nzchar(vapply(hits, `[`, character(1), 5L)),
    se = parse_se_token(vapply(hits, `[`, character(1), 6L)),
    outfit_msq = as.numeric(vapply(hits, `[`, character(1), 7L)),
    infit_msq = as.numeric(vapply(hits, `[`, character(1), 8L)),
    stringsAsFactors = FALSE
  )
}

parse_rsm_step_shw <- function(block_lines, K) {
  rx <- "^\\s*(\\d+)\\s+([-+]?\\d+\\.\\d+)(\\*?)(?:\\s+(\\d+\\.\\d+|\\.))?\\s+(\\d+\\.\\d+)"
  hits <- regmatches(block_lines, regexec(rx, block_lines, perl = TRUE))
  hits <- Filter(function(m) length(m) > 0L && as.integer(m[2L]) >= 1L, hits)
  if (length(hits) != (K - 1L)) {
    stop(sprintf("Expected %d RSM steps in ConQuest table, parsed %d", K - 1L, length(hits)))
  }
  data.frame(
    step = as.integer(vapply(hits, `[`, character(1), 2L)),
    shw_est = as.numeric(vapply(hits, `[`, character(1), 3L)),
    se = parse_se_token(vapply(hits, `[`, character(1), 5L)),
    stringsAsFactors = FALSE
  )
}

parse_pcm_step_shw <- function(block_lines, K_vec) {
  rx <- "^\\s*(\\d+)\\s+(\\S+)\\s+(\\d+)\\s+([-+]?\\d+\\.\\d+)(\\*?)(?:\\s+(\\d+\\.\\d+|\\.))?\\s+(\\d+\\.\\d+)"
  hits <- regmatches(block_lines, regexec(rx, block_lines, perl = TRUE))
  hits <- Filter(function(m) length(m) > 0L && as.integer(m[4L]) >= 1L, hits)
  df <- data.frame(
    item_idx = as.integer(vapply(hits, `[`, character(1), 2L)),
    step = as.integer(vapply(hits, `[`, character(1), 4L)),
    shw_est = as.numeric(vapply(hits, `[`, character(1), 5L)),
    se = parse_se_token(vapply(hits, `[`, character(1), 7L)),
    stringsAsFactors = FALSE
  )
  expected_rows <- sum(ifelse(K_vec >= 3L, K_vec - 1L, 0L))
  if (nrow(df) != expected_rows) {
    stop(sprintf("Expected %d PCM step rows in ConQuest table, parsed %d", expected_rows, nrow(df)))
  }
  df
}

parse_bias_shw <- function(block_lines, I, R) {
  rx <- "^\\s*(\\d+)\\s+(\\S+)\\s+(\\d+)\\s+(\\S+)\\s+([-+]?\\d+\\.\\d+)(\\*?)(?:\\s+(\\d+\\.\\d+|\\.))?\\s+(\\d+\\.\\d+)\\s+\\([^)]+\\)\\s+[-+]?\\d+\\.\\d+\\s+(\\d+\\.\\d+)"
  hits <- regmatches(block_lines, regexec(rx, block_lines, perl = TRUE))
  hits <- Filter(function(m) length(m) > 0L, hits)
  if (length(hits) != I * R) {
    stop(sprintf("Expected %d interaction cells in ConQuest table, parsed %d", I * R, length(hits)))
  }
  data.frame(
    item_idx = as.integer(vapply(hits, `[`, character(1), 2L)),
    rater_idx = as.integer(vapply(hits, `[`, character(1), 4L)),
    shw_est = as.numeric(vapply(hits, `[`, character(1), 6L)),
    se = parse_se_token(vapply(hits, `[`, character(1), 8L)),
    outfit_msq = as.numeric(vapply(hits, `[`, character(1), 9L)),
    infit_msq = as.numeric(vapply(hits, `[`, character(1), 10L)),
    stringsAsFactors = FALSE
  )
}

run_conquest_mfrm <- function(df, P, I, R, K_vec,
                              model_type = c("rsm", "pcm", "hetero_pcm", "rsm_bias"),
                              cq_seed = 101L,
                              conquest_exe = conquest_bin) {
  model_type <- match.arg(model_type)
  if (length(K_vec) == 1L) K_vec <- rep(as.integer(K_vec), I)
  K_max <- max(K_vec)
  cq_seed_safe <- as.integer(((cq_seed - 1L) %% 900L) + 10L)

  work <- tempfile("cq_run_")
  dir.create(work, recursive = TRUE)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)

  # Convert long (person, item, rater, score in 1..K) to ConQuest wide (pid, rater, item_1..item_I in 0..K-1)
  p_num <- as.integer(sub("^P", "", df$person))
  i_num <- as.integer(sub("^I", "", df$item))
  r_num <- as.integer(sub("^R", "", df$rater))
  y0 <- as.integer(df$score - 1L)

  pr_pairs <- unique(data.frame(p = p_num, r = r_num))
  pr_pairs <- pr_pairs[order(pr_pairs$p, pr_pairs$r), , drop = FALSE]
  resp_mat <- matrix(" ", nrow = nrow(pr_pairs), ncol = I)
  row_map <- match(paste(p_num, r_num, sep = ":"), paste(pr_pairs$p, pr_pairs$r, sep = ":"))
  for (k in seq_len(nrow(df))) {
    resp_mat[row_map[k], i_num[k]] <- as.character(y0[k])
  }
  dat_lines <- sprintf(
    "%04d%02d%s",
    pr_pairs$p, pr_pairs$r,
    apply(resp_mat, 1L, paste0, collapse = "")
  )
  writeLines(dat_lines, file.path(work, "dat.txt"))

  cq_model_str <- switch(
    model_type,
    rsm        = "item + rater + step",
    pcm        = "item + rater + item*step",
    hetero_pcm = "item + rater + item*step",
    rsm_bias   = "item + rater - item*rater + step"
  )
  codes_str <- paste(0:(K_max - 1L), collapse = ",")

  exec_once <- function(se_opt) {
    cqc <- c(
      "datafile dat.txt;",
      sprintf("format pid 1-4 rater 5-6 responses 7-%d;", 6L + I),
      sprintf("codes %s;", codes_str),
      sprintf("model %s;", cq_model_str),
      "set constraints=items;",
      sprintf("set seed=%d;", cq_seed_safe),
      "set n_plausible=50;",
      "set fitdraws=50;",
      "set warnings=no;",
      sprintf("estimate ! method=gauss, nodes=25, converge=0.0001, deviancechange=0.0001, iter=2000, stderr=%s;", se_opt),
      "export par >> par.txt;",
      "show >> show.shw ! estimates=latent, tables=1:2:3;",
      "show cases >> eap.out ! estimates=eap;",
      "quit;"
    )
    cqc_path <- file.path(work, paste0("run_", se_opt, ".cqc"))
    writeLines(cqc, cqc_path)
    t0 <- proc.time()[["elapsed"]]
    old_wd <- getwd()
    on.exit(setwd(old_wd), add = TRUE)
    setwd(work)
    log_out <- suppressWarnings(system2(conquest_exe, stdin = basename(cqc_path), stdout = TRUE, stderr = TRUE))
    elapsed <- proc.time()[["elapsed"]] - t0
    ok_shw <- file.exists(file.path(work, "show.shw.txt"))
    shw_chk <- if (ok_shw) readLines(file.path(work, "show.shw.txt"), warn = FALSE) else character(0)
    se_failed <- any(grepl(
      "Covariance matrix is not positive definite|calculation of variance-covariance matrix failed|Empirical standard errors not available|fatal",
      c(log_out, shw_chk), ignore.case = TRUE
    ))
    list(log = log_out, elapsed = elapsed, ok = ok_shw, se_failed = se_failed)
  }

  r1 <- exec_once("empirical")
  empirical_se_failed <- FALSE
  se_mode <- "empirical"
  run_used <- r1
  if (!r1$ok || isTRUE(r1$se_failed)) {
    empirical_se_failed <- TRUE
    se_mode <- "quick"
    unlink(file.path(work, c("show.shw.txt", "eap.out.txt", "par.txt")))
    run_used <- exec_once("quick")
  }
  if (!run_used$ok) {
    stop("ConQuest failed to produce show.shw.txt:\n", paste(utils::tail(run_used$log, 30L), collapse = "\n"))
  }

  shw_lines <- readLines(file.path(work, "show.shw.txt"), warn = FALSE)
  par_lines <- readLines(file.path(work, "par.txt"), warn = FALSE)

  # Convergence & deviance metadata
  conv_line <- grep("^Iterations terminated because", shw_lines, value = TRUE)
  converged <- length(conv_line) > 0L && grepl("convergence criteria", conv_line[1L], ignore.case = TRUE)
  it_line <- grep("^The number of iterations:", shw_lines, value = TRUE)
  n_iter <- if (length(it_line)) as.integer(sub(".*:\\s*(\\d+).*", "\\1", it_line[1L])) else NA_integer_
  dev_line <- grep("^Final Deviance:", shw_lines, value = TRUE)
  deviance <- if (length(dev_line)) as.numeric(sub(".*:\\s*([-+]?\\d+\\.\\d+).*", "\\1", dev_line[1L])) else NA_real_

  # Population mean & SD
  int_line <- grep("^Intercept\\s+", shw_lines, value = TRUE)
  mu_theta <- if (length(int_line)) as.numeric(strsplit(trimws(int_line[1L]), "\\s+")[[1L]][2L]) else 0
  var_line <- grep("^Variance\\s+", shw_lines, value = TRUE)
  var_theta <- if (length(var_line)) as.numeric(strsplit(trimws(var_line[1L]), "\\s+")[[1L]][2L]) else NA_real_
  sigma_theta <- sqrt(max(var_theta, 0))

  # Parse 1D item & rater tables from show.shw.txt
  item_shw <- parse_1d_shw_table(parse_conquest_term_block(shw_lines, "^TERM 1:\\s+item\\b"), I)
  rater_shw <- parse_1d_shw_table(parse_conquest_term_block(shw_lines, "^TERM 2:\\s+rater\\b"), R)

  # Parse 5-decimal free parameters from par.txt and reconstruct exact sum-to-zero vectors!
  par_vals <- as.numeric(vapply(strsplit(trimws(par_lines[nzchar(trimws(par_lines))]), "\\s+"), `[`, character(1), 2L))
  pos <- 0L
  take_par <- function(n) {
    if (n <= 0L) return(numeric(0))
    idx <- pos + seq_len(n)
    pos <<- pos + n
    par_vals[idx]
  }
  item_free <- take_par(I - 1L)
  item_est <- c(item_free, -sum(item_free))
  rater_free <- take_par(R - 1L)
  rater_est <- c(rater_free, -sum(rater_free))

  # Self-check 5-decimal reconstruction against ConQuest's printed 3-decimal table
  stopifnot(max(abs(item_est - item_shw$shw_est)) < 0.002)
  stopifnot(max(abs(rater_est - rater_shw$shw_est)) < 0.002)
  item_shw$est <- item_est
  rater_shw$est <- rater_est

  bias_mat <- NULL
  bias_se_mat <- NULL
  tau_list <- vector("list", I)
  tau_se_list <- vector("list", I)

  if (model_type == "rsm") {
    step_shw <- parse_rsm_step_shw(parse_conquest_term_block(shw_lines, "^TERM 3:\\s+step\\b"), K_max)
    step_free <- take_par(K_max - 2L)
    step_est <- c(step_free, -sum(step_free))
    stopifnot(max(abs(step_est - step_shw$shw_est)) < 0.002)
    step_shw$est <- step_est
    for (i in seq_len(I)) {
      tau_list[[i]] <- step_est
      tau_se_list[[i]] <- step_shw$se
    }
  } else if (model_type %in% c("pcm", "hetero_pcm")) {
    pstep_shw <- parse_pcm_step_shw(parse_conquest_term_block(shw_lines, "^TERM 3:\\s+item\\*step\\b"), K_vec)
    for (i in seq_len(I)) {
      Ki <- K_vec[i]
      if (Ki == 2L) {
        tau_list[[i]] <- 0
        tau_se_list[[i]] <- NA_real_
      } else {
        sf <- take_par(Ki - 2L)
        se_i <- c(sf, -sum(sf))
        sub_df <- pstep_shw[pstep_shw$item_idx == i, , drop = FALSE]
        stopifnot(max(abs(se_i - sub_df$shw_est)) < 0.002)
        tau_list[[i]] <- se_i
        tau_se_list[[i]] <- sub_df$se
      }
    }
  } else if (model_type == "rsm_bias") {
    b_shw <- parse_bias_shw(parse_conquest_term_block(shw_lines, "^TERM 3:\\s+\\(-\\)item\\*rater\\b"), I, R)
    # ConQuest stores free interaction parameters in rater-major order matching TERM 3 listing:
    # r = 1..(R-1), i = 1..(I-1) -> verify orientation against b_shw$shw_est!
    b_free <- take_par((I - 1L) * (R - 1L))
    complete_2d <- function(mat_free) {
      m <- matrix(0, nrow = I, ncol = R)
      m[seq_len(I - 1L), seq_len(R - 1L)] <- mat_free
      m[I, seq_len(R - 1L)] <- -colSums(m[seq_len(I - 1L), seq_len(R - 1L), drop = FALSE])
      m[, R] <- -rowSums(m[, seq_len(R - 1L), drop = FALSE])
      m
    }
    cand_colmajor <- complete_2d(matrix(b_free, nrow = I - 1L, ncol = R - 1L, byrow = FALSE))
    cand_rowmajor <- complete_2d(matrix(b_free, nrow = I - 1L, ncol = R - 1L, byrow = TRUE))
    shw_mat <- matrix(NA_real_, nrow = I, ncol = R)
    se_mat <- matrix(NA_real_, nrow = I, ncol = R)
    for (row_k in seq_len(nrow(b_shw))) {
      shw_mat[b_shw$item_idx[row_k], b_shw$rater_idx[row_k]] <- b_shw$shw_est[row_k]
      se_mat[b_shw$item_idx[row_k], b_shw$rater_idx[row_k]] <- b_shw$se[row_k]
    }
    err_col <- max(abs(cand_colmajor - shw_mat))
    err_row <- max(abs(cand_rowmajor - shw_mat))
    bias_mat <- if (err_col <= err_row && err_col < 0.002) {
      cand_colmajor
    } else if (err_row < 0.002) {
      cand_rowmajor
    } else {
      shw_mat
    }
    bias_se_mat <- se_mat

    step_shw <- parse_rsm_step_shw(parse_conquest_term_block(shw_lines, "^TERM 4:\\s+step\\b"), K_max)
    step_free <- take_par(K_max - 2L)
    step_est <- c(step_free, -sum(step_free))
    stopifnot(max(abs(step_est - step_shw$shw_est)) < 0.002)
    step_shw$est <- step_est
    for (i in seq_len(I)) {
      tau_list[[i]] <- step_est
      tau_se_list[[i]] <- step_shw$se
    }
  }

  # Parse case EAP estimates from eap.out.txt
  eap_tab <- utils::read.table(file.path(work, "eap.out.txt"), header = FALSE, stringsAsFactors = FALSE)
  nc_eap <- ncol(eap_tab)
  pid_int <- as.integer(eap_tab[[2L]])
  eap_est <- rep(NA_real_, P)
  eap_sd <- rep(NA_real_, P)
  eap_est[pid_int] <- as.numeric(eap_tab[[nc_eap - 2L]])
  eap_sd[pid_int] <- as.numeric(eap_tab[[nc_eap - 1L]])

  fair_est <- compute_point_fair_scores(eap_est, item_shw$est, tau_list)

  list(
    converged = converged,
    n_iter = n_iter,
    deviance = deviance,
    empirical_se_failed = empirical_se_failed,
    se_mode = se_mode,
    elapsed_sec = run_used$elapsed,
    mu_theta = mu_theta,
    sigma_theta = sigma_theta,
    item = item_shw,
    rater = rater_shw,
    tau = tau_list,
    tau_se = tau_se_list,
    step_table = if (exists("step_shw", inherits = FALSE)) step_shw else pstep_shw,
    bias = bias_mat,
    bias_se = bias_se_mat,
    eap = eap_est,
    eap_sd = eap_sd,
    fair_score = fair_est
  )
}

# Helper to extract BayesMFRM posterior summaries aligned 1-to-1 with ground truth & ConQuest
extract_bayes_estimates <- function(fit, P, I, R, K_vec) {
  if (length(K_vec) == 1L) K_vec <- rep(as.integer(K_vec), I)
  p_tab <- facet_summary(fit, "person")
  i_tab <- facet_summary(fit, "item")
  r_tab <- facet_summary(fit, "rater")
  i_fit <- facet_fit(fit, "item")
  r_fit <- facet_fit(fit, "rater")
  fs_tab <- suppressMessages(fair_scores(fit, adjust_facets = "rater", standardize = TRUE))

  p_ord <- match(sprintf("P%03d", seq_len(P)), p_tab$level)
  i_ord <- match(sprintf("I%02d", seq_len(I)), i_tab$level)
  r_ord <- match(sprintf("R%02d", seq_len(R)), r_tab$level)
  if_ord <- match(sprintf("I%02d", seq_len(I)), i_fit$label)
  rf_ord <- match(sprintf("R%02d", seq_len(R)), r_fit$label)
  fs_ord <- match(sprintf("P%03d", seq_len(P)), fs_tab$person)

  mu_d <- as.numeric(unclass(fit$fit$draws("mu_theta", format = "matrix")))
  sig_d <- as.numeric(unclass(fit$fit$draws("sigma_theta", format = "matrix")))

  tau_draws <- unclass(fit$fit$draws("tau", format = "matrix"))
  tau_mean_list <- vector("list", I)
  tau_sd_list <- vector("list", I)
  tau_lo_list <- vector("list", I)
  tau_hi_list <- vector("list", I)

  if (identical(fit$spec$family, "rsm")) {
    K_u <- K_vec[1L]
    cols <- paste0("tau[", seq_len(K_u - 1L), "]")
    tm <- unname(colMeans(tau_draws[, cols, drop = FALSE]))
    ts <- unname(apply(tau_draws[, cols, drop = FALSE], 2L, stats::sd))
    tl <- unname(apply(tau_draws[, cols, drop = FALSE], 2L, stats::quantile, 0.025))
    th <- unname(apply(tau_draws[, cols, drop = FALSE], 2L, stats::quantile, 0.975))
    for (i in seq_len(I)) {
      tau_mean_list[[i]] <- tm; tau_sd_list[[i]] <- ts
      tau_lo_list[[i]] <- tl; tau_hi_list[[i]] <- th
    }
  } else if (isTRUE(fit$spec$K_heterogeneous)) {
    pos <- 1L
    for (i in seq_len(I)) {
      Ki <- K_vec[i]
      len_i <- Ki - 1L
      cols <- paste0("tau[", pos:(pos + len_i - 1L), "]")
      pos <- pos + len_i
      tau_mean_list[[i]] <- unname(colMeans(tau_draws[, cols, drop = FALSE]))
      tau_sd_list[[i]] <- unname(apply(tau_draws[, cols, drop = FALSE], 2L, stats::sd))
      tau_lo_list[[i]] <- unname(apply(tau_draws[, cols, drop = FALSE], 2L, stats::quantile, 0.025))
      tau_hi_list[[i]] <- unname(apply(tau_draws[, cols, drop = FALSE], 2L, stats::quantile, 0.975))
    }
  } else {
    K_u <- K_vec[1L]
    for (i in seq_len(I)) {
      cols <- paste0("tau[", i, ",", seq_len(K_u - 1L), "]")
      tau_mean_list[[i]] <- unname(colMeans(tau_draws[, cols, drop = FALSE]))
      tau_sd_list[[i]] <- unname(apply(tau_draws[, cols, drop = FALSE], 2L, stats::sd))
      tau_lo_list[[i]] <- unname(apply(tau_draws[, cols, drop = FALSE], 2L, stats::quantile, 0.025))
      tau_hi_list[[i]] <- unname(apply(tau_draws[, cols, drop = FALSE], 2L, stats::quantile, 0.975))
    }
  }

  bias_mean <- bias_sd <- bias_lo <- bias_hi <- NULL
  if (length(fit$spec$facets_bias) > 0L) {
    bt <- summarise_bias(fit, "item:rater", diagnostics = TRUE)
    bias_mean <- bias_sd <- bias_lo <- bias_hi <- matrix(NA_real_, nrow = I, ncol = R)
    for (k in seq_len(nrow(bt))) {
      ii <- match(bt$item[k], sprintf("I%02d", seq_len(I)))
      rr <- match(bt$rater[k], sprintf("R%02d", seq_len(R)))
      bias_mean[ii, rr] <- bt$bias_mean[k]
      bias_sd[ii, rr] <- bt$bias_sd[k]
      bias_lo[ii, rr] <- bt$lower_95[k]
      bias_hi[ii, rr] <- bt$upper_95[k]
    }
  }

  list(
    converged = isTRUE(fit$convergence$converged),
    max_rhat = fit$convergence$max_rhat,
    min_ess_bulk = fit$convergence$min_ess_bulk,
    min_ess_tail = fit$convergence$min_ess_tail,
    n_divergent = fit$convergence$n_divergent,
    n_max_treedepth = fit$convergence$n_max_treedepth,
    min_ebfmi = fit$convergence$min_ebfmi,
    max_mcse_mean = fit$convergence$max_mcse_mean,
    mu_theta = mean(mu_d),
    sigma_theta = mean(sig_d),
    theta = p_tab$mean[p_ord],
    theta_sd = p_tab$sd[p_ord],
    theta_lo = p_tab$q2.5[p_ord],
    theta_hi = p_tab$q97.5[p_ord],
    item = i_tab$mean[i_ord],
    item_sd = i_tab$sd[i_ord],
    item_lo = i_tab$q2.5[i_ord],
    item_hi = i_tab$q97.5[i_ord],
    rater = r_tab$mean[r_ord],
    rater_sd = r_tab$sd[r_ord],
    rater_lo = r_tab$q2.5[r_ord],
    rater_hi = r_tab$q97.5[r_ord],
    item_infit_msq = i_fit$infit_msq[if_ord],
    item_infit_plugin = i_fit$infit_msq_plugin[if_ord],
    item_outfit_msq = i_fit$outfit_msq[if_ord],
    item_outfit_plugin = i_fit$outfit_msq_plugin[if_ord],
    rater_infit_msq = r_fit$infit_msq[rf_ord],
    rater_infit_plugin = r_fit$infit_msq_plugin[rf_ord],
    rater_outfit_msq = r_fit$outfit_msq[rf_ord],
    rater_outfit_plugin = r_fit$outfit_msq_plugin[rf_ord],
    tau = tau_mean_list,
    tau_sd = tau_sd_list,
    tau_lo = tau_lo_list,
    tau_hi = tau_hi_list,
    bias = bias_mean,
    bias_sd = bias_sd,
    bias_lo = bias_lo,
    bias_hi = bias_hi,
    fair_score = fs_tab$mean_fair_score[fs_ord],
    fair_lo = fs_tab$lower_95[fs_ord],
    fair_hi = fs_tab$upper_95[fs_ord]
  )
}

# =============================================================================
# EXPERIMENT 1: Balanced Connected Parity Across All 4 Model Families (P = P_exp1)
# =============================================================================
cat("=== EXPERIMENT 1: 4-Model Cross-Software Parity (P =", P_exp1, ") ===\n")
I1 <- 5L; R1 <- 4L

# Balanced cyclic rater bridge: each person rated by 2 adjacent raters on all 5 items (2,000 ratings)
build_exp1_design <- function(P, I, R) {
  rows <- lapply(seq_len(P), function(p) {
    r1 <- ((p - 1L) %% R) + 1L
    r2 <- (p %% R) + 1L
    expand.grid(
      person = sprintf("P%03d", p),
      item = sprintf("I%02d", seq_len(I)),
      rater = sprintf("R%02d", c(r1, r2)),
      stringsAsFactors = FALSE
    )
  })
  d <- do.call(rbind, rows)
  rownames(d) <- NULL
  d
}

des_exp1 <- build_exp1_design(P_exp1, I1, R1)
set.seed(20260100L)
theta_exp1 <- stats::rnorm(P_exp1, mean = 0.10, sd = 1.05)
item_exp1 <- c(-0.90, -0.35, 0.05, 0.45, 0.75)
rater_exp1 <- c(-0.60, -0.15, 0.20, 0.55)

exp1_specs <- list(
  M1_RSM = list(
    label = "1. RSM (K = 5)",
    cq_type = "rsm",
    formula = score ~ person + item + rater,
    family = "rsm",
    threshold_facet = NULL,
    K_spec = 5L,
    K_vec = rep(5L, I1),
    tau_true = replicate(I1, c(-1.20, -0.40, 0.40, 1.20), simplify = FALSE),
    bias_true = NULL
  ),
  M2_PCM_Uniform = list(
    label = "2. Uniform PCM (K = 4)",
    cq_type = "pcm",
    formula = score ~ person + item + rater,
    family = "pcm",
    threshold_facet = "item",
    K_spec = 4L,
    K_vec = rep(4L, I1),
    tau_true = list(
      c(-1.10, 0.00, 1.10),
      c(-0.45, 0.00, 0.45),
      c(-0.90, 0.25, 0.65),
      c( 0.20, -0.50, 0.30), # reversed step reversal (valid in PCM!)
      c(-0.70, -0.10, 0.80)
    ),
    bias_true = NULL
  ),
  M3_PCM_Hetero = list(
    label = "3. Hetero PCM (K = 2,3,4,5,4)",
    cq_type = "hetero_pcm",
    formula = score ~ person + item + rater,
    family = "pcm",
    threshold_facet = "item",
    K_spec = c(I01 = 2L, I02 = 3L, I03 = 4L, I04 = 5L, I05 = 4L),
    K_vec = c(2L, 3L, 4L, 5L, 4L),
    tau_true = list(
      0,
      c(-0.75, 0.75),
      c(-0.95, 0.05, 0.90),
      c(-1.25, -0.35, 0.45, 1.15),
      c(-0.60, -0.10, 0.70)
    ),
    bias_true = NULL
  ),
  M4_RSM_Bias = list(
    label = "4. Interaction RSM (item:rater)",
    cq_type = "rsm_bias",
    formula = score ~ person + item + rater + item:rater,
    family = "rsm",
    threshold_facet = NULL,
    K_spec = 4L,
    K_vec = rep(4L, I1),
    tau_true = replicate(I1, c(-0.95, 0.00, 0.95), simplify = FALSE),
    bias_true = matrix(c(
       0.55, -0.55,  0.00,  0.00,
      -0.55,  0.55,  0.00,  0.00,
       0.00,  0.00,  0.45, -0.45,
       0.00,  0.00, -0.45,  0.45,
       0.00,  0.00,  0.00,  0.00
    ), nrow = I1, ncol = R1, byrow = TRUE)
  )
)

make_exp1_priors <- function(has_bias, mode = c("weak_marginal", "default_hier")) {
  mode <- match.arg(mode)
  if (mode == "weak_marginal") {
    ps <- c(
      prior("normal(0, 2)", "mean", facet = "person"),
      prior("lognormal(0, 0.5)", "sd", facet = "person"),
      prior("normal(0, 2)", "item", scale = "marginal"),
      prior("normal(0, 2)", "rater", scale = "marginal"),
      prior("normal(0, 2)", "tau", scale = "marginal")
    )
    hier <- "person"
    if (has_bias) {
      ps <- c(ps, list(prior("normal(0, 2)", "bias", facet = "item:rater", scale = "marginal")))
    }
  } else {
    ps <- c(
      prior("normal(0, 2)", "mean", facet = "person"),
      prior("lognormal(0, 0.5)", "sd", facet = "person"),
      prior("normal(0, 1)", "item", scale = "marginal"),
      prior("normal(0, 1)", "sd", facet = "rater", scale = "marginal"),
      prior("normal(0, 1)", "tau", scale = "marginal")
    )
    hier <- c("person", "rater")
    if (has_bias) {
      ps <- c(ps, list(prior("normal(0, 0.6)", "sd", facet = "item:rater", scale = "marginal")))
      hier <- c(hier, "item:rater")
    }
  }
  list(priors = ps, hierarchical = hier)
}

exp1_param_rows <- list()
exp1_person_rows <- list()
exp1_fitstat_rows <- list()
exp1_summary_rows <- list()

for (m_idx in seq_along(exp1_specs)) {
  m_id <- names(exp1_specs)[m_idx]
  sp <- exp1_specs[[m_id]]
  cat(sprintf("[Exp 1 | %s] Simulating & running ConQuest + BayesMFRM...\n", sp$label))

  dat_m <- simulate_mfrm_dataset(
    des_exp1, theta_exp1, item_exp1, rater_exp1, sp$tau_true, sp$bias_true,
    seed = 20260100L + m_idx
  )
  true_fair <- compute_point_fair_scores(theta_exp1, item_exp1, sp$tau_true)

  cq_res <- run_conquest_mfrm(
    dat_m, P = P_exp1, I = I1, R = R1, K_vec = sp$K_vec,
    model_type = sp$cq_type, cq_seed = 100L + m_idx
  )

  for (p_mode in c("weak_marginal", "default_hier")) {
    pr_cfg <- make_exp1_priors(!is.null(sp$bias_true), p_mode)
    t0 <- proc.time()[["elapsed"]]
    bfit <- suppressMessages(bmfrm(
      formula = sp$formula,
      data = dat_m,
      K = sp$K_spec,
      family = sp$family,
      threshold_facet = sp$threshold_facet,
      hierarchical = pr_cfg$hierarchical,
      priors = pr_cfg$priors,
      save_gq = "minimal",
      cache_dir = stan_cache_dir,
      model_name = paste0("exp1_", m_id, "_", p_mode),
      design_check = "error",
      chains = 4L, cores = 4L, iter = 2000L, warmup = 1000L,
      adapt_delta = if (p_mode == "default_hier") 0.99 else 0.95, max_treedepth = 12L,
      seed = 71000L + m_idx * 10L + match(p_mode, c("weak_marginal", "default_hier")),
      refresh = 0L, show_messages = FALSE, show_exceptions = FALSE,
      cpp_options = list(PRECOMPILED_HEADERS = "false")
    ))
    b_elapsed <- proc.time()[["elapsed"]] - t0
    bres <- extract_bayes_estimates(bfit, P = P_exp1, I = I1, R = R1, K_vec = sp$K_vec)

    # Flatten structural parameters (item, rater, free/all steps, optional bias)
    add_p <- function(p_class, p_names, truth_v, cq_v, cq_se_v, by_v, by_sd_v, by_lo_v, by_hi_v) {
      data.frame(
        model_id = m_id, model_label = sp$label, prior_mode = p_mode,
        param_class = p_class, param = p_names,
        truth = truth_v, cq_est = cq_v, cq_se = cq_se_v,
        bayes_mean = by_v, bayes_sd = by_sd_v, bayes_lo = by_lo_v, bayes_hi = by_hi_v,
        stringsAsFactors = FALSE
      )
    }
    p_blocks <- list(
      add_p("item", sprintf("item[%d]", 1:I1), item_exp1, cq_res$item$est, cq_res$item$se,
            bres$item, bres$item_sd, bres$item_lo, bres$item_hi),
      add_p("rater", sprintf("rater[%d]", 1:R1), rater_exp1, cq_res$rater$est, cq_res$rater$se,
            bres$rater, bres$rater_sd, bres$rater_lo, bres$rater_hi)
    )
    if (sp$family == "rsm") {
      Km1 <- sp$K_vec[1L] - 1L
      p_blocks[[length(p_blocks) + 1L]] <- add_p(
        "step", sprintf("tau[%d]", seq_len(Km1)),
        sp$tau_true[[1L]], cq_res$tau[[1L]], cq_res$tau_se[[1L]],
        bres$tau[[1L]], bres$tau_sd[[1L]], bres$tau_lo[[1L]], bres$tau_hi[[1L]]
      )
    } else {
      for (ii in seq_len(I1)) {
        if (sp$K_vec[ii] >= 3L) {
          Km1 <- sp$K_vec[ii] - 1L
          p_blocks[[length(p_blocks) + 1L]] <- add_p(
            "step", sprintf("tau[%d,%d]", ii, seq_len(Km1)),
            sp$tau_true[[ii]], cq_res$tau[[ii]], cq_res$tau_se[[ii]],
            bres$tau[[ii]], bres$tau_sd[[ii]], bres$tau_lo[[ii]], bres$tau_hi[[ii]]
          )
        }
      }
    }
    if (!is.null(sp$bias_true)) {
      grid_ir <- expand.grid(i = seq_len(I1), r = seq_len(R1))
      p_blocks[[length(p_blocks) + 1L]] <- add_p(
        "interaction", sprintf("bias[%d,%d]", grid_ir$i, grid_ir$r),
        as.vector(sp$bias_true), as.vector(cq_res$bias), as.vector(cq_res$bias_se),
        as.vector(bres$bias), as.vector(bres$bias_sd), as.vector(bres$bias_lo), as.vector(bres$bias_hi)
      )
    }
    p_df <- do.call(rbind, p_blocks)
    exp1_param_rows[[length(exp1_param_rows) + 1L]] <- p_df

    if (p_mode == "weak_marginal") {
      exp1_person_rows[[m_id]] <- data.frame(
        model_id = m_id, model_label = sp$label, person = sprintf("P%03d", seq_len(P_exp1)),
        theta_true = theta_exp1, cq_eap = cq_res$eap, cq_eap_sd = cq_res$eap_sd,
        bayes_theta = bres$theta, bayes_theta_sd = bres$theta_sd,
        true_fair = true_fair, cq_fair = cq_res$fair_score, bayes_fair = bres$fair_score,
        stringsAsFactors = FALSE
      )
      exp1_fitstat_rows[[m_id]] <- rbind(
        data.frame(
          model_id = m_id, model_label = sp$label, facet = "item", level = sprintf("I%02d", 1:I1),
          cq_infit = cq_res$item$infit_msq, bayes_infit_post = bres$item_infit_msq, bayes_infit_plugin = bres$item_infit_plugin,
          cq_outfit = cq_res$item$outfit_msq, bayes_outfit_post = bres$item_outfit_msq, bayes_outfit_plugin = bres$item_outfit_plugin,
          stringsAsFactors = FALSE
        ),
        data.frame(
          model_id = m_id, model_label = sp$label, facet = "rater", level = sprintf("R%02d", 1:R1),
          cq_infit = cq_res$rater$infit_msq, bayes_infit_post = bres$rater_infit_msq, bayes_infit_plugin = bres$rater_infit_plugin,
          cq_outfit = cq_res$rater$outfit_msq, bayes_outfit_post = bres$rater_outfit_msq, bayes_outfit_plugin = bres$rater_outfit_plugin,
          stringsAsFactors = FALSE
        )
      )
    }

    exp1_summary_rows[[length(exp1_summary_rows) + 1L]] <- data.frame(
      model_id = m_id,
      model_label = sp$label,
      prior_mode = p_mode,
      cq_converged = cq_res$converged,
      cq_se_mode = cq_res$se_mode,
      cq_iter = cq_res$n_iter,
      cq_sec = round(cq_res$elapsed_sec, 2),
      bayes_converged = bres$converged,
      max_rhat = bres$max_rhat,
      min_ess_bulk = bres$min_ess_bulk,
      min_ess_tail = bres$min_ess_tail,
      divergences = bres$n_divergent,
      bayes_sec = round(b_elapsed, 2),
      param_cor = stats::cor(p_df$cq_est, p_df$bayes_mean),
      param_rmsd_cq_vs_bayes = sqrt(mean((p_df$cq_est - p_df$bayes_mean)^2)),
      param_max_abs_diff = max(abs(p_df$cq_est - p_df$bayes_mean)),
      cq_rmse_truth = sqrt(mean((p_df$cq_est - p_df$truth)^2)),
      bayes_rmse_truth = sqrt(mean((p_df$bayes_mean - p_df$truth)^2)),
      person_eap_cor = stats::cor(cq_res$eap, bres$theta),
      person_eap_rmsd = sqrt(mean((cq_res$eap - bres$theta)^2)),
      fair_score_cor = stats::cor(cq_res$fair_score, bres$fair_score),
      fair_score_rmsd = sqrt(mean((cq_res$fair_score - bres$fair_score)^2)),
      cq_sigma_theta = cq_res$sigma_theta,
      bayes_sigma_theta = bres$sigma_theta,
      infit_rmsd_post = sqrt(mean((c(cq_res$item$infit_msq, cq_res$rater$infit_msq) -
                                   c(bres$item_infit_msq, bres$rater_infit_msq))^2)),
      infit_rmsd_plugin = sqrt(mean((c(cq_res$item$infit_msq, cq_res$rater$infit_msq) -
                                     c(bres$item_infit_plugin, bres$rater_infit_plugin))^2)),
      outfit_rmsd_post = sqrt(mean((c(cq_res$item$outfit_msq, cq_res$rater$outfit_msq) -
                                    c(bres$item_outfit_msq, bres$rater_outfit_msq))^2)),
      stringsAsFactors = FALSE
    )
  }
}

exp1_params_df <- do.call(rbind, exp1_param_rows)
exp1_persons_df <- do.call(rbind, exp1_person_rows)
exp1_fitstats_df <- do.call(rbind, exp1_fitstat_rows)
exp1_summary_df <- do.call(rbind, exp1_summary_rows)

utils::write.csv(exp1_params_df, file.path(out_dir, "exp1_parameters.csv"), row.names = FALSE)
utils::write.csv(exp1_persons_df, file.path(out_dir, "exp1_persons.csv"), row.names = FALSE)
utils::write.csv(exp1_fitstats_df, file.path(out_dir, "exp1_fitstats.csv"), row.names = FALSE)
utils::write.csv(exp1_summary_df, file.path(out_dir, "exp1_summary.csv"), row.names = FALSE)

cat("\n--- Experiment 1 Summary (ConQuest MML vs. BayesMFRM) ---\n")
print(exp1_summary_df[, c("model_id", "prior_mode", "cq_converged", "bayes_converged",
                          "max_rhat", "min_ess_bulk", "param_cor", "param_rmsd_cq_vs_bayes",
                          "person_eap_cor", "fair_score_cor", "infit_rmsd_post", "infit_rmsd_plugin")],
      digits = 4, row.names = FALSE)

# =============================================================================
# EXPERIMENT 2: Sparse Unbalanced Design + Rater x Item Bias (B = B_reps)
# =============================================================================
cat(sprintf("\n=== EXPERIMENT 2: Sparse Unbalanced Bias Benchmark (P = %d, B = %d) ===\n", P_exp2, B_reps))
I2 <- 5L; R2 <- 6L; K2 <- 5L

item_true2 <- c(-0.90, -0.40, 0.00, 0.50, 0.80)
rater_true2 <- c(-0.60, -0.25, -0.05, 0.10, 0.20, 0.60)
tau_true2 <- replicate(I2, c(-1.15, -0.35, 0.35, 1.15), simplify = FALSE)

# Localized sparse bias matrix (exact row & column sum-to-zero):
#   - High-volume rater pair (R01, R02) on I01..I02: +/- 0.70
#   - Low-volume rater pair (R05, R06) on I03..I04: +/- 0.65
bias_true2 <- matrix(0, nrow = I2, ncol = R2)
bias_true2[1, 1] <-  0.70; bias_true2[1, 2] <- -0.70
bias_true2[2, 1] <- -0.70; bias_true2[2, 2] <-  0.70
bias_true2[3, 5] <-  0.65; bias_true2[3, 6] <- -0.65
bias_true2[4, 5] <- -0.65; bias_true2[4, 6] <-  0.65
stopifnot(max(abs(c(rowSums(bias_true2), colSums(bias_true2)))) < 1e-12)

# Connected sparse & unbalanced rating design:
#   - First P_dense = round(P/3) anchor examinees: 2 raters x all 5 items (10 ratings)
#   - Remaining P_sparse examinees: 1 rater x 3 items (3 ratings -> frequent extreme scores!)
#   - Raters R01..R04 handle the bulk; R05 & R06 are low-volume raters (~6-8 ratings per item-rater cell)
build_exp2_sparse_design <- function(P, I, R, seed) {
  set.seed(seed)
  n_dense <- as.integer(round(P / 3))
  rows <- list()
  # Ensure R05 and R06 are bridged cleanly through subsets of anchor examinees
  r_pool_dense <- list(
    c(1L, 2L), c(2L, 3L), c(3L, 4L), c(4L, 1L),
    c(1L, 5L), c(2L, 6L), c(3L, 5L), c(4L, 6L)
  )
  # Weight main raters higher so R05/R06 stay genuinely sparse (~7-8 obs per cell)
  pair_idx <- c(rep(1:4, length.out = max(1L, n_dense - 10L)), rep(5:8, length.out = min(10L, n_dense)))
  for (p in seq_len(n_dense)) {
    rp <- r_pool_dense[[pair_idx[p]]]
    rows[[length(rows) + 1L]] <- expand.grid(
      person = sprintf("P%03d", p),
      item = sprintf("I%02d", seq_len(I)),
      rater = sprintf("R%02d", rp),
      stringsAsFactors = FALSE
    )
  }
  sparse_raters <- sample(c(rep(1:4, each = 18L), 5L, 5L, 6L, 6L), P - n_dense, replace = TRUE)
  for (k in seq_len(P - n_dense)) {
    p <- n_dense + k
    r_one <- sparse_raters[k]
    items_sub <- sort(sample.int(I, 3L, replace = FALSE))
    rows[[length(rows) + 1L]] <- data.frame(
      person = sprintf("P%03d", p),
      item = sprintf("I%02d", items_sub),
      rater = sprintf("R%02d", r_one),
      stringsAsFactors = FALSE
    )
  }
  d <- do.call(rbind, rows)
  rownames(d) <- NULL
  d
}

# Calibrate Horseshoe global scale so induced prior cell RMS matches Gaussian/Student-t (s = 0.50, marginal)
calibrate_hs_scale <- function(I, R, K, seed = 91001L) {
  tpl <- data.frame(
    person = sprintf("P%03d", ((seq_len(max(I, R)) - 1L) %% 2L) + 1L),
    item   = sprintf("I%02d", ((seq_len(max(I, R)) - 1L) %% I) + 1L),
    rater  = sprintf("R%02d", ((seq_len(max(I, R)) - 1L) %% R) + 1L),
    stringsAsFactors = FALSE
  )
  fml <- score ~ person + item + rater + item:rater
  p_iso <- c(
    prior("normal(0, 2)", "mean", facet = "person"),
    prior("lognormal(0, 0.5)", "sd", facet = "person"),
    prior("normal(0, 1)", "item", scale = "marginal"),
    prior("normal(0, 1)", "sd", facet = "rater", scale = "marginal"),
    prior("normal(0, 1)", "tau", scale = "marginal"),
    prior("normal(0, 0.5)", "sd", facet = "item:rater", scale = "marginal")
  )
  pp_iso <- prior_predict_bmfrm(fml, tpl, K = K, hierarchical = c("person", "rater", "item:rater"),
                                priors = p_iso, ndraws = 4000L, seed = seed)
  rms_target <- sqrt(mean(vapply(pp_iso$parameters, function(p) as.vector(p$itemrater), numeric(I * R))^2))
  f_obj <- function(s_g) {
    p_hs <- p_iso
    p_hs[[6L]] <- prior(sprintf("normal(0, %.6f)", s_g), "sd", facet = "item:rater",
                        scale = "raw", local_shrinkage = "horseshoe")
    pp_h <- prior_predict_bmfrm(fml, tpl, K = K, hierarchical = c("person", "rater", "item:rater"),
                                priors = p_hs, ndraws = 3000L, seed = seed)
    sqrt(mean(vapply(pp_h$parameters, function(p) as.vector(p$itemrater), numeric(I * R))^2)) - rms_target
  }
  round(stats::uniroot(f_obj, interval = c(0.05, 1.0), tol = 1e-3)$root, 4L)
}

hs_s_cal <- calibrate_hs_scale(I2, R2, K2, seed = 91001L)
cat(sprintf("Calibrated Horseshoe global scale for %dx%d interaction (target prior cell RMS ~0.50): %.4f\n",
            I2, R2, hs_s_cal))

fml_exp2 <- score ~ person + item + rater + item:rater
hier_exp2 <- c("person", "rater", "item:rater")
base_priors_exp2 <- c(
  prior("normal(0, 2)", "mean", facet = "person"),
  prior("lognormal(0, 0.5)", "sd", facet = "person"),
  prior("normal(0, 1)", "item", scale = "marginal"),
  prior("normal(0, 1)", "sd", facet = "rater", scale = "marginal"),
  prior("normal(0, 1)", "tau", scale = "marginal")
)
bayes_methods_exp2 <- list(
  Bayes_Gaussian = list(
    priors = c(base_priors_exp2, list(prior("normal(0, 0.5)", "sd", facet = "item:rater", scale = "marginal", local_shrinkage = "none"))),
    adapt_delta = 0.96, seed_offset = 1000L
  ),
  Bayes_StudentT = list(
    priors = c(base_priors_exp2, list(prior("normal(0, 0.5)", "sd", facet = "item:rater", scale = "marginal", local_shrinkage = "student_t"))),
    adapt_delta = 0.995, seed_offset = 2000L
  ),
  Bayes_Horseshoe = list(
    priors = c(base_priors_exp2, list(prior(sprintf("normal(0, %.4f)", hs_s_cal), "sd", facet = "item:rater", scale = "raw", local_shrinkage = "horseshoe"))),
    adapt_delta = 0.98, seed_offset = 3000L
  )
)

exp2_rep_rows <- list()
exp2_rep1_cells <- NULL

for (b in seq_len(B_reps)) {
  cat(sprintf("[Exp 2 | Replication %02d / %02d] ", b, B_reps))
  des_b <- build_exp2_sparse_design(P_exp2, I2, R2, seed = 20260500L + b)
  dc_b <- check_design_bmfrm(fml_exp2, des_b, K = K2)
  # Ensure full structural rank for fair unpenalized ConQuest comparison; retry deterministically if needed
  attempt <- 0L
  while (!dc_b$identifiable && attempt < 20L) {
    attempt <- attempt + 1L
    des_b <- build_exp2_sparse_design(P_exp2, I2, R2, seed = 20260500L + b * 100L + attempt)
    dc_b <- check_design_bmfrm(fml_exp2, des_b, K = K2)
  }
  stopifnot(dc_b$identifiable)

  set.seed(20261127L + b)
  theta_b <- stats::rnorm(P_exp2, mean = 0, sd = 1.0)
  dat_b <- simulate_mfrm_dataset(des_b, theta_b, item_true2, rater_true2, tau_true2, bias_true2, seed = 20262127L + b)
  true_fair_b <- compute_point_fair_scores(theta_b, item_true2, tau_true2)

  # Cell observation counts & examinee subgroups
  i_id_b <- as.integer(sub("^I", "", dat_b$item))
  r_id_b <- as.integer(sub("^R", "", dat_b$rater))
  cell_counts <- matrix(tabulate(i_id_b + (r_id_b - 1L) * I2, I2 * R2), nrow = I2, ncol = R2)
  low_cell_mask <- cell_counts <= 10L
  nonzero_mask <- abs(bias_true2) > 1e-9

  n_per_p <- as.integer(table(factor(dat_b$person, levels = sprintf("P%03d", seq_len(P_exp2)))))
  is_sparse_p <- n_per_p <= 3L
  is_extreme_p <- vapply(split(dat_b$score, factor(dat_b$person, levels = sprintf("P%03d", seq_len(P_exp2)))),
                         function(y) all(y == 1L) || all(y == K2), logical(1))

  score_method_metrics <- function(method_name, converged_ok, empirical_se_fail, elapsed_s,
                                   max_rhat, min_ess, n_div,
                                   item_est, item_lo, item_hi,
                                   rater_est, rater_lo, rater_hi,
                                   step_est, bias_est, bias_lo, bias_hi,
                                   th_est, th_lo, th_hi, fair_est) {
    cov_fn <- function(tr, lo, hi) {
      ok <- is.finite(lo) & is.finite(hi)
      if (!any(ok)) return(NA_real_)
      mean(tr[ok] >= lo[ok] & tr[ok] <= hi[ok])
    }
    data.frame(
      replication = b,
      method = method_name,
      converged = converged_ok,
      empirical_se_failed = empirical_se_fail,
      elapsed_sec = elapsed_s,
      max_rhat = max_rhat,
      min_ess_bulk = min_ess,
      divergences = n_div,
      n_ratings = nrow(dat_b),
      n_low_cells = sum(low_cell_mask),
      n_extreme_persons = sum(is_extreme_p),
      bias_rmse_all = sqrt(mean((bias_est - bias_true2)^2)),
      bias_rmse_low_cells = sqrt(mean((bias_est[low_cell_mask] - bias_true2[low_cell_mask])^2)),
      bias_rmse_null_cells = sqrt(mean((bias_est[!nonzero_mask] - bias_true2[!nonzero_mask])^2)),
      bias_rmse_active_cells = sqrt(mean((bias_est[nonzero_mask] - bias_true2[nonzero_mask])^2)),
      bias_max_abs_err = max(abs(bias_est - bias_true2)),
      rater_rmse = sqrt(mean((rater_est - rater_true2)^2)),
      item_rmse = sqrt(mean((item_est - item_true2)^2)),
      step_rmse = sqrt(mean((step_est - tau_true2[[1L]])^2)),
      theta_rmse_all = sqrt(mean((th_est - theta_b)^2)),
      theta_rmse_sparse3 = sqrt(mean((th_est[is_sparse_p] - theta_b[is_sparse_p])^2)),
      theta_rmse_dense10 = sqrt(mean((th_est[!is_sparse_p] - theta_b[!is_sparse_p])^2)),
      theta_rmse_extreme = if (any(is_extreme_p)) sqrt(mean((th_est[is_extreme_p] - theta_b[is_extreme_p])^2)) else NA_real_,
      fair_rmse_all = sqrt(mean((fair_est - true_fair_b)^2)),
      fair_rmse_sparse3 = sqrt(mean((fair_est[is_sparse_p] - true_fair_b[is_sparse_p])^2)),
      bias_cov95 = cov_fn(as.vector(bias_true2), as.vector(bias_lo), as.vector(bias_hi)),
      theta_cov95 = cov_fn(theta_b, th_lo, th_hi),
      bias_ci_width = mean(bias_hi - bias_lo, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }

  # 1. ACER ConQuest MML
  cq_b <- run_conquest_mfrm(
    dat_b, P = P_exp2, I = I2, R = R2, K_vec = K2,
    model_type = "rsm_bias", cq_seed = 200L + b
  )
  cat(sprintf("ConQuest(conv=%s,se=%s) ", cq_b$converged, cq_b$se_mode))
  exp2_rep_rows[[length(exp2_rep_rows) + 1L]] <- score_method_metrics(
    method_name = "ConQuest_MML",
    converged_ok = cq_b$converged,
    empirical_se_fail = cq_b$empirical_se_failed,
    elapsed_s = cq_b$elapsed_sec,
    max_rhat = NA_real_, min_ess = NA_real_, n_div = NA_integer_,
    item_est = cq_b$item$est,
    item_lo = cq_b$item$est - 1.96 * cq_b$item$se,
    item_hi = cq_b$item$est + 1.96 * cq_b$item$se,
    rater_est = cq_b$rater$est,
    rater_lo = cq_b$rater$est - 1.96 * cq_b$rater$se,
    rater_hi = cq_b$rater$est + 1.96 * cq_b$rater$se,
    step_est = cq_b$tau[[1L]],
    bias_est = cq_b$bias,
    bias_lo = cq_b$bias - 1.96 * cq_b$bias_se,
    bias_hi = cq_b$bias + 1.96 * cq_b$bias_se,
    th_est = cq_b$eap,
    th_lo = cq_b$eap - 1.96 * cq_b$eap_sd,
    th_hi = cq_b$eap + 1.96 * cq_b$eap_sd,
    fair_est = cq_b$fair_score
  )

  rep1_list <- list()
  if (b == 1L) {
    rep1_list[["ConQuest_MML"]] <- data.frame(
      method = "ConQuest_MML",
      cell = sprintf("I%02d:R%02d", rep(1:I2, R2), rep(1:R2, each = I2)),
      n_obs = as.vector(cell_counts),
      truth = as.vector(bias_true2),
      est = as.vector(cq_b$bias),
      lo = as.vector(cq_b$bias - 1.96 * cq_b$bias_se),
      hi = as.vector(cq_b$bias + 1.96 * cq_b$bias_se),
      stringsAsFactors = FALSE
    )
  }

  # 2. BayesMFRM Shrinkage Suite
  for (m_nm in names(bayes_methods_exp2)) {
    m_cfg <- bayes_methods_exp2[[m_nm]]
    t0 <- proc.time()[["elapsed"]]
    bf_b <- suppressWarnings(suppressMessages(bmfrm(
      formula = fml_exp2,
      data = dat_b,
      K = K2,
      family = "rsm",
      hierarchical = hier_exp2,
      priors = m_cfg$priors,
      save_gq = "minimal",
      cache_dir = stan_cache_dir,
      model_name = paste0("exp2_", m_nm),
      design_check = "none",
      chains = 4L, cores = 4L, iter = 2000L, warmup = 1000L,
      adapt_delta = m_cfg$adapt_delta, max_treedepth = 12L,
      seed = 81000L + m_cfg$seed_offset + b,
      refresh = 0L, show_messages = FALSE, show_exceptions = FALSE,
      cpp_options = list(PRECOMPILED_HEADERS = "false")
    )))
    el_b <- proc.time()[["elapsed"]] - t0
    br_b <- extract_bayes_estimates(bf_b, P = P_exp2, I = I2, R = R2, K_vec = K2)
    cat(sprintf("%s(conv=%s,Rhat=%.3f) ", m_nm, br_b$converged, br_b$max_rhat))

    exp2_rep_rows[[length(exp2_rep_rows) + 1L]] <- score_method_metrics(
      method_name = m_nm,
      converged_ok = br_b$converged,
      empirical_se_fail = FALSE,
      elapsed_s = el_b,
      max_rhat = br_b$max_rhat,
      min_ess = br_b$min_ess_bulk,
      n_div = br_b$n_divergent,
      item_est = br_b$item, item_lo = br_b$item_lo, item_hi = br_b$item_hi,
      rater_est = br_b$rater, rater_lo = br_b$rater_lo, rater_hi = br_b$rater_hi,
      step_est = br_b$tau[[1L]],
      bias_est = br_b$bias, bias_lo = br_b$bias_lo, bias_hi = br_b$bias_hi,
      th_est = br_b$theta, th_lo = br_b$theta_lo, th_hi = br_b$theta_hi,
      fair_est = br_b$fair_score
    )
    if (b == 1L) {
      rep1_list[[m_nm]] <- data.frame(
        method = m_nm,
        cell = sprintf("I%02d:R%02d", rep(1:I2, R2), rep(1:R2, each = I2)),
        n_obs = as.vector(cell_counts),
        truth = as.vector(bias_true2),
        est = as.vector(br_b$bias),
        lo = as.vector(br_b$bias_lo),
        hi = as.vector(br_b$bias_hi),
        stringsAsFactors = FALSE
      )
    }
    unlink(bf_b$fit$output_files())
  }
  if (b == 1L) exp2_rep1_cells <- do.call(rbind, rep1_list)
  cat("\n")
}

exp2_reps_df <- do.call(rbind, exp2_rep_rows)
utils::write.csv(exp2_reps_df, file.path(out_dir, "exp2_replications.csv"), row.names = FALSE)
utils::write.csv(exp2_rep1_cells, file.path(out_dir, "exp2_rep1_cells.csv"), row.names = FALSE)

# Aggregate strictly over converged runs AND save all-completed table separately
summarise_exp2 <- function(df_in, converged_only = TRUE) {
  methods <- c("ConQuest_MML", names(bayes_methods_exp2))
  do.call(rbind, lapply(methods, function(m) {
    all_m <- df_in[df_in$method == m, , drop = FALSE]
    sub_m <- if (converged_only) all_m[all_m$converged, , drop = FALSE] else all_m
    mcse_fn <- function(v) {
      v <- v[is.finite(v)]
      if (length(v) <= 1L) return(NA_real_)
      stats::sd(v) / sqrt(length(v))
    }
    data.frame(
      method = m,
      n_total = nrow(all_m),
      n_converged = sum(all_m$converged),
      n_empirical_se_failed = sum(all_m$empirical_se_failed),
      mean_sec = mean(sub_m$elapsed_sec),
      bias_rmse_all = mean(sub_m$bias_rmse_all),
      bias_rmse_all_mcse = mcse_fn(sub_m$bias_rmse_all),
      bias_rmse_low_cells = mean(sub_m$bias_rmse_low_cells),
      bias_rmse_low_mcse = mcse_fn(sub_m$bias_rmse_low_cells),
      bias_rmse_null_cells = mean(sub_m$bias_rmse_null_cells),
      bias_rmse_active_cells = mean(sub_m$bias_rmse_active_cells),
      bias_max_abs_err = mean(sub_m$bias_max_abs_err),
      rater_rmse = mean(sub_m$rater_rmse),
      item_rmse = mean(sub_m$item_rmse),
      theta_rmse_all = mean(sub_m$theta_rmse_all),
      theta_rmse_sparse3 = mean(sub_m$theta_rmse_sparse3),
      theta_rmse_extreme = mean(sub_m$theta_rmse_extreme, na.rm = TRUE),
      fair_rmse_all = mean(sub_m$fair_rmse_all),
      fair_rmse_sparse3 = mean(sub_m$fair_rmse_sparse3),
      bias_cov95 = mean(sub_m$bias_cov95, na.rm = TRUE),
      theta_cov95 = mean(sub_m$theta_cov95, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
}

exp2_summary_conv <- summarise_exp2(exp2_reps_df, converged_only = TRUE)
exp2_summary_all <- summarise_exp2(exp2_reps_df, converged_only = FALSE)
utils::write.csv(exp2_summary_conv, file.path(out_dir, "exp2_summary_converged.csv"), row.names = FALSE)
utils::write.csv(exp2_summary_all, file.path(out_dir, "exp2_summary_all.csv"), row.names = FALSE)

cat("\n--- Experiment 2 Summary (Converged Runs Only; note Null-Cell Noise Reduction vs. Active-Cell Shrinkage Attenuation) ---\n")
print(exp2_summary_conv[, c("method", "n_converged", "n_empirical_se_failed",
                            "bias_rmse_all", "bias_rmse_null_cells", "bias_rmse_active_cells",
                            "bias_rmse_low_cells", "rater_rmse",
                            "theta_rmse_sparse3", "fair_rmse_sparse3", "bias_cov95", "theta_cov95")],
      digits = 4, row.names = FALSE)

# =============================================================================
# 3. Publication-Grade Visualizations (PNG + PDF)
# =============================================================================
render_dual <- function(stem, width_in, height_in, plot_fn) {
  grDevices::png(file.path(out_dir, paste0(stem, ".png")),
                 width = width_in * 180, height = height_in * 180, res = 180)
  plot_fn()
  grDevices::dev.off()
  grDevices::pdf(file.path(out_dir, paste0(stem, ".pdf")),
                 width = width_in, height = height_in)
  plot_fn()
  grDevices::dev.off()
}

# Figure 1: 4-Model Parameter Parity (ConQuest MML vs BayesMFRM 95% CrI)
render_dual("fig1_exp1_parameter_parity", 10.5, 8.5, function() {
  graphics::par(mfrow = c(2, 2), mar = c(4.2, 4.4, 3.0, 1.2), oma = c(0, 0, 2.2, 0))
  cols_map <- c(item = "#1b6ca8", rater = "#d95f02", step = "#2a9d8f", interaction = "#7b2cbf")
  pch_map <- c(item = 16L, rater = 17L, step = 15L, interaction = 18L)

  for (m_id in names(exp1_specs)) {
    sub_w <- exp1_params_df[exp1_params_df$model_id == m_id & exp1_params_df$prior_mode == "weak_marginal", ]
    lim <- range(c(sub_w$cq_est, sub_w$bayes_lo, sub_w$bayes_hi), na.rm = TRUE)
    pad <- 0.08 * diff(lim)
    lim <- c(lim[1L] - pad, lim[2L] + pad)

    graphics::plot(NA, xlim = lim, ylim = lim,
                   xlab = "ACER ConQuest 5.47.5 Estimate (logits)",
                   ylab = "BayesMFRM Posterior Mean & 95% CrI (logits)",
                   main = exp1_specs[[m_id]]$label, las = 1, cex.main = 1.05)
    graphics::grid(col = "#e6e6e6", lty = 1)
    graphics::abline(0, 1, col = "#444444", lty = 2, lwd = 1.3)
    graphics::segments(sub_w$cq_est, sub_w$bayes_lo, sub_w$cq_est, sub_w$bayes_hi,
                       col = grDevices::adjustcolor(cols_map[sub_w$param_class], alpha.f = 0.45), lwd = 1.4)
    graphics::points(sub_w$cq_est, sub_w$bayes_mean,
                     col = cols_map[sub_w$param_class], pch = pch_map[sub_w$param_class], cex = 1.15)

    r_val <- stats::cor(sub_w$cq_est, sub_w$bayes_mean)
    rmsd_val <- sqrt(mean((sub_w$cq_est - sub_w$bayes_mean)^2))
    graphics::legend("topleft",
                     legend = c(sprintf("r = %.4f", r_val), sprintf("RMSD = %.3f", rmsd_val)),
                     bty = "n", cex = 0.88, text.col = "#222222")
    present_cls <- unique(sub_w$param_class)
    graphics::legend("bottomright", legend = present_cls,
                     col = cols_map[present_cls], pch = pch_map[present_cls],
                     bty = "o", bg = "white", box.col = "#cccccc", cex = 0.82)
  }
  graphics::mtext(sprintf("Experiment 1: Cross-Software Parameter Parity (P = %d, Weakly-Informative Marginal Priors)", P_exp1),
                  outer = TRUE, cex = 1.15, font = 2, line = 0.5)
})

# Figure 2: Person EAP, Uncertainty, & Fit-Statistic Parity
render_dual("fig2_exp1_eap_and_fit_parity", 11.5, 4.0, function() {
  graphics::par(mfrow = c(1, 3), mar = c(4.3, 4.3, 3.0, 1.0), oma = c(0, 0, 2.0, 0))
  m_cols <- stats::setNames(c("#1b6ca8", "#2a9d8f", "#e76f51", "#7b2cbf"), names(exp1_specs))

  # Panel A: Person EAP vs BayesMFRM theta posterior mean
  lim_a <- range(c(exp1_persons_df$cq_eap, exp1_persons_df$bayes_theta), na.rm = TRUE)
  graphics::plot(exp1_persons_df$cq_eap, exp1_persons_df$bayes_theta,
                 col = grDevices::adjustcolor(m_cols[exp1_persons_df$model_id], alpha.f = 0.45),
                 pch = 16, cex = 0.75, xlim = lim_a, ylim = lim_a, las = 1,
                 xlab = "ConQuest Person EAP (logits)",
                 ylab = "BayesMFRM Person Posterior Mean (logits)",
                 main = "A. Examinee Ability (EAP vs. Bayes)")
  graphics::grid(col = "#e6e6e6", lty = 1)
  graphics::abline(0, 1, lty = 2, lwd = 1.3, col = "#333333")
  graphics::legend("topleft",
                   legend = sprintf("Overall r = %.4f\nRMSD = %.3f",
                                    stats::cor(exp1_persons_df$cq_eap, exp1_persons_df$bayes_theta),
                                    sqrt(mean((exp1_persons_df$cq_eap - exp1_persons_df$bayes_theta)^2))),
                   bty = "n", cex = 0.88)

  # Panel B: Standardized Fair Scores (1..K rubric scale)
  lim_b <- range(c(exp1_persons_df$cq_fair, exp1_persons_df$bayes_fair), na.rm = TRUE)
  graphics::plot(exp1_persons_df$cq_fair, exp1_persons_df$bayes_fair,
                 col = grDevices::adjustcolor(m_cols[exp1_persons_df$model_id], alpha.f = 0.45),
                 pch = 16, cex = 0.75, xlim = lim_b, ylim = lim_b, las = 1,
                 xlab = "ConQuest Standardized Fair Score",
                 ylab = "BayesMFRM fair_scores(standardize = TRUE)",
                 main = "B. Standardized Fair Scores")
  graphics::grid(col = "#e6e6e6", lty = 1)
  graphics::abline(0, 1, lty = 2, lwd = 1.3, col = "#333333")
  graphics::legend("topleft",
                   legend = c("M1 RSM", "M2 PCM", "M3 Hetero PCM", "M4 Bias RSM"),
                   col = m_cols, pch = 16, bty = "o", bg = "white", box.col = "#cccccc", cex = 0.78)

  # Panel C: Weighted & Unweighted MNSQ (Posterior-averaged vs. Plug-in EAP shrinkage)
  lim_c <- range(c(exp1_fitstats_df$cq_infit, exp1_fitstats_df$bayes_infit_post,
                   exp1_fitstats_df$bayes_infit_plugin, exp1_fitstats_df$cq_outfit,
                   exp1_fitstats_df$bayes_outfit_post), na.rm = TRUE)
  pad_c <- max(0.04, 0.12 * diff(lim_c))
  lim_c <- c(lim_c[1L] - pad_c, lim_c[2L] + pad_c)
  graphics::plot(NA, xlim = lim_c, ylim = lim_c, las = 1,
                 xlab = "ConQuest PV-Averaged MNSQ (50 PVs)",
                 ylab = "BayesMFRM facet_fit() MNSQ",
                 main = "C. Infit & Outfit MNSQ (Posterior vs. Plug-in)")
  graphics::grid(col = "#e6e6e6", lty = 1)
  graphics::abline(0, 1, lty = 2, lwd = 1.3, col = "#333333")
  graphics::points(exp1_fitstats_df$cq_infit, exp1_fitstats_df$bayes_infit_plugin,
                   col = "#8d99ae", pch = 1, cex = 1.0, lwd = 1.3)
  graphics::points(exp1_fitstats_df$cq_infit, exp1_fitstats_df$bayes_infit_post,
                   col = "#1b6ca8", pch = 16, cex = 1.1)
  graphics::points(exp1_fitstats_df$cq_outfit, exp1_fitstats_df$bayes_outfit_post,
                   col = "#d95f02", pch = 17, cex = 1.1)
  rmsd_in_post <- sqrt(mean((exp1_fitstats_df$cq_infit - exp1_fitstats_df$bayes_infit_post)^2))
  rmsd_out_post <- sqrt(mean((exp1_fitstats_df$cq_outfit - exp1_fitstats_df$bayes_outfit_post)^2))
  rmsd_in_plug <- sqrt(mean((exp1_fitstats_df$cq_infit - exp1_fitstats_df$bayes_infit_plugin)^2))
  graphics::legend("topleft",
                   legend = c(sprintf("Posterior Infit (RMSD = %.3f)", rmsd_in_post),
                              sprintf("Posterior Outfit (RMSD = %.3f)", rmsd_out_post),
                              sprintf("Plug-in Infit (RMSD = %.3f)", rmsd_in_plug)),
                   col = c("#1b6ca8", "#d95f02", "#8d99ae"), pch = c(16, 17, 1),
                   bty = "o", bg = "white", box.col = "#cccccc", cex = 0.78)
  graphics::mtext("Experiment 1: Person EAP, Fair Score, and Plausible-Value / Posterior MNSQ Agreement",
                  outer = TRUE, cex = 1.12, font = 2, line = 0.4)
})

# Figure 3: Experiment 2 Sparse & Unbalanced Bias Benchmark
render_dual("fig3_exp2_sparse_bias_benchmark", 12.0, 4.3, function() {
  graphics::par(mfrow = c(1, 3), mar = c(5.2, 4.3, 3.0, 1.0), oma = c(0, 0, 2.0, 0))
  m_order <- c("ConQuest_MML", "Bayes_Gaussian", "Bayes_StudentT", "Bayes_Horseshoe")
  m_short <- c("ConQuest\nMML", "Bayes\nGaussian", "Bayes\nStudent-t", "Bayes\nHorseshoe")
  m_pal <- c("#8d99ae", "#1b6ca8", "#2a9d8f", "#7b2cbf")

  sub_conv <- exp2_summary_conv[match(m_order, exp2_summary_conv$method), ]

  # Panel A: Interaction Bias RMSE (All cells vs. Low-count sparse cells)
  mat_a <- rbind(sub_conv$bias_rmse_all, sub_conv$bias_rmse_low_cells)
  bp <- graphics::barplot(
    mat_a, beside = TRUE, names.arg = m_short,
    col = c("#3a86ff", "#fb5607"), border = NA, las = 1,
    ylim = c(0, max(mat_a, na.rm = TRUE) * 1.28),
    ylab = "Interaction Bias RMSE (logits, lower is better)",
    main = "A. Bias Recovery (All vs. Low-Count Cells)"
  )
  graphics::grid(nx = NA, ny = NULL, col = "#e6e6e6", lty = 1)
  graphics::barplot(mat_a, beside = TRUE, names.arg = m_short,
                    col = c("#3a86ff", "#fb5607"), border = NA, las = 1, add = TRUE)
  se_a <- rbind(sub_conv$bias_rmse_all_mcse, sub_conv$bias_rmse_low_mcse)
  ok_se <- is.finite(se_a) & se_a > 0
  graphics::segments(bp[ok_se], mat_a[ok_se] - se_a[ok_se], bp[ok_se], mat_a[ok_se] + se_a[ok_se], lwd = 1.5)
  graphics::legend("topright", legend = c("All 30 cells", "Low-count cells (n <= 10)"),
                   fill = c("#3a86ff", "#fb5607"), border = NA, bty = "o", bg = "white", box.col = "#cccccc", cex = 0.80)

  # Panel B: Replication 1 Shrinkage Profile (ConQuest MML vs Bayes Student-t)
  c_cq <- exp2_rep1_cells[exp2_rep1_cells$method == "ConQuest_MML", ]
  c_st <- exp2_rep1_cells[exp2_rep1_cells$method == "Bayes_StudentT", ]
  ord_c <- order(c_cq$truth, c_cq$est)
  idx_x <- seq_along(ord_c)
  lim_y <- range(c(c_cq$est, c_st$lo, c_st$hi, c_cq$truth), na.rm = TRUE)
  graphics::plot(idx_x, c_cq$truth[ord_c], type = "n", ylim = lim_y, las = 1,
                 xlab = "Interaction Cells (sorted by true bias)",
                 ylab = "Estimated Bias Effect (logits)",
                 main = "B. Cell-Level Noise vs. Shrinkage (Rep 1)")
  graphics::grid(col = "#e6e6e6", lty = 1)
  graphics::abline(h = 0, col = "#888888", lty = 3)
  graphics::segments(idx_x, c_st$lo[ord_c], idx_x, c_st$hi[ord_c],
                     col = grDevices::adjustcolor("#2a9d8f", 0.35), lwd = 1.8)
  graphics::points(idx_x, c_cq$est[ord_c], col = "#d90429", pch = 4, lwd = 1.5, cex = 0.95)
  graphics::points(idx_x, c_st$est[ord_c], col = "#2a9d8f", pch = 16, cex = 0.95)
  graphics::points(idx_x, c_cq$truth[ord_c], col = "#111111", pch = 1, lwd = 1.4, cex = 1.05)
  graphics::legend("topleft",
                   legend = c("Ground truth", "ConQuest MML (unpenalized)", "Bayes Student-t + 95% CrI"),
                   col = c("#111111", "#d90429", "#2a9d8f"), pch = c(1, 4, 16),
                   bty = "o", bg = "white", box.col = "#cccccc", cex = 0.78)

  # Panel C: Sparse-Examinee Fair-Score RMSE & Nominal 95% Bias Coverage
  mat_c <- rbind(sub_conv$fair_rmse_all, sub_conv$fair_rmse_sparse3)
  graphics::barplot(
    mat_c, beside = TRUE, names.arg = m_short,
    col = c("#457b9d", "#e76f51"), border = NA, las = 1,
    ylim = c(0, max(mat_c, na.rm = TRUE) * 1.32),
    ylab = "Standardized Fair-Score RMSE (1..5 scale)",
    main = "C. Fair-Score Error & 95% Bias Coverage"
  )
  graphics::grid(nx = NA, ny = NULL, col = "#e6e6e6", lty = 1)
  bp_c <- graphics::barplot(mat_c, beside = TRUE, names.arg = m_short,
                            col = c("#457b9d", "#e76f51"), border = NA, las = 1, add = TRUE)
  cov_lbl <- sprintf("Cov:%.0f%%", 100 * sub_conv$bias_cov95)
  graphics::text(colMeans(bp_c), apply(mat_c, 2L, max) + 0.03 * max(mat_c),
                 labels = cov_lbl, cex = 0.78, font = 2, col = "#222222")
  graphics::legend("topright", legend = c("All examinees", "3-rating sparse examinees"),
                   fill = c("#457b9d", "#e76f51"), border = NA, bty = "o", bg = "white", box.col = "#cccccc", cex = 0.78)

  graphics::mtext(sprintf("Experiment 2: Sparse Unbalanced Design Benchmark (%d Replications, P = %d, I = %d, R = %d)",
                          B_reps, P_exp2, I2, R2),
                  outer = TRUE, cex = 1.12, font = 2, line = 0.4)
})

cat("\nSaved all CSV summaries, provenance manifest, and figures (PNG + PDF) to:", out_dir, "\n")
