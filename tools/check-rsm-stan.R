# Deterministic density/analytic-gradient parity against the compiled Stan target.
# This command compiles one model and evaluates fixed points; it runs no MCMC.
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1L])
root <- dirname(dirname(normalizePath(script)))
suppressPackageStartupMessages(pkgload::load_all(root, quiet = TRUE))
out <- file.path(root, ".cache", "isometric")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

# Independent closed-form Helmert construction, category-index softmax, and
# analytical derivative. Deliberately do not call the package probability/basis code.
basis <- function(n, scale) {
  ans <- matrix(0, n, max(n - 1L, 0L))
  if (n > 1L) for (j in seq_len(n - 1L)) {
    ans[seq_len(j), j] <- -1 / sqrt(j * (j + 1))
    ans[j + 1L, j] <- j / sqrt(j * (j + 1))
  }
  if (n > 1L && scale == "marginal") ans <- ans * sqrt(n / (n - 1L))
  ans
}
reference <- function(q, d, scales) {
  pos <- 0L
  take <- function(n) {
    idx <- pos + seq_len(n)
    pos <<- pos + n
    idx
  }
  iz <- take(d$J_person); im <- take(d$estimate_mu); il <- take(d$estimate_sigma)
  facets <- setdiff(sub("^J_", "", grep("^J_", names(d), value = TRUE)), "person")
  indices <- lapply(facets, function(f) take(d[[paste0("J_", f)]] - 1L))
  it <- take(d$K - 2L)
  stopifnot(pos == length(q))
  mu <- if (length(im)) d$ability_mean_location + d$ability_mean_scale * q[im] else d$ability_mean_location
  sigma <- if (length(il)) exp(d$ability_sd_location + d$ability_sd_scale * q[il]) else d$ability_sd_location
  bases <- lapply(seq_along(facets), function(j) basis(d[[paste0("J_", facets[j])]], scales[j]))
  Bt <- basis(d$K - 1L, scales[length(scales)])
  tau <- as.vector(d$prior_sd_tau * Bt %*% q[it])
  eta <- mu + sigma * q[iz][d$person_id]
  for (j in seq_along(facets)) {
    f <- facets[j]
    effect <- as.vector(d[[paste0("prior_sd_", f)]] * bases[[j]] %*% q[indices[[j]]])
    eta <- eta - effect[d[[paste0(f, "_id")]]]
  }
  k <- 0:(d$K - 1L)
  logits <- outer(eta, k) - matrix(c(0, cumsum(tau)), d$N, d$K, byrow = TRUE)
  logits <- logits - apply(logits, 1L, max)
  pr <- exp(logits) / rowSums(exp(logits))
  ll <- logits[cbind(seq_len(d$N), d$X)] - log(rowSums(exp(logits)))
  score <- d$X - 1L - as.vector(pr %*% k)
  group_sum <- function(ids, n) vapply(seq_len(n), function(j) sum(score[ids == j]), numeric(1))
  g <- -q
  g[iz] <- g[iz] + sigma * group_sum(d$person_id, d$J_person)
  if (length(im)) g[im] <- g[im] + d$ability_mean_scale * sum(score)
  if (length(il)) g[il] <- g[il] + d$ability_sd_scale * sigma * sum(q[iz][d$person_id] * score)
  for (j in seq_along(facets)) {
    f <- facets[j]; idx <- indices[[j]]
    g[idx] <- g[idx] - as.vector(d[[paste0("prior_sd_", f)]] * crossprod(bases[[j]],
      group_sum(d[[paste0(f, "_id")]], d[[paste0("J_", f)]])))
  }
  gt <- vapply(seq_len(d$K - 1L), function(h) sum(rowSums(pr[, (h + 1L):d$K, drop = FALSE]) - (d$X - 1L >= h)), numeric(1))
  g[it] <- g[it] + as.vector(d$prior_sd_tau * crossprod(Bt, gt))
  list(lp = sum(ll) - sum(q^2) / 2, gradient = g)
}

cases <- list(
  hierarchical = list(K = 4L, items = 3L, raters = 2L, sds = c(.8, 1.3, .6)),
  fixed_binary = list(K = 2L, items = 1L, raters = 1L, mean = 1.2, sd = .7),
  mean_only = list(K = 3L, items = 2L, raters = 1L, sd = .9),
  sd_only = list(K = 5L, items = 1L, raters = 3L, mean = -.4),
  raw_hierarchical = list(K = 4L, items = 3L, raters = 2L, scales = rep("raw", 3)),
  mixed_independent = list(K = 4L, items = 3L, raters = 2L, hierarchical = character(),
                           scales = c("marginal", "raw", "marginal")),
  ten_facets = list(K = 4L, items = 3L, raters = 2L, extras = paste0("f", 4:10)))
rows <- list()
for (name in names(cases)) {
  cs <- cases[[name]]
  hierarchical <- cs[["hierarchical"]] %||% "person"
  facets <- c("item", "rater", cs[["extras"]])
  scales <- cs[["scales"]] %||% rep("marginal", length(facets) + 1L)
  sds <- cs[["sds"]] %||% rep(1, length(scales))
  ps <- if ("person" %in% hierarchical) list(prior(cs[["mean"]] %||% "normal(0, 2.5)", "mean", facet = "person"),
    prior(cs[["sd"]] %||% "lognormal(0, 0.5)", "sd", facet = "person")) else list(prior("normal(.4, .8)", "theta"))
  classes <- c(ifelse(facets %in% c("item", "rater"), facets, paste0("facet_", facets)), "tau")
  for (j in seq_along(classes)) ps <- c(ps, list(prior(sprintf("normal(0, %.17g)", sds[j]), classes[j], scale = scales[j])))
  dat <- expand.grid(person = 1:12, item = seq_len(cs[["items"]]), rater = seq_len(cs[["raters"]]))
  for (j in seq_along(cs[["extras"]])) dat[[cs[["extras"]][j]]] <- (seq_len(nrow(dat)) %/% 2^(j - 1L)) %% 2L + 1L
  dat$score <- (seq_len(nrow(dat)) - 1L) %% cs[["K"]] + 1L
  spec <- parse_bmfrm_formula(stats::reformulate(c("person", facets), response = "score"), dat, K = cs[["K"]], hierarchical = hierarchical)
  d <- prepare_data_bmfrm(spec, dat, ps)$stan_data
  compiled <- compile_bmfrm(build_stan_code(spec, ps), "isometric", out,
                           cpp_options = list(PRECOMPILED_HEADERS = "false"))
  dim <- d$J_person + d$estimate_mu + d$estimate_sigma + sum(unlist(d[paste0("J_", facets)]) - 1L) + d$K - 2L
  points <- rbind(rep(0, dim), .3 * sin(seq_len(dim)), .7 * cos(seq_len(dim)), seq(-.8, .8, length.out = dim))
  data_path <- file.path(out, paste0(name, "-data.json"))
  cmdstanr::write_stan_json(d, data_path)
  # One flat vector per invocation avoids a multi-point JSON memory-order ambiguity.
  actual <- do.call(rbind, lapply(seq_len(nrow(points)), function(j) {
    point_path <- file.path(out, paste0(name, "-point-", j, ".json"))
    result_path <- file.path(out, paste0(name, "-point-", j, ".csv"))
    cmdstanr::write_stan_json(list(params_r = as.numeric(points[j, ])), point_path)
    status <- system2(compiled$model$exe_file(), c("log_prob", paste0("unconstrained_params=", shQuote(point_path)),
      "data", paste0("file=", shQuote(data_path)), "output", paste0("file=", shQuote(result_path)), "sig_figs=18"),
      stdout = file.path(out, paste0(name, "-point-", j, ".log")), stderr = "")
    stopifnot(status == 0L)
    utils::read.csv(result_path, comment.char = "#", check.names = FALSE)
  }))
  expected <- lapply(seq_len(nrow(points)), function(j) reference(points[j, ], d, scales))
  lp <- vapply(expected, `[[`, numeric(1), "lp")
  # CmdStan can omit parameter-independent normalizing constants.
  lp_error <- max(abs((actual[[1L]] - actual[[1L]][1L]) - (lp - lp[1L])))
  grad_error <- max(abs(as.matrix(actual[, -1L]) - do.call(rbind, lapply(expected, `[[`, "gradient"))))
  cat(name, "density difference error:", lp_error, "gradient error:", grad_error, "\n")
  stopifnot(lp_error < 1e-10, grad_error < 1e-10)
  rows[[name]] <- data.frame(case = name, points = nrow(points), dimensions = dim,
                             density_difference_error = lp_error, gradient_error = grad_error)
}
result <- do.call(rbind, rows)
utils::write.csv(result, file.path(out, "parity.csv"), row.names = FALSE)
print(result, row.names = FALSE)
