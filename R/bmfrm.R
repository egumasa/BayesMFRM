#' Fit a Bayesian Many-Facet Rasch Model
#'
#' Fits Rating Scale Models (RSM) or Partial Credit Models (PCM) with arbitrary
#' main facets, hierarchical partial pooling, and two-way interactions (including
#' Student-t and regularized horseshoe local shrinkage) via CmdStan.
#'
#' @param formula R formula specifying response, ability facet, main facets, and
#'   optional two-way interactions (e.g., `score ~ person + item + rater + item:rater`).
#' @param data data.frame containing rating observations.
#' @param K Number of rating categories (`1..K`): scalar integer `>= 2` for RSM
#'   or uniform PCM, a named integer vector matching `threshold_facet` levels for
#'   heterogeneous-category PCM, or `NULL` (inferred from `max(score)`).
#' @param priors Optional list of [prior()] specifications. Unspecified blocks
#'   receive [get_default_priors()].
#' @param stan_file Optional path to a custom `.stan` file (bypasses Stan code
#'   generation; `priors` and `hierarchical` cannot be supplied).
#' @param family Model family: `"rsm"` (shared step thresholds) or `"pcm"`
#'   (level-specific step thresholds on `threshold_facet`).
#' @param threshold_facet Required for `family = "pcm"`: main facet name whose
#'   levels have distinct step thresholds.
#' @param person_facet Name of the unconstrained examinee / ability main facet in
#'   `formula` (default `"person"`).
#' @param adjust_facets Optional character vector of nuisance main facets (e.g.,
#'   `"rater"`, `"judge"`) to adjust in fair-score calculations. When `NULL`,
#'   auto-detects `"rater"` and `"interlocutor"` if present.
#' @param model_name Optional base name for compiled Stan binary caching.
#' @param cache_dir Directory for compiled Stan binaries and fit caches. Defaults
#'   to `tempdir()` (CRAN-compliant, never pollutes the working directory); pass
#'   `tools::R_user_dir("BayesMFRM", "cache")` for persistent cross-session caching.
#' @param save_csvs Logical; must remain `TRUE` because `cmdstanr` draw accessors
#'   read output CSVs. Use `save_gq = "minimal"` (default) to keep CSVs compact.
#' @param save_gq Either `"minimal"` (default; omits observation-level Stan
#'   `generated quantities` and streams all post-processing in R) or `"full"`
#'   (writes `log_lik`, `y_rep`, `mu`, `sigma2`, `mu_fair` to Stan CSVs).
#' @param design_check Structural rank check mode before fitting: `"warn"`
#'   (default), `"error"`, or `"none"`.
#' @param cpp_options Named list of C++ options passed to `cmdstanr::cmdstan_model()`.
#' @param refit Cache policy: `"on_change"` (default; reuses cache when a fixed
#'   `seed` is supplied and code/data/settings match), `"never"`, or `"always"`.
#' @param iter Total iterations per chain (warmup + sampling; default `2000`).
#' @param warmup Warmup iterations per chain (default `floor(iter / 2)`).
#' @param chains Number of Markov chains (default `4`).
#' @param cores Number of parallel chains (default `getOption("mc.cores", 1L)`).
#' @param seed Optional integer RNG seed (`NULL` runs fresh random chains).
#' @param hierarchical Character vector of hierarchical facets/terms. Defaults to
#'   `person_facet` (empirical/Bayesian partial pooling of examinee abilities,
#'   matching MML-EAP); pass `character(0)` for unpooled fixed-scale priors.
#' @param obs_id Optional column name or length-`nrow(data)` vector of unique
#'   observation IDs for model-comparison alignment across missing-score filtering.
#' @param rhat_max Convergence diagnostic threshold for split R-hat (default `1.01`).
#' @param ess_min Convergence diagnostic minimum bulk/tail ESS (default `400`).
#' @param ... Additional named arguments passed to `CmdStanModel$sample()` (e.g.,
#'   `adapt_delta = 0.95`, `max_treedepth = 12`).
#'
#' @return A `bmfrm_fit` object.
#' @export
bmfrm <- function(
  formula,
  data,
  K = NULL,
  priors = NULL,
  stan_file = NULL,
  family = c("rsm", "pcm"),
  threshold_facet = NULL,
  person_facet = "person",
  adjust_facets = NULL,
  model_name = NULL,
  cache_dir = tempdir(),
  save_csvs = TRUE,
  save_gq = c("minimal", "full"),
  design_check = c("warn", "error", "none"),
  cpp_options = list(),
  refit = c("on_change", "never", "always"),
  iter = 2000,
  warmup = floor(iter / 2),
  chains = 4,
  cores = getOption("mc.cores", 1L),
  seed = NULL,
  hierarchical = person_facet,
  obs_id = NULL,
  rhat_max = 1.01,
  ess_min = 400,
  ...
) {
  family <- match.arg(family)
  refit <- match.arg(refit)
  save_gq_explicit <- !missing(save_gq)
  hierarchical_explicit <- !missing(hierarchical)
  save_gq <- match.arg(save_gq)
  design_check <- match.arg(design_check)

  if (!isTRUE(save_csvs)) {
    stop("save_csvs must be TRUE: cached CmdStan fits require CSV files; use save_gq = 'minimal' to reduce output")
  }
  if (!is.null(stan_file)) {
    if (!is.null(priors)) {
      stop("priors cannot be supplied with stan_file; specify priors inside the Stan file")
    }
    if (!hierarchical_explicit) {
      hierarchical <- character()
    }
    if (!save_gq_explicit) {
      save_gq <- "full"
    } else if (identical(save_gq, "minimal")) {
      stop("save_gq = 'minimal' requires generated Stan code; custom stan_file models define their own generated quantities")
    }
  }
  if ("model" %in% names(sys.call())) stop("The development model argument was removed; use hierarchical and priors")

  if (!is.numeric(iter) || length(iter) != 1L || !is.finite(iter) || iter != floor(iter) || iter < 2L) {
    stop("iter must be a finite integer >= 2")
  }
  if (!is.numeric(warmup) || length(warmup) != 1L || !is.finite(warmup) || warmup != floor(warmup) ||
      warmup < 0L || warmup >= iter) {
    stop("warmup must be a non-negative integer strictly less than iter")
  }
  if (!is.numeric(chains) || length(chains) != 1L || !is.finite(chains) || chains != floor(chains) || chains < 1L) {
    stop("chains must be a positive integer")
  }
  if (is.null(cores) || (length(cores) == 1L && is.na(cores) && !("cores" %in% names(sys.call())))) {
    cores <- 1L
  }
  if (!is.numeric(cores) || length(cores) != 1L || !is.finite(cores) || cores != floor(cores) || cores < 1L) {
    stop("cores must be a positive integer")
  }
  dots <- list(...)
  if (length(dots) > 0L) {
    if (is.null(names(dots)) || any(!nzchar(names(dots)))) {
      stop("Additional sampling arguments in '...' must be named")
    }
    reserved_dots <- intersect(
      names(dots),
      c("data", "seed", "chains", "parallel_chains", "output_dir",
        "iter_sampling", "iter_warmup", "save_warmup")
    )
    if (length(reserved_dots) > 0L) {
      stop("Do not pass reserved cmdstanr argument(s) (",
           paste(reserved_dots, collapse = ", "),
           ") in '...'; use bmfrm()'s iter, warmup, chains, cores, and seed arguments instead")
    }
  }

  if (!check_cmdstan()) {
    stop("CmdStan not found. Please install it using install_cmdstan_if_needed()")
  }

  message("Parsing model formula...")
  spec <- parse_bmfrm_formula(
    formula = formula,
    data = data,
    K = K,
    family = family,
    threshold_facet = threshold_facet,
    person_facet = person_facet,
    adjust_facets = adjust_facets,
    model_name = model_name,
    cache_dir = cache_dir,
    hierarchical = hierarchical,
    save_gq = save_gq
  )
  spec$hierarchical_explicit <- hierarchical_explicit
  spec$custom_stan <- !is.null(stan_file)
  spec <- resolve_bmfrm_spec(spec, priors)
  priors <- spec$priors
  if (!is.null(stan_file) && !is.null(spec$gaussian_model)) {
    stop("Hierarchical ability and marginal scaling require generated Stan code")
  }

  message("Preparing data...")
  data_info <- prepare_data_bmfrm(spec, data, obs_id = obs_id)
  if (!is.null(data_info$spec)) spec <- data_info$spec

  design_info <- NULL
  if (design_check != "none" && is.null(stan_file)) {
    design_K <- if (isTRUE(spec$K_heterogeneous)) data_info$K_by_level else spec$K
    design_info <- check_design_bmfrm(
      formula, data_info$data_clean, K = design_K, family = family,
      threshold_facet = spec$threshold_facet, person_facet = spec$person_facet
    )
    if (!design_info$identifiable) {
      msg <- paste0("Rank-deficient ", toupper(family), " design: ",
                    design_info$adjacent_logit_parameters - design_info$adjacent_logit_rank,
                    " unidentified direction(s). Proper priors do not establish likelihood identification.")
      if (design_check == "error") stop(msg) else warning(msg, call. = FALSE)
    }
  }

  if (is.null(stan_file)) {
    message("Generating Stan code...")
    stan_code <- build_stan_code(spec, priors)
  } else {
    message("Using hardcoded Stan file: ", basename(stan_file))
    stan_code <- paste(readLines(stan_file, warn = FALSE), collapse = "\n")
  }
  spec$eval_contract <- resolve_eval_contract(spec, stan_code)

  message("Compiling model...")
  comp <- compile_bmfrm(
    stan_code = stan_code,
    model_name = spec$model_name,
    cache_dir = spec$cache_dir,
    cpp_options = cpp_options
  )

  message("Fitting model...")
  fit_code <- if (!is.null(comp$model) && is.function(comp$model$code)) comp$model$code() else stan_code
  fit_cache_file <- resolve_fit_cache_file(
    base_cache_file = comp$cache_file,
    code = fit_code,
    data = data_info$stan_data,
    seed = seed,
    chains = chains,
    parallel_chains = min(chains, cores),
    extra_sampling = c(
      list(iter_sampling = iter - warmup, iter_warmup = warmup, save_warmup = FALSE),
      dots
    )
  )
  fit <- fit_cmdstan_cached(
    model = comp$model,
    data = data_info$stan_data,
    file = fit_cache_file,
    refit = refit,
    seed = seed,
    chains = chains,
    parallel_chains = min(chains, cores),
    iter_sampling = iter - warmup,
    iter_warmup = warmup,
    save_warmup = FALSE,
    ...
  )

  message("Post-processing results...")
  out <- postprocess_bmfrm(
    fit = fit,
    spec = spec,
    data_info = data_info,
    priors = priors,
    stan_code = stan_code,
    stan_code_hash = comp$hash,
    rhat_max = rhat_max,
    ess_min = ess_min
  )

  out$design <- design_info
  out$fit_args <- list(
    formula = formula,
    K = spec$K,
    priors = priors,
    stan_file = stan_file,
    family = family,
    threshold_facet = spec$threshold_facet,
    person_facet = spec$person_facet,
    adjust_facets = if (isTRUE(spec$adjust_facets_explicit)) spec$adjust_facets else NULL,
    hierarchical = spec$hierarchical,
    save_gq = spec$save_gq,
    cache_dir = cache_dir,
    cpp_options = cpp_options,
    iter = iter,
    warmup = warmup,
    chains = chains,
    cores = cores,
    seed = seed,
    obs_id = obs_id,
    rhat_max = rhat_max,
    ess_min = ess_min,
    sampler_args = list(...)
  )
  message("Model fitting complete!")
  out
}

#' Post-process bmfrm fit results
#' @keywords internal
postprocess_bmfrm <- function(fit, spec, data_info, priors, stan_code, stan_code_hash,
                              rhat_max = 1.01, ess_min = 400, ebfmi_min = 0.3) {
  convergence_info <- check_convergence(
    fit,
    variables = model_diagnostic_variables(spec, data_info),
    rhat_max = rhat_max,
    ess_min = ess_min,
    ebfmi_min = ebfmi_min
  )

  out <- structure(
    list(
      fit = fit,
      spec = spec,
      data_info = data_info,
      stan_code = stan_code,
      stan_code_hash = stan_code_hash,
      priors = priors,
      hierarchical = spec$hierarchical,
      convergence = convergence_info,
      timestamp = Sys.time(),
      bmfrm_version = utils::packageVersion("BayesMFRM"),
      cmdstan_version = cmdstanr::cmdstan_version()
    ),
    class = "bmfrm_fit"
  )

  if (!convergence_info$converged) {
    warning("Model may not have converged. Check diagnostics with summary(fit).")
  } else if (isTRUE(convergence_info$efficiency_warning)) {
    message("Note: ", convergence_info$n_max_treedepth,
            " transition(s) hit max_treedepth (sampler efficiency note; consider raising max_treedepth).")
  }

  out
}

#' Check MCMC convergence and sampler diagnostics
#'
#' Evaluates split R-hat, bulk/tail effective sample sizes, divergent
#' transitions, and energy Bayesian fraction of missing information (E-BFMI).
#' Maximum treedepth hits are recorded as an efficiency warning
#' (`efficiency_warning = TRUE`) rather than non-convergence (`converged = FALSE`),
#' following Stan diagnostic conventions.
#'
#' @param fit CmdStanMCMC object.
#' @param variables Optional parameter names.
#' @param rhat_max Upper threshold for split R-hat (default `1.01`).
#' @param ess_min Minimum bulk and tail ESS (default `400`).
#' @param ebfmi_min Minimum per-chain E-BFMI (default `0.3`).
#' @return Named list of convergence and efficiency diagnostics.
#' @export
check_convergence <- function(fit, variables = NULL, rhat_max = 1.01, ess_min = 400, ebfmi_min = 0.3) {
  if (!is.numeric(rhat_max) || length(rhat_max) != 1L || !is.finite(rhat_max) || rhat_max <= 1) {
    stop("rhat_max must be a finite number > 1")
  }
  if (!is.numeric(ess_min) || length(ess_min) != 1L || !is.finite(ess_min) || ess_min <= 0) {
    stop("ess_min must be a finite positive number")
  }
  if (!is.numeric(ebfmi_min) || length(ebfmi_min) != 1L || !is.finite(ebfmi_min) || ebfmi_min <= 0) {
    stop("ebfmi_min must be a finite positive number")
  }
  if (is.null(variables)) {
    metadata <- fit$metadata()
    parameter_names <- metadata$stan_variables %||% metadata$model_params
    parameter_names <- unique(sub("\\[.*$", "", parameter_names))
    raw <- parameter_names[grepl("_raw$", parameter_names)]
    if (length(raw)) variables <- raw
  }
  summ <- fit$summary(variables, "mean", "sd", "rhat", "ess_bulk", "ess_tail", "mcse_mean")
  sampler <- fit$diagnostic_summary(quiet = TRUE)
  chains <- fit$num_chains()
  required <- c("rhat", "ess_bulk", "ess_tail", "mcse_mean")
  complete <- all(required %in% names(summ)) && nrow(summ) > 0L &&
    all(vapply(summ[required], function(x) all(is.finite(x)), logical(1))) &&
    length(sampler$ebfmi) == chains && all(is.finite(sampler$ebfmi)) &&
    length(sampler$num_divergent) == chains && all(is.finite(sampler$num_divergent)) &&
    length(sampler$num_max_treedepth) == chains && all(is.finite(sampler$num_max_treedepth))
  finite_min <- function(x) if (length(x) && all(is.finite(x))) min(x) else NA_real_
  finite_max <- function(x) if (length(x) && all(is.finite(x))) max(x) else NA_real_

  ok <- complete && chains >= 2L &&
    all(summ$rhat < rhat_max) &&
    all(summ$ess_bulk >= ess_min) &&
    all(summ$ess_tail >= ess_min) &&
    all(sampler$num_divergent == 0) &&
    all(sampler$ebfmi >= ebfmi_min)

  n_td <- sum(sampler$num_max_treedepth)
  eff_warn <- complete && (n_td > 0L)

  list(
    converged = ok,
    status = if (!complete || chains < 2L) "incomplete" else if (!ok) "warning" else if (eff_warn) "ok_treedepth_hit" else "ok",
    efficiency_warning = eff_warn,
    thresholds = list(rhat_max = rhat_max, ess_min = ess_min, ebfmi_min = ebfmi_min),
    n_chains = chains,
    parameters = summ,
    per_chain = sampler,
    max_rhat = finite_max(summ$rhat),
    rhat_problems = sum(summ$rhat >= rhat_max, na.rm = TRUE),
    min_ess_bulk = finite_min(summ$ess_bulk),
    min_ess_tail = finite_min(summ$ess_tail),
    low_ess_bulk = sum(summ$ess_bulk < ess_min, na.rm = TRUE),
    low_ess_tail = sum(summ$ess_tail < ess_min, na.rm = TRUE),
    n_divergent = sum(sampler$num_divergent),
    n_max_treedepth = n_td,
    min_ebfmi = finite_min(sampler$ebfmi),
    max_mcse_mean = finite_max(summ$mcse_mean)
  )
}
