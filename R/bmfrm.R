#' Fit Bayesian Many-Facet Rasch Models
#'
#' This is the main function for fitting Bayesian many-facet Rasch models using Stan.
#' It provides a user-friendly interface for specifying models with multiple facets
#' and optional bias terms.
#'
#' @param formula R formula specifying the model (e.g., score ~ person + item + rater + rater:item)
#' @param data data.frame containing all variables in the formula
#' @param K number of rating categories (if NULL, inferred from max(score))
#' @param priors list of prior() objects (if NULL, uses defaults)
#' @param stan_file optional path to a hardcoded Stan file (bypasses code generation)
#' @param family model family: "rsm" (Rating Scale Model with shared thresholds) or
#'   "pcm" (Partial Credit Model with facet-specific thresholds)
#' @param threshold_facet for PCM only: which facet gets item-specific thresholds (e.g., "item")
#' @param model_name optional name for the model (auto-generated if NULL)
#' @param cache_dir directory for caching compiled models and fits
#' @param design_check How to handle a rank-deficient generated additive RSM: warn, error, or none.
#' @param cpp_options Named compilation options passed to cmdstanr.
#' @param save_csvs whether to save MCMC output as CSV files
#' @param refit refit behavior: "never", "always", or "on_change"
#' @param iter number of MCMC iterations
#' @param warmup number of warmup iterations
#' @param chains number of MCMC chains
#' @param cores number of CPU cores to use
#' @param seed random seed
#' @param ... additional arguments passed to cmdstanr::sample()
#'
#' @return bmfrm_fit object containing the fitted model and metadata
#' @export
#'
#' @examples
#' \dontrun{
#' # Basic 3-facet RSM model (shared thresholds)
#' fit <- bmfrm(score ~ person + item + rater, data = ratings, K = 6, family = "rsm")
#'
#' # PCM model with item-specific thresholds
#' fit <- bmfrm(score ~ person + item + rater, data = ratings, K = 6,
#'              family = "pcm", threshold_facet = "item")
#'
#' # Model with bias term
#' fit <- bmfrm(score ~ person + item + rater + rater:item, data = ratings)
#'
#' # Custom priors
#' priors <- c(
#'   prior("normal(0, 2)", class = "theta"),
#'   prior("normal(0, 1)", class = "item"),
#'   prior("student_t(3, 0, 1)", class = "rater"),
#'   prior("normal(0, 0.5)", class = "bias", facet = "rater:item"),
#'   prior("normal(0, 3)", class = "tau")
#' )
#' fit <- bmfrm(score ~ person + item + rater + rater:item,
#'              data = ratings, priors = priors)
#' }
bmfrm <- function(
  formula,
  data,
  K = NULL,
  priors = NULL,
  stan_file = NULL,
  family = c("rsm", "pcm"),
  threshold_facet = NULL,
  model_name = NULL,
  cache_dir = "stan_cache",
  save_csvs = TRUE,
  design_check = c("warn", "error", "none"),
  cpp_options = list(),
  refit = c("on_change", "never", "always"),
  iter = 2000,
  warmup = floor(iter / 2),
  chains = 4,
  cores = getOption("mc.cores", 1L),
  seed = NULL,
  ...
) {

  # Argument validation
  family <- match.arg(family)
  refit <- match.arg(refit)
  design_check <- match.arg(design_check)
  
  # Check cmdstan availability
  if (!check_cmdstan()) {
    stop("CmdStan not found. Please install it using install_cmdstan_if_needed()")
  }
  
  # 1. Parse formula → spec
  message("Parsing model formula...")
  spec <- parse_bmfrm_formula(
    formula = formula,
    data = data,
    K = K,
    family = family,
    threshold_facet = threshold_facet,
    model_name = model_name,
    cache_dir = cache_dir
  )
  
  # 2. Prepare data → stan_data + facet labels  
  message("Preparing data...")
  data_info <- prepare_data_bmfrm(spec, data)
  
  design_info <- NULL
  if (design_check != "none" && is.null(stan_file) && family == "rsm" && !length(spec$facets_bias)) {
    design_info <- check_design_bmfrm(formula, data_info$data_clean, K = spec$K)
    if (!design_info$identifiable) {
      msg <- paste0("Rank-deficient RSM design: ", design_info$nullity,
                    " unidentified direction(s). Proper priors do not establish likelihood identification.")
      if (design_check == "error") stop(msg) else warning(msg, call. = FALSE)
    }
  }

  # 3. Validate/set priors
  message("Setting up priors...")
  priors <- validate_priors(priors, spec$facets_main, spec$facets_bias, spec$family)
  
  # 4. Generate or use provided Stan code
  if (is.null(stan_file)) {
    message("Generating Stan code...")
    stan_code <- build_stan_code(spec, priors)
  } else {
    message("Using hardcoded Stan file: ", basename(stan_file))
    stan_code <- paste(readLines(stan_file), collapse = "\n")
  }
  
  # 5. Compile Stan model (with caching based on code hash)
  message("Compiling model...")
  comp <- compile_bmfrm(
    stan_code = stan_code,
    model_name = spec$model_name,
    cache_dir = spec$cache_dir,
    cpp_options = cpp_options
  )
  
  # 6. Fit using existing fit_cmdstan_cached()
  message("Fitting model...")
  fit <- fit_cmdstan_cached(
    model = comp$model,
    data = data_info$stan_data,
    file = comp$cache_file,
    refit = refit,
    seed = seed,
    chains = chains,
    parallel_chains = min(chains, cores),
    iter_sampling = iter - warmup,
    iter_warmup = warmup,
    save_warmup = FALSE,
    ...
  )
  
  # 7. Wrap as bmfrm_fit object
  message("Post-processing results...")
  out <- postprocess_bmfrm(
    fit = fit,
    spec = spec,
    data_info = data_info,
    priors = priors,
    stan_code = stan_code,
    stan_code_hash = comp$hash
  )
  
  out$design <- design_info
  message("Model fitting complete!")
  out
}

#' Post-process bmfrm fit results
#'
#' @param fit cmdstanr fit object
#' @param spec bmfrm_spec object
#' @param data_info data preparation results
#' @param priors list of priors used
#' @param stan_code Stan code string
#' @param stan_code_hash hash of Stan code
#'
#' @return bmfrm_fit object
#' @keywords internal
postprocess_bmfrm <- function(fit, spec, data_info, priors, stan_code, stan_code_hash) {
  
  # Check convergence
  convergence_info <- check_convergence(fit)
  
  # Create bmfrm_fit object
  out <- structure(
    list(
      # Core fitting results
      fit = fit,
      
      # Model specification
      spec = spec,
      
      # Data and labels
      data_info = data_info,
      
      # Stan code and compilation
      stan_code = stan_code,
      stan_code_hash = stan_code_hash,
      
      # Priors used
      priors = priors,
      
      # Convergence diagnostics
      convergence = convergence_info,
      
      # Metadata
      timestamp = Sys.time(),
      bmfrm_version = utils::packageVersion("BayesMFRM"),
      cmdstan_version = cmdstanr::cmdstan_version()
    ),
    class = "bmfrm_fit"
  )
  
  # Warn about convergence issues
  if (!convergence_info$converged) {
    warning("Model may not have converged. Check diagnostics with summary(fit).")
  }
  
  out
}

#' Check MCMC convergence
#' @param fit CmdStanMCMC object.
#' @param variables Optional parameter names. By default monitor raw model parameters
#'   where available, avoiding deterministic generated quantities.
#' @return Individual diagnostics and a conservative status; this is not a proof of convergence.
#' @keywords internal
check_convergence <- function(fit, variables = NULL) {
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
  ok <- complete && chains >= 2L && all(summ$rhat < 1.01) &&
    all(summ$ess_bulk >= 400) && all(summ$ess_tail >= 400) &&
    all(sampler$num_divergent == 0) && all(sampler$num_max_treedepth == 0) &&
    all(sampler$ebfmi >= 0.3)
  list(converged = ok, status = if (!complete || chains < 2L) "incomplete" else if (ok) "ok" else "warning",
       n_chains = chains, parameters = summ, per_chain = sampler,
       max_rhat = finite_max(summ$rhat), rhat_problems = sum(summ$rhat >= 1.01, na.rm = TRUE),
       min_ess_bulk = finite_min(summ$ess_bulk), min_ess_tail = finite_min(summ$ess_tail),
       low_ess_bulk = sum(summ$ess_bulk < 400, na.rm = TRUE),
       low_ess_tail = sum(summ$ess_tail < 400, na.rm = TRUE),
       n_divergent = sum(sampler$num_divergent), n_max_treedepth = sum(sampler$num_max_treedepth),
       min_ebfmi = finite_min(sampler$ebfmi), max_mcse_mean = finite_max(summ$mcse_mean))
}
