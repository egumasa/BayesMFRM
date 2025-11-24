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
#' @param family model family ("rating_scale" for v0.1)
#' @param model_name optional name for the model (auto-generated if NULL)
#' @param cache_dir directory for caching compiled models and fits
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
#' # Basic 3-facet model
#' fit <- bmfrm(score ~ person + item + rater, data = ratings, K = 6)
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
  family = c("rating_scale", "partial_credit"),
  model_name = NULL,
  cache_dir = "stan_cache", 
  save_csvs = TRUE,
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
  
  if (family != "rating_scale") {
    stop("Only 'rating_scale' family supported in v0.1")
  }
  
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
    model_name = model_name,
    cache_dir = cache_dir
  )
  
  # 2. Prepare data → stan_data + facet labels  
  message("Preparing data...")
  data_info <- prepare_data_bmfrm(spec, data)
  
  # 3. Validate/set priors
  message("Setting up priors...")
  priors <- validate_priors(priors, spec$facets_main, spec$facets_bias)
  
  # 4. Generate Stan code string
  message("Generating Stan code...")
  stan_code <- build_stan_code(spec, priors)
  
  # 5. Compile Stan model (with caching based on code hash)
  message("Compiling model...")
  comp <- compile_bmfrm(
    stan_code = stan_code,
    model_name = spec$model_name,
    cache_dir = spec$cache_dir
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
#'
#' @param fit cmdstanr fit object
#'
#' @return list with convergence information
#' @keywords internal
check_convergence <- function(fit) {
  
  # Get summary for key parameters
  summ <- fit$summary()
  
  # Check R-hat
  max_rhat <- max(summ$rhat, na.rm = TRUE)
  rhat_problems <- sum(summ$rhat > 1.01, na.rm = TRUE)
  
  # Check effective sample size
  min_ess_bulk <- min(summ$ess_bulk, na.rm = TRUE)
  min_ess_tail <- min(summ$ess_tail, na.rm = TRUE)
  low_ess_bulk <- sum(summ$ess_bulk < 400, na.rm = TRUE)
  low_ess_tail <- sum(summ$ess_tail < 400, na.rm = TRUE)
  
  # Check divergences and other diagnostics
  diagnostics <- fit$sampler_diagnostics()
  n_divergent <- sum(diagnostics[,,"divergent__"], na.rm = TRUE)
  n_max_treedepth <- sum(diagnostics[,,"treedepth__"] >= 10, na.rm = TRUE)
  
  # Overall convergence assessment
  converged <- (
    max_rhat < 1.01 &&
    rhat_problems == 0 &&
    min_ess_bulk >= 400 &&
    min_ess_tail >= 400 &&
    n_divergent == 0 &&
    n_max_treedepth < 0.01 * nrow(summ)
  )
  
  list(
    converged = converged,
    max_rhat = max_rhat,
    rhat_problems = rhat_problems,
    min_ess_bulk = min_ess_bulk,
    min_ess_tail = min_ess_tail,
    low_ess_bulk = low_ess_bulk,
    low_ess_tail = low_ess_tail,
    n_divergent = n_divergent,
    n_max_treedepth = n_max_treedepth
  )
}
