#' Cached cmdstanr fitting with hash-based validation
#'
#' This function provides a robust caching mechanism for cmdstanr fits, with
#' hash-based validation to ensure the cached fit matches the current model and data.
#'
#' @param model compiled cmdstan_model
#' @param data list passed to model$sample()
#' @param file path to .rds cache file
#' @param refit one of "never", "always", "on_change"
#' @param seed random seed
#' @param chains number of chains
#' @param parallel_chains number of parallel chains
#' @param ... additional cmdstanr::sample() args
#'
#' @return cmdstanr fit object
#' @export
fit_cmdstan_cached <- function(model, data, file, refit = c("never", "always", "on_change"),
                               seed = NULL, chains = 4, parallel_chains = chains, ...) {
  refit <- match.arg(refit)
  dots <- list(...)
  if (length(dots) && (is.null(names(dots)) || any(names(dots) == ""))) stop("Sampling arguments must be named")
  if (any(c("data", "seed", "chains", "parallel_chains", "output_dir") %in% names(dots))) stop("Reserved sampling argument supplied twice")
  ordered_dots <- if (length(dots)) dots[order(names(dots))] else list()
  sampling <- c(list(seed = seed, chains = chains, parallel_chains = parallel_chains), ordered_dots)
  signature <- list(schema = 2L, code = model$code(), data = data, sampling = sampling,
                    cmdstan = tryCatch(as.character(cmdstanr::cmdstan_version()), error = function(e) NA_character_))
  current_hash <- digest::digest(signature)
  if (refit != "always") {
    if (file.exists(file)) {
      obj <- tryCatch(readRDS(file), error = function(e) NULL)
      matches <- is.list(obj) && identical(obj$meta$schema, 2L) && identical(obj$meta$hash, current_hash)
      files <- if (matches) tryCatch(obj$fit$output_files(), error = function(e) character()) else character()
      available <- length(files) > 0L && all(file.exists(files))
      if (matches && available) {
        message("Loading cached fit (model, data, and sampling settings unchanged).")
        return(obj$fit)
      }
    }
    if (refit == "never") stop("Cached fit cannot be verified or its output files are missing; use refit = 'on_change' or 'always'")
  }
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  csv_dir <- file.path(dirname(file), "csv")
  dir.create(csv_dir, recursive = TRUE, showWarnings = FALSE)
  message("Fitting model and saving verified cache...")
  fit <- do.call(model$sample, c(list(data = data, output_dir = csv_dir), sampling))
  meta <- list(schema = 2L, hash = current_hash, signature = signature)
  saveRDS(list(fit = fit, meta = meta), file = file)
  fit
}

#' Resolve a signature-specific .rds fit cache path while honoring matching legacy caches
#' @keywords internal
resolve_fit_cache_file <- function(base_cache_file, code, data, seed, chains,
                                   parallel_chains, extra_sampling = list()) {
  ordered_extra <- if (length(extra_sampling)) extra_sampling[order(names(extra_sampling))] else list()
  sampling <- c(list(seed = seed, chains = chains, parallel_chains = parallel_chains), ordered_extra)
  signature <- list(
    schema = 2L,
    code = code,
    data = data,
    sampling = sampling,
    cmdstan = tryCatch(as.character(cmdstanr::cmdstan_version()), error = function(e) NA_character_)
  )
  sig_hash <- digest::digest(signature)
  stem <- tools::file_path_sans_ext(base_cache_file)
  specific_file <- paste0(stem, "_", substr(sig_hash, 1L, 12L), ".rds")
  if (!file.exists(specific_file) && file.exists(base_cache_file)) {
    legacy_obj <- tryCatch(readRDS(base_cache_file), error = function(e) NULL)
    if (is.list(legacy_obj) && identical(legacy_obj$meta$schema, 2L) &&
        identical(legacy_obj$meta$hash, sig_hash)) {
      return(base_cache_file)
    }
  }
  specific_file
}

`%||%` <- function(x, y) if (is.null(x)) y else x
