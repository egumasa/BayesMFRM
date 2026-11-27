#' Cached cmdstanr fitting with portable relative CSV resolution
#'
#' Caches compiled CmdStan fits with a full model/data/sampler signature hash.
#' Stores chain CSV paths both directly and relative to `dirname(file)` so
#' moving or renaming a project/cache directory does not invalidate saved caches.
#' When `seed = NULL` and `refit = "on_change"`, fresh chains are sampled rather
#' than silently returning a previous unseeded run; pass an explicit integer
#' `seed` (or `refit = "never"`) to reuse cached draws deterministically.
#'
#' @param model compiled cmdstan_model
#' @param data list passed to model$sample()
#' @param file path to .rds cache file
#' @param refit one of "never", "always", "on_change"
#' @param seed random seed (`NULL` draws fresh chains under `refit = "on_change"`)
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
  signature <- list(schema = 3L, code = model$code(), data = data, sampling = sampling,
                    cmdstan = tryCatch(as.character(cmdstanr::cmdstan_version()), error = function(e) NA_character_))
  current_hash <- digest::digest(signature)

  can_reuse_cache <- (refit == "never") || (refit == "on_change" && !is.null(seed))
  if (can_reuse_cache) {
    if (file.exists(file)) {
      obj <- tryCatch(readRDS(file), error = function(e) NULL)
      schema_ok <- is.list(obj) && isTRUE(obj$meta$schema %in% c(2L, 3L))
      # Schema 2/3 share identical payload fields; verify signature match
      sig_ok <- schema_ok && (identical(obj$meta$hash, current_hash) ||
        (identical(obj$meta$schema, 2L) &&
         identical(obj$meta$signature[c("code", "data", "sampling", "cmdstan")],
                   signature[c("code", "data", "sampling", "cmdstan")])))
      if (sig_ok) {
        files <- tryCatch(obj$fit$output_files(), error = function(e) character())
        if (length(files) > 0L && all(file.exists(files))) {
          message("Loading cached fit (model, data, and sampling settings unchanged).")
          return(obj$fit)
        }
        # Relocate portable relative CSV paths when the cache directory moved
        if (!is.null(obj$meta$csv_rel) && length(obj$meta$csv_rel) > 0L) {
          relocated <- file.path(dirname(file), obj$meta$csv_rel)
          if (all(file.exists(relocated))) {
            message("Loading relocated cached fit from relative CSV paths.")
            return(cmdstanr::as_cmdstan_fit(relocated))
          }
        }
      }
    }
    if (refit == "never") {
      stop("Cached fit cannot be verified or its output files are missing; use refit = 'on_change' or 'always'")
    }
  }

  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  csv_dir <- file.path(dirname(file), "csv")
  dir.create(csv_dir, recursive = TRUE, showWarnings = FALSE)
  message("Fitting model and saving verified cache...")
  fit <- do.call(model$sample, c(list(data = data, output_dir = csv_dir), sampling))
  out_files <- tryCatch(fit$output_files(), error = function(e) character())
  csv_rel <- if (length(out_files)) file.path("csv", basename(out_files)) else NULL
  meta <- list(schema = 3L, hash = current_hash, signature = signature, csv_rel = csv_rel)
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
    schema = 3L,
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
    if (is.list(legacy_obj) && isTRUE(legacy_obj$meta$schema %in% c(2L, 3L)) &&
        identical(legacy_obj$meta$signature[c("code", "data", "sampling", "cmdstan")],
                  signature[c("code", "data", "sampling", "cmdstan")])) {
      return(base_cache_file)
    }
  }
  specific_file
}

`%||%` <- function(x, y) if (is.null(x)) y else x
