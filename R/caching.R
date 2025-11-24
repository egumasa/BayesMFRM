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
fit_cmdstan_cached <- function(
  model,           # compiled cmdstan_model
  data,            # list passed to model$sample()
  file,            # path to .rds cache file
  refit = c("never", "always", "on_change"),
  seed = NULL,
  chains = 4,
  parallel_chains = chains,
  ...              # additional cmdstanr::sample() args
) {
  refit <- match.arg(refit)

  # ---- helpers ----
  compute_hash <- function(model, data) {
    digest::digest(list(code = model$code(), data = data))
  }

  get_model_name_safe <- function(model) {
    md <- try(model$metadata(), silent = TRUE)
    if (inherits(md, "try-error") || is.null(md$model_name)) {
      return(NA_character_)
    } else {
      return(md$model_name)
    }
  }

  # ---- ensure parent folder exists ----
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)

  # ---- persistent CSV directory ----
  csv_dir <- file.path(dirname(file), "csv")
  dir.create(csv_dir, recursive = TRUE, showWarnings = FALSE)

  # ---- current model signature ----
  current_hash       <- compute_hash(model, data)
  current_model_name <- get_model_name_safe(model)

  # -------------------------------------------------------------------
  # Case 1: Explicit refit = "always"
  # -------------------------------------------------------------------
  if (refit == "always") {
    message("Refitting model (refit = 'always')...")
    fit <- model$sample(
      data = data,
      seed = seed,
      output_dir = csv_dir,
      chains = chains,
      parallel_chains = parallel_chains,
      ...
    )

    meta <- list(
      hash       = current_hash,
      model_name = current_model_name
    )

    saveRDS(list(fit = fit, meta = meta), file = file)
    return(fit)
  }

  # -------------------------------------------------------------------
  # Case 2: Cached file exists → validate before reuse
  # -------------------------------------------------------------------
  if (file.exists(file)) {
    obj <- readRDS(file)

    # New-style object with meta
    if (is.list(obj) && !is.null(obj$fit) && !is.null(obj$meta$hash)) {

      cached_hash       <- obj$meta$hash
      cached_model_name <- obj$meta$model_name %||% NA_character_

      # --- safety check: model name mismatch ---
      if (!is.na(current_model_name) &&
          !is.na(cached_model_name) &&
          !identical(current_model_name, cached_model_name)) {

        stop(
          "Cached file '", file, "' was created for model '", cached_model_name,
          "', but current model is '", current_model_name, "'.\n",
          "To avoid contamination, please use a different 'file' path or delete the old cache."
        )
      }

      # --- safety check: hash mismatch ---
      if (!identical(cached_hash, current_hash)) {
        if (refit == "never") {
          stop(
            "Cached file '", file, "' does not match current model+data hash, ",
            "but refit = 'never'.\n",
            "Refusing to reuse potentially contaminated cache. ",
            "Use refit = 'on_change' or 'always', or delete the file."
          )
        }

        # refit == "on_change"
        message("Model or data changed → refitting (refit = 'on_change').")
        fit <- model$sample(
          data = data,
          seed = seed,
          output_dir = csv_dir,
          chains = chains,
          parallel_chains = parallel_chains,
          ...
        )

        meta <- list(
          hash       = current_hash,
          model_name = current_model_name
        )

        saveRDS(list(fit = fit, meta = meta), file = file)
        return(fit)
      }

      # --- hashes match → safe to reuse ---
      message("Loading cached fit (model & data unchanged).")
      return(obj$fit)
    }

    # -----------------------------------------------------------------
    # Old-style RDS: just a CmdStanMCMC object, no meta
    # -----------------------------------------------------------------
    if (refit == "never") {
      message("Loading old-style cached fit (no meta; refit = 'never'). ",
              "Cannot fully verify contamination, use with care.")
      return(obj)
    }

    if (refit == "on_change") {
      message("Migrating old cached fit (no meta) → wrapping & adding hash (no refit).")

      # Ensure CSVs exist
      missing <- !file.exists(obj$output_files())
      if (any(missing)) {
        stop("Old cached fit refers to missing CSV files. You must refit the model.")
      }

      meta <- list(
        hash       = current_hash,
        model_name = current_model_name
      )

      obj_new <- list(fit = obj, meta = meta)
      saveRDS(obj_new, file = file)
      return(obj_new$fit)
    }
  } else {
    message("No cached file found → fitting model.")
  }

  # -------------------------------------------------------------------
  # Case 3: Must refit (no cache, or migrated/invalid)
  # -------------------------------------------------------------------
  message("Refitting model...")
  fit <- model$sample(
    data = data,
    seed = seed,
    output_dir = csv_dir,
    chains = chains,
    parallel_chains = parallel_chains,
    ...
  )

  meta <- list(
    hash       = current_hash,
    model_name = current_model_name
  )

  saveRDS(list(fit = fit, meta = meta), file = file)
  return(fit)
}

# Helper for null-coalescing operator
`%||%` <- function(x, y) if (is.null(x)) y else x
