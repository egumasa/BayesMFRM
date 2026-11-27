#' Compile bmfrm Stan model with caching
#'
#' @param stan_code character string with Stan code
#' @param model_name name for the compiled model
#' @param cache_dir directory for caching compiled models
#' @param cpp_options Named C++ compilation options passed to cmdstanr.
#'
#' @return list with model, hash, stan_file, and cache_file
#' @export
compile_bmfrm <- function(stan_code, model_name, cache_dir, cpp_options = list()) {
  
  # Ensure cache directory exists
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Compute hash of Stan code
  code_hash <- digest::digest(list(stan_code, cpp_options, cmdstanr::cmdstan_version()))
  
  # Build file paths
  stan_file <- file.path(cache_dir, paste0(model_name, "_", code_hash, ".stan"))
  cache_file <- file.path(cache_dir, paste0(model_name, "_", code_hash, ".rds"))
  
  # Write Stan code to file if it doesn't exist
  if (!file.exists(stan_file)) {
    writeLines(stan_code, stan_file)
  }
  
  # Check if model is already compiled (Windows .exe or extensionless Unix binary)
  exe_candidates <- c(
    paste0(stan_file, ".exe"),
    paste0(tools::file_path_sans_ext(stan_file), ".exe"),
    tools::file_path_sans_ext(stan_file)
  )
  if (any(file.exists(exe_candidates) & !dir.exists(exe_candidates))) {
    message("Loading cached compiled model...")
    model <- cmdstanr::cmdstan_model(stan_file, cpp_options = cpp_options)
  } else {
    # Compile the model
    message("Compiling Stan model...")
    model <- cmdstanr::cmdstan_model(stan_file, cpp_options = cpp_options)
  }
  
  list(
    model      = model,
    hash       = code_hash,
    stan_file  = stan_file,
    cache_file = cache_file
  )
}

#' Check if CmdStan is available
#'
#' @return TRUE if CmdStan is installed, FALSE otherwise
#' @export
check_cmdstan <- function() {
  tryCatch({
    cmdstanr::cmdstan_version()
    TRUE
  }, error = function(e) {
    FALSE
  })
}

#' Install CmdStan if needed (wrapper around cmdstanr)
#'
#' @param cores number of cores to use for compilation
#' @param quiet suppress installation messages
#'
#' @return invisible TRUE if successful
#' @export
install_cmdstan_if_needed <- function(cores = 2, quiet = FALSE) {
  if (!check_cmdstan()) {
    if (!quiet) {
      message("CmdStan not found. Installing...")
    }
    cmdstanr::install_cmdstan(cores = cores, quiet = quiet)
  } else {
    if (!quiet) {
      message("CmdStan is already installed.")
    }
  }
  invisible(TRUE)
}

#' Get model compilation info
#'
#' @param model_path path to compiled Stan model
#'
#' @return list with compilation details
#' @export
get_model_info <- function(model_path) {
  if (!file.exists(model_path)) {
    stop("Model file not found: ", model_path)
  }
  
  # Try to create model object to get info
  tryCatch({
    model <- cmdstanr::cmdstan_model(model_path)
    list(
      path = model_path,
      exists = TRUE,
      exe_time = file.mtime(model$exe_file()),
      code_hash = digest::digest(readLines(model_path)),
      cmdstan_version = cmdstanr::cmdstan_version()
    )
  }, error = function(e) {
    list(
      path = model_path,
      exists = FALSE,
      error = e$message
    )
  })
}

#' Clean old compiled models from cache
#'
#' @param cache_dir cache directory
#' @param max_age_days maximum age in days before cleaning
#' @param dry_run if TRUE, just show what would be deleted
#' @param include_csv if TRUE, also delete expired chain CSV files under
#'   `cache_dir/csv/` and `cache_dir/cv_folds/csv/` that are not referenced by a
#'   non-expired cache `.rds` in `cache_dir`. Defaults to `FALSE` so fitted
#'   `bmfrm_fit` objects saved outside `cache_dir` or held in memory retain
#'   their underlying CmdStan CSV outputs.
#'
#' @return character vector of files that were (or would be) deleted
#' @export
clean_model_cache <- function(cache_dir, max_age_days = 30, dry_run = FALSE, include_csv = FALSE) {
  if (!dir.exists(cache_dir)) {
    message("Cache directory does not exist: ", cache_dir)
    return(character(0))
  }
  
  # Find all Stan files, executables, and fit caches (including CV fold caches)
  stan_files <- list.files(cache_dir, pattern = "\\.stan$", full.names = TRUE)
  exe_files <- list.files(cache_dir, pattern = "\\.exe$", full.names = TRUE)
  rds_files <- list.files(cache_dir, pattern = "\\.rds$", full.names = TRUE)
  cv_dir <- file.path(cache_dir, "cv_folds")
  cv_rds_files <- if (dir.exists(cv_dir)) {
    list.files(cv_dir, pattern = "\\.rds$", full.names = TRUE)
  } else {
    character(0)
  }
  all_rds_files <- c(rds_files, cv_rds_files)
  
  # CmdStan executables on Unix have no extension. Include binaries paired with
  # a Stan source plus any orphaned <name>_<32-hex-hash> ELF/Mach-O executable,
  # so unrelated extensionless user files are never removed.
  is_unix_binary <- function(f) {
    if (!file.exists(f) || dir.exists(f)) return(FALSE)
    hdr <- tryCatch(readBin(f, what = "raw", n = 4L), error = function(e) raw(0))
    if (length(hdr) < 4L) return(FALSE)
    # ELF (\x7fELF) or Mach-O 32/64/fat magic bytes
    identical(hdr, as.raw(c(0x7f, 0x45, 0x4c, 0x46))) ||
      identical(hdr, as.raw(c(0xcf, 0xfa, 0xed, 0xfe))) ||
      identical(hdr, as.raw(c(0xfe, 0xed, 0xfa, 0xcf))) ||
      identical(hdr, as.raw(c(0xce, 0xfa, 0xed, 0xfe))) ||
      identical(hdr, as.raw(c(0xfe, 0xed, 0xfa, 0xce))) ||
      identical(hdr, as.raw(c(0xca, 0xfe, 0xba, 0xbe)))
  }
  paired_unix_exe <- tools::file_path_sans_ext(stan_files)
  cand_orphan <- list.files(cache_dir, pattern = "_[0-9a-f]{32}$", full.names = TRUE)
  orphan_unix_exe <- cand_orphan[vapply(cand_orphan, is_unix_binary, logical(1))]
  unix_exe_files <- unique(c(paired_unix_exe, orphan_unix_exe))
  cutoff_time <- Sys.time() - (max_age_days * 24 * 60 * 60)

  managed_csv_files <- character(0)
  if (isTRUE(include_csv)) {
    managed_csv_dirs <- c(file.path(cache_dir, "csv"), file.path(cv_dir, "csv"))
    managed_csv_dirs <- managed_csv_dirs[dir.exists(managed_csv_dirs)]
    managed_csv_files <- unlist(
      lapply(managed_csv_dirs, list.files, pattern = "\\.csv$", full.names = TRUE),
      use.names = FALSE
    )
    # Preserve any managed CSV still referenced by a non-expired .rds fit cache
    kept_rds <- if (length(all_rds_files)) {
      rt <- file.mtime(all_rds_files)
      all_rds_files[is.na(rt) | rt >= cutoff_time]
    } else {
      character(0)
    }
    active_csvs <- unique(unlist(lapply(kept_rds, function(rf) {
      obj <- tryCatch(readRDS(rf), error = function(e) NULL)
      if (is.list(obj) && !is.null(obj$fit)) {
        tryCatch(normalizePath(obj$fit$output_files(), mustWork = FALSE), error = function(e) character(0))
      } else {
        character(0)
      }
    }), use.names = FALSE))
    if (length(managed_csv_files) && length(active_csvs)) {
      norm_csv <- normalizePath(managed_csv_files, mustWork = FALSE)
      managed_csv_files <- managed_csv_files[!(norm_csv %in% active_csvs)]
    }
  }

  all_files <- unique(c(stan_files, exe_files, unix_exe_files, all_rds_files, managed_csv_files))
  all_files <- all_files[file.exists(all_files) & !dir.exists(all_files)]
  
  if (length(all_files) == 0) {
    message("No files found in cache directory.")
    return(character(0))
  }
  
  # Check file ages (keep .stan paired with its compiled binary so .stan is not
  # deleted while its newer binary remains in cache_dir)
  file_times <- file.mtime(all_files)
  names(file_times) <- all_files
  for (sf in intersect(stan_files, all_files)) {
    stem <- tools::file_path_sans_ext(sf)
    bins <- intersect(c(stem, paste0(stem, ".exe"), paste0(sf, ".exe")), all_files)
    if (length(bins)) {
      bt <- file_times[bins]
      if (any(!is.na(bt))) file_times[sf] <- max(c(file_times[sf], bt), na.rm = TRUE)
    }
  }
  old_files <- all_files[!is.na(file_times) & file_times < cutoff_time]
  
  if (length(old_files) == 0) {
    message("No old files found to clean.")
    return(character(0))
  }
  
  if (dry_run) {
    message("Would delete ", length(old_files), " old files:")
    for (f in old_files) {
      message("  ", f)
    }
    return(old_files)
  }
  
  # Actually delete files
  deleted <- character(0)
  for (f in old_files) {
    tryCatch({
      if (file.remove(f)) deleted <- c(deleted, f)
    }, error = function(e) {
      warning("Could not delete ", f, ": ", e$message)
    })
  }
  
  message("Deleted ", length(deleted), " old cache files.")
  deleted
}

#' Validate Stan code syntax
#'
#' @param stan_code character string with Stan code
#'
#' @return list with is_valid (logical) and any error messages
#' @export
validate_stan_code <- function(stan_code) {
  # Write to temporary file
  temp_file <- tempfile(fileext = ".stan")
  on.exit(unlink(temp_file))
  
  writeLines(stan_code, temp_file)
  diagnostics <- character()
  
  # Try to parse with cmdstan
  tryCatch({
    model <- cmdstanr::cmdstan_model(temp_file, compile = FALSE)
    withCallingHandlers(model$check_syntax(quiet = TRUE), message = function(m) {
      diagnostics <<- c(diagnostics, conditionMessage(m))
      invokeRestart("muffleMessage")
    })
    list(is_valid = TRUE, message = "Stan code syntax is valid")
  }, error = function(e) {
    list(is_valid = FALSE, error = paste(c(diagnostics, conditionMessage(e)), collapse = "\n"))
  })
}
