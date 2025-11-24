#' Compile bmfrm Stan model with caching
#'
#' @param stan_code character string with Stan code
#' @param model_name name for the compiled model
#' @param cache_dir directory for caching compiled models
#'
#' @return list with model, hash, stan_file, and cache_file
#' @export
compile_bmfrm <- function(stan_code, model_name, cache_dir) {
  
  # Ensure cache directory exists
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Compute hash of Stan code
  code_hash <- digest::digest(stan_code)
  
  # Build file paths
  stan_file <- file.path(cache_dir, paste0(model_name, "_", code_hash, ".stan"))
  cache_file <- file.path(cache_dir, paste0(model_name, "_", code_hash, ".rds"))
  
  # Write Stan code to file if it doesn't exist
  if (!file.exists(stan_file)) {
    writeLines(stan_code, stan_file)
  }
  
  # Check if model is already compiled
  if (file.exists(paste0(stan_file, ".exe")) || file.exists(paste0(tools::file_path_sans_ext(stan_file), ".exe"))) {
    message("Loading cached compiled model...")
    model <- cmdstanr::cmdstan_model(stan_file)
  } else {
    # Compile the model
    message("Compiling Stan model...")
    model <- cmdstanr::cmdstan_model(stan_file)
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
#'
#' @return character vector of files that were (or would be) deleted
#' @export
clean_model_cache <- function(cache_dir, max_age_days = 30, dry_run = FALSE) {
  if (!dir.exists(cache_dir)) {
    message("Cache directory does not exist: ", cache_dir)
    return(character(0))
  }
  
  # Find all Stan files and executables
  stan_files <- list.files(cache_dir, pattern = "\\.stan$", full.names = TRUE)
  exe_files <- list.files(cache_dir, pattern = "\\.exe$", full.names = TRUE)
  rds_files <- list.files(cache_dir, pattern = "\\.rds$", full.names = TRUE)
  
  all_files <- c(stan_files, exe_files, rds_files)
  
  if (length(all_files) == 0) {
    message("No files found in cache directory.")
    return(character(0))
  }
  
  # Check file ages
  cutoff_time <- Sys.time() - (max_age_days * 24 * 60 * 60)
  file_times <- file.mtime(all_files)
  old_files <- all_files[file_times < cutoff_time]
  
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
      file.remove(f)
      deleted <- c(deleted, f)
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
  
  # Try to parse with cmdstan
  tryCatch({
    model <- cmdstanr::cmdstan_model(temp_file, compile = FALSE)
    # If we get here, syntax is valid
    list(is_valid = TRUE, message = "Stan code syntax is valid")
  }, error = function(e) {
    list(is_valid = FALSE, error = e$message)
  })
}
