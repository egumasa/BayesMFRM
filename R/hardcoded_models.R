#' Get path to a hardcoded Stan model
#'
#' @param name Model name (e.g., "three-facet/person-rater-scenario")
#' @return Absolute path to the .stan file
#' @keywords internal
get_hardcoded_model_path <- function(name) {
  # Normalize: add .stan if missing
  if (!grepl("\\.stan$", name)) {
    name <- paste0(name, ".stan")
  }

  # Try installed package location
  pkg_path <- system.file("stan", "models", name, package = "BayesMFRM")
  if (pkg_path != "" && file.exists(pkg_path)) {
    return(pkg_path)
  }

  # Try development location
  dev_path <- file.path("inst", "stan", "models", name)
  if (file.exists(dev_path)) {
    return(normalizePath(dev_path))
  }

  stop("Model '", name, "' not found.\n",
       "Use bmfrm_models() to see available models.", call. = FALSE)
}

#' List available hardcoded Stan models
#'
#' Lists all pre-defined Stan models available in the BayesMFRM package.
#' These models can be used with bmfrm() via the stan_file argument.
#'
#' @return Character vector of model names (invisibly)
#' @export
#'
#' @examples
#' \dontrun{
#' bmfrm_models()
#' }
bmfrm_models <- function() {
  # Try installed package location
  pkg_dir <- system.file("stan", "models", package = "BayesMFRM")

  if (pkg_dir == "" || !dir.exists(pkg_dir)) {
    # Try development location
    pkg_dir <- file.path("inst", "stan", "models")
  }

  if (!dir.exists(pkg_dir)) {
    message("No hardcoded models directory found.")
    return(character(0))
  }

  # Find all .stan files recursively
  files <- list.files(pkg_dir, pattern = "\\.stan$", recursive = TRUE)

  # Remove .stan extension for display
  models <- sub("\\.stan$", "", files)

  if (length(models) == 0) {
    message("No hardcoded models found.")
    return(character(0))
  }

  cat("Available hardcoded models:\n")
  for (m in models) {
    cat("  ", m, "\n")
  }
  invisible(models)
}

#' Load a hardcoded Stan model
#'
#' Loads a pre-defined Stan model from the BayesMFRM package.
#' The returned object contains the path to the Stan file, which can be
#' passed to bmfrm() via the stan_file argument.
#'
#' @param name Model name (e.g., "three-facet/person-rater-scenario-main")
#' @param compile Logical; compile the model immediately? Default FALSE.
#' @param cache_dir Directory for compiled model cache.
#'
#' @return A list with:
#'   \itemize{
#'     \item \code{name}: Model name
#'     \item \code{path}: Absolute path to the .stan file
#'     \item \code{code}: Stan code as character string
#'     \item \code{stanmodel}: Compiled CmdStanModel (if compile = TRUE)
#'   }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Get model path for use with bmfrm()
#' model <- bmfrm_model("four-facet/person-item-rater-scenario-main")
#'
#' fit <- bmfrm(
#'   score ~ person + item + rater + scenario,
#'   data = ratings,
#'   stan_file = model$path
#' )
#' }
bmfrm_model <- function(name, compile = FALSE, cache_dir = "stan_cache") {
  path <- get_hardcoded_model_path(name)
  code <- paste(readLines(path), collapse = "\n")

  result <- list(
    name = name,
    path = path,
    code = code,
    stanmodel = NULL
  )

  if (compile) {
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    result$stanmodel <- cmdstanr::cmdstan_model(path, dir = cache_dir)
  }

  class(result) <- "bmfrm_hardcoded_model"
  result
}

#' Internal registry of verified bundled Stan models and their R-engine compatibility
#'
#' Maps each bundled model under `inst/stan/models/` to its verified likelihood
#' structure, whether the R posterior-quantity engine matches its parameter names
#' and linear predictor, and which facets its built-in `mu_fair` adjusts.
#'
#' @keywords internal
bundled_stan_model_registry <- function() {
  list(
    "two-facet/person-rater-main" = list(
      family = "rsm",
      facets_main = c("person", "rater"),
      facets_bias = character(0),
      param_map = c(person = "theta", rater = "rater"),
      r_recompute = TRUE,
      stan_mu_fair_adjust = "rater",
      supported_ops = c("residuals", "facet_fit", "fair_scores", "standardized_fair_scores", "pp_check", "predict", "loo")
    ),
    "three-facet/person-rater-scenario-main" = list(
      family = "rsm",
      facets_main = c("person", "rater", "scenario"),
      facets_bias = character(0),
      param_map = c(person = "theta", rater = "rater", scenario = "scenario"),
      r_recompute = TRUE,
      stan_mu_fair_adjust = c("rater", "scenario"),
      supported_ops = c("residuals", "facet_fit", "fair_scores", "standardized_fair_scores", "pp_check", "predict", "loo")
    ),
    "four-facet/person-item-rater-scenario-main" = list(
      family = "rsm",
      facets_main = c("person", "item", "rater", "scenario"),
      facets_bias = character(0),
      param_map = c(person = "theta", item = "beta", rater = "rater", scenario = "scenario"),
      r_recompute = FALSE,
      stan_mu_fair_adjust = c("rater", "scenario"),
      supported_ops = c("stan_gq_only")
    )
  )
}

#' Match Stan source code exactly against verified bundled Stan models
#' @keywords internal
match_bundled_stan_model <- function(stan_code) {
  if (is.null(stan_code) || !is.character(stan_code) || length(stan_code) != 1L || !nzchar(stan_code)) {
    return(NULL)
  }
  norm_code <- sub("\n+$", "", gsub("\r\n", "\n", stan_code, fixed = TRUE))
  reg <- bundled_stan_model_registry()
  for (name in names(reg)) {
    path <- tryCatch(get_hardcoded_model_path(name), error = function(e) NULL)
    if (is.null(path) || !file.exists(path)) next
    ref_code <- paste(readLines(path, warn = FALSE), collapse = "\n")
    norm_ref <- sub("\n+$", "", gsub("\r\n", "\n", ref_code, fixed = TRUE))
    if (identical(norm_code, norm_ref)) {
      return(c(list(model_id = name), reg[[name]]))
    }
  }
  NULL
}

#' Print method for bmfrm_hardcoded_model objects
#' @param x bmfrm_hardcoded_model object
#' @param ... unused
#' @export
print.bmfrm_hardcoded_model <- function(x, ...) {
  cat("BayesMFRM Hardcoded Model\n")
  cat("  Name:", x$name, "\n")
  cat("  Path:", x$path, "\n")
  cat("  Compiled:", !is.null(x$stanmodel), "\n")
  invisible(x)
}
