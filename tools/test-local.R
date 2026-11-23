# Run from any directory with Rscript --vanilla path/to/tools/test-local.R.
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
root <- dirname(dirname(normalizePath(sub("^--file=", "", script_arg[1]))))
setwd(root)
suppressPackageStartupMessages(pkgload::load_all(root, quiet = TRUE))
# The existing integration tests are opt-in; the replication project provides
# the validated small RSM workflow using the real CmdStan sampler.
Sys.setenv(NOT_CRAN = if (identical(Sys.getenv("BAYESMFRM_RUN_STAN"), "true")) "true" else "false")
x <- testthat::test_dir(file.path(root, "tests", "testthat"), reporter = "summary", stop_on_failure = TRUE)
df <- as.data.frame(x)
print(colSums(df[, c("failed", "error", "warning", "skipped", "passed")]))
