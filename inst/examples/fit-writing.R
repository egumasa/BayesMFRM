# Shared fitted model for the help examples. Requires CmdStan.
data(bmfrm_example, package = "BayesMFRM")
ratings <- bmfrm_example$ratings
model <- score ~ person + item + rater + item:rater
fit <- bmfrm(model, ratings, K = 4, design_check = "error",
             seed = 2409, chains = 4, cores = 2, iter = 2000, adapt_delta = 0.95,
             cpp_options = list(PRECOMPILED_HEADERS = "false"))
