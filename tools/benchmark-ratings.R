# Source after loading BayesMFRM and fitting the bundled writing example.
benchmark_ratings <- function(fit, rounds = 9L, calls = 5L) {
  data("bmfrm_example", package = "BayesMFRM", envir = environment())
  ratings <- bmfrm_example$ratings
  model <- score ~ person + item + rater
  fixed <- simulate_bmfrm(model, ratings, 4, seed = 47)$parameters[[1L]]
  rubrics <- c(Content = 2, Evidence = 4, Structure = 6)
  jobs <- list(
    agreement = function() rater_agreement(fit),
    fixed_rsm = function() simulate_bmfrm(model, ratings, 4, parameters = fixed, nsim = 200, seed = 48),
    prior_rsm = function() pp_bmfrm(model, ratings, 4, ndraws = 200, seed = 48),
    prior_pcm = function() pp_bmfrm(model, ratings, rubrics, family = "pcm", threshold_facet = "item",
                                    ndraws = 200, seed = 48)
  )
  vapply(jobs, function(job) {
    invisible(job())
    elapsed <- replicate(rounds, system.time(for (i in seq_len(calls)) job())[["elapsed"]] / calls)
    median(elapsed)
  }, numeric(1L))
}
