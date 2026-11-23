bias_fixture <- function(h = c("person", "rater:item"), scale = "marginal", sd = "normal(0, .5)") {
  d <- expand.grid(person = paste0("P", 1:6), item = paste0("I", 1:3), rater = c("A", "B"))
  d$score <- 1L
  f <- score ~ person + item + rater + rater:item
  s <- parse_bmfrm_formula(f, d, K = 4, hierarchical = h)
  ps <- get_default_priors(s$facets_main, s$facets_bias, hierarchical = h)
  at <- which(vapply(ps, function(p) !is.null(p$facet) && grepl(":", p$facet), logical(1)))
  ps[[at]] <- prior(sd, if (length(intersect(h, c("item:rater", "rater:item")))) "sd" else "bias", facet = "rater:item", scale = scale)
  list(data = d, spec = s, priors = ps)
}

test_that("interaction names are matched independent of component order", {
  x <- bias_fixture()
  s <- resolve_bmfrm_spec(x$spec, x$priors)
  expect_identical(s$hierarchical, c("person", "item:rater"))
  ps <- x$priors
  at <- which(vapply(ps, function(p) identical(p$facet, "rater:item"), logical(1)))
  ps[[at]]$facet <- "item:rater"
  expect_identical(build_stan_code(s), build_stan_code(x$spec, ps))
  expect_error(resolve_bmfrm_spec(x$spec, c(ps, list(x$priors[[at]]))), "Duplicate")
  expect_error(validate_hierarchical(c("item:rater", "rater:item"), s$facets_main, s$facets_bias, "rsm"), "unique")
  expect_error(parse_bmfrm_formula(score ~ person + item + rater + person:item:rater, x$data, 4), "two-way")
  expect_error(validate_interaction_terms(c("person", "item", "rater", "itemrater"), "item:rater"), "collide")
  expect_error(prior("normal(0,1)", "mean", facet = "item:rater"), "no population mean")
  expect_error(prior("normal(1,1)", "sd", facet = "item:rater"), "requires")
  di <- prepare_data_bmfrm(s,x$data)
  expect_identical(get_facet_labels(di,"rater:item"),paste("item",rep(paste0("I",1:3),2),"rater",rep(c("A","B"),each=3),sep=":"))
})

test_that("interaction covariance and hierarchy respect the declared scale", {
  for (scale in c("raw", "marginal")) {
    x <- bias_fixture(sd = .7, scale = scale)
    s <- resolve_bmfrm_spec(x$spec, x$priors)
    d <- prepare_data_bmfrm(s, x$data)$stan_data
    M <- kronecker(d$B_itemrater_2, d$B_itemrater_1)
    V <- .7^2 * tcrossprod(M)
    expected <- .7^2 * if (scale == "raw") (1-1/3)*(1-1/2) else 1
    expect_equal(diag(V), rep(expected, 6))
    expect_equal(V[1, 2], -expected/2)
    expect_equal(V[1, 4], -expected)
    expect_equal(V[1, 5], expected/2)
    expect_equal(qr(M)$rank, 2)
  }
  x <- bias_fixture()
  sim <- prior_predict_bmfrm(x$spec$formula, x$data, 4, priors = x$priors,
    hierarchical = x$spec$hierarchical, ndraws = 2500, seed = 96)
  sigma <- vapply(sim$parameters, `[[`, numeric(1), "sigma_itemrater")
  b <- vapply(sim$parameters, function(p) p$itemrater[1, 1], numeric(1))
  expect_equal(mean(sigma^2), .25, tolerance = .025)
  expect_equal(stats::var(b), .25, tolerance = .035)
  expect_equal(stats::var(b / sigma), 1, tolerance = .08)
  expect_true(all(vapply(sim$parameters, function(p) max(abs(c(rowSums(p$itemrater), colSums(p$itemrater)))) < 1e-12, logical(1))))
  s <- sim$spec; n <- c(3L, 2L); Q <- rsm_prior_basis(3, "marginal"); T <- rsm_prior_basis(2, "marginal")
  z <- matrix(c(.3, -.8), 2, 1); b <- Q %*% z %*% t(T)
  perm <- c(3,1,2); other <- 2:1
  transported <- crossprod(Q, b[perm, other]) %*% T * (2/3) * (1/2)
  expect_equal(sum(z^2), sum(transported^2), tolerance = 1e-14)
  expect_equal(Q %*% transported %*% t(T), b[perm, other], tolerance = 1e-14)
})

test_that("interaction scale priors and singleton blocks are explicit", {
  for (p in c("normal(0,.5)", "student_t(4,0,.5)", "cauchy(0,.5)", "lognormal(-1,.3)")) {
    x <- bias_fixture(sd = p)
    sim <- simulate_bmfrm(x$spec$formula, x$data, 4, priors = x$priors, hierarchical = x$spec$hierarchical, seed = 12)
    expect_gt(sim$parameters[[1]]$sigma_itemrater, 0)
    expect_match(build_stan_code(sim$spec), paste0("sigma_itemrater_raw ~ ", p, ";"), fixed = TRUE)
  }
  x <- bias_fixture(sd = .7)
  x$data$rater <- "only"
  sim <- simulate_bmfrm(x$spec$formula, x$data, 2, priors = x$priors, hierarchical = x$spec$hierarchical, seed = 12)
  expect_equal(sim$parameters[[1]]$itemrater, matrix(0, 3, 1))
  expect_equal(sim$parameters[[1]]$sigma_itemrater, .7)
  vars <- model_diagnostic_variables(sim$spec, prepare_data_bmfrm(sim$spec, sim$data))
  expect_false(any(c("itemrater", "itemrater_raw", "sigma_itemrater", "sigma_itemrater_raw") %in% vars))
})

test_that("fixed bias matrices preserve labels and respect sign and constraints", {
  x <- bias_fixture(h = character(), scale = "raw", sd = "normal(0,.5)")
  b <- cbind(c(-.7,.2,.5), c(.7,-.2,-.5))
  pars <- list(theta = seq(-1,1,length.out=6), item=c(-.2,0,.2), rater=c(-.1,.1), tau=c(-.4,0,.4), itemrater=b)
  sim <- simulate_bmfrm(x$spec$formula,x$data,4,parameters=pars,nsim=2,seed=18)
  # Independent categorical probability at one planned cell.
  eta <- pars$theta[1]-pars$item[1]-pars$rater[1]+b[1,1]
  logits <- c(0,cumsum(eta-pars$tau)); pr <- exp(logits)/sum(exp(logits))
  expect_equal(as.vector(rsm_probabilities(eta,pars$tau)),pr)
  dimnames(pars$itemrater) <- list(paste0("I",1:3),c("A","B"))
  pars$itemrater <- pars$itemrater[3:1,2:1]
  expect_identical(sim$y_rep, simulate_bmfrm(x$spec$formula,x$data,4,parameters=pars,nsim=2,seed=18)$y_rep)
  pars$itemrater[1,1] <- pars$itemrater[1,1]+1
  expect_error(simulate_bmfrm(x$spec$formula,x$data,4,parameters=pars),"each sum to zero")
})

test_that("rank checks include interactions and preserve unobserved cells", {
  x <- bias_fixture()
  ch <- check_design_bmfrm(x$spec$formula,x$data,4)
  expect_equal(ch$n_observations,36)
  expect_equal(ch$n_parameters,11)
  expect_true(ch$identifiable)
  d <- subset(x$data,!(item=="I3" & rater=="B"))
  bad <- check_design_bmfrm(x$spec$formula,d,4)
  expect_true(bad$connected)
  expect_equal(bad$nullity,1)
  expect_equal(sum(bad$interaction_levels[[1]]$n_observations==0),1)
  expect_lt(max(abs(bad$design_matrix %*% bad$null_space)),1e-12)
  # Standard treatment-coded model has the same column space on these rows.
  ref <- stats::model.matrix(~ factor(person)+factor(item)*factor(rater),d)
  expect_equal(qr(ref)$rank,bad$rank)
  expect_equal(qr(cbind(ref,bad$design_matrix))$rank,bad$rank)
  pcm <- check_design_bmfrm(x$spec$formula,d,4,family="pcm",threshold_facet="item")
  expect_equal(pcm$adjacent_logit_rank,bad$rank+6)
  expect_error(check_design_bmfrm(x$spec$formula,x$data,4,max_cells=10),"limit")
})

test_that("bias summaries use generated names, original labels and practical thresholds", {
  x <- bias_fixture()
  ch <- check_design_bmfrm(x$spec$formula,x$data,4)
  b <- c(-.5,0,.5,.5,0,-.5)
  mat <- matrix(rep(b,each=50),50,6,dimnames=list(NULL,paste0("itemrater[",rep(1:3,2),",",rep(1:2,each=3),"]")))
  obj <- structure(list(spec=x$spec,design=ch,
    data_info=list(facet_labels=list(person=paste0("P",1:6),item=paste0("I",1:3),rater=c("A","B")),data_clean=x$data),
    fit=list(draws=function(...) mat[,6:1])),class="bmfrm_fit")
  a <- summarise_bias(obj,"rater:item",delta=.2)
  expect_identical(a$item,rep(paste0("I",1:3),2))
  expect_identical(a$rater,rep(c("A","B"),each=3))
  expect_equal(unname(a$bias_mean),b)
  expect_equal(unname(a$prob_outside_delta),as.numeric(abs(b)>.2))
  expect_true(all(a$likelihood_identified))
  expect_true(all(a$n_observations==6))
  expect_equal(sum(a$flag_bias=="positive_bias"),2)
  expect_equal(sum(summarise_bias(obj,"item:rater",delta=.6)$flag_bias=="no_clear_bias"),6)
  obj$data_info$data_clean <- subset(x$data,!(item=="I3" & rater=="B"))
  obj$design <- check_design_bmfrm(x$spec$formula,obj$data_info$data_clean,4)
  a <- summarise_bias(obj,"item:rater",delta=.2)
  expect_true(all(!a$likelihood_identified))
  expect_true(all(a$flag_bias=="not_likelihood_identified"))
  expect_equal(a$n_observations[6],0)
  expect_error(summarise_bias(obj,delta=-1),"nonnegative")
})

test_that("fair predictions retain interactions among unadjusted facets", {
  d <- expand.grid(person=1:4,item=1:3,rater=1:2,task=1:2); d$score <- 1
  s <- parse_bmfrm_formula(score~person+item+rater+task+item:rater+item:task,d,3)
  code <- generate_eta_sum(s,fair=TRUE)
  expect_match(code,"itemtask[",fixed=TRUE)
  expect_false(grepl("itemrater[",code,fixed=TRUE))
  expect_match(build_stan_code(s),"matrix[J_item, J_rater] itemrater;",fixed=TRUE)
})
