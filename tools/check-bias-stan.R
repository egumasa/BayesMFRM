# Fixed-point density and analytic gradient checks. No MCMC or package probability
# / basis helpers in the reference. Includes the positive-SD Jacobian explicitly.
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value=TRUE)[1L])
root <- dirname(dirname(normalizePath(script)))
suppressPackageStartupMessages(pkgload::load_all(root,quiet=TRUE))
out <- file.path(root,".cache","bias-parity")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
H <- function(n, scale="raw") {
  b <- matrix(0,n,max(0L,n-1L))
  if (n>1L) for (j in seq_len(n-1L)) {
    b[seq_len(j),j] <- -1/sqrt(j*(j+1)); b[j+1L,j] <- j/sqrt(j*(j+1))
  }
  if (scale=="marginal" && n>1L) b <- b*sqrt(n/(n-1))
  b
}
scale_density <- function(h, kind) {
  s <- exp(h)
  switch(kind,
    normal=list(lp=-.5*(s/.6)^2+h,g=1-(s/.6)^2),
    student_t=list(lp=stats::dt(s/.6,4,log=TRUE)-log(.6)+h,g=1-5*s^2/(4*.6^2+s^2)),
    cauchy=list(lp=stats::dcauchy(s,0,.6,log=TRUE)+h,g=1-2*s^2/(.6^2+s^2)),
    lognormal=list(lp=stats::dnorm(h,-.4,.5,log=TRUE),g=-(h+.4)/.5^2))
}
reference <- function(q,d,cs) {
  raw <- isTRUE(cs$raw)
  pos <- 0L
  take <- function(n) { i <- pos+seq_len(n); pos <<- pos+n; i }
  facets <- c("item","rater",if (isTRUE(cs$multiple)) "task")
  ip <- take(d$J_person)
  im <- if (!raw) take(d$estimate_mu) else integer()
  il <- if (!raw) take(d$estimate_sigma) else integer()
  ids <- lapply(facets,function(f) take(d[[paste0("J_",f)]]-as.integer(!raw))); names(ids) <- facets
  Js <- if (cs$family=="pcm") d$J_item else 1L
  it <- take(Js*(d$K-if(raw) 1L else 2L))
  pairs <- list(c("item","rater")); if (isTRUE(cs$multiple)) pairs[[2L]] <- c("rater","task")
  terms <- lapply(seq_along(pairs),function(j) {
    f <- pairs[[j]]; n <- unlist(d[paste0("J_",f)])
    list(f=f,n=n,i=take(prod(n-as.integer(!raw))),
      ih=if (!raw && j==1L && cs$sd!="fixed") take(1L) else integer())
  })
  stopifnot(length(q)==pos)
  mu <- if (length(im)) d$ability_mean_location+d$ability_mean_scale*q[im] else if(raw) .4 else d$ability_mean_location
  sigma <- if (length(il)) exp(d$ability_sd_location+d$ability_sd_scale*q[il]) else if(raw) .8 else d$ability_sd_location
  theta <- if(raw) q[ip] else mu+sigma*q[ip]
  eta <- theta[d$person_id]
  bases <- lapply(facets,function(f) H(d[[paste0("J_",f)]],cs$scale)); names(bases) <- facets
  for (f in facets) {
    b <- if(raw) q[ids[[f]]]-mean(q[ids[[f]]]) else as.vector(bases[[f]]%*%q[ids[[f]]])
    eta <- eta-b[d[[paste0(f,"_id")]]]
  }
  Bt <- H(d$K-1L,cs$scale)
  tau <- if(raw) {
    z <- matrix(q[it],Js,d$K-1L); z-rowMeans(z)
  } else matrix(q[it],Js,d$K-2L)%*%t(Bt)
  lp_prior <- -sum(q^2)/2; g <- -q
  if(raw) {
    lp_prior <- sum(stats::dnorm(q[ip],.4,.8,log=TRUE))+sum(stats::dnorm(q[unlist(ids)],log=TRUE))+
      sum(stats::dt(q[it]/.5,4,log=TRUE)-log(.5))
    g[ip] <- -(q[ip]-.4)/.8^2
    g[it] <- -5*q[it]/(4*.5^2+q[it]^2)
  }
  for (j in seq_along(terms)) {
    b <- terms[[j]]
    b$left <- H(b$n[1L],cs$scale); b$right <- H(b$n[2L],cs$scale)
    b$s <- if(length(b$ih)) exp(q[b$ih]) else .6
    z <- matrix(q[b$i],b$n[1L]-as.integer(!raw),b$n[2L]-as.integer(!raw))
    if(raw) {
      b$effect <- sweep(z-rowMeans(z),2L,colMeans(z-rowMeans(z)),"-")
      if(cs$dist=="normal") {
        lp_prior <- lp_prior+sum(stats::dnorm(q[b$i],0,.6,log=TRUE)); g[b$i] <- -q[b$i]/.6^2
      } else if(cs$dist=="student_t") {
        lp_prior <- lp_prior+sum(stats::dt(q[b$i]/.6,4,log=TRUE)-log(.6)); g[b$i] <- -5*q[b$i]/(4*.6^2+q[b$i]^2)
      } else {
        lp_prior <- lp_prior-sum(abs(q[b$i])/.6)-length(b$i)*log(1.2); g[b$i] <- -sign(q[b$i])/.6
      }
    } else {
      b$effect <- b$s*b$left%*%z%*%t(b$right)
      if(length(b$ih)) {
        density <- scale_density(q[b$ih],cs$sd)
        lp_prior <- lp_prior+q[b$ih]^2/2+density$lp; g[b$ih] <- density$g
      }
    }
    b$row <- d[[paste0(b$f[1L],"_id")]]; b$col <- d[[paste0(b$f[2L],"_id")]]
    eta <- eta+b$effect[cbind(b$row,b$col)]; terms[[j]] <- b
  }
  si <- if(cs$family=="pcm") d$item_id else rep(1L,d$N)
  logits <- outer(eta,0:(d$K-1L))-t(vapply(si,function(j) c(0,cumsum(tau[j,])),numeric(d$K)))
  logits <- logits-apply(logits,1L,max)
  pr <- exp(logits)/rowSums(exp(logits))
  ll <- logits[cbind(seq_len(d$N),d$X)]-log(rowSums(exp(logits)))
  score <- d$X-1L-as.vector(pr%*%(0:(d$K-1L)))
  sums <- function(id,n) vapply(seq_len(n),function(j) sum(score[id==j]),numeric(1))
  g[ip] <- g[ip]+if(raw) sums(d$person_id,d$J_person) else sigma*sums(d$person_id,d$J_person)
  if(length(im)) g[im] <- g[im]+d$ability_mean_scale*sum(score)
  if(length(il)) g[il] <- g[il]+d$ability_sd_scale*sigma*sum(score*q[ip][d$person_id])
  for(f in facets) {
    gs <- sums(d[[paste0(f,"_id")]],d[[paste0("J_",f)]])
    g[ids[[f]]] <- g[ids[[f]]]-if(raw) gs-mean(gs) else as.vector(crossprod(bases[[f]],gs))
  }
  gt <- matrix(0,Js,d$K-1L)
  for(h in seq_len(d$K-1L)) for(j in seq_len(Js)) {
    gt[j,h] <- sum((rowSums(pr[,(h+1L):d$K,drop=FALSE])-(d$X>h))[si==j])
  }
  g[it] <- g[it]+if(raw) as.vector(gt-rowMeans(gt)) else as.vector(gt%*%Bt)
  for(b in terms) {
    gb <- matrix(0,b$n[1L],b$n[2L])
    for(n in seq_len(d$N)) gb[b$row[n],b$col[n]] <- gb[b$row[n],b$col[n]]+score[n]
    g[b$i] <- g[b$i]+if(raw) as.vector(sweep(gb-rowMeans(gb),2L,colMeans(gb-rowMeans(gb)),"-")) else
      as.vector(b$s*crossprod(b$left,gb)%*%b$right)
    if(length(b$ih)) g[b$ih] <- g[b$ih]+sum(gb*b$effect)
  }
  list(lp=sum(ll)+lp_prior,gradient=g)
}
cases <- list(
  rsm_half_normal=list(sd="normal"), pcm_half_normal=list(sd="normal",family="pcm"),
  rsm_half_t=list(sd="student_t"), rsm_half_cauchy=list(sd="cauchy"), rsm_lognormal=list(sd="lognormal"),
  rsm_fixed_sd=list(sd="fixed"), rsm_raw_scale=list(sd="normal",scale="raw"),
  pcm_singleton_binary=list(sd="normal",family="pcm",Jr=1L,K=2L),
  rsm_multiple=list(sd="normal",multiple=TRUE),
  legacy_normal=list(raw=TRUE,dist="normal"), legacy_student_t=list(raw=TRUE,dist="student_t",family="pcm"),
  legacy_laplace=list(raw=TRUE,dist="double_exponential"))
rows <- list()
for(name in names(cases)) {
  cs <- utils::modifyList(list(sd="fixed",raw=FALSE,scale="marginal",family="rsm",Jr=2L,K=4L,multiple=FALSE),cases[[name]])
  d <- expand.grid(person=1:6,item=1:3,rater=seq_len(cs$Jr))
  if(cs$multiple) d <- merge(d,data.frame(task=1:2))
  facets <- c("item","rater",if(cs$multiple) "task")
  terms <- c("item:rater",if(cs$multiple) "rater:task")
  d$score <- (seq_len(nrow(d))-1L)%%cs$K+1L
  f <- stats::reformulate(c("person",facets,terms),response="score")
  hier <- if(cs$raw) character() else c("person","rater:item")
  ps <- if(cs$raw) c(prior("normal(.4,.8)","theta")) else c(prior("normal(0,2.5)","mean",facet="person"),prior("lognormal(0,.5)","sd",facet="person"))
  for(facet in facets) ps <- c(ps,list(prior("normal(0,1)",if(facet%in%c("item","rater")) facet else paste0("facet_",facet),scale=if(cs$raw) "raw" else cs$scale)))
  ps <- c(ps,list(prior(if(cs$raw) "student_t(4,0,.5)" else "normal(0,1)","tau",scale=if(cs$raw) "raw" else cs$scale)))
  sd <- switch(cs$sd,normal="normal(0,.6)",student_t="student_t(4,0,.6)",cauchy="cauchy(0,.6)",lognormal="lognormal(-.4,.5)",fixed=.6)
  bp <- if(cs$raw) prior(switch(cs$dist,normal="normal(0,.6)",student_t="student_t(4,0,.6)",double_exponential="double_exponential(0,.6)"),"bias",facet="rater:item") else
    prior(sd,"sd",facet="rater:item",scale=cs$scale)
  ps <- c(ps,list(bp))
  if(cs$multiple) ps <- c(ps,list(prior("normal(0,.6)","bias",facet="task:rater",scale=cs$scale)))
  spec <- parse_bmfrm_formula(f,d,cs$K,family=cs$family,threshold_facet=if(cs$family=="pcm") "item" else NULL,hierarchical=hier)
  data <- prepare_data_bmfrm(spec,d,ps)$stan_data
  compiled <- compile_bmfrm(build_stan_code(spec,ps),"bias",out,cpp_options=list(PRECOMPILED_HEADERS="false"))
  np <- data$J_person+(if(cs$raw) 0L else data$estimate_mu+data$estimate_sigma)+
    sum(unlist(data[paste0("J_",facets)])-as.integer(!cs$raw))+
    (if(cs$family=="pcm") data$J_item else 1L)*(data$K-if(cs$raw) 1L else 2L)+
    (data$J_item-as.integer(!cs$raw))*(data$J_rater-as.integer(!cs$raw))+
    as.integer(!cs$raw&&cs$sd!="fixed")+
    if(cs$multiple) (data$J_rater-1L)*(data$J_task-1L) else 0L
  # Avoid Laplace cusps for derivative checks, including the first reference point.
  points <- rbind(rep(.11,np),.3*sin(seq_len(np)),.7*cos(seq_len(np)),seq(-.8,.8,length.out=np)+.013)
  dp <- file.path(out,paste0(name,"-data.json")); cmdstanr::write_stan_json(data,dp)
  actual <- do.call(rbind,lapply(seq_len(nrow(points)),function(j) {
    qp <- file.path(out,paste0(name,"-",j,".json")); op <- file.path(out,paste0(name,"-",j,".csv"))
    cmdstanr::write_stan_json(list(params_r=as.numeric(points[j,])),qp)
    code <- system2(compiled$model$exe_file(),c("log_prob","jacobian=1",paste0("unconstrained_params=",shQuote(qp)),
      "data",paste0("file=",shQuote(dp)),"output",paste0("file=",shQuote(op)),"sig_figs=18"),stdout=file.path(out,paste0(name,"-",j,".log")),stderr="")
    stopifnot(code==0L); utils::read.csv(op,comment.char="#",check.names=FALSE)
  }))
  ref <- lapply(seq_len(nrow(points)),function(j) reference(points[j,],data,cs))
  lp <- vapply(ref,`[[`,numeric(1),"lp")
  de <- max(abs((actual[[1]]-actual[[1]][1])-(lp-lp[1])))
  ge <- max(abs(as.matrix(actual[,-1])-do.call(rbind,lapply(ref,`[[`,"gradient"))))
  cat(name,"density difference error:",de,"gradient error:",ge,"\n")
  stopifnot(de<1e-10,ge<1e-10)
  rows[[name]] <- data.frame(case=name,dimensions=np,points=4,density_difference_error=de,gradient_error=ge)
}
utils::write.csv(do.call(rbind,rows),file.path(out,"parity.csv"),row.names=FALSE)
