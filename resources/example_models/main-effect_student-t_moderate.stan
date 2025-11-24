// Rating-scale MFRM with 4 facets:
// Examinee, Item(Criterion), Rater, Interlocutor
data {
  int<lower=1> J;                // # examinees (units at the person level)
  int<lower=1> I;                // # items / criteria (units at item level)
  int<lower=1> R;                // # raters (units at rater level)
  int<lower=1> L;                // # interlocutors (units at interlocutor level)
  int<lower=2> K;                // # rating categories
  int<lower=1> N;                // # observations (person×item×rater×interlocutor entries)

  array[N] int<lower=1, upper=J> ExamineeID;
  array[N] int<lower=1, upper=I> ItemID;
  array[N] int<lower=1, upper=R> RaterID;
  array[N] int<lower=1, upper=L> InterlocutorID;
  array[N] int<lower=1, upper=K> X;   // observed category (1..K)
}

parameters {
  // ─────────────────────────────────────────────────────────────────────
  // Individual‐level parameters (one per unit)
  vector[J]       theta_raw;        // examinee ability z-scores → individual persons
  vector[I]       beta_raw;         // item difficulty z-scores → individual items
  vector[R]       rater_raw;        // rater severity z-scores → individual raters
  vector[L]       interlocutor_raw; // interlocutor effect z-scores → individual interlocutors
  vector[K - 1]   tau_raw;          // rating‐scale thresholds z-scores → thresholds (common across units)

  // ─────────────────────────────────────────────────────────────────────
  // Population‐level (hyper)parameters (modeling the distribution of those individual units)
  real<lower=0> sigma_theta;            // SD of examinee abilities → population spread of persons
  real<lower=0> sigma_beta;             // SD of item difficulties → population spread of items
  real<lower=0> sigma_rater;            // SD of rater severities → population spread of raters
  real<lower=0> sigma_interlocutor;     // SD of interlocutor effects → population spread of interlocutors
  real<lower=0> sigma_tau;              // SD of thresholds → spread of threshold parameters across categories
}

transformed parameters {
  // ─────────────────────────────────────────────────────────────────────
  // These map the raw individual‐level parameters into “actual” scale using the population SDs
  vector[J]     theta;          // actual person abilities
  vector[I]     beta;           // actual item difficulties
  vector[R]     rater;          // actual rater severities
  vector[L]     interlocutor;   // actual interlocutor effects
  vector[K - 1] tau;            // actual thresholds

  // Persons: non-center facet (we leave mean free; population SD matters)
  theta = sigma_theta * theta_raw;

  // Other facets: centered (mean zero) then scaled by population SDs
  beta         = sigma_beta         * (beta_raw         - mean(beta_raw));
  rater        = sigma_rater        * (rater_raw        - mean(rater_raw));
  interlocutor = sigma_interlocutor * (interlocutor_raw - mean(interlocutor_raw));
  tau          = sigma_tau          * (tau_raw          - mean(tau_raw));
}

model {
  // ─────────────────────────────────────────────────────────────────────
  // Population‐level priors (hyper-priors) – they govern the distributions of the individual units
  sigma_theta        ~ normal(0, 3);    // prior on person‐ability SD
  sigma_beta         ~ normal(0, 2);    // prior on item‐difficulty SD
  sigma_rater        ~ normal(0, 2);    // prior on rater‐severity SD
  sigma_interlocutor ~ normal(0, 2);    // prior on interlocutor SD
  sigma_tau          ~ normal(0, 2);    // prior on threshold spread

  // Individual‐level priors – each individual unit parameter is drawn from a common (standardized) distribution
  theta_raw         ~ student_t(3, 0, 1);     // each person’s z‐score
  beta_raw          ~ student_t(3, 0, 1);     // each item’s z‐score
  rater_raw         ~ student_t(3, 0, 1);     // each rater’s z‐score
  interlocutor_raw  ~ student_t(3, 0, 1);     // each interlocutor’s z‐score
  tau_raw           ~ student_t(3, 0, 1);     // each threshold’s z‐score

  // ─────────────────────────────────────────────────────────────────────
  // Likelihood: data model at the observation level (lowest level)
  for (n in 1:N) {
    int j = ExamineeID[n];
    int i = ItemID[n];
    int r = RaterID[n];
    int l = InterlocutorID[n];
    int x = X[n];

    real eta = theta[j] - beta[i] - rater[r] - interlocutor[l];

    vector[K] logits;
    logits[1] = 0;
    for (k in 2:K) {
      logits[k] = logits[k - 1] + eta - tau[k - 1];
    }
    X[n] ~ categorical_logit(logits);
  }
}

generated quantities {
  vector[N] log_lik;
  array[N] int<lower=1, upper=K> y_rep;
  vector[N] mu;        // expected score E[Y_n] under actual facets
  vector[N] sigma2;    // model variance Var(Y_n) under actual facets
  vector[N] mu_fair;   // expected score with "fair" rater & interlocutor

  for (n in 1:N) {
    int j = ExamineeID[n];
    int i = ItemID[n];
    int r = RaterID[n];
    int l = InterlocutorID[n];
    int x = X[n];

    // ---- 1) Actual facets ----
    real eta = theta[j] - beta[i] - rater[r] - interlocutor[l];

    vector[K] logits;
    vector[K] p;
    real ex;
    real ex2;

    // build cumulative logits for rating-scale model
    logits[1] = 0;
    for (k in 2:K) {
      logits[k] = logits[k - 1] + eta - tau[k - 1];
    }

    // probabilities
    p = softmax(logits);

    // expected score and variance on 1..K scale
    ex  = 0;
    ex2 = 0;
    for (k in 1:K) {
      ex  += k * p[k];
      ex2 += k * k * p[k];
    }
    mu[n]     = ex;
    sigma2[n] = ex2 - ex * ex;

    // log-likelihood & posterior predictive
    log_lik[n] = categorical_logit_lpmf(x | logits);
    y_rep[n]   = categorical_logit_rng(logits);

    // ---- 2) "Fair" expected score: rater = 0, interlocutor = 0 ----
    {
      real eta_fair = theta[j] - beta[i];  // average rater & interlocutor = 0
      vector[K] logits_fair;
      vector[K] p_fair;
      real ex_fair;

      logits_fair[1] = 0;
      for (k in 2:K) {
        logits_fair[k] = logits_fair[k - 1] + eta_fair - tau[k - 1];
      }

      p_fair = softmax(logits_fair);

      ex_fair = 0;
      for (k in 1:K) {
        ex_fair += k * p_fair[k];
      }
      mu_fair[n] = ex_fair;
    }
  }
}