// Rating-scale MFRM with 3 facets:
// person, rater, scenario
data {
  int<lower=1> J_person;         // # persons
  int<lower=1> J_rater;          // # raters
  int<lower=1> J_scenario;       // # scenarios
  int<lower=2> K;                // # rating categories
  int<lower=1> N;                // # observations

  array[N] int<lower=1, upper=J_person> person_id;
  array[N] int<lower=1, upper=J_rater> rater_id;
  array[N] int<lower=1, upper=J_scenario> scenario_id;
  array[N] int<lower=1, upper=K> X;   // observed category (1..K)
}

parameters {
  vector[J_person] theta_raw;        // person ability
  vector[J_rater] rater_raw;         // rater severity
  vector[J_scenario] scenario_raw;   // scenario effect
  vector[K - 1] tau_raw;             // rating-scale thresholds (unconstrained)
}

transformed parameters {
  vector[J_person] theta;
  vector[J_rater] rater;
  vector[J_scenario] scenario;
  vector[K - 1] tau;

  // sum-to-zero identification (Eckes-style)
  theta        = theta_raw;        // To mimic FACET behavior this should not be centered  - mean(theta_raw);
  rater        = rater_raw        - mean(rater_raw);
  scenario     = scenario_raw     - mean(scenario_raw);
  tau          = tau_raw          - mean(tau_raw);
}

model {
  // Priors
  theta_raw        ~ student_t(3, 0, 4);
  rater_raw        ~ student_t(3, 0, 2);
  scenario_raw     ~ student_t(3, 0, 2);
  tau_raw          ~ student_t(3, 0, 2);

  // Likelihood: Many-Facet Rating Scale Model
  for (n in 1:N) {
    int j = person_id[n];
    int r = rater_id[n];
    int s = scenario_id[n];
    int x = X[n];

    // core Rasch/MFRM equation
    real eta = theta[j] - rater[r] - scenario[s];

    vector[K] logits;
    logits[1] = 0; // base category

    // successive sums (Andrich rating-scale structure)
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
  vector[N] mu_fair;   // expected score with "fair" rater & scenario

  for (n in 1:N) {
    int j = person_id[n];
    int r = rater_id[n];
    int s = scenario_id[n];
    int x = X[n];

    // ---- 1) Actual facets ----
    real eta = theta[j] - rater[r] - scenario[s];

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

    // ---- 2) "Fair" expected score: rater = 0, scenario = 0 ----
    {
      real eta_fair = theta[j];  // average rater & scenario = 0
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
