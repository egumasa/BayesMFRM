// Rating-scale MFRM with 4 facets:
// Examinee, Item(Criterion), Rater, Interlocutor
data {
  int<lower=1> J;                // # examinees
  int<lower=1> I;                // # items / criteria
  int<lower=1> R;                // # raters
  int<lower=1> L;                // # interlocutors
  int<lower=2> K;                // # rating categories
  int<lower=1> N;                // # observations

  array[N] int<lower=1, upper=J> ExamineeID;
  array[N] int<lower=1, upper=I> ItemID;
  array[N] int<lower=1, upper=R> RaterID;
  array[N] int<lower=1, upper=L> InterlocutorID;
  array[N] int<lower=1, upper=K> X;   // observed category (1..K)
}

parameters {
  vector[J] theta_raw;        // examinee ability
  vector[I] beta_raw;         // item difficulty
  vector[R] rater_raw;        // rater severity
  vector[L] interlocutor_raw; // interlocutor effect
  vector[K - 1] tau_raw;      // rating-scale thresholds (unconstrained)
}

transformed parameters {
  vector[J] theta;
  vector[I] beta;
  vector[R] rater;
  vector[L] interlocutor;
  vector[K - 1] tau;

  // sum-to-zero identification (Eckes-style)
  theta        = theta_raw;        // To mimic FACET behavior this should not be centered  - mean(theta_raw);
  beta         = beta_raw         - mean(beta_raw);
  rater        = rater_raw        - mean(rater_raw);
  interlocutor = interlocutor_raw - mean(interlocutor_raw);
  tau          = tau_raw          - mean(tau_raw);
}

model {
  // Priors
  theta_raw        ~ student_t(3, 0, 4);
  beta_raw         ~ student_t(3, 0, 2);
  rater_raw        ~ student_t(3, 0, 2);
  interlocutor_raw ~ student_t(3, 0, 2);
  tau_raw          ~ student_t(3, 0, 2);

  // Likelihood: Many-Facet Rating Scale Model
  for (n in 1:N) {
    int j = ExamineeID[n];
    int i = ItemID[n];
    int r = RaterID[n];
    int l = InterlocutorID[n];
    int x = X[n];

    // core Rasch/MFRM equation
    real eta = theta[j] - beta[i] - rater[r] - interlocutor[l];

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