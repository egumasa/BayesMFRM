// four-facets-main-effect-hierarchical-model5-full.stan
// 4-facet hierarchical rating-scale MFRM with:
//   - Main effects: Examinee, Item, Rater, Interlocutor, Thresholds
//   - Two-way interactions (bias):
//       delta_er[j, r] : Examinee × Rater
//       delta_ei[j, l] : Examinee × Interlocutor
//       delta_ri[r, l] : Rater × Interlocutor
//   - Three-way interaction (bias):
//       delta_eri[jl, r] : Examinee × Rater × Interlocutor
//         where jl = (j - 1) * L + l

data {
  int<lower=1> J;                // # examinees
  int<lower=1> I;                // # items / criteria
  int<lower=1> R;                // # raters
  int<lower=1> L;                // # interlocutors
  int<lower=2> K;                // # categories
  int<lower=1> N;                // # observations

  array[N] int<lower=1, upper=J> ExamineeID;
  array[N] int<lower=1, upper=I> ItemID;
  array[N] int<lower=1, upper=R> RaterID;
  array[N] int<lower=1, upper=L> InterlocutorID;
  array[N] int<lower=1, upper=K> X;   // observed category (1..K)
}

parameters {
  // ─────────────────────────────────────────────────────────────
  // Individual-level (raw) parameters
  vector[J]       theta_raw;        // examinee ability (z)
  vector[I]       beta_raw;         // item difficulty (z)
  vector[R]       rater_raw;        // rater severity (z)
  vector[L]       interlocutor_raw; // interlocutor effect (z)
  vector[K - 1]   tau_raw;          // rating-scale thresholds (z)

  // Two-way interactions (bias) – raw z-scores
  vector[J * R]   delta_er_raw;     // Examinee × Rater
  vector[J * L]   delta_ei_raw;     // Examinee × Interlocutor
  vector[R * L]   delta_ri_raw;     // Rater × Interlocutor

  // Three-way interaction (bias) – raw z-scores
  // stored as (J * L) × R after reshaping
  vector[J * R * L] delta_eri_raw;  // Examinee × Rater × Interlocutor

  // ─────────────────────────────────────────────────────────────
  // Hyperparameters (population-level SDs)
  real<lower=0> sigma_theta;
  real<lower=0> sigma_beta;
  real<lower=0> sigma_rater;
  real<lower=0> sigma_interlocutor;
  real<lower=0> sigma_tau;

  real<lower=0> sigma_delta_er;     // SD of Examinee × Rater bias
  real<lower=0> sigma_delta_ei;     // SD of Examinee × Interlocutor bias
  real<lower=0> sigma_delta_ri;     // SD of Rater × Interlocutor bias
  real<lower=0> sigma_delta_eri;    // SD of 3-way bias
}

transformed parameters {
  // ─────────────────────────────────────────────────────────────
  // Main effects on actual scale
  vector[J]     theta;
  vector[I]     beta;
  vector[R]     rater;
  vector[L]     interlocutor;
  vector[K - 1] tau;

  // Interaction terms on actual scale
  matrix[J, R]       delta_er;      // Examinee × Rater bias
  matrix[J, L]       delta_ei;      // Examinee × Interlocutor bias
  matrix[R, L]       delta_ri;      // Rater × Interlocutor bias
  matrix[J * L, R]   delta_eri;     // 3-way bias: row jl = (j-1)*L + l, col r

  // Main effects:
  theta = sigma_theta * theta_raw;   // persons free to drift; SD identifies scale
  beta  = sigma_beta  * (beta_raw         - mean(beta_raw));
  rater = sigma_rater * (rater_raw        - mean(rater_raw));
  interlocutor = sigma_interlocutor * (interlocutor_raw - mean(interlocutor_raw));
  tau   = sigma_tau   * (tau_raw          - mean(tau_raw));

  // ─────────────────────────────────────────────────────────────
  // Two-way interactions: centering → sum-to-zero by row

  {
    // Examinee × Rater: matrix[J, R]
    matrix[J, R] tmp_er = to_matrix(delta_er_raw, J, R);
    for (j in 1:J) {
      real mn = mean(tmp_er[j]);
      for (r in 1:R) {
        delta_er[j, r] = sigma_delta_er * (tmp_er[j, r] - mn);
      }
    }
  }

  {
    // Examinee × Interlocutor: matrix[J, L]
    matrix[J, L] tmp_ei = to_matrix(delta_ei_raw, J, L);
    for (j in 1:J) {
      real mn = mean(tmp_ei[j]);
      for (l in 1:L) {
        delta_ei[j, l] = sigma_delta_ei * (tmp_ei[j, l] - mn);
      }
    }
  }

  {
    // Rater × Interlocutor: matrix[R, L]
    matrix[R, L] tmp_ri = to_matrix(delta_ri_raw, R, L);
    for (r in 1:R) {
      real mn = mean(tmp_ri[r]);
      for (l in 1:L) {
        delta_ri[r, l] = sigma_delta_ri * (tmp_ri[r, l] - mn);
      }
    }
  }

  {
    // Three-way: store as matrix[J*L, R]
    // Each row s = (j-1)*L + l corresponds to (j, l), cols = raters
    matrix[J * L, R] tmp_eri = to_matrix(delta_eri_raw, J * L, R);

    for (s in 1:(J * L)) {
      real mn = mean(tmp_eri[s]);
      for (r in 1:R) {
        delta_eri[s, r] = sigma_delta_eri * (tmp_eri[s, r] - mn);
      }
    }
    // Interpretation: for each (j, l), sum_r delta_eri[row(s), r] = 0
  }
}

model {
  // ─────────────────────────────────────────────────────────────
  // Hyperpriors
  sigma_theta        ~ normal(0, 3);
  sigma_beta         ~ normal(0, 2);
  sigma_rater        ~ normal(0, 2);
  sigma_interlocutor ~ normal(0, 2);
  sigma_tau          ~ normal(0, 2);

  sigma_delta_er     ~ normal(0, 1);
  sigma_delta_ei     ~ normal(0, 1);
  sigma_delta_ri     ~ normal(0, 1);
  sigma_delta_eri    ~ normal(0, 1);

  // Priors on standardized individual-level parameters
  theta_raw         ~ student_t(3, 0, 1);
  beta_raw          ~ student_t(3, 0, 1);
  rater_raw         ~ student_t(3, 0, 1);
  interlocutor_raw  ~ student_t(3, 0, 1);
  tau_raw           ~ student_t(3, 0, 1);

  delta_er_raw      ~ student_t(3, 0, 1);
  delta_ei_raw      ~ student_t(3, 0, 1);
  delta_ri_raw      ~ student_t(3, 0, 1);
  delta_eri_raw     ~ student_t(3, 0, 1);

  // ─────────────────────────────────────────────────────────────
  // Likelihood
  for (n in 1:N) {
    int j = ExamineeID[n];
    int i = ItemID[n];
    int r = RaterID[n];
    int l = InterlocutorID[n];
    int x = X[n];

    int jl = (j - 1) * L + l;  // row index for (j, l) combination in delta_eri

    // Linear predictor with all main + interaction effects
    real eta = theta[j]
               - beta[i]
               - rater[r]
               - interlocutor[l]
               + delta_er[j, r]
               + delta_ei[j, l]
               + delta_ri[r, l]
               + delta_eri[jl, r];

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
  vector[N] mu;        // expected score under actual facets + bias
  vector[N] sigma2;    // model variance
  vector[N] mu_fair;   // "fair" expected score (no rater/interlocutor/bias)

  for (n in 1:N) {
    int j = ExamineeID[n];
    int i = ItemID[n];
    int r = RaterID[n];
    int l = InterlocutorID[n];
    int x = X[n];

    int jl = (j - 1) * L + l;

    // ── 1) Actual facets + all bias terms
    real eta = theta[j]
               - beta[i]
               - rater[r]
               - interlocutor[l]
               + delta_er[j, r]
               + delta_ei[j, l]
               + delta_ri[r, l]
               + delta_eri[jl, r];

    vector[K] logits;
    vector[K] p;
    real ex;
    real ex2;

    logits[1] = 0;
    for (k in 2:K) {
      logits[k] = logits[k - 1] + eta - tau[k - 1];
    }

    p = softmax(logits);

    ex  = 0;
    ex2 = 0;
    for (k in 1:K) {
      ex  += k * p[k];
      ex2 += k * k * p[k];
    }
    mu[n]     = ex;
    sigma2[n] = ex2 - ex * ex;

    log_lik[n] = categorical_logit_lpmf(x | logits);
    y_rep[n]   = categorical_logit_rng(logits);

    // ── 2) "Fair" expected score:
    //     average rater & interlocutor (0), no bias terms
    {
      real eta_fair = theta[j] - beta[i];
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