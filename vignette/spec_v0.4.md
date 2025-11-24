# BayesMFRM – Design Specification (v0.3)

1. Big picture

Goal:
BayesMFRM is an R package for Bayesian many-facet Rasch models using Stan (via cmdstanr).

It should:
	•	Feel as easy as brms, but specialized for MFRM.
	•	Hide most Stan details, but still allow advanced users to inspect Stan code.
	•	Produce FACETS-style outputs: rater severity, item difficulty, fair averages, infit/outfit, etc.

We start simple:
	•	Rating-scale Rasch model (shared thresholds).
	•	Main facet effects (person, item, rater, etc.).
	•	Optional 2-way bias facets (e.g., rater × item).
	•	A clean, simple prior system.

Later, we can extend to partial credit, centrality, drift, etc.

⸻

1. User-facing workflow

A typical user should write code like:


```r
library(BayesMFRM)

priors <- c(
  prior("normal(0, 2)",       class = "theta"),        # person
  prior("normal(0, 2)",       class = "item"),         # item difficulty
  prior("student_t(3, 0, 1)", class = "rater"),        # rater severity
  prior("normal(0, 1)",       class = "bias", facet = "rater:item"),
  prior("normal(0, 3)",       class = "tau")           # thresholds
)

fit <- bmfrm(
  score ~ person + item + rater + (rater:item),
  data       = ratings,
  K          = 6,
  priors     = priors,
  family     = "rating_scale",
  model_name = "demo_rater_item",
  cache_dir  = "stan_cache",
  refit      = "on_change",
  iter       = 2000,
  warmup     = 1000,
  chains     = 4,
  cores      = 4,
  seed       = 1234
)

summary(fit, facets = c("person", "rater"))
rater_tab  <- facet_summary(fit, "rater")
fair_scores <- fair_scores(fit)

bias    <- summarise_bias(fit) # when bias term is estimated
resid_tbl  <- extract_residuals(fit) #return FACET like residuals
crit_cor   <- residual_cor_criteria(resid_tbl)
stancode(fit)  # inspect underlying Stan code
```

Key design idea:

User only thinks in terms of facets & Rasch.
The package handles: data preparation → Stan code → sampling → Rasch tables.

⸻

2. Core concepts explained

2.1 What is a “facet”?

In MFRM, facets are “dimensions” of the assessment design:
	•	person / examinee
	•	item or criterion
	•	rater
	•	task
	•	interlocutor
	•	etc.

Each facet is:
	•	a categorical variable in the data (IDs),
	•	with a latent measure on the Rasch logit scale (ability, difficulty, severity).

In R, a facet is just a column:

```
ratings$person     # examinee ID
ratings$item       # item or criterion ID
ratings$rater      # rater ID
ratings$score      # rating 1…K
```

In Stan, we convert these into integer indices 1…J_f (e.g. 1…J_person).

2.2 Main facets vs bias facets
	•	Main facets: person, item, rater, etc.
	•	Bias facets (interactions): combinations like rater:item that represent specific patterns (rater-item bias).

In model formula:

```
score ~ person + item + rater + (rater:item)
```

	•	person, item, rater are main facets.
	•	(rater:item) is an interaction facet (bias).

We treat these differently in Stan:
	•	main facets → vectors of latent parameters.
	•	bias facets → matrices or arrays of parameters.

⸻

3. Statistical model (v0.1)

We start with a rating-scale Rasch bmfrm.

For rating n:
	•	Person: p_n
	•	Item: i_n
	•	Rater: r_n
	•	Score: X_n \in \{1, \dots, K\}

Linear predictor:

\eta_n
= \theta_{p_n}
- \beta_{i_n}
- \rho_{r_n}
+ \text{(optional bias terms)}

Rating-scale structure:

\log\frac{P(X_n \ge k)}{P(X_n < k)} = \eta_n - \tau_{k-1}, \quad k=2,\dots,K
	•	\theta: person ability
	•	\beta: item difficulty
	•	\rho: rater severity
	•	\tau_k: rating thresholds (shared across items/raters in v0.1)

We assume:
	•	equal discrimination (Rasch).
	•	Optional additional facets in the same additive way (e.g., task difficulty).

⸻

4. Top-level API design

4.1 bmfrm() function

```r
bmfrm <- function(
  formula,
  data,
  K,
  priors      = NULL,
  family      = c("rating_scale", "partial_credit"),
  model_name  = NULL,
  cache_dir   = "stan_cache",
  save_csvs = TRUE,         # Always save CSVs for posterior access
  refit       = c("on_change", "never", "always"),
  iter        = 2000,
  warmup      = floor(iter / 2),
  chains      = 4,
  cores       = getOption("mc.cores", 1L),
  seed        = NULL,
  ...
)
```

Arguments
	•	formula
	•	R formula specifying the model:

score ~ person + item + rater + (rater:item)


	•	LHS: score column (integer 1…K).
	•	RHS: facet columns; interaction terms in parentheses for bias facets.

	•	data
	•	Data frame with all variables named in formula.
	•	K
	•	Number of score categories. If NULL, infer from max(score).
	•	priors
	•	Vector of prior() objects (see below). If NULL, use safe defaults.
	•	family
	•	"rating_scale" (v0.1) and later "partial_credit".
	•	model_name
	•	Optional string used to label compilation and caching.
	•	If NULL, automatically generate something like "bmfrm_<hash>".
	•	cache_dir
	•	Directory to store compiled Stan models and optionally cached fits.
	•	refit
	•	"on_change": refit only if Stan code or data structure changed.
	•	"never": never refit (just reuse existing cached fit if present).
	•	"always": always refit.
	•	iter, warmup, chains, cores, seed
	•	Standard cmdstanr::sample() arguments.

4.2 prior() object

We design a small, Rasch-focused prior system.

Usage:

```r
priors <- c(
  prior("normal(0, 2)",       class = "theta"),        # person
  prior("normal(0, 2)",       class = "item"),         # item difficulty
  prior("student_t(3, 0, 1)", class = "rater"),        # rater severity
  prior("normal(0, 1)",       class = "bias", facet = "rater:item"),
  prior("normal(0, 3)",       class = "tau")           # thresholds
)
```

Definition:

```r
prior <- function(spec, class, facet = NULL) {
  structure(
    list(
      spec  = spec,    # e.g. "normal(0,2)"
      class = class,   # e.g. "theta", "item", "rater", "bias", "tau"
      facet = facet    # e.g. "rater:item" (only for bias)
    ),
    class = "bmfrm_prior"
  )
}
```

We restrict class to a small set to keep things simple:
	•	"theta" – person abilities
	•	"item" – item difficulties
	•	"rater" – rater severities
	•	"facet_<name>" – for extra main facets (e.g. "facet_task")
	•	"bias" – interaction (bias) terms; must specify facet = "rater:item" etc.
	•	"tau" – rating thresholds

Design choice (why this is simple):
	•	For a beginner, it’s easier to say “I want a Normal(0, 2) prior on person ability” than to think about abstract “b”/“sd”/“cor” classes.
	•	We define a limited vocabulary that maps directly to psychometric concepts they know.

If priors = NULL, we use defaults such as:
	•	theta ~ normal(0, 2)
	•	item ~ normal(0, 1)
	•	rater ~ normal(0, 1)
	•	bias ~ normal(0, 0.5)
	•	tau ~ normal(0, 3) on ordered thresholds (or on spacings).

pipeline for the bmfrm function

```r
bmfrm <- function(formula, data, K, priors = NULL, ...) {
  # 1. Parse formula
  spec <- parse_bmfrm_formula(formula, data, K)
  
  # 2. Prepare data
  data_info <- prepare_data_bmfrm(spec, data)
  
  # 3. Generate Stan code
  stan_code <- build_stan_code(spec, priors)
  
  # 4. Compile model
  model <- compile_stan_model(stan_code, cache_dir)
  
  # 5. Fit using your existing function (see below)
  fit <- fit_cmdstan_cached(
    model = model,
    data = data_info$stan_data,
    file = cache_file,
    ...
  )
  
  # 6. Return bmfrm_fit object
  structure(
    list(
      fit = fit,
      spec = spec,
      data_info = data_info,
      stan_code = stan_code
    ),
    class = "bmfrm_fit"
  )
}

```

⸻

4.5 Existing Working Components (v0.3)

As of v0.3, we have robust working code for critical components:

### **4.5.1 fit_cmdstan_cached() - Complete caching solution**

Handles cmdstanr's CSV storage properly with:
- Hash-based validation of model + data changes
- Persistent CSV directory structure (`cache_dir/model_name/csv/`)
- Smart refit logic: `"never"`, `"always"`, `"on_change"`
- Model metadata tracking and contamination prevention
- Automatic migration from old cache formats

### **4.5.2 extract_residuals() - Comprehensive residual analysis**

Extracts observation-level information with:
- Posterior summaries for `mu[n]` and `sigma2[n]` from generated quantities
- Standardized residuals: `z = (x - mu_hat) / sqrt(sigma2_hat)`
- Optional posterior draws stored as list-columns for Bayesian uncertainty
- Automatic facet detection and attachment
- Model comparison support

### **4.5.3 facet_infit_outfit() - FACETS-style fit statistics**

Computes fit statistics for any facet:
- Infit and outfit mean squares per facet level
- Supports any facet column name
- Bayesian credible intervals via `add_fit_cri_bayes()`
- Z-standardized fit statistics via `add_zstd()`

### **4.5.4 summarise_bias_ir() - Interaction analysis**

Analyzes rater × item bias effects:
- Posterior summaries for bias matrices
- Probability-based flagging of significant bias
- Automatic label mapping from integer indices
- Credible intervals and effect directions

### **4.5.5 run_ppc() - Posterior predictive checks**

Generates model validation plots:
- PPC bar charts comparing observed vs replicated data
- Distribution checks for mean and variance
- Flexible draw sampling for computational efficiency

### **4.5.6 Revised Package Architecture**

Given these working components, the new architecture becomes:

```
BayesMFRM/
├── R/
│   ├── bmfrm.R              # Main function (NEW)
│   ├── formula_parsing.R   # Parse formula to spec (NEW)
│   ├── data_prep.R         # Facet to integer conversion (NEW)
│   ├── stan_generation.R   # Build Stan code from spec (NEW)
│   ├── prior.R             # Prior specification (NEW)
│   ├── caching.R           # fit_cmdstan_cached() (EXISTING ✓)
│   ├── fit_stats.R         # extract_residuals, facet_infit_outfit (EXISTING ✓)
│   ├── bias_analysis.R     # summarise_bias_ir (EXISTING ✓)
│   ├── ppc.R               # run_ppc (EXISTING ✓)
│   └── methods.R           # S3 methods integrating existing functions (NEW)
└── inst/stan/              # Stan template files (NEW)
```

⸻

5. Internal architecture

For a beginner, it’s useful to see the pipeline:

```
bmfrm()
  ├─ parse_bmfrm_formula()  # understand formula
  ├─ prepare_data_bmfrm()   # build stan_data list
  ├─ build_stan_code()     # create Stan model string
  ├─ compile_or_reuse()    # compiled Stan model with cmdstanr
  ├─ sample()              # run MCMC (cmdstanr::sample)
  └─ postprocess_bmfrm()    # wrap & summarize results
```

We’ll define each piece.

⸻

5.1 Parsing the formula

Function: parse_bmfrm_formula(formula, data, K, family)

Input:
	•	formula like score ~ person + item + rater + (rater:item).
	•	data data.frame.
	•	K (optional).
	•	family ("rating_scale" now).

Output (a “spec” object):

```r
list(
  response     = "score",
  facets_main  = c("person", "item", "rater"),
  facets_bias  = c("rater:item"),
  K            = K_inferred,
  family       = "rating_scale",
  formula      = formula
)
```

How:
	1.	Use terms() or lme4-style parsing to get RHS terms.
	2.	Any RHS term without : becomes a main facet.
	3.	Any RHS term with : becomes a bias facet.
	•	We require interactions in parentheses, e.g. (rater:item) just for clarity, but technically rater:item is enough.
	4.	We check that each variable name exists in data.
	5.	Infer K if not given:

```r
if (is.null(K)) K <- max(data[[response]], na.rm = TRUE)
```


Design choice:
	•	We interpret every RHS variable as a facet.
This is different from general regression (where a variable could be continuous). This makes everything conceptually simpler for MFRM.

⸻

5.2 Data preparation (prepare_data_bmfrm)

Function: prepare_data_bmfrm(spec, data)

Goal:
	•	Convert data frame into a stan_data list with integer indices and counts.
	•	Keep mapping between original IDs and integer indices.

Steps:
	1.	Drop missing scores:

```r
df <- data[!is.na(data[[spec$response]]), , drop = FALSE]
```

	2.	For each main facet in spec$facets_main, create integer indices:

```r
facet_names <- spec$facets_main

idx_list <- list()
J_list   <- list()
label_map <- list()

for (f in facet_names) {
  fac_vec <- df[[f]]
  fac_factor <- factor(fac_vec)           # ensures levels correspond to unique IDs
  idx <- as.integer(fac_factor)          # 1...J_f
  idx_list[[f]] <- idx
  J_list[[f]]   <- nlevels(fac_factor)
  label_map[[f]] <- levels(fac_factor)   # store original labels
}
```

	3.	Build stan_data:
```r
stan_data <- list(
  N = nrow(df),
  K = spec$K,
  X = as.integer(df[[spec$response]])
)

for (f in facet_names) {
  stan_data[[f]] <- idx_list[[f]]
  stan_data[[paste0("J_", f)]] <- J_list[[f]]
}

```
	4.	Return both stan_data and meta-information (labels):

```r
list(
  stan_data = stan_data,
  data_clean = df,
  facet_labels = label_map
)
```


Design choice:
	•	Using factor() + as.integer() is the simplest way to map arbitrary IDs (like "R001", "A", etc.) to consecutive integers 1…J.
	•	We store facet_labels to re-attach human-readable IDs when summarizing.

⸻

5.3 Stan template and code generation

We keep one main Stan template with placeholders.

5.3.1 Template skeleton (simplified)

```stan
data {
  int<lower=1> N;
  int<lower=2> K;

  // {{FACET_DATA_DECL}}

  array[N] int<lower=1, upper=K> X;
}

parameters {
  // person abilities
  vector[J_person] theta_raw;

  // {{FACET_PARAM_DECL}}

  // thresholds
  ordered[K-1] tau_raw;

  // {{BIAS_PARAM_DECL}}
}

transformed parameters {
  vector[J_person] theta;
  // {{FACET_TRANSFORM}}

  vector[K-1] tau = tau_raw;
}

model {
  // {{PRIOR_BLOCK}}

  for (n in 1:N) {
    int person_n = person[n];
    // {{INDEX_EXTRACT}}

    real eta = theta[person_n]
               // {{ETA_SUM}};

    vector[K] logits;
    logits[1] = 0;
    for (k in 2:K)
      logits[k] = logits[k-1] + eta - tau[k-1];

    X[n] ~ categorical_logit(logits);
  }
}

generated quantities {
  vector[N] log_lik;
  array[N] int y_rep;
  vector[N] mu;
  vector[N] sigma2;
  vector[N] mu_fair;

  for (n in 1:N) {
    int person_n = person[n];
    // {{INDEX_EXTRACT}}

    real eta_full = theta[person_n]
                    // {{ETA_SUM}};
    real eta_fair = theta[person_n]
                    // {{ETA_SUM_FAIR}};

    vector[K] logits_full;
    vector[K] logits_fair;
    logits_full[1] = 0;
    logits_fair[1] = 0;
    for (k in 2:K) {
      logits_full[k] = logits_full[k-1] + eta_full - tau[k-1];
      logits_fair[k] = logits_fair[k-1] + eta_fair - tau[k-1];
    }

    log_lik[n] = categorical_logit_lpmf(X[n] | logits_full);
    y_rep[n]   = categorical_logit_rng(logits_full);

    real mu_n = 0;
    real mu_f = 0;
    {
      vector[K] p_full = softmax(logits_full);
      vector[K] p_fair = softmax(logits_fair);
      for (k in 1:K) {
        mu_n += k * p_full[k];
        mu_f += k * p_fair[k];
      }
      mu[n]      = mu_n;
      mu_fair[n] = mu_f;

      real var_n = 0;
      for (k in 1:K)
        var_n += square(k - mu_n) * p_full[k];
      sigma2[n] = var_n;
    }
  }
}
```


5.3.2 Filling placeholders
Given spec$facets_main = c("person", "item", "rater"), the code generator would insert:
	•	{{FACET_DATA_DECL}}:

```stan
array[N] int<lower=1, upper=J_person> person;
array[N] int<lower=1, upper=J_item>   item;
array[N] int<lower=1, upper=J_rater>  rater;

int<lower=1> J_person;
int<lower=1> J_item;
int<lower=1> J_rater;


	•	{{FACET_PARAM_DECL}}:

vector[J_item]  item_raw;
vector[J_rater] rater_raw;


	•	{{FACET_TRANSFORM}}:

{
  real mean_theta = mean(theta_raw);
  theta = theta_raw - mean_theta;
}
{
  real mean_item = mean(item_raw);
  item = item_raw - mean_item;
}
{
  real mean_rater = mean(rater_raw);
  rater = rater_raw - mean_rater;
}

(the actual code might declare vector[J_item] item; in transformed parameters too.)

	•	{{BIAS_PARAM_DECL}} (if facets_bias includes "rater:item"):

matrix[J_rater, J_item] bias_raw;


	•	{{INDEX_EXTRACT}}:

int item_n  = item[n];
int rater_n = rater[n];


	•	{{ETA_SUM}}:
For score ~ person + item + rater + (rater:item):

- item[item_n]
- rater[rater_n]
- bias[rater_n, item_n]


	•	{{ETA_SUM_FAIR}}:
For fair scores, remove rater & bias:

- item[item_n]
// (no rater, no bias)


	•	{{PRIOR_BLOCK}}:
Use priors to generate lines like:

theta_raw  ~ normal(0, 2);
item_raw   ~ normal(0, 2);
rater_raw  ~ student_t(3, 0, 1);
to_vector(bias_raw) ~ normal(0, 1);
tau_raw    ~ normal(0, 3);

```

Design choice:
	•	Using string templates with placeholders keeps Stan code readable and maintainable.
	•	The generator is just replacing certain tokens—simple string gsub() or glue-style expansions.

⸻

5.4 Compilation & caching

Function: compile_bmfrm(stan_code, model_name, cache_dir)

Steps:
	1.	Compute hash of stan_code (e.g., using digest::digest()).
	2.	Build file names:

```r
stan_file   <- file.path(cache_dir, paste0(model_name, "_", hash, ".stan"))
model_rds   <- file.path(cache_dir, paste0(model_name, "_", hash, ".rds"))
```

	3.	If stan_file does not exist, write it.
	4.	If compiled model not yet in memory, call:

```r
model <- cmdstanr::cmdstan_model(stan_file)
```

	5.	Return compiled model and the hash.

Design choice:
	•	Hashing the full Stan code ensures we recompile only when the model structure changes (formula, family, priors, K, etc.).
	•	This is inspired by how brms and your fit_cmdstan_cached() work.

For fitting:

```r
fit_cmdstan_cached(
  model,
  data,
  file,
  refit = c("never", "always", "on_change"),
  seed  = NULL,
  chains = 4,
  parallel_chains = chains,
  ...
)
```
Here is the current implementation

```r
fit_cmdstan_cached <- function(
  model,           # compiled cmdstan_model
  data,            # list passed to model$sample()
  file,            # path to .rds cache file
  refit = c("never", "always", "on_change"),
  seed = NULL,
  chains = 4,
  parallel_chains = chains,
  ...              # additional cmdstanr::sample() args
) {
  refit <- match.arg(refit)

  # ---- helpers ----
  compute_hash <- function(model, data) {
    digest::digest(list(code = model$code(), data = data))
  }

  get_model_name_safe <- function(model) {
    md <- try(model$metadata(), silent = TRUE)
    if (inherits(md, "try-error") || is.null(md$model_name)) {
      return(NA_character_)
    } else {
      return(md$model_name)
    }
  }

  # ---- ensure parent folder exists ----
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)

  # ---- persistent CSV directory ----
  csv_dir <- file.path(dirname(file), "csv")
  dir.create(csv_dir, recursive = TRUE, showWarnings = FALSE)

  # ---- current model signature ----
  current_hash       <- compute_hash(model, data)
  current_model_name <- get_model_name_safe(model)

  # -------------------------------------------------------------------
  # Case 1: Explicit refit = "always"
  # -------------------------------------------------------------------
  if (refit == "always") {
    message("Refitting model (refit = 'always')...")
    fit <- model$sample(
      data = data,
      seed = seed,
      output_dir = csv_dir,
      chains = chains,
      parallel_chains = parallel_chains,
      ...
    )

    meta <- list(
      hash       = current_hash,
      model_name = current_model_name
    )

    saveRDS(list(fit = fit, meta = meta), file = file)
    return(fit)
  }

  # -------------------------------------------------------------------
  # Case 2: Cached file exists → validate before reuse
  # -------------------------------------------------------------------
  if (file.exists(file)) {
    obj <- readRDS(file)

    # New-style object with meta
    if (is.list(obj) && !is.null(obj$fit) && !is.null(obj$meta$hash)) {

      cached_hash       <- obj$meta$hash
      cached_model_name <- obj$meta$model_name %||% NA_character_

      # --- safety check: model name mismatch ---
      if (!is.na(current_model_name) &&
          !is.na(cached_model_name) &&
          !identical(current_model_name, cached_model_name)) {

        stop(
          "Cached file '", file, "' was created for model '", cached_model_name,
          "', but current model is '", current_model_name, "'.\n",
          "To avoid contamination, please use a different 'file' path or delete the old cache."
        )
      }

      # --- safety check: hash mismatch ---
      if (!identical(cached_hash, current_hash)) {
        if (refit == "never") {
          stop(
            "Cached file '", file, "' does not match current model+data hash, ",
            "but refit = 'never'.\n",
            "Refusing to reuse potentially contaminated cache. ",
            "Use refit = 'on_change' or 'always', or delete the file."
          )
        }

        # refit == "on_change"
        message("Model or data changed → refitting (refit = 'on_change').")
        fit <- model$sample(
          data = data,
          seed = seed,
          output_dir = csv_dir,
          chains = chains,
          parallel_chains = parallel_chains,
          ...
        )

        meta <- list(
          hash       = current_hash,
          model_name = current_model_name
        )

        saveRDS(list(fit = fit, meta = meta), file = file)
        return(fit)
      }

      # --- hashes match → safe to reuse ---
      message("Loading cached fit (model & data unchanged).")
      return(obj$fit)
    }

    # -----------------------------------------------------------------
    # Old-style RDS: just a CmdStanMCMC object, no meta
    # -----------------------------------------------------------------
    if (refit == "never") {
      message("Loading old-style cached fit (no meta; refit = 'never'). ",
              "Cannot fully verify contamination, use with care.")
      return(obj)
    }

    if (refit == "on_change") {
      message("Migrating old cached fit (no meta) → wrapping & adding hash (no refit).")

      # Ensure CSVs exist
      missing <- !file.exists(obj$output_files())
      if (any(missing)) {
        stop("Old cached fit refers to missing CSV files. You must refit the model.")
      }

      meta <- list(
        hash       = current_hash,
        model_name = current_model_name
      )

      obj_new <- list(fit = obj, meta = meta)
      saveRDS(obj_new, file = file)
      return(obj_new$fit)
    }
  } else {
    message("No cached file found → fitting model.")
  }

  # -------------------------------------------------------------------
  # Case 3: Must refit (no cache, or migrated/invalid)
  # -------------------------------------------------------------------
  message("Refitting model...")
  fit <- model$sample(
    data = data,
    seed = seed,
    output_dir = csv_dir,
    chains = chains,
    parallel_chains = parallel_chains,
    ...
  )

  meta <- list(
    hash       = current_hash,
    model_name = current_model_name
  )

  saveRDS(list(fit = fit, meta = meta), file = file)
  return(fit)
}

```

⸻

5.5 Post-processing

Function: postprocess_bmfrm(fit, spec, data_info, priors, stan_code_hash)

It should build an object, say class "bmfrm_fit":

```r
out <- list(
  fit          = fit,              # cmdstanr fit
  spec         = spec,             # formula, facets, family, K
  data_info    = data_info,        # cleaned data, facet label maps
  priors       = priors,
  stan_code    = stan_code_string,
  stan_hash    = stan_code_hash
)

class(out) <- "bmfrm_fit"
out
```

We then define methods:
	•	print.bmfrm_fit()
	•	summary.bmfrm_fit()
	•	stancode.bmfrm_fit() – returns Stan code string.
	•	standata.bmfrm_fit() – returns stan_data list.
	•	facet_summary(fit, facet = "rater")
	•	fair_scores(fit, newdata = NULL)
	•	extract_residuals(fit, ...)

Design choice:
	•	For a beginner, having standard S3 methods means they can call summary(fit) and plot(fit) like in base R / brms.
	•	Exposing stancode() and standata() makes advanced usage and debugging easier.

⸻

6. Utility layers

6.1 facet_summary()

Goal: FACETS-style facet tables.

```r
facet_summary <- function(object, facet, probs = c(0.025, 0.5, 0.975)) {
  # 1. Extract posterior draws for that facet’s parameters.
  # 2. Summarize: mean, sd, quantiles.
  # 3. Attach original labels from `object$data_info$facet_labels[[facet]]`.
}
```

This gives a data frame like:

facet	label	mean	sd	q2.5	q50	q97.5
rater	R01	…	..	…	…	…

6.2 fair_scores()

Uses mu_fair[n] from generated quantities:

```r
fair_scores <- function(object, newdata = NULL, summary = TRUE) {
  # For now, just return a data.frame with:
  # person, item, rater, score, mu_fair, etc.
}
```

6.3 extract_residuals()

You already have a design: use mu, sigma2 to compute standardized residuals and then add facet columns.
Here is the current idea. 

```r
#' Extract residual info from a cmdstanr fit (with optional posterior draws)
#'
#' @param fit         cmdstanr fit object (must have mu[n], sigma2[n] in GQ)
#' @param stan_data   list with X and one or more facet columns
#' @param facet_cols  which columns from stan_data are facets
#' @param save_draws  if TRUE, save posterior draws of mu and sigma2 per obs
#' @param model_name  optional model label
#'
#' @return tibble with:
#'   n, x, mu_hat, sigma2_hat, resid, z, weight, model,
#'   facet columns,
#'   and (if save_draws=TRUE):
#'       mu_draws[n]      list-column of length-D vectors
#'       sigma2_draws[n]  list-column of length-D vectors
extract_residuals <- function(
    fit,
    stan_data,
    facet_cols = NULL,
    save_draws = TRUE,
    model_name = deparse(substitute(fit))
) {
  # ───────────────────────────────────────────────────────────────
  # 1. Posterior summaries for mu and sigma2
  # ───────────────────────────────────────────────────────────────
  mu_summ  <- fit$summary("mu")
  sig_summ <- fit$summary("sigma2")
  
  get_index <- function(x) as.integer(stringr::str_extract(x, "(?<=\\[)\\d+(?=\\])"))
  
  mu_df <- mu_summ %>%
    transmute(
      n      = get_index(variable),
      mu_hat = mean
    )
  
  sig_df <- sig_summ %>%
    transmute(
      n          = get_index(variable),
      sigma2_hat = mean
    )
  
  # ───────────────────────────────────────────────────────────────
  # 2. Determine facets
  # ───────────────────────────────────────────────────────────────
  N <- length(stan_data$X)
  
  if (is.null(facet_cols)) {
    lens <- vapply(stan_data, length, integer(1))
    facet_cols <- names(stan_data)[lens == N & names(stan_data) != "X"]
  }
  
  # ───────────────────────────────────────────────────────────────
  # 3. Build observation-level table
  # ───────────────────────────────────────────────────────────────
  obs_df <- tibble(
    n = seq_len(N),
    x = stan_data$X
  )
  
  for (nm in facet_cols) {
    obs_df[[nm]] <- stan_data[[nm]]
  }
  
  # ───────────────────────────────────────────────────────────────
  # 4. Add point summaries
  # ───────────────────────────────────────────────────────────────
  out <- obs_df %>%
    left_join(mu_df,  by = "n") %>%
    left_join(sig_df, by = "n") %>%
    mutate(
      resid  = x - mu_hat,
      weight = sigma2_hat,
      z      = resid / sqrt(sigma2_hat),
      model  = model_name
    )
  
  # ───────────────────────────────────────────────────────────────
  # 5. Optional: Attach posterior draws (compact list-column)
  # ───────────────────────────────────────────────────────────────
  if (save_draws) {
    # raw draws (iterations x N)
    mu_mat  <- fit$draws("mu",     format = "matrix")
    sig_mat <- fit$draws("sigma2", format = "matrix")
    
    # sort columns by [n]
    mu_mat  <- mu_mat[,  order(get_index(colnames(mu_mat))), drop=FALSE]
    sig_mat <- sig_mat[, order(get_index(colnames(sig_mat))), drop=FALSE]
    
    # attach list columns
    out$mu_draws    <- split(as.data.frame(t(mu_mat)),    seq_len(N))
    out$sigma2_draws <- split(as.data.frame(t(sig_mat)), seq_len(N))
    
    # each list element will be a numeric vector of draws
    out <- out %>%
      mutate(
        mu_draws    = purrr::map(mu_draws,    unlist),
        sigma2_draws = purrr::map(sigma2_draws, unlist)
      )
  }
  
  # attach metadata: facet columns
  attr(out, "facet_cols") <- facet_cols
  
  out
}
```

6.4 fit statistics

```r
# In methods.R
fit_statistics <- function(object, facet) {
  res_df <- extract_residuals(
    fit = object$fit,
    stan_data = object$data_info$stan_data,
    facet_cols = object$spec$facets_main
  )
  
  facet_infit_outfit(res_df, facet)
}
```

```r
#' Compute infit and outfit mean square for any facet
#'
#' @param res_df tibble from extract_residuals()
#' @param facet  facet column to group by; can be a bare name or a string
#'               e.g. RaterID, "RaterID", ExamineeID, "TaskID", "CriterionID", ...
#'
#' @return tibble with facet_value, model, n_obs, infit_msq, outfit_msq
facet_infit_outfit <- function(res_df, facet) {
  # allow both strings and bare names
  facet_sym <- if (is.character(facet)) sym(facet) else enquo(facet)
  
  res_df %>%
    group_by(!!facet_sym, model) %>%
    summarise(
      n_obs      = n(),
      outfit_msq = mean(z^2, na.rm = TRUE),
      infit_msq  = sum(weight * z^2, na.rm = TRUE) / sum(weight, na.rm = TRUE),
      .groups    = "drop"
    ) %>%
    rename(facet_value = !!facet_sym)
}

```

```r
#' Bayesian CrIs for infit & outfit using draws stored in residuals
add_fit_cri_bayes <- function(fit_stats,
                              res_df,
                              facet,
                              prob = 0.95) {
  facet_name <- if (is.character(facet)) facet else rlang::as_name(rlang::enquo(facet))
  
  prob_lo <- (1 - prob) / 2
  prob_hi <- 1 - prob_lo
  
  facet_levels <- sort(unique(res_df[[facet_name]]))
  
  cri_df <- purrr::map_dfr(
    facet_levels,
    function(f_val) {
      idx <- which(res_df[[facet_name]] == f_val)
      
      mu_list  <- res_df$mu_draws[idx]      # list of numeric vectors
      sig_list <- res_df$sigma2_draws[idx]
      
      # matrix: draws x obs_in_facet
      mu_mat  <- do.call(cbind, mu_list)
      sig_mat <- do.call(cbind, sig_list)
      
      D <- nrow(mu_mat)
      n_obs <- length(idx)
      
      X_sub <- matrix(res_df$x[idx], nrow=D, ncol=n_obs, byrow=TRUE)
      
      z2 <- (X_sub - mu_mat)^2 / sig_mat
      w  <- sig_mat
      
      outfit <- rowMeans(z2)
      infit  <- rowSums(w * z2) / rowSums(w)
      
      tibble(
        facet_value = f_val,
        outfit_lo = quantile(outfit, prob_lo),
        outfit_hi = quantile(outfit, prob_hi),
        infit_lo  = quantile(infit,  prob_lo),
        infit_hi  = quantile(infit,  prob_hi)
      )
    }
  )
  
  fit_stats %>%
    left_join(cri_df, by = "facet_value")
}

```

```r
facet_infit_outfit_all <- function(res_df) {
  facet_cols <- attr(res_df, "facet_cols")
  
  map_dfr(
    facet_cols,
    ~ facet_infit_outfit(res_df, .x) %>%
      mutate(facet = .x),
    .id = NULL
  ) %>%
    relocate(facet, facet_value)
}
```

⸻

7. Design decisions and rationale (for beginners)
	1.	Why formula syntax?
	•	R users are familiar with y ~ x1 + x2. We reuse that mental model, but reinterpret x’s as facets, not predictors.
	2.	Why integer indices?
	•	Stan likes integer indices for categorical things. Converting "R01" → 1, "R02" → 2 is standard practice in hierarchical models.
	3.	Why sum-to-zero constraints?
	•	Rasch/MFRM parameters are only identified up to a linear shift (we can add a constant to all person abilities and subtract it from item difficulties).
	•	Enforcing that a facet’s parameters sum to zero anchors that facet’s scale and avoids non-identifiability.
  •	Row-wise centering for theta to retain global meaning beyond specific raters.
	4.	Why only rating-scale Rasch first?
	•	Implementing everything (partial credit, GRM, drift, multidimensional) at once is complicated.
	•	Starting with one clean, well-understood model lets us:
	•	test infrastructure,
	•	debug Stan code generation,
	•	design a stable R API.
	5.	Why a small prior vocabulary?
	•	brms supports many classes, but that’s overwhelming early on.
	•	MFRM use cases mostly need priors on:
	•	person ability,
	•	item difficulty,
	•	rater severity,
	•	thresholds,
	•	bias terms.
	•	So a small, Rasch-specific dictionary is enough.
	6.	Why caching?
	•	Compiling Stan models is expensive.
	•	In research, you often re-run the same model with minor data changes.
	•	Caching + hashing the Stan code saves time and frustration.


⸻

8. Implementation roadmap (v0.3 - leveraging existing code)

You (or a beginner dev) can follow this order:
	1.	Package skeleton
	•	usethis::create_package("BayesMFRM")
	•	Add R/, inst/stan/, DESCRIPTION with imports: cmdstanr, posterior, dplyr, tibble, rlang, glue, digest.
	2.	Core helpers
	•	Implement prior() and basic validation.
	•	Implement parse_bmfrm_formula() for simple cases.
	3.	Data prep
	•	Implement prepare_data_bmfrm():
	•	stan_data, facet_labels.
	4.	Stan template
	•	Save base template as a string (or .stan file with tokens).
	•	Implement build_stan_code(spec, priors) that fills placeholders for a 3-facet model.
	5.	Compile & fit
	•	Implement compile_bmfrm().
	•	Reuse fit_cmdstan_cached().
	•	Implement bmfrm() itself (wire everything together).
	6.	Post-processing
	•	Define bmfrm_fit class.
	•	Implement summary.bmfrm_fit(), stancode.bmfrm_fit(), standata.bmfrm_fit().
	•	Implement facet_summary(), fair_scores().
	7.	Testing
	•	Simulate small Rasch data sets, fit the model, and verify:
	•	Parameter recovery.
	•	Reasonable diagnostics (Rhat, ESS).
	•	Reasonable fair-score behavior.
	8.	Documentation
	•	Write vignettes:
	•	“Basic 3-facet MFRM”
	•	“Adding a rater×item bias facet”
	•	“Inspecting Stan code and priors”

9. Example usecase


