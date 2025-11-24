# BayesMFRM – Design Specification (v0.4)

## 1. Big picture

**Goal:**  
BayesMFRM is an R package for Bayesian many-facet Rasch models using Stan (via cmdstanr).

It should:

- Feel as easy as **brms**, but specialized for MFRM.
- Hide most Stan details, but still allow advanced users to inspect Stan code.
- Produce **FACETS-style** outputs: rater severity, item difficulty, fair averages, infit/outfit, etc.

We start simple (v0.1):

- Rating-scale Rasch model (shared thresholds).
- Main facet effects (person, item, rater, etc.).
- Optional 2-way bias facets (e.g., rater × item).
- A clean, simple prior system.

Later, we can extend to:

- Partial credit / GRM
- Centrality / drift
- Multidimensional abilities
- More complex bias structures (rater × task, etc.)

---

## 2. User-facing workflow

### 2.1 Data Requirements

Before calling `bmfrm()`, ensure your data meets these requirements:

- **Score column**: Integer values from 1 to K (not 0-indexed)
- **Facet columns**: Any values that can be converted to factors (strings, numbers, etc.)
- **Missing values**: 
  - NA scores will be dropped with a message
  - NA in facets will cause an error (clean your data first)
  
**Example data formats:**
```r
# Good data (ready for bmfrm)
ratings$score   # 1, 2, 3, 4, 5, 6
ratings$person  # "P001", "P002", ...
ratings$rater   # "R01", "R02", ...

# Also acceptable (will be auto-converted)
ratings$person  # 1, 2, 3, ... (numeric IDs)
ratings$rater   # 101, 102, 103, ... (numeric IDs)

# Will cause an error - fix before calling bmfrm()
ratings$score   # 0, 1, 2, 3, 4, 5  (0-indexed - add 1 first)
ratings$rater   # "R01", NA, "R03"  (NA in facets not allowed)
```

### 2.2 Core Workflow

The ideal user workflow follows four simple steps:

1. **Load data** → 2. **Fit with `bmfrm()`** → 3. **Check model with `pp_check()`** → 4. **Analyze residuals**

```r
library(BayesMFRM)

# Step 1: Load data (ensure it meets requirements above)
ratings <- readr::read_csv("ratings.csv")

# Step 2: Fit the model (minimal specification)
fit <- bmfrm(
  score ~ person + item + rater + rater:item,
  data = ratings,
  K = 6  # Optional if max(score) = 6
)

# Step 3: Model checking
summary(fit)
pp_check(fit)                    # Default posterior predictive check
pp_check(fit, type = "stat")     # Check distributional properties

# Step 4: Residual analysis and facet summaries
resid_tbl <- residuals(fit)      # Extract residuals
rater_fit <- facet_fit(fit, "rater")   # Rater fit statistics
item_fit  <- facet_fit(fit, "item")    # Item fit statistics
bias_tab  <- summarise_bias(fit)       # Bias analysis
crit_cor  <- residual_cor_criteria(resid_tbl)  # Cross-criteria correlations

# Additional summaries
facet_summary(fit, "rater")      # Detailed rater severity table
fair_scores(fit)                 # Fair scores (rater-bias adjusted)
```

### 2.3 Advanced Options (Optional)

For more control, you can specify priors, MCMC settings, and caching:

```r
# Custom priors
priors <- c(
  prior("normal(0, 2)",       class = "theta"),        # person
  prior("normal(0, 2)",       class = "item"),         # item difficulty
  prior("student_t(3, 0, 1)", class = "rater"),        # rater severity
  prior("normal(0, 1)",       class = "bias", facet = "rater:item"),
  prior("normal(0, 3)",       class = "tau")           # thresholds
)

fit <- bmfrm(
  score ~ person + item + rater + rater:item,
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

# For debugging/inspection (advanced users)
stancode(fit)   # Inspect generated Stan code
standata(fit)   # Inspect Stan data structure
```

Key design idea:

**The user thinks in terms of facets & Rasch. The package handles everything else:**
- Data preparation → Stan code generation → MCMC sampling → Rasch-style diagnostics
- Simple S3 methods hide complexity while preserving access to underlying details

⸻

## 3. S3 Methods Reference

BayesMFRM provides a clean S3 interface that hides complexity while preserving access to details.

### 3.1 Primary User Interface

These methods support the core workflow and should be used by typical users:

```r
# Model fitting
fit <- bmfrm(formula, data, ...)         # Main fitting function

# Model summaries  
print(fit)                               # Basic model info
summary(fit)                             # Comprehensive summary
facet_summary(fit, "rater")              # Facet-level summaries
fair_scores(fit)                         # Fair score extraction

# Model checking
pp_check(fit)                            # Posterior predictive checks
pp_check(fit, type = "stat")            # Distributional checks

# Residual analysis  
residuals(fit)                           # Extract residuals (auto-handles stan_data)
facet_fit(fit, "rater")                  # Facet fit statistics (infit/outfit)
summarise_bias(fit)                      # Bias analysis
```

### 3.2 Utility Functions

These functions work with outputs from the primary interface:

```r
# Data processing utilities
residual_cor_criteria(resid_tbl)         # Cross-criteria correlations

# Convenience functions  
facet_fit_all(fit)                       # Fit statistics for all facets
```

### 3.3 Advanced/Developer Interface

These methods are for debugging, inspection, or advanced customization:

```r
# Stan inspection
stancode(fit)                            # View generated Stan code
standata(fit)                            # View Stan data structure

# Low-level data processing (normally called internally)
extract_residuals(fit, stan_data, ...)   # Direct residual extraction
run_ppc(fit, stan_data, ...)            # Direct PPC computation
prepare_data_bmfrm(spec, data)           # Data preparation
parse_bmfrm_formula(formula, data, ...)  # Formula parsing
```

### 3.4 Method Signatures

**Key principle**: User-facing methods should not require `stan_data` or other internal objects.

```r
# User-facing signatures (simple)
residuals.bmfrm_fit <- function(object, save_draws = FALSE, ...)
pp_check.bmfrm_fit <- function(object, type = "bars", ...)  
facet_fit.bmfrm_fit <- function(object, facet, ...)

# Internal signatures (complex - for advanced users only)
extract_residuals <- function(fit, stan_data, facet_cols, save_draws = TRUE, ...)
run_ppc <- function(fit, stan_data, type = "bars", ...)
```

This separation ensures beginners can use simple commands while experts can access full control when needed.

⸻

## 4. Core concepts

4.1 What is a "facet"?

In MFRM, facets are "dimensions" of the assessment design:
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
ratings$person   # examinee ID
ratings$item     # item or criterion ID
ratings$rater    # rater ID
ratings$score    # rating 1…K
```


In Stan, we convert these into integer indices 1…J_f (e.g., 1…J_person).

4.2 Main facets vs bias facets
	•	Main facets: person, item, rater, etc.
	•	Bias facets (interactions): combinations like rater:item that represent specific patterns (rater–item bias).

In the model formula:

```
score ~ person + item + rater + rater:item
```

	•	person, item, rater are main facets.
	•	(rater:item) is a bias facet.

In Stan:
	•	main facets → vectors of latent parameters.
	•	bias facets → matrices or arrays of parameters.

⸻

5. Statistical model (v0.1)

We start with a rating-scale Rasch bmfrm.

For rating (n):
	•	Person: (p_n)
	•	Item: (i_n)
	•	Rater: (r_n)
	•	Score: (X_n \in {1, \dots, K})

Linear predictor

[
\eta_n
= \theta_{p_n}
	•	\beta_{i_n}
	•	\rho_{r_n}

	•	\text{(optional bias terms)}
]

Rating-scale structure

[
\log\frac{P(X_n \ge k)}{P(X_n < k)} = \eta_n - \tau_{k-1}, \quad k=2,\dots,K
]
	•	(\theta): person ability
	•	(\beta): item difficulty
	•	(\rho): rater severity
	•	(\tau_k): rating thresholds (shared across items/raters in v0.1)

Assumptions:
	•	Equal discrimination (Rasch).
	•	Optional additional facets enter additively (e.g., task difficulty).
	•	Fair scores are derived from a version of (\eta_n) that removes selected nuisance facets (e.g., rater, bias).

⸻

6. Top-level API design

6.1 bmfrm() function

```r
bmfrm <- function(
  formula,
  data,
  K,
  priors      = NULL,
  family      = c("rating_scale", "partial_credit"),
  model_name  = NULL,
  cache_dir   = "stan_cache",
  save_csvs   = TRUE,         # Always save CSVs for posterior access
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
R formula specifying the model, e.g.:

score ~ person + item + rater + (rater:item)

	•	LHS: score column (integer 1…K).
	•	RHS: facet columns; interaction terms (optionally in parentheses) for bias facets.

	•	data
Data frame with all variables named in formula.
	•	K
Number of score categories. If NULL, infer from max(score).
	•	priors
Vector of prior() objects (see below). If NULL, use safe defaults.
	•	family
"rating_scale" (v0.1 only). Later: "partial_credit" etc.
	•	model_name
Optional string used to label compilation and caching.
If NULL, automatically generate something like "bmfrm_<hash>".
	•	cache_dir
Directory to store compiled Stan models and cached fits.
	•	refit
	•	"on_change": refit only if Stan code or data structure changed.
	•	"never": never refit (just reuse existing cached fit if present).
	•	"always": always refit.
	•	iter, warmup, chains, cores, seed
Standard cmdstanr::sample() arguments.

5.2 prior() object

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

Allowed classes (v0.1)
	•	"theta" – person abilities
	•	"item" – item difficulties
	•	"rater" – rater severities
	•	"facet_<name>" – for extra main facets (e.g. "facet_task")
	•	"bias" – interaction (bias) terms; must specify facet = "rater:item" etc.
	•	"tau" – rating thresholds

Design choice:
	•	For a beginner, it’s easier to say “Normal(0, 2) prior on person ability” than to reason about abstract b / sd / cor classes.
	•	We define a limited vocabulary that maps directly to psychometric concepts.

Defaults when priors = NULL (v0.1 idea):
	•	theta ~ normal(0, 2)
	•	item  ~ normal(0, 1)
	•	rater ~ normal(0, 1)
	•	bias  ~ normal(0, 0.5)
	•	tau   ~ normal(0, 3) on ordered thresholds (or on spacings)

⸻

5.3 bmfrm() pipeline (internal)

Conceptual implementation:

```r
bmfrm <- function(formula, data, K, priors = NULL, family = "rating_scale", ...) {
  # 1. Parse formula → spec
  spec <- parse_bmfrm_formula(formula, data, K, family = family)
  
  # 2. Prepare data → stan_data + facet labels
  data_info <- prepare_data_bmfrm(spec, data)
  
  # 3. Generate Stan code string
  stan_code <- build_stan_code(spec, priors)
  
  # 4. Compile Stan model (with caching based on code hash)
  comp <- compile_bmfrm(stan_code, model_name = spec$model_name, cache_dir = spec$cache_dir)
  model <- comp$model
  cache_file <- comp$cache_file
  
  # 5. Fit using existing fit_cmdstan_cached()
  fit <- fit_cmdstan_cached(
    model = model,
    data  = data_info$stan_data,
    file  = cache_file,
    ...
  )
  
  # 6. Wrap as bmfrm_fit object
  out <- postprocess_bmfrm(
    fit          = fit,
    spec         = spec,
    data_info    = data_info,
    priors       = priors,
    stan_code    = stan_code,
    stan_code_hash = comp$hash
  )
  
  out
}

```

⸻

6. Existing working components (v0.3/v0.4)

As of now, we have robust working code for critical components (in helper.R):

6.1 fit_cmdstan_cached() – complete caching solution

Handles cmdstanr’s CSV storage properly with:
	•	Hash-based validation of model + data changes
	•	Persistent CSV directory structure (cache_dir/model_name/csv/)
	•	Smart refit logic: "never", "always", "on_change"
	•	Model metadata tracking and contamination prevention
	•	Automatic migration from old cache formats

This function is reused by bmfrm() and should not be newly reimplemented, but should be reused.

⸻

6.2 extract_residuals() – comprehensive residual analysis

Extracts observation-level information with:
	•	Posterior summaries for mu[n] and sigma2[n] from generated quantities
	•	Standardized residuals: z = (x - mu_hat) / sqrt(sigma2_hat)
	•	Optional posterior draws stored as list-columns for Bayesian uncertainty
	•	Automatic facet detection and attachment
	•	Support for model comparison (via model label)

bmfrm will call extract_residuals() as-is, with appropriate stan_data and facet columns.

⸻

6.3 facet_infit_outfit() / facet_infit_outfit_all() – FACETS-style fit statistics

Compute fit statistics for any facet:
	•	Infit and outfit mean squares per facet level
	•	Works with any facet column name
	•	Bayesian credible intervals via add_fit_cri_bayes()
	•	Z-standardized fit statistics via add_zstd()

These functions will be used in S3 methods like:
	•	fit_statistics.bmfrm_fit() (future)
	•	or helpers such as facet_fit(fit, facet = "rater")

⸻

6.4 summarise_bias() 

Analyzes bias effects:
	•	Posterior summaries for bias matrices
	•	Probability-based flagging of “significant bias”
	•	Automatic label mapping from integer indices
	•	Credible intervals and effect directions

In v0.4, the user-facing function is:

summarise_bias(fit, facet = "rater:item")


⸻

6.5 run_ppc() – posterior predictive checks

Generates model validation summaries/plots:
	•	PPC bar charts comparing observed vs replicated data
	•	Distribution checks for mean and variance
	•	Flexible draw sampling for computational efficiency

bmfrm will expose PPC via a helper (e.g., pp_check.bmfrm_fit() or ppc_bmfrm()).

⸻

6.6 Revised package architecture

Given these working components, the package layout is:

```
BayesMFRM/
├── R/
│   ├── bmfrm.R              # Main function (NEW)
│   ├── formula_parsing.R    # parse_bmfrm_formula() (NEW)
│   ├── data_prep.R          # prepare_data_bmfrm() (NEW)
│   ├── stan_generation.R    # build_stan_code() (NEW)
│   ├── prior.R              # prior() + validation (NEW)
│   ├── compile.R            # compile_bmfrm() (NEW)
│   ├── caching.R            # fit_cmdstan_cached() (EXISTING ✓)
│   ├── fit_stats.R          # extract_residuals, facet_infit_outfit (EXISTING ✓)
│   ├── bias_analysis.R      # summarise_bias_ir, summarise_bias() (EXISTING + NEW)
│   ├── ppc.R                # run_ppc (EXISTING ✓)
│   ├── residual_cor.R       # residual_cor_criteria() (NEW)
│   └── methods.R            # S3 methods for bmfrm_fit (NEW)
└── inst/stan/
    ├── bmfrm_template.stan  # Base template with placeholders (NEW)
    └── ...
```

⸻

1. Internal architecture

For a beginner, it’s useful to see the pipeline:

```
bmfrm()
  ├─ parse_bmfrm_formula()   # understand formula → spec
  ├─ prepare_data_bmfrm()    # build stan_data list + labels
  ├─ build_stan_code()       # create Stan model string from spec + priors
  ├─ compile_bmfrm()         # compiled Stan model with cmdstanr (and hash)
  ├─ fit_cmdstan_cached()    # run MCMC with caching
  └─ postprocess_bmfrm()     # wrap & summarize results as bmfrm_fit
```

We define each piece below.

⸻

7.1 Parsing the formula – parse_bmfrm_formula()

Function:

parse_bmfrm_formula <- function(formula, data, K = NULL, family = "rating_scale") { ... }

Input:
	•	formula like score ~ person + item + rater + (rater:item)
	•	data data.frame
	•	K (optional)
	•	family (currently "rating_scale")

Output (a “spec” object):

```r
list(
  response     = "score",
  facets_main  = c("person", "item", "rater"),
  facets_bias  = c("rater:item"),
  K            = K_inferred,
  family       = "rating_scale",
  formula      = formula,
  model_name   = model_name_or_default,
  cache_dir    = cache_dir
)
```

Steps:
	1.	Use terms() or similar to parse the formula.
	2.	Any RHS term without : → main facet.
	3.	Any RHS term with : → bias facet.
	•	Parentheses are allowed but not required:
	•	(rater:item) or rater:item are both accepted.
	4.	Check that each variable exists in data.
	5.	Infer K if not given:

```r
if (is.null(K)) {
  resp <- data[[response]]
  if (!all(resp == floor(resp))) stop("Response must be integer for rating-scale MFRM.")
  K <- max(resp, na.rm = TRUE)
}
```

Design choice:
	•	Every RHS variable is interpreted as a facet (categorical).
	•	This simplifies the mental model for MFRM users.

⸻

7.2 Data preparation – prepare_data_bmfrm()

Function:

```r
prepare_data_bmfrm <- function(spec, data) { ... }
```

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

idx_list   <- list()
J_list     <- list()
label_map  <- list()

for (f in facet_names) {
  fac_vec    <- df[[f]]
  fac_factor <- factor(fac_vec)
  idx        <- as.integer(fac_factor)          # 1...J_f
  idx_list[[f]]   <- idx
  J_list[[f]]     <- nlevels(fac_factor)
  label_map[[f]]  <- levels(fac_factor)         # store original labels
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
  stan_data[[f]]               <- idx_list[[f]]
  stan_data[[paste0("J_", f)]] <- J_list[[f]]
}
```

	4.	Return both stan_data and meta-information:

```r
list(
  stan_data    = stan_data,
  data_clean   = df,
  facet_labels = label_map
)
```

⸻

7.3 Stan template and code generation – build_stan_code()

We keep one main Stan template with placeholders.

7.3.1 Template skeleton (simplified)

```
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


7.3.2 Filling placeholders
Given spec$facets_main = c("person", "item", "rater"), the code generator inserts:
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
```

	•	{{FACET_TRANSFORM}}: sum-to-zero constraints per facet:

```
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
```

	•	{{BIAS_PARAM_DECL}} (if facets_bias includes "rater:item"):

```
matrix[J_rater, J_item] bias_raw;
```
plus a centered bias matrix if desired.

	•	{{INDEX_EXTRACT}}:
```
int item_n  = item[n];
int rater_n = rater[n];
```

	•	{{ETA_SUM}}
For score ~ person + item + rater + (rater:item):

```
+ (- item[item_n])
+ (- rater[rater_n])
+ bias[rater_n, item_n]
```

	•	{{ETA_SUM_FAIR}}
For fair scores, drop rater and bias, keeping item:
```
+ (- item[item_n])
// no rater, no bias
```

	•	{{PRIOR_BLOCK}}:
uses priors to generate lines like:
```
theta_raw  ~ normal(0, 2);
item_raw   ~ normal(0, 2);
rater_raw  ~ student_t(3, 0, 1);
to_vector(bias_raw) ~ normal(0, 1);
tau_raw    ~ normal(0, 3);
```


Design choice:
	•	Using string templates with placeholders keeps Stan code readable and maintainable.
	•	The generator simply replaces tokens using glue or gsub.

⸻

7.4 Compilation & caching – compile_bmfrm()

Function:

```r
compile_bmfrm <- function(stan_code, model_name, cache_dir) { ... }
```

Steps:
	1.	Compute hash of stan_code (e.g., digest::digest(stan_code)).
	2.	Build file names:

```r
stan_file <- file.path(cache_dir, paste0(model_name, "_", hash, ".stan"))
cache_file <- file.path(cache_dir, paste0(model_name, "_", hash, ".rds"))
```

	3.	If stan_file does not exist, write it.
	4.	Compile:

```r
model <- cmdstanr::cmdstan_model(stan_file)
```

	5.	Return list:

```r
list(
  model      = model,
  hash       = hash,
  stan_file  = stan_file,
  cache_file = cache_file
)

```

Design choice:
	•	Hashing the full Stan code ensures we recompile only when the model structure changes.
	•	bmfrm() uses fit_cmdstan_cached() for sampling + reuse of fits.

⸻

7.5 Post-processing – postprocess_bmfrm()

Function:

```r
postprocess_bmfrm <- function(fit, spec, data_info, priors, stan_code, stan_code_hash) { ... }
```

Build a "bmfrm_fit" S3 object:

```r
out <- list(
  fit           = fit,              # cmdstanr fit
  spec          = spec,             # formula, facets, family, K
  data_info     = data_info,        # cleaned data, facet label maps, stan_data
  priors        = priors,
  stan_code     = stan_code,
  stan_hash     = stan_code_hash
)

class(out) <- "bmfrm_fit"
out
```

S3 methods:
	•	print.bmfrm_fit()
	•	summary.bmfrm_fit()
	•	stancode.bmfrm_fit() – returns Stan code string.
	•	standata.bmfrm_fit() – returns stan_data list.
	•	facet_summary() – facet-level summary.
	•	fair_scores() – fair averages / scores.
	•	extract_residuals() – wraps existing helper using fit + data_info$stan_data.

⸻

8. Utility layers

8.1 facet_summary()

Goal: FACETS-style facet tables.

```r
facet_summary <- function(object, facet, probs = c(0.025, 0.5, 0.975)) {
  # 1. Extract posterior draws for that facet's parameters from object$fit.
  # 2. Summarise: mean, sd, quantiles.
  # 3. Attach original labels from object$data_info$facet_labels[[facet]].
}
```

Output:

facet	label	mean	sd	q2.5	q50	q97.5
rater	R01	…	..	…	…	…


⸻

8.2 fair_scores()

Uses mu_fair[n] from generated quantities.

Function:

```r
fair_scores <- function(object, summary = TRUE) {
  # Extract mu_fair[n] draws or summaries from object$fit
  # Join with person/item/rater IDs from object$data_info
}
```

v0.1 proposal:
	•	If summary = FALSE: return an observation-level data.frame:
	•	columns: person, item, rater, score, mu, mu_fair, etc.
	•	If summary = TRUE: additionally return/attach:
	•	person-level fair scores (e.g., mean mu_fair by person).

Exact format can evolve, but spec emphasizes:
	•	mu_fair removes rater/bias effects, preserving examinee and item structure.

⸻

8.3 extract_residuals() (existing) – contract

We reuse the existing implementation. The contract for the Stan model:
	•	generated quantities must contain:
	•	mu[N]
	•	sigma2[N]
	•	stan_data must contain:
	•	X (observed scores)
	•	facet integer indices (person, item, rater, etc.)

Then:

```r
resid_tbl <- extract_residuals(
  fit        = object$fit,
  stan_data  = object$data_info$stan_data,
  facet_cols = object$spec$facets_main,
  save_draws = TRUE
)
```

returns a tibble with:
	•	n, x, mu_hat, sigma2_hat, resid, z, weight, model
	•	facet columns (person, item, rater, …)
	•	optional mu_draws, sigma2_draws list-columns.

attr(resid_tbl, "facet_cols") lists facet columns for facet_infit_outfit_all().

⸻

8.4 Fit statistics – facet_infit_outfit() / facet_infit_outfit_all()

We reuse existing helpers.

Example in methods.R:

```r
fit_statistics <- function(object, facet) {
  res_df <- extract_residuals(
    fit        = object$fit,
    stan_data  = object$data_info$stan_data,
    facet_cols = object$spec$facets_main
  )
  
  facet_infit_outfit(res_df, facet)
}
```

facet_infit_outfit_all(res_df) computes infit/outfit across all facets.

add_fit_cri_bayes() and add_zstd() add Bayesian CrIs and z-standardized fit stats.

⸻

8.5 Bias summaries – summarise_bias() (wrapper) + summarise_bias_ir() (existing)

User-facing API:

```r
summarise_bias <- function(object, facet = NULL, ...) {
  # If facet is NULL and spec has only one bias facet, use that.
  # If facet == "rater:item", call summarise_bias_ir(...)
  # Future: extend to other bias facets.
}
```

Internally, for "rater:item":

```r
bias_rc <- summarise_bias_ir(
  object$fit,
  item_labels  = object$data_info$facet_labels$item,
  rater_labels = object$data_info$facet_labels$rater,
  prob_thresh  = 0.95
)
```

⸻

8.6 Residual correlations – residual_cor_criteria()

New helper for analyzing correlations of residuals across criteria/items.

Function:

```r
residual_cor_criteria <- function(resid_tbl) {
  # 1. Pivot to wide: one row per person (or person–rater), one column per item.
  # 2. Compute correlation matrix across item columns.
}
```


Basic implementation idea:

```r
crit_wide <- resid_tbl %>%
  dplyr::select(person, item, resid) %>%
  dplyr::distinct() %>%
  dplyr::mutate(item = paste0("item_", item)) %>%
  tidyr::pivot_wider(
    names_from  = item,
    values_from = resid
  )

crit_cor <- stats::cor(
  dplyr::select(crit_wide, -person),
  use = "pairwise.complete.obs"
)
```

residual_cor_criteria() returns crit_cor (matrix) or a tidy form.

⸻

9. Design decisions and rationale (for beginners)
	1.	Why formula syntax?
R users are familiar with y ~ x1 + x2. We reuse that mental model but reinterpret RHS terms as facets, not continuous predictors.
	2.	Why integer indices?
Stan likes integer indices for categorical facets. Mapping "R01" → 1, "R02" → 2 is standard in hierarchical models.
	3.	Why sum-to-zero constraints?
Rasch/MFRM parameters are only identified up to a linear shift (we can add a constant to all person abilities and subtract it from item difficulties). We therefore:
	•	Mean-center each main facet’s parameters (sum-to-zero across levels).
	•	This makes scales comparable and avoids non-identifiability.
	•	The person scale can still be interpreted relative to the centered item/rater scales.
	4.	Why only rating-scale Rasch first?
Implementing everything (partial credit, GRM, drift, multidimensional) at once is complicated. Starting with one clean, well-understood model lets us:
	•	Test infrastructure
	•	Debug Stan code generation
	•	Design a stable R API
	5.	Why a small prior vocabulary?
brms supports many classes, but that’s overwhelming early on. MFRM use cases mostly need priors on:
	•	person ability
	•	item difficulty
	•	rater severity
	•	thresholds
	•	bias terms
So a small, Rasch-specific dictionary is enough.
	6.	Why caching?
Compiling Stan models is expensive. In research, you often re-run the same model with minor data changes. Caching + hashing the Stan code saves time and frustration.
	7.	Why S3 classes (bmfrm_fit)?
S3 is simple and familiar:
	•	fit <- bmfrm(...)
	•	summary(fit), plot(fit), stancode(fit), standata(fit)
Internally, bmfrm_fit stores a cmdstanr fit (R6), but the user interacts through simple S3 generics.

⸻

10. Implementation roadmap (v0.4 – leveraging existing code)

A developer can follow this order:
	1.	Package skeleton
	•	usethis::create_package("BayesMFRM")
	•	Add imports: cmdstanr, posterior, dplyr, tibble, rlang, glue, digest, tidyr, purrr.
	2.	Core helpers
	•	Implement prior() and basic validation.
	•	Implement parse_bmfrm_formula() for simple cases (person + item + rater + one bias).
	3.	Data prep
	•	Implement prepare_data_bmfrm():
	•	stan_data, facet_labels, data_clean.
	4.	Stan template
	•	Save base template as a string or .stan file with placeholders.
	•	Implement build_stan_code(spec, priors) for a 3-facet model.
	5.	Compile & fit
	•	Implement compile_bmfrm().
	•	Reuse fit_cmdstan_cached().
	•	Implement bmfrm() itself (wire everything together).
	•	Implement postprocess_bmfrm() to construct bmfrm_fit.
	6.	S3 methods & utilities
	•	Define "bmfrm_fit" class.
	•	Implement print.bmfrm_fit(), summary.bmfrm_fit(), stancode.bmfrm_fit(), standata.bmfrm_fit().
	•	Implement facet_summary(), fair_scores().
	•	Implement simple wrapper summarise_bias() around summarise_bias_ir().
	•	Implement residual_cor_criteria().
	7.	Testing
	•	Simulate small Rasch/MFRM data sets.
	•	Fit the model with bmfrm() and verify:
	•	Parameter recovery.
	•	Diagnostics (Rhat, ESS).
	•	Reasonable fair-score behavior.
	•	Reasonable infit/outfit.
	8.	Documentation
	•	Write vignettes:
	•	“Basic 3-facet MFRM”
	•	“Adding a rater×item bias facet”
	•	“Inspecting Stan code and priors”
	9.	Future extensions
	•	Partial-credit / GRM family.
	•	Additional facets (task, interlocutor).
	•	Additional bias facets (rater:task, rater:interlocutor).
	•	Multidimensional theta.



```r
ratings <- readr::read_csv("ratings.csv")

fit <- bmfrm(
  score ~ person + item + rater + rater:item,
  data = ratings,
  K    = 6
)

summary(fit)
facet_summary(fit, "rater")
fair_scores(fit)

# ppcheck
pp_check(fit)

# residual analysis
resid_tbl <- residuals(fit)
rater_fit <- facet_fit(fit, "rater")
item_fit  <- facet_fit(fit, "item")
bias_tab  <- summarise_bias(fit)
crit_cor  <- residual_cor_criteria(resid_tbl)
```
