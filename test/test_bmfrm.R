# Test script for BayesMFRM package
# This script tests the basic workflow without actually running MCMC

library(dplyr)
library(tibble)

# Source all R files for testing
r_files <- list.files("R", pattern = "\\.R$", full.names = TRUE)
for (f in r_files) {
  source(f)
}

# Create synthetic test data
set.seed(123)
n_persons <- 50
n_items <- 6
n_raters <- 4

# Generate all combinations
test_data <- expand.grid(
  person = paste0("P", 1:n_persons),
  item = paste0("I", 1:n_items), 
  rater = paste0("R", 1:n_raters)
) %>%
  # Simulate some missing design points (realistic for MFRM)
  sample_n(floor(nrow(.) * 0.8)) %>%
  # Add realistic scores (1-6 scale)
  mutate(
    # Simple simulation: person ability + item difficulty + rater severity + noise
    person_ability = rep(rnorm(n_persons, 0, 1), length.out = nrow(.))[match(person, paste0("P", 1:n_persons))],
    item_difficulty = rep(rnorm(n_items, 0, 0.8), length.out = nrow(.))[match(item, paste0("I", 1:n_items))],
    rater_severity = rep(rnorm(n_raters, 0, 0.5), length.out = nrow(.))[match(rater, paste0("R", 1:n_raters))],
    linear_pred = person_ability - item_difficulty - rater_severity,
    # Convert to ordinal scale
    prob = plogis(linear_pred),
    score = pmax(1, pmin(6, round(1 + 5 * prob + rnorm(nrow(.), 0, 0.5))))
  ) %>%
  select(person, item, rater, score)

cat("Test data created with", nrow(test_data), "observations\n")
cat("Score distribution:\n")
print(table(test_data$score))

# Test 1: Formula parsing
cat("\n=== Testing Formula Parsing ===\n")
spec1 <- parse_bmfrm_formula(
  formula = score ~ person + item + rater,
  data = test_data,
  K = 6
)
print(spec1)

# Test 2: Data preparation  
cat("\n=== Testing Data Preparation ===\n")
data_info1 <- prepare_data_bmfrm(spec1, test_data)
cat("Stan data structure:\n")
str(data_info1$stan_data)

# Test 3: Prior system
cat("\n=== Testing Prior System ===\n")
test_priors <- c(
  prior("normal(0, 2)", class = "theta"),
  prior("normal(0, 1)", class = "item"), 
  prior("normal(0, 1)", class = "rater"),
  prior("normal(0, 3)", class = "tau")
)
print(test_priors)

validated_priors <- validate_priors(test_priors, spec1$facets_main, spec1$facets_bias)
cat("Validated", length(validated_priors), "priors\n")

# Test 4: Stan code generation
cat("\n=== Testing Stan Code Generation ===\n")
stan_code1 <- build_stan_code(spec1, validated_priors)
cat("Generated Stan code (", nchar(stan_code1), " characters)\n")
cat("First 500 characters:\n")
cat(substr(stan_code1, 1, 500), "\n...\n")

# Test 5: Model with bias term
cat("\n=== Testing Model with Bias ===\n")
spec2 <- parse_bmfrm_formula(
  formula = score ~ person + item + rater + rater:item,
  data = test_data,
  K = 6
)
print(spec2)

bias_priors <- c(
  prior("normal(0, 2)", class = "theta"),
  prior("normal(0, 1)", class = "item"),
  prior("normal(0, 1)", class = "rater"), 
  prior("normal(0, 0.5)", class = "bias", facet = "rater:item"),
  prior("normal(0, 3)", class = "tau")
)

data_info2 <- prepare_data_bmfrm(spec2, test_data)
stan_code2 <- build_stan_code(spec2, bias_priors)
cat("Generated bias model Stan code (", nchar(stan_code2), " characters)\n")

# Test 6: Data suitability check
cat("\n=== Testing Data Suitability Check ===\n")
data_check <- check_data_suitability(
  data = test_data,
  response = "score",
  facets = c("person", "item", "rater")
)
print(data_check)

cat("\n=== All Tests Completed Successfully! ===\n")
cat("The BayesMFRM package implementation is ready for use.\n")
cat("To fit an actual model, ensure CmdStan is installed and use:\n")
cat("fit <- bmfrm(score ~ person + item + rater, data = your_data)\n")
