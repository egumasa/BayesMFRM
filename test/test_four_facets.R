# Test script to reproduce four-facets-main-effect_without_mono-moderate.stan
# Tests the formula: score ~ person + criteria + rater + interlocutor

library(dplyr)
library(tibble)

# Source all R files for testing
r_files <- list.files("R", pattern = "\\.R$", full.names = TRUE)
for (f in r_files) {
  source(f)
}

cat("=== Testing Four-Facet Model Generation ===\n")
cat("Target: Reproduce four-facets-main-effect_without_mono-moderate.stan\n")
cat("Formula: score ~ person + criteria + rater + interlocutor\n\n")

# 1. Create synthetic test data with four facets
set.seed(42)  # for reproducibility
n_persons <- 30
n_criteria <- 8
n_raters <- 5
n_interlocutors <- 4

cat("Creating synthetic data with:\n")
cat("  - Persons:", n_persons, "\n")
cat("  - Criteria:", n_criteria, "\n") 
cat("  - Raters:", n_raters, "\n")
cat("  - Interlocutors:", n_interlocutors, "\n\n")

# Generate all combinations but sample to create realistic sparsity
full_design <- expand.grid(
  person = paste0("P", sprintf("%02d", 1:n_persons)),
  criteria = paste0("C", sprintf("%02d", 1:n_criteria)), 
  rater = paste0("R", sprintf("%02d", 1:n_raters)),
  interlocutor = paste0("I", sprintf("%02d", 1:n_interlocutors))
)

# Sample about 60% of possible combinations (realistic for MFRM)
test_data <- full_design %>%
  sample_n(floor(nrow(.) * 0.6)) %>%
  # Add realistic scores (1-6 scale) based on Many-Facet Rasch Model
  mutate(
    # Simulate underlying abilities and difficulties
    person_ability = rep(rnorm(n_persons, 0, 1), length.out = nrow(.))[match(person, paste0("P", sprintf("%02d", 1:n_persons)))],
    criteria_difficulty = rep(rnorm(n_criteria, 0, 0.8), length.out = nrow(.))[match(criteria, paste0("C", sprintf("%02d", 1:n_criteria)))],
    rater_severity = rep(rnorm(n_raters, 0, 0.5), length.out = nrow(.))[match(rater, paste0("R", sprintf("%02d", 1:n_raters)))],
    interlocutor_effect = rep(rnorm(n_interlocutors, 0, 0.3), length.out = nrow(.))[match(interlocutor, paste0("I", sprintf("%02d", 1:n_interlocutors)))],
    
    # MFRM linear predictor: ability - difficulty - severity - effect
    linear_pred = person_ability - criteria_difficulty - rater_severity - interlocutor_effect,
    
    # Convert to 6-point ordinal scale with some noise
    prob = plogis(linear_pred),
    score = pmax(1, pmin(6, round(1 + 5 * prob + rnorm(nrow(.), 0, 0.4))))
  ) %>%
  select(person, criteria, rater, interlocutor, score)

cat("Generated", nrow(test_data), "observations\n")
cat("Score distribution:\n")
print(table(test_data$score))
cat("\n")

# 2. Test formula parsing for four facets
cat("=== Testing Formula Parsing ===\n")
spec_four_facets <- parse_bmfrm_formula(
  formula = score ~ person + criteria + rater + interlocutor,
  data = test_data,
  K = 6
)

print(spec_four_facets)
cat("\n")

# 3. Test data preparation
cat("=== Testing Data Preparation ===\n")
data_info <- prepare_data_bmfrm(spec_four_facets, test_data)
cat("Stan data structure:\n")
str(data_info$stan_data)
cat("\n")

# 4. Configure priors to match the example model
cat("=== Configuring Priors (matching example model) ===\n")
four_facet_priors <- list(
  # theta_raw ~ student_t(3, 0, 4)
  prior("student_t(3, 0, 4)", class = "theta"),
  
  # beta_raw ~ student_t(3, 0, 2) (criteria = item difficulty)
  prior("student_t(3, 0, 2)", class = "facet_criteria"), 
  
  # rater_raw ~ student_t(3, 0, 2)
  prior("student_t(3, 0, 2)", class = "rater"),
  
  # interlocutor_raw ~ student_t(3, 0, 2)
  prior("student_t(3, 0, 2)", class = "facet_interlocutor"),
  
  # tau_raw ~ student_t(3, 0, 2)
  prior("student_t(3, 0, 2)", class = "tau")
)

cat("Configured priors:\n")
for (p in four_facet_priors) {
  print(p)
}
cat("\n")

# Validate priors
validated_priors <- validate_priors(four_facet_priors, 
                                   spec_four_facets$facets_main, 
                                   spec_four_facets$facets_bias)
cat("Validated", length(validated_priors), "priors successfully\n\n")

# 5. Generate Stan code
cat("=== Generating Stan Code ===\n")
generated_stan_code <- build_stan_code(spec_four_facets, validated_priors)

cat("Generated Stan code (", nchar(generated_stan_code), " characters)\n")

# 6. Save generated code for comparison
output_file <- "test/generated_four_facets.stan"
writeLines(generated_stan_code, output_file)
cat("Saved generated code to:", output_file, "\n\n")

# 7. Basic validation checks
cat("=== Validation Checks ===\n")

# Check for key structural elements
checks <- list(
  "Four facets in data block" = all(c(
    grepl("J_person", generated_stan_code, fixed = TRUE),
    grepl("J_criteria", generated_stan_code, fixed = TRUE), 
    grepl("J_rater", generated_stan_code, fixed = TRUE),
    grepl("J_interlocutor", generated_stan_code, fixed = TRUE)
  )),
  
  "Parameter declarations present" = all(c(
    grepl("theta_raw", generated_stan_code, fixed = TRUE),
    grepl("criteria_raw", generated_stan_code, fixed = TRUE),
    grepl("rater_raw", generated_stan_code, fixed = TRUE),
    grepl("interlocutor_raw", generated_stan_code, fixed = TRUE),
    grepl("tau_raw", generated_stan_code, fixed = TRUE)
  )),
  
  "Student-t priors applied" = all(c(
    grepl("theta_raw ~ student_t(3, 0, 4)", generated_stan_code, fixed = TRUE),
    grepl("criteria_raw ~ student_t(3, 0, 2)", generated_stan_code, fixed = TRUE),
    grepl("rater_raw ~ student_t(3, 0, 2)", generated_stan_code, fixed = TRUE),
    grepl("interlocutor_raw ~ student_t(3, 0, 2)", generated_stan_code, fixed = TRUE),
    grepl("tau_raw ~ student_t(3, 0, 2)", generated_stan_code, fixed = TRUE)
  )),
  
  "Rating-scale model structure" = all(c(
    grepl("categorical_logit", generated_stan_code, fixed = TRUE),
    grepl("logits\\[k-1\\] \\+ eta - tau\\[k-1\\]", generated_stan_code, fixed = FALSE),
    grepl("log_lik", generated_stan_code, fixed = TRUE),
    grepl("y_rep", generated_stan_code, fixed = TRUE)
  )),
  
  "Fair scores calculation" = all(c(
    grepl("mu_fair", generated_stan_code, fixed = TRUE),
    grepl("eta_fair", generated_stan_code, fixed = TRUE)
  ))
)

# Report validation results
for (check_name in names(checks)) {
  status <- if (checks[[check_name]]) "✓ PASS" else "✗ FAIL"
  cat(sprintf("%-35s %s\n", check_name, status))
}

cat("\n=== Test Summary ===\n")
all_passed <- all(unlist(checks))
if (all_passed) {
  cat("🎉 All validation checks PASSED!\n")
  cat("The BayesMFRM package successfully generates Stan code\n")
  cat("that matches the structure of the four-facet example model.\n")
} else {
  cat("⚠️  Some validation checks FAILED.\n")
  cat("Review the generated code for issues.\n")
}

cat("\nGenerated file:", output_file, "\n")
cat("Compare with: resources/example_models/four-facets-main-effect_without_mono-moderate.stan\n")
