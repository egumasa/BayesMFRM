# Example usage of BayesMFRM package
# This demonstrates the intended workflow from the specification

# Note: This requires CmdStan to be installed
# install_cmdstan_if_needed()

library(BayesMFRM)

# Create example data (normally you would load your own data)
set.seed(42)
ratings <- data.frame(
  person = rep(paste0("P", 1:20), each = 12),
  item = rep(paste0("Item", 1:4), times = 60), 
  rater = rep(paste0("Rater", 1:3), times = 80),
  score = sample(1:6, 240, replace = TRUE, prob = c(0.05, 0.1, 0.2, 0.3, 0.25, 0.1))
)

# Check data suitability
data_check <- check_data_suitability(ratings, "score", c("person", "item", "rater"))
print(data_check)

# Example 1: Basic 3-facet model
cat("\\n=== Example 1: Basic Model ===\\n")
fit1 <- bmfrm(
  score ~ person + item + rater,
  data = ratings,
  K = 6,
  chains = 2,  # Use fewer chains for demo
  iter = 1000   # Use fewer iterations for demo
)

# View results
print(fit1)
summary(fit1)

# Example 2: Model with bias term
cat("\\n=== Example 2: Model with Bias ===\\n") 
fit2 <- bmfrm(
  score ~ person + item + rater + rater:item,
  data = ratings,
  chains = 2,
  iter = 1000
)

# Model checking
pp_check(fit2)

# Residual analysis
resid_tbl <- residuals(fit2)
rater_fit <- facet_fit(fit2, "rater")
item_fit <- facet_fit(fit2, "item")

print("Rater fit statistics:")
print(rater_fit)

# Bias analysis
bias_tab <- summarise_bias(fit2)
print("Bias analysis:")
print(bias_tab)

# Fair scores
fair_scores_result <- fair_scores(fit2)
print("Fair scores (person-level summary):")
print(head(fair_scores_result))

# Cross-criteria correlations
crit_cor <- residual_cor_criteria(resid_tbl)
print(crit_cor)

cat("\\nWorkflow complete! Check results above.\\n")
