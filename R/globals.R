# Declare column symbols used in dplyr/tidyr non-standard evaluation so
# R CMD check's codetools pass emits zero global-variable notes.
utils::globalVariables(c(
  "bias_mean", "bias_sd", "ci_lower", "ci_upper", "correlation", "facet",
  "facet_value", "flag_bias", "i", "idx", "infit_hi", "infit_lo", "infit_msq",
  "infit_ppp", "infit_rep_hi", "infit_rep_lo", "infit_rep_mean", "infit_zstd",
  "item", "item1", "item2", "label", "lower_95", "model", "mu_hat", "n_obs",
  "outfit_hi", "outfit_lo", "outfit_msq", "outfit_ppp", "outfit_rep_hi",
  "outfit_rep_lo", "outfit_rep_mean", "outfit_zstd", "param", "prob_neg",
  "prob_pos", "r", "rater", "resid", "sigma2_hat", "upper_95", "value",
  "variable", "weight", "x", "z"
))
