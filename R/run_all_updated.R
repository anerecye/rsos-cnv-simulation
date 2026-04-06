set.seed(2025)

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

source("scripts/run_scenario_single.R")
source("scripts/run_scenario_A.R")
source("scripts/run_scenario_B.R")

cat("Running single-chromosome simulation...\n")
single_res <- run_scenario_single()
write.csv(single_res, "results/single_chr_results.csv", row.names = FALSE)

cat("Running Scenario A (uniform)...\n")
scenario_A_res <- run_scenario_A()
write.csv(scenario_A_res, "results/scenario_A_results.csv", row.names = FALSE)

cat("Running Scenario B (single affected)...\n")
scenario_B_res <- run_scenario_B()
write.csv(scenario_B_res, "results/scenario_B_results.csv", row.names = FALSE)

single_subset <- single_res[single_res$cov == 10 & single_res$m %in% c(0, 0.15, 0.20), ]
single_normal <- single_subset$prop_normal

table5 <- data.frame(
  M = c(0, 0.15, 0.20),
  single_chr_normal = single_normal,
  scenario_A_normal = scenario_A_res$prop_normal,
  scenario_B_normal = scenario_B_res$prop_normal
)

write.csv(table5, "results/table5_comparison.csv", row.names = FALSE)

# ── Supplementary sensitivity analyses ─────────────────────────────────────

cat("\nRunning supplementary: binomial cell sampling sensitivity...\n")
source("scripts/sensitivity_binom_cells.R")

cat("\nRunning supplementary: PCR multiplier (eta) sensitivity...\n")
source("scripts/sensitivity_eta.R")

cat("\nRunning supplementary: Scenario B filter factor sensitivity...\n")
source("scripts/sensitivity_filter_factor.R")

# ── Bayesian classifier (Sections 2.5, 2.6, 3.4) ───────────────────────────
# NOTE: this step is the slowest (~10–20 min at n_replicates = 5000)
# because posterior_prob_mosaic() calls integrate() once per embryo.
# To run a quick check first, reduce n_replicates to 500 interactively:
#   source("R/bayes_classifier.R"); build_table3(n_replicates = 500)

cat("\nRunning Bayesian classifier analysis...\n")
source("scripts/run_bayes.R")

cat("\n✅ All simulations complete!\n")
cat("Main results:          results/\n")
cat("Supplementary results: results/binom_sensitivity.csv\n")
cat("                       results/sensitivity_eta.csv\n")
cat("                       results/filter_factor_sensitivity.csv\n")
cat("Bayesian results:      results/table3_bayes_priors.csv\n")
cat("                       results/posterior_probs_embryos.csv\n")
cat("                       results/table3_age_stratified.csv\n")
cat("Figures:               figures/\n")
