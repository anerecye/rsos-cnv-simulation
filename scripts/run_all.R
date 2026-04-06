
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

cat("\n✅ All simulations complete!\n")
cat("Results saved in results/\n")

