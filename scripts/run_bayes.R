# ============================================================
# scripts/run_bayes.R
#
# Runs the Bayesian classifier analysis (Sections 2.5, 2.6, 3.4).
#
# Outputs:
#   results/table3_bayes_priors.csv       — Table 3 (sensitivity/FPR by prior)
#   results/posterior_probs_embryos.csv   — per-embryo posterior probabilities
#   results/table3_age_stratified.csv     — age-stratified prior results
#   figures/fig_bayes_sensitivity_fpr.png — sensitivity vs FPR figure
#
# Called from run_all.R or sourced interactively.
# ============================================================

source("R/bayes_classifier.R")

set.seed(2025)
dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)


# ── 1. Table 3: sensitivity / FPR across all priors ───────────────────────

cat("\n[1/4] Building Table 3 (Bayesian vs fixed threshold, all priors)...\n")

table3 <- build_table3(
  true_m       = 0.15,
  coverage     = 10,
  n_replicates = 5000
)

# Round for readability
table3_rounded <- table3
table3_rounded[, c("sensitivity_fixed", "fpr_fixed",
                    "sensitivity_bayes", "fpr_bayes")] <-
  round(table3_rounded[, c("sensitivity_fixed", "fpr_fixed",
                            "sensitivity_bayes", "fpr_bayes")], 3)

write.csv(table3_rounded, "results/table3_bayes_priors.csv", row.names = FALSE)

cat("\nTable 3:\n")
print(table3_rounded[, c("prior_label", "mean_M",
                          "sensitivity_bayes", "fpr_bayes",
                          "sensitivity_fixed", "fpr_fixed")])


# ── 2. Per-embryo posterior probabilities ─────────────────────────────────

cat("\n[2/4] Computing per-embryo posterior probabilities...\n")

# Simulate a representative set of embryos at several true M levels
# using main prior Beta(2, 10) and default coverage 10x

m_levels    <- c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30)
n_per_level <- 500   # 500 per M level: tractable, enough for CSV illustration

sigma_default <- 0.05 * 1.2 / sqrt(10 / 10)   # 0.06

embryo_rows <- lapply(m_levels, function(m) {
  r_obs <- rnorm(n_per_level, mean = log2(1 + m / 2), sd = sigma_default)
  data.frame(
    true_m  = m,
    r_obs   = r_obs,
    stringsAsFactors = FALSE
  )
})
embryo_df <- do.call(rbind, embryo_rows)

cat("  Computing posteriors for", nrow(embryo_df), "embryos...\n")

bayes_out <- classify_bayes(
  r_obs_vec       = embryo_df$r_obs,
  alpha           = 2,
  beta            = 10,
  m_threshold     = 0.10,
  posterior_cutoff = 0.5,
  coverage        = 10
)

embryo_df$posterior_prob  <- bayes_out$posterior_prob
embryo_df$bayes_call      <- bayes_out$bayes_call
embryo_df$fixed_call      <- ifelse(abs(embryo_df$r_obs) >= 0.15,
                                     "mosaic-relevant", "normal")

write.csv(embryo_df, "results/posterior_probs_embryos.csv", row.names = FALSE)
cat("  Saved: results/posterior_probs_embryos.csv\n")


# ── 3. Age-stratified priors (Section 2.6) ────────────────────────────────

cat("\n[3/4] Age-stratified priors (Section 2.6)...\n")

age_groups <- c("<35", "35-40", ">40")

age_results <- lapply(age_groups, function(ag) {
  p <- make_age_prior(ag)
  cat("  Age group:", p$label, "| Beta(", p$alpha, ",", p$beta, ")\n")
  res <- eval_classifier(
    alpha        = p$alpha,
    beta         = p$beta,
    prior_label  = p$label,
    true_m       = 0.15,
    coverage     = 10,
    n_replicates = 5000
  )
  res
})

age_table <- do.call(rbind, age_results)
age_table[, c("sensitivity_fixed", "fpr_fixed",
              "sensitivity_bayes", "fpr_bayes")] <-
  round(age_table[, c("sensitivity_fixed", "fpr_fixed",
                       "sensitivity_bayes", "fpr_bayes")], 3)

write.csv(age_table, "results/table3_age_stratified.csv", row.names = FALSE)

cat("\nAge-stratified results:\n")
print(age_table[, c("prior_label", "mean_M",
                     "sensitivity_bayes", "fpr_bayes")])


# ── 4. Figure: sensitivity vs FPR ─────────────────────────────────────────

cat("\n[4/4] Generating sensitivity vs FPR figure...\n")

# Build ROC-style curve by sweeping posterior_cutoff from 0.05 to 0.95
# for the main prior (Beta(2,10), M = 0.15, coverage = 10x)

n_curve     <- 2000   # embryos per class for the curve
r_diploid_c <- rnorm(n_curve, mean = 0,              sd = sigma_default)
r_mosaic_c  <- rnorm(n_curve, mean = log2(1 + 0.15 / 2), sd = sigma_default)

# Pre-compute posteriors once (expensive), then sweep cutoffs
post_diploid <- vapply(r_diploid_c, posterior_prob_mosaic,
                       numeric(1), alpha = 2, beta = 10,
                       m_threshold = 0.10, coverage = 10)

post_mosaic  <- vapply(r_mosaic_c, posterior_prob_mosaic,
                       numeric(1), alpha = 2, beta = 10,
                       m_threshold = 0.10, coverage = 10)

cutoffs <- seq(0.05, 0.95, by = 0.025)

roc_rows <- lapply(cutoffs, function(co) {
  data.frame(
    cutoff      = co,
    sensitivity = mean(post_mosaic  > co, na.rm = TRUE),
    fpr         = mean(post_diploid > co, na.rm = TRUE),
    method      = "Bayesian Beta(2,10)"
  )
})
roc_df <- do.call(rbind, roc_rows)

# Fixed threshold "curve" — single point for each low_thr value
fixed_pts <- data.frame(
  cutoff      = c(0.10, 0.15, 0.20),
  sensitivity = c(
    mean(abs(r_mosaic_c)  >= 0.10),
    mean(abs(r_mosaic_c)  >= 0.15),
    mean(abs(r_mosaic_c)  >= 0.20)
  ),
  fpr = c(
    mean(abs(r_diploid_c) >= 0.10),
    mean(abs(r_diploid_c) >= 0.15),
    mean(abs(r_diploid_c) >= 0.20)
  ),
  method = "Fixed threshold"
)

# ---- plot with base R (no ggplot2 dependency) ----

png("figures/fig_bayes_sensitivity_fpr.png",
    width = 700, height = 600, res = 120)

par(mar = c(5, 5, 4, 2))

plot(roc_df$fpr, roc_df$sensitivity,
     type = "l", lwd = 2.5, col = "#2166AC",
     xlim = c(0, 0.5), ylim = c(0, 1),
     xlab = "False positive rate (diploid embryos called mosaic-relevant)",
     ylab = "Sensitivity (mosaic embryos correctly flagged)",
     main = "Bayesian vs fixed threshold classifier\nM = 15%, coverage 10×, prior Beta(2,10)",
     cex.main = 0.95, cex.lab = 0.9)

grid(col = "grey85", lty = 1)

# Annotate the operating point used in the paper (cutoff = 0.5)
op_idx <- which.min(abs(roc_df$cutoff - 0.5))
points(roc_df$fpr[op_idx], roc_df$sensitivity[op_idx],
       pch = 19, cex = 1.5, col = "#2166AC")
text(roc_df$fpr[op_idx] + 0.02, roc_df$sensitivity[op_idx] - 0.04,
     "p > 0.5\n(paper)", cex = 0.75, col = "#2166AC", adj = 0)

# Fixed threshold points
points(fixed_pts$fpr, fixed_pts$sensitivity,
       pch = 17, cex = 1.4, col = "#D6604D")
text(fixed_pts$fpr + 0.015, fixed_pts$sensitivity + 0.03,
     paste0("thr=", fixed_pts$cutoff),
     cex = 0.72, col = "#D6604D", adj = 0)

# Diagonal reference
abline(a = 0, b = 1, lty = 2, col = "grey60")

legend("bottomright",
       legend = c("Bayesian Beta(2,10)", "Fixed threshold points"),
       col    = c("#2166AC", "#D6604D"),
       lty    = c(1, NA), pch = c(19, 17),
       lwd    = c(2.5, NA), pt.cex = 1.3,
       bty    = "n", cex = 0.85)

dev.off()
cat("  Saved: figures/fig_bayes_sensitivity_fpr.png\n")


# ── Summary ───────────────────────────────────────────────────────────────

cat("\n✅ Bayesian analysis complete.\n")
cat("  results/table3_bayes_priors.csv\n")
cat("  results/posterior_probs_embryos.csv\n")
cat("  results/table3_age_stratified.csv\n")
cat("  figures/fig_bayes_sensitivity_fpr.png\n\n")

cat("Key result (main prior Beta(2,10), M=15%, 10x):\n")
main_row <- table3_rounded[table3_rounded$prior_label == "Beta(2, 10) [main]", ]
cat(sprintf("  Sensitivity: %.1f%%  |  FPR: %.1f%%\n",
            main_row$sensitivity_bayes * 100,
            main_row$fpr_bayes * 100))
fixed_row_out <- table3_rounded[table3_rounded$prior_label == "Fixed threshold 0.15/0.30", ]
cat(sprintf("  vs fixed:    %.1f%%  |  FPR: %.1f%%\n",
            fixed_row_out$sensitivity_bayes * 100,
            fixed_row_out$fpr_bayes * 100))
