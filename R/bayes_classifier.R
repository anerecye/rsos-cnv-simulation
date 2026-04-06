# ============================================================
# R/bayes_classifier.R
#
# Bayesian posterior classifier for mosaic CNV detection in PGT-A.
#
# Implements Section 2.5 and 2.6 of the paper:
#   - Beta prior on true mosaicism M
#   - Gaussian likelihood from the sequencing noise model
#   - Posterior P(M >= threshold | R_obs) via adaptive quadrature
#   - Age-stratified priors (Section 2.6)
#
# Main exported functions:
#   posterior_prob_mosaic()   — core single-observation classifier
#   classify_bayes()          — vectorised wrapper
#   make_age_prior()          — age-stratified Beta prior parameters
#   eval_classifier()         — sensitivity / FPR for one prior
#   build_table3()            — reproduce Table 3 across all priors
#
# Dependencies: base R only (uses integrate()).
# ============================================================


# ------------------------------------------------------------------
# 1. Noise model (mirrors simulate_cnv.R)
# ------------------------------------------------------------------

#' Log2-ratio likelihood given true mosaicism M
#'
#' The noise model: observed log2(CR) ~ N(log2(1 + M/2), sigma^2)
#' where sigma = noise_sd_base * pcr_duplicates / sqrt(coverage / 10)
#'
#' @param r_obs  Observed log2 ratio (scalar)
#' @param M      True mosaicism (scalar in [0,1])
#' @param coverage Sequencing coverage
#' @param noise_sd_base Baseline SD at 10x (default 0.05)
#' @param pcr_duplicates PCR multiplier (default 1.2)
#' @return Gaussian likelihood value (unnormalised over M)

likelihood <- function(M, r_obs,
                       coverage      = 10,
                       noise_sd_base = 0.05,
                       pcr_duplicates = 1.2) {
  true_log2r <- log2(1 + M / 2)
  sigma <- noise_sd_base * pcr_duplicates / sqrt(coverage / 10)
  dnorm(r_obs, mean = true_log2r, sd = sigma)
}


# ------------------------------------------------------------------
# 2. Prior specification
# ------------------------------------------------------------------

#' Beta prior density on M
#'
#' @param M     True mosaicism (scalar or vector in [0,1])
#' @param alpha Beta shape parameter alpha
#' @param beta  Beta shape parameter beta
#' @return Prior density

prior_beta <- function(M, alpha, beta) {
  dbeta(M, shape1 = alpha, shape2 = beta)
}


#' Build age-stratified Beta prior parameters (Section 2.6)
#'
#' Priors derived from large PGT-A cohorts (Viotti et al. 2021;
#' Capalbo et al. 2021) via method of moments with concentration
#' alpha + beta = 20. These should be treated as lower bounds on
#' true prevalence (see Section 2.6 caveat).
#'
#' @param age_group One of "<35", "35-40", ">40"
#' @return Named list with alpha, beta, mean_M, label

make_age_prior <- function(age_group = c("<35", "35-40", ">40")) {
  age_group <- match.arg(age_group)

  params <- list(
    "<35"   = list(alpha = 1.2,  beta = 18.8, mean_M = 0.06,  label = "Age <35"),
    "35-40" = list(alpha = 2.4,  beta = 17.6, mean_M = 0.12,  label = "Age 35-40"),
    ">40"   = list(alpha = 4.0,  beta = 16.0, mean_M = 0.20,  label = "Age >40")
  )
  params[[age_group]]
}


# ------------------------------------------------------------------
# 3. Core posterior computation
# ------------------------------------------------------------------

#' Posterior probability that M >= m_threshold given R_obs
#'
#' P(M >= m_threshold | R_obs) = integral_{m_threshold}^{1} L(R|M) * pi(M) dM
#'                               / integral_{0}^{1}           L(R|M) * pi(M) dM
#'
#' Integration via adaptive quadrature (base::integrate, tol = 1e-6).
#' Returns NA with a warning if either integral fails to converge.
#'
#' @param r_obs         Observed log2 ratio
#' @param alpha         Beta prior alpha
#' @param beta          Beta prior beta
#' @param m_threshold   Clinical relevance threshold (default 0.10)
#' @param coverage      Sequencing coverage
#' @param noise_sd_base Baseline noise SD
#' @param pcr_duplicates PCR multiplier
#' @return Posterior probability in [0, 1]

posterior_prob_mosaic <- function(r_obs,
                                  alpha          = 2,
                                  beta           = 10,
                                  m_threshold    = 0.10,
                                  coverage       = 10,
                                  noise_sd_base  = 0.05,
                                  pcr_duplicates = 1.2) {

  integrand <- function(M) {
    likelihood(M, r_obs,
               coverage       = coverage,
               noise_sd_base  = noise_sd_base,
               pcr_duplicates = pcr_duplicates) *
      prior_beta(M, alpha, beta)
  }

  # Marginal (denominator)
  denom <- tryCatch(
    integrate(integrand, lower = 0, upper = 1,
              rel.tol = 1e-6, abs.tol = 1e-9)$value,
    error = function(e) NA_real_
  )

  if (is.na(denom) || denom <= 0) return(NA_real_)

  # Numerator: P(M >= m_threshold | data)
  numer <- tryCatch(
    integrate(integrand, lower = m_threshold, upper = 1,
              rel.tol = 1e-6, abs.tol = 1e-9)$value,
    error = function(e) NA_real_
  )

  if (is.na(numer)) return(NA_real_)

  min(max(numer / denom, 0), 1)   # clamp to [0, 1]
}


# ------------------------------------------------------------------
# 4. Vectorised classifier
# ------------------------------------------------------------------

#' Classify a vector of observed log2 ratios using the Bayesian rule
#'
#' An embryo is flagged as "mosaic-relevant" if
#' P(M >= m_threshold | R_obs) > posterior_cutoff (default 0.5).
#'
#' @param r_obs_vec      Numeric vector of observed log2 ratios
#' @param alpha          Beta prior alpha
#' @param beta           Beta prior beta
#' @param m_threshold    Clinical relevance threshold
#' @param posterior_cutoff Decision boundary (default 0.5)
#' @param coverage       Sequencing coverage
#' @param noise_sd_base  Baseline noise SD
#' @param pcr_duplicates PCR multiplier
#' @return Data frame: r_obs, posterior_prob, bayes_call ("mosaic-relevant" / "normal")

classify_bayes <- function(r_obs_vec,
                           alpha           = 2,
                           beta            = 10,
                           m_threshold     = 0.10,
                           posterior_cutoff = 0.5,
                           coverage        = 10,
                           noise_sd_base   = 0.05,
                           pcr_duplicates  = 1.2) {

  probs <- vapply(r_obs_vec, posterior_prob_mosaic,
                  numeric(1),
                  alpha          = alpha,
                  beta           = beta,
                  m_threshold    = m_threshold,
                  coverage       = coverage,
                  noise_sd_base  = noise_sd_base,
                  pcr_duplicates = pcr_duplicates)

  data.frame(
    r_obs         = r_obs_vec,
    posterior_prob = probs,
    bayes_call    = ifelse(probs > posterior_cutoff,
                           "mosaic-relevant", "normal"),
    stringsAsFactors = FALSE
  )
}


# ------------------------------------------------------------------
# 5. Evaluate classifier performance for one prior
# ------------------------------------------------------------------

#' Compute sensitivity and FPR for one prior specification
#'
#' Simulates n_replicates embryos at M = true_m (mosaic) and M = 0
#' (diploid), classifies with both the fixed threshold and the
#' Bayesian rule, returns a one-row summary data frame.
#'
#' @param alpha           Beta prior alpha
#' @param beta            Beta prior beta
#' @param prior_label     Human-readable label for the prior
#' @param true_m          True mosaicism level to test sensitivity
#' @param coverage        Sequencing coverage
#' @param low_thr         Fixed lower threshold (for comparison)
#' @param n_replicates    Number of simulated embryos per class
#' @param noise_sd_base   Baseline noise SD
#' @param pcr_duplicates  PCR multiplier
#' @param posterior_cutoff Bayesian decision boundary
#' @param m_threshold     Clinical relevance threshold for posterior
#' @return One-row data frame

eval_classifier <- function(alpha,
                            beta,
                            prior_label,
                            true_m          = 0.15,
                            coverage        = 10,
                            low_thr         = 0.15,
                            n_replicates    = 5000,
                            noise_sd_base   = 0.05,
                            pcr_duplicates  = 1.2,
                            posterior_cutoff = 0.5,
                            m_threshold     = 0.10) {

  sigma <- noise_sd_base * pcr_duplicates / sqrt(coverage / 10)

  # --- Simulate log2 ratios ---
  # Diploid embryos (M = 0): true log2 CR = 0
  r_diploid <- rnorm(n_replicates, mean = 0, sd = sigma)

  # Mosaic embryos (M = true_m)
  r_mosaic  <- rnorm(n_replicates,
                     mean = log2(1 + true_m / 2),
                     sd   = sigma)

  # --- Fixed threshold classification ---
  fixed_fp  <- mean(abs(r_diploid) >= low_thr)
  fixed_sens <- mean(abs(r_mosaic) >= low_thr)

  # --- Bayesian classification ---
  bayes_diploid <- classify_bayes(r_diploid,
                                  alpha            = alpha,
                                  beta             = beta,
                                  m_threshold      = m_threshold,
                                  posterior_cutoff = posterior_cutoff,
                                  coverage         = coverage,
                                  noise_sd_base    = noise_sd_base,
                                  pcr_duplicates   = pcr_duplicates)

  bayes_mosaic  <- classify_bayes(r_mosaic,
                                  alpha            = alpha,
                                  beta             = beta,
                                  m_threshold      = m_threshold,
                                  posterior_cutoff = posterior_cutoff,
                                  coverage         = coverage,
                                  noise_sd_base    = noise_sd_base,
                                  pcr_duplicates   = pcr_duplicates)

  bayes_fp   <- mean(bayes_diploid$bayes_call == "mosaic-relevant", na.rm = TRUE)
  bayes_sens <- mean(bayes_mosaic$bayes_call  == "mosaic-relevant", na.rm = TRUE)

  data.frame(
    prior_label   = prior_label,
    alpha         = alpha,
    beta          = beta,
    mean_M        = alpha / (alpha + beta),
    sensitivity_fixed = fixed_sens,
    fpr_fixed         = fixed_fp,
    sensitivity_bayes = bayes_sens,
    fpr_bayes         = bayes_fp,
    stringsAsFactors  = FALSE
  )
}


# ------------------------------------------------------------------
# 6. Reproduce Table 3 across all priors from the paper
# ------------------------------------------------------------------

#' Build Table 3: Bayesian classifier performance by prior
#'
#' Priors tested (Section 2.5):
#'   Sceptical  Beta(1, 49)  mean ≈ 0.02
#'   Beta(0.5, 10)           mean ≈ 0.048
#'   Uniform    Beta(1, 1)   mean = 0.50
#'   Main       Beta(2, 10)  mean ≈ 0.167
#'   Beta(5, 20)             mean = 0.20
#'
#' @param true_m       True mosaicism level (default 0.15)
#' @param coverage     Sequencing coverage (default 10)
#' @param n_replicates Replicates per class
#' @return Data frame matching Table 3 structure

build_table3 <- function(true_m       = 0.15,
                         coverage     = 10,
                         n_replicates = 5000) {

  priors <- list(
    list(alpha = 1,   beta = 49, label = "Sceptical (mean=0.02)"),
    list(alpha = 0.5, beta = 10, label = "Beta(0.5, 10)"),
    list(alpha = 1,   beta = 1,  label = "Beta(1,1) uniform"),
    list(alpha = 2,   beta = 10, label = "Beta(2, 10) [main]"),
    list(alpha = 5,   beta = 20, label = "Beta(5, 20)")
  )

  results <- lapply(priors, function(p) {
    cat("  Evaluating prior:", p$label, "\n")
    eval_classifier(
      alpha        = p$alpha,
      beta         = p$beta,
      prior_label  = p$label,
      true_m       = true_m,
      coverage     = coverage,
      n_replicates = n_replicates
    )
  })

  table3 <- do.call(rbind, results)

  # Add fixed threshold row for direct comparison (same format)
  # Already computed inside eval_classifier; pull from first row
  fixed_row <- data.frame(
    prior_label       = "Fixed threshold 0.15/0.30",
    alpha             = NA,
    beta              = NA,
    mean_M            = NA,
    sensitivity_fixed = table3$sensitivity_fixed[1],
    fpr_fixed         = table3$fpr_fixed[1],
    sensitivity_bayes = table3$sensitivity_fixed[1],  # same value
    fpr_bayes         = table3$fpr_fixed[1],
    stringsAsFactors  = FALSE
  )

  table3_full <- rbind(table3, fixed_row)
  rownames(table3_full) <- NULL
  table3_full
}
