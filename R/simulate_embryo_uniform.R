
R version 4.5.3 (2026-03-11 ucrt) -- "Reassured Reassurer"
Copyright (C) 2026 The R Foundation for Statistical Computing
Platform: x86_64-w64-mingw32/x64

R is free software and comes with ABSOLUTELY NO WARRANTY.
You are welcome to redistribute it under certain conditions.
Type 'license()' or 'licence()' for distribution details.

R is a collaborative project with many contributors.
Type 'contributors()' for more information and
'citation()' on how to cite R or R packages in publications.

Type 'demo()' for some demos, 'help()' for on-line help, or
'help.start()' for an HTML browser interface to help.
Type 'q()' to quit R.

> # simulate_embryo_uniform.R
> # Scenario A: all chromosomes have the same true mosaicism level
> 
> simulate_embryo_uniform <- function(true_mosaicism, n_chr = 24, coverage = 10,
+                                     noise_sd_base = 0.05, pcr_duplicates = 1.2,
+                                     batch_sd = 0.05, correlated_noise = TRUE) {
+   
+   # True log2 ratio (same for all chromosomes)
+   true_log2 <- log2(1 + true_mosaicism / 2)
+   
+   # Independent noise per chromosome
+   chr_noise <- rnorm(n_chr, mean = 0, 
+                      sd = noise_sd_base * pcr_duplicates / sqrt(coverage / 10))
+   
+   # Optional shared batch effect
+   batch_effect <- if (correlated_noise) rnorm(1, mean = 0, sd = batch_sd) else 0
+   
+   # Observed log2 ratios
+   observed <- rep(true_log2, n_chr) + chr_noise + batch_effect
+   
+   return(observed)
+ }
