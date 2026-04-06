
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

> # simulate_cnv.R
> # Simulate observed log2 ratio for one CNV on one chromosome
> 
> simulate_cnv <- function(true_mosaicism, cnv_type = "dup",
+                          coverage = 10, noise_sd_base = 0.05,
+                          pcr_duplicates = 1.2) {
+   
+   # True copy ratio based on mosaicism level and CNV type
+   if (cnv_type == "dup") {
+     true_cr <- 1 + true_mosaicism / 2
+   } else if (cnv_type == "del") {
+     true_cr <- 1 - true_mosaicism / 2
+   } else {
+     stop("cnv_type must be 'dup' or 'del'")
+   }
+   
+   # Noise decreases with higher coverage
+   effective_noise_sd <- noise_sd_base * pcr_duplicates / sqrt(coverage / 10)
+   
+   # Add Gaussian noise
+   observed_cr <- true_cr + rnorm(1, mean = 0, sd = effective_noise_sd)
+   
+   # Log2 transformation for threshold application
+   log2_ratio <- log2(observed_cr)
+   
+   return(log2_ratio)
+ }
