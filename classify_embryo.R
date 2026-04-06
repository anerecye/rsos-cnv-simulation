
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

> # classify_embryo.R
> # Classify one chromosome based on log2 ratio and thresholds
> 
> classify_embryo <- function(log2r, low_thr, high_thr) {
+   ifelse(abs(log2r) < low_thr, "normal",
+          ifelse(abs(log2r) < high_thr, "mosaic", "aneuploid"))
+ }
