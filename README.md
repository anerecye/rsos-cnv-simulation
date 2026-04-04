# Borderline CNVs in PGT-A: a simulation study

## Overview

This repository contains the complete reproducible simulation code for the paper:

> **Title**: *Borderline copy number variations in preimplantation genetic testing: when signal becomes noise*  
> **Authors**: Anere Cye
> **Journal**: Royal Society Open Science (2026)

The code simulates embryo-level CNV detection under different mosaicism levels (0–30%), sequencing coverages (5×, 10×, 20×, 50×), and classification thresholds (0.10–0.20 low, 0.30–0.40 high). It quantifies the proportion of embryos misclassified as "normal" despite clinically significant mosaicism.

## Requirements

- R (≥4.0)
- R packages: `tidyverse`, `patchwork`, `scales`

Install packages with:
```r
install.packages(c("tidyverse", "patchwork", "scales"))
