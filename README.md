# Borderline CNVs in PGT-A: a simulation study

## Overview
This repository contains reproducible simulation code for the paper:

Title: Borderline copy number variations in preimplantation genetic testing: when signal becomes noise
Authors: Anere Cye
Journal: Royal Society Open Science (2026)
DOI: 10.5281/zenodo.19421746

## Repository structure

- R/ - simulation functions
- scripts/ - run_all.R and scenario scripts  
- results/ - generated automatically

## How to run

Open R in the project root and run:
source(scripts/run_all.R)

All simulations use set.seed(2025). 5000 replicates per condition.

## Key results (Scenario B, n=5000)

False normal rate M=15%: 34.5% [33.2-35.9%]
False positive rate diploid: 54.4% [53.0-55.8%]
