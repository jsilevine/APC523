#!/usr/bin/env bash
echo "running poisson"

## run .cpp program and generate final_phi.csv
./poisson 128 1e-6

## run p1a.R to generate density plot
Rscript p1a.R
