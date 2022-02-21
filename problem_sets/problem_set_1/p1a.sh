#!/usr/bin/env bash
echo "running poisson"

## run .cpp program and generate final_phi.csv
./poisson 128 1e-6

## run fig_gen.R to generate density plot
Rscript p1a_fig_gen.R
