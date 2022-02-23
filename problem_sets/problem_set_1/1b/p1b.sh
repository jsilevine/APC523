#!/usr/bin/env bash

## run script for different core numbers
for nthreads in 1 2 4 8 16
do
    ./a.out 128 1e-6 $nthreads
done

## generate timing plot
Rscript p1b_fig_gen.R

