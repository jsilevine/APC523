#!/usr/bin/env bash                                                                                                                                                                                                                                                             

## run script for different core numbers                                                                                                                                                                                                                                        
for h in 0.1 0.05 0.01 0.005 0.001
do
    ./a.out $h
done

## generate timing plot                                                                                                                                                                                                                                                         
Rscript p2b.R



