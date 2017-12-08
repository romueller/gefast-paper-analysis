#!/bin/bash

# Runs the complete analysis for the paper: 
# downloading and preparing the data,
# running the tools, and visualising the results.



# ===== Preparations =====

bash scripts/prepare.sh


# ===== Clustering quality =====

# even data set, threshold 1 - 10
bash scripts/analysis_quality.sh even 1 10

# uneven data set, threshold 1 - 10
bash scripts/analysis_quality.sh uneven 1 10


# ===== Performance (runtime & memory) =====

# eldermet data set, threshold 1 - 10
bash scripts/analysis_performance.sh eldermet 1 10

# subsets of eldermet data set, threshold 1 - 2
bash scripts/analysis_subsampling.sh eldermet 1 2

